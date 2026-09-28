#!/usr/bin/env python3
"""
Python module for reading/writing HDF5 files in KID format (fast readout).

:note: Module is inspired on the tesfdm.fdmlib.hdffile module.
"""
from   dataclasses      import dataclass
from   datetime         import datetime
import getpass
import logging
import socket
import time
import matplotlib.pyplot as plt

import h5py
import numpy            as     np

logger = logging.getLogger(__name__)

# Exported symbols.
__all__ = ( 'NoiseFile')


# function copied from the module kids.lib.comb

def binToFreq(bins, framelen, sampleRateMHz):
    """Convert bin number to frequency in MHz.

    :param bins:            bin number(s) as single integer or Numpy array
    :param framelen:        2-log of number of points per frame
    :param sampleRateMHz:   sample rate in MHz
    :return:                frequency(s) as integer or Numpy array

    :note: Negative frequencies may be specified as negative bin numbers or
           as bin numbers above 2^(framelen-1).
    """
    # pylint: disable=invalid-name
    k = 2**framelen
    b = (bins + k//2) % k - k//2
    return b * sampleRateMHz / float(k)



# ##########################################################
#
#  class HdfFile
#

class HdfFile:
    """An open HDF5 file in KID format (read-only). """

    def __init__(self, filename, mode='r'):
        """Construct a new HdfFile object referring to an existing HDF5 file on disk.

        :param str  filename: Full path name to the HDF5 file on disk
        """
        self.filename = filename
        # Open HDF5 file.
        self.hdf = h5py.File(filename, mode,
                             track_order=True # This is needed to preserve ordering of the datasets as in the file
                             )

        # Check file format version number.
        version = self.hdf.attrs.get('format_version', 'n/a')
        version = version.decode() if isinstance(version, bytes) else version

        # Retrieve DAC and ADC sample rates for later use in conversions
        channelAttrs = self.hdf['rack_00']['channel_00'].attrs
        self.dacSampleRateMhz = self._getBoardSampleRate(channelAttrs, 'DAC')
        self.adcSampleRateMhz = self._getBoardSampleRate(channelAttrs, 'ADC')

        # LO frequency could be unavailable, e.g. when using a test system with only an IF-loopback
        self.loFrequency      = channelAttrs.get('LO frequency', 0.0) # used to display pixels with their RF frequency 
        if not self.loFrequency:
            logger.warning('No LO frequency available, pixels are displayed with IF frequencies')

    @staticmethod
    def _getBoardSampleRate(channelAttrs: dict, boardType: str):
        """ Retrieves the board sample rate from the HDF5 file (attribute).

        In case the attribute isn't present in the file, it is derived from the board type. """ 
        sampleRate = channelAttrs.get(f'{boardType} sample rate MHz', None)

        if sampleRate is None:
            # This fall back is for when the sample rate is not in the HDF5 file yet (was added later)
            logger.debug("File doesn't contain %s board sample rate, defaulting to 2000 MHz", boardType)
            sampleRate = 2000.0
        return float(sampleRate)

    def __del__(self):
        """ Close file handle, if it was still open """
        if hasattr(self, 'hdf') and self.hdf:
            self.close()

    def close(self):
        """Close underlying HDF5 file. This object and its associated
        Dataset objects become invalid."""

        if self.hdf:
            self.hdf.close()
        self.hdf = None

    def getFileInfo(self):
        """Return dictionary of meta-data attributes at root level."""

        return {(key, value.decode() if isinstance(value, bytes) else value) \
                for key, value in self.hdf.attrs.items()}

    @property
    def startTime(self):
        """ Return datetime object with start time of the measurement. """
        return datetime.fromtimestamp(self.hdf.attrs['time_start_utc'] / 1000)

    def __iter__(self):
        """ Make this class iterable (not an iterator itself). """
        yield from self.hdf.__iter__()

    def countDisabledItems(self):
        """ Return number of disabled items in pulse or noise file """

        count = 0
        total = 0
        for pixelItem in self:
            for sample in pixelItem:
                total += 1
                if not sample.enabled:
                    count += 1
        return count, total


@dataclass
class PixelInfo:
    """ This contains the parameters of a pixel.
    By using this class, parameters no longer have to be passed on individually. """
    freq: float
    dacbin: int
    pixelNr: int
    centerI: float
    centerQ: float
    baselineAngle: float


class NoiseSample:
    """ Single virtual noise sample (part of a longer time-series) """

    def __init__(self, pixelInfo: PixelInfo, dataset, start, stop):
        self.pixelInfo = pixelInfo
        self.dacbin = pixelInfo.dacbin
        self.freq   = pixelInfo.freq
        self._start = start
        self._stop = stop
        self.enabled = True # Sample can be excluded, e.g. when it contains a cosmic ray

        # Pre-process data (from HDF5 format to real KID phases)
        # Scale digital 16 bit values
        # 3 extra is needed for the moving-average filter (as done in FW for pulse-data)
        # dataIs = dataset[start:stop + 3,0].copy() * (1.0 / 2**15)
        # dataQs = dataset[start:stop + 3,1].copy() * (1.0 / 2**15)

        dataIs = dataset[start:stop,0].copy() * (1.0 / 2**15)
        dataQs = dataset[start:stop,1].copy() * (1.0 / 2**15)

        # Apply KID center (no -= here to avoid unmutable data issues)
        dataI = dataIs - pixelInfo.centerI
        dataQ = dataQs - pixelInfo.centerQ

        # Now calculate phase
        phase = np.arctan2(dataQ, dataI) # Note: Q is the intended first parameter of the arctan function

        # Unwrap (raw phase can jitter around +/- Pi)
        phase = np.unwrap(phase)

        # Correct for baseline angle
        self._phase = phase - pixelInfo.baselineAngle

        test = phase - pixelInfo.baselineAngle

        # if np.max(test > np.pi):
        #     plt.figure()
        #     plt.plot(test), plt.title(f'dacbin {self.dacbin}')
        #     plt.plot(phase)
        #     plt.show()

        # Now perform moving average (also done in FW for pulse-data, before pulse-detection)
        # The only difference is that the FW still has integers at this point, and we are already using floating point
        # numbers (needed for the arctan2 for example)
        # self._phase = self._moving_average(self._phase, 4) # This loses the 3 extra samples taken above

    @staticmethod
    def _moving_average(x, w):
        # below is a quick floating-point moving average
        # taken from https://stackoverflow.com/questions/14313510/how-to-calculate-rolling-moving-average-using-python-numpy-scipy#54628145
        return np.convolve(x, np.ones(w), 'valid') / w

    def __str__(self):
        return f"NoiseSample dacbin={self.dacbin}, data[{self._start}:{self._stop}] (shape={self._phase.shape})"

    @property
    def phase(self):
        """ return the KID phase with respect to the baseline. """
        return self._phase.copy()



class PixelNoise:
    """ All virtual noise samples of a pixel. """

    def __init__(self, pixelInfo: PixelInfo, dataset, windowLength):
        self.dacbin       = pixelInfo.dacbin
        self.pixelNr      = pixelInfo.pixelNr
        self.freq         = pixelInfo.freq
        self.dataset      = dataset
        self.windowLength = windowLength
        assert self.dataset.shape[1] == 2, f"Expecting IQ data, shape was {self.dataset.shape}"

        self.samples : NoiseSample = []

        # construct virtual samples from the noise time-series
        start = 0
        stop  = windowLength

        while stop <= len(self.dataset):
            # Don't slice the data itself immediately to reduce the number of operations during the data-review phase
            # Only slice upon use
            self.samples.append(NoiseSample(pixelInfo, self.dataset, start, stop))
            start += self.windowLength
            stop  += self.windowLength

    def __iter__(self) -> NoiseSample:
        """ Make this class iterable (not an iterator itself). """
        yield from self.samples.__iter__()

    def __len__(self) -> int:
        """ Gives the number of enabled samples. Is used to check whether at least one sample is enabled. """
        return sum(sample.enabled for sample in self.samples)

    def __getitem__(self, sampleId) -> NoiseSample:
        return self.samples[sampleId]

    def __str__(self) -> str:
        return f"Pixel noise wrapper (dacbin={self.dacbin}, {len(self)} samples)"

    @property
    def allPhaseData(self) -> np.ndarray:
        """ Returns a contiguous array with the noise phase-data for processing. """
        return np.concatenate([sample.phase.copy() for sample in self if sample.enabled])




class NoiseFile(HdfFile):
    """An open HDF5 file in KID format containing noise data (read-only). """

    #: Noise sample length, must be equivalent to the pulse window length
    DEFAULT_WINDOW_LENGTH = 1024

    def __init__(self, filename, windowLength=DEFAULT_WINDOW_LENGTH):
        """Construct a new object referring to an existing HDF5 file on disk containing noise-data.

        :param str filename:     Full path name to the HDF5 file on disk.
        :param int windowLength: Noise sample length, must be equivalent to the pulse window length.
        """
        super().__init__(filename)

        self.pixels = {}

        # Check that it's really a file with noise
        # TODO: Note, this assumes 1 channel, perhaps use h5py's .visit() function, or iterate over the H5Groups
        channelData    = self.hdf['rack_00']['channel_00']
        self.numPixelGroups = 0
        self.numPixels      = 0

        # Each 8-pixel group has 2 datasets, 1 containing the pixel configuration, 1 with the noise data
        for dataSetName, dataSet in channelData.items():
            if not dataSetName.startswith('continuous_group'):
                continue
            # Fetch the configuration dataset also
            pixelGroup = int(dataSetName[16:])
            try:
                configSet = channelData[f'configuration_group{pixelGroup}']
            except KeyError as ex:
                raise ValueError(
                    f"Noise file contents are inconsistent (missing pixel configuration data): {ex}") from None

            # pixel numbers were added to the continuous files later, so check whether they're present
            # This is the number that is used by the ADC board firmware to identify the pixel
            if 'pixel' in configSet.dtype.fields.keys():
                pixelNumbers = configSet['pixel']
            else:
                pixelNumbers = [None] * len(configSet)

            for indexInGroup, (dacbin, pixelNr, centerI, centerQ, baselineAngle) in enumerate(
                zip(configSet['dacbin'], pixelNumbers,
                    configSet['center_i'], configSet['center_q'], configSet['baseline_angle'])):
                freq = binToFreq(dacbin, 19, self.dacSampleRateMhz) + self.loFrequency
                pixelInfo = PixelInfo(freq, dacbin, pixelNr, centerI, centerQ, baselineAngle)
                self.pixels[dacbin] = PixelNoise(pixelInfo, dataSet[:,indexInGroup*2:indexInGroup*2 + 2],
                                                 windowLength)
                print(f"\rloading noise data for pixel {self.numPixels + indexInGroup}", end="")

            self.numPixelGroups += 1
            self.numPixels = len(self.pixels)
        print("\r--- done loading noise data ---          ")

        logger.info("Found %s pixel groups with %s pixels in total", self.numPixelGroups, self.numPixels)

    @property
    def allPixels(self):
        """ Get list of DAC bins. """
        return list(self.pixels.keys())

    def __len__(self):
        """ Returns the number of pixels for which there is noise data in this file. """
        return len(self.pixels)

    def __getitem__(self, dacbin) -> PixelNoise:
        """ Acces pixel noise-data. """
        assert self.hdf
        return self.pixels[dacbin]

    def __iter__(self) -> PixelNoise:
        """ Make this class iterable (not an iterator itself). """
        yield from self.pixels.values().__iter__()



if __name__=='__main__':
    def main():
        """ Test function for the module. """
        # Don't pollute global module
        from argparse import ArgumentParser # pylint: disable=import-outside-toplevel
        from pprint   import pprint         # pylint: disable=import-outside-toplevel
        parser = ArgumentParser(description=__doc__)
        parser.add_argument('--file', type=str, default='noise.h5',
                            help='HDF5 noise file for testing (default: %(default)s)')
        args = parser.parse_args()
        # No logging config is loaded from the Egress configuration file, so apply one here
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s %(module)25s:%(lineno)-4d : %(levelname)-7s: %(message)s')

        nFile = NoiseFile(args.file)
        print(f"File was recorded on: {nFile.startTime}")
        print("File metadata:")
        pprint(nFile.getFileInfo())
        print()
        print(f"{len(nFile.allPixels)} pixels (indexed per DAC bin):")
        pprint(nFile.allPixels)
        print()
        print("First pixel time series")
        pixelNoise = nFile[nFile.allPixels[0]].allPhaseData
        print(f"{pixelNoise} ({pixelNoise.shape} samples)")

    main()

