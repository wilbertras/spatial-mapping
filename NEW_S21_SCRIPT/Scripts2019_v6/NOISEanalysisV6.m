function [NOISE,NOISEParameters] = NOISEanalysisV6
% This function will search in the specified directory with noise data for
% all (KID,Pread) combinations and read the S21(T)/FFT(T)/TimeDomain(T)
% sweeps. It will then process all S21(F) sweeps to obtain Fres,Q and other
% KID specific parameters, as well as calibrate the resonance circle. It
% will then combine different powers (at base [lowest] temperature) to
% determine the optimum power to operate KID and compares the noise to the
% Jiansong Gao line and slope. Finally all information is stored to figure
% and csv files. Only base temperature power dependence is shown in figures.
%
% 
%
% NOTE: This script does not do the NEP calculation.
%
% NOTE: 2D (power and temperature) or 1D power data is handled by treating
%   each (KID,power) combination as a new KID.
%
% NOTE: Function assumes External subroutines are located in the
%   ..\subroutines subdirectory, while this function is in the main
%   directory.
%
% NOTE: Before processing checks the size of the FFT noise files. If a file
% is regularly large or small it is disregarded. This could have
% significant effects if you are processing different (KID,P) combinations
% at once that each have a different number of temperatures measured.
%
% NOTE: in time domain files it is assumed that a single timestep dt is
% used for the sweeps at all temperatures. Hence dt of each (KID,P)
% combination is only obtained from the first subheader.
%
% replaced for most purposes by 2D version, except if you have 2D data and
% want to read in only the first T then use this one.
%
%Reads in all noise (*FFT.dat) files in specified dir. and analyzes them.
%NOT WORKING WITH 2D DATA
%NB: will check the first file read in to get # powers and # temperatures. If any of the two gets larger in other files it is ignored
%this makes the program resistant against bla fuckup from setup!
%
% Current VERSION: 5.2
%   Based upon the original "NOISEANALYSIS_HybridsV2" program
%   V1.0 original script by Stephen Yates & Jochem Baselmans
%   V1.1 (2010-09-10,JB/SY):  improved fitting of fnoise data -used to
%                             throw away one point for fit, now uses all
%                             data. Also sorts data in power, not Pint
%                             (which could go wrong with overdriven KIDS)
%   V1.2 (2012-03-03,JB/SY): Added plotlimrange option
%   V2.0(2012-10-24,RJ): Full revamp of the code
%   V2.1(2012-12-11,RJ): Based on TUD converted data: Build in checksize
%   switch, and MeanFreqNoise interpolation was given 1 as default value
%   for extrapolation. (TUD data measures from 5 Hz and up).
%   V5.0 (2012-12-xx,JB): Fits now also the frequency noise vs P at
%   1,10 and 100 Hz (in addition to the 1 kHz).
%   V5.1 (2013-02-14,RJ): Now uses FitS21main5.
%   V5.2 (2014-02-14,JB): Now uses FitS21main5.
%   V6 (2014-02-14 JB): Allows to be ran again with Popt.csv used for
%   finding Popt. this allows you to get the SF fits ok. Should be done as
%   standard procedure. Also contains phase and amplitude noise wrt complex
%   plane (will be plotted in NoisecompareV8 and higher)
%
%INPUT:
%   SPECIFIED INTERNALLY
%
% FFTsubdir == subdirectory specifying the location of the S21 files.
%               Starts with "\", ends without "\"
%
% ChipInfo == a struct containing information about the measured chip in
%           the following elements
%   ChipInfo.path == Root path containing the measurement data.
%   ChipInfo.INdexPref == Index of the power that is used as a reference
%                         for Popt finding. 1=lowest power, 2= second
%                         lowest power, etc. If the index is higher than
%                         the number of powers, the highest power is taken.
%
% EXTRA: BWrange, FITF0, T0only (See below).
%==========================================================================
% Set ChipInfo inside function.
%==========================================================================
%sub-dir inside ChipPath where the Noise Data [FFT(T),S21(T),ft(T)] is located. Starts with \ ends without \.
ChipInfo.path = [cd '/..']; %root path where data is, one higher than the scripts
FFTsubdir='/Noise100mK_CPW/FFT/Power';
%Note: pwd is current matlab directory. Usually directory where this m-file
%is located.
    ChipInfo.IndexPref = 1;                                     %Index of the power that is used as reference for Popt finding.
%==========================================================================
% Internal variables which may be user modified.
%==========================================================================
global FITF0
ReadPoptfile=0; %default=0. If set to 1 optimum power takemn from Popt file, that is created by this script. Usefull to be abel to rerun the sript for the correct Pot.
BWrange=1;      %Width aroung Fres of the reference power that is considered to find Popt. Default 1
%small enough an freq shift to be allowed to give Popt based %upon noise level requirement.
FITF0=3;          %Switch used to determine which method is used for determination of F0.
%Please see the FitS21main5 routine for details. FITF0 = 0 is recommended due to the presence of overdriven
%resonators in Power sweeps (it will only update Q, not F0. F0 mot very important here).
T0only=1;         %Switch, if true (1) all figures except base temperature are closed.
%For all temperatures a figure will be made and saved, but keeping them
%open can crash your pc.
checksize = 1;    %Switch, if true (1) the FFT files are checked in advance
%to see if their size (in bits) is not deviating from the others. 
MaxnT = 30 ;    %def = 30; max nuymber of temperatures looked at. Must always be set (set to very large if you want to get all for sure)
                % prohram only reads in the firsn MaxnT temperatures
%==========================================================================
%==========================================================================
%OUTPUT:
%   NOISE(N). == struct array of N = # KIDs (actually # (KID,P) combinations),
%                which for all temperatures (nT = # temperatures)
%NB:any loop over the length of NOISE has index p and length(NOISE)
%NB: anyt loop over KIDnumber has index kidn
%       NOISE(N).KIDnumber = KIDnumber (ID)
%       NOISE(N).filename = a single string containing the full
%                           path+filename that contains the S21datafile.
%       NOISE(N).Temperature(nT,1) = vector containing the chip temperature [K]
%       NOISE(N).ReadPower(nT,1) = vector containing the used read power [dBm]
%       NOISE(n).S21Parameters(nT,1) == [Fres,S21min,Ql,Qi,Qc] as
%                       determined by labview and stored in datafile header.
%       NOISE(N).S21_MPplane{nT,1}(:,1:5) = cell array containing the S21data at each Temperature
%              Each cell contains a double array with the columns:
%              [F (GHz), |S21|, phase (rad)]
%       NOISE(N).S21_IQplane{nT,1}(:,1:5) = cell array containing the S21data at each Temperature
%              Each cell contains a double array with the columns:
%              [F (GHz), I, Q]
%       NOISE(n).Ql(nT,1) == Measured Loaded Q at each temperature [From S21 fit]
%       NOISE(n).Qi(nT,1) == Measured Internal Q at each temperature [From S21 fit]
%       NOISE(n).Qc(nT,1) == Measured Coupling (External) Q at each temperature [From S21 fit]
%       NOISE(n).Fres(nT,1) == Measured Resonance Frequency [From S21 fit, subsampling resolution]
%       NOISE(n).S21min(nT,1) == Minimum value S21 in dB [From S21 fit, subsampling resolution]
%       NOISE(n).S21Fit{nT,1} == Fit to the S21 (MPplane) evaluated. [Freq,|S21|]
%       NOISE(n).InternalPower(nT,1) == KID internal power [dBm]
%       NOISE(n).Fread(nT,1) == Frequency used to read the noise at.
%       NOISE(n).FFTnoise{nT,1} == cell array containing the measured noise at each temperature.
%                   Each cell contains a double array with the columns
%                   [1 F(Hz),2 Phase Noise dBc/Hz,3 Amp Noise dBc/Hz, 4 SF/F^2 (Gao, NOT log scale), 5 Freq Noise (Sx) (Mazin) NOT log scale]
%                   and col 6 col 7 is phase noise and ampl noise the noise
%                   resp. wrt the complex plane dBc/Hz (Set in NoisecompareV8)
%       NOISE(n).MeanNoise(nT,3) == array containing for each temperature
%                   [S_A(100Hz),S_A(10Hz),S_ph(10Hz),S_A(100Hz),S_ph(100Hz)]
%                   1) the mean amplitude noise between 20 and 500 Hz (required
%                   for Popt determination)
%                   2) the mean amplitude noise at 10 Hz
%                   3) the mean phase noise at 10 Hz
%                   4) the mean amplitude noise at 100 Hz
%                   5) the mean phase noise at 100 Hz
%       NOISE(n).MeanFreqNoise(nT,4) == array containing for each
%                   temperature the mean frequency noise (Gao definition)
%                   at: [1 Hz, 10 Hz, 100 Hz, 1 kHz]
%       NOISE(n).OptimumPower(nT*,4) == An array that is essentially a
%                   look-up table. In this table (which is KID precific,
%                   rather that (KID,P) combination specific) there are 4
%                   columns(dim2): [Temperature,Tindex,Optimum Power,Popt Index].
%                   The array has nT* rows (dim1) where nT* is not
%                   neccessarily equal to nT and the values in Temperature
%                   are not neccessarily equal to NOISE(n).Temperature.
%                   The look-up table works as follows:
%                   For the KID NOISE(q).KIDnumber at temperature
%                   NOISE(q).OptimumPower(r,1) in column 1, the optimum 
%                   power is NOISE(q).OptimumPower(r,3). The data (for
%                   example the resonance frequency)accompanying this 
%                   temperature and power can be found in
%                   NOISE(NOISE(q).OptimumPower(r,4)).Fres(NOISE(q).OptimumPower(r,2),1)
%
%
%   NOISEParameters == Nx16 array contains the most important information
%                  about all N=#(KID,P) combinations at base temperature
%                  It contains the parameters (Dim 2):
%               1 = KIDID
%               2 = Tbase [K]
%               3 = Fres(Tbase) [GHz]
%               4 =
%               5 = Q
%               6 = Qi
%               7 = Qc
%               8 = Popt [dBm]
%               9 = Pread [dBm]
%               10 = Pint [dBm]
%               11 = mean freq noise (SF/F^2 in dBc/Hz) at 1 Hz
%               12 = mean freq noise (SF/F^2 in dBc/Hz) at 10 Hz
%               13 = mean freq noise (SF/F^2 in dBc/Hz) at 100 Hz
%               14 = mean freq noise (SF/F^2 in dBc/Hz) at 1 kHz
%
%   In addition all results are output as 3 types of files (Described below):
%   In the FFTsubdir where all measurements are:
%   One overview file with information of all KIDs at the base temperature.
%   One file per KID containing S21(F) at base temperature
%   One file per KID containing analysis results at all measurement temperatures.
%   In the measurement main directory:
%   A summary file containing xx columns with the values of NOISEParameters.
%
%
%SUBROUTINES:
% import_data (External)
% FitS21main3 (External)
% WriteSRONcsv (External)
% GenerateColorsFromMap (External)
% 	colormapstorage.mat (External)
%
%
%DATE: October 24, 2012
%AUTHOR: Reinier Janssen, Jochem Baselmans
%==========================================================================

%==========================================================================
% Setting some routine default values.
%==========================================================================
global Noisepath
format('long','g'); %Set display format of numbers to 7 digits
%Enable subroutines by adding path in search path.
addpath([pwd,filesep,'subroutines']);
close all; %close all open plots to remove clutter.
%==========================================================================
% In the datapath find all files and identify KIDs.
%==========================================================================
Noisepath = [ChipInfo.path,FFTsubdir,filesep]; %Path containing the raw noise data
if ReadPoptfile==1
    PoptFile = [Noisepath,'Popt.csv'];
    fid = fopen(PoptFile);
    if fid == -1
        disp([Noisepath,'Popt.csv'])
        error('Error: NOISEAnalysis: Cannot find Popt.csv, set ReadPoptfile to 0.')
        
    else
        fclose(fid);
        %There is a Popt file. Read it and use it.
        fprintf('Found Popt.csv, this will be used to indicate Popt.\n')
        [~,PoptData] = ReadSRONcsvV2(PoptFile,'',0);
        %PoptData will be Nx3. Where dim 2 contains [T,KIDid,|Popt [dBm]|]
    end
end
%Search the noisepath for all FFT files and filter out all files that are
%significantly longer or shorter than the mean. (This usually means
%something went wrong during measurements)
RawFFTfiles = dir([Noisepath,'KID*FFT*.dat']);

if checksize
    datasize = zeros(1,length(RawFFTfiles));
    for p=1:length(RawFFTfiles)
        datasize(p)=RawFFTfiles(p).bytes;
    end
    FFT2read = 0.5*mean(datasize) < datasize & datasize < 1.5*mean(datasize);
    RawFFTfiles = RawFFTfiles(FFT2read);
    clear datasize
end

%Determine the KID numbers (IDs) from the names of the good files.
KIDnumbers = zeros(1,length(RawFFTfiles));
for p=1:length(RawFFTfiles) %Loop over all FFT files
    %Determine the KIDid (its number, XX) and power from each filename.
    KIDnumbers(p) = cell2mat(textscan(RawFFTfiles(p).name,'%*3s %f %*s'));
    NOISE(p).KIDnumber = KIDnumbers(p);
end
%Determine all unique KIDs.
KIDnumbers = unique(KIDnumbers);

%Print to screen some of the reading information.
fprintf('Search for FFT data performed in:\n')
disp(Noisepath)
fprintf(['Inside the path the following KIDs are available: ',num2str(KIDnumbers),'\n'])
fprintf(['Inside this path a total number of ',num2str(length(RawFFTfiles)),' files were found.\n'])

%==========================================================================
% Read in all the data files
%==========================================================================
%Predefine some variables
NOISEParameters = zeros(length(RawFFTfiles),10); %A quick access reference array.

for p=1:length(RawFFTfiles) %LOOP OVER ALL FILES (aka KID-P-combinations)
    %Each (KID,P) combination has 4 different files in which all
    %temperatures are concatenated. These 4 files all have the same name
    %with the exception of the text just before the .dat extension. Hence a
    %common section of the name can be constructed from the RawFFTfile
    %name.
    
    %Construct the root full path name.
    LocEndRoot = strfind(RawFFTfiles(p).name,'FFT.dat');
    if isempty(LocEndRoot)
        fprintf('ERROR Noise Analysis: Cannot find FFT.dat in the name of the noise file.\n')
        fprintf('Most likely text has been placed between FFT and .dat extension.\n')
        error('Cannot create RootName for file detection.\n')
    end
    RootName = [Noisepath,RawFFTfiles(p).name(1:LocEndRoot-1)];
    
    %Store the RootName
    NOISE(p).filename = RootName;
    
    %======================================================================
    %One by one read in the data files
    %======================================================================
    
    %Reading FFT file
    FFTfile = [RootName,'FFT.dat'];
    [Data,Temperature,Power,FFTheader] = import_data(FFTfile);
    NOISE(p).FFTnoise = Data;
    if length(Temperature) > MaxnT
        Temperature(MaxnT+1:end)=[];
    else
        MaxnT = length(Temperature);
    end
    NOISE(p).Temperature = Temperature;
    NOISE(p).ReadPower = -1*Power*ones(length(Temperature),1);
    
    %Find in FFTheader the frequency of the tone used to read the noise.
    for hl = 1:length(FFTheader)
        Freadindex = cell2mat(strfind(FFTheader(hl,1),'F used')); %Try to find F used in this line
        ColonIndex = cell2mat(strfind(FFTheader(hl,1),':')); %Check for : in this line
        if isempty(Freadindex)
            %Not in this line
        else
            %If it is in the line, find the first : after F used.
            Freadstart = find(ColonIndex>Freadindex,1);
            Freadstart = ColonIndex(Freadstart);
            Fread = cell2mat(textscan(FFTheader{hl,1}(Freadstart+1:end),'%f'));
        end
    end
    NOISE(p).Fread = Fread*ones(length(Temperature),1);
        
    %======================================================================
    %Reading S21 file (GHz,Re,Im)
    S21file = [RootName,'S21.dat'];
    [Data,Temperature,Power,~] = import_data(S21file);
    %Check for equality of temperature and power
    if abs(-1*Power-NOISE(p).ReadPower(1)) > 0.5
        disp('Warning: Read Power does not match between FFT and S21 file.')
    end
    if length(Temperature) < length(NOISE(p).Temperature)
        (Temperature) 
        (NOISE(p).Temperature)
        error('ERROR: Number of Temperatures in S21 file too small wrt FFT file')
    end
    Tcheck = sum(abs(Temperature(1:MaxnT)-NOISE(p).Temperature) > 0.001*NOISE(p).Temperature);
    if Tcheck > 0
        disp('Warning: ',num2str(Tcheck),' Temperatures does not match between FFT and S21 Re Im file.')
    end
    %Copy data for struct
    NOISE(p).S21_IQplane = Data;
    
    %======================================================================
    %Reading S21 file (GHz,dB,rad)
    S21DBfile = [RootName,'S21dB.dat'];
    [Data,Temperature,Power,S21dBheader] = import_data(S21DBfile);
    %Check for equality of temperature and power
    if abs(-1*Power-NOISE(p).ReadPower(1)) > 0.5
        disp('Warning: Read Power does not match between FFT and S21dB file.')
    end
    if length(Temperature) < length(NOISE(p).Temperature)
        (Temperature) 
        (NOISE(p).Temperature)
        error('ERROR: Number of Temperatures in S21dB file too small wrt FFT file')
    end
    Tcheck = sum(abs(Temperature(1:MaxnT)-NOISE(p).Temperature) > 0.001*NOISE(p).Temperature);
    if Tcheck > 0
        disp('Warning: ',num2str(Tcheck),' Temperatures do not match between FFT and S21dB file.')
    end
    %Copy data for struct
    NOISE(p).S21_MPplane = Data;
    
    %======================================================================
    %Reading time-domain file
    TDfile = [RootName,'td.dat'];
    [Data,Temperature,Power,TDHeader] = import_data(TDfile);
    %Check for equality of temperature and power
    if abs(-1*Power-NOISE(p).ReadPower(1)) > 0.5
        disp('Warning: Read Power does not match between FFT and TD file.')
    end
    if length(Temperature) < length(NOISE(p).Temperature)
        error('ERROR: Number of Temperatures in TD file too small')
    end
    Tcheck = sum(abs(Temperature(1:MaxnT)-NOISE(p).Temperature) > 0.001*NOISE(p).Temperature);
    if Tcheck > 0
        disp('Warning: ',num2str(Tcheck),' Temperatures does not match between FFT and TD file.')
    end
    
    %Extract the timestep from the first subheader.
    fid=fopen(TDfile);
    dt = cell2mat(textscan(fid,'%*s%*s%*s%f',1,'headerlines',length(TDHeader)+1));
    fclose(fid);
    
    %Copy data for struct
    for nT=1:MaxnT
        CompleteData = zeros(size(Data{nT},1),size(Data{nT},2)+1);
        CompleteData(:,1) = dt*(1:size(Data{nT},1))';
        CompleteData(:,2:end) = Data{nT}(:,:);
        NOISE(p).TDIQ{nT,1} = CompleteData;
    end
    clear Temperature
    %======================================================================
    %Analyse the S21 data measured in [mag,phase] space for each temperature
    %======================================================================
    %Initialize some temporary storage variables.
    % limit amount of temperatures

    Ntemperatures = length(NOISE(p).Temperature);
    fres = zeros(Ntemperatures,1); % resonance frequency
    Q = zeros(Ntemperatures,3); %[Ql,Qi,Qc]
    S21min = zeros(Ntemperatures,1); %|S21(Fres)| (dB)
    FitResult = cell(Ntemperatures,1); %[F(GHz),S21(dB)] Evaluation of the fit
    
    for nT=1:length(NOISE(p).Temperature)
        %Normalize the S21 data measured in magnitude plane
        filtered=smooth(NOISE(p).S21_MPplane{nT}(:,2),3); %smoothing the |S21| data in dB space
        %normalise S21 in log space to the max(|S21|)
        NOISE(p).S21_MPplane{nT}(:,2)=NOISE(p).S21_MPplane{nT}(:,2)-max(filtered);
        %Convert dB to magnitude
        NOISE(p).S21_MPplane{nT}(:,2) = 10.^(NOISE(p).S21_MPplane{nT}(:,2)/20);
        
        %Perform fit to obtain resonator parameters. Note the FITF0 is
        %recommended for overdriven resonators.
        [fres(nT),Q(nT,:),S21min(nT),FitResult{nT}] = FitS21main5(NOISE(p).S21_MPplane{nT}(:,1:3),FITF0);
    end
    %Put the temporary storage variables into the NOISE struct
    NOISE(p).Ql=Q(:,1);
    NOISE(p).Qi=Q(:,2);
    NOISE(p).Qc=Q(:,3);
    NOISE(p).Fres=fres;
    NOISE(p).S21min=S21min;%in dB!
    NOISE(p).S21fit = FitResult;
    %Calculate the Internal Power
    NOISE(p).InternalPower = 10*log10((2/pi)*10.^(NOISE(p).ReadPower/10).*(NOISE(p).Ql.^2./NOISE(p).Qc));
    
    %Calculate the resonator ring time
    NOISE(p).TauRes = NOISE(p).Ql./(pi*NOISE(p).Fres);
    %Calculate the resonator bandwidth (used later in Popt determination)
    NOISE(p).Bandwidth = NOISE(p).Fres./NOISE(p).Ql;
    
    
    NOISE(p).MeanNoise = zeros(Ntemperatures,3);
    NOISE(p).MeanFreqNoise = zeros(Ntemperatures,4);
    %Calculate frequency noise
    for nT=1:Ntemperatures
        %Normalized Frequency Noise (Sf/F^2) [1/Hz] (defined by J. Gao)
        NOISE(p).FFTnoise{nT}(:,4) = (10.^(NOISE(p).FFTnoise{nT}(:,2)/10))*(1/(4*NOISE(p).Ql(nT,1))^2);
        %Frequency Noise (defined by B. Mazin)
        NOISE(p).FFTnoise{nT}(:,5) = NOISE(p).FFTnoise{nT}(:,4)*(NOISE(p).Fres(nT,1)*1e9)^2;
        %Phase noise WRT complex plane
        S21min_a=10^(( NOISE(p).S21min(nT) )/20);%from dB to magnitude, not stored
        R=20*log10((1-S21min_a)/2);%correction due to circl radius vs complex plane radius in dB
        NOISE(p).FFTnoise{nT}(:,6)=NOISE(p).FFTnoise{nT}(:,2)+R;
        %Amplitude noise WRT complex plane
        NOISE(p).FFTnoise{nT}(:,7)=NOISE(p).FFTnoise{nT}(:,3)+R;
        
        %Calculate Mean amplitude noise between 20 and 500 Hz (Required for
        %Popt determination)
        Frange = 20 <= NOISE(p).FFTnoise{nT,1}(:,1) & NOISE(p).FFTnoise{nT,1}(:,1) <= 500;
        NOISE(p).MeanNoise(nT,1) = mean(NOISE(p).FFTnoise{nT,1}(Frange,3));
        
        %Calculate Mean amplitude and phase noise around 10 Hz
        NOISE(p).MeanNoise(nT,2) = ... % Amp @ 10 Hz
            mean(interp1(NOISE(p).FFTnoise{nT,1}(:,1),NOISE(p).FFTnoise{nT,1}(:,3),(9:0.1:11),'linear',1));
        NOISE(p).MeanNoise(nT,3) = ... % Phase @ 10 Hz
            mean(interp1(NOISE(p).FFTnoise{nT,1}(:,1),NOISE(p).FFTnoise{nT,1}(:,2),(9:0.1:11),'linear',1));
        
        %Calculate Mean amplitude and phase noise around 100 Hz
        NOISE(p).MeanNoise(nT,4) = ... % Amp @ 100 Hz
            mean(interp1(NOISE(p).FFTnoise{nT,1}(:,1),NOISE(p).FFTnoise{nT,1}(:,3),(90:1:110),'linear',1));
        NOISE(p).MeanNoise(nT,5) = ... % Phase @ 100 Hz
            mean(interp1(NOISE(p).FFTnoise{nT,1}(:,1),NOISE(p).FFTnoise{nT,1}(:,2),(90:1:110),'linear',1));
        
        %Calculate the mean frequency noise (Gao definition) at:
        NOISE(p).MeanFreqNoise(nT,1) = ... % 1 Hz
            mean(interp1(NOISE(p).FFTnoise{nT,1}(:,1),NOISE(p).FFTnoise{nT,1}(:,4),(0.8:0.1:1.2),'linear',1));
        NOISE(p).MeanFreqNoise(nT,2) = ... % 10 Hz
            mean(interp1(NOISE(p).FFTnoise{nT,1}(:,1),NOISE(p).FFTnoise{nT,1}(:,4),(9:0.1:11),'linear',1));
        NOISE(p).MeanFreqNoise(nT,3) = ... % 100 Hz
            mean(interp1(NOISE(p).FFTnoise{nT,1}(:,1),NOISE(p).FFTnoise{nT,1}(:,4),(90:1:110),'linear',1));
        NOISE(p).MeanFreqNoise(nT,4) = ... % 1 kHz
            mean(interp1(NOISE(p).FFTnoise{nT,1}(:,1),NOISE(p).FFTnoise{nT,1}(:,4),(900:10:1100),'linear',1));
    end
    
    %FILL THE NOISEParameters ARRAY WITH BASE TEMPERATURE INFORMATION
    [~,BaseTindex] = min(NOISE(p).Temperature);
    
    NOISEParameters(p,1) = NOISE(p).KIDnumber; %KID ID
    NOISEParameters(p,2) = NOISE(p).Temperature(BaseTindex); %Tbase
    NOISEParameters(p,3) = NOISE(p).Fres(BaseTindex); %Fres
    NOISEParameters(p,5) = NOISE(p).Ql(BaseTindex); %Q
    NOISEParameters(p,6) = NOISE(p).Qi(BaseTindex); %Qi
    NOISEParameters(p,7) = NOISE(p).Qc(BaseTindex); %Qc
    NOISEParameters(p,9) = NOISE(p).ReadPower(BaseTindex); %Pread
    NOISEParameters(p,10) = NOISE(p).InternalPower(BaseTindex); %Pint
    NOISEParameters(p,11:14) = NOISE(p).MeanFreqNoise(BaseTindex,:); %Frequency Noise
    
end %END OF LOOP OVER ALL FILES (aka KID-Pbb-combinations)

%==========================================================================
% Now that all files have been read we switch to looking per (KID,T)
% combination for all powers and determine Popt.
%==========================================================================
IndexPsort = cell(length(KIDnumbers),2);
%initialize cell array to store for each KID an array that contains:
% 1 = contains vector of length NP. This vector holds the indexes in NOISE
% that contain information about this resonator. The indexes are given in
% order of increasing power. Hence for this resonator there are NP
% different powers.
% BV: NOISE(IndexPsort{kidn,1}(nP,1)) indexes power 'nP' of KID 'kidn'
% 2 = Stores for each KID the PTindex array (defined at line 411).
LookUpTables = cell(length(KIDnumbers),1);
%Initialize a cell array to store for each KID 2 lookup table array that
%contains in each cell (dim2=1) a double array with
%   [T,Tindex(in NOISE struct elements),Popt,PoptIndex(in NOISE)]

for kidn=1:length(KIDnumbers) % LOOP OVER ALL UNIQUE KIDS
    %======================================================================
    % For each KIDs sort their stuff with increasing power
    %======================================================================
    KIDlocations = find(NOISEParameters(:,1) == KIDnumbers(kidn));
    [~,IndexPsort{kidn,1}]=sort(NOISEParameters(KIDlocations,9)); %Sort on Pread and save the indexes.
    IndexPsort{kidn,1}(:,1) = KIDlocations(IndexPsort{kidn,1}(:)); %Recalibrate the indexes to the full KID(Parameters) variable
    
    %Check if the desired reference index is within the number of powers
    %measured for this KID
    if ChipInfo.IndexPref < 1
        disp('WARNING NoiseAnalysis: IndexPref must be >0. Assuming minimum Power to be reference.')
        ChipInfo.IndexPref = 1;
    elseif ChipInfo.IndexPref > length(IndexPsort{kidn,1}(:,1))
        disp('WARNING NoiseAnalysis: IndexPref exceeds number of Powers. Assuming maximum Power to be reference.')
        ChipInfo.IndexPref = length(IndexPsort{kidn,1}(:,1));
    end
    %Shorthand variable for the reference power index
    IndexPref = IndexPsort{kidn}(ChipInfo.IndexPref,1);
    %Find the number of temperatures at the reference power.
    Ntemperatures = length(NOISE(IndexPref).Temperature);
    
    %Put some user info to screen.
    disp(['Searchin optimum power for KID ',num2str(KIDnumbers(kidn)),...
        ' in Power range: ',num2str(NOISE(IndexPsort{kidn,1}(1,1)).ReadPower(1)),...
        ' - ',num2str(NOISE(IndexPsort{kidn,1}(end,1)).ReadPower(1)),' dBm'])
    
    %======================================================================
    % Determine for each (KID,T) combination the indexes of all (P,T)
    % combinations matching closest in T.
    %======================================================================
    PTindex = zeros(length(IndexPsort{kidn,1}(:,1)),Ntemperatures,2,2);
    %For KID under consideration note for each (P,T) combination the index
    %in NOISE struct (dim 3=1) and temperature index (dim3 = 2), which gives
    %the minimal difference between requested Power and Temperature. If
    %dim4 = 1 indexes are given. If dim4= 2 actual values are given of P
    %and T.
    
    PT_AmpNoise = zeros(length(IndexPsort{kidn,1}(:,1)),Ntemperatures);
    PT_Fres = zeros(length(IndexPsort{kidn,1}(:,1)),Ntemperatures);
    %Two Array storing for each Power,Temperature either the mean amplitude
    %noise between 20 Hz and 500 Hz or the resonance frequenecy. These are
    %required for Popt finding.
    
    TPopt = zeros(Ntemperatures,4);
    %For the active KID a list of temperatures (dim2=1) with corresponding
    %Tindices (dim2=2) for the T index inside NOISE struct fields. It also
    %gives the Popt values (dim2=3) and indices (dim2=4) for NOISE struct.
    
    for nP=1:length(IndexPsort{kidn,1}(:,1)) %LOOP OVER ALL POWERS this kID
        %FIND FOR EACH POWER THE TEMPERATURES MATCHING CLOSEST TO THE
        %TEMPERATURES OF THE REFERENCE POWER
        PTindex(nP,:,1,1) = IndexPsort{kidn,1}(nP,1);%P value
        PTindex(nP,:,1,2) = NOISE(PTindex(nP,1,1,1)).ReadPower(1); %Power always the same in one element of the NOISE struct
        Tref = ones(length(NOISE(IndexPsort{kidn,1}(nP,1)).Temperature),1)*NOISE(IndexPref).Temperature'; %Temperatures @Pref
        Tmeas = (ones(length(NOISE(IndexPref).Temperature),1)*NOISE(IndexPsort{kidn,1}(nP,1)).Temperature')'; %Temperatures @P under investigation
        [dTmin,PTindex(nP,:,2,1)]=min(abs(Tmeas - Tref));
        PTindex(nP,:,2,2) = NOISE(IndexPsort{kidn,1}(nP,1)).Temperature(PTindex(nP,:,2,1));
        if sum(dTmin' > 0.05*NOISE(IndexPref).Temperature) > 0
            disp(['Warning: For readpower ',num2str(NOISE(PTindex(nP,1,1,1)).ReadPower(1)), 'dBm ',...
                num2str(sum(dTmin' > 0.05*NOISE(IndexPref).Temperature)),...
                ' temperatures have more that 5% difference with T@Pref (',num2str(NOISE(IndexPref).ReadPower(1)),' dBm).'])
        end
        %Store the Mean Amplitude noise and resonance frequency in
        %convenient arrays for later Popt finding.
        PT_AmpNoise(nP,:) = NOISE(IndexPsort{kidn,1}(nP,1)).MeanNoise(PTindex(nP,:,2,1),1);%is 50-300Hz amp noise
        PT_Fres(nP,:) = NOISE(IndexPsort{kidn,1}(nP,1)).Fres(PTindex(nP,:,2,1),1);
    end
    IndexPsort{kidn,2} = PTindex; %Store for later checks
    
    for nT=1:Ntemperatures %LOOP OVER ALL TEMPERATURES AVAILABLE TO THE KID AT REFERENCE POWER
        %==================================================================
        % DETERMINE Popt
        %==================================================================
        %Determine which powers have to big a resonance frequency shift
        FresPref = NOISE(IndexPref).Fres(nT);
        Frange = BWrange*NOISE(IndexPref).Bandwidth(nT);
        AllowedPindices = FresPref-Frange <= PT_Fres(:,nT) & ...
            PT_Fres(:,nT) <= FresPref+Frange;
        SkippedPindices = ~AllowedPindices;
        AllowedPindices = find(AllowedPindices);
        
        %Display information on temperature and skipped powers to screen
        disp(['Searchin optimum power at T: ',num2str(min(PTindex(:,nT,2,2))),' - ',num2str(max(PTindex(:,nT,2,2))),' K'])
        if sum(SkippedPindices) > 0
            disp(['Powers excluded due to KID dip frequency drift: ',num2str(PTindex(SkippedPindices,nT,1,2)'),' dBm'])
        end
        
        %Determine Popt in the allowed P range using AmpNoise.
        [~,PoptIndex] = min(PT_AmpNoise(AllowedPindices,nT));
        PoptIndex = AllowedPindices(PoptIndex);
        
        %store in PT index the Popt refs from Popt.csv
        if ReadPoptfile==1
            Poptdata_Tindices= PoptData(:,1)>=0.9*NOISE(IndexPref).Temperature(nT,1) & PoptData(:,1)<=1.1*NOISE(IndexPref).Temperature(nT,1);%Note Ntemperatures = length(NOISE(IndexPref).Temperature);
            Poptdata_thisKID_thisT=find(PoptData(Poptdata_Tindices,2)==NOISE(IndexPref).KIDnumber);
            bla=find(Poptdata_Tindices);%this is an array of indiced in PoptData that are allowed wrt temperature
            Popt_fromfile=PoptData(bla(Poptdata_thisKID_thisT),3);
            PoptIndex=find(PTindex(:,nT,1,2)==Popt_fromfile);
        end
        
        %Construct a Temperature Popt lookup table.
        TPopt(nT,1) = PTindex(PoptIndex,nT,2,2); %T value
        TPopt(nT,2) = PTindex(PoptIndex,nT,2,1); %T index
        TPopt(nT,3) = PTindex(PoptIndex,nT,1,2); %P value
        TPopt(nT,4) = PTindex(PoptIndex,nT,1,1); %P index
        
            
            
            
    end %END OF LOOP OVER ALL TEMPERATURES
    
    for nP=1:length(IndexPsort{kidn,1}(:,1)) %LOOP OVER ALL POWERS
        NOISE(IndexPsort{kidn,1}(nP,1)).OptimumPower = TPopt;
        [~,BaseTindex] = min(NOISE(IndexPsort{kidn,1}(nP,1)).Temperature);
        NOISEParameters(IndexPsort{kidn,1}(nP,1),8) = TPopt(BaseTindex,3);
    end
    %Store TPopt for easy writing to Popt file
    LookUpTables{kidn,1}=TPopt;
    
    %======================================================================
    % Fit the 1 kHz data to obtain the noise at Pint = -40 dBm 
    %======================================================================
    SfPi = zeros(length(IndexPsort{kidn,1}(:,1)),3,Ntemperatures);
    %Stores Pi (dim2=1), Sf/f2 (dim2=2) and if P<Popt (dim2=3 1 or 0)
    %This is stored for all temperatures
    for nP=1:length(IndexPsort{kidn,1}(:,1)) %create arrays vs power for all temperatures (note we are in a KIDn loop)
        SfPi(nP,1,:) = NOISE(IndexPsort{kidn,1}(nP,1)).InternalPower(PTindex(nP,:,2,1),1);
        SfPi(nP,2,:) = NOISE(IndexPsort{kidn,1}(nP,1)).ReadPower(PTindex(nP,:,2,1),1) <= TPopt(:,3);
        SfPi(nP,3,:) = NOISE(IndexPsort{kidn,1}(nP,1)).MeanFreqNoise(PTindex(nP,:,2,1),4);%1kHz
        SfPi(nP,4,:) = NOISE(IndexPsort{kidn,1}(nP,1)).MeanFreqNoise(PTindex(nP,:,2,1),3);%100Hz
        SfPi(nP,5,:) = NOISE(IndexPsort{kidn,1}(nP,1)).MeanFreqNoise(PTindex(nP,:,2,1),2);%10Hz
        SfPi(nP,6,:) = NOISE(IndexPsort{kidn,1}(nP,1)).MeanFreqNoise(PTindex(nP,:,2,1),1);%1Hz
    end
    
    FitOptionsS1kHz = fitoptions('Method','NonlinearLeastSquares','StartPoint',[-0.5 -100]);
    FitTypeS1kHz = fittype('a*(x+40)+b','options',FitOptionsS1kHz);
    
    FitS1kHz = cell(Ntemperatures,1);
    %Stores the fit to the 1 kHz frequency noise as a function of internal
    %power.
    FreqNoise1kHz = zeros(Ntemperatures,2);
    FreqNoise100Hz = zeros(Ntemperatures,2);
    FreqNoise10Hz = zeros(Ntemperatures,2);
    FreqNoise1Hz = zeros(Ntemperatures,2);
    for nT=1:Ntemperatures
        PsubPopt = SfPi(:,2,nT) == 1;
        if sum(PsubPopt) > 1
            FitS1kHz{nT} = fit(SfPi(PsubPopt,1,nT),10*log10(SfPi(PsubPopt,3,nT)),FitTypeS1kHz);
            FitS100Hz{nT} = fit(SfPi(PsubPopt,1,nT),10*log10(SfPi(PsubPopt,4,nT)),FitTypeS1kHz);
            FitS10Hz{nT} = fit(SfPi(PsubPopt,1,nT),10*log10(SfPi(PsubPopt,5,nT)),FitTypeS1kHz);
            FitS1Hz{nT} = fit(SfPi(PsubPopt,1,nT),10*log10(SfPi(PsubPopt,6,nT)),FitTypeS1kHz);
            FreqNoise1kHz(nT,1) = FitS1kHz{nT,1}.b;
            FreqNoise1kHz(nT,2) = FitS1kHz{nT,1}.a;
        else
            %Not enought points for a linear fit. Leave all empty.
            FitS1kHz{nT} = [];
            FitS100Hz{nT} = [];
            FitS10Hz{nT} = [];
            FitS1Hz{nT} = [];
            FreqNoise1kHz(nT,1) = 0;
            FreqNoise1kHz(nT,2) = 0;
        end
    end
    
    for nP=1:length(IndexPsort{kidn,1}(:,1))% strtoes the fit; note that the parameter is the same for all powers 
        NOISE(IndexPsort{kidn,1}(nP,1)).FitS1kHz = FitS1kHz;
        NOISE(IndexPsort{kidn,1}(nP,1)).FitS100Hz = FitS100Hz;
        NOISE(IndexPsort{kidn,1}(nP,1)).FitS10Hz = FitS10Hz;
        NOISE(IndexPsort{kidn,1}(nP,1)).FitS1Hz = FitS1Hz;
    end
    
    %======================================================================
    % Plot for each (KID,T) combination the Power dependencies
    %======================================================================
    [~,BaseTindex] = min(NOISE(IndexPref).Temperature);
    %Make an array with all variables to make it convenient for plotting.
    %(Cannot call multiple indexes in struct array. For example
    %NOISE(1:2).Temperatures(1) is impossible).
    
    DataArray2Plot = zeros(length(IndexPsort{kidn,1}(:,1)),5,Ntemperatures);
    for nP=1:length(IndexPsort{kidn,1}(:,1))
        DataArray2Plot(nP,1,:) = NOISE(IndexPsort{kidn,1}(nP,1)).ReadPower(PTindex(nP,:,2,1),1) <= TPopt(:,3); %P above/below Popt
        DataArray2Plot(nP,2,:) = NOISE(IndexPsort{kidn,1}(nP,1)).InternalPower(PTindex(nP,:,2,1),1); %Pint
        DataArray2Plot(nP,3,:) = NOISE(IndexPsort{kidn,1}(nP,1)).Fres(PTindex(nP,:,2,1),1); %Fres
        DataArray2Plot(nP,4,:) = NOISE(IndexPsort{kidn,1}(nP,1)).Qi(PTindex(nP,:,2,1),1); %Qi
        DataArray2Plot(nP,5,:) = NOISE(IndexPsort{kidn,1}(nP,1)).MeanFreqNoise(PTindex(nP,:,2,1),4); %Frequency Noise @ 1kHz
        DataArray2Plot(nP,6,:) = NOISE(IndexPsort{kidn,1}(nP,1)).ReadPower(PTindex(nP,:,2,1),1) <= TPopt(:,3)+4; %P above/below Popt+4dB
    end
    
    for nT=1:Ntemperatures
        %==================================================================
        % Figure one: mainly checks for the quality of the analysis
        % routine
        %==================================================================
        [Pcolors,~] = GenerateColorsFromMap(length(IndexPsort{kidn,1}(:,1)),'RainbowReinier');
        
        figurehandle1 = 10000*KIDnumbers(kidn)+nT*10+1;
        figure(figurehandle1)
        clf
        
        clear PowerLegend
        PowerLegend = cell(length(IndexPsort{kidn,1}(:,1))+1,2);
        PowerLegend{1,1} = 'P_{opt}';
        PowerLegend{1,2} = 'P^{int}_{opt}';
        %Resonance circle as a function of power (incl noise blobs)
        subplot(2,3,1)
        hold on
        for nP=1:length(IndexPsort{kidn,1}(:,1))
            plot(NOISE(IndexPsort{kidn,1}(nP,1)).S21_IQplane{PTindex(nP,nT,2,1)}(:,2),...
                NOISE(IndexPsort{kidn,1}(nP,1)).S21_IQplane{PTindex(nP,nT,2,1)}(:,3),...
                '-','color',Pcolors(nP,:),'LineWidth',1)
            PowerLegend{nP+1,1} = num2str(NOISE(IndexPsort{kidn,1}(nP,1)).ReadPower(PTindex(nP,nT,2,1),1),'%2.0f');
            PowerLegend{nP+1,2} = num2str(NOISE(IndexPsort{kidn,1}(nP,1)).InternalPower(PTindex(nP,nT,2,1),1),'%2.0f');
        end
        for nP=1:length(IndexPsort{kidn,1}(:,1)) %Second Loop to get the legend correct.
           plot(NOISE(IndexPsort{kidn,1}(nP,1)).TDIQ{PTindex(nP,nT,2,1)}(:,2),...
                NOISE(IndexPsort{kidn,1}(nP,1)).TDIQ{PTindex(nP,nT,2,1)}(:,3),...
                '.','color',Pcolors(nP,:),'MarkerSize',6)
        end
        legend(PowerLegend(2:end,1))
        xlabel('Re')
        ylabel('Im')
        title(['KID ',num2str(KIDnumbers(kidn),'%.0f'),' @T=',num2str(NOISE(IndexPref).Temperature(nT),'%.3g'),' K'])
        hold off
        
        %Resonance Dip as a function of power, Incl reference lines around
        %reference power
        subplot(2,3,2)
        hold on
        for nP=1:length(IndexPsort{kidn,1}(:,1))
            plot(NOISE(IndexPsort{kidn,1}(nP,1)).S21_MPplane{PTindex(nP,nT,2,1)}(:,1),...
                20*log10(NOISE(IndexPsort{kidn,1}(nP,1)).S21_MPplane{PTindex(nP,nT,2,1)}(:,2)),...
                '-','color',Pcolors(nP,:),'LineWidth',1)
        end
        Flow = NOISE(IndexPref).Fres(nT)-BWrange*NOISE(IndexPref).Bandwidth(nT);
        plot([Flow,Flow],[0,-20],'r-','LineWidth',2)
        Fhigh = NOISE(IndexPref).Fres(nT)+BWrange*NOISE(IndexPref).Bandwidth(nT);
        plot([Fhigh,Fhigh],[0,-20],'r-','LineWidth',2)
        xlabel('F [GHz]')
        ylabel('|S21| [dB]')
        title(['KID ',num2str(KIDnumbers(kidn),'%.0f'),' P_{opt}=',num2str(TPopt(nT,3)),' dBm'])
        axis tight;
        hold off
       
        %Time Domain Trace as a function of time
        subplot(2,3,3)
        hold on
        for nP=1:length(IndexPsort{kidn,1}(:,1))
            plot(NOISE(IndexPsort{kidn,1}(nP,1)).TDIQ{PTindex(nP,nT,2,1)}(:,1),...
                NOISE(IndexPsort{kidn,1}(nP,1)).TDIQ{PTindex(nP,nT,2,1)}(:,2),...
                'o','color',Pcolors(nP,:),'MarkerSize',3,'MarkerFaceColor',Pcolors(nP,:))
            plot(NOISE(IndexPsort{kidn,1}(nP,1)).TDIQ{PTindex(nP,nT,2,1)}(:,1),...
                NOISE(IndexPsort{kidn,1}(nP,1)).TDIQ{PTindex(nP,nT,2,1)}(:,3),...
                'd','color',Pcolors(nP,:),'MarkerSize',3,'MarkerFaceColor',Pcolors(nP,:))
        end
        xlabel('t [sec]')
        ylabel('Re or Im')
        legend('Re(S21)','Im(S21)')
        hold off
        
        %Frequency Noise
        subplot(2,3,4)
        semilogx(NOISE(TPopt(nT,4)).FFTnoise{TPopt(nT,2)}(:,1),...
                10*log10(NOISE(TPopt(nT,4)).FFTnoise{TPopt(nT,2)}(:,4)),...
                '-','color','k','LineWidth',3)
        hold on
        for nP=1:length(IndexPsort{kidn,1}(:,1))
            semilogx(NOISE(IndexPsort{kidn,1}(nP,1)).FFTnoise{PTindex(nP,nT,2,1)}(:,1),...
                10*log10(NOISE(IndexPsort{kidn,1}(nP,1)).FFTnoise{PTindex(nP,nT,2,1)}(:,4)),...
                '-','color',Pcolors(nP,:),'LineWidth',1)
        end
        xlabel('F [Hz]')
        ylabel('S_F/F^2 [dBc/Hz]')
        legend(PowerLegend(:,2))
        title(['KID ',num2str(KIDnumbers(kidn),'%.0f'),' P^{int}_{opt}=',num2str(NOISE(TPopt(nT,4)).InternalPower(TPopt(nT,2),1)),' dBm'])
        xlim([0.1,1e6])
        ylim([-260,-140])
        hold off
        
        
        %Amplitude Noise
        subplot(2,3,5)
        semilogx(NOISE(TPopt(nT,4)).FFTnoise{TPopt(nT,2)}(:,1),...
                NOISE(TPopt(nT,4)).FFTnoise{TPopt(nT,2)}(:,3),...
                '-','color','k','LineWidth',3)
        hold on
        for nP=1:length(IndexPsort{kidn,1}(:,1))
            semilogx(NOISE(IndexPsort{kidn,1}(nP,1)).FFTnoise{PTindex(nP,nT,2,1)}(:,1),...
                NOISE(IndexPsort{kidn,1}(nP,1)).FFTnoise{PTindex(nP,nT,2,1)}(:,3),...
                '-','color',Pcolors(nP,:),'LineWidth',1)
            semilogx(100,NOISE(IndexPsort{kidn,1}(nP,1)).MeanNoise(PTindex(nP,nT,2,1),1),...
                'kd','MarkerSize',7,'MarkerFaceColor',Pcolors(nP,:))
        end
        xlabel('F [Hz]')
        ylabel('S_A [dBc/Hz]')
        xlim([0.1,1e6])
        ylim([-120,-30])
        hold off
        
        %Phase Noise (and Amp Noise)
        subplot(2,3,6)
        semilogx(NOISE(IndexPsort{kidn,1}(nP,1)).FFTnoise{PTindex(nP,nT,2,1)}(:,1),...
            NOISE(TPopt(nT,4)).FFTnoise{TPopt(nT,2)}(:,3),...
            '--','color','k','LineWidth',3)
        hold on
        semilogx(NOISE(TPopt(nT,4)).FFTnoise{TPopt(nT,2)}(:,1),...
            NOISE(TPopt(nT,4)).FFTnoise{TPopt(nT,2)}(:,2),...
            '-','color','k','LineWidth',3)
        for nP=1:length(IndexPsort{kidn,1}(:,1))
            semilogx(NOISE(IndexPsort{kidn,1}(nP,1)).FFTnoise{PTindex(nP,nT,2,1)}(:,1),...
                NOISE(IndexPsort{kidn,1}(nP,1)).FFTnoise{PTindex(nP,nT,2,1)}(:,3),...
                '--','color',Pcolors(nP,:),'LineWidth',1)
            semilogx(NOISE(IndexPsort{kidn,1}(nP,1)).FFTnoise{PTindex(nP,nT,2,1)}(:,1),...
                NOISE(IndexPsort{kidn,1}(nP,1)).FFTnoise{PTindex(nP,nT,2,1)}(:,2),...
                '-','color',Pcolors(nP,:),'LineWidth',1)
        end
        xlabel('F [Hz]')
        ylabel('S_x [dBc/Hz]')
        legend('S_A','S_{\theta}')
        xlim([0.1,1e6])
        ylim([-120,-30])
        hold off
        
        %SAVE the figure
        
        Figfile=[Noisepath,'KID',num2str(KIDnumbers(kidn),'%.0f'),'_',...
            num2str(NOISE(IndexPref).Temperature(nT),'%.3g'),'K_NOISE1'];
        MakeGoodFigure(15,12,14,Figfile)
        %saveas(gcf,Figfile,'fig')

        %Clean up if desired.
        if T0only && nT ~= BaseTindex  
            %close the figure
            close(figurehandle1);
        end
        
        %==================================================================
        % Figure two: Overview of resonator parameters as a function of
        % power.
        %==================================================================
        figurehandle2 = 10000*KIDnumbers(kidn)+nT*10+2;
        figure(figurehandle2)
        clf
        
        %Frequency Noise as a function of internal Power
        PsubPopt = find(DataArray2Plot(:,1,nT) == 1);
        PsuperPopt = find(DataArray2Plot(:,1,nT) == 0);
        Py =(-70:-30);
        if isempty(FitS1kHz{nT,1})
            FitS1kHzEval = [];
        else
            FitS1kHzEval = feval(FitS1kHz{nT,1},Py);
        end
        
        %Make the legend for the noise figure
        clear NoiseLegend
        NoiseLegend{1,1} = 'P<P_{opt}'; %Since there is always at least 1 measurement, which is Popt
        NoiseLegend{2,1} = 'P=P_{opt}'; %Since there is always at least 1 measurement, which is Popt
        NLindex = 3;
        if isempty(PsuperPopt) == 0
            NoiseLegend{NLindex,1} = 'P>P_{opt}';
            NLindex = NLindex+1;
        end 
        if isempty(FitS1kHzEval) == 0
            NoiseLegend{NLindex,1} = 'Fit';
            NLindex = NLindex+1;
        end
        NoiseLegend{NLindex,1} = 'Gao Line';
        
        subplot(2,2,1) %Frequency Noise at 1 kHz as a function of Pint
        %Measurement Points
        semilogy(DataArray2Plot(PsubPopt,2,nT),DataArray2Plot(PsubPopt,5,nT),'ro','MarkerSize',6) %Below Popt
        hold on
        semilogy(DataArray2Plot(PsubPopt(end),2,nT),DataArray2Plot(PsubPopt(end),5,nT),'ro','MarkerSize',6,'MarkerFaceColor','r') %At Popt
        if isempty(PsuperPopt) == 0
            semilogy(DataArray2Plot(PsuperPopt,2,nT),DataArray2Plot(PsuperPopt,5,nT),'kd','MarkerSize',6) %Above Popt
        end
        %Lines
        if isempty(FitS1kHzEval) == 0 %Only if fit has been made
            semilogy(Py,10.^(FitS1kHzEval/10),'r-','LineWidth',1) %Fit to P<Popt
        end
        semilogy(Py,10.^((-189-0.5*(Py+40))/10),'b-','LineWidth',1) %Jiansong Line
        %make a nice plot.
        axis tight
        title(['KID ',num2str(KIDnumbers(kidn),'%.0f'),' @T=',num2str(NOISE(IndexPref).Temperature(nT),'%.3g'),' K'])
        legend(NoiseLegend)
        xlabel('P_{int} (dBm)')
        ylabel('S_f/f^2 @1kHz [dB]')
        hold off
        
        subplot(2,2,2) %Resonance Frequency as a function of Pint
        hold on
        %Measurement Points
        plot(DataArray2Plot(PsubPopt,2,nT),DataArray2Plot(PsubPopt,3,nT),'ro','MarkerSize',6) %Below Popt
        plot(DataArray2Plot(PsuperPopt,2,nT),DataArray2Plot(PsuperPopt,3,nT),'kd','MarkerSize',6) %Above Popt
        plot(DataArray2Plot(PsubPopt(end),2,nT),DataArray2Plot(PsubPopt(end),3,nT),'ro','MarkerSize',6,'MarkerFaceColor','r') %At Popt
        xlabel('P_{int} (dBm)')
        ylabel('f_{res} (GHz)')
        title(['KID ',num2str(KIDnumbers(kidn),'%.0f'),' P_{opt}=',num2str(TPopt(nT,3)),' dBm'])
        hold off
        
        subplot(2,2,3) %Qi as a function of Pint
        hold on
        %Measurement Points
        plot(DataArray2Plot(PsubPopt,2,nT),DataArray2Plot(PsubPopt,4,nT),'ro','MarkerSize',6) %Below Popt
        plot(DataArray2Plot(PsuperPopt,2,nT),DataArray2Plot(PsuperPopt,4,nT),'kd','MarkerSize',6) %Above Popt
        plot(DataArray2Plot(PsubPopt(end),2,nT),DataArray2Plot(PsubPopt(end),4,nT),'ro','MarkerSize',6,'MarkerFaceColor','r') %At Popt
        xlabel('P_{int} (dBm)')
        ylabel('Q_{i} (dBm)')
        hold off
        
        subplot(2,2,4) %S21 resonance fit at Popt
        
        %Measurement Points
        plot(NOISE(TPopt(nT,4)).S21_MPplane{TPopt(nT,2),1}(:,1),NOISE(TPopt(nT,4)).S21_MPplane{TPopt(nT,2),1}(:,2),'b.','MarkerSize',5)
        axis tight;hold on
        %Fit
        plot(NOISE(TPopt(nT,4)).S21fit{TPopt(nT,2),1}(:,1),10.^(NOISE(TPopt(nT,4)).S21fit{TPopt(nT,2),1}(:,2)/20),'r-','LineWidth',1)
        plot(NOISE(TPopt(nT,4)).Fres(TPopt(nT,2),1),10.^(NOISE(TPopt(nT,4)).S21min(TPopt(nT,2),1)/20),'kd',...
            'MarkerFaceColor','r','MarkerSize',6)%The determined Fres and |S21(Fres)|
        %make the plot nice
        xlabel('F (GHz)')
        ylabel('|S21|')
        title('Resonance @Popt')
        hold off
        
        %SAVE the figure
        Figfile=[Noisepath,'KID',num2str(KIDnumbers(kidn),'%.0f'),'_',...
            num2str(NOISE(IndexPref).Temperature(nT),'%.3g'),'K_NOISE2.fig'];
        saveas(gcf,Figfile,'fig')

        %Clean up if desired.
        if T0only && nT ~= BaseTindex  
            %close the figure
            close(figurehandle2);
        end
        
    end
    
end %END OF LOOP OVER ALL UNIQUE KIDS

%==========================================================================
% Write the Popt file
%==========================================================================
PoptArrayLength = 0;
for kidn=1:length(KIDnumbers)
    PoptArrayLength = PoptArrayLength + size(LookUpTables{kidn,1}(:,:),1);
end
PoptData = zeros(PoptArrayLength,3);
LengthStart = 1;
LengthStop = 0;
for kidn=1:length(KIDnumbers)
    LengthStop = LengthStop + size(LookUpTables{kidn,1}(:,:),1);
    PoptData(LengthStart:LengthStop,1) = LookUpTables{kidn,1}(:,1);
    PoptData(LengthStart:LengthStop,2) = KIDnumbers(kidn);
    PoptData(LengthStart:LengthStop,3) = LookUpTables{kidn,1}(:,3);
    LengthStart = LengthStart + size(LookUpTables{kidn,1}(:,:),1);
end

PoptHeader = {'T (K)','KIDID','Popt (dBm)'}';
WriteSRONcsv([Noisepath,'Popt.csv'],PoptData,PoptHeader,'%.3g')
%==========================================================================
% Wrap up and routine closure
%==========================================================================
clear AllowedPindices BaseTindex CompleteData Data FFT2read FFTfile Fhigh Flow
clear fres FreqNoise1kHz FresPref Header IndexPref KIDlocations LocEndRoot NLindex
clear NoiseLegend Ntemperatures PT_AmpNoise PT_Fres Pcolors PoptIndex Power PowerLegend
clear PsubPopt PsuperPopt DataArray2Plot S21DBfile S21file S21min SfPi SkippedPindices
clear Figfile Frange FitResult FitS1kHz nT nP filtered TDfile TPopt Tcheck Temperature
clear figurehandle1 figurehandle2 fid dt ans Tmeas Tref dTmin 
clear LengthStart LengthStop PoptArrayLength
clear kidn q Q Py
rmpath([pwd,filesep,'subroutines']);
save([Noisepath,'Noise.mat'])
end