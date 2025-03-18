function ConvertNOISE_TUD2SRON

%This function searches the user specified directory for Noise measurements
%done at Delft. It then sorts these by KIDnumber and power. It then creates
%the 4 noise files made at SRON for each (#,P) combination. A single file 
%of each kind is created containing all temperatures in the style of SRON.
%These output files are then good to be analysed by the NOISEanalysis
%function. Note that not all data required for the SRON analysis is
%measured at Delft. Missing data is instead set to 0.
%
%
%
%THIS ROUTINE MAKE A LOT OFF ASSUMPTIONS. THESE ARE GENERALLY TRUE BUT
%SHOULD BE CHECKED.
%ASSUMPTIONS:
%1: TUD name format is assumed to be: 
%       [#,'-',T[mK],'-',|CH1Atten[dB]|,'-S21.txt'] (S21 for VNA)
%       [#,'-',|P[dBm]|,'-',T[mK],'-S21.txt'] (S21 for IQ mixer)
%       [#,'all','-',|P[dBm]|,'-',T[mK],'.txt'] (FFT of noise)
%               6 columns: [F[Hz],Amp Off res, Phase Off res,Amp On res,Phase On res,Freq. noise] all noises in dB
%       [#,'-',|P[dBm]|,'-',T[mK],'-',X,'-IQ.txt'] 
%               Where X=ON/OFF Contains IQ mixer settings
%2: SRON name format is assumed to be:
%       ['KID',#,'_',|P|,'dBm_',ChipID,'_FFT.dat'] (FFT of noise giving F,Phase,Amp)
%       ['KID',#,'_',|P|,'dBm_',ChipID,'_S21.dat'] (S21 scan in IQ space)
%       ['KID',#,'_',|P|,'dBm_',ChipID,'_S21dB.dat'] (S21 scan in magnitude phase space)
%       ['KID',#,'_',|P|,'dBm_',ChipID,'_td.dat'] (part of the IQ time domain trace)
%3: Function assumes it is located in the ..\TUD2SRON subdirectory while
%   External subroutines are located in the ..\subroutines subdirectory
%   of the same main directory.
%
% NOTE: Each noise measurement has a S21 scan. The S21 scan is used as the
%   reference file (file that issearched for).
% NOTE: This function automatically recognizes if the data was taken with
%   VNA or the IQ-mixer.
% NOTE: ReadSRONcsvV2 is used to read all 3 types of data files (S21,FFT
% and IQ mixer settings). However, because the mixer settings does not have
% actual data, the result of the read in is different. In the IQ mixer
% settings the numeric values are seen as data (instead of part of the
% header in the case of FFT and S21).
%
%
%SUBROUTINES:
% ReadSRONcsvV2 [External]
% ReadHeader
% WriteSRONS21
% WriteSRONFFT
%
%REQUIRED FILES:
%
%VERSION: 1.1
%   V1.0(2012-12-04,RJ) Start writing this code
%   V1.1(2012-12-11,RJ) Build in error / warning messages for missing
%   files.
%
%DATE: December 4, 2012
%AUTHOR: Reinier Janssen
%
%==========================================================================
% EXAMPLE INPUT FILES: A single file per KIDnr-Power-Temperature combination.
%==========================================================================
%   EXAMPLE LAYOUT DELFT S21 FILE: VNA (starts with empty line)
% 
% S21 TRANS FWD
% Center frequency:3.930112212E+9
% Span frequency:9.113461882E+6
% Power (dBm):0.0
% Step Attenuator (dB):10
% IF Bandwidth:3 kHz
% Date:vr 29 jun 2012
% Time:17:38
% Ql:2.986006678E+3
% Qi:355.720108665E+3
% Qc:3.011284177E+3
% Resonance frequency:3.930084871E+9
% Power at chip:-63.007776742E+0
% Pinternal:-30.254677440E+0
% Temperatures:522.000000E-3
% f(Hz) mag(dB) angle(deg)
% 
% 
% 3.9255554807E+9	-8.3828399369E-1	-1.6243411355E+1
% 3.9255600374E+9	-8.3876372587E-1	-1.6270945705E+1
%   END EXAMPLE LAYOUT DELFT S21 FILE: VNA
%
%   EXAMPLE LAYOUT DELFT S21 FILE: IQ-Mixer (starts with empty line)
% 
% Center frequency:4.362423763E+9
% Span frequency:378.773812747E+3
% Number of points:1.500000000E+3
% Power (dBm):16.000000000E+0
% Ql:292.105397524E+3
% Qi:354.358254420E+3
% Qc:1.662734273E+6
% Resonance frequency:4.362366277E+9
% Power at chip:-49.258441039E+0
% Pinternal:-4.117076659E+0
% Temperatures:320.000000E-3
% 
% 
% 4.3622343757E+9	-1.1420761163E-1	-2.3870144019E+0
% 4.3622346283E+9	-1.1193444929E-1	-2.4189451236E+0
%   END EXAMPLE LAYOUT DELFT S21 FILE: IQ-Mixer
%
%   EXAMPLE LAYOUT DELFT FFT NOISE FILE (starts with empty line)
% 
% Ql:13.291482181E+3
% Qi:352.663432376E+3
% Qc:13.812042288E+3
% Resonance frequency:3.975142405E+9
% Pinternal:-23.202118619E+0
% Power at chip:-62.309809123E+0
% Temperatures:350.225756600E-3
% 
% 
% 5.0000000000E+0	-1.1924889606E+1	-1.9372944740E+1	-1.1936411367E+1	-1.8338544189E+0	-9.6346522513E+1
% 1.0000000000E+1	-9.0462258263E+1	-9.6055955359E+1	-1.0637221506E+2	-7.9614731549E+1	-1.7412739964E+2
% 1.5000000000E+1	-9.1209616611E+1	-9.6300931410E+1	-1.0628068423E+2	-7.9197543206E+1	-1.7371021130E+2
% 2.0000000000E+1	-9.1906410708E+1	-9.7536486325E+1	-1.0780093710E+2	-8.1215989181E+1	-1.7572865727E+2
% 2.5000000000E+1	-9.1630991865E+1	-1.0136858377E+2	-1.0782375632E+2	-8.4010670292E+1	-1.7852333839E+2
% 3.0000000000E+1	-9.2663607437E+1	-1.0139562021E+2	-1.0785786785E+2	-8.3342961419E+1	-1.7785562951E+2
%   END EXAMPLE LAYOUT DELFT FFT NOISE FILE
%
%   EXAMPLE LAYOUT DELFT IQ FILE (starts with empty line)
% 
% Frequency (Hz):3.969245022E+9
% Power at chip:-82.304466322E+0
% Sample frequency:250.000000000E+3
% Record length:50.000000000E+3
% Averages or steps:10.000000000E+0
% Pinternal:-43.560767421E+0
% Temperatures:347.533382491E-3
%   END EXAMPLE LAYOUT DELFT IQ FILE
%==========================================================================

%==========================================================================
%USER SPECIFIED INPUT
%==========================================================================
%TUDdatadir == Where all the noise data taken with TUD systems are located.
TUDdatadir = 'K:\ns\nf\nf-shared\Cristina\cristina 2\Cristina\130304_SiN4um\315mK';
%Outputdir == Directory where you want the new files to be written.
%A separate folder is required for the on and off resonance noise data
OutputOndir = 'K:\ns\nf\nf-shared\Cristina\cristina 2\Cristina\130304_SiN4um\315mK\SRON ON';
OutputOffdir = 'K:\ns\nf\nf-shared\Cristina\cristina 2\Cristina\130304_SiN4um\315mK\SRON OFF';
%[Optional,Default=''] ChipID: string ID of the chip. Will be added in the SRON name format.
ChipID = '';

%Minimum size of the jump (in dBm) between 2 different powers at chip level
sigmaP = 0.5; %(dBm).
%Typically anything around 1 is good due to minimum step of 2dB by the step
%attenuator.

%==========================================================================
% Set some internal parameters
%==========================================================================
format('long','e'); %Set display format of numbers to 7 digits

%Enable subroutines by adding path in search path.
PATHcurrent = pwd; %current directory
slashlocations = strfind(PATHcurrent,filesep); %locate the last \ in current directory.
PATHmain = PATHcurrent(1:slashlocations(end));%main routine path will be
addpath([PATHmain,'subroutines']); %add subroutines path

%==========================================================================
% Find all files in the input dir and extract from each file #,P,T,S21data
%==========================================================================
TUDfiles = dir([TUDdatadir,filesep,'*-S21.txt']);

%Preallocate some variables
%The NPT struct containing various things that need to be known for later
%use.
NPT.KIDnumber = zeros(length(TUDfiles),1); %KID number
NPT.ReadPower = zeros(length(TUDfiles),2); %Read power at chip level (dBm)
%Dim 2: 
%   1 == as written in Delft measurement file
%   2 == Mean Read power at chip level rounded to nearest dBm. Here the mean
%           is taken over all powers with similar power level and KID number.
NPT.Temperature = zeros(length(TUDfiles),1); %Temperature (K)
NPT.Fres = zeros(length(TUDfiles),1); %Resonance frequency (Hz)
NPT.Q = zeros(length(TUDfiles),3); %[Q,Qi,Qc]
NPT.S21min = zeros(length(TUDfiles),1); %S21min (dB)
NPT.Fsample = zeros(length(TUDfiles),2); %Sampling frequency of time domain [On,Off]
NPT.Fon = zeros(length(TUDfiles),1); %Frequency at which ON resonance noise is measured
NPT.Foff = zeros(length(TUDfiles),1); %Frequency at which OFF resonance noise is measured
NPT.Poff = zeros(length(TUDfiles),2); %Power [read,internal] during OFF resonance noise measurement

S21data = cell(length(TUDfiles),1);
%cell array to be filled with Nx3 arrays containing S21 transmission.
%Dim 2 of all these arrays will be [F(GHz), |S21|(dB), Phase(rad)]
FFTdata = cell(length(TUDfiles),1);
%cell array to be filled with Nx6 arrays containing FFT noise spectra
%Dim 2 of all these arrays will be (noises in dB)
%[F(Hz), Amp Noise Off res, Phase Noise Off res, Amp Noise On res, Phase Noise On res,Frequency Noise]

%Copies to split in on and off resonance data.
FFTONdata = FFTdata;
FFTOFFdata = FFTdata; 

%Dummy cell arrays (will be filled with zeros) to amend the lack of the
%resonance circle and time domain data.
S21ReIm = cell(length(TUDfiles),1);
TDdata = cell(length(TUDfiles),1);


for p=1:length(TUDfiles) %loop over all files
    %Define the files to be read in:
    S21file = [TUDdatadir,filesep,TUDfiles(p).name];
    dashlocation = strfind(TUDfiles(p).name,'-');
    FFTfile = [TUDdatadir,filesep,TUDfiles(p).name(1:(dashlocation(1)-1)),'all',...
        TUDfiles(p).name(dashlocation(1):(dashlocation(end)-1)),'.txt'];
    IQONfile = [TUDdatadir,filesep,TUDfiles(p).name(1:dashlocation(end)),'ON-IQ.txt'];
    IQOFFfile = [TUDdatadir,filesep,TUDfiles(p).name(1:dashlocation(end)),'OFF-IQ.txt'];
    
    %======================================================================
    %Read in the information from the S21 file
    disp(['Reading file: ',S21file])
    [headerS21,S21data{p,1}]=ReadSRONcsvV2(S21file,'',0);
    %Extract (Fres,Q),P,T from the header
    [NPT.Fres(p),NPT.Q(p,:),NPT.ReadPower(p,1),~,TempT] = ReadS21Header(headerS21);
    if isempty(TempT)
        NPT.Temperature(p) = str2double(TUDfiles(p).name(dashlocation(1)+1:dashlocation(2)-1))/1000;
    else
        NPT.Temperature(p) = TempT;
    end
    %Extract KID number from filename
    %KID number is located before first '-'
    NPT.KIDnumber(p) = str2num(TUDfiles(p).name(1:(dashlocation(1)-1)));
    %Calculate S21min (dB) for later use in SRON script print
    NPT.S21min(p) = 20*log10(abs(NPT.Q(p,1)/NPT.Q(p,2)));
    %Convert S21data first column from Hz to GHz
    S21data{p,1}(:,1) = S21data{p,1}(:,1)/1e9;
    %Convert S21data thrid column from degree to radians
    S21data{p,1}(:,3) = S21data{p,1}(:,3)*pi/180.0;
    
    %Create the resonance circle dummy data
    S21ReIm{p,1} = zeros(size(S21data{p,1}));
    
    %======================================================================
    %Read in the information from the IQ ON mixer file
    disp(['Reading file: ',IQONfile])
    %Check if file exists
    fid = fopen(IQONfile);
    if fid == -1
        %There is no file
        NPT.Fsample(p,1) = 1; %Dummy
        %Print warning
        disp('Warning: This IQ ON txt file has not been found. Using dummy sample frequency.')
    else
        %There is a file
        fclose(fid);
        [headerIQON,IQONdata]=ReadSRONcsvV2(IQONfile,'',0);
        for q=1:length(headerIQON)
            FsampleHere = regexp(headerIQON{q,1},'Sample frequency');
            if isempty(FsampleHere)
            else
                NPT.Fsample(p,1) = IQONdata(q);
            end
            %Rest off the properties should be know from S21 header.
        end
    end
    
    %======================================================================
    %Read in the information from the IQ OFF mixer file
    disp(['Reading file: ',IQOFFfile])
    %Check if file exists
    fid = fopen(IQOFFfile);
    if fid == -1
        %There is no file
        %DUMMIES
        NPT.Fsample(p,2) = 1;
        NPT.Poff(p,1) = 1;
        NPT.Poff(p,2) = 1;
        NPT.Poff(p,1) = 0;
        %Print warning
        disp('Warning: This IQ OFF txt file has not been found. Using dummy frequencies and powers.')
    else
        %The file is there
        fclose(fid);
        [headerIQOFF,IQOFFdata]=ReadSRONcsvV2(IQOFFfile,'',0);
        for q=1:length(headerIQOFF)
            %Sample Frequency
            FsampleHere = regexp(headerIQOFF{q,1},'Sample frequency');
            if isempty(FsampleHere)
            else
                NPT.Fsample(p,2) = IQOFFdata(q);
            end
            
            %Read Power
            PoffHere = regexp(headerIQOFF{q,1},'Power at chip');
            if isempty(PoffHere)
            else
                NPT.Poff(p,1) = IQOFFdata(q);
            end
            
            %Internal Power
            PoffHere = regexp(headerIQOFF{q,1},'Pinternal');
            if isempty(PoffHere)
            else
                NPT.Poff(p,2) = IQOFFdata(q);
            end
            
            %Frequency
            FoffHere = regexp(headerIQOFF{q,1},'Frequency');
            if isempty(FoffHere)
            else
                NPT.Foff(p,1) = IQOFFdata(q);
            end
        end
    end
    %======================================================================
    %Read in the information from the all file (noise spectra)
    disp(['Reading file: ',FFTfile])
    fid = fopen(FFTfile);
    if fid == -1
        %error, cannot proceed without this file
        error('Missing Noise Spectra. Cannot proceed without this file.')
    else
        fclose(fid);
        [~,FFTdata{p,1}]=ReadSRONcsvV2(FFTfile,'',0);
        %We don't need the header (all this info is in S21 header)
        
        %Now split the FFTdata in On and Off resonance
        FFTONdata{p,1} = FFTdata{p,1}(:,[1,5,4]);
        FFTOFFdata{p,1} = FFTdata{p,1}(:,[1,3,2]);
        
        %Some dummy time domain trace
        TDdata{p,1} = zeros(10,2);
    end
end

%==========================================================================
% Identify each unique combination of (#,P)
%==========================================================================
KIDnumbers = unique(NPT.KIDnumber);
for n = 1:length(KIDnumbers) %loop over kid number
    KIDlocation = find(NPT.KIDnumber == KIDnumbers(n)); %find which data is of this KID
    %Sort all powers of this KID and use large jumps to identify unique
    %powers.
    [Psorted,Porder]=sort(NPT.ReadPower(KIDlocation,1)); %sort for increasing P
    dP = Psorted(2:end) - Psorted(1:end-1);
    Pjump = find(dP > sigmaP); %These indexes are the last of the lower power.
    if isempty(Pjump)
        %Only a single power
        Powers = zeros(1,1); %Hard reset to the Powers variable
        Powers(1,1) = round(mean(Psorted));
        NPT.ReadPower(KIDlocation,2) = Powers(1,1);
    else
        %Multiple powers
        %Number of power is number of jumps+1
        Powers = zeros(length(Pjump)+1,1); %Hard reset to the Powers variable
        for p = 1:length(Pjump)+1
            if p == 1
                Powers(p,1) = round(mean(Psorted(1:Pjump(p))));
                NPT.ReadPower(KIDlocation(Porder(1:Pjump(p))),2) = Powers(p,1);
            elseif p == length(Pjump)+1
                Powers(p,1) = round(mean(Psorted(Pjump(p-1)+1:end)));
                NPT.ReadPower(KIDlocation(Porder(Pjump(p-1)+1:end)),2) = Powers(p,1);
            else
                Powers(p,1) = round(mean(Psorted(Pjump(p-1)+1:Pjump(p))));
                NPT.ReadPower(KIDlocation(Porder(Pjump(p-1)+1:Pjump(p))),2) = Powers(p,1);
            end
        end
    end
    
    %ALTERNATIVE: There is an alternative to this uniqueness finding using
    %a combination of round/ceil/floor and unique. However, this leaves
    %potential errors when the read power is floating around x.5 of x.0
    %dBm. Hence it is done in this ugly (but working) manner.
    
    %We now know that KID number "KIDnumbers(n)" has the unique powers "Powers" and
    %that these can be easily reference to NPT(:,4).
    for p = 1:length(Powers)
        %==================================================================
        % For each unique combination of (#,P) write the SRON style files
        %==================================================================
        UniqueNP = find(NPT.KIDnumber == KIDnumbers(n) & NPT.ReadPower(:,2) == Powers(p));
        [~,IndexMinT] = min(NPT.Temperature(UniqueNP));
        IndexMinT = UniqueNP(IndexMinT);
        
        %==================================================================
        %Write the S21dB file
        %Predefine the variables required for WriteSRONS21 subroutine
        OutputS21dBFile = [OutputOndir,filesep,'KID',num2str(KIDnumbers(n)),'_',...
            num2str(abs(Powers(p))),'dBm_',ChipID,'_S21dB.dat'];
        headerstring = [TUDdatadir,filesep,num2str(KIDnumbers(n)),'-*-S21.txt'];

        headerinfo = zeros(6,1);
        headerinfo(1,1) = abs(Powers(p));
        headerinfo(2,1) = NPT.Fres(IndexMinT)/1e9; %Hz to GHz
        headerinfo(3:5,1) = NPT.Q(IndexMinT,[1,3,2]);
        headerinfo(6,1) = NPT.S21min(IndexMinT);
        %Call the subroutine to write the actual output file we are after.
        WriteSRONS21(OutputS21dBFile,headerstring,headerinfo,NPT.Temperature(UniqueNP),S21data(UniqueNP,1),0)
        
        %==================================================================
        %Write the S21 file (with complex circle)
        %Predefine the variables required for WriteSRONS21 subroutine
        OutputS21File = [OutputOndir,filesep,'KID',num2str(KIDnumbers(n)),'_',...
            num2str(abs(Powers(p))),'dBm_',ChipID,'_S21.dat'];
        
        %Call the subroutine to write the actual output file we are after.
        WriteSRONS21(OutputS21File,headerstring,headerinfo,NPT.Temperature(UniqueNP),S21ReIm(UniqueNP,1),1)
        
        %==================================================================
        %Write the FFT file (with the noise spectra) [ON RESONANCE]
        OutputFFTfile = [OutputOndir,filesep,'KID',num2str(KIDnumbers(n)),'_',...
            num2str(abs(Powers(p))),'dBm_',ChipID,'_FFT.dat'];
        FFTheaderstring = [TUDdatadir,filesep,num2str(KIDnumbers(n)),'all-*.txt'];
        
        FFTheaderinfo = headerinfo;
        FFTheaderinfo(6,1) = NPT.Fon(IndexMinT);
        
        %Call the subroutine to write the actual output file we are after.
        WriteSRONFFT(OutputFFTfile,FFTheaderstring,FFTheaderinfo,...
            NPT.Temperature(UniqueNP),FFTONdata(UniqueNP,1))
        
        %==================================================================
        %Write the FFT file (with the noise spectra) [OFF RESONANCE]
        OutputFFTOFF = [OutputOffdir,filesep,'KID',num2str(KIDnumbers(n)),'_',...
            num2str(abs(NPT.Poff(IndexMinT,1)),'%2.0f'),'dBm_',ChipID,'_FFT.dat'];
        
        FFTOFFheaderinfo = FFTheaderinfo;
        FFTOFFheaderinfo(1,1) = abs(NPT.Poff(IndexMinT,1));
        FFTOFFheaderinfo(6,1) = NPT.Foff(IndexMinT);
        
        %Call the subroutine to write the actual output file we are after.
        WriteSRONFFT(OutputFFTOFF,FFTheaderstring,FFTOFFheaderinfo,...
            NPT.Temperature(UniqueNP),FFTOFFdata(UniqueNP,1))
        
        %==================================================================
        %Write the td file (with the time domain traces) [ON RESONANCE]
        OutputTDON = [OutputOndir,filesep,'KID',num2str(KIDnumbers(n)),'_',...
            num2str(abs(Powers(p))),'dBm_',ChipID,'_td.dat'];
        TDONheaderstring = [TUDdatadir,filesep,num2str(KIDnumbers(n)),'-*-ON-IQ.txt'];
        
        TDONheader = zeros(7,1);
        TDONheader(1:5,1) = FFTheaderinfo(1:5,1);
        TDONheader(6,1) = headerinfo(6,1);
        TDONheader(7,1) = 1/NPT.Fsample(IndexMinT,1);
        
        WriteSRONTD(OutputTDON,TDONheaderstring,TDONheader,...
            NPT.Temperature(UniqueNP),TDdata(UniqueNP,1))
        
        %==================================================================
        %Write the td file (with the time domain traces) [OFF RESONANCE]
        OutputTDOFF = [OutputOffdir,filesep,'KID',num2str(KIDnumbers(n)),'_',...
            num2str(abs(Powers(p))),'dBm_',ChipID,'_td.dat'];
        TDOFFheaderstring = [TUDdatadir,filesep,num2str(KIDnumbers(n)),'-*-OFF-IQ.txt'];
        
        TDOFFheader = zeros(7,1);
        TDOFFheader(1:5,1) = FFTOFFheaderinfo(1:5,1);
        TDOFFheader(6,1) = headerinfo(6,1);
        TDOFFheader(7,1) = 1/NPT.Fsample(IndexMinT,2);
        
        WriteSRONTD(OutputTDOFF,TDOFFheaderstring,TDOFFheader,...
            NPT.Temperature(UniqueNP),TDdata(UniqueNP,1))
    end
end

%==========================================================================
% Save and clean up
%==========================================================================
%Save the complete dataspace
save([OutputOndir,filesep,'TUD2SRONnoise.mat'])

%Clean up subroutines path
rmpath([PATHmain,'subroutines']); %add subroutines path
%==========================================================================
% END OF MAIN ROUTINE
%==========================================================================
end

function [Fres,Q,Pread,Pint,T] = ReadS21Header(header)
%this function extracts the internal power, read-out power and temperature
%from the header of the S21 datafiles. This can be done for various different
%headers based on the different devices used in the Delft setup:
%DVNA : Delft VNA (Using ZVMwithControl.vi)
%DIQM : Delft IQmixer (Using S21scopefast.vi)

%NOTE: FOR the noise files use the NoiseHeader function.
%NOTE: Last modification does not using device input anymore. Looks itself
%which device is used.

% a matrix with the correct position of each thing in header
% NOTE: These locations might only be correct after the importing by
% ReadSRONcsv.m
%[Fres,Ql,Qi,Qc,Pread,Pint,T]
locations = [13 10 11 12 14 15 16;  %DVNA
             9 6 7 8 10 11 12]'; %DIQM

%switch device
%    case 'DVNA'
%        device = 1;
%    case 'DIQM'
%        device = 2;
%end

%Use the location of the 'Center frequency' string to identify the device.
CFL = 0;
for p=1:length(header)
    if length(header{p}) >= 16
        if strcmpi(header{p}(1:16),'Center frequency')
            CFL = p;
        end
    end
end

switch CFL
    case 2
        device = 2;
    case 3
        device = 1;
    otherwise
        error('Header from unknown device')
end
Q = zeros(3,1); %[Ql,Qi,Qc]

%==========================================================================
% Resonance Frequency
%==========================================================================
PTstring = header{locations(1,device),1};
locdd = find(PTstring == ':');
Fres = str2num(PTstring(locdd+1:end));

%==========================================================================
% loaded Q
%==========================================================================
PTstring = header{locations(2,device),1};
locdd = find(PTstring == ':');
Q(1,1) = str2num(PTstring(locdd+1:end));

%==========================================================================
% internal Q
%==========================================================================
PTstring = header{locations(3,device),1};
locdd = find(PTstring == ':');
Q(2,1) = str2num(PTstring(locdd+1:end));

%==========================================================================
% coupling Q
%==========================================================================
PTstring = header{locations(4,device),1};
locdd = find(PTstring == ':');
Q(3,1) = str2num(PTstring(locdd+1:end));

%==========================================================================
% Read Out Power
%==========================================================================
PTstring = header{locations(5,device),1};
locdd = find(PTstring == ':');
Pread = str2num(PTstring(locdd+1:end));

%==========================================================================
% Internal Power
%==========================================================================
PTstring = header{locations(6,device),1};
locdd = find(PTstring == ':');
Pint = str2num(PTstring(locdd+1:end));
%==========================================================================
% Temperature
%==========================================================================
PTstring = header{locations(7,device),1};
locdd = find(PTstring == ':');
T = str2num(PTstring(locdd+1:end));
%==========================================================================
% END OF READHEADER SUBROUTINE
%==========================================================================
end

function WriteSRONS21(fullpath,headerstring,headerinfo,T,S21data,IQcircle)
% This function writes to file identified by fullpath a SRON-style S21
% file. This starts with the headerstring (a path), followed by the
% information from the headerinfo. 
% headerinfo is a 6 element vector containing: 
% [|Pread[dBm]|, Fres[GHz], Q, Qc, Qi, S21min[dB]]
% IQcircle = 1 true if IQ data (F,Re,Im) instead of data in the magnitude
% phase  plane.
%
% This is then followed by all temperatures T with adjoining S21data.
% S21data must be a cellarray with the same length as T, which contains Nx3
% double arrays that have the S21 transmission as [F(GHz),|S21|(dB),phase(rad)]
%
%OUTPUT (IQcircle = 0): THE RESULTING FILE SHOULD LOOK SOMETHING LIKE THE TEMPLATE BELOW
%==========================================================================
%   EXAMPLE LAYOUT SRON S21 FILE
%C:\KID Metingen\ADR\B6Ch7__7_11_06 16_02\S21\2D\KID41_73dBm_JB1.dat
%Power at KID:73dBm
%resonance Frequency in GHz :4.893858
%Q=30785.808661, Qc=31780.243467, Qi=983855.842304, S21min=-30.091618
%
%
%Temperature in K:0.099953
%GHz	dB	Rad
%4.893115000	-6.722063065	2.685341760
%4.893118870	-6.721852779	2.685013925
%4.894663000	-6.540477276	2.794459459
%
%Temperature in K:0.105161
%GHz	dB	Rad
%4.893118740	-6.718620300	2.685038958
%4.893122571	-6.718233109	2.684799274
%   END EXAMPLE LAYOUT SRON S21 FILE
%==========================================================================
%OUTPUT (IQcircle = 1): THE RESULTING FILE SHOULD LOOK SOMETHING LIKE THE TEMPLATE BELOW
%==========================================================================
%   EXAMPLE LAYOUT SRON S21 FILE
%C:\KID Metingen\ADR\B6Ch7__7_11_06 16_02\S21\2D\KID41_73dBm_JB1.dat
%Power at KID:73dBm
%resonance Frequency in GHz :4.893858
%Q=30785.808661, Qc=31780.243467, Qi=983855.842304, S21min=-30.091618
%
%
%Temperature in K:0.099953
%GHz	Re	Im
%
%4.893115000	-6.722063065	2.685341760
%4.893118870	-6.721852779	2.685013925
%4.894663000	-6.540477276	2.794459459
%
%Temperature in K:0.105161
%GHz	Re	Im
%
%4.893118740	-6.718620300	2.685038958
%4.893122571	-6.718233109	2.684799274
%   END EXAMPLE LAYOUT SRON S21 FILE
%==========================================================================
fid = fopen(fullpath,'w'); %OVERWRITES EXISTING FILES

%Writing the header string
%This is a bit complicated because matlab recognises '\' in a string as an
%identifier. '\\' is the actual backslash. Since the headerstring is a
%complete filepath, containing a lot off '\', we need to do some clever
%tricks.

%First replace all '\' by '\\' by calling:
headerstring = regexprep(headerstring,'\\','\\\\');
%Note that you ask the replacement of '\' (for matlab '\\') by a double
%slash (for matlab ['\\','\\'])
%This modified headerstring can now be printed using fprintf. fprintf will
%recognise all double \ as identifiers for \ and will only print a single \
%in the file.
fprintf(fid,headerstring);
fprintf(fid,'\n');
%Still don't understand? Play around with strings containing '\' and the
%functions regexprep and fprintf in your maltab window.

%Write the second header line [containing power]
Line2Format = 'Power at KID:%2.0fdBm\n';
fprintf(fid,Line2Format,headerinfo(1));
%Write thrid header line [containing resonance frequency]
Line3Format = 'resonance Frequency in GHz :%11.9f\n';
fprintf(fid,Line3Format,headerinfo(2));
%Write the fourth header line [containing Q's and S21min]
Line4Format = 'Q=%7.6f, Qc=%7.6f, Qi=%7.6f, S21min=%7.5f\n';
fprintf(fid,Line4Format,headerinfo(3),headerinfo(4),headerinfo(5),headerinfo(6));

%NOTE: During HybridS21 analysis only header line 2 is important the rest
%will be skipped. The number of lines for each part is important though.
fprintf(fid,'\n'); %Hence print first empty line

%Print temperatures in increasing order
[Tsorted,Tindex] = sort(T);
for p = 1:length(Tsorted)
    %For each element in T print:
    fprintf(fid,'\n'); %an empty line
    LineTFormat = 'Temperature in K:%10.6f\n';
    fprintf(fid,LineTFormat,Tsorted(p)); %A temperature line
    if IQcircle
        fprintf(fid,'GHz	Re	Im\n\n'); %A column header line for complex circle
    else
        fprintf(fid,'GHz	dB	Rad\n'); %A column header line for magnitude phase
    end
    %And now the actual data
    DataFormat = '%10.9f\t%10.9f\t%10.9f\n';
    fprintf(fid,DataFormat,S21data{Tindex(p),1}(:,:)');    
end

%Close the file
fclose(fid);

%==========================================================================
% END OF WriteSRONS21 SUBROUTINE
%==========================================================================
end

function WriteSRONFFT(fullpath,headerstring,headerinfo,T,FFTdata)
% This function writes to file identified by fullpath a SRON-style FFT
% file. This starts with the headerstring (a path), followed by the
% information from the headerinfo. 
% headerinfo is a 6 element vector containing: 
% [|Pread[dBm]|, Fres[GHz], Q, Qc, Qi, Fmeas[GHz]]
% Here Fmeas is the frequency at which the noise is measured.
%
% This is then followed by all temperatures T with adjoining FFTdata.
% FFTdata must be a cellarray with the same length as T, which contains Nx3
% double arrays that have the Noise spectra as [F(Hz),Phase(dB),Amp(rad)]
%
%OUTPUT: THE RESULTING FILE SHOULD LOOK SOMETHING LIKE THE TEMPLATE BELOW
%==========================================================================
%   EXAMPLE LAYOUT SRON FFT FILE
% C:\KID metingen\ADR Box in Box\H20_1 APEX 13_2_12\FFT\Power\KID3_92dBm__FFT.dat
% Power at KID:92dBm
% resonance Frequency in GHz :5.755258665F used in GHz :5.755258665
% Q=31099.980177, Qc=32808.413552, Qi=597237.812105
% KID_LNA -4, LNA 4, GainAMP2 35, W.ch1 10, sys Gain KID_IQin 83, Tsys 292, NoiseFloor -80.2
% 
% 
% Temperature in K:0.100018
% Hz	Phase noise	Amp noise
% 0.000000000	-42.183568303	-0.585862560
% 0.762939453	-45.093054177	-3.596210471
% 1.525878906	-66.620912819	-88.229268984
%
%Temperature in K:0.105161
%GHz	dB	Rad
% 0.000000000	-42.183568303	-0.585862560
% 0.762939453	-45.093054177	-3.596210471
% 1.525878906	-66.620912819	-88.229268984
%   END EXAMPLE LAYOUT SRON FFT FILE
%==========================================================================
fid = fopen(fullpath,'w'); %OVERWRITES EXISTING FILES

%Writing the header string
%This is a bit complicated because matlab recognises '\' in a string as an
%identifier. '\\' is the actual backslash. Since the headerstring is a
%complete filepath, containing a lot off '\', we need to do some clever
%tricks.

%First replace all '\' by '\\' by calling:
headerstring = regexprep(headerstring,'\\','\\\\');
%Note that you ask the replacement of '\' (for matlab '\\') by a double
%slash (for matlab ['\\','\\'])
%This modified headerstring can now be printed using fprintf. fprintf will
%recognise all double \ as identifiers for \ and will only print a single \
%in the file.
fprintf(fid,headerstring);
fprintf(fid,'\n');
%Still don't understand? Play around with strings containing '\' and the
%functions regexprep and fprintf in your maltab window.

%Write the second header line [containing power]
Line2Format = 'Power at KID:%2.0fdBm\n';
fprintf(fid,Line2Format,headerinfo(1));
%Write thrid header line [containing resonance frequency]
Line3Format = 'resonance Frequency in GHz :%11.9f F used in GHz:%11.9f\n';
fprintf(fid,Line3Format,headerinfo(2),headerinfo(6));
%Write the fourth header line [containing Q's]
Line4Format = 'Q=%7.6f, Qc=%7.6f, Qi=%7.6f\n';
fprintf(fid,Line4Format,headerinfo(3),headerinfo(4),headerinfo(5));
Line5Format = 'KID_LNA 0, LNA 0, GainAMP2 00, W.ch1 00, sys Gain KID_IQin 00, Tsys 000, NoiseFloor 00.0\n';
fprintf(fid,Line5Format);

%NOTE: During Noise analysis only header line 2 and 3 are important the rest
%will be skipped. The number of lines for each part is important though.
%This might change in the future
fprintf(fid,'\n'); %Hence print first empty line 

%Print temperatures in increasing order
[Tsorted,Tindex] = sort(T);
for p = 1:length(Tsorted)
    %For each element in T print:
    fprintf(fid,'\n'); %an empty line
    LineTFormat = 'Temperature in K:%10.6f\n';
    fprintf(fid,LineTFormat,Tsorted(p)); %A temperature line
    fprintf(fid,'Hz	Phase noise	Amp noise\n'); %A column header line
    %And now the actual data
    DataFormat = '%10.9f\t%10.9f\t%10.9f\n';
    fprintf(fid,DataFormat,FFTdata{Tindex(p),1}(:,:)');    
end

%Close the file
fclose(fid);

%==========================================================================
% END OF WriteSRONFFT SUBROUTINE
%==========================================================================
end

function WriteSRONTD(fullpath,headerstring,headerinfo,T,TDdata)
% This function writes to file identified by fullpath a SRON-style time
% domain file. This starts with the headerstring (a path), followed by the
% information from the headerinfo. 
% headerinfo is a 7 element vector containing: 
% [|Pread[dBm]|, Fres[GHz], Q, Qc, Qi, S21min[dB], dt [sec]]
% Here dt is the timestep between sequential (I,Q) measurement points.
%
% This is then followed by all temperatures T with adjoining time domain 
% trace. The TDdata must be a cellarray with the same length as T, which
% contains Nx2 double arrays that have the I and Q as a function of time
% [I,Q]
%
%OUTPUT: THE RESULTING FILE SHOULD LOOK SOMETHING LIKE THE TEMPLATE BELOW
%==========================================================================
%   EXAMPLE LAYOUT SRON TD FILE
% C:\KID metingen\ADR Box in Box\k98w1_1_6THz_19-9-2012\twoSelectedKIDs\FFT\2D\KID1_76dBm__td.dat
% Power at KID:76dBm
% resonance Frequency in GHz :5.318665
% Q=50907.280548, Qc=151393.852988, Qi=76697.305522  S21min=-3.560004
% 
% 
% Temperature in K:0.116421
% I	Q	dt= 0.000020
% -0.256526093	0.331262266
% -0.255443434	0.331289519
% -0.257608966	0.331132319
%
% Temperature in K:0.135205
% I	Q	dt= 0.000020
% -0.256014504	0.332769778
% -0.257341317	0.332118123
% -0.260047902	0.332367196
%   END EXAMPLE LAYOUT SRON FFT FILE
%==========================================================================
fid = fopen(fullpath,'w'); %OVERWRITES EXISTING FILES

%Writing the header string
%This is a bit complicated because matlab recognises '\' in a string as an
%identifier. '\\' is the actual backslash. Since the headerstring is a
%complete filepath, containing a lot off '\', we need to do some clever
%tricks.

%First replace all '\' by '\\' by calling:
headerstring = regexprep(headerstring,'\\','\\\\');
%Note that you ask the replacement of '\' (for matlab '\\') by a double
%slash (for matlab ['\\','\\'])
%This modified headerstring can now be printed using fprintf. fprintf will
%recognise all double \ as identifiers for \ and will only print a single \
%in the file.
fprintf(fid,headerstring);
fprintf(fid,'\n');
%Still don't understand? Play around with strings containing '\' and the
%functions regexprep and fprintf in your maltab window.

%Write the second header line [containing power]
Line2Format = 'Power at KID:%2.0fdBm\n';
fprintf(fid,Line2Format,headerinfo(1));
%Write thrid header line [containing resonance frequency]
Line3Format = 'resonance Frequency in GHz :%11.9f\n';
fprintf(fid,Line3Format,headerinfo(2));
%Write the fourth header line [containing Q's]
Line4Format = 'Q=%7.6f, Qc=%7.6f, Qi=%7.6f, S21min=%7.5f\n';
fprintf(fid,Line4Format,headerinfo(3),headerinfo(4),headerinfo(5),headerinfo(6));

%NOTE: During Noise analysis only header line 2 and 3 are important the rest
%will be skipped. The number of lines for each part is important though.
%This might change in the future
fprintf(fid,'\n'); %Hence print first empty line 

%Print temperatures in increasing order
[Tsorted,Tindex] = sort(T);
for p = 1:length(Tsorted)
    %For each element in T print:
    fprintf(fid,'\n'); %an empty line
    LineTFormat = 'Temperature in K:%10.6f\n';
    fprintf(fid,LineTFormat,Tsorted(p)); %A temperature line
    LineColFormat = 'I	Q	dt= %7.6f\n';
    fprintf(fid,LineColFormat,headerinfo(7)); %A column header line
    %And now the actual data
    DataFormat = '%10.9f\t%10.9f\n';
    fprintf(fid,DataFormat,TDdata{Tindex(p),1}(:,:)');    
end

%Close the file
fclose(fid);

%==========================================================================
% END OF WriteSRONFFT SUBROUTINE
%==========================================================================
end