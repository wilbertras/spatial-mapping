function ConvertS21_TUD2SRON

%This function searches the user specified directory for S21 measurements
%done at Delft. It then sorts these by KIDnumber and power. And for each
%(#,P) combination creates a single file containing all temperatures in the
%style of SRON. This output file is then good to be analysed by the 
%HybridsS21 function
%
%
%
%THIS ROUTINE MAKE A LOT OFF ASSUMPTIONS. THESE ARE GENERALLY TRUE BUT
%SHOULD BE CHECKED.
%ASSUMPTIONS:
%1: TUD name format is assumed to be: 
%       [#,'-',T[mK],'-',|CH1Atten[dB]|,'-S21.txt'] (for VNA)
%       [#,'-',|P[dBm]|,'-',T[mK],'-S21.txt'] (for IQ mixer)
%2: SRON name format is assumed to be:
%       ['KID',#,'_',|P|,'dBm_',ChipID,'.dat']
%3: Function assumes it is located in the ..\TUD2SRON subdirectory while
%   External subroutines are located in the ..\subroutines subdirectory
%   of the same main directory.
%
% NOTE: This function automatically recognizes if the data was taken with
% VNA or the IQ-mixer.
%
%SUBROUTINES:
% ReadSRONcsvV2 [External]
% ReadHeader
% WriteSRONS21
%
%REQUIRED FILES:
%
%VERSION: 1.0
%   V1.0(2012-08-29,RJ) Start writing this code
%   V1.0(2012-08-30,RJ) Added a bit comments
%
%DATE: August 30, 2012
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
%==========================================================================

%==========================================================================
%USER SPECIFIED INPUT
%==========================================================================
%TUDdatadir == Where all the S21 data taken with TUD systems are located.
TUDdatadir = 'E:\My Documents\PhD\Matlab\Data Analysis\TUD-SRON KID analysis\Dummy Data\TUD\S21 orig';
%Outputdir == Directory where you want the new files to be written.
Outputdir = 'E:\My Documents\PhD\Matlab\Data Analysis\TUD-SRON KID analysis\Dummy Data\TUD\S21 conv';
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

S21data = cell(length(TUDfiles),1);
%cell array to be filled with Nx3 arrays containing S21 transmission.
%Dim 2 of all these arrays will be [F(GHz), |S21|(dB), Phase(rad)]

for p=1:length(TUDfiles) %loop over all files
    %Read in the information from the file
    disp(['Reading file: ',TUDdatadir,filesep,TUDfiles(p).name])
    [header,S21data{p,1}]=ReadSRONcsvV2([TUDdatadir,filesep],TUDfiles(p).name,0);
    %Extract (Fres,Q),P,T from the header
    [NPT.Fres(p),NPT.Q(p,:),NPT.ReadPower(p,1),~,TempT] = ReadHeader(header);
    if isempty(TempT)
        dashlocation = strfind(TUDfiles(p).name,'-');
        NPT.Temperature(p) = str2double(TUDfiles(p).name(dashlocation(1)+1:dashlocation(2)-1))/1000;
    else
        NPT.Temperature(p) = TempT;
    end
    %Extract KID number from filename
    dashlocation = strfind(TUDfiles(p).name,'-');
    %KID number is located before first '-'
    NPT.KIDnumber(p) = str2num(TUDfiles(p).name(1:(dashlocation(1)-1)));
    %Calculate S21min (dB) for later use in SRON script print
    NPT.S21min(p) = 20*log10(abs(NPT.Q(p,1)/NPT.Q(p,2)));
    %Convert S21data first column from Hz to GHz
    S21data{p,1}(:,1) = S21data{p,1}(:,1)/1e9;
    %Convert S21data thrid column from degree to radians
    S21data{p,1}(:,3) = S21data{p,1}(:,3)*pi/180.0;
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
        % For each unique combination of (#,P) write the SRON style file
        %==================================================================
        UniqueNP = find(NPT.KIDnumber == KIDnumbers(n) & NPT.ReadPower(:,2) == Powers(p));
        [~,IndexMinT] = min(NPT.Temperature(UniqueNP));
        IndexMinT = UniqueNP(IndexMinT);
        %Predefine the variables required for WriteSRONS21 subroutine
        OutputS21File = [Outputdir,filesep,'KID',num2str(KIDnumbers(n)),'_',...
            num2str(abs(Powers(p))),'dBm_',ChipID,'.dat'];
        headerstring = [TUDdatadir,filesep,num2str(KIDnumbers(n)),'-*-S21.txt'];

        headerinfo = zeros(6,1);
        headerinfo(1,1) = abs(Powers(p));
        headerinfo(2,1) = NPT.Fres(IndexMinT)/1e9; %Hz to GHz
        headerinfo(3:5,1) = NPT.Q(IndexMinT,[1,3,2]);
        headerinfo(6,1) = NPT.S21min(IndexMinT);
        %Call the subroutine to write the actual output file we are after.
        WriteSRONS21(OutputS21File,headerstring,headerinfo,NPT.Temperature(UniqueNP),S21data(UniqueNP,1))
    end
end

%==========================================================================
% Save and clean up
%==========================================================================
%Save the complete dataspace
save([Outputdir,filesep,'TUD2SRONs21.mat'])

%Clean up subroutines path
rmpath([PATHmain,'subroutines']); %add subroutines path
%==========================================================================
% END OF MAIN ROUTINE
%==========================================================================
end

function [Fres,Q,Pread,Pint,T] = ReadHeader(header)
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

function WriteSRONS21(fullpath,headerstring,headerinfo,T,S21data)
% This function writes to file identified by fullpath a SRON-style S21
% file. This starts with the headerstring (a path), followed by the
% information from the headerinfo. 
% headerinfo is a 6 element vector containing: 
% [|Pread[dBm]|, Fres[GHz], Q, Qc, Qi, S21min[dB]]
%
% This is then followed by all temperatures T with adjoining S21data.
% S21data must be a cellarray with the same length as T, which contains Nx3
% double arrays that have the S21 transmission as [F(GHz),|S21|(dB),phase(rad)]
%
%OUTPUT: THE RESULTING FILE SHOULD LOOK SOMETHING LIKE THE TEMPLATE BELOW
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
%(somehow '\n\n' is needed, RJ does not completely understand)

%Print temperatures in increasing order
[Tsorted,Tindex] = sort(T);
for p = 1:length(Tsorted)
    %For each element in T print:
    fprintf(fid,'\n'); %an empty line
    LineTFormat = 'Temperature in K:%10.6f\n';
    fprintf(fid,LineTFormat,Tsorted(p)); %A temperature line
    fprintf(fid,'GHz	dB	Rad\n'); %A column header line
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