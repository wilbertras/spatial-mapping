function [KID] = BB_S21analysisV5
% Jochem 12-12-2015. Taken from S21analysis V5 this script reads and
% processes S21(Tbb) data.
% This function will search in the specified directory for all (KID,Pread)
% combinations and read the S21(TBB) sweep. It will then process all S21(F)
% sweeps to obtain Fres,Q and other KID specific parameters, as well as
% calibrate the resonance circle. After it uses blackbody_11.m to convert T into
% power @ detector front. Finally all information is stored to figure and csv files.
%
% NOTE: 2D (power and temperature) or 1D power data is handled by treating 
%   each (KID,power) combination as a new KID. Only the T dependence of
%   (KID,min(P)) is shown to screen unless 'showallplots' is set to 1 on
%   line 73. Further Power Dependence can be analysed using
%   PowerDependenceS21.m
%
% NOTE: Function assumes External subroutines are located in the 
%   ..\subroutines subdirectory, while this function is in the main
%   directory.
%
% NOTE: (TODO) These modifications could be made in the future:
%   *stopoffset variable for each resonator
%   *replace subplot(3,3,4) or (3,3,5) with a graph of Fres(T)
%
%INPUT:
%   SPECIFIED INTERNALLY
%
% S21subdir == subdirectory specifying the location of the S21 files.
%               Starts with "\", ends without "\"
% 
% ChipInfo == a struct containing information about the measured chip in 
%           the following elements
%   ChipInfo.path == Root path containing the measurement data.
% EXTRA: stopoffset, showallplots, TCfrac, FITF0. (See below).
%==========================================================================
% Set ChipInfo inside function.
%==========================================================================
%==========================================================================
%==========================================================================
%OUTPUT:
%   KID(n). == struct array of N = # KIDs, which for all 
%                     temperatures (nT = # temperatures).
%       KID(n).KIDnumber = KIDnumber (ID)
%       KID(n).filename = a single string containing the full
%                         path+filename that contains the S21datafile.
%       KID(n).Temperature(nT,1) = vector containing the chip temperature [K]
%       KID(n).ReadPower(nT,1) = vector containing the used read power [dBm]
%       KID(n).S21data{nT,1}(:,1:5) = cell array containing the S21data at each Temperature
%              Each cell contains a double array with the columns:
%              [F (GHz), |S21|, phase (rad), Re S21 (normalized), Im S21 (normalized)]
%       KID(n).Ql(nT,1) == Measured Loaded Q at each temperature [From S21 fit]
%       KID(n).Qi(nT,1) == Measured Internal Q at each temperature [From S21 fit]
%       KID(n).Qc(nT,1) == Measured Coupling (External) Q at each temperature [From S21 fit]
%       KID(n).Fres(nT,1) == Measured Resonance Frequency [From S21 fit, subsampling resolution]
%       KID(n).S21min(nT,1) == Minimum value S21 [From S21 fit, subsampling resolution]
%       KID(n).Paramerrors(nT,1:5) == errors in [Fres,S21min(dB),Ql,Qi,Qc] obtained from S21 fitting.
%       KID(n).InternalPower(nT,1) == KID internal power [dBm]
%       KID(n).Pbb    == BB power for each TBB
%       KID(n).ReImF0(:,1:2) == The value of S21(Re(1),Im(2)) in the calibrated
%                               resonance circle at Fres(T=0[K])
%       KID(n).Response(:,1:2) == The (change in) phase (2) and radius (1) at the 
%                                 base temperature resonance frequency wrt
%                                 the circle at T == 0 [K]
%       KID(n).responsivityM2(:,1:2) == phase (1) and radius (2) responsivity
%                                       as a function of temperature as 
%                                       determined by Method 2 (see above)
%   KIDParameters == Nx16 array contains the most important information 
%                  about all N=#(KID,P) combinations at base temperature
%                  It contains the parameters (Dim 2):
%               1 = KIDID
%               2 = Tbase [K]
%               3 = Fres(Tbase) [GHz]
%               4 = Fdesign [GHz]
%               5  = Q
%               6 = Qi
%               7 = Qc
%               8 = Qdesign
%               9 = Pread [dBm]
%               10 = Pint [dBm]
%               11 = ddx/dNqp
%               12 = dR/dNqp'
%               13 = Tc [K]
%               14 = Delta [J]
%               15 = Film Thickness [um]
%               16 = Active Area [um^2]
%
%   In addition all results are output as 4 types of files (Described below):
%   In the S21subdir where all measurements are:
%   One overview file with information of all KIDs at the base temperature.
%   One file per KID containing S21(F) at base temperature
%   One file per KID containing analysis results at all measurement temperatures.
%   In the measurement main directory:
%   A summary file containing 16 columns with the values of KIDParameters.
%
%
%SUBROUTINES:
% import_S21data
%   import_data(External)
% findparameters
% FitS21main5 (External)
% blackbody_11.m (External)
% findresponse
% savethedata
% 	WriteSRONcsv (External)
% plotTdependence
%
%REQUIRED FILES:
%
%VERSION: 4.3
%   Based upon the original "process_S21_2" program
%   V1.0 original script by Stephen Yates & Jochem Baselmans
%   V1.1(2011-03-16,JB) Added ALPascale option, using delta/lTc = 2.10
%   V1.2(2012-08-03,RJ) clean-up and addition of comments.
%                       Changed code (not effect) in which unique KIDnumbers
%                       are determined.
%   V2.0(2012-08-21,RJ) Addition of more comments, switch to FitS21_2
%                       routine,
%   V2.0(2012-09-05,RJ) Fixed a 2011b compatibility issue of matlabs unique
%                       function
%   V2.1 (2012-10-18,RJ) Made is mac/unix compatible by using filesep
%                       instead of /. Output Fres now in Hz (ipv GHz) for
%                       backwards compatibility.
%   V2.2 (2012-10-26,RJ) Updated some comments
%   V4.0 (2012-10-30,RJ) Modified the import_S21data function to use the
%                        general import_data subroutine. 
%   V4.1 (2012-11-07,RJ) Renamed to from HybridsS21 to RESPONSEanalysis.
%                       Added a clear statement to remove useless variables
%                       from the workspace before saving.
%   V4.2 (2013-02-14,RJ) Now uses FitS21main3
%   V4.3 (2013-03-06,RJ) Now made robust against poor data in which no KID
%                       can be found. (It will give warning but will not
%                       crash the code anymore).
%   V5.0 (2013-11-20,RJ) Now uses FitS21main5 en stores error estimates
%DATE: September 5, 2012
%AUTHOR: Reinier Janssen
%
%==========================================================================
% Layout of user INPUT defined files
%==========================================================================
% ResonatorInfoFile: A four column - tab separated - file without header
% lines containing for each KID the following information.
% ID    Fres[GHz]   Qdesign     Active Area [um^2]
% 17	4.2         14000   	42156
%
% KID S21 measurements (SRON): A single file per KIDnr-Power combination.
% Temperatures are added sequentially in these files.
%   EXAMPLE LAYOUT SRON S21 FILE
%C:\KID Metingen\ADR\B6Ch7__7_11_06 16_02\S21\2D\KID41_73dBm_JB1.dat
%Power at KID:73dBm
%resonance Frequency in GHz :4.893858
%Q=40792.348421, Qc=56509.887020, Qi=146662.340671  
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
% Layout of OUTPUT files
%==========================================================================
% KID_KIDid_PdBm_Tdep.csv: A 26 column - comma separated - file with header line
% containing for a KID at a single power the following parameters:
% 1=KIDID,    2=Temperature,	3=Q,	4=Qc,	5=Qi
% 6=Fres,   7=S21@Fres,	8=power,	9=Pinternal,	10=response
% 11=dtheta/dNqp FIT,	12=dtheta/dNqp numerical,     13=Nqp,	14=dx/dN,	15=V[um3]
% 16=Delta [J],	 17=Re point,	18=Im point,	19=drdnqp=ampresp,	20=dRdB
% 21=grad 1/Qi vs nqp, 22=Tc,   23=Fres des, 24=Q des, 25=Area des,
% 26=thickness
%
% (response = )
% (dRdB = numerical variation of radius with 1/Qi, normalized to lowest T)
%
% 
% KID_KIDid_PdBm_S21lowT.csv: A 5 column - comma separated - file with
% header line containing for a KID at a single power the transmission as a
% function of frequency at the lowest available temperature
% Freq[GHz],S21[dB],theta[rad],Re Norm, Im Norm
% 5.754142630000e+00,-3.573115700000e-02,-6.639363800000e-02,4.806554635962e-01,1.857603917372e-02
%
% allTdep.csv: A 26 column - comma separated - file containing the same
% parameters as KID_KIDid_PdBm_Tdep.csv. However, this file contains information
% of all analysed KIDs but only at base temperature.
%==========================================================================
%sub-dir inside ChipPath where the S21 T dependence is located. Starts with \ ends without \.
S21subdir='/S21/GOOD';                                          
if nargin==0
    ChipInfo.path = ...      
        [cd '/../']; %root path where data is, one higher than the scripts
end
%==========================================================================
% Internal variables which may be user modified.
%==========================================================================
global FITF0 InterpolatPhaseref
%BB settings
method.filter='850 GHz';    %Cardiff filter
method.tp='lambda^2';       %'lambda^2' or 'Geometrical'
method.eta_c = 1.0;         %optical efficiency CST
method.pol=1;               % one polarization
stopoffset=0;    %Removes the highest temperatures from the data analysis.
                  %Used to stop before the measurement has lost any KID due to dipjumping, extremely shallow dips, or similar.
showallplots = 0; %IF 1 the T dependences will be shown in plot for all (KID,P) combinations.
                  %IF 0 [default], only for the minimum read-out power the T dependence will be shown on screen.
                  %The other figures will be save and closed.
FITF0=0;          %Switch used to determine which method is used for determination of F0. 
                  %Please see the FitS21main5 routine for details.
% 0 == use minimum value of measured |S21| to determine fres and S21min.
%       Fit is only used to determine Q.
% 1 == use a logspace Lorentzian Fit to obtain Q, 
%       Fres and S21min. (FitS21_3 by Pieter de Visser). 
% 2 == use the function by Khalil et.al. to fit to |S21|(dB) (FitS21_2)
InterpolatPhaseref=1;    %1=find phaseref by interpolartion (default), best for VNA deep KIDs. 0 takes phase at measured datapoint closes to Fres.

%==========================================================================
% Setting some routine default values
%==========================================================================
global S21path
format('long','e'); %Set display format of numbers to 7 digits
%Enable subroutines by adding path in search path.
addpath([pwd,filesep,'subroutines']);
%close all; %close all open plots to remove clutter.
%==========================================================================
% In the datapath find all files and identify KIDs.
%==========================================================================
S21path=[ChipInfo.path S21subdir filesep] %Path containing raw S21 data
RawS21files=dir([S21path 'KID*.dat']); %Reading the names of all the raw S21 files into an array of cells
NumberOfS21Files = length(RawS21files);

KIDnumbers = zeros(1,NumberOfS21Files);
for n=1:NumberOfS21Files %Loop over all S21 files
    %Determine the KIDid (its number, XX) from each filename.
    KIDnumbers(n) = cell2mat(textscan(RawS21files(n).name,'%*3s %f %*s'));
end
%Determine all unique KIDs.
KIDnumbers = unique(KIDnumbers);

%Print to screen some of the reading information.
fprintf('Search for S21 data performed in:\n')
disp(S21path)
fprintf(['Inside the path the following KIDs are available: ',num2str(KIDnumbers),'\n'])
fprintf(['Inside this path a total number of ',num2str(NumberOfS21Files),' files were found.\n'])

%==========================================================================
% Analyse each S21 file. For each raw file a results file is
% created with a similar name to the S21data file, but with the
% ''extension.
% All results save in KID struct array.
% A graph is only given if a single power per KID is measured. Graphs are
% saved to disk in any case.
%==========================================================================

%Predefine some parameters:
transformation = cell(NumberOfS21Files,1);
%Stores for each (KID#,P) combination the calibration transformation of the 
%resonance circle in the complex plane
fitTbase = cell(NumberOfS21Files,1);
%Stores for each (KID#,P) combination the [F(GHz),|S21|(dB)] that is the
%evaluated version of the fit to the base temperature S21 transmission.
BBcal.T=0; 
for n=1:NumberOfS21Files %Loop over all S21 files
    %Extract from the filename piece of information: KIDnumber
    info=textscan(RawS21files(n).name,'%*3s %f %*1s %*f %*s'); %Get KIDnumber
    %Note: second float is the power [dBm]. Replace %*f by %f to get Pread = info{2}
    KID(n).KIDnumber = info{1}; %Storing KIDnumber
    KID(n).filename = [S21path,RawS21files(n).name]; %Storing original data file
    %======================================================================
    
    %Import the S21 data from the specified file.
    [KID(n).S21data,KID(n).Temperature,KID(n).ReadPower] = import_S21data(KID(n).filename,stopoffset);
    % Assigns:
    %   KID(n).S21data{:,1}(:,1:3) Frequency Magn (dB) Phase
    %   KID(n).Temperature(:,1)
    %   KID(n).ReadPower(:,1)
    %======================================================================
    
    %Analyse the S21data to obtain:
    [KID,transformation{n},fitTbase{n}]=findparameters(n,KID);
    % Requires: KID(n).Temperature, KID(n).S21data{:,1}(:,1:3), KID(n).ReadPower
    % Assigns: 
    %   KID(n).S21data{:,1}(:,2) [update, normalized & converted to magnitude space]
    %   KID(n).S21data{:,1}(:,4:5) 
    %   KID(n).S21data{:,1}  =  Freq magnitude(raw) phase(raw) Re(kidcircle) Im(Kidcircle)  
    %   KID(n).Ql(:,1)
    %   KID(n).Qi(:,1)
    %   KID(n).Qc(:,1)
    %   KID(n).Fres(:,1) [From S21 fit, subsampling resolution]
    %   KID(n).S21min(:,1) [From S21 fit, subsampling resolution]
    %   KID(n).Paramerrors(:,1:5)
    %   KID(n).InternalPower
    %   
    %   transformation{n} = 2 element vector containing the mean (over T) phase
    %   rotationa and Real axis translation applied to the resonance circle
    %   for calibration.
    %   fitTbase{n} = FitResult of the fit performed by FitS21 subroutine.
    %======================================================================  
    
    % Calculated the number of quasi-particles in the film at each temperature.
    % KID = QPnumber(n,KID,ChipInfo,designvalues);
    [KID(n).Pbb(:,1),~,BBcal,method]=blackbody_int(KID(n).Temperature(:,1),BBcal,method);
    % Requires: KID(n).Temperature, KID(n).KIDnumber
    % Assigns: 
    %   KID(n).Pbb(:,1)
    %SAVE the figure
Figfile=[S21path,'BBcallibration.fig'];
saveas(gcf,Figfile,'fig')
close(gcf)
    %======================================================================
    
    if length(KID(n).Temperature) < 4 %Value 4 is arbitrary. Can be anything >= 2
        %If the number of temperatures is too low reponse measurements are useless
        
        %Must make some assignments to avoid problems later.
        KID(n).ddxdNqp = zeros(1,2);
        KID(n).dBdNqp = zeros(1,2);
        KID(n).ResponsivityM1 = zeros(length(KID(n).Temperature),2);
        KID(n).ResponsivityM2 = zeros(length(KID(n).Temperature),2);
        KID(n).Response = zeros(length(KID(n).Temperature),2);
        KID(n).ReImF0 = zeros(length(KID(n).Temperature),2);
    else
        %If the number of temperatures is high enough reponse measurements
        %are useful.
        % Calculate the temperature respons of the KID in various ways.
        KID = findBBresponse(n,KID);
        % Requires: KID(n).Temperature, KID(n).Fres, KID(n).Pbb, KID(n).Ql,
        %           KID(n).Qi, KID(n).Qc,KID(n).S21data
        % Assigns:
        %   KID(n).ReImF0(:,1:2)
        %   KID(n).Response(:,1:2)
        %   KID(n).ResponsivityM2(:,1:2)
    end
    
    %Plot the data to screen at least once. Such that it will be saved for
    %later use. Close it to avoid clutter.
    plotTdependence(n,KID,fitTbase{n},transformation{n},-1*KID(n).KIDnumber)
    disp(['KID ' num2str(KID(n).KIDnumber) ' done'])
end
%==========================================================================
% Plot the data to screen. Crashes KIDParameters does not exist
%==========================================================================
% %If there are multiple powers per KID
% if showallplots
%     %User requests all plots are shown
%     for n=1:NumberOfS21Files %Loop over all S21 files (all KIDs)
%         plotTdependence(n,KID,fitTbase{n}(:,:),transformation{n}(:,:),n)
%     end
% elseif NumberOfS21Files == length(KIDnumbers)
%     %No mixing with Power sweep.
%     %Plot the results of all resonators to file and to screen to show user.
%     for n=1:NumberOfS21Files %Loop over all S21 files (all KIDs)
%         plotTdependence(n,KID,fitTbase{n}(:,:),transformation{n}(:,:),KID(n).KIDnumber)
%     end
%     
% else
%     save([S21path,'HybridsS21.mat'])
%     %Mixing with Power sweep. A slightly different method will be used to
%     %plot only the lowest power to screen for each resonator.
%     for n=1:length(KIDnumbers)
%         KIDid = KIDnumbers(n);
%         KIDlocations = find(KIDParameters(:,1) == KIDid);
%         [~,IndexPmin]=min(KIDParameters(KIDlocations,9));
%         IndexPmin = KIDlocations(IndexPmin);
%         plotTdependence(IndexPmin,KID,fitTbase{IndexPmin}(:,:),transformation{IndexPmin}(:,:),KID(n).KIDnumber)
%     end
% end

%==========================================================================
% Workspace saving and wrap-up
%==========================================================================
% clear up nonsense variables.
clear NumberOfS21Files RawS21files S21subdir ans fid info n showallplots
% Save the matlab workspace
save([S21path,'ResponseS21.mat'])
%Remove Path containing Subroutines from pathlist.
rmpath([pwd,filesep,'subroutines']);
%==========================================================================
%END OF HybridsS21_1RJ MAIN FUNCTION
%==========================================================================
end


function[powersum,NEP,BBcal,method]=blackbody_int(TTT,BBcal,method)
% calls external function blackbody_2.m if BBcal.T=0 (first KID in main
% program). backbody.m Calculates the powers vs black body T at many
% temperatures and puts the results in BBcal struct.
% afterwards, and also if BBcal.T~=0 the script inperpoaltes the BBcal
% struct ti get parameters.

% TTT col of temperatures
% BBcal struct with tempertaure, power and photon noise NEP
% pol polarisation 1 or 2
% method: string describing the optics of setup (see blackbody.m)
if isequal(BBcal.T,0) % if BB.cal is empty call blackbody.m
    BBcal.T=[3:0.05:100]';
    [BBcal.power,~,BBcal.NEP,method]=blackbody_14(BBcal.T,method,1); % old slow script
end

% interpolate BBcal struct, log space
powersum=10.^(interp1(log10(BBcal.T),log10(BBcal.power),log10(TTT)));

NEP.poisson=10.^(interp1(log10(BBcal.T),log10(BBcal.NEP.poisson),log10(TTT),'spline'));%2phF
NEP.g_r=10.^(interp1(log10(BBcal.T),log10(BBcal.NEP.g_r),log10(TTT),'spline'));%2PDelta/eta
NEP.wave=10.^(interp1(log10(BBcal.T),log10(BBcal.NEP.wave),log10(TTT),'spline'));%wave term
NEP.totphoton=10.^(interp1(log10(BBcal.T),log10(BBcal.NEP.totphoton),log10(TTT),'spline'));
end

function [data,Temperature,Power] = import_S21data(file,stopoffset)
%This function imports the S21 data from specified file. 
%From the main header only the power information is obtained.
%From the subheaders the temperature information is obtained.
%Always outputs the data, Temperatures and Power in order of unique
%monotonically increasing temperature
%
%INPUT:
% file == full file path
% stopoffset == number of temperatures at the end of the measurement disregarded
%OUTPUT:
% FOR N the number of temperatures, the following elements are read
% data == Nx1 cell array. Each cell contains a Mx3 array with [F (GHz), |S21| (dB), Phase (Rad)]
% Temperatures == Nx1 vector containing at the temperatures at which S21 is taken
% Power == Nx1 vector containing at the power (dBm) at which S21 is taken
%
%==========================================================================
% Example S21 data format inside the file.
%==========================================================================
%C:\KID Metingen\ADR\B6Ch7__7_11_06 16_02\S21\2D\KID41_73dBm_JB1.dat
%Power at KID:73dBm
%resonance Frequency in GHz :4.893858
%Q=40792.348421, Qc=56509.887020, Qi=146662.340671  
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
%==========================================================================

%Matlab textscan routine hard to understand (Reinier Janssen, 6-8-2012).
%Leave the code as is as it works. It just reads the data.

%==========================================================================
[tempdata,Temperature,Power,~] = import_data(file);

%Sort the temperatures and cut off the predefined number from the high T end.
%Making monotonic increasing Temperature in which all Temperature are
%unique. Non-unique values give problems later during interpolation steps.
[tempT,IncreasingTindex] = unique(Temperature(:,1),'first'); %'sorted' means picks the first of each encountered T duplicate.

%Number of temperatures corrected for the predefined high temperature cut-off
if stopoffset >= length(IncreasingTindex)
    %Would reduce the number of temperatures to zero or lower.
    %Give warning an ignore stopoffset
    fprintf('WARING [import_data]: Variable stopoffset is larger than the number of temperatures.\n')
    fprintf('Assuming stopoffset = 0 for this resonator.\n')
    Nbbtemperatures = length(IncreasingTindex);
else
    %Normal operation
    Nbbtemperatures = length(IncreasingTindex)-stopoffset; 
end

%Fill the output data in the correct way (monotically increasing T).
Temperature = tempT(1:Nbbtemperatures,1); %Temperature
data = tempdata(IncreasingTindex(1:Nbbtemperatures),1); %S21data
%General Matlab Note: tempdata and data are cell arrays usually to be
%addressed by {} brackets. However when copying cellcc to cell () brackets
%are to be used.
Power = -1*Power*ones(Nbbtemperatures,1); %Power

%Print information to screen
fprintf(['Selected temperatures range from ' num2str(Temperature(1,1)) ' to ' num2str(Temperature(end,1)) ' K\n']);
%==========================================================================
%END OF import_data SUBROUTINE
%==========================================================================
end

function [KID,transformation,fitTbase] = findparameters (n,KID)
%This function analyses the S21 data of all temperatures individually. For
%each temperature it does the following:
%1) Normalize |S21|(dB) to zero and convert it to magnitude space
%2) Fit the |S21|curve using one of the FitS21_x subroutine as selected by 
% FITF0 to obtain accurate values of Fres,Ql,Qc,Qi,|S21(Fres)|
%3) Calibrate the S21 resonance circle.

%INPUT:
% n == index in KIDstruct of KID considered.
% KID == large struct array containing all data.
%INPUT VIA GLOBAL:
% FITF0 == Switch index to determine the exact fitting routine used to 
% determine Fres and other parameters more accurately (To sub-frequency
% sampling resolution).
% The options for FITF0 are:
% 0 == use minimum value of measured |S21|. Q if fitted using FitS21_1
% 1 == use the original fit routine by Yates and Baselmans. (FitS21_1)
% 2 == use the function by Khalil et.al. to fit to |S21|(dB) (FitS21_2)
% 3 == [FUTURE]use the function by Khalil et.al. to fit to S21 (complex number, FitS21_3)
% 4 == [FUTURE] make a point-by-point fit in the complex plane. AND/OR make
%       use of theta(F) to determine bias point.
%
%==========================================================================
global FITF0 InterpolatPhaseref
Nbbtemperatures=length(KID(n).Temperature); %Find the number of temperatures.
S21data=KID(n).S21data; %Assign shorthand for S21data.

%Predefine a number of vectors with the number of temperatures
fres = zeros(Nbbtemperatures,1); % resonance frequency
Q = zeros(Nbbtemperatures,3); %[Ql,Qi,Qc]
S21min = zeros(Nbbtemperatures,1); %|S21(Fres)| (dB)
Paramerrors = zeros(Nbbtemperatures,5); %errors in [Fres,S21min(dB),Ql,Qi,Qc]

%==========================================================================
% Calibrate and convert |S21|
%==========================================================================
for p=1:Nbbtemperatures %Loop over all measured Temperatures
    %Smooth and normalize the |S21| data in dB space
    filtered=smooth(S21data{p}(:,2),3); %smoothing the |S21| data in dB space
    %normalise S21 in log space to the max(|S21|)
    S21data{p}(:,2)=S21data{p}(:,2)-max(filtered);
    
    %Convert to magnitude space
    S21data{p}(:,2)=10.^(S21data{p}(:,2)/20); %from log(magn) to magnitude
    %Division by 20 because dB is always defined in power. S21 is a
    %measurement in voltage. Since P \propto V^2: 
    %|S21| [dB] = 10*log10(|S21|^2) = 20*log10(|S21|)
    
    %Estimate the resonance parameters using fit and estimation routines.
    %Actual mode used is determined by FITF0 parameter.
    if p==1
        %At base temperature we want the evaluated fit function for future
        %plotting.
        [fres(p),Q(p,:),S21min(p),fitTbase,FitErrors] = FitS21main5(S21data{p}(:,1:3),FITF0);
    else
        [fres(p),Q(p,:),S21min(p),~,FitErrors] = FitS21main5(S21data{p}(:,1:3),FITF0);
    end
    Paramerrors(p,:) = FitErrors';
end

%SAVE THE ESTIMATES IN THE KID STRUCT
KID(n).Ql=Q(:,1);
KID(n).Qi=Q(:,2);
KID(n).Qc=Q(:,3);
KID(n).Fres=fres;
KID(n).S21min=S21min;
KID(n).S21data = S21data;
KID(n).Paramerrors = Paramerrors;
%Calculate the Internal Power
KID(n).InternalPower = 10*log10((2/pi)*10.^(KID(n).ReadPower/10).*(KID(n).Ql.^2./KID(n).Qc));

%==========================================================================
% Rotate the phase of the resonance circle such that:
% * Fres(T) is located in the direction of (-0.5,0)
% * Circle at T0 is centered at (0,0)
% * Points far off-resonance are, for all temperatures, located at (0.5,0)
%==========================================================================
%Determine the phase at resonance for each KID.
phaseres = zeros(Nbbtemperatures,1); %Phase at resonance
for p=1:Nbbtemperatures
    if InterpolatPhaseref==1
        %Interpolate Phase data to get on resonance phase.
        phaseres(p) = interp1(S21data{p}(:,1),unwrap(S21data{p}(:,3)),fres(p),'linear');
        
    else
        %ALTERNATIVE (ORIGINAL)
        [~,F0index] = min(S21data{p}(:,2));
        phaseres(p) = S21data{p}(F0index,3);
    end
    % Figure present for testing purposes
    if FITF0 == -1
        figure(KID(n).KIDnumber*100+p)
        clf
        hold on
        plot(S21data{p}(:,1),S21data{p}(:,3),'b-')
        plot(S21data{p}(:,1),unwrap(S21data{p}(:,3)),'k-')
        plot(fres(p),phaseres(p),'ro')
        hold off
    end
end
% Use the mean rotation to avoid problems for deepest KIDs
rotation = mean(phaseres);

%Use the fact the data is sorted that the lowest temperature is at index 1.
%Here the maximum value of Fres will be located. This is the resonance
%circle used as reference.

%rotating the circle
newphase=cell(Nbbtemperatures,1);
S21Real=cell(Nbbtemperatures,1);
S21Imag=cell(Nbbtemperatures,1);
OffsetReal=zeros(Nbbtemperatures,1);
for p=1:Nbbtemperatures
    newphase{p} = unwrap(S21data{p}(:,3)) - rotation; %Correct the phases of the S21 data
    S21Real{p} = S21data{p}(:,2).*cos(newphase{p}(:)); %Calculate Re{S21} with the calibrated phases
    S21Imag{p} = S21data{p}(:,2).*sin(newphase{p}(:)); %Calculate Im{S21} with the calibrated phases
    
     %Shift the S21 circle to center of the complex plane.
     OffsetReal(p) = (max(S21Real{p}(:))+10.^(S21min(1)/20))/2; %Q: Why not min(S21Real{1}(:)) ??
     %The use of S21min(1) means points far off resonance end up in the same
     %location (0.5,0.0)
     %ALTERNATIVE (Q: Why not this?):
     %OffsetReal(p)=(max(S21Real{p}(:))+min(S21Real{p}(:)))/2;
     %This means all circles are centered at (0,0) but not necessarily with
     %points far off-resonance in the same location (radius different, phase is equal)
     S21Real{p}(:)=S21Real{p}(:)-OffsetReal(p);
    
    %Update the S21data with calibrated S21 circle.
    S21data{p}(:,4)=S21Real{p}(:);
    S21data{p}(:,5)=S21Imag{p}(:);
end

%SAVE THE UPDATED S21data IN THE KID STRUCT
KID(n).S21data = S21data;

%Save the transformation for the use in plotting.
transformation(1)=rotation; %Mean phase correction [radians]
transformation(2)=mean(OffsetReal); %Mean offset correction.
%transformation(2)=OffsetReal(end); %Highest T used as offset correction.

%==========================================================================
%END OF findparameters SUBROUTINE
%==========================================================================
end

function  [KID] = findBBresponse(n,KID)
%This function takes the KID struct and using the information inside
%determines:
%   Response:
%   The (change in) phase and radius at the base temperature resonance
%   frequency is determined.
%
%   Responsivity Method 2:
%   Using change in the phase and radius at the base temperature resonance
%   frequency a numerical difference method is used to determine the
%   responsivity for thermal quasiparticles.
%
%NOTE: Makes use of the fact the lowest (base) temperature is always at index 1.
%And that from there the temperature is monotonically increasing.
%==========================================================================
global FITF0
Nbbtemperatures=size(KID(n).Temperature,1); %Number of temperatures / Length of most vectors

%Define Some Variables for easy use.
fres=KID(n).Fres; %resonant frequency
T=KID(n).Temperature; %Temperature
Pbb=KID(n).Pbb; %Quasiparticle Number
dx=KID(n).Fres/KID(n).Fres(1)-1; %dx = (F0(T)-F0(Tbase))/F0(Tbase)
Q=KID(n).Ql; %Q
Qc=KID(n).Qc; %Qc
Qi=KID(n).Qi; %Qi
S21data = KID(n).S21data; %S21data

%==========================================================================
% Phase Response && Radius Response
% As a function of temperature determine the value of R and Theta at
% frequency = f0(Tbase)
%==========================================================================
Response = zeros(Nbbtemperatures,8);
%dim 2:
%   1 == Re{S21(T)} @f0(T0) on callibrated MKID circle EXPORTED KID(n).ReImF0
%   2 == Im{S21(T)} @f0(T0) on callibrated MKID circle EXPORTED KID(n).ReImF0
%   3 == R{S21(T)} @f0(T0) on callibrated MKID circle
%   4 == Theta{S21(T)} @f0(T0) on callibrated MKID circle
%   5 == R{S21(T)} - R{S21(T0)} @f0(T0)   EXPORTED  KID(n).Response
%   6 == Theta{S21(T)} - Theta{S21(T)} @f0(T0)   EXPORTED KID(n).Response
%   7 == Smoothed version of (:,3) for responsivity calc
%   8 == Smoothed version of (:,4) for responsivity calc
S21Response = zeros(Nbbtemperatures,4);
%   1 == Magn(raw)(T) @f0(T0) EXPORTED KID(n).S21Response
%   2 == Phase(raw)(T) @f0(T0) EXPORTED KID(n).S21Response
%   1 == Magn(raw)(T) @f0(T0) - Magn(raw)(T=0) @f0(T0)EXPORTED  KID(n).S21Response
%   2 == Phase(raw)(T) @f0(T0) - Phase(raw)(T=0) @f0(T0) EXPORTED  KID(n).S21Response

for p=1:Nbbtemperatures
    %Use linear interpolation between the data points to obtain the value at f0(T0)
    S21Response(p,1) = interp1(S21data{p}(:,1),S21data{p}(:,2),fres(1),'spline');%magn raw
    S21Response(p,2) = interp1(S21data{p}(:,1),S21data{p}(:,3),fres(1),'spline');%phase raw
    Response(p,1) = interp1(S21data{p}(:,1),S21data{p}(:,4),fres(1),'spline'); %Real in calibrated circle
    Response(p,2) = interp1(S21data{p}(:,1),S21data{p}(:,5),fres(1),'spline'); %Imag in calibrated circle
end
Response(:,3) = abs(Response(:,1)+1i*Response(:,2)); %Radius
Response(:,4) = atan2(Response(:,2),Response(:,1)); %Theta [rad]
%Correction to remove 2pi phase jumps
Response(:,4) = unwrap(Response(:,4));
%OLD METHOD of removing jumps
%PhaseJumps = Response(:,4) <-0.75*pi;
%Response(PhaseJumps,4) = Response(PhaseJumps,4)+2*pi;

%Calculate the response to temperature (change in R and Theta)
Response(:,5) = Response(:,3) - Response(1,3); %R{S21(T)} - R{S21(T0)} @f0(T0)
Response(:,6) = Response(:,4) - Response(1,4); %Theta{S21(T)} - Theta{S21(T)} @f0(T0)
Response(:,7) = Response(:,5)/(1-10^(KID(n).S21min(1)/20))/2; %R{S21(T)} - R{S21(T0)} @f0(T0) normalised 
S21Response(:,3) = S21Response(:,1) - S21Response(1,1); %
S21Response(:,4) = S21Response(:,2) - S21Response(1,2); %


% Figure present for testing purposes
if FITF0 == -1
    figure(KID(n).KIDnumber*100)
    clf
    hold on
    plot(T,Response(:,4),'b-')
    plot(T,Response(:,6),'b--')
    xlabel('T [K]')
    ylabel('Phase [rad]')
    legend('\theta(S21(T)) @F0(T0)','\delta\theta(S21(T)) @F0(T0)')
    hold off
end

%Save the interesting values to the KID struct
KID(n).ReImF0 = Response(:,1:2);    %Real and Im point on callibrated circle @ Fres(T0)
KID(n).Response = Response(:,5:7);
KID(n).S21Response = S21Response;

%==========================================================================
% Method 2: Direct numerical differentiation (after smoothing) to obtain responsivity.
% Phase Responsivity: [Theta(Ti+1) - Theta(Ti)]/[Nqp(Ti+1) - Nqp(Ti)]
%   (Suffers from backbending due to TLS and the fact that at high T you
%   reach the other side of the circle)
% Amplitude Responsivity: [R(Ti+1) - R(Ti)]/[Nqp(Ti+1) - Nqp(Ti)]
%==========================================================================

%Smooth the responses for a nicer responsivity
Response(:,7) = smooth(Response(:,3),3); %Smooth R
Response(:,8) = smooth(Response(:,4),3); %Smooth Theta

%Determine the responsivity by numerical differentiation.
ResponsivityM2 = zeros(Nbbtemperatures,2);
ResponsivityM2(2:end,1) = (Response(2:end,8)-Response(1:end-1,8))./(Pbb(2:end)-Pbb(1:end-1));%Phase Responsivity.
ResponsivityM2(2:end,2) = (Response(2:end,7)-Response(1:end-1,7))./(Pbb(2:end)-Pbb(1:end-1));%Radius Responsivity.

%Save the method 2 responsivity to KID struct
KID(n).ResponsivityM2 = ResponsivityM2;

%==========================================================================
%END OF findBBresponse SUBROUTINE
%==========================================================================
end
    
function plotTdependence(n,KID,fitTbase,transformation,figurehandle)
% This function takes KID struct array and plots all the available
% information inside. This should give insight in the analysis done by
% other subroutines in the HybridS21 analysis program.
%
% INPUT
% n == number of the KID in the KID struct
% KID == KID struct containing all the data
% fitTbase == Nx2 array [F,|S21|(dB)] that is the evaluated function of the
%             fit result given by the FitS21_X routine.
% transformation == mean (in T) transformation applied to calibrate the
%                   resonance circle
% figurehandle == absolute value of this number will become the number of
%                   the figure. If negative, the figure will be closed
%                   after saving. 0 is not allowed.
%
%==========================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%modified slightly to show three methods of responsivity 20-2-2007
%Modified aug 28 2006 JochemB; minor cosmetic changes, 
%addition of saturation point in 3.3.6; resp. points in 3.3.3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set format for printing to screen
format('long','e');
%Generate a color rainbow for temperature sweeps
Nbbtemperatures = length(KID(n).Temperature);
[Tcolors,~] = GenerateColorsFromMap(Nbbtemperatures,'RainbowReinier');

%==========================================================================
% Figure == KIDnumber: Overview of all temperature dependencies of the KID
%==========================================================================
figure(abs(figurehandle));
clf %clear current figure

%SUBPLOT 3.3.1 KID RESONANCE at base temperature
subplot(3,3,1)
hold on
%Plot the data, fit and resonance frequency
plot(KID(n).S21data{1}(:,1),KID(n).S21data{1}(:,2),'bx','MarkerSize',5) %The measured resonance
plot(fitTbase(:,1),10.^(fitTbase(:,2)/20),'k-','LineWidth',1) %Plot the base Temperature fit
plot(KID(n).Fres(1),10^(KID(n).S21min(1)/20),'kd',...
    'MarkerFaceColor','r','MarkerSize',6)%The determined Fres(0) and |S21(Fres(0))|
%Make a nice layout.
xlim([min(KID(n).S21data{1}(:,1)),max(KID(n).S21data{1}(:,1))])
ylim([0,1.1]) %Due to being in magnitude 0<|S21|<= 1
legend('hide') %Otherwise an automatic legend pops up from the plot fitresult combination.
xlabel('F [GHz]')
ylabel('|S21| [mag]')
title(['KID ',num2str(KID(n).KIDnumber),' @ ',num2str(KID(n).ReadPower(1)),' [dBm]'])
hold off

%SUBPLOT 3.3.2 KID RESONANCES as a function of temperature
subplot(3,3,2)
hold on
for p=1:Nbbtemperatures
    plot(KID(n).S21data{p}(:,1),20*log10(KID(n).S21data{p}(:,2)),'.-',...
        'color',Tcolors(p,:),'MarkerSize',2) %Plot the resonances in dB space
end
%Two stage loop so that resonance freq dots are on top.
for p=1:Nbbtemperatures
    plot(KID(n).Fres(p),KID(n).S21min(p),'ko','MarkerSize',4,'MarkerFaceColor','k') %Plot Fres
end
%Make a nice layout
xlabel('F [GHz]')
ylabel('S21 [dB]')
title('KID resonances')
axis tight;
hold off

%SUBPLOT 3.3.3 Q FACTORS
subplot(3,3,3)
semilogy(KID(n).Temperature,KID(n).Ql,'ko-',... %Loaded Q
    'MarkerSize',4,'LineWidth',1)
hold on %Hold on after first semilogy to keep y-axis in log scale
semilogy(KID(n).Temperature,KID(n).Qc,'ro-',... %Coupling Q
    'MarkerSize',4,'LineWidth',1)
semilogy(KID(n).Temperature,KID(n).Qi,'bo-',... %Internal Q
    'MarkerSize',4,'LineWidth',1)
%Make a nice layout
legend('Q_l','Q_c','Q_i')
xlabel('T [K]')
ylabel('Q')
title('Q factors')
hold off

%SUBPLOT 3.3.4 KID RESONANCES in complex plane (As Measured)
subplot(3,3,4)
hold on
for p=1:Nbbtemperatures
    %Plot the resonance as measured
    MeasuredRe = KID(n).S21data{p}(:,2) .* cos(KID(n).S21data{p}(:,3));
    MeasuredIm = KID(n).S21data{p}(:,2) .* sin(KID(n).S21data{p}(:,3));
    plot(MeasuredRe,MeasuredIm,'o-',...
        'color',Tcolors(p,:),'MarkerSize',2);
end
%Two stage loop to have points on top.
for p=1:Nbbtemperatures
    %Recalculate the resonance as measured
    MeasuredRe = KID(n).S21data{p}(:,2) .* cos(KID(n).S21data{p}(:,3));
    MeasuredIm = KID(n).S21data{p}(:,2) .* sin(KID(n).S21data{p}(:,3));
    %Determine the measurement point closest to resonance
    [~,F0index] = min(abs(KID(n).Fres(p) - KID(n).S21data{p}(:,1))); 
    plot(MeasuredRe(F0index),MeasuredIm(F0index),'ko',...
        'MarkerSize',5,'MarkerFaceColor','k')
    %Apply the inverse transformation to the (Re,Im)@F0(T0) value.
    RealF0 = KID(n).ReImF0(p,1)+transformation(2); %Reverse transformation on Real axies
    Phase = atan2(KID(n).ReImF0(p,2),RealF0)+transformation(1); %Calculate Phase & Reverse transformation theta
    Radius = abs(RealF0+i*KID(n).ReImF0(p,2)); %Calculate new radius.
    %Plot the trace of F0(T0) in the complex plane
    plot(Radius*cos(Phase),Radius*sin(Phase),'kd',...
        'MarkerSize',5,'MarkerFaceColor','r')
end
%Make a nice layout
grid on
xlabel('Re')
ylabel('Im')
title('Measured KID resonances in complex plane')
axis([-1 1 -1 1])
hold off

%SUBPLOT 3.3.5 KID RESONANCES in complex plane (Calibrated)
subplot(3,3,5)
hold on
for p=1:Nbbtemperatures
    %Plot the calibrated resonance circles
    plot(KID(n).S21data{p}(:,4),KID(n).S21data{p}(:,5),'.-',...
        'color',Tcolors(p,:),'MarkerSize',1); %Calibrated Measurent in Circle
end
for p=1:Nbbtemperatures
    %Determine the measurement point closest to resonance
    [~,F0index] = min(abs(KID(n).Fres(p) - KID(n).S21data{p}(:,1))); 
    plot(KID(n).S21data{p}(F0index,4),KID(n).S21data{p}(F0index,5),'ko','MarkerSize',5,'MarkerFaceColor','k'); %Location Fres(T)
    %Plot the F0(T0) trace in the complex plane
    plot(KID(n).ReImF0(p,1),KID(n).ReImF0(p,2),'kd','MarkerSize',5,'MarkerFaceColor','r'); %Location Fres(Tbase)
end
%Make a nice layout
grid on
xlabel('Re')
ylabel('Im')
title('Calibrated KID resonances in complex plane')
axis([-1 1 -1 1]*0.5)
hold off

%SUBPLOT 3.3.7 Frequency versus Pbb
subplot(3,3,6)
hold on
DX = (KID(n).Fres - KID(n).Fres(1))/KID(n).Fres(1);
plot(KID(n).Pbb*1e12,DX,'k.-',... %(F-F0)/F0 as function of Nqp
    'MarkerSize',4,'LineWidth',1)
%Make a nice layout
xlabel('Pbb (pW)')
ylabel('dx=(F-F0)/F0')
title('Frequency vs Nqp');

%SUBPLOT 3.3.8 KID RESPONSE
subplot(3,3,7)
plot(KID(n).Pbb*1e12,KID(n).Response(:,2),'-or');hold on; %Phase response
plot(KID(n).Pbb*1e12,10*KID(n).Response(:,3),'-ob'); %Normaised Amplitude response
%General Stuff
xlabel('Pbb [pW]');ylabel('Response');
legend('Phase angle (Rad)','10*(1-Radius/,max(Radius))')
title('KID response')
hold off




%SAVE the figure
Figfile=[KID(n).filename(1:end-3),'fig'];
saveas(gcf,Figfile,'fig')

%Clean up if desired.
if figurehandle > 0
    %keep plot on the screen
else
    %close the figure
    close(abs(figurehandle));
end

%==========================================================================
%END OF plotTdependence SUBROUTINE
%==========================================================================
end