function [NOISE]= NOISEcompareV8

% This function can be run after NOISEanalysis has finished analysing the
% data. It load the NOISEanalysis workspace from file.
%
%INPUT: All user defined in this function (see below)
%
%OUTPUT:
%
% note on output
% FitSf{p,1}.b has the Gao F noise at 1 Hz and -40 dBm power at temperatue
% set in the script
% FitSf{p,2}.b has the Gao F noise at 10 Hz and -40 dBm power at temperatue
% set in the script
% FitSf{p,3}.b has the Gao F noise at 100 Hz and -40 dBm power at temperatue
% set in the script
% FitSf{p,4}.b has the Gao F noise at 1000 Hz and -40 dBm power at temperatue
% set in the script
%
% NOTE: NOISE is updated with offres data
%
% The script uses an array with length KIDnumbers that has all the data.
% KID numbers (with typically index kidn to adress one) has for each KID,
% the KID ID, readout power and Temperature to use. (T is set in this script, P=Popt from Pop.csv file)
% Lookuptables is a cell array length KIDnumbers with in each sell a px4 array with
% [T value,T index,P value,P index], all T and P=Ppt we want to use.
% p=index that is length(NOISE), i.e. all P and T combinations long
%
% later in the plots the script uses PoptT0index, KIDnumbers liong
% that has for each kidn the indices PoptT0index(kidn,1) for P and
% PoptT0index(kidn,2) for T for the desired T and P to use
%
%SUBROUTINES:
% GenerateColorsFromMap (External)
% 	colormapstorage.mat (External)
%
%REQUIRED FILES:
% Noise.mat (Produced by NOISEanalysis.m)
% Popt.csv (optional, produced by NOISEanalysis.m)
%
%VERSION: 5.0
%   Contains ideas of the original NoiseCompare V3 by JB
%   V4.0 (2012/11/16, RJ): Start writing of this script.
%   V5.0 (18-12-2012, JB): small changes and consistency check
%   V7.0 (2012-12-xx, JB): Temperature under investigation can now be
%                           selected. New NOISEparameters is created at
%                           desired temperature. Minor changes in what is
%                           saved.
%DATE: November 16, 2012
%AUTHOR: Reinier Janssen
% V8 Jochem 2014-3-3)   Displays also the noise wrt the complex plane
%==========================================================================
% User defined input
%==========================================================================
% Full directory containing the data and analysis output that has been
% processed by NOISEanalysis. It is the combination of the following
% variables in the NOISEanalysis space: Noisepath = [ChipInfo.path,FFTsubdir,'\']
Path = [cd '/..']; %root path where data is, one higher than the scripts
FFT_subdir='/FFT/Power_test';
FFTOffRes_subdir='';  % local path with ofres noisedata
% Note that this data has not been analyzed.%
BADKIDS = []; %Numbers of MKIDs to be excluded of plots that use the fit of noise vs Pi

Noisepath = [Path,FFT_subdir,filesep];
OffresNoisepath = [Path,FFTOffRes_subdir,filesep];%Leave empty if not available
T_to_use=0.10; %in K; put in temperature you want to analyze. for only one T taken will allways work

%System Noise Temperature
Tsystem = 7; %K (Default values: ADR = 7K)

%If the chip under analysis is the Akira Test Chip with varying Al length
%please set this to 1. It will convert KID number (KID ID) to aluminum
%percentage (% of length).
Akiratestchip=0;                        %set to 1 if you want to convert KID ID to %Al
if Akiratestchip
    %alperc=[0 0 20 20 20 20 40 40 60 60 60 60 80 80 100 100 100 100 110 110]; %design data Akira test chip
    alperc=[0 0 4.14 4.26 11.2 11.6 18 18.6 31.4 32.7 64.4 67.4 100 100 110 110 ]; %design data Akira test chip
    
end

%==========================================================================
%Standard settings to determine when to give warnings. Can be user modified
%==========================================================================

%==========================================================================
% Setting some routine default values
%==========================================================================
format('long','g'); %Set display format of numbers to 7 digits
%Enable subroutines by adding path in search path.
addpath([pwd,filesep,'subroutines']);
close all; %close all open plots to remove clutter.

%==========================================================================
% Load the workspace of HybridS21_2
%==========================================================================
NoiseVariables = {'NOISE','NOISEParameters','KIDnumbers','ChipInfo','IndexPsort'};
load([Noisepath,'Noise.mat'],NoiseVariables{:})
%This matlab workspace storage should at least contain the variables:
% NOISE, NOISEParameters,KIDnumbers,Noisepath (overwrites current)
% PoptData (will be overwritten), ChipInfo, IndexPsort

disp(['Noise will be compared for KIDS: ',num2str(KIDnumbers)])
%==========================================================================
% Read the Popt file (Or use the data found in the environment of
% NOISEanalysis)
%==========================================================================
PoptFile = [Noisepath,'Popt.csv'];
fid = fopen(PoptFile);
if fid == -1
    disp('Warning NOISEcompare: Cannot find Popt.csv. Using PoptData from NOISEanalysis environment instead.')
    load([Noisepath,'Noise.mat'],'PoptData')
else
    fclose(fid);
    %There is a Popt file. Read it and use it.
    fprintf('Found Popt.csv, this will be used to indicate Popt.\n')
    [~,PoptData] = ReadSRONcsvV2(PoptFile,'',0);
    %PoptData will be Nx3. Where dim 2 contains [T,KIDid,|Popt [dBm]|]
end
%JB12-1-13 now create NOISEparameters at the desired temperature
for p=1:length(NOISE)
    [~,BaseTindex] = min(abs(T_to_use-NOISE(p).Temperature));
    
    NOISEParameters(p,1) = NOISE(p).KIDnumber; %KID ID
    NOISEParameters(p,2) = NOISE(p).Temperature(BaseTindex); %correct temperature
    NOISEParameters(p,3) = NOISE(p).Fres(BaseTindex); %Fres
    NOISEParameters(p,5) = NOISE(p).Ql(BaseTindex); %Q
    NOISEParameters(p,6) = NOISE(p).Qi(BaseTindex); %Qi
    NOISEParameters(p,7) = NOISE(p).Qc(BaseTindex); %Qc
    NOISEParameters(p,9) = NOISE(p).ReadPower(BaseTindex); %Pread
    NOISEParameters(p,10) = NOISE(p).InternalPower(BaseTindex); %Pint
    NOISEParameters(p,11:14) = NOISE(p).MeanFreqNoise(BaseTindex,:); %Frequency Noise
    
end


%Convert the PoptData back into the LookUpTables (PoptData sorted into a
%cell with a cell element per KID).
LookUpTables = cell(length(KIDnumbers),1);
for kidn=1:length(KIDnumbers)
    Match = PoptData(:,2) == KIDnumbers(kidn);
    LookUpTables{kidn,1} = zeros(sum(Match),4); %will be [T value,T index,Popt value,Popt index]. index=1foer the T or P top use
    LookUpTables{kidn,1}(:,1) = PoptData(Match,1);
    LookUpTables{kidn,1}(:,3) = PoptData(Match,3);
    %Find the indexes corresponding to these powers and temperatures
    for n=1:sum(Match)
        Pindex = find(NOISEParameters(:,1) == KIDnumbers(kidn) & ...
            NOISEParameters(:,9) == LookUpTables{kidn,1}(n,3));
        LookUpTables{kidn,1}(n,4) = Pindex;
        [~,Tindex] = min(abs(NOISE(Pindex).Temperature(:) - LookUpTables{kidn,1}(n,1)));
        LookUpTables{kidn,1}(n,2) = Tindex;
    end
end
clear Match Pindex Tindex
Offres=zeros(1,length(NOISE))+10;%if equals 1 we have offres noise, initialise with 10 to easy spot errors
%==========================================================================
% Import system noise measured off resonance if available
%==========================================================================
if isempty(OffresNoisepath)
    %No off resonance noise data available
    disp('No off resonance noise measurements found.')
    for p=1:length(NOISE)
        NOISE(p).Offresfile = '';
        NOISE(p).OffresFFTnoise = cell(size(NOISE(p).Temperature));
        NOISE(p).OffresFread = [];
        NOISE(p).OffresTDIQ = cell(size(NOISE(p).Temperature));
        Offres(p)=0;
    end
else
    for p=1:length(NOISE) %NB: NOISE goes through all KIDs and all noise powers
        %Check each (KID,P) combination to see if there are off resonance noise
        %files.
        OffresFFT = dir([OffresNoisepath,'KID',num2str(NOISE(p).KIDnumber),...
            '_',num2str(abs(NOISE(p).ReadPower(1,1))),'dBm_*_FFT.dat']);
        %OffresTD = dir([OffresNoisepath,'KID',num2str(NOISE(p).KIDnumber),...
        %   '_',num2str(abs(NOISE(p).Pread(1,1))),'dBm_*_td.dat']);
        
        if isempty(OffresFFT)
            %disp(['No off resonance noise measurements found for KID ',num2str(NOISE(p).KIDnumber),' at ',num2str(NOISE(p).ReadPower(1,1)),' dBm']);
            Offres(p)=0;
            
        else
            if length(OffresFFT) >= 2
                disp(['Warning: Multiple off resonance noise data files found. Using first only.'])
            end
            disp(['Found off resonance noise measurements for KID ',num2str(NOISE(p).KIDnumber),...
                ' at ',num2str(NOISE(p).ReadPower(1,1)),' dBm'])
            Offres(p)=1;
            
            %Reading Off resonance FFT file
            [Data,Temperature,Power,FFTheader] = import_data([OffresNoisepath,OffresFFT(1).name]);
            
            %Check for equality of power
            if abs(-1*Power-NOISE(p).ReadPower(1)) > 0.5
                disp('Warning: Read Power does not match between FFT on and off resonance.')
            end
            
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
            
            %Check matching temperatures and copy variables into the right
            %places
            NOISE(p).OffresFFTnoise = cell(length(NOISE(p).Temperature),1);
            NOISE(p).OffresFread = zeros(length(NOISE(p).Temperature),1);
            for nT=1:length(Temperature)
                [Tdiff,Tindex] = min(abs(NOISE(p).Temperature - Temperature(nT)));
                NOISE(p).OffresFFTnoise(Tindex,1) = Data(nT,1);
                if Tdiff > 10e-3
                    disp('Warning: More that 10 mK difference between on and off resonance noise data')
                end
                NOISE(p).OffresFread(Tindex,1) = Fread;
            end
        end %END OF IF NO FILES FOUND
    end %END OF LOOP OVER ALL KID,POWER COMBINATIONS
end

%==========================================================================
% Calculate setup noise
%==========================================================================
for p=1:length(NOISE) %LOOP OVER ALL KID,POWER COMBINATIONS
    S21minMag = 10.^(NOISE(p).S21min/20);
    NOISE(p).SystemNoiseT = Tsystem; %Store system noise Temperature
    NOISE(p).SetupNoise = 10*log10(1.38e-23*Tsystem*(2./(1-S21minMag)).^2)-(NOISE(p).ReadPower-30);
end %END OF LOOP OVER ALL KID,POWER COMBINATIONS

%==========================================================================
% Determine frequency noise at -40 dBm
%==========================================================================
PoptT0index = zeros(length(KIDnumbers),2); %[Pindex,T0index]
PsubPopt = cell(length(KIDnumbers),1);
FitSf = cell(length(KIDnumbers),4);
%Stores the fit to the frequency noise as a function of internal power for
%1 Hz, 10 Hz, 100 Hz and 1 kHz

NOISEParameters(p,11:14) = NOISE(p).MeanFreqNoise(BaseTindex,:); %Frequency Noise
%Why ???? RJ (2013-02-14)

for kidn=1:length(KIDnumbers)
    %For each KID we determine the frequency noise at -40 dBm of internal
    %kidnower using a linear fit to all powers equal to or below Popt.
    % and at the desired temperature
    
    [~,T0LUT] = min(abs(T_to_use-LookUpTables{kidn,1}(:,1)));%find desired T to use
    PoptT0index(kidn,1) = LookUpTables{kidn,1}(T0LUT,4);% correct power index for kidn in the NOISE array.
    PoptT0index(kidn,2) = LookUpTables{kidn,1}(T0LUT,2);% correct T index
    
    TempIndex = NOISEParameters(IndexPsort{kidn,1}(:,1),9) <= NOISEParameters(PoptT0index(kidn,1),9);
    PsubPopt{kidn,1} = IndexPsort{kidn,1}(TempIndex,1);
    
    %Put the Sf/F^2 fits in nice arays
    
    FitSf{kidn,1} = NOISE(PoptT0index(kidn,1)).FitS1Hz(PoptT0index(kidn,2));
    FitSf{kidn,2} = NOISE(PoptT0index(kidn,1)).FitS10Hz(PoptT0index(kidn,2));
    FitSf{kidn,3} = NOISE(PoptT0index(kidn,1)).FitS100Hz(PoptT0index(kidn,2));
    FitSf{kidn,4} = NOISE(PoptT0index(kidn,1)).FitS1kHz(PoptT0index(kidn,2));
    
    
end %END OF LOOP OVER ALL KIDS
%==========================================================================
% Create plots
%==========================================================================

%Generate Colors for different KIDS
[KIDcolors,~] = GenerateColorsFromMap(length(KIDnumbers),'RainbowReinier');

%Generate KID ID (here Al% can be taken into account)
if Akiratestchip
    ID = alperc(NOISEParameters(PoptT0index(:,1),1));
else
    ID = NOISEParameters(PoptT0index(:,1),1);
end

%Generate the colormap to view Q distribution
QCmap = colormap(jet(41));
minQ = floor(min(NOISEParameters(PoptT0index(:,1),5))/1000);
Qrange = ceil(max(NOISEParameters(PoptT0index(:,1),5))/1000)-minQ;
QCindex = floor(((NOISEParameters(PoptT0index(:,1),5)/1000-minQ)/Qrange)*(length(QCmap)-1))+1;
QCtick = [1; 11; 21; 31; 41];
tick2plot = round(((QCtick-1)*Qrange+(length(QCmap)-1)*minQ)/(length(QCmap)-1));

Py = -70:-20; %Internal powers over which fits are evaluated.

% [~,Tindex] = min(abs(NOISE(Pindex).Temperature(:) - LookUpTables{p,1}(n,1)));
%       LookUpTables{p,1}(n,2) = Tindex;

figure(1)

subplot(2,3,1) %At Popt Amp,Phase and Off resonance Noise
offreslegend=0;
for kidn=1:length(KIDnumbers)
    semilogx(NOISE(PoptT0index(kidn,1)).FFTnoise{PoptT0index(kidn,2),1}(:,1),...
        NOISE(PoptT0index(kidn,1)).FFTnoise{PoptT0index(kidn,2),1}(:,2),...
        '-','color',KIDcolors(kidn,:),'LineWidth',2) %Phase Noise
    hold on
    semilogx(NOISE(PoptT0index(kidn,1)).FFTnoise{PoptT0index(kidn,2),1}(:,1),...
        NOISE(PoptT0index(kidn,1)).FFTnoise{PoptT0index(kidn,2),1}(:,3),...
        '--','color',KIDcolors(kidn,:),'LineWidth',2) %Amp Noise
    
    %if %isempty(NOISE(PoptT0index(kidn,1)).OffresFread)NOISE(kidn).OffresFread(Tindex,1)
    if   Offres(PoptT0index(kidn,1))==1
        %Off resonance noise known at optimum power
        offreslegend=1;
        semilogx(NOISE(PoptT0index(kidn,1)).OffresFFTnoise{PoptT0index(kidn,2),1}(:,1),...
            NOISE(PoptT0index(kidn,1)).OffresFFTnoise{PoptT0index(kidn,2),1}(:,2),...
            '-.','color',KIDcolors(kidn,:)) %Phase Noise
        semilogx(NOISE(PoptT0index(kidn,1)).OffresFFTnoise{PoptT0index(kidn,2),1}(:,1),...
            NOISE(PoptT0index(kidn,1)).OffresFFTnoise{PoptT0index(kidn,2),1}(:,3),...
            ':','color',KIDcolors(kidn,:)) %Amp Noise
    end
    
end
xlim([1 1e6])
xlabel('F [Hz]')
ylabel('S_x')
title(['Noise at P_{opt} and ', num2str(NOISE(PoptT0index(1,1)).Temperature(PoptT0index(1,2)),2) ' K'])
if offreslegend==0
    legend('S_{\theta} On','S_A On')
else
    legend('S_{\theta} On','S_A On','S_{\theta} Off','S_A Off')
end
hold off


subplot(2,3,2) %Frequency noise at Popt
for kidn=1:length(KIDnumbers)
    semilogx(NOISE(PoptT0index(kidn,1)).FFTnoise{PoptT0index(kidn,2),1}(:,1),...
        10*log10(NOISE(PoptT0index(kidn,1)).FFTnoise{PoptT0index(kidn,2),1}(:,4)),...
        '-','color',KIDcolors(kidn,:),'LineWidth',2) %Frequency Noise
    hold on
end
xlim([1 1e6])
ylim([-210 -150])
xlabel('F [Hz]')
ylabel('S_{f}/f^2')
legend(num2str(KIDnumbers'))
title(['Frequency Noise at P_{opt} and ', num2str(NOISE(PoptT0index(1,1)).Temperature(PoptT0index(1,2)),2) ' K'])
grid on
hold off

subplot(2,3,3) %System noise
for kidn=1:length(KIDnumbers)
    semilogx(NOISE(PoptT0index(kidn,1)).FFTnoise{PoptT0index(kidn,2),1}(:,1),...
        NOISE(PoptT0index(kidn,1)).FFTnoise{PoptT0index(kidn,2),1}(:,3),...
        '--','color',KIDcolors(kidn,:),'LineWidth',2) %Amp Noise
    hold on
    if Offres(PoptT0index(kidn,1))==1
        %off resonance noise known
        
        semilogx(NOISE(PoptT0index(kidn,1)).OffresFFTnoise{PoptT0index(kidn,2),1}(:,1),...
            NOISE(PoptT0index(kidn,1)).OffresFFTnoise{PoptT0index(kidn,2),1}(:,3),...
            ':','color',KIDcolors(kidn,:)) %Amp Noise
    end
    semilogx(1e5,...
        NOISE(PoptT0index(kidn,1)).SetupNoise(PoptT0index(kidn,2),1),...
        'ko','MarkerFaceColor',KIDcolors(kidn,:),'LineWidth',2,'MarkerSize',8)%Setup Noise
end
xlim([1 1e6])
xlabel('F [Hz]')
ylabel('S_A')
title(['System Noise at P_{opt}. T_{sys} = ',num2str(Tsystem),' K'])
if offreslegend==0;
    legend('S_A On','Setup @ T_{sys}')
else
    legend('S_A On','S_A Off','Setup @ T_{sys}')
end
hold off

subplot(2,3,4) %Freq Noise as a function of power (1 kHz)
for kidn=1:length(KIDnumbers)%NOISE(IndexPsort{p,1}(nP,1)).MeanFreqNoise(PTindex(nP,:,2,1),4);%1kHz
    if isempty(find(kidn==BADKIDS,1))
        plot(NOISEParameters(PsubPopt{kidn,1}(:),10),10*log10(NOISEParameters(PsubPopt{kidn,1}(:),14)),...
            'o','Color',KIDcolors(kidn,:),'MarkerFaceColor',KIDcolors(kidn,:),'MarkerSize',8) %Measurements
        hold on
        plot(FitSf{kidn,4}{1},'k') %Fit FitSf{kidn,4}
    end
end
plot(Py,-189-0.5*(Py+40),'k-','LineWidth',2) %Gao Line
axis tight
legend('off')
xlabel('P_{int} [dBm]')
ylabel('S_F/F^2(1 kHz)')
hold off

subplot(2,3,5) %Freq Noise as a function of power (10 Hz)
for kidn=1:length(KIDnumbers)
    if isempty(find(kidn==BADKIDS,1))
        plot(NOISEParameters(PsubPopt{kidn,1}(:),10),10*log10(NOISEParameters(PsubPopt{kidn,1}(:),12)),...
            'o','Color',KIDcolors(kidn,:),'MarkerFaceColor',KIDcolors(kidn,:),'MarkerSize',8) %Measurements
        hold on
        plot(FitSf{kidn,2}{1},'k') %Fit
    end
end
plot(Py,-189-0.5*(Py+40)+10,'k-','LineWidth',2) %Translated Gao Line
axis tight
legend('off')
xlabel('P_{int} [dBm]')
ylabel('S_F/F^2(10 Hz)')
title('Jiansong line 10 dB increase (S^2 \propto \sqrt{F})')
hold off

subplot(2,3,6) %Noise WRT complex plane
for kidn=1:length(KIDnumbers)
    semilogx(NOISE(PoptT0index(kidn,1)).FFTnoise{PoptT0index(kidn,2),1}(:,1),...
        NOISE(PoptT0index(kidn,1)).FFTnoise{PoptT0index(kidn,2),1}(:,6),...
        '-','color',KIDcolors(kidn,:),'LineWidth',2) %Phase Noise
    hold on
    semilogx(NOISE(PoptT0index(kidn,1)).FFTnoise{PoptT0index(kidn,2),1}(:,1),...
        NOISE(PoptT0index(kidn,1)).FFTnoise{PoptT0index(kidn,2),1}(:,7),...
        '--','color',KIDcolors(kidn,:),'LineWidth',2) %Amp Noise
    
    
end
xlim([1 1e6])
xlabel('F [Hz]')
ylabel('S_x')
title(['Noise at P_{opt} and ', num2str(NOISE(PoptT0index(1,1)).Temperature(PoptT0index(1,2)),2) ' K'])
legend('S_{\theta} wrt complex plane','S_A wrt complex plane')

hold off

%SAVE the figure
Figfile=[Noisepath,'KIDs@_',num2str(NOISE(PoptT0index(1,1)).Temperature(PoptT0index(1,2)),2),'K_NOISEcompare1'];
MakeGoodFigure(15,12,14,Figfile)
% saveas(gcf,Figfile,'fig')

figure(2)
QCmap = colormap(jet(41));
subplot(2,3,1) %Popt as a function of KID ID (and Q)
for kidn=1:length(KIDnumbers)
    plot(NOISEParameters(PoptT0index(kidn,1),1),NOISEParameters(PoptT0index(kidn,1),10),...
        'ko','MarkerSize',8,'LineWidth',2,'MarkerFaceColor',QCmap(QCindex(kidn),:))
    hold on
end
axis tight
hh=colorbar;
set(hh,'YTick',QCtick)
set(hh,'YTickLabel',num2str(tick2plot));
xlabel('KID ID')
ylabel('P^{int}_{opt} (dBm)')
title('color = Q/1000')
hold off

subplot(2,3,2) %Q as a function of KID ID
semilogy(NOISEParameters(PoptT0index(:,1),1),NOISEParameters(PoptT0index(:,1),5),...
    'ko','MarkerSize',5,'LineWidth',2,'MarkerFaceColor','k')
hold on
semilogy(NOISEParameters(PoptT0index(:,1),1),NOISEParameters(PoptT0index(:,1),6),...
    'bo','MarkerSize',5,'LineWidth',2,'MarkerFaceColor','b')
semilogy(NOISEParameters(PoptT0index(:,1),1),NOISEParameters(PoptT0index(:,1),7),...
    'ro','MarkerSize',5,'LineWidth',2,'MarkerFaceColor','r')
axis tight
xlabel('KID ID')
ylabel('Q')
legend('Q','Qi','Qc')
hold off
title(['All data at P_{opt} and T_{sys} = ',num2str(Tsystem),' K'])



subplot(2,3,3) %Freq Noise at Popt,10Hz as a function of KID (and Q)
for kidn=1:length(KIDnumbers)
    plot(NOISE(PoptT0index(kidn,1)).KIDnumber,10*log10(NOISE(PoptT0index(kidn,1)).MeanFreqNoise(PoptT0index(kidn,2),2)),...
        'ko','MarkerSize',8,'LineWidth',2,'MarkerFaceColor',QCmap(QCindex(kidn),:))
    hold on
end
axis tight
hh=colorbar;
set(hh,'YTick',QCtick)
set(hh,'YTickLabel',num2str(tick2plot));
xlabel('KID ID')
ylabel('S_F/F^2 (dBk/Hz)')
title('S_F/F^2(10 Hz) @P_{opt}')
hold off


subplot(2,3,4) %Freq Noise as a function of KID ID and modulation Freq at -40dBm internal power
for kidn=1:length(KIDnumbers)
    if isempty(find(kidn==BADKIDS,1))
        if ~isempty(FitSf{kidn,1}{1})
            plot(ID(kidn),FitSf{kidn,1}{1}.b,'go','MarkerSize',5,'MarkerFaceColor','g');hold on
        end
        if ~isempty(FitSf{kidn,2}{1})
            plot(ID(kidn),FitSf{kidn,2}{1}.b,'ro','MarkerSize',5,'MarkerFaceColor','r');hold on
        end
        if  ~isempty(FitSf{kidn,3}{1})
            plot(ID(kidn),FitSf{kidn,3}{1}.b,'bo','MarkerSize',5,'MarkerFaceColor','b');hold on
        end
        if ~isempty(FitSf{kidn,4}{1})
            plot(ID(kidn),FitSf{kidn,4}{1}.b,'ko','MarkerSize',5,'MarkerFaceColor','k');hold on
        end
        if kidn==1
            plot([-10,110],[-189,-189],'g-')
        end
    end
end
if Akiratestchip
    xlabel('Al percentage')
else
    xlabel('KID ID')
end
ylabel('S_F/F^2 [dBk/Hz] @-40dBm internal power')
legend('1 Hz','10 Hz','100 Hz','1 kHz','JSG-line 1 kHz 120mK')
xlim([min(ID)-0.1*max(ID),1.1*max(ID)])
ylim([-200,-160])
grid on
title('5 dB noise increase for 10x higher F = \sqrt(F) dependence.')
hold off

subplot(2,3,5) %Phase and Amp Noise at 10 Hz & Popt as function of KID ID.
for kidn=1:length(KIDnumbers)
    plot(NOISE(PoptT0index(kidn,1)).KIDnumber,NOISE(PoptT0index(kidn,1)).MeanNoise(PoptT0index(kidn,2),2),...
        'ro','MarkerSize',5,'MarkerFaceColor','r') %Amplitude
    hold on
    plot(NOISE(PoptT0index(kidn,1)).KIDnumber,NOISE(PoptT0index(kidn,1)).MeanNoise(PoptT0index(kidn,2),3),...
        'bo','MarkerSize',5,'MarkerFaceColor','b') %Phase
end
legend('S_A','S_{\theta}')
xlabel('KID ID')
ylabel('S_x (dBk/Hz)')
title('Noise at 10 Hz at P_{opt}')
hold off


subplot(2,3,6) %Freq Noise as a function of Fres and modulation Freq at -40dBm internal power
for kidn=1:length(KIDnumbers)
    if isempty(find(kidn==BADKIDS,1))
        if ~isempty(FitSf{kidn,1}{1})
            plot(NOISEParameters(PoptT0index(1,1),3)...
                ,FitSf{kidn,1}{1}.b,'go','MarkerSize',5,'MarkerFaceColor','g');hold on
        end
        if ~isempty(FitSf{kidn,2}{1})
            plot(NOISEParameters(PoptT0index(1,1),3)...
                ,FitSf{kidn,2}{1}.b,'ro','MarkerSize',5,'MarkerFaceColor','r');hold on
        end
        if  ~isempty(FitSf{kidn,3}{1})
            plot(NOISEParameters(PoptT0index(1,1),3)...
                ,FitSf{kidn,3}{1}.b,'bo','MarkerSize',5,'MarkerFaceColor','b');hold on
        end
        if ~isempty(FitSf{kidn,4}{1})
            plot(NOISEParameters(PoptT0index(1,1),3)...
                ,FitSf{kidn,4}{1}.b,'ko','MarkerSize',5,'MarkerFaceColor','k');hold on
        end
        if kidn==1
            plot([-10,110],[-189,-189],'g-')
        end
    end
end
if Akiratestchip
    xlabel('Al percentage')
else
    xlabel('KID Fres')
end
ylabel('S_F/F^2 [dBk/Hz] @-40dBm internal power')
legend('1 Hz','10 Hz','100 Hz','1 kHz','JSG-line 1 kHz 120mK')

ylim([-200,-160])
grid on
title('5 dB noise increase for 10x higher F = \sqrt(F) dependence.')
hold off


%SAVE the figure
Figfile=[Noisepath,'KIDs@_',num2str(NOISE(PoptT0index(1,1)).Temperature(PoptT0index(1,2)),2),'K_NOISEcompare2'];
MakeGoodFigure(15,12,14,Figfile)

%==========================================================================
% Wrap up
%==========================================================================
clear kidn
save([Noisepath,'Fits_at',num2str(NOISE(PoptT0index(1,1)).Temperature(PoptT0index(1,2)),2),'K.mat'],'FitSf');
save([Noisepath,'NOISEParameters_at',num2str(NOISE(PoptT0index(1,1)).Temperature(PoptT0index(1,2)),2),'K.mat'],'NOISEParameters');
save([Noisepath,'Noise.mat']);
rmpath([pwd,filesep,'subroutines']);
end