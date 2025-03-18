% NOISEanalysis_PdepV7
%
%DATE: October 24, 2012 - Aug 20 2019
%AUTHOR: Reinier Janssen, Jochem Baselmans
%==========================================================================
close all
clear all
global FITF0
%==========================================================================
% Inputs
%==========================================================================
%ChipInfo.path = [cd filesep '..']; %root path where data is, one higher than the scripts
ChipInfo.path = [ '..' filesep '..'];%without filesep at end
FFTsubdir = [ filesep 'FFT' filesep 'Power'];     %
ChipInfo.IndexPref = 3;     %Index of the power that is used as reference for Popt finding.
ReadPoptfile = 1;             %default=0. If set to 1 optimum power takemn from Popt file, that is created by this script. Usefull to be abel to rerun the sript for the correct Pot.
numlowPremoved = 1;         %for noise plots(P) - removes lowest power for less clutter
SaveStuff = 1;              %0 to not save files (faster) 1 to save
Tsystem = 7;                %estimated system noise, 8K oid
BWrange=1;                  %Width aroung Fres of the reference power that is considered to find Popt. Default 1
FITF0=0;                    %Switch used to determine which method is used for determination of F0.
%Please see the FitS21main5 routine for details. FITF0 = 0 is recommended due to the presence of overdriven
%resonators in Power sweeps (it will only update Q, not F0. F0 mot very important here).

%==========================================================================

%==========================================================================
% Start by setting some bla
%==========================================================================
format('long','g');     %Set display format of numbers to 7 digits
addpath([pwd,filesep,'..',filesep,'subroutines']);
close all;              %close all open plots to remove clutter.

%==========================================================================
% read Poptfile.
%==========================================================================
Noisepath = [ChipInfo.path,FFTsubdir,filesep]; %Path containing the raw noise data
if ReadPoptfile==1
    PoptFile = [Noisepath,'Popt.csv'];
    fid = fopen(PoptFile);
    if fid == -1 %no Popt file
        disp([Noisepath,'Popt.csv not found and ignored'])
        ReadPoptfile = 0;
    else %There is a Popt file. Read it and use it.
        fclose(fid);
        fprintf('Found Popt.csv, this will be used to indicate Popt.\n')
        [~,PoptData] = ReadSRONcsvV2(PoptFile,'',0);    %PoptData will be Nx3, the cols are [T,KIDid,|Popt [dBm]|]
    end
end
clear PoptFile FFTsubdir
%==========================================================================
%Search the noisepath for all FFT files and filter out all files that are
%significantly longer or shorter than the mean.
%==========================================================================

RawFFTfiles = dir([Noisepath,'KID*FFT*.dat']);
if isempty(RawFFTfiles)
    error(['No data available in' Noisepath]);
end

% %check FFT size for compatibility NOT NEEDED ANYMORE
% datasize = zeros(1,length(RawFFTfiles));
% for p=1:length(RawFFTfiles)
%     datasize(p)=RawFFTfiles(p).bytes;
% end
% FFT2read = 0.5*mean(datasize) < datasize & datasize < 1.5*mean(datasize);
% %RawFFTfiles = RawFFTfiles(FFT2read);
% clear datasize FFT2read 

%Determine the KID numbers (IDs) from the names of the good files.
KIDnumbers = zeros(1,length(RawFFTfiles));
NOISE(length(RawFFTfiles)).KIDnumber = ''; %pre allocate
for p=1:length(RawFFTfiles) %Loop over all FFT files
    KIDnumbers(p) = cell2mat(textscan(RawFFTfiles(p).name,'%*3s %f %*s'));
    NOISE(p).KIDnumber = KIDnumbers(p);
end
KIDnumbers = unique(KIDnumbers); %Determine all unique KIDs.

%Print to screen some of the reading information.
fprintf('Search for FFT data performed in:\n')
disp(Noisepath)
fprintf(['Inside the path the following KIDs are available: ',num2str(KIDnumbers),'\n'])
fprintf(['Inside this path a total number of ',num2str(length(RawFFTfiles)),' files were found.\n'])

%==========================================================================
% Read in all the data files
%==========================================================================
for p=1:length(RawFFTfiles) %LOOP OVER ALL FILES (aka KID-P-combinations)
    LocEndRoot = strfind(RawFFTfiles(p).name,'FFT.dat');
    if isempty(LocEndRoot)
        fprintf('ERROR Noise Analysis: Cannot find FFT.dat in the name of the noise file.\n')
        fprintf('Most likely text has been placed between FFT and .dat extension.\n')
        error('Cannot create RootName for file detection.\n')
    end
    
    RootName = [Noisepath,RawFFTfiles(p).name(1:LocEndRoot-1)];
    NOISE(p).filename = RootName;
    clear LocEndRoot RootName
    %======================================================================
    %One by one read in the data files
    %======================================================================
    
    %======================================================================
    %Reading FFT file
    FFTfile = [NOISE(p).filename,'FFT.dat'];
    [Data,Temperature,Power,FFTheader] = import_data(FFTfile);%Data = cell array
    % limit T range to 1
    if length(Temperature) > 1
        Temperature(MaxnT+1:end)=[];
    end
    NOISE(p).Temperature = Temperature;
    Data(2:end)=[];
    NOISE(p).FFTnoise = Data;
    clear Data;
    NOISE(p).ReadPower = -1*Power;
    
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
    NOISE(p).Fread = Fread;
    clear ColonIndex Fread Freadindex Freadstart FFTfile hl FFTheader
    %======================================================================
    %Reading S21 file (GHz,Re,Im)
    S21file = [NOISE(p).filename,'S21.dat'];
    [Data,~,Power,~] = import_data(S21file);
    %Check for equality of temperature and power
    if abs(-1*Power-NOISE(p).ReadPower) > 0.5
        disp('Warning: Read Power does not match between FFT and S21 file.')
    end
    %Copy data for struct
    Data(2:end)=[];
    NOISE(p).S21_IQplane = Data;
    clear S21file Data
    %======================================================================
    %Reading S21 file (GHz,dB,rad)
    S21DBfile = [NOISE(p).filename,'S21dB.dat'];
    [Data,~,Power] = import_data(S21DBfile);
    %Check for equality of temperature and power
    if abs(-1*Power-NOISE(p).ReadPower(1)) > 0.5
        disp('Warning: Read Power does not match between FFT and S21dB file.')
    end
    %Copy data for struct
    Data(2:end)=[];
    NOISE(p).S21_MPplane = Data;
    clear Data S21DBfile
    %======================================================================
    %Reading time-domain file
    TDfile = [NOISE(p).filename,'td.dat'];
    [Data,~,Power,TDHeader] = import_data(TDfile);
    %Check for equality of temperature and power
    if abs(-1*Power-NOISE(p).ReadPower(1)) > 0.5
        disp('Warning: Read Power does not match between FFT and TD file.')
    end
    %Extract the timestep from the first subheader.
    fid=fopen(TDfile);
    dt = cell2mat(textscan(fid,'%*s%*s%*s%f',1,'headerlines',length(TDHeader)+1));
    fclose(fid);
    %Copy data for struct
    CompleteData = zeros(size(Data{1},1),size(Data{1},2)+1);
    CompleteData(:,1) = dt*(1:size(Data{1},1))';
    CompleteData(:,2:end) = Data{1}(:,:);
    NOISE(p).TDIQ{1,1} = CompleteData;
    
    disp(['Highest Read: ' num2str(NOISE(p).Temperature(end))]);
    clear Temperature Data dt TDfile CompleteData fid TDHeader
    
    %======================================================================
    %Analyse the S21 data measured in [mag,phase] space for each temperature
    %======================================================================
    %Normalize the S21 data measured in magnitude plane
    filtered=smooth(NOISE(p).S21_MPplane{1}(:,2),3); %smoothing the |S21| data in dB space
    %normalise S21 in log space to the max(|S21|)
    NOISE(p).S21_MPplane{1}(:,2)=NOISE(p).S21_MPplane{1}(:,2)-max(filtered);
    %Convert dB to magnitude
    NOISE(p).S21_MPplane{1}(:,2) = 10.^(NOISE(p).S21_MPplane{1}(:,2)/20);
    clear filtered
    %Perform fit to obtain resonator parameters. Note the FITF0 is
    %recommended for overdriven resonators.
    [fres,Q,S21min,FitResult] = FitS21main5(NOISE(p).S21_MPplane{1}(:,1:3),FITF0);
    
    %Put the temporary storage variables into the NOISE struct
    NOISE(p).Ql=Q(1);
    NOISE(p).Qi=Q(2);
    NOISE(p).Qc=Q(3);
    NOISE(p).Fres=fres;
    NOISE(p).S21min=S21min;%in dB!
    NOISE(p).S21fit = FitResult;
    %Calculate the Internal Power
    NOISE(p).InternalPower = 10*log10((2/pi)*10.^(NOISE(p).ReadPower/10).*(NOISE(p).Ql.^2./NOISE(p).Qc));
    %Calculate the resonator ring time
    NOISE(p).TauRes = NOISE(p).Ql./(pi*NOISE(p).Fres);
    %Calculate the resonator bandwidth (used later in Popt determination)
    NOISE(p).Bandwidth = NOISE(p).Fres./NOISE(p).Ql;
    
    %Normalized Frequency Noise (Sf/F^2) [1/Hz] (defined by J. Gao) - setup corrected %
    tempnoise = (10.^(NOISE(p).FFTnoise{1}(:,2)/10))*(1/(4*NOISE(p).Ql(1,1))^2);
    NOISE(p).Setup_level = mean(tempnoise(end-10:end-5));
    NOISE(p).FFTnoise{1}(:,4) = tempnoise - NOISE(p).Setup_level;
    clear tempnoise Setup_level fres Q S21min FitResult
    
    %Frequency Noise (defined by B. Mazin) not used
    %NOISE(p).FFTnoise{1}(:,5) = NOISE(p).FFTnoise{1}(:,4)*(NOISE(p).Fres(1,1)*1e9)^2;
    
    %Phase noise WRT complex plane
    S21min_a=10^(( NOISE(p).S21min(1) )/20);%from dB to magnitude, not stored
    R=20*log10((1-S21min_a)/2);%correction due to circl radius vs complex plane radius in dB
    NOISE(p).FFTnoise{1}(:,6)=NOISE(p).FFTnoise{1}(:,2)+R;
    
    %Amplitude noise WRT complex plane
    NOISE(p).FFTnoise{1}(:,7)=NOISE(p).FFTnoise{1}(:,3)+R;
    
    %Calculate Mean amplitude noise between 20 and 500 Hz (Required for Popt determination)%
    Frange = 20 <= NOISE(p).FFTnoise{1,1}(:,1) & NOISE(p).FFTnoise{1,1}(:,1) <= 500;
    NOISE(p).MeanNoiseforPopt = mean(NOISE(p).FFTnoise{1,1}(Frange,3));
    
    %Calculate the mean frequency noise (Gao definition), with setup
    %contribution subtracted (was done already in defining Gao's freq
    %noise)
    if min(NOISE(p).FFTnoise{1,1}(:,1)) < 1000 %data exists (no cosmiuc ray proble,)
        NOISE(p).MeanFreqNoise1Hz = ... % 1 Hz
            mean(interp1(NOISE(p).FFTnoise{1,1}(:,1),NOISE(p).FFTnoise{1,1}(:,4),(0.8:0.1:1.2),'linear',1));
        NOISE(p).MeanFreqNoise10Hz = ... % 10 Hz
            mean(interp1(NOISE(p).FFTnoise{1,1}(:,1),NOISE(p).FFTnoise{1,1}(:,4),(9:0.1:11),'linear',1));
        NOISE(p).MeanFreqNoise100Hz = ... % 100 Hz
            mean(interp1(NOISE(p).FFTnoise{1,1}(:,1),NOISE(p).FFTnoise{1,1}(:,4),(90:1:110),'linear',1));
        NOISE(p).MeanFreqNoise1000Hz = ... % 1 kHz
            mean(interp1(NOISE(p).FFTnoise{1,1}(:,1),NOISE(p).FFTnoise{1,1}(:,4),(900:10:1100),'linear',1));
    else
        NOISE(p).MeanFreqNoise1Hz = NaN;
        NOISE(p).MeanFreqNoise01Hz = NaN;
        NOISE(p).MeanFreqNoise100Hz = NaN;
        NOISE(p).MeanFreqNoise1000Hz = NaN;
    end
    % system noise
    S21minMag = 10.^(NOISE(p).S21min/20);
    NOISE(p).SystemNoiseT = Tsystem; %Store system noise Temperature
    NOISE(p).SetupNoise = 10*log10(1.38e-23*Tsystem*(2./(1-S21minMag)).^2)-(NOISE(p).ReadPower-30);
    
end %END OF LOOP OVER ALL FILES (aka KID-Pread-combinations)
clear Power R S21min_a Data Frange RawFFTfiles FITF0

%==========================================================================
% Now that all files have been read we determine Popt.
%==========================================================================
IndexPsort = cell(length(KIDnumbers),1);
% initialize cell array to store for each KID an array that contains:
% contains vector of length = number of readout powers (NP). This vector holds the indexes in NOISE
% that contain information about this resonator. The indexes are given in
% order of increasing power.
% IndexPsort{kidn} = array all indices for kidn in NOISE

NOISEParams = cell(length(KIDnumbers),1);IndexP_sub_opt = NOISEParams;IndexP_super_opt = NOISEParams;
IndexPopt = zeros(length(KIDnumbers),1);
for kidn=1:length(KIDnumbers) % LOOP OVER ALL UNIQUE KIDS,
    %======================================================================
    % create IndexPsort
    %======================================================================
    KIDlocations = find([NOISE(:).KIDnumber] == KIDnumbers(kidn));  %gets the indices in NOISE that contains data on this KID (i.e which files this KID is descrived in)
    [~,blaindices]=sort([NOISE(KIDlocations).ReadPower]);           %Sort on Pread and save the indexes.
    IndexPsort{kidn} = KIDlocations(blaindices);                    %Recalibrate the indexes to the NOISE indices. This vector can be used to get the sorted P read for this KID, valid for all temperaturesd
    clear KIDlocations blaindices
    %Check if the desired reference index is within the number of powers measured for this KID%
    if ChipInfo.IndexPref > length(IndexPsort{kidn})
        disp('WARNING NoiseAnalysis: IndexPref exceeds number of Powers. Assuming maximum Power to be reference.')
        ChipInfo.IndexPref = length(IndexPsort{kidn});
    end
    
    %local Shorthand variable for the reference power index (valid locally for this kid as we loop over KIDs here)%
    IndexPref = IndexPsort{kidn}(ChipInfo.IndexPref);
    
    %==================================================================
    % DETERMINE Popt
    %==================================================================
    %Determine which powers have to big a resonance frequency shift
    if ReadPoptfile == 0
        %Put some user info to screen. ReadPower(1) refers to lowest T
        disp(['Searchin optimum power for KID ',num2str(KIDnumbers(kidn)),...
            ' in Power range: ',num2str(NOISE(IndexPsort{kidn}(1)).ReadPower(1)),...
            ' - ',num2str(NOISE(IndexPsort{kidn}(end)).ReadPower(1)),' dBm']);
        %
        FresPref = NOISE(IndexPref).Fres;
        Frange = BWrange*NOISE(IndexPref).Bandwidth;
        AllowedPindices = FresPref-Frange <= [NOISE(IndexPsort{kidn}).Fres] & ...
            [NOISE(IndexPsort{kidn}).Fres] <= FresPref+Frange;
        AllowedPindices = find(AllowedPindices);%convert logical to index array. valid in PT_Fres (i.e. for this KID) == valid inside IndexPsort{kidn})%
        
        %Determine Popt in the allowed P range using AmpNoise.
        [~,PoptIndex] = min([NOISE(IndexPsort{kidn}(AllowedPindices)).MeanNoiseforPopt]);
        tempPoptIndex = AllowedPindices(PoptIndex);
        
        % optionally
    elseif ReadPoptfile==1
         %Put some user info to screen. ReadPower(1) refers to lowest T
        disp(['Reading optimum power for KID ',num2str(KIDnumbers(kidn))]);
        %
        Poptdata_Tindices= PoptData(:,1)>=0.9*NOISE(IndexPref).Temperature(1,1) & PoptData(:,1)<=1.1*NOISE(IndexPref).Temperature(1,1);%Note MaxnT = length(NOISE(IndexPref).Temperature);
        Poptdata_thisKID_thisT=PoptData(Poptdata_Tindices,2)==NOISE(IndexPref).KIDnumber;
        bla=find(Poptdata_Tindices);%this is an array of indiced in PoptData that are allowed wrt temperature
        Popt_fromfile=PoptData(bla(Poptdata_thisKID_thisT),3);
        tempPoptIndex=find([NOISE(IndexPsort{kidn}).ReadPower]==Popt_fromfile);
    else
        error('ReadPoptfile must be 0 or 1')
    end
    %==================================================================
    % New indices for NOISE struct: One for Popt, one for all powers below Popt%
    %==================================================================
    % ******************* IndexPopt *******************
    IndexPopt(kidn) = IndexPsort{kidn}(tempPoptIndex); %Index of Popt in full NOise struct (for each T)  IndexPopt(kidn)
    
    % find the indices above and below in IndexPsort{kidn}
    temp_sub_opt    = ([NOISE(IndexPsort{kidn}).ReadPower] <= NOISE(IndexPopt(kidn)).ReadPower); %index in Indexpsort
    temp_super_opt  = ([NOISE(IndexPsort{kidn}).ReadPower] > NOISE(IndexPopt(kidn)).ReadPower);
    
    % ******************* indices above and below popt *******************
    IndexP_sub_opt{kidn}    = IndexPsort{kidn}(temp_sub_opt);%Index of all powers below Popt in full NOise struct (for each T)  IndexPopt(kidn)
    IndexP_super_opt{kidn}  = IndexPsort{kidn}(temp_super_opt);
    
    % create Noise paraemeters struct to allow eacy KID to KID plots
    NOISEParams{kidn}.Popt = NOISE(IndexPopt(kidn)).ReadPower;
    NOISEParams{kidn}.Popt_int = NOISE(IndexPopt(kidn)).InternalPower;
    clear FresPref AllowedPindices PoptIndex tempPoptIndex Poptdata_Tindices Poptdata_thisKID_thisT bla Popt_fromfile temp_sub_opt temp_super_opt
    
    %======================================================================
    % Fit all the Pint(noise @ 1 F) data to obtain the noise at Pint = -40 dBm for all T
    %======================================================================
    
    FitOptionsS1kHz = fitoptions('Method','NonlinearLeastSquares','StartPoint',[-0.5 -100]);
    FitTypeS1kHz = fittype('a*(x+50)+b','options',FitOptionsS1kHz);
    
    % swicthh off fit warning
    warning('off', 'curvefit:fit:complexYusingOnlyReal');
    %Stores the fit to the 1 kHz frequency noise as a function of internal power.%
    
        tfp     = [NOISE(IndexP_sub_opt{kidn}).InternalPower]';
        tdfata  = 10*log10([NOISE(IndexP_sub_opt{kidn}).MeanFreqNoise1000Hz]');
     if length(tdfata(~isnan(tdfata))) > 1
        FitS1kHz    = fit( tfp(~isnan(tdfata)) , tdfata(~isnan(tdfata)) , FitTypeS1kHz);
        
        tdfata  = 10*log10([NOISE(IndexP_sub_opt{kidn}).MeanFreqNoise100Hz]');
        FitS100Hz   = fit( tfp(~isnan(tdfata)) , tdfata(~isnan(tdfata)) , FitTypeS1kHz);
        
        tdfata  = 10*log10([NOISE(IndexP_sub_opt{kidn}).MeanFreqNoise10Hz]');
        FitS10Hz	= fit( tfp(~isnan(tdfata)) , tdfata(~isnan(tdfata)) , FitTypeS1kHz);
        
        tdfata  = 10*log10([NOISE(IndexP_sub_opt{kidn}).MeanFreqNoise1Hz]');
        FitS1Hz     = fit( tfp(~isnan(tdfata)) , tdfata(~isnan(tdfata)) , FitTypeS1kHz);
    else %Not enought points for a linear fit. Leave all empty.
        FitS1kHz    = [];
        FitS100Hz	= [];
        FitS10Hz	= [];
        FitS1Hz     = [];
    end
    NOISEParams{kidn}.FitS1kHz    = FitS1kHz;
    NOISEParams{kidn}.FitS100Hz   = FitS100Hz;
    NOISEParams{kidn}.FitS10Hz    = FitS10Hz;
    NOISEParams{kidn}.FitS1Hz     = FitS1Hz;
    clear FitS1kHz FitS100Hz FitS10Hz FitS1Hz FitOptionsS1kHz FitTypeS1kHz
    
    %==================================================================
    % Figure one: mainly checks for the quality of the analysis
    % routine
    %==================================================================
    Pcolors = colormapJetJB(length(IndexPsort{kidn}));
    
    figure(1000*KIDnumbers(kidn)+1)
    clf
    
    clear PowerLegend
    PowerLegend = cell(length(IndexPsort{kidn})+1,2);
    PowerLegend{1,1} = 'P_{opt}';
    PowerLegend{1,2} = 'P^{int}_{opt}';
    %==============================================================
    %Resonance circle as a function of power (incl noise blobs)
    %==============================================================
    subplot(2,3,1)
    hold on
    for nP=1:length(IndexPsort{kidn})
        plot(NOISE(IndexPsort{kidn}(nP)).S21_IQplane{1}(:,2) , NOISE(IndexPsort{kidn}(nP)).S21_IQplane{1}(:,3),'-','color',Pcolors(nP,:),'LineWidth',1)
        PowerLegend{nP+1,1} = num2str(NOISE(IndexPsort{kidn}(nP)).ReadPower,'%.1f');
        PowerLegend{nP+1,2} = num2str(NOISE(IndexPsort{kidn}(nP)).InternalPower,'%.1f');
    end
    for nP=1:length(IndexPsort{kidn}) %Second Loop to get the legend correct.
        plot(NOISE(IndexPsort{kidn}(nP)).TDIQ{1}(:,2),NOISE(IndexPsort{kidn}(nP)).TDIQ{1}(:,3),'.','color',Pcolors(nP,:),'MarkerSize',6)
    end
    legend(PowerLegend(2:end,1))
    xlabel('Re');ylabel('Im')
    title(['KID ',num2str(KIDnumbers(kidn),'%.0f'),' @T=',num2str(NOISE(IndexPref).Temperature,'%.3g'),' K'])
    box on;grid on;
    hold off
    
    %==============================================================
    %Resonance Dip as a function of power, Incl reference lines around reference power%
    %==============================================================
    subplot(2,3,2)
    hold on
    for nP=1:length(IndexPsort{kidn})
        plot(NOISE(IndexPsort{kidn}(nP)).S21_MPplane{1}(:,1),20*log10(NOISE(IndexPsort{kidn}(nP)).S21_MPplane{1}(:,2)),'-','color',Pcolors(nP,:),'LineWidth',1)
    end
    Flow = NOISE(IndexPopt(kidn)).Fres-BWrange*NOISE(IndexPopt(kidn)).Bandwidth;
    plot([Flow,Flow],[0,-20],'r-','LineWidth',2)
    Fhigh = NOISE(IndexPopt(kidn)).Fres+BWrange*NOISE(IndexPopt(kidn)).Bandwidth;
    plot([Fhigh,Fhigh],[0,-20],'r-','LineWidth',2)
    xlabel('F [GHz]');ylabel('|S21| [dB]')
    title(['KID ',num2str(KIDnumbers(kidn),'%.0f'),' P_{opt}=',num2str(NOISE(IndexPopt(kidn)).ReadPower),' dBm'])
    axis tight;
    xlim([NOISE(IndexPopt(kidn)).Fres - 2*BWrange*NOISE(IndexPopt(kidn)).Bandwidth...
        NOISE(IndexPopt(kidn)).Fres   + 2*BWrange*NOISE(IndexPopt(kidn)).Bandwidth])
    box on;grid on;
    hold off
    
    %==============================================================
    %Time Domain Trace as a function of time
    %==============================================================
    subplot(2,3,3)
    for nP=1:length(IndexPsort{kidn})
        plot(NOISE(IndexPsort{kidn}(nP)).TDIQ{1}(:,1),NOISE(IndexPsort{kidn}(nP)).TDIQ{1}(:,2),...
            'o','color',Pcolors(nP,:),'MarkerSize',3,'MarkerFaceColor',Pcolors(nP,:));
        hold on
        plot(NOISE(IndexPsort{kidn}(nP)).TDIQ{1}(:,1),NOISE(IndexPsort{kidn}(nP)).TDIQ{1}(:,3),...
            '.','color',Pcolors(nP,:),'MarkerSize',3);
    end
    xlabel('t [sec]'); ylabel('Re or Im')
    legend('Re(S21)','Im(S21)');
    box on;
    hold off
    
    %==============================================================
    %Frequency Noise
    %==============================================================
    subplot(2,3,4)
    warning('off', 'MATLAB:plot:IgnoreImaginaryXYPart');
    toplot = NOISE(IndexPopt(kidn)).FFTnoise{1}(:,4) > 0;
    semilogx(NOISE(IndexPopt(kidn)).FFTnoise{1}(toplot,1),10*log10(NOISE(IndexPopt(kidn)).FFTnoise{1}(toplot,4)),...
        '-','color','k','LineWidth',3)
    hold on
    for nP=1:length(IndexPsort{kidn})
        toplot = NOISE(IndexPsort{kidn}(nP)).FFTnoise{1}(:,4) > 0;
        semilogx(NOISE(IndexPsort{kidn}(nP)).FFTnoise{1}(toplot,1),10*log10(NOISE(IndexPsort{kidn}(nP)).FFTnoise{1}(toplot,4)),...
            '-','color',Pcolors(nP,:),'LineWidth',1)
    end
    xlabel('F [Hz]');ylabel('S_F/F^2 [dBc/Hz]')
    legend(PowerLegend(:,2))
    title(['KID ',num2str(KIDnumbers(kidn),'%.0f'),' P^{int}_{opt}=',num2str(NOISE(IndexPopt(kidn)).InternalPower,'%.1f'),' dBm'])
    xlim([0.5,1e4]);grid on;ylim([-220,-140])
    hold off
    
    %==============================================================
    %Phase Noise (and Amp Noise)
    %==============================================================
    subplot(2,3,[5,6])
    semilogx(NOISE(IndexPopt(kidn)).FFTnoise{1}(:,1),NOISE(IndexPopt(kidn)).FFTnoise{1}(:,3),...
        '-','color','k','LineWidth',3)
    hold on
    semilogx(NOISE(IndexPopt(kidn)).FFTnoise{1}(:,1),NOISE(IndexPopt(kidn)).FFTnoise{1}(:,2),...
        '-','color','k','LineWidth',4)
    maxnp=zeros(length(IndexPsort{kidn})-numlowPremoved,1);minnp=maxnp;
    for nP=numlowPremoved+1:length(IndexPsort{kidn})
        % finding range
        maxnp(nP-numlowPremoved) = max(NOISE(IndexPsort{kidn}(nP)).FFTnoise{1}(NOISE(IndexPsort{kidn}(nP)).FFTnoise{1}(:,1) > 10,2));
        minnp(nP-numlowPremoved) = mean(NOISE(IndexPsort{kidn}(nP)).FFTnoise{1}(20:end-5,3));
        semilogx(NOISE(IndexPsort{kidn}(nP)).FFTnoise{1}(:,1),NOISE(IndexPsort{kidn}(nP)).FFTnoise{1}(:,3),...
            '-','color',Pcolors(nP,:),'LineWidth',1)
        semilogx(NOISE(IndexPsort{kidn}(nP)).FFTnoise{1}(:,1),NOISE(IndexPsort{kidn}(nP)).FFTnoise{1}(:,2),...
            '-','color',Pcolors(nP,:),'LineWidth',2)
    end
    grid on;
    xlabel('F [Hz]');ylabel('S_x [dBc/Hz]')
    legend('S_A','S_{\theta}')
    xlim([10,0.3e6]);ylim([10*round(min(0.1*minnp))-5,10*round(max(0.1*maxnp))+5]);
    hold off
    
    %==============================================================
    %SAVE the figure
    %==============================================================
    Figfile=[Noisepath,'KID',num2str(KIDnumbers(kidn),'%.0f'),'_',...
        num2str(NOISE(IndexPref).Temperature,'%.3g'),'K_NOISE1'];
    if SaveStuff == 1
        MakeGoodFigure(14,13,12,Figfile)
    else
        MakeGoodFigure(14,13,12)
    end
    clear Fhigh Flow maxnp minnp nP PowerLegend toplot Figfile
    
    %==================================================================
    %==================================================================
    % Figure two: Overview of resonator parameters as a function of
    % power.
    %==================================================================
    %==================================================================
    figure(1000*KIDnumbers(kidn)+2)
    clf
    % evaluate the noise fit function
    Py =(-80:-20);
    if isempty(NOISEParams{kidn}.FitS1kHz)
        FitS1kHzEval = [];
        NOISEParams{kidn}.FitS1kHzEval = [];
    else
        FitS1kHzEval = feval(NOISEParams{kidn}.FitS1kHz,Py);
        NOISEParams{kidn}.FitS1kHzEval(:,2) = 10.^(FitS1kHzEval/10);
        NOISEParams{kidn}.FitS1kHzEval(:,1) = Py;
    end
    
    %Make the legend for the noise figure
    NoiseLegend{1} = 'P<P_{opt}'; %Since there is always at least 1 measurement, which is Popt
    NoiseLegend{2} = 'P=P_{opt}'; %Since there is always at least 1 measurement, which is Popt
    NoiseLegend{3} = 'Gao Line';
    if length(IndexPsort{kidn}) > length(IndexP_sub_opt{kidn})
        NoiseLegend{4} = 'P>P_{opt}';
        if ~isempty(FitS1kHzEval)
            NoiseLegend{5} = 'Fit';
        end
    else
        if ~isempty(FitS1kHzEval)
            NoiseLegend{4} = 'Fit';
        end
    end
    %==============================================================
    %Frequency Noise at 1 kHz as a function of Pint
    %==============================================================
    subplot(2,2,1)
    warning('off','MATLAB:Axes:NegativeDataInLogAxis');
    %Measurement Points
    semilogy([NOISE(IndexPsort{kidn}).InternalPower],[NOISE(IndexPsort{kidn}).MeanFreqNoise1000Hz],'ro','MarkerSize',6) %all P
    hold on
    semilogy([NOISE(IndexPopt(kidn)).InternalPower],[NOISE(IndexPopt(kidn)).MeanFreqNoise1000Hz],'ro','MarkerSize',6,'MarkerFaceColor','r') %At Popt
    semilogy(Py,10.^((-189-0.5*(Py+40))/10),'b-','LineWidth',1) %Jiansong Line
    if ~isempty(IndexP_super_opt{kidn})
        semilogy([NOISE(IndexP_super_opt{kidn}).InternalPower],[NOISE(IndexP_super_opt{kidn}).MeanFreqNoise1000Hz],'kd','MarkerSize',6) %Above Popt
    end
    %Lines
    if isempty(FitS1kHzEval) == 0 %Only if fit has been made
        semilogy(Py,10.^(FitS1kHzEval/10),'r-','LineWidth',1) %Fit to P<Popt
    end
    %make a nice plot.
    axis tight
    title(['KID ',num2str(KIDnumbers(kidn),'%.0f'),' @T=',num2str(NOISE(IndexPref).Temperature,'%.3g'),' K'])
    legend(NoiseLegend)
    xlabel('P_{int} (dBm)')
    ylabel('S_f/f^2 @1kHz [dB]')
    box on;grid on;
    hold off
    
    %==============================================================
    %Resonance Frequency as a function of Pint
    %==============================================================
    subplot(2,2,2)
    hold on
    semilogy([NOISE(IndexPsort{kidn}).InternalPower],[NOISE(IndexPsort{kidn}).Fres],'ro','MarkerSize',6) %all P
    hold on
    semilogy([NOISE(IndexPopt(kidn)).InternalPower],[NOISE(IndexPopt(kidn)).Fres],'ro','MarkerSize',6,'MarkerFaceColor','r') %At Popt
    if ~isempty(IndexP_super_opt{kidn})
        semilogy([NOISE(IndexP_super_opt{kidn}).InternalPower],[NOISE(IndexP_super_opt{kidn}).Fres],'kd','MarkerSize',6) %Above Popt
    end
    xlabel('P_{int} (dBm)')
    ylabel('f_{res} (GHz)')
    title(['KID ',num2str(KIDnumbers(kidn),'%.0f'),' P_{opt}=',num2str(NOISE(IndexPopt(kidn)).ReadPower),' dBm'])
    box on;grid on;
    hold off
    
    %==============================================================
    %Qi as a function of Pint
    %==============================================================
    subplot(2,2,3)
    hold on
    %Measurement Points
    semilogy([NOISE(IndexPsort{kidn}).InternalPower],[NOISE(IndexPsort{kidn}).Qi],'ro','MarkerSize',6); %Below Popt
    hold on
    semilogy([NOISE(IndexPopt(kidn)).InternalPower],[NOISE(IndexPopt(kidn)).Qi],'ro','MarkerSize',6,'MarkerFaceColor','r');%At Popt
    
    if ~isempty(IndexP_super_opt{kidn})
        semilogy([NOISE(IndexP_super_opt{kidn}).InternalPower],[NOISE(IndexP_super_opt{kidn}).Qi],'kd','MarkerSize',6) %Above Popt
    end
    xlabel('P_{int} (dBm)')
    ylabel('Q_{i} (dBm)')
    box on;grid on;
    hold off
    
    %==============================================================
    %S21 resonance fit at Popt
    %==============================================================
    subplot(2,2,4)
    %Measurement Points
    plot(NOISE(IndexPopt(kidn)).S21_MPplane{1}(:,1),NOISE(IndexPopt(kidn)).S21_MPplane{1}(:,2),'b.','MarkerSize',5)
    axis tight;hold on
    %Fit
    plot(NOISE(IndexPopt(kidn)).S21fit(:,1),10.^(NOISE(IndexPopt(kidn)).S21fit(:,2)/20),'r-','LineWidth',1)
    plot(NOISE(IndexPopt(kidn)).Fres,10.^(NOISE(IndexPopt(kidn)).S21min/20),'kd',...
        'MarkerFaceColor','r','MarkerSize',6)%The determined Fres and |S21(Fres)|
    
    %make the plot nice
    xlabel('F (GHz)')
    ylabel('|S21|')
    title('Resonance @Popt')
    hold off
    %==============================================================
    %SAVE the figure
    %==============================================================
    Figfile=[Noisepath,'KID',num2str(KIDnumbers(kidn),'%.0f'),'_',...
        num2str(NOISE(IndexPref).Temperature,'%.3g'),'K_NOISE2'];
    if SaveStuff == 1
        MakeGoodFigure(12,8,14,Figfile)
    else
        MakeGoodFigure(12,8,14)
    end
    
    clear Pcolors Py FitS1kHzEval NoiseLegend NLindex Figfile
    
end %END OF LOOP OVER ALL UNIQUE KIDS

%==========================================================================
% Write the Popt file
%==========================================================================
PoptData = zeros(length(KIDnumbers),3);
PoutData = zeros(length(KIDnumbers),9);
for kidn=1:length(KIDnumbers)
    PoptData(kidn,1) = NOISE(IndexPopt(kidn)).Temperature;
    PoptData(kidn,2) = NOISE(IndexPopt(kidn)).KIDnumber;
    PoptData(kidn,3) = NOISEParams{kidn}.Popt;
    PoutData(kidn,4) = round(NOISE(IndexPopt(kidn)).InternalPower);
    PoutData(kidn,5) = NOISE(IndexPopt(kidn)).Fres;
    PoutData(kidn,6) = 0.1*round(NOISE(IndexPopt(kidn)).S21min*10);
    PoutData(kidn,7) = round(NOISE(IndexPopt(kidn)).Ql/1e3);
    PoutData(kidn,8) = round(NOISE(IndexPopt(kidn)).Qc/1e3);
    PoutData(kidn,9) = round(NOISE(IndexPopt(kidn)).Qi/1e3);
end
PoutData(:,1:3) = PoptData;
PoptHeader = {'T (K)','KIDID','Popt (dBm)'}';
WriteSRONcsv([Noisepath,'Popt.csv'],PoptData,PoptHeader,'%.3g')

PoutHeader = {'T (K)','KIDID','Popt (dBm)',...
    'Pint_opt (dBm)','Fres (GHz)','S21 (dB)','Q','Qc','Qi'}';
WriteSRONcsv([Noisepath,'Pout.csv'],PoutData,PoutHeader,'%.6g')

clear PoptData Pcolors BWrange p PoptHeader kidn Chipinfo

%==========================================================================
% Do KID by KID analysis - this was originally in Noisecompare
%==========================================================================

if ReadPoptfile == 1
    figure(1)
    %simple colors
    kidcolors = colormapJetJB(length(KIDnumbers));
    %simple legend
    lstr = cell(length(KIDnumbers),1);
    %Generate the colormap to view Q distribution
    QCmap = colormapJetJB(31);%range from 3:0.1:6 = log(Qc)
    QCrange = 3:0.1:6;
    QCtick = [3,4,5,6];
    %tick2plot = num2cell(10.^ QCtick);%could be used, but less nice
    PR = [min([NOISE(IndexPopt).KIDnumber])-1, max([NOISE(IndexPopt).KIDnumber])+1];
    minp = zeros(length(KIDnumbers),1);maxnp = minp;
    for kidn=1:length(KIDnumbers)
        %===============================
        %noise spectra @ Popt
        %===============================
        subplot(2,3,1)
        semilogx(NOISE(IndexPopt(kidn)).FFTnoise{1}(:,1),NOISE(IndexPopt(kidn)).FFTnoise{1}(:,3),...
            '-','color',kidcolors(kidn,:),'LineWidth',1)
        grid on;
        %setup noise
        semilogx(1e5,NOISE(IndexPopt(kidn)).SetupNoise,'ko','MarkerFaceColor',kidcolors(kidn,:),'LineWidth',2,'MarkerSize',8)
        xlabel('F [Hz]');ylabel('S_x [dBc/Hz]')
        xlim([10,0.5e6]);
        semilogx(NOISE(IndexPopt(kidn)).FFTnoise{1}(:,1),NOISE(IndexPopt(kidn)).FFTnoise{1}(:,2),...
            '-','color',kidcolors(kidn,:),'LineWidth',2)
        hold on;
        
        %===============================
        %Sf's
        %===============================
        subplot(2,3,2)
        semilogx(NOISE(IndexPopt(kidn)).FFTnoise{1}(:,1),10*log10(NOISE(IndexPopt(kidn)).FFTnoise{1}(:,4)),...
            '-','color',kidcolors(kidn,:),'LineWidth',2)
        hold on;box on;
        xlabel('F  (Hz)');ylabel('S_F/F^2 [dBk/Hz] ')
        xlim([5,2e4]);grid on;ylim([-220,-140])
        lstr{kidn} = ['KID ' num2str(NOISE(IndexPopt(kidn)).KIDnumber)];
        minp(kidn) = min(10*log10(abs(NOISE(IndexPopt(kidn)).FFTnoise{1}(20:end-25,4))));
        maxnp(kidn) = max(10*log10(abs(NOISE(IndexPopt(kidn)).FFTnoise{1}(20:end-25,4))));%ylim([floor(min(minp)) ,ceil(max(maxp))])
        %===============================
        %1 kHz Jiansong data + fit (P<=Popt)
        %===============================
        subplot(2,3,3)
        semilogy([NOISE(IndexP_sub_opt{kidn}).InternalPower],[NOISE(IndexP_sub_opt{kidn}).MeanFreqNoise1000Hz],...
            'o','color',kidcolors(kidn,:),'MarkerSize',6);
        hold on
        %Lines
        if ~isempty(NOISEParams{kidn}.FitS1kHzEval)  %Only if fit has been made
            semilogy(NOISEParams{kidn}.FitS1kHzEval(:,1),NOISEParams{kidn}.FitS1kHzEval(:,2),'-','color',kidcolors(kidn,:),'LineWidth',1) %Fit to P<=Popt
            semilogy(NOISEParams{kidn}.FitS1kHzEval(:,1),10.^((-189-0.5*(NOISEParams{kidn}.FitS1kHzEval(:,1)+40))/10),'k-','LineWidth',2) %Jiansong Line
            semilogy(-50,10.^(NOISEParams{kidn}.FitS1kHz.b/10),'o','color',kidcolors(kidn,:),'MarkerSize',6,'MarkerFaceColor',kidcolors(kidn,:));
        end
        axis tight
        if ~isempty(NOISEParams{kidn}.FitS1kHzEval)  %Only if fit has been made
            legend('data','Fit','Jiansong','-50dBm value')
        else
            legend('data')
        end
        xlabel('P_{int} (dBm)');ylabel('S_f/f^2 @1kHz [dB]')
        box on;grid on;
        title('Setup subtracted')
        
        
        %===============================
        %Pint (for P<=Popt)
        %===============================
        subplot(2,3,4) %Popt as a function of KID ID (and Q)
        [~,QCindex] = min(abs(log10(NOISE(IndexPopt(kidn)).Qc)-QCrange));%mapping Qc on the Qrange
        plot(NOISE(IndexPopt(kidn)).KIDnumber,NOISE(IndexPopt(kidn)).InternalPower,...
            'ko','MarkerSize',12,'LineWidth',2,'MarkerFaceColor',QCmap(QCindex,:))
        hold on
        axis tight
        caxis([3, 6]);
        colorbar('Ticks',QCtick,'TickLabels',{'10^3','10^4','10^5','10^6'});
        %colorbar('Ticks',QCtick,'TickLabels',tick2plot); %also works but looks less nice%
        xlabel('KID ID');ylabel('P^{int}_{opt} (dBm)')
        xlim(PR)
        title('color = Q');
        grid on;
        
        
        %===============================
        %Q's (for P=Popt)
        %===============================
        subplot(2,3,5) %Popt as a function of KID ID (and Q)
        semilogy(NOISE(IndexPopt(kidn)).KIDnumber,NOISE(IndexPopt(kidn)).Ql,'ko','MarkerSize',8,'MarkerFaceColor','k');
        hold on;
        semilogy(NOISE(IndexPopt(kidn)).KIDnumber,NOISE(IndexPopt(kidn)).Qi,'pb','MarkerSize',8,'MarkerFaceColor','b');
        semilogy(NOISE(IndexPopt(kidn)).KIDnumber,NOISE(IndexPopt(kidn)).Qc,'rs','MarkerSize',8,'MarkerFaceColor','r');
        xlabel('KID ID');ylabel('Q factors');
        grid on
        legend('Q','Qi','Qc');
        xlim(PR)
        
        %===============================
        %Sf's (for P=-50dBm)
        %===============================
        if ~isempty(NOISEParams{kidn}.FitS1kHz)
            subplot(2,3,6)
            plot(NOISE(IndexPopt(kidn)).KIDnumber,NOISEParams{kidn}.FitS1Hz.b,'go','MarkerSize',8,'MarkerFaceColor','g');
            hold on;
            plot(NOISE(IndexPopt(kidn)).KIDnumber,NOISEParams{kidn}.FitS10Hz.b,'ro','MarkerSize',8,'MarkerFaceColor','r');
            plot(NOISE(IndexPopt(kidn)).KIDnumber,NOISEParams{kidn}.FitS100Hz.b,'bo','MarkerSize',8,'MarkerFaceColor','b');
            plot(NOISE(IndexPopt(kidn)).KIDnumber,NOISEParams{kidn}.FitS1kHz.b,'ko','MarkerSize',8,'MarkerFaceColor','k');
            if kidn==length(KIDnumbers)
                plot(PR,[-189+5,-189+5],'k-','linewidth',2);
            end
            warning('off', 'MATLAB:legend:IgnoringExtraEntries' );
            legend('1 Hz','10 Hz','100 Hz','1kHz')
            xlabel('KID ID');ylabel('S_F/F^2 [dBk/Hz] @-50dBm internal power')
            title('Black JSG-line 1 kHz 120mK -50dBm')
            grid on
            xlim(PR)
        end
    end
    subplot(2,3,2)
    ylim([floor(min(minp)) ,ceil(max(maxnp))])
    xlim([5 0.5e5])
    legend(lstr);
    
    Figfile=[Noisepath,'All_KIDS_',...
        num2str(NOISE(IndexPref).Temperature,'%.3g'),'K'];
    if SaveStuff == 1
        MakeGoodFigure(20,13,12,Figfile)
    else
        MakeGoodFigure(20,13,12)
    end
    
    clear S21minMag tick2plot QCtick QCrange QCmap QCindex numlowPremoved lstr kidn kidcolors
end

rmpath([pwd,filesep,'..',filesep,'subroutines']);
if SaveStuff == 1
    save([Noisepath,'Noise_P.mat'])
end
