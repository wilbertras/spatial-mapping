% NOISEanalysis_TdepV7
close all
clear all

%sub-dir inside ChipPath where the Noise Data [FFT(T),S21(T),ft(T)] is located. Starts with \ ends without \.
% ChipInfo.path = [cd filesep '/..']; %root path where data is, one higher than the scripts
ChipInfo.path = [ '..' filesep '..'];%without filesep at end
FFTsubdir = [filesep 'Noise_Temperature' filesep 'FFT' filesep 'Temp'];     %
%Note: pwd is current matlab directory. Usually directory where this m-file
%is located.
%==========================================================================
% Internal variables which may be user modified.
%==========================================================================
global FITF0

FITF0=0;          %Switch used to determine which method is used for determination of F0.
%Please see the FitS21main5 routine for details. FITF0 = 0 is recommended due to the presence of overdriven
%resonators in Power sweeps (it will only update Q, not F0. F0 mot very important here).
SaveStuff = 1; 	%0 to not save files (faster) 1 to save
MaxnT = 24;     %max. #temperatures to be read.   
%==========================================================================
% Setting some routine default values.
%==========================================================================
global Noisepath
format('long','g'); 
addpath([pwd,filesep,'..',filesep,'subroutines']);
close all; 
%==========================================================================
% In the datapath find all files and identify KIDs.
%==========================================================================
Noisepath = [ChipInfo.path,FFTsubdir,filesep]; %Path containing the raw noise data

%Search the noisepath for all FFT files and filter out all files that are
%significantly longer or shorter than the mean. (This usually means
%something went wrong during measurements)
RawFFTfiles = dir([Noisepath,'KID*FFT*.dat']);
if isempty(RawFFTfiles)
    error([ 'No files found in ' [Noisepath,'KID*FFT*.dat'] ]);
end

% %check FFT size for compatibility
% datasize = zeros(1,length(RawFFTfiles));
% for p=1:length(RawFFTfiles)
%     datasize(p)=RawFFTfiles(p).bytes;
% end
% FFT2read = 0.5*mean(datasize) < datasize & datasize < 1.5*mean(datasize);
% RawFFTfiles = RawFFTfiles(FFT2read);
% clear datasize

%Determine the KID numbers (IDs) from the names of the good files.
KIDnumbers = zeros(1,length(RawFFTfiles));
NOISE(length(RawFFTfiles)).KIDnumbers = ''; %pre allocate
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
    clear LocEndRoot
    %======================================================================
    %One by one read in the data files
    %======================================================================
    
    %Reading FFT file
    %======================================================================
    FFTfile = [NOISE(p).filename,'FFT.dat'];
    [Data,Temperature,Power,FFTheader] = import_data(FFTfile);%Data = cell array
    % limit T range
    if p == 1
        disp('Temperatures:')
        disp(Temperature);
    end
    if length(Temperature) > MaxnT
        Temperature(MaxnT+1:end)=[];
    end
    MaxnT = length(Temperature);
    NOISE(p).Temperature = Temperature;%        
    Data(MaxnT+1:end)=[];
    NOISE(p).FFTnoise = Data;
    clear Data;
    NOISE(p).ReadPower = -1*Power;
    
    %Find in FFTheader the frequency of the tone used to read the noise (at nT=1).
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
    clear Data;
    
    %======================================================================
    %Reading S21 file (GHz,Re,Im)
    S21file = [RootName,'S21.dat'];
    [Data,~,Power,~] = import_data(S21file);
    %Check for equality of temperature and power
    if abs(-1*Power-NOISE(p).ReadPower(1)) > 0.5
        disp('Warning: Read Power does not match between FFT and S21 file.')
    end
    
    %Copy data for struct
    Data(MaxnT+1:end)=[];
    NOISE(p).S21_IQplane = Data;
    clear Data S21file;
    %======================================================================
    %Reading S21 file (GHz,dB,rad)
    S21DBfile = [RootName,'S21dB.dat'];
    [Data,~,Power] = import_data(S21DBfile);
    %Check for equality of temperature and power
    if abs(-1*Power-NOISE(p).ReadPower(1)) > 0.5
        disp('Warning: Read Power does not match between FFT and S21dB file.')
    end
    %Copy data for struct
    Data(MaxnT+1:end)=[];
    NOISE(p).S21_MPplane = Data;
    clear Data S21DBfile
    %======================================================================
    %Reading time-domain file
    TDfile = [RootName,'td.dat'];
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
    Data(MaxnT+1:end)=[];
    for nT=1:MaxnT
        CompleteData = zeros(size(Data{1},1),size(Data{1},2)+1);
        CompleteData(:,1) = dt*(1:size(Data{nT},1))';
        CompleteData(:,2:3) = Data{nT}(:,:);
        NOISE(p).TDIQ{nT,1} = CompleteData;
    end
    disp(['Highest Read: ' num2str(NOISE(p).Temperature(end))]);
    clear Temperature Data dt TDfile CompleteData fid TDHeader
    
    %======================================================================
    %Analyse the S21 data measured in [mag,phase] space for each temperature
    %======================================================================
    
    for nT=1:length(NOISE(p).Temperature)
        %Normalize the S21 data measured in magnitude plane
        filtered=smooth(NOISE(p).S21_MPplane{1}(:,2),3); %smoothing the |S21| data in dB space
        %normalise S21 in log space to the max(|S21|)
        NOISE(p).S21_MPplane{nT}(:,2)=NOISE(p).S21_MPplane{nT}(:,2)-max(filtered);
        %Convert dB to magnitude
        NOISE(p).S21_MPplane{nT}(:,2) = 10.^(NOISE(p).S21_MPplane{nT}(:,2)/20);
        %Perform fit to obtain resonator parameters. Note the FITF0 is
        %recommended for overdriven resonators.
        [fres,Q,S21min,FitResult] = FitS21main5(NOISE(p).S21_MPplane{nT}(:,1:3),FITF0);
        %Put the temporary storage variables into the NOISE struct
        NOISE(p).Ql(nT)=Q(1);
        NOISE(p).Qi(nT)=Q(2);
        NOISE(p).Qc(nT)=Q(3);
        NOISE(p).Fres(nT)=fres;
        NOISE(p).S21min(nT)=S21min;%in dB!
        NOISE(p).S21fit{nT,1} = FitResult;
        %Calculate the Internal Power
        NOISE(p).InternalPower(nT) = 10*log10((2/pi)*10.^(NOISE(p).ReadPower/10).*(NOISE(p).Ql(nT).^2./NOISE(p).Qc(nT)));
    end
    %Calculate the resonator ring time
    NOISE(p).TauRes = NOISE(p).Ql(nT)./(pi*NOISE(p).Fres);
    %Calculate the resonator bandwidth (used later in Popt determination)
    NOISE(p).Bandwidth = NOISE(p).Fres./NOISE(p).Ql;
    NOISE(p).MeanNoise = zeros(MaxnT,3);
    NOISE(p).MeanFreqNoise = zeros(MaxnT,4);
    
    for nT=1:length(NOISE(p).Temperature)
        %Calculate frequency noise
        %Normalized Frequency Noise (Sf/F^2) [1/Hz] (defined by J. Gao) - setup
        %corrected
        tempnoise = (10.^(NOISE(p).FFTnoise{nT}(:,2)/10))*(1/(4*NOISE(p).Ql(nT))^2);
        Setup_level = mean(tempnoise(end-10:end-5));
        NOISE(p).FFTnoise{nT}(:,4) = tempnoise - Setup_level;
        %NOISE(p).FFTnoise{1}(NOISE(p).FFTnoise{1}(:,4)<0,4) = NaN;
        clear tempnoise Setup_level fres Q S21min FitResult
        
        %Frequency Noise (defined by B. Mazin) not used
        %NOISE(p).FFTnoise{nT}(:,5) = NOISE(p).FFTnoise{nT}(:,4)*(NOISE(p).Fres(nT)*1e9)^2;
        
        %Phase noise WRT complex plane
        S21min_a=10^(( NOISE(p).S21min(nT) )/20);%from dB to magnitude, not stored
        R=20*log10((1-S21min_a)/2);%correction due to circl radius vs complex plane radius in dB
        NOISE(p).FFTnoise{nT}(:,6)=NOISE(p).FFTnoise{nT}(:,2)+R;
        
        %Amplitude noise WRT complex plane
        NOISE(p).FFTnoise{nT}(:,7)=NOISE(p).FFTnoise{nT}(:,3)+R;
        
        %Calculate Mean amplitude noise between 20 and 500 Hz (Required for
        %Popt determination)
        Frange = 20 <= NOISE(p).FFTnoise{nT}(:,1) & NOISE(p).FFTnoise{nT}(:,1) <= 500;
        NOISE(p).MeanNoiseforPopt(nT) = mean(NOISE(p).FFTnoise{nT}(Frange,3));
        
        %Calculate the mean frequency noise (Gao definition), with setup
        %contribution subtracted (was done already in defining Gao's freq
        %noise)
        
        NOISE(p).MeanFreqNoise1Hz(nT) = ... % 1 Hz
            mean(interp1(NOISE(p).FFTnoise{(nT)}(:,1),NOISE(p).FFTnoise{(nT)}(:,4),(0.8:0.1:1.2),'linear',1));
        NOISE(p).MeanFreqNoise10Hz(nT) = ... % 10 Hz
            mean(interp1(NOISE(p).FFTnoise{(nT)}(:,1),NOISE(p).FFTnoise{(nT)}(:,4),(9:0.1:11),'linear',1));
        NOISE(p).MeanFreqNoise100Hz(nT) = ... % 100 Hz
            mean(interp1(NOISE(p).FFTnoise{(nT)}(:,1),NOISE(p).FFTnoise{(nT)}(:,4),(90:1:110),'linear',1));
        NOISE(p).MeanFreqNoise1000Hz(nT) = ... % 1 kHz
            mean(interp1(NOISE(p).FFTnoise{(nT)}(:,1),NOISE(p).FFTnoise{(nT)}(:,4),(900:10:1100),'linear',1));
    end
    
end %END OF LOOP OVER ALL FILES (aka KID-Pread-combinations)
clear Power R S21min_a Data Frange RawFFTfiles FITF0 Freadstart Freadindex hl FFT2read ColonIndex filtered nT p FFTheader FFTsubdir FFTfile MaxnT
%==========================================================================
% Now that all files have been read we plot
%==========================================================================
for kidn=1:length(KIDnumbers) % LOOP OVER ALL UNIQUE KIDS,
    %==================================================================
    % Figure one: mainly checks for the quality of the analysis
    % routine
    %==================================================================
    figure(1000*NOISE(kidn).KIDnumber(1)+1)
    clf
    
    Tcolors = colormapJetJB(length(NOISE(kidn).Temperature));
    clear TempLegend
    TempLegend = cell(length(NOISE(kidn).Temperature)+1,1);
    TempLegend{1} = 'P_{opt}';
    
    %Resonance circle as a function of T (incl noise blobs)
    subplot(2,3,1)
    hold on
    for nT=1:length(NOISE(kidn).Temperature)
        plot(NOISE(kidn).S21_IQplane{nT}(:,2) , NOISE(kidn).S21_IQplane{nT}(:,3),'-','color',Tcolors(nT,:),'LineWidth',1)
        TempLegend{nT+1} = num2str(NOISE(kidn).Temperature(nT),'%.3f');
    end
    for nT=1:length(NOISE(kidn).Temperature) %Second Loop to get the legend correct.
        plot(NOISE(kidn).TDIQ{nT}(:,2),NOISE(kidn).TDIQ{nT}(:,3),'.','color',Tcolors(nT,:),'MarkerSize',6)
    end
    legend(TempLegend(2:end))
    xlabel('Re');ylabel('Im')
    title(['KID ',num2str(NOISE(kidn).KIDnumber(1),'%.0f'),' @P_{read}=',num2str(NOISE(kidn).ReadPower,'%.0f'),' dBm'])
    box on;grid on;
    hold off
    
    %Resonance Dip as a function of T, Incl reference lines around reference power%
    subplot(2,3,2)
    hold on
    for nT=1:length(NOISE(kidn).Temperature)
        plot(NOISE(kidn).S21_MPplane{nT}(:,1),20*log10(NOISE(kidn).S21_MPplane{nT}(:,2)),'-','color',Tcolors(nT,:),'LineWidth',1)
    end
    xlabel('F [GHz]');ylabel('|S21| [dB]')
    title(['KID ',num2str(NOISE(kidn).KIDnumber(1),'%.0f'),' P_{opt}=',num2str(NOISE(kidn).ReadPower),' dBm'])
    axis tight;
    box on;grid on;
    hold off
    
    %Time Domain Trace as a function of time
    subplot(2,3,3)
    for nT=1:length(NOISE(kidn).Temperature)
        plot(NOISE(kidn).TDIQ{nT}(:,1),NOISE(kidn).TDIQ{nT}(:,2),...
            'o','color',Tcolors(nT,:),'MarkerSize',3,'MarkerFaceColor',Tcolors(nT,:));
        hold on
        plot(NOISE(kidn).TDIQ{nT}(:,1),NOISE(kidn).TDIQ{nT}(:,3),...
            '.','color',Tcolors(nT,:),'MarkerSize',3);
    end
    xlabel('t [sec]'); ylabel('Re or Im')
    legend('Re(S21)','Im(S21)');
    box on;
    hold off
    
    %Frequency Noise
    subplot(2,3,4)
    warning('off', 'MATLAB:plot:IgnoreImaginaryXYPart');
    toplot = NOISE(kidn).FFTnoise{1}(:,4) > 0;
    semilogx(NOISE(kidn).FFTnoise{1}(toplot,1),10*log10(NOISE(kidn).FFTnoise{1}(toplot,4)),...
        '-','color','k','LineWidth',2)
    hold on
    for nT=1:length(NOISE(kidn).Temperature)
        toplot = NOISE(kidn).FFTnoise{nT}(:,4) > 0;
        semilogx(NOISE(kidn).FFTnoise{nT}(toplot,1),10*log10(NOISE(kidn).FFTnoise{nT}(toplot,4)),...
            '-','color',Tcolors(nT,:),'LineWidth',1)
    end
    xlabel('F [Hz]');ylabel('S_F/F^2 [dBc/Hz]')
    xlim([0.5,1e5]);grid on;ylim([-220,-140])
    hold off
    
    %Phase Noise (and Amp Noise)
    subplot(2,3,[5,6])
    semilogx(NOISE(kidn).FFTnoise{1}(:,1),NOISE(kidn).FFTnoise{1}(:,3),...
        '-','color','k','LineWidth',3)
    hold on
    semilogx(NOISE(kidn).FFTnoise{1}(:,1),NOISE(kidn).FFTnoise{1}(:,2),...
        '-','color','k','LineWidth',4)
    minnp = zeros(length(NOISE(kidn).Temperature),1);maxnp = minnp;
    for nT=1:length(NOISE(kidn).Temperature)
        % finding range
        maxnp(nT) = max(NOISE(kidn).FFTnoise{nT}(NOISE(kidn).FFTnoise{nT}(:,1) > 10,2));
        minnp(nT) = mean(NOISE(kidn).FFTnoise{nT}(20:end-5,3));
        semilogx(NOISE(kidn).FFTnoise{nT}(:,1),NOISE(kidn).FFTnoise{nT}(:,3),...
            '-','color',Tcolors(nT,:),'LineWidth',1)
        semilogx(NOISE(kidn).FFTnoise{nT}(:,1),NOISE(kidn).FFTnoise{nT}(:,2),...
            '-','color',Tcolors(nT,:),'LineWidth',2)
    end
    grid on;
    xlabel('F [Hz]');ylabel('S_x [dBc/Hz]')
    legend('S_A','S_{\theta}')
    xlim([3,0.3e6]);ylim([10*floor(min(0.1*minnp))-5,10*ceil(max(0.1*maxnp))+5]);
    hold off
    
    %SAVE the figure
    Figfile=[Noisepath,'KID',num2str(NOISE(kidn).KIDnumber(1),'%.0f'),'_P_R_',...
        num2str(NOISE(kidn).ReadPower,'%.3g'),'dBm_NOISE_vsT1'];
    if SaveStuff == 1
        MakeGoodFigure(15,15,14,Figfile)
    else
        MakeGoodFigure(15,15,14)
    end
    
    %==================================================================
    % Figure two: Overview of resonator parameters as a function of
    % power.
    %==================================================================
    figure(1000*NOISE(kidn).KIDnumber(1)+2)
    clf
    %Frequency Noise at 1 kHz as a function of Pint
    subplot(2,2,1)
    warning('off','MATLAB:Axes:NegativeDataInLogAxis');
    %Measurement Points
    semilogy(NOISE(kidn).Temperature,NOISE(kidn).MeanFreqNoise1000Hz,'ro','MarkerSize',6,'MarkerFaceColor','r') ;
    title(['KID ',num2str(NOISE(kidn).KIDnumber(1),'%.0f')])
    xlabel('T  (K)');ylabel('S_f/f^2 @1kHz (dB)')
    box on;grid on;
    hold off
    
    subplot(2,2,2) %Resonance Frequency as a function of Pint
    hold on
    semilogy(NOISE(kidn).Temperature,NOISE(kidn).Fres,'ro','MarkerSize',6,'MarkerFaceColor','r') %all P
    xlabel('T  (K)');ylabel('f_{res} (GHz)')
    box on;grid on;
    hold off
    
    subplot(2,2,3) %Qi as a function of Pint
    semilogy(NOISE(kidn).Temperature,NOISE(kidn).Qi,'ro','MarkerSize',6,'MarkerFaceColor','r');%At Popt
    xlabel('T  (K)');ylabel('Q_{i} (dBm)')
    box on;grid on;
    hold off
    
    subplot(2,2,4) %S21 resonance fit at Tbase
    plot(NOISE(kidn).S21_MPplane{1}(:,1),NOISE(kidn).S21_MPplane{1}(:,2),'b.','MarkerSize',5)
    axis tight;hold on
    %Fit
    plot(NOISE(kidn).S21fit{1}(:,1),10.^(NOISE(kidn).S21fit{1}(:,2)/20),'r-','LineWidth',1)
    plot(NOISE(kidn).Fres(1),10.^(NOISE(kidn).S21min(1)/20),'kd','MarkerFaceColor','r','MarkerSize',6)%The determined Fres and |S21(Fres)|
    %make the plot nice
    xlabel('F (GHz)')
    ylabel('|S21|')
    title('Resonance @Popt')
    hold off
    
    %SAVE the figure
    Figfile=[Noisepath,'KID',num2str(NOISE(kidn).KIDnumber(1),'%.0f'),'_P_R_',...
        num2str(NOISE(kidn).ReadPower,'%.3g'),'dBm_NOISE_vsT2'];
    if SaveStuff == 1
        MakeGoodFigure(13,9,14,Figfile)
    else
        MakeGoodFigure(13,9,14)
    end
    
end %END OF LOOP OVER ALL UNIQUE KIDS

PoutData = zeros(length(KIDnumbers),9);
for kidn=1:length(KIDnumbers)
    PoutData(kidn,1) = NOISE(kidn).Temperature(1);
    PoutData(kidn,2) = NOISE(kidn).KIDnumber(1);
    PoutData(kidn,3) = NOISE(kidn).ReadPower(1);
    PoutData(kidn,4) = round(NOISE(kidn).InternalPower(1));
    PoutData(kidn,5) = NOISE(kidn).Fres(1);
    PoutData(kidn,6) = 0.1*round(NOISE(kidn).S21min(1)*10);
    PoutData(kidn,7) = round(NOISE(kidn).Ql(1)/1e3);
    PoutData(kidn,8) = round(NOISE(kidn).Qc(1)/1e3);
    PoutData(kidn,9) = round(NOISE(kidn).Qi(1)/1e3);
end


PoutHeader = {'T_{min}','KIDID','Pread (dBm)','Pint (dBm)','Fres (GHz)','S21 (dB)','Q','Qc','Qi'}';
WriteSRONcsv([Noisepath,'Pout.csv'],PoutData,PoutHeader,'%.6g')

clear PoptData Pcolors BWrange p PoptHeader kidn Chipinfo PoutData

%==========================================================================
% Wrap up and routine closure
%==========================================================================
clear Figfile TempLegend Tcolors RootName maxnp minnp kidn  Fread 

rmpath([pwd,filesep,'..',filesep,'subroutines']);
if SaveStuff == 1
    save([Noisepath,'Noise_T.mat'])
end
