function TDanalysisV2
% reads td binary data, calulates and plots the PSD's
% saves the data into the NOISE struct and appends to the exiting .mat file 
%
clear variables;close all
% set by HW VI's in labview
freqmed = 50e3;              %standard is 50e3;
freqfast = 1e6;              %standard is 1e6;
timefast = 0.2;              %standard is 0.2;

%================================================================================
% Input
%================================================================================
ChipInfo_path = ['..' filesep '..']; %root path where data is, one higher than the scripts %ChipInfo_path = [cd filesep '..'];
Pdep = 1;           %=1 for Pdep analysis, =0 for Tdep analysis, error otheriwse
timemed = 120;       % duration of td data @ 50 ksample/sec; is varied so check the files, 40 is old default
TDoption.nrsigma = 8;                 % 8 peaks>TDoption.nrsigma*sigma are rejected (based on filtered TD-stream) %6 seems to be really minimum, 8 good value also for smaller peaks
TDoption.smoothtime = 1e-3;           % 1e-3quasiparticle lifetime guess used for fit (note it will be updated if there is a resultfile present)
average_med = 64;                       % 32,  (not more) Nr of TDoption.averages for PSD, 
                                        % should be a multiple of 2 AND nrpoints/TDoption.average/2 should be integer. 
                                      % Average is also the number of pieces over which peak rejection is done, ie if TDoption.average=32 and a peak is found, 
                                      % that 1/32th piece is thrown away. 
                                      % scipt uses TDp[tion.average, is set to 32  for the fast data and to average for med data%
TDoption.ppd = 30;                    % 30 points per decade in frequency for logsmooth after calculating PSD, 
                                      % 30 (20) will give spectrum of 158 (105) points. 0 will give you the full spectrum, usually a million points
                                      % 30 is used standard in the SRON labview    
                                      % NOTE: if you want to use fitranges from earlier processing, make sure ppd is the same!!!
TDoption.savex = 0;                   % 0 save filtered time domain data in csv-file (note some 90 Mb per file); Saves in data direcotory
TDoption.corroff = 1;                 % 1 = correct each piece of length/average with linear offset (to correct for very slow drifts, only relevant for 50 kHz files

average_fast = 32; %do not touch

if rem( freqmed * timemed / average_med , 1)~=0
    error('TDoption.average = wrong for med data')
end

%================================================================================

if Pdep == 1
    FFTsubsubdir=[ 'FFT' filesep 'Power'];                   %FFTsubdir = [filesep 'Noise_Powers_165mK' filesep 'FFT' filesep 'Power'];     %
    TDsubsubdir=[ 'TD_Power'];
    Pdependence = 1;
    matfile = 'Noise_P.mat';
    load([ChipInfo_path,filesep,FFTsubsubdir,filesep,matfile],'NOISE','IndexP_sub_opt','KIDnumbers','IndexPopt');
elseif Pdep ==0
    %Tdep
    FFTsubsubdir=[  'FFT' filesep 'Temp'];                   %FFTsubdir = [filesep 'Noise_Powers_165mK' filesep 'FFT' filesep 'Power'];     %
    TDsubsubdir=[   'TD_Temp'];
    Pdependence = 0;
    matfile = 'Noise_T.mat';
    load([ChipInfo_path,filesep,FFTsubsubdir,filesep,matfile],'NOISE','KIDnumbers');
else
    error('Pdep not defined')
end

addpath([pwd,filesep,'..',filesep,'subroutines']);
addpath([pwd,filesep,'..',filesep,'Time Domain analysis']);

ChipInfo.path = ChipInfo_path;clear ChipInfo_path;

TDpath = [ChipInfo.path,filesep,TDsubsubdir,filesep]; %end with filesep
disp(['TDpath: ',TDpath]);


for kidn=1:length(KIDnumbers) % LOOP OVER ALL UNIQUE KIDS,
    %construct filename
    if Pdependence == 1
        ID = num2str(NOISE(IndexPopt(kidn)).KIDnumber);
        for p=1:length(IndexP_sub_opt{kidn})% over Power
            % construct filename
            Pr = num2str(-1*NOISE(IndexP_sub_opt{kidn}(p)).ReadPower);
            Tk = num2str(round(1e3*NOISE(IndexP_sub_opt{kidn}(p)).Temperature),'%.0f');
            fnmed = ['KID' ID '_' Pr 'dBm__TDmed_TmK' Tk '.bin' ];  %bin filename - 50ksample/sec
            fnfast = ['KID' ID '_' Pr 'dBm__TDfast_TmK' Tk '.bin' ];%bin filename 2Msample/sec
            figname = [TDpath 'KID' ID '_' Pr 'dBm__TD_TmK' Tk];    %figure name for export
            %check if files are there
            fid = fopen([TDpath fnmed ]);fid2 = fopen([TDpath fnfast ]);
            if fid == -1 || fid2 == -1
                error(['ERROR import_data: Cannot find ' fnmed ])
            else
                disp(['Reading: ' fnmed])
            end
            fclose(fid);fclose(fid2);
            %reading data; figure name defined as binfilefunction4 makes
            %plots
            figure(100000*NOISE(IndexP_sub_opt{kidn}(p)).KIDnumber + abs(NOISE(IndexP_sub_opt{kidn}(p)).ReadPower));
            [fmfast,SRRrawfast,SPPrawfast,SPRrawfast]   = binfilefunction4([TDpath fnfast ],timefast,freqfast,average_fast,TDoption);
            [fmmed,SRRrawmed,SPPrawmed,SPRrawmed]       = binfilefunction4([TDpath fnmed ],timemed,freqmed,average_med,TDoption);
            %correcting and patching spectral data med to fast
            fmmed(end-5:end) = [];SRRrawmed(end-5:end) = [];SPPrawmed(end-5:end) = [];SPRrawmed(end-5:end) = [];
            fmfast(1:39) = [];SRRrawfast(1:39) = [];SPPrawfast(1:39) = [];SPRrawfast(1:39) = [];
            fm = [fmmed fmfast];SRRraw = [SRRrawmed SRRrawfast];SPPraw = [SPPrawmed SPPrawfast];SPRraw = [SPRrawmed SPRrawfast];
            %Put data in noise struct=
            NOISE(IndexP_sub_opt{kidn}(p)).CrossPSD(:,4) = SRRraw;
            NOISE(IndexP_sub_opt{kidn}(p)).CrossPSD(:,3) = SPPraw;
            NOISE(IndexP_sub_opt{kidn}(p)).CrossPSD(:,2) = SPRraw;
            NOISE(IndexP_sub_opt{kidn}(p)).CrossPSD(:,1) = fm;
            
            %===================================================================
            %Plotting (cntnd from binfilefunction4)
            %===================================================================
            
            subplot(2,2,3)%plot R and theta labview
            semilogx(NOISE(IndexP_sub_opt{kidn}(p)).FFTnoise{1}(:,1),NOISE(IndexP_sub_opt{kidn}(p)).FFTnoise{1}(:,2),...
                'r-','LineWidth',2);hold on;
            semilogx(fm,SPPraw,'b-','LineWidth',2);
            %add R
            semilogx(NOISE(IndexP_sub_opt{kidn}(p)).FFTnoise{1}(:,1),NOISE(IndexP_sub_opt{kidn}(p)).FFTnoise{1}(:,3),'r-');
            semilogx(fm,SRRraw,'b-');
            grid on;axis tight;xlim([10,0.3e6]);
            xlabel('F [Hz]');ylabel('S_x [dBc/Hz]');
            legend('\theta Labview','\theta Matlab','R Labview','R Matlab')
            
            subplot(2,2,4)
            semilogx(NOISE(IndexP_sub_opt{kidn}(p)).CrossPSD(:,1),-1*real(NOISE(IndexP_sub_opt{kidn}(p)).CrossPSD(:,2)),'b-');grid on;
            xlabel('F [Hz]');ylabel('S_{cross} [dBc/Hz]');
            grid on;axis tight;xlim([10,0.1e6]);
            title(['KID ' num2str(NOISE(IndexP_sub_opt{kidn}(p)).KIDnumber) ' @Pread = ' num2str(NOISE(IndexP_sub_opt{kidn}(p)).ReadPower)...
                ' dBm, T = ' num2str(NOISE(IndexP_sub_opt{kidn}(p)).Temperature)]);
            MakeGoodFigure(8,6,8,figname,1); %save png
            close(gcf)
        end
    elseif Pdependence == 0
        ID = num2str(NOISE(kidn).KIDnumber);
        for nT=1:length(NOISE(kidn).Temperature) % over T
            % construct filename
            Pr = num2str(-1*NOISE(kidn).ReadPower);
            Tk = num2str(round(1e3*NOISE(kidn).Temperature(nT)),'%.0f');
            fnmed = ['KID' ID '_' Pr 'dBm__TDmed_TmK' Tk '.bin' ];
            fnfast = ['KID' ID '_' Pr 'dBm__TDfast_TmK' Tk '.bin' ];
            figname = [TDpath 'KID' ID '_' Pr 'dBm__TD_TmK' Tk];
            %check if files are there
            fid = fopen([TDpath fnmed ]);fid2 = fopen([TDpath fnfast ]);
            if fid == -1 || fid2 == -1
                error(['ERROR import_data: Cannot find ' fnmed ])
            else
                disp(['Reading: ' fnmed])
            end
            fclose(fid);fclose(fid2);
            %start reading
            figure(100000*NOISE(kidn).KIDnumber + round(1000*NOISE(kidn).Temperature(nT)));
            [fmfast,SRRrawfast,SPPrawfast,SPRrawfast]   = binfilefunction4([TDpath fnfast ],timefast,freqfast,average_fast,TDoption);
            [fmmed,SRRrawmed,SPPrawmed,SPRrawmed]       = binfilefunction4([TDpath fnmed ],timemed,freqmed,average_med,TDoption);
            %correcting and patching data
            fmmed(end-5:end) = [];SRRrawmed(end-5:end) = [];SPPrawmed(end-5:end) = [];SPRrawmed(end-5:end) = [];
            fmfast(1:39) = [];SRRrawfast(1:39) = [];SPPrawfast(1:39) = [];SPRrawfast(1:39) = [];
            fm = [fmmed fmfast];SRRraw = [SRRrawmed SRRrawfast];SPPraw = [SPPrawmed SPPrawfast];SPRraw = [SPRrawmed SPRrawfast];
            %Put data in noise struct=
            NOISE(kidn).CrossPSD{nT}(:,4) = SRRraw;
            NOISE(kidn).CrossPSD{nT}(:,3) = SPPraw;
            NOISE(kidn).CrossPSD{nT}(:,2) = SPRraw;
            NOISE(kidn).CrossPSD{nT}(:,1) = fm;
            
            %===================================================================
            %Plotting (cntnd from binfilefunction4)
            %===================================================================
            
            subplot(2,2,3)%phase
            semilogx(NOISE(kidn).FFTnoise{nT}(:,1),NOISE(kidn).FFTnoise{nT}(:,2),...
                'r-','LineWidth',2);hold on;
            semilogx(fm,SPPraw,'b-','LineWidth',2);
            %add R
            semilogx(NOISE(kidn).FFTnoise{nT}(:,1),NOISE(kidn).FFTnoise{nT}(:,3),'r-');
            semilogx(fm,SRRraw,'b-');
            grid on;axis tight;xlim([10,0.3e6]);
            xlabel('F [Hz]');ylabel('S_x [dBc/Hz]');
            legend('\theta Labview','\theta Matlab','R Labview','R Matlab')
            
            subplot(2,2,4)
            semilogx(NOISE(kidn).CrossPSD{nT}(:,1),-1*real(NOISE(kidn).CrossPSD{nT}(:,2)),'b-');grid on;
            xlabel('F [Hz]');ylabel('S_{cross} [dBc/Hz]');
            grid on;axis tight;xlim([10,0.1e6]);
            title(['KID ' num2str(NOISE(kidn).KIDnumber) ' @Pread = ' num2str(NOISE(kidn).ReadPower)...
                ' dBm, T = ' num2str(NOISE(kidn).Temperature(nT))]);
            MakeGoodFigure(8,6,8,figname,1); %save png
            close(gcf);
        end
    end
end

save([ChipInfo.path,filesep,FFTsubsubdir,filesep,matfile],'NOISE','-append');