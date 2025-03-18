function CrossPSDJBV1
% reads td binary data, calulates and plots the PSD's
% saves the data into the relevant struct
%no further analysis
%
clear variables;close all

%================================================================================
% Input
%================================================================================
ChipInfo_path = [ '..' filesep '..'];%without filesep at end %ChipInfo_path = [cd filesep '..'];
Pdep = 1;           %=1 for Pdep analysis, =0 for Tdep analysis, error otheriwse
lowFupper = 90;     %upper bound of low F range (level determined by quasiparticles) in Hz
lowFlower = 40;     %upper bound of low F range (level determined by quasiparticles)
nsigfitc = 1.0;   	%1.2 # sigma in the lowF noise that is condsidered enough to start a fit
plotall = 0;        %if =1 a plot of each fit is made and stored as png
tauinmax = 1.5e-3;    %maximum possible lifetime
highFpt = 1.3e4;    %1.3e4: high frequency reference pt for the setiup level for P dependence. ange will be 0.7 .. 1x this number%
%================================================================================

if Pdep == 1
    %FFTsubsubdir=['FFT' filesep 'Power'];                   %no filesep, name of the noise data in /FFT (== TD data dir name)
    FFTsubsubdir = [ 'FFT' filesep 'Power'];
    Pdependence = 1;
    matfile = 'Noise_P.mat';
elseif Pdep ==0
    %Tdep
    %FFTsubsubdir=['FFT' filesep 'Temp'];                   %no filesep, name of the noise data in /FFT (== TD data dir name)
    FFTsubsubdir = [filesep 'Noise_Temperature' filesep 'FFT' filesep 'Temp']; 
    Pdependence = 0;
    matfile = 'Noise_T.mat';
else
    error('Pdep not defined')
end

addpath([pwd,filesep,'..',filesep,'subroutines']);

%load relevant data
if Pdependence == 1
    load([ChipInfo_path,filesep,FFTsubsubdir,filesep,matfile],'NOISE','IndexP_sub_opt','KIDnumbers');
else
    load([ChipInfo_path,filesep,FFTsubsubdir,filesep,matfile],'NOISE','KIDnumbers');
end
ChipInfo.path = ChipInfo_path;clear ChipInfo_path; %needed to keep struct ok wrt path names when data is moved in between analysis ruyns
warning('off', 'MATLAB:Axes:NegativeDataInLogAxis')

for kidn=1:length(KIDnumbers) % LOOP OVER ALL UNIQUE KIDS,
    %construct filename
    legc = 1;%legend counter
    clear lstr
    if Pdependence == 1
        kleur = colormapJetJB(length(IndexP_sub_opt{kidn}));
        for p=1:length(IndexP_sub_opt{kidn})% over Power
            %quick variables
            CPSD = -1*real(NOISE(IndexP_sub_opt{kidn}(p)).CrossPSD(:,2));
            FCPSD = NOISE(IndexP_sub_opt{kidn}(p)).CrossPSD(:,1);
            %auto condition to allow for fit
            lowF = mean(CPSD(lowFlower <= FCPSD & FCPSD <= lowFupper));
            s_lowF = std(CPSD(lowFlower < FCPSD & FCPSD <= lowFupper));
            highF = mean(CPSD(FCPSD > 0.7*highFpt & FCPSD < highFpt));
            %s_highF = std(CPSD(end-10:end));%ignored as being much smaller
            if plotall == 1
                %plot data for each
                figure(100000*NOISE(IndexP_sub_opt{kidn}(p)).KIDnumber + abs(round(1000*NOISE(IndexP_sub_opt{kidn}(p)).ReadPower)));
                subplot(1,2,1)
                semilogx(FCPSD,CPSD,'b-');hold on;
                semilogx([lowFlower lowFupper],[lowF lowF],'r-','LineWidth',2);
                semilogx([lowFlower lowFupper],[lowF lowF]-s_lowF,'r--','LineWidth',1);semilogx([lowFlower lowFupper],[lowF lowF]+s_lowF,'r--','LineWidth',1);
                semilogx([1e4 2e4],[highF highF],'g-','LineWidth',4);
                subplot(1,2,2)
                loglog(FCPSD,CPSD,'b-');hold on;
                loglog([lowFlower lowFupper],[lowF lowF],'r-','LineWidth',2);
                loglog([lowFlower lowFupper],[lowF lowF]-s_lowF,'r--','LineWidth',1);semilogx([lowFlower lowFupper],[lowF lowF]+s_lowF,'r--','LineWidth',1);
                loglog([FCPSD(end-10) FCPSD(end)],[highF highF],'g-','LineWidth',4);
            end
            % see if we can fit and then do so
            if lowF - highF > nsigfitc*s_lowF
                % we fit starting at the low F range
                [tau,level,setupnoise,taumin,taumax] = crossfit(FCPSD(FCPSD > lowFlower),CPSD(FCPSD > lowFlower),lowF,tauinmax);
                NOISE(IndexP_sub_opt{kidn}(p)).CrossPSD(:,5) = level./(1 + (2*pi*tau*FCPSD).^2)+ setupnoise;
                if plotall == 1
                    subplot(1,2,1)
                    semilogx(FCPSD,NOISE(IndexP_sub_opt{kidn}(p)).CrossPSD(:,5),'k-');
                    legend('Data','GR level','GR-\sigma','GR+\sigma','highF level',['Fit, tau = ' num2str(tau*1000,'%0.3f') 'msec'])
                    subplot(1,2,2)
                    loglog(FCPSD,NOISE(IndexP_sub_opt{kidn}(p)).CrossPSD(:,5),'k-');
                end
            else
                tau = NaN;
                level = 0;
                setupnoise = 0;
                taumin = 0;taumax = 0;
            end
            
            %store results
            NOISE(IndexP_sub_opt{kidn}(p)).tau = tau;
            NOISE(IndexP_sub_opt{kidn}(p)).taumin = taumin;NOISE(IndexP_sub_opt{kidn}(p)).taumax = taumax;
            NOISE(IndexP_sub_opt{kidn}(p)).crosslevel = level;
            NOISE(IndexP_sub_opt{kidn}(p)).setupnoise = setupnoise;
            if plotall == 1
                %finish plot
                subplot(1,2,1)
                xlabel('F [Hz]');ylabel('S_{cross} [dBc/Hz]');
                grid on;axis tight;
                xlim([10,0.5e6]);
                subplot(1,2,2)
                xlabel('F [Hz]');ylabel('S_{cross} [dBc/Hz]');
                grid on;axis tight;
                xlim([10,0.5e6]);
                title(['KID ' num2str(NOISE(IndexP_sub_opt{kidn}(p)).KIDnumber) ' @Pread = ' num2str(NOISE(IndexP_sub_opt{kidn}(p)).ReadPower)...
                    ' dBm, T = ' num2str(NOISE(IndexP_sub_opt{kidn}(p)).Temperature)]);
                % save figure
                Figfile = [ChipInfo.path,filesep,FFTsubsubdir,filesep,...
                    'KID',num2str(NOISE(IndexP_sub_opt{kidn}(p)).KIDnumber,'%.0f'),'_',num2str(NOISE(IndexP_sub_opt{kidn}(p)).ReadPower,'%.3g'),'dBm_',...
                    num2str(round(1000*NOISE(IndexP_sub_opt{kidn}(p)).Temperature)) 'mK_CrossPSD'];
                MakeGoodFigure(14,6,12,Figfile,1); %save png
                close gcf
            end
            
            %nice final figure
            figure(IndexP_sub_opt{kidn}(1))
            subplot(1,3,1) %Cross PSD linear
            semilogx(NOISE(IndexP_sub_opt{kidn}(p)).CrossPSD(:,1),-1*real(NOISE(IndexP_sub_opt{kidn}(p)).CrossPSD(:,2)),'-',...
                'linewidth',2,'color',kleur(p,:));hold on;
            if ~isnan(NOISE(IndexP_sub_opt{kidn}(p)).tau)%add fit
                semilogx(NOISE(IndexP_sub_opt{kidn}(p)).CrossPSD(:,1),NOISE(IndexP_sub_opt{kidn}(p)).CrossPSD(:,5),'-k');
                lstr{legc}=['P_{int} = ' num2str(round(NOISE(IndexP_sub_opt{kidn}(p)).InternalPower)) ' dBm'];
                lstr{legc+1}=['\tau = ' num2str(round(1e6*NOISE(IndexP_sub_opt{kidn}(p)).tau)/1e3) ' msec'];
                legc = legc+2;%legend counter
            else
                lstr{legc}=['P_{int} = ' num2str(round(NOISE(IndexP_sub_opt{kidn}(p)).InternalPower)) ' dBm'];
                legc = legc+1;%legend counter
            end
            
            subplot(1,3,2) %Cross PSD log
            loglog(NOISE(IndexP_sub_opt{kidn}(p)).CrossPSD(:,1),-1*real(NOISE(IndexP_sub_opt{kidn}(p)).CrossPSD(:,2)),'-',...
                'linewidth',2,'color',kleur(p,:));hold on;
            if ~isnan(NOISE(IndexP_sub_opt{kidn}(p)).tau)%add fit
                loglog(NOISE(IndexP_sub_opt{kidn}(p)).CrossPSD(:,1),NOISE(IndexP_sub_opt{kidn}(p)).CrossPSD(:,5),'-k');
            end
            
            subplot(1,3,3)%lifetime vs P
            semilogy(NOISE(IndexP_sub_opt{kidn}(p)).InternalPower,NOISE(IndexP_sub_opt{kidn}(p)).tau*1e3,'o',...
                'color',kleur(p,:),'MarkerFaceColor',kleur(p,:));hold on
            semilogy(NOISE(IndexP_sub_opt{kidn}(p)).InternalPower,NOISE(IndexP_sub_opt{kidn}(p)).taumin*1e3,'o','color',kleur(p,:));
            semilogy(NOISE(IndexP_sub_opt{kidn}(p)).InternalPower,NOISE(IndexP_sub_opt{kidn}(p)).taumax*1e3,'o','color',kleur(p,:));
            
        end
        subplot(1,3,1)
        title(['KID ' num2str(NOISE(IndexP_sub_opt{kidn}(1)).KIDnumber) ' @ T = ' num2str(NOISE(IndexP_sub_opt{kidn}(1)).Temperature) ' K']);
        xlabel('F [Hz]');ylabel('S_{cross} [dBc/Hz]');
        grid on;axis tight;
        xlim([10,0.5e6]);
        legend(lstr);
        
        subplot(1,3,2)
        xlabel('F [Hz]');ylabel('S_{cross} [dBc/Hz]');
        grid on;axis tight;
        xlim([10,1e5]);
        
        subplot(1,3,3)
        xlabel('P_{internal} (dBm)');ylabel('\tau (msec)');grid on;ylim([0.01 tauinmax*1000]);
        Figfile = [ChipInfo.path,filesep,FFTsubsubdir,filesep,...
            'KID',num2str(num2str(NOISE(IndexP_sub_opt{kidn}(1)).KIDnumber),'%.0f'),'_crossPSD_Tdep'];
        MakeGoodFigure(18,6,12,Figfile); %save png and fig
        
        
        %==================================================================================================================================
        %==================================================================================================================================
    elseif Pdependence == 0 %Tdep
        %==================================================================================================================================
        %==================================================================================================================================
        kleur = colormapJetJB(length(NOISE(kidn).Temperature));
        for nT=1:length(NOISE(kidn).Temperature) % over T
            
            %quick variables
            CPSD = -1*real(NOISE(kidn).CrossPSD{nT}(:,2));
            FCPSD = NOISE(kidn).CrossPSD{nT}(:,1);
            %auto condition to allow for fit
            lowF = mean(CPSD(lowFlower <= FCPSD & FCPSD <= lowFupper));
            s_lowF = std(CPSD(lowFlower < FCPSD & FCPSD <= lowFupper));
            highF = mean(CPSD(end-10:end));
            %s_highF = std(CPSD(end-10:end));%ignored as being much smaller
            if plotall == 1
                %plot data for each
                figure(100000*NOISE(kidn).KIDnumber + round(1000*NOISE(kidn).Temperature(nT)));
                subplot(1,2,1)
                semilogx(FCPSD,CPSD,'b-');hold on;
                semilogx([lowFlower lowFupper],[lowF lowF],'r-','LineWidth',2);
                semilogx([lowFlower lowFupper],[lowF lowF]-s_lowF,'r--','LineWidth',1);semilogx([lowFlower lowFupper],[lowF lowF]+s_lowF,'r--','LineWidth',1);
                semilogx([FCPSD(end-10) FCPSD(end)],[highF highF],'g-','LineWidth',4);
                subplot(1,2,2)
                loglog(FCPSD,CPSD,'b-');hold on;
                loglog([lowFlower lowFupper],[lowF lowF],'r-','LineWidth',2);
                loglog([lowFlower lowFupper],[lowF lowF]-s_lowF,'r--','LineWidth',1);semilogx([lowFlower lowFupper],[lowF lowF]+s_lowF,'r--','LineWidth',1);
                loglog([FCPSD(end-10) FCPSD(end)],[highF highF],'g-','LineWidth',4);
            end
            % see if we can fit and then do so
            if lowF - highF > nsigfitc*s_lowF
                % we fit starting at the low F range
                [tau,level,setupnoise,taumin,taumax] = crossfit(FCPSD(FCPSD > lowFlower),CPSD(FCPSD > lowFlower),lowF,tauinmax);
                NOISE(kidn).CrossPSD{nT}(:,5) = level./(1 + (2*pi*tau*FCPSD).^2)+ setupnoise;
                if plotall == 1
                    subplot(1,2,1)
                    semilogx(FCPSD,NOISE(kidn).CrossPSD{nT}(:,5),'k-');
                    legend('Data','GR level','GR-range','GR+range','highF level',['Fit, tau = ' num2str(tau*1000,'%0.3f') 'msec'])
                    subplot(1,2,2)
                    loglog(FCPSD,NOISE(IndexP_sub_opt{kidn}(p)).CrossPSD(:,5),'k-');
                end
            else
                tau = NaN;
                level = 0;
                setupnoise = 0;
                taumin = 0;taumax = 0;
            end
            
            %store results
            NOISE(kidn).tau(nT) = tau;
            NOISE(kidn).taumin(nT) = taumin;NOISE(kidn).taumax(nT) = taumax;
            NOISE(kidn).crosslevel(nT) = level;
            NOISE(kidn).setupnoise(nT) = setupnoise;
            if plotall == 1
                %finish plot
                subplot(1,2,1)
                xlabel('F [Hz]');ylabel('S_{cross} [dBc/Hz]');
                grid on;axis tight;
                xlim([10,0.5e6]);
                subplot(1,2,2)
                xlabel('F [Hz]');ylabel('S_{cross} [dBc/Hz]');
                grid on;axis tight;
                xlim([10,0.5e6]);
                title(['KID ' num2str(NOISE(kidn).KIDnumber) ' @Pread = ' num2str(NOISE(kidn).ReadPower)...
                    ' dBm, T = ' num2str(NOISE(kidn).Temperature(nT))]);
                % save figure
                Figfile = [ChipInfo.path,filesep,FFTsubsubdir,filesep,...
                    'KID',num2str(NOISE(kidn).KIDnumber,'%.0f'),'_',num2str(NOISE(kidn).ReadPower,'%.3g'),'dBm_',...
                    num2str(round(1000*NOISE(kidn).Temperature(nT))) 'mK_CrossPSD'];
                MakeGoodFigure(14,6,12,Figfile,1); %save png
                close gcf
            end
            
            %nice final figure
            figure(kidn)
            subplot(1,3,1)%noise + fit lin
            semilogx(NOISE(kidn).CrossPSD{nT}(:,1),-1*real(NOISE(kidn).CrossPSD{nT}(:,2)),'-',...
                'linewidth',2,'color',kleur(nT,:));hold on;
            if ~isnan(NOISE(kidn).tau(nT))
                semilogx(NOISE(kidn).CrossPSD{nT}(:,1),NOISE(kidn).CrossPSD{nT}(:,5),'-k');
                lstr{legc}=['T = ' num2str(round(1000*NOISE(kidn).Temperature(nT))) 'mK'];
                lstr{legc+1}=['\tau = ' num2str(round(1e6*NOISE(kidn).tau(nT))/1e3) 'msec'];
                legc = legc + 2;
            else
                lstr{legc}=['T = ' num2str(round(1000*NOISE(kidn).Temperature(nT))) 'mK'];
                legc = legc + 1;
            end
            
            subplot(1,3,2)%noise + fit log
            loglog(NOISE(kidn).CrossPSD{nT}(:,1),-1*real(NOISE(kidn).CrossPSD{nT}(:,2)),'-',...
                'linewidth',2,'color',kleur(nT,:));hold on;
            if ~isnan(NOISE(kidn).tau(nT))
                loglog(NOISE(kidn).CrossPSD{nT}(:,1),NOISE(kidn).CrossPSD{nT}(:,5),'-k');
            else
                lstr{legc}=['T = ' num2str(round(1000*NOISE(kidn).Temperature(nT))) 'mK'];
            end
            
            subplot(1,3,2)
            xlabel('F [Hz]');ylabel('S_{cross} [dBc/Hz]');
            grid on;axis tight;
            xlim([10,1e5]);
            
            subplot(1,3,3)
            semilogy(NOISE(kidn).Temperature(nT),NOISE(kidn).tau(nT)*1e3,'o',...
                'color',kleur(nT,:),'MarkerFaceColor',kleur(nT,:));hold on
            semilogy(NOISE(kidn).Temperature(nT),NOISE(kidn).taumin(nT)*1e3,'o','color',kleur(nT,:));
            semilogy(NOISE(kidn).Temperature(nT),NOISE(kidn).taumax(nT)*1e3,'o','color',kleur(nT,:));
        end
        subplot(1,3,1)
        xlabel('F [Hz]');ylabel('S_{cross} [dBc/Hz]');
        grid on;axis tight;
        xlim([10,0.5e6]);
        title(['KID ' num2str(NOISE(kidn).KIDnumber) ' @Pread = ' num2str(NOISE(kidn).ReadPower) ' dBm']);
        legend(lstr);
        
        subplot(1,3,3)
        xlabel('T (K)');ylabel('\tau (msec)');grid on;ylim([0.01 round(1.5*tauinmax*1000)]);xlim([0.1 0.35])
        Figfile = [ChipInfo.path,filesep,FFTsubsubdir,filesep,...
            'KID',num2str(NOISE(kidn).KIDnumber,'%.0f'),'_crossPSD_Tdep'];
        MakeGoodFigure(18,6,12,Figfile); %save png and fig
    end
end
save([ChipInfo.path,filesep,FFTsubsubdir,filesep,matfile],'NOISE','-append');
rmpath([pwd,filesep,'..',filesep,'subroutines']);
end