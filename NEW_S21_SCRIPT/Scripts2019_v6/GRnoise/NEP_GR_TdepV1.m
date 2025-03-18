% NEP_GR_TdepV1
% reads td binary data, calulates and plots the PSD's
% saves the data into the relevant struct
%

clear variables;close all;clc

%================================================================================
% Input
%================================================================================
ChipInfo_path = '../..';%'/Volumes/kid/KIDonSun/experiments/Entropy ADR/LT169_chip5'; %root path where data is, one higher than the scripts %ChipInfo_path = [cd filesep '..'];
S21file = [filesep 'S21' filesep 'Temp' filesep 'ResponseS21.mat'];
eta_pb = 0.4;
lowT = 0.2;     %Tbelow tau is ~ constant. 
FixTc = 1;  %1 = default. Tc not fitted in Kaplan fit. set 0 to fit it
%================================================================================

%Tdep
%FFTsubsubdir=['FFT' filesep 'Temp'];                   %no filesep, name of the noise data in /FFT (== TD data dir name)
FFTsubsubdir = [filesep 'FFT' filesep 'Temp'];     %

matfile = 'Noise_T.mat';

addpath([pwd,filesep,'..',filesep,'subroutines']);

%load relevant data S21
load([ChipInfo_path,filesep,S21file],'KID');

%load relevant data noise files (level, lifetime etc)
load([ChipInfo_path,filesep,FFTsubsubdir,filesep,matfile],'NOISE','KIDnumbers');
ChipInfo.path = ChipInfo_path;clear ChipInfo_path; %needed to keep struct ok wrt path names when data is moved in between analysis ruyns

%define NEP folder
NEPfolder = [ChipInfo.path, filesep, 'NEP'] ;
if ~isfolder(NEPfolder)
    mkdir(NEPfolder);
end
%define output array
PoutData = zeros(length(KIDnumbers),11);

for kidn=1:length(KIDnumbers) % LOOP OVER ALL UNIQUE KIDS,
    %construct filename
    legc = 1;%legend counter
    clear lstr
    
    %indexing in response array and getting some parameters for easy
    %working
    KIDID = NOISE(kidn).KIDnumber;
    bla = [KID(:).KIDnumber];
    KIDindresp = (bla == KIDID);
    Fdesign = KID(KIDindresp).Fdesign;
    Delta = KID(KIDindresp).Delta;%not T dependent in the data!
    Tc = KID(KIDindresp).Tc;%not T dependent in the data!
    KIDID_resp = KID(KIDindresp).KIDnumber;%not T dependent in the data!
    Area = KID(KIDindresp).Area;%not T dependent in the data!
    Volume = KID(KIDindresp).Volume;%not T dependent in the data!
    Fres_resp = KID(KIDindresp).Fres(1);%Fres, not T dependent in the data!
    
    %==================================================================================================================================
    %PLOT
    %==================================================================================================================================
    figure(NOISE(kidn).KIDnumber)
    kleur = colormapJetJB(length(NOISE(kidn).Temperature));
    ki = 1;%index for color plots
    notauforthiskid = 1;
    for nT=1:length(NOISE(kidn).Temperature) % over T
        %check if we have a lifetome, othewsie we do not bother
        if NOISE(kidn).crosslevel(nT) ~= 0 %if this is 0. no tau fit, so no NEP
            % create local variables with relevant data
            % Cross PSD results: NOte we also use the SR and Stheta
            % from the cross PSD program - the matlab is better than
            % labview
            notauforthiskid = 0;
            F	 = NOISE(kidn).CrossPSD{nT}(:,1); %frequency cross PSD
            CPSD = 10*log10(abs(NOISE(kidn).CrossPSD{nT}(:,2)));      %cross PSD
            SR = NOISE(kidn).CrossPSD{nT}(:,4);             %in dBc/Hz
            Stheta = NOISE(kidn).CrossPSD{nT}(:,3);         %in dBc/Hz
            %cross fit results
            tau = NOISE(kidn).tau(nT);
            crosslevel = NOISE(kidn).crosslevel(nT);
            crosssetup = NOISE(kidn).setupnoise(nT);
            %KID ids
            KIDID = NOISE(kidn).KIDnumber;
            Pint = NOISE(kidn).InternalPower(nT);
            Temperature = NOISE(kidn).Temperature(nT);
            Fres = NOISE(kidn).Fres(nT);
            % response
            KIDtau_T = 1e-9*KID(KIDindresp).Ql./(pi*KID(KIDindresp).Fres); %tau res vs T in sec
            
            %get reponse and other params @ correct T by interpolation
            if min(KID(KIDindresp).Temperature) <= Temperature && max(KID(KIDindresp).Temperature) > Temperature
                dthetadN =  interp1(KID(KIDindresp).Temperature, KID(KIDindresp).ResponsivityM1(:,1), Temperature);
                dRdN =      interp1(KID(KIDindresp).Temperature, KID(KIDindresp).ResponsivityM1(:,2), Temperature);
                taures =    interp1(KID(KIDindresp).Temperature, KIDtau_T, Temperature);
            elseif min(KID(KIDindresp).Temperature) > Temperature && max(KID(KIDindresp).Temperature) > Temperature
                dthetadN =  KID(KIDindresp).ResponsivityM1(1,1);
                dRdN =      KID(KIDindresp).ResponsivityM1(1,2);
                taures =    KIDtau_T(1);
            else
                error('Repsonse T range terrible, cannot consruct NEP')
            end
            
            % Get the NEP
            [NEPR,NEPtheta,NEPcross,Nqp,nqp,NEPGR,tau0] = getDARKNEP(dthetadN, dRdN, Stheta, SR, CPSD, F, eta_pb, tau, taures, crosslevel, Volume, Delta, Tc);
            
            %store NEP spectra
            NOISE(kidn).NEP{nT}(:,1) = NOISE(kidn).CrossPSD{nT}(:,1); % == F
            NOISE(kidn).NEP{nT}(:,2) = NEPR;       %Radius NEP
            NOISE(kidn).NEP{nT}(:,3) = NEPtheta;   %Phase NEP
            NOISE(kidn).NEP{nT}(:,4) = NEPcross;   %cross NEP
            NOISE(kidn).Nqp(nT)      = Nqp;
            NOISE(kidn).nqp(nT)      = nqp;
            NOISE(kidn).NEPGR(nT)    = NEPGR;
            NOISE(kidn).tau0_PSD(nT) = tau0;
            NOISE(kidn).dthetadN(nT) = dthetadN;
            NOISE(kidn).dRdN(nT)     = dRdN;
            
            % get NEP values smoothed above 10 Hz and below 1 kHz
            % NOISE(kidn).NEPRmin(nT,1) = F,
            % NOISE(kidn).NEPRmin(nT,2) = NEP value
            FR = find(NOISE(kidn).NEP{nT}(:,1) < 1e3 & NOISE(kidn).NEP{nT}(:,1) > 15);
            [NOISE(kidn).NEPRmin(nT,2),indi] = min(smooth(NEPR(FR),5));
            NOISE(kidn).NEPRmin(nT,1) = NOISE(kidn).NEP{nT}(FR(indi),1);
            [NOISE(kidn).NEPthetamin(nT,2),indi] = min(smooth(NEPtheta(FR),5));
            NOISE(kidn).NEPthetamin(nT,1) = NOISE(kidn).NEP{nT}(FR(indi),1);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            subplot(2,3,1) %Noise
            semilogx(NOISE(kidn).CrossPSD{nT}(:,1),NOISE(kidn).CrossPSD{nT}(:,4),...
                '-','color',kleur(ki,:),'LineWidth',1);%R
            hold on;grid on;xlim([1,1e5]);
            semilogx(NOISE(kidn).CrossPSD{nT}(:,1),NOISE(kidn).CrossPSD{nT}(:,3),...
                '-','color',kleur(ki,:),'LineWidth',2);%Theta
            
            subplot(2,3,2)% NEP
            loglog(NOISE(kidn).NEP{nT}(:,1),NOISE(kidn).NEP{nT}(:,2),'-','Linewidth',1,'color',kleur(ki,:));%R
            grid on;hold on;
            loglog(NOISE(kidn).NEP{nT}(:,1),NOISE(kidn).NEP{nT}(:,3),'-','Linewidth',2,'color',kleur(ki,:));%theta
            loglog(NOISE(kidn).NEPRmin(nT,1),NOISE(kidn).NEPRmin(nT,2),'o','color',kleur(ki,:),'MarkerSize',8);
            loglog(NOISE(kidn).NEPthetamin(nT,1),NOISE(kidn).NEPthetamin(nT,2),'o','color',kleur(ki,:),'markerfacecolor',kleur(ki,:),'MarkerSize',8);
            
            subplot(2,3,3) %Cross PSD
            semilogx(NOISE(kidn).CrossPSD{nT}(:,1),-1*real(NOISE(kidn).CrossPSD{nT}(:,2)),'-',...
                'linewidth',2,'color',kleur(ki,:));hold on;grid on
            semilogx(NOISE(kidn).CrossPSD{nT}(:,1),NOISE(kidn).CrossPSD{nT}(:,5),'-k');
            
            subplot(2,3,4)%lifetime vs P
            semilogy(NOISE(kidn).Temperature(nT),NOISE(kidn).tau(nT)*1e3,'o',...
                'color',kleur(ki,:),'MarkerFaceColor',kleur(ki,:));hold on;grid on;
            semilogy(NOISE(kidn).Temperature(nT),NOISE(kidn).taumin(nT)*1e3,'o','color',kleur(ki,:));
            semilogy(NOISE(kidn).Temperature(nT),NOISE(kidn).taumax(nT)*1e3,'o','color',kleur(ki,:));
            
            subplot(2,3,5)%nqp vs P Not made here
            
            subplot(2,3,6)%NEP vs P not made here
            
            legc = legc+2;%legend counter
            ki = ki +1; %kleur index
        else
            %store NEP values as NaN
            NOISE(kidn).Nqp(nT)      = NaN;
            NOISE(kidn).nqp(nT)      = NaN;
            NOISE(kidn).NEPGR(nT)    = NaN;
            NOISE(kidn).tau0_PSD(nT) = NaN;
            NOISE(kidn).dthetadN(nT) = NaN;
            NOISE(kidn).dRdN(nT)     = NaN;
        end %if statement if there is a fit
    end %end loop over all T
    
    %Now do T dependent stuff: Kaplan fit and Tau0 estimate
    
    if notauforthiskid == 0 % will nly be accessed for at least 1 fitted NEP
        %Kaplan fit
        if FixTc == 1
            [result,tau_fitted, nqpfullfit, nqp] = Fit_Kaplan2(NOISE(kidn).Temperature,NOISE(kidn).tau',lowT,Tc);%Tc
        elseif FixTc == 0
            [result,tau_fitted, nqpfullfit, nqp] = Fit_Kaplan2(NOISE(kidn).Temperature,NOISE(kidn).tau',lowT);%Tc
        else
            error('FixTc set wrong')
        end
        % I use nqp obtained from the Kaplan fit result tau_0 and measured lifetime.  
        %rest plot
        NOISE(kidn).tau_0_Kaplan = result.tau_0;
        NOISE(kidn).NEPGR_Kaplan = 2*Delta/eta_pb * (Volume*nqp./(2*NOISE(kidn).tau')).^0.5;     %PdV 2.41; using the CORRECT single particle tau = 2*tau measured%
        NOISE(kidn).taufitted = tau_fitted;
        NOISE(kidn).nqp_Kaplan = nqp;
        NOISE(kidn).Tc_Kaplan = result.Tc;
        
        subplot(2,3,1)
        legend('S_A','S_{theta}');
        xlabel('F (Hz)');ylabel('S_x (dBc/Hz)')
        title(['KID ' num2str(NOISE(kidn).KIDnumber) ]);
        
        subplot(2,3,2)
        legend('NEP_A','NEP_{\theta}','Location','SouthEast');
        xlim([1,1e4]);ylim([1e-20 3e-17]);
        xlabel('F (Hz)');ylabel('NEP (W/\surd Hz)')
        title(['V = ' num2str(Volume) '\mum^3 , Area = ' num2str(Area) '\mum^2']);
        
        subplot(2,3,3)
        xlabel('F [Hz]');ylabel('S_{cross} [dBc/Hz]');
        grid on;axis tight;
        xlim([10,1e5]);
        legend('Data','Fit');
        title('T dependence @ P_{opt}')
        
        subplot(2,3,4)
        semilogy(NOISE(kidn).Temperature,1e3*NOISE(kidn).taufitted,'k-');
        xlabel('T  (K)');ylabel('\tau (msec)');grid on;ylim([0.01 3]);xlim([0.1 0.3])
        title(['Kaplan \tau_0 = ' num2str(1e9*NOISE(kidn).tau_0_Kaplan,'%.0f') ' nsec.']);
        
        subplot(2,3,5)
        semilogy(NOISE(kidn).Temperature,NOISE(kidn).nqp,'kx','MarkerSize',9);hold on;grid on;
        semilogy(NOISE(kidn).Temperature,NOISE(kidn).nqp_Kaplan,'ks','MarkerSize',9);%nqpfullfit
        semilogy(NOISE(kidn).Temperature,nqpfullfit,'k-');
        xlabel('T  (K)');ylabel('n_{qp} (\mum^{-3})');grid on;axis tight;
        xlim([0.1 0.3]);ylim([1 max(NOISE(kidn).nqp)*1.1]);
        legend(['n_{qp} from PSD_{cross}, Tc=' num2str(Tc ,'%0.2f') 'K'],...
            ['n_{qp} from \tau_0 from Kaplan fit and \tau_m. Tc=' num2str(NOISE(kidn).Tc_Kaplan ,'%0.2f') ' K'],...
            ['n_{qp} from \tau_0 and \tau from Kaplan fit. Tc=' num2str(NOISE(kidn).Tc_Kaplan ,'%0.2f') ' K']...
            ,'location','SouthEast')
        
        subplot(2,3,6)
        semilogy(NOISE(kidn).Temperature, NOISE(kidn).NEPGR,'kx','MarkerSize',9);hold on
        semilogy(NOISE(kidn).Temperature, NOISE(kidn).NEPGR_Kaplan,'ks','MarkerSize',9);
        ki = 1;
        for nT=1:length(NOISE(kidn).Temperature) % over T
            %check if we have a lifetome, othewsie we do not bother
            if NOISE(kidn).crosslevel(nT) ~= 0 %if this is 0. no tau fit, so no NEP
                semilogy(NOISE(kidn).Temperature(nT),NOISE(kidn).NEPRmin(nT,2),'o',...
                    'color',kleur(ki,:));grid on;
                semilogy(NOISE(kidn).Temperature(nT),NOISE(kidn).NEPthetamin(nT,2),'o',...
                    'color',kleur(ki,:),'MarkerFaceColor',kleur(ki,:));hold on;grid on;
                ki = ki +1; %kleur index
            end
        end
        legend('NEP from from PSD_{cross}','NEP from FItted Kaplan \tau_0 +measured \tau',...
            'NEP_{R}','NEP_{\theta}','Location','NorthWest');
        xlabel('T  (K)');ylabel('NEP (W /\surd Hz)');grid on;ylim([1e-20 1e-16]);xlim([0.1 0.3])
        title(['Fres = ' num2str(NOISE(kidn).Fres(1)) ' GHz']);
        
       
        Figfile = [NEPfolder,filesep,...
            'KID',num2str(NOISE(kidn).KIDnumber,'%.0f'),'_NEP_Tdep1'];
        MakeGoodFigure(16,15,9,Figfile); %save png and fig
    end % if tau is found from cross PSD
end %KID loop

save([NEPfolder,filesep,'NEP_T.mat'],'NOISE','KIDnumbers');

rmpath([pwd,filesep,'..',filesep,'subroutines']);
