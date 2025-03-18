% NEP_GR_PdepV1
% reads td binary data, calulates and plots the PSD's
% saves the data into the relevant struct
%

clear variables;close all;clc

%================================================================================
% Input
%================================================================================
ChipInfo_path = [ '..' filesep '..'];%without filesep at end %'/Volumes/kid/KIDonSun/experiments/Entropy ADR/LT169_chip5'; %root path where data is, one higher than the scripts %ChipInfo_path = [cd filesep '..'];
S21file = [filesep 'S21' filesep 'Temp' filesep 'ResponseS21.mat'];

eta_pb = 0.4;
tau_0 = 800e-9; %in sec, estimated or obtained from T dep data. Could be KID:KID coded (data is present) but getting good data is sketchy and this parameter shold be the same for all devices.
%================================================================================

%FFTsubsubdir=['FFT' filesep 'Power'];                   %no filesep, name of the noise data in /FFT (== TD data dir name)
FFTsubsubdir = [filesep 'FFT' filesep 'Power'];     %

matfile = 'Noise_P.mat';

addpath([pwd,filesep,'..',filesep,'subroutines']);

%load relevant data S21
load([ChipInfo_path,filesep,S21file],'KID');

%load relevant data noise files (level, lifetime etc)
load([ChipInfo_path,filesep,FFTsubsubdir,filesep,matfile],'NOISE','IndexP_sub_opt','KIDnumbers','IndexPopt','IndexPsort');
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
    
    %indexing in response array and getting some parameters
    KIDID = NOISE(IndexPopt(kidn)).KIDnumber;
    bla = [KID(:).KIDnumber];
    KIDindresp = (bla == KIDID);
    Fdesign = KID(KIDindresp).Fdesign;
    Delta = KID(KIDindresp).Delta;%not T dependent in the data!
    Tc = KID(KIDindresp).Tc;%not T dependent in the data!
    KIDID_resp = KID(KIDindresp).KIDnumber;%not T dependent in the data!
    Area = KID(KIDindresp).Area;%not T dependent in the data!
    Volume = KID(KIDindresp).Volume;%not T dependent in the data!
    Fres_resp = KID(KIDindresp).Fres(1);%Fres, not T dependent in the data!
    
    figure(NOISE(IndexP_sub_opt{kidn}(1)).KIDnumber)
    kleur = colormapJetJB(sum([NOISE(IndexP_sub_opt{kidn}).crosslevel]' ~= 0));
    ki = 1;%index for color plots
    notauforthiskid = 1;
    
    for p=1:length(IndexP_sub_opt{kidn})% over Power
        %check if we have a lifetome, othewsie we do not bother
        if NOISE(IndexP_sub_opt{kidn}(p)).crosslevel ~= 0 %if this is 0. no tau fit, so no NEP
            % create local variables with relevant data
            % Cross PSD results: NOte we also use the SR and Stheta
            % from the cross PSD program - the matlab is better than
            % labview
            notauforthiskid = 0;
            F	 = NOISE(IndexP_sub_opt{kidn}(p)).CrossPSD(:,1); %frequency cross PSD
            CPSD = 10*log10(abs(NOISE(IndexP_sub_opt{kidn}(p)).CrossPSD(:,2))); %cross PSD
            SR = NOISE(IndexP_sub_opt{kidn}(p)).CrossPSD(:,4);
            Stheta = NOISE(IndexP_sub_opt{kidn}(p)).CrossPSD(:,3);
            %cross fit results
            tau = NOISE(IndexP_sub_opt{kidn}(p)).tau;
            crosslevel = NOISE(IndexP_sub_opt{kidn}(p)).crosslevel;
            crosssetup = NOISE(IndexP_sub_opt{kidn}(p)).setupnoise;
            %KID ids
            KIDID = NOISE(IndexP_sub_opt{kidn}(p)).KIDnumber;
            Pint = NOISE(IndexP_sub_opt{kidn}(p)).InternalPower;
            Temperature = NOISE(IndexP_sub_opt{kidn}(p)).Temperature;
            Fres = NOISE(IndexP_sub_opt{kidn}(p)).Fres;
            % response
            KIDtau_T = 1e-9*KID(KIDindresp).Ql./(pi*KID(KIDindresp).Fres); %tau res vs T in sec
            
            %disp(['Check Fres: ' num2str([Fres Fres_resp])]);
            %disp(['Check ID: ' num2str([KIDID KIDID_resp])]);
            
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
            
            NOISE(IndexP_sub_opt{kidn}(p)).NEP(:,1) = NOISE(IndexP_sub_opt{kidn}(p)).CrossPSD(:,1); % == F
            NOISE(IndexP_sub_opt{kidn}(p)).NEP(:,2) = NEPR;
            NOISE(IndexP_sub_opt{kidn}(p)).NEP(:,3) = NEPtheta;
            NOISE(IndexP_sub_opt{kidn}(p)).NEP(:,4) = NEPcross;
            NOISE(IndexP_sub_opt{kidn}(p)).Nqp = Nqp;
            NOISE(IndexP_sub_opt{kidn}(p)).nqp = nqp;
            NOISE(IndexP_sub_opt{kidn}(p)).NEPGR = NEPGR;
            NOISE(IndexP_sub_opt{kidn}(p)).tau0 = tau0;
            NOISE(IndexP_sub_opt{kidn}(p)).dthetadN = dthetadN;
            NOISE(IndexP_sub_opt{kidn}(p)).dRdN = dRdN;
            NOISE(IndexP_sub_opt{kidn}(p)).Fdesign = Fdesign;
            
            % get NEP values smoothed above 10 Hz and below 1 kHz
            % NOISE(kidn).NEPRmin(nT,1) = F,
            % NOISE(kidn).NEPRmin(nT,2) = NEP value
            FR = find(NOISE(IndexP_sub_opt{kidn}(p)).NEP(:,1) < 1e3 & NOISE(IndexP_sub_opt{kidn}(p)).NEP(:,1) > 15);
            [NOISE(IndexP_sub_opt{kidn}(p)).NEPRmin(1,2),indi]      = min(smooth(NEPR(FR),5));
            NOISE(IndexP_sub_opt{kidn}(p)).NEPRmin(1,1)            = NOISE(IndexP_sub_opt{kidn}(p)).NEP(FR(indi),1);
            [NOISE(IndexP_sub_opt{kidn}(p)).NEPthetamin(1,2),indi]  = min(smooth(NEPtheta(FR),5));
            NOISE(IndexP_sub_opt{kidn}(p)).NEPthetamin(1,1)        = NOISE(IndexP_sub_opt{kidn}(p)).NEP(FR(indi),1);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            subplot(2,3,1) %Noise
            semilogx(NOISE(IndexP_sub_opt{kidn}(p)).CrossPSD(:,1),NOISE(IndexP_sub_opt{kidn}(p)).CrossPSD(:,4),...
                '-','color',kleur(ki,:),'LineWidth',1);%R
            hold on;grid on;xlim([1,1e5]);
            semilogx(NOISE(IndexP_sub_opt{kidn}(p)).CrossPSD(:,1),NOISE(IndexP_sub_opt{kidn}(p)).CrossPSD(:,3),...
                '-','color',kleur(ki,:),'LineWidth',2);%Theta
            
            subplot(2,3,2)% NEP
            loglog(NOISE(IndexP_sub_opt{kidn}(p)).NEP(:,1),         NOISE(IndexP_sub_opt{kidn}(p)).NEP(:,2),'-','Linewidth',1,'color',kleur(ki,:));%R
            grid on;hold on;
            loglog(NOISE(IndexP_sub_opt{kidn}(p)).NEP(:,1),         NOISE(IndexP_sub_opt{kidn}(p)).NEP(:,3),'-','Linewidth',2,'color',kleur(ki,:));%theta
            loglog(NOISE(IndexP_sub_opt{kidn}(p)).NEPRmin(1),       NOISE(IndexP_sub_opt{kidn}(p)).NEPRmin(2),'o','color',kleur(ki,:),'MarkerSize',8);
            loglog(NOISE(IndexP_sub_opt{kidn}(p)).NEPthetamin(1),   NOISE(IndexP_sub_opt{kidn}(p)).NEPthetamin(2),'o','color',kleur(ki,:),'markerfacecolor',kleur(ki,:),'MarkerSize',8);
            
            subplot(2,3,3) %Cross PSD
            semilogx(NOISE(IndexP_sub_opt{kidn}(p)).CrossPSD(:,1),-1*real(NOISE(IndexP_sub_opt{kidn}(p)).CrossPSD(:,2)),'-',...
                'linewidth',2,'color',kleur(ki,:));hold on;grid on
            semilogx(NOISE(IndexP_sub_opt{kidn}(p)).CrossPSD(:,1),NOISE(IndexP_sub_opt{kidn}(p)).CrossPSD(:,5),'-k');
            
            subplot(2,3,4)%lifetime vs P
            semilogy(NOISE(IndexP_sub_opt{kidn}(p)).InternalPower,NOISE(IndexP_sub_opt{kidn}(p)).tau*1e3,'o',...
                'color',kleur(ki,:),'MarkerFaceColor',kleur(ki,:));hold on;grid on;
            semilogy(NOISE(IndexP_sub_opt{kidn}(p)).InternalPower,NOISE(IndexP_sub_opt{kidn}(p)).taumin*1e3,'o','color',kleur(ki,:));
            semilogy(NOISE(IndexP_sub_opt{kidn}(p)).InternalPower,NOISE(IndexP_sub_opt{kidn}(p)).taumax*1e3,'o','color',kleur(ki,:));
            
            subplot(2,3,5)%nqp vs P
            semilogy(NOISE(IndexP_sub_opt{kidn}(p)).InternalPower,NOISE(IndexP_sub_opt{kidn}(p)).nqp,'o',...
                'color',kleur(ki,:),'MarkerFaceColor',kleur(ki,:));hold on;grid on;
            
            subplot(2,3,6)%NEP vs P
            semilogy(NOISE(IndexP_sub_opt{kidn}(p)).InternalPower,NOISE(IndexP_sub_opt{kidn}(p)).NEPRmin(2),'o',...
                'color',kleur(ki,:));hold on;grid on;
            semilogy(NOISE(IndexP_sub_opt{kidn}(p)).InternalPower,NOISE(IndexP_sub_opt{kidn}(p)).NEPthetamin(2),'o',...
                'color',kleur(ki,:),'MarkerFaceColor',kleur(ki,:));hold on;grid on;
            
            legc = legc+2;%legend counter
            ki = ki +1; %kleur index
        end %conditional oop if a tau fit is there
        
    end %power loop (p)
    
    
    if notauforthiskid == 0
        subplot(2,3,1)
        xlabel('F (Hz)');ylabel('S_x (dBc/Hz)')
        title(['KID ' num2str(NOISE(IndexP_sub_opt{kidn}(1)).KIDnumber) ]);
        
        subplot(2,3,2)
        legend('NEP_A','NEP_{\theta}','Location','SouthEast');
        xlim([10,1e3]);ylim([1e-20 5e-18]);
        xlabel('F (Hz)');ylabel('NEP (W/\surd Hz)')
        title(['V = ' num2str(Volume) '\mum^3, A = ' num2str(Area) '\mum^2']);
        
        subplot(2,3,3)
        xlabel('F [Hz]');ylabel('S_{cross} [dBc/Hz]');
        grid on;axis tight;
        xlim([10,1e5]);
        title(['Power Dependence @ ' num2str(round(1000*(NOISE(IndexP_sub_opt{kidn}(1)).Temperature))) ' mK'] );
        
        subplot(2,3,4)
        xlabel('P  (dBm)');ylabel('\tau (msec)');grid on;ylim([0.01 3]);
        title(['F_{res} = ' num2str(NOISE(IndexP_sub_opt{kidn}(1)).Fres) ' GHz'])
        
        subplot(2,3,5)
        xlabel('T  (K)');ylabel('n_{qp} (\mum^{-3})');grid on;axis tight;
        if ~isnan(max([NOISE(IndexP_sub_opt{kidn}).nqp]))
            ylim([1 max([NOISE(IndexP_sub_opt{kidn}).nqp])*1.1]);
        end
        title('n_{qp} from PSD_{cross}, responsivity and \tau')
        
        subplot(2,3,6)%NEP vs P
        legend('NEP_{R}','NEP_{\theta}','Location','NorthWest');
        xlabel('P  (dBm)');ylabel('NEP (W /\surd Hz)');grid on;ylim([1e-20 1e-17]);
        
    end
    
    Figfile = [NEPfolder,filesep,...
        'KID',num2str(NOISE(IndexP_sub_opt{kidn}(1)).KIDnumber,'%.0f'),'_NEP_Pdep1'];
    MakeGoodFigure(16,15,9,Figfile); %save png and fig
    disp(num2str(kidn));
    %for output csv, all params at Max power (pOpt according to noise)
    PoutData(kidn,1) = NOISE(IndexPopt(kidn)).Temperature;
    PoutData(kidn,2) = NOISE(IndexPopt(kidn)).KIDnumber;
    PoutData(kidn,3) = NOISE(IndexPopt(kidn)).InternalPower;
    PoutData(kidn,4) = NOISE(IndexPopt(kidn)).Fres;
    PoutData(kidn,5) = 0.1*round(NOISE(IndexPopt(kidn)).S21min*10);
    PoutData(kidn,6) = round(NOISE(IndexPopt(kidn)).Ql/1e3);
    PoutData(kidn,7) = round(NOISE(IndexPopt(kidn)).Qc/1e3);
    PoutData(kidn,8) = round(NOISE(IndexPopt(kidn)).Qi/1e3);
    PoutData(kidn,9) = Fdesign;
    PoutData(kidn,10) = Area;
    PoutData(kidn,11) = NOISE(IndexPopt(kidn)).MeanFreqNoise1000Hz;
    if notauforthiskid == 0
        %NB: the values below are only avilable for some powers, not
        %stored in an index. The raay is initialized with 0's. So
        %min/max works
        PoutData(kidn,12) = min([NOISE(IndexP_sub_opt{kidn}).NEPthetamin]);
        PoutData(kidn,13) = 1e3*max([NOISE(IndexP_sub_opt{kidn}).tau]);
    end
    %==================================================================================================================================
    %==================================================================================================================================
    
end %KID loop

save([NEPfolder,filesep,'NEP_P.mat'],'NOISE','IndexP_sub_opt','KIDnumbers','IndexPopt','IndexPsort');
PoutHeader = {'T (K)','KIDID','Pint_opt (dBm)','Fres (GHz)','S21 (dB)','Q','Qc','Qi','Fdesign','Area (um^2)','1kHzFnoise','NEPtheta','tau (msec)'}';
WriteSRONcsv([NEPfolder,filesep,'NEP_P.csv'],PoutData,PoutHeader,'%.6g')
rmpath([pwd,filesep,'..',filesep,'subroutines']);
