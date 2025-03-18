function NEP_vs_loadingv15C
% run after NEP_vs_loadingv#A !!!!!!!!!!!!
% updated wrt V14:
% error estimate optical NEP changed inot a weighted mean.
close all;
clear all;
clc
addpath([pwd,filesep,'subroutines']);                           %Enable subroutines by adding path in search path.

path=[cd '/../']; %root path where data is, one higher than the scripts
resppathy_C = [path '/2D_BB/2D_BB/'];
pbbtplot        = []; %for last figure, give array of PBB to plot (in sorted order), leave [] as default
KillPlots       = 0; %kills the plots after creation; required for more than 10 KIDs
% optical powers to remove for optical efficiency calc
startoptcut     = 5; %   start at BB power 1+ this value for opt eff calc (removes lowest Pbb)
endoptcut       = 1; %   end at end BB power range - this for opt eff calc (removes highest Pbb)
% settings to get the efficiency for each Pbb
usephase_eta    = 1; % =1 uses the opt_eff from the phase for final NE{ plot, for 0  uses amplitude
gettau          = 1; % gets the tau automatically, is wil output maxtau if fit fails
usephasefortau  = 0; % 0 for radius 1 for phase. If radius is present it works better (use 0)
maxtau          = 3e-3; % maximum realistic lifetime
plottaufit      = 0; % set to 1 only for debugging - creates the tau fit + data, but overplots in its frame
nepplotrange    = [1e-18 1e-15];% ylim for all NEP wrt source power figures
etasetup        = 1.4462e-08/((2.99e8/870e9)^2); % coupling efficiency due to apertures etc in the setup. 
                % only affects the last opt eff plot, NOT the outputted values.%
                % data from BP integrations @ 850 GHz (lambda^2 TP)
                % 1mmlens_3mm_17mm_850_GHz has TP = 1.4462e-08 @ 870 GHz,=> eta = 0.122%
                % 1mmlens_3mm_15mm_850_GHz has TP = 1.7301e-08; @ 870 GHz,=> eta = 0.146%%
                % NB: \eta = TP/lambda^2 for single mode

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(groot,'defaultLegendAutoUpdate','off');
%%%%%%%%%%%%%%%%%%%%%% Read in KIDparam.mat %%%%%%%%%%%%%%%%%%%%%%
load([resppathy_C 'KIDparam.mat'])
resppathy=resppathy_C;%catches issues with windows PC; resppathy is saved also in previous matlab.mat that we just loaded.
clear resppathy_C;
method.etasetup = etasetup;
%%%%%%%%%%%%%%%%%%%%%% Read in Popt.csv %%%%%%%%%%%%%%%%%%%%%%
% the optimum power values from Popt.csv are inported into the KIDparam.mat

[~,PoptData] = ReadSRONcsvV2([resppathy 'Popt.csv'],'',0);
rowi=1;
for nKID=1:nokids
    for tbb=1:length(KIDparam(nKID).Popt) %one Popt per BB temperature
        KIDparam(nKID).Popt(tbb)=PoptData(rowi,3); %Storing Popt
        %find Poptindex and store
        KIDparam(nKID).Poptindex(tbb)=find(KIDparam(nKID).Pread(:,tbb)==PoptData(rowi,3));
        
        rowi=rowi+1;
    end
end
if nokids>10
    KillPlots=1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% calculate optical efficiency, where NEPopt is used for every TBB

if noBBTS<=startoptcut-1
    error('optical efficiency not calculated, reduce startoptcut');
end

for nKID=1:nokids
    disp(['start KID no ' num2str(nKID) ' with ID: ' num2str(KIDparam(nKID).KIDid(1))]);
    %fill also Poptindexmatrix in loop
    Pindexmatrix=zeros(size(KIDparam(nKID).Pread));
    for PBB_n=1:noBBTS
        Pindexmatrix(KIDparam(nKID).Poptindex(PBB_n),PBB_n)=1;%ewhich power we need to dake from all readout powers (Popt as defined in Popt.cvs after script ...B)
    end
    Pindexmatrix=logical(Pindexmatrix); % logical to grab the optimal readout power for each Pbb
    
    %%%%%%%%%%%%Make range for opteff fit%%%%%%%%%%%%%%%%
    [KIDparam(nKID).Pbbnoise_Popt,Pbb_SI]=sort(KIDparam(nKID).Pbbnoise(Pindexmatrix));
    takenum=Pbb_SI(1+startoptcut:end-endoptcut);% takenum is the indices we want to use in all data that is vs Pbb (@ Popt) (UNSORTED)%
    % KIDparam(nKID).Pbbnoise_Popt is sorted PObb @ each Popt readout
    
    %%%%%%%%%%%%%%%% calculation of the optical efficiency in Amplitude @ Poptimal readout %%%%%%%%%%%%%%%%
    KIDparam(nKID).effradNEPdet=KIDparam(nKID).radiusNEPamplifier(Pindexmatrix);%LNA contribution
    % efficiency for each datapoint: Eff = (NEP_Poisson^2 + NEP_GR^2) /(NEP_exp_blip^2 - NEP_wave^2 )
    % with NEP_exp_blip^2 = NEP_measured^2 - NEP_detector^2 %
    KIDparam(nKID).radopteffvsP=(KIDparam(nKID).poisson(Pindexmatrix).^2 + KIDparam(nKID).g_r(Pindexmatrix).^2)./...
        (KIDparam(nKID).radiusNEPfref(Pindexmatrix).^2-KIDparam(nKID).effradNEPdet.^2-KIDparam(nKID).wave(Pindexmatrix).^2);
    % catching negatives
    KIDparam(nKID).radopteffvsP(KIDparam(nKID).radopteffvsP<0) = NaN;%catching neg values
    takenum(isnan(KIDparam(nKID).radopteffvsP(takenum)))=[];
    
    % error in efficiency per datapoint: 
    % Eff = (NEP_Poisson^2 + NEP_GR^2) /(NEP_exp^2 - NEP_wave^2 ). So, sigma(Eff) = sigma(NEP_exp) dEff/dNEP_exp %
    % rewrite Eff = (A+B)/(x^2-C); dEff/dx = 2x(A+B)/(x^2-C^2) and use
    % wolfram differentiator
    % sigma_Eff = sigma_NEP_exp * 2*NEP_exp*(NEP_Poisson^2 + NEP_GR^2)/(NEP_exp^2 - NEP_wave^2 )^2%
    KIDparam(nKID).stdradopteffvsP = ...
        KIDparam(nKID).stdradiusNEPfref(Pindexmatrix) .* 2 .* KIDparam(nKID).radiusNEPfref(Pindexmatrix) ...
        .* (KIDparam(nKID).poisson(Pindexmatrix).^2 + KIDparam(nKID).g_r(Pindexmatrix).^2)  ./  ...
        ((KIDparam(nKID).radiusNEPfref(Pindexmatrix).^2-KIDparam(nKID).effradNEPdet.^2-KIDparam(nKID).wave(Pindexmatrix).^2).^2);
    % old eqn (wrong) = KIDparam(nKID).radopteffvsP .* (2*KIDparam(nKID).stdradiusNEPfref(Pindexmatrix) ./ KIDparam(nKID).radiusNEPfref(Pindexmatrix));%error
 
    % Optical Efficiency: weighted mean of the optical efficiency radius readout
    KIDparam(nKID).optradeff = wmean(KIDparam(nKID).radopteffvsP(takenum),1./KIDparam(nKID).stdradopteffvsP(takenum)); %weighted mean of allowed values with 1/errors as weights
    % Error Opt. eff: weighted mean of the error of the optical efficiency radius readout, devided by sqrt(# points) %  
    mean_error = wmean(KIDparam(nKID).stdradopteffvsP(takenum),1./KIDparam(nKID).stdradopteffvsP(takenum));
    KIDparam(nKID).stdradeff = mean_error/sqrt(length(takenum));
    clear mean_error;
    % Recalculating the measured NEP from obtained optical efficiency and amplifier contribution EACH POINT INDIVIDUALLY %
    NEPradfit(:,nKID)=sqrt((KIDparam(nKID).poisson(Pindexmatrix).^2+KIDparam(nKID).g_r(Pindexmatrix).^2)./KIDparam(nKID).radopteffvsP...
        +KIDparam(nKID).wave(Pindexmatrix).^2+KIDparam(nKID).effradNEPdet.^2); %The fitted value of the NEP. i.e. recalculating the measured NEP including optical efficiency and amplifier contribution
    KIDparam(nKID).takenumrad = takenum;
    clear takenum
    
    %%%%%%%%%%%%%%%% calculation of the optical efficiency Phase (using subtraction of high post detect. F NEP as correction (=same as for Amplitude)%%%%%%%%%%%%%%%%
    % not commented, see amplitude for comments)
    takenum=Pbb_SI(1+startoptcut:end-endoptcut);                                    %indices sorted Pbb, At optimal readout power for each Pbb
    KIDparam(nKID).effPhaseNEPdet=KIDparam(nKID).phaseNEPamplifier(Pindexmatrix);   %LNA contribution
    %efficiency for each datapoint
    KIDparam(nKID).phaseopteffvsP=(KIDparam(nKID).poisson(Pindexmatrix).^2+KIDparam(nKID).g_r(Pindexmatrix).^2)./...
        (KIDparam(nKID).phaseNEPfref(Pindexmatrix).^2-KIDparam(nKID).effPhaseNEPdet.^2-KIDparam(nKID).wave(Pindexmatrix).^2); 
    % catching negatives
    KIDparam(nKID).phaseopteffvsP(KIDparam(nKID).phaseopteffvsP<0) = NaN;
    takenum(isnan(KIDparam(nKID).phaseopteffvsP(takenum)))=[];
    % error in efficiency per datapoint:
    KIDparam(nKID).stdphaseopteffvsP = ...
        KIDparam(nKID).stdphaseNEPfref(Pindexmatrix) .* 2 .* KIDparam(nKID).phaseNEPfref(Pindexmatrix) ...
        .* (KIDparam(nKID).poisson(Pindexmatrix).^2 + KIDparam(nKID).g_r(Pindexmatrix).^2)  ./  ...
        ((KIDparam(nKID).radiusNEPfref(Pindexmatrix).^2-KIDparam(nKID).effPhaseNEPdet.^2-KIDparam(nKID).wave(Pindexmatrix).^2).^2);
    
    % Optical Efficiency: weighted mean of the optical efficiency phase readout
    KIDparam(nKID).optphaseeff = wmean(KIDparam(nKID).phaseopteffvsP(takenum),1./KIDparam(nKID).stdphaseopteffvsP(takenum)); %weighted mean of allowed values with 1/errors as weights
   % Error Opt. eff: weighted mean of the error of the optical efficiency radius readout %  
    mean_error = wmean(KIDparam(nKID).stdphaseopteffvsP(takenum),1./KIDparam(nKID).stdphaseopteffvsP(takenum));
    KIDparam(nKID).stdphaseeff = mean_error/sqrt(length(takenum));
    clear mean_error;
    % Recalculating the measured NEP from obtained optical efficiency and amplifier contribution EACH POINT INDIVIDUALLY %
    NEPphasefit(:,nKID)=sqrt((KIDparam(nKID).poisson(Pindexmatrix).^2+KIDparam(nKID).g_r(Pindexmatrix).^2)./KIDparam(nKID).phaseopteffvsP...
        +KIDparam(nKID).wave(Pindexmatrix).^2+KIDparam(nKID).effPhaseNEPdet.^2);
    KIDparam(nKID).takenumphase = takenum;
    clear takenum
    
    
    %%% Now we calculate the NEP(Pabs) using the effciency from specified method
    if usephase_eta == 1
        method.eta_c = KIDparam(nKID).optphaseeff;
    elseif usephase_eta == 0
        method.eta_c = KIDparam(nKID).optradeff;
    else
        error('usephase_eta not 1 or 0')
    end
    
    %%%%%%%%%%%%%%%%%%%%%% call new BB script to find Pabs and the NEP's %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    BBcal.T=0; %to make sure we call the BB interation again, and for each MKID since eta changes.
    [KIDparam(nKID).Pbbnoise_abs,poo,BBcal]=blackbody_int(KIDparam(nKID).Tbbnoise(Pindexmatrix),BBcal,method);
    if nKID ~= 1
       close 123451;
    end
    set(gcf,'Color','White')
    Figfile=[resppathy 'KID_' num2str(respkids(nKID)) '_' num2str(KIDparam(nKID).Tchip(1,1),'%.2g') 'BBfigure.fig'];
    saveas(gcf,Figfile,'fig')
    KIDparam(nKID).wave_abs=poo.wave;
    KIDparam(nKID).g_r_abs=poo.g_r;
    KIDparam(nKID).poisson_abs=poo.poisson;
    KIDparam(nKID).totphoton_abs=poo.totphoton;
    KIDparam(nKID).phaseNEPfref_abs=method.eta_c*KIDparam(nKID).phaseNEPfref(Pindexmatrix);%NEP_exp = \eta*NEP_exp(eta=1) (only response is affected now!)%
    KIDparam(nKID).stdphaseNEPfref_abs=method.eta_c*KIDparam(nKID).stdphaseNEPfref(Pindexmatrix);%NEP_exp = \eta*NEP_exp(eta=1) (only response is affected now!)%
    
    KIDparam(nKID).radiusNEPfref_abs=method.eta_c*KIDparam(nKID).radiusNEPfref(Pindexmatrix);%see word file
    KIDparam(nKID).stdradiusNEPfref_abs=method.eta_c*KIDparam(nKID).stdradiusNEPfref(Pindexmatrix);%NEP_exp = \eta*NEP_exp(eta=1) (only response is affected now!)%
    
    clear meh poo
    %now to be able to plot full ranges
    [KIDparam(nKID).Pbbnoise_absFR,poo,BBcal]=blackbody_int([2:0.5:30 35:5:300],BBcal,method);
    KIDparam(nKID).wave_absFR=poo.wave;
    KIDparam(nKID).g_r_absFR=poo.g_r;
    KIDparam(nKID).poisson_absFR=poo.poisson;
    KIDparam(nKID).totphoton_absFR=poo.totphoton;
    KIDparam(nKID).phaseNEPfref_absFR=method.eta_c*KIDparam(nKID).phaseNEPfref(Pindexmatrix);%NEP_exp = \eta*NEP_exp(eta=1) (only response is affected now!)%
    KIDparam(nKID).radiusNEPfref_absFR=method.eta_c*KIDparam(nKID).radiusNEPfref(Pindexmatrix);%see word file
    clear meh poo
    %%%%%%%%%%%%%%%%%%%
    
    %get the tau from noise rolloff - a value is allwasy given to tau, out
    %of rang eif fit not enables
    
    if gettau == 1
        for PBB_n=1:noBBTS
            if usephasefortau==1;
                noisedata = KIDparam(nKID).phasenoise{KIDparam(nKID).Poptindex(PBB_n),PBB_n};
            else
                noisedata = KIDparam(nKID).ampnoise{KIDparam(nKID).Poptindex(PBB_n),PBB_n};
            end
            fdata = KIDparam(nKID).f_noise{KIDparam(nKID).Poptindex(PBB_n),PBB_n};
            KIDparam(nKID).tau{KIDparam(nKID).Poptindex(PBB_n),PBB_n} = ...
                gettau_noise(fdata,noisedata,indfref,blawindow,maxtau,plottaufit);
            
        end
        
    else
        for PBB_n=1:noBBTS
            KIDparam(nKID).tau{KIDparam(nKID).Poptindex(PBB_n),PBB_n} = 1e-6;%set out of range value
            
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FIRST FIGURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% vs  Psource %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    figure(respkids(nKID));
    set(gcf,'Position', [100, 0, 1200, 800]);
    set(gcf,'Color','White')
    colorsPBB = colormap(jet(noBBTS));
    for PBB_nunsorted=1:noBBTS
        PBB_n = Pbb_SI(PBB_nunsorted);%is increasing in power
        subplot(2,3,1);
        plot(KIDparam(nKID).S21Real{KIDparam(nKID).Poptindex(PBB_n),PBB_n},KIDparam(nKID).S21Imag{KIDparam(nKID).Poptindex(PBB_n),PBB_n},'-k','LineWidth',1);hold on;%volt S21
        xlabel('Re');ylabel('Im ');
        axis tight;
        plot(KIDparam(nKID).Reresp{KIDparam(nKID).Poptindex(PBB_n),PBB_n},KIDparam(nKID).Imresp{KIDparam(nKID).Poptindex(PBB_n),PBB_n},...
            '.','MarkerSize',8,'color',colorsPBB(PBB_nunsorted,:))
        
        legend('Fscan','TD data')
        if PBB_n==1 % low P Q factor
            title(['KID ' num2str(respkids(nKID),'%.0f') ',' 'Q= ' num2str(KIDparam(nKID).Q(PBB_n)/1e5,'%.3g') 'e5' ', Qi= '  num2str(KIDparam(nKID).Qi(PBB_n)/1e5,'%.3g') 'e5, at Popt']);
        end
        
        subplot(2,3,2);
        plot(KIDparam(nKID).Pbb{KIDparam(nKID).Poptindex(PBB_n),PBB_n}*1E15,KIDparam(nKID).Phaseresp{KIDparam(nKID).Poptindex(PBB_n),PBB_n},...
            'o','MarkerSize',4,'color',colorsPBB(PBB_nunsorted,:),'MarkerFaceColor',colorsPBB(PBB_nunsorted,:) );hold on;
        plot(KIDparam(nKID).Pbb{KIDparam(nKID).Poptindex(PBB_n),PBB_n}*1E15,KIDparam(nKID).Radiusresp{KIDparam(nKID).Poptindex(PBB_n),PBB_n}-1,...
            'o','MarkerSize',4,'color',colorsPBB(PBB_nunsorted,:))
        
        legend('Phase','Radius-1')
        ylabel('Reponse')
        xlabel('P_{s} (fW)')
        plot(KIDparam(nKID).Pbb{KIDparam(nKID).Poptindex(PBB_n),PBB_n}*1E15,KIDparam(nKID).respfitRup{KIDparam(nKID).Poptindex(PBB_n),PBB_n},'k');% fits
        plot(KIDparam(nKID).Pbb{KIDparam(nKID).Poptindex(PBB_n),PBB_n}*1E15,KIDparam(nKID).respfitthup{KIDparam(nKID).Poptindex(PBB_n),PBB_n},'k');
        if KIDparam(nKID).plotthrespdown == 1
            plot(KIDparam(nKID).Pbb{KIDparam(nKID).Poptindex(PBB_n),PBB_n}*1E15,KIDparam(nKID).respfitRdown{KIDparam(nKID).Poptindex(PBB_n),PBB_n},'--k');
        end
        if KIDparam(nKID).plotRrespdown == 1
            plot(KIDparam(nKID).Pbb{KIDparam(nKID).Poptindex(PBB_n),PBB_n}*1E15,KIDparam(nKID).respfitthdown{KIDparam(nKID).Poptindex(PBB_n),PBB_n},'--k');
        end
        title('Old figure vs source power')
        
        subplot(2,3,4);
        semilogx(KIDparam(nKID).f_noise{KIDparam(nKID).Poptindex(PBB_n),PBB_n},KIDparam(nKID).phasenoise{KIDparam(nKID).Poptindex(PBB_n),PBB_n},'-','color',colorsPBB(PBB_nunsorted,:)) % ,'r');xlim([0.5 20000]);grid on;hold on;%ylim([-100 -50]);
        xlim([0.5 200000]);grid on;hold on;
        %semilogx(KIDparam(nKID).f_noise{KIDparam(nKID).Poptindex(PBB_n),PBB_n},KIDparam(nKID).ampnoise{KIDparam(nKID).Poptindex(PBB_n),PBB_n},'--','color',colorsPBB(PBB_n,:))%,'b','LineWidth',1);%xlim([1 5000]);grid on;ylim([10^-19 10^-16]);hold on;
        ylabel('Noise (dBk/Hz)');xlabel('F [Hz]')
        %title(['Noise at ' num2str(qrt(1,1)) ' K and @Popt= ' num2str(Popt(nKID))])
        
        subplot(2,3,5);%method.eta_c*KIDparam(nKID).phaseNEPfref(Pindexmatrix);
        loglog(KIDparam(nKID).f_noise{KIDparam(nKID).Poptindex(PBB_n),PBB_n},KIDparam(nKID).PhaseoptNEP{KIDparam(nKID).Poptindex(PBB_n),PBB_n},'-','color',colorsPBB(PBB_nunsorted,:))
        xlim([0.5 200000]);grid on;hold on;%ylim([10^-19 10^-16]);
        loglog(KIDparam(nKID).f_noise{KIDparam(nKID).Poptindex(PBB_n),PBB_n},KIDparam(nKID).RadiusoptNEP{KIDparam(nKID).Poptindex(PBB_n),PBB_n},'--','color',colorsPBB(PBB_nunsorted,:)) % 'b','LineWidth',1);%xlim([1 5000]);grid on;ylim([10^-19 10^-16]);hold on;
        ylabel('NEP(P_s) [W\surd Hz]');xlabel('F [Hz]')
        ylim(nepplotrange)
        legendst{PBB_nunsorted}=['P_s = ' num2str(KIDparam(nKID).Pbbnoise((PBB_n))*1e15,'%.3g') ' fW'];
    end
    subplot(2,3,4);
    legend(legendst);
    for PBB_nunsorted=1:noBBTS
        PBB_n = Pbb_SI(PBB_nunsorted);%is increasing in power
        semilogx(KIDparam(nKID).f_noise{KIDparam(nKID).Poptindex(PBB_n),PBB_n},KIDparam(nKID).ampnoise{KIDparam(nKID).Poptindex(PBB_n),PBB_n},'--','color',colorsPBB(PBB_n,:))%,'b','LineWidth',1);%xlim([1 5000]);grid on;ylim([10^-19 10^-16]);hold on;
    end
    
    
    subplot(2,3,3)
    semilogx(KIDparam(nKID).Pbbnoise(Pindexmatrix)*1e15,KIDparam(nKID).dthdP(Pindexmatrix),'-.r');hold on;
    semilogx(KIDparam(nKID).Pbbnoise(Pindexmatrix)*1e15,-KIDparam(nKID).dRdP(Pindexmatrix),'-+k')
    errorbar(KIDparam(nKID).Pbbnoise(Pindexmatrix)*1e15,KIDparam(nKID).dthdP(Pindexmatrix),KIDparam(nKID).dthdPstd(Pindexmatrix),'-.r');
    errorbar(KIDparam(nKID).Pbbnoise(Pindexmatrix)*1e15,-KIDparam(nKID).dRdP(Pindexmatrix),KIDparam(nKID).dRdPstd(Pindexmatrix),'-+k');
    ylabel('Responsivity (/W)')
    xlabel('P_{s} (fW)')
    legend('phase','radius')
    axis tight
    
    subplot(2,3,6)
    [Psorted{nKID},PsI{nKID}]=sort(KIDparam(nKID).Pbbnoise(Pindexmatrix)*1e15);
    loglog(KIDparam(nKID).Pbbnoise(Pindexmatrix)*1e15,KIDparam(nKID).phaseNEPfref(Pindexmatrix),'or','MarkerSize',6,'MarkerFaceColor','r');hold on;
    loglog(KIDparam(nKID).Pbbnoise(Pindexmatrix)*1e15,KIDparam(nKID).radiusNEPfref(Pindexmatrix),'ok','MarkerSize',6,'MarkerFaceColor','k');
    photonlimit{nKID}=KIDparam(nKID).totphoton(Pindexmatrix);
    loglog(Psorted{nKID}, photonlimit{nKID}(PsI{nKID}),'-b')
    loglog(KIDparam(nKID).Pbbnoise(Pindexmatrix)*1e15,KIDparam(nKID).effPhaseNEPdet,'sr')
    loglog(KIDparam(nKID).Pbbnoise(Pindexmatrix)*1e15,KIDparam(nKID).effradNEPdet,'sk')
    loglog(Psorted{nKID},NEPphasefit(PsI{nKID},nKID),'--r')
    loglog(Psorted{nKID},NEPradfit(PsI{nKID},nKID),'--k')
    %loglog(KIDparam(nKID).Pbbnoise(Pindexmatrix)*1e15,sqrt((KIDparam(nKID).poisson(Pindexmatrix).^2+KIDparam(nKID).g_r(Pindexmatrix).^2)./KIDparam(nKID).optradeff+KIDparam(nKID).wave(Pindexmatrix).^2),'--b')
    %loglog(KIDparam(nKID).Pbbnoise(Pindexmatrix)*1e15,sqrt((KIDparam(nKID).poisson(Pindexmatrix).^2+KIDparam(nKID).g_r(Pindexmatrix).^2)./KIDparam(nKID).optphaseeff+KIDparam(nKID).wave(Pindexmatrix).^2),'-.b')
    ylabel('NEP(P_s) @ ref. f')
    xlabel('P_s (fW)')
    legend('phase','radius','photon','Setup NEP phase','Setup NEP amplitute','fit phase','fit radius')%,'photon/\surd \eta','photon/\surd \eta')
    title(['Optical eff radius: ' num2str(KIDparam(nKID).optradeff,'%.3g') ' phase: ' num2str(KIDparam(nKID).optphaseeff,'%.3g') ' @' num2str(fref) 'Hz'])
    ylim(nepplotrange)
    grid on;
    Figfile=[resppathy 'KID_' num2str(respkids(nKID)) '_' num2str(KIDparam(nKID).Tchip(1,1),'%.2g') 'K_vsP_s.fig'];
    
    saveas(gcf,Figfile,'fig')
    if KillPlots==1
        close(gcf);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FIGURE 1000 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PHASE and AMPLITUDE vs Pabs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clear legendst
    figure(respkids(nKID)+1000) %phase plot
    set(gcf,'Position', [100, 0, 1200, 800]);
    set(gcf,'Color','White')
    nepplotrangePabs = nepplotrange * 10^floor(log10(KIDparam(nKID).optphaseeff));
    
    if ~isempty(pbbtplot)
        noBBtoplot = length(pbbtplot);
    else
        noBBtoplot = noBBTS;
    end
    for PBB_nunsorted=1:noBBtoplot
        if isempty(pbbtplot)
            PBB_n = Pbb_SI(PBB_nunsorted);%is increasing in power
        else
            PBB_ni = Pbb_SI(PBB_nunsorted);%is increasing in power
            PBB_n = PBB_n(pbbtplot);clear PBB_ni
        end
        
        %%%% NEP %%%%
        subplot(2,2,2)
        %phase
        TrueNEP = method.eta_c*KIDparam(nKID).PhaseoptNEP{KIDparam(nKID).Poptindex(PBB_n),PBB_n} .* ...
            (1 + (2*pi*KIDparam(nKID).f_noise{KIDparam(nKID).Poptindex(PBB_n),PBB_n} * KIDparam(nKID).tau{KIDparam(nKID).Poptindex(PBB_n),PBB_n}).^2).^0.5;
        loglog(KIDparam(nKID).f_noise{KIDparam(nKID).Poptindex(PBB_n),PBB_n},TrueNEP,'color',colorsPBB(PBB_nunsorted,:))%,'b','LineWidth',1);%xlim([1 5000]);grid on;ylim([10^-19 10^-16]);hold on;
        hold on;
        clear TrueNEP;
        %%%% Amplitude NEP, corrected for roll off
        TrueNEP = method.eta_c*KIDparam(nKID).RadiusoptNEP{KIDparam(nKID).Poptindex(PBB_n),PBB_n} .* ...
            (1 + (2*pi*KIDparam(nKID).f_noise{KIDparam(nKID).Poptindex(PBB_n),PBB_n} * KIDparam(nKID).tau{KIDparam(nKID).Poptindex(PBB_n),PBB_n}).^2).^0.5;
        loglog(KIDparam(nKID).f_noise{KIDparam(nKID).Poptindex(PBB_n),PBB_n},TrueNEP,'--','color',colorsPBB(PBB_nunsorted,:))%,'b','LineWidth',1);%xlim([1 5000]);grid on;ylim([10^-19 10^-16]);hold on;
        clear TrueNEP;
        
        %phase NEP points
        loglog(fref,KIDparam(nKID).phaseNEPfref_abs(PBB_n),'o','color',colorsPBB(PBB_nunsorted,:),'MarkerSize',7,'MarkerFaceColor',colorsPBB(PBB_nunsorted,:)); %plot fref points on the line
        %Amp NEP
        loglog(fref,KIDparam(nKID).radiusNEPfref_absFR(PBB_n),'s','color',colorsPBB(PBB_nunsorted,:),'MarkerSize',7,'MarkerFaceColor',colorsPBB(PBB_nunsorted,:)); %plot fref points on the line
        
        grid on; grid minor;
        ylabel('NEP(P_{abs}) (W/\surd Hz)');
        xlabel('F (Hz)');xlim([0.5 10000]);ylim(nepplotrangePabs)
        legend('Phase Readout','Amplitude Readout')
        
        %Noise spectra
        subplot(2,2,1) %phase
        semilogx(KIDparam(nKID).f_noise{KIDparam(nKID).Poptindex(PBB_n),PBB_n},KIDparam(nKID).phasenoise{KIDparam(nKID).Poptindex(PBB_n),PBB_n},'color',colorsPBB(PBB_nunsorted,:))%,'b','LineWidth',1);%xlim([1 5000]);grid on;ylim([10^-19 10^-16]);hold on;
        hold on;
        if KIDparam(nKID).Pbbnoise_abs(PsI{nKID}(PBB_n))<1e-15
            legendst{PBB_nunsorted}=['P = ' num2str(KIDparam(nKID).Pbbnoise_abs(PsI{nKID}(PBB_n))*1e18,'%.3g') ' aW'];
        elseif KIDparam(nKID).Pbbnoise_abs(PsI{nKID}(PBB_n))<1e-12
            legendst{PBB_nunsorted}=['P = ' num2str(KIDparam(nKID).Pbbnoise_abs(PsI{nKID}(PBB_n))*1e15,'%.3g') ' fW'];
        elseif KIDparam(nKID).Pbbnoise_abs(PsI{nKID}(PBB_n))<1e-9
            legendst{PBB_nunsorted}=['P = ' num2str(KIDparam(nKID).Pbbnoise_abs(PsI{nKID}(PBB_n))*1e12,'%.3g') ' pW'];
        else
            legendst{PBB_nunsorted}=['P = ' num2str(KIDparam(nKID).Pbbnoise_abs(PsI{nKID}(PBB_n))*1e9,'%.3g') ' nW'];
        end
    end
    legend(legendst);
    
    %amp noise, to allow correct legend
    for PBB_nunsorted=1:noBBtoplot
        if isempty(pbbtplot)
            PBB_n = Pbb_SI(PBB_nunsorted);%is increasing in power
        else
            PBB_ni = Pbb_SI(PBB_nunsorted);%is increasing in power
            PBB_n = PBB_n(pbbtplot);clear PBB_ni
        end
        subplot(2,2,1)
        semilogx(KIDparam(nKID).f_noise{KIDparam(nKID).Poptindex(PBB_n),PBB_n},KIDparam(nKID).ampnoise{KIDparam(nKID).Poptindex(PBB_n),PBB_n},'--','color',colorsPBB(PBB_nunsorted,:))%,'b','LineWidth',1);%xlim([1 5000]);grid on;ylim([10^-19 10^-16]);hold on;
        xlim([0.5 350000]);grid on;hold on;grid minor;
        ylabel('S_{R} (dBc/Hz)');xlabel('F (Hz)')
    end
    
    %NEP vs Pabs
    PPlot = KIDparam(nKID).Pbbnoise_abs*1e15; %unsorted Pbb
    subplot(2,2,3)
    loglog(PPlot(PsI{nKID}),KIDparam(nKID).phaseNEPfref_abs(PsI{nKID}),'ob','MarkerSize',7,'MarkerFaceColor','b'); hold on;%Measured NEP
    loglog(PPlot(PsI{nKID}),KIDparam(nKID).radiusNEPfref_abs(PsI{nKID}),'sr','MarkerSize',7,'MarkerFaceColor','r'); hold on;%Measured NEP
    legend('NEP_{\theta}(P_{abs})','NEP_{Amp}(P_{abs})','NEP_{BLIP}(P_{abs})','Setup contribution','Location','Best')
    loglog(PPlot(PsI{nKID}),KIDparam(nKID).totphoton_abs(PsI{nKID}),'o-k','MarkerSize',4,'MarkerFaceColor','k')
    ylabel('NEP at f_{ref} (W/\surd Hz) ')
    xlabel('P_{abs} (fW)')
    if usephase_eta==1
        title(['NEP vs Pabs abotained using \eta_{phase} = ' num2str(KIDparam(nKID).optphaseeff,'%.3g')]);
    else
        title(['NEP vs Pabs abotained using \eta_{R} = ' num2str(KIDparam(nKID).optradeff,'%.3g')]);
    end
    grid on;grid minor;
    ylim(nepplotrangePabs)
    
    %Efficiency KIDparam(nKID).takenumphase
    subplot(2,2,4)
    semilogx(PPlot(KIDparam(nKID).takenumphase),...
        KIDparam(nKID).phaseopteffvsP(KIDparam(nKID).takenumphase)/method.etasetup,'ob','MarkerFaceColor','b','MarkerSize',10);hold on
    semilogx(PPlot(KIDparam(nKID).takenumrad),...
        KIDparam(nKID).radopteffvsP(KIDparam(nKID).takenumrad)/method.etasetup,'rs','MarkerFaceColor','r','MarkerSize',10);
    legend(['\eta_{theta, setup corrected} = ' num2str(KIDparam(nKID).optphaseeff/method.etasetup,'%.2g') ' \pm ' num2str(KIDparam(nKID).stdphaseeff/method.etasetup,'%.2f')],...
        ['\eta_{R, setup corrected} = ' num2str(KIDparam(nKID).optradeff/method.etasetup,'%.2g') ' \pm ' num2str(KIDparam(nKID).stdradeff/method.etasetup,'%.2f')],...
        'AutoUpdate','off');
    %all points
    semilogx(PPlot(PsI{nKID}),KIDparam(nKID).phaseopteffvsP(PsI{nKID})/method.etasetup,'ob','MarkerSize',10);hold on;
    semilogx(PPlot(PsI{nKID}),KIDparam(nKID).radopteffvsP(PsI{nKID})/method.etasetup,'rs','MarkerSize',10);
    %error bars
    errorbar(PPlot(PsI{nKID}),KIDparam(nKID).phaseopteffvsP(PsI{nKID})/method.etasetup, KIDparam(nKID).stdphaseopteffvsP(PsI{nKID})/method.etasetup,'ob','MarkerSize',10);
    errorbar(PPlot(PsI{nKID}),KIDparam(nKID).radopteffvsP(PsI{nKID})  /method.etasetup, KIDparam(nKID).stdradopteffvsP(PsI{nKID})  /method.etasetup,'or','MarkerSize',10);
    %lines
    
    semilogx([min(PPlot) max(PPlot)],[1 1]*KIDparam(nKID).optphaseeff/method.etasetup,'b-');
    semilogx([min(PPlot) max(PPlot)],[1 1]*(KIDparam(nKID).optphaseeff+KIDparam(nKID).stdphaseeff)/method.etasetup,'b--');
    semilogx([min(PPlot) max(PPlot)],[1 1]*(KIDparam(nKID).optphaseeff-KIDparam(nKID).stdphaseeff)/method.etasetup,'b--');
    semilogx([min(PPlot) max(PPlot)],[1 1]*KIDparam(nKID).optradeff/method.etasetup,'r-');
    semilogx([min(PPlot) max(PPlot)],[1 1]*(KIDparam(nKID).optradeff+KIDparam(nKID).stdradeff)/method.etasetup,'r--');
    semilogx([min(PPlot) max(PPlot)],[1 1]*(KIDparam(nKID).optradeff-KIDparam(nKID).stdradeff)/method.etasetup,'r--');
    
    xlabel('P_{abs}');
    ylabel('\eta_{opt}');
    title(['Setup corrected efficiency @ ' num2str(fref) ' \pm ' num2str(frefrange) ' Hz'])
    grid on;grid minor;
    ylim([0 1]);
    
    
    Figfile=[resppathy 'KID_' num2str(respkids(nKID)) '_vsPabs.fig'];
    saveas(gcf,Figfile,'fig')
    if KillPlots==1
        close(gcf);
    end
    %%%%%%%%end combined figure%%%%%%%%%%
    
    disp(['Done KID no ' num2str(nKID) ' with ID: ' num2str(KIDparam(nKID).KIDid(1))])
end



clear pbbtplot etasetup legendstr legendstr2 ftype endoptcut Figfile gettau ind1pW KillPlots legendst path Pbb_SI  ...
    plottaufit Poptfile rowi startoptcut usephase_eta usephasefortau maxtau plottaufit

save([resppathy 'KIDparam.mat']) %

fprintf('NEP_vs_loadingv12C is finished \n')

rmpath([pwd,filesep,'subroutines']);                           %Disable subroutines by removing path in search path.

%%%%%%%%%%%%%OUTPUT: this describes only what is added to KIDparam in C,
%%%%%%%%%%%%%the rest is described at the end of function A
% KIDparam(nKID).effradNEPdet                       %Detector NEP from tail of amplitude NEP, 1D double
% KIDparam(nKID).radopteffvsP                       %optical efficiency for different PBB in amplitude, 1D double
% KIDparam(nKID).optradeff                          %optical efficiency from range in PBB in amplitude, 1 number
% KIDparam(nKID).erroptradeff                       %std in optical efficiency, 1 number
% KIDparam(nKID).optphaseeff                        %phase optical efficiency 1 number
% KIDparam(nKID).effphaseNEPdet                     %Detector NEP from tail of phase NEP, 1D double
% KIDparam(nKID).phaseopteffvsP                     %optical efficiency for different PBB in phase, 1D double
% KIDparam(nKID).RadiusoptNEProlloff{mpread,PBB_n}  %(optional) amplitude NEP with the roll-off term, cell
% KIDparam(nKID).PhaseoptNEProlloff{mpread,PBB_n}   %(optional) phase NEP with the roll-off term, cell

end %of main function