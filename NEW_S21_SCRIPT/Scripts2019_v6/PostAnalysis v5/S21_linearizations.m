close all
clear all
ChipInfo.path = ...
    [cd '/../../LyotStopLT072_1']; %root path where data is, one higher than the scripts
S21subdir='/S21/1D';
S21path=[ChipInfo.path S21subdir filesep] %Path containing raw S21 data

load([S21path,'ResponseS21.mat'])
for n=1:length(KIDnumbers)
    figure(n)
    
    Nbbtemperatures=length(KID(n).Temperature);
    Tcolors=colormap(jet(Nbbtemperatures));
    %---------------- KID CIRCLE ----------------%
    subplot(2,2,1)
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
        %Plot the F0(T0) points in the complex plane
        plot(KID(n).ReImF0(p,1),KID(n).ReImF0(p,2),'kd','MarkerSize',5,'MarkerFaceColor','r'); %Location Fres(Tbase)
    end
    %Make a nice layout
    grid on
    xlabel('Re')
    ylabel('Im')
    title(['KID ' num2str(KID(n).KIDnumber) ' @ Pread= ' num2str(KID(n).ReadPower(1))])
    axis([-1 1 -1 1]*0.5)
    hold off
    %% ---------------- making laura plot ----------------%
    close all
    p=1;
    ReS21=KID(n).S21data{p}(:,4);
    ImS21=KID(n).S21data{p}(:,5);
    FreqS21=KID(n).S21data{p}(:,1);
    angleS21=(atan2(ImS21,ReS21));
    
    for p=1:Nbbtemperatures
    ReBB(p)=KID(n).ReImF0(p,1)
    ImBB(p)=KID(n).ReImF0(p,2)
    Fresbb(p)=KID(n).Fres(1);
    end
    angleBB=atan2(ImBB,ReBB);
    %making 0 angle as the reference
    angleS21=flip(unwrap(flip(angleS21-angleBB(1))));%flip prevents wrong jumps of unwrap
    angleBB=(angleBB-angleBB(1));
    
    Fres_Laura = 2*KID(n).Fres(1)-interp1(angleS21,FreqS21,angleBB,'spline');%F_laura=Fres0-dF and dF=Fi(interp)-Fres0
    
    figure(1000)
    subplot(2,1,1)
    R=1:length(ReS21);
    plot(ReS21(R),ImS21(R),'bx');hold on
    plot(ReBB,ImBB,'rx');hold on
    subplot(2,1,2)
    plot(FreqS21(R),angleS21(R),'bx');hold on
    plot(Fresbb,angleBB,'rx');hold on
    %% ---------------- KID Dip ----------------%
    subplot(2,2,2)
    hold on
    
    for p=1:Nbbtemperatures
        plot(KID(n).S21data{p}(:,1),20*log10(KID(n).S21data{p}(:,2)),'-',...
            'color',Tcolors(p,:),'MarkerSize',1) %Plot the resonances in dB space
    end
    p=1;
    plot(KID(n).S21data{p}(:,1),20*log10(KID(n).S21data{p}(:,2)),'-k','LineWidth',2) %Plot the first resonance again
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
    
    %---------------- Laura callibration ----------------%
    subplot(2,2,3)
    hold on
    p=1;
    DX = (KID(n).Fres - KID(n).Fres(1))/KID(n).Fres(1);%response/BW BW=Fres/Q  (Fres in GHz)
    %add the response data
    plot(DX,phaseBBresponse,'or') %
    
    dxLauraResponse = interp1(KIDphase,dF_F0,KID(n).Response(:,2),'spline');%convert theta response to F response
    
    
    %---------------- Various responses ----------------%
    
%     %find the Fres from max (dtheta/dF) assuming dF is the same between
%     %each point
%     SmoothKIDphase=smooth(KIDphase,5);
%     DiffPhase=diff(SmoothKIDphase);
%     tofit=DiffPhase((max(DiffPhase)-min(DiffPhase))/2);
%     figure(2000)
%     plot(dF_F0(tofit),DiffPhase(tofit),'o');hold on
%     
%     g = fittype('a*(x+b)^2+c',...
%             'dependent',{'y'},'independent',{'x'},'coefficients',{'a','b','c'});
%     myfit = fit(dF_F0(tofit),DiffPhase(tofit),g);    
%     
%     plot(myfit)
%     
%     figure(n)
    subplot(2,2,4)
    %fit first 10 pts
    DXfit= polyfit(KID(n).Pbb(1:10)*1e12,DX(1:10),1);
    phasefit= polyfit(KID(n).Pbb(1:10)*1e12,KID(n).Response(1:10,2),1);
    Rfit=polyfit(KID(n).Pbb(1:10)*1e12,KID(n).Response(1:10,3),1);
    dxLaurafit= polyfit(KID(n).Pbb(1:10)*1e12,dxLauraResponse(1:10,1),1);
    
    Norm_dx         =   -DX;%Frequency response in units of MKID BW
    Norm_phaseresp  =   -DXfit(1) * KID(n).Response(:,2)/phasefit(1);
    Norm_Rresp      =   -DXfit(1) * KID(n).Response(:,3)/Rfit(1);
    Norm_Lauraresp  =   -DXfit(1) * dxLauraResponse/dxLaurafit(1);
    
    hold on
    
    plot(KID(n).Pbb*1e12,Norm_dx,'-ok','MarkerFaceColor','k'); 
    plot(KID(n).Pbb*1e12,Norm_Lauraresp,'-ob','MarkerFaceColor','b'); 
    plot(KID(n).Pbb*1e12,Norm_phaseresp,'-or','MarkerFaceColor','r'); 
    plot(KID(n).Pbb*1e12,Norm_Rresp,'-og','MarkerFaceColor','g'); 
    plot(KID(n).Pbb*1e12,-1*polyval(DF_BWfit,KID(n).Pbb*1e12),'-k','LineWidth',2)
    xlabel('Power (pW)')
    ylabel('dFres/Fres ')
    legend('Measured dF','Laura response','phase response','amplitude response','small signal dx extrapolated')
    
    axis tight
    
    Figfile=[S21path,'S21analysis_KID ' num2str(KID(n).KIDnumber) '_' num2str(KID(n).ReadPower(1)) 'dBm.fig'];
    saveas(gcf,Figfile,'fig')

end