%after running the pulseanalysis script all info is in the stdout ans. This
%outpout has the same stuff in it as drdtheta.csv. 
%it is used here tp make nice plots


%% NB: rename the output of pulsanaysis to drdtheta
close all;clc
KIDnumbers=[2 3 7];
Kcol=colormap(lines(length(KIDnumbers)));
for m=1:length(KIDnumbers)
KID_id=find(drdtheta(:,1)==KIDnumbers(m));
powers=unique(drdtheta(KID_id,5));
npowers=length(powers);
Pcol=colormap(lines(npowers));

figure(m)
subplot(2,2,1)
for n=1:npowers
   thisP=find(powers(n)==drdtheta(KID_id,5)) ;
   plot(drdtheta(KID_id(thisP),4),1e3*drdtheta(KID_id(thisP),11),'--o',...
       'Color',Pcol(n,:),'MarkerFaceColor',Pcol(n,:),'MarkerEdgeColor',Pcol(n,:));hold on;
   Lstr{n}=['P= -' num2str(powers(n))];
   %get value @ Tmin
   [~,indTmin]=min(drdtheta(KID_id(thisP),4));
   tau_P_Tmin(n,2)=drdtheta(KID_id(thisP(indTmin)),11);
   tau_P_Tmin(n,1)=drdtheta(KID_id(thisP(indTmin)),4);
end
legend(Lstr)
xlabel('T [K]');ylabel('phase lifetime  [msec]')
title(['KID' num2str(KIDnumbers(m))])

subplot(2,2,2)
for n=1:npowers
   thisP=find(powers(n)==drdtheta(KID_id,5)) ;
   plot(drdtheta(KID_id(thisP),4),1e3*drdtheta(KID_id(thisP),13),'--s',...
       'Color',Pcol(n,:),'MarkerFaceColor',Pcol(n,:),'MarkerEdgeColor',Pcol(n,:));hold on;
end
legend(Lstr)
xlabel('T [K]');ylabel('amplitude lifetime  [msec]')


subplot(2,2,3)
plot(-powers,1e-3*tau_P_Tmin(:,2),'--o',...
    'Color',Kcol(m,:),'MarkerFaceColor',Kcol(m,:),'MarkerEdgeColor',Kcol(m,:));hold on;
ylabel('tau [msec]');xlabel('P readout  [dBm]')

subplot(2,2,4)
plot(powers,tau_P_Tmin(:,1),'--o',...
    'Color',Kcol(m,:),'MarkerFaceColor',Kcol(m,:),'MarkerEdgeColor',Kcol(m,:));hold on;
ylabel('T [K] of tau measurement');xlabel('P readout  [dBm]')

end
%% 
