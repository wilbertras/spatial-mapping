%after running the S21 script inport 'allTdep.csv'. 
%This script plots the power
%dependence of the response and Q


%% NB: rename the output of pulsanaysis to drdtheta
close all;clc
KIDnumbers=[3,14];
Kcol=colormap(lines(length(KIDnumbers)));
for m=1:length(KIDnumbers)
KID_id=find(allTdep(:,1)==KIDnumbers(m));
powers=(allTdep(KID_id,8));
Q=(allTdep(KID_id,5));
dxdN=(allTdep(KID_id,14));
npowers=length(powers);

figure(1)
subplot(2,1,1)
plot(powers,Q,'--o',...
    'Color',Kcol(m,:),'MarkerFaceColor',Kcol(m,:),'MarkerEdgeColor',Kcol(m,:));hold on;
xlabel('Pint  [dBm]');ylabel('Q');
Lstr{m}=['KID ' num2str(KIDnumbers(m))];

subplot(2,1,2)
plot(powers,-dxdN,'--o',...
    'Color',Kcol(m,:),'MarkerFaceColor',Kcol(m,:),'MarkerEdgeColor',Kcol(m,:));hold on;
xlabel('Pint  [dBm]');ylabel('|dXdN|');



end
legend(Lstr)