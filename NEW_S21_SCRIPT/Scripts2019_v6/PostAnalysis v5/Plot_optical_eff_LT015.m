%This script plots the optical efficiency (from phase readout) per KID
clear all;
resppathyn=[cd filesep '..' filesep '..' filesep '2D_BB' filesep '2D_BB_lambda_2' filesep];
load([resppathyn 'KIDparam.mat'])

figure(1)
for n=1:length(KIDparam)
    opticalEff(n)=KIDparam(n).optphaseeff;
    ID(n)=KIDparam(n).KIDid(1);
    fres(n)=KIDparam(n).fres(1);
end
close all;
subplot(2,1,1)
plot(ID,100*opticalEff,'--o');axis tight
xlabel('KID ID');ylabel('phase optical eficiency %')
hold on

subplot(2,1,2)
[fressorted,fresSI]=sort(fres)
plot(fressorted,100*opticalEff(fresSI),'--o');axis tight
xlabel('fres [GHz]');ylabel('phase optical eficiency %')
hold on


for n=1:length(KIDparam)
    opticalEff(n)=KIDparam(n).optradeff;
    ID(n)=KIDparam(n).KIDid(1);
    fres(n)=KIDparam(n).fres(1);
end

subplot(2,1,1)
plot(ID,100*opticalEff,'--or');axis tight
xlabel('KID ID');ylabel('optical eficiency %')
legend('Phase','Radius')
ylim([0 100])
subplot(2,1,2)
[fressorted,fresSI]=sort(fres)
plot(fressorted,100*opticalEff(fresSI),'--or');axis tight
xlabel('fres [GHz]');ylabel('optical eficiency %')
ylim([0 100])