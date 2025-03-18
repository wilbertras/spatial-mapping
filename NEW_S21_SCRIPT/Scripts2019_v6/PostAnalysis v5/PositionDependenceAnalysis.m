%% ini
clc
close all
clear all
datapath='/Volumes/KID/KIDonSun/experiments/ADR general/BoxinBoxsetup/Hybrids/LT013';
%% Read chip design
%Positions.mat has 'Positions' matrix, 4 cols: Fres, X pos, Ypos, index (not used). The
%original chip writing script must be updated to get this correctly. See
%The LT013 directory for example

Positiondata='Design and fab/Positions.mat';
load([datapath filesep Positiondata]);
clear Positiondata 
Positions=sortrows(Positions,1)
%% Read Fsweep
Fsweepfile='Fsweep/S21_130mK_3_5KBBworkspace.mat' ;%.mat after Fsweep analysis
load([datapath filesep Fsweepfile],'freq','magn','magn_syscal','KIDparam');
clear Fsweepfile 

%% plot Fsweep
close all;
figure(1)
subplot(2,1,1)
plot(freq,20*log10(magn_syscal),'-b');hold on
subplot(2,1,2)
plot(freq,20*log10(magn_syscal),'-b');hold on
xlabel('F [GHz]');ylabel('S21');title('green: Fdesign, Blue: Fsweep, Red: mfile')

figure(2);subplot(1,2,1)
plot(KIDparam.fres,'-x');hold on;
subplot(1,2,2)
plot(KIDparam.fres,'-x');hold on;
%% read m file that was used
mfile='selected_many.m' ;%.mat after Fsweep analysis
Mfile=dlmread([datapath filesep mfile]);
clear mfile 
Mfile(1:2,:)=[];
figure(1);
for mmm=1:2
    subplot(2,1,mmm)
    for n=1:length(Mfile(:,1))
        plot([Mfile(n,2) Mfile(n,2)],get(gca,'Ylim'),'-r')
        text(Mfile(n,2),2,num2str(Mfile(n,1)),'HorizontalAlignment','center','color','red');hold on;
    end
end
%% plot Fsweep and original design, with dF factor
dF=1.002
figure(1);subplot(2,1,1)
for n=1:length(Positions(:,3))
    Positions(n,5)=n;%giving positions an increasing index
    plot(dF*[Positions(n,1) Positions(n,1)],get(gca,'Ylim'),'g');hold on;
    text(dF*Positions(n,1),4,num2str(Positions(n,5)),'HorizontalAlignment','center')
end
figure(2);subplot(1,2,1)
hold on;
plot(Positions(:,5),dF*Positions(:,1),'og');hold on;
xlabel('KID number');ylabel('Fres');legend('Measured','Design')
%% Trying to adjust the Spositions col 4 to get 1:1 correspondence with M file
Spositions=Positions;
%OPTIONAL: KID at Fres=5 GHz is due to error. Design positions modified
Spositions(5:end,5)=Spositions(5:end,5)+1;%KID @ 5 GHz is not a KID, we add one in the design
Spositions(31,:)=[];
Spositions(31:end,5)=Spositions(31:end,5)-1;%KID with index 31 is missing most likely

figure(1);subplot(2,1,2)
for n=1:length(Spositions(:,3))
    
    plot(dF*[Spositions(n,1) Spositions(n,1)],get(gca,'Ylim'),'g');hold on;
    text(dF*Spositions(n,1),4,num2str(Spositions(n,5)),'HorizontalAlignment','center')
end
figure(2);subplot(1,2,2)
hold on;
plot(Spositions(:,5),dF*Spositions(:,1),'or');hold on;
xlabel('KID number');ylabel('Fres');legend('Measured','Modelled Reality (for positions)')
%% Read simulation results antenna and plot
%loads in the distance effect .mat file from Reinier and plots the distance
%dependent coupling efficiency using colormap(Jet)
figure(3)
SIMpath='/Volumes/KID/KIDonSun/experiments/ADR general/BoxinBoxsetup/Hybrids/LT020Chip8/LT020 Chip 8 Run Alicia/Scripts 2014_8_v3/CSTbeampatterns/LT013 1mm lens 850 GHz 20x12 holder';
fn='1mmlens_3mmAp_15_4mmdist_OffsetEffect_850GHz.mat';
load([SIMpath filesep fn]);
offeff=Offset(:,:,3);
contourf(OffsetLine*1e3,OffsetLine*1e3,offeff'/max(max(offeff)),[0:0.05:1]);colormap(jet(128));%in mm
colorbar('YTick',[0:0.05:1])%,'TickLength',0.08)
hold on;
% adding the data to the plot
%Read optical efficiency from 2dBB
NEPdata='2D_BB/2D_BB/KIDparam.mat';
load([datapath filesep NEPdata],'KIDparam');
clear NEPdata 
ShiftX=-1.5;ShiftY=-0.8;
posX=Spositions(:,2)/1e3+ShiftX;
posY=Spositions(:,3)/1e3+ShiftY;



bla=colormap(jet(128));
for n=1:length(KIDparam)
    ptu=KIDparam(n).KIDid(1)==Spositions(:,5);
plot(posX(ptu),posY(ptu),'o','MarkerSize',20,...
'MarkerFaceColor',bla(round(128*KIDparam(n).optphaseeff/max([KIDparam.optphaseeff])),:));%,... 
%'MarkerEdgeColor',bla(round(128*KIDparam(n).optphaseeff/max([KIDparam.optphaseeff])),:));hold on;
%text(Positions(n,2)-0.1,Positions(n,3),num2str(KIDparam(n).optphaseeff/max([KIDparam.optphaseeff]),3));
title(['Peak efficiency: ' num2str(max([KIDparam.optphaseeff]))]);
end