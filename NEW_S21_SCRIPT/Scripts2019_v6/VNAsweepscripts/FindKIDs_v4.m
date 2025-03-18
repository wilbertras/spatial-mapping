function [blindrange,fres_ind,fres] = FindKIDs_v4
% From V4 of the original Fscan analysis script. Simplified for MUX people
% blindrange = logical on the original frequency array where blind MKIDs can be placed%
% fres_ind = logical on the original frequency array where the MKID resonances are (tones must be placed here)%
% fres = frequency array where the MKID resonances are (tones must be placed here)%
%
% script uses 2 thesholds, one for the MKIDs and one for the blind tone
% placement that can be set manual or automatic. 
%% SET VARIABLES
close all;clear variables;

%%%%%%%%%%% give file and dir for data import %%%%%%%%%%%
filename='Fsweep_86dBm.dat';
dir='/Volumes/KID/KIDonSun/experiments/Entropy ADR/LT102MIXED/singletone/';%directory incl filesep at end

%%%%%%%%%%%  variables to find the dips  %%%%%%%%%%%
diff2treshold=10e8;         %Treshold in second derivative to find the peaks. Set to 0 for automatic
diff2treshold_blinds= 1e8;  %Treshold in second derivative to find the peaks. Set to 0 for automatic
BW_blinds = 1;              %BW in MHz around each found resonance that we do not want to place blind MKIDs
%variables that are normaly ok (for dip finding) but must be settable
difsmooth=3;        %Typ. 10. span in smooth function typ=10-100
maxnumpeaks=5000;   %DEFAULT=5000. max number of peaks that are allowed to be detected (to limit memory use)
thres_sample=2;     %DEFAULT=2. Thres_sample is number of samples above threshold in second der to make datarange qualify as having a KID resonance

%% data read. Fy expects 3 cols without headers: F[Ghz] S21 [dB] rad (not used)

format long g;
datafile=[dir filename];
global S21data;

fid=fopen(datafile);
if fid==-1 %%check for exists
    error('tsss, no file given in fy');
else
    S21data=cell2mat(textscan(fid,'%f%f%f','headerlines',0));%read in KID number and dRdtheta
    fclose(fid);
end

%% Data initialization

[~,b]=unique(S21data(:,1)); %GHz
good= S21data(b,:);
freq  = good(:,1);
magn  = 10.^(good(:,2)/20); %data in dB
clear S21data b good;
dF=freq(2)-freq(1);%frequency difference adjecent points (GHz)
BW_blinds=round((BW_blinds*1e-3)/dF); %convert BW blinds into points

%numerical differentiation of data
der2 = numddif([freq 20*log10(magn)],difsmooth);

%find diff2treshold if set to 0
if diff2treshold == 0
    diff2treshold = 20*std(der2(1:100,2));    %20 sigma treshold
    disp(['diff2treshold set to: ' num2str(diff2treshold,3)])
else
    disp(['diff2treshold used from input: ' num2str(diff2treshold,3)])
end
%blinds
if diff2treshold_blinds == 0
    diff2treshold_blinds = 5*std(der2(1:100,2));    %15 sigma treshold
    disp(['diff2treshold_blinds set to: ' num2str(diff2treshold_blinds,3)])
else
    disp(['diff2treshold_blinds used from input: ' num2str(diff2treshold_blinds,3)])
end


%% logic to find peaks
c=0;cactive=0;j=1;  %c counts peaks, cactive is switch, j counts # points above threshold
foundpeak=0;        %remebers if there is 1 peak or more found
peak=cell(1,floor(length(der2(:,2))/(thres_sample+1)));     %ini of peak to maximum possible number
for n=1+floor(2.5*difsmooth):length(der2(:,2))-floor(2.5*difsmooth) % we do not look at the extreme start and end of der2 as it spikes there
    if der2(n,2)>diff2treshold
        foundpeak=1;%peak found
        if cactive==0  % only update if new peak
            c=c+1;
            cactive=1;  % now we know we have detected 1 point in 1 peak
        end
        peak{c}(j)=n;%contains a list of indices around each peak
        j=j+1; %numpoints(c)=j;   %
    else %data smaller than threshold
        if cactive == 1 && j <= thres_sample %if last peak found has less than thres_sample datapoints -> peak invalid
            c = c-1; % overwrite data
        end
        cactive=0;j=1;  %if peak is valid enable detection new peak
    end
end

if foundpeak
    numpeaks = c;
    peak(numpeaks+1:end)=[];
    if numpeaks>maxnumpeaks
        error([ num2str(numpeaks) ' > ' num2str(maxnumpeaks) ' peaks found, adjust diff2treshold parameter if too low, or adjust difsmooth (in function call) or maxnumpeaks (in program)']);
    end
else
    error('No peaks found, adjust diff2treshold parameter');
end

%% Finding KID Resonance frequencies from the minimum value within the found 'peak' ranges range
ind=zeros(0,numpeaks);S21minRAW=zeros(0,numpeaks);fres=zeros(0,numpeaks);fres_ind=zeros(0,numpeaks);
blindrange1 = ones(length(freq),1);
for c=1:numpeaks %finding parameters
    [S21minRAW(c),ind(c)]     = min(magn(peak{c}));    % finds dip in in peak magn data
    fres_ind(c) = peak{c}(ind(c));
    fres(c)    = freq(fres_ind(c));
    blindrange1(fres_ind(c)-BW_blinds:fres_ind(c)+BW_blinds)=0;
end

blindrange2 = diff2treshold_blinds > abs(der2(:,2));

blindrange = logical(blindrange2) & logical(blindrange1);

%% plot
%plotting F sweep
figure(1)
subplot(2,1,1)
hold on;
plot(freq,20*log10(magn),'k');	%raw data
plot(fres,20*log10(S21minRAW),'or','MarkerFaceColor','r','MarkerSize',4);% found peaks
plot(freq(logical(blindrange)),20*log10(magn(logical(blindrange))),'.b','MarkerSize',9);
legend('S21  (dB)','KID resonances','blind tone area')
title([num2str(numpeaks) ' dips found'])
grid on
xlabel('Frequency (GHz)')
ylabel('Magnitude (dB)')

%plotting 2nd derivative and range;
subplot(2,1,2)
plot(der2(:,1),der2(:,2),'k')
grid on; hold on
plot(der2(:,1),ones(length(der2(:,1)),1)*diff2treshold,'r')
plot(der2(:,1),ones(length(der2(:,1)),1)*diff2treshold_blinds,'b')
%plot(der2(logical(blindrange1),1),der2(logical(blindrange1),2),'.k','MarkerSize',9);
%plot(der2(logical(blindrange2),1),der2(logical(blindrange2),2),'.r','MarkerSize',9);
plot(der2(logical(blindrange),1),der2(logical(blindrange),2),'.b','MarkerSize',9);
plot(der2(:,1),ones(length(der2(:,1)),1)*-1*diff2treshold_blinds,'b')
xlabel('Frequency [GHz]')
ylabel('Second Derivative')
legend('Second Derivative',['diff2treshold: ' num2str(diff2treshold,3)],['diff2treshold blinds: ' num2str(diff2treshold_blinds,3)],'blind tone area' )

end


function ddif=numddif(raw,span)
%raw is colx,coly. returns num differential in dif (x, dy/dx)
%span is usuall span parameter in smooth function
dif=zeros(length(raw),3);ddif=zeros(length(raw),3);
%
y=raw(:,2);
dx=raw(2,1)-raw(1,1);
dif(:,1)=raw(:,1)-dx/2;%correct x value
dif(2:end,3)=diff(y);
%'invent' first point
dif(1,3)=dif(2,3);
%devide by dx
dif(:,3)=dif(:,3)/(dx);
%smooth
dif(:,2)=smooth(dif(:,3),span);

%double der
ddif(:,1)=dif(:,1)-dx/2;%correct x value
ddif(2:end,3)=diff(dif(:,2));
%'invent' first point
ddif(1,3)=ddif(2,3);
%devide by dx
ddif(:,3)=ddif(:,3)/(dx);
ddif(:,2)=smooth(ddif(:,3),span);
disp('bla')
end


