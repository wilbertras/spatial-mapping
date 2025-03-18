function [arrayparam,KIDparam] = FscanAnalysis_v1(inputmat,fileparam)
% NewestVersion JB 22-5-13
% Figure 1 is not saved and not overwritten, usefull to combine bla
% Function finds all KID dips within scan range, determines the parameters
% of all KIDs, but does not fit the dips (for speed). Outputs some plots and
% structs with parameters. Plots are saved automatically.

% runs on F[Ghz] Re Im or F [Ghz] dBmagn, and must be told which format is used
%
% INPUT
% 2 options:
% 1 - use the arguments of the function
% inputmat has 3 columns: Freq [GHz] magn(dB) phase (which is not used) or
% Freq [GHz] re Im.
% fileparam is a text string for easy reference
% of the figures that are saved. data is saved in program root
%
% 2 - load file with coirrect format (and no header lines) in directly
% filename and dir specified in fisrt cell
%
% OUTPUT
%
% arrayparam is a struct with all parameters of the array (mean Q, F spacing and so on)%
%
% KIDparam is a 2D struct with all paramaeters of each individual KID (less
% usefull in this context)%
%
% fig file and 2 .csv files containing all info from arrayparam and
% KIDparam structs. Saved into the script root (called with arguments)
% or in the specified dir (when called without arguments)
%
% NB: after running [data,power]=importFscanGroningen you can analyze for
% example power 5 like FscanAnalysis_magndBonly_1(data{5},'power5')
%% SET VARIABLES
close all;
logdB=1;    %=1 (default) if supplied with Freq [GHz] magn(dB) phase (phase col not used and does not need to be there).
% use 0 for Groningen data as this data is Freq [GHz] re Im given
BWcon=0.5;  % fraction of BW that results in bad KID if spacing is smaller than this value

% give file and dir for firect data import, ignored if fy is called with arguments
filename='S21_130mK_3_5KBB.dat';
dir='/Volumes/KID/KIDonSun/experiments/ADR general/BoxinBoxsetup/Hybrids/LT013/Fsweep/';%directory incl \ at end

%variables for the hystograms
histdf=0.1;         %DEFAULT=0.5. df for Histgram plot in MHz. 0.1 for high Q>100.000
histmax=30;         %DEFAULT=20 maximum df value in the histogram (in MHz)

%variables to find the dips (might need some playing with)
diff2treshold=40e6;  %80for hot ones use 0.4e8 100 mK use 1e8;  %Typ. 0.5-5e8. diff2treshold is treshold value in second derivation to find KIDs. 5e8 or so %

%fort poltting
plotdots=1;

%variables that are normaly ok
difsmooth=10;       %Typ. 10. difsmooth is a smoothing in the derivative fundtion used to find the KID peaks, typ=10-100%
smp=0;              %DEFAULT=0. optional smoothing in the amplitude ripple removal. must be odd or 0. If zero no smoothing
rangefactor=8;      %DEFAULT=12. sets addn width of KID peaks to be reomved from data prior to S21 cal
maxnumpeaks=5000;   %DEFAULT=5000. max number of peaks that are allowed to be detected
thres_sample=4;     %DEFAULT=4. Thres_sample is number of samples needed in dip to make it valid set to default 4
%close all;
%% data read. Fy expects 3 cols without headers F[Ghz] S21 [dB] rad (not used)

format long g;

datafile=[dir filename]
global S21data;

if nargin==2
    S21data=inputmat;
    dir=cd('../');%data outputted to current Dir of fy
    
elseif nargin==0
    fid=fopen(datafile);
    if fid==-1 %%check for exists
        error('tsss, no file given in fy');
    else
        S21data=cell2mat(textscan(fid,'%f%f%f','headerlines',0));%read in KID number and dRdtheta
        fclose(fid);
    end
    fileparam=filename(1:end-4);
else
    error('not correct # iputs')
    
end
%% Data initialization
%converting to logdBif required
if logdB==0
    S21data(:,2)=10*log10(S21data(:,2).^2+S21data(:,3).^2);
    disp('data converted from re Im to log dB')
elseif logdB==1
    disp('data in F[GHz] dB magn phase')
end



[bla,b]=unique(S21data(:,1));
good= S21data(b,:);
freq  = good(:,1);
magn  = 10.^(good(:,2)/20);
clear S21data bla b good;

%numerical differentiation of data
der=numdif([freq 20*log10(magn)],difsmooth);
der2=numdif(der,difsmooth);
clear der;

%plotting 2nd derivative and range;

figure(2)
%subplot(1,2,2)
hold on;
plot(der2(:,1),der2(:,2))
grid on
hold on
plot(der2(:,1),ones(length(der2(:,1)),1)*diff2treshold,'r')
xlabel('Frequency [GHz]')
ylabel('Second Derivative [dB]')
title('Second Derivative [dB] and threshold (thres\_diff2)')
legend('Second Derivative [dB]','diff2treshold')
figure(1)
%subplot(1,2,1)
plot(freq,20*log10(magn),'k');hold on;
figure(10)
%subplot(1,2,1)
hold on;
plot(freq,20*log10(magn),'k');hold on;

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


%% Use rangefactor to remove data where KID dips are and store the indices in range
% NB: Now we use double derivative part*rangefactor to get the range over
% which a KID is present. No KID BW is used because this makes software
% slow
range    = [];  % define empty
numpeaks
sample_up = 0;
sample_down = 0;
ind=zeros(numpeaks,1);peaks_10BW=cell(numpeaks,1);listFres=zeros(numpeaks,1);BW=zeros(numpeaks,1);
%Removing data with double der above threshold (i.e. where the KIDs are) to be able to smooth the S21 data
for c=1:numpeaks
    [S21minRAW(c),index]  = min(magn(peak{c}));
    Windex = length(peak{c});    % #points within KID dip as found by der2
    % multiply width left and right rangefactor times bandwidth
    range_start      =  peak{c}(index)   - round(rangefactor*Windex);%beginning of KID expanded downwards
    range_end        =  peak{c}(index)   + round(rangefactor*Windex);
    
    if (range_start <1)   % make sure index higher than 1
        range_start = 1;
    end
    if (range_end > length(magn))   % make sure index max equal to length
        range_end = length(magn);
    end
    
    peaks_10BW{c} = (range_start:range_end);%range within rangefactor*KIDBW, i.e. where KID modifies system transmission
    if (c == 1)  % range does not yet exits
        range(1:length(peaks_10BW{c})) =  peaks_10BW{c};
    else
        range(end+1:end+length(peaks_10BW{c})) =  peaks_10BW{c};  % all indices in data where KIDs are (sum of all peaks_10BW)
    end
end

%% Smoothes S21 data (removes ripple)
%
% range gives the indices containing KID -> to be thrown away

index_cal        = (1:length(magn));%all datapoints
if (not(isempty(range)))  % peaks found
    index_cal(range) = [];% remove data containing KID
end


if length(index_cal)<2
    error(['No data without KIDs, please reduce rangefactor (in program) below '  num2str(rangefactor)])
end

%amp correction: creating smoothed data and devide by it (in normal
if smp>0
    rippletemp = smooth(freq(index_cal),magn(index_cal),smp); %data smoothed
    ripple=interp1(freq(index_cal),rippletemp,freq,'linear','extrap');%smoothed data interpolated (to bridge gaps were KIDs ar removed)%
    magn_syscal  = magn  ./ ripple;
else %no smoothing
    ripple=interp1(freq(index_cal),magn(index_cal),freq,'linear','extrap');%data interpolated (to bridge gaps were KIDs ar removed)%
    magn_syscal  = magn  ./ ripple;
end
clear ripple rippletemp;
%% Finding KID parameters and the rest of the callibration parameters
%NB: First we get 3dB BW and the rest using system callibrated data
S21_3dBpt=zeros(0,numpeaks);
BWindex=zeros(0,numpeaks);
indices3dB=zeros(numpeaks,2);
for c=1:numpeaks %finding parameters
    [S21min,ind(c)]     = min(magn_syscal(peak{c}));    % finds dip in in peak magn data
    S21max              = mean(magn_syscal(index_cal)); %S21 off resonance is mean of S21 data that does not have KIDs in it.
    Fres                = freq(peak{c}(ind(c)));
    listFres(c)         = Fres;
    
    %S21_3dBpt(c)        =
    %20*log10((10^(S21min/10)+10^(S21max/10))/2);%finds half of KID dip (3dB for deep KID). 3dB=0.5 in power
    S21_3dBpt(c)        = sqrt((S21min^2+S21max^2)/2);%finds half of KID dip (3dB for deep KID). 3dB=0.5 in power
    
    %finding KID BW here lower F edge KID 3dB pts
    for jj=1:length(magn_syscal)-peak{c}(ind(c))
        if magn_syscal(peak{c}(ind(c))+jj)>S21_3dBpt(c)
            sample_up = jj; %is used to obtain ranges, highest index NOT kid affected, started from Dip center
            break
        end
    end
    %finding high f point dip
    for jj=1:peak{c}(ind(c))-1
        if magn_syscal(peak{c}(ind(c))-jj)>S21_3dBpt(c)
            sample_down = jj;
            break
        end
    end
    
    BWindex(c) = (sample_up + sample_down);           %full KID dip width in points
    BW(c) = BWindex(c) * (freq(2)-freq(1));             %the total KID dip width in frequency (points * dfper point), this is the bandwidth
    indices3dB(c,1)=peak{c}(ind(c))-sample_down;        %is the index at low F 3dB point of dip (c), used for plotting
    indices3dB(c,2)=peak{c}(ind(c))+sample_up;          %is the index at high F 3dB point of dip (c). used for plotting
end

for c=1:numpeaks
    if foundpeak
        KIDparam.BW(c)    =  BW(c);
        KIDparam.S21min(c)=  magn_syscal(peak{c}(ind(c)));  %S21 min amplitude
        KIDparam.S21minRAW(c)=S21minRAW(c);
        KIDparam.fres(c)  =  freq(peak{c}(ind(c)));
        if c>1
            KIDparam.df(c)=abs(KIDparam.fres(c)-KIDparam.fres(c-1));
        else
            KIDparam.df(c)=0;
        end
        KIDparam.Q(c)  =  KIDparam.fres(c)/KIDparam.BW(c);
        KIDparam.ind(c)   =  ind(c); % index in dataset where closest to true resonance frequency of the KID in scan center
        KIDparam.Qi(c)=KIDparam.Q(c)/KIDparam.S21min(c);
        KIDparam.Qc(c)=KIDparam.Qi(c)*KIDparam.Q(c)/(KIDparam.Qi(c)-KIDparam.Q(c));  %Get Qc=QiQ/(Qi-Q)
        
    else
        error('scan without peaks');
    end
end
%calculates means
[bla,blaind]=find(KIDparam.Qc<0);%invalid Qc
crazykids=length(bla);
if crazykids>0
    disp(['negative Qc found: number= ' num2str(length(bla))])
end
[~,blaind]=find(KIDparam.Qc>0);%invalid Qc
arrayparam.meanQ=mean(KIDparam.Q(blaind));
arrayparam.meanBW=mean(KIDparam.BW(blaind));
arrayparam.stdQ=std(KIDparam.Q(blaind));
arrayparam.meanQc=mean(KIDparam.Qc(blaind));
arrayparam.stdQc=std(KIDparam.Qc(blaind));
arrayparam.meanQi=mean(KIDparam.Qi(blaind));
arrayparam.stdQi=std(KIDparam.Qi(blaind));
arrayparam.meanS21dB=20*log10(mean(KIDparam.S21min(blaind)));
arrayparam.stdS21dB=20*log10(std(KIDparam.S21min(blaind)));

[~,blaind]=find(KIDparam.Qc(2:end)>0);%invalid Qc 2:end for dF
arrayparam.meandf=mean(KIDparam.df(blaind));
arrayparam.stddf=std(KIDparam.df(blaind));
arrayparam.df10MHz=length(find(KIDparam.df(blaind)<0.01))-1;%number dF below 5 MHz. -1 beacuse df first KID with itself=0.
arrayparam.df5MHz=length(find(KIDparam.df(blaind)<0.005))-1;%number dF below 5 MHz. -1 beacuse df first KID with itself=0.
arrayparam.df1MHz=length(find(KIDparam.df(blaind)<0.001))-1;%number dF below 1 MHz. -1 beacuse df first KID with itself=0.
arrayparam.df02MHz=length(find(KIDparam.df(blaind)<0.0002))-1;%number dF below 0.2 MHz. -1 beacuse df first KID with itself=0.

%dead pixel analysis
st=KIDparam.df(2:end)<BWcon*KIDparam.BW(2:end);
m=0;cm=0;
for nn=1:length(st)
    if st(nn)==1
        cm=cm+1;m=1;
        KIDparam.goodKID(nn)=0;
    elseif st(nn)==0 && m==1
        cm=cm+1;m=0;
        KIDparam.goodKID(nn)=0;
    elseif st(nn)==0 && m==0
        KIDparam.goodKID(nn)=1;
    end
end
if st(nn)==1 %add 1 if last KID was too close
    cm=cm+1;
    KIDparam.goodKID(nn+1)=0;
else
    KIDparam.goodKID(nn+1)=1;
end


disp(['# pixels within ' num2str(BWcon) ' x BW = ' num2str(cm) ' for n= ' num2str(numpeaks)])


%format short G
fprintf('# usable tonesfound = %.0f  (not too close)\n', numpeaks-crazykids-cm);
fprintf('#tonesfound = %.0f \n', numpeaks);
fprintf('df: %.1f MHz pm %.1f \n',arrayparam.meandf*1000,arrayparam.stddf*1000)
fprintf('Q: %.3f pm %.3f \n', arrayparam.meanQ/1e3, arrayparam.stdQ/1e3)
fprintf('BW: %.3f \n', arrayparam.meanBW*1e9)
fprintf('Qi: %.3f pm %.3f \n', arrayparam.meanQi/1e3, arrayparam.stdQi/1e3)
fprintf('Qc: %.3f pm %.3f \n', arrayparam.meanQc/1e3, arrayparam.stdQc/1e3)
fprintf('S21[dB]: %.1f pm %.1f \n', arrayparam.meanS21dB, arrayparam.stdS21dB)
fprintf('# tones with df < 1MHz: %.0f \n', arrayparam.df1MHz)
fprintf('# tones with df < 5MHz: %.0f \n', arrayparam.df5MHz)
fprintf('# tones with df < 10MHz: %.0f \n', arrayparam.df10MHz)


%% plot rest Fig 1
%if saveall > 0
figure(10)
hold on
%subplot(1,2,1)
%for legend
plot(freq(peaks_10BW{1}),20*log10(magn(peaks_10BW{1})),'b');hold on;%data to be fitted   
plot(KIDparam.fres(1),20*log10(KIDparam.S21minRAW(1)),'or','MarkerFaceColor','r','MarkerSize',4);hold on;
plot(KIDparam.fres(1),20*log10(KIDparam.S21minRAW(1)),'ok','MarkerFaceColor','k','MarkerSize',4);hold on;
legend('Measured magnitude','KIDS found','Bad KID','Usable KID')
for c=1:numpeaks
    plot(freq(peaks_10BW{c}),20*log10(magn(peaks_10BW{c})),'b');hold on;%data to be fitted   
end

for c=1:numpeaks
    if KIDparam.goodKID(c)==0
        plot(KIDparam.fres(c),20*log10(KIDparam.S21minRAW(c)),'or','MarkerFaceColor','r','MarkerSize',4);hold on;
    else
        plot(KIDparam.fres(c),20*log10(KIDparam.S21minRAW(c)),'ok','MarkerFaceColor','k','MarkerSize',4);hold on;
    end
end
title([num2str(numpeaks) ' dips found, of which ' num2str(numpeaks-crazykids-cm) 'usable, rest df <' num2str(BWcon) ' bandwidths'])
grid on
xlabel('Frequency [GHz]')
ylabel('Magnitude [dB]')



figure(2)
for c=1:numpeaks
    plot(der2(peaks_10BW{c},1),der2(peaks_10BW{c},2),'.b');hold on;%data to be fitted)
    
   
    
end


figure(4)
subplot(2,2,1)
semilogy(KIDparam.Qi,'bx');
hold on;
semilogy(KIDparam.Q,'ro');
semilogy(KIDparam.Qc,'k.');
xlabel('index')
ylabel('Q factors')
legend('Qi','Q','Qc')
subplot(2,2,2)
hist(KIDparam.df*1000,0:histdf:histmax)
xlim([0 histmax])
xlabel('df [MHz]')
ylabel('count')
title(['mean ' num2str(arrayparam.meandf*1000)  'std '  num2str(arrayparam.stddf*1000)])
subplot(2,2,3)
hist((KIDparam.Qc),histmax)
xlabel('Qc ')
ylabel('count')
title(['mean ' num2str(arrayparam.meanQc)  'std '  num2str(arrayparam.stdQc)])
subplot(2,2,4)
hist(20*log10(KIDparam.S21min),-40:2:0)
xlim([-40 0])
xlabel('S21min [dB]')
ylabel('count')
title(['mean ' num2str(arrayparam.meanS21dB)  'std '  num2str(arrayparam.stdS21dB)])

figure(3)
plot(freq,20*log10(magn_syscal),'-b');hold on
if plotdots==1
for c=1:numpeaks
    if KIDparam.goodKID(c)==0
        plot(KIDparam.fres(c),20*log10(KIDparam.S21min(c)),'or','MarkerFaceColor','r','MarkerSize',4);hold on;
    else
        plot(KIDparam.fres(c),20*log10(KIDparam.S21min(c)),'ok','MarkerFaceColor','k','MarkerSize',4);hold on;
    end
    %plot(freq(indices3dB(c,:)),20*log10(S21_3dBpt(c)),'go','MarkerFaceColor','g','MarkerSize',2);%
end
end
title(['# peaks: ' num2str(numpeaks) '. OK KIDs have black dot'  ]);

saveas(figure(10),[dir fileparam 'Figure1'  '.fig']);
saveas(figure(2),[dir fileparam 'Figure2'  '.fig']);
saveas(figure(3),[dir fileparam 'Figure3'  '.fig']);
saveas(figure(4),[dir fileparam 'Figure4'  '.fig']);
dlmwrite([dir fileparam 'KIDparam.csv'], cell2mat(struct2cell(KIDparam)), 'precision', '%.9g', 'delimiter', ',', 'newline', 'pc');
keeswrite('Q, stdQ, Qc, stdQc, Qi, stdQi, S21, stdS21, df, stddf, df10MHz, df5 MHz, df1 MHz, df 0.2MHz',[dir fileparam 'arrayparam.csv']);
dlmwrite([dir fileparam 'arrayparam.csv'], cell2mat(struct2cell(arrayparam))', '-append', 'precision', '%.9g', 'delimiter', ',', 'newline', 'pc');

fid=fopen([dir fileparam 'textoutput.txt'],'wt');


fprintf(fid,'#tonesfound = %.0f \n', numpeaks);
fprintf(fid,'# usable tonesfound = %.0f  (not too close)\n', numpeaks-crazykids-cm);
fprintf(fid,'# dead pixels within %0.f BW =%.0f \n', BWcon,cm);
fprintf(fid,'df: %.1f MHz pm %.1f \n',arrayparam.meandf*1000,arrayparam.stddf*1000);
fprintf(fid,'Q: %.3f pm %.3f \n', arrayparam.meanQ/1e3, arrayparam.stdQ/1e3);
fprintf(fid,'BW: %.3f \n', arrayparam.meanBW*1e9);
fprintf(fid,'Qi: %.3f pm %.3f \n', arrayparam.meanQi/1e3, arrayparam.stdQi/1e3);
fprintf(fid,'Qc: %.3f pm %.3f \n', arrayparam.meanQc/1e3, arrayparam.stdQc/1e3);
fprintf(fid,'S21[dB]: %.0f pm %.0f \n', arrayparam.meanS21dB, arrayparam.stdS21dB);
fprintf(fid,'# tones with df < 0.2MHz: %.0f \n', arrayparam.df02MHz);
fprintf(fid,'# tones with df < 1MHz: %.0f \n', arrayparam.df1MHz);
fprintf(fid,'# tones with df < 5MHz: %.0f \n', arrayparam.df5MHz);
fprintf(fid,'# tones with df < 10MHz: %.0f \n', arrayparam.df10MHz);
fclose(fid);
save([dir fileparam 'workspace'])
end

function dif=numdif(raw,span)
%raw is colx,coly. returns num differential in dif (x, dy/dx)
%span is usuall span parameter in smooth function

y=raw(:,2);
dx=raw(2,1)-raw(1,1);
dif=zeros(length(y),1);
dif(:,1)=raw(:,1);
dif(2:end,2)=diff(y);
%equal ends
dif(1,2)=dif(2,2);dif(length(y),2)=dif(length(y)-1,2);
%devide by dx
dif(:,2)=dif(:,2)/(dx);
%smooth
dif(:,2)=smooth(dif(:,2),span);

end


function keeswrite(strmat,filename)
%CASEWRITE Writes casenames from a string matrix to a file.
%   CASEWRITE(STRMAT,FILENAME) writes a list of names to a file, one per line.
%   FILENAME is the complete path to the desired file. If FILENAME does not
%   include directory information, the file will appear in the current directory.
%
%   CASEWRITE with just one input displays the File Open dialog box allowing
%   interactive naming of the file.

%   Copyright 1993-2004 The MathWorks, Inc.
%   $Revision: 2.9.2.1 $  $Date: 2003/11/01 04:25:22 $

if (nargin == 0)
    error('stats:casewrite:TooFewInputs',...
        'CASEWRITE requires at least one argument.');
end
if nargin == 1
    [F,P]=uiputfile('*');
    filename = [P,F];
end
fid = fopen(filename,'wt');

if fid == -1
    disp('Unable to open file.');
    return
end

if strcmp(computer,'MAC2')
    lf = setstr(13);
else
    lf = setstr(10);
end

lf = lf(ones(size(strmat,1),1),:);
lines  = [strmat lf]';

fprintf(fid,'%s',lines);
fclose(fid);
end
