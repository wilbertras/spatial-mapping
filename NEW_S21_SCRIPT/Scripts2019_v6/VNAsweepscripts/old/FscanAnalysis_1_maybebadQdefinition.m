function [arrayparam,KIDparam] = FscanAnalysis_1(S21data,diff2treshold,difsmooth,fileparam)

% Function finds all KID dips within scan range, determines the parameters
% of all KIDs (using estimation of BW, not a fit) and outputs some plots and 
% structs with parameters. Plots are saved automatically.
% Program is used to analysis a single F scan of a MKID array.
% version 2 has much simpler 3dB BW finding (all error catching removed)
% NB: KID param found from data, so too course scan gives bad results. Addn.figure 3 should have a flat baseline otherwise
% KID parameters not valid as the mean baseline without KIDs is used as
% reference to get KID dip depth and BW.
%
% INPUT
%
% S21data has 3 columns: Freq [GHz] Re Im
%
% diff2treshold is treshold value in second derivation to find KIDs. 5e8 or so %
%
% difsmooth is a smoothing in the derivative fundtion used to find the KID peaks, typ=10-100%
%
% fileparam is a text string (might contain the chip ID) for easy reference
% of the figures that are saved
%
% OUTPUT
%
%arrayparam is a struct with all parameters of the array (mean Q, F spacing and so on)%
%
%KIDparam is a 2D struct with all paramaeters of each individual KID (less usefull in this context)%

%% Set variables
smp=0;              %best=0, otherwise program slow. optional smoothing in the amplitude ripple removal. must be odd or 0. If zero no smoothing
rangefactor=12;     %=12. total area around 1 KID marked as KID modified is set by this parameter. Not critical 
maxnumpeaks=5000;   %max number of peaks that are allowed to be detected
thres_sample=4;     %=4Thres_sample is number of samples needed in 1 dip to mark it as a KID.

%% Data initialization
[bla,b]=unique(S21data(:,1));
good= S21data(b,:);
freq  = good(:,1);
Z=complex(good(:,2),good(:,3));
magn  = abs(Z);
clear S21data bla b good;

%numerical differentiation of data
der=numdif([freq 20*log10(magn)],difsmooth);
der2=numdif(der,difsmooth);
clear der;

%plotting 2nd derivative and range;

figure(1)
subplot(2,1,2)
hold on;
plot(der2(:,1),der2(:,2))
grid on
hold on
plot(der2(:,1),ones(length(der2(:,1)),1)*diff2treshold,'r')
xlabel('Frequency [GHz]')
ylabel('Second Derivative [dB]')
title('Second Derivative [dB] and threshold (thres\_diff2)')
legend('Second Derivative [dB]','diff2treshold')
hold on
subplot(2,1,1)
plot(freq,20*log10(magn),'b');


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
%Removing data with double der above threshold to be able to smooth the S21 data
for c=1:numpeaks
    [S21min,index]  = min(magn(peak{c}));
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
    
    BWindex(c) = (sample_up + sample_down);           %full of total KID dip width in points
    BW(c) = BWindex(c) * (freq(2)-freq(1));             %half of the total KID dip width in frequency, this is the bandwidth
    indices3dB(c,1)=peak{c}(ind(c))-sample_down;        %is the index at low F 3dB point of dip (c), used for plotting
    indices3dB(c,2)=peak{c}(ind(c))+sample_up;          %is the index at high F 3dB point of dip (c). used for plotting
end

for c=1:numpeaks
    if foundpeak
        KIDparam.BW(c)    =  BW(c);
        KIDparam.S21min(c)=  magn_syscal(peak{c}(ind(c)));  %S21 min amplitude
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
arrayparam.meanQi=mean(KIDparam.Qi);
arrayparam.stdQi=std(KIDparam.Qi);
arrayparam.meanS21dB=mean(20*log10(KIDparam.S21min));
arrayparam.stdS21dB=std(20*log10(KIDparam.S21min));
arrayparam.meandf=mean(KIDparam.df(2:end));
arrayparam.stddf=std(KIDparam.df(2:end));
arrayparam.df1MHz=length(find(KIDparam.df<0.001))-1;%number dF below 1 MHz. -1 beacuse df first KID with itself=0.
arrayparam.df02MHz=length(find(KIDparam.df<0.0002))-1;%number dF below 0.2 MHz. -1 beacuse df first KID with itself=0.
disp(['#tonesfound = ' num2str(numpeaks)]);
disp(['df: ' num2str(arrayparam.meandf*1000)  'pm '  num2str(arrayparam.stddf*1000)])
disp(['Qi: ' num2str(arrayparam.meanQi/1e6)  'pm '  num2str(arrayparam.stdQi/1e6)])
disp(['S21[dB]: ' num2str(arrayparam.meanS21dB)  'pm '  num2str(arrayparam.stdS21dB)])
disp(['# tones with df < 1MHz: ' num2str(arrayparam.df1MHz)])
disp(['# tones with df < 0.2MHz: ' num2str(arrayparam.df02MHz)])



%% plot rest Fig 1
%if saveall > 0
    figure(1)
    hold on
    subplot(2,1,1)
    for c=1:numpeaks
        plot(freq(peaks_10BW{c}),20*log10(magn(peaks_10BW{c})),'.r');hold on;%data to be fitted
    end
    
    grid on
    xlabel('Frequency [GHz]')
    ylabel('Magnitude [dB]')
    title('Normalized magnitude (dBm) and KIDS found 10BW')
    legend('Measured magnitude','KIDS found')
   
    subplot(2,1,2)
    for c=1:numpeaks
        plot(der2(peaks_10BW{c},1),der2(peaks_10BW{c},2),'.b');hold on;%data to be fitted)
        
    end


    figure (2)
    subplot(2,2,1)
    semilogy(KIDparam.Qi,'bx');
    hold on;
    semilogy(KIDparam.Q,'ro');
    semilogy(KIDparam.Qc,'k.');
    xlabel('index')
    ylabel('Q factors')
    legend('Qi','Q','Qc')
    subplot(2,2,2)
    hist(KIDparam.df*1000,0:0.5:20)
    xlim([0 20])
    xlabel('df [MHz]')
    ylabel('count')
    title(['mean ' num2str(arrayparam.meandf*1000)  'std '  num2str(arrayparam.stddf*1000)])
    subplot(2,2,3)
    hist(KIDparam.Qi/1e6,20)
    xlabel('Qi [10^6]')
    ylabel('count')
    title(['mean ' num2str(arrayparam.meanQi/1e6)  'std '  num2str(arrayparam.stdQi/1e6)])
    subplot(2,2,4)
    hist(20*log10(KIDparam.S21min),-40:2:0)
    xlim([-40 0])
    xlabel('S21min [dB]')
    ylabel('count')
    title(['mean ' num2str(arrayparam.meanS21dB)  'std '  num2str(arrayparam.stdS21dB)])

    figure(3)
    plot(freq,20*log10(magn_syscal),'-b');hold on
    for c=1:numpeaks
        plot(KIDparam.fres(c),20*log10(KIDparam.S21min(c)),'or','MarkerFaceColor','b','MarkerSize',2);
        plot(freq(indices3dB(c,:)),20*log10(S21_3dBpt(c)),'go','MarkerFaceColor','g','MarkerSize',2);%
    end
    title(['# peaks: ' num2str(numpeaks)  ]);
    %if saveall==2
        saveas(figure(1),['Figure1' fileparam '.fig']);
        saveas(figure(2),['Figure2' fileparam '.fig']);
        saveas(figure(3),['Figure3' fileparam '.fig']);
    %end
%end
dlmwrite([fileparam 'KIDparam.dat'], cell2mat(struct2cell(KIDparam)), 'precision', '%.9g', 'delimiter', '\t', 'newline', 'pc');
dlmwrite([fileparam 'arrayparam.dat'], cell2mat(struct2cell(arrayparam)), 'precision', '%.9g', 'delimiter', '\t', 'newline', 'pc');

end

function dif=numdif(raw,span)
%raw is colx,coly. returns num differential in dif (x, dy/dx)
%span is usuall span parameter in smooth function

y=raw(:,2);
dx=raw(2,1)-raw(1,1);
dif=zeros(length(y),1);
dif(:,1)=raw(:,1);
for i=2:length(y)-1
    dif(i,2)=(y(i+1)-y(i-1))/(2*dx);%numerical differentiating of y
end
dif(1,2)=dif(2,2);
dif(length(y),2)=dif(length(y)-1,2);
dif(:,2)=smooth(dif(:,2),span);
%dif(:,2)= sgolayfilt(dif(:,2),3,span); % 3rd order
%dif(:,2)= dif(:,2); % 3rd order
end