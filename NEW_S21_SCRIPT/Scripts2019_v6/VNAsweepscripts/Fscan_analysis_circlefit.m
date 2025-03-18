function [arrayparam,KIDparam,procdata]=Fscan_analysis_circlefit(data,fileparam,ReadP,cryostat)

% 6/5/2013, analysis f-scan using S21phase to find KIDs, plus use KID phase
% to calculate Qs
% Power is the readout power of VNA
if nargin==2
    ReadP=0;
    cryostat='He10';
end
%%

method='S21logmagn'%method='S21phase' NEVER USE IT
findmethod='stdev';%'diff2treshold'

smoothwindow=10 % for smoothing dphi/df

dthetawindow=4
nbstdevthreshold=4 % in KID finding wrt stdeviation, number above stdev to be flagged KID
minlength=2 % number of points required per dip
plotty=0;
nptstofit=15;
diff2treshold=1e4; % if diff2threshold used
filteronQC=1; %can throw away bad Qc dips
nbstev=4%2.5 % nb of std deviation in log(Qc) to reject KIDs for noisy data

MUXv4dir='/Users/stepheny/Documents/Localwork/matlabprogram/MUX v4' % directory for scripts
cleanpath=1 % remove from path

if isempty(strfind(path,MUXv4dir))
    addpath(MUXv4dir)
    addpath([MUXv4dir '/general analysis'])%,[MUXv4dir '/utrecht'],[MUXv4dir '/utrecht/filterfiles'],[MUXv4dir '/utrecht/subroutines'])    
end

dir='';%directory incl \ at end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
comdata=complex(data(:,2),data(:,3));

%% 
% remove phase delay
rangetofit=1:length(data);
phasedelay=polyfit(data(rangetofit,1),unwrap(angle(smooth(comdata(rangetofit),1000))),9)

maxrawdata=max(abs(comdata));
cor_comdata=exp(-1i*polyval(phasedelay,data(rangetofit,1))).*comdata/maxrawdata;

%% find KID dips

if strcmp(method,'S21phase')
    ddata_df=numdif([data(:,1),(angle(cor_comdata))],smoothwindow);
%clf
%plot(ddata_df(:,1),ddata_df(:,2));
elseif strcmp(method,'S21logmagn') 
    ddata_df=numdif([data(:,1),20*log10(abs(cor_comdata))],smoothwindow);
    ddata_df=numdif(ddata_df,smoothwindow); % double derivative
end

if strcmp(findmethod,'stdev')
    stdddata_df=std(ddata_df(20:1000,2));
    KIDregion=find(ddata_df(smoothwindow:end-smoothwindow,2)>stdddata_df*nbstdevthreshold)+smoothwindow-1;
elseif strcmp(findmethod,'diff2treshold')
    KIDregion=find(ddata_df(:,2)>diff2treshold);
end
%%
clf
plot(ddata_df(:,1),ddata_df(:,2));
hold on
plot(ddata_df(KIDregion,1),ddata_df(KIDregion,2),'+r')
%%
% Determine KIDs from breaks in indices

nonsequentialregion=zeros(length(KIDregion),2);
nonsequentialregion(:,1)=(1:length(KIDregion))';
for i=(2:length(KIDregion))
    nonsequentialregion(i,2)=(KIDregion(i)-KIDregion(i-1));%numerical differentiating of y
end

discontinuity=find(nonsequentialregion(:,2)>1);
KIDind=[];
KIDind(1,:)=[1 discontinuity(1)-1]';
KIDind(2:(length(discontinuity)),:)=[discontinuity((2:(length(discontinuity)))-1) discontinuity((2:(length(discontinuity))))-1];
KIDind(length(discontinuity)+1,:)=[discontinuity(end) length(KIDregion)];
KIDind=KIDregion(KIDind);



KIDind((KIDind(:,2)-KIDind(:,1)<minlength),:)=[];
KIDsfound=size(KIDind,1);

%
%clf
%plot(data(:,1),20*log10(abs(cor_comdata))); hold on
%colors = CFmymap(KIDsfound);
%%
padding=nptstofit*5;% extra around KIDs
noKIDind=[1:(KIDind(1,1)-padding)];
for fuck =1:KIDsfound
    indextoplot=KIDind(fuck,1):KIDind(fuck,2);
    if fuck~=KIDsfound
    noKIDind=[noKIDind (KIDind(fuck,2)+padding):(KIDind(fuck+1,1)-padding)];
    end
   % plot(data(indextoplot,1),20*log10(abs(cor_comdata(indextoplot))),'r')%,'Color',colors(fuck,:));hold on
    [~,minind(fuck)]=min(abs(cor_comdata(indextoplot)));
    minind(fuck)=indextoplot(minind(fuck));
   % plot(data(minind(fuck),1),20*log10(abs(cor_comdata(minind(fuck)))),'og')%,'Color',colors(fuck,:));hold on
end

noKIDind=[noKIDind (KIDind(fuck,2)+padding):length(ddata_df)];
noKIDind=makeunique(noKIDind);
noKIDind=sort(noKIDind);
%title(['Number of KIDs found ' num2str(KIDsfound)])
%
%
% ripple correction
smp=20
rippletemp = smooth(abs(cor_comdata(noKIDind)),smp); %data smoothed
ripple=interp1(data(noKIDind,1),rippletemp,data(:,1),'linear','extrap');%smoothed data interpolated (to bridge gaps were KIDs ar removed)%
cor_comdata  = cor_comdata  ./ ripple;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% circle fit over range >5stdev of noise

for bob=1:KIDsfound
    
    fitrange=(minind(bob)-nptstofit):(minind(bob)+nptstofit);%KIDind(bob,1):KIDind(bob,2);
    fitrange(fitrange<1)=[];
    fitrange(fitrange>length(cor_comdata))=[];
    
    
    [xfit,yfit,Rfit,a] = circfit(real(cor_comdata(fitrange)),imag(cor_comdata(fitrange)),1);
    procdata(bob).cor_comdata=cor_comdata(fitrange);
    reelshift=sqrt(xfit^2+yfit^2);
    rotation=atan2(yfit,xfit);
    if plotty==1
        figure(1001)
        plot(cor_comdata(KIDind(bob,1):KIDind(bob,2))); hold on
        
        plot(cor_comdata(fitrange),'+k') 
        rectangle('position',[xfit-Rfit,yfit-Rfit,Rfit*2,Rfit*2],...
        'curvature',[1,1],'linestyle','-','edgecolor','r');%plots circle 8-)
        axis square
    end
    CALparam(bob).KIDcalphase=rotation;
    CALparam(bob).KIDcalreal=reelshift;
    CALparam(bob).KIDcalR=Rfit;
    
    procdata(bob).KIDcalS21data=KIDcalNEW([procdata(bob).cor_comdata],1,CALparam(bob),1);
    procdata(bob).f=data(fitrange,1);
    KIDparam((bob)).fres=data(minind(bob),1);
%    usedind=minind(bob)-KIDind(bob,1)+1;
    
    KIDparam(bob).S21min=min(abs(procdata(bob).cor_comdata));
 %   KIDparam(bob).S21used=abs((procdata(bob).cor_comdata(usedind)));
    

    KIDparam(bob).diffed=numdif([procdata(bob).f*1e6 angle(exp(-1i*pi)*procdata(bob).KIDcalS21data)],dthetawindow);
    KIDparam(bob).dthetadfmin=min(KIDparam(bob).diffed(:,2));
    KIDparam(bob).Q=- KIDparam(bob).dthetadfmin*KIDparam(bob).fres*1e6/4; %get Q from dtheta/df max
%    KIDparam(bob).dthetadf_used=(KIDparam(bob).diffed(usedind,2));
    KIDparam(bob).Qi=KIDparam(bob).Q/KIDparam(bob).S21min;
    KIDparam(bob).Qc=KIDparam(bob).Qi*KIDparam(bob).Q/(KIDparam(bob).Qi-KIDparam(bob).Q);
end
%
% reject non KIDs
if filteronQC==1
pQc=polyfit(1:KIDsfound,log10([KIDparam.Qc]),1);
badlist=(abs(log10(([KIDparam.Qc]))-(polyval(pQc,1:KIDsfound)))>std(log10([KIDparam.Qc]))*nbstev);
disp(num2str(find(badlist)))
disp(num2str([KIDparam(find(badlist)).Qc]))

pQc=polyfit(find(~badlist),log10([KIDparam(~badlist).Qc]),1);
badlist=(abs(log10(([KIDparam.Qc]))-(polyval(pQc,1:KIDsfound)))>std(log10([KIDparam(~badlist).Qc]))*nbstev);


KIDsfound=sum(~badlist);
KIDparam(badlist)=[];
procdata(badlist)=[];
KIDind(badlist,:)=[];
minind(badlist)=[];
end
%

%% Calculate Readout Power and internal power
if strcmp(cryostat,'He10')
   load('kingcryo_thrucals.mat') 
else %Unknown cryo
    f=1000:100:10000;
    Cryoin_to_KID=ones(size(f));
end
% calculate 
for bob=1:KIDsfound
    if bob>1
        KIDparam(bob).df=KIDparam(bob).fres-KIDparam(bob-1).fres;
    else
        KIDparam(bob).df=0;
    end
    KIDparam(bob).ReadP=interp1(f,Cryoin_to_KID,KIDparam(bob).fres)+ReadP;
    KIDparam(bob).Pint=10*log10(10.^(KIDparam(bob).ReadP/10)*(2/pi)*(KIDparam(bob).Q^2/KIDparam(bob).Qc));
    KIDparam(bob).BW=KIDparam(bob).fres/KIDparam(bob).Q;
    KIDparam(bob).dfBW=KIDparam(bob).df/KIDparam(bob).BW;   
end


%% 
figure(1011)
subplot(2,2,1)
semilogy([KIDparam.Q],'b.')
hold on
semilogy([KIDparam.Qc],'r.')
semilogy([KIDparam.Qi],'k.')
legend('Q','Qc','Qi')
xlabel('KID #')
ylabel('Q factor')
ylim([1e4 2e6])
title(['KIDs found ' num2str(KIDsfound) ', Readout P ' num2str(mean([KIDparam.ReadP])) 'dBm'])
subplot(2,2,2)
hist([KIDparam.df]*1e3,0:0.5:10)
xlim([0 10])
xlabel('\Delta f (MHz)')

title(['\Delta f mean ' num2str(mean([KIDparam.df]*1e3),'%.1f') ' std ' num2str(std([KIDparam.df]*1e3),'%.1f')])
subplot(2,2,3)
hist([KIDparam.Qc]/1e3,(10:10:100))
title(['Qc mean ' num2str(mean([KIDparam.Qc]/1e3),'%.1f') ' std ' num2str(std([KIDparam.Qc]/1e3),'%.1f') 'e3'])
xlabel('Qc/1e3')
xlim([0 100])
subplot(2,2,4)
hist(20*log10([KIDparam.S21min]))
xlabel('S21min (dB)')
title(['S21min mean ' num2str(mean(20*log10([KIDparam.S21min])),'%.1f') ' std ' num2str(std(20*log10([KIDparam.S21min])),'%.1f')])

%hist([KIDparam.dfBW],1:1:100)
%xlim([0 100])
figure(999)
clf
subplot(2,1,1)

plot(data(:,1),20*log10(abs(cor_comdata))); hold on
subplot(2,1,2)
plot(ddata_df(:,1),ddata_df(:,2)); hold on
colors = CFmymap(KIDsfound);
goodlist=1:KIDsfound;%find(~badlist);
for fuck =1:length(goodlist)
    subplot(2,1,1)
    indextoplot=KIDind(goodlist(fuck),1):KIDind(goodlist(fuck),2);
    plot(data(indextoplot,1),20*log10(abs(cor_comdata(indextoplot))),'r')%,'Color',colors(fuck,:));hold on
 %   [~,minind(fuck)]=min(abs(cor_comdata(indextoplot)));
%    minind(fuck)=indextoplot(minind(fuck));
    plot(data(minind(goodlist(fuck)),1),20*log10(abs(cor_comdata(minind(goodlist(fuck))))),'og')%,'Color',colors(fuck,:));hold on

    subplot(2,1,2)
    plot(ddata_df(indextoplot,1),ddata_df(indextoplot,2),'r')%'Color',colors(fuck,:)); hold on;
end
title(['Number of KIDs found ' num2str(KIDsfound)])

%hist([KIDparam(~badlist).Qc]);
%% statistics
arrayparam.meanQ=mean([KIDparam.Q]);
arrayparam.meanBW=mean([KIDparam.BW]);
arrayparam.stdQ=std([KIDparam.Q]);
arrayparam.meanQc=mean([KIDparam.Qc]);
arrayparam.stdQc=std([KIDparam.Qc]);
arrayparam.meanQi=mean([KIDparam.Qi]);
arrayparam.stdQi=std([KIDparam.Qi]);
arrayparam.meanS21dB=mean(20*log10([KIDparam.S21min]));
arrayparam.stdS21dB=std(20*log10([KIDparam.S21min]));

arrayparam.meandf=mean([KIDparam.df]);
arrayparam.stddf=std([KIDparam.df]);
arrayparam.df10MHz=length(find([KIDparam.df]<0.01))-1;%number dF below 5 MHz. -1 beacuse df first KID with itself=0.
arrayparam.df5MHz=length(find([KIDparam.df]<0.005))-1;%number dF below 5 MHz. -1 beacuse df first KID with itself=0.
arrayparam.df1MHz=length(find([KIDparam.df]<0.001))-1;%number dF below 1 MHz. -1 beacuse df first KID with itself=0.
arrayparam.df02MHz=length(find([KIDparam.df]<0.0002))-1;%number dF below 0.2 MHz. -1 beacuse df first KID with itself=0.
arrayparam.dfBW_1=length(find([KIDparam.dfBW]<1))-1;
arrayparam.dfBW=mean([KIDparam.dfBW]);
arrayparam.dfBWstd=std([KIDparam.dfBW]);

arrayparam.ReadPstart=[KIDparam(1).ReadP];
arrayparam.ReadPstop=[KIDparam(end).ReadP];
arrayparam.Pintmean=mean([KIDparam.Pint]);
arrayparam.Pintstdev=std([KIDparam.Pint]);
arrayparam.AMKIDnbpixels=length(find(([KIDparam.dfBW]>1)&(([KIDparam.Qc]<100e3)&((([KIDparam.fres]>4.26)&([KIDparam.fres]<6.2))|(([KIDparam.fres]>6.3)&([KIDparam.fres]<8.2))))));

% display
fprintf('#tonesfound = %.0f \n', KIDsfound);
fprintf('df: %.1f MHz pm %.1f \n',arrayparam.meandf*1000,arrayparam.stddf*1000)
fprintf('# BW separation: %.0f pm %.0f \n', arrayparam.dfBW, arrayparam.dfBWstd)
fprintf('Q/1e3: %.3f pm %.3f \n', arrayparam.meanQ/1e3, arrayparam.stdQ/1e3)
fprintf('BW: %.0f \n', arrayparam.meanBW*1e9)
fprintf('Qi/1e3: %.3f pm %.3f \n', arrayparam.meanQi/1e3, arrayparam.stdQi/1e3)
fprintf('Qc/1e3: %.3f pm %.3f \n', arrayparam.meanQc/1e3, arrayparam.stdQc/1e3)
fprintf('S21[dB]: %.0f pm %.0f \n', arrayparam.meanS21dB, arrayparam.stdS21dB)
fprintf('# tones with df < 1MHz: %.0f \n', arrayparam.df1MHz)
fprintf('# tones with df < 5MHz: %.0f \n', arrayparam.df5MHz)
fprintf('# tones with df < 10MHz: %.0f \n', arrayparam.df10MHz)
fprintf('# tones within 1BW %.1f \n', arrayparam.dfBW_1)


saveas(figure(999),[dir fileparam 'Figure1'  '.fig']);
saveas(figure(1011),[dir fileparam 'Figure2'  '.fig']);
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

function   [xc,yc,R,a] = circfit(x,y,seedR)
%
%   [xc yx R] = circfit(x,y)
%
%   fits a circle  in x,y plane in a more accurate
%   (less prone to ill condition )
%  procedure than circfit2 but using more memory
%  x,y are column vector where (x(i),y(i)) is a measured point
%
%  result is center point (yc,xc) and radius R
%  an optional output is the vector of coeficient a
% describing the circle's equation
%
%   x^2+y^2+a(1)*x+a(2)*y+a(3)=0
%
%  By:  Izhak bucher 25/oct /1991, 
   x=x(:); y=y(:);
   a=[x y seedR*ones(size(x))]\(-(x.^2+y.^2));%square bracjedts used before here
   xc = -.5*a(1);
   yc = -.5*a(2);
   R  =  sqrt((a(1)^2+a(2)^2)/4-a(3));
end