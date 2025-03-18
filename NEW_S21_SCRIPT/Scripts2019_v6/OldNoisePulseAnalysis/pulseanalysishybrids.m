function [dRdthetaresult]=pulseanalysishybrids
clear all;
close all;
% pulseanalysis.m: Loads pulse data, S21 card data. Fits s21 to get Q factor and hence
% internal power. Throughs away beginning of pulse (given by ring) to remove 
% excessive ringing. 
% Calculates standard deviation (s) of pulse phase data, uses that to
% calculate data range to fit: from nostd*s to sizestd*s in phase, uses
% this range to fit dRdtheta
% Also fits lifetime, from the same start point as for dRdtheta, but to
% nlife*lifetime for stop point.

% important variables for fit:
% fitty =1 , fit compared to standard deviation
% seedlife - seed lifetime
% minsizestd : min amount of data to fit
% sizestd: max amount of data to fit
% norings: use Q factor to calculate ring time, don't fit this amount of
% data, not robust, so use 
% ring: start at this fraction of data - short t constant this might cause
% problems
% lifeiterations: multiple fits to converge data over:
% nlife, time of data to fit to
% windowsize: moving average filter over data, set to 1 if not wanted

% program history:
% 13/3/10 cleaned up and made stand alone for new hybrid software, removed
% qp calculations
% 25/11/08 modified bounding problem, changed fit to exp((x-t0)/tau),
% offset to data so pulse start is same for fit, improves seeding of
% amp/offset to fit.
%
% Modified from ampresp7.m 
%Then fits pulses to extract lifetime in phase and amplitude
%%%% version history:
% 2/9/08 renamed to pulseanalysis
%1/9/08 amprespv7 cleaned up program, removed binning in phase, reorganised
% Also: reworked std deviation fit range, works now
% added exp fit to pulse data for lifetime, using similar fit range
% added output file at each power JB


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
path = [cd '/..']; %root path where data is, one higher than the scripts
lifetimeDir='/Session2/Lifetime_120mK';
% 
% path='/Volumes/KID/KIDonSun/experiments/ADR general/BoxinBoxsetup/Hybrids/LT021 Chip7/Selected KIDs';%Directory of M file
% lifetimeDir='Lifetime';%can be only 'Lifetime' 
makeplots=1;%set to 0 if more than 100 files to be analyzed; no plots wil be saved, 1 otherwise

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fitting parameters
% select (pulse) powers to fit
powerstopoffset=0% throw away last couple of powers
powerstartoffset=0 % throw away start powers
% fit
seedlife=1e-3 %% seed lifetime s

% fitting parameters, range in phase cf to peak to peak deviation of noise
% Small values give longer lifetime, but can be more sensitive to other
% reponses (thermal/2LS) and so give inconsistent lifetime and dRdTheta vs
% pulse power
% data fitted is a pulse that starts with an amplitude sizstd and is fitted over a time range
% by given # of tau (see below) fit is also stopped
%
sizestd=10 %signal size where the fit starts (i.e. fit pulse start) in st.dev. of noise, NORMAL= 20
minsizestd=5 %signal size to reject pulses with amplitude smaller than this one in st.dev. of noise, NORMAL=5

% set offset to remove at start of pulse, remove ringing etc.
norings=5; %5 offset at  beginning of data by this number of KID response times
ring=0.12% fraction of data to start. Pulse at 0.1, removes ringing if higher, 0.11, not critical

% fit done iteratively, fitting over a number of lifetime
lifeiterations=3; % number of iterations to fit lifetime. NORMAL=3 
nlife=5 %number of timeconstants to fit over, NORMAL=5
% smooth data
windowsize=25 %filter data if nequal to 1 with this moving average window, NORMAL=25

% fit to maximium phase (1 rad normal)
maxallowphase=1;  % maximum allowed phase to fit NORMAL=1 

figures=1; %just one figure

format('long','g');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%Finding files %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lifetimepath=[path filesep lifetimeDir filesep];%contains files
Exportpath=lifetimepath;
lifetimefiles=dir([lifetimepath '*mK.dat']);

disp(['number of files to be processed: ' num2str(length(lifetimefiles))]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
seedlife1=seedlife;

for n=1:length(lifetimefiles)%every file
    % Read lifetime files in
    
    disp(['processing ' lifetimefiles(n).name]);
    %reading file. T dependence relies on sorting of the fils by dir command, so that all files of each kid are done before another KID is done%
    fid=fopen([lifetimepath lifetimefiles(n).name],'r');
    KIDID=cell2mat(textscan(lifetimefiles(n).name,'%*3s %f %*s'));%gets KIDID
    readpower(n)=cell2mat(textscan(fid, '%*s %*s %*4s%f ',1, 'headerLines' ,1));
    T=cell2mat(textscan(fid, '%*s %*s %*2s%f ',1, 'headerLines' ,5)); %% get T
    qrt(n,1)=T;
    C=cell2mat(textscan(fid,'',1, 'headerLines' ,9));
    nopowers=size(C,2)/4; % calculate no. of columns
    fseek(fid,0,'bof'); %reset pointer
    %format read strings
    %readstringy=['%*s ' '%*s ' '%*s' '%*s']; % remove pulse powers at KID, old
    readstringy=['%*s ' '%*s ' '%*s' '%*s' '%*s']; % remove "LED time: 100, LED currents:,"
    readstringy2=[];
    for i=1:nopowers
        readstringy=[readstringy '%f '];  %read powers 
        readstringy=[readstringy '%*s ']; % ignore , between powers
        readstringy2=[readstringy2 '%f %f %f %f '];  % read string for file, 4 columns per power
    end
    nopowers=nopowers-powerstopoffset; % remove last couple of powers
    pulsepowers=cell2mat(textscan(fid,readstringy,1, 'headerLines' ,8));
    pulsedata=cell2mat(textscan(fid,readstringy2,'headerlines',1));
    fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %%% Load S21 file in, get seed Q factors
    S21filename=[lifetimepath lifetimefiles(n).name(1:length(lifetimefiles(n).name)-4) 'S21.dat'];
    fid=fopen(S21filename,'r');
    %%% load data, get seed values
    temp=textscan(fid,'%*s%*s%*s%*s%*1c%f',1,'headerlines',2);%getting F0 from header
    qrt(n,9)=temp{1};%Fo
    temp=textscan(fid,'%*2c%f%*5c%f%*5c%f%*7c%f',1,'headerlines',0);%getting Q's from header
    qrt(n,2)=temp{1};%Q
    qrt(n,3)=temp{2};%Qc
    qrt(n,4)=temp{3};%Qi
    qrt(n,10)=10^(temp{4}/20);%S21 magnitude for fit
    qrt(n,13)=-readpower(n);
    textscan(fid,'%*[^\n]',2);%scanning  empty lines to skip rest of header
    fseek(fid,0,'bof');
    S21data{n}=cell2mat(textscan(fid,'%f%f%f','headerlines',8)); %error in programs?
    fclose(fid);
 
    %%% fit S21, data not properly normalised/scaled, so lost absolute,

    %scale to remove arbitary gain factor,
    fuckybla=sqrt(mean(S21data{n}(:,2).^2+S21data{n}(:,3).^2));
    % shift real back to 1, scale from 1 to 0, Q factor independent of dip
    % depth so correct
    data(:,1)=S21data{n}(:,1); % f
    data(:,2)=(((1+0.5)/2+S21data{n}(:,2)/fuckybla/2/2).^2+(S21data{n}(:,3)/fuckybla/2/2).^2).^(1/2); % magn S21 for fit
    %data(:,2)=S21data{n}(:,2);
    %%% Fit S21
    [a(n),b(n),cc(n),l(n),stretch(n),tempresult{n}]=FitS21_1(data(:,:),qrt(n,9),qrt(n,2),0.5);   %% fres, Qf, S21min fitted, assymetric
    qrt(n,9)=a(n); %fres obtained from skewed Lorenztian, S21min
    qrt(n,2)=b(n);  %Q factor Obtained from Fit of parabolic
    qrt(n,4)=qrt(n,2)./qrt(n,10);%Get Qi=Q/S21min, S21min from data header
    qrt(n,3)=(qrt(n,4).*qrt(n,2))/(qrt(n,4)-qrt(n,2));%Get Qc=QiQ/(Qi-Q)
    qrt(n,18)=10*log10((2/pi)*(10.^(qrt(n,13)./10).*(qrt(n,2).^2./qrt(n,3))));%internal power
    pinternal(n)=qrt(n,18);
    ringtime=qrt(n,2)/qrt(n,9)/1e9; % KID ringtime= Q/f0
    %assymetry parameters
    %qrt(n,18)=l(n); 
    qrt(n,19)=stretch(n);
     %assigns qrt col.12 qrt 1,2,3,4,9,10,11,12,17 known%
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Main data proccesing
    %find the correct part to fit; this will be the ringing on of the KIDs
    
    pcount=0;
    nopoints=length(pulsedata(:,1));
    pulseindex(1:nopowers)=round(nopoints*ring); %seed to start above pulse, plus ringing
    phaselife=[]; % empty lifetime array
    amplife=[];
    for p=(1+powerstartoffset):nopowers
%%%%%%% Load data into new cells for future reference, take only data after
%%%%%%% cut for pulse beginning and ringing
        offset=pulsedata(round(0.1*nopoints),4*p-3); % pulse time
        time{p}=pulsedata(pulseindex(p):nopoints,4*p-3);
        Re{p}=pulsedata(pulseindex(p):nopoints,4*p-2);
        Im{p}=pulsedata(pulseindex(p):nopoints,4*p-1);
        phase{p}=pulsedata(pulseindex(p):nopoints,4*p);
        newlength=length(phase{p});
        
        time{p}=time{p}; % t=0 is the pulse point
        timingbla=time{p}(2)-time{p}(1);
      
        ringindex=ceil(norings*ringtime/timingbla); %round up no. of points in ring time, throw away this data

        phase{p}=phase{p}(ringindex:newlength)-mean(phase{p}((round(newlength*0.95)):newlength)); % take last points for removal of offset
        Rnorm{p}=(pulsedata(pulseindex(p):nopoints,4*p-1).^2+pulsedata(pulseindex(p):nopoints,4*p-2).^2).^0.5;
        Rnorm{p}=Rnorm{p}(ringindex:newlength)/mean(Rnorm{p}(round(newlength*0.95):newlength));
        newlength=length(phase{p});
        
%%%%%%% standard deviation of noise, calculate before any data sorting
        stdphase=std(phase{p}(round(newlength*8/10):newlength)); % noise peak to peak deviation
%%%%%%% get maximum of data
        maxphase(p)=max(phase{p}(:));

        %%%%%%% catch saturated response

        if maxphase(p)>2.5
            phasemaximumindex=find(phase{p}(:)>2.5,1,'last');
            maxphase(p)=phase{p}(phasemaximumindex);
%%%%%%% remove ringing
            phase{p}=phase{p}(phasemaximumindex:newlength);
            Rnorm{p}=Rnorm{p}(phasemaximumindex:newlength);
            time{p}=time{p}(phasemaximumindex:newlength);
            phasemaximumindex=1;
            newlength=length(phase{p});

            disp(['KID ringing KID' num2str(p) ]);
        else
            phasemaximumindex=find(phase{p}(:)==max(phase{p}(:)),1,'first');
        end

        
%%%%%%% Sort data 
     %   [fuckedd,fuck]=sort(phase{p},'descend');
    %    phasesort{p}=phase{p}(fuck);

        %if fitty==1 % fit cf std of data, or fit cf to phase conditions
            startphase=sizestd*stdphase ;                       % set phase conditions rel to std
            minphasecondition=minsizestd*stdphase;
       % end 
        if startphase>maxallowphase
            startphase=maxallowphase; %% limit range to fit
        end


%%%%%%% find fitting range points in data, work on sorted data
        if isempty(find(phase{p}(1:length(phase{p}))>startphase,1,'last') )
            startindex(p)=1;  
        else
            startindex(p)=find(phase{p}(1:length(phase{p}))>startphase,1,'last') ;
        end
        maxphase(p)=phase{p}(startindex(p));
      %  if relstopphase>maxphase(p);relstopphase=maxphase;end 
%%%%% stop index no longer used except for data rejection
        if (maxphase(p)>minphasecondition)% large enough response to fit
% stop point in fitting given by iterative finding of renge cf to tau
            stopindex(p)=newlength;
        else %needed for pulse without any response
            %stopindex(p)=stoptindex(p);
            stopindex(p)=startindex(p);
        end
        %don't fit if not enough data
        if ((phase{p}(startindex(p))-phase{p}(stopindex(p))))<(minsizestd*stdphase)
            startindex(p)=newlength ;
        end
        thisPhasdrdtheta=0;
        if (stopindex(p)-startindex(p)>(2*windowsize+12)) %enough point (3) and low enough max phase
                pcount=pcount+1;
                thisPhasdrdtheta=1;
% filter data with moving average filter
                if windowsize >1
                    phase{p}=filter(ones(1,windowsize)/windowsize,1,phase{p});
                    Rnorm{p}=filter(ones(1,windowsize)/windowsize,1,Rnorm{p});
                end

%%%%% fit lifetime, similar ranges

                phaselife(p)=seedlife1;
                amplife(p)=seedlife1;
                for fuckybla=1:lifeiterations
%%%%% for lifetime, set stop for fit to end, or if multiple fits, to the
%%%%% lifetime
                    if lifeiterations>1
                        stopindex(p)=startindex(p)+round(nlife*phaselife(p)/timingbla);
                        if stopindex(p)>newlength
                            stopindex(p)=newlength;
                        end
                    else
                        stopindex(p)=newlength;
                    end
                    if (stopindex(p)-2*windowsize)<12 % catch broken fit
                        stopindex(p)=newlength
                    end
                    if (stopindex(p)-windowsize-2-startindex(p)-windowsize-2)<10
                        startindex(p)=stopindex(p)-2*windowsize-8 %shouldn't be happening!
                        stopindex
                    end
                    poorange2=(startindex(p)+windowsize+2):(stopindex(p)-windowsize-2);
                    if size(poorange2,2)<5
                        poorange2
                    end
                    seedphase=phase{p}(startindex(p)+windowsize+2)/exp(-(time{p}(startindex(p)+windowsize+2)-offset)/phaselife(p));
%                    s = fitoptions('Method','NonlinearLeastSquares',...
%                       'Startpoint',[1 seedlife1 0.001]);
                    s = fitoptions('Method','NonlinearLeastSquares',...
                       'Startpoint',[seedphase seedlife1 0.001],...
                    'Upper',[100,0.1,1],...
                    'Lower',[-100,1e-6,-1]);
      
                    f = fittype('a*exp(-(x-xoffset)/b)+c','problem','xoffset','options',s);
                    [lifefit,bler,output] = fit(time{p}(poorange2),phase{p}(poorange2),f,'problem',offset);
                    phaselife(p)=lifefit.b;
                    seedamp=-(1-Rnorm{p}(startindex(p)+windowsize+2))/exp(-(time{p}(startindex(p)+windowsize+2)-offset)/amplife(p));

                   s = fitoptions('Method','NonlinearLeastSquares',...
                       'Startpoint',[seedamp seedlife1 1],...
                       'Upper',[10,0.1,10],...
                       'Lower',[-100,1e-6,-10]);

                    f2=fittype('a*exp(-(x-xoffset)/b)+c','problem','xoffset','options',s);
                    ampfit= fit(time{p}(poorange2),Rnorm{p}(poorange2),f2,'problem',offset);
                    amplife(p)=ampfit.b;
                end
                seedlife1=phaselife(p)
% fit amplitude responsivity, same range as for lifetime 
                s = fitoptions('Method','NonlinearLeastSquares',...
                   'Upper',[0.1,2],...
                   'Lower',[-2,0],...
                   'Startpoint',[-0.25 1]);
                f = fittype('a*x+b','options',s);
                %poorange=(startindex(p)+windowsize+2):(stopindex(p)-windowsize-2);
                fitresult = fit(phase{p}(poorange2),Rnorm{p}(poorange2),f);

                dRdtheta{n}(pcount)=fitresult.a;
                Pdone{n}(pcount)=pulsepowers(p);
%%%% plot
            if makeplots   
                if (figures==1)&(pcount==1)
                    clf;
                elseif figures~=1
                    figure(n); 
                else
                    figure(gcf)
                end
                subplot(2,2,1); plot(time{p}(windowsize:(newlength-windowsize)),phase{p}(windowsize:(newlength-windowsize)),'.b','MarkerSize',1);ylabel('phase [rad]');xlabel('time [msec.]');title(['KID ' num2str(KIDID) ' @ ' num2str(T) ' @ ' num2str(readpower(n))]);hold on;
                plot(time{p}(poorange2),lifefit.a*exp(-(time{p}(poorange2)-offset)/(lifefit.b))+lifefit.c,'r-');%fit
                subplot(2,2,2); plot(time{p}(windowsize:(newlength-windowsize)),Rnorm{p}(windowsize:(newlength-windowsize)),'.b','MarkerSize',1);ylabel('R/<R> []');xlabel('time [msec.]');title(['P ' num2str(pulsepowers)]);hold on;
                plot(time{p}(poorange2),ampfit.a*exp(-(time{p}(poorange2)-offset)/(ampfit.b))+ampfit.c,'r-');%fit
                subplot(2,2,4); plot(Re{p},Im{p});axis tight;hold on;plot(S21data{n}(:,2),S21data{n}(:,3),'b-');
                subplot(2,2,3); 
                plot(phase{p}(windowsize:(newlength-windowsize)),Rnorm{p}(windowsize:(newlength-windowsize)),'.g','MarkerSize',1);xlabel('phase [rad]');ylabel('R/<R> []');hold on;%data
                plot(phase{p}(poorange2),(fitresult.b+fitresult.a*phase{p}(poorange2)),'--','linewidth',2);%fit

            end
        else
            if pcount==0
                Pdone{n}(1)=0;dRdtheta{n}(1)=0;
            end
        end
        %output arrays
        arraydata(p,1)=p;
        
        if thisPhasdrdtheta
            arraydata(p,2)=phaselife(p);
            arraydata(p,3)=amplife(p);
            arraydata(p,4)=dRdtheta{n}(pcount);
        else
            arraydata(p,2)=0;
            arraydata(p,3)=0;
            arraydata(p,4)=0;
        end
            
    end
    
    if isempty(phaselife)~=1
    path=[Exportpath 'dRdtheta_KID' num2str(KIDID) '.csv'];
    keeswrite('Power, phaselifetime, amlifetime, drdtheta',path);
    dlmwrite(path,arraydata ,'-append','newline', 'pc', 'precision', '%.6g','delimiter', ',');
    Figfile=[Exportpath 'K' num2str(KIDID) '_' num2str(readpower(n)) 'dBm' '_dRdtheta_' num2str(T,'%.3f') 'K.fig'];
    end
   % [Pdone{n}',dRdtheta{n}']
   
    %ignore zero values!
    dRdthetaresult(n,1)=KIDID;
    dRdthetaresult(n,2)=mean(nonzeros(dRdtheta{n}));
    dRdthetaresult(n,3)=std(nonzeros(dRdtheta{n}));
    dRdthetaresult(n,4)=T;
    dRdthetaresult(n,5)=readpower(n);
    dRdthetaresult(n,6)=pinternal(n);
    dRdthetaresult(n,7)=qrt(n,2); %Q
    dRdthetaresult(n,8)=qrt(n,3); %Qc
    dRdthetaresult(n,9)=qrt(n,4); %Qi
    dRdthetaresult(n,10)=qrt(n,9); %f0
    dRdthetaresult(n,11)=mean(nonzeros(phaselife(:))); %lifetime
    dRdthetaresult(n,12)=std(nonzeros(phaselife(:))); %lifetime
    dRdthetaresult(n,13)=mean(nonzeros(amplife(:))); %lifetime
    dRdthetaresult(n,14)=std(nonzeros(amplife(:))); %lifetime
%    dRdthetaresult(n,15)=qrt(n,12); %noqps
%    dRdthetaresult(n,16)=dRdthetaresult(n,15)*delta*1.602e-19/(dRdthetaresult(n,11)*0.57); %eff sky power    
    
    saveas(gcf,Figfile,'fig');
end
format short g;
disp(num2str(dRdthetaresult));
disp(num2str([dRdthetaresult(:,1) dRdthetaresult(:,11) dRdthetaresult(:,13)]));
path=[Exportpath 'dRdtheta.csv'];
%keeswrite('KID,  dRdtheta, std, T[K], readpower, pinternal, Q, Qc, Qi, f0, phase life, std, amp life, std, Nqp, Psky',path);
keeswrite('KID,  dRdtheta, std, T[K], readpower, pinternal, Q, Qc, Qi, f0, phase life, std, amp life, std',path);
dlmwrite(path,dRdthetaresult ,'-append','newline', 'pc', 'precision', '%.6g','delimiter', ',');
end

function  [Fresres,Qfres,S21minres,l,stretch,result] = FitS21_1(data,Fres,Qf,S21min)
%version 13/3/10 copied here
%modified 24/6/09 - modified to fit Q factor over 2 bandwidths in log space
%fit f0/smin from narrower range

%perfoms skewd Lorenztain fit for Q, parabolic fit for Fres and S21min
%Data should contain F[GHz], S21 magn., (phase optional). Fres is res freq
%in GHz, Qf is initial Q factor and S21min is original S21min. These values
%are used as guess values for a fit to Bens eqn. 
%07-03-07: modified for bad ini guesses bandwidth
format('long','e');

no=2; %width in resonator bandwidths to fit over
datasize=size(data(:,2),1);
%limit width of fit
bandwidth=Fres/Qf; %resonator bandwidth
minfitindex=find(data(:,1)>(Fres-no*bandwidth/2),1);
maxfitindex=find(data(:,1)>(Fres+no*bandwidth/2),1);
if minfitindex+4<maxfitindex
    x=data(minfitindex:maxfitindex,1);
    %y=data(minfitindex:maxfitindex,2).^2;
    y=20*log10(data(minfitindex:maxfitindex,2));
else
    x=data(:,1);%y=data(:,2).^2;
    y=20*log10(data(:,2));
end
% estimate linear background
l=(data(1,2)^2-data(datasize,2)^2)/(data(1,1)-data(datasize,1)); %%linear background drift

%skewterm

% catch error in finding guess parameter
low3dbpoint=find(data(:,1)>(Fres-bandwidth/2),1);
if (isempty(low3dbpoint)==1)
    low3dbpoint=1
end
high3dbpoint=find(data(:,1)>(Fres+bandwidth/2),1);
if (isempty(high3dbpoint)==1)
    high3dbpoint=size(data(:,1),1)
end

stretch=(data(low3dbpoint,2)^2-data(high3dbpoint,2)^2)/(bandwidth); % estimate skew parameter

if Qf<100; Qf=100; end; if Qf>2e6; Qf=2e6; end

s = fitoptions('Method','NonlinearLeastSquares', 'Startpoint',[Fres Qf S21min l stretch 1],'Lower',[Fres*0.9 Qf/2 S21min/4 -25*abs(l) -25*abs(stretch) 0],...
    'Upper',[Fres*1.1 Qf*2 S21min*2 abs(l)*25 +25*abs(stretch) 10^(y(1)/10)*5],'MaxFunEvals',1000);

%Fit Lorentzian like in power space with linear background and stretch
ftype = fittype('10*log10(abs(y0*(1-(1-Smin^2)*(1+stretch*(x-Fr)/Fr)/(1+(2*Q*(x-Fr)/Fr)^2)+l*(x-Fr)/Fr)))','options', s);  
[result,gof,output]=fit(x,y,ftype);
%figure(1);clf;hold on;plot(x,y,'+r');plot(result,'r');
Qfres=result.Q;
l=result.l;
stretch=result.stretch;
%Ffit=result.Fr;
%Smin=result.Smin;

%%new: a simple lorentzian in log space fit to find resonance
%%frequency/S21min. Is done in power space
so=0.4;

minfitindex=find(data(:,1)>(Fres-so*bandwidth/2),1);
maxfitindex=find(data(:,1)>(Fres+so*bandwidth/2),1);
if maxfitindex-minfitindex>10
    x=(data(minfitindex:maxfitindex,1));
    y=20*log10(data(minfitindex:maxfitindex,2));
else
    x=(data(:,1));
    y=20*log10(data(:,2));
    disp('bandwidthguess for fit far off, whole S21 curve used for Fo and S21');
end

s = fitoptions('Method','NonlinearLeastSquares', 'Startpoint',[Fres Qfres S21min], 'MaxFunEvals',300);

ftype = fittype('10*log10(1-(1-Smin^2)/(1+(2*Q*(x-Fr)/Fr)^2))','options', s);

[resultFres]=fit(x,y,ftype);

Fresres=resultFres.Fr;

S21minres=resultFres.Smin;


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