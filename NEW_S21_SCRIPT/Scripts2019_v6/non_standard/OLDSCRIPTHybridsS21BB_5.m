function [qrt,data]=HybridsS21BB_5
% 18/4/2012 modified to be compatible with blacbody_2.m 
% 10/2/2012 Added full large aperture calculation of power and NEP.
% Includes correction for beam pattern on aperture efficiency (cf to
% lambda^2), photon NEP corrected for wave+g-r noise 


% analyses BB scan measured with VNA. Either SPICA or 650GHz BB, set via
% setup string
% file KIDdata.dat has: 
% 17	4.2         14000       42.156
% ID    Fres[Ghz]   Qdesign     A in mm^2
% OUTPUT FILES
% 1: Tdep.csv (for each KID)
% 1KIDID	2T	3Q	4Qc	5Qi	6F0[]	7S21@res	8power	9Pinternal	10response	11dtheta/dN FIT	12dtheta/dN numerical	
% 13qpnumber	14 dx/dN	 15 V[um3]	 16Delta [J]	 17Re point	 18Im point	19drdnqp=ampresp	20dRdB	21grad 1/Qi vs nqp 22 Tc %
% 23 Fres des 24 Q des 25 Area des 26=thickness 
%
% 2: low T S21 for each KID
% Freq	S21dB	theta	Re Norm	 Im Norm
% 
% 3: all T dep parameters at the lowest temperature the program has used for all KIDs in 1 file
% =same data as in 1; Tdep.csv but than all kifs at base T in 1 file
% drdb= numerical variation of radius with 1/Qi,m normalised to lowest T (from SY)
%ADD PLOT OF THE CIRCLE USING RFIT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%PARAMETERS TO SET%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
path='C:\Documents and Settings\Jochem\My Documents\iFolder\KID\Technology_Experiment\DATA\H20_1 APEX 13_2_12';%Directory of M file
localDir='TempBB';  %TempBB
%%%%%%%%% metod%%%%%%%%%%
%method.filter='325 GHz AB'         %AB old BP filter for A-MKID
method.filter='350 GHz WT1275';     %Cardiff filter
%method.filter='325 GHz He7'        % He7 system filtering used in A-MKID L band He7 cooler
%
method.lensdiameter=2;              %1 or 2 %lens size
%
%method.aperture='He7'                      %He7 system
method.aperture='ADR small';                %small aperture in the ADR, 0.071rad is the angle to 1 mm radius hole, theta=atan(1mm/14mm), angle to normal, TP=2pi(1-costheta)
%method.aperture='ADR small Geometrical';   %small 2 mm aperture usinng geometrical calculation
%method.aperture='ADR full'                 %full aperture in the ADR
                                            %NB: 2 mm lenses only implemented for ADR small and ADR small Geometrical
%%%%%%%%%%%%%%%%%%%%%%
startoffset=0;
stopoffset=0;       %NB program does filter for too shallow dips
Pstart=5E-15;       %power below which there is no response measureable. 5e-15 for 350 Ghz small, 2.1e-13 350 Ghz full thruput, 1e-16 SPICA CHECK WITH IT!!!!
FITF0=1;            %Does not use fitted S21 data for F0, handy for very skewed KIDs. NB: But fit still done for Q factor
maxnumkids=100;     %maximum number of KIDs that can be analyzed (used for allocating variables)
smprm=3;            %2x this number+1 is window used for differentiating dx. use 2-4 here
dipmin=-0.5;        %dip depth (in dB, negaitve value!)below which data is thropwn away

pol=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
format('long','e');
clear data qrt;
S21path=[path '\S21\' localDir '\'];%Contains Raw S21 data
disp(S21path);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Finding files %%%%%%%%%%%%%%%%%%%%%%%%%%%%
files=dir([S21path 'KID*.dat']);%%%%Reading in all the raw S21 filenames in array of cells
m=1;c=zeros(1,maxnumkids);
for n=1:length(files)%for loop used to display KID numbers and not all KID files
    temp=cell2mat(textscan(files(n).name,'%*3s %f %*s'));
    if m==1
        c(m)=temp; 
        m=m+1;
    elseif n~=1 && c(m-1)~=temp 
        c(m)=temp; 
        m=m+1;
    end
end
clear temp;
c(m:end)=[];
numberofscans=length(files);
disp(['Kids availabe: ' num2str(c) ', total number of files ' num2str(numberofscans)]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Working each file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n=1:length(files)%all files = all KIDs
    clear data qrt;
    file=files(n).name
    info=textscan(file,'%*3s %f %*1s %f %s');
    ChipID=char(info{3}(1));
    ChipID=ChipID(4:cell2mat(strfind(info{3},'.dat'))-1);%automatically finds Chip ID in Filename
    place=regexpi(file,'KID','end');%after which number follows
    info=textscan(file,['%*' num2str(place) 's %f %*1s %f %*s']);
    kidnumber=info{1};
    kidpower=info{2};
    %%%%%%%%%%%%%%%%%%reading in and processing KID data%%%%%%%%%%%%%%%%%%%%%%%
    clear data qrt;
    filetoread=[S21path file];
    if n==1
        plaatje=1;
    else
        plaatje=0;
    end
    [data,qrt]=import_data(filetoread,startoffset,stopoffset,method,pol,plaatje);%column 1,5, 12 and 13 of qrt assigned
    [data,index,qrt,transformation,fitresult]=findparameters(data,qrt,FITF0,dipmin);%assigns qrt columns 2 3 4 9 10 11 18 and col. 4,5 of data{} index=minS21(all KIDs)
    [qrt]=findresponse(data,qrt,Pstart,smprm,transformation);%rest of qrt assigned
    qrt(:,15)=kidnumber;
    savethedata(S21path,qrt,kidnumber,kidpower,ChipID);
    plotdata(data,qrt,kidnumber,kidpower,S21path,fitresult,Pstart,index,transformation(:,3))%last is R of circle
end
end

function [data,qrt]=import_data(file,startoffset,stopoffset,method,pol,plt)
%imports datafile and kills the header, except for the power information.
%does invert order in case of decreasing T. Added JB 15-08-2008
%From the subheaders in the file the Temp is read.
% File format
%
%C:\KID Metingen\ADR\B6Ch7__7_11_06 16_02\S21\2D\KID41_73dBm_JB1.dat
%Power at KID:73dBm
%resonance Frequency in GHz :4.893858
%Q=40792.348421, Qc=56509.887020, Qi=146662.340671  
%
%
%Temperature in K:0.099953
%GHz	dB	Rad
%4.893115000	-6.722063065	2.685341760
%4.893118870	-6.721852779	2.685013925
%4.894663000	-6.540477276	2.794459459
%
%Temperature in K:0.105161
%GHz	dB	Rad
%4.893118740	-6.718620300	2.685038958
%4.893122571	-6.718233109	2.684799274
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

format('long','e');
fid=fopen(file);
temp=textscan(fid,'%*s%*s%*c%*c%*c%*c%f%*s',1,'headerlines',1);
power=-1*temp{1};
textscan(fid,'%*[^\r\n]',2);

n=1;% goes through all the temperatures
tempqrt=zeros(1000,1);tempdata=cell(1000,1);
while n<1000
    i=1;bla=[1,1,1];
    temp=textscan(fid,'%*s%*s%*c%*c%f%*s',1,'headerlines',0);
    if isempty(temp{1})
        break;
    end
    tempqrt(n,1)=temp{1};
    textscan(fid,'%*s%*s%*s',1,'headerlines',1);
    while ~isempty(bla)%bla is not empty
        bla=cell2mat(textscan(fid, '%f%f%f', 1));
        if ~isempty(bla)
            tempdata{n}(i,:)=bla;
        end
        i=i+1;
    end
    n=n+1;
end
clear temp;
tempqrt(n:end,:)=[];tempdata(n:end)=[];
lastdata=length(tempdata)-stopoffset;
data=cell(1,lastdata-startoffset);qrt=zeros(lastdata-startoffset,24);
for bla=1:lastdata-startoffset
    data{bla}=tempdata{bla+startoffset};
    qrt(bla,1)=tempqrt(bla+startoffset,:);
end    
clear tempdata;
if qrt(1,1)>qrt(lastdata-startoffset,1) %decreasing T, not normal situation
    qrt=flipud(qrt);
    data=flipud(data);
end
qrt(:,13)=power;
qrt(:,5)=blackbody_2(qrt(:,1),pol,method,plt);%
%qrt(:,12)=meh.totphoton;

fclose(fid);
end

function [data,index,qrt,transformation,fitresult] = findparameters(data,qrt,FITF0,dipmin)
%finds parameters of the KID data read in using getdata functions
%this version uses automatically individual phase rotation from the center
%of lowest T res. circle
%FITFO=1 is standard; uses fits to get Fo and S21, 0 handy bad KIDs

stop=size(data,2);
Rfit=zeros(stop,1);maximum=zeros(stop,1);S21res=zeros(stop,1);index=zeros(stop,1);fres=zeros(stop,1);phaseres=zeros(stop,1);Q=zeros(stop,1);reelshift=zeros(stop,1);
for n=1:stop%all measured Temperatures
   filtered=smooth(data{n}(:,2),3); %smoothing
   maximum(n)=max(filtered); %find max value
   data{n}(:,2)=data{n}(:,2)-maximum(n);  %normalise S21 in log space
   data{n}(:,2)=10.^(data{n}(:,2)/20); %from log(magn) to magnitude
   [fres(n),S21res(n),index(n),Q(n),phaseres(n),reelshift(n),Rfit(n)]=kidcal(data{n});%finding cal param
   [data{n}(:,4), data{n}(:,5)]=applykidcal(data{n},phaseres(1),reelshift(1));%callibrate with lowest T circle = first one as reference
   
end

%removing data with shallow dip
for n=1:stop
    if 20*log10(S21res(n))>-1*abs(dipmin)
        disp(['datasets ' num2str(n) ' till ' num2str(stop) ' thrown away because dip shallower than ' num2str(dipmin) 'dB'])
        break
    end
end
qrt(n+1:end,:)=[];
data(n+1:end)=[];
Q(n+1:end)=[];
fres(n+1:end)=[];
S21res(n+1:end)=[];
phaseres(n+1:end)=[];
reelshift(n+1:end)=[];
phaseres(n+1:end)=[];
Rfit(n+1:end)=[];
index(n+1:end)=[];
stop=n;

transformation(:,1)=phaseres;%in radians
transformation(:,2)=reelshift;
transformation(:,3)=Rfit;
qrt(:,2)=Q';
qrt(:,9)=fres';
qrt(:,10)=S21res';
qrt(:,11)=phaseres';

%%%Assymetric fit to data using FitS21_12 Fits F vs Magnitude (col 1 and 2
%%%of data)
%%%Assymetry used to get Q; resonance requency obtained from parabolic fit
a=zeros(1,stop);b=zeros(1,stop);c=zeros(1,stop);l=zeros(1,stop);stretch=zeros(1,stop);tempresult=cell(1,stop);
for n=1:stop%all measured Ts FITTING
    if FITF0
        [a(n),b(n),c(n),l(n),stretch(n),tempresult{n}]=FitS21(data{n}(:,:),qrt(n,9),qrt(n,2),qrt(n,10));   %% fres, Qf, S21min fitted, assymetric
        qrt(n,9)=a(n);qrt(n,10)=c(n); %fres obtained from skewed Lorenztian, S21min
        qrt(n,2)=b(n);  %Q factor Obtained from Fit of parabolic
    end
    qrt(n,4)=qrt(n,2)./qrt(n,10);%Get Qi=Q/S21min
    qrt(n,3)=(qrt(n,4).*qrt(n,2))/(qrt(n,4)-qrt(n,2));%Get Qc=QiQ/(Qi-Q)
    qrt(n,18)=10*log10((2/pi)*(10.^(qrt(n,13)./10).*(qrt(n,2).^2./qrt(n,3))));%internal power

end
%size(qrt)
fitresult=tempresult{1};
end

function  [qrt] = findresponse(data,qrt,Pstart,smprm,transformation)
%
%modified 23/5/07 - amplitude responsivity, from numerical relation of Q
%factor
stop=size(data,2);
[T,index] = min(qrt(:,1));%lowest T data
fres=qrt(:,9); % resonant frequency
%%%%%%%%%%%%%%%%%%%%%%%Finding response %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%loop goes from low T to high T
Rfit=transformation(:,3);
response2=zeros(stop,4);
for n=1:1:stop-1
    response2(n+1,2)=interp1(data{n+1}(:,1)',data{n+1}(:,4)',fres(1),'linear');%real in corrected circle
    response2(n+1,3)=interp1(data{n+1}(:,1)',data{n+1}(:,5)',fres(1),'linear');% Im in corrected circle
    response2(n+1,1)=atan2(response2(n+1,3),-response2(n+1,2));%phase response in rad
    response2(n+1,4)=1/Rfit(1)*(response2(n+1,3).^2+response2(n+1,2).^2).^0.5;%amplitude response, note normalisation
    if response2(n+1,1)<-0.5*pi %to removes steps
        response2(n+1,1)=response2(n+1,1)+2*pi;
    end
end
response2(1,:)=response2(2,:);
qrt(:,14)=qrt(:,9)/qrt(index,9)-1;%dF/F0 = (Fres-Fres0)/Fres0=Fres/Fres0-1 = normalized f response
qrt(:,6)=response2(:,1);%phase response in rad
qrt(:,16)=response2(:,4);%R response
qrt(:,7)=response2(:,2);%real pt
qrt(:,8)=response2(:,3);%imaginary pt
clear response;
%%%%% new - sort for hysteritic data
qrt=sortrows(qrt,5);%
%%%%%%%%%%%%%%%%%%%%%%%% Finfing responsivity by local lineair fit to the data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% finding start point for response measurements %% 
noresponsedata=qrt(:,5)<Pstart;     %low powers is no response
if sum(noresponsedata)~=0
    stdf0=std(qrt(noresponsedata,14));  %std of F0 below where response happens
    meanf0=mean(qrt(noresponsedata,14)); 
    %respdata=-1*qrt(:,14)>2*stdf0+meanf0;      %logical array with data to be analyzed at 1.5 sigma and higher response
    %respdata=-1*qrt(:,14)>(2*stdf0+meanf0); 
    respdata=qrt(:,14)<(meanf0-2*stdf0);
    for n=length(respdata):-1:1         %adjusting logical array, starting from high response find first point that has too small response
        if respdata(n)==0
            break
        end
    end
    minpowerindex=n+1;                  %index where response starts
    respdata(1:minpowerindex-1)=0;      %makes sure that we look at continues data
    if minpowerindex-smprm<1
        error('you try to calculate derivatives outside datarange, reduce smprm value or add more data at low BB power or increase  Pstart');
    end
else
    minpowerindex=1+smprm;
    respdata=zeros(1+smprm+1,1); 
end
%finding derivative for dx, R and theta  by derivative with window of measured data = tru responivity%
%indexing minindex-smprm makes sure that numdiffit gives first result for
%minpowerindex
qrt(minpowerindex-smprm:end,17)=1E15*numdiffit([qrt(minpowerindex-smprm:end,5)*1E15 qrt(minpowerindex-smprm:end,14)],smprm);    %derivative of F0 response
qrt(minpowerindex-smprm:end,19)=1E15*numdiffit([qrt(minpowerindex-smprm:end,5)*1E15 qrt(minpowerindex-smprm:end,6)],smprm);    %derivative of phase response
%derivative of R response. Needs renormalisation as the reference circle is
%not at n=1 but at n=minpowerindex.
qrt(minpowerindex-smprm:end,20)=1E15*Rfit(1)/Rfit(minpowerindex)*numdiffit([qrt(minpowerindex-smprm:end,5)*1E15 qrt(minpowerindex-smprm:end,16)],smprm);    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% find derivative for R, theta response vs PBB. dx by recallibrating data at each tempearture and getting local response  %%%%%
%%% this is the responsivity vs loading if we can retune the KID readout  %%
%res=polyfit(data(n-window:n+window,1),data(n-window:n+window,2),1);
    
%respdata(end-nptsforfit+1:end)=0;      %adjusting logical array to only plot data that has correct resp.


for n=minpowerindex:length(respdata)-smprm
    phaseres=transformation(n,1);   %getting circle cal from dataset {n}
    reelshift=transformation(n,2);
    re=zeros(2*smprm+1,1);im=zeros(2*smprm+1,1);
    datare=cell(2*smprm+1,1);dataim=cell(2*smprm+1,1);
    mnm=1;
    for m=n-smprm:n+smprm
        [datare{mnm}, dataim{mnm}]=applykidcal(data{m},phaseres,reelshift);%callibrate data at for circle at n
        re(mnm)=interp1(data{m}(:,1),datare{mnm},fres(n),'linear');%real in corrected circle
        im(mnm)=interp1(data{m}(:,1),dataim{mnm},fres(n),'linear');% Im in corrected circle
        mnm=mnm+1;
    end
    prsp=atan2(im,-re);        %local phase response in rad
    Rrsp=1/Rfit(n)*(re.^2+im.^2).^0.5;   %local R response
    
%     %     %testplot
%     close all;figure(123);subplot(2,2,1);clear ree iim;
%     for mnm=1:length(re)
%         plot(datare{mnm},dataim{mnm},'-k');hold on;
%         if mnm==smprm+1
%             plot(datare{mnm},dataim{mnm},'-k','LineWidth',2)
%         end
%     end
%     plot(re,im,'ro','MarkerFaceColor','r');
%     subplot(2,2,2);plot(qrt(n-smprm:n+smprm,5),prsp,'bo');
%     subplot(2,2,3);plot(qrt(n-smprm:n+smprm,5),Rrsp,'bo');
%     %end test plot     %

    bla=polyfit(qrt(n-smprm:n+smprm,5)*1E15,prsp,1);    
    qrt(n,21)=bla(1)*1E15;   %slope of the line = dtheta/dPBB at each T
    bla=polyfit(qrt(n-smprm:n+smprm,5)*1E15,Rrsp,1);    
    qrt(n,22)=bla(1)*1E15;   %slope of the line = dR/dPBB at each T
    bla=polyfit(prsp,Rrsp,1)*1E15;    
    qrt(n,25)=bla(1)*1E15; %dR/dtheta
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% getting dRdtheta%%
tofit=qrt(:,6)<pi/4;
for n=1:length(qrt(:,6))
    if tofit(n)==0
        break
    end
end
tofit(n+1:end)=0;
tofit(1:minpowerindex)=0;
drdthetafit=polyfit(qrt(tofit,6),qrt(tofit,16),1);   %drdtheta 
qrt(:,23)=drdthetafit(1);
qrt(:,24)=drdthetafit(1)*qrt(:,6)+drdthetafit(2);%fitted equivalent of qrt(:,6) for drdtheta plot

end

function savethedata(path,qrt,kidnumber,kidpower,ChipID)

resultfile=[path 'KID_' num2str(kidnumber) '_' num2str(abs(kidpower)) 'dBm' ChipID  'Tdep.csv'];
keeswrite('T,Q,Qc,Qi,Pbb,Theta,Rept,Impt,Fres,S21@Fres,theta@res,photonNEP,P,dx,ID,R,dx/dP,Pi,dtheta/dPbase,dR/dPbase,dtheta/dPlocal,dR/dPlocal,dR/dthetabase,dRdthetabasFIT,dR/dthetalocal',resultfile);
dlmwrite(resultfile, qrt, '-append','newline', 'pc', 'precision', '%.12e','delimiter', ',');

end 

function plotdata(data,qrt,kidnumber,kidpower,path,fitresult,Pstart,index,Rfit)
format('long','e');
%h=gcf;
figure(kidnumber);hold off;
clf
 
stop=size(data,2);
%START plotting SUBPLOT 3.3.1 KID RESONANCE
if qrt(1,1)>qrt(stop,1) %temperature
    BLAH=stop;%low T index
else
    BLAH=1;
end
%Single dip squared, squared t show fitresults
subplot(3,3,1);hold off;
xlim([min(data{BLAH}(:,1)) max(data{BLAH}(:,1))]);
plot(fitresult);hold on;ylim([0 1.2]);legend('hide');
plot(data{BLAH}(:,1),(data{BLAH}(:,2)).^2,'.','LineWidth',1);xlim([data{BLAH}(index(BLAH),1)*(1-5/qrt(BLAH,2)) data{BLAH}(index(BLAH),1)*(1+5/qrt(BLAH,2))]);
xlabel('F [GHz]');ylabel('S21 [mag.]');
plot(qrt(BLAH,9),(qrt(BLAH,10)).^2,'ro','MarkerSize',6,'MarkerFaceColor','r');%S21min and Fres used and saved
title(['KID ' num2str(kidnumber) ' @-' num2str(kidpower) ' dBm'])

%SUBPLOT 3.3.2 KID Resonances
subplot(3,3,2);
hold on;
for n=1:stop 
    plot(data{n}(:,1),20*log10(data{n}(:,2)),'-o','MarkerSize',1);
%    plot(data{n}(index(n),1)',20*log10(data{n}(index(n),2))','go','MarkerSize',5,'MarkerFaceColor','g');
    plot(qrt(n,9),20*log10(qrt(n,10)),'ro','MarkerSize',4,'MarkerFaceColor','r');
end

xlabel('F [GHz]');ylabel('S21 [dB]');title('KID resonances');axis tight;%S21 vs freq

%SUBPLOT 3.3.3 KID Resonances in complex plane
subplot(3,3,3);hold off;
xlabel('Re');ylabel('Im');title('KID resonances in complex plane');hold on;%legend(num2str(qrt(:,1)),2);
hold on;
for n=1:stop 
    plot(data{n}(:,4)',data{n}(:,5)','-o','MarkerSize',1);   
    plot(data{n}(index(n),4)',data{n}(index(n),5)','go','MarkerSize',4,'MarkerFaceColor','g');
    plot(qrt(n,7),qrt(n,8),'ro','MarkerSize',4,'MarkerFaceColor','r');
end
rectangle('position',[-Rfit(1),-Rfit(1),Rfit(1)*2,Rfit(1)*2],'curvature',[1,1],'linestyle','-','edgecolor','r');%plots circle 8-)
axis tight;
 
%SUBPLOT 3.3.4 KID response
subplot(3,3,4);hold off;
semilogx(qrt(:,5),qrt(:,6),'-xr');xlabel('Radiated power [W]');ylabel('normalised signal');title('KID response');axis tight;hold on;
semilogx(qrt(:,5),qrt(:,16)-1,'-xb','MarkerSize',3);axis tight;
xlim([2e-17 max(qrt(:,5))]); ylim([-1 2])%SPICA range
legend('phase','R-1')
ptp=logical(qrt(:,17));%points used for derivative
semilogx(qrt(ptp,5),qrt(ptp,6),'ok');
semilogx(qrt(ptp,5),qrt(ptp,16)-1,'ok');

%SUBPLOT 3.3.5 phase and R resp.
subplot(3,3,5);
ptp=logical(qrt(:,19));
semilogx(qrt(ptp,5),qrt(ptp,19)/pi,'-ro');hold on
semilogx(qrt(ptp,5),qrt(ptp,20),'-bx')
ptp=logical(qrt(:,21));
semilogx(qrt(ptp,5),qrt(ptp,21)/pi,'.--r');hold on
semilogx(qrt(ptp,5),qrt(ptp,22),'.--b');axis tight
legend('phase/pi','radius','recentered','recentered');
semilogx(qrt(:,5),zeros(length(qrt(:,1)),1),'-k');
xlabel('P_{BB}');ylabel('normalised power responsivities');


%SUBPLOT 3.3.6 Q Factors
subplot(3,3,6);hold off; %Q data
loglog(qrt(:,5),qrt(:,2),'-ok','MarkerSize',3);
hold on;loglog(qrt(:,5),qrt(:,3),'-or','MarkerSize',3);
hold on;loglog(qrt(:,5),qrt(:,4),'-og','MarkerSize',3);legend('Q','Qc','Qi');
xlabel('P_{BB}  [W]');ylabel('Q ');title('Q factors');xlim([Pstart/10 max(qrt(:,5))])

%SUBPLOT 3.3.7 dx response
subplot(3,3,7);hold off;bla=qrt(:,14)<0;
loglog(qrt(bla,5),abs(qrt(bla,14)),'--xr');xlabel('Radiated power [W]');ylabel('dx=(f-fo)/fo');axis tight;hold on;
ptp=logical(qrt(:,17));%only where we calcualted F responsivity
loglog(qrt(ptp,5),-1*qrt(ptp,14),'ok');
xlim([Pstart max(qrt(:,5))])

%SUBPLOT 3.3.8 KID dx responsivity
subplot(3,3,8);hold off;
ptp=logical(qrt(:,17));
loglog(qrt(ptp,5),-1*qrt(ptp,17),'-x')
xlabel('P_{BB}');ylabel('\delta x/\delta P');

%SUBPLOT 3.3.9 drdtheta
subplot(3,3,9);hold off;
plot(qrt(:,6),qrt(:,16),'--xr');xlabel('\theta');ylabel('R');axis tight;hold on;xlim([-0.3 pi+0.3])
plot(qrt(:,6),qrt(:,24),'-k')
legend('data','fit');
title(['dRdtheta= ' num2str(qrt(1,23))]);

Figfile=[path 'KID_' num2str(kidnumber) '_' num2str(abs(kidpower)) 'dBm_VNABB.Fig'];
saveas(gcf,Figfile,'fig')

end

function [Fresres,Qfres,S21minres,l,stretch,result] = FitS21(data,Fres,Qf,S21min)
%perfoms skewd Lorenztain fit for Q, parabolic fit for Fres and S21min
%Data should contain F[GHz], S21 magn., (phase optional). Fres is res freq
%in GHz, Qf is initial Q factor and S21min is original S21min. These values
%are used as guess values for a fit to Bens eqn. 
%07-03-07: modified for bad ini guesses bandwidth
format('long','e');

no=4; %width in resonator bandwidths to fit over
datasize=size(data(:,2),1);
%limit width of fit
bandwidth=Fres/Qf; %resonator bandwidth
minfitindex=find(data(:,1)>(Fres-no*bandwidth/2),1);
maxfitindex=find(data(:,1)>(Fres+no*bandwidth/2),1);
if minfitindex+4<maxfitindex
    x=data(minfitindex:maxfitindex,1);y=data(minfitindex:maxfitindex,2).^2;
else
    x=data(:,1);y=data(:,2).^2;
end
% estimate linear background
l=(data(1,2)^2-data(datasize,2)^2)/(data(1,1)-data(datasize,1)); %%linear background drift

%skewterm

% catch error in finding guess parameter
low3dbpoint=find(data(:,1)>(Fres-bandwidth/2),1);
if (isempty(low3dbpoint)==1)
    low3dbpoint=1;
end
high3dbpoint=find(data(:,1)>(Fres+bandwidth/2),1);
if (isempty(high3dbpoint)==1)
    high3dbpoint=size(data(:,1),1);
end

stretch=(data(low3dbpoint,2)^2-data(high3dbpoint,2)^2)/(bandwidth); % estimate skew parameter

s = fitoptions('Method','NonlinearLeastSquares', 'Startpoint',[Fres Qf S21min l stretch],'Lower',[Fres*0.9 Qf/2 S21min/4 -25*abs(l) -25*abs(stretch)],...
    'Upper',[Fres*1.1 Qf*2 S21min*2 abs(l)*25 +25*abs(stretch)],'MaxFunEvals',1000);
% figure(1);plot(x,y,'.');hold on;
%Fit Lorentzian like in power space with linear background and stretch
ftype = fittype('(1-(1-Smin^2)*(1+stretch*(x-Fr)/Fr)/(1+(2*Q*(x-Fr)/Fr)^2)+l*(x-Fr)/Fr)','options', s);  
%[result,gof,output]=fit(x,y,ftype);
result=fit(x,y,ftype);
% figure(1);plot(result);hold off;
Qfres=result.Q;
l=result.l;
stretch=result.stretch;
%%new: a simple lorentzian in log space fit to find resonance
%%frequency/S21min. Is done in power space
so=0.3;

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

s = fitoptions('Method','NonlinearLeastSquares', 'Startpoint',[Fres Qf S21min], 'MaxFunEvals',100);

ftype = fittype('10*log10(1-(1-Smin^2)/(1+(2*Q*(x-Fr)/Fr)^2))','options', s);
[resultFres]=fit(x,y,ftype);
Fresres=resultFres.Fr;
S21minres=resultFres.Smin;
end


%old bb scripts
function[powersum,NEP,freq,filter]=blackbody650GHz(TTT,plotdata) 
% 9/9/10, Modified from old blackbodyTtoP of 10/6/09, compatible with SPICA
% L Band, need to set throughput in this program!
% calculates power emitted from blackbody as function of T
% combined main program, calculates blackbody emmission (2 polarisations) for given T array, can
% save data or (savedata=1) and plot (plotdata=1)
% throughput input by etendue, set!
bwidth=0.05;
savedata=0;

centrefreq=650e+9; %% 650Ghz
nopts=200;
nowave=1000;
minwavelength=0.9e-4;
%minwavelength=0.05e-4; %wideband
 
maxwavelength=10e-3;
maxinvwave=1/minwavelength;
mininvwave=1/maxwavelength;

difinvwave=maxinvwave-mininvwave;

h = 6.6262e-34;		% J.s
c = 2.9998e8;		% m/s
k = 1.3806e-23;		% J/K

centrelambda=c/centrefreq;

%Etendue=0;p=0;dp=0;sangle=0;theta=0;TotIrad=0;dEtendue=0;
%etendue=(centrelambda)^2
%etendue=2.299e-3*(0.52*centrelambda)^2*pi % Nyquist sampled airy disc
%etendue=2.299e-3*pi*(3.4e-3/2)^2; %lens 1.5mm, 2.3e-3m^2.str/m^2 from blackbody horn at detector plain
%sangle=0;
%etendue=(1.5e-3)^2*pi*pi
%etendue=2e-7; %3.4mm elliptical lens, 650GHz
%Throughput for 3.4mm diameter lens, in ADR blackbody, multiple energy
%density by lens area.
%etendue=2.299e-3*pi*((3.4e-3)/2)^2; %lens dia 3.4mm, 2.3e-3m^2.str/m^2 from blackbody horn at detector plain
% 1.9 mm lens, 
etendue=2.299e-3*pi*((1.9e-3)/2)^2; %lens dia 1.9mm, 2.3e-3m^2.str/m^2 from blackbody horn at detector plain
bob(:,1)=TTT(:);

%Load filters, in increaseing wavelength
nofilters=2;
bba=flipdim(dlmread('B768 450umBP1K5.txt','\t'),1);
filterform{1}(:,1)=[1/maxinvwave; 1./bba(:,1)/100; 1/mininvwave]; 
bbc=size(bba(:,2),2);
filterform{1}(:,2)=[bba(1,2); bba(:,2);bba(bbc,2)] ;
clear bbc bba;
bba=flipdim(dlmread('W969 37cmLPESCUBAII.txt','\t'),1);
filterform{2}(:,1)=[1/maxinvwave; 1./bba(:,1)/100; 1/mininvwave];
bbc=size(bba(:,2),2);
filterform{2}(:,2)=[bba(1,2); bba(:,2);bba(bbc,2)] ;
clear bbc bba;
for I=1:nowave
    freq(I)=c*((1-I)*difinvwave/(nowave-1)+maxinvwave); % frequency spaced equally decreasing in f space
    filter(I)=1;
    for pee=1:nofilters
        filter(I)=filter(I)*interp1(filterform{pee}(:,1),filterform{pee}(:,2),c/freq(I));
    end
end

for i=1:size(TTT,1)
    for I=1:nowave
        brill(I)=(2*h*freq(I)^3)/(c^2)/(exp(h*freq(I)/(k*TTT(i)))-1); %%brillance W.m^2.str/Hz
        filterirad(I)=0;
        %filter calculations
        irad=brill*etendue;  
        filterirad(I)=irad(I).*filter(I);

    end
    %numerical integration of brilliance
    powersum(i)=sum(filterirad,2)*c*difinvwave/(nowave); 
    NEP(i)=sqrt(2*powersum(i)*centrefreq*h); 
end
bob(:,2)=powersum;
bob(:,3)=NEP(:);
if plotdata==1
    subplot(3,1,1)
    plot(bob(:,1),bob(:,3))
    ylabel('Photon NEP W/Hz^1/2')
    xlabel('T (K)')
    subplot(3,1,2)
    plot(bob(:,1),bob(:,2))
    ylabel('Power in W')
    xlabel('T (K)')
    subplot(3,1,3)
    plot(c./filterform{1}(:,1)/1e12,filterform{1}(:,2),c./filterform{2}(:,1)/1e12,filterform{2}(:,2),freq(:)/1e12,filter(:))
    ylabel('Transmission')
    xlabel('Frequency (THz)')
    axis([0 4 0 1])
end
if savedata==1
    poo=what;
    pathy=poo.path;
    result=[pathy '\bbscan.csv'];
    saveas(gcf,[pathy '\result' '.fig'],'fig')

    casewrite(['centre f=' num2str(centrefreq) 'bwidth' num2str(bwidth) 'throughput=' num2str(etendue) 'f,Power,NEP,Equivilant T stability'],result)
    dlmwrite(result, bob, '-append','newline', 'pc','precision', '%9e', 'delimiter', ',')
end
end

function [Pbb,PhotonNEP,freq,filterT]=SPICA_Lband(TTT,plaatje) 
% 10/6/09
% calculates power emitted from blackbody as function of T
% Integrates over the entire range of filter data that is loaded in,
% assumes 0 transmission outside filter band
% bob= cols P and photon NEP

h = 6.6262e-34;		% J.s
c = 2.9998e8;		% m/s
k = 1.3806e-23;		% J/K

filter=dlmread('totalLbandilters.txt','\t');
if filter(2,1)<filter(1,1) %DEcreasing in frequency
    filter=flipdim(filter,1);
end
filter(:,1)=100*filter(:,1)*c; %converting to f space from cm-1

%making sure array has equal f spacing
npts=length(filter);
df=(filter(end,1)-filter(1,1))/npts; %mean df between points
freq=zeros(npts,1);filterT=zeros(npts,1);
for n=1:npts
    freq(n)=(df*(n-1)+filter(1,1));%equally spaced f array
    filterT(n)=interp1(filter(:,1),filter(:,2),freq(n));
end

[maxFt,maxFti]=max(filterT); 
centrefreq=filter(maxFti,1);
ntemperatures=length(TTT);
powersum=zeros(ntemperatures,1);NEP=zeros(ntemperatures,1);
for mm=1:ntemperatures                                           %all temperatures
    filterirad=zeros(1,size(filter,1));
    for n=1:npts
        brill=(2*h*freq(n)^3)/(c^2)/(exp(h*freq(n)/(k*TTT(mm)))-1); %brillance W.m^2.str/Hz
        irad=brill*(c/freq(n)).^2;                                  %radiated power from BB using lambda^2 throughput
        filterirad(n)=irad.*filterT(n);                                %power transmitted through filter
    end
    %numerical integration of brilliance
    powersum(mm)=sum(filterirad)*df; 
    NEP(mm)=sqrt(2*powersum(mm)*centrefreq*h); 
end

Pbb=powersum;
PhotonNEP=NEP;
if plaatje==1
    figure(2)
    subplot(1,3,1)
    semilogy(TTT,PhotonNEP)
    ylabel('Photon NEP W/Hz^1/2')
    xlabel('T (K)')
    subplot(1,3,2)
    semilogy(TTT,Pbb)
    ylabel('Power in W')
    xlabel('T (K)')
    subplot(1,3,3)
    plot(freq,filterT)
    ylabel('Transmission')
    xlabel('Frequency (Hz)')
end

end

function [der,offs]=numdiffit(data,window)
%data is 2 cols x and y
%function will fit linear over 2*window+1 number of points, no padding at the
%ends. saves slope and offset of fitted line

npts=length(data);
der=zeros(npts,1);offs=zeros(npts,1);
for n=1+window:npts-window
   res=polyfit(data(n-window:n+window,1),data(n-window:n+window,2),1);%lineair fit
   der(n)=res(1);
   offs(n)=der(2);
end
end

function [fres,S21min,indexfres,Q,rotation,reelshift,Rfit]=kidcal(data)
%data is F GHz magnitude (not in dB) rad (theta)
%range sets #pts for lineair fit
%function finds all callibration parameters but does not apply them

F=data(:,1);magn=data(:,2);phase=data(:,3);
S21max=(magn(1)+magn(end))/2;
[S21min,indexfres]=min(magn);
fres=F(indexfres);
dBpt=sqrt((S21min^2+S21max^2)/2);
%find 3 dB point indices
[temp,bwindexdown]=min(abs(data(1:indexfres,2).^2-dBpt^2));
[temp,bwindexup]=min(abs(data(indexfres:end,2).^2-dBpt^2));
bwindexup=bwindexup+indexfres-1;
if bwindexdown<2
    bwindexdown=2*indexfres-bwindexup;
else if bwindexup>length(F)
        bwindexup=2*indexfres-bwindexdown;
    end
end

Q=fres/(F(bwindexup)-F(bwindexdown));

%now finding the rotation parameters for the circle
reel=magn.*cos(phase);
imagi=magn.*sin(phase);
%finding shifts by fitting circle to the data within 3dB pts
x=reel(bwindexdown:bwindexup);y=imagi(bwindexdown:bwindexup);
[xfit,yfit,Rfit]=circfit(x,y);
reelshift=sqrt(xfit^2+yfit^2);
rotation=atan2(yfit,xfit);
%getting Q factor

% %test plot
% figure(122334)
% plot(reel,imagi,'b.');hold on;
% rectangle('position',[xfit-Rfit,yfit-Rfit,Rfit*2,Rfit*2],...
%     'curvature',[1,1],'linestyle','-','edgecolor','r');%plots circle 8-)
% [recal,imcal]=applykidcal(data,rotation,reelshift);
% hold on;plot(recal,imcal,'.');hold on
% rectangle('position',[-Rfit,-Rfit,Rfit*2,Rfit*2],...
%     'curvature',[1,1],'linestyle','-','edgecolor','g');%plots circle 8-)

end

function [reelcal, imcal]=applykidcal(data,rotation,reelshift)
%data = F[GHz] magn (NOT dB) phase rad. rotation and reelshioft are
%obtained using kidcal on same data

magn=data(:,2);phase=data(:,3);

newphase=phase-rotation;
reelcal=magn.*cos(newphase)-reelshift;
imcal=magn.*sin(newphase);

%plot(reel,im);
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
   lf = char(13);%were both setstr before
else
   lf = char(10);
end

lf = lf(ones(size(strmat,1),1),:);
lines  = [strmat lf]';

fprintf(fid,'%s',lines);
fclose(fid);
end

function   [xc,yc,R,a] = circfit(x,y)
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
   a=[x y ones(size(x))]\(-(x.^2+y.^2));%square bracjedts used before here
   xc = -.5*a(1);
   yc = -.5*a(2);
   R  =  sqrt((a(1)^2+a(2)^2)/4-a(3));
end

% 1	T
% 2	Q
% 3	Qc
% 4	Qi
% 5	PBB in W
% 6	phase response in rad wrt lowest T circle
% 7	real point wrt lowst T circle
% 8	im point wrt lowst T circle
% 9	Fres
% 10	S21@ res (magn, not dB)
% 11	phase at resonance
% 12	photon noise NEP
% 13	Preadout
% 14	f_{res}/F_{res,0}-1
% 15
% 16	R response wrt lowest T circle
% 17	normalised Fo responsivity (derivative of 14 vs PBB)
% 18	Pi in dBm
% 19	normalised phase responsivity (derivative of 6 vs PBB)
% 20	normalised R responsivity (derivative of 16 vs PBB)
% 21	phase responsivity at temperature T = der of dtheta/dPBB around T(n)
% 22	R responsivity at temperature T = der of dR/dPBB around T(n)
% 23	slope of dR/dtheta at base temperature (local dR/dtheta)
% 24	fitted line of dR/dtheta (=dtheta from fit =dR/dtheta(fit) * theta)
% 25 	dR/dtheta at temperature T = der of dR/dtheta around T(n)
% 26	estimate in P where phase response is 10% of max: instantaneaous dyn range

%data has F magnitude(normalised, NOT in log) theta in rad Re(rotated) Im(rotated)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
