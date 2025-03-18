function KIDparam=Load_ResponseV1
% Hacked from the NEP-vs_Loadig to plot the local response of the KIDs

% NOTE(PdV Jan2013): NEP_vs_loading has been splitted in A,B and C as of
% vs9. While editing please make sure changes in A,B or C do not conflict
% with the other functions, since you need all of them for a complete analysis
% The output workspace contains now a lot of variables to carry on to C,
% but only KIDparam is the real output.
% OUTPUT is described at the end of the function
%
% runs on the directories created by the optical NEP labview script
%
%
% method. A struct containing information about the used setup including:
%         filters, lens size, and aperture determination. The following
%         struct elements are required:
%
%   method.filter='325 GHz AB'      %AB old BP filter for A-MKID
%   method.filter='350 GHz WT1275'  %Cardiff filters (only 2)
%   method.filter='325 GHz He7'     %He7 system filtering used in A-MKID L band He7 cooler
%   method.filter='1_6THz'          %1.6 THz SPICA SAFARI stack
%
%   method.lensdiameter=1 %lens diameter in millimeters
%   method.lensdiameter=2 %lens diameter in millimeters
%
%   method.aperture='He7'                   %He7 system
%   method.aperture='ADR small'             %350 GHz small 2 mm diameter aperture in the ADR at 14 mm above the chip.
%                                               Using beam pattern calculation values for throughput determination.
%   method.aperture='ADR small Geometrical' %small 2 mm diameter aperture in the ADR at 14 mm above the chip.
%                                               Using geometrical calculation for throughput determination.
%   method.aperture='ADR full'              %full aperture in the ADR
%
%   method.aperture='ADR full 1.6 THz'      %full aperture ADR with 1.6 Thz SAFARI stack [only for 2 mm lens] (JB,PdV 21-9-2012)
%
%   method.freqresolution = 2 %[Optional, Default = 2] frequency resolution in GHz used
%                               to integrate the blackbody spectrum with all filters.
%                               NOTE: for the 350 GHz and 325 GHz filters 2 GHz is
%                                       adviced as minimum.
%   method.Delta = 45.6 %[Optional, Default = 45.6=Al] value of Delta (half
%                        the superconducting gap) in GHz. Only used for G-R
%                        NEP estimate by Pblackbody_4.m
%   method.eta_pb = 0.57 %[Optional, Default = 0.57] value of the photon
%                        pair breaking efficiency. Default from Kozorezov
%                        et.al. PRB 2000
%   NOTE: 2 mm lenses only implemented for ADR small and ADR small Geometrical


%%%%%%%%%%%%%%%%%%%%%%%%%Parameters to Set%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%==========================================================================
% Setting some routine default values
%==========================================================================
scheiding=filesep;                                              %for poor Mac users use '/' 

path=[cd '/..']; %root path where data is, one higher than the scripts
resperpath=[filesep 'Response' filesep 'Time' filesep];
Plottype=1;%plots vs measured variable, 1 plots vs index
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%NB call global resppathy in main%%%
closefigures=0;         %close figures after saving
readallKIDs=1;          % 0 does not read all KIDs, then KIDstr=[12 14] or so should give the KIDs to read. 
%NB: use readallKIDs=1 before proceeding to the next scripts
KIDstr=[9];         % [3 7 9]if not all KIDs are to be read these ones are
pol=1;                  % 1 calc for 1 polarisation, 2= for 2

method.filter='850 GHz'; %Cardiff filter
method.lensdiameter=1; %in mm
method.aperture='ADR small Geometrical'; %small aperture in the ADR, 0.071rad is the angle to 1 mm radius hole, theta=atan(1mm/14mm), angle to normal, TP=2pi(1-costheta)
method.eta_pb = 0.57; %pairbreaking efficiency.
method.EpsTrans = 1.0; %optical efficiency due to translation wrt pinhole
%Known EpsTrans values
%First measurement H20 chip, KID16 0.79\pm 0.05
%First measurement H20 chip, KID3 0.71\pm 0.05

kb = 1.38e-23; %J/K
h = 6.626e-34; %J/Hz
Tc = 1.283; %K
method.Delta = 1.76*kb*Tc/h/1e9; %Half gap in GHz. 

PdVtoys=0; %some optional stuff for further analysis (under construction, PdV) %now: fit responsivity vs Pread and draw GRnoiseline

            
Prange = 0.1; %Halfwidth of the range given 
% EITHER (RRT=1 or 4) absolute value in [W] OR (RRT=2 or 5) fraction of central optical power

%The power range over which the phase response is fitted can be adjusted
%wrt to power range over which the amplitude reponse fits happens.
notheta=1;  %(Power range of phase fit)/(power range of amplitude fit)
%In general notheta <= 1, because the phase response is stronger than the
%amplitude response.

%%%%%%%%%%%%%%%%%%%%%%%%%  Defaults Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
maxnumkids=100;         % maximum number of KIDs that can be analyzed (used fro allocating variables)
maxnopowers=20;         % maximum number P 
maxpts=5000;            %max # points in any datafile at 1 power and 1 T for 1 ID

format('long','e');                                             %Set display format of numbers to 7 digits
set(0,'defaultAxesFontName', 'Arial');                          %readable axis fonts on Mac of JB 8-)
set(0,'defaultAxesFontSize', 13);                               %poor PdV does not have such an enormous screen :(
set(0,'defaultTextFontName', 'Arial');
set(0,'defaultTextFontSize', 13);

addpath([pwd,filesep,'subroutines']);                           %Enable subroutines by adding path in search path.

kleur='RainbowReinier';
% The available colormaps are:
% *RainbowReinier [Default]: Blue-Cyan-Green-Yellow-Red. Often found in MSc
%                            Thesis R.M.J. Janssen.
% *RBG_Miranda: Red-Black-Green. Created to make a chemical heat map for
%                            M.G.M. Kok
% *RWB_Miranda: Red-White-Blue. Created to make a chemical heat map for
%                            M.G.M. Kok
% *PGO_Miranda: Purple-Green-Orange. Created to make a chemical heat map for
%                            M.G.M. Kok
% *PieterdV: Stolen from Pieter de Visser
%This routine requires the 'colormapstorage.mat' file to be present.

%==========================================================================
% Warning shutdown
%==========================================================================
%In this program it often happens that data containing negative values is
%plotted in log scales. This creates matlab warnings. Usually the negative
%data means something silly has gone wrong or a default, do not look at
%this values has been set. Thus the warning can generally be ignored.
w(1).identifier = 'MATLAB:Axes:NegativeDataInLogAxis';
for nw=1:length(w)
    warning('off',w(nw).identifier)
end
%To add warnings, use the line
%   w = warning('query','last')
%directly after the warning is generated. Then take the identifier and add
%it to the structlist above (before loop)
%   w(nw).identifier = 'MATLAB:AnotherUselessWarning';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  START  REAL  PROGRAM  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;     %close all open plots to remove clutter.
clc
resppathy=[path resperpath];                                    %contains directories with the data. 1 dir for each background
format('long','g');
dirs=dir(resppathy);

nn=1;
dirnum=zeros(1,length(dirs));
for n=1:length(dirs)
    %temp=str2double(dirs(n).name);
    if dirs(n).isdir==1
        if ~isnan(str2double(dirs(n).name))
            dirnum(nn)= str2double(dirs(n).name);
            nn=nn+1;
        end
    end
end
noBBTS=nn-1;
dirnum(noBBTS+1:end)=[];
dirnum=sort(dirnum);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%              Reading in dir structure                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

colorsPBB = GenerateColorsFromMap(noBBTS,kleur);   %colorsPBB is for BB power
BBcal.T=0; 
kidsP=zeros(maxnumkids,maxnopowers,noBBTS);
Temperaturenoisefile=zeros(noBBTS,maxnumkids,maxnopowers);

    
for PBB_n=1:noBBTS                                  %PBB_n runs over all PBB, i.e. all dirs. This loop is done dir for dir
    resppath=[path resperpath num2str(dirnum(PBB_n)) scheiding];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %%%%%%%%%%%%%%     Finding response files      %%%%%%%%%%%%%% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
    temprespfiles=dir([resppath '*ReImT*.dat']);
    
    %removing all files that are crap, i.e. either longer or shorter than the mean%
    seis=zeros(1,length(temprespfiles));
    for n=1:length(temprespfiles)
        seis(n)=temprespfiles(n).bytes;
    end
    datatoberead=seis>0.5*mean(seis) & seis<1.5*mean(seis);
    respfiles=temprespfiles(datatoberead);
    clear seis;
    
    trespkids=zeros(length(respfiles),2);
    for nn=1:length(respfiles)
        if isempty(respfiles(nn).name)==1
            error('file empty');
        end
        temp=cell2mat(textscan(respfiles(nn).name,'%*3s %f %*c %f %*s'));
        trespkids(nn,1)=temp(1);
        trespkids(nn,2)=temp(2);
        
    end
    allkidfiles=sortrows(trespkids,[1 2]);%col1 has KID ID, col 2 has KID P, all sorted, i.e. increasing KIDID and per ID increasing P
    
    %get all KID ID's in 1 array and put all powers in a cell array for each KID
    respkids=zeros(maxnumkids,1);
        
    bla=1;bob=0; %size(allkidfiles);
    respkids(bla)=allkidfiles(bla,1);
    
    for n=1:size(allkidfiles,1)
        if respkids(bla)==allkidfiles(n,1);
            bob=bob+1;
            kidsP(bla,bob,PBB_n)=allkidfiles(n,2);
        else
            bla=bla+1;
            bob=1;
            respkids(bla)=allkidfiles(n,1);
            kidsP(bla,bob,PBB_n)=allkidfiles(n,2);
        end
    end
    respkids(n+1:end)=[];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    %%%%%%%%%%%%%%     Reading in data + process KID by KID      %%%%%%%%%%%%%% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    
    if readallKIDs==0
        KIDitoread=zeros(length(KIDstr),1);
        for n=1:length(KIDstr)
            KIDitoread(n)=find(KIDstr(n)==respkids);%index of the KID in respKIDs that we want to read
        end
    else %DO NOT CHANGE THIS PART IT IS WORKING
        KIDitoread=zeros(length(respkids),1);
        for n=1:length(respkids)
            if respkids(n)~=0
                KIDitoread(n)=n;%index of the KID in respKIDs that we want to read
            else
                KIDitoread(n)=0;
            end
        end
    end
    
    
    nokids=length(find(KIDitoread~=0));
    for n=1:nokids                                       % everyKID in the particular PBB folder
        %DO NOT CHANGE THIS PART (net 3 lines) IT IS WORKING for many
        %readout powers
        nKID=KIDitoread(n);%index of the KID we read in in the respkids array (and all other arrays)
        nopread=length(find(kidsP(nKID,:,PBB_n)~=0));
        colorPreadout=GenerateColorsFromMap(nopread,'RainbowReinier');
        for m=1:nopread                                     %m=every readout power for that KID at that PBB
 
            clear tempqrt S21Re S21Im response S21F S21dB S21rad data;
            tempqrt=zeros(maxpts,22);

            %%%%%%%%%%%   Reading S21 file GHz Re Im  that was measured together with the response file  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            filetjes=dir([resppath 'KID' num2str(respkids(nKID)) '_' num2str(kidsP(nKID,m,PBB_n)) '*_S21.dat']);%
            if length(filetjes)==1
                S21file=[resppath filetjes.name];
            else
                error('identical S21 files found')
            end
            fid=fopen(S21file);                             %GHz Re Im of KID circle that is normalized
            if fid==-1 %%check for exists
                error('S21 file not found');
            end
            textscan(fid,'%*[^\r\n]',6);                    %scanning 4 empty lines to skip header
            KIDparam(nKID).S21_IQplane{m,PBB_n}=cell2mat(textscan(fid, '%f%f%f'));
            S21Re=KIDparam(nKID).S21_IQplane{m,PBB_n}(:,2);
            S21Im=KIDparam(nKID).S21_IQplane{m,PBB_n}(:,3);
            KIDparam(nKID).S21Real{m,PBB_n}=S21Re;
            KIDparam(nKID).S21Imag{m,PBB_n}=S21Im;
            fclose(fid);
            [~,~,KIDparam(nKID).Rfit(m,PBB_n)] = circfit(S21Re,S21Im);
            KIDparam(nKID).KIDid=respkids(nKID);
            %%%%%%%%%%% Reading response file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            filetjes=dir([resppath 'KID' num2str(respkids(nKID)) '_' num2str(kidsP(nKID,m,PBB_n)) '*_ReImT.dat']);%
            if length(filetjes)==1
                respfilename=[resppath filetjes.name];
            else
                error('identical response files found')
            end
            disp(respfilename);
            fid=fopen(respfilename);
            if fid==-1                                          %check for exists
                disp('no datafile line 231')
            end
            textscan(fid,'%*[^\r\n]',4);                        %scanning empty lines to skip rest of header
            temp=textscan(fid,'%*s%*s%*2c%f%*s',1);
            ChipTemperature=temp{1};%chip temperature
            textscan(fid,'%*[^\r\n]',1);
            response=cell2mat(textscan(fid, '%f%f%f%f%f%f'));
            endcut=1;                                           % 1 cut end of response data, fraction
            tempqrt(1:round(length(response(:,1))*endcut),1)=response(1:round(length(response(:,1))*endcut),1);%TBB
            tempqrt(round(length(response(:,1))*endcut)+1:end,:)=[];
            tempqrt(:,11)=ChipTemperature;
            tempqrt(:,13)=-1*kidsP(nKID,m,PBB_n);               %readout power from response
            tempqrt(1:round(length(response(:,1))*endcut),7)=response(1:round(length(response(:,1))*endcut),2);%re on cirvle
            tempqrt(1:round(length(response(:,1))*endcut),8)=response(1:round(length(response(:,1))*endcut),3);%Im on cirle

            fclose(fid);
                                                                %cleaning up qrt: make unique, sort and remove hysteresis
            mmm=1;clear b;b=zeros(length(tempqrt(:,1))-1,1);
            for nn=1:length(tempqrt(:,1))-1                     %cleaning up orignal TD data
                if abs(tempqrt(nn,1)-tempqrt(nn+1))/tempqrt(nn,1)>0.2 || abs(tempqrt(nn,1)-tempqrt(nn+1))/tempqrt(nn,1)<1E-4 %T setp too small or too large
                    b(mmm)=nn;
                    mmm=mmm+1;
                end
            end
            b(mmm:end)=[];
            tempqrt(b,:)=[];
            KIDparam(nKID).Tbb{m,PBB_n}=tempqrt(:,1);
            KIDparam(nKID).Tchip(m,PBB_n)=tempqrt(1,11);
            KIDparam(nKID).Pread(m,PBB_n)=tempqrt(1,13);
            KIDparam(nKID).Reresp{m,PBB_n}=tempqrt(:,7);
            KIDparam(nKID).Imresp{m,PBB_n}=tempqrt(:,8);

            % Call new bb script NOT!
            if Plottype==1 %Change X axis of plot into index (for Time measuremets)
                for nn=1:length(KIDparam(nKID).Tbb{m,PBB_n})
                    KIDparam(nKID).Tbb{m,PBB_n}(nn)=nn;
                end
            end
            KIDparam(nKID).Radiusresp{m,PBB_n}=sqrt(KIDparam(nKID).Reresp{m,PBB_n}.^2+KIDparam(nKID).Imresp{m,PBB_n}.^2)/KIDparam(nKID).Rfit(m,PBB_n);
            KIDparam(nKID).Phaseresp{m,PBB_n}=atan2(KIDparam(nKID).Imresp{m,PBB_n},-1*KIDparam(nKID).Reresp{m,PBB_n});
            
            %%%%%%%%%%%Reading S21 file GHz dB rad  that was measured together with the KID resp%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            filetjes=dir([resppath 'KID' num2str(respkids(nKID)) '_' num2str(kidsP(nKID,m,PBB_n)) '*_S21dB.dat']);%
            if length(filetjes)==1
                S21file=[resppath filetjes.name];
            else
                error('identical S21dB files found')
            end
            fid=fopen(S21file);
            if fid==-1
                error('file not found line 122');
            end
            textscan(fid,'%*[^\r\n]',2);
            textscan(fid,'%*s%*s%*s%*s%*1c%*f',1);              %last%f has Fres, not  used anynmore
            temp=textscan(fid,'%*2c%f%*s%*s%*s%*7c%f',1);
            KIDparam(nKID).Q(m,PBB_n)=temp{1};
            KIDparam(nKID).S21min(m,PBB_n)=10^(temp{2}/20);     %in magnitude space

            textscan(fid,'%*[^\r\n]',2);
            KIDparam(nKID).S21_MPplane{m,PBB_n} = cell2mat(textscan(fid, '%f%f%f'));
            fclose(fid);
            %correcting length
            S21F=KIDparam(nKID).S21_MPplane{m,PBB_n}(:,1);%GHz
            S21dB=KIDparam(nKID).S21_MPplane{m,PBB_n}(:,2);
            
            %==============================================================
            %Normalize the S21 data measured in magnitude plane
            filtered=smooth(KIDparam(nKID).S21_MPplane{m,PBB_n}(:,2),3); %smoothing the |S21| data in dB space
            %normalise S21 in log space to the max(|S21|)
            KIDparam(nKID).S21_MPplane{m,PBB_n}(:,2)=KIDparam(nKID).S21_MPplane{m,PBB_n}(:,2)-max(filtered);
            %Convert dB to magnitude
            KIDparam(nKID).S21_MPplane{m,PBB_n}(:,2) = 10.^(KIDparam(nKID).S21_MPplane{m,PBB_n}(:,2)/20);

            %%%%%%%%%%%%%% Getting KID parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            KIDparam(nKID).Qi(m,PBB_n)=KIDparam(nKID).Q(m,PBB_n)/KIDparam(nKID).S21min(m,PBB_n);
            KIDparam(nKID).Qc(m,PBB_n)=(KIDparam(nKID).Qi(m,PBB_n)*KIDparam(nKID).Q(m,PBB_n))/(KIDparam(nKID).Qi(m,PBB_n)-KIDparam(nKID).Q(m,PBB_n));
            [KIDparam(nKID).S21min(m,PBB_n),s21minindex]=min(smooth(S21dB,5));%NO FITTING TO FIND S21min (dB) and F0
            KIDparam(nKID).fres(m,PBB_n)=S21F(s21minindex);
            KIDparam(nKID).Pint(m,PBB_n)=10*log10(10.^(KIDparam(nKID).Pread(m,PBB_n)/10)*(2/pi)*(KIDparam(nKID).Q(m,PBB_n)^2/KIDparam(nKID).Qc(m,PBB_n)));%9=Pi10*log10((2/pi)*(10.^(qrt(nKID,13)./10).*(qrt(nKID,2).^2./qrt(nKID,3))))
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Getting the noise file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %%%%%%%%%%%%%%%%%%%%%% call new BB script to find the background power %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             [KIDparam(nKID).Pbbnoise(m,PBB_n),poo,BBcal]=blackbody_int(Temperaturenoisefile(PBB_n,nKID,m),BBcal,pol,method);
%             KIDparam(nKID).Tbbnoise(m,PBB_n)=Temperaturenoisefile(PBB_n,nKID,m);
%             KIDparam(nKID).wave(m,PBB_n)=poo.wave;
%             KIDparam(nKID).g_r(m,PBB_n)=poo.g_r;
%             KIDparam(nKID).poisson(m,PBB_n)=poo.poisson;
%             KIDparam(nKID).totphoton(m,PBB_n)=poo.totphoton;
%             clear meh poo
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%Processing data: finding the phase radius and dx response  %%%%%%%%%%%%%%
            KIDparam(nKID).df{m,PBB_n}=-1*KIDparam(nKID).Phaseresp{m,PBB_n}/(4*KIDparam(nKID).Q(m,PBB_n)); %estimation of F response
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%now finding response at base =  used for the NEP calcul%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            %%%%%%%%%%%%%%%%%%%%%%%% Finding responsivity by local lineair fit to the data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %save response file with Tbb,Pbb,radius resp, phase resp, frequency resp
%             respresultsfile = [resppath 'KID' num2str(respkids(nKID)) '_' num2str(kidsP(nKID,m,PBB_n)) 'dBm_Responses.csv'];
%             WriteSRONcsv(respresultsfile,[KIDparam(nKID).Tbb{m,PBB_n} KIDparam(nKID).Pbb{m,PBB_n} KIDparam(nKID).Radiusresp{m,PBB_n} KIDparam(nKID).Phaseresp{m,PBB_n} KIDparam(nKID).df{m,PBB_n}],'Tbb, Pbb, radiusresp, phaseresp, freqresp \nKID','%.12e');
    
        end
        
        
%        disp(['*********** index= ' num2str(PBB_n) ', Pbb= ' num2str(KIDparam(nKID).Pbbnoise(m,PBB_n)) ', TBB= ' num2str(KIDparam(nKID).Tbbnoise(m,PBB_n)) '***********'])
        
        %%%%%%%%%%%%%%%%Calculate Popt for every PBB, ie minimise in NEP
       
        
        %radresponse=figure(respkids(nKID)*10000000+PBB_n);
        %%%%%%%%%%%%%%%Make plots that will be generated for every PBB,KID to show readout power dependences%%%%%
        preadfig=figure(respkids(nKID)*10000+PBB_n);
                
        subplot(2,2,1)
        semilogy(-nonzeros(kidsP(nKID,:,PBB_n)),nonzeros(KIDparam(nKID).Q(:,PBB_n)),'s',-nonzeros(kidsP(nKID,:,PBB_n)),nonzeros(KIDparam(nKID).Qi(:,PBB_n)),'d',-nonzeros(kidsP(nKID,:,PBB_n)),nonzeros(KIDparam(nKID).Qc(:,PBB_n)),'^')
        xlabel('P_{read}(dBm)');ylabel('quality factor from S21');legend('Q','Qi','Qc')
        
        for mpread=1:length(nonzeros(kidsP(nKID,:,PBB_n)))
            
            subplot(2,2,2)
            plot(KIDparam(nKID).Tbb{m,PBB_n},KIDparam(nKID).Phaseresp{mpread,PBB_n},'s','color',colorPreadout(nopread+1-mpread,:));hold on;
            plot(KIDparam(nKID).Tbb{mpread,PBB_n},KIDparam(nKID).Radiusresp{mpread,PBB_n}-1,'+','color',colorPreadout(nopread+1-mpread,:))
            title('Response + fit')
            legend('Phase','Radius-1')
            ylabel('Reponse')
            xlabel('Scanned parameter')
            title(num2str(KIDparam(nKID).KIDid))
%            plot(KIDparam(nKID).Pbb{mpread,PBB_n}*1E15,KIDparam(nKID).respfitR{mpread,PBB_n},'color',colorPreadout(nopread+1-mpread,:));% fits
 %           plot(KIDparam(nKID).Pbb{mpread,PBB_n}*1E15,KIDparam(nKID).respfitth{mpread,PBB_n},':','color',colorPreadout(nopread+1-mpread,:));
            
           
            subplot(2,2,4)
            pbla(mpread)=plot(KIDparam(nKID).S21Real{mpread,PBB_n},KIDparam(nKID).S21Imag{mpread,PBB_n},'-k','LineWidth',1,'color',colorPreadout(nopread+1-mpread,:));hold on;%volt S21
            xlabel('Re');ylabel('Im ');
            axis tight;
            plot(KIDparam(nKID).Reresp{mpread,PBB_n},KIDparam(nKID).Imresp{mpread,PBB_n},'.','color',colorPreadout(nopread+1-mpread,:))
            %legend('Fscan','TD data')
            blalegend{mpread,1}=[num2str(-kidsP(nKID,mpread,PBB_n)) ' dBm']; %T in mK
            
        end
        legend(pbla, blalegend); %legend to the circle figure, subplot(3,3,9), wich shows the colour coding for the readout powers which appears in all figures.
        clear pbla blalegend 
        Figfile=[resppath 'KID_' num2str(respkids(nKID)) '_' num2str(KIDparam(nKID).Tchip(1,PBB_n),'%.3f') 'K_Preaddependence.fig'];
        saveas(preadfig,Figfile,'fig');
        if closefigures==1
            close(preadfig);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Creating plots and saving figures%%%%%%%%%%%%%%%%%%%%%5%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end
end

%clear the variables that are not needed for function C
clear BBcal blalegend ChipTemperature Figfile Pofset S21F S21dB S21Im S21Re S21file TDpreadindex TDtempindex Temperaturenoisefile allkidfiles...
    ans b bla bob dRdPfit datatoberead dirnum dirs dthetadPfit endcut factorbound fid filetjes fithighbound fitlowbound...
    m maxnopowers maxnumkids maxpts mmm mpread n nKID nn noisefilename nopread phasefitresult phaseresult radiusfitresult radiusresult...
    respfiles response resppath result2 result2a result2p s s21minindex searchtemp stddRdPfit stddthetadPfit temp tempqrt temprespfiles trespkids...
    weightsdR

save([resppathy 'KIDparam.mat']) %Need to save whole workspace, not only KIDparam, for use in B and C (PdV)

rmpath([pwd,filesep,'subroutines']);

fprintf('`load Response Finished \n')

%%%%%%%%%%%%%%%%%OUTPUT (only KIDparam is useful, the rest are dummy variables that need to be carried to function C)
% KIDparam(nKID).S21Real{mpread,PBB_n}                      %Data of resonance circle, real part, cell
% KIDparam(nKID).S21Imag{mpread,PBB_n}                      %Data of resonance circle, imag part, cell
% KIDparam(nKID).S21_MPplane{mpread,PBB_n}(:,1:3) = cell array. Each cell contains a double array with the columns:[F (GHz), |S21|, phase (rad)]
% KIDparam(nKID).S21_IQplane{mpread,PBB_n}(:,1:3) = cell array. Each cell contains a double array with the columns:[F (GHz), I, Q]
% KIDparam(nKID).Rfit[mpread,PBB_n]                         %Radius from circfit, double
% KIDparam(nKID).Tbb{mpread,PBB_n}                          %Blackbody temperature of response measurement, corresponds to Reresp,Imresp, cell
% KIDparam(nKID).Tchip[mpread,PBB_n]                        %Measured chip temperature (is rounded to only one digit for optical measurements), double
% KIDparam(nKID).Pread[mpread,PBB_n]                        %Microwave readout power, double
% KIDparam(nKID).Reresp{mpread,PBB_n}                       %Measured optical response upon BB sweep, real part, cell
% KIDparam(nKID).Imresp{mpread,PBB_n}                       %Measured optical response upon BB sweep, imag part, cell
% KIDparam(nKID).Pbb{mpread,PBB_n}                          %Blackbody power of response measurement, corresponds to Tbb, cell             
% KIDparam(nKID).Radiusresp{mpread,PBB_n}                   %Measured optical response upon BB sweep, amplitude, cell
% KIDparam(nKID).Phaseresp{mpread,PBB_n}                    %Measured optical response upon BB sweep, phase, cell
% KIDparam(nKID).Q[mpread,PBB_n]                            %Loaded Q as fitted from S21, double
% KIDparam(nKID).S21min[mpread,PBB_n]                       %S21min as fitted from S21, double
% KIDparam(nKID).Qi[mpread,PBB_n]                           %Internal Q as fitted from S21, double
% KIDparam(nKID).Qc[mpread,PBB_n]                           %Coupling Q as fitted from S21, double
% KIDparam(nKID).fres[mpread,PBB_n]                         %Resonant frequency as fitted from S21, double
% KIDparam(nKID).Pint[mpread,PBB_n]                         %Internal power, double
% KIDparam(nKID).f_noise{mpread,PBB_n}                      %PSD frequency, comes with phasenoise and ampnoise, cell
% KIDparam(nKID).phasenoise{mpread,PBB_n}                   %Phase PSD, cell
% KIDparam(nKID).ampnoise{mpread,PBB_n}                     %Amplitude PSD, cell
% KIDparam(nKID).Pbbnoise[mpread,PBB_n]                     %Blackbody power where noise is measured, double
% KIDparam(nKID).Tbbnoise[mpread,PBB_n]                     %Blackbody temperature where noise is measured, double
% KIDparam(nKID).wave[mpread,PBB_n]                         %NEP-wave (photon-bunching term), determined by blackbody_4, double
% KIDparam(nKID).g_r[mpread,PBB_n]                          %NEP-gr, (recombination term) determined by blackbody_4, double
% KIDparam(nKID).poisson[mpread,PBB_n]                      %NEP-poisson (photon-poisson term), determined by blackbody_4, double
% KIDparam(nKID).totphoton[mpread,PBB_n]                    %Total photon noise, sum of NEP-wave, NEP-gr, NEP-poisson, double
% KIDparam(nKID).df{mpread,PBB_n}                           %phase response devided by 4Q, cell
% KIDparam(nKID).tofitR{mpread,PBB_n}                       %Range to fit amplitude responsivity, cell
% KIDparam(nKID).tofittheta{mpread,PBB_n}                   %Range to fit phase responsivity, cell
% KIDparam(nKID).dthdP[mpread,PBB_n]                        %Phase responsivity from fit, double
% KIDparam(nKID).dthdPstd[mpread,PBB_n]                     %Std in Phase responsivity from fit, double
% KIDparam(nKID).dRdP[mpread,PBB_n]                         %Amplitude responsivity from fit, double
% KIDparam(nKID).dRdPstd[mpread,PBB_n]                      %Std in Amplitude responsivity from fit, double
% KIDparam(nKID).ddxdP[mpread,PBB_n]                        %X-responsivity, ie fit to df=phase/4Q, double
% KIDparam(nKID).ddx2dP[mpread,PBB_n]                       %X2-responsivity, ie fit to amplitude/4Q, double
% KIDparam(nKID).respfitR{mpread,PBB_n}                     %The amplitude responsivity fit (to plot), cell
% KIDparam(nKID).respfitth{mpread,PBB_n}                    %The phase responsivity fit (to plot), cell
% KIDparam(nKID).KIDid[mpread,PBB_n]                        %KID number as given by Labview, double
% KIDparam(nKID).PhaseoptNEP{mpread,PBB_n}                  %Phase NEP spectrum, cell
% KIDparam(nKID).RadiusoptNEP{mpread,PBB_n}                 %Amplitude NEP spectrum, cell
% KIDparam(nKID).phaseNEPfref[mpread,PBB_n]                 %Phase NEP at fref, double
% KIDparam(nKID).radiusNEPfref[mpread,PBB_n]                %Amplitude NEP at fref, double
% KIDparam(nKID).phaserootNoiselevelfref[mpread,PBB_n]      %sqrt(phase_noise_level_@fref), used to calulated the std in for the NEP at fref, double
% KIDparam(nKID).radiusrootNoiselevelfref[mpread,PBB_n]     %sqrt(amplitude_noise_level_@fref), used to calulated the std in for the NEP at fref, double
% KIDparam(nKID).stdphaserootNoiselevelfref[mpread,PBB_n]   %std of phaserootNoiselevelfref, double
% KIDparam(nKID).stdradiusrootNoiselevelfref[mpread,PBB_n]  %std of radiusrootNoiselevelfref, double
% KIDparam(nKID).stdphaseNEPfref[mpread,PBB_n]              %std in phase NEP, double
% KIDparam(nKID).stdradiusNEPfref[mpread,PBB_n]             %std in amplitude NEP, double
% KIDparam(nKID).radiusNEPamplifier[mpread,PBB_n]           %amplifier NEP determined from tail of amplitude PSD, double
% KIDparam(nKID).radiusNEPPopt[PBB_n]                       %amplitude NEP at optimum readout power, 1D double
% KIDparam(nKID).Poptindex[PBB_n]                           %index of optimum readout power, 1D double
% KIDparam(nKID).phaseNEPPopt[PBB_n]                        %phase NEP at optimum readout power, 1D double
% KIDparam(nKID).Popt[PBB_n]                                %optimum readout power, 1D double
% KIDparam(nKID).preadcoefficientR[PBB_n]                   %(optional) coefficient of fit to the readou-power dependent amplitude responsivity, 1D double
% KIDparam(nKID).preadcoefficientth[PBB_n]                  %(optional) coefficient of fit to the readou-power dependent phase responsivity, 1D double
% KIDparam(nKID).stdpreadcoefficientR[PBB_n]                %(optional) std of preadcoefficientR, 1D double
% KIDparam(nKID).stdpreadcoefficientth[PBB_n]               %(optional) std of preadcoefficientth, 1D double
        
end %of the main function

function[powersum,NEP,BBcal]=blackbody_int(TTT,BBcal,pol,method)
% calls external function blackbody_2.m if BBcal.T=0 (first KID in main
% program). backbody.m Calculates the powers vs black body T at many
% temperatures and puts the results in BBcal struct.
% afterwards, and also if BBcal.T~=0 the script inperpoaltes the BBcal
% struct ti get parameters.

% TTT col of temperatures
% BBcal struct with tempertaure, power and photon noise NEP
% pol polarisation 1 or 2
% method: string describing the optics of setup (see blackbody.m)
if isequal(BBcal.T,0) % if BB.cal is empty call blackbody.m
    BBcal.T=[3:0.05:100]';
    [BBcal.power,~,~,BBcal.NEP]=blackbody_8(BBcal.T,pol,method,1); % old slow script
end

% interpolate BBcal struct, log space
powersum=10.^(interp1(log10(BBcal.T),log10(BBcal.power),log10(TTT)));

NEP.poisson=10.^(interp1(log10(BBcal.T),log10(BBcal.NEP.poisson),log10(TTT)));%2phF
NEP.g_r=10.^(interp1(log10(BBcal.T),log10(BBcal.NEP.g_r),log10(TTT)));%2PDelta/eta
NEP.wave=10.^(interp1(log10(BBcal.T),log10(BBcal.NEP.wave),log10(TTT)));%wave term
NEP.totphoton=10.^(interp1(log10(BBcal.T),log10(BBcal.NEP.totphoton),log10(TTT)));
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
