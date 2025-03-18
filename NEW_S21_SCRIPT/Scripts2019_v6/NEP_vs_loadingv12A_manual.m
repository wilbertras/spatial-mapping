function resppathy=NEP_vs_loadingv12A_manual

% Comments by Juan
%
% This script is identical to NEP_vs_loadingv12A but does the fitting
% manually. It is especially useful when you want to fit the responsivity
% at the lowest powers.
%
% The script will ask you to fill up a value for the 'fitlowpower' (a good
% first guess for the lowest optical powers is ~2e-19). It will fit the
% responsivity automatically and plot it. You will be asked if it is a good
% plot. Answer '0' or '1' accordingly. If you don't like it, a new value for
% fitlowpower needs to be filled. This process repeats itself until a good
% fit is obtained. Once you think that the data has been fitted properly,
% the script will save this fit and move to the next power point. You will
% be asked again for a new 'fitlowpower' value. The process repeats itself
% until all the power points have been fitted.
%


%global resppathy
 
% v12A: updated with new TP calcualtions LT20 and LT013 1mm chips @ 850
% GHz. runs withblackbody_12
% included Deshima

% updated to catch missing pts in FFT spectra (can happen with deglitch
% option)
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
%   method.filter='350 GHz WT1275'  %Cardiff filter
%   method.filter='325 GHz He7'     %He7 system filtering used in A-MKID L band He7 cooler
%   method.filter='1_6THz'          % 1.6 THz SPICA SAFARI stack
%   method.filter='350 GHz 4Filters' %350 Cardiff bandpass with additional LPF
%   method.filter='850 GHz'         %850 GHz filterstack used first with LT010
%
%   method.lensdiameter=1 %lens diameter in millimeters
%   method.lensdiameter=2 %lens diameter in millimeters
%
%   method.aperture='He7'                   %He7 system
%   method.aperture='ADR small'             %small 2 mm diameter aperture in the ADR at 15.05 mm above the chip.
%                                               Using beam pattern calculation values for throughput determination.
%   method.aperture='ADR small Geometrical' %small 2 mm diameter aperture in the ADR at 15.05 mm above the chip.
%                                               Using geometrical calculation for throughput determination.
%   method.aperture='ADR full'              %full aperture in the ADR, 0.62 truncation (H10) @350 GHz%
%   method.aperture='ADR full 350 GHz'           %full aperture without any truncation
%
%   method.aperture='ADR full 1.6 THz'          %full aperture ADR with 1.6 Thz SAFARI stack [only for 2 mm lens] (JB,PdV 21-9-2012)
%   method.aperture='ADR full 850 GHz'          %full aperture ADR with 850 GHz Filter stack [only for 2 mm lens] (RJ 20-2-2014)

% (method.aperture,'1mmlens_3mm_15mm_850_GHz') %1 mm lens, 3 mm aperture, 15.04 mm distance, 850 GHz: old LT holder 20x12 used for LT013.
%       %berekening jochem jan 2015%  TP=1.7301e-08;

% (method.aperture,'1mmlens_3mm_17mm_850_GHz') %1 mm lens, 3 mm aperture, 17 mm distance, 850 GHz: beam pattern method.
%       % berekening jochem dec 2014 TP=1.4462e-08;

% strcmp(method.aperture,'1mmlens_2mm_17mm_850_GHz') %1 mm lens, 2 mm aperture, 17 mm distance, 850 GHz: Spacekids holder.
%       %berekening jochem jan 2015%    TP=6.7185e-09;

% (method.aperture,'1mmlens_2mm_15mm_850_GHz') %1 mm lens, 2 mm aperture, 15.7 mm distance, 850 GHz: Spacekids holder.
        %berekening jochem jan 2015%     TP=7.8363e-09 ; 
%(method.aperture,'2mmlens_3mm_17mm_850_GHz') %2 mm lens, 3 mm aperture, 17 mm distance, 850 GHz: Spacekids holder.
        %berekening jochem may 2015%     TP=4.2651e-08 ; 
        


%   method.freqresolution = 2 %[Optional, Default = 2] frequency resolution in GHz used
%                               to integrate the blackbody spectrum with all filters.
%                               NOTE: for the 350 GHz and 325 GHz filters 2 GHz is
%                                       adviced as minimum.
%   method.Delta = 45.6 %[Optional, Default = 45.6] value of Delta (half
%                        the superconducting gap) in GHz
%   method.eta_pb = 0.57 %[Optional, Default = 0.57] value of the photon
%                        pair breaking efficiency. Default from Kozorezov
%                        et.al. PRB 2000
%   method.EpsTrans = 1.0 %[Optional, Default = 1.00] fraction of the
%                         received photon power with respect to co-aligned 
%                         source and lens optical axis
%   NOTE: 2 mm lenses only implemented for ADR small and ADR small Geometrical


%%%%%%%%%%%%%%%%%%%%%%%%%Parameters to Set%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%==========================================================================
% Setting some routine default values
%==========================================================================
                                         %for poor Mac users use '/' 

path=[cd '/../']; %root path where data is, one higher than the scripts
resperpath=[filesep '2D_BB' filesep '2D_BB' filesep];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TDspectra = 0;          % option to use PSD's from TDanalysis (1), you have to run TDanalysis before! If 0 it will use the *FFT.dat files as usual
if TDspectra ==1
    TDpath = [path filesep 'TD_optNEP2D_BB' filesep]; %path where the bin-files are and the result-matlab-workspace from TDanalysis
    load ([TDpath 'TDresults_5.mat'], 'TDparam','FITparam'); %load TD results into workspace, make sure you take the right version of your results
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%NB call global resppathy in main%%%
closefigures=0;         %close figures after saving
readallKIDs=0;          % 0 does not read all KIDs, then KIDstr=[12 14] or so should give the KIDs to read. 
%NB: use readallKIDs=1 before proceeding to the next scripts
KIDstr=[3];         % [3]if not all KIDs are to be read these ones are
fref=100;                % ref post detection frequency to calulate NEP. 
frefSETUPin=0; %default 100.000! If sett o 0 the program, takes the minimum Amplitude and phasenoise value below 300.000 Hz
methodPopt='phase';    % determine optimum readout power from 'radius' NEP or 'phase' NEP
fitlowpower = 2e-19;      % low power limit, can be adjusted to give good fits of responsivity
                        % if background < this value resp fit from 0-5x this value%
                        % if background >100x this value will fit over 90-110% background power%
                        % in between it will fit from 50% to 2x the background power%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

method.filter='1_6THz';
method.tp='lambda^2' ;
method.eta_c = 0.13 ;%total optical efficiency
method.pol = 1;
method.eta_pb = 0.57;
method.EpsTrans = 1.0; %optical efficiency due to translation wrt pinhole
%Known EpsTrans values
%First measurement H20 chip, KID16 0.79\pm 0.05
%First measurement H20 chip, KID3 0.71\pm 0.05

kb = 1.38e-23; %J/K
h = 6.626e-34; %J/Hz
Tc = 1.283; %K
method.Delta = 1.76*kb*Tc/h/1e9; %Half gap in GHz. 



PdVtoys=0; %some optional stuff for further analysis (under construction, PdV) %now: fit responsivity vs Pread and draw GRnoiseline
if PdVtoys==1
    %some bla for GR to photon NEP crossover (under construction, PdV)
    lifefactor=440e-9*190*1.72e10/1.602e-19*(1.38e-23*1.283)^3/2/(1.76*1.38e-23*1.283)^2;
    taumax=150e-6; %the 2.5 is just bla to get it matched
    Nqpmax=lifefactor/taumax;
end

%==========================================================================
% Settings for range over which response fit is done
%==========================================================================
ResponseRangeType = 6;%5, 6 or 3. Best is 6
%Value determines the method used to determine the range over which the
%(amplitude) response is fitted. The range is always centered (in optical
%power) at the optical power at which the noise spectrum (and S21 circle)
%are measured. The options for ResponseRangeType are:
%   1 = (RJ) Fixed optical power range around the base optical power in fW.
%   2 = (PdV/RJ) Fixed percentage of the base optical power
%   3 = (PdV) At low powers where GR-noise limited a percentage of the base
%       power / low power limit, and at high powers a power dependent
%       (empirical equation) fraction of the base optical power is taken.
%   4 = (RJ) as 1, but with a catch is Pbb sweep does not
%       contain Pbb noise (due to for example thermal drift). It takes
%       reference power at the power nearest to Pbbnoise.
%   5 = (RJ) as 2, but with a catch is Pbb sweep does not
%       contain Pbb noise (due to for example thermal drift). It takes
%       reference power at the power nearest to Pbbnoise.
%   6 = (PdV/RJ/JB) Idenical to 3, but if a hysteresis  is present at the
%   range selected by option 3 it will fit the only incrasing popwer part of the response.            

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

colorsPBB = colormap(jet(noBBTS));   %colorsPBB is for BB power
BBcal.T=0; 
kidsP=zeros(maxnumkids,maxnopowers,noBBTS);
Temperaturenoisefile=zeros(noBBTS,maxnumkids,maxnopowers);

    
for PBB_n=1:noBBTS                                  %PBB_n runs over all PBB, i.e. all dirs. This loop is done dir for dir
    resppath=[path resperpath num2str(dirnum(PBB_n)) filesep];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %%%%%%%%%%%%%%     Finding response files      %%%%%%%%%%%%%% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
    temprespfiles=dir([resppath '*ReImT*.dat']);
    if isempty(temprespfiles)
        error(['No files in: ' resppath])
    end
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
        if nopread == 1
            cPr=[ 0 0 1];   %blue 
        else
            cPr=colormap(jet(nopread));
        end
        for m=1:nopread                                     %m=every readout power for that KID at that PBB
 
            clear tempqrt S21Re S21Im response S21F S21dB S21rad data;
            tempqrt=zeros(maxpts,22);

            %%%%%%%%%%%   Reading S21 file GHz Re Im  that was measured together with the response file  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            filetjes=dir([resppath 'KID' num2str(respkids(nKID)) '_' num2str(kidsP(nKID,m,PBB_n)) '*_S21.dat']);%
            if length(filetjes)==1
                S21file=[resppath filetjes.name];
            elseif length(filetjes)>1
                error('identical S21 files found')
            else
                error(['KID' num2str(respkids(nKID)) '_' num2str(kidsP(nKID,m,PBB_n)) '*_S21.dat DOes not exist' ])
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

            % Call new bb script
            [meh,~,BBcal,method]=blackbody_int(KIDparam(nKID).Tbb{m,PBB_n},BBcal,method);
            KIDparam(nKID).Pbb{m,PBB_n}=meh;                    %done this way because of cell bla
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

            noisefilename=dir([resppath 'KID' num2str(respkids(nKID)) '_' num2str(kidsP(nKID,m,PBB_n)) 'dBm_*_FFT.dat']);
            disp(noisefilename.name);
            fid=fopen([resppath noisefilename.name]);
            if fid==-1%%check for exists
                disp('no noise datafile')
            end
            textscan(fid,'%*s%*s%*4c%f%*s',1,'headerlines',1);          %NOT getting power from header
            %temp=textscan(fid,'%*s%*s%*4c%f%*s',1,'headerlines',1);    %getting power from header
            %Powernoisefile=-1*temp{1};                                 %power is correct (negative) whereas kidsP{nKID}(:) is positve since taken from filename
            textscan(fid,'%*[^\r\n]',3);                                %scanning more empty lines to skip rest of header
            temp=textscan(fid,'%*s%*s%*2c%f%*s',1);
            Temperaturenoisefile(PBB_n,nKID,m)=temp{1};                 %ie get the blackbody temperature

            textscan(fid,'%*[^\r\n]',1);
            bla=cell2mat(textscan(fid, '%f%f%f'));
            fclose(fid);
            
            if TDspectra==1 % take the spectra out of the TDparam struct
                searchtemp=round(1000*Temperaturenoisefile(PBB_n,nKID,m)); %BB temp in mK rounded to integer (as it appears in binfile filename)
                if respkids(nKID)~=TDparam(nKID).kidnr; error('KIDnr of TDresults does not match'); end
                TDpreadindex=find(TDparam(nKID).Pread(:,1)==kidsP(nKID,m,PBB_n));
                TDtempindex=find(TDparam(nKID).Temp(TDpreadindex,:)==searchtemp);
                KIDparam(nKID).f_noise{m,PBB_n}=TDparam(nKID).fmtotal{TDpreadindex,TDtempindex};
                KIDparam(nKID).phasenoise{m,PBB_n}=TDparam(nKID).SPPtotalraw{TDpreadindex,TDtempindex}; %in dB
                KIDparam(nKID).ampnoise{m,PBB_n}=TDparam(nKID).SRRtotalraw{TDpreadindex,TDtempindex};   %in dB
            else %take the spectra from the FFT.dat file
                KIDparam(nKID).f_noise{m,PBB_n}=bla(:,1);
                KIDparam(nKID).phasenoise{m,PBB_n}=bla(:,2);                %in dB
                KIDparam(nKID).ampnoise{m,PBB_n}=bla(:,3);                  %in dB
            end
            
            %%%%%%%%%%%%%%%%%%%%%% call new BB script to find the background power %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [KIDparam(nKID).Pbbnoise(m,PBB_n),poo,BBcal]=blackbody_int(Temperaturenoisefile(PBB_n,nKID,m),BBcal,method);
            KIDparam(nKID).Tbbnoise(m,PBB_n)=Temperaturenoisefile(PBB_n,nKID,m);
            KIDparam(nKID).wave(m,PBB_n)=poo.wave;
            KIDparam(nKID).g_r(m,PBB_n)=poo.g_r;
            KIDparam(nKID).poisson(m,PBB_n)=poo.poisson;
            KIDparam(nKID).totphoton(m,PBB_n)=poo.totphoton;
            clear meh poo
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%Processing data: finding the phase radius and dx response  %%%%%%%%%%%%%%
            KIDparam(nKID).df{m,PBB_n}=-1*KIDparam(nKID).Phaseresp{m,PBB_n}/(4*KIDparam(nKID).Q(m,PBB_n)); %estimation of F response
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%now finding response at base =  used for the NEP calcul%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            
% It fits manually the responsivity for each blackbody power

            fitagain = 0

            while fitagain == 0
            
                prompt = (['Change the low power limit for the responsivity fit. Old value equals: ', num2str(fitlowpower),'\n Enter new value: ']);
                fitlowpower = input(prompt);
                
                %This piece of code determines the blackbody power range over
                %which the response is fitted. This generally needs to be
                %adjusted for each sample type and observed effects. As a
                %reference the base optical power is taken, which is the
                %optical power at which the noise spectrum and S21 circle are
                %measured. This value is stored in KIDparam(nKID).Pbbnoise(m,PBB_n)
                %ResponseRangeType:
                if ResponseRangeType == 1
                    if (KIDparam(nKID).Pbbnoise(m,PBB_n) < min(KIDparam(nKID).Pbb{m,PBB_n})) || (KIDparam(nKID).Pbbnoise(m,PBB_n) > max(KIDparam(nKID).Pbb{m,PBB_n}))
                        disp('WARNING: Pbb @ Noise outside response sweep.')
                        disp('Investigate data and try a different ResponseRangeType')
                    end
                    fitlowbound = KIDparam(nKID).Pbbnoise(m,PBB_n) - Prange;
                    fithighbound = KIDparam(nKID).Pbbnoise(m,PBB_n) + Prange;
                elseif ResponseRangeType == 2
                    fitlowbound = (1-Prange)*KIDparam(nKID).Pbbnoise(m,PBB_n);
                    fithighbound = (1+Prange)*KIDparam(nKID).Pbbnoise(m,PBB_n);
                elseif ResponseRangeType == 3 || ResponseRangeType == 6 %Both same ranges
                    if KIDparam(nKID).Pbbnoise(m,PBB_n) < 10*fitlowpower  
                        fitlowbound=max([fitlowpower 1.2*KIDparam(nKID).Pbbnoise(m,PBB_n)]); %the 1.2 is to not fit the noise at PBBnoise
                        fithighbound=18*max([fitlowpower KIDparam(nKID).Pbbnoise(m,PBB_n)]);
                        %             elseif (KIDparam(nKID).Pbbnoise(m,PBB_n) > 10*fitlowpower)&&(KIDparam(nKID).Pbbnoise(m,PBB_n) < 2000*fitlowpower)
                        %                 fitlowbound=KIDparam(nKID).Pbbnoise(m,PBB_n)/2;
                        %                 fithighbound=KIDparam(nKID).Pbbnoise(m,PBB_n)*1.5;
                        %             else
                        %                 fitlowbound=.9*KIDparam(nKID).Pbbnoise(m,PBB_n);
                        %                 fithighbound=1.1*KIDparam(nKID).Pbbnoise(m,PBB_n);
                        %             end
                    else
                        factorbound=1/(KIDparam(nKID).Pbbnoise(m,PBB_n)/fitlowpower)^(1/4);
                        if factorbound>=1; factorbound=1e-1; end
                        fitlowbound=(1-factorbound)*KIDparam(nKID).Pbbnoise(m,PBB_n);
                        fithighbound=(1+factorbound)*KIDparam(nKID).Pbbnoise(m,PBB_n);
                    end
                elseif ResponseRangeType == 4
                    if KIDparam(nKID).Pbbnoise(m,PBB_n) < min(KIDparam(nKID).Pbb{m,PBB_n})
                        %due to thermal drift it can happen that the power at
                        %which the noise is measured is not within the power
                        %range in this case the nearest (lowest) power in the 
                        %sweep is taken as reference power.
                        disp('WARNING: Pbb @ Noise outside (and below) response sweep')
                        fitlowbound = min(KIDparam(nKID).Pbb{m,PBB_n});
                        fithighbound = fitlowbound + Prange;
                    elseif KIDparam(nKID).Pbbnoise(m,PBB_n) > max(KIDparam(nKID).Pbb{m,PBB_n})
                        %due to thermal drift it can happen that the power at
                        %which the noise is measured is not within the power
                        %range in this case the nearest (highest) power in the 
                        %sweep is taken as reference power.
                        disp('WARNING: Pbb @ Noise outside (and above) response sweep')
                        fithighbound = max(KIDparam(nKID).Pbb{m,PBB_n});
                        fitlowbound = fithighbound - Prange;
                    else
                        %regular range as per scenario 1.
                        fitlowbound = KIDparam(nKID).Pbbnoise(m,PBB_n) - Prange;
                        fithighbound = KIDparam(nKID).Pbbnoise(m,PBB_n) + Prange;
                    end
                elseif ResponseRangeType == 5
                    if KIDparam(nKID).Pbbnoise(m,PBB_n) < min(KIDparam(nKID).Pbb{m,PBB_n})
                        %due to thermal drift it can happen that the power at
                        %which the noise is measured is not within the power
                        %range in this case the nearest (lowest) power in the
                        %sweep is taken as reference power.
                        disp('WARNING: Pbb @ Noise outside (and below) response sweep')
                        fitlowbound = min(KIDparam(nKID).Pbb{m,PBB_n});
                        fithighbound = fitlowbound*(1 + Prange);
                    elseif KIDparam(nKID).Pbbnoise(m,PBB_n) > max(KIDparam(nKID).Pbb{m,PBB_n})
                        %due to thermal drift it can happen that the power at
                        %which the noise is measured is not within the power
                        %range in this case the nearest (highest) power in the
                        %sweep is taken as reference power.
                        disp('WARNING: Pbb @ Noise outside (and above) response sweep')
                        fithighbound = max(KIDparam(nKID).Pbb{m,PBB_n});
                        fitlowbound = fithighbound*(1 - Prange);
                    else
                        %regular range as per scenario 1.
                        fitlowbound = KIDparam(nKID).Pbbnoise(m,PBB_n)*(1 - Prange);
                        fithighbound = KIDparam(nKID).Pbbnoise(m,PBB_n)*(1 + Prange);
                    end

                else
                    error('NEPvsPbb: Unknown ResponseRangeType. This must be a value from [1,2,3,4,5]')
                end

                %Use the calculated upper and lower bounds in power to
                %determine which Pbb of the response measurement fall within
                %this range.
                KIDparam(nKID).tofitR{m,PBB_n} = (KIDparam(nKID).Pbb{m,PBB_n}>fitlowbound) & (KIDparam(nKID).Pbb{m,PBB_n}<fithighbound) & (KIDparam(nKID).Phaseresp{m,PBB_n}<1);
                KIDparam(nKID).tofittheta{m,PBB_n} = (KIDparam(nKID).Pbb{m,PBB_n}>fitlowbound/notheta) & (KIDparam(nKID).Pbb{m,PBB_n}<notheta*fithighbound) & (KIDparam(nKID).Phaseresp{m,PBB_n}<1);

                if ResponseRangeType == 6 %we now want to fit onl;y up and down in P of BB
                    dP=[0 ; sign(KIDparam(nKID).Pbb{m,PBB_n}(2:end)-KIDparam(nKID).Pbb{m,PBB_n}(1:end-1))];
                    KIDparam(nKID).tofitR{m,PBB_n} = (KIDparam(nKID).Pbb{m,PBB_n}>fitlowbound) & (KIDparam(nKID).Pbb{m,PBB_n}<fithighbound) & ...
                        (KIDparam(nKID).Phaseresp{m,PBB_n}<1 & dP>0);
                    KIDparam(nKID).tofittheta{m,PBB_n} = (KIDparam(nKID).Pbb{m,PBB_n}>fitlowbound/notheta) & (KIDparam(nKID).Pbb{m,PBB_n}<notheta*fithighbound) & ...
                        (KIDparam(nKID).Phaseresp{m,PBB_n}<1 & dP>0);

                end
                %%%%%%%%%%%%%%%%%%%%%%%% Finding responsivity by local lineair fit to the data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                %radius 
                if nnz(KIDparam(nKID).tofitR{m,PBB_n})>=3
                    %Minimum of 3 data points is required to be able to do a
                    %constrained linear fit for which the confidence intervals
                    %can be calculated.
                    %This catches the situations in which either the fit range
                    %is ill defined or the resonator is overdriven and phase>1.
                    if nnz(KIDparam(nKID).tofitR{m,PBB_n})<=10
                        disp(['WARNING: Fitting radius response to ',num2str(nnz(KIDparam(nKID).tofitR{m,PBB_n})),' datapoints.\n'])
                        disp(['Check KID ',num2str(respkids(nKID)),' at ',num2str(KIDparam(nKID).Pread(m,PBB_n)),'dBm'])
                        disp(['Main problems are typically: ill defined power range, thermal drift or overdriving.'])
                    end

                    s = fitoptions('Method','LinearLeastSquares'); 
                    ftype = fittype({'x',num2str(1)},'options', s); % performs a linear 'a.*x+b' fit in fW (reduces errors). corrected later 
                    [radiusresult]=fit(KIDparam(nKID).Pbb{m,PBB_n}(KIDparam(nKID).tofitR{m,PBB_n})*1E15,KIDparam(nKID).Radiusresp{m,PBB_n}(KIDparam(nKID).tofitR{m,PBB_n})-1,ftype);
                    dRdPfit=radiusresult.a; 
                    astart=radiusresult.b; 
                    result2=confint(radiusresult);                  %confint extracts the 95% confidence intervals from result 
                    stddRdPfit=(dRdPfit-result2(1,1))/1.96;         %convert 95% interval to standard deviation 
                else 
                    dRdPfit=1e-20; %very small number to make the NEP for this power useless (very high)
                    astart=0; 
                    stddRdPfit=1e-10; 
                end 
                %phase 
                if nnz(KIDparam(nKID).tofittheta{m,PBB_n})>=3
                    %Minimum of 3 data points is required to be able to do a
                    %constrained linear fit for which the confidence intervals
                    %can be calculated.
                    %This catches the situations in which either the fit range
                    %is ill defined or the resonator is overdriven and phase>1.
                    if nnz(KIDparam(nKID).tofittheta{m,PBB_n})<=10
                        disp(['WARNING: Fitting phase response to ',num2str(nnz(KIDparam(nKID).tofitR{m,PBB_n})),' datapoints.\n'])
                        disp(['Check KID ',num2str(respkids(nKID)),' at ',num2str(KIDparam(nKID).Pread(m,PBB_n)),'dBm'])
                        disp(['Main problems are typically: ill defined power range, thermal drift or overdriving.'])
                    end

                    [phaseresult]=fit(KIDparam(nKID).Pbb{m,PBB_n}(KIDparam(nKID).tofittheta{m,PBB_n})*1E15,KIDparam(nKID).Phaseresp{m,PBB_n}(KIDparam(nKID).tofittheta{m,PBB_n}),ftype);
                    dthetadPfit=phaseresult.a;
                    pstart=phaseresult.b;
                    result2=confint(phaseresult);                   %confint extracts the 95% confidence intervals from result
                    stddthetadPfit=(dthetadPfit-result2(1,1))/1.96; %convert 95% interval to standard deviation 
                else 
                    dthetadPfit=1e-20; 
                    pstart=0; 
                    stddthetadPfit=1e-10; 
                end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                KIDparam(nKID).dthdP(m,PBB_n)=dthetadPfit(1)*1e15;  %response corrected back to W
                KIDparam(nKID).dthdPoffset(m,PBB_n)=pstart;
                KIDparam(nKID).dthdPstd(m,PBB_n)=stddthetadPfit*1e15;
                KIDparam(nKID).dRdP(m,PBB_n)=dRdPfit(1)*1e15;
                KIDparam(nKID).dRdPoffset(m,PBB_n)=astart;
                KIDparam(nKID).dRdPstd(m,PBB_n)=stddRdPfit*1e15;
                KIDparam(nKID).ddxdP(m,PBB_n)=KIDparam(nKID).dthdP(m,PBB_n)./KIDparam(nKID).Q(m,PBB_n)/4; %% ddx/dP
                KIDparam(nKID).ddx2dP(m,PBB_n)=-KIDparam(nKID).dRdP(m,PBB_n)./KIDparam(nKID).Q(m,PBB_n)/4; %% ddx2/dP, ie radius response devided by Q
                KIDparam(nKID).respfitR{m,PBB_n}=astart+KIDparam(nKID).Pbb{m,PBB_n}*KIDparam(nKID).dRdP(m,PBB_n);
                KIDparam(nKID).respfitth{m,PBB_n}=pstart+KIDparam(nKID).Pbb{m,PBB_n}*KIDparam(nKID).dthdP(m,PBB_n);
            
            
            
                for mpread=1:length(nonzeros(kidsP(nKID,:,PBB_n)))
                    figure(111)
                    semilogx(KIDparam(nKID).Pbb{mpread,PBB_n},KIDparam(nKID).Phaseresp{mpread,PBB_n},'ob');
                    hold on;
%                     plot(KIDparam(nKID).Pbb{mpread,PBB_n},KIDparam(nKID).Radiusresp{mpread,PBB_n}-1,'or')
                    title('Response + fit')
                    legend('Phase','Radius-1')
                    ylabel('Reponse')
                    xlabel('Loading (W)')
                    semilogx(KIDparam(nKID).Pbb{mpread,PBB_n},KIDparam(nKID).respfitth{mpread,PBB_n},'-r','Linewidth',4);% fits
%                     plot(KIDparam(nKID).Pbb{mpread,PBB_n},KIDparam(nKID).respfitR{mpread,PBB_n},'-r');
                    hold off;
                end
                
                prompt = ('Are you happy with the fit? \n If yes, type 1.  If not, type 0. \n');
                fitagain = input(prompt);
                
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%% getting optical NEP data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            KIDparam(nKID).KIDid(m,PBB_n)=respkids(nKID);
            KIDparam(nKID).PhaseoptNEP{m,PBB_n}=10.^(KIDparam(nKID).phasenoise{m,PBB_n}/20)/KIDparam(nKID).dthdP(m,PBB_n);%phase NEP
            KIDparam(nKID).RadiusoptNEP{m,PBB_n}=10.^(KIDparam(nKID).ampnoise{m,PBB_n}/20)/(-KIDparam(nKID).dRdP(m,PBB_n));%R nep
            indfref=find(KIDparam(nKID).f_noise{m,PBB_n}>fref,1);
            blawindow=2;
            if ~isempty(indfref) && indfref-blawindow>0
                KIDparam(nKID).phaseNEPfref(m,PBB_n)=mean(KIDparam(nKID).PhaseoptNEP{m,PBB_n}((indfref-blawindow):(indfref+blawindow))); % ref Hz value
                KIDparam(nKID).radiusNEPfref(m,PBB_n)=mean(KIDparam(nKID).RadiusoptNEP{m,PBB_n}((indfref-blawindow):(indfref+blawindow))); % ref Hz value
            else
                KIDparam(nKID).phaseNEPfref(m,PBB_n)=NaN; % ref Hz value
                KIDparam(nKID).radiusNEPfref(m,PBB_n)=NaN; % ref Hz value
            end
            %%%%%%%%%%%%%%%%%%%%%%%% calculated error bar on optical NEP from std_response and std_noiselevel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %STD of NEP has to be calculated by combining std of noise level and std of reponsivity. If you only take the std of the
            %NEP over 'blawindow' you neglect the uncertainty in responsivity
            %%%rootNoiselevel is chosen because it goes directly into NEP
            if ~isempty(indfref) && indfref-blawindow>0 %catching issue with the FFT not beign complete sometimes
                KIDparam(nKID).phaserootNoiselevelfref(m,PBB_n)= mean(10.^(KIDparam(nKID).phasenoise{m,PBB_n}((indfref-blawindow):(indfref+blawindow))/20));
                KIDparam(nKID).radiusrootNoiselevelfref(m,PBB_n)= mean(10.^(KIDparam(nKID).ampnoise{m,PBB_n}((indfref-blawindow):(indfref+blawindow))/20));
                KIDparam(nKID).stdphaserootNoiselevelfref(m,PBB_n)= std(10.^(KIDparam(nKID).phasenoise{m,PBB_n}((indfref-blawindow):(indfref+blawindow))/20));
                KIDparam(nKID).stdradiusrootNoiselevelfref(m,PBB_n)= std(10.^(KIDparam(nKID).ampnoise{m,PBB_n}((indfref-blawindow):(indfref+blawindow))/20));
                KIDparam(nKID).stdphaseNEPfref(m,PBB_n)=sqrt((KIDparam(nKID).dthdPstd(m,PBB_n)/KIDparam(nKID).dthdP(m,PBB_n))^2+(KIDparam(nKID).stdphaserootNoiselevelfref(m,PBB_n)/KIDparam(nKID).phaserootNoiselevelfref(m,PBB_n))^2)*KIDparam(nKID).phaseNEPfref(m,PBB_n); % 20Hz value
                KIDparam(nKID).stdradiusNEPfref(m,PBB_n)=sqrt((KIDparam(nKID).dRdPstd(m,PBB_n)/KIDparam(nKID).dRdP(m,PBB_n))^2+(KIDparam(nKID).stdradiusrootNoiselevelfref(m,PBB_n)/KIDparam(nKID).radiusrootNoiselevelfref(m,PBB_n))^2)*KIDparam(nKID).radiusNEPfref(m,PBB_n); % 20Hz value
            else
                KIDparam(nKID).phaserootNoiselevelfref(m,PBB_n)= NaN;
                KIDparam(nKID).radiusrootNoiselevelfref(m,PBB_n)= NaN;
                KIDparam(nKID).stdphaserootNoiselevelfref(m,PBB_n)= NaN;
                KIDparam(nKID).stdradiusrootNoiselevelfref(m,PBB_n)= NaN;
                KIDparam(nKID).stdphaseNEPfref(m,PBB_n)=NaN; % 20Hz value
                KIDparam(nKID).stdradiusNEPfref(m,PBB_n)=NaN; % 20Hz value
                
            end
            %%%%%%%%%%%%%%%%% calculate amplifier limited optNEP
            if frefSETUPin ~=0
                highfind=find(KIDparam(nKID).f_noise{m,PBB_n}>frefSETUPin,1);
                KIDparam(nKID).frefSETUPphase(m,PBB_n)=frefSETUPin;
                KIDparam(nKID).frefSETUPradius(m,PBB_n)=frefSETUPin;
            else 
                rtuF=find(KIDparam(nKID).f_noise{m,PBB_n}<3e5 & KIDparam(nKID).f_noise{m,PBB_n} > 1000);
                [~,thf]=min(KIDparam(nKID).ampnoise{m,PBB_n}(rtuF));%frefSETUPin=0, we take the miniumu value up to 3e5Hz
                highfind=rtuF(thf);
                KIDparam(nKID).frefSETUPradius(m,PBB_n)=KIDparam(nKID).f_noise{m,PBB_n}(highfind);KIDparam(nKID).f_noise{m,PBB_n}(highfind)
                %phase
                [~,thf2]=min(KIDparam(nKID).phasenoise{m,PBB_n}(rtuF));%frefSETUPin=0, we take the miniumu value up to 3e5Hz
                highfind2=rtuF(thf2);
                KIDparam(nKID).frefSETUPphase(m,PBB_n)=KIDparam(nKID).f_noise{m,PBB_n}(highfind2);KIDparam(nKID).f_noise{m,PBB_n}(highfind2)
                %
            end
            if ~isempty(highfind) && highfind-blawindow-blawindow>0
                KIDparam(nKID).radiusNEPamplifier(m,PBB_n)=mean(KIDparam(nKID).RadiusoptNEP{m,PBB_n}((highfind-blawindow):(highfind+blawindow)));
                KIDparam(nKID).radiusnoiseamplifier(m,PBB_n)=mean(KIDparam(nKID).ampnoise{m,PBB_n}((highfind-blawindow):(highfind+blawindow)));
            else
                KIDparam(nKID).radiusNEPamplifier(m,PBBfrefSETUPfrefSETUP_n)=NaN;
            end
            %phase
            if ~isempty(highfind2) && highfind2-blawindow-blawindow>0
                KIDparam(nKID).phaseNEPamplifier(m,PBB_n)=mean(KIDparam(nKID).PhaseoptNEP{m,PBB_n}((highfind2-blawindow):(highfind2+blawindow)));
                KIDparam(nKID).phasenoiseamplifier(m,PBB_n)=mean(KIDparam(nKID).phasenoise{m,PBB_n}((highfind2-blawindow):(highfind2+blawindow)));   
            else
                KIDparam(nKID).phaseNEPamplifier(m,PBB_n)=NaN;
            end
            
            %save response file with Tbb,Pbb,radius resp, phase resp, frequency resp
            respresultsfile = [resppath 'KID' num2str(respkids(nKID)) '_' num2str(kidsP(nKID,m,PBB_n)) 'dBm_Responses.csv'];
            WriteSRONcsv(respresultsfile,[KIDparam(nKID).Tbb{m,PBB_n} KIDparam(nKID).Pbb{m,PBB_n} KIDparam(nKID).Radiusresp{m,PBB_n} KIDparam(nKID).Phaseresp{m,PBB_n} KIDparam(nKID).df{m,PBB_n}],'Tbb, Pbb, radiusresp, phaseresp, freqresp \nKID','%.12e');
    
        end %loop over P readout
        
        
        disp(['*********** index= ' num2str(PBB_n) ', Pbb= ' num2str(KIDparam(nKID).Pbbnoise(m,PBB_n)) ', TBB= ' num2str(KIDparam(nKID).Tbbnoise(m,PBB_n)) '***********Done'])
        
        
            %%%%%%%%%%%%%%%%Calculate Popt for every PBB, ie minimise in NEP
            radtemp=KIDparam(nKID).radiusNEPfref(:,PBB_n);
            if strcmp(methodPopt,'radius')
                radtemp(radtemp<0)=1e99; %to catch negative values, can occur where KID is switching in S21
                [KIDparam(nKID).radiusNEPPopt(PBB_n), KIDparam(nKID).Poptindex(PBB_n)]=min(nonzeros(radtemp));
                KIDparam(nKID).phaseNEPPopt(PBB_n)=KIDparam(nKID).phaseNEPfref(KIDparam(nKID).Poptindex(PBB_n),PBB_n);
                clear radtemp
            elseif strcmp(methodPopt,'phase')
                phasetemp=KIDparam(nKID).phaseNEPfref(:,PBB_n);phasetemp(radtemp<0)=1e99;%to catch negative values, can occur where KID is switching in S21
                [KIDparam(nKID).phaseNEPPopt(PBB_n), KIDparam(nKID).Poptindex(PBB_n)]=min(nonzeros(phasetemp));
                KIDparam(nKID).radiusNEPPopt(PBB_n)=KIDparam(nKID).radiusNEPfref(KIDparam(nKID).Poptindex(PBB_n),PBB_n);
                clear phasetemp
            end
            KIDparam(nKID).Popt(PBB_n)=kidsP(nKID,KIDparam(nKID).Poptindex(PBB_n),PBB_n);
        
            
            
        %%%%%%%%%%%%Optional: Fit a power law to the responsivity versus microwave power
        if PdVtoys==1 %fit responsivity vs Pread, power law coefficient
            if nnz(KIDparam(nKID).Pread(:,PBB_n)) >= 2
            weightsdR=1./nonzeros(KIDparam(nKID).dRdPstd(:,PBB_n));
            s = fitoptions('Method','LinearLeastSquares','Weights',weightsdR);
            ftype = fittype({'x',num2str(1)},'options', s);
            [radiusfitresult]=fit(-nonzeros(kidsP(nKID,:,PBB_n)),log10(nonzeros(-KIDparam(nKID).dRdP(:,PBB_n))),ftype);
            [phasefitresult]=fit(-nonzeros(kidsP(nKID,:,PBB_n)),log10(nonzeros(KIDparam(nKID).dthdP(:,PBB_n))),ftype);
            KIDparam(nKID).preadcoefficientR(PBB_n)=radiusfitresult.a;
            KIDparam(nKID).preadcoefficientth(PBB_n)=phasefitresult.a;
            KIDparam(nKID).preadoffsetR(PBB_n)=radiusfitresult.b;   %needed for plot
            KIDparam(nKID).preadoffsetth(PBB_n)=phasefitresult.b;    %needed for plot
            result2a=confint(radiusfitresult);
            result2p=confint(phasefitresult);
            KIDparam(nKID).stdpreadcoefficientR(PBB_n)=(KIDparam(nKID).preadcoefficientR(PBB_n)-result2a(1,1))/1.96;
            KIDparam(nKID).stdpreadcoefficientth(PBB_n)=(KIDparam(nKID).preadcoefficientth(PBB_n)-result2p(1,1))/1.96;
            end
        end
        
        %radresponse=figure(respkids(nKID)*10000000+PBB_n);
        %%%%%%%%%%%%%%%Make plots that will be generated for every PBB,KID to show readout power dependences%%%%%
        preadfig=figure(respkids(nKID)*10000+PBB_n);
        subplot(3,3,1)
        semilogy(-nonzeros(kidsP(nKID,:,PBB_n)),nonzeros(KIDparam(nKID).phaseNEPfref(:,PBB_n)),'ro',...
                -nonzeros(kidsP(nKID,:,PBB_n)),nonzeros(KIDparam(nKID).radiusNEPfref(:,PBB_n)),'bo',...
                -nonzeros(kidsP(nKID,:,PBB_n)),nonzeros(KIDparam(nKID).phaseNEPamplifier(:,PBB_n)),'rx',...
                -nonzeros(kidsP(nKID,:,PBB_n)),nonzeros(KIDparam(nKID).radiusNEPamplifier(:,PBB_n)),'bx');hold on
            xlabel('P_{read}(dBm)');
            ylabel(['optical NEP at ' num2str(fref) ' Hz (W/Hz^{1/2})']);
            legend('phase','radius','radius amplifier (>350 kHz)','phase amplifier')
            errorbar(-nonzeros(kidsP(nKID,:,PBB_n)),nonzeros(KIDparam(nKID).phaseNEPfref(:,PBB_n)),nonzeros(KIDparam(nKID).stdphaseNEPfref(:,PBB_n)),'ro');
            errorbar(-nonzeros(kidsP(nKID,:,PBB_n)),nonzeros(KIDparam(nKID).radiusNEPfref(:,PBB_n)),nonzeros(KIDparam(nKID).stdradiusNEPfref(:,PBB_n)),'bo');
            title(['KID ' num2str(respkids(nKID)) ', ' num2str(KIDparam(nKID).Pbbnoise(1,PBB_n)*1e15) ' fW, P_{read,opt} = -' num2str(KIDparam(nKID).Popt(PBB_n)) 'dBm']);
        
        subplot(3,3,2)
        semilogy(-nonzeros(kidsP(nKID,:,PBB_n)),nonzeros(KIDparam(nKID).dthdP(:,PBB_n)),'ro',-nonzeros(kidsP(nKID,:,PBB_n)),nonzeros(-KIDparam(nKID).dRdP(:,PBB_n)),'bx');hold on
        errorbar(-nonzeros(kidsP(nKID,:,PBB_n)),nonzeros(KIDparam(nKID).dthdP(:,PBB_n)),nonzeros(KIDparam(nKID).dthdPstd(:,PBB_n)),'ro');
        errorbar(-nonzeros(kidsP(nKID,:,PBB_n)),-nonzeros(KIDparam(nKID).dRdP(:,PBB_n)),nonzeros(KIDparam(nKID).dRdPstd(:,PBB_n)),'bx');
        xlabel('P_{read}(dBm)');ylabel('Responsivity (/W)'); legend('phase','-radius')
        if PdVtoys==1
            if nnz(KIDparam(nKID).Pread(:,PBB_n)) >= 2
            semilogy(-nonzeros(kidsP(nKID,:,PBB_n)),10.^(KIDparam(nKID).preadoffsetR(PBB_n)-nonzeros(kidsP(nKID,:,PBB_n)).*KIDparam(nKID).preadcoefficientR(PBB_n)),'b');
            semilogy(-nonzeros(kidsP(nKID,:,PBB_n)),10.^(KIDparam(nKID).preadoffsetth(PBB_n)-nonzeros(kidsP(nKID,:,PBB_n)).*KIDparam(nKID).preadcoefficientth(PBB_n)),'r');
            if strcmp(methodPopt,'radius')
                title(['fit coefficient R: ' num2str(KIDparam(nKID).preadcoefficientR(PBB_n)*10,3) '\pm' num2str(KIDparam(nKID).stdpreadcoefficientR(PBB_n)*10,2) ])
            elseif strcmp(methodPopt,'phase')
                title(['fit coefficient \theta: ' num2str(KIDparam(nKID).preadcoefficientth(PBB_n)*10,3) '\pm' num2str(KIDparam(nKID).stdpreadcoefficientth(PBB_n)*10,2) ])
            end
            end
        end
            
        subplot(3,3,3)
        semilogy(-nonzeros(kidsP(nKID,:,PBB_n)),nonzeros(KIDparam(nKID).ddxdP(:,PBB_n)),'ro');hold on;
        semilogy(-nonzeros(kidsP(nKID,:,PBB_n)),nonzeros(KIDparam(nKID).ddx2dP(:,PBB_n)),'bx');
        xlabel('P_{read}(dBm)');ylabel('dxdP (/W)');legend('dxdP=d\theta/dP*1/4Q','dx2dP=-dR/dP*1/4Q')
        
        subplot(3,3,4)
        semilogy(-nonzeros(kidsP(nKID,:,PBB_n)),nonzeros(KIDparam(nKID).Q(:,PBB_n)),'s',-nonzeros(kidsP(nKID,:,PBB_n)),nonzeros(KIDparam(nKID).Qi(:,PBB_n)),'d',-nonzeros(kidsP(nKID,:,PBB_n)),nonzeros(KIDparam(nKID).Qc(:,PBB_n)),'^')
        xlabel('P_{read}(dBm)');ylabel('quality factor from S21');legend('Q','Qi','Qc')
        
        for mpread=1:length(nonzeros(kidsP(nKID,:,PBB_n)))
            subplot(3,3,5)
            loglog(KIDparam(nKID).f_noise{mpread,PBB_n},KIDparam(nKID).PhaseoptNEP{mpread,PBB_n},'--','color',cPr(nopread+1-mpread,:))
            xlim([0.5 320000]);grid on;hold on;
            loglog(KIDparam(nKID).f_noise{mpread,PBB_n},KIDparam(nKID).RadiusoptNEP{mpread,PBB_n},'color',cPr(nopread+1-mpread,:))
            ylabel('NEP (W/Hz)');xlabel('Frequency (Hz)');legend('Phase','Radius');
            title(['P_{read,opt} = -' num2str(KIDparam(nKID).Popt(PBB_n)) 'dBm']);
            loglog(fref,KIDparam(nKID).radiusNEPfref(mpread,PBB_n),'o','color',cPr(nopread+1-mpread,:),'MarkerFaceColor',cPr(nopread+1-mpread,:))
            loglog(fref,KIDparam(nKID).phaseNEPfref(mpread,PBB_n),'o','color',cPr(nopread+1-mpread,:))
            loglog(KIDparam(nKID).frefSETUPradius(m,PBB_n),KIDparam(nKID).radiusNEPamplifier(mpread,PBB_n),'s','color',cPr(nopread+1-mpread,:),'MarkerFaceColor',cPr(nopread+1-mpread,:))
            loglog(KIDparam(nKID).frefSETUPphase(m,PBB_n),KIDparam(nKID).phaseNEPamplifier(mpread,PBB_n),'s','color',cPr(nopread+1-mpread,:))
            
            subplot(3,3,6)
            semilogx(KIDparam(nKID).f_noise{mpread,PBB_n},KIDparam(nKID).phasenoise{mpread,PBB_n},'color',cPr(nopread+1-mpread,:))
            xlim([0.5 300000]);grid on;hold on;
            semilogx(KIDparam(nKID).f_noise{mpread,PBB_n},KIDparam(nKID).ampnoise{mpread,PBB_n},'color',cPr(nopread+1-mpread,:))  
            semilogx(KIDparam(nKID).frefSETUPradius(m,PBB_n),KIDparam(nKID).radiusnoiseamplifier(mpread,PBB_n),'s','color',cPr(nopread+1-mpread,:),'MarkerFaceColor',cPr(nopread+1-mpread,:))
            semilogx(KIDparam(nKID).frefSETUPphase(m,PBB_n),KIDparam(nKID).phasenoiseamplifier(mpread,PBB_n),'s','color',cPr(nopread+1-mpread,:))
            ylabel('Noise spectrum (dBc/Hz)');xlabel('Frequency (Hz)');legend('Phase','Radius');title('Noise')
            
            
            subplot(3,3,7)
            plot(KIDparam(nKID).Pbb{mpread,PBB_n}*1E15,KIDparam(nKID).Phaseresp{mpread,PBB_n},'s','color',cPr(nopread+1-mpread,:));hold on;
            plot(KIDparam(nKID).Pbb{mpread,PBB_n}*1E15,KIDparam(nKID).Radiusresp{mpread,PBB_n}-1,'o','color',cPr(nopread+1-mpread,:))
            title('Response + fit')
            ylabel('Reponse')
            xlabel('Loading (fW)')
             %Overplot part of response that was used for the fit (closed symbols)%
            plot(KIDparam(nKID).Pbb{mpread,PBB_n}(KIDparam(nKID).tofitR{mpread,PBB_n})*1E15,... 
                KIDparam(nKID).Phaseresp{mpread,PBB_n}(KIDparam(nKID).tofitR{mpread,PBB_n}),...
                's','MarkerFaceColor',cPr(nopread+1-mpread,:),'color',cPr(nopread+1-mpread,:));hold on;
            plot(KIDparam(nKID).Pbb{mpread,PBB_n}(KIDparam(nKID).tofitR{mpread,PBB_n})*1E15,...
                KIDparam(nKID).Radiusresp{mpread,PBB_n}(KIDparam(nKID).tofitR{mpread,PBB_n})-1,...
                'o','MarkerFaceColor',cPr(nopread+1-mpread,:),'color',cPr(nopread+1-mpread,:));hold on;
            %fits
             plot(KIDparam(nKID).Pbb{mpread,PBB_n}*1E15,KIDparam(nKID).respfitR{mpread,PBB_n},'k','Linewidth',2);% fits
            plot(KIDparam(nKID).Pbb{mpread,PBB_n}*1E15,KIDparam(nKID).respfitth{mpread,PBB_n},'k','Linewidth',2);% fits
            %P point
            pt1=min([min(KIDparam(nKID).Phaseresp{mpread,PBB_n}) min(KIDparam(nKID).Radiusresp{mpread,PBB_n}-1)]);
            pt2=max([max(KIDparam(nKID).Phaseresp{mpread,PBB_n}) max(KIDparam(nKID).Radiusresp{mpread,PBB_n}-1)]);
            plot([KIDparam(nKID).Pbbnoise(1,PBB_n)*1e15 KIDparam(nKID).Pbbnoise(1,PBB_n)*1e15],[pt1 pt2],'-r')
            clear pt1 pt2
            
            
            subplot(3,3,8)
            Pofset=max(KIDparam(nKID).Pbb{mpread,PBB_n}*1E15)-min(KIDparam(nKID).Pbb{mpread,PBB_n}*1E15);
            plot(KIDparam(nKID).Pbb{mpread,PBB_n}*1E15+(mpread-1)*Pofset,KIDparam(nKID).Phaseresp{mpread,PBB_n},'s','color',cPr(nopread+1-mpread,:));hold on;
            plot(KIDparam(nKID).Pbb{mpread,PBB_n}*1E15+(mpread-1)*Pofset,KIDparam(nKID).Radiusresp{mpread,PBB_n}-1,'o','color',cPr(nopread+1-mpread,:))
            title('Response+fit+range, offs. for clarity')
            ylabel('Reponse')
            xlabel('Loading (fW)')
            %Overplot part of response that was used for the fit (closed symbols)%
            plot(KIDparam(nKID).Pbb{mpread,PBB_n}(KIDparam(nKID).tofitR{mpread,PBB_n})*1E15+(mpread-1)*Pofset,... 
                KIDparam(nKID).Phaseresp{mpread,PBB_n}(KIDparam(nKID).tofitR{mpread,PBB_n}),...
                's','MarkerFaceColor',cPr(nopread+1-mpread,:),'color',cPr(nopread+1-mpread,:));hold on;
            plot(KIDparam(nKID).Pbb{mpread,PBB_n}(KIDparam(nKID).tofitR{mpread,PBB_n})*1E15+(mpread-1)*Pofset,...
                KIDparam(nKID).Radiusresp{mpread,PBB_n}(KIDparam(nKID).tofitR{mpread,PBB_n})-1,...
                'o','MarkerFaceColor',cPr(nopread+1-mpread,:),'color',cPr(nopread+1-mpread,:));hold on;
            %Overplot with fit
            plot(KIDparam(nKID).Pbb{mpread,PBB_n}*1E15+(mpread-1)*Pofset,KIDparam(nKID).respfitR{mpread,PBB_n},'k','Linewidth',2);% fits
            plot(KIDparam(nKID).Pbb{mpread,PBB_n}*1E15+(mpread-1)*Pofset,KIDparam(nKID).respfitth{mpread,PBB_n},'k','Linewidth',2);% fits
           
            subplot(3,3,9)
            pbla(mpread)=plot(KIDparam(nKID).S21Real{mpread,PBB_n},KIDparam(nKID).S21Imag{mpread,PBB_n},'-k','LineWidth',1,'color',cPr(nopread+1-mpread,:));hold on;%volt S21
            xlabel('Re');ylabel('Im ');
            axis tight;
            plot(KIDparam(nKID).Reresp{mpread,PBB_n},KIDparam(nKID).Imresp{mpread,PBB_n},'.','color',cPr(nopread+1-mpread,:))
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
%Store Popt .csv
rowi=1;
for nKID=1:nokids
    for tbb=1:length(KIDparam(nKID).Popt) %one Popt per BB temperature
        OptP(rowi,2)=KIDparam(nKID).KIDid(1,1);%KID id
        OptP(rowi,1)=KIDparam(nKID).Tbbnoise(KIDparam(nKID).Poptindex(tbb),tbb);%emperature
        OptP(rowi,3)=KIDparam(nKID).Pread(KIDparam(nKID).Poptindex(tbb),tbb);%Popt check
        OptP(rowi,4)=KIDparam(nKID).Popt(tbb); %Activate if you want to
        %check if all is fine; this line will not be used
        rowi=rowi+1;
    end
end
PoptFile = [resppathy,'Popt.csv'];
PoptHeader = {'T (K)','KIDID','Popt (dBm)'}';
WriteSRONcsv(PoptFile,OptP(:,1:3),PoptHeader,'%.3g')
%

rmpath([pwd,filesep,'subroutines']);

fprintf('NEP_vs_loadingv12B is finished \n')

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
% KIDparam(nKID).phaseNEPamplifier[mpread,PBB_n]                         %SAME FOR PHASE
% KIDparam(nKID).radiusNEPPopt[PBB_n]                       %amplitude NEP at optimum readout power, 1D double
% KIDparam(nKID).Poptindex[PBB_n]                           %index of optimum readout power, 1D double
% KIDparam(nKID).phaseNEPPopt[PBB_n]                        %phase NEP at optimum readout power, 1D double
% KIDparam(nKID).Popt[PBB_n]                                %optimum readout power, 1D double
% KIDparam(nKID).preadcoefficientR[PBB_n]                   %(optional) coefficient of fit to the readou-power dependent amplitude responsivity, 1D double
% KIDparam(nKID).preadcoefficientth[PBB_n]                  %(optional) coefficient of fit to the readou-power dependent phase responsivity, 1D double
% KIDparam(nKID).stdpreadcoefficientR[PBB_n]                %(optional) std of preadcoefficientR, 1D double
% KIDparam(nKID).stdpreadcoefficientth[PBB_n]               %(optional) std of preadcoefficientth, 1D double
        
end %of the main function

function[powersum,NEP,BBcal,method]=blackbody_int(TTT,BBcal,method)
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
    BBcal.T=[2:0.05:100]';
    [BBcal.power,~,BBcal.NEP,method]=blackbody_12(BBcal.T,method,1); % [TotalPbb,Filtertransmission,NEP,method]
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
