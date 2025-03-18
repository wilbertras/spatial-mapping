%%
%function [TDparam, TDoption, FITparam, allkidfiles] = noise_for_Jochem
function [allkidfiles] = noise_for_Jochem

close all
path = [cd '/../../TD_Power/'];        %Directory of bin-files
pathcal = 'empty'                   ;%Directory of off-res bin-files
saveindex = '1';                     %integer number, used for labeling your stored Matlab workspace (prevents overwriting of previous results)
                                    %filename will be: 'TDresults_saveindex.mat'
oldsaveindex = '1';                   %earlier results workspace to be loaded 'TDresults_oldsaveindex.mat' (to retrieve earlier defined fit ranges)


%This TDoption struct is convenient to transport the options through the functions

TDoption.nrsigma = 7;                 % peaks>TDoption.nrsigma*sigma are rejected (based on filtered TD-stream) %6 seems to be really minimum, 8 good value also for smaller peaks


TDoption.saveplot = 1;                %save plot with spectra
TDoption.usecalibration = 2;          %0-no calibration; 1-use offres-binfiles (requires off-res bin files!); 2-subtract amplifier noise floor from tail of PSD

%define acquisition frequency and size of the files for both bin-files
TDoption.freqmed = 50e3;              %standard is 50e3;
TDoption.rowsmed = 2e6;               %standard is 2e6; (ie sampling for 40 s at 50 kHz)
TDoption.freqfast = 1e6;              %standard is 1e6;
TDoption.rowsfast = 2e5;              %standard is 2e5; (ie sampling for 200 ms at 1 MHz)

TDoption.average = 32;                %Nr of TDoption.averages for PSD, should be a multiple of 2 and nrpoints/TDoption.average should be integer. 
                                      %Average is also the number of pieces over which peak rejection is done, ie if TDoption.average=32 and a peak is found, 
                                      %that 1/32th piece is thrown away
TDoption.ppd = 30;                    %points per decade in frequency for logsmooth after calculating PSD, 
                                      %30 (20) will give spectrum of 158 (105) points. 0 will give you the full spectrum, usually a million points
                                      %30 is used standard in the SRON labview    
                                      %NOTE: if you want to use fitranges from earlier processing, make sure ppd is the same!!!
TDoption.method = 0;                  %0=pwelch and cpsd (needs the Matlab signal processing toolbox!!) 1=with correlationfunction from Rami
TDoption.plotx = 1;                   %1: plot and save time domain plot for every temperature, useful to look for peaks in TD-stream 
                                      %(note, gives some 100Mb data per step); Saves in data directory; 2: plot PSD's for the different TDoption.methods, 'correlation' and 'cpsd', just for testing
TDoption.savex = 0;                   %save filtered time domain data in csv-file (note some 90 Mb per file); Saves in data direcotory
TDoption.rejectpeak = 1;              %1 only available option. If you want no rejection, make TDoption.nrsigma high
TDoption.twotimefit = 0;              %do also the fit to crossPSD with 2 timescales
TDoption.corroff = 1;                 %correct each piece of length/average with linear offset (to correct for very slow drifts, only done for 50 kHz files
TDoption.crosstime = 0;               %computes cross-correlation of amp and phase in TIME, takes quite some calculation time, because of >1M points

%check for binfiles in directory (only the 50 kHz sampled ones)%%%%%%%%%%%%%%%%%%%%%
noisefiles_med = dir([path '*TDmed*.bin']);
noisefiles_fast = dir([path '*TDfast*.bin']); 
number_files = length(noisefiles_med)+length(noisefiles_fast);
disp([num2str(number_files) ' files found in ' path ' '])
maxnumkids = 20;
maxnumpower = 20;
maxnumT = 50;       %used for initial array size

%get kidnr,power,temperature from filenames
tnoisekids = zeros(length(noisefiles_med),3); %this is a temporary variable
for nn=1:length(noisefiles_med)
    %error for textscan
    if isempty(noisefiles_med(nn).name) == 1
        error('file empty');
    end
    %This alternative for the textscan method works for both GR noise (dark) and photon noise (optical) bin-files
    %filename structure is always: 'KIDn_mdBm_Tchipy_TD*_TmKl.bin' where n=KID number, m=readout power, 
    %* = either 'med'(dark), 'fastt'(optical) for the 50 kHz file and 'fast' (dark) and 'slow' (optical) for the 1 MHz file
    %'Tchipy' is only there for optical files for dark files empty, l=temperature of bath(dark), blackbody (optical)
    [ligstreep] = regexpi(noisefiles_med(nn).name,'_');   %identify underscore positions in filename
    [ndBm] = regexpi(noisefiles_med(nn).name,'dBm');      %where 'dBm' starts, ie the position of the d
    [~,nTmK] = regexpi(noisefiles_med(nn).name,'TmK');    %end of 'TmK', ie position of the K
    [ndotbin] = regexpi(noisefiles_med(nn).name,'.bin');  %position of the dot
    tnoisekids(nn,1) = str2double(noisefiles_med(nn).name(4:ligstreep(1)-1));         %kid number
    tnoisekids(nn,2) = str2double(noisefiles_med(nn).name(ligstreep(1)+1:ndBm-1));    %readout power
    tnoisekids(nn,3) = str2double(noisefiles_med(nn).name(nTmK+1:ndotbin-1));         %temperature in mK, is either base temperature (dark) or blackbody temperature (optical)
end
%sort on kidnr 1,power 2,temperature 3
allkidfiles = sortrows(tnoisekids,[1 2 3]);
clear ligstreep ndBm nTmK ndotbin tnoisekids

%get all KID ID's in 1 array and put all powers in a cell array for each KID
n = 1;      % n=kids
m = 1;      % m=powers
k = 1;      % k=temperatures
maxm = 1;
maxk = 1;   % allows different #powers and #temp for each kid/power
noisekids = zeros(maxnumkids,1);
noisekids(1) = allkidfiles(1,1);
kidsP = cell(maxnumkids,1);
kidsPtemp = zeros(maxnumpower,1);
kidsPtemp(1) = allkidfiles(1,2);
kidsT = cell(maxnumkids,maxnumpower);
kidsTtemp = zeros(maxnumT,1);
kidsTtemp(1) = allkidfiles(1,3);
for nn=1:length(allkidfiles)
    if m > maxm; 
        maxm = m;
    end
    if noisekids(n) ~= allkidfiles(nn,1) %if new KID name found
        kidsPtemp(m+1:end) = [];
        kidsP{n} = kidsPtemp;
        kidsPtemp = zeros(maxnumpower,1);
        kidsTtemp(k:end) = [];
        kidsT{n,m} = kidsTtemp;
        kidsTtemp = zeros(maxnumT,1);
        n = n+1;
        kidsPtemp(1) = allkidfiles(nn,2);%this is to skip the 'new Power' if there is a new KID
        m = 1;
        k = 1;
        noisekids(n) = allkidfiles(nn,1);
    end
    if kidsPtemp(m)~=allkidfiles(nn,2) %if new Power found
        kidsTtemp(k:end)=[];kidsT{n,m}=kidsTtemp;kidsTtemp=zeros(maxnumT,1);
        m=m+1;
        if k>maxk; maxk=k-1; end
        k=1;
        kidsPtemp(m)=allkidfiles(nn,2);
    end
    kidsTtemp(k)=allkidfiles(nn,3);
    k=k+1;
end
%write into array for last step
kidsPtemp(m+1:end) = [];
kidsP{n} = kidsPtemp;
kidsTtemp(k:end) = [];
kidsT{n,m} = kidsTtemp;
clear kidsPtemp kidsTtemp nn
%delete some unneeded extra cells
noisekids(n+1:end) = [];
kidsP(n+1:end) = [];
kidsT(n+1:end,:) = [];
kidsT(:,m+1:end) = [];
%%%%%%here ends the filenamestuff%%%%%%%%%%%%%%%%%%%%%%%%%

kidnrstart = 1;
kidnrstop = length(noisekids);
powernrstart = 1;
powernrstop = length(kidsP{n});
tempnrstart=1;
tempnrstop=length(kidsT{n,m});
%loop KIDS
for n = kidnrstart:kidnrstop
    %loop powers
    powerarray=zeros(maxm,1);
    for m=powernrstart:powernrstop
        powerarray(m,1)=kidsP{n}(m);
        %loop temperatures
        for k=tempnrstart:tempnrstop
            %display current status
            disp(['Now processing: KID ' num2str(noisekids(n)) ' | -' num2str(kidsP{n}(m)) 'dBm | ' num2str(kidsT{n,m}(k)) ' mK' ]);
            
            TDparam(n).kidnr=noisekids(n);          %Put kid number in output struct
            TDparam(n).Pread(m,k)=kidsP{n}(m);      %Put readout power in output struct
            TDparam(n).Temp(m,k)=kidsT{n,m}(k);     %Put temperature in output struct
            
            %filename stuff of bin-files
            %med
            filetjes_med = dir([path 'KID' num2str(noisekids(n)) '_' num2str(kidsP{n}(m)) 'dBm__TDmed_TmK' num2str(kidsT{n,m}(k)) '.bin']);%
            if length(filetjes_med) == 1
                TDparam(n).filemed{m,k} = [path filetjes_med.name];
                if TDoption.usecalibration == 1
                    TDparam(n).filemedcal{m,k} = [pathcal filetjes_med.name];
                else
                    TDparam(n).filemedcal{m,k} = 'none';
                end
            else
                error('identical med files found')
            end
            
            [fm_med,SRRcal,SPPcal,SPRcal,SRRraw_med,SPPraw_med,SPRraw_med,nrrejected,sigmax5] = binfilefunction2(filetjes_med,TDoption.rowsmed,TDoption.freqmed,TDoption);
            
            SPRraw_med = 10*log10(real(SPRraw_med));
            figure(100)
            semilogx(fm_med,SPPraw_med,'b')
            hold
            semilogx(fm_med,SPRraw_med,'k')
            semilogx(fm_med,SRRraw_med,'r')
            
            %fast
            filetjes_fast = dir([path 'KID' num2str(noisekids(n)) '_' num2str(kidsP{n}(m)) 'dBm__TDfast_TmK' num2str(kidsT{n,m}(k)) '.bin']);%
            if length(filetjes_fast) == 1
                TDparam(n).filefast{m,k} = [path filetjes_fast.name];
                if TDoption.usecalibration == 1
                    TDparam(n).filefastcal{m,k} = [pathcal filetjes_fast.name];
                else
                    TDparam(n).filefastcal{m,k} = 'none';
                end
            else
                error('identical fast files found')
            end
            
            [fm_fast,SRRcal,SPPcal,SPRcal,SRRraw_fast,SPPraw_fast,SPRraw_fast,nrrejected,sigmax5] = binfilefunction2(filetjes_fast,TDoption.rowsfast,TDoption.freqfast,TDoption);
            
            SPRraw_fast = 10*log10(real(SPRraw_fast));
            figure(100)
            semilogx(fm_fast,SPPraw_fast,'--b')
            semilogx(fm_fast,SPRraw_fast,'--k')
            semilogx(fm_fast,SRRraw_fast,'--r')

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end %loop temperature
    end %loop powers
end %loop KIDS