%%
function [fm,SRRraw,SPPraw,SPRraw,nrrejected,sigmaRx5,sigmaPx5] = binfilefunction4(filename,meastime,freq,aver,TDoption)
% simplified from 2 by Jochem Aug 2019
% this is done in a separate function since it has to be done twice (50kHz and 1MHz sampled bin-files) for every noise spectrum
% 
% vs2 is with subtraction of calibration file (in FFT)
% 
% opens a binary file made by labview and reads 64 bit floating point data.
% labview 'header': First 4 bytes: rows in uint, second 4 bytes: columns in uint. 
% Than data row by row in BIG ENDIAN
% PWELCH = PSD * 2/samplefreq!!!
% 
% Variable names for the spectra:
% SRR: PSD of Radius
% SPP: PSD of Phase
% SPR: cross-PSD of R and P
% 
% %%%%%%%%%%%%%%%%%%%%%%%INPUT of binfilefucntion2%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%first four variables are different for med/fast bin-file, TDoption has the common options
% filename: name of the on-resonance bin-file to be processed
% meastime: time data is taken (if cal is used they have to be equal size).
% is the last item in xxdBm__all_td_averaged.dat file
% freq: sampling frequency of data-acquisition
%
% %%%%%%%%%%%%%%%%%%%%%%%%OUTPUT of binfilefunction2%%%%%%%%%%%%%%%%%%%%%%%%%
% fm: frequency of the following PSDs (x-axis) 
% SRRraw: radius PSD of the on-resonance file 
% SPPraw: phase PSD of the on-resonance file
% SPRraw: cross PSD of the on-resonance file; NOTE: phase and radius PSD come out in log-space, cross-PSD comes out in linear space and it can be a complex number
% nrrejected: number of rejected time-domain pieces due to peaks
% sigmax5: the actual radius-level corresponding to the rejection threshold
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% This TDoption struct is convenient to transport the options through the functions
% define acquisition frequency and size of the files for both bin-files
%%%%%%%%%%%FOR DATA-FILE ON RESONANCE
TDoption.average = aver;

[fid,message] = fopen(filename,'r','ieee-be.l64'); %open in Big endian 64 bit floating point number
if ~isempty(message)
    fprintf([message '\n']);
end

columns = 2;

%reduces memory
rows = meastime * freq;
readlength = 10000; %read readlength x columns 
pp = ceil(rows/readlength);
I(1:rows) = 0; 
Q(1:rows) = 0;

idx = 1;
for tel = 1:pp
    rlen = rows-(tel-1)*readlength;
    if rlen > readlength 
        rlen = readlength;
    end
    Matrix = fread(fid,[columns,rlen],'float64')';
    I2 = Matrix(:,1)';
    Q2 = Matrix(:,2)';
    I( idx : idx+rlen-1 ) = I2;
    Q( idx : idx+rlen-1 ) = Q2;
    idx=idx+rlen;
end

fclose(fid);
clear Matrix;

%renormalize the radius to 1 and calculate the phase with respect to the correct axis
makeIQrp = 1; %leave this option 1!!!%here variables I and Q are translated in R and P
if makeIQrp == 1
    r = sqrt(I.^2+Q.^2);
    R = r/mean(r);      %normalize radius to 1
    p = atan2(Q,I);     %This phase is with respect to the axis in the direction (I,Q)=(1,0) in stead of (-1,0) therefore next line
    P = pi-mod(p,2*pi);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% END OF READING THE 2 FILES %%%% left with I,Q,Ical,Qcal %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t = linspace(0,1,rows)*rows/freq;
tplot = linspace(0,t(end),10);%to plot straight lines based on a few points
R_reject = zeros(1,rows)*rows/freq;
P_reject = zeros(1,rows)*rows/freq;
Rmmsmooth_reject = zeros(1,rows);Pmmsmooth_reject = Rmmsmooth_reject; 
n = length(R);
telgood = 0;
while telgood < 2  %this is to catch that all data is filtered, it iteratively increases TDoption.nrsigma until some data gets through  
    %Filter (ie moving average) with the quasiparticle lifetime, to
    %detect in a better way the cosmic ray hits (and any other
    %quasiparticle creating hits)
    %NOW (6-12-10) use smoothtime,qplifetime /2
    nrpointsfilter = round(freq*TDoption.smoothtime/2); %ie freq/filterfreq
    evencheck = nrpointsfilter/2-floor(nrpointsfilter/2);
    if evencheck == 0
        nrpointsfilter = nrpointsfilter+1;
    end
    
    %subtract the mean from both variables, so we are calculating covariances below
    %for the PSD it does not matter, the possible DC signal we subtract here would only imply a peak at f=0
    Rmm = R-mean(R);
    Rmmsmooth = smooth(Rmm,nrpointsfilter)';
    Rmmsmoothorig = Rmmsmooth;
    stdRmm = std(Rmm);
    stdRmmsmooth1 = std(Rmmsmooth);
    
    Pmm = P-mean(P);
    Pmmsmooth = smooth(Pmm,nrpointsfilter)';
    Pmmsmoothorig = Pmmsmooth;
    stdPmm = std(Pmm);
    stdPmmsmooth1 = std(Pmmsmooth);

    sigmaarrayR = zeros(1,TDoption.average);
    sigmaarrayP = zeros(1,TDoption.average);
    for k = 1:TDoption.average
        sigmaarrayR(k) = std(Rmmsmooth(1,(k-1)*n/TDoption.average+1:k*n/TDoption.average));
        sigmaarrayP(k) = std(Pmmsmooth(1,(k-1)*n/TDoption.average+1:k*n/TDoption.average));
    end
    stdRmmsmooth = min(sigmaarrayR)*1.2;
    stdPmmsmooth = min(sigmaarrayP)*1.2;
    sigmaRx5 = stdRmmsmooth*TDoption.nrsigma;%for output
    sigmaPx5 = stdPmmsmooth*TDoption.nrsigma;%for output
    
    srrcal = zeros(n/TDoption.average/2+1,1);
    sppcal = zeros(n/TDoption.average/2+1,1);
    sprcal = zeros(n/TDoption.average/2+1,1);
    srrraw = zeros(n/TDoption.average/2+1,1);
    sprraw = zeros(n/TDoption.average/2+1,1);
    sppraw = zeros(n/TDoption.average/2+1,1);
    
    telgood = 0;
    for k = 1:TDoption.average
        %reset std of parts to std of total time series
        %refine std, if std of a part is LOWER than std of total 
        %(ie, make the criterion stricter)
        % Radius
        if stdRmm > std(Rmm(1,(k-1)*n/TDoption.average+1:k*n/TDoption.average))
            stdRmm = std(Rmm(1,(k-1)*n/TDoption.average+1:k*n/TDoption.average));
        end
        % Phase
        if stdPmm > std(Pmm(1,(k-1)*n/TDoption.average+1:k*n/TDoption.average))
            stdPmm = std(Pmm(1,(k-1)*n/TDoption.average+1:k*n/TDoption.average));
        end        
    end
    %continue with strictest possible stdRmm and stdPmm
    for k = 1:TDoption.average
        %first check for peaks higher than 10x the threshold in non-smoothed data
        % Radius
        stdcheckR = abs(Rmm(1,(k-1)*n/TDoption.average+1:k*n/TDoption.average)) > (TDoption.nrsigma*stdRmm*10);    %factor of 10 because criterion is refined later
        %stdcheckcalR = 0; %if no cal-file
         % Phase
        stdcheckP = abs(Pmm(1,(k-1)*n/TDoption.average+1:k*n/TDoption.average)) > (TDoption.nrsigma*stdPmm*10);    %factor of 10 because criterion is refined later
        %stdcheckcalP = 0; %if no cal-file
        
        %correct for low frequency ripples due to cooler, linear correction (only in on-resonance)
        %The rows criterion includes that the linear correction term will not be applied in the 'fast' (1MHz sampled) data
        if (TDoption.corroff == 1) && (TDoption.average > 6) && (rows > 1e6) && (TDoption.average < 129) && (sum(stdcheckR) == 0) && (sum(stdcheckP) == 0)
            avPmm1 = mean(Pmm(1,(k-1)*n/TDoption.average+1:(k-1)*n/TDoption.average+1+100));
            avPmm2 = mean(Pmm(1,k*n/TDoption.average-100:k*n/TDoption.average));
            avRmm1 = mean(Rmm(1,(k-1)*n/TDoption.average+1:(k-1)*n/TDoption.average+1+100));
            avRmm2 = mean(Rmm(1,k*n/TDoption.average-100:k*n/TDoption.average));
            corrPmm = linspace(avPmm1,avPmm2,rows/TDoption.average);
            corrRmm = linspace(avRmm1,avRmm2,rows/TDoption.average);
            %do the correction
            Rmm(1,(k-1)*n/TDoption.average+1:k*n/TDoption.average) = Rmm(1,(k-1)*n/TDoption.average+1:k*n/TDoption.average)-corrRmm;
            Pmm(1,(k-1)*n/TDoption.average+1:k*n/TDoption.average) = Pmm(1,(k-1)*n/TDoption.average+1:k*n/TDoption.average)-corrPmm;
            Rmmsmooth(1,(k-1)*n/TDoption.average+1:k*n/TDoption.average) = Rmmsmooth(1,(k-1)*n/TDoption.average+1:k*n/TDoption.average)-corrRmm;
            Pmmsmooth(1,(k-1)*n/TDoption.average+1:k*n/TDoption.average) = Pmmsmooth(1,(k-1)*n/TDoption.average+1:k*n/TDoption.average)-corrPmm;
        end
        
        %determine on the basis of the smoothed data whether the piece will be rejected
        % Radius
        stdchecksmoothR = abs(Rmmsmooth(1,(k-1)*n/TDoption.average+1:k*n/TDoption.average)) > (TDoption.nrsigma*stdRmmsmooth);
        %stdcheckcalsmoothR = 0;
        stdcheckR = abs(Rmm(1,(k-1)*n/TDoption.average+1:k*n/TDoption.average)) > (TDoption.nrsigma*stdRmm);   %redefine to filter out single point peaks
        % Phase = ignored here
        stdchecksmoothP = abs(Pmmsmooth(1,(k-1)*n/TDoption.average+1:k*n/TDoption.average)) > (TDoption.nrsigma*stdPmmsmooth);
        %stdcheckcalsmoothP = 0;
        stdcheckP = abs(Pmm(1,(k-1)*n/TDoption.average+1:k*n/TDoption.average)) > (TDoption.nrsigma*stdPmm);   %redefine to filter out single point peaks
        
        %this 'if' is the peak REJECTION
        %ie if either the data or the calibration file has a peak data is rejected
        if (sum(stdcheckR) > 0) || (sum(stdchecksmoothR) > 0) || (sum(stdcheckP) > 0) || (sum(stdchecksmoothP) > 0)  
            %REJECT
            %set values in rejected range to 0 just for plotting, to show
            %which parts are rejected.
            index = find(stdcheckR);
            R_reject(1,(k-1)*n/TDoption.average+1:k*n/TDoption.average) = Rmm(1,(k-1)*n/TDoption.average+1:k*n/TDoption.average);
            Rmmsmooth_reject(1,(k-1)*n/TDoption.average+1:k*n/TDoption.average) = Rmmsmooth(1,(k-1)*n/TDoption.average+1:k*n/TDoption.average);
            Rmm(1,(k-1)*n/TDoption.average+1:k*n/TDoption.average) = 0;
            P_reject(1,(k-1)*n/TDoption.average+1:k*n/TDoption.average) = Pmm(1,(k-1)*n/TDoption.average+1:k*n/TDoption.average);
            Pmmsmooth_reject(1,(k-1)*n/TDoption.average+1:k*n/TDoption.average) = Pmmsmooth(1,(k-1)*n/TDoption.average+1:k*n/TDoption.average);
            Pmm(1,(k-1)*n/TDoption.average+1:k*n/TDoption.average) = 0;
        else
            %fft is done for pieces that are not rejected (without peaks)
            %for method==0 you need the Matlab signal processing toolbox
            [srrpartraw,fmb] = pwelch(Rmm(1,(k-1)*n/TDoption.average+1:k*n/TDoption.average),n/TDoption.average,[],n/TDoption.average,freq,'onesided');
            [spppartraw,fmb] = pwelch(Pmm(1,(k-1)*n/TDoption.average+1:k*n/TDoption.average),n/TDoption.average,[],n/TDoption.average,freq,'onesided');
            [sprpartraw,fm2b] = cpsd(Rmm(1,(k-1)*n/TDoption.average+1:k*n/TDoption.average),Pmm(1,(k-1)*n/TDoption.average+1:k*n/TDoption.average),n/TDoption.average,[],n/TDoption.average,freq,'onesided');
            %do sum over the non-rejected pieces, store averaged raw
            srrraw = srrraw + srrpartraw;
            sprraw = sprraw + spppartraw;
            sppraw = sppraw + sprpartraw;
            
            telgood = telgood+1; %count how many good datapieces there were
        end

    end
    min5stdRmm = stdRmmsmooth*TDoption.nrsigma;%value used for plots
    min5stdPmm = stdPmmsmooth*TDoption.nrsigma;%value used for plots
    
    %devide sum by nr of averages that had no peaks (telgood)
    srrcal = srrcal/telgood;
    sppcal = sppcal/telgood;
    sprcal = sprcal/telgood;
    srrraw = srrraw/telgood;
    sprraw = sprraw/telgood;
    sppraw = sppraw/telgood;
    nrrejected = TDoption.average-telgood;
    if telgood < 2; 
        TDoption.nrsigma = 2*TDoption.nrsigma;
    end %if there were no pieces without peak, make nrsigma 2x larger and do it again (happes sometimes for too much readout power).
end

%plot time domain of smoothed data

minstdplus = zeros(size(tplot));
minstdplus(:,:) = min5stdRmm;
minstdmin = zeros(size(tplot));
minstdmin(:,:) = -min5stdRmm;

%time domain of smoothed data
subplot(2,2,1)
plot(t,Rmmsmoothorig','b');% data
hold on;
plot(t,Rmmsmooth_reject,'r-')
plot(tplot,minstdplus,'-g',tplot,minstdmin,'-g');%
title('R_{smooth} - mean, \tau_{qp} filtered')
legend('smoothed data','Rejected points','Location','South')
xlabel('time  (sec)');ylabel('R response');

minstdplus = zeros(size(tplot));
minstdplus(:,:) = min5stdPmm;
minstdmin = zeros(size(tplot));
minstdmin(:,:) = -min5stdPmm;

subplot(2,2,2)
plot(t,Pmmsmoothorig,'b');hold on
plot(t,Pmmsmooth_reject,'r-');
plot(tplot,minstdplus,'-g',tplot,minstdmin,'-g');
title(['\theta_{smooth} - mean, ' num2str(TDoption.nrsigma) '\sigma filtered with \tau_{qp}'])
legend('smoothed data','Rejected points','Location','North')
xlabel('time  (sec)');ylabel('\theta response');

%saveas(figh,[Figfiletd])
% close(figh) %these fig-files are so large that they need to be closed immediately

%logsmooth smoothes the data with an equal amount of points per decade in frequency
%do the logsmooth before converting to dBc!!!
if TDoption.ppd>0
        [fm2,sppraw] = logsmooth(fm2b,sppraw,TDoption.ppd);
        [fm,srrraw] = logsmooth(fmb,srrraw,TDoption.ppd);
        [fm,sprraw] = logsmooth(fmb,sprraw,TDoption.ppd);
        
        [fm2,sprcal] = logsmooth(fm2b,sprcal,TDoption.ppd);
        [fm,srrcal] = logsmooth(fmb,srrcal,TDoption.ppd);
        [fm,sppcal] = logsmooth(fmb,sppcal,TDoption.ppd);
else
    fm = fmb;
    fm2 = fm2b;
end

%these are the output spectra
%NOTE: CROSS COMES OUT OF BINFILEFUNCTION2 NON-LOGARITMIC
SRRraw = 10*log10(srrraw);
SPPraw = 10*log10(sprraw);
SPRraw = sppraw;%ie output is now (possibly, and often) complex!

if TDoption.savex == 1 %Save filtered time domain
    Filetdfilter = [filename(1:end-4) '-tdfiltered.csv'];
    keeswrite('Rfilter,Pfilter',Filetdfilter);
    dlmwrite(Filetdfilter, [Rmmsmoothorig' Pmmsmoothorig'],'-append','newline', 'pc', 'precision', '%.12e','delimiter', ',');
end

end%function binfilefunction2