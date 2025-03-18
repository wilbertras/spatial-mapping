function [Ich,Qch]=binfile
%R. Barends, Delft, May 2010

%opens a binary file made by labview and reads 64 bit floating point data.
%labview 'header': First 4 bytes: rows in uint, second 4 bytes: columns in uint. 
%Than data row by row in BIG ENDIAN
%PWELCH = PSD * 2/samplefreq!!!

%filename='D:\pieterdevisser\Desktop\dataKID19\KID19_75dBm_Tdep_TDmed_TmK217.bin'
filename='Q:\k30_Pieter__14_10_10 12_16\TD_2D\KID9_79dBm__TDmed_TmK120.bin';

%filename='D:\pieterdevisser\Desktop\dataKID19\KID19_75dBm_Tdep_TDfast_TmK217.bin'
%filename='D:\pieterdevisser\Desktop\D19C3\D19C3noise\4-74-350-ON-IQ.bin'


%filename='d:\kids\matlab\ananoise\Ktest-77-0-On-IQ.bin';
%filename='d:\kids\matlab\ananoise\3-87-334-ON-IQ.bin';

%filename='O:\MeasPC_BackUp100716\F\data\D14C3-3rdrun\noise\4-83-312-OFF-IQ.bin';
%filename='O:\MeasPC_BackUp100716\F\data\B4-C60-SiOx-2\LEDnoise\K39lednoisenicest\nK39-1k-70-337-off-IQ.bin'; % PRL   PRL   PRL
%filename='e:\kids-measurements\TaNb\150nmTa-B4-C60-SiOx2\LEDnoise\k39lednoisenicest\nK39-noled-70-337-on-IQ.bin';

%filename='e:\kids-measurements\TaNb\150nmTa-B4-C60-SiOx2\LEDnoise\k39lednoisenicest\nK39-1k-70-749-on-IQ.bin';
%filename='e:\kids-measurements\TaNb\150nmTa-B4-C60-SiOx2\LEDnoise\k39lednoisenicest\nK39-1k-10m-70-500-on-IQ.bin';

%filename='e:\kids-measurements\TaNb\150nmTa-B4-C60-SiOx2\LEDnoise\k39lednoisenicest\nK39-noled-70-750-on-IQ.bin';
%filename='e:\kids-measurements\TaNb\150nmTa-B4-C60-SiOx2\LEDnoise\k39lednoisenicest\nK39-noled-10m-70-500-on-IQ.bin';
%filename='d:\kids-measurements\TaNb\k44-lowfreq\K44-77-331-on-iq.bin';


%filename='e:\kids-measurements\TaNb\150nmTa-B4-C60-SiOx2\LEDnoise\k44lednoisemooi\K44-1k-77-350-on-IQ.bin';
%filename='e:\kids-measurements\TaNb\150nmTa-B4-C60-SiOx2\LEDnoise\k44lednoisemooi\K44-1k-77-350-off-IQ.bin';
%filename='e:\kids-measurements\TaNb\150nmTa-B4-C60-SiOx2\LEDnoise\k44lednoisemooi\K44-77-335-on-IQ.bin';

%filename='e:\kids-jochem\taled\K9\KID9_64dBm_LED0mA_TDfastT0.099996.bin';
%filename='e:\kids-jochem\taled\K9\KID9_64dBm_LED1mA_TDmedT0.099980.bin';
%filename='e:\kids-jochem\taled\K9\KID9_64dBm_LED2mA_TDmedT0.099996.bin';
%filename='e:\kids-jochem\taled\K9\KID9_64dBm_LED4mA_TDfastT0.099980.bin';
%filename='e:\kids-jochem\taled\K9\KID9_64dBm_LED0_125mA_TDmedT0.099996.bin';

%filename='E:\KIDs-Jochem\measurements\B6C58K13_cont_ill\TD_0D\kid13_92dbm_led0mA_TDmed_TmK100.bin';

%filename='E:\KIDs-Jochem\measurements\B9C37K39\TD LEDFFT\KID39_61dBm_LED32mA_T100m.bin';


%filename='e:\aplkid-fig2data\d8\d8c3noise\3-77-799-on-iq.bin'; %

[fid,message]=fopen(filename,'r','ieee-be.l64'); %open in Big endian 64 bit floating point number
message
filenameforsize=dir(filename);

%rows=fread(fid,1,'uint')
%columns=fread(fid,1,'uint')

%fseek(fid,8,-1);

%Matrix=fread(fid,[columns,rows],'float64')';

jochem=0;
if jochem==1
    rows=filenameforsize.bytes/8/3;
else
    rows=filenameforsize.bytes/8/2;
end

rows

columns=2
%med
rows=2e6;
freq=50.000000E+3;
filterfreq=1/2e-3;
%fast
%rows=2e5;
%freq=1e6;

%
%rows=50e3;
%freq=250e3;

%M1=fread(fid,[columns,rows],'float64')'
%M2=fread(fid,[columns,rows],'float64')'

%blabla




%reduces memory
readlength=20000; %read readlength columns at a time
pp=ceil(rows/readlength);
I(1:rows)=0;
Q=I;

idx=1;
for tel=1:pp
    rlen=rows-(tel-1)*readlength;
    if rlen>readlength rlen=readlength;end
    if jochem==1
        Matrix=fread(fid,[columns+1,rlen],'float64')';
        I2=Matrix(:,2)';
        Q2=Matrix(:,3)';
        freq=round(1/abs(Matrix(1,1)-Matrix(2,1)));
    else
        Matrix=fread(fid,[columns,rlen],'float64')';
        I2=Matrix(:,1)';
        Q2=Matrix(:,2)';
    end
    I( idx : idx+rlen-1 )=I2;
    Q( idx : idx+rlen-1 )=Q2;
    idx=idx+rlen;
end

fclose(fid);
clear Matrix;

freq



jochem=1;
if jochem==1
    %renormalize IQ
    r=sqrt(I.^2+Q.^2);
    p=atan2(Q,I);
    r=r/mean(r);%set amplitude to 1 (in stead of some value like 0.3
    I=r.*cos(p);
    Q=r.*sin(p);
end


%resample to conserve memory
%newfreq=freq/1; %new frequency of resampled data
%I=resample(I,newfreq,freq);
%Q=resample(Q,newfreq,freq);
%rows=length(I);
%freq=newfreq;




t=linspace(0,1,rows)*rows/freq;
tplot=linspace(0,t(end),10);
%ghguy

figh=figure;
subplot(2,2,1)
stdI=std(I)

stdplus=zeros(size(tplot));stdplus(:,:)=mean(I)+5*stdI;
stdmin=zeros(size(tplot));stdmin(:,:)=mean(I)-5*stdI;
plot(t,I,t,Q,tplot,stdplus,'-r',tplot,stdmin,'-r')%plot(t,I,t,Q)
title('I(t) & Q(t)')
legend('I','Q','I+5\sigma','I-5\sigma')

nrpointsfilter=freq/filterfreq
evencheck=nrpointsfilter/2-floor(nrpointsfilter/2)
if evencheck==0
    nrpointsfilter=nrpointsfilter+1
end
Ifilter=smooth(I,nrpointsfilter)';
Qfilter=smooth(Q,nrpointsfilter)';


r=sqrt(I.^2+Q.^2);
p=atan2(Q,I);
p=mod(p,2*pi);

n=length(I);
p=p-mean(p);
r=r-mean(r);


skip=0;
if skip==0
average=32;

%[ppp,f]=pwelch(p,hann(n),[],n,freq,'onesided');
[ppp,fm]=pwelch(I-mean(I),n/average,[],n/average,freq,'onesided');
%[ppp,fm]=pwelch(I,n/average,[],n/average,freq,'onesided');
[prr,fm]=pwelch(Q-mean(Q),n/average,[],n/average,freq,'onesided');
%[prr,fm]=pwelch(Q,n/average,[],n/average,freq,'onesided');
%[prr,fm]=psd(Q-mean(Q),n/average,freq,n/average);%prr=prr*2/freq;

[ppr,fm2]=csd(I-mean(I),Q-mean(Q),n/average,freq,n/average);ppr=ppr*2/freq;  % csd is differently defined than in the book, angle is negative. I use: R_{xy}(tau)=E{ x(t)y(t-tau) }
[pprnew,fm3]=cpsd(I-mean(I),Q-mean(Q),n/average,[],n/average,freq,'onesided');
%[ppr,fm2]=csd(I,Q,n/average,freq,n/average);ppr=ppr*2/freq;  % csd is differently defined than in the book, angle is negative. I use: R_{xy}(tau)=E{ x(t)y(t-tau) }

%csd(I-mean(I),Q-mean(Q),n/average,freq,n/average);
else
    fm=0;
    ppp=0;
    prr=0;
    fm2=0;
    ppr=0;    
end



%[ppp,f]=pwelch(p);%;,[],[],[],freq,'onesided');
%[prr,f]=pwelch(r);%,[],[],[],freq,'onesided');
%[pxx,f]=pwelch(I,[],[],[],freq,'onesided');
subplot(2,2,2)

PPP=10*log10(ppp);
PRR=10*log10(prr);
PPR=10*log10(abs(ppr)); %- 10log10 (freq/2)
PPRnew=10*log10(abs(pprnew));
%[ff,pp]=correlation(p,p,freq,5);
figure(figh);

semilogx(fm,PPP,fm,PRR,fm2,PPR)
title('psd I, Q and csd IQ using pwelch')
legend('radius','phase','cross')


%average=1



makeIQrp=1;
if makeIQrp==1
    r=sqrt(I.^2+Q.^2);
    p=atan2(I,Q);
    I=r;
    Q=p;
end


%Q(1:100)-mean(Q)
%b46hnuj


[f,cpsdII,ph]=correlation(I,I,freq,average,1);

[f,cpsdQQ,ph]=correlation(Q,Q,freq,average,1);

[f2,cpsdIQ,cpsdIQphase]=correlation(I,Q,freq,average,1);
cpsdIQphase=unwrap(cpsdIQphase,pi/1.5);


if f~=f2
    disp('f and f2 not equal')
end
fig3=2;
if fig3==1
    figure(figh)
    subplot(2,2,3)
    semilogx(f,cpsdII,f,cpsdQQ,f,cpsdIQ,fm,PPP,fm,PRR,fm2,PPR,fm3,PPRnew);
    legend('cpsdII','cpsdQQ','cpsdIQ','radius','phase','cross','crossnew')

    title('psd I psd Q and mag. of crosspsd IQ')
    xlabel('f(Hz)')
    ylabel('dBc/Hz')
elseif fig3==2
    figure(figh)
    subplot(2,2,3)
    stdIfilter=std(Ifilter)
    stdplusfilter=zeros(size(tplot));stdplusfilter(:,:)=mean(Ifilter)+5*stdIfilter;
    stdminfilter=zeros(size(tplot));stdminfilter(:,:)=mean(Ifilter)-5*stdIfilter;
    plot(t,Ifilter,t,Qfilter,tplot,stdplusfilter,'-r',tplot,stdminfilter,'-r')%plot(t,I,t,Q)
    title('I(t) & Q(t) filtered at qp-lifetime')
    legend('I filtered','Q filtered','I filtered+5\sigma','I filtered-5\sigma')
end

subplot(2,2,4)
semilogx(f,cpsdIQphase);
title('phase of crosspsd IQ')

ppd=30;
[a,b]=logsmooth(fm2,PPR,ppd);
[c,d]=logsmooth(f2,cpsdIQ,ppd);
[e,f]=logsmooth(fm3,PPRnew,ppd);
figure
semilogx(a,b,c,d,'--',e,f,':')
legend('csd','correlation','cpsd')

writegraph([f;cpsdII;cpsdQQ;cpsdIQ;cpsdIQphase],'f AA pp ApMag Apph','psd')

%secondspectrum(Q,freq,100)