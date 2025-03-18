%%
figure(2)
KID=KID11_80_LED;%from direct import of raw file, has F, phase noise, R noise
close all;
semilogx(KID(3:end,1),KID(3:end,3));axis tight;hold on;
%%
photonnoise=-98.3;
ampnoise=-101.3;
photonnoise=10*log10(10^(photonnoise/10)-10^(ampnoise/10));
tau=110e-6;
roll_off=10*log10(10^(ampnoise/10)+10^(photonnoise/10)./(1+(KID(:,1)*2*pi).^2*tau.^2));
semilogx(KID(:,1),roll_off,'g');
%%
corr_R_data=10*log10(10.^(KID(3:end,3)/10)-10^(ampnoise/10));
semilogx(KID(3:end,1),corr_R_data,'k');
roll_off2=10*log10(10^(photonnoise/10)./(1+(KID(:,1)*2*pi).^2*tau.^2));
semilogx(KID(:,1),roll_off2,'r');
%%
legend('data',['fit for tau= ' num2str(tau*1e3) ' msec'],'data-noisefloor',['fit for tau= ' num2str(tau*1e3)]);