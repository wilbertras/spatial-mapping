Ktd = 3;    %Kid to plot

nPbb = length(KIDparam(Ktd).Q);

kl = colormap(jet(nPbb));
for n = 1:nPbb
    %semilogx(KIDparam(Ktd).f_noise{n},KIDparam(Ktd).phasenoise{n},'color',kl(n,:));hold on
    KIDparam(Ktd).freqnoise{n} = 10*log10(10.^(KIDparam(Ktd).phasenoise{n}/10)  * (1/(4*KIDparam(Ktd).Q(n))^2));
    semilogx(KIDparam(Ktd).f_noise{n},KIDparam(Ktd).freqnoise{n},'color',kl(n,:));hold on
    xlabel('Frequency  (Hz)');
    ylabel('Normalised Frequency noise  (Hz^2/Hz)')
    xlim([0.5 0.5e6]);ylim([-220 -160])
    grid on;
    legendst{n}=['P = ' num2str(KIDparam(nKID).Pbbnoise(n)*1e15,'%.3g') ' fW'];
end
legend(legendst,'Location','Best');