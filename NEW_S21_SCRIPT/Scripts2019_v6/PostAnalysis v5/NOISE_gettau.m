function tau = NOISE_gettau
%
% This function can be run after NOISEanalysis has finished analysing the
% data. It load the NOISEanalysis workspace from file and fits the spectrum
% to get tau
% only works for data with  1 readout power to get Tdep.
%
%INPUT: All user defined in this function (see below)
usephasefortau  = 1; % 0 for radius 1 for phase. If radius is present it works better (use 0)
maxtau          = 5e-3; % maximum realistic lifetime
fref            = 80; %Ref frequecny wheret he noise is FOR SURE flat
blawindow       = 5; %amount of pts around fref used to get n oise level. default = 5
taustart        = 0.1e-3; % start value fit
Path            = [cd '/..']; %root path where data is, one higher than the scripts
FFT_subdir      = '/FFT/Power';

Noisepath = [Path,FFT_subdir,filesep];
format('long','g'); %Set display format of numbers to 7 digits
addpath([pwd,filesep,'subroutines']);
close all; %close all open plots to remove clutter.

%==========================================================================
% Load the workspace
%==========================================================================
NoiseVariables = {'NOISE','KIDnumbers','IndexPsort','ChipInfo'};
load([Noisepath,'Noise.mat'],NoiseVariables{:})
for kidn=1:length(KIDnumbers) % LOOP OVER ALL UNIQUE KID
    figure(kidn)
    PTindex = IndexPsort{kidn,2};
    IndexPref = IndexPsort{kidn}(ChipInfo.IndexPref,1);%index of optimum power in noise array
    Ntemperatures = length(NOISE(IndexPref).Temperature);
    Tcolors = colormapJetJB(Ntemperatures);
    for nT=1:Ntemperatures
        nP = 1; %look only at first readout power (no idea what all these indices mean anymore)
        semilogx(NOISE(IndexPref).FFTnoise{PTindex(nP,nT,2,1)}(:,1),...
            NOISE(IndexPref).FFTnoise{PTindex(nP,nT,2,1)}(:,3),...
            '--','color',Tcolors(nT,:),'LineWidth',2);hold on
        % fit
        if usephasefortau==1
            noisedata = NOISE(IndexPsort{kidn,1}(nP,1)).FFTnoise{PTindex(nP,nT,2,1)}(:,2);
        else
            noisedata = NOISE(IndexPsort{kidn,1}(nP,1)).FFTnoise{PTindex(nP,nT,2,1)}(:,3);
        end
        fdata = NOISE(IndexPsort{kidn,1}(nP,1)).FFTnoise{PTindex(nP,nT,2,1)}(:,1);
        indfref =find(fdata>fref,1);
        [tau(kidn,nP,nT), noiselevel(kidn,nP,nT), setupnoise(kidn,nP,nT)] = gettau_noise(fdata,noisedata,indfref,blawindow,maxtau,0,taustart);
        tau_out{kidn}(:,1) =  tau(kidn,nP,nT);
        tau_out{kidn}(:,2) =  NOISE(IndexPref).Temperature(nT);
        legendstr{nT} = ['\tau = ' num2str(1e3*tau(kidn,nP,nT),'%.3g') ' msec'];
    end
    legend(legendstr,'location','best')
    for nT=1:Ntemperatures
        for nP=1:length(IndexPsort{kidn,1}(:,1))
            semilogx(NOISE(IndexPsort{kidn,1}(nP,1)).FFTnoise{PTindex(nP,nT,2,1)}(:,1),...
                NOISE(IndexPsort{kidn,1}(nP,1)).FFTnoise{PTindex(nP,nT,2,1)}(:,2),...
                '-','color',Tcolors(nT,:),'LineWidth',2)
            %plot fit
            semilogx(fdata,10*log10(noiselevel(kidn,nP,nT)./(1 + (2*pi*tau(kidn,nP,nT)*fdata).^2)+setupnoise(kidn,nP,nT)),'k--');
        end
        xlabel('F [Hz]')
        ylabel('S_x [dBc/Hz]')
        xlim([10,0.2e6])
        ylim([-100,-40])
        grid on
        
        
    end
    title(['KID' num2str(KIDnumbers(kidn),'%.0f')]);
    MakeGoodFigure(15,9,11,[Noisepath 'KID' num2str(KIDnumbers(kidn),'%.0f') '_tau']);
end


end