% Tbb = A 1xN vector of temperatures of the blackbody.
% pad = path where the filterfiles are, can be left out, pad = pwd in that case. it will look in [pad,filesep,'filterfiles',filesep];%
% method. A struct containing information about the used setup including:
%         filters, lens size, and aperture determination. The following
%         struct elements are required:
%
%   method.pol = Number of polarizations used (equal to either 1 or 2)
%
%   method.filter='1_6THz'              % 1.6 THz SPICA SAFARI stack
%   method.filter='350 GHz 4Filters'    %350 Cardiff bandpass with additional LPF
%   method.filter='350 GHz Deshima'     %350 GHz LPF (2x) on cold box, and K1817BPF and K1785 LPG on chip and 1 Thz LPF on BB%
%   method.filter='850 GHz'             %850 GHz filterstack used first with LT010
%   method.filter='650 GHz'             %650 GHz filterstack with B768.dat
%   method.tp='lambda^2'                %using lambda^2 throughput over the entire filter band%
%   method.tp='Geometrical'             %Using geometrical calculation for throughput determination.
%   
%   ONLY FOR method.tp='Geometrical'%
%   method.lensdiameter=1 %lens diameter in millimeters. 
%   method.opening_angle=10 %angle wrt optical propagation (i.e. total angle is 2x larger) in degrees%
%
%   ONLY FOR method.tp='lambda^2'%
%   method.eta_c = 1.0              %[Optional, Default = 1.00] total SO couplig (CST), only for lambda^2%
%       
%   OPTIONAL    
%   method.GR = 4;                  % prefactor in the GR noise term of the NEP. default = 4
%                                   % following Flanigan et al (APL 108, 083502 (2016). %
%                                   % We used GR=2 in the past (in Pieter/Reineir/SY work%
%   method.freqresolution = 2       %[Optional, Default = 2] frequency resolution in GHz used
%                                       to integrate the blackbody spectrum with all filters.
%                                       NOTE: for the 350 GHz and 325 GHz filters 2 GHz is
%                                       adviced as minimum.
%   method.MaxFreq = 5e12           %[Optional, Default= 5e12] max Freq. of integration%
%   method.Delta = 45.6             %[Optional, Default = 45.6] value of Delta (half
%                                       the superconducting gap) in GHz
%   method.eta_pb = 0.57            %[Optional, Default = 0.57] value of the photon
%                                       pair breaking efficiency. Default from Kozorezov
%                                       et.al. PRB 2000
%
%  plotdata. (either 1 [true] or 0 [false]) If true, the routine plots the
%            photon noise NEP's of Al resonators. 
% 
%OUTPUT:
% TotalPbb.             % 1xN vector of blackbody power received by the lens/antenna ateach blackbody temperature%
% method.Etendue.       % Optical throughput of the system at the center frequency of the band %
% method.centrefreq     % central frequency
% method.filterBW       % effective BW
% Filtertransmission.   %2xM vector of total filter transmission [frequency (Hz),transmission].
%                       frequency is given with the resolution specified in method.freqresolution  
% NEP.                  %struct containing four 1xN vectors that for each blackbody temperature give the expected NEP due to:%
%   NEP.g_r         %Generation Recombination Noise (using method.Delta)
%   NEP.poisson     %Poissonian photon noise
%   NEP.wave        %Wave Bunching Photon Noise
%   NEP.totphoton   %Total NEP due to Photon, GR and wavebunching noise.
%
%SUBROUTINES:
% getfilterform (included below)
% plotresults (included below)
%
%REQUIRED FILES (filters, located in /filterfiles subdir):
% W969 37cmLPESCUBAII.txt
% W1275_350GHz.dat
% totalLbandfilters.txt
% B386 18cm LPE.dat
% W1052 14cm LPE SCUBAII.dat
% H20coupling efficiency upto1_3THz.txt
% H10coupling efficiency.txt

%close all
clear all
clc
plotfiguresBBscript=1;
method.filter='350 GHz 4Filters';    %Cardiff filter
method.tp='lambda^2';       %'lambda^2' 'Geometrical'
%only Lambda^2
method.eta_c = 1;        %optical efficiency CST
%only geometrical
method.lensdiameter=1;%lens diameter in millimeters. 
method.opening_angle=10 ;

method.pol=1;               % one polarization
method.GR = 4;              %R noise prefactor (4)

Tbb = [2.2:0.2:10 11:1:40 ]';
DesiredPowers = [1e-6;2e-6;5e-6;1e-5;2e-5;5e-5;1e-4;2e-4;5e-4;0.001;0.002;0.005;0.01;0.02;0.05;0.1;0.2;0.5;1;2;5;10;20;50;100;200;500;1000;2000;5000]*1e-15;  %to get the Tbb for these Powers
pad=(pwd);
addpath([pwd,filesep,'subroutines']);
%%
close all
disp(method);
[TotalPbb,FilterTransmission,NEP,method]=blackbody_14(Tbb,method,plotfiguresBBscript,pad);
%%
DesiredTemperatures = interp1(log10(TotalPbb),Tbb,log10(DesiredPowers),'linear');
bla = isnan(DesiredTemperatures);
DesiredTemperatures(bla)=[];
DesiredPowers(bla)=[];
TParray = [DesiredTemperatures 1e15*DesiredPowers];
clear bla DesiredTemperatures DesiredPowers
%% power figure
figure(1)
semilogy(Tbb,TotalPbb*1e15,'k-','LineWidth',2); hold on
semilogy(TParray(:,1),TParray(:,2),'ro','MarkerSize',6,'MarkerFaceColor','r')
xlabel('T [K]')
ylabel('P [fW]')
axis tight
hold off
grid on
title([method.filter ', \eta_c = ' num2str(method.eta_c)])

%% extra figure
figure(2)
axes('LineWidth',2,'FontSize',16)
subplot(2,2,1)

plot(FilterTransmission(1,:)/1e9,FilterTransmission(2,:),'k-');hold on
%
plot([method.centrefreq-method.filterBW/2 method.centrefreq-method.filterBW/2 ...
    method.centrefreq+method.filterBW/2 method.centrefreq+method.filterBW/2]/1e9,...
    [0 max(FilterTransmission(2,:)) max(FilterTransmission(2,:)) 0],'-b')
plot([method.centrefreq method.centrefreq]/1e9,...
    [0 1.2*max(FilterTransmission(2,:))],'-r')
xlim([min(FilterTransmission(1,:))/1e9,1e3]);
ylim([0 1.2*max(FilterTransmission(2,:))])
xlabel('Frequency [GHz]')
ylabel('Filter Transmission')

title(['F0= ' num2str(round(method.centrefreq/1e9)) ' BW= ' num2str(round(method.filterBW/1e9)) ...
    ' GHz, T= ' num2str(0.001*max(round(FilterTransmission(2,:)*1000)))])
hold off
%
subplot(2,2,2)

semilogy(Tbb,TotalPbb*1e15,'k-','LineWidth',2)
xlabel('T [K]')
ylabel('P [fW]')
hold off
grid on

subplot(2,2,3)
LH = '';%num2str(round(method.centrefreq/1e9));
%Overview of the calculated NEP's use to photon induced pair breaking
loglog(TotalPbb*1e15,NEP.poisson,['b'],'MarkerSize',8)
hold on
loglog(TotalPbb*1e15,NEP.wave,['g'],'MarkerSize',8)
loglog(TotalPbb*1e15,NEP.g_r,['r'],'MarkerSize',8)
loglog(TotalPbb*1e15,NEP.totphoton, ['k'],'MarkerSize',8)
legend([LH 'Poisson'],[LH 'Bunching'],[LH 'Recombination'],[LH 'NEP_{BLIP}'],'Location','Best')
ylabel('NEP [W/\surd{Hz}Hz]')
xlabel('P_{BB} (fW)')
title('NEP calc. overview @ detector')
axis tight;
grid on


subplot(2,2,4)
LH = num2str(round(method.centrefreq/1e9));
%Overview of the calculated NEP's use to photon induced pair breaking
hold on
NEP_nowave=NEP.totphoton./(NEP.g_r.^2+NEP.poisson.^2).^0.5;
NEP_noGR=NEP.totphoton./(NEP.wave.^2+NEP.poisson.^2).^0.5;
%NEP_nopoisson=(sqrt(NEP.wave.^2+NEP.g_r.^2)./NEP.totphoton)'
semilogx(TotalPbb*1e15,NEP_nowave',['g'],'MarkerSize',8)
semilogx(TotalPbb*1e15,NEP_noGR',['r'],'MarkerSize',8)
semilogx(TotalPbb*1e15,NEP.wave./NEP.totphoton,['--g'],'MarkerSize',8)
semilogx(TotalPbb*1e15,NEP.g_r./NEP.totphoton,['--r'],'MarkerSize',8)
%semilogx(TotalPbb*1e15,NEP_nopoisson, ['k'],'MarkerSize',8)

legend([LH ' NEP without wave/NEPtot'],[LH ' NEP without gr/NEPtot'],[LH ' NEPwave/NEPtot'],[LH ' NEPgr/NEPtot'],'Location','Best')
ylabel('Photon NEP [W/Hz^1/2]')
xlabel('P_{BB} (fW)')
title('NEP calc. overview @ detector')

axis tight;
xlim([1e-3 max(TotalPbb*1e15)])
grid on

