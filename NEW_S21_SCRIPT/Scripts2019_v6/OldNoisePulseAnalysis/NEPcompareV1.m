function [KID,NOISE,NEP]=NEPcompareV1
% calcul;ates NEP and so on. Note that drdtheta.csv should only contain
% data at 1 temperature and noise and Tau should be taken at same T. No T
% dependence taken into account here. 2D version should be written if
% needed.
%NB: for large Q variation: The response arrays are resized such that the minimum
%array size in T is the only response remaining, all other T's are
%discarded.

%%% modifications:
% 15/9/10 changed so can run on one KID
% 16-3-11 JB: addad few paramaters to output file

%REQUIRED INPUT FILES:
%THESE MUST BE COPIED TO THE ROOT DIIRECTORY
% Popt.csv
% dRdtheta.csv

% Noise.mat
% NoiseCompare.mat
% ResponseS21.mat

%==========================================================================
% Main Input Parameters
%==========================================================================
Mainpath= [cd '/../Session2']; %Main path containing all data subdirectories of the measurement 
         
NoiseDir = '/FFT/Power'; %Subdirectory path to the noise measurement (and noise analysis output)
ResponseDir = '/S21/2D'; %Subdirectory path to the response measuremetn (and response analysis output)

etapb = 0.57; %Assumed pairbreaking efficiency

%==========================================================================
% Optional input parameters
%==========================================================================

%==========================================================================
% Initialize
%==========================================================================
close all; %to remove clutter
addpath([pwd,filesep,'subroutines']); %enable subroutines
format('long','g'); %set screen printing format

NEPheader = {'F [Hz]','NEPphase','NEPamp'}';
%==========================================================================
% Load the required variables from the noise and response analysis
% environment
%==========================================================================
Noisefile = [Mainpath,NoiseDir,filesep,'Noise.mat'];
Noisevariables = {'NOISE','KIDnumbers','IndexPsort','NOISEParameters'};
load(Noisefile,Noisevariables{:})

%NoiseCfile = [Mainpath,NoiseDir,filesep,'NoiseCompare.mat'];
%NoiseCvariables = {'NOISE'};
%load(NoiseCfile,NoiseCvariables{:})

Responsefile = [Mainpath,ResponseDir,filesep,'ResponseS21.mat'];
Responsevariables = {'KID','ChipInfo','KIDParameters'};
load(Responsefile,Responsevariables{:})

%==========================================================================
% Read Popt file (or load from noise.mat)
%==========================================================================
PoptFile = [Mainpath,filesep,'Popt.csv'];
fid = fopen(PoptFile);
if fid == -1
    disp('Warning NEPcompare: Cannot find Popt.csv. Using PoptData from NOISEanalysis environment instead.')
    load(Noisefile,'PoptData')
else
    fclose(fid);
    %There is a Popt file. Read it and use it.
    fprintf('Found Popt.csv, this will be used to indicate Popt.\n')
    [~,PoptData] = ReadSRONcsvV2(PoptFile,'',0);
    %PoptData will be Nx3. Where dim 2 contains [T,KIDid,|Popt [dBm]|]
end

%Convert the PoptData back into the LookUpTables (PoptData sorted into a
%cell with a cell element per KID). This takes into account any potential
%user modifications made in the csv file.
LookUpTables = cell(length(KIDnumbers),1);
for p=1:length(KIDnumbers)
    Match = PoptData(:,2) == KIDnumbers(p);
    LookUpTables{p,1} = zeros(sum(Match),4); %[T value,T index,P value,P index]
    LookUpTables{p,1}(:,1) = PoptData(Match,1);
    LookUpTables{p,1}(:,3) = PoptData(Match,3);
    %Find the indexes corresponding to these powers and temperatures
    for n=1:sum(Match)
        Pindex = find(NOISEParameters(:,1) == KIDnumbers(p) & ...
            NOISEParameters(:,9) == LookUpTables{p,1}(n,3));
        LookUpTables{p,1}(n,4) = Pindex;
        
        [~,Tindex] = min(abs(NOISE(Pindex).Temperature(:) - LookUpTables{p,1}(n,1)));
        LookUpTables{p,1}(n,2) = Tindex;
    end
end
clear Match Pindex Tindex

%==========================================================================
% Read dRdtheta file
%==========================================================================
dRdthetafile = [Mainpath,filesep,'dRdtheta.csv'];
fid = fopen(dRdthetafile);
if fid == -1
    error('NEPcompare: No dRdtheta.csv file found.')
else
    fclose(fid)
    fprintf('Found dRdtheta.csv.\n')
    [~,LifetimeData] = ReadSRONcsvV2(dRdthetafile,'',0);
end

%==========================================================================
% Loop over all the KIDs to calculate the NEP
%==========================================================================
%Initialize some variables
PoptT0index = zeros(length(KIDnumbers),2,3); 
%dim 2: [Pindex,T0index]
%dim 3: [NOISE,KID(Response),Lifetime]
Popt = zeros(length(KIDnumbers),1);
T0 = zeros(length(KIDnumbers),1);

NEPcounter = 0; %Counts the number of times there is full noise,response,lifetime matching

for n=1:length(KIDnumbers) %Loop over all KIDs with noise data
    %These KIDs will also have noise data for Popt as given by the LoopUpTable
    [~,T0LUT] = min(LookUpTables{n,1}(:,1));
    PoptT0index(n,1,1) = LookUpTables{n,1}(T0LUT,4);
    PoptT0index(n,2,1) = LookUpTables{n,1}(T0LUT,2);
    %Check for matching KID data (with as close as possible P and T0) in
    %the response
    KIDmatch = find(KIDParameters(:,1) == KIDnumbers(n));
    LTmatch = find(LifetimeData(:,1) == KIDnumbers(n));
    if KIDmatch(1) >= 1 && LTmatch(1) >=1
        %For the response data: seek minimum difference in read power
        [~,Pmatch] = min(abs(KIDParameters(KIDmatch,9) - NOISE(PoptT0index(n,1,1)).ReadPower(PoptT0index(n,2,1))));
        PoptT0index(n,1,2) = KIDmatch(Pmatch);
        %At minimum difference in read power seek minimum T difference
        [~,Tmatch] = min(abs(KID(PoptT0index(n,1,2)).Temperature(:) - NOISE(PoptT0index(n,1,1)).Temperature(PoptT0index(n,2,1))));
        PoptT0index(n,2,2) = Tmatch;
        
        if length(LTmatch) == 1
            %Only one matching KID index. Use this.
            PoptT0index(n,:,3) = LTmatch;
        else
            %For lifetime data: seek minimum difference in read power
            dPs = abs(-1*LifetimeData(LTmatch,5) - NOISE(PoptT0index(n,1,1)).ReadPower(PoptT0index(n,2,1)));
            Pmatches = find(dPs == min(dPs));
            if length(Pmatches) == 1
                %use this match
                PoptT0index(n,:,3) = LTmatch(Pmatch);
            else
                %find the minimum temperature difference
                [~,Tmatch] = min(abs(LifetimeData(LTmatch(Pmatch),4) - NOISE(PoptT0index(n,1,1)).Temperature(PoptT0index(n,2,1))));
                PoptT0index(n,:,3) = LTmatch(Pmatch(Tmatch));
            end
        end
        
        %Since there is noise, response and lifetime data for the KIDs we
        %extract all required data from NEP calculation
        NEPcounter = NEPcounter +1;
        NEP(NEPcounter).KIDnumber = KIDnumbers(n);
        NEP(NEPcounter).Temperature = ...
            [NOISE(PoptT0index(n,1,1)).Temperature(PoptT0index(n,2,1)),...
            KID(PoptT0index(n,1,2)).Temperature(PoptT0index(n,2,2)),...
            LifetimeData(PoptT0index(n,1,3),4)];
        NEP(NEPcounter).ReadPower = ...
            [NOISE(PoptT0index(n,1,1)).ReadPower(PoptT0index(n,2,1)),...
            KID(PoptT0index(n,1,2)).ReadPower(PoptT0index(n,2,2)),...
            -1*LifetimeData(PoptT0index(n,1,3),5)];
        %Obtain data from noise struct
        NEP(NEPcounter).FFTnoise = NOISE(PoptT0index(n,1,1)).FFTnoise{PoptT0index(n,2,1),1}(:,1:3);
        NEP(NEPcounter).MeanNoise = NOISE(PoptT0index(n,1,1)).MeanNoise(PoptT0index(n,2,1),:);
        NEP(NEPcounter).MeanFreqNoise = NOISE(PoptT0index(n,1,1)).MeanFreqNoise(PoptT0index(n,2,1),:);
        %Obtain data from response struct
        NEP(NEPcounter).Fres = KID(PoptT0index(n,1,2)).Fres(PoptT0index(n,2,1),1);
        NEP(NEPcounter).Ql(1) = KID(PoptT0index(n,1,2)).Ql(PoptT0index(n,2,1),1); %at T0
        NEP(NEPcounter).Ql(2) = max(KID(PoptT0index(n,1,2)).Ql(:,1)); %max(at any T)
        NEP(NEPcounter).Qi(1) = KID(PoptT0index(n,1,2)).Qi(PoptT0index(n,2,1),1); %at T0
        NEP(NEPcounter).Qi(2) = max(KID(PoptT0index(n,1,2)).Qi(:,1)); %max(at any T)
        NEP(NEPcounter).Qc(1) = KID(PoptT0index(n,1,2)).Qc(PoptT0index(n,2,1),1); %at T0
        NEP(NEPcounter).Qc(2) = max(KID(PoptT0index(n,1,2)).Qc(:,1)); %max(at any T)
        NEP(NEPcounter).Delta = KID(PoptT0index(n,1,2)).Delta;
        NEP(NEPcounter).Responsivity = KID(PoptT0index(n,1,2)).ResponsivityM1(PoptT0index(n,2,1),1:2); %[phase,radius] from fit and method 1
        NEP(NEPcounter).ddxdNqp = KID(PoptT0index(n,1,2)).ddxdNqp(1); %from fit of dx(Nqp)
        NEP(NEPcounter).etapb = etapb;
        %Obtain data from lifetime data array
        NEP(NEPcounter).dRdtheta = LifetimeData(PoptT0index(n,1,3),2);
        NEP(NEPcounter).tau = [LifetimeData(PoptT0index(n,1,3),11),LifetimeData(PoptT0index(n,1,3),11)]; %[tau(phase),tau(R)]
        
        
        %Display the match between noise, response and lifetime for user verification
        disp(['NEP calculation for KID ' num2str(KIDnumbers(n))])
        disp(['Noise data taken at P_{opt} = ',num2str(NEP(NEPcounter).ReadPower(1)),...
            ' (dBm) and T = ',num2str(NEP(NEPcounter).Temperature(1)),' (K)'])
        disp(['Response data taken at P = ',num2str(NEP(NEPcounter).ReadPower(2)),...
            ' (dBm) and T = ',num2str(NEP(NEPcounter).Temperature(2)),' (K)'])
        disp(['Lifetime data taken at P = ',num2str(NEP(NEPcounter).ReadPower(3)),...
            ' (dBm) and T = ',num2str(NEP(NEPcounter).Temperature(3)),' (K)'])
        
        %NEP Calculation
        NEP(NEPcounter).InternalPower = 10*log10(2/pi*(NEP(NEPcounter).Ql(1)^2/NEP(NEPcounter).Qc(1))*10.^(NEP(NEPcounter).ReadPower/10));
        NEP(NEPcounter).ResRingTime = NEP(NEPcounter).Ql(1)/pi/(NEP(NEPcounter).Fres*1e9);
        IFfreq2 = (2*pi*NEP(NEPcounter).FFTnoise(:,1)).^2;
        NEP(NEPcounter).RollOff = sqrt(1+IFfreq2*NEP(NEPcounter).tau(1)^2).*sqrt(1+IFfreq2*NEP(NEPcounter).ResRingTime^2);
        NEP(NEPcounter).dydPdark(1) = NEP(NEPcounter).Responsivity(1)*etapb*NEP(NEPcounter).tau(1)/NEP(NEPcounter).Delta; %(dtheta/dP) dark
        NEP(NEPcounter).dydPdark(2) = NEP(NEPcounter).dydPdark(1)*NEP(NEPcounter).dRdtheta; %(dR/dP) dark
        NEP(NEPcounter).darkNEP = NEP(NEPcounter).FFTnoise;
        NEP(NEPcounter).darkNEP(:,2) = 10.^(NEP(NEPcounter).darkNEP(:,2)/20)/abs(NEP(NEPcounter).dydPdark(1)).*NEP(NEPcounter).RollOff;
        NEP(NEPcounter).darkNEP(:,3) = 10.^(NEP(NEPcounter).darkNEP(:,3)/20)/abs(NEP(NEPcounter).dydPdark(2)).*NEP(NEPcounter).RollOff;
        
        %Saving the optimum NEP to file for each KID.
        NEPfile=[Mainpath,filesep,'KID_' num2str(KIDnumbers(n)),'_optimumNEP.csv'];
        WriteSRONcsv(NEPfile,NEP(NEPcounter).darkNEP,NEPheader,'%.6g');
    end
end

%==========================================================================
% Estimate alpha (kinetic inductance fraction) by fitting Fdesign vs Fres
% (measured)
%==========================================================================
if size(KIDParameters,1) > 1
    AlphaFitOptions = fitoptions('Method','NonlinearLeastSquares','Startpoint',[1,0]);
    AlphaFitType = fittype('a*x+b','options',AlphaFitOptions);  %%%fit algorithm
    AlphaFit = fit(KIDParameters(:,4)*1e9,KIDParameters(:,3),AlphaFitType); %F0=a*Fdes
    alpha = 1-sqrt(AlphaFit.a);
    disp(['The kinetic induction fraction of this chip is estimated to be: ',num2str(alpha)])
else
    AlphaFit=[];
    disp('Not enough KIDs to estimate the kinetic induction fraction')
end

%==========================================================================
% Export the most important parameters to csv file
%==========================================================================
NEPParameters = zeros(length(NEP),22);
for p=1:length(NEP)
    KIDindex = find(KIDnumbers == NEP(p).KIDnumber);
    
    NEPParameters(p,1) = NEP(p).KIDnumber;
    NEPParameters(p,2) = NEP(p).Temperature(1);
    NEPParameters(p,3) = NEP(p).Fres;
    NEPParameters(p,4) = KID(PoptT0index(KIDindex,1,2)).Fdesign;
    NEPParameters(p,5) = NEP(p).Ql(1);
    NEPParameters(p,6) = NEP(p).Qi(1);
    NEPParameters(p,7) = NEP(p).Qc(1);
    NEPParameters(p,8) = KID(PoptT0index(KIDindex,1,2)).Qdesign;
    NEPParameters(p,9) = NEP(p).ReadPower(1);
    NEPParameters(p,10) = NEP(p).InternalPower(1);
    
    NEPParameters(p,11) = NEP(p).tau(1);
    NEPParameters(p,12) = NEP(p).ddxdNqp;
    NEPParameters(p,13) = NEP(p).dRdtheta;
    NEPParameters(p,14) = (-1*NEP(p).ddxdNqp*NEP(p).etapb*NEP(p).tau(1)/NEP(p).Delta)^(-1); %dF0/dPdark*1/F0
    NEPParameters(p,15) = KID(PoptT0index(KIDindex,1,2)).Tc;
    NEPParameters(p,16) = NEP(p).Delta;
    NEPParameters(p,17) = KID(PoptT0index(KIDindex,1,2)).Thickness;
    NEPParameters(p,18) = KID(PoptT0index(KIDindex,1,2)).Area;
    NEPParameters(p,19) = NEP(p).etapb;
    NEPParameters(p,20) = feval(NOISE(PoptT0index(KIDindex,1,1)).FitS1kHz{PoptT0index(KIDindex,2,1),1},-40);
    
    NEPParameters(p,21) = min(NEP(p).darkNEP(:,2));
    NEPParameters(p,22) = min(NEP(p).darkNEP(:,3));
end

NEPparfile = [Mainpath,filesep,'NEPsummary.csv'];
NEPparheader = {'KIDID','T0','F0','Fdesign','Q','Qi','Qc','Qdesign',...
    'Popt [dBm]','Piopt [dBm]','tau [msec]',...
    'ddx/dNqp','dR/dtheta','ddx/dPdark',...
    'Tc [K]','Delta [J]','Film Thickness [um]','Active Area [um^2]','eta_pb'...
    'Sf1kHz40dBm','minNEPphase','minNEPamp'}';
WriteSRONcsv(NEPparfile,NEPParameters,NEPparheader,'%.6g');

%==========================================================================
% Plotting
%==========================================================================

NEPkids = zeros(NEPcounter,1);
for p=1:NEPcounter
    NEPkids(p)=NEP(p).KIDnumber;
end
[NEPcolors,~] = GenerateColorsFromMap(NEPcounter,'RainbowReinier');

%FIGURE 1: Shows the Phase and Amplitude noise as optimum power for each of
%the KIDs for which complete data was available.
figure(1)
clf

subplot(1,2,1) %Phase NEP
for p=1:NEPcounter
    loglog(NEP(p).darkNEP(:,1),NEP(p).darkNEP(:,2),'-','color',NEPcolors(p,:),'LineWidth',2)
    hold on
end
axis([1 5000 1e-19 1e-16])
grid on
xlabel('F [Hz]')
ylabel('NEP_{\theta} [W Hz^{-0.5}]')
hold off

subplot(1,2,2) %Amplitude NEP
for p=1:NEPcounter
    loglog(NEP(p).darkNEP(:,1),NEP(p).darkNEP(:,3),'-','color',NEPcolors(p,:),'LineWidth',2)
    hold on
end
axis([1 5000 1e-19 1e-16])
grid on
xlabel('F [Hz]')
ylabel('NEP_{A} [W Hz^{-0.5}]')
legend(num2str(NEPkids(:)),'Location','Best');
hold off

%Save the figure
Figfile=[Mainpath,filesep,'NEPcurves.fig'];
saveas(gcf,Figfile,'fig')

%==========================================================================
%FIGURE 2: Overview of KID parameters for all KIDs (function of KID ID)
figure(2)
clf

subplot(2,3,1); %Minimum NEP (phase and amplitude)
semilogy(NEPParameters(:,1),NEPParameters(:,22),'ko','MarkerSize',5,'MarkerFaceColor','k')
hold on
semilogy(NEPParameters(:,1),NEPParameters(:,21),'ro','MarkerSize',5,'MarkerFaceColor','r')
axis tight
xlabel('KID ID')
ylabel('NEP_x [W/Hz^{0.5}]')
legend('NEP_{A}','NEP_{\theta}')
hold off

subplot(2,3,2) %Lifetime
for p=1:NEPcounter
    plot(NEP(p).KIDnumber,NEP(p).tau(1)*1000,'ro','MarkerSize',5,'MarkerFaceColor','r')
    hold on
    plot(NEP(p).KIDnumber,NEP(p).tau(2)*1000,'bo','MarkerSize',5,'MarkerFaceColor','b')
end
axis tight
xlabel('KID ID')
ylabel('Lifetime [ms]')
legend('\tau_{\theta}','\tau_{A}')
hold off

subplot(2,3,3) %Noise at Popt and F=100Hz
for p=1:NEPcounter
    plot(NEP(p).KIDnumber,10*log10(NEP(p).MeanFreqNoise(3)),'ro','MarkerSize',5,'MarkerFaceColor','r')
    hold on
    plot(NEP(p).KIDnumber,NEP(p).MeanNoise(4),'ko','MarkerSize',5,'MarkerFaceColor','k')
end
axis tight
xlabel('KID ID')
ylabel('Noise [dB]')
title('Noise at P_{opt} at F = 100 Hz')
legend('S_f/f^2','S_A')
hold off

subplot(2,3,4) %Frequency Response
plot(NEPParameters(:,1),NEPParameters(:,12),'ko','MarkerSize',5,'MarkerFaceColor','k')
axis tight
xlabel('KID ID')
ylabel('ddx/dN_{qp}')
title('QP responsivity')
hold off

subplot(2,3,5) %dRdtheta
plot(NEPParameters(:,1),NEPParameters(:,13),'ko','MarkerSize',5,'MarkerFaceColor','k')
axis tight
xlabel('KID ID')
ylabel('dR/d\theta')
hold off

subplot(2,3,6) %Internal Power at Popt.
plot(NEPParameters(:,1),NEPParameters(:,10),'ko','MarkerSize',5,'MarkerFaceColor','k')
axis tight
xlabel('KID ID')
ylabel('P_{int,opt}')
hold off

%Save the figure
Figfile=[Mainpath,filesep,'NEPoverview1.fig'];
saveas(gcf,Figfile,'fig')

%==========================================================================
%FIGURE 3: Overview of KID relations
figure(3)
clf

subplot(2,3,1); %Minumum NEP as function of maximum Q
for p=1:NEPcounter
    loglog(NEP(p).Ql(2),min(NEP(p).darkNEP(:,3)),'ko','MarkerSize',5,'MarkerFaceColor','k')
    hold on
    loglog(NEP(p).Ql(2),min(NEP(p).darkNEP(:,2)),'ro','MarkerSize',5,'MarkerFaceColor','r')
end
axis tight
xlabel('max(Q_{meas})')
ylabel('NEP_x [W/Hz^{0.5}]')
legend('NEP_{A}','NEP_{\theta}')
hold off

subplot(2,3,2); %maximum Qi distribution
for p=1:NEPcounter
    semilogy(NEP(p).KIDnumber,NEP(p).Qi(2),'ko','MarkerSize',5,'MarkerFaceColor','k')
    hold on
end
axis tight
xlabel('KID ID')
ylabel('max(Q_i)')
hold off

subplot(2,3,4); %Optimum internal power as a function of Q
for p=1:NEPcounter
    semilogx(NEP(p).Ql(2),NEP(p).InternalPower,'ko','MarkerSize',5,'MarkerFaceColor','k')
    hold on
end
axis tight
xlabel('max(Q_{meas})')
ylabel('P_{int,opt}')
hold off

subplot(2,3,5); %Scatter in coupling Q
plot(NEPParameters(:,8),NEPParameters(:,7)./NEPParameters(:,8),'ko','MarkerSize',5,'MarkerFaceColor','k')
axis tight
xlabel('Q_{design}')
ylabel('Q_{c}/Q_{design}')
hold off

subplot(2,3,3); %Frequency change wrt design
if size(KIDParameters,1) > 1
    plot(KIDParameters(:,4)*1e9,KIDParameters(:,3),'ko','MarkerSize',5,'MarkerFaceColor','k')
    hold on
    plot(AlphaFit)
    axis tight
    xlabel('F_{design} [Hz]')
    ylabel('F_{res} [Hz]')
    title(['\alpha = ',num2str(alpha)])
    hold off
end

subplot(2,3,6); %Frequency change wrt design
if size(KIDParameters,1) > 1
    plot(KIDParameters(:,4)*1e9,(KIDParameters(:,3)-(AlphaFit.a*KIDParameters(:,4)*1e9+AlphaFit.b))/1e6,'ko','MarkerSize',5,'MarkerFaceColor','k')
    axis tight
    xlabel('F_{design} [Hz]')
    ylabel('F_{res} - F_{design}*\delta F_{res}/\delta F_{design} [MHz]')
    title('Scatter in F_{res}')
    hold off
end

%Save the figure
Figfile=[Mainpath,filesep,'NEPoverview2.fig'];
saveas(gcf,Figfile,'fig')

%==========================================================================
% Clean up, save workspace and close
%==========================================================================
save([Mainpath,filesep,'NEPcompare.mat'])

rmpath([pwd,filesep,'subroutines']);
end
