function [IndexPopt] = PowerDependenceS21

% This function can be run after RESPONSEanalysisV5 has finished analysing the
% data. It load the HybridS21_2 workspace from file.
% It then does 2 things:
% 1) ONLY IF A 'Popt.csv' IS PRESENT IN THE S21path: for each KID the index
%    in the KID struct array and KIDParameters array (these 2 indexes are
%    identical) will be determined.
%   [Optional:] if 'Popt.csv' is not found or if the KID has no listed Popt
%   in it, the minimum available power is assumed to be Popt.
%   NOTE: The 'Popt.csv' is produced by the NoiseAnalysisHybrids routine.
%   It has to be manually copied to the directory containing the S21 data
%   on which this script runs. The format of 'Popt.csv' is given below
%
% 2) Make a number of plots showing the power dependence of Q's and Fres.
%    As well as the (power,temperature) dependence of Qi and Fres if 2D
%    data is present.
%
%INPUT: All user defined in this function (see below)
%
%OUTPUT:
%   IndexPopt == A Nx3 array. Here N is the number of unique KIDs in the
%                specified directory. It is given by the KIDid/KIDnumber in
%                IndexPopt(:,1). In IndexPopt(:,2) the index is given
%                where the optimum power parameters are located for this
%                KID in the "KID" and "KIDParameters" variables produced by
%                HybridsS21_2.m. Finally IndexPopt(:,3) gives the optimum
%                power in dBm.
%
%SUBROUTINES:
% GenerateColorsFromMap (External)
% 	colormapstorage.mat (External)
%
%REQUIRED FILES:
% HybridsS21.mat (Produced by HybridsS21_2.m)
% Popt.csv (optional, produced by NoiseAnalysisHybrids.m)
%
%VERSION: 1.2
%   Contains ideas of the original Q_Power1_1.m
%   V1.0 (2012/08/30, RJ): Start writing of this script.
%   V1.2 (2012/11/07, RJ): Added clear statements to remove some useless
%                           variables from the stored mat file.
%
%DATE: August 30, 2012
%AUTHOR: Reinier Janssen

%==========================================================================
% User defined input
%==========================================================================
% Full directory containing the data and analysis output that has been
% processed by HybridS21_2. It is the combination of the following
% variables in the HybridS21_2 space: S21path = [ChipInfo.path,S21subdir,'\']
S21path = ...  %NOTE: Should end with '\'ChipInfo.path = [cd '/..']; %root path where data is, one higher than the scripts
    [cd '/../S21/Power/'];

%==========================================================================
%Standard settings to determine when to give warnings. Can be user modified
%==========================================================================
dT0limit = 0.1; %[K](Default = 0.1) For each KID at all powers the base temperature is looked up.
           %If the maximum and minimum deviate with more than dT0 a warning is given.
dPlimit = 0.5; %[dBm](Default = 0.5). If the optimum power found in Popt deviates 
           %more than dPlimit from the nearest power in the S21(T,P) sweep a warning is given.
PoptPmin = 0; %if set to 1, then Pmin with be taken as Popt if Popt is not found
           %(either by absence of or absence in 'Popt.csv'-file.)
%==========================================================================
% Setting some routine default values
%==========================================================================
format('long','e'); %Set display format of numbers to 7 digits
%Enable subroutines by adding path in search path.
addpath([pwd,filesep,'subroutines']);
close all; %close all open plots to remove clutter.

%==========================================================================
% Load the workspace of HybridS21_2
%==========================================================================
S21variables={'KID','KIDParameters','KIDnumbers','ChipInfo'};
load([S21path,'ResponseS21.mat'],S21variables{:}) 
%This matlab workspace storage should at least contain the variables:
% KID, KIDParameters,KIDnumbers,S21path (overwrites current), ChipInfo


%==========================================================================
% Check if a Popt.csv file is present
%==========================================================================
PoptFile = [ChipInfo.path,filesep,'Popt.csv'];
fid = fopen(PoptFile);
if fid == -1
    %There is no Popt file. Do own Popt analysis internally
    DataPopt = [];
else
    fclose(fid);
    %There is a Popt file. Read it and use it.
    fprintf('Found Popt.csv, this will be used to indicate Popt.\n')
    [~,DataPopt] = ReadSRONcsvV2(PoptFile,'',0);
    %DataPopt will be Nx3. Where dim 2 contains [T,KIDid,|Popt [dBm]|]
end

%==========================================================================
% For all KIDs sort their stuff with increasing power
%==========================================================================
IndexPsort = cell(length(KIDnumbers),1); %initialize cell array.
for p=1:length(KIDnumbers)
    KIDlocations = find(KIDParameters(:,1) == KIDnumbers(p));
    [~,IndexPsort{p,1}]=sort(KIDParameters(KIDlocations,9)); %Sort on Pread and save the indexes.
    IndexPsort{p,1} = KIDlocations(IndexPsort{p,1}(:)); %Recalibrate the indexes to the full KID(Parameters) variable
end

%==========================================================================
% Loop over all KIDs to find the index of Popt
%==========================================================================
%Initialize the output variable
IndexPopt = zeros(length(KIDnumbers),3) -1; %set to -1 for error catching.
IndexPopt(:,3) = 1e6+1; %set actual Popt to 1.000.001 for error catching.

%Initialize some additional parameters
T0bounds = zeros(length(KIDnumbers),2); % [Min,Max] base temperature found in the power sweep

for p=1:length(KIDnumbers)
    KIDid = KIDnumbers(p); %Determine KIDid for this loop
    IndexPopt(p,1) = KIDid; %Put KID number in the output variable
    KIDlocations = find(KIDParameters(:,1) == KIDid); %find location of this KIDid in the larger variables.
    
    %Calculate min and max of Base Temperature
    T0bounds(p,1) = min(KIDParameters(KIDlocations,2));
    T0bounds(p,2) = max(KIDParameters(KIDlocations,2));
    
    if T0bounds(p,2)-T0bounds(p,1) > dT0limit %Large fluctuation in base temperature
        %Print a warning
        fprintf(['WARNING [plotPdependence]: For KID ',num2str(KIDid),...
            ' a large variation in base temperatures is detected\n'])
    end
    
    %======================================================================
    % Determination of the Popt and its index in KID
    %======================================================================

    %USING Popt FILE
    %See if Popt can be found in the Popt file.
    if isempty(DataPopt) == 0
        %There is Popt data from a file. Search it to find an identical
        %KIDid and a base temperature between 0.9*T0min and 1.1*T0max
        KIDinDataPopt = find(DataPopt(:,2) == KIDid & ...
            DataPopt(:,1) <= 1.1*T0bounds(p,2) & DataPopt(:,1) >= 0.9*T0bounds(p,1));
        if isempty(KIDinDataPopt)
            %KID (at these base temperatures) is not in DataPopt
        else
            if length(KIDinDataPopt) >= 2
                %KID (at these base temperatures) is more then once in
                %DataPopt. This is strange. Give warning and pick the lowest
                %temperature
                fprintf(['WARNING [plotPdependence]: For KID ',num2str(KIDid),...
                    ' multiple Popt found within Popt.csv. Taking lowest T.\n'])
                [~,IndexTlow]=min(DataPopt(KIDinDataPopt,1));
                IndexPopt(p,3) = -1*DataPopt(KIDinDataPopt(IndexTlow),3);
            else
                %As regular (expected) only a single match.
                IndexPopt(p,3) = -1*DataPopt(KIDinDataPopt,3);
            end
            [dP,IndexPopt(p,2)] = min(abs(IndexPopt(p,3)-KIDParameters(KIDlocations,9)));
            IndexPopt(p,2) = KIDlocations(IndexPopt(p,2));
            if dP > dPlimit %If the difference is bigger that specified
                fprintf(['WARNING [plotPdependence]: For KID ',num2str(KIDid),...
                    ' Popt from Popt.csv is more that 0.5dB away from nearest measured value.\n'])
            end
        end
    end

    %INTERNAL Popt DETERMINATION (IF PoptPmin == 1)
    %If Popt has not been found in the Popt file, determine Popt internally
    if IndexPopt(p,2) < 0 && PoptPmin == 1
        %This is a very dumb statement, but it works for now. (Basically Popt = Pmin)
        fprintf(['WARNING [plotPdependence]: For KID ',num2str(KIDid),...
            ' taking minimum power as Popt.\n'])
        IndexPopt(p,2) = IndexPsort{p}(1); %Lowest Power will be at the first sorted index.
        IndexPopt(p,3) = KIDParameters(IndexPopt(p,2),9);
    end
end

%==========================================================================
%Create the legend of the following figures
%==========================================================================
LegendText = cell(length(KIDnumbers),1); %Initialize the legend
for p=1:length(KIDnumbers)
    KIDid = KIDnumbers(p);
    %Construct the legend text
    LegendText{p} = ['KID ',num2str(KIDid),...
        ' @ ',num2str(T0bounds(p,1)),' \leq T_{base} \leq ',num2str(T0bounds(p,2)),' [K]'];
end

%==========================================================================
% Make a figure showing the power dependence of each KID
%==========================================================================
for p=1:length(KIDnumbers)
    [Pcolors,~]=GenerateColorsFromMap(length(IndexPsort{p,1}(:)),'RainbowReinier');
    
    figure(KIDnumbers(p)*10) %multiply by 10 so any figures of HybridsS21 are not overwritten.
    clf
    
    %SUBPLOT(2,2,1) |S21| [magnitude] as a function of Power
    subplot(2,2,1)
    hold on
    for q=1:length(IndexPsort{p,1}(:))
        plot(KID(IndexPsort{p,1}(q)).S21data{1}(:,1),KID(IndexPsort{p,1}(q)).S21data{1}(:,2),...
            '.-','color',Pcolors(q,:),'MarkerSize',4)
    end
    if IndexPopt(p,2) > 0 %Only if Popt has been identified.
        plot(KID(IndexPopt(p,2)).Fres(1),10.^(KID(IndexPopt(p,2)).S21min(1)/20),'kd',...
            'MarkerSize',8,'LineWidth',1,'MarkerFaceColor','r')
    end
    axis tight
    %make legend with title
    h=legend(num2str(KIDParameters(IndexPsort{p,1}(:),9)));
    v = get(h,'title');
    set(v,'string','P_{read} [dBm]');
    %Set general title and axis labels
    title('|S21|(F,P)')
    xlabel('F [GHz]')
    ylabel('|S21| [magn.]')
    hold off
    
    %SUBPLOT(2,2,2) Q at Tbase as function of power
    subplot(2,2,2)
    semilogy(KIDParameters(IndexPsort{p,1}(:),9),KIDParameters(IndexPsort{p,1}(:),5),'k.-',... %Loaded Q
        'MarkerSize',6,'LineWidth',1)
    hold on %Hold on after first semilogy to keep y-axis in log scale
    semilogy(KIDParameters(IndexPsort{p,1}(:),9),KIDParameters(IndexPsort{p,1}(:),7),'r.-',... %Coupling Q
        'MarkerSize',6,'LineWidth',1)
    semilogy(KIDParameters(IndexPsort{p,1}(:),9),KIDParameters(IndexPsort{p,1}(:),6),'b.-',... %Internal Q
        'MarkerSize',6,'LineWidth',1)
    if IndexPopt(p,2) > 0 %Only if Popt has been identified.
        %Popt POINTS
        semilogy(IndexPopt(p,3),KIDParameters(IndexPopt(p,2),5),'ko',... %Loaded Q @Popt
            'MarkerSize',8,'LineWidth',1)
        semilogy(IndexPopt(p,3),KIDParameters(IndexPopt(p,2),7),'ro',... %Coupling Q @Popt
            'MarkerSize',8,'LineWidth',1)
        semilogy(IndexPopt(p,3),KIDParameters(IndexPopt(p,2),6),'bo',... %Internal Q @Popt
            'MarkerSize',8,'LineWidth',1)
    end
    %Make a nice layout
    legend('Q_l','Q_c','Q_i')
    xlabel('P @ KID [dBm]')
    ylabel('Q (diamonds at Popt)')
    title(LegendText{p})
    hold off
    
    %SUBPLOT(2,2,3) Fres as a function of temperature (x-axis) and power (color)
    subplot(2,2,3)
    hold on
    for q=1:length(IndexPsort{p,1}(:))
        plot(KID(IndexPsort{p,1}(q)).Temperature,KID(IndexPsort{p,1}(q)).Fres,... %Fres
            '.-','color',Pcolors(q,:),'MarkerSize',6)
    end
    if IndexPopt(p,2) > 0 %Only if Popt has been identified.
        plot(KID(IndexPopt(p,2)).Temperature,KID(IndexPopt(p,2)).Fres,'k.-',... %Fres
            'MarkerSize',10,'LineWidth',2)
        title(['Popt of ',num2str(IndexPopt(p,3)),' [dBm] given by black line'])
    end
    %make legend with title
    h=legend(num2str(KIDParameters(IndexPsort{p,1}(:),9)));
    v = get(h,'title');
    set(v,'string','P_{read} [dBm]');
    %Set axis labels
    xlabel('T [K]')
    ylabel('F_{res} [GHz]')
    hold off
    
    %SUBPLOT(2,2,4) Qi as a function of temperature (x-axis) and power (color)
    AX4 = subplot(2,2,4);
    set(AX4,'Box','on','Yscale','log')
    hold on
    for q=1:length(IndexPsort{p,1}(:))
        plot(KID(IndexPsort{p,1}(q)).Temperature,KID(IndexPsort{p,1}(q)).Ql,... %Qi
            '.-','color',Pcolors(q,:),'MarkerSize',6)
    end
    if IndexPopt(p,2) > 0 %Only if Popt has been identified.
        plot(KID(IndexPopt(p,2)).Temperature,KID(IndexPopt(p,2)).Qi,'k.-',... %Qi
            'MarkerSize',10,'LineWidth',2)
        title(['Popt of ',num2str(IndexPopt(p,3)),' [dBm] given by black line'])
    end
    %make legend with title
    legend(num2str(KIDParameters(IndexPsort{p,1}(:),9)));
    legendtitle('P_{read} [dBm]')
    %v = get(h,'title');
    %set(v,'string','P_{read} [dBm]');
    %Set axis labels
    xlabel('T [K]')
    ylabel('Qi')
    hold off
    
    %SAVE the figure
    Figfile=[S21path,'KID',num2str(KIDnumbers(p)),'_PowerDependence.fig'];
    saveas(gcf,Figfile,'fig')
end

%==========================================================================
%Make an overview plot of Q(P) and Qi(P) for all KIDs in KIDnumbers
%==========================================================================
[KIDcolors,~] = GenerateColorsFromMap(length(KIDnumbers),'RainbowReinier'); %Generate some colors

%Determine the limits for the Y axes
Qlims(1,1) = min(KIDParameters(:,5));
Qlims(2,1) = max(KIDParameters(:,5));
Qlims(1,2) = min(KIDParameters(:,6));
Qlims(2,2) = max(KIDParameters(:,6));
Ylims(1,1) = min(Qlims(1,:));
Ylims(2,1) = max(Qlims(2,:));

figure(1e6)
clf
%SUBPLOT(1,2,1): Q as function of Pint
AX1 = subplot(1,2,1);
set(AX1,'Box','on','Yscale','log','Ylim',Ylims)
hold on
for p=1:length(KIDnumbers)
    KIDid = KIDnumbers(p);
    KIDlocation = find(KIDParameters(:,1) == KIDid);
    %Plot the data
    plot(KIDParameters(KIDlocation,10),KIDParameters(KIDlocation,5),'.',...
        'color',KIDcolors(p,:),'MarkerSize',15)
end
ylabel('Q')
xlabel('P_{int} [dBm]')
legend(LegendText)
hold off

%SUBPLOT(1,2,2): Q as function of Pint
AX2 = subplot(1,2,2);
set(AX2,'Box','on','Yscale','log','Ylim',Ylims)
hold on
for p=1:length(KIDnumbers)
    KIDid = KIDnumbers(p);
    KIDlocation = find(KIDParameters(:,1) == KIDid);
    plot(KIDParameters(KIDlocation,10),KIDParameters(KIDlocation,6),'.',...
        'color',KIDcolors(p,:),'MarkerSize',15)
end
ylabel('Q_{i}')
xlabel('P_{int} [dBm]')
hold off

%SAVE the figure
Figfile=[S21path,'CompareQvsPint.fig'];
saveas(gcf,Figfile,'fig');

%==========================================================================
% Workspace saving and wrap-up
%==========================================================================
%Clean up nonsense variables
clear AX1 AX2 AX4 Figfile KIDcolors KIDid KIDlocation KIDlocations LegendText
clear Pcolors Qlims T0bounds Ylims ans fid h p q v
% Save the matlab workspace
save([S21path,'PowerDependenceS21.mat']);
%Remove Path containing Subroutines from pathlist.
rmpath([pwd,filesep,'subroutines']);
%==========================================================================
%END OF plotPdependence SUBROUTINE
%==========================================================================
end