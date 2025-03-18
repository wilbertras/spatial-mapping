%Plots the optical efficinecy vs position of the MKIDs
%Positions is a 2xn matrix with xcoordinate, ycoordinate of the KIDs on the
%chip, referred to the chip center and sorted in Fres



addpath([pwd,filesep,'subroutines']);                           %Enable subroutines by adding path in search path.

hold on;
bla=colormap(cool(100));
bla=GenerateColorsFromMap(100,'RainbowReinier')
for n=1:length(KIDparam)
plot(Positions(n,2),Positions(n,3),'o','MarkerSize',50,...
'MarkerFaceColor',bla(round(100*KIDparam(n).optphaseeff/max([KIDparam.optphaseeff])),:),...
'MarkerEdgeColor',bla(round(100*KIDparam(n).optphaseeff/max([KIDparam.optphaseeff])),:));hold on;
text(Positions(n,2)-0.1,Positions(n,3),num2str(KIDparam(n).optphaseeff/max([KIDparam.optphaseeff]),3));
title(['Peak efficiency: ' num2str(max([KIDparam.optphaseeff]))]);
end


rmpath([pwd,filesep,'subroutines']);              