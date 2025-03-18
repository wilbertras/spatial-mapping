
function  [result,fitted] = FitS21simple(fdata,Fres,Qf,S21min,no)
% perfoms Lorenztain fit in log space over no bandwidths
% Data should contain 2 cols F[Hz], S21 dB.
% Fres is res freq in GHz, 
% Qf is initial Q factor and 
% S21min is original S21min in magn. (NOT in dB)%
% no = width in resonator bandwidths to fit over
% result is a struct with fields 
% result.Q;
% result.Fr;
% result.Smin;
% this is a much more simple version of the more general function
% FitS21main5.m

%format('long','e');
%fdata(:,2)=10.^(fdata(:,2)/20);%converting data to magnitude. NB S21 is in magn!!!
bandwidth=Fres/Qf; %resonator bandwidth

%a simple lorentzian in log space
%find range
minfitindex=find(fdata(:,1)>(Fres-no*bandwidth/2),1);
maxfitindex=find(fdata(:,1)>(Fres+no*bandwidth/2),1);
if maxfitindex-minfitindex>10
    x=fdata(minfitindex:maxfitindex,1);
    y=fdata(minfitindex:maxfitindex,2);
else
    x=(fdata(:,1));
    y=fdata(:,2);
    disp('bandwidthguess for fit far off, whole S21 curve used for Fo and S21');
end

s = fitoptions('Method','NonlinearLeastSquares', 'Startpoint',[Fres Qf S21min], 'MaxFunEvals',100);

%Fit the log10 of a  Lorentzian in power space 
ftype = fittype('10*log10(1-(1-Smin^2)/(1+(2*Q*(x-Fr)/Fr)^2))','options', s);
[result]=fit(x,y,ftype);
%  figure(100);
%  hold on;
%  plot(x,y,'o');
%  plot(x,10*log10(1-(1-S21min.^2)./(1+(2*Qf*(x-Fres)/Fres).^2)),'--');
%  plot(result);
%  legend('data','guess','fit')
%  
%  close(100);

fitted=10*log10(1-(1-result.Smin^2)./(1+(2*result.Q*((fdata(:,1))-result.Fr)/result.Fr).^2));

end