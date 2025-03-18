close all;
%First Run PostAbalysisLifetime to be able to grab the data
KIDnumber=3;
power=87;

KID_id=find(drdtheta(:,1)==KIDnumber);
powers=unique(drdtheta(KID_id,5));
npowers=length(powers);
thisP=find(power==drdtheta(KID_id,5)) ;
%%Fit only 1 power
%T=drdtheta(KID_id(thisP),4);%T in K
%tau_m=1e3*drdtheta(KID_id(thisP),11);%tau in msec
%Fit all powers at once
T=drdtheta(KID_id(:),4);%T in K
tau_m=1e3*drdtheta(KID_id(:),11);%tau in msec

T(isnan(tau_m))=[];
tau_m(isnan(tau_m))=[];
[T,ind]=sort(T);
tau_m=tau_m(ind);

%material and other params
Tcstart=1.25;%Al
maxlifetime=0.8;%limiting tau and lowest T (saturation regime)
tau0start=400e-9;%Kaplan Value
plot(T,tau_m,'o')


s = fitoptions('Method','NonlinearLeastSquares',...
    'Startpoint',[ Tcstart tau0start],...
    'Upper',[1.2*Tcstart 3*tau0start],...
    'Lower',[0.8*Tcstart 0.3*tau0start],...
    'MaxFunEvals',1000);
%Specify function to fit
%ftype = fittype('(maxtau*1000*(tau_0/41.2)*(Tc/x).^0.5.*exp(1.76*Tc/T))/((maxtau+1000*(tau_0/41.2)*(Tc/x).^0.5.*exp(1.76*Tc/x)))','options', s);

ftype = fittype('(maxtau*1000*(tau_0/41.2)*(Tc./x).^0.5.*exp(1.76*Tc./x))./(maxtau+1000*(tau_0/41.2)*(Tc./x).^0.5.*exp(1.76*Tc./x))','problem','maxtau','options', s);
%Perform the actual fit
[result]=fit(T,tau_m,ftype,'problem',maxlifetime)


s = fitoptions('Method','NonlinearLeastSquares',...
    'Startpoint',[ Tcstart tau0start],...
    'Upper',[1.1*Tcstart 2*tau0start],...
    'Lower',[0.9*Tcstart 0.5*tau0start],...
    'MaxFunEvals',1000);
ftype = fittype('(maxtau*1000*(tau_0/41.2)*(Tc./x).^0.5.*exp(1.76*Tc./x))./(maxtau+1000*(tau_0/41.2)*(Tc./x).^0.5.*exp(1.76*Tc./x))','problem','maxtau','options', s);
[result2]=fit(T,tau_m,ftype,'problem',maxlifetime)


%tau_real=1000*(result.tau_0/41.2)*(result.Tc./T).^0.5.*exp(1.76*result.Tc./T)
%tau=1000*(result.tau_0/41.2)*(result.Tc./T).^0.5.*exp(1.76*result.Tc./T)
%tau=1000*(resulttau.tau_0/41.2)*(resulttau.Tc./T).^0.5.*exp(1.76*resulttau.Tc./T);%Kaplan, valid for BCS SC: Delta/k_BTc=1.76
%tau_real=(maxtau.*tau)./(maxtau+tau);
tau_real=(maxlifetime*1000*(tau0start/41.2)*(Tcstart./T).^0.5.*exp(1.76*Tcstart./T))./(maxlifetime+1000*(tau0start/41.2)*(Tcstart./T).^0.5.*exp(1.76*Tcstart./T));
tau_realfitted=(result.maxtau*1000*(result.tau_0/41.2)*(result.Tc./T).^0.5.*exp(1.76*result.Tc./T))./...
    (result.maxtau+1000*(result.tau_0/41.2)*(result.Tc./T).^0.5.*exp(1.76*result.Tc./T));
tau_realfitted2=(result2.maxtau*1000*(result2.tau_0/41.2)*(result2.Tc./T).^0.5.*exp(1.76*result2.Tc./T))./...
    (result2.maxtau+1000*(result2.tau_0/41.2)*(result2.Tc./T).^0.5.*exp(1.76*result2.Tc./T));
%figure(123)
hold on;
plot(T,tau_real,'r');hold on;
plot(T,tau_realfitted,'b');hold on;
plot(T,tau_realfitted2,'k');hold on;
legend('data','start values',['fit tau_0=' num2str(result.tau_0) ' Tc=' num2str(result.Tc)],...
   ['fit tau_0=' num2str(result2.tau_0) ' Tc=' num2str(result2.Tc)] )