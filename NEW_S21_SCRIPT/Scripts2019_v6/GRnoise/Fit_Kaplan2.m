function [result,tau_fitted, nqp_fit, nqp_tau0] = Fit_Kaplan2(T, tau_m, lowT, Tc)
% fits Kaplan theory to lifetime (T) in log space to get fit reliable for
% the T dependence
% uses Delta = 1.76 kTc throughout
% also calcualtes the GR linited NEP based upon the resulting Tc and
% measured lifetime, if V is given as input
% INPUT
% tau:      qp lifetime in (sec)
% T:        T in K
% optional:
% V:        Alu volme (um^3). NEP calcualted if provided%
% eta_pb:	pair breaking efficiency, set to 0.4 if not given
% lowT:     T below which lifetime gets saturated, best to take a big margin (defauls = 0.2K)%
%           I found that a lowT range + fitting only the highT range, but including a
%           maximum tau value gives the most consistent results.
% Tc        Measured Tc, if supplied, Tc will NOT be fitted
% OUTPUT
% result        fitobject with fields:
%  .tau_0       correct Kaplan tau_0 (sec)
%               note that this refers to the single particle interaction
%               time, NOT our old wrong method. See Pieter's note%
%  .Tc          Criticat T in (K)
%  .maxtau      (sec) limiting lifetime at low T
% tau_fitted    fitted lifetime (sec) with same size as T
% nqp_fitted       	resulting qp number / um^3 based upon fitted tau, and fitted params tau_0, maxtau and fitted Tc, assuming Delta = 1.76 kbTc%
% nqp_tau0      	resulting qp number  based upon measured tau_m and fitted params tau_0 maxtau and fitted Tc, assuming Delta = 1.76 kbTc%%
% NEPGR         NEP based upon measured tau and fitted Tc, assuming Delta = 1.76 kbTc, with same size as T%


%%%%%%%% constants %%%%%%%%
e_c=1.602e-19;          %single electron charge
N0 = 1.72e10/e_c;       %density of states in /J/um^3 of aluminium
kb=1.3806503e-23;       %boltzmann contstant J/K
%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
    error('Fit_Kaplan: not enough input arguments')
end

if nargin < 3
    lowT = 0.2;
end

if length(T) ~= length(tau_m)
    error('Fit_Kaplan: unequal length T and tau');
elseif length(T) < 3
    disp('Warning: Kaplan fit ignored due to lack of data')
    result.tau_0 = NaN;
    result.Tc = NaN;
    result.maxtau = NaN;
    tau_fitted = zeros(size(T));
end

% Fit F range
Ri = T > lowT;
% catch missing values
gff = ~isnan(tau_m);
%defin fit data  
FRi = Ri & gff;

% start values for the fit
tau0start   = 300e-9;%Kaplan Value
if ~isempty(tau_m(T<lowT & gff))
    maxtaustart = mean(tau_m(T<lowT & gff));
else
    maxtaustart = 1e-3;
end

if nargin < 4 %allow Tc fit
    Tcstart     = 1.25;%Al
    % fitoptions
    s = fitoptions('Method','NonlinearLeastSquares',...
        'Startpoint',   [ Tcstart       maxtaustart     tau0start],...
        'Upper',        [1.2*Tcstart    2*maxtaustart   3*tau0start],...
        'Lower',        [0.8*Tcstart    0.5*maxtaustart 0.33*tau0start],...
        'MaxFunEvals',1000);
else
    s = fitoptions('Method','NonlinearLeastSquares',...
        'Startpoint',   [ Tc       maxtaustart     tau0start],...
        'Upper',        [1.0*Tc    2*maxtaustart   5*tau0start],...
        'Lower',        [1.0*Tc    0.5*maxtaustart 0.2*tau0start],...
        'MaxFunEvals',1000); %not using problem parameter here as this will make code more bulky
end
% Eqn 5.2 PdV: using Delta = 1.76 kTc becomes: 
% tau(sec) = tau_0 * 0.0243 * sqrt(Tc/T) * exp(1.76*Tc/T) with tau the
% SINGLE particle recombination time, tau=2*tau_measured!
% and adding taumax: tau = tau*maxtau/(maxtau + tau)

ftype = fittype('log10(maxtau*0.0243*tau_0*(Tc./x).^0.5.*exp(1.76*Tc./x) ./ (maxtau + 0.0243*tau_0*(Tc./x).^0.5.*exp(1.76*Tc./x))) ','options', s);
%coeffnames(ftype) %uncomment to get coefficient order

%Perform the actual fit using SINGLE particle tau as input
% SINGLE particle recombination time, tau=2*tau_measured!
Y = log10(2*tau_m(FRi));
[result]=fit(T(FRi),Y,ftype);

% get data for easy plotting in sec
% tau_fitted is experimental one, =tau/2
tau_fitted  = 0.5*(result.maxtau * 0.0243*result.tau_0  * (result.Tc./T).^0.5 .*exp(1.76*result.Tc./T)   ./ (result.maxtau + 0.0243*result.tau_0   * (result.Tc./T).^0.5 .*exp(1.76*result.Tc./T)));

Delta = 1.76*kb*result.Tc;
nqp_fit = result.tau_0 * N0 * (kb*result.Tc)^3 ./(tau_fitted*2*Delta^2); %qp concentration per um^3
nqp_tau0 = result.tau_0 * N0 * (kb*result.Tc)^3 ./(tau_m*2*Delta^2); %qp concentration per um^3


% %help figure in msec
% figure(123)
% subplot(1,2,1) %raw fit
% plot(T(FRi),Y,'bo','MarkerFaceColor','b');hold on
% plot(result,'b');
% 
% subplot(1,2,2) %result plot
% semilogy(T,tau_m,'bo');hold on
% semilogy(T(FRi),tau_m(FRi),'bo','MarkerFaceColor','b');hold on
% plot(T,tau_fitted,'b');axis tight
% legend('data','fitted',['fit tau_0=' num2str(result.tau_0) ' Tc=' num2str(result.Tc)] )
end