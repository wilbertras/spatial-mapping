function [NEPR,NEPtheta,NEPcross,Nqp,nqp,NEPGR,tau0] = getDARKNEP(dthetadN, dRdN, Stheta, SR, Scross, F, eta_pb, tau, taures, crosslevel, V, Delta, Tc)
% Calculates the experimental dark NEP and other params from the cross PSD
% Not sure if implemented ok wrt factor 2 Pieter
%
% INPUT
% single values unless otherwise specified
% dthetadN, dRdN:   phase and amplitude responsivities p[er qp (obtained using S21 analysis)%
% Stheta, SR:       Noise  in phase and amplitude (dBc/Hz), vector %
% F:                Frequency vector corresponding to the noise data (Hz)
% eta_pb:           pair breaking efficiency
% tau:              experimental qp lifetime (sec)
% taures:           resonator ring time (sec)
% crosslevel:       PSD noise level of the cross PSD (1/Hz)%
% V:                Alu volme (um^3)
% Delta:            Gap in (J)
% Tc:               critical temperature (K) - optional
%
% OUTPUT
% NEPR, NEPtheta    NEP in amplitude, phase (W/rt(Hz) 
% NEPcross          NEP cross estimate (W/rt(Hz) 
% Nqp               qp number
% nqp               qp density (um^-3)
% NEPGR             GR noise limited NEP base dupon tau and Nqp (W\rt(Hz))
% tau0              kaplan tau0 for experimental data, consistent with
% earlier publications. To compare with Kaplan theory, this value is wrong
% by factor 2

%%%%%%%% constants %%%%%%%%
e_c=1.602e-19;          %single electron charge
N0 = 1.72e10/e_c;       %density of states in /J/um^3 of aluminium
kb=1.3806503e-23;       %boltzmann contstant J/K
%%%%%%%%%%%%%%%%%%%%%%%%%
    
%since we need Delta and Tc I want to check consistency D=1.76*kTc
if nargin == 12
    Tc_check = Delta/(1.76*kb);
    if round(Tc*1e3) ~= round(1e3*Tc_check)
        error('Delta and Tc inconsistent')
    end  
end
%calcualte responsivities
dthetaPdark = dthetadN * eta_pb * tau / Delta;  %=dtheta/dN * dN/dP: dN/dP => PdV 2.48. 
% Note that this eqn. is correct fully when using the MEASURED lifetime (2 particle)%
dRPdark = dRdN * eta_pb * tau / Delta;
dcrossPdark = (dRdN*dthetadN)^0.5 * eta_pb * tau / Delta;

%calcualte NEP's (these are arrays)             PdV 2.47
NEPtheta    = (10.^(Stheta/20) / dthetaPdark) .* (1+(2*pi*F*tau).^2).^0.5 .* (1+(2*pi*F*taures).^2).^0.5;
NEPR        = (10.^(SR/20) / dRPdark)         .* (1+(2*pi*F*tau).^2).^0.5 .* (1+(2*pi*F*taures).^2).^0.5;
NEPcross    = (10.^(Scross/20) / dcrossPdark)     .* (1+(2*pi*F*tau).^2).^0.5 .* (1+(2*pi*F*taures).^2).^0.5;

% get Nqp and nqp using cross PSD level
%S_level,cross = 4 N tau dA/dN dtheta/dN            PdV 6.1 and 6.2   
Nqp   = crosslevel /(4*tau*dRdN*dthetadN);      %tau is the 2 particle MEASURED lifetime (checked)
nqp = Nqp/V; %qp per um^3

% get tau_0
% nqq tauqp = tau0 N0 (kbTc)^3/(2 D^2)          PdV 2.29
tau0 = nqp*2*tau*2*Delta^2/(N0*(kb*Tc)^3);  %tau_0 MUST be realted to the SINGLE particle lifetime; tau_sp=2xtau exp.%

% GR limited NEP
NEPGR     = 2*Delta/eta_pb * sqrt(Nqp/(2*tau));     %tau MUST be  the SINGLE particle lifetime; tau_sp=2xtau exp.%
% Note: By implementing this last 2 eqns correcfly, nothing has changed
% compared to our old method as the factors 2 compensate.
end