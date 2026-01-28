%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Système stationnaire                                                    %
%                                                                         %
% Cécile DURIEU, 25 janvier 2025                                          %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Système analogique
Ka = 1; tau = 0.1;
p = tf('s'); Ha = Ka/(1+tau*p);

% Système échantillonné
Te = 0.1; Fe = 1/Te;
a = -exp(-Te/tau); b = Ka*(1+a);
z = tf('z',Te); sys = b*z^-1/(1+a*z^-1);
Hn = c2d(Ha,Te);

% Instants d'acquisition
t = (0:N-1)*Te; t = t(:);