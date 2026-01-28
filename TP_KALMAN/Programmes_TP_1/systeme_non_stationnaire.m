%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Système non stationnaire                                                %
%                                                                         %
% Cécile DURIEU, 25 janvier 2025                                          %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Système analogique
Ka1 = 1; Ka2 = 0.8; % valeurs extrêmes du gain
tau1 = 0.1; tau2 = 0.08; % valeurs extrêmes de la constante de temps

% Variation linéaire du gain et de la constante de temps
% Ka = [Ka1*ones(N/2,1); Ka2*ones(N/2,1)]; 
% Ka(N/2-100:N/2+100) = Ka1:(Ka2-Ka1)/200:K2;
% tau = [tau1*ones(N/2,1); tau2*ones(N/2,1)];
% tau(N/2-100:N/2+100) = tau1:(tau2-tau1)/200:tau2;

% Saut du gain et de la constante de temps
Ka = [Ka1*ones(N/2,1); Ka2*ones(N/2,1)];
tau = [tau1*ones(N/2,1); tau2*ones(N/2,1)];

% Système échantillonné bloqué : évolution des paramètres en fonction du temps
Te = 0.1;  Fe = 1/Te; z = tf('z',Te);
a = -exp(-Te./tau); b = Ka.*(1+a); 

% Instants d'acquisition
t = (0:N-1)*Te; t = t(:);