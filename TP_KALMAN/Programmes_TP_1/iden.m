%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         % 
% Identification d'un système du 1er ordre échantillonné avec :           %
%   - la méthode du modèle (parallèle) / méthode à erreur de sortie,      %
%     implémentée sous Matlab avec la fonction oe (output error),         %
%   - la méthode des moindres carrés (MC) / méthode du modèle             %
%     série-parallèle / méthode arx, implémentée sous Matlab avec la      %
%     fonction oe.                                                        %
%                                                                         %
% Les points étudiés sont :                                               %
%   - comparaison des résultats obtenus avec l'implémentation de la       %
%     méthode des MC et la fonction méthode arx,                          %
%   - comparaison des résultats obtenus avec la méthode du modèle         %
%     parallèle et la méthode du modèle série-parallèle,                  %
%   - effet du bruit de mesure : étude statistique,                       %
%   - suivi de paramètres avec la méthode des MCR et facteur d'oubli.     % 
%                                                                         %  
% Cécile Durieu, 28 février 2018                                          % 
% modifié le 16 février 2019, 28 janvier 2024 et 22 janvier 2025          %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
clear all; close all; clc;
choix = menu('Choisir la partie du programme à exécuter',...
    '1/ Comparaison méthode des MC et modèle ARX',...
    '2/ Comparaison méthode du modèle et méthode des MC : étude statistique',...
    '3/ MCR : illustration du fonctionnement',...
    '4/ MCR avec facteur d''oubli : illustration du fonctionnement',...
    '5/ MCR avec facteur d''oubli : étude statistique avec système stationnaire',...
    '6/ MCR avec facteur d''oubli : étude statistique avec système non stationnaire');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         % 
% Comparaison des résultats obtenus avec la méthode des MC et la fonction %
% ARX de Matlab                                                           %
%                                                                         % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if choix == 1;
N = 100; % nombre d'échantillons pour l'estimation
 
% Système
systeme_stationnaire;
% Signal d'entrée
entree;
% Signal de sortie du système
RSB = input(' rapport signal sur bruit (RSB) : '); if RSB == 0; RSB = 1e-16; end;
s = lsim(sys,e); s = s(:); Ps = sum(s.^2)/N; % sortie non bruitée
Pb = Ps/RSB; sigma_b = sqrt(Pb); % bruit de sortie
s_b = s + sigma_b*randn(N,1); % sortie bruitée
% Tracé des signaux
figure; stairs(t,e,'r'); hold on; plot(t,s,'b--',t,s_b,'b'); title('Signaux'); 
legend('entrée','sortie','sortie bruitée'); xlabel('temps (s)');
axis([0 N*Te -Ka*max(abs(e))-0.5 Ka*max(abs(e))+0.5]); grid on; pause
 
% Estimation avec la méthode des MC
x_mc = s_b(2:N); C = [-s_b(1:N-1) e(1:N-1)];
est = inv(C.'*C)*C.'*x_mc; % ou est = (C.'*C)\C.'*x_mc;
mod_mc = est(2)*z^-1/(1+est(1)*z^-1);
a_est_mc = est(1); b_est_mc = est(2);

% Estimation avec la méthode ARX
mod_arx = arx([s_b e],[1 1 1]);
a_est_arx = mod_arx.a(2); b_est_arx = mod_arx.b(2);

% Affichage des résultats
text       = [' système :          a = ',num2str(a),...
                           '       b = ',num2str(b)];  
text_mc    = [' méthode des MC :   a = ',num2str(a_est_mc),...
                           '       b = ',num2str(b_est_mc)];
text_arx   = [' méthode ARX :      a = ',num2str(a_est_arx),...
                           '       b = ',num2str(b_est_arx)];
text_delta = [' écart :          d_a = ',num2str(a_est_mc-a_est_arx),...
                              '  d_b = ',num2str(b_est_mc-b_est_arx)];
disp(' '); disp(text); disp(text_mc); disp(text_arx); disp(text_delta);

end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         % 
% Comparaison des résultats obtenus avec la méthode du modèle et la       %
% méthode des MC : étude statistique                                      %
%                                                                         % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if choix == 2;
N = 100; % nombre d'échantillons pour l'estimation
I = 1000; % nombre de simulations

% Système
systeme_stationnaire;
% Signal d'entrée
entree;
% Signal de sortie du système
RSB = input(' rapport signal sur bruit (RSB) : '); if RSB == 0; RSB = 1e-16; end;
s = lsim(sys,e); s = s(:); Ps = sum(s.^2)/N; % sortie non bruitée
Pb = Ps/RSB; sigma_b = sqrt(Pb); % bruit de sortie
s_b = s + sigma_b*randn(N,1); % sortie bruitée
figure; stairs(t,e,'r'); hold on; plot(t,s,'b--',t,s_b,'b');  title('Signaux');
legend('entrée','sortie','sortie bruitée'); xlabel('temps (s)');
axis([0 t(end)+Te -Ka*max(abs(e))-0.5 Ka*max(abs(e))+0.5]);  grid on; pause
     
A_est_oe = []; B_est_oe = []; Crit_oe = []; Pui_er_oe = [];
A_est_mc = []; B_est_mc = []; Crit_mc = []; Pui_er_mc = [];
for i = 1:I
    s_b = s + sigma_b*randn(N,1); % sortie bruitée
    % Estimation avec la méthode du modèle
    mod_oe = oe([s_b e],[1 1 1]); s_oe = lsim(mod_oe,e);
    A_est_oe = [A_est_oe; mod_oe.f(2)]; B_est_oe = [B_est_oe; mod_oe.b(2)];
    Crit_oe = [Crit_oe; sum((s_oe-s_b).^2)];
    Pui_er_oe = [Pui_er_oe; sum((s_oe-s).^2)];
    % Estimation avec la méthode des MC
    X_mc = s_b(2:N); C = [-s_b(1:N-1) e(1:N-1) ]; est = inv(C.'*C)*C.'*X_mc;
    mod_mc = est(2)*z^-1/(1+est(1)*z^-1); s_mc = lsim(mod_mc,e);
    s_pred_mc = [0; C*est];
    A_est_mc = [A_est_mc; est(1)]; B_est_mc = [B_est_mc; est(2)];
    Crit_mc = [Crit_mc; sum((s_pred_mc-s_b).^2)];
    Pui_er_mc = [Pui_er_mc; sum((s_mc-s).^2)];
end;
 
% Visualisation de la valeur estimée des paramètres
figure;
subplot(121); plot(1:I,a*ones(I,1),'r',1:I,A_est_mc,'g',1:I,A_est_oe,'y');
title('a'); legend('valeur exacte','mc','oe'); xlabel('# de la simulation');
axis([0 I min([A_est_mc; A_est_oe]) max([A_est_mc; A_est_oe])]); grid on;
subplot(122); plot(1:I,b*ones(I,1),'r',1:I,B_est_mc,'g',1:I,B_est_oe,'y');
title('b'); legend('valeur exacte','mc','oe'); xlabel('# de la simulation');
axis([0 I min([B_est_mc; B_est_oe]) max([B_est_mc; B_est_oe])]); grid on; pause;

moy_a_oe = sum(A_est_oe)/I; biais_a_oe = sum(A_est_oe)/I - a; 
std_a_oe = sqrt(sum((A_est_oe-moy_a_oe).^2)/I);
eqm_a_oe = sum((A_est_oe-a).^2)/I;
moy_b_oe = sum(B_est_oe)/I; biais_b_oe = sum(B_est_oe)/I - b; 
std_b_oe = sqrt(sum((B_est_oe-moy_b_oe).^2)/I);
eqm_b_oe = sum((B_est_oe-b).^2)/I;
crit_oe = sum(Crit_oe)/(N*I);
pui_er_oe = sum((Pui_er_oe))/(N*I);

moy_a_mc = sum(A_est_mc)/I; biais_a_mc = sum(A_est_mc)/I - a; 
std_a_mc = sqrt(sum((A_est_mc-moy_a_mc).^2)/I);
eqm_a_mc = sum((A_est_mc-a).^2)/I;
moy_b_mc = sum(B_est_mc)/I; biais_b_mc = sum(B_est_mc)/I - b; 
std_b_mc = sqrt(sum((B_est_mc-moy_b_mc).^2)/I);
eqm_b_mc = sum((B_est_mc-b).^2)/I;
crit_mc = sum(Crit_mc)/(N*I);
pui_er_mc = sum((Pui_er_mc))/(N*I);

% Visualisation du dernier relevé et de la sortie estimée
figure; plot(t,s,'b--',t,s_b,'b',t,s_mc,'g',t,s_oe,'y');
title('Comparaison des sorties pour le dernier relevé');
legend('sortie','sortie bruitée','sortie mc','sortie oe'); xlabel('temps (s)'); 
axis([0 t(end)+Te -Ka*max(abs(e))-0.5 Ka*max(abs(e))+0.5]); grid on; pause

% Visualisation du dernier relevé, de la sortie prédite et de la sortie prédite
s_pred_mc = [0; C*est];
figure; plot(t,s,'b--',t,s_b,'b',t,s_mc,'g',t,s_pred_mc,'y');
title('Comparaison des sorties pour le dernier relevé');
legend('sortie','sortie bruitée','sortie mc','sortie prédite'); xlabel('temps (s)');
axis([0 t(end)++Te -Ka*max(abs(e))-0.5 Ka*max(abs(e))+0.5]); grid on; pause

% Réponse unitaire avec les paramètres du dernier relevé
N = 21; t = (0:N-1)*Te;
e = ones(N,1); ru = lsim(sys,e);
ru_oe = lsim(mod_oe,e); ru_mc = lsim(mod_mc,e);
figure; plot(t,ru,'r',t,ru_mc,'g',t,ru_oe,'y');
legend('ru','ru mc','ru oe');
xlabel('temps (s)'); title('Comparaison des RU pour le dernier relevé');
axis([0 t(end)+Te 0 Ka*max(abs(e))+0.1]); grid on; pause

% Affichage des résultats
text_mc_a  = [' méthode des MC :    a biais = ',num2str(biais_a_mc,4),...
                           '     écart type = ',num2str(std_a_mc,4)];
text_mc_b  = ['                     b biais = ',num2str(biais_b_mc,4),...
                           '     écart type = ',num2str(std_b_mc,4)];
text_oe_a = [' méthode du modèle : a biais = ',num2str(biais_a_oe,4),...
                          '     écart type = ',num2str(std_a_oe,4)];
text_oe_b = ['                     b biais = ',num2str(biais_b_oe,4),...
                          '     écart type = ',num2str(std_b_oe,4)];
disp(' '); disp(text_mc_a); disp(text_mc_b); disp(text_oe_a); disp(text_oe_b);

end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         % 
% Méthode des MCR : illustration du fonctionnement                        % 
%                                                                         % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if choix == 3;   
N = 100; % nombre d'échantillons

% Système
systeme_stationnaire;
% Signal d'entrée
entree;
% Signal de sortie du système
RSB = input(' rapport signal sur bruit (RSB) : '); if RSB == 0; RSB = 1e-16; end;
s = lsim(sys,e); s = s(:); Ps = sum(s.^2)/N; % sortie non bruitée
Pb = Ps/RSB; sigma_b = sqrt(Pb); s_b = s+sigma_b*randn(N,1); % sortie bruitée

% Initialisation
N1 = 1;
a_est(1) = -1; b_est(1) = 1; est = [a_est(1); b_est(1)]; 
P = diag([1^2,20^2]); P11(N1) = P(1,1); P22(N1) = P(2,2); 
% ou
N1 = 4; 
x = s_b(2:N1); C = [-s_b(1:N1-1) e(1:N1-1)];
est = inv(C.'*C)*C.'*x; a_est(N1) = est(1); b_est(N1) = est(2); 
mod_mc = est(2)*z^-1/(1+est(1)*z^-1);
P = inv(C'*C);
s_est = lsim(mod_mc,e); s_pred = s_est;
P = diag([1^2,20^2]); P11(N1) = P(1,1); P22(N1) = P(2,2); 

% Estimation récursive
for n = N1+1:N;
    c = [-s_b(n-1) e(n-1)]';
    K = P*c/(c'*P*c+1); K1(n) = K(1); K2(n) = K(2);
    P = (P-K*c.'*P); P11(n) = P(1,1); P22(n) = P(2,2); 
    est = est+K*(s_b(n)-c.'*est); a_est(n) = est(1); b_est(n) = est(2);
    s_pred(n) = b_est(n-1)*e(n-1)-a_est(n-1)*s_b(n-1);
    s_est(n)= b_est(n)*e(n-1)-a_est(n)*s_est(n-1);    
end;

figure;
subplot(211); stairs(t,e,'r'); hold on; plot(t,s,'b--',t,s_b,'b'); title('Signaux');
legend('entrée','sortie','sortie bruitée'); xlabel('temps (s)'); 
axis([0 N*Te -Ka*max(abs(e))-0.5 Ka*max(abs(e))+0.5]); grid on;
subplot(212); plot(t,a*ones(1,N),'r',t(N1:N),a_est(N1:N),'g',...
                   t,b*ones(1,N),'-.r',t(N1:N),b_est(N1:N),'-.g');
title('Paramètres estimés par MCR');
legend('a : exact','a : estimé','b : exact','b : estimé'); xlabel('temps (s)');
axis([0 N*Te min([a_est b_est])-0.1 max([a_est b_est])+0.1]); grid; pause;

figure;
subplot(221); stairs(t,e,'r'); hold on; plot(t,s,'b--',t,s_b,'b'); title('Signaux');
legend('entrée','sortie','sortie bruitée'); xlabel('temps (s)');
axis([0 N*Te -Ka*max(abs(e))-0.5 Ka*max(abs(e))+0.5]); grid on;
subplot(222); stairs(t,e,'r'); hold on; plot(t,s,'b--',t,s_b,'b'); title('Signaux');
legend('entrée','sortie','sortie bruitée');
xlabel('temps (s)'); axis([0 N*Te -Ka*max(abs(e))-0.5 Ka*max(abs(e))+0.5]); grid on;
subplot(223); plot(t(N1:N),P11(N1:N).^0.5,'g',t(N1:N),P22(N1:N).^0.5,'-.g'); 
title('Matrice P'); legend('a : sqrt(P(1,1))','b : sprt(P(2,2)'); 
xlabel('temps (s)'); axis([0 N*Te 0 2]); grid on;
subplot(224); plot(t(N1+1:N),K1(N1+1:N),'g',t(N1+1:N),K2(N1+1:N),'-.g'); 
title('Gain K'); legend('K(1)','K(2)');
xlabel('temps (s)'); axis([0 N*Te -0.5 0.5]); grid on;

end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         % 
% Méthode des MCR avec facteur d'oubli : illustration du fonctionnement   % 
% avec un système dont les paramètres varient (saut ou variaion linéaire) %
% au milieu des acquisitions
%                                                                         % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if choix == 4;   
N = 500; % nombre d'échantillons

% Système
systeme_non_stationnaire;
% Signal d'entrée
entree;
% Signal de sortie du système
RSB = input(' rapport signal sur bruit (RSB) au début : '); if RSB == 0; RSB = 1e-16; end;
s(1) = 0;
for n = 2:N; 
    s(n) = b(n)*e(n-1)-a(n)*s(n-1); end; s = s(:); Ps = sum(s(1:N/2-100).^2)/(N/2-100);
Pb = Ps/RSB; sigma_b = sqrt(Pb); s_b = s+sigma_b*randn(N,1); % sortie bruitée

% Initialisation
N1 = 4; 
x = s_b(2:N1); C = [-s_b(1:N1-1) e(1:N1-1)];
est = inv(C.'*C)*C.'*x; a_est(N1) = est(1); b_est(N1) = est(2); 
mod_mc = est(2)*z^-1/(1+est(1)*z^-1);
P = inv(C'*C);
s_est = lsim(mod_mc,e); s_pred = s_est;
P = diag([1^2,20^2]); P11(N1) = P(1,1); P22(N1) = P(2,2); 

% Estimation récursive
lambda = input(' facteur d''oubli : ');
for n = N1:N;
    c = [-s_b(n-1) e(n-1)]';
    K = P*c/(c'*P*c+lambda); K1(n) = K(1); K2(n) = K(2);
    P = 1/lambda*(P-K*c.'*P);     P11(n) = P(1,1); P22(n) = P(2,2); 
    est = est+K*(s_b(n)-c.'*est); a_est(n) = est(1); b_est(n) = est(2);
    s_pred(n) = b_est(n-1)*e(n-1)-a_est(n-1)*s_b(n-1);
    s_est(n)= b_est(n)*e(n-1)-a_est(n)*s_est(n-1);    
end;

figure;
subplot(211); stairs(t,e,'r'); hold on; plot(t,s,'b--',t,s_b,'b'); 
title('Signaux'); legend('entrée','sortie','sortie bruitée');
xlabel('temps (s)'); axis([0 N*Te -max(abs(e))-0.5 max(abs(e))+0.5]); grid on;
subplot(212); plot(t,a,'r',t(N1:N),a_est(N1:N),'g',t,b,'-.r',t(N1:N),b_est(N1:N),'-.g');               
title('Paramètres estimés par MCR');
legend('a : exact','a : estimé','b : exact','b : estimé');
xlabel('temps (s)'); axis([0 N*Te min([a_est b_est])-0.1 max([a_est b_est])+0.1]); grid; pause;

figure;
subplot(221); stairs(t,e,'r'); hold on; plot(t,s,'b--',t,s_b,'b'); title('Signaux');
legend('entrée','sortie','sortie bruitée'); xlabel('temps (s)'); 
axis([0 N*Te -max(abs(e))-0.5 max(abs(e))+0.5]); grid on;
subplot(222); stairs(t,e,'r'); hold on; plot(t,s,'b--',t,s_b,'b'); title('Signaux');
legend('entrée','sortie','sortie '); xlabel('temps (s)');
axis([0 N*Te -max(abs(e))-0.5 max(abs(e))+0.5]); grid on;
subplot(223); plot(t(N1:N),P11(N1:N).^0.5,'g',t(N1:N),P22(N1:N).^0.5,'-.g'); 
legend('a : sqrt(P(1,1))','b : sprt(P(2,2)'); 
title('Matrice P'); xlabel('temps (s)'); axis([0 N*Te 0 2]); grid on;
subplot(224); plot(t(N1+1:N),K1(N1+1:N),'g',t(N1+1:N),K2(N1+1:N),'-.g'); 
title('Gain K'); legend('K(1)','K(2)');
xlabel('temps (s)'); axis([0 N*Te -0.5 0.5]); grid on; pause;

% Réponse unitaire
e = ones(N,1); ru(N1) = 0; ru_est(N1) = 0;
for n = N1+1:N;
    ru(n) = b(n)*e(n-1)-a(n)*ru(n-1);
    ru_est(n) = b_est(n)*e(n-1)-a_est(n)*ru_est(n-1);
end;
ru = ru(:); ru_est = ru_est(:);
figure;
plot(t(N1:N),ru(N1:N),'b',t(N1:N),ru_est(N1:N),'g'); title('Comparaison des RU');
legend('ru','ru estimée mcr'); xlabel('temps (s)');  grid on;

end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         % 
% Méthode des MCR avec facteur d'oubli : étude statistique avec un        % 
% système stationnaire                                                    %
%                                                                         % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if choix == 5;   
N = 500; % nombre d'échantillons
I = 100; % nombre de simulations

% Système
systeme_stationnaire;
% Signal d'entrée
entree;
% Signal de sortie du système
RSB = input(' rapport signal sur bruit (RSB) : '); if RSB == 0; RSB = 1e-16; end;
s(1) = 0;
s = lsim(sys,e); s = s(:); Ps = sum(s.^2)/N; % sortie non bruitée
Ps = sum(s(1:(N/2-100)).^2)/(N/2-100);
Pb = Ps/RSB; sigma_b = sqrt(Pb); s_b = s+sigma_b*randn(N,1); % sortie bruitée

% Initialisation par MC
N1 = 4; 
x = s_b(2:N1); C = [-s_b(1:N1-1) e(1:N1-1)];
est = inv(C.'*C)*C.'*x; a_est(N1) = est(1); b_est(N1) = est(2); 
P = inv(C'*C);
% Estimation récursive avec facteur d'oubli
lambda = input(' choix de la valeur du facteur d''oubli : ');
for n = N1+1:N;
    c = [-s_b(n-1) e(n-1)]';
    K = P*c/(lambda+c'*P*c);
    P = 1/lambda*(P-K*c.'*P);
    est = est+K*(s_b(n)-c.'*est); a_est(n) = est(1); b_est(n) = est(2);
end;
figure;
subplot(121); plot(t,a*ones(1,N),'r',t,a_est,'g'); title('a');
xlabel('temps (s)'); axis([0 N*Te min(a)-0.1 max(a)+0.1]); grid on; hold on; 
subplot(122); plot(t,b*ones(1,N),'r',t,b_est,'g'); title('b');
xlabel('temps (s)'); axis([0 N*Te min(b)-0.1 max(b)+0.1]); grid on; hold on; 
for i = 2:I;
    s_b = s+sigma_b*randn(N,1);
    % Initialisation par MC
    x = s_b(2:N1); C = [-s_b(1:N1-1) e(1:N1-1)];
    est = inv(C.'*C)*C.'*x; a_est(N1) = est(1); b_est(N1) = est(2); 
    P = inv(C'*C);
    for n = N1+1:N;
        c = [-s_b(n-1) e(n-1)]';
        K = P*c/(lambda+c'*P*c);
        P = 1/lambda*(P-K*c.'*P);
        est = est+K*(s_b(n)-c.'*est); a_est(n) = est(1); b_est(n) = est(2);
    end;
    subplot(121); plot(t,a_est,'g');
    subplot(122); plot(t,b_est,'g');
end;
subplot(121); plot(t,a*ones(1,N),'r',t,a_est,'g'); legend('valeur exacte','valeur estimée');
subplot(122); plot(t,b*ones(1,N),'r',t,b_est,'g'); legend('valeur exacte','valeur estimée');

end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         % 
% Méthode des MCR avec facteur d'oubli : étude statistique avec un        % 
% système dont les paramètres varient                                     %
%                                                                         % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if choix == 6;   
N = 1000; % nombre d'échantillons
I = 100; % nombre de simulations

% Système
systeme_non_stationnaire;
% Signal d'entrée
entree;
% Signal de sortie du système
RSB = input(' rapport signal sur bruit (RSB) au début : '); if RSB == 0; RSB = 1e-16; end;
s(1) = 0;
for n = 2:N; 
    s(n) = b(n)*e(n-1)-a(n)*s(n-1); end; s = s(:); Ps = sum(s(1:(N/2-100)).^2)/(N/2-100);
Pb = Ps/RSB; sigma_b = sqrt(Pb); s_b = s+sigma_b*randn(N,1); % sortie bruitée

% Initialisation par MC
N1 = 4; 
x = s_b(2:N1); C = [-s_b(1:N1-1) e(1:N1-1)];
est = inv(C.'*C)*C.'*x; a_est(N1) = est(1); b_est(N1) = est(2); 
P = inv(C'*C);
% Estimation récursive avec facteur d'oubli
lambda = input(' choix de la valeur du facteur d''oubli : ');
for n = N1+1:N;
    c = [-s_b(n-1) e(n-1)]';
    K = P*c/(lambda+c'*P*c);
    P = 1/lambda*(P-K*c.'*P);
    est = est+K*(s_b(n)-c.'*est); a_est(n) = est(1); b_est(n) = est(2);
end;
figure;
subplot(121); plot(t,a,'r',t,a_est,'g'); title('a');
xlabel('temps (s)'); axis([0 N*Te min(a)-0.1 max(a)+0.1]); grid on; hold on; 
subplot(122); plot(t,b,'r',t,b_est,'g'); title('b');
xlabel('temps (s)'); axis([0 N*Te min(b)-0.1 max(b)+0.1]); grid on; hold on; 
for i = 2:I;
    s_b = s+sigma_b*randn(N,1);
    % Initialisation par MC
    x = s_b(2:N1); C = [-s_b(1:N1-1) e(1:N1-1)];
    est = inv(C.'*C)*C.'*x; a_est(N1) = est(1); b_est(N1) = est(2); 
    P = inv(C'*C);
    for n = N1+1:N;
        c = [-s_b(n-1) e(n-1)]';
        K = P*c/(lambda+c'*P*c);
        P = 1/lambda*(P-K*c.'*P);
        est = est+K*(s_b(n)-c.'*est); a_est(n) = est(1); b_est(n) = est(2);
    end;
    subplot(121); plot(t,a_est,'g');
    subplot(122); plot(t,b_est,'g');
end;
subplot(121); plot(t,a,'r',t,a_est,'g'); legend('valeur exacte','valeur estimée');
subplot(122); plot(t,b,'r',t,b_est,'g'); legend('valeur exacte','valeur estimée');

end;