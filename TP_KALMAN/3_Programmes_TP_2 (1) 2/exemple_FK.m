%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         % 
% Estimation de l'état d'un système par filtrage de Kalman : véhicule en  %
% déplacement rectiligne quasiuniforme avec mesure la distance parcourue  %
% et état = (position,vitesse).                                           %
%                                                                         %
% Les points étudiée sont :                                               %
%   - illustration du filtre de Kalman,                                   %
%   - étude statistique,                                                  %
%   - effet des erreurs de modèle.                                        % 
%                                                                         %  
% Cécile Durieu                                                           %
% 4 avril 2018, modifié 2 avril 2019 et 11 février 2025                   % 
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear all; clc;

choix = menu('Filtre de Kalman',...
       '1/ Illustration du fonctionnement du filtre de Kalman',...
       '2/ Etude statistique avec des bruits gaussiens',...
       '3/ Effet de la ddp : étude statistique avec des bruits uniformes',...
       '4/ Effet de l''état initial',...
       '5/ Effet de la covariance initiale',...
       '6/ Effet de l''erreur sur la covariance du bruit d''état',...
       '7/ Effet de l''erreur sur la covariance du bruit d''observation',...
       '8/ Test de cohérence');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         % 
% Modèle du système                                                       %
%                                                                         % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% Hypothèses : 
%    - accélération aléatoire et constante entre les instants de mesure
%    - bruit d'observation additif
%    - bruit d'état et bruit d'observation blancs et décorrélés entre eux 
%      et décorrélés de l'état initial
%
% Notations
% T : écart entre les instants de mesure
% g : accélération
% v : vitesse
% d : distance parcourue
% y : observation
% w : bruit d'observation
%
% Equations d'évolution / d'état :
% v(k+1) = v(k) + g(k)*T;
% d(k+1) = d(k) + v(k)*T + g(k)*T^2/2;
% état = (d;v);
% x(k+1) = A(k)*x(k) + b(k);
%        = A(k)*x(k) + B*g(k)
%        = | 1 T |*x(k) + | T^2/2 |*g(k)
%          | 0 1 |        |   T   |    
%    g(k) étant une variabla aléatoire centrée d'écart type std_g ainsi
%    b(k) est une variable aléatoire centrée de matrice de covariance
%    Q = std_g^2*B*B';
%
% Equations d'observation :
% y(k) = d(k) + w(k);       
%      = C(k)*x(k) + w(k);      
%      = |1 0| * x(k) + w(k);
%    w(k) étant une variabla aléatoire centrée d'écart type std_w ou encore
%    de covariance R = std_w^2

I = eye(2);
T = 1; % cadence des mesures
v_0 = 1; m_0 = [0;v_0]; P_0 = 100*I; % état initial

A = [1 T; 0 1]; % matrice d'état
V = [T^2/2; T]; % matrice du bruit d'état
std_g = 0.01; Q = std_g^2*V*V'; % matrice de covariance du bruit d'état
C = [1 0]; % matrice d'observationx
std_w = 0.1*v_0*T; R = std_w^2; % matrice de covariance du bruit d'observation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         % 
% Illustration du filtre de Kalman                                        %
%                                                                         % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if choix == 1;
    
% Simulation d'une trajectoire  
k_max = 20; % nombre d'itérations / d'observations
rng('default'); % initialisation du générateur de nombres aléatoires
g = std_g*randn(1,k_max); % bruit d'état
w = std_w*randn(1,k_max); % bruit d'observation
x(:,1) = m_0; y(1) = C*x(:,1)+w(1);
for k = 1:k_max-1; % simulation des observations
    x(:,k+1) = A*x(:,k)+V*g(k);
    y(k+1) = C*x(:,k+1)+w(k+1);
end
% estimation de la vitesse par dérivation de la position (dérivée arrière)
x2_d(2:k_max) = (y(2:k_max)-y(1:k_max-1))/T; 
figure; plot(x(1,:),x(2,:),'r');
        title('évolution de l''état'); legend('état exact');
        xlabel('position'); ylabel('vitesse');
        m1 = min(x(1,:)); M1 = max(x(1,:)); m2 = min(x(2,:)); M2 = max(x(2,:));
        axis([m1-1 M1+1 m2-0.3 M2+0.3]); pause;
figure;
subplot(211); plot(1:k_max,x(1,:),'r');
              title('évolution de la position'); legend('position exacte');
              xlabel('#'); ylabel('position'); axis([1 k_max m1-1 M1+1]);
subplot(212); plot(1:k_max,x(2,:),'r'); legend('vitesse exacte');
              title('évolution de la vitesse'); xlabel('#'); ylabel('vitesse');
              axis([1 k_max m2-0.3 M2+0.3]); pause;
figure;
subplot(211); plot(1:k_max,x(1,:),'r',1:k_max,y(1,:),'m');
              legend('position exacte','position observée'); 
              title('évolution de la position');
              xlabel('#'); ylabel('position'); axis([0 k_max m1-1 M1+1]);
subplot(212); plot(1:k_max,x(2,:),'r',2:k_max,x2_d(1,2:k_max),'m');
              legend('vitesse exacte','dérivée de la position observée');
              title('évolution de la vitesse');
              xlabel('#'); ylabel('vitesse'); axis([0 k_max m2-0.3 M2+0.3]); pause;

% Estimation de la trajectoire
x_est(:,1) = m_0; P = P_0;
for k = 1:k_max-1;
    x_pred(:,k+1) = A*x_est(:,k);
    P = A*P*A'+Q; 
    K = P*C'*inv(C*P*C'+R);
    P = (I-K*C)*P;
    y_pred(k+1) = C*x_pred(:,k+1); 
    x_est(:,k+1) = x_pred(:,k+1)+K*(y(k+1)-y_pred(k+1));
end
figure; plot(x(1,:),x(2,:),'r',x_est(1,:),x_est(2,:),'b');
        legend('état exact','état estimé par FK');
        title('évolution de l''état'); xlabel('position'); ylabel('vitesse');
        axis([m1-1 M1+1 m2-0.3 M2+0.3]); pause;
figure; 
subplot(211); plot(1:k_max,x(1,:),'r',1:k_max,x_est(1,:),'b');
              legend('position exacte','position estimée par FK');
              title('évolution de la position'); xlabel('#'); ylabel('position');
              axis([0 k_max m1-1 M1+1]);
subplot(212); plot(1:k_max,x(2,:),'r',1:k_max,x_est(2,:),'b');
              legend('vitesse exacte','vitesse estimée par FK');
              title('évolution de la vitesse'); xlabel('#'); ylabel('vitesse');
              axis([0 k_max m2-0.3 M2+0.3]); pause;
figure; 
subplot(211); plot(1:k_max,x(1,:),'r',1:k_max,x_est(1,:),'b',1:k_max,y(1,:),'m')
              legend('position exacte','position estimée par FK','position observée');
              title('évolution de la position'); xlabel('#'); ylabel('position');
              axis([0 k_max m1-1 M1+1]);
subplot(212); plot(1:k_max,x(2,:),'r',1:k_max,x_est(2,:),'b',2:k_max,x2_d(1,2:k_max),'m');
              legend('vitesse exacte','vitesse estimée par FK','dérivée de la position observée');
              title('évolution de la vitesse'); xlabel('#'); ylabel('vitesse');
               axis([0 k_max m2-0.3 M2+0.3]); pause;
figure; 
subplot(211); plot(1:k_max,x_est(1,:)-x(1,:),'b',1:k_max,y-x(1,:),'m')
              legend('estimation par FK','position observée');
              title('évolution de l''erreur de position'); xlabel('#'); ylabel('position');
              axis([0 k_max -3*std_w 3*std_w]);
subplot(212); plot(1:k_max,x_est(2,:)-x(2,:),'b',2:k_max,x2_d(2:k_max)-x(2,2:k_max),'m');
              legend('estimation par FK','dérivée de la position observée');
              title('évolution de l''erreur de vitesse'); xlabel('#'); ylabel('vitesse');
              axis([0 k_max -3*sqrt(2)*std_w/T 3*sqrt(2)*std_w/T]); pause;

% Estimation d'une trajectoire avec sauvegarde de grandeurs et affichage
x_est(:,1) = m_0; P = P_0; P1_est(:,1) = P(:,1); P2_est(:,1) = P(:,2);
figure; hold on
for k = 1:k_max-1;
    rectangle('Position',[y(k+1)-3*R^0.5,-10,6*R^0.5,20],'LineStyle','none','FaceColor',[1 1 0.8]);
end;
trace_ellipse_P(x_est(:,1),9*P,'b');
for k = 1:k_max-1;
    x_pred(:,k+1) = A*x_est(:,k);
    e_pred(:,k+1) = x_pred(:,k+1)-x(:,k+1);
    P = A*P*A';
    trace_ellipse_P(x_pred(:,k+1),9*P,':g');
    P = P+Q;  P1_pred(:,k+1) = P(:,1); P2_pred(:,k+1) = P(:,2); 
    trace_ellipse_P(x_pred(:,k+1),9*P,'g');
    K = P*C'*inv(C*P*C'+R);
    P = (I-K*C)*P; P1_est(:,k+1) = P(:,1); P2_est(:,k+1) = P(:,2); 
    y_pred(k+1) = C*x_pred(:,k+1); 
    x_est(:,k+1) = x_pred(:,k+1)+K*(y(k+1)-y_pred(k+1));
    e_est(:,k+1) = x_est(:,k+1)-x(:,k+1);  
    trace_ellipse_P(x_est(:,k+1),9*P,'b');
end
title('évolution de l''état et précision'); xlabel('position'); ylabel('vitesse');
legend('estimation','prédiction sans bruit d''état','prédiction','correction');
axis([m1-1 M1+1 m2-0.1 M2+0.1]); pause

kz = input('préciser le numéro de l''observation pour le zoom : ');
x_est(:,1) = m_0; P = P_0;
for k = 1:kz-1;
    x_pred(:,k+1) = A*x_est(:,k);
    P = A*P*A'+Q;
    K = P*C'*inv(C*P*C'+R); P = (I-K*C)*P;
    y_pred(k+1) = C*x_pred(:,k+1); 
    x_est(:,k+1) = x_pred(:,k+1)+K*(y(k+1)-y_pred(k+1));    
end;
figure;
plot(x_est(1,k+1),x_est(2,k+1),'*b'); hold on; 
trace_ellipse_P(x_est(:,k+1),9*P,'b'); 
legend('état estimé à (k-1)','précision à (k-1)');
title('évolution de l''état et précision'); xlabel('position'); ylabel('vitesse');
axis([x_est(1,k+1)-1 x_est(1,k+1)+2 x_est(2,k)-0.3 x_est(2,k)+0.3]); pause;
k = k+1; x_pred(:,k+1) = A*x_est(:,k); P = A*P*A';
plot(x_pred(1,k+1),x_pred(2,k+1),'*g'); trace_ellipse_P(x_pred(:,k+1),9*P,':g'); 
legend('état estimé à (k-1)','précision estimation à (k-1)',...
       'état prédit à k','précision prédiction sans bruit d''état à k'); pause;
P = P+Q; trace_ellipse_P(x_pred(:,k+1),9*P,'g'); 
legend('état estimé à (k-1)','précision estimation à (k-1)',...
       'état prédit à k','précision prédiction sans bruit d''état à k',...
       'précision prédiction à k'); pause;
uistack(rectangle('Position',[y(k+1)-3*R^0.5,-10,6*R^0.5,20],...
                  'LineStyle','none','FaceColor',[1 1 0.8]),'bottom'); pause;
K = P*C'*inv(C*P*C'+R); P = (I-K*C)*P;
y_pred(k+1) = C*x_pred(:,k+1); x_est(:,k+1) = x_pred(:,k+1)+K*(y(k+1)-y_pred(k+1));
plot(x_est(1,k+1),x_est(2,k+1),'*b'); trace_ellipse_P(x_est(:,k+1),9*P,'b'); 
legend('état estimé à (k-1)','précision estimation à (k-1)',...
       'état prédit à k','précision prédiction sans bruit d''état à k',...
       'précision prédiction à k','état estimé à k','précision estimation à k'); pause;
   
figure;
subplot(211); plot(1:k_max,P1_est(1,:).^0.5,'b');
              legend('std position calculée par FK');
              title('précision de l''estimation de la position');
              xlabel('#'); ylabel('position');
              axis([0 k_max 0 3*std_w]);
subplot(212); plot(1:k_max,P2_est(2,:).^0.5,'b');
              legend('std vitesse calculée par FK');
              title('précision de l''estimation de la vitesse');
              xlabel('#'); ylabel('vitesse');
              axis([0 k_max 0 3*sqrt(2)*std_w/T]); pause;
figure;
subplot(211); plot(1:k_max,P1_est(1,:).^0.5,'b',1:k_max,std_w*ones(1,k_max),'m');
              legend('std position calculée par FK','std observation de la position'); 
              title('précision de l''estimation de la position'); xlabel('#'); ylabel('position');
              axis([0 k_max 0 3*std_w]);
subplot(212); plot(1:k_max,P2_est(2,:).^0.5,'b',1:k_max,std_w*ones(1,k_max),'m');
              legend('std vitesse calculée par FK','std dérivée de la position observée');
              title('précision de l''estimation de la vitesse'); xlabel('#'); ylabel('vitesse');
              axis([0 k_max 0 3*sqrt(2)*std_w/T]); pause;
figure;
subplot(211); plot(1:k_max,x(1,:),'r',1:k_max,x_est(1,:),'b',...
                   1:k_max,x_est(1,:)-3*P1_est(1,:).^0.5,':b',...
                   1:k_max,x_est(1,:)+3*P1_est(1,:).^0.5,':b',...
                   1:k_max,y,'g');
              legend('position exacte','position estimée par FK',....
                     'position estimée -3*std','position estimée +3*std',...
                     'position observée')
              title('estimation avec précision associée');  
              xlabel('#'); ylabel('position');  axis([0 k_max m1-1 M1+1]);
subplot(212); plot(1:k_max,x(2,:),'r',1:k_max,x_est(2,:),'b',...
                   1:k_max,x_est(2,:)-3*P2_est(2,:).^0.5,':b',...
                   1:k_max,x_est(2,:)+3*P2_est(2,:).^0.5,':b',....
                   2:k_max,x2_d(2:k_max),'g');
              legend('vitesse exacte','vitesse estimée par FK',...
                     'vitesse estimée -3*std','vitesse estimée +3*std',...
                     'dérivée de la position observée')
              title('estimation avec précision associée'); xlabel('#'); ylabel('vitesse');
              axis([0 k_max m2-1 M2+1]); pause;
figure;
subplot(211); plot(1:k_max,x_est(1,:)-x(1,:),'b',...
                   1:k_max,-3*P1_est(1,:).^0.5,':b',1:k_max,3*P1_est(1,:).^0.5,':b');
              legend('erreur FK','-3*std','+3*std'); title('erreur d''estimation de la position'); 
              xlabel('#'); ylabel('position'); axis([0 k_max -5*std_w 5*std_w]);
subplot(212); plot(1:k_max,x_est(2,:)-x(2,:),'b',...
                   1:k_max,-3*P2_est(2,:).^0.5,':b',1:k_max,+3*P2_est(2,:).^0.5,':b');
              legend('erreur FK','-3*std','+3*std'); title('erreur d''estimation de la vitesse');
              xlabel('#'); ylabel('vitesse'); axis([0 k_max -5*sqrt(2)*std_w/T 5*sqrt(2)*std_w/T]);
              
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         % 
% Etude statistique avec des bruits gaussiens                             %
%                                                                         % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if choix == 2;
 
k_max = 20; % nombe d'itérations / d'observations
N = 10000; % nombre de trajectoires
nuage_n(A,V,C,std_g,std_w,m_0,P_0,Q,R,k_max,N,1);
title('erreur d''estimation et ellipsoide de confiance avec des bruits gaussiens');

end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         % 
%  Effet de la ddp : étude statistique avec des bruits uniformes          %
%                                                                         % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if choix == 3;
  
k_max = 20; % nombe d'itérations / d'observations
N = 10000; % nombre de trajectoire
nuage_u(A,V,C,std_g,std_w,m_0,P_0,Q,R,k_max,N,1);
title('erreur d''estimation et ellipsoide de confiance avec des bruits uniformes');

end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         % 
% Effet de l'initialisation : état                                        %
%                                                                         % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if choix == 4;

k_max = 20; % nombe d'itérations / d'observations
rng('default');
g = std_g*randn(1,k_max); % bruit d'état
w = std_w*randn(1,k_max); % bruit d'observation
x(:,1) = m_0; y(1) = C*x(:,1)+w(1);
for k = 1:k_max-1;
    x(:,k+1) = A*x(:,k)+V*g(k);
    y(k+1) = C*x(:,k+1)+w(k+1);
end

% 1re initialisation
x_est(:,1) = m_0; P = P_0;
for k = 1:k_max-1;
    x_pred(:,k+1) = A*x_est(:,k);
    P = A*P*A'+Q; 
    K = P*C'*inv(C*P*C'+R);
    P = (I-K*C)*P;
    y_pred(k+1) = C*x_pred(:,k+1); 
    x_est(:,k+1) = x_pred(:,k+1)+K*(y(k+1)-y_pred(k+1));
end
m1 = min(x(1,:)); M1 = max(x(1,:)); m2 = min(x(2,:)); M2 = max(x(2,:));
figure; subplot(211); plot(1:k_max,x(1,:),'r',1:k_max,x_est(1,:),'b');
        subplot(212); plot(1:k_max,x(2,:),'r',1:k_max,x_est(2,:),'b');

% 2e initialisation
x_est(:,1) = [0; 0]; P = P_0;
for k = 1:k_max-1;
    x_pred(:,k+1) = A*x_est(:,k);
    P = A*P*A'+Q; 
    K = P*C'*inv(C*P*C'+R);
    P = (I-K*C)*P;
    y_pred(k+1) = C*x_pred(:,k+1); 
    x_est(:,k+1) = x_pred(:,k+1)+K*(y(k+1)-y_pred(k+1));
end
subplot(211); hold on; plot(1:k_max,x_est(1,:),'m');
              legend('position exacte','position estimée par FK 1','position estimée par FK 2');
              title('évolution de la position'); xlabel('#'); ylabel('position');
              axis([0 k_max m1-1 M1+1]);
subplot(212); hold on; plot(1:k_max,x_est(2,:),'m');
              legend('vitesse exacte','vitesse estimée par FK 1','vitesse estimée par FK 2');
              title('évolution de la vitesse'); xlabel('#'); ylabel('vitesse');
              axis([0 k_max m2-1 M2+1]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         % 
% Effet de l'initialisation : covariance                                  %
%                                                                         % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if choix == 5;

k_max = 50; % nombe d'itérations / d'observations
rng('default');
g = std_g*randn(1,k_max); % bruit d'état
w = std_w*randn(1,k_max); % bruit d'observation
x(:,1) = m_0; y(1) = C*x(:,1)+w(1);
for k = 1:k_max-1;
    x(:,k+1) = A*x(:,k)+V*g(k);
    y(k+1) = C*x(:,k+1)+w(k+1);
end

% 1re initialisation
x_est(:,1) = [0; 0]; P = P_0; P1_est(:,1) = P(:,1); P2_est(:,1) = P(:,2);
for k = 1:k_max-1;
    x_pred(:,k+1) = A*x_est(:,k);
    P = A*P*A'+Q; 
    K = P*C'*inv(C*P*C'+R);
    P = (I-K*C)*P; P1_est(:,k+1) = P(:,1); P2_est(:,k+1) = P(:,2);
    y_pred(k+1) = C*x_pred(:,k+1); 
    x_est(:,k+1) = x_pred(:,k+1)+K*(y(k+1)-y_pred(k+1));
end
figure(1); subplot(211); plot(1:k_max,x(1,:),'r',1:k_max,x_est(1,:),'-+b');
           subplot(212); plot(1:k_max,x(2,:),'r',1:k_max,x_est(2,:),'-+b');
figure(2); subplot(211); plot(1:k_max,P1_est(1,:).^0.5,'-+b');
           subplot(212); plot(1:k_max,P2_est(2,:).^0.5,'-+b');

% 2e initialisation
x_est(:,1) = [0; 0]; P = 100*P_0; P1_est(:,1) = P(:,1); P2_est(:,1) = P(:,2);
for k = 1:k_max-1;
    x_pred(:,k+1) = A*x_est(:,k);
    P = A*P*A'+Q; 
    K = P*C'*inv(C*P*C'+R);
    P = (I-K*C)*P; P1_est(:,k+1) = P(:,1); P2_est(:,k+1) = P(:,2);
    y_pred(k+1) = C*x_pred(:,k+1); 
    x_est(:,k+1) = x_pred(:,k+1)+K*(y(k+1)-y_pred(k+1));
end
figure(1); subplot(211); hold on; plot(1:k_max,x_est(1,:),'-om');
           subplot(212); hold on; plot(1:k_max,x_est(2,:),'-om');
figure(2); subplot(211); hold on; plot(1:k_max,P1_est(1,:).^0.5,'-om');
           subplot(212); hold on; plot(1:k_max,P2_est(2,:).^0.5,'-om');

% 3ème initialisation
x_est(:,1) = [0; 0]; P = 0.01*P_0; P1_est(:,1) = P(:,1); P2_est(:,1) = P(:,2);
for k = 1:k_max-1;
    x_pred(:,k+1) = A*x_est(:,k);
    P = A*P*A'+Q; 
    K = P*C'*inv(C*P*C'+R);
    P = (I-K*C)*P; P1_est(:,k+1) = P(:,1); P2_est(:,k+1) = P(:,2);
    y_pred(k+1) = C*x_pred(:,k+1); 
    x_est(:,k+1) = x_pred(:,k+1)+K*(y(k+1)-y_pred(k+1));
end
m1 = min(x(1,:)); M1 = max(x(1,:)); m2 = min(x(2,:)); M2 = max(x(2,:));
figure(1)
subplot(211); hold on; plot(1:k_max,x_est(1,:),'g');
              legend('position exacte','position estimée par FK : P_0',...
                     'position estimée par FK : 100*P_0','position estimée par FK : 0.01*P_0');
              title('évolution de la position'); xlabel('#'); ylabel('position');
              axis([0 k_max m1-1 M1+1]);
subplot(212); hold on; plot(1:k_max,x_est(2,:),'g');
              legend('vitesse exacte','vitesse estimée par FK : P_0',...
                      'vitesse estimée par FK : 100*P_0','vitesse estimée par FK : 0.01*P_0');
              title('évolution de la vitesse'); xlabel('#'); ylabel('vitesse');
              axis([0 k_max m2-1 M2+1]); pause
figure(2);
subplot(211); plot(1:k_max,P1_est(1,:).^0.5,'g');
              legend('std position calculée par FK : P_0',...
                     'std position calculée par FK : 100*P_0',...
                     'std position calculée par FK : 0.01*P_0');
              title('précision de l''estimation de la position');
              xlabel('#'); ylabel('position'); axis( [0 k_max 0 5*std_w]);
subplot(212);plot(1:k_max,P2_est(2,:).^0.5,'g');
             legend('std vitesse calculée par FK : P_0',...
                    'std vitesse calculée par FK : 100*P_0',...
                    'std vitesse calculée par FK : 0.01*P_0');
             title('précision de l''estimation de la vitesse')
             xlabel('#'); ylabel('vitesse'); axis([0 k_max 0 5*sqrt(2)*std_w/T]);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         % 
% Effet de la covriance du bruit d'état                                   %
%                                                                         % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if choix == 6;

k_max = 50; % nombe d'itérations / d'observations
rng('default');
g = std_g*randn(1,k_max); % bruit d'état
w = std_w*randn(1,k_max); % bruit d'observation
x(:,1) = m_0; y(1) = C*x(:,1)+w(1);
for k = 1:k_max-1;
    x(:,k+1) = A*x(:,k)+V*g(k);
    y(k+1) = C*x(:,k+1)+w(k+1);
end

% 1re initialisation
x_est(:,1) = [0; 0]; P = P_0; P1_est(:,1) = P(:,1); P2_est(:,1) = P(:,2);
for k = 1:k_max-1;
    x_pred(:,k+1) = A*x_est(:,k);
    P = A*P*A'+Q; 
    K = P*C'*inv(C*P*C'+R);
    P = (I-K*C)*P; P1_est(:,k+1) = P(:,1); P2_est(:,k+1) = P(:,2);
    y_pred(k+1) = C*x_pred(:,k+1); 
    x_est(:,k+1) = x_pred(:,k+1)+K*(y(k+1)-y_pred(k+1));
end
figure(1); subplot(211); plot(1:k_max,x(1,:),'r',1:k_max,x_est(1,:),'b');
           subplot(212); plot(1:k_max,x(2,:),'r',1:k_max,x_est(2,:),'b');
figure(2); subplot(211); plot(1:k_max,P1_est(1,:).^0.5,'b');
           subplot(212); plot(1:k_max,P2_est(2,:).^0.5,'b');

% 2e initialisation
x_est(:,1) = [0; 0]; P = P_0; P1_est(:,1) = P(:,1); P2_est(:,1) = P(:,2);
for k = 1:k_max-1;
    x_pred(:,k+1) = A*x_est(:,k);
    P = A*P*A'+10*Q; 
    K = P*C'*inv(C*P*C'+R);
    P = (I-K*C)*P; P1_est(:,k+1) = P(:,1); P2_est(:,k+1) = P(:,2);
    y_pred(k+1) = C*x_pred(:,k+1); 
    x_est(:,k+1) = x_pred(:,k+1)+K*(y(k+1)-y_pred(k+1));
end
figure(1); subplot(211); hold on; plot(1:k_max,x_est(1,:),'m');
           subplot(212); hold on; plot(1:k_max,x_est(2,:),'m');
figure(2); subplot(211); hold on; plot(1:k_max,P1_est(1,:).^0.5,'m');
           subplot(212); hold on; plot(1:k_max,P2_est(2,:).^0.5,'m');

% 3e initialisation
x_est(:,1) = [0; 0]; P = P_0; P1_est(:,1) = P(:,1); P2_est(:,1) = P(:,2);
for k = 1:k_max-1;
    x_pred(:,k+1) = A*x_est(:,k);
    P = A*P*A'+0.1*Q; 
    K = P*C'*inv(C*P*C'+R);
    P = (I-K*C)*P; P1_est(:,k+1) = P(:,1); P2_est(:,k+1) = P(:,2);
    y_pred(k+1) = C*x_pred(:,k+1); 
    x_est(:,k+1) = x_pred(:,k+1)+K*(y(k+1)-y_pred(k+1));
end
m1 = min(x(1,:)); M1 = max(x(1,:)); m2 = min(x(2,:)); M2 = max(x(2,:));
figure(1);
subplot(211); hold on; plot(1:k_max,x_est(1,:),'g');
              legend('position exacte','position estimée par FK : Q',...
                     'position estimée par FK : 10*Q','position estimée par FK : 0.1*Q');
              title('évolution de la position'); xlabel('#'); ylabel('position'); 
              axis([0 k_max m1-1 M1+1]);
subplot(212); hold on; plot(1:k_max,x_est(2,:),'g');
              legend('vitesse exacte','vitesse estimée par FK : Q',...
                     'vitesse estimée par FK : 10*Q','vitesse estimée par FK : 0.1*Q');
              title('évolution de la vitesse'); xlabel('#'); ylabel('vitesse');
              axis([0 k_max m2-1 M2+1]);
figure(2);
subplot(211); plot(1:k_max,P1_est(1,:).^0.5,'g');
              legend('std position calculée par FK : Q',...
                     'std position calculée par FK : 10*Q',...
                     'std position calculée par FK : 0.1*Q');
              title('précision de l''estimation de la position'); xlabel('#'); ylabel('position');
               axis( [0 k_max 0 5*std_w]);
subplot(212); plot(1:k_max,P2_est(2,:).^0.5,'g');
              legend('std vitesse calculée par FK : Q',...
                     'std vitesse calculée par FK : 10*Q',...
                     'std vitesse calculée par FK : 0.1*Q');
              title('précision de l''estimation de la vitesse'); xlabel('#'); ylabel('vitesse');
              axis([0 k_max 0 5*sqrt(2)*std_w/T]); pause    
              
N = 10000;              
nuage_n(A,V,C,std_g,std_w,m_0,P_0,Q,R,k_max,N,3);
title('erreur d''estimation et ellipsoide de confiance avec Q');
axis([-0.25 0.25 -0.25 0.25]); axis('square');
nuage_n(A,V,C,std_g,std_w,m_0,P_0,10*Q,R,k_max,N,4);
title('erreur d''estimation et ellipsoide de confiance avec 10*Q');
axis([-0.25 0.25 -0.25 0.25]); axis('square');
nuage_n(A,V,C,std_g,std_w,m_0,P_0,0.1*Q,R,k_max,N,5);
title('erreur d''estimation et ellipsoide de confiance avec 0.1*Q');
axis([-0.25 0.25 -0.25 0.25]); axis('square');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         % 
% Effet de la covariance du bruit d'observation                           %
%                                                                         % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if choix == 7;

k_max = 50; % nombe d'itérations / d'observations
rng('default');
g = std_g*randn(1,k_max); % bruit d'état
w = std_w*randn(1,k_max); % bruit d'observation
x(:,1) = m_0; y(1) = C*x(:,1)+w(1);
for k = 1:k_max-1;
    x(:,k+1) = A*x(:,k)+V*g(k);
    y(k+1) = C*x(:,k+1)+w(k+1);
end

% 1re simulation : R
x_est(:,1) = [0; 0]; P = P_0; P1_est(:,1) = P(:,1); P2_est(:,1) = P(:,2);
for k = 1:k_max-1;
    x_pred(:,k+1) = A*x_est(:,k);
    P = A*P*A'+Q; 
    K = P*C'*inv(C*P*C'+R);
    P = (I-K*C)*P; P1_est(:,k+1) = P(:,1); P2_est(:,k+1) = P(:,2);
    y_pred(k+1) = C*x_pred(:,k+1);
    x_est(:,k+1) = x_pred(:,k+1)+K*(y(k+1)-y_pred(k+1));
end
figure(1); subplot(211); plot(1:k_max,x(1,:),'r',1:k_max,x_est(1,:),'b');
           subplot(212); plot(1:k_max,x(2,:),'r',1:k_max,x_est(2,:),'b');
figure(2); subplot(211); plot(1:k_max,P1_est(1,:).^0.5,'b');
           subplot(212); plot(1:k_max,P2_est(2,:).^0.5,'b');

% 2e simulation  10*R
x_est(:,1) = [0; 0]; P = P_0; P1_est(:,1) = P(:,1); P2_est(:,1) = P(:,2);
for k = 1:k_max-1;
    x_pred(:,k+1) = A*x_est(:,k);
    P = A*P*A'+Q; 
    K = P*C'*inv(C*P*C'+10*R);
    P = (I-K*C)*P; P1_est(:,k+1) = P(:,1); P2_est(:,k+1) = P(:,2);
    y_pred(k+1) = C*x_pred(:,k+1); 
    x_est(:,k+1) = x_pred(:,k+1)+K*(y(k+1)-y_pred(k+1));
end
figure(1); subplot(211); hold on; plot(1:k_max,x_est(1,:),'m');
           subplot(212); hold on; plot(1:k_max,x_est(2,:),'m');
figure(2); subplot(211); hold on; plot(1:k_max,P1_est(1,:).^0.5,'m');
           subplot(212); hold on; plot(1:k_max,P2_est(2,:).^0.5,'m');

% 3e simulation : 0.1*R
x_est(:,1) = [0; 0]; P = P_0; P1_est(:,1) = P(:,1); P2_est(:,1) = P(:,2);
for k = 1:k_max-1;
    x_pred(:,k+1) = A*x_est(:,k);
    P = A*P*A'+Q; 
    K = P*C'*inv(C*P*C'+0.1*R);
    P = (I-K*C)*P; P1_est(:,k+1) = P(:,1); P2_est(:,k+1) = P(:,2);
    y_pred(k+1) = C*x_pred(:,k+1); 
    x_est(:,k+1) = x_pred(:,k+1)+K*(y(k+1)-y_pred(k+1));
end
m1 = min(x(1,:)); M1 = max(x(1,:)); m2 = min(x(2,:)); M2 = max(x(2,:));
figure(1); 
subplot(211); hold on; plot(1:k_max,x_est(1,:),'g');
              legend('position exacte','position estimée par FK : R',...
                     'position estimée par FK : 10*R','position estimée par FK : 0.1*R');
              title('évolution de la position'); xlabel('#'); ylabel('position'); 
              axis([0 k_max m1-1 M1+1]);
subplot(212); hold on; plot(1:k_max,x_est(2,:),'g');
              legend('vitesse exacte','vitesse estimée par FK : R',...
                     'vitesse estimée par FK : 10*R','vitesse estimée par FK : 0.1*R');
              title('évolution de la vitesse'); xlabel('#'); ylabel('vitesse');
              axis([0 k_max m2-1 M2+1]);
figure(2); 
subplot(211); plot(1:k_max,P1_est(1,:).^0.5,'g');
              legend('std position calculée par FK : R',...
                     'std position calculée par FK : 10*R',...
                     'std position calculée par FK : 0.1*R');
              title('précision de l''estimation de la position');
              xlabel('#'); ylabel('position'); axis( [0 k_max 0 5*std_w]);
subplot(212); plot(1:k_max,P2_est(2,:).^0.5,'g');
              legend('std vitesse calculée par FK : R',...
                     'std vitesse calculée par FK : 10*R',...
                     'std vitesse calculée par FK : 0.1*R');
              title('précision de l''estimation de la vitesse');
              xlabel('#'); ylabel('vitesse'); axis([0 k_max 0 5*sqrt(2)*std_w/T]); pause
              
N = 10000;              
nuage_n(A,V,C,std_g,std_w,m_0,P_0,Q,R,k_max,N,3);
title('erreur d''estimation et ellipsoide de confiance avec R');
axis([-0.25 0.25 -0.25 0.25]); axis('square');
nuage_n(A,V,C,std_g,std_w,m_0,P_0,Q,10*R,k_max,N,4);
title('erreur d''estimation et ellipsoide de confiance avec 10*R');
axis([-0.25 0.25 -0.25 0.25]); axis('square');
nuage_n(A,V,C,std_g,std_w,m_0,P_0,0.1*Q,0.1*R,k_max,N,5);
title('erreur d''estimation et ellipsoide de confiance avec 0.1*R');
axis([-0.25 0.25 -0.25 0.25]); axis('square');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         % 
% Test de cohérence                                                       %
%                                                                         % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if choix == 8;

k_max = 100; % nombe d'itérations / d'observations
rng('default');
g = std_g*randn(1,k_max); % bruit d'état
w = std_w*randn(1,k_max); % bruit d'observation
x(:,1) = m_0; y(1) = C*x(:,1)+w(1);
for k = 1:k_max-1;
    x(:,k+1) = A*x(:,k)+V*g(k);
    y(k+1) = C*x(:,k+1)+w(k+1);
end

% Etude de l'effet de Q
% 1re simulation : Q
x_est(:,1) = [0; 0]; P = P_0; P1_est(:,1) = P(:,1); P2_est(:,1) = P(:,2);
for k = 1:k_max-1;
    x_pred(:,k+1) = A*x_est(:,k);
    P = A*P*A'+Q; 
    M = inv(C*P*C'+R);
    K = P*C'*M;
    P = (I-K*C)*P; P1_est(:,k+1) = P(:,1); P2_est(:,k+1) = P(:,2);
    y_pred(k+1) = C*x_pred(:,k+1);
    delta(k+1) = y(k+1)-y_pred(k+1);
    dist(k+1) = sqrt(delta(k+1)*M*delta(k+1));
    x_est(:,k+1) = x_pred(:,k+1)+K*delta(k+1);
end
figure(1); hold on;
plot(2:k_max,dist(2:k_max),'b');

% 2e simulation : 10*Q
x_est(:,1) = [0; 0]; P = P_0; P1_est(:,1) = P(:,1); P2_est(:,1) = P(:,2);
for k = 1:k_max-1;
 x_pred(:,k+1) = A*x_est(:,k);
    P = A*P*A'+10*Q; 
    M = inv(C*P*C'+R);
    K = P*C'*M;
    P = (I-K*C)*P; P1_est(:,k+1) = P(:,1); P2_est(:,k+1) = P(:,2);
    y_pred(k+1) = C*x_pred(:,k+1);
    delta(k+1) = y(k+1)-y_pred(k+1);
    dist(k+1) = sqrt(delta(k+1)*M*delta(k+1));
    x_est(:,k+1) = x_pred(:,k+1)+K*delta(k+1);
end
plot(2:k_max,dist(2:k_max),'m');

% 3e simulation : 0.1*Q
x_est(:,1) = [0; 0]; P = P_0; P1_est(:,1) = P(:,1); P2_est(:,1) = P(:,2);
for k = 1:k_max-1;
    x_pred(:,k+1) = A*x_est(:,k);
    P = A*P*A'+0.1*Q; 
    M = inv(C*P*C'+R);
    K = P*C'*M;
    P = (I-K*C)*P; P1_est(:,k+1) = P(:,1); P2_est(:,k+1) = P(:,2);
    y_pred(k+1) = C*x_pred(:,k+1);
    delta(k+1) = y(k+1)-y_pred(k+1);
    dist(k+1) = sqrt(delta(k+1)*M*delta(k+1));
    x_est(:,k+1) = x_pred(:,k+1)+K*delta(k+1);
end
m1 = min(x(1,:)); M1 = max(x(1,:)); m2 = min(x(2,:)); M2 = max(x(2,:));
plot(2:k_max,dist(2:k_max),'g');
legend('distance : Q','distance : 10*Q','distance : 0.1*Q');
title('évolution de la distance Mahalanobis'); xlabel('#'); ylabel('distance'); 
axis([0 k_max 0 5]);

% Etude de l'effet de R
% 1re simulation : R
x_est(:,1) = [0; 0]; P = P_0; P1_est(:,1) = P(:,1); P2_est(:,1) = P(:,2);
for k = 1:k_max-1;
    x_pred(:,k+1) = A*x_est(:,k);
    P = A*P*A'+Q; 
    M = inv(C*P*C'+R);
    K = P*C'*M;
    P = (I-K*C)*P; P1_est(:,k+1) = P(:,1); P2_est(:,k+1) = P(:,2);
    y_pred(k+1) = C*x_pred(:,k+1);
    delta(k+1) = y(k+1)-y_pred(k+1);
    dist(k+1) = sqrt(delta(k+1)*M*delta(k+1));
    x_est(:,k+1) = x_pred(:,k+1)+K*delta(k+1);
end
figure(2); hold on;
plot(2:k_max,dist(2:k_max),'b');

% 2e simulation : 10*R
x_est(:,1) = [0; 0]; P = P_0; P1_est(:,1) = P(:,1); P2_est(:,1) = P(:,2);
for k = 1:k_max-1;
 x_pred(:,k+1) = A*x_est(:,k);
    P = A*P*A'+Q; 
    M = inv(C*P*C'+10*R);
    K = P*C'*M;
    P = (I-K*C)*P; P1_est(:,k+1) = P(:,1); P2_est(:,k+1) = P(:,2);
    y_pred(k+1) = C*x_pred(:,k+1);
    delta(k+1) = y(k+1)-y_pred(k+1);
    dist(k+1) = sqrt(delta(k+1)*M*delta(k+1));
    x_est(:,k+1) = x_pred(:,k+1)+K*delta(k+1);
end
plot(2:k_max,dist(2:k_max),'m');

% 3e simulation : 0.1*R
x_est(:,1) = [0; 0]; P = P_0; P1_est(:,1) = P(:,1); P2_est(:,1) = P(:,2);
for k = 1:k_max-1;
    x_pred(:,k+1) = A*x_est(:,k);
    P = A*P*A'+Q; 
    M = inv(C*P*C'+0.1*R);
    K = P*C'*M;
    P = (I-K*C)*P; P1_est(:,k+1) = P(:,1); P2_est(:,k+1) = P(:,2);
    y_pred(k+1) = C*x_pred(:,k+1);
    delta(k+1) = y(k+1)-y_pred(k+1);
    dist(k+1) = sqrt(delta(k+1)*M*delta(k+1));
    x_est(:,k+1) = x_pred(:,k+1)+K*delta(k+1);
end
m1 = min(x(1,:)); M1 = max(x(1,:)); m2 = min(x(2,:)); M2 = max(x(2,:));
plot(2:k_max,dist(2:k_max),'g');
legend('distance : R','distance : 10*R','distance : 0.1*R');
title('évolution de la distance de Mahalanobis'); xlabel('#'); ylabel('distance'); 
axis([0 k_max 0 10]);

end