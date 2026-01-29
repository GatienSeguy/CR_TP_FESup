% Script pour executer l'option 1 du filtre de Kalman
% et sauvegarder les figures

close all; clear all; clc;

% Dossier de sauvegarde des figures
output_dir = '../LATEX2/figures';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modele du systeme
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

I = eye(2);
T = 1; % cadence des mesures
v_0 = 1; m_0 = [0;v_0]; P_0 = 100*I; % etat initial

A = [1 T; 0 1]; % matrice d'etat
V = [T^2/2; T]; % matrice du bruit d'etat
std_g = 0.01; Q = std_g^2*V*V'; % matrice de covariance du bruit d'etat
C = [1 0]; % matrice d'observation
std_w = 0.1*v_0*T; R = std_w^2; % matrice de covariance du bruit d'observation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Option 1 : Illustration du filtre de Kalman
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Simulation d'une trajectoire
k_max = 20; % nombre d'iterations / d'observations
rng('default'); % initialisation du generateur de nombres aleatoires
g = std_g*randn(1,k_max); % bruit d'etat
w = std_w*randn(1,k_max); % bruit d'observation
x(:,1) = m_0; y(1) = C*x(:,1)+w(1);
for k = 1:k_max-1
    x(:,k+1) = A*x(:,k)+V*g(k);
    y(k+1) = C*x(:,k+1)+w(k+1);
end

% estimation de la vitesse par derivation de la position (derivee arriere)
x2_d(2:k_max) = (y(2:k_max)-y(1:k_max-1))/T;

m1 = min(x(1,:)); M1 = max(x(1,:)); m2 = min(x(2,:)); M2 = max(x(2,:));

% Figure 1 : Evolution de l'etat (position vs vitesse)
figure('Position', [100 100 800 600]);
plot(x(1,:),x(2,:),'r', 'LineWidth', 1.5);
title('Evolution de l''etat');
legend('Etat exact');
xlabel('Position (m)'); ylabel('Vitesse (m/s)');
axis([m1-1 M1+1 m2-0.3 M2+0.3]);
grid on;
saveas(gcf, fullfile(output_dir, 'fig1_etat_exact.png'));

% Figure 2 : Evolution position et vitesse en fonction du temps
figure('Position', [100 100 800 600]);
subplot(211);
plot(1:k_max,x(1,:),'r', 'LineWidth', 1.5);
title('Evolution de la position');
legend('Position exacte');
xlabel('k'); ylabel('Position (m)');
axis([1 k_max m1-1 M1+1]);
grid on;
subplot(212);
plot(1:k_max,x(2,:),'r', 'LineWidth', 1.5);
legend('Vitesse exacte');
title('Evolution de la vitesse');
xlabel('k'); ylabel('Vitesse (m/s)');
axis([1 k_max m2-0.3 M2+0.3]);
grid on;
saveas(gcf, fullfile(output_dir, 'fig2_position_vitesse_exactes.png'));

% Figure 3 : Comparaison etat exact vs observation
figure('Position', [100 100 800 600]);
subplot(211);
plot(1:k_max,x(1,:),'r', 1:k_max,y(1,:),'m', 'LineWidth', 1.5);
legend('Position exacte','Position observee');
title('Evolution de la position');
xlabel('k'); ylabel('Position (m)');
axis([0 k_max m1-1 M1+1]);
grid on;
subplot(212);
plot(1:k_max,x(2,:),'r', 2:k_max,x2_d(1,2:k_max),'m', 'LineWidth', 1.5);
legend('Vitesse exacte','Derivee de la position observee');
title('Evolution de la vitesse');
xlabel('k'); ylabel('Vitesse (m/s)');
axis([0 k_max m2-0.3 M2+0.3]);
grid on;
saveas(gcf, fullfile(output_dir, 'fig3_exact_vs_observation.png'));

% Estimation de la trajectoire par Filtre de Kalman
x_est(:,1) = m_0; P = P_0;
for k = 1:k_max-1
    x_pred(:,k+1) = A*x_est(:,k);
    P = A*P*A'+Q;
    K = P*C'*inv(C*P*C'+R);
    P = (I-K*C)*P;
    y_pred(k+1) = C*x_pred(:,k+1);
    x_est(:,k+1) = x_pred(:,k+1)+K*(y(k+1)-y_pred(k+1));
end

% Figure 4 : Etat exact vs estime par FK
figure('Position', [100 100 800 600]);
plot(x(1,:),x(2,:),'r', x_est(1,:),x_est(2,:),'b', 'LineWidth', 1.5);
legend('Etat exact','Etat estime par FK');
title('Evolution de l''etat');
xlabel('Position (m)'); ylabel('Vitesse (m/s)');
axis([m1-1 M1+1 m2-0.3 M2+0.3]);
grid on;
saveas(gcf, fullfile(output_dir, 'fig4_etat_exact_vs_FK.png'));

% Figure 5 : Position et vitesse estimees par FK
figure('Position', [100 100 800 600]);
subplot(211);
plot(1:k_max,x(1,:),'r', 1:k_max,x_est(1,:),'b', 'LineWidth', 1.5);
legend('Position exacte','Position estimee par FK');
title('Evolution de la position');
xlabel('k'); ylabel('Position (m)');
axis([0 k_max m1-1 M1+1]);
grid on;
subplot(212);
plot(1:k_max,x(2,:),'r', 1:k_max,x_est(2,:),'b', 'LineWidth', 1.5);
legend('Vitesse exacte','Vitesse estimee par FK');
title('Evolution de la vitesse');
xlabel('k'); ylabel('Vitesse (m/s)');
axis([0 k_max m2-0.3 M2+0.3]);
grid on;
saveas(gcf, fullfile(output_dir, 'fig5_FK_estimation.png'));

% Figure 6 : Comparaison FK vs observation vs exact
figure('Position', [100 100 800 600]);
subplot(211);
plot(1:k_max,x(1,:),'r', 1:k_max,x_est(1,:),'b', 1:k_max,y(1,:),'m', 'LineWidth', 1.5);
legend('Position exacte','Position estimee par FK','Position observee');
title('Evolution de la position');
xlabel('k'); ylabel('Position (m)');
axis([0 k_max m1-1 M1+1]);
grid on;
subplot(212);
plot(1:k_max,x(2,:),'r', 1:k_max,x_est(2,:),'b', 2:k_max,x2_d(1,2:k_max),'m', 'LineWidth', 1.5);
legend('Vitesse exacte','Vitesse estimee par FK','Derivee de la position observee');
title('Evolution de la vitesse');
xlabel('k'); ylabel('Vitesse (m/s)');
axis([0 k_max m2-0.3 M2+0.3]);
grid on;
saveas(gcf, fullfile(output_dir, 'fig6_comparaison_complete.png'));

% Figure 7 : Erreurs d'estimation
figure('Position', [100 100 800 600]);
subplot(211);
plot(1:k_max,x_est(1,:)-x(1,:),'b', 1:k_max,y-x(1,:),'m', 'LineWidth', 1.5);
legend('Erreur FK','Erreur observation');
title('Evolution de l''erreur de position');
xlabel('k'); ylabel('Erreur position (m)');
axis([0 k_max -3*std_w 3*std_w]);
grid on;
subplot(212);
plot(1:k_max,x_est(2,:)-x(2,:),'b', 2:k_max,x2_d(2:k_max)-x(2,2:k_max),'m', 'LineWidth', 1.5);
legend('Erreur FK','Erreur derivee observation');
title('Evolution de l''erreur de vitesse');
xlabel('k'); ylabel('Erreur vitesse (m/s)');
axis([0 k_max -3*sqrt(2)*std_w/T 3*sqrt(2)*std_w/T]);
grid on;
saveas(gcf, fullfile(output_dir, 'fig7_erreurs_estimation.png'));

% Estimation avec sauvegarde de P pour les ellipses
x_est(:,1) = m_0; P = P_0; P1_est(:,1) = P(:,1); P2_est(:,1) = P(:,2);
for k = 1:k_max-1
    x_pred(:,k+1) = A*x_est(:,k);
    P = A*P*A'+Q;
    K = P*C'*inv(C*P*C'+R);
    P = (I-K*C)*P; P1_est(:,k+1) = P(:,1); P2_est(:,k+1) = P(:,2);
    y_pred(k+1) = C*x_pred(:,k+1);
    x_est(:,k+1) = x_pred(:,k+1)+K*(y(k+1)-y_pred(k+1));
end

% Figure 8 : Evolution de l'ecart-type
figure('Position', [100 100 800 600]);
subplot(211);
plot(1:k_max,P1_est(1,:).^0.5,'b', 1:k_max,std_w*ones(1,k_max),'m', 'LineWidth', 1.5);
legend('Std position FK','Std observation position');
title('Precision de l''estimation de la position');
xlabel('k'); ylabel('Ecart-type (m)');
axis([0 k_max 0 3*std_w]);
grid on;
subplot(212);
plot(1:k_max,P2_est(2,:).^0.5,'b', 1:k_max,sqrt(2)*std_w/T*ones(1,k_max),'m', 'LineWidth', 1.5);
legend('Std vitesse FK','Std derivee observation');
title('Precision de l''estimation de la vitesse');
xlabel('k'); ylabel('Ecart-type (m/s)');
axis([0 k_max 0 3*sqrt(2)*std_w/T]);
grid on;
saveas(gcf, fullfile(output_dir, 'fig8_precision_estimation.png'));

% Figure 9 : Estimation avec intervalle de confiance a 3 sigma
figure('Position', [100 100 800 600]);
subplot(211);
plot(1:k_max,x(1,:),'r', 1:k_max,x_est(1,:),'b',...
     1:k_max,x_est(1,:)-3*P1_est(1,:).^0.5,':b',...
     1:k_max,x_est(1,:)+3*P1_est(1,:).^0.5,':b',...
     1:k_max,y,'g', 'LineWidth', 1.2);
legend('Position exacte','Position estimee FK',...
       'Estimation -3\sigma','Estimation +3\sigma',...
       'Position observee');
title('Estimation avec precision associee');
xlabel('k'); ylabel('Position (m)');
axis([0 k_max m1-1 M1+1]);
grid on;
subplot(212);
plot(1:k_max,x(2,:),'r', 1:k_max,x_est(2,:),'b',...
     1:k_max,x_est(2,:)-3*P2_est(2,:).^0.5,':b',...
     1:k_max,x_est(2,:)+3*P2_est(2,:).^0.5,':b',...
     2:k_max,x2_d(2:k_max),'g', 'LineWidth', 1.2);
legend('Vitesse exacte','Vitesse estimee FK',...
       'Estimation -3\sigma','Estimation +3\sigma',...
       'Derivee observation');
title('Estimation avec precision associee');
xlabel('k'); ylabel('Vitesse (m/s)');
axis([0 k_max m2-1 M2+1]);
grid on;
saveas(gcf, fullfile(output_dir, 'fig9_estimation_intervalle_confiance.png'));

% Figure 10 : Erreur avec bornes a 3 sigma
figure('Position', [100 100 800 600]);
subplot(211);
plot(1:k_max,x_est(1,:)-x(1,:),'b',...
     1:k_max,-3*P1_est(1,:).^0.5,':r', 1:k_max,3*P1_est(1,:).^0.5,':r', 'LineWidth', 1.5);
legend('Erreur FK','-3\sigma','+3\sigma');
title('Erreur d''estimation de la position');
xlabel('k'); ylabel('Erreur (m)');
axis([0 k_max -5*std_w 5*std_w]);
grid on;
subplot(212);
plot(1:k_max,x_est(2,:)-x(2,:),'b',...
     1:k_max,-3*P2_est(2,:).^0.5,':r', 1:k_max,+3*P2_est(2,:).^0.5,':r', 'LineWidth', 1.5);
legend('Erreur FK','-3\sigma','+3\sigma');
title('Erreur d''estimation de la vitesse');
xlabel('k'); ylabel('Erreur (m/s)');
axis([0 k_max -5*sqrt(2)*std_w/T 5*sqrt(2)*std_w/T]);
grid on;
saveas(gcf, fullfile(output_dir, 'fig10_erreur_bornes_3sigma.png'));

fprintf('\n=== Figures sauvegardees dans %s ===\n', output_dir);
fprintf('10 figures generees pour l''option 1 (Illustration du FK)\n');

% Fermer toutes les figures
close all;
