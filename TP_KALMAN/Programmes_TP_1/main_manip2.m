%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MANIPULATION 2 : MCR avec facteur d'oubli
% Script principal (version sans toolbox)
% Sauvegarde les figures pour le rapport LaTeX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc;

% Chemin pour sauvegarder les figures
figPath = '../LATEX/figures/';
if ~exist(figPath, 'dir')
    mkdir(figPath);
end

%% ========================================================================
% PARTIE 1 : MCR sans facteur d'oubli - Système stationnaire
% =========================================================================
fprintf('\n=== PARTIE 1 : MCR sur système stationnaire ===\n');

N = 100;
RSB = 10;

% Paramètres système stationnaire
Ka = 1; tau = 0.1;
Te = 0.1; Fe = 1/Te;
a = -exp(-Te/tau); b = Ka*(1+a);
t = (0:N-1)*Te; t = t(:);

% Entrée chirp
f_min = 0; f_max = 0.25*Fe;
e = chirp(t, f_min, max(t), f_max); e = e(:);

% Sortie (simulation manuelle du filtre)
s = filter([0 b], [1 a], e);
Ps = sum(s.^2)/N;
Pb = Ps/RSB; sigma_b = sqrt(Pb);
s_b = s + sigma_b*randn(N,1);

% Initialisation MCR
N1 = 4;
x = s_b(2:N1); C = [-s_b(1:N1-1) e(1:N1-1)];
est = (C.'*C)\C.'*x;
a_est = zeros(1,N); b_est = zeros(1,N);
a_est(N1) = est(1); b_est(N1) = est(2);
P = inv(C'*C);
P = diag([1^2, 20^2]);
P11 = zeros(1,N); P22 = zeros(1,N);
K1 = zeros(1,N); K2 = zeros(1,N);
P11(N1) = P(1,1); P22(N1) = P(2,2);

% Estimation récursive (lambda = 1, pas d'oubli)
for n = N1+1:N
    c = [-s_b(n-1) e(n-1)]';
    K = P*c/(c'*P*c+1);
    K1(n) = K(1); K2(n) = K(2);
    P = P - K*c.'*P;
    P11(n) = P(1,1); P22(n) = P(2,2);
    est = est + K*(s_b(n)-c.'*est);
    a_est(n) = est(1); b_est(n) = est(2);
end

fig = figure('Position', [100 100 900 600], 'Visible', 'off');
subplot(2,1,1);
stairs(t, e, 'r', 'LineWidth', 1); hold on;
plot(t, s, 'b--', 'LineWidth', 1.5);
plot(t, s_b, 'b', 'LineWidth', 0.5);
title('Signaux', 'FontSize', 12);
legend('Entrée', 'Sortie', 'Sortie bruitée', 'Location', 'best');
xlabel('Temps (s)'); grid on;

subplot(2,1,2);
plot(t, a*ones(1,N), 'r', 'LineWidth', 2); hold on;
plot(t(N1:N), a_est(N1:N), 'g', 'LineWidth', 1.5);
plot(t, b*ones(1,N), 'r--', 'LineWidth', 2);
plot(t(N1:N), b_est(N1:N), 'g--', 'LineWidth', 1.5);
title('Paramètres estimés par MCR (sans facteur d''oubli)', 'FontSize', 12);
legend('a exact', 'a estimé', 'b exact', 'b estimé', 'Location', 'best');
xlabel('Temps (s)'); grid on;
saveas(fig, [figPath 'manip2_MCR_sans_oubli.png']);
close(fig);

% Figure supplémentaire : évolution de P et K
fig2 = figure('Position', [100 100 1000 700], 'Visible', 'off');
subplot(2,2,1);
stairs(t, e, 'r', 'LineWidth', 1); hold on;
plot(t, s, 'b--', 'LineWidth', 1.5);
plot(t, s_b, 'b', 'LineWidth', 0.5);
title('Signaux'); legend('Entrée', 'Sortie', 'Sortie bruitée');
xlabel('Temps (s)'); grid on;

subplot(2,2,2);
plot(t, a*ones(1,N), 'r', 'LineWidth', 2); hold on;
plot(t(N1:N), a_est(N1:N), 'g', 'LineWidth', 1.5);
plot(t, b*ones(1,N), 'r--', 'LineWidth', 2);
plot(t(N1:N), b_est(N1:N), 'g--', 'LineWidth', 1.5);
title('Paramètres estimés');
legend('a exact', 'a estimé', 'b exact', 'b estimé');
xlabel('Temps (s)'); grid on;

subplot(2,2,3);
plot(t(N1:N), sqrt(P11(N1:N)), 'g', 'LineWidth', 1.5); hold on;
plot(t(N1:N), sqrt(P22(N1:N)), 'g--', 'LineWidth', 1.5);
title('Matrice P (écart-type)');
legend('\sigma_a = sqrt(P_{11})', '\sigma_b = sqrt(P_{22})');
xlabel('Temps (s)'); grid on;

subplot(2,2,4);
plot(t(N1+1:N), K1(N1+1:N), 'g', 'LineWidth', 1.5); hold on;
plot(t(N1+1:N), K2(N1+1:N), 'g--', 'LineWidth', 1.5);
title('Gain de Kalman K');
legend('K_1', 'K_2');
xlabel('Temps (s)'); grid on;

sgtitle('MCR : Évolution des paramètres et gains', 'FontSize', 14);
saveas(fig2, [figPath 'manip2_MCR_details.png']);
close(fig2);

%% ========================================================================
% PARTIE 2 : MCR avec facteur d'oubli - Système non stationnaire
% =========================================================================
fprintf('\n=== PARTIE 2 : MCR avec facteur d''oubli (système non stationnaire) ===\n');

N = 500;
RSB = 10;

% Système non stationnaire (saut de paramètres)
Ka1 = 1; Ka2 = 0.8;
tau1 = 0.1; tau2 = 0.08;
Ka_vec = [Ka1*ones(N/2,1); Ka2*ones(N/2,1)];
tau_vec = [tau1*ones(N/2,1); tau2*ones(N/2,1)];
Te = 0.1; Fe = 1/Te;
a_true = -exp(-Te./tau_vec);
b_true = Ka_vec.*(1+a_true);
t = (0:N-1)*Te; t = t(:);

% Entrée chirp
f_min = 0; f_max = 0.25*Fe;
e = chirp(t, f_min, max(t), f_max); e = e(:);

% Sortie du système non stationnaire
s = zeros(N,1);
for n = 2:N
    s(n) = b_true(n)*e(n-1) - a_true(n)*s(n-1);
end
Ps = sum(s(1:N/2-100).^2)/(N/2-100);
Pb = Ps/RSB; sigma_b = sqrt(Pb);
s_b = s + sigma_b*randn(N,1);

% Test avec différentes valeurs de lambda
lambda_values = [1, 0.98, 0.95, 0.9];

for idx = 1:length(lambda_values)
    lambda = lambda_values(idx);
    fprintf('Lambda = %.2f\n', lambda);

    % Initialisation
    N1 = 4;
    x = s_b(2:N1); C = [-s_b(1:N1-1) e(1:N1-1)];
    est = (C.'*C)\C.'*x;
    a_est = zeros(1,N); b_est = zeros(1,N);
    a_est(N1) = est(1); b_est(N1) = est(2);
    P = diag([1^2, 20^2]);

    % Estimation récursive avec facteur d'oubli
    for n = N1+1:N
        c = [-s_b(n-1) e(n-1)]';
        K = P*c/(c'*P*c + lambda);
        P = 1/lambda*(P - K*c.'*P);
        est = est + K*(s_b(n) - c.'*est);
        a_est(n) = est(1); b_est(n) = est(2);
    end

    fig = figure('Position', [100 100 900 600], 'Visible', 'off');
    subplot(2,1,1);
    stairs(t, e, 'r', 'LineWidth', 0.5); hold on;
    plot(t, s, 'b--', 'LineWidth', 1);
    plot(t, s_b, 'b', 'LineWidth', 0.3);
    title(sprintf('Signaux (\\lambda = %.2f)', lambda), 'FontSize', 12);
    legend('Entrée', 'Sortie', 'Sortie bruitée', 'Location', 'best');
    xlabel('Temps (s)'); grid on;

    subplot(2,1,2);
    plot(t, a_true, 'r', 'LineWidth', 2); hold on;
    plot(t(N1:N), a_est(N1:N), 'g', 'LineWidth', 1.5);
    plot(t, b_true, 'r--', 'LineWidth', 2);
    plot(t(N1:N), b_est(N1:N), 'g--', 'LineWidth', 1.5);
    xline(t(N/2), 'k--', 'LineWidth', 1);
    title(sprintf('Paramètres estimés par MCR (\\lambda = %.2f)', lambda), 'FontSize', 12);
    legend('a exact', 'a estimé', 'b exact', 'b estimé', 'Saut', 'Location', 'best');
    xlabel('Temps (s)'); grid on;

    lambda_str = strrep(num2str(lambda), '.', '_');
    saveas(fig, [figPath 'manip2_MCR_lambda' lambda_str '.png']);
    close(fig);
end

%% ========================================================================
% PARTIE 3 : Comparaison des différentes valeurs de lambda
% =========================================================================
fprintf('\n=== PARTIE 3 : Comparaison des valeurs de lambda ===\n');

fig = figure('Position', [100 100 1000 700], 'Visible', 'off');
colors = {'b', [0 0.5 0], 'm', 'c'};
legend_entries_a = {'Exact'};
legend_entries_b = {'Exact'};

for idx = 1:length(lambda_values)
    lambda = lambda_values(idx);

    % Initialisation
    N1 = 4;
    x = s_b(2:N1); C = [-s_b(1:N1-1) e(1:N1-1)];
    est = (C.'*C)\C.'*x;
    a_est = zeros(1,N); b_est = zeros(1,N);
    a_est(N1) = est(1); b_est(N1) = est(2);
    P = diag([1^2, 20^2]);

    for n = N1+1:N
        c = [-s_b(n-1) e(n-1)]';
        K = P*c/(c'*P*c + lambda);
        P = 1/lambda*(P - K*c.'*P);
        est = est + K*(s_b(n) - c.'*est);
        a_est(n) = est(1); b_est(n) = est(2);
    end

    subplot(2,1,1);
    if idx == 1
        plot(t, a_true, 'r', 'LineWidth', 2); hold on;
    end
    plot(t(N1:N), a_est(N1:N), 'Color', colors{idx}, 'LineWidth', 1.2);
    legend_entries_a{end+1} = sprintf('\\lambda=%.2f', lambda);

    subplot(2,1,2);
    if idx == 1
        plot(t, b_true, 'r', 'LineWidth', 2); hold on;
    end
    plot(t(N1:N), b_est(N1:N), 'Color', colors{idx}, 'LineWidth', 1.2);
    legend_entries_b{end+1} = sprintf('\\lambda=%.2f', lambda);
end

subplot(2,1,1);
xline(t(N/2), 'k--', 'LineWidth', 1);
title('Comparaison de l''estimation de a pour différentes valeurs de \lambda', 'FontSize', 12);
legend(legend_entries_a, 'Location', 'best');
xlabel('Temps (s)'); ylabel('a'); grid on;

subplot(2,1,2);
xline(t(N/2), 'k--', 'LineWidth', 1);
title('Comparaison de l''estimation de b pour différentes valeurs de \lambda', 'FontSize', 12);
legend(legend_entries_b, 'Location', 'best');
xlabel('Temps (s)'); ylabel('b'); grid on;

sgtitle('Influence du facteur d''oubli \lambda sur le suivi de paramètres', 'FontSize', 14);
saveas(fig, [figPath 'manip2_comparaison_lambda.png']);
close(fig);

%% ========================================================================
% PARTIE 4 : Étude statistique MCR avec facteur d'oubli
% =========================================================================
fprintf('\n=== PARTIE 4 : Étude statistique MCR ===\n');

I = 50; % nombre de simulations
lambda = 0.95;

fig = figure('Position', [100 100 1000 500], 'Visible', 'off');

subplot(1,2,1);
plot(t, a_true, 'r', 'LineWidth', 2); hold on;
title(sprintf('Estimation de a (\\lambda = %.2f)', lambda));
xlabel('Temps (s)'); ylabel('a'); grid on;

subplot(1,2,2);
plot(t, b_true, 'r', 'LineWidth', 2); hold on;
title(sprintf('Estimation de b (\\lambda = %.2f)', lambda));
xlabel('Temps (s)'); ylabel('b'); grid on;

for i = 1:I
    s_b_i = s + sigma_b*randn(N,1);

    % Initialisation
    N1 = 4;
    x = s_b_i(2:N1); C = [-s_b_i(1:N1-1) e(1:N1-1)];
    est = (C.'*C)\C.'*x;
    a_est = zeros(1,N); b_est = zeros(1,N);
    a_est(N1) = est(1); b_est(N1) = est(2);
    P = diag([1^2, 20^2]);

    for n = N1+1:N
        c = [-s_b_i(n-1) e(n-1)]';
        K = P*c/(c'*P*c + lambda);
        P = 1/lambda*(P - K*c.'*P);
        est = est + K*(s_b_i(n) - c.'*est);
        a_est(n) = est(1); b_est(n) = est(2);
    end

    subplot(1,2,1);
    plot(t(N1:N), a_est(N1:N), 'g', 'LineWidth', 0.3);

    subplot(1,2,2);
    plot(t(N1:N), b_est(N1:N), 'g', 'LineWidth', 0.3);
end

subplot(1,2,1);
plot(t, a_true, 'r', 'LineWidth', 2);
legend('Valeur exacte', 'Estimations', 'Location', 'best');

subplot(1,2,2);
plot(t, b_true, 'r', 'LineWidth', 2);
legend('Valeur exacte', 'Estimations', 'Location', 'best');

sgtitle(sprintf('MCR avec facteur d''oubli : %d simulations', I), 'FontSize', 14);
saveas(fig, [figPath 'manip2_etude_statistique.png']);
close(fig);

fprintf('\nFigures sauvegardées dans %s\n', figPath);
fprintf('=== FIN MANIPULATION 2 ===\n');
