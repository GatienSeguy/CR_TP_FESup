%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MANIPULATION 1 : Méthode des Moindres Carrés
% Script principal (version sans toolbox)
% Sauvegarde les figures pour le rapport LaTeX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc;

% Chemin pour sauvegarder les figures
figPath = '../LATEX/figures/';
if ~exist(figPath, 'dir')
    mkdir(figPath);
end

%% Paramètres du système (équivalent à systeme_stationnaire.m)
Ka = 1; tau = 0.1;
Te = 0.1; Fe = 1/Te;
a = -exp(-Te/tau);  % a = -0.3679
b = Ka*(1+a);       % b = 0.6321

% Fonction pour simuler le système H(z) = b*z^-1 / (1 + a*z^-1)
simulate_sys = @(a_val, b_val, e_val) filter([0 b_val], [1 a_val], e_val);

%% ========================================================================
% PARTIE 1 : Comparaison méthode des MC pour différents RSB
% =========================================================================
fprintf('\n=== PARTIE 1 : Méthode des Moindres Carrés ===\n');

N = 100;
t = (0:N-1)*Te; t = t(:);

% Signal d'entrée : chirp
f_min = 0; f_max = 0.25*Fe;
e = chirp(t, f_min, max(t), f_max); e = e(:);

% Test avec différents RSB
RSB_values = [100, 10, 1];

for idx = 1:length(RSB_values)
    RSB = RSB_values(idx);

    % Signal de sortie
    s = simulate_sys(a, b, e);
    Ps = sum(s.^2)/N;
    Pb = Ps/RSB; sigma_b = sqrt(Pb);
    s_b = s + sigma_b*randn(N,1);

    % Estimation avec la méthode des MC
    x_mc = s_b(2:N);
    C = [-s_b(1:N-1) e(1:N-1)];
    est = (C.'*C)\C.'*x_mc;
    a_est_mc = est(1);
    b_est_mc = est(2);

    % Sortie simulée avec les paramètres estimés
    s_est = simulate_sys(a_est_mc, b_est_mc, e);

    % Erreur d'estimation
    err_a = abs(a_est_mc - a);
    err_b = abs(b_est_mc - b);

    % Affichage
    fprintf('\nRSB = %d\n', RSB);
    fprintf('  Système réel :    a = %.4f, b = %.4f\n', a, b);
    fprintf('  Méthode MC :      a = %.4f, b = %.4f\n', a_est_mc, b_est_mc);
    fprintf('  Erreur :        |da| = %.4f, |db| = %.4f\n', err_a, err_b);

    % Figure des signaux
    fig = figure('Position', [100 100 900 600], 'Visible', 'off');

    subplot(2,1,1);
    stairs(t, e, 'r', 'LineWidth', 1); hold on;
    plot(t, s, 'b--', 'LineWidth', 1.5);
    plot(t, s_b, 'b', 'LineWidth', 0.5);
    title(sprintf('Signaux d''entrée et de sortie (RSB = %d)', RSB), 'FontSize', 12);
    legend('Entrée (chirp)', 'Sortie non bruitée', 'Sortie bruitée', 'Location', 'best');
    xlabel('Temps (s)'); ylabel('Amplitude');
    grid on;

    subplot(2,1,2);
    plot(t, s, 'b--', 'LineWidth', 1.5); hold on;
    plot(t, s_b, 'b.', 'MarkerSize', 4);
    plot(t, s_est, 'g', 'LineWidth', 1.5);
    title(sprintf('Comparaison sortie réelle et estimée (RSB = %d)', RSB), 'FontSize', 12);
    legend('Sortie exacte', 'Sortie bruitée', 'Sortie estimée MC', 'Location', 'best');
    xlabel('Temps (s)'); ylabel('Amplitude');
    grid on;

    saveas(fig, [figPath 'manip1_signaux_RSB' num2str(RSB) '.png']);
    close(fig);
end

%% ========================================================================
% PARTIE 2 : Étude statistique - Effet du bruit sur l'estimateur MC
% =========================================================================
fprintf('\n=== PARTIE 2 : Étude statistique de l''effet du bruit ===\n');

N = 100;
I = 500; % nombre de simulations

t = (0:N-1)*Te; t = t(:);

% Entrée chirp
f_min = 0; f_max = 0.25*Fe;
e = chirp(t, f_min, max(t), f_max); e = e(:);

% Sortie non bruitée
s = simulate_sys(a, b, e);

% Test pour différents RSB
RSB_test = [100, 10, 5, 2];

fig_stat = figure('Position', [100 100 1200 800], 'Visible', 'off');

for k = 1:length(RSB_test)
    RSB = RSB_test(k);
    Ps = sum(s.^2)/N;
    Pb = Ps/RSB; sigma_b = sqrt(Pb);

    % Stockage des résultats
    A_est_mc = zeros(I,1);
    B_est_mc = zeros(I,1);

    for i = 1:I
        s_b = s + sigma_b*randn(N,1);

        % Méthode des MC
        X_mc = s_b(2:N);
        C = [-s_b(1:N-1) e(1:N-1)];
        est = (C.'*C)\C.'*X_mc;
        A_est_mc(i) = est(1);
        B_est_mc(i) = est(2);
    end

    % Statistiques
    biais_a = mean(A_est_mc) - a;
    biais_b = mean(B_est_mc) - b;
    std_a = std(A_est_mc);
    std_b = std(B_est_mc);

    fprintf('\nRSB = %d (%d simulations)\n', RSB, I);
    fprintf('  Paramètre a : biais = %.4f, écart-type = %.4f\n', biais_a, std_a);
    fprintf('  Paramètre b : biais = %.4f, écart-type = %.4f\n', biais_b, std_b);

    % Histogrammes
    subplot(2, length(RSB_test), k);
    histogram(A_est_mc, 30, 'FaceColor', 'g', 'FaceAlpha', 0.7);
    xline(a, 'r', 'LineWidth', 2);
    xline(mean(A_est_mc), 'b--', 'LineWidth', 1.5);
    title(sprintf('RSB=%d : a', RSB));
    xlabel('a'); ylabel('Fréquence');
    legend('Estimations', 'Valeur vraie', 'Moyenne', 'Location', 'best');

    subplot(2, length(RSB_test), k + length(RSB_test));
    histogram(B_est_mc, 30, 'FaceColor', 'b', 'FaceAlpha', 0.7);
    xline(b, 'r', 'LineWidth', 2);
    xline(mean(B_est_mc), 'g--', 'LineWidth', 1.5);
    title(sprintf('RSB=%d : b', RSB));
    xlabel('b'); ylabel('Fréquence');
end

sgtitle('Distribution des estimations MC pour différents RSB', 'FontSize', 14);
saveas(fig_stat, [figPath 'manip1_etude_stat_RSB.png']);
close(fig_stat);

%% ========================================================================
% PARTIE 3 : Évolution des estimations au cours des simulations
% =========================================================================
fprintf('\n=== PARTIE 3 : Évolution des estimations ===\n');

RSB = 10;
Ps = sum(s.^2)/N;
Pb = Ps/RSB; sigma_b = sqrt(Pb);

A_est_mc = zeros(I,1);
B_est_mc = zeros(I,1);

for i = 1:I
    s_b = s + sigma_b*randn(N,1);
    X_mc = s_b(2:N);
    C = [-s_b(1:N-1) e(1:N-1)];
    est = (C.'*C)\C.'*X_mc;
    A_est_mc(i) = est(1);
    B_est_mc(i) = est(2);
end

fig = figure('Position', [100 100 1000 500], 'Visible', 'off');
subplot(1,2,1);
plot(1:I, a*ones(I,1), 'r', 'LineWidth', 2); hold on;
plot(1:I, A_est_mc, 'g.', 'MarkerSize', 4);
plot(1:I, cumsum(A_est_mc)./(1:I)', 'b', 'LineWidth', 1.5);
title('Estimation du paramètre a');
legend('Valeur exacte', 'Estimations', 'Moyenne cumulée');
xlabel('Numéro de simulation'); ylabel('a');
grid on;

subplot(1,2,2);
plot(1:I, b*ones(I,1), 'r', 'LineWidth', 2); hold on;
plot(1:I, B_est_mc, 'g.', 'MarkerSize', 4);
plot(1:I, cumsum(B_est_mc)./(1:I)', 'b', 'LineWidth', 1.5);
title('Estimation du paramètre b');
legend('Valeur exacte', 'Estimations', 'Moyenne cumulée');
xlabel('Numéro de simulation'); ylabel('b');
grid on;

sgtitle(sprintf('Méthode MC : %d simulations (RSB = %d)', I, RSB), 'FontSize', 14);
saveas(fig, [figPath 'manip1_evolution_estimations.png']);
close(fig);

fprintf('\nFigures sauvegardées dans %s\n', figPath);
fprintf('=== FIN MANIPULATION 1 ===\n');
