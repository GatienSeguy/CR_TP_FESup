function nuage_u(A,V,C,std_g,std_w,m_0,P_0,Q,R,k_max,N,fig);

% nuage_u(A,V,C,std_g,std_w,m_0,P_0,Q,R,k_max,N, fig);
%
% Cette fonction effectue une étude statistique pour la kième itération
% de la trajectoire du robot avec des bruits uniformes.
%
% Entrées : A     matrice d'état
%           V     matrice du bruit d'état
%           C     matrice d'observation
%           std_g écart type du bruit d'état
%           std_w écart type du bruit d'observation
%           m_0   état initial
%           P_0   matrice de covariance de l'état initial
%           Q     matrice de covariance du bruit d'état pour le FK
%           R     matrice de covariance du bruit d'observation pour le FK
%           k_max itération de la trajectoire pour l'étude statistique
%           N     nombre de simulations pour la simulation statistique
%           fig   numéro de la figure pour les tracés
%
% Cécile Durieu, 2 avril 2019

I = eye(2);
rng('default');
for n = 1:N;
    g = std_g*(rand(1,k_max)-0.5)*2*sqrt(3); % bruit d'état
    w = std_w*(rand(1,k_max)-0.5)*2*sqrt(3); % bruit d'observation
    x(:,1) = m_0; y(1) = C*x(:,1)+w(1);
    for k = 1:k_max-1;
        x(:,k+1) = A*x(:,k)+V*g(k);
        y(k+1) = C*x(:,k+1)+w(k+1);
    end
    x_est(:,1) = m_0; P = P_0;
    for k = 1:k_max-1;
        x_pred(:,k+1) = A*x_est(:,k);
        P = A*P*A'+Q;
        K = P*C'*inv(C*P*C'+R);
        P = (I-K*C)*P;
        y_pred(k+1) = C*x_pred(:,k+1);
        x_est(:,k+1) = x_pred(:,k+1)+K*(y(k+1)-y_pred(k+1));
    end
    e_est(:,n) = x_est(:,k+1)-x(:,k+1);
    eqm1(n) = (x_est(1,k+1)-x(1,k+1)).^2;
    eqm2(n) = (x_est(2,k+1)-x(2,k+1)).^2;
end;
figure(fig); plot(e_est(1,:),e_est(2,:),'.b');
title('erreur d''estimation'); xlabel('position'); ylabel('vitesse')
hold on; trace_ellipse_P([0 0],9*P,'r');
axis([-5*P(1,1).^0.5 5*P(1,1).^0.5 -5*P(2,2).^0.5 5*P(2,2).^0.5]);

moy_e = mean(e_est'); cov_e = cov(e_est');
text = [' covariance initiale :',num2str(mean(P))];
text_moy = [' moyenne de l''erreur :       ',num2str(mean(moy_e(1))),'    ',num2str(mean(moy_e(2)))];
text_C1  = [' covariance de l''erreur :     ',num2str(cov_e(1,1)),'    ',num2str(cov_e(1,2))];
text_C2  = ['                              ',num2str(cov_e(2,1)),'    ',num2str(cov_e(2,2))];
text_eqm = [' eqm :                        ',num2str(mean(eqm1)),'    ',num2str(mean(eqm2))];
text_P1  = [' covariance calculée par FK : ',num2str(P(1,1)),'    ',num2str(P(1,2))];
text_P2  = ['                              ' ,num2str(P(2,1)),'    ',num2str(P(2,2))];
disp(' '); disp(text_moy); disp(text_C1);disp(text_C2); disp(text_eqm); disp(text_P1);disp(text_P2);