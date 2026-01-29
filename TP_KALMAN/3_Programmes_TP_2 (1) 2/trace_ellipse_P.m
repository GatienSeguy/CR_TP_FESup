function trace_ellipse_P(Centre,P,couleur)
  
% trace_ellipse_P(Centre,P,couleur)
% 
% Cette fonction trace une ellipse caracterisée par :
%    - son centre Centre = [xc yc]'
%    - sa matrice (covariance) P.
% L'équation d'un point (x,y) appartenant à l'ellipse est :
%      [x-xc y-yc]'*inv(P)*[x-xc y-yc] = 1
% La couleur et le style du trait sont mposéss par couleur.
%
% Cécile Durieu, 18 avril 1995

% Centre doit être un vecteur colonne
s = size(Centre); if s(1) < s(2); Centre = Centre'; end

N = 100;                        % Nombre de points pour tracer l'ellipse

% Changement de repere 
[V,D] = eig(P);                 % Recherche des valeurs propres et des
                                % vecteurs propres de la matrice P
teta_0 = atan2(V(2,1),V(1,1));  % Une des directions principales de l'ellipse
R = [cos(teta_0) -sin(teta_0)
     sin(teta_0) cos(teta_0)];  % Matrice de rotation % axe principal retenu

% Trace de l'ellipse
teta = (0:N)*2*pi/N;        
pts= Centre*ones(1,N+1) + R*sqrt(D)*[cos(teta);sin(teta)];
plot(pts(1,:),pts(2,:),eval('couleur'));
