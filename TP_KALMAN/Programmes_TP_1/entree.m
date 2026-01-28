%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Signal d'entrée                                                         %
%                                                                         %
% Cécile DURIEU, 25 janvier 2025                                          %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

e = input(' choix du signal d''entrée : carré (1), chirp (2), SBPA (3) : ');
if e == 1;
   e = square(2*pi*Fe/25*t);
elseif e ==2;
   f_min = 0; f_max = 0.25*Fe; e = chirp(t,f_min,max(t),f_max);
else e = idinput(N); end;
e = e(:);