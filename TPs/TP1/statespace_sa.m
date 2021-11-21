function xdot=statespace_sa(x,ux,uy,t)
%% Calcula as variaveis de estado para o atuador solenoide
% Codigo adaptado de Luis Antonio Aguirre (https://www.researchgate.net/project/Scripts-on-Nonlinear-Dynamics)
% Documentacao original:
% function xdot=dvCord(x,ux,uy,t)
% implements the Cord system
% x state vector
% the system has no inputs so ux=uy=0
% xd time derivative of x (vector field at x)


% x1 = I(corrente), x2 = z (posicao da armadura), x3 = dz/dt (velocidade da armadura)

% Parametros do circuito equivalente do solenoide
R = 3; % resistencia da armadura[ohm]
L = 5e-3; % Indutancia da bobina [H]
K = 6; % Permeabilidade magnetica [N/A^2]

% Parametros mecanicos do sistema
m = 0.03; % massa da valvula-armadura [kg]
b = 12; % atrito viscoso [Ns/m]
k = 6e3; % constante elastica [N/m]

xd(1) = -(R/L*x(1)) - (K/L*x(1)*x(3)) + (1/L*ux);
xd(2) = x(3);
xd(3) = -(k/m*x(2)) - (b/m*x(3)) + (K/(2*m)*(x(1)^2));

xdot=xd';