function x=rk4(x0,ux,uy,h,t)
%% Implementa Runge-Kutta de 4° ordem
% Autor: Luis Antonio Aguirre (https://www.researchgate.net/project/Scripts-on-Nonlinear-Dynamics)
% 
% This function implements a numerical integration algorithm known as 4th-order Runge-Kutta  
% x0 is the state vector BEFORE calling the function (i.e. the initial condition at each integration step) 
% ux and uy (if different from zero) are external forces ("control actions") added to the first and second
% state equations (in dvExample.m), respectively. It is assumed that such
% control actions do NOT change during the integrarion period h. If
% ux=uy=0, the autonomous case is simulated.
% h integration interval.
% t the time BEFORE calling the function.
% The vector field (the equations) are in dvCord.m

% LAA 15/8/18
%%

% 1st evaluation
xd=statespace_sa(x0,ux,uy,t);
savex0=x0;
phi=xd;
x0=savex0+0.5*h*xd;

% 2nd evaluation
xd=statespace_sa(x0,ux,uy,t+0.5*h);
phi=phi+2*xd;
x0=savex0+0.5*h*xd;

% 3rd evaluation
xd=statespace_sa(x0,ux,uy,t+0.5*h);
phi=phi+2*xd;
x0=savex0+h*xd;

% 4th evaluation
xd=statespace_sa(x0,ux,uy,t+h);
x=savex0+(phi+xd)*h/6;
