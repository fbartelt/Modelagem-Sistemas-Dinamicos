function xdot=dvCord(x,ux,uy,t)
%
% function xdot=dvCord(x,ux,uy,t)
% implements the Cord system
% x state vector
% the system has no inputs so ux=uy=0
% xd time derivative of x (vector field at x)

% The cord system

%
% First published in
%
% Aguirre, L.A., Letellier, C., Investigating observability properties from data in nonlinear dynamics�, 
% Phisical Reivew E, 83:066209, 2011. DOI: 10.1103/PhysRevE.83.066209.
%
% Analysed in
%
% Letellier, C., Aguirre, L.A., Required criteria for recognizing new types of chaos: Application to the �cord� attractor�, 
% Phisical Reivew E, 85:036204, 2012. DOI: 10.1103/PhysRevE.85.036204.

a=0.258;
b=4.033;
F=8;
G=1;

xd(1)=-x(2)-x(3)-a*x(1)+a*F;
xd(2)=x(1)*x(2)-b*x(1)*x(3)-x(2)+G;
xd(3)=b*x(1)*x(2)+x(1)*x(3)-x(3);

xdot=xd';