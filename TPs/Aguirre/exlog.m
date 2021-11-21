% Exemplo identificação da eq. logisitica
% (c) Luis Aguirre, 22/01/99

% Carregar dados
load exlog.dat;
y=exlog;
clear exlog;

% Montar matriz de regressores
reg1=ones(8,1);		% constante
reg2=y(12:19);		% y(k-1)
reg3=y(11:18);		% y(k-2)
reg4=reg2.^2;		% y(k-1)^2
reg5=reg2.*reg3;	% y(k-1)y(k-2)
reg6=reg3.^2;		% y(k-2)^2

psi=[reg1 reg2 reg3 reg4 reg5 reg6];
vec=y(13:20);
Psi=[psi vec];

load exlog.dat;
y=exlog;
clear exlog;

% Montar matriz de regressores
reg1=ones(8,1);		% constante
reg2=y(12:19);		% y(k-1)
reg3=y(11:18);		% y(k-2)
reg4=reg2.^2;		% y(k-1)^2
reg5=reg2.*reg3;	% y(k-1)y(k-2)
reg6=reg3.^2;		% y(k-2)^2

psi=[reg1 reg2 reg3 reg4 reg5 reg6];
vec=y(13:20);
Psi=[psi vec];

[A,err,piv]=myhouse(Psi,6);
% numero de parametros no modelo final
np=3;
Psit=Psi(:,piv(1:np));
teta=inv(Psit'*Psit)*Psit'*vec
% numero de parametros no modelo final
np=3;
Psit=Psi(:,piv(1:np));
teta=inv(Psit'*Psit)*Psit'*vec

pause

%eliminando o regressor "constante" 

Psit=Psi(:,piv(2:3));
teta=inv(Psit'*Psit)*Psit'*vec

	