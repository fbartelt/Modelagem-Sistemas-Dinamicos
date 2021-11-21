% Exemplo Chua
% (c) Luis Aguirre, 26/01/99

% Carregar dados
load dsvc1.dat;
y=dsvc1(2001:3000-87,2);
y=y./max(y);
clear dsvc1;

% Montar matriz de regressores
reg1=ones(998-87,1);	% constante
reg2=y(2:999-87);	% y(k-1)
reg3=y(1:998-87);	% y(k-2)
reg4=reg2.^2;		% y(k-1)^2
reg5=reg2.*reg3;	% y(k-1)y(k-2)
reg6=reg3.^2;		% y(k-2)^2
reg7=reg4.*reg2;  	% y(k-1)^3
reg8=reg6.*reg3;  	% y(k-2)^3
reg9=reg4.*reg3;  	% y(k-1)^2y(k-2)
reg10=reg2.*reg6; 	% y(k-1)y(k-2)^2

psi=[reg1 reg2 reg3 reg4 reg5 reg6 reg7 reg8 reg9 reg10];
vec=y(3:1000-87);
Psi=[psi vec];

[A,err,piv]=myhouse(Psi,10);
% numero de parametros no modelo final
np=10;
Psit=Psi(:,piv(1:np));
teta=inv(Psit'*Psit)*Psit'*vec


if 1>2
% Analise de agrupamentos de termos para Dupla volta
% Agrupamento constante
S0=teta(find (piv==1));
% Agrupamento linear
S1=sum(teta(find (piv>1&piv<4)));
% Agrupamento quadratico
S2=sum(teta(find (piv>3&piv<7)));
% Agrupamento cubico
S3=sum(teta(find (piv>6&piv<11)));
end;
