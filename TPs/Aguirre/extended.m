% Exemplo Chua
% (c) Luis Aguirre, 12/10/99

% Simular o sistema
e=randn(500,1)*0.05;
u=randn(500,1);
y(1)=0;
y(2)=0;
for k=3:500
  y(k)=0.4*y(k-1)-0.3*exp(-u(k-2)^2)*y(k-1)+0.5*y(k-2)+u(k-1)+e(k);
end;

% Montar matriz de regressores
reg1=ones(488,1);	% constante
reg2=y(12:499);		% y(k-1)
reg3=y(11:498);		% y(k-2)
reg4=reg2.^2;		% y(k-1)^2
reg5=reg2.*reg3;	% y(k-1)y(k-2)
reg6=reg3.^2;		% y(k-2)^2
reg7=u(12:499);		% u(k-1)
reg8=u(11:498);		% u(k-2)
reg9=reg7.^2;		% u(k-1)^2
reg10=reg7.*reg8;	% u(k-1)u(k-2)
reg11=reg8.^2;		% u(k-2)^2
reg12=exp(-reg4).*reg2; % exp(-y(k-1)^2)y(k-1)
reg13=exp(-reg4).*reg3; % exp(-y(k-1)^2)y(k-2)
reg14=exp(-reg6).*reg2; % exp(-y(k-2)^2)y(k-1)
reg15=exp(-reg6).*reg3; % exp(-y(k-2)^2)y(k-2)
reg16=exp(-reg9').*reg2; % exp(-u(k-1)^2)y(k-1)
reg17=exp(-reg9').*reg3; % exp(-u(k-1)^2)y(k-2)
reg18=exp(-reg11').*reg2; % exp(-u(k-2)^2)y(k-1)
reg19=exp(-reg11').*reg3; % exp(-u(k-2)^2)y(k-2)
reg20=exp(-reg10').*reg2; % exp(-u(k-1)u(k-2))y(k-1)
reg21=exp(-reg10').*reg3; % exp(-u(k-1)u(k-2))y(k-2)

psi=[reg1 reg2' reg3' reg4' reg5' reg6' reg7 reg8 reg9 reg10 reg11];
psi=[psi reg12' reg13' reg14' reg15' reg16' reg17' reg18' reg19' reg20' reg21'];
vec=y(13:500)';
Psi=[psi vec];

[A,err,piv]=myhouse(Psi,21);
% numero de parametros no modelo final
np=4;
Psit=Psi(:,piv(1:np));
teta=inv(Psit'*Psit)*Psit'*vec

