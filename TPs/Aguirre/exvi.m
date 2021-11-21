% investiga o problema de polarização do estimador VI
% 
% (c) Luis Antonio Aguirre, 1999


% Escolha o tipo do erro na equacao de regressao atribuindo 1
MA=0;
AR=1;
% ruido de observacao
RO=0;

% Escolha do tipo de instrumentos. P1=predicoes de um passo aa frente
% e U=sinais de entrada
P1=0;
U=1;

% condicoes iniciais
y=zeros(2,1);
teta=zeros(2,1);


% determina o numero de realizacoes
for j=1:400
   
% nivel de ruído usado 
e=0.2*randn(502,1);
nu=0.2*randn(502,1);

if MA==1
% geração de ruído MA
for k=2:length(nu)
  e(k)=0.8*nu(k-1)+nu(k);
end;
end;
if AR==1
% geração de ruído AR
for k=2:length(nu)
  e(k)=0.8*e(k-1)+nu(k);
end;
end

% entrada aleatoria
u=randn(502,1);

% sistema com erro de modelagem
if RO==0
for k=2:502
 y(k) = 0.7*y(k-1)+0.5*u(k-1)+e(k);
end;
end;


% sistema com ruido de observacao
if RO==1
for k=2:502
 y(k) = 0.7*y(k-1)+0.5*u(k-1);
end;
y=y+nu;
end;

% monta matriz de regressores e vetor de observacoes
Psi=[y(1:501) u(1:501)];
V=y(2:502);
tetapol(:,j)=inv(Psi'*Psi)*Psi'*V;

% monta matriz de instrumentos 
if P1==1
 y1pf=Psi*tetapol(:,j);
 Z=[y1pf(1:499) u(2:500)];
elseif U==1
 Z=[u(2:500) u(1:499)];
end;

teta(:,j)=inv(Z'*Psi(2:500,:))*Z'*V(2:500);

end;
 
subplot(221)
hist(tetapol(1,:),100);
title('(a)');
subplot(222)
hist(teta(1,:),100);
title('(b)');
subplot(223)
hist(tetapol(2,:),100);
axis([0.45 0.55 0 15])
title('(c)')
subplot(224)
hist(teta(2,:),100);
axis([0.45 0.55 0 15])
title('(d)');
