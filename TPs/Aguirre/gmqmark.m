% investiga o problema de polarização do estimador GMQ
% 
% (c) Luis Antonio Aguirre, 1999

% Escolha o estimador atribuindo 1
GMQ=0;
MARKOV=1;

% Escolha o tipo do erro na equacao de regressao atribuindo 1
MA=1;
AR=0;

% condicoes iniciais
y=zeros(2,1);
if MARKOV==1
  teta=zeros(2,1);
else
  teta=zeros(3,1);
end;
tetapol=zeros(2,1);


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
for k=2:502
 y(k) = 0.7*y(k-1)+0.5*u(k-1)+e(k);
end;

% Estimador GMQ
if GMQ==1
% estimação de parametros (estimador GMQ)
% na primeira iteracao usa-se o estimador MQ

Psi=[y(1:501) u(1:501)];
V=y(2:502);
tetapol(:,j)=inv(Psi'*Psi)*Psi'*V;
xi=V-Psi*tetapol(:,j);


yfilt=y;
ufilt=u;

% iteracao para estimacao do modelo do ruído
for i=1:4
PsiRes=xi(1:499);
VRes=xi(2:500);
tetaRes=inv(PsiRes'*PsiRes)*PsiRes'*VRes;

% gera novas series filtradas
yfilt(2:500)=y(2:500)-tetaRes*y(1:499);
ufilt(2:500)=u(2:500)-tetaRes*u(1:499);

% Compoe matriz de regressores e faz estimacao com novos sinais filtrados
Psifilt=[yfilt(1:501) ufilt(1:501)];
Vfilt=yfilt(2:502);
tetaux=inv(Psifilt'*Psifilt)*Psifilt'*Vfilt;
xi=V-Psi*tetaux;
end
teta(:,j)=[tetaux; tetaRes];
end;

% estimador de MARKOV
if MARKOV==1
Psi=[y(1:100) u(1:100)];
V=y(2:101);

R=diag(ones(100,1)*(1+0.8^2));
for c=2:100
  R(c,c-1)=0.8;
  R(c-1,c)=0.8;
end;
R=std(nu)^2*R;
M1=Psi'*inv(R);
teta(:,j)=inv(M1*Psi)*M1*V;
end;


end

if GMQ==1
subplot(221)
hist(tetapol(1,:),100);
title('(a)');
subplot(222)
hist(teta(1,:),100);
title('(b)');
subplot(223)
hist(tetapol(2,:),100);
title('(c)')
subplot(224)
hist(teta(2,:),100);
title('(d)');


figure(2)
hist(teta(3,:),100);
pause
end;


if MARKOV==1
hist(teta(1,:),100);
axis([0.55 0.85 0 18])
figure(2)
hist(teta(2,:),100); 
pause
end;



