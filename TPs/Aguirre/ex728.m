% Gera figura do Exemplo 7.2.8

% LAA 13/05/2014
%%
clear
close all

N=1502;

% inicializacao
xi1=zeros(N,1);
xi2=zeros(N,1);

% condicoes iniciais
y1=zeros(2,1);

% entradas aleatorias
u1=randn(N,1);

% coeficientes do polinomio C
a1=-0.3;
a2=0.6;

for n=1:400
% ruido 
nu=randn(N,1)*0.3;

% sistemas
% sistema que sera usado com erro de observacao em y
for k=3:N
   y1(k) = a1*y1(k-1)+a2*y1(k-2)+ 0.5*u1(k-1);
end;
y1=y1+nu;


% estimacao de parametros
% primeira iteracao, MQ
Psi1=[y1(2:N-1) y1(1:N-2) u1(2:N-1)];
teta1=inv(Psi1'*Psi1)*Psi1'*y1(3:N);
xi1(3:N)=[y1(3:N)-Psi1*teta1];

% parametros estimados por MQ classico
Teta1(n)=teta1(1,end);
Teta2(n)=teta1(2,end);
Teta3(n)=teta1(3,end);

Psi2=[y1(2:N-2) y1(1:N-3) u1(2:N-2) xi1(2:N-2) xi1(1:N-3)];
teta2=inv(Psi2'*Psi2)*Psi2'*y1(3:N-1);
xi2(3:N-1)=y1(3:N-1)-Psi2*teta2;
%xi2=xi2';

for k=1:10
    Psi3=[y1(2:N-3-k) y1(1:N-4-k) u1(2:N-3-k) xi2(2:N-3-k) xi2(1:N-4-k)];
    teta3(:,k)=inv(Psi3'*Psi3)*Psi3'*y1(3:N-2-k);
    xi2(3:N-2-k)=y1(3:N-2-k)-Psi3*teta3(:,k);
end;

% parametros estimados por EMQ
Teta1(n)=teta3(1,end);
Teta2(n)=teta3(2,end);
Teta3(n)=teta3(3,end);
Teta4(n)=teta3(4,end);
Teta5(n)=teta3(5,end);
end

figure(1)
subplot(231)
set(gca,'FontSize',14)
hist(-Teta1,50)
axis([0.2 0.4 0 25])
ylabel('(a)');
subplot(232)
set(gca,'FontSize',14)
hist(-Teta2,50)
axis([-0.8 -0.4 0 25])
ylabel('(b)');
subplot(233)
set(gca,'FontSize',14)
hist(Teta3,50)
axis([0.45 0.55 0 25])
ylabel('(c)');
subplot(234)
set(gca,'FontSize',14)
hist(Teta4,50)
axis([0.2 0.4 0 25])
ylabel('(d)');
subplot(235)
set(gca,'FontSize',14)
hist(Teta5,50)
axis([-0.8 -0.4 0 25])
ylabel('(e)');


%%
% Verificacao da condicao de "positividade real" Eq. 7.17

% vetor de frequencias
w=logspace(-2,1);
% polinomio C
c1=-a1;
c2=-a2;
Cw=1+c1*exp(j.*w)+c2*exp(j*2.*w);
% teste
figure(2)
set(gca,'FontSize',14)
semilogx(w,real(1./Cw),'k',w,0.5*ones(size(w)),'k--');
xlabel('\omega')
ylabel('Re[1/C(e^{j\omega} )]');


