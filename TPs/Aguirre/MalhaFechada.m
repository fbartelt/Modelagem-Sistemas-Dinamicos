% Script de apoio aa aplicacao do metodo de Yuwana e Seborg (1982)
% aos dados bfg3.dat
%
% LAA 07/04/2020

clear
close all

% carrega dados
load mf643.dat
U=mf643(:,2);
Y=mf643(:,1);
clear mf643

figure(1)
plot(1:length(U),U,'k-.',1:length(Y),Y,'b-');
set(gca,'FontSize',16)
axis([0 length(U) 220 280])
xlabel('amostras')
ylabel('sinais medidos')


% print -djpeg mf643.jpg



% ajustar os dados, consideraremos que a entrada inicia na amostra 22.

% variacao de saida
Delta_Y=mean(Y(end-15:end))-mean(Y(1:21));
% variacao da entrada
Delta_U=mean(U(end-15:end))-mean(U(1:21));


% fazendo todos os ajuste, lembrando de retirar o offset tem-se
u=(U(21:end)-mean(U(1:21)));
y=(Y(21:end)-mean(Y(1:21)));


figure(2)
plot(1:length(u),u,'k-.',1:length(y),y,'b-');
set(gca,'FontSize',16)
axis([0 length(y) -60 10])
xlabel('amostras')
ylabel('sinais ajustados')

% print -djpeg mf643a.jpg


figure(3)
plot(1:50,y(1:50),'k');
set(gca,'FontSize',16)
axis([0 50 -60 10])
xlabel('amostras')
ylabel('saida ajustada')
grid


% print -djpeg mf643a_zoom.jpg

% y_infinity
y_inf=mean(y(end-5:end))
% yp1, se a resposta for positiva, trocar min por max
yp1=min(y);
% ym, se a resposta for positiva, trocar max por min. O ponto 10 foi obtido
% graficamente e corresponde a yp1.
ym=max(y(10:end));
% yp2, se a resposta for positiva, trocar min por max. O ponto 17 foi obtido
% graficamente e corresponde a ym.
yp2=min(y(17:end));
% Delta_t, obtido graficamente
Delta_t=13-6.5;

% Quem realizou o teste informou que o ganho proporcional do controlador
% era Kc=2
Kc=2;

%%%%% Calculo dos parametros do modelo de processo em malha aberta
K=y_inf/(Kc*(Delta_U-y_inf)); % Eq. 3.30
Kf=Kc*K; % relacao dada pelo metodo

% variaveis intermediarias
aux1=log((y_inf-ym)/(yp1-y_inf));
aux2=log((yp2-y_inf)/(yp1-y_inf));

zeta1=-aux1/(sqrt(pi^2+aux1^2)); % Eq. 3.33
zeta2=-aux2/(sqrt(4*pi^2+aux2^2)); % Eq. 3.34
zeta=(zeta1+zeta2)/2;

% mais variaveis intermediarias
aux3=zeta*sqrt(Kf+1);
aux4=sqrt(zeta^2*(Kf+1)+Kf);
aux5=sqrt((1-zeta^2)*(Kf+1));

tau=Delta_t*(aux3+aux4)*aux5/pi; % Eq. 3.31
taud=2*Delta_t*aux5/(pi*(aux3+aux4)); % Eq. 3.32

% para fins de simulacao aproximaremos o atraso puro de tempo pela sua
% aproximacao de Pade de 3a ordem

[numa,dena]=pade(taud,3);
atraso=tf(numa,dena);

% o modelo FOPDT eh
Planta=series(tf(K,[tau 1]),atraso);

step(Planta);
% print -djpeg step.jpg


% simulacao em malha fechada para validar
MF=feedback(Kc*Planta,1);
ymodelo=lsim(MF,u,1:length(u));


figure(4)
subplot(211)
plot(1:length(u),u,'k-.',1:length(y),y,'b-');
set(gca,'FontSize',16)
axis([0 length(y) -60 10])
ylabel('sinais ajustados')

subplot(212)
plot(1:length(u),u,'k-.',1:length(ymodelo),ymodelo,'r-');
set(gca,'FontSize',16)
axis([0 length(y) -60 10])
xlabel('amostras')
ylabel('simulacao')

% print -djpeg validacao.jpg






