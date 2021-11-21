% Modelagem de sitemas com modelos FOPDT (First-Order Plus Dead Time)
%
% LAA 31/3/2020

% sistema de alta ordem SEM atraso de tempo
num=50;
% dois polos em s=-1
den1=conv([1 1],[1 1]);

% dois polos em s=-2
den2=conv([1 2],[1 2]);

% dois polos em s=-3
den3=conv([1 3],[1 3]);

% sistema com denominador de 6a ordem
den=conv(den1,den2);
den=conv(den,den3);

t=0:0.1:15;
y=step(num,den,t);

figure(1)
plot(t,y);
set(gca,'FontSize',16)
grid
xlabel('tempo')

% atraso puro de tempo aprox. 1,5
Taud=1.6;
% ganho, como a entrada foi um degrau unitario, basta pegar o valor da
% saida em estado estacionario
K=y(end);
% estimando a constante de tempo dominante. 63% do valor final eh:
y63=y(end)*0.63;
% busca os indices do vetor y para os quais y>y63
I=find(y>y63);
% o primeiro indice corresponde aproximadamente a y=y63
% o tempo correspondente eh
t63=t(I(1));
% a constante de tempo deve descontar o tempo morto, logo
tau=t63-Taud;

% simulação do modelo ajustado dessa forma SEM atraso
ym=step(K,[tau 1],t);

% para comparar, vamos ajustar o atraso na mao:
figure(2)
plot(t,y,'k',t+Taud,ym,'r');
axis([0 15 0 1.5])
set(gca,'FontSize',16)
xlabel('tempo')

% ajuste fino. Aumentamos um pouco o atraso
Tauda=Taud+0.4;
% reduzimos a constante de tempo
taua=tau*0.7
% simulação do modelo ajustado dessa forma SEM atraso
yma=step(K,[taua 1],t);

figure(3)
plot(t,y,'k',t+Tauda,yma,'r');
axis([0 15 0 1.5])
set(gca,'FontSize',16)
xlabel('tempo')

% margem de fase do sistema
SYS=tf(num,den);
[Gm,Pm,Wcg,Wcp] = margin(SYS);

% fase devida ao atraso puro de tempo em Wcp
j=sqrt(-1);
fase_atraso_em_Wcp=phase(exp(-j*Tauda*Wcp))*180/pi;
% fase devida ao modelo de 1a ordem sem atraso, em Wcp
MOD=tf(K,[taua 1]);
[mag,pha]=bode(MOD,Wcp);
% margem de fase do modelo FOPDT
180+pha+fase_atraso_em_Wcp

figure(4)
bode(SYS,MOD,logspace(-1,1))
xlabel('rad/s')
grid

% 
w=logspace(-0.5,0.5);
pha=(phase(1./(j.*w*taua +1))+phase(exp(-j.*w*Tauda)))*180/pi;
% busca os indices do vetor pha para os quais pha<-180
I=find(pha<-180);
% o primeiro indice corresponde aproximadamente a pha=-180
% a frequencia onde leremos a margem de ganho eh aproximadamente
wcg=w(I(1));
% o ganho do modelo nessa frequencia eh
mag=abs(1./(j.*wcg*taua +1))
% portanto a margem de ganho eh
MG=1/mag

% print -djpeg fopdt.jpg
% print -djpeg fopdt2.jpg
% print -djpeg fopdt3.jpg
% print -djpeg fopdt4.jpg
 