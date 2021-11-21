% Alguns exemplos do capitulo 4

% (c) Luis Aguirre, 1999, 2011, 2012

%%
clear
close all

% O sistema
u=prbs(378,6,1);
u=u-0.5;
lu=length(u);
y=dlsim([1 0.5],[1 -1.5 0.7],u);
e=randn(378,1);
% ajustar variancia para dar aproximadamente 0,1
e=0.3*(e-mean(e));
yi=dimpulse([1 0.5],[1 -1.5 0.7],lu);

figure(1)
subplot(211)
plot(1:100,u(1:100),'k-',1:100,u(1:100),'ko');
axis([0 100 -0.7 0.7]);
title('(a)');
subplot(212)
plot(1:100,y(1:100),'k-',1:100,y(1:100),'ko');
title('(b)');
xlabel('amostras')
axis([0 100 -6 6]);



[t,ruy,l,B1]=myccf([y(378-314:378)+e(378-314:378) u(378-314:378)'],40,0,0,'k');
%assumindo que a entrada eh perfeitamente aleatoria, ou seja, Ruu diagonal.
h=ruy/var(u);

figure(2)
set(gca,'FontSize',18)
plot(1:41,yi(1:41),'k-',1:41,yi(1:41),'ko',1:41,h(1:41)*B1,'k+');
%title('(a)');

axis([0 40 -1 3]);
xlabel('amostras');

[t,ruu,l,B]=myccf([u(378-314:378)' u(378-314:378)'],80,1,0,'k');
% sem assumir que Ruu eh diagonal
Ruu=B*ruu(41:end)';
for i=1:length(ruy)-1
    Ruu=[Ruu B*ruu(41-i:end-i)'];
end
h2=inv(Ruu)*ruy';

figure(3)
set(gca,'FontSize',18)
plot(1:41,yi(1:41),'k-',1:41,h(1:41)*B1,'ko',1:41,h2(1:41),'k+');
axis([0 40 -1 3]);
xlabel('amostras');
% print -deps2c /Users/aguirre/mydocuments/latex/book2/figures/iddimp4.eps



[t,ruu,l,B]=myccf([u(378-314:378)' u(378-314:378)'],20,0,0,'k');
figure(4)
subplot(211)
set(gca,'FontSize',18)
plot([0 20],[l l],'k--',[0 20],[-l -l],'k--');
axis([-0.5 20.5 -1.1 1.1]);
hold on
stem(t,ruu,'k');
hold off
xlabel('atraso');
title('(a)');

[t,rue,l,B]=myccf([u(378-314:378)' e(378-314:378)],40,1,0,'k');
subplot(212)
set(gca,'FontSize',18)
plot([0 20],[l l],'k--',[0 20],[-l -l],'k--');
axis([-0.5 20.5 -1.1 1.1]);
hold on
stem(t,rue,'k');
hold off
xlabel('atraso');
title('(b)');
% print -deps2c /Users/aguirre/mydocuments/latex/book2/figures/ccfs2.eps

%%
%%%% Dominio da frequencia %%%%%
N=512;
u=prbs(N,11,1);
u=u-0.5;
lu=length(u);
y=dlsim([1 0.5],[1 -1.5 0.7],u);
e=randn(N,1);
Y=fft(y+e);
U=fft(u');
H=Y./U;

freq=1/(length(y))*(0:length(y)/2); 
semilogx(freq,20*log10(abs(H(1:length(freq))))); 

w=logspace(-1.4,0.4,100);
mw=max(w);
I=find(2*pi*freq<=mw);
freq=freq(I);
[mag,pha]=dbode([1 0.5],[1 -1.5 0.7],1,w);
subplot(211)
set(gca,'FontSize',14)
semilogx(w,20*log10(mag),'k-',2*pi*freq,20*log10(abs(H(1:length(freq)))),'k--');
title('(a)');
ylabel('|H(j\omega)| (dB)');
axis([0.04 2.5 -20 30]);
subplot(212)
set(gca,'FontSize',14)
semilogx(w,pha,'k-',2*pi*freq,angle(H(1:length(freq)))*180/pi,'k--');
title('(b)');
ylabel('fase de H(j\omega) (graus)');
xlabel('frequencia (rad/s)')
axis([0.04 2.5 -200 50]);



[t,ruy,l,B]=myccf([y+e u'],N,0,0,'k');
Ruy=fft(ruy*B);
[t,ruu,l,B]=myccf([u' u'],N,0,0,'k');
Ruu=fft(ruu*B);
Hr=Ruy./Ruu;

figure(2)
w=logspace(-1.4,0.4,100);
mw=max(w);
I=find(2*pi*freq<=mw);
freq=freq(I);
[mag,pha]=dbode([1 0.5],[1 -1.5 0.7],1,w);
subplot(311)
set(gca,'FontSize',14)
semilogx(w,20*log10(mag),'k-',2*pi*freq,20*log10(abs(Hr(1:length(freq)))),'k--');
title('(a)');
ylabel('|H(j\omega)| (dB)');
axis([0.04 2.5 -20 30]);
subplot(312)
set(gca,'FontSize',14)
semilogx(w,pha,'k-',2*pi*freq,angle(Hr(1:length(freq)))*180/pi,'k--');
title('(b)');
ylabel('fase de H(j\omega) (graus)');
axis([0.04 2.5 -200 50]);


% calculo da funcao de coerencia (18/08/2011)

[t,ryy,l,B]=myccf([y+e y+e],N,0,0,'k');
Ryy=fft(ryy*B);
guy=sqrt((abs(Ruy).^2)./(abs(Ruu).*abs(Ryy)));
subplot(313)
set(gca,'FontSize',14)
semilogx(2*pi*freq,guy(1:length(freq)),'k',[freq(2) max(2*pi*freq)],[0.6 0.6],'k');
axis([0.04 2.5 0 3]);
xlabel('frequencia (rad/s)')
ylabel('\gamma_{uy}')
title('(c)');

