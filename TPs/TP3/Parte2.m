%Felipe Bartelt de Assis Pessoa - 2016026841
% Script adaptado das rotinas fornecidas por Luis Antonio Aguirre

close all
clear

%% a)
tau = 0.0165;
u = prbs(1024, 10, ceil(tau));
u = u - mean(u);

lu=length(u);

figure(1)
subplot(211)
plot(1:lu,u(1:lu));
axis([0 lu -0.6 0.6]);
title('Sinal PRBS');
xlabel('Amostras')

lag = 30;
[t,ruu,l,~]=myccf([u' u'],lag,0,0,'k');
subplot(212)
plot([0 lag],[l l],'b--',[0 lag],[-l -l],'b--');
axis([-0.5 lag+0.5 -0.2 1.1]);
hold on
stem(t,ruu, 'k');
hold off
xlabel('Atraso');
title('Funcao de autocorrelacao');

%% b)

t0 = 0; % Tempo inicial
tff = 10; % Tempo final
h = 0.0001; % Intervalo de integracao
t = t0:1:1033; % Vetor de tempo
x0 = [0;0;0]; %cond iniciais
% x1 - corrente, x2 - posicao armadura, x3 - velocidade armadura
x = [x0 zeros(length(x0),length(t)-1)];
xt = [x0 zeros(length(x0),length(t)-1)];
e_in = zeros(length(t), 1);
w = (pi/2 + 2*pi)/0.5; % 5pi rad (permite 2 ciclos ate t=0.5s)
e_in(t>=10) = 10*u';  % degrau em t=0.5, amplitude=10
u = e_in;
ut = ones(size(u))*5;
ut(t>1) = 0;

for k=2:length(t)
    x(:,k)=rk4(x(:,k-1),u(k),u(k),h,t(k));
    xt(:,k)=rk4(xt(:,k-1),ut(k),ut(k),h,t(k));
end

y = x(2,:);
noise = 0.05 * y .* randn(size(y));
yn = y + noise;
figure(2)
myccf([yn' u],100, 1, 1,'k');
axis([-51 51 -0.15 0.15])
title('FCC')

%% c)

[tt,ruy,l,B1] = myccf([yn' u],50, 0, 0,'k');
h=ruy/var(u);

figure(3)
stem(1:length(h), h*B1)
hold on
stem(1:length(h), xt(2, 1:length(h)))
hold off
%title('(a)');
xlim([0 length(h)+0.5])
title('Aproximacao do sistema')
xlabel('Amostras');
ylabel('h(k)');
legend('Estimado', 'Real')

% figure(4)
% plot(t,y*1e3);
% set(gca,'FontSize',16)
% xlabel('Tempo [s]')
% ylabel('x_2: Posição da armadura [mm]')
% grid on