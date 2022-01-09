% Felipe Bartelt  de Assis Pessoa - 2016026841
% 
% Script baseado nos codigos fornecidos por Bruno Otávio Soares Teixeira e
% Luis Antonio Aguirre
%

clear; close all
t0 = 0; % Tempo inicial
tff = 10; % Tempo final
h = 0.0001; % Intervalo de integracao
t = t0:h:tff; % Vetor de tempo

%% Atuador Solenoide
x0 = [0;0;0]; %cond iniciais
% x1 - corrente, x2 - posicao armadura, x3 - velocidade armadura
x = [x0 zeros(length(x0),length(t)-1)];
e_in = zeros(length(t), 1);
w = (pi/2 + 2*pi)/0.5; % 5pi rad (permite 2 ciclos ate t=0.5s)
e_in(t>=0.5) = 10;  % degrau em t=0.5, amplitude=10
e_in(t>=5.5) = 15;
u = e_in;

for k=2:length(t)
    x(:,k)=rk4(x(:,k-1),u(k),u(k),h,t(k));
end
%% Pontos de operacao
%u1,t1,y -> pontode operacao 0
%u2,t2,y2 -> pontode operacao = estacionario do anterior
change = find(t==5.5)
change2 = find(t==5)
y = x(2,1:change-1);
y2 = x(2,change2:end);
u1 = u(1:change-1);
u2 = u(change2:end);
t1 = t(1:change-1);
t2 = t(change2:end);

%% Modelo 1° ordem
K = (y(end) - y(1))/(u1(end) - u1(1))
tau_d = 0;
tau = 0.5165 - 0.5;
sys = tf(K, [tau 1])
[yt1, ~] = lsim(sys, u1, t1);
[yt2, ~] = lsim(sys, u2, t2);

%% Modelo 2° ordem
step_time = find(t==0.5)
K2 = (y(end) - y(1))/(u1(end) - u1(1))
ya = y(step_time:end);
ua = u1(step_time:end)./max(u1);
ta = t1(step_time:end);

%K2 = (y(end) - y(1))/(ua(end) - ua(1))
w = log(abs(1 - ya./(K2*ua.')));
coef = polyfit(ta(20000:end),w(20000:end),1);
tau1a = -1/coef(1)*-1e-18; %aprox 0
%tau1a = 1e64;

v = log(abs(exp(coef(2))*exp(-(ta)./tau1a) - (1 - ya./(K2*ua.'))));
coef2 = polyfit(ta(1:60),v(1:60),1);

tau2c = -1/coef2(1)*0.665e-4;

%K2 = K2/3e303
G2a = tf(K2, [tau1a*tau2c  tau1a+tau2c  1]);
y22 = lsim(G2a, u1, t1);
y22v = lsim(G2a, u2, t2);

figure(3)
subplot(2,1,1)
plot(ta, w)
xlabel('t [s]'); 
ylabel('w(t)');
xlim([ta(1) ta(end)]);
hold on
plot(ta, coef(1)*ta + coef(2), 'm-.', 'LineWidth', 2)
hold off
title('Passo a Passo (Resposta Complementar)')
ylim([w(1)-0.1 w(end)+0.1])
subplot(2,1,2)
plot(ta, v, 'LineWidth', 2)
hold on
xlabel('t [s]')
ylabel('v(t)')
plot(ta, coef2(1)*ta + coef2(2), 'm-.', 'LineWidth', 2)
hold off

%% Visualizacao
% Primeira ordem
colororder({'#0072BD','#D95319'})
figure(1)
subplot(2,1,1)
yyaxis left
plot(t1,y,'-','color','m', 'LineWidth', 3);
hold on 
plot(t1,yt1,'-.','color','#0072BD', 'LineWidth', 2);
hold off
ylabel('x_2: Posição da armadura [m]')
yyaxis right
plot(t1,u1, 'LineWidth', 1.5, 'color', '#D95319');
ylabel('Amplitude do Degrau')
ylim([0 11])
xlim([t1(1) t1(end)])
xlabel('Tempo [s]')
legend('Sinal original', 'Resposta do modelo', 'Degrau', 'location', 'se')
title('1° ordem - Primeiro Ponto de Operacao')
grid on

subplot(2,1,2)
yyaxis left
plot(t2,y2,'-','color','m', 'LineWidth', 3);
hold on 
plot(t2,yt2,'-.','color','#0072BD', 'LineWidth', 2);
hold off
ylabel('x_2: Posição da armadura [m]')
yyaxis right
plot(t2,u2, 'LineWidth', 1.5, 'color', '#D95319');
ylabel('Amplitude do Degrau')
ylim([0 u2(end)+1])
xlim([t2(1) t2(end)])
xlabel('Tempo [s]')
legend('Sinal original', 'Resposta do modelo', 'Degrau', 'location', 'se')
title('Segundo Ponto de Operacao')
grid on


% Segunda ordem
colororder({'#0072BD','#77AC30'})
figure(2)
subplot(2,1,1)
yyaxis left
plot(t1,y,'-','color','m', 'LineWidth', 3);
hold on 
plot(t1,y22,'-.','color','#0072BD', 'LineWidth', 2);
hold off
ylabel('x_2: Posição da armadura [m]')
ylim([-0.1e-3 6e-3])
yyaxis right
plot(t1,u1, 'LineWidth', 1.5);
ylabel('Amplitude do Degrau')
ylim([0 11])
xlim([t1(1) t1(end)])
xlabel('Tempo [s]')
legend('Sinal original', 'Resposta do modelo', 'Degrau', 'location', 'se')
title('2° ordem - Primeiro Ponto de Operacao')
grid on

subplot(2,1,2)
yyaxis left
plot(t2,y2,'-','color','m', 'LineWidth', 3);
hold on 
plot(t2,y22v,'-.','color','#0072BD', 'LineWidth', 2);
hold off
ylim([-0.1e-3 2e-2])
ylabel('x_2: Posição da armadura [m]')
yyaxis right
plot(t2,u2, 'LineWidth', 1.5);
ylabel('Amplitude do Degrau')
ylim([0 u2(end)+1])
xlim([t2(1) t2(end)])
xlabel('Tempo [s]')
legend('Sinal original', 'Resposta do modelo', 'Degrau', 'location', 'se')
title('Segundo Ponto de Operacao')
grid on