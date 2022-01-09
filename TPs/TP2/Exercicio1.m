% Felipe Bartelt  de Assis Pessoa - 2016026841
% 
% Script baseado nos codigos fornecidos por Bruno Otávio Soares Teixeira e
% Luis Antonio Aguirre
%
clear
close all
addpath(genpath('./Tarefa2_DadosExperimentais'))
S = load('ENS_25.DAT');
V = load('ENS_26.DAT'); % para validacao

t = S(:, 1);
y = S(:, 2);

tv = V(:, 1); % validation data
yv = V(:, 2); % validation data

%% Parametros
to =0;		% tempo inicial em segundos
tff = t(end); 	% tempo de simulacao em segundos
tff2 = tv(end);% tempo de simulacao em segundos

TD = 100;	% tempo em que ocorre o degrau
u25 = 16.34*ones(size(t)); % 16.34 a 17.05 mA é aplicado
u25(t > TD) = 17.05;
u26 = 17.05*ones(size(t));% 17.05 a 16.34 mA
u26(t > TD) = 16.34;

%% Modelo metodo 1
K = (mean(y(end-20:end)) - mean(y(1:20)))/(u25(end) - u25(1))/24;
% tau_d = 126; %s
% tau = 370 - tau_d;
tau_d = 120;
%tau = 414.8;
tau = 370; % 90% do anterior
%tau = 244
sys = tf(K, [tau 1])
s = tf('s');
delay = pade(exp(-tau_d * s), 5);
sys = sys * delay
[yt, ~] = lsim(sys, u25, t); % sinal resultante da aproximacao
offset = y(1);
yt = yt + offset;
[ytv, ~] = lsim(sys, u26, tv); % sinal aproximado pelo modelo na validacao
ytv = ytv + offset;

%% Modelo metodo das areas
Ts = mean(t(2:end) - t(1:end-1))
y0 = mean(y(1:20));
yff = mean(y(end-20:end));
K_a = (yff - y0)/(u25(end) - u25(1));
fator = 1/yff;

% normalizacao dos dados de entrada e saida
yn = y - y0;
yn = yn./fator;
un = u25 - u25(1);
un = un./fator;

area = sum(Ts*(un - yn));
K_a = K_a*0.042;

tau_a = exp(1)*sum(Ts*yn(1:find(t==t(round(area/Ts + 1)))));
tau_a = tau_a*19
theta = area - tau_a;

sys_a = tf(K_a, [tau_a 1]);
delay_a = pade(exp(-theta * s), 5);
sys_a = sys_a * delay_a

[y_a, ~] = lsim(sys_a, u25, t);
y_a = y_a + y0;

[ytv_a, ~] = lsim(sys_a, u26, tv); % sinal aproximado pelo modelo na validacao
ytv_a = ytv_a + y0;

%% Modelo 2°ordem (Resposta Complementar)
y0 = mean(y(1:20));
yff = mean(y(end-20:end));
K2 = (yff - y0)/(u25(end) - u25(1));
%find(abs(t-theta) <= 0.8) % 116
ya = y(116:end);
ua = u25(116:end)./max(u25);
%K2 = K2 * max(u25)
ta = t(1:length(t) - theta/Ts);
w = log(abs(1 - ya./(K2*ua)));
coef = polyfit(ta(1500:3000),w(1500:3000),1);

tau1a = -1/coef(1)*0.005;

v = log(abs(exp(coef(2))*exp(-(ta)./tau1a) - (1 - ya./(K2*ua))));
coef2 = polyfit(ta(1:900),v(1:900),1);

tau2c = -1/coef2(1)*0.7;

K2 = K2/24
G2a = tf(K2, [tau1a*tau2c  tau1a+tau2c  1], 'ioDelay', theta);
y2 = lsim(G2a, u25, t) + y0;
y2v = lsim(G2a, u26, tv) + y0;

figure(1)
subplot(2,1,1)
plot(ta, w)
xlabel('t [s]'); 
ylabel('w(t)');
xlim([to tff]);
hold on
plot(ta, coef(1)*ta + coef(2), 'm-.', 'LineWidth', 2)
hold off
subplot(2,1,2)
plot(ta, v)
xlabel('t [s]')
ylabel('v(t)')
title('Passo a Passo (Resposta Complementar)')
hold on
plot(ta, coef2(1)*ta + coef2(2), 'm-.', 'LineWidth', 2)
hold off


%% Visualizacao
figure(3)
colororder({'k','#77AC30'})
subplot(2,1,1)
yyaxis left
plot(t,y, 'c:')
xlabel('time[s]')
ylabel('Saida')
hold on
plot(t, yt, '-','color','#D95319', 'LineWidth',2)
plot(t, y_a, '--', 'color', '#0072BD', 'LineWidth',2)
plot(t, y2, ':', 'color', '#EDB120', 'LineWidth',3)
hold off
title('Resultado para modelagem')

xlim([to tff])
ylim([-0.05, 0.4])
yyaxis right
plot(t, u25, 'color', '#77AC30')
legend('Sinal original', '1°ordem (Metodo 1)', '1°ordem (Metodo das areas)', '2°ordem (Resposta Complementar)', 'Degrau','Location', 'southeast')
xlim([to tff])
ylabel('Entrada')
grid on

subplot(2,1,2)
yyaxis left
plot(tv, yv, 'c:')
xlabel('time[s]')
ylabel('Saida')
hold on
plot(tv, ytv, '-','color','#D95319', 'LineWidth',2)
plot(tv, ytv_a, '--', 'color', '#0072BD', 'LineWidth',2)
plot(t, y2v, ':','color', '#EDB120', 'LineWidth',3)
hold off
title('Resultado para validacao')

xlim([to tff2])
ylim([-0.2, 0.38])
yyaxis right
plot(tv, u26, 'color', '#77AC30')
legend('Sinal validacao', '1°ordem (Metodo 1)', '1°ordem (Metodo das areas)', '2°ordem (Resposta Complementar)', 'Degrau','Location', 'southeast')
xlim([to tff])
ylabel('Entrada')
grid on
