clear
%close all
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
%u25 = zeros(size(t)); % 16.34 a 17.05 mA é aplicado
%u25(t > TD) = 0.3;
u26 = 17.05*ones(size(t));% 17.05 a 16.34 mA
u26(t > TD) = 16.34;

%% Modelo metodo 1
K = (mean(y(end-20:end)) - mean(y(1:20)))/(u25(end) - u25(1))/24;
% tau_d = 26
% tau_d = 126; %s
% tau = 370 - tau_d;
tau_d = 120;
%tau = 414.8;
tau = 373.32; % 90% do anterior
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

%% Visualizacao

figure(1)
plot(t,y, 'c:')
xlabel('time[s]')
ylabel('Saida')
hold on
plot(t, yt, '-.', 'LineWidth',2)
plot(t, y_a, '--', 'color', '#0072BD', 'LineWidth',2)
hold off
title('Resultado para modelagem')
legend('Sinal original', 'Aproximacao (Metodo 1)', 'Aproximacao (Metodo das areas)', 'Location', 'southeast')
xlim([to tff])
grid on

figure(2)
plot(tv, yv, 'c:')
xlabel('time[s]')
ylabel('Saida')
hold on
plot(tv, ytv, '-.', 'LineWidth',2)
plot(tv, ytv_a, '--', 'color', '#0072BD', 'LineWidth',2)
hold off
title('Resultado para validacao')
legend('Sinal de validacao', 'Aproximacao (Metodo 1)', 'Aproximacao (Metodo das areas)', 'Location', 'southeast')
xlim([to tff2])
grid on
