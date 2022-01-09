% Felipe Bartelt  de Assis Pessoa - 2016026841
% 
% Script baseado nos codigos fornecidos por Bruno Otávio Soares Teixeira e
% Luis Antonio Aguirre
%
clear
close all
addpath(genpath('./Tarefa2_DadosExperimentais'))
S = load('teste1.txt');
V = load('teste2.txt'); % para validacao

ts = S(:, 1);
tv = V(:, 1);
ys = S(:, 2);
yv = V(:, 2);

t0 = find(ts==0); %tempo de aplicacao do degrau
t1 = find(tv==0);
t = ts(t0:end);
tv = tv(t1:end);
y = ys(t0:end);
yv = yv(t1:end)

N_ciclos = 61 %52 61
zeta = 0.6/N_ciclos
yf = mean(y(end-20:end));
K = (yf - y(1)) % supondo degrau unitario
[~, l] = findpeaks(-y(247:500)); % para aproximar o periodo de um ciclo
tx = t(l); % tempos onde terminam alguns ciclos
T = mean(tx(2:end)-tx(1:end-1)) % periodo dos ciclos
wn = 2*pi/T

G = tf(K*wn^2, [1 2*zeta*wn wn^2])
yt = step(G, t) + y(1);


subplot(2,1,1)
plot(t, y, 'LineWidth',1.5)
hold on
plot(t, yt, 'LineWidth',1.5)
hold off
legend('Sinal original', 'Modelo 2°ordem','Location', 'southeast')
xlim([t(1) t(end)])
ylim([min(y)-0.05 max(y)+0.05])
title('Resposta da Aproximaçao')
ylabel('Saida')
xlabel('Tempo [s]')
subplot(2,1,2)
plot(tv, yv, 'LineWidth',1.5)
hold on
plot(t, yt, 'LineWidth',1.5)
hold off
xlim([tv(1) tv(end)])
ylim([min(yv)-0.05 max(yv)+0.05])
ylabel('Saida')
xlabel('Tempo [s]')
title('Resposta da Validacao')
legend('Sinal validacao', 'Modelo 2°ordem','Location', 'southeast')