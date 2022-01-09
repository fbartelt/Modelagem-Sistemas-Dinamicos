%close all; 
clear all; clc

% Sistema de primeira ordem com atraso puro de tempo:
tau = 4;
theta = 3;
K = 0.75;
G1 = tf(K, [tau 1], 'ioDelay', theta);

Ts = 0.01;
u0 = 0; y0 = 0;
t = [0:Ts:30-Ts]';
u = [0; ones(length(t)-1,1)] + u0;
y = lsim(G1, u, t) + y0;
y = y + 0.02*randn(length(y),1);

% figure(1); 
% subplot(211); stairs(t, u,'LineWidth',2); xlabel('t (s)'); ylabel('u(t)');
% subplot(212); plot(t, y); xlabel('t (s)'); ylabel('y(t)');

figure(3); stairs(t, u, 'LineWidth',2); xlabel('t (s)'); ylabel('u(t)');
figure(1); plot(t, y); xlabel('t (s)'); ylabel('y(t)');

figure(4); stairs(t, u,'k-.','LineWidth',2); xlabel('t (s)'); ylabel('u(t), y(t)');
hold on; plot(t, y); 

% Metodo das areas:
K1 = (mean(y(end-20:end)) - mean(y(1:20)))/(u(end) - u(1));
yn = y - y0;
yn = yn./K1;
area = sum(Ts*(u - yn));

tau1 = exp(1)*sum(Ts*yn(1:find(t==round(area))));
theta1 = area - tau1;

G1a = tf(K1, [tau1 1], 'ioDelay', theta1+0.3);
y1 = lsim(G1a, u, t) + y0;

% figure(1); subplot(212); hold on; plot(t, y1, 'r-.','LineWidth',2); xlabel('t (s)'); ylabel('y(t)');
figure(1); hold on; plot(t, y1, 'r-.','LineWidth',2); xlabel('t (s)'); ylabel('y(t)');

% Resposta complementar:
yy = log(1 - y./(K1*u));
figure(2); plot(t, yy); xlabel('t (s)'); ylabel('ln(y(t)/(K*u(t)))');
coef = polyfit(t(300:1000),yy(300:1000),1);
hold on; plot(t, coef(1)*t + coef(2), 'm-.', 'LineWidth', 2);

tau2 = -1/coef(1);
G2a = tf(K1, [tau2 1], 'ioDelay', theta1);
y2 = lsim(G2a, u, t) + y0;

% figure(1); subplot(212); hold on; plot(t, y2, 'm--','LineWidth',2); xlabel('t (s)'); ylabel('y(t)');
figure(1);  hold on; plot(t, y2, 'm--','LineWidth',2); xlabel('t (s)'); ylabel('y(t)');


