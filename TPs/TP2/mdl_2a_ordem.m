% close all; 
clear all; clc; close all

% Sistema de segunda ordem com atraso puro de tempo:
tau1 = 3; tau2 = 2;
theta = 3;
K = 0.75;
G1 = tf(K, [tau1*tau2 tau1+tau2 1], 'ioDelay', theta);

Ts = 0.01;
u0 = 0; y0 = 0;
t = [0:Ts:40-Ts]';
u = [0; ones(length(t)-1,1)] + u0;
y = lsim(G1, u, t) + y0;
y = y + 1*0.02*randn(length(y),1);

% figure(1); 
% subplot(211); stairs(t, u,'LineWidth',2); xlabel('t (s)'); ylabel('u(t)');
% subplot(212); plot(t, y); xlabel('t (s)'); ylabel('y(t)');
figure(3); stairs(t, u,'LineWidth',2); xlabel('t (s)'); ylabel('u(t)');
figure(1); plot(t, y); xlabel('t (s)'); ylabel('y(t)');

%% Metodo das areas (aproximacao de primeira ordem):
K1 = (mean(y(end-20:end)) - mean(y(1:20)))/(u(end) - u(1));
yn = y./K1;
area = sum(Ts*(u - yn));

tau1a = exp([1])*sum(Ts*yn(1:find(t==round(area))));
theta1a = area - tau1a;

G1a = tf(K1, [tau1a 1], 'ioDelay', theta1a);
y1 = lsim(G1a, u, t) + y0;

% figure(1); subplot(212); hold on; plot(t, y1, 'r-.','LineWidth',2); xlabel('t (s)'); ylabel('y(t)');
figure(1);  hold on; plot(t, y1, 'r-.','LineWidth',2); xlabel('t (s)'); ylabel('y(t)');


%% Resposta complementar:
yy = log(abs(1 - y./(K1*u)));
% figure(2); subplot(211); plot(t, yy); xlabel('t (s)'); ylabel('ln(1 - y(t)/(K*u(t)))');
figure(2); plot(t, yy); xlabel('t (s)'); ylabel('w(t)');
coef = polyfit(t(900:1600),yy(900:1600),1);
hold on; plot(t, coef(1)*t + coef(2), 'm-.', 'LineWidth', 2); axis([0 18 -8 2]);
tau2a = -1/coef(1);

yy2 = log(abs(exp(coef(2))*exp(-(t)./tau2a) - (1 - y./(K1*u))));
coef2 = polyfit(t(10:700),yy2(10:700),1);
% figure(2); subplot(212); plot(t, yy2); xlabel('t (s)'); ylabel('ln(\tau_1/(\tau_1 - \tau_2)e^{t/\tau_1} - (1 - y(t)/(K*u(t))))');
figure(5);  plot(t, yy2); xlabel('t (s)'); ylabel('v(t)');
hold on; plot(t, coef2(1)*t + coef2(2), 'm-.', 'LineWidth', 2); axis([0 18 -8 2]);
tau2b = -1/coef2(1);

G2a = tf(K1, [tau2a*tau2b  tau2a+tau2b  1], 'ioDelay', theta);
y2 = lsim(G2a, u, t) + y0;

% figure(1); subplot(212); hold on; plot(t, y2, 'm--','LineWidth',2); xlabel('t (s)'); ylabel('y(t)');
figure(1); hold on; plot(t, y2, 'm--','LineWidth',2); xlabel('t (s)'); ylabel('y(t)');

