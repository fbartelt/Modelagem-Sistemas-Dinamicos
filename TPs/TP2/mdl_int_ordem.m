close all; clear all; clc

% Sistema de primeira ordem com atraso puro de tempo:
theta = 3;
K = 0.75;
G1 = tf(K, [1 0], 'ioDelay', theta);

Ts = 0.01;
u0 = 0; y0 = 0;
t = [0:Ts:30-Ts]';
u = [0; ones(length(t)/3,1); zeros(2*length(t)/3-1,1)] + u0;
y = lsim(G1, u, t) + y0;
y = y + 0.1*randn(length(y),1);

% figure(1); 
% subplot(211); stairs(t, u,'LineWidth',2); xlabel('t (s)'); ylabel('u(t)');
% subplot(212); plot(t, y); xlabel('t (s)'); ylabel('y(t)');
figure(2); stairs(t, u, 'k-.', 'LineWidth',2); 
hold on; plot(t, y); xlabel('t (s)'); ylabel('u(t), y(t)');
figure(1); plot(t, y); xlabel('t (s)'); ylabel('y(t)');

% Metodo das areas:
theta1 = 3.2;
tw = 10;
deltaU = u(500) - u(1);
deltaY = mean(y(end-300:end)) - mean(y(1:300));
K1 = (deltaY)/(tw*(deltaU));

G1a = tf(K1, [1 0], 'ioDelay', theta1);
y1 = lsim(G1a, u, t) + y0;

% figure(1); subplot(212); hold on; plot(t, y1, 'r-.','LineWidth',2); xlabel('t (s)'); ylabel('y(t)');
figure(1); hold on; plot(t, y1, 'g-.','LineWidth',2); xlabel('t (s)'); ylabel('y(t)');

