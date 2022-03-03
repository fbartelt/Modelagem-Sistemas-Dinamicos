%Felipe Bartelt de Assis Pessoa - 2016026841
% Script adaptado dos algoritmos de Aguirre, L. A.
format long;clear;clc;close all;

%% Projeto de teste
Hs = tf([0.1 1], [0.5 0.3 0.1]);
H = c2d(Hs,1,'zoh');
theta_real = [H.Denominator{:}(2) H.Denominator{:}(3) H.Numerator{:}(2) H.Numerator{:}(3)];

% figure(1);step(Hs); hold on; step(H); hold off; legend('H(s)', 'H(z)')
% tau = 3.71

T_b = 1;
t = 0:700;
u = prbs(length(t),9,T_b);
u = u - mean(u);
y = lsim(H, u, t);
sigma_e = std(y)/10;
noise = sigma_e * randn(size(y));
ym = y + noise;
[y_step, t_step] = step(H);
u_step = ones(length(y_step), 1);

figure(1)
subplot(211);
plot(t, u);
xlabel('Tempo [s]')
ylabel('Amplitude')
title('Entrada')
ylim([min(u)-0.1 max(u)+0.1]);

subplot(212);
plot(t, y, '-', 'LineWidth', 2);
hold on 
plot(t, ym, '-.', 'LineWidth', 2);
legend('y', 'y_m')
xlabel('Tempo [s]')
ylabel('Amplitude')
title('Saída')
%ylim([-0.2 10]);

%% Estimação de parâmetros

for k=1:length(y)-2
   psi(k,:) = [flip(-y(k:k+1).') flip(u(k:k+1))];
   psim(k,:) = [flip(-ym(k:k+1).') flip(u(k:k+1))];
end

theta_real
theta = pinv(psi) * y(3:end)
thetam = pinv(psim) * ym(3:end)

psi_t = [-y_step(2:end-1) -y_step(1:end-2) u_step(2:end-1) u_step(1:end-2)];
yhat = psi_t * theta;
yhatm = psi_t * thetam;

figure(2)
plot(t_step, y_step, '-', 'LineWidth', 2)
hold on
plot(t_step(3:end), yhat, '--', 'LineWidth', 2)
plot(t_step(3:end), yhatm, '-.', 'LineWidth', 2)
hold off
legend('Resposta ao degrau', 'Resposta estimada sem ruido(2° ordem)', 'Resposta estimada com ruido(2° ordem)', 'Location', 'SE')
xlabel('Tempo[s]')
ylabel('Amplitude')
axis([-0.05 24 0 11])

%% Estrutura do modelo

%1ordem
for k=1:length(y)-1
   psi1(k,:) = [flip(-y(k).') flip(u(k))];
   psi1m(k,:) = [flip(-ym(k).') flip(u(k))];
end

theta_real
theta1 = pinv(psi1) * y(2:end)
theta1m = pinv(psi1m) * ym(2:end)

psi_t1 = [-y_step(1:end-1) u_step(1:end-1)];
yhat1 = psi_t1 * theta1;
yhatm1 = psi_t1 * theta1m;

figure(3)
plot(t_step, y_step, '-', 'LineWidth', 2)
hold on
plot(t_step(2:end), yhat1, '--', 'LineWidth', 2)
plot(t_step(2:end), yhatm1, '-.', 'LineWidth', 2)
hold off
legend('Resposta ao degrau', 'Resposta estimada sem ruido(1° ordem)', 'Resposta estimada com ruido(1° ordem)', 'Location', 'SE')
xlabel('Tempo[s]')
ylabel('Amplitude')
axis([-0.05 24 0 11])

%3ordem
for k=1:length(y)-3
   psi3(k,:) = [flip(-y(k:k+2).') flip(u(k:k+2))];
   psi3m(k,:) = [flip(-ym(k:k+2).') flip(u(k:k+2))];
end

theta_real
theta3 = pinv(psi3) * y(4:end)
theta3m = pinv(psi3m) * ym(4:end)

psi_t3 = [-y_step(3:end-1) -y_step(2:end-2) -y_step(1:end-3) u_step(3:end-1) u_step(2:end-2) u_step(1:end-3)];
yhat3 = psi_t3 * theta3;
yhatm3 = psi_t3 * theta3m;

figure(4)
plot(t_step, y_step, '-', 'LineWidth', 2)
hold on
plot(t_step(4:end), yhat3, '--', 'LineWidth', 2)
plot(t_step(4:end), yhatm3, '-.', 'LineWidth', 2)
hold off
legend('Resposta ao degrau', 'Resposta estimada sem ruido(3° ordem)', 'Resposta estimada com ruido(3° ordem)', 'Location', 'SE')
xlabel('Tempo[s]')
ylabel('Amplitude')
axis([-0.05 24 0 11])
