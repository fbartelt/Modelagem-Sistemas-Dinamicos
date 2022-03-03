%Felipe Bartelt de Assis Pessoa - 2016026841
% Script adaptado dos algoritmos de Aguirre, L. A. 
clear
clc
close all
%% Separacao dos dados
data=load('tmsd_tarefa5.txt');

t = data(:, 1);
u = data(:, 2);
u = u - mean(u);

y = data(:, 3);
y = y - mean(y);
[N, m] = size(data);
N_train = ceil(0.7*N);

%Construçao
t_train = t(1:N_train);
u_train = u(1:N_train);
y_train = y(1:N_train);

%Validacao
t_val = t(N_train:end);
u_val = u(N_train:end);
y_val = y(N_train:end);


%ENTRADA
figure(1);
subplot(2, 1, 1)
plot(t_train, u_train)
hold on 
plot(t_val, u_val)
hold off
xlabel('Tempo [s]')
ylabel('Amplitude')
title('Dados de entrada')
axis([0 200 -1 1]);
legend('Identificaçao', 'Validaçao')

%SAIDA
subplot(2, 1, 2)
plot(t_train, y_train)
hold on 
plot(t_val, y_val)
hold off
xlabel('Tempo [s]')
ylabel('Amplitude')
title('Dados de saída')
axis([0 200 -1 1]);
legend('Identificaçao', 'Validaçao')

%% Tempo de amostragem 

figure(2)
subplot(211)
[~,ry,~,~]=myccf([y_train y_train],300,0,1,'k');
title('Funçao de autocorrelaçao da saida')

subplot(212)
y_train_2 = y_train.^2 - mean(y_train.^2);
[~,ry2,~,~]=myccf([y_train_2 y_train_2],300,0,1,'k');
title('Funçao de autocorrelaçao da saida ao quadrado')

tau_m1 = min(find(ry==min(ry)), find(ry2==min(ry2)));
% tau_m de 33. Decimacao de 2 --> 10 <= tau_m' <= 20

y_train_dec = y_train(1:2:end);
u_train_dec = u_train(1:2:end);
t_train_dec = 1:length(y_train_dec);

figure(3)
subplot(211)
[~,ry,~,~]=myccf([y_train_dec y_train_dec],300,0,1,'k');
title('Funçao de autocorrelaçao da saida decimada')
subplot(212)
y_train_dec2 = y_train_dec.^2 - mean(y_train_dec.^2);
[~,ry2,~,~]=myccf([y_train_dec2 y_train_dec2],300,0,1,'k');
title('Funçao de autocorrelaçao da saida decimada ao quadrado')

tau_m = min(find(ry==min(ry)), find(ry2==min(ry2)));

figure(4)
subplot(211)
plot(t_train_dec, u_train_dec)
title('Sinal de entrada decimado')
subplot(212)
plot(t_train_dec, y_train_dec)
title('Sinal de saida decimado')
xlabel('Amostras')
ylabel('Amplitude')

% Correlaçao entrada/saida

figure(5)
[~,ry,~,~]=myccf([y_train_dec u_train_dec],300,0,1,'k');
title('Funçao de correlaçao entre entrada e saida decimadas')

%% Selecao de estrutura
% Criterio de Akaike
degree = 6;
num_delay = degree / 2;
Y = y_train_dec(num_delay+1:end);
idx_aic = 1;

psi = [];
AIC = zeros(degree, 1);

for delay = 1:num_delay
    
    psi = [psi y_train_dec(num_delay-delay+1:end-delay)];
    theta_hat = pinv(psi)*Y;
    residual = Y-psi*theta_hat;

    AIC(idx_aic) = length(residual)*log(var(residual))+2*length(theta_hat);
    
    psi = [psi u_train_dec(num_delay-delay+1:end-delay)];
    theta_hat = pinv(psi)*Y;
    residual = Y-psi*theta_hat;

    AIC(idx_aic+1) = length(residual)*log(var(residual))+2*length(theta_hat);
    
    idx_aic = idx_aic + 2;
end

figure(6);
plot(1:length(AIC), AIC, '-o');
xlabel('Número de regressores');
ylabel('AIC');
xlim([0 degree+0.2]);


%Agrupamento de termos
idx_sigma = 1;

psi = [];
sigma = zeros(degree, 1);

yt=zeros(length(y_train_dec), 1);
for k=2:length(y_train_dec)
    yt(k)=yt(k-1)+y_train_dec(k);
end

Y = yt(num_delay+1:end);

for delay = 1:num_delay
    
    psi = [psi yt(num_delay-delay+1:end-delay)];
    theta_hat = pinv(psi)*Y;
    residual = Y-psi*theta_hat;

    sigma(idx_sigma) = sum(theta_hat(1:2:end));
    
    psi = [psi u_train_dec(num_delay-delay+1:end-delay)];
    theta_hat = pinv(psi)*Y;
    residual = Y-psi*theta_hat;

    sigma(idx_sigma+1) = sum(theta_hat(1:2:end));
    
    idx_sigma = idx_sigma + 2;
end

figure(7);
plot(1:length(sigma), abs(sigma-1), '-o');
xlabel('Número de regressores + 1');
ylabel('|\Sigma_y - 1|');
xlim([0 degree+1]);

%% Estimacao de parâmetros

for k=1:length(y_train_dec)-3
   PSI(k,:) = [flip(y_train_dec(k:k+2).') flip(u_train_dec(k:k+1).')];
end

Y = y_train_dec(4:end);
Theta = pinv(PSI) * Y
%res = Y-PSI*Theta;
% figure(8)
% subplot(211)
% [~,r_rr,~,~]=myccf([res res],1000,0,1,'k');
% axis([-50 1050 -1.1 1.1]);
% title('Função de autocorrelação do resíduos');
% subplot(212)
% [~,r_ur,~,~]=myccf([u_train_dec(6:705) res],1000,0,1,'k');
% axis([-50 1050 -1.1 1.1]);
% title('Função de correlação cruzada entre a entrada e os resíduos');



%% Validacao
% Simulacao um passo a frente
PSI = [y_val(3:end-1) y_val(2:end-2) y_val(1:end-3) u_val(2:end-2) u_val(1:end-3)];
y_hat = PSI*Theta;
res = y_val(4:end)-y_hat;
RMSE = sqrt((y_val(4:end) - y_hat).' * (y_val(4:end) - y_hat)/length(y_hat))

figure(9)
plot(t_val, y_val, '.', 'LineWidth', 2.5);
hold on;
plot(t_val(4:end), y_hat, '--', 'LineWidth', 2.5);
axis([t_val(1) t_val(end)+0.1 min(y_val)-0.2 max(y_val)+0.2]);
legend('Sinal de Validaçao', 'Sinal aproximado')
title('Simulaçao um passo a frente')
hold off

figure(10)
subplot(211)
[~,r_rr,~,~]=myccf([res res],500,0,1,'k');
axis([-0.1 501 -1.1 1.1]);
title('Função de autocorrelação do resíduos');
subplot(212)
[~,r_ur,~,~]=myccf([u_val(4:end) res],500,0,1,'k');
axis([-0.1 501 -1.1 1.1]);
title('Função de correlação cruzada entre a entrada e os resíduos');

% Simulacao livre
delay = 3;

y_hat = zeros(length(y_val)-delay, 1);
y_hat(1) = [y_val(delay) y_val(delay-1) y_val(delay-2) u_val(delay)  u_val(delay-1)]*Theta;
y_hat(2) = [y_hat(1) y_val(delay) y_val(delay-1) u_val(delay+1)  u_val(delay)]*Theta;
y_hat(3) = [y_hat(2) y_hat(1) y_val(delay) u_val(delay+2)  u_val(delay+1)]*Theta;
y_hat(4) = [y_hat(3) y_hat(2) y_hat(1) u_val(delay+3) u_val(delay+2)]*Theta;

for i = 5:length(y_val)-delay
    y_hat(i) = [y_hat(i-1) y_hat(i-2) y_hat(i-3) u_val(delay+i-1) u_val(delay+i-2)]*Theta;
end
res = y_val(4:end) - y_hat;
RMSE = sqrt((y_val(4:end) - y_hat).' * (y_val(4:end) - y_hat)/length(y_hat))

figure(12);
plot(t_val, y_val, '.');
hold on;
plot(t_val(4:end), y_hat(1:end), '--');
axis([t_val(1) t_val(end)+0.1 min(y_val)-0.2 max(y_val)+0.2]);
hold off
legend('Sinal de Validaçao', 'Sinal aproximado')
title('Simulaçao livre')

figure(13)
subplot(211)
[~,r_rr,~,~]=myccf([res res],500,0,1,'k');
axis([-0.1 501 -1.1 1.1]);
title('Função de autocorrelação do resíduos');
subplot(212)
[~,r_ur,~,~]=myccf([u_val(4:end) res],500,0,1,'k');
axis([-0.1 501 -1.1 1.1]);
title('Função de correlação cruzada entre a entrada e os resíduos');