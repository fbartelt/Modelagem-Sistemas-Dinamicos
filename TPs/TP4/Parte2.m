%Felipe Bartelt de Assis Pessoa - 2016026841
% Script adaptado dos algoritmos de Aguirre, L. A. (ermqfe.m)
close all
clear

%% A
K = 2;
tau = 5;
theta = 4;
u = prbs(400, 6, round(tau/3));
u = u - mean(u);
ts = 1;
t = 0:ts:399;
sys = tf(K, [tau 1], 'InputDelay', theta)

[y, ~] = lsim(sys, u, t);
figure(1)
lsim(sys, u, t)

%% B
lag = 20;
lu = round(length(y)/2);
figure(2)
[tx,ruy,l,B1]=myccf([y u'],lag,0,1,'b');
theta_hat = 4;

%% C
init = 28;
y_ = y(init:end);
y_ = y_ - mean(y_);
t_ = t(init:end);
u_ = u(init-theta_hat:end-theta_hat);
u_ = u_ - mean(u_);

P=eye(2)*10^6;
teta(:,init-1)=[1; 0.86];
lambda=0.99;

% Algoritmo recursivo
for k=init:length(y_)
   psi_k=[y_(k-1);u_(k-1)];
   K_k = (P*psi_k)/(psi_k'*P*psi_k+lambda);
   teta(:,k)=teta(:,k-1)+K_k*(y_(k)-psi_k'*teta(:,k-1));
   P=(P-((P * (psi_k * psi_k')*P)/(psi_k'*P*psi_k+lambda)))/lambda;
end

[a b]=size(teta);

% constante de tempo e ganho
for k=1:b   
   tau(k)=-ts/(teta(1,k)-1);
   ganho(k)=tau(k)*teta(2,k)/ts;
end

figure(4)
subplot(211);
plot(ts:ts:(b-2*init+2)*ts,tau(init-1:b-init),'k-');
ylabel('constante de tempo (s)');
title('(a)');
subplot(212)
plot(ts:ts:(b-2*init+2)*ts,ganho(init-1:b-init),'k-');
ylabel('ganho');
xlabel('tempo (s)');
title('(b)');

sys_hat = tf(2, [5.5 1], 'InputDelay', theta_hat)

figure(5);
subplot(211)
step(sys, 'r-');
hold on;
step(sys_hat, 'b-.');
hold off
legend('Real', 'Estimado')
ylim([-0.2 2.25])
xlim([-0.2 40])
subplot(212)
impulse(sys, 'r-');
hold on;
impulse(sys_hat, 'b-.');
legend('Real', 'Estimado')
hold off
xlim([-0.2 40])
ylim([-0.02 0.4])