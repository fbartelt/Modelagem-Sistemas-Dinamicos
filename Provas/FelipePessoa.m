%%Felipe Bartelt de Assis Pessoa - 2016026841
%Codigo conta com bastante reutilizaçao de scripts desenvolvidos nas
%Tarefas do semestre e scripts disponibilizados por Aguirre L.A.

clear; clc; close all
%numero 11
data = load('data_prova_tmsd_11.dat');

t = data(:, 1); %ts do enunciado =1 -> despreza-se esse vetor
u = data(:, 2);
y = data(:, 3);

ts = 1; %enunciado
t = 1:length(y);

%% Separaçao em amostras de Identificaçao e Validaçao
% Separa-se os dados em 70% de identificaçao e 30% de validaçao
N_id = ceil(0.7*length(y))

t_id = t(1:N_id);
u_id = u(1:N_id);
y_id = y(1:N_id);

t_val = t(N_id:end);
u_val = u(N_id:end);
y_val = y(N_id:end);

% Os dados de identificaçao e Validaçao sao mostrados na figura a seguir
%ENTRADA
figure(1);
subplot(211)
plot(t_id, u_id)
hold on 
plot(t_val, u_val)
hold off
xlabel('Tempo [s]')
ylabel('Amplitude')
title('Dados de entrada')
axis([min(t)-0.1 max(t)+0.1 min(u)-0.1 max(u)+0.1]);
legend('Identificaçao', 'Validaçao')

%SAIDA
subplot(212)
plot(t_id, y_id)
hold on 
plot(t_val, y_val)
hold off
xlabel('Tempo [s]')
ylabel('Amplitude')
title('Dados de saída')
axis([min(t)-0.1 max(t)+0.1 min(y)-0.1 max(y)+0.1]);
legend('Identificaçao', 'Validaçao')

%% Determinaçao de tempo morto
% Para determinar o tempo morto, pode-se utilizar a funçao de
% correlaçao cruzada. Toma-se como tempo morto como a quantidade de segundos
% iniciais ate que a autocorrelaçao tenha valor significativo
% Tomando-se a figura abaixo, pode-se notar que na amostra 1 o sistema
% esboça um pico de correlaçao significativo e, uma vez que esse grafico
% esta no dominio discreto, pode-se concluir que o sistema respondeu a um
% estimulo anterior. Portanto, tomar-se-a o atraso como 0.

figure(2)
[~,ry,~,~]=myccf([y_id u_id],100,0,1,'k'); %essa funçao originalmente tem nome de myccf2 nos scripts do Aguirre

theta_hat = 0

%% Determinaçao de ganho e constante de tempo
% Para determinar o ganho 'K' e constante de tempo 'tau', utiliza-se o
% metodo recursivo de minimos quadrados apresentado por Aguirre L.A.,
% utilizado durante a tarefa 4. Esse metodo eh baseado na aproximaçao da
% derivada temporal da saida. O script aqui utilizado foi adaptado do
% material disponibilizado por Aguirre L.A na plataforma ResearchGate.

% Primeiramente, despreza-se as 3 primeiras amostras, de forma a se 
% desprezar qualquer transiente inicial. Alem disso, remove-se os valores 
% medios do sinal de entrada e saida

init = 4;
y_ = y_id(init : end);
y_ = y_ - mean(y_);
t_ = t_id(init : end);
u_ = u_id(init - theta_hat : end - theta_hat);
u_ = u_ - mean(u_);

P=eye(2)*10^6;
teta(:,init-1)=[2; 2];
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

figure(3)
subplot(211);
plot(ts:ts:(b-2*init+2)*ts,tau(init-1:b-init),'k-');
ylabel('constante de tempo (s)');
title('(a)');
subplot(212)
plot(ts:ts:(b-2*init+2)*ts,ganho(init-1:b-init),'k-');
ylabel('ganho');
xlabel('tempo (s)');
title('(b)');

% Nota-se que a constante de tempo e ganho parecem convergir a valores
% proximos de 5.4 e 2.68, respectivamente. Assim, pode-se definir o sistema
% aproximado 'sys_hat'
K_hat = 2.68
tau_hat = 5.4
sys_hat = tf(K_hat, [tau_hat 1], 'InputDelay', theta_hat)

%% Validaçao
% Para validar o sistema aproximado, basta analisar a resposta do modelo
% obtido 'sys_hat' para os dados de entrada de validaçao 'u_val' e
% compara-los com as saidas de validaçao 'y_val'

y_hat = lsim(sys_hat, u_val, t_val);

figure(4)
plot(t_val, y_val)
hold on
plot(t_val, y_hat)
hold off
legend('Y_{val}', 'Y aproximado')

% Nota-se que o modelo obtido consegue aproximar satisfatoriamente o modelo
% real. O modelo tem resposta bastante similar aos valores reais a partir de
% 155 segundos, requisitando talvez um ajuste fino de ganho e constante de
% tempo para deixa-lo mais preciso. Especula-se que o fato do modelo nao 
% aproximar bem os valores iniciais eh devido a falta de condiçoes iniciais
% corretas na simulaçao, uma vez que, passado pouco mais de 10 segundos, o
% modelo passa a exibir uma resposta satisfatoria.