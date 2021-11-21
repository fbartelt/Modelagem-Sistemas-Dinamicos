%% Simulacao do atuador solenoide
% Felipe Bartelt de Assis Pessoa - 2016026841
%
% Codigo principal e dependencias adaptados dos codigos de Luis Antonio 
% Aguirre (https://www.researchgate.net/project/Scripts-on-Nonlinear-Dynamics)
%
% O script nao e automatico, de forma que se deve descomentar diversas
% linhas para que as variaveis de entrada e os plots coincidam com o
% objetivo.

clear; close all

t0 = 0; % Tempo inicial
tf = 8; % Tempo final
h = 0.0001; % Intervalo de integracao
t = t0:h:tf; % Vetor de tempo

%% Condicoes iniciais
x0 = [0;0;0];

%% Inicializacoes
x = [x0 zeros(length(x0),length(t)-1)];
e_in = zeros(length(t), 1);
w = (pi/2 + 2*pi)/0.5; % 5pi rad (permite 2 ciclos ate t=0.5s)

%% Entradas
%e_in(t>=0.5) = 10;                   % degrau em t=0.5, amplitude=10
%e_in = e_in + 0.3*randn(size(e_in)); % degrau com ruido

%e_in(5000:5010) = 10;                % pulso de duracao 1ms

%e_in = 10*sin(w*t);                  % sinal senoidal

%% Superposicao
% x1 = x;
% x2 = x;
% e_in1 = e_in;
% e_in2 = e_in;
% e_in(t>=0.5) = 10;                   % degrau em t=0.5, amplitude=10
% e_in1(t>=0.5) = 7;                   % degrau em t=0.5, amplitude=7
% e_in2(t>=0.5) = 3;                   % degrau em t=0.5, amplitude=3
% u1 = e_in1;
% u2 = e_in2;

%% Diversos pontos de operacao
niveis = 0:1:15;
largura = round(length(t)/length(niveis));
aux = [];
for j=1:length(niveis)-1
   aux = [aux niveis(j)*ones(1, largura)];
end
aux = [aux niveis(j+1)*ones(1, length(t)-length(aux))];
e_in = aux;


u = e_in;
for k=2:length(t)
    x(:,k)=rk4(x(:,k-1),u(k),u(k),h,t(k));
    %x1(:,k)=rk4(x1(:,k-1),u1(k),u1(k),h,t(k)); % para superposicao
    %x2(:,k)=rk4(x2(:,k-1),u2(k),u2(k),h,t(k)); % para superposicao
end

%% Plots para respostas de cada estado
% figure(1)
% plot(t,x(1,:));
% set(gca,'FontSize',16)
% xlabel('Tempo [s]')
% ylabel('x_1: Corrente [A]')
% grid on
% 
% figure(2)
% plot(t,x(2,:)*1e3);
% set(gca,'FontSize',16)
% xlabel('Tempo [s]')
% ylabel('x_2: Posição da armadura [mm]')
% grid on
% 
% figure(3)
% plot(t,x(3,:));
% set(gca,'FontSize',16)
% xlabel('Tempo [s]')
% ylabel('x_3 Velocidade da armadura [m/s]')
% grid on
% 
% % figure(4)
% % plot3(x(1,:),x(2,:),x(3,:));
% % set(gca,'FontSize',16)
% % xlabel('x_1')
% % ylabel('x_2')
% % zlabel('x_3')
% % grid on

%% Plots para testar superposicao
% figure(1)
% yyaxis left
% plot(t,x(2,:)*1e3,'LineWidth',2);
% hold on
% plot(t,(x1(2,:)+x2(2,:))*1e3,'LineWidth',2);
% plot(t, x1(2,:)*1e3,'LineWidth',2);
% plot(t, x2(2,:)*1e3,'LineWidth',2);
% ylabel('x_2: Posição da armadura [mm]')
% yyaxis right
% plot(t, u,'LineWidth',1);
% plot(t, u1+u2,'LineWidth',2);
% plot(t, u1,'LineWidth',2);
% plot(t, u2,'LineWidth',2);
% ylim([0 11])
% ylabel('Amplitude do Degrau')
% set(gca,'FontSize',16)
% xlabel('Tempo [s]')
% legend('x_2^{10}', 'x_2^7+x_2^3', 'x_2^7', 'x_2^3','u_{10}', 'u_7+u_3', 'u_7', 'u_3', 'Location','northwest', 'NumColumns',2)
% grid on
% hold off

%% Plot para diversos pontos de operacao
yyaxis left
plot(t, x(2,:)*1e3,'LineWidth',2)
ylabel('x_2: Posição da armadura [mm]')
yyaxis right
plot(t, u,'LineWidth',2)
set(gca,'FontSize',14)
ylabel('Entrada [V]')
xlabel('Tempo [s]')
grid on

