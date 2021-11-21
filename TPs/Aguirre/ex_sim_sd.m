% Exemplo de como simular um sistema dinamico sem e com ruido, incluindo
% caso de aproximacao discreta

close all; clear all;

% Sistema verdadeiro G(s) = K/(tau*s +1):
tau = 20;
K = 5;
theta = 10; % Atraso puro de tempo em segundos
Gs = tf([K], [tau 1], 'ioDelay', theta);

% Sistema discretizado equivalente:
% y(k) = a*y(k-1) + b*u(k-1)
Ts = tau/10 % Intervalo de amostragem Ts = 2s
d = theta/Ts % Atraso puro de tempo em amostras d = 5 amostras
a = 1 - (Ts/tau); % Tarefa: Como obter esses valores? Ver Exemplo 8.5.1 do livro
b = (Ts*K)/tau;

%% Simulacao sem ruido:

% Dados de entrada:
N = 400;
t = [0:Ts:(N-1)*Ts]';
Tb = tau/5 % Tempo de permanencia em s Tb = 4s
m = Tb/Ts  % Tempo de permanencia em amostras m = 2 amostras
% m = 1;
u = prbs(N,10,m); u = u(:);

stairs(t, u); xlabel('t (segundos)'); ylabel('u(kT_s)');

% Simulacao 1:
y1 = lsim(Gs, u, t);

% Simulacao 2:
y2(1:d+1) = zeros(d+1,1);
for k = 2+d:N
    y2(k) = a*y2(k-1) + b*u(k-1-d);
end
y2 = y2(:);

figure;
stairs(t, y1, 'b'); hold on;
stairs(t, y2, 'r'); 
xlabel('t (segundos)'); ylabel('y_1(t), y_2(kT_s)');

%% Simulacao com ruido na equacao (de processo):

w = 0.01*randn(length(u),1);

% Simulacao 1:
yw1 = lsim(Gs, u + w, t);
% Ou seja, conhecemos u, mas na pratica o sistema e excitado por u + w

% Simulacao 2:
yw2(1:d+1) = zeros(d+1,1);
for k = 2+d:N
    yw2(k) = a*yw2(k-1) + b*(u(k-1-d) + w(k-1));
end
yw2 = yw2(:);


hold on;
stairs(t, yw1, 'g'); hold on;
stairs(t, yw2, 'k'); 
xlabel('t (segundos)'); ylabel('y_1(t), y_2(kT_s)');


%% Simulacao com ruido na saida (de medicao):

r = 0.05*randn(length(y1),1); % Como o ganho do sistema eh  K = 5, escolhi um ruido 5 vezes maior que o de processo

% Simulacao 1:
y1 = lsim(Gs, u, t);
yr1 = y1 + r; % Observe que ruido afeta apenas a medicao da saida

% Simulacao 2:
y2(1:d+1) = zeros(d+1,1);
for k = 2+d:N
    y2(k) = a*y2(k-1) + b*u(k-1-d);
end
y2 = y2(:);
yr2 = y2 + r;

hold on;
stairs(t, yr1, 'c'); hold on;
stairs(t, yr2, 'm'); 
xlabel('t (segundos)'); ylabel('y_1(t), y_2(kT_s)');


%% Estimacao de atraso puro de tempo:
figure; myccf([u, y2],100,1,1,'k'); % Dados sem ru?do
hold on; myccf([u, yr2],100,1,1,'m'); % Dados com ru?do

% Observar que o m?ximo tem ocorrido em torno de lag = -7. Sempre haver? um
% atraso m?nimo de -1 pelo fato do sistema ser discreto no tempo. Portanto,
% podemos estimar aqui dest = 6, valor pr?ximo do verdadeiro d = 5.
