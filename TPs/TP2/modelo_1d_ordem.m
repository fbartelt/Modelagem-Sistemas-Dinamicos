addpath(genpath('./Tarefa2_DadosExperimentais'))
S = load('ENS_25.DAT');
V = load('ENS_26.DAT'); % para validacao

t = S(:, 1);
y = S(:, 2);

tv = S(:, 1); % validation data
yv = S(:, 2); % validation data
%% Visualizacao


%% Parametros
to =0;		% tempo inicial em segundos
tff = t(end); 	% tempo de simulacao em segundos
TD = 100;	% tempo em que ocorre o degrau
u25 = 16.34*ones(size(t)); % 16.34 a 17.05 mA é aplicado
u25(t > TD) = 17.05;
%u25 = zeros(size(t)); % 16.34 a 17.05 mA é aplicado
%u25(t > TD) = 0.3;
u26 = 17.05*ones(size(t));% 17.05 a 16.34 mA
u26(t > TD) = 16.34;


K = (y(end) - y(1))/(u25(end) - u25(1))/24;
% tau_d = 26
% tau_d = 126; %s
% tau = 370 - tau_d;
tau_d = 120;
%tau = 414.8;
tau = 373.32; % 90% do anterior
sys = tf(K, [tau 1])
s = tf('s');
delay = pade(exp(-tau_d * s), 5);
sys = sys * delay
[yt, ~] = lsim(sys, u25, t);
offset = y(1);
yt = yt + offset;
figure(1)

%yyaxis left
plot(t,y)
xlabel('time[s]')
ylabel('Saida')
hold on
plot(t, yt)
hold off
%yyaxis right
%plot(t, u25)
%ylabel('Corrente [mA]')
xlim([to tff])
grid on