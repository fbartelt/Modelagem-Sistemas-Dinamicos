% Simulação da planta de bombeamento modelada na secao 1.4
% Esta rotina chama a funcao mypant.m onde esta codificado o modelo
% Copyright (c) 1998 por Luis A. Aguirre. Todos os direitos reservados.
% Versao original codificada por Rodolfo Siqueira Santana, 1998.

% Declaracao de variaveis globais
global Us
global k

% Definicao de parametros de simulacao 
to =0;		% tempo inicial em segundos
tff = 1500; 	% tempo de simulacao em segundos
TD = 100;	% tempo em que ocorre o degrau
degrau=1;	% execute esta linha se o degrau for positivo
%degrau=-1;	% execute esta linha se o degrau for negativo

%%%%%%%%%%%  Inicio da secao que gera o sinal de entrada %%%%%%%%%
%
if degrau>0
  % Para degraus positivos
  SP1 = 16.34;	% patamar inferior do degrau de entrada
  % A correção no ganho DC é feita aqui conforme descrito na secao 1.4.4
  SP2 = 17.05*1.012; % patamar superior
  % amplitude do degrau original= 17.05-16.34=0.71
  % amplitude do degrau (ganho corrigido) = 17.05*1.012-16.34=0.91
  % Condicao inicial de nivel
  xo = 0.075;	% metros
else
  % Para degraus negativos
  SP2 = 16.34/1.012;
  SP1 = 17.05;
  % amplitude do degrau original= -17.05+16.34=-0.71
  % amplitude do degrau (ganho corrigido) = -17.05+16.34/1.012=-0.91
  % Condicao inicial de nivel e tempo
  xo = 0.27;	
end

k = 0;
% criacao do vetor de entrada Us em degrau
for i = to:tff
	k = k + 1;
	if k < TD
	  Us(k) = SP1;
	else 
	  Us(k) = SP2;
    end
end
%
%%%%%%%%%%%  Fim da secao que gera o sinal de entrada %%%%%%%%%


% Simula o modelo codificado em myplant 
[t,x] = ode23('eq_bomb',[to tff],xo);

% Carrega os dados medidos para comparacao
if degrau>0
  % Dados com o degrau positivo (figuras 1.9/1.10)
  load ens_25.dat;
  h_ens=ens_25(:,2);
  t_ens=ens_25(:,1);
else
  % Dados com o degrau negativo (figura 1.11)
  load ens_26.dat;
  h_ens=ens_26(:,2);
  t_ens=ens_26(:,1);
end


% gera figura 1.10 ou figura 1.11
plot(t_ens(1:5:1500),h_ens(1:5:1500),'r:',t,[ones(1,1)*x(1); x(1:length(t)-1)],'b-');
axis([0 1500 0 0.4]);
xlabel('tempo (s)');
ylabel('Nível em m');
