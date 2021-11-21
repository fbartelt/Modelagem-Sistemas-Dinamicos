% Exemplo 5.5.2

% (c) Luis Aguirre, 2002

% representacao do circuito em espaco de estados
% Escolha de parametros eletricos hipoteticos
L=0.5;
c=0.8;
R=0.2;

% Matrizes de estado
A=[-R/L -1/L;1/c 0];
B=[1/L;0];
C=[1 0;0 1];
D=[0;0];

% Gera entrada
v=prbs(2000,7,8);

% Simula
inc=0.05;
t=0:inc:length(v)*inc-inc;
[y,x]=lsim(A,B,C,D,v',t',[0;0]);

% Atribuicao de nomes
i=y(:,1);
vc=y(:,2);
vc_ponto=i./c;
i_ponto=-R*i./L-vc./L+v'./L;
vc_dois_pontos=i_ponto./c;

% Montar vetor e matriz de regressores
% indices onde serao tomadas as restricoes:
In=10:250:2000;
Vc=vc(In);
V=v(In)';
Vcp=vc_ponto(In);
Vc2p=vc_dois_pontos(In);

Psi=[V -Vcp -Vc2p];
Teta=inv(Psi'*Psi)*Psi'*(Vc./c);
