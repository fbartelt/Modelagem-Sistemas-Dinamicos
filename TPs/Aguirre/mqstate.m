% Estimacao recursiva das matrizes de estado de um sistema discreto

% (c) Luis A Aguirre, BH 03/06/99

phi=[0.8 0.1 0.005;0 0.7 0.1;0 0 0.9];
gama=[1 0.05 0.0025;0 1 0.05;0 0 1];

% vetor de entradas
u=randn(3,300);

% inicializacao
x=[0;0;0];

% simulacao
for k=2:300
  x(:,k)=phi*x(:,k-1)+gama*u(:,k-1)+0.2*randn(3,1);
end;

% algoritmo recursivo

% inicializacao
delta=zeros(6,3);
P=1000*eye(6);

% desprezar os primeiros 9 conjuntos de dados
for k=10:200
m=[x(:,k-1)' u(:,k-1)'];
K=P*m'/(1+m*P*m');
delta=delta+K*(x(:,k)'-m*delta);
P=P-(P*m'*m*P')/(1+m*P*m');
end;

Phi=delta(1:3,1:3)';
Gama=delta(4:6,1:3)';

% simulacao do modelo estimado
X=[0;0;0];
for k=2:300
  X(:,k)=Phi*X(:,k-1)+Gama*u(:,k-1);
end;

% Comparando a saída, supondo matriz de observacao [1 1 1]
plot(1:100,[1 1 1]*x(:,201:300),'k-',1:100,[1 1 1]*X(:,201:300),'k-.');
xlabel('número de amostras');



