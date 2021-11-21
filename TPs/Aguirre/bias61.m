% investiga o problema de polarização em modelos de erro na
% saida e ARX, cujos parâmetros são estimados por MQ

% Luis Aguirre, 1999

% condicoes iniciais
y1=zeros(2,1);
y2=zeros(2,1);

% entradas aleatorias
u1=randn(502,1);

% sistema que será usado com erro de observacao em y
for k=3:502
   y1(k) = 1.5*y1(k-1)-0.7*y1(k-2)+ 0.5*u1(k-1);
end;


for i=1:200
   
% nivel de ruído usado em observacoes de y
e=0.2*randn(502,1);

% erro de observacao em y
Y1=y1+e;

% sistema com erro de modelagem
for k=3:502
 y2(k) = 1.5*y2(k-1)-0.7*y2(k-2)+0.5*u1(k-1)+e(k);
end;


% estimação de parametros
Psi1=[Y1(2:501) Y1(1:500) u1(2:501)];
V1=Y1(3:502);
teta1(:,i)=inv(Psi1'*Psi1)*Psi1'*V1;

Psi2=[y2(2:501) y2(1:500) u1(2:501)];
V2=y2(3:502);
teta2(:,i)=inv(Psi2'*Psi2)*Psi2'*V2;
end


subplot(321)
hist(teta1(1,:),100);
ylabel('(a)')
axis([1.3 1.55 0 10]);
subplot(323)
hist(teta1(2,:),100);
ylabel('(b)')
axis([-0.75 -0.45 0 10]);
subplot(325)
hist(teta1(3,:),100);
ylabel('(c)')
axis([0.45 0.55 0 10]);

subplot(322)
hist(teta2(1,:),100);
ylabel('(d)');
axis([1.3 1.55 0 10]);
subplot(324)
hist(teta2(2,:),100);
ylabel('(e)')
axis([-0.75 -0.45 0 10]);
subplot(326)
hist(teta2(3,:),100);
ylabel('(f)')
axis([0.45 0.55 0 10]);


