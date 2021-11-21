% investiga o problema de erro nas variaveis


% condicoes iniciais
y=zeros(2,1);

% entradas aleatorias
ui=randn(502,1);


% sistema que será usado com erro de observacao em y
for k=3:502
   y(k) = 1.5*y(k-1)-0.7*y(k-2)+ 0.5*ui(k-1);
end;
yi=y;

for i=1:400
   
% nivel de ruído usado em observacoes de y
  ey=0.2*randn(502,1);
% nivel de ruído usado em observacoes de u
  eu=0.2*randn(502,1);


% erro de observacao em y
y=yi+ey;
% erro de observacao em u
u=ui+eu;


Psi=[y(2:501) y(1:500) u(2:501)];
V=y(3:502);
teta(:,i)=inv(Psi'*Psi)*Psi'*V;

[U,S,V] = svd([Psi V]);
% Duas formas diferentes de calcular a solucao
% tetas(:,i)=inv(Psi'*Psi-S(4,4)^2*eye(3))*Psi'*V;
tetas(:,i)=-V(1:3,4)/V(4,4);

end


subplot(321)
hist(teta(1,:),100);
%axis([1.3 1.55 0 20]);
ylabel('(a)');
subplot(322)
hist(tetas(1,:),100);
%axis([1.3 1.55 0 20]);
ylabel('(b)');
subplot(323)
hist(teta(2,:),100);
%axis([-0.75 -0.45 0 20]);
ylabel('(c)');
subplot(324)
hist(tetas(2,:),100);
%axis([-0.75 -0.45 0 20]);
ylabel('(d)');
subplot(325)
hist(teta(3,:),100);
%axis([0.4 0.6 0 20]);
ylabel('(e)');
subplot(326)
hist(tetas(3,:),100);
%axis([0.4 0.6 0 20]);
ylabel('(f)');

  