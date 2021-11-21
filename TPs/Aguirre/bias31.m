% investiga o problema de polarização

% condicoes iniciais
y3=zeros(2,1);

% entradas aleatorias
u1=randn(502,1);

for i=1:400
   
% ruído 
nu=0.2*randn(502,1);


% geração de ruído MA
for k=2:length(nu)
  e(k)=0.8*nu(k-1)+nu(k);
end;


% outro modelo com erro de modelagem MA
for k=2:502
 y3(k) =0.7*y3(k-1)+0.5*u1(k-1)+e(k);
end;

% estimação de parametros usando estimador MQ

Psi3=[y3(2:501) u1(2:501)];
V3=y3(3:502);
teta3(:,i)=inv(Psi3'*Psi3)*Psi3'*V3;
xi=V3-Psi3*teta3(:,i);

%estimação usando o estimador EMQ
%usamos a ultima matriz de regressores do estimador MQ como ponto de partida
xi2=[0; xi(1:499)];
for k=1:10 % numero de iteracoes precisa ser suficiente para convergencia
 PsiE=[Psi3(3:500,:) xi2(2:499) xi2(1:498)];
 tetaE(:,i)=inv(PsiE'*PsiE)*PsiE'*V3(3:500);
 xi2=[0;0; V3(3:500)-PsiE*tetaE(:,i)];
% PsiE=[Psi3(2:499,:) xi2(1:498)];
% tetaE(:,i)=inv(PsiE'*PsiE)*PsiE'*V3(2:499);
% xi2=[0; V3(2:499)-PsiE*tetaE(:,i)];
end;
end;

figure(1)
hist(teta3(1,:),100);
figure(2)
hist(teta3(2,:),100);
pause;


[t,r,l]=myccf([e' e'],10,0,0,'k');
[t,rx,l]=myccf([xi xi],10,0,0,'k');
plot(t,r*std(e)^2,'ko',t,rx*std(xi)^2,'k+',[t(1) t(length(t))],[l l]*std(xi)^2,'k:',[t(1) t(length(t))],[-l -l]*std(xi)^2,'k:');
xlabel('atrasos (lags)')
hold on;
plot(t,r*std(e)^2,'k-',t,rx*std(xi)^2,'k--',[t(1) t(length(t))],[l l]*std(xi)^2,'k:',[t(1) t(length(t))],[-l -l]*std(xi)^2,'k:');
hold off
pause;


%gera ruido estimado (xi) a partir de residuos (xi2) de modelo nao polarizado
Xi(1:2)=xi2(1:2);
for k=3:499
 Xi(k)=mean(tetaE(3,:))*xi2(k-1)+mean(tetaE(4,:))*xi2(k-2)+xi2(k);
end;

[t,r,l]=myccf([e' e'],10,0,0,'k');
[t,rx,l]=myccf([Xi' Xi'],10,0,0,'k');
plot(t,r*std(e)^2,'ko',t,rx*std(Xi)^2,'k+',[t(1) t(length(t))],[l l]*std(Xi)^2,'k:',[t(1) t(length(t))],[-l -l]*std(Xi)^2,'k:');
xlabel('atrasos (lags)')
hold on;
plot(t,r*std(e)^2,'k-',t,rx*std(Xi)^2,'k--',[t(1) t(length(t))],[l l]*std(Xi)^2,'k:',[t(1) t(length(t))],[-l -l]*std(Xi)^2,'k:');
hold off
pause;


subplot(221)
hist(tetaE(1,:),100);
title('(a)');
subplot(222)
hist(tetaE(2,:),100);
title('(b)');
axis([0.46 0.54 0 15]);
subplot(223)
hist(tetaE(3,:),100);
title('(c)')
subplot(224)
hist(tetaE(4,:),100);
title('(d)');



