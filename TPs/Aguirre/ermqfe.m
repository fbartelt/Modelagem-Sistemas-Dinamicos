% Estimação recursiva de parâmetros variantes no tempo
%
% (c) Luis Antonio Aguirre, BH 17/06/99
clear
% carregar dados
load f0707.dat;

% fator de dizimação a ser usado
fd=1;
ini = 20;
teta(:,ini-1)=[(0.9+0.94)/2; 2.5/13];
% período de amostragem em segundos
Ts=fd*70;

% referir os dados ao ponto de operação em que foi feito o teste
[a b]=size(f0707);
u=f0707(1:fd:a,1);
y=f0707(1:fd:a,2);

%[a b]=size(teta);
figure(1)
subplot(211);
plot(Ts:Ts:(b-2*ini+2)*Ts,u(ini-1:b-ini),'k-');
ylabel('razão cíclica (%)');
title('(a)');
%axis([0 15000 2 14])
subplot(212)
plot(Ts:Ts:(b-2*ini+2)*Ts,y(ini-1:b-ini),'k-');
ylabel('temperatura (%)');
xlabel('tempo (s)');
title('(b)');
%axis([0 15000 25 60])

u=f0707(1:fd:a,1)-f0707(1,1);
y=f0707(1:fd:a,2)-f0707(1,2);
clear f0707.dat;

% condições iniciais
P=eye(2)*10^6;
ini=20;
% estimativa de constante de tempo inicial approx 900s
% estimativa de ganho inicial 25/10=2,5
teta(:,ini-1)=[(0.9+0.94)/2; 2.5/13];

% fator de esquecimento
lambda=0.99;

% Algoritmo recursivo
for k=ini:length(y);
   psi_k=[y(k-1);u(k-1)];
   K_k = (P*psi_k)/(psi_k'*P*psi_k+lambda);
   teta(:,k)=teta(:,k-1)+K_k*(y(k)-psi_k'*teta(:,k-1));
   P=(P-(P*psi_k*psi_k'*P)/(psi_k'*P*psi_k+lambda))/lambda;
end;

[a b]=size(teta);

% constante de tempo e ganho
for k=1:b   
   tau(k)=-Ts/(teta(1,k)-1);
   ganho(k)=tau(k)*teta(2,k)/Ts;
end;
figure(2)
subplot(211);
plot(Ts:Ts:(b-2*ini+2)*Ts,tau(ini-1:b-ini),'k-');
ylabel('constante de tempo (s)');
title('(a)');
subplot(212)
plot(Ts:Ts:(b-2*ini+2)*Ts,ganho(ini-1:b-ini),'k-');
ylabel('ganho');
xlabel('tempo (s)');
title('(b)');
axis([0 15000 2 3]);



tau99=tau;

% condições iniciais
P=eye(2)*10^6;
ini=20;
teta(:,ini-1)=[(0.9+0.94)/2; 2.5/13];
% fator de esquecimento
lambda=1;

% Algoritmo recursivo
for k=ini:length(y);
   psi_k=[y(k-1);u(k-1)];
   K_k = (P*psi_k)/(psi_k'*P*psi_k+lambda);
   teta(:,k)=teta(:,k-1)+K_k*(y(k)-psi_k'*teta(:,k-1));
   P=(P-(P*psi_k*psi_k'*P)/(psi_k'*P*psi_k+lambda))/lambda;
end;

% constante de tempo e ganho
for k=1:b   
   tau1(k)=-Ts/(teta(1,k)-1);
end;

% condições iniciais
P=eye(2)*10^6;
ini=20;
teta(:,ini-1)=[(0.9+0.94)/2; 2.5/13];

% fator de esquecimento
lambda=0.97;

% Algoritmo recursivo
for k=ini:length(y);
   psi_k=[y(k-1);u(k-1)];
   K_k = (P*psi_k)/(psi_k'*P*psi_k+lambda);
   teta(:,k)=teta(:,k-1)+K_k*(y(k)-psi_k'*teta(:,k-1));
   P=(P-(P*psi_k*psi_k'*P)/(psi_k'*P*psi_k+lambda))/lambda;
end;

% constante de tempo e ganho
for k=1:b   
   tau97(k)=-Ts/(teta(1,k)-1);
end;
figure(3)
plot(Ts:Ts:(b-2*ini+2)*Ts,tau1(ini-1:b-ini),'k:',Ts:Ts:(b-2*ini+2)*Ts,tau99(ini-1:b-ini),'k--',Ts:Ts:(b-2*ini+2)*Ts,tau97(ini-1:b-ini),'k-');
ylabel('constante de tempo (s)');
xlabel('tempo (s)');
axis([0 15000 400 1400]);
text(11000,1000,'lambda=1');
text(11000,1130,'lambda=0,99');
text(8500,1250,'lambda=0,97');
hold on
plot([8300 8300],[400 1400],'k:');
hold off


