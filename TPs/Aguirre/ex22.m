% Exemplos do capítulo 2
% (c) Luis Antonio Aguirre 

% Exemplo de decomposição
num=[5 3];
den=poly([-1 -2 -3]);
num1=-1;
den1=[1 1];
num2=7;
den2=[1 2];
num3=-6;
den3=[1 3];

t=0:0.025:3;
y=impulse(num,den,t);
y1=impulse(num1,den1,t);
y2=impulse(num2,den2,t);
y3=impulse(num3,den3,t);
plot(t,y,'w-',t,y1,'w-.',t,y2,'w--',t,y3,'w:');
xlabel('t');
pause;

% Exemplo de redução usando MDI
% Função de transferência original
num=[7.526 10.357 0.92051 0.63827];
den=[20.0141 11.749 16.987 1.1602 1];

% Funções de transferência reduzidas
num1=[2.9668e-2 4.1527e-2];
den1=[1 2.8134e-2 6.5061e-2];
num2=[0.49496 0.49017];
den2=[1 0.55890 0.76797];

t=0:0.1:50;
y=step(num,den,t);
y1=step(num1,den1,t);
y2=step(num2,den2,t);
plot(t,y,'w-',t,y1,'w-.',t,y2,'w--');
xlabel('t');
pause

w=logspace(-1,1,100);
[mag,phase]=bode(num,den,w);
subplot(211)
semilogx(w,20*log10(mag),'w-');
ylabel('20 log |H(jw)|');
subplot(212)
semilogx(w,phase,'w-');
ylabel('graus');
xlabel('freqüência (Hz)');



