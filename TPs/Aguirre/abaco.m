% Gera ábacos e gráficos correlatos do capítulo 2
% (c) Luis Antonio Aguirre

% Resposta ao degrau de sistemas de primeira ordem
t=0:0.1:40;
y=step(1,[10 1],t);
plot(t,y,'w-');
xlabel('t');
ylabel('y(t)');
axis([0 40 0 1.2]);
text(10,0.55,'tau=10s');
text(1,0.95,'tau=2s');
hold on
for tau=2:2:8
  y=step(1,[tau 1],t);
  plot(t,y,'w-');
end;
hold off;
pause


% Resposta ao degrau de sistemas de segunda ordem
t=0:0.1:15;
[n,d]=seg_ord(1,0.1,1);
y=step(n,d,t);
plot(t,y,'w-');
xlabel('t');
ylabel('y(t)');
axis([0 15 0 1.8]);
text(2.5,0.55,'zeta=1,3');
text(0.6,1.65,'zeta=0,1');
hold on
for zeta=0.4:0.3:1.3
 [n,d]=seg_ord(1,zeta,1);
 y=step(n,d,t);
 plot(t,y,'w-');
end;
hold off
pause


% Resposta em frequência de k, 1/jw e jw
w=logspace(-1,1);
[mag1,pha1]=bode(1,1,w);
[mag2,pha2]=bode(1,[1 0],w);
[mag3,pha3]=bode([1 0],1,w);
subplot(211)
semilogx(w,20*log10(mag1),'w-',w,20*log10(mag2),'w--',w,20*log10(mag3),'w-.');
ylabel('|H(jw)| (dB)');
title('(a)');
subplot(212)
semilogx(w,pha1,'w-',w,pha2,'w--',w,pha3,'w-.');
xlabel('freqüência')
ylabel('fase de H(jw) (graus)');
title('(b)')
pause


% Resposta em frequencia de sistemas de primeira ordem
w=logspace(-2,1);
[mag,pha]=bode(1,[10 1],w);
subplot(211)
semilogx(w,20*log10(mag),'w-');
ylabel('|H(jw)| (dB)');
title('(a)');
axis([0.01 10 -40 10]);
text(0.3,-20,'tau=10s');
text(1,-2,'tau=2s');
hold on
subplot(212)
semilogx(w,pha,'w-');
xlabel('freqüência')
ylabel('fase de H(jw) (graus)');
title('(b)')
text(0.08,-70,'tau=10s');
text(0.8,-45,'tau=2s');
hold on
for tau=2:2:8
  [mag,pha]=bode(1,[tau 1],w);
  subplot(211)
  semilogx(w,20*log10(mag),'w-');
  subplot(212)
  semilogx(w,pha,'w-');
end;
hold off;
pause


% Resposta em frequencia de 1+jw\tau
w=logspace(-2,1);
[mag,pha]=bode([10 1],1,w);
subplot(211)
hold off
semilogx(w,20*log10(mag),'w-');
ylabel('|H(jw)| (dB)');
title('(a)');
axis([0.01 10 -12 40]);
text(0.3,20,'tau=10s');
text(1,2,'tau=2s');
hold on
subplot(212)
semilogx(w,pha,'w-');
xlabel('freqüência')
ylabel('fase de H(jw) (graus)');
title('(b)')
text(0.08,70,'tau=10s');
text(0.8,45,'tau=2s');
hold on
for tau=2:2:8
  [mag,pha]=bode([tau 1],1,w);
  subplot(211)
  semilogx(w,20*log10(mag),'w-');
  subplot(212)
  semilogx(w,pha,'w-');
end;
hold off;
pause


% Resposta em frequencia de sistemas de segunda ordem
w=logspace(-1,1,100);
[n,d]=seg_ord(1,0.1,1);
[mag,pha,w]=bode(n,d,w);
subplot(211)
hold off
semilogx(w,20*log10(mag),'w-');
ylabel('|H(jw)| (dB)');
title('(a)');
axis([0.1 10 -40 20]);
text(1.2,14,'zeta=0,1');
text(0.4,-12,'zeta=1,3');
hold on
subplot(212)
semilogx(w,pha,'w-');
xlabel('freqüência')
ylabel('fase de H(jw) (graus)');
title('(b)')
axis([0.1 10 -200 10]);
text(0.9,-30,'zeta=0,1');
text(0.2,-70,'zeta=1,3');
hold on;
for zeta=0.4:0.3:1.3
 [n,d]=seg_ord(1,zeta,1);
 [mag,pha,w]=bode(n,d,w);
 subplot(211)
 semilogx(w,20*log10(mag),'w-');
 subplot(212)
 semilogx(w,pha,'w-');
end
hold off

