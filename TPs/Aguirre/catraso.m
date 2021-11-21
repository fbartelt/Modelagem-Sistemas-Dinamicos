% gera grafico com atraso puro de tempo

t=0:0.01:4;
[nd1,dd1]=delay(1,1);
[nd2,dd2]=delay(1,2);
[nd3,dd3]=delay(1,3);
num=1;
den=[1 1];
y1=step(conv(num,nd1),conv(den,dd1),t);
y2=step(conv(num,nd2),conv(den,dd2),t);
y3=step(conv(num,nd3),conv(den,dd3),t);
y=step(num,den,t);
lt=length(t);

plot(t,y1,'w:',t,y2,'w-.',t,y3,'w--',t,[zeros(100,1); y(1:lt-100)],'w-');
xlabel('t');
pause

w=logspace(-1,1.3,100);
[mag1,pha1]=bode(nd1,dd1,w);
[mag2,pha2]=bode(nd2,dd2,w);
[mag3,pha3]=bode(nd3,dd3,w);
pha4=-w*57.3;
subplot(211)
semilogx(w,20*log10(mag1),'w-',w,20*log10(mag2),'w--',w,20*log10(mag3),'w-.');
axis([0.1 20 -1 1]);
ylabel('|H(jw)| (dB)');
title('(a)');
subplot(212)
semilogx(w,pha4,'w-',w,pha1,'w:',w,pha2,'w-.',w,pha3,'w--');
axis([0.1 20 -800 0]);
xlabel('freqüência')
ylabel('fase de H(jw) (graus)');
title('(b)')
pause
print -deps2 c:\pctex\book2\figures\fatraso.eps
