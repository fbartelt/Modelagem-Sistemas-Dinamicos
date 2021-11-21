% exemplo secao 3.2.2
% (c) Luis Antonio Aguirre

t=0:0.1:35;
[n1,d1]=seg_ord(1,0.15,1);
[n2,d2]=seg_ord(1,0.4,1);
y1=step(n1,d1,t);
y2=step(n2,d2,t);
plot(t,y1,'w-',t,y2,'w--');
xlabel('t');
ylabel('y(t)');

