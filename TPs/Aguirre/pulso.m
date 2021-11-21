% modela a planta de bombeamento a partir da resposta ao pulso
% (c) Luis Aguirre BH 2/03/99

% Carrega os dados
load pulso10.dat;
t=pulso10(:,1);
vazao=pulso10(:,2);
nivel=pulso10(:,3);

% Tempo de amostragem Ts=1.0440s

% filtragem
for i=1:50
   v(i)=mean(vazao(i:i+8));
end;
for i=51:60
   v(i)=vazao(i);
end;
for i=61:length(vazao)-8
   v(i)=mean(vazao(i:i+8));
end;

% sinal de nivel
% Na realidade este deve ser o sinal que foi para o inversor
% a divisao por 250 muda a escala de 1 a 5V para 4 a 20mA.
n=nivel(2:731)*1000/250;
plot(t(2:731),vazao(2:731),t(2:731),v(1:730));
pause
subplot(211)
plot(t(1:730),n(1:730),'b');
axis([1 730 14 21]);
ylabel('mA');
title('(a)')
subplot(212)
plot(t(2:731),v(1:730),'b');
ylabel('volts');
xlabel('t(s)');
title('(b)')
axis([1 730 3.6 4.0]);
pause

h=1.0440; % tempo de amostragem
d=1; % fator de dizimacao
fim=2^14-730+49;
wv=[v(50:730) ones(1,fim)*3.83];
wn=[n(47:727)' ones(1,fim)*3.78*1000/250];


V=fft(wv);
N=fft(wn);
t=0:d*1.0440:1.0440*(length(wv)-1);
f = 0:1/t(length(t)):1/(d*h);
f=2*pi*f;

figure(1);
subplot(211);
semilogx(f,abs(N),'b-');
title('(a)');
ylabel('|U(jw)|');
subplot(212);
semilogx(f,abs(V),'b-');
title('(b)');
ylabel('|Y(jw)|');
xlabel('w (rad/s)');
pause;

%so interessa pontos onde os espectros são significativos

%fmax=input('frequencia max: ');
fmax=0.015;
ind=find(f>=fmax);
xind=ind(1);
fft_v=V(1:xind);
fft_n=N(1:xind);
f=f(1:xind);

fft_tf = fft_v ./ fft_n;
fft_tfa = abs(fft_tf);
fft_tff = angle(fft_tf);

figure(1);
subplot(211);
semilogx(f, 20*log10(fft_tfa), 'b');
ylabel('|H(jw)| (dB)');
title('(a)');
subplot(212);
semilogx(f, (180/pi)*fft_tff, 'b');
xlabel('w (rad/s)');
title('(b)');
ylabel('fase de H(jw) (graus)');
pause;

% Metodo de Levy
lf=length(f);
l0=lf;
l2=sum(f.^2);
l4=sum(f.^4);

T1=imag(fft_tf)*f';
T3=imag(fft_tf)*(f.^3)';
T5=imag(fft_tf)*(f.^5)';

S0=real(fft_tf)*(f.^0)';
S2=real(fft_tf)*(f.^2)';
S4=real(fft_tf)*(f.^4)';
S6=real(fft_tf)*(f.^6)';

U2=(real(fft_tf).^2+imag(fft_tf).^2)*(f.^2)';
U4=(real(fft_tf).^2+imag(fft_tf).^2)*(f.^4)';
U6=(real(fft_tf).^2+imag(fft_tf).^2)*(f.^6)';

% primeira ordem
M=[l0 T1;T1 U2];
vec=[S0;0];
teta1=inv(M)*vec;
num1=teta1(1);
den1=[teta1(2) 1];
[mag1 phase1]=bode(num1,den1,f);

% segunda ordem
M=[l0 0 T1 S2;0 l2 -S2 T3;T1 -S2 U2 0; S2 T3 0 U4];
vec=[S0;T1; 0;U2];
teta2=inv(M)*vec;
num2=teta2(2:-1:1)';
den2=[teta2(4:-1:3)' 1];
[mag2 phase2]=bode(num2,den2,f);

% terceira ordem
M=[l0 0 -l2 T1 S2 -T3;0 l2 0 -S2 T3 S4;l2 0 -l4 T3 S4 -T5;T1 -S2 -T3 U2 0 -U4];
M=[M; S2 T3 -S4 0 U4 0;T3 -S4 T5 U4 0 -U6];
vec=[S0;T1; S2; 0;U2; 0];
teta3=inv(M)*vec;
num3=teta3(3:-1:1)';
den3=[teta3(6:-1:4)' 1];
[mag3 phase3]=bode(num3,den3,f);


figure(2);
subplot(211);
semilogx(f, 20*log10(fft_tfa), 'b',f, 20*log10(mag1), 'b--',f, 20*log10(mag2), 'b-.',f, 20*log10(mag3), 'b:');
ylabel('|H(jw)| (dB)');
title('(a)');
subplot(212);
semilogx(f, (180/pi)*fft_tff, 'b',f, phase1, 'b--',f, phase2, 'b-.',f, phase3, 'b:');
xlabel('w (rad/s)');
ylabel('fase de H(jw) (graus)');
title('(b)');



 