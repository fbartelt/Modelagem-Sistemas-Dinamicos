%Felipe Bartelt de Assis Pessoa - 2016026841
% Script adaptado das rotinas fornecidas por Luis Antonio Aguirre

clear
close all
%% E 3.10
u=randn(1,100);      % para a), c)
%u = [1 zeros(1,99)]; % para b), d)
lu=length(u);
y=dlsim(1,[1 0.2 0.8],u);
e= 0.05*max(u)*randn(size(u')); % para c), d)
%e=zeros(size(u)); % para a), b)

U=u';
for i=1:99
U=[U [zeros(i,1); u(1:lu-i)']];
end

h=inv(U)*(y+e);

figure(2)
subplot(311)
plot(1:lu, u, 'LineWidth', 1.5)
title('Entrada');
xlim([1 lu])
subplot(312)
plot(1:lu,y,'-', 'LineWidth', 1.5);
hold on
plot(1:lu,U*h,'m--', 'LineWidth', 1.5)
hold off
title('Saida');
legend('Real', 'Estimada')
xlim([0 lu]);
xlabel('amostras');
subplot(313)
stem(0:9, h(1:10))
xlim([-0.2, 9.2])
title('H estimado')
sgtitle('Exercicio 3.10')
a = funcmassa(3)
%% E 4.15
u=prbs(378,6,1);
lu=length(u);

figure(2)
subplot(211)
plot(1:190,u(1:190),'k-');
axis([0 190 -0.2 1.2]);
title('Sinal PRBS');
xlabel('Amostras')

lag = 40;
[t,ruu,l,~]=myccf([u' u'],lag,0,0,'k');
subplot(212)
set(gca,'FontSize',18)
plot([0 lag],[l l],'b--',[0 lag],[-l -l],'b--');
axis([-0.5 lag+0.5 -0.2 1.1]);
hold on
stem(t,ruu,'k');
hold off
xlabel('Atraso');
title('Funcao de autocorrelacao');
sgtitle('Exercicio 4.15')
%% E 4.16
tb_list = 6*[1, 2, 3, 4];

for i=1:length(tb_list)
    tb = tb_list(i)
    u=prbs(378,6, tb);
    lu=length(u);
    figure(2+i)
    subplot(311)
    plot(1:lu,u,'k-');
    axis([0 lu -0.2 1.2]);
    title(['Sinal PRBS (Tb = ', num2str(tb), ')']);
    xlabel('Amostras')

    subplot(312)
    lag = 40;
    [t,ruu,l,B1]=myccf([u' u'],lag,0,0,'k');
    set(gca,'FontSize',18)
    plot([0 lag],[l l],'b--',[0 lag],[-l -l],'b--');
    xlim([-0.5 lag+0.5]);
    hold on
    stem(t,ruu,'k');
    hold off
    xlabel('Atraso');
    title('Funcao de autocorrelacao');
    sgtitle('Exercicio 4.16')

    subplot(313)
    U=fft(u');
    freq=1/(length(u))*(0:length(u)/2); 
    semilogx(freq,20*log10(abs(U(1:length(freq))))); 
    title('Espectro');
    xlim([min(freq) max(freq)])
    ylabel('|U(j\omega)| (dB)');
end
%% E 4.20
tb_list = [1, 100, 1000, 10000];
sys = tf(1, [1000 1])

for i=1:length(tb_list)
    figure(3+i)
    tb = tb_list(i)
    u=prbs(1e5, 6, tb);
    lu=length(u);
    y = lsim(sys, u, 0:lu-1);
    subplot(2, 1, 1)
    plot(1:lu,u,'-');
    axis([0 lu -0.2 1.1]);
    title(['PRBS (Tb =' num2str(tb) ')']);
    subplot(2, 1, 2)
    plot(1:lu,y,'-');
    title(['y (Tb =' num2str(tb) ')']);
    xlabel('amostras')
    xlim([-0.5 lu])
    sgtitle('Exercicio 4.20')
end