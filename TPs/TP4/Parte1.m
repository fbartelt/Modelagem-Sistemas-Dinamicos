%Felipe Bartelt de Assis Pessoa - 2016026841
% Script adaptado dos algoritmos de Aguirre, L. A. 
clear 
close all

%% Q 5.18
u = prbs(500, 6, 1);
y(1) = 0;
y(2) = 0;
ye(1) = 0;
ye(2) = 0;

for k = 3:length(u)
   y(k) = 1.5 * y(k - 1) - 0.7 * y(k-2) + u(k - 1) + 0.5 * u(k-2);
end
for k = 3:length(u)
   ye(k) = 1.5 * ye(k - 1) - 0.7 * ye(k-2) + u(k - 1) + 0.5 * u(k-2) + 0.05 * std(y) * randn();
end

figure(1)
plot(1:500, y, '--', 'LineWidth', 1.8)
hold on
plot(1:500, ye, '-.', 'LineWidth', 1.5)
hold off
title('Resposta do sistema a entrada PRBS')
legend('Sistema sem ruido', 'Sistema com ruido')
xlabel('Amostras')
ylabel('Amplitude')

% MQ 4 x 4

Y1 = y(5:8).';
Y2 = y(205:208).';
Psi1 = [flip(y(3:4)) flip(u(3:4)); 
        flip(y(4:5)) flip(u(4:5)); 
        flip(y(5:6)) flip(u(5:6)); 
        flip(y(6:7)) flip(u(6:7))]
Psi2 = [flip(y(203:204)) flip(u(203:204)); 
        flip(y(204:205)) flip(u(204:205)); 
        flip(y(205:206)) flip(u(205:206)); 
        flip(y(206:207)) flip(u(206:207))]
theta1 = pinv(Psi1) * Y1
theta2 = pinv(Psi2) * Y2
%xi1 = Y1 - Psi1 * theta1
%xi2 = Y2 - Psi2 * theta2

Y1e = ye(5:8).';
Y2e = ye(205:208).';
Psi1e = [flip(ye(3:4)) flip(u(3:4)); 
        flip(ye(4:5)) flip(u(4:5)); 
        flip(ye(5:6)) flip(u(5:6)); 
        flip(ye(6:7)) flip(u(6:7))]
Psi2e = [flip(ye(203:204)) flip(u(203:204)); 
        flip(ye(204:205)) flip(u(204:205)); 
        flip(ye(205:206)) flip(u(205:206)); 
        flip(ye(206:207)) flip(u(206:207))]
theta1e = pinv(Psi1e) * Y1e
theta2e = pinv(Psi2e) * Y2e
%xi1e = Y1e - Psi1e * theta1e
%xi2e = Y2e - Psi2e * theta2e

% MQ 498 x 4

for k=1:length(u)-2
   PSI(k,:) = [flip(y(k:k+1)) flip(u(k:k+1))];
   PSIe(k, :) = [flip(ye(k:k+1)) flip(u(k:k+1))];
end
Y = y(3:end)';
Ye = ye(3:end)';
Theta = pinv(PSI) * Y
Thetae = pinv(PSIe) * Ye
