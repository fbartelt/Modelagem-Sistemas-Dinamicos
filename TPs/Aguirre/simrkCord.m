% Simulates the Cord system in dvCord.m by performing numerical integration
% using the algorithm implemented in rkCord.m 

% LAA 15/8/18

clear; close all


t0=0; % initial time
tf=150; % final time
% integration interval
h=0.01;
t=t0:h:tf; % time vector


% initial conditions
x0=[0.1;0.1;0.1];

% initialization
x=[x0 zeros(length(x0),length(t)-1)];
% no inputs (no external force)
u=zeros(length(t),1);


for k=2:length(t)
    x(:,k)=rkCord(x(:,k-1),u(k),u(k),h,t(k));
end

figure(1)
plot(t,x(1,:));
set(gca,'FontSize',18)
xlabel('time')
ylabel('x_1')
figure(2)
plot(t,x(2,:));
set(gca,'FontSize',18)
xlabel('time')
ylabel('x_2')
figure(3)
plot3(x(1,1000:end),x(2,1000:end),x(3,1000:end));
set(gca,'FontSize',18)
xlabel('x_1')
ylabel('x_2')


