% This code generates the early and late time plots 
% for the piecewise and Gaussian initial data:
% N = 12,48,100, S = 5, m1 = 0.5, m2 = 0.5, b1 = 0.5, b2 =
% 0.5, d1 = 0.5, d2 = 0.5, c11 = 0.5, c22 = 0.5, c12 = 0, c21 = 0.5. 

clc
clear all
close all

n1 = 3; % early
n2 = 100; % late

% piecewise continuos
load('piecewise2sNH.mat')

% early time
Col = ['g','c','b'];
Col1 = ['m','k','r']; 
figure
for i = 1:cc
    hold on
    plot(x,densn(n1,:,i),'color',Col(i),'linewidth',3)
    plot(x,densm(n1,:,i),'color',Col1(i),'linewidth',3)
end
plot(x,y(1:mx,n1),'b-*','linewidth',3)
plot(x,y(1+mx:2*mx,n1),'r-*','linewidth',3)
xlabel('x','FontSize',12,'FontWeight','bold') 
ylabel('\phi,\psi','FontSize',12,'FontWeight','bold') 
title(['Population densities at time t=' num2str(n1)])
legend({'simulations A N=12',...
        'simulations B N=12',...
        'simulations A N=48',...
        'simulations B N=48',...
        'simulations A N=100',...
        'simulations B N=100',...
        'mean field theory A',...
        'mean field theory B'},...
 'FontSize',12,'FontWeight','bold')
 set(gca,'FontSize',12,'FontWeight','bold')

axis([0 100 0 0.9])
 
% late time
figure 
for i = 1:cc
    hold on
    plot(x,densn(n2,:,i),'color',Col(i),'linewidth',3)
    plot(x,densm(n2,:,i),'color',Col1(i),'linewidth',3)
end
plot(x,y(1:mx,n2),'b-*','linewidth',3)
hold on
plot(x,y(1+mx:2*mx,n2),'r-*','linewidth',3)

xlabel('x','FontSize',12,'FontWeight','bold') 
ylabel('\phi,\psi','FontSize',12,'FontWeight','bold') 
title(['Population densities at time t=' num2str(n2)])
legend({'simulations A N=12',...
        'simulations B N=12',...
        'simulations A N=48',...
        'simulations B N=48',...
        'simulations A N=100',...
        'simulations B N=100',...
        'mean field theory A',...
        'mean field theory B'},...
 'FontSize',12,'FontWeight','bold')
 set(gca,'FontSize',12,'FontWeight','bold')

axis([0 100 0 0.6])

% Gaussian 
load('Gaussian2sNH.mat')

% early time 
figure
for i = 1:cc
    hold on
    plot(x,densn(n1,:,i),'color',Col(i),'linewidth',3)
    plot(x,densm(n1,:,i),'color',Col1(i),'linewidth',3)
end
plot(x,y(1:mx,n1),'b-*','linewidth',3)
plot(x,y(1+mx:2*mx,n1),'r-*','linewidth',3)
xlabel('x','FontSize',12,'FontWeight','bold') 
ylabel('\phi,\psi','FontSize',12,'FontWeight','bold') 
title(['Population densities at time t=' num2str(n1)])
legend({'simulations A N=12',...
        'simulations B N=12',...
        'simulations A N=48',...
        'simulations B N=48',...
        'simulations A N=100',...
        'simulations B N=100',...
        'mean field theory A',...
        'mean field theory B'},...
 'FontSize',12,'FontWeight','bold')
 set(gca,'FontSize',12,'FontWeight','bold')

axis([0 100 0 0.9])
 
% late time
figure 
for i = 1:cc
    hold on
    plot(x,densn(n2,:,i),'color',Col(i),'linewidth',3)
    plot(x,densm(n2,:,i),'color',Col1(i),'linewidth',3)
end
plot(x,y(1:mx,n2),'b-*','linewidth',3)
hold on
plot(x,y(1+mx:2*mx,n2),'r-*','linewidth',3)

xlabel('x','FontSize',12,'FontWeight','bold') 
ylabel('\phi,\psi','FontSize',12,'FontWeight','bold') 
title(['Population densities at time t=' num2str(n2)])
legend({'simulations A N=12',...
        'simulations B N=12',...
        'simulations A N=48',...
        'simulations B N=48',...
        'simulations A N=100',...
        'simulations B N=100',...
        'mean field theory A',...
        'mean field theory B'},...
 'FontSize',12,'FontWeight','bold')
 set(gca,'FontSize',12,'FontWeight','bold')

axis([0 100 0 0.6])

