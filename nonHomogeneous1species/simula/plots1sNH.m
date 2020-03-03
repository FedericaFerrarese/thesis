% This code generates the early and late time plots 
% for the piecewise continuos and Gaussian initial data:
% N = 10, 50, 100. Parameters: b = c = d = m = 0.5.

clc
clear all
close all
n1 = 3; % early
n2 = 100; % late

% piecewise continuos
load('piecewise1sNH.mat')

% early time
% figure
Col=['k','g','b'];
for i = 1:cc
    hold on
    plot(x,densn(n1,:,i),'color',Col(i),'linewidth',3)
end
 plot(x,y(:,n1),'r','linewidth',3)
xlabel('x','FontSize',12,'FontWeight','bold') 
ylabel('\phi,\phi^N','FontSize',12,'FontWeight','bold') 
title(['Population density at time t=' num2str(n1)])
legend({'simulations N=10',...
         'simulations N=50',...
         'simulations N=100',...
         'mean field theory'},...
 'FontSize',12,'FontWeight','bold')
 set(gca,'FontSize',12,'FontWeight','bold')

axis([0 100 0 0.6])
 
% late time
figure 
for i = 1:cc
    hold on
    plot(x,densn(n2,:,i),'color',Col(i),'linewidth',3)
end
plot(x,y(:,n2),'r','linewidth',3)
xlabel('x','FontSize',12,'FontWeight','bold') 
ylabel('\phi,\phi^N','FontSize',12,'FontWeight','bold') 
title(['Population density at time t=' num2str(n2)])
legend({'simulations N=10',...
        'simulations N=50',...
         'simulations N=100',...
         'mean field theory'},...
 'FontSize',12,'FontWeight','bold')
 set(gca,'FontSize',12,'FontWeight','bold')

axis([0 100 0 0.6])



% Gaussian 
load('Gaussian1sNH.mat')

% early time
figure
for i = 1:cc
    hold on
    plot(x,densn(n1,:,i),'color',Col(i),'linewidth',3)
end
plot(x,y(:,n1),'r','linewidth',3)
xlabel('x','FontSize',12,'FontWeight','bold') 
ylabel('\phi,\phi^N','FontSize',12,'FontWeight','bold') 
title(['Population density at time t=' num2str(n1)])
legend({'simulations N=10',...
         'simulations N=50',...
         'simulations N=100',...
         'mean field theory'},...
 'FontSize',12,'FontWeight','bold')
 set(gca,'FontSize',12,'FontWeight','bold')

axis([0 100 0 0.6])
 
% late time
figure 
for i = 1:cc
    hold on
    plot(x,densn(n2,:,i),'color',Col(i),'linewidth',3)
end
plot(x,y(:,n2),'r','linewidth',3)
xlabel('x','FontSize',12,'FontWeight','bold') 
ylabel('\phi,\phi^N','FontSize',12,'FontWeight','bold') 
title(['Population density at time t=' num2str(n2)])
legend({'simulations N=10',...
         'simulations N=50',...
         'simulations N=100',...
         'mean field theory'},...
 'FontSize',12,'FontWeight','bold')
 set(gca,'FontSize',12,'FontWeight','bold')

axis([0 100 0 0.6])


