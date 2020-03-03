% This code plots the regularity values for predators and preys populations
% as the parameters p_1 and p_2 vary. 
clc
clear all
close all

load('regularity.mat')

% regularity predators and labels 
figure
contourf(P1,P2,RN,100,'LineStyle','none')
shading flat
title('Regularity values predators','FontSize',12,'FontWeight','bold')
xlabel('p1','FontSize',12,'FontWeight','bold')
ylabel('p2','FontSize',12,'FontWeight','bold')
set(gca,'FontSize',12,'FontWeight','bold')
colorbar
colormap(jet)
axis([0.1 0.6 0.1 0.6])
caxis([0 15]) 
 
 
% regularity preys 
figure
contourf(P1,P2,RM,100,'LineStyle','none')
shading flat
title('Regularity values preys','FontSize',12,'FontWeight','bold')
xlabel('p1','FontSize',12,'FontWeight','bold')
ylabel('p2','FontSize',12,'FontWeight','bold')
set(gca,'FontSize',12,'FontWeight','bold')
colorbar
colormap(jet)
axis([0.1 0.6 0.1 0.6])
% caxis([0 15])