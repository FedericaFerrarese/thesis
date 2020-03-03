% MONTE CARLO ALGORITHM, MEAN FIELD SOLUTIONS AND ERROR 
% NON HOMOGENEOUS PREDATORS PREYS MODEL 
% This code implements the Monte Carlo algorithm. The numerical solutions
% to the mean field equations are compared with the simulated ones. 

clc
clear all
close all

% parameters
N = 3200;
% transition rates 
mu = 0.5;
b = 0.1; % preys birth rate 
d1 = 0.1; % predators death rate 
d2 = 0.; % preys death rate 
p1 = 0.25; % competition rate (+1 predator -1 prey)
p2 = 0.05; % competition rate (-1 prey) 

% scaled parameters 
bt = mu*b;
d1t = (1-mu)*d1;
d2t = (1-mu)*d2;
p1t = mu*p1;
p2t = mu*p2;

r = 2*bt-d2t;
K = 1-d2t/(2*bt);


% time interval 
t = linspace(0,1200,1200);
T = length(t);


% mean field equations 
f = @(t,x)[x(1)*(2*p1t*x(2)-d1t);...
    x(2)*(r*(1-(x(2))/K)-2*(p1t+p2t+bt)*x(1))];


n0 = N/4; % initial predators population 
m0 = N/4; % initial preys population 

y0 = [n0/N;m0/N]; % initial densities 

% mean field solutions 
[tout,YOUT] = ode45(f,t,y0);

S = 10; % number of simulations 

for s = 1:S
    
    M = zeros(1,N);
    I = randperm(N);
    pn = n0; % initial predators population 
    pm = m0; % initial preys population 
    M(I(1:pn)) = 1;
    M(I(pn+1:pn+pm)) = 2;
   
    dens1N(1,s) = pn/N; % initial predators density
    dens1M(1,s) = pm/N; % initial preys density
 
    for i=1:T-1
        
        I = randperm(N);
        J = randperm(N);
        
        Mmu1 = M(I(1:round(N*mu)));
        Mmu2 = M(J(1:round(N*mu)));
        
        M1mmu = M(I(round(N*mu)+1:N));
        
        % preys B birth
        % BE in BB 
        p22 = sum(Mmu1==2 & Mmu2==0);
        r2 = rand(1,round(p22))<b;
        r = sum(r2);
        pm = pm+r;

        % EB in BB
        p33 = sum(Mmu1==0 & Mmu2==2);
        r3 = rand(1,round(p33))<b;
        r = sum(r3);
        pm = pm+r;
        
        % predators A and preys B death 
        % A in E
        p4A = round(sum(M1mmu==1));
        r4A = rand(1,p4A)<d1;
        r = sum(r4A);
        pn = pn-r;
     
        % B in E
        p44 = sum(M1mmu==2);
        r4 = rand(1,p44)<d2;
        r = sum(r4);
        pm = pm-r;
        
        % Competition
        % AB in AA
        p55 = sum(Mmu1==1 & Mmu2==2);
        r5 = rand(1,round(p55))<p1;
        r = round(sum(r5));
        pn = pn+r;
        pm = pm-r;
        
        % BA in AA
        p55 = sum(Mmu1==2 & Mmu2==1);
        r5 = rand(1,round(p55))<p1;
        r = round(sum(r5));
        pn = pn+r;
        pm = pm-r;
        % Competition
        % AB in AE
        p55 = sum(Mmu1==1 & Mmu2==2);
        r5 = rand(1,round(p55))<p2;
        r = round(sum(r5));
        pm = pm-r;
        
        % BA in EA
        p55 = sum(Mmu1==2 & Mmu2==1);
        r5 = rand(1,round(p55))<p2;
        r = round(sum(r5));
        pm = pm-r;
        
       
        I1 = randperm(N);
        M = zeros(1,N);
        M(I1(1:pn)) = 1;
        M(I1(pn+1:pn+pm)) = 2;
       
        dens1N(i+1,s) = pn/N;
        dens1M(i+1,s) = pm/N;
    end
end
% mean densities
densN = mean(dens1N,2);
densM = mean(dens1M,2);
   
% comparison between numerical and simulated solutions predators
figure
plot(t,densN,'b','linewidth',3)
hold on
plot(t,YOUT(:,1),'r','linewidth',3)
xlabel('t','FontSize',12,'FontWeight','bold') 
ylabel('\psi^N,\psi','FontSize',12,'FontWeight','bold') 
title('Evolution of prey population')
legend({'simulations N=3200',...
        'mean field theory'},...
        'FontSize',12,'FontWeight','bold')
set(gca,'FontSize',12,'FontWeight','bold')

% comparison between numerical and simulated solutions preys
figure
plot(t,densM,'b','linewidth',3)
hold on
plot(t,YOUT(:,2),'r','linewidth',3)
xlabel('t','FontSize',12,'FontWeight','bold') 
ylabel('\psi^N,\psi','FontSize',12,'FontWeight','bold') 
title('Evolution of prey population')
legend({'simulations N=3200',...
        'mean field theory'},...
        'FontSize',12,'FontWeight','bold')
set(gca,'FontSize',12,'FontWeight','bold')
