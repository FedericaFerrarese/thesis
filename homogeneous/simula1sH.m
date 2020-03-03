% MONTE CARLO ALGORITHM, MEAN FIELD SOLUTIONS AND ERROR
% NON HOMOGENEOUS ONE SPECIES MODEL 
% This code implements the Monte Carlo algorithm. The numerical solutions
% to the mean field equation are compared with the simulated ones. 
% The error is proportional to 1/sqrt(N) 
clc
clear all
close all

% parameters 
Nrange = [50,100,200,1000]; % population size

% time 
T = 500;
t = 0:T;

% transition rates 
mu = 0.5;
b = 0.5; % birth rate
c = 0.5; % competition rate
d = 0.5; % death rate

% scaled parameters
bt = mu*b;
ct = mu*c;
dt = (1-mu)*d;

% number of simulations
S = 20;
    
cc = 0; % counter 
for N=Nrange
    cc=cc+1;
    
    % initial population
    n0 = N/2;
    % intial density
    y0 = n0/N; 
    
    % mean field equation
    f = @(t,phi)2*bt*phi*(1-phi)-ct*phi^2-dt*phi;
    
    % mean field solution
    [tout,yout] = ode45(f,t,y0);
    
    for s = 1:S
        % initial random distribution of the population 
        M = zeros(1,N);
        I = randperm(N);
        p = n0;  % initial population
        M(I(1:p)) = 1;
        dens1(1,s) = sum(M)/N; % initial density

        for i=1:T
            I = randperm(N);
            J = randperm(N);

            Mmu1 = M(I(1:round(N*mu)));
            Mmu2 = M(J(1:round(N*mu)));

            M1mmu = M(I(round(N*mu)+1:N));


            % AA in AE
            p11 = sum(Mmu1==1 & Mmu2==1);
            r1 = rand(1,p11)<c;
            r = sum(r1);
            p = p-r;

            % AE in AA
            p2 = sum(Mmu1==1 & Mmu2==0);
            r2 = rand(1,p2)<b;
            r = sum(r2);
            p = p+r;

            % EA in AA
            p33 = sum(Mmu1==0 & Mmu2==1);
            r3 = rand(1,p33)<b;
            r = sum(r3);
            p = p+r;

            % A in E
            p44 = round(sum(M1mmu==1));
            r4 = rand(1,p44)<d;
            r = sum(r4);
            p = p-r;
            
            % density at time i+1
            dens1(i+1,s) = p/N; 

            M = zeros(1,N);
            M(I(1:p)) = 1;

       end
    end
     % density at each time step (mean of s simulations) 
     dens(:,cc) = mean(dens1,2);
     % error between simulated and numerical solutions 
     err(cc) = norm(dens(:,cc)-yout,inf);
end

% comparison between numerical and simulated solutions 
figure
for i=1:cc
    hold on 
    plot(t,dens(:,i),'linewidth',3)
end
plot (t, yout,'r','linewidth',3)
title('Time evolution of population density')
xlabel('t')
ylabel('\phi^N,\phi')
legend({'simulations N=50',...
         'simulations N=100',...
         'simulations N=200',...
         'simulations N=1000',...
         'mean field theory'},...
'FontSize',12,'FontWeight','bold')
 set(gca,'FontSize',12,'FontWeight','bold')
axis([0 500 0 1])

% error
figure
loglog(Nrange,err,'b*-','linewidth',3)
hold on
loglog(Nrange, err(end).*(Nrange/Nrange(end)).^(-1/2),'k','linewidth',3)
xlabel('N')
ylabel('Error')
legend({'Monte Carlo error',...
'Real error as N^{(-1/2)}'},...
'FontSize',12,'FontWeight','bold')
title('Error of the Monte Carlo method')
set(gca,'FontSize',12,'FontWeight','bold')
