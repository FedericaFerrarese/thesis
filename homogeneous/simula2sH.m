% MONTE CARLO ALGORITHM, MEAN FIELD SOLUTIONS AND ERROR 
% NON HOMOGENEOUS TWO SPECIES MODEL 
% This code implements the Monte Carlo algorithm. The numerical solutions
% to the mean field equations are compared with the simulated ones. 
% The error is proportional to 1/sqrt(N)
clc
clear all
close all

% parameters 
 Nrange = [50,100,200,1000]; % population sizes 
 
% time 
T = 500;
t = 0:T;

% transition rates 
mu = 0.5;   
b1 = 0.5;   % birth rate population A
b2 = 0.5;   % birth rate population B
d1 = 0.5;   % death rate population A 
d2 = 0.5;   % death rate population B
c11 = 0.5;  % intraspecies competition rate A
c22 = 0.5;  % intraspecies competition rate B
c12 = 0.0;  % interspecific competition rate AB
c21 = 0.;   % interspecific competition rate BA 

% scaled parameters 
b1t = mu*b1;
b2t = mu*b2;
c11t = mu*c11;
c12t = mu*c12;
c22t = mu*c22;
c21t = mu*c21;
d1t = (1-mu)*d1;
d2t = (1-mu)*d2;

% number of simulations 
S = 20;

cc = 0; % counter 
for N = Nrange
    cc = cc+1;
    
    % initial populations
    n0 = N/2;
    m0 = N/2;
    % initial densities
    y0 = [n0/N;m0/N];
    
    % mean field equations
    f = @(t,y)[2*b1t*y(1)*(1-y(1)-y(2))-(c11t*y(1)^2+d1t*y(1)+2*c12t*y(1)*y(2));...
    2*b2t*y(2)*(1-y(1)-y(2))-(c22t*y(2)^2+d2t*y(2)+2*c21t*y(1)*y(2))];
    
    % mean field solutions 
    [TOUT,YOUT] = ode45(f,t,y0);
    
    for s = 1:S
        % initial random distribution of the populations
        M = zeros(1,N);
        I = randperm(N);
        
        % initial populations 
        pn = n0;
        pm = m0;
        M(I(1:pn)) = 1;
        M(I(pn+1:pm+pn)) = 2;
        % initial densities 
        dens1N(1,s) = pn/N;
        dens1M(1,s) = pm/N;

        for i = 1:T

            I = randperm(N);
            J = randperm(N);

            Mmu1 = M(I(1:round(N*mu)));
            Mmu2 = M(J(1:round(N*mu)));

            M1mmu = M(I(round(N*mu)+1:N));

            % population A
            
            % AA in AE
            p11 = sum(Mmu1==1 & Mmu2==1);
            r1 = rand(1,p11)<c11;
            r = sum(r1);
            pn = pn-r;

            % AE in AA
            p2 = sum(Mmu1==1 & Mmu2==0);
            r2 = rand(1,p2)<b1;
            r = sum(r2);
            pn = pn+r;

            % EA in AA
            p33 = sum(Mmu1==0 & Mmu2==1);
            r3 = rand(1,p33)<b1;
            r = sum(r3);
            pn = pn+r;

            % A in E
            p44 = round(sum(M1mmu==1));
            r4 = rand(1,p44)<d1;
            r = sum(r4);
            pn = pn-r;

            % Population B

            % BB in BE
            p11 = sum(Mmu1==2 & Mmu2==2);
            r1 = rand(1,p11)<c22;
            r = sum(r1);
            pm = pm-r;

            % BE in BB
            p2 = sum(Mmu1==2 & Mmu2==0);
            r2 = rand(1,p2)<b2;
            r = sum(r2);
            pm = pm+r;

            % EB in BB
            p33 = sum(Mmu1==0 & Mmu2==2);
            r3 = rand(1,p33)<b2;
            r = sum(r3);
            pm = pm+r;

            % B in E
            p44 = round(sum(M1mmu==2));
            r4 = rand(1,p44)<d2;
            r = sum(r4);
            pm=pm-r;

            % Competition AB/BA
            % AB in AE 
            pMN1 = sum(Mmu1==1 & Mmu2==2);
            r5 = rand(1,round(pMN1))<c21;
            r = sum(r5);
            pm = pm-r;
 
            % AB in EB 
            r51 = rand(1,pMN1)<c12;
            r = sum(r51);
            pn = pn-r;

            % BA in EB
            pMN2 = sum(Mmu1==2 & Mmu2==1);
            r6 = rand(1,round(pMN2))<c12;
            r1 = sum(r6);
            pn = pn-r1;

            % BA in AE 
            r61 = rand(1,pMN2)<c21;
            r = sum(r61);
            pm = pm-r;
            
            
            I = randperm(N);
            M = zeros(1,N);
            M(I(1:pn)) = 1;
            M(I(pn+1:pm+pn)) = 2;
            
            % densities at time i+1
            dens1N(i+1,s) = pn/N;
            dens1M(i+1,s) = pm/N;
        end
        
        % errors at each simulation
        errN1(s) = norm(dens1N(:,s)-YOUT(:,1),inf);
        errM1(s) = norm(dens1M(:,s)-YOUT(:,2),inf);
    end
    % mean densities
    densN(:,cc) = mean(dens1N,2);
    densM(:,cc) = mean(dens1M,2);
     
    % mean errors
    errN(cc) = mean(errN1);
    errM(cc) = mean(errM1);
 
end



% comparison between numerical and simulated solutions population A
figure
for i=1:cc
    hold on 
    plot(t,densN(:,i),'linewidth',3)
end
plot(t,YOUT(:,1),'b','linewidth',3)
xlabel('t','FontSize',12,'FontWeight','bold') 
ylabel('\phi^N,\phi','FontSize',12,'FontWeight','bold') 
title('Time evolution of population A density')
legend({'simulations N=50',...
        'simulations N=100',...
        'simulations N=200',...
        'simulations N=1000',...
        'mean field theory'},...
        'FontSize',12,'FontWeight','bold')
set(gca,'FontSize',12,'FontWeight','bold')

% comparison between numerical and simulated solutions population B
figure
for i=1:cc
    hold on 
    plot(t,densM(:,i),'linewidth',3)
end
plot(t,YOUT(:,2),'r','linewidth',3)
xlabel('t','FontSize',12,'FontWeight','bold') 
ylabel('\psi^N,\psi','FontSize',12,'FontWeight','bold') 
title('Time evolution of population B density')
legend({'simulations N=50',...
         'simulations N=100',...
         'simulations N=200',...
         'simulations N=1000',...
         'mean field theory'},...
         'FontSize',12,'FontWeight','bold')
set(gca,'FontSize',12,'FontWeight','bold')

 
% error population A 
figure
loglog(Nrange,errN,'b*-','linewidth',3)
hold on
loglog(Nrange,errN(end).*(Nrange/Nrange(end)).^(-1/2),'k','linewidth',3)
xlabel('N')
ylabel('Error')
legend({'Monte Carlo error',...
'Real error as N^{(-1/2)}'},...
'FontSize',12,'FontWeight','bold')
title('Error of the Monte Carlo method (population A), S=20')
set(gca,'FontSize',12,'FontWeight','bold')

% error population B
figure
loglog(Nrange,errM,'r*-','linewidth',3)
hold on
loglog(Nrange,errM(end).*(Nrange/Nrange(end)).^(-1/2),'k','linewidth',3)
xlabel('N')
ylabel('Error')
legend({'Monte Carlo error',...
'Real error as N^{(-1/2)}'},...
'FontSize',12,'FontWeight','bold')
title('Error of the Monte Carlo method (population B), S=20')
set(gca,'FontSize',12,'FontWeight','bold')