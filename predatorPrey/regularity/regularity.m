% REGULARITY: this code computes the regularity of the period via
% periodogram and autocorrelation method, both for predators and preys
% density
clc
clear all
close all

% parameters
N = 12000; % populations size
n0 = N/4;  % initial predators population 
m0 = N/4;  % initial preys population 

y0 = [n0/N;m0/N]; % initial densities

% transition rates 
mu = 0.5;
b = 0.1; % preys birth 
d1 = 0.1; % predators death 
d2 = 0.; % preys death 


y1 = 0; % counter p1 

% time 
T = 10000;
t = linspace(0,T/10,T);
ht = (t(end)-t(1))/T;


for p1 = 0.1:0.05:0.6 % competition rate 
    y2 = 0; % counter p2
    y1 = y1+1;
    for p2 = 0.1:0.05:0.6 % competition rate 
        y2 = y2+1;
        rng('default')
        S = 1; % number of simulation 
        for s = 1:S
            M = zeros(1,N);
            I = randperm(N);
            pn = n0; % initial predators population 
            pm = m0; % initial preys population 
           
            M(I(1:pn)) = 1;
            M(I(pn+1:pn+pm)) = 2;

            dens1N(1,s) = pn/N; % initial predators density 
            dens1M(1,s) = pm/N; % initial preys density 

            for i = 1:T-1

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
       
        dens1N(i+1,s)=pn/N;
        dens1M(i+1,s)=pm/N;
            end

        end
        densN = mean(dens1N,2);
        densM = mean(dens1M,2);


        % REGULARITY
        % compute the signal
        xN = densN-mean(densN);
        xM = densM-mean(densM);

        % periodogram
        [pxxN,fN] = periodogram(xN,[],[],1);
        [pxxM,fM] = periodogram(xM,[],[],1);

        % find values of the periodogram over a threshold value 
        ptN = max(pxxN)/2;
        j1N = find(pxxN>ptN);
        ptM = max(pxxM)/2;
        j1M = find(pxxM>ptM);

        % autocorelation function 
        [cN,lagsN] = xcorr(xN);
        cN = cN(T:end);
        lagsN = lagsN(T:end);
        [cM,lagsM] = xcorr(xM);
        cM = cM(T:end);
        lagsM = lagsM(T:end);
        
        % verify which intervals are acceptable (predators) 
        j = 1;
        for i1 = 1:length(j1N)
           c1N(i1) = 1/fN(j1N(i1));
           if cN(round(c1N(i1)))+10 < cN(round(c1N(i1))+1)+10
                intN(j) = c1N(i1);
                j = j+1;
           end
        end
        
        % verify which intervals are acceptable (preys)
        j = 1;
        for i2 = 1:length(j1M)
           c1M(i2) = 1/fM(j1M(i2));
           if cM(round(c1M(i2)))+10 < cM(round(c1M(i2))+1)+10
                intM(j) = c1M(i2);
                j = j+1;
           end
        end
        
        % compute the regularity 
        sigmaN = std(intN);
        mN = mean(intN);
        RN(y1,y2) = mN/sigmaN;
        
        sigmaM = std(intM);
        mM = mean(intM);
        RM(y1,y2) = mM/sigmaM;
        P1(y1,y2) = p1;
        P2(y1,y2) = p2;
    end
end

% Signal (predators population)
figure
plot(xN,'b')
title('Signal (predators population)')
xlabel('t','FontSize',12,'FontWeight','bold')
ylabel('x_A','FontSize',12,'FontWeight','bold')
set(gca,'FontSize',12,'FontWeight','bold')

% Signal (preys population)
figure
plot(xM,'b')
title('Signal (preys population)','FontSize',12,'FontWeight','bold')
xlabel('t','FontSize',12,'FontWeight','bold')
ylabel('x_B','FontSize',12,'FontWeight','bold')
set(gca,'FontSize',12,'FontWeight','bold')

% Periodogram (predators population)
figure
plot(fN,pxxN,'b')
title('Periodogram (predators population)','FontSize',12,'FontWeight','bold')
xlabel('\omega','FontSize',12,'FontWeight','bold')
ylabel('P(\omega)','FontSize',12,'FontWeight','bold')
set(gca,'FontSize',12,'FontWeight','bold')
axis([0 0.05 0 2])

% Periodogram (preys population)
figure
plot(fM,pxxM,'b')
title('Periodogram (preys population)','FontSize',12,'FontWeight','bold')
xlabel('\omega','FontSize',12,'FontWeight','bold')
ylabel('P(\omega)','FontSize',12,'FontWeight','bold')
set(gca,'FontSize',12,'FontWeight','bold')
axis([0 0.05 0 2])

% Regularity values predators
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

% Regularity values preys
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

