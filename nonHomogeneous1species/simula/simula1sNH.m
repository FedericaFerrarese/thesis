% MONTE CARLO ALGORITHM AND COMPARISON WITH NUMERICAL SOLUTION
% NON HOMOGENEOUS 1 SPECIES 1 DIMENSION
% This code implements the Monte Carlo algorithm. Simulated and numerical
% solutions to the mean field solutions are compared.

clc
clear all
close all 

% parameters
global N mx m b c d n0 q1 q2

Nrange = [10,50,100]; % population size 
q1 = 1/3;
q2 = 1/3; 

% transition rates  
b = 0.5; % birth 
d = 0.5; % death
c = 0.5; % competition
m = 0.5; % migration

% spatial discretization
mx = 100;
ax = 0;
bx = 100;
hx = (bx-ax) /(mx-1); % dx
x = linspace(ax,bx,mx);
x = x';

% time discretization
mt = 100;
at1 = 0;
bt1 = 30;
ht1 = (bt1-at1)/(mt-1);%dt
t1 = linspace(at1,bt1,mt);

% number of simulations
S = 5;

% selection of the initial data 
disp('[1] piecewise continuous initial data'); 
disp('[2]  Gaussian initial data'); 
k = input('initial data 1-2: ');

cc = 0; % counter 
for N = Nrange
    cc = cc+1;
    
    n0 = N/2; % initial population in each cell 
    phi0 = n0/N; % initial density in each cell 

    for s=1:S
        
        % initial population distribution 
        [M1,M2,pn] = initialData1sNH(k);

        % initial density 
        dens1n(1,:,s) = pn/N;

        for i = 1:mt-1 % time
            j1 = 1;

            for  j = 0:N:N*mx-N % cells 
                
                % population in each cell
                pn = sol1sNH(M1,M2,pn,j,j1);
                j1=j1+1;
                M1=zeros(1,N*mx);
                M2=zeros(1,N*mx);


                j2=1;
                % distribution of population in cells
                for jk=0:N:N*mx-N 
                    I=randperm(N);
                    J=randperm(N);
                    M=M1(jk+1:jk+N);
                    Mt=M2(jk+1:jk+N);
                    M(I(1:pn(j2)))=1;

                    Mt(J(1:pn(j2)))=1;

                    M1(jk+1:jk+N)=M;
                    M2(jk+1:jk+N)=Mt;
                    j2=j2+1;
                end      


            end 
        
            % population in each cell at time i+1, simulation s
            dens1n(i+1,:,s)=pn/N;

        end

    end
    % mean of the s simulations of population at each time and in each cell
        densn(:,:,cc)=mean(dens1n,3);
end


% finite differences method with explicit Euler

% finite differences matrix 
C = toeplitz(sparse([1,1],[1,2],[-2,1]/hx^2,1,mx));


% reaction term
g = @(y) 2*b*y.*(1-y)-c*y.^2-d*y;

if k == 1 
    % Neumann boudary conditions
    C(1,1:2) = [-2,2]/hx^2;
    C(mx,mx-1:mx) = [2,-2]/hx^2;
 
    % piecewise continuous initial condition
    y0 = phi0* (x < 50)+  0* (x >=50);
end

if k == 2 
    % Dirichelet boudary conditions
    C(1,1:2) = [1,0]/hx^2;
    C(mx,mx-1:mx) = [0,1]/hx^2;
    
    % Gaussian initial condition
    y0 = zeros(mx,1);
    jn = mx/4; 
    y0(jn) = n0/N;

end
C1 = C*m;

% Explicit Euler
y(:,1) = y0;

figure
for n = 1:mt
    % comparison between simulated and numerical solutions
    plot(x,densn(n,:,cc),'b',x,y(:,n),'r','linewidth',3)
    xlabel('x','FontSize',12,'FontWeight','bold') 
    ylabel('\phi,\phi^N','FontSize',12,'FontWeight','bold') 
    title(['Population density for N=100 at time t=' num2str(n)])
    set(gca,'FontSize',12,'FontWeight','bold')
    axis([0 100 0 0.6])
    pause(0.1)
    y(:,n+1) = y(:,n)+ht1*(C1*y(:,n)+g(y(:,n)));
end

