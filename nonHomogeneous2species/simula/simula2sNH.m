% MONTE CARLO ALGORITHM AND COMPARISON WITH NUMERICAL SOLUTION
% NON HOMOGENEOUS 2 SPECIES 1 DIMENSION
% This code implements the Monte Carlo algorithm. Simulated and numerical
% solutions to the mean field solutions are compared.

clc
clear all
close all 
% parameters
global  N n0 m0 q1 q2 m1 m2 b1 b2 d1 d2 c11 c12 c22 c21 mx

Nrange = [12,48,100];

q1 = 1/3;
q2 = 1/3; 

% transition rates  
m1 = 0.5;
m2 = 0.5;
b1 = 0.5;
b2 = 0.5;
d1 = 0.5;
d2 = 0.5;
c11 = 0.5;
c22 = 0.5;
c12 = 0.;
c21 = 0.5;


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
    
    n0 = N/4; % initial population A in each cell 
    m0 = 3*N/4; % initial population B in each cell 
    phi0 = n0/N; % initial density A in each cell 
    psi0 = m0/N;  % initial density B in each cell
    
    k1 = 30; 
    
    for s = 1:S    
        % initial population distribution 
        [M1,M2,pn,pm] = initialData2sNH(k,k1);

        % initial densities 
        dens1n(1,:,s) = pn/N;
        dens1m(1,:,s) = pm/N;
        
        for i = 1:mt-1 % time
            j1 = 1;

            for  j = 0:N:N*mx-N % cells 
                
                % population in each cell
                [pn,pm] = sol2sNH(M1,M2,pn,pm,j,j1);
                j1 = j1+1;
                M1 = zeros(1,N*mx);
                M2 = zeros(1,N*mx);


                j2=1;
                % distribution of populations in cells
               for jk = 0:N:N*mx-N
                    I = randperm(N);
                    J = randperm(N);
                    M = M1(jk+1:jk+N);
                    Mt = M2(jk+1:jk+N);
                    M(I(1:pn(j2))) = 1;
                    M(I(pn(j2)+1:pn(j2)+pm(j2))) = 2;
                    Mt(J(1:pn(j2))) = 1;
                    Mt(J(pn(j2)+1:pn(j2)+pm(j2))) = 2;
                    M1(jk+1:jk+N) = M;
                    M2(jk+1:jk+N) = Mt;
                    j2 = j2+1;
                end        


            end 
        
            % populations in each cell at time i+1, simulation s
            dens1n(i+1,:,s) = pn/N;
            dens1m(i+1,:,s) = pm/N;
        end
      
    end
    % mean of the s simulations of populations at each time and in each cell
        densn(:,:,cc) = mean(dens1n,3);
        densm(:,:,cc) = mean(dens1m,3);
end



% finite difference matrix without bc 
C = toeplitz(sparse([1,1],[1,2],[-2,1]/hx^2,1,mx));

if k == 1 
    % Neumann boundary conditions 
    C(1,1:2) = [-2,2]/hx^2;
    C(mx,mx-1:mx) = [2,-2]/hx^2;
    
    % piecewise continuos initial data
    y0phi = phi0* (x < 50)+  0* (x >=50);
    y0psi = 0* (x < 50)+  psi0* (x >=50);
    y0 = [y0phi;y0psi];

end

if k == 2 
    % Dirichelet boundary conditions 
    C(1,1:2) = [1,0]/hx^2;
    C(mx,mx-1:mx) = [0,1]/hx^2;
    
   % Gaussian initial data
    jn = N*mx/2;
    y0phi = zeros(mx,1);
    y0psi = zeros(mx,1);
    y0phi(jn/N) = n0/N;
    y0psi(jn/N) = (m0/(2*k1))/N;
    for i = 1:k1
        y0psi(jn/N-i) = (m0/(2*k1))/N;
        y0psi(jn/N+i) = (m0/(2*k1))/N;
    end
    y0 = [y0phi;y0psi];
end

% finite differences matrix with bc
C1 = [m1*C,sparse(zeros(mx));sparse(zeros(mx)),m2*C];


% cross diffusion term 
C2 = @(y)[y(1:mx);y(mx+1:2*mx)]'.*C1*[y(mx+1:2*mx);y(1:mx)]-...
    [y(mx+1:2*mx);y(1:mx)]'.*C1*[y(1:mx);y(mx+1:2*mx)];

% reaction term 
g = @(y) [2*b1*y(1:mx).*(1-y(1:mx)-y(mx+1:2*mx))-c11*y(1:mx).^2-...
    2*c12*y(1:mx).*y(mx+1:2*mx)-d1*y(1:mx);...
    2*b2*y(mx+1:2*mx).*(1-y(1:mx)-y(mx+1:2*mx))-c22*y(mx+1:2*mx).^2-...
    2*c21*y(1:mx).*y(mx+1:2*mx)-d2*y(mx+1:2*mx)];

% function 
f = @(y)C1*y+C2(y)+g(y);

y(:,1) = y0;

figure
for n = 1:mt
    plot(x,y(1:mx,n),'b',x,y(mx+1:2*mx,n),'r',...
    x,densn(n,:,cc),'b*-',x,densm(n,:,cc),'r*-','linewidth',3)
    xlabel('x','FontSize',12,'FontWeight','bold') 
    ylabel('\phi,\psi','FontSize',12,'FontWeight','bold') 
    title(['Population density for N = 100 at time t = ' num2str(n)])
    set(gca,'FontSize',12,'FontWeight','bold')
    axis([0 100 0 0.8])
    pause(0.1)
    y(:,n+1) = y(:,n)+ht1*(f(y(:,n)));
end
