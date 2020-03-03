% FINITE DIFFERENCE METHOD AND IMPLICIT/SEMIMPLICIT/EXPLICIT METHOD 
% NON HOMOGENEOUS ONE SPECIES TWO DIMENSIONAL CASE 
% This code implements the finite diffence method and implicit Euler method
% for the two dimensional one species case. 
clc
clear all
close all

% parameters
N = 120; % population size 

n0 = N/4; % initial population
phi0 = n0/N; % initial density 

% parameters 
mu = 0.5;

b = 0.5; % birth 
d = 0.5; % death 
c = 0.5; % competition
m = 0.5; % migration 

% scaled parameters 
b = b*mu; 
d = d*mu; 
c = c*mu;
m = m*mu; 

% spatial and time discretizations 
mx = 50; % number of spatial nodes x-axis
my = 50; % number of spatial nodes y-axis
mt = 100; % number of time nodes 

% spatial interval 
ax = 0; 
bx = 100;
hx = (bx-ax) /(mx-1); % dx
x = linspace(ax,bx,mx);
x = x';
ay = 0; 
by = 50;
hy = (by-ay) /(my-1); % dy
y1 = linspace(ay,by,my);
y1 = y1';

[X,Y]=ndgrid(x,y1);

% time interval
at1 = 0;
bt1 = 30;
ht1 = (bt1-at1)/(mt-1); % dt
t1 = linspace(at1,bt1,mt);


% finite differences
Cx = toeplitz(sparse([1,1],[1,2],[-2,1]/hx^2,1,mx));
Cy = toeplitz(sparse([1,1],[1,2],[-2,1]/hy^2,1,mx));
idx = eye(mx);
idy = eye(my);

% Dirichelet boundary condition 
Cx(1,1:2) = [1,0]/hx^2;
Cx(mx,mx-1:mx) = [0,1]/hx^2;
Cy(1,1:2) = [1,0]/hy^2;
Cy(my,my-1:my) = [0,1]/hy^2;

C = kron(idy,Cx)+kron(Cy,idy);
C1 = m*C;


% reaction term3
g = @(y) 2*b*y.*(1-y)-c*y.^2-d*y;
dg = @(y) 2*b-(4*b+2*c)*y-d; 

% Gaussian initial data  
Z=[X(:) Y(:)];
mu = [50 25];
sigma = [5 0; 0 5];
y0phi = 10*mvnpdf(Z,mu,sigma);
y0=y0phi(:);

% selection of the method
disp('[1] Implicit Euler method');
disp('[2]  Explicit Euler method');   
k = input('initial data 1-2: ');

f = @(y)C1*y+g(y);
jf = @(y)C1+diag(dg(y));


if k == 1 % IMPLICIT METHOD
    % implicit Euler
    F = @(y,yn,k) y-k*feval(f,y)-yn;
    JF = @(y,k) eye(length(y))-k*feval(jf,y);
    tol = ht1/100; 
    y(:,1) = y0;
    yphi(:,:,1) = reshape(y0,[mx,my]);
    % Newton method
    for n = 1:mt-1
        y(:,n+1) = y(:,n);

        errest = -JF(y(:,n+1),ht1)\F(y(:,n+1),y(:,n),ht1);
        while norm(errest,inf)> tol
            y(:,n+1) = y(:,n+1)+errest;
            errest = -JF(y(:,n+1),ht1)\F(y(:,n+1),y(:,n),ht1);
        end
        y(:,n+1) = y(:,n+1)+errest;

        yphi(:,:,n+1) = reshape(y(:,n+1),[mx,my]);
    end

end


if k == 2 % EXPLICIT METHOD

    y(:,1) = y0;
    yphi(:,:,1) = reshape(y0,[mx,my]);
    for n = 1:mt-1
        y(:,n+1) = y(:,n)+ht1*f(y(:,n));
        yphi(:,:,n+1) = reshape(y(:,n+1),[mx,my]);
    end
end

n1 = 0; % early time
n2 = mt/2-1; % middle time
n = mt-1; % late time

% early time
figure
contourf(X,Y,yphi(:,:,n1+1),100,'LineStyle','none')
shading flat
title(['Population density at time t=' num2str(n1+1)])
xlabel('x','FontSize',12,'FontWeight','bold') 
ylabel('y','FontSize',12,'FontWeight','bold') 
axis equal
axis([0 100 0 50])   
colormap(jet)
colorbar
set(gca,'FontSize',12,'FontWeight','bold')
caxis([0 0.3])


% middle time 
figure
contourf(X,Y,yphi(:,:,n2+1),100,'LineStyle','none')
shading flat
title(['Population densitiy at time t=' num2str(n2+1)])
xlabel('x','FontSize',12,'FontWeight','bold') 
ylabel('y','FontSize',12,'FontWeight','bold') 
axis equal
axis([0 100 0 50])   
colormap(jet)
colorbar
set(gca,'FontSize',12,'FontWeight','bold')
caxis([0 0.3])

% late time 
figure
contourf(X,Y,yphi(:,:,n+1),100,'LineStyle','none')
shading flat
title(['Population density at time t=' num2str(n+1)])
xlabel('x','FontSize',12,'FontWeight','bold') 
ylabel('y','FontSize',12,'FontWeight','bold') 
axis equal
axis([0 100 0 50])   
colormap(jet)
colorbar
set(gca,'FontSize',12,'FontWeight','bold')
caxis([0 0.3])
