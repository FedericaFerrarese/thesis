% FINITE DIFFERENCE METHOD AND IMPLICIT/SEMIMPLICIT/EXPLICIT METHOD 
% NON HOMOGENEOUS TWO SPECIES TWO DIMENSIONAL CASE 
% This code implements the finite diffence method and implicit Euler method
% for the two dimensional two species case. 
clc
clear all
close all

% parameters
N = 100;  % population sizes 

n0 = N/4; % initial population A 
m0 = N/4; % initial population B

phi0 = n0/N; % initial density A
psi0 = m0/N; % initial density B


% transition rates 
m1 = 0.5; % migration rate population A 
m2 = 0.5; % migration rate population B 
b1 = 0.5; % birth rate population A 
b2 = 0.5; % birth rate population B 
d1 = 0.5; % death rate population A 
d2 = 0.5; % death rate population B 
c11 = 0.5; % competition rate population A 
c22 = 0.5; % competition rate population B
c12 = 0.5; % competition rate population A B
c21 = 0.; % competition rate population B A 

% spatial and time discretizations 
mx = 50; % number of spatial nodes x-axis
my = 50; % number of spatial nodes y-axis
mt = 50; % number of time nodes 

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
bt1 = 15;
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

C1 = [m1*C,sparse(zeros(mx*my));...
    sparse(zeros(mx*my)),m2*C];

% cross diffusion term
C2 = @(y)[y(1:mx*my);y(mx*my+1:2*mx*my)]'*C1*[y(mx*my+1:2*mx*my);y(1:mx*my)]-...
    [y(mx*my+1:2*mx*my);y(1:mx*my)]'*C1*[y(1:mx*my);y(mx*my+1:2*mx*my)];

% term without diffusion
g = @(y) [2*b1*y(1:mx*my).*(1-y(1:mx*my)-y(mx*my+1:2*mx*my))-c11*y(1:mx*my).^2+...
    2*c12*y(1:mx*my).*y(mx*my+1:2*mx*my)-d1*y(1:mx*my);...
    2*b2*y(mx*my+1:2*mx*my).*(1-y(1:mx*my)-y(mx*my+1:2*mx*my))-c22*y(mx*my+1:2*mx*my).^2+...
    2*c21*y(1:mx*my).*y(mx*my+1:2*mx*my)-d2*y(mx*my+1:2*mx*my)];

% derivatives 
dC2 = @(y) diag(m1*C(1:mx*my,1:mx*my)*y(mx*my+1:2*mx*my)-m1*(y(mx*my+1:2*mx*my)'*C(1:mx*my,1:mx*my))');
dC21 = @(y) diag(m1*(y(1:mx*my)'*C(1:mx*my,1:mx*my))'-m1*C(1:mx*my,1:mx*my)*y(1:mx*my));
dC22 = @(y) diag(m2*C(1:mx*my,1:mx*my)*y(1:mx*my)-m2*(y(1:mx*my)'*C(1:mx*my,1:mx*my))');
dC23 = @(y)diag(m2*(y(mx*my+1:2*mx*my)'*C(1:mx*my,1:mx*my))'-m1*C(1:mx*my,1:mx*my)*y(mx*my+1:2*mx*my));

dg = @(y)diag(2*b1*(1-2*y(1:mx*my)-y(mx*my+1:2*mx*my))-...
    2*c11*y(1:mx*my)+2*c12*y(mx*my+1:2*mx*my)-d1); % phi
dg1 = @(y)diag((-2*b1*y(1:mx*my)+2*c12*y(1:mx*my))); %psi
dg2 = @(y)diag(-2*b2*y(mx*my+1:2*mx*my)+2*c21*y(mx*my+1:2*mx*my)); % phi
dg3 = @(y)diag(2*b2*(1-y(1:mx*my)-2*y(mx*my+1:2*mx*my))-...
    2*c22*y(mx*my+1:2*mx*my)-d2+2*c21*y(1:mx*my)); %psi



dC = @(y)[dC2(y),dC21(y);dC22(y),dC23(y)];
dG = @(y)[dg(y),dg1(y);dg2(y),dg3(y)];

% fuction and derivative
f = @(y)C1*y+C2(y)+g(y);
jf = @(y)C1+dC(y)+dG(y);

% selection of the method
disp('[1] Implicit Euler method');
disp('[2]  Explicit Euler method');   
k = input('initial data 1-2: ');

% initial data
Z=[X(:) Y(:)];
sigma1n=3;
sigma1m=6;
sigma2n=1/((2*pi*phi0)*sigma1n);
sigma2m=1/((2*pi*psi0)*sigma1m);
muA = [Z(mx*my/2-mx/2,1) Z(mx*my/2-mx/2,2)];
sigmaA = [sigma1n^2 0;0 sigma2n^2];
muB = [Z(mx*my/2-mx/2,1) Z(mx*my/2-mx/2,2)];
sigmaB = [sigma1m^2 0;0 sigma2m^2];
y0phi = mvnpdf(Z,muA,sigmaA);
y0psi = mvnpdf(Z,muB,sigmaB);
y0 = [y0phi(:);y0psi(:)];

if k == 1 % implicit Euler
    F = @(y,yn,k) y-k*feval(f,y)-yn;
    JF = @(y,k) eye(length(y))-k*feval(jf,y);

    % Newton method
    tol = ht1/10; 
    y(:,1) = y0;
    yphi(:,:,1) = reshape(y0(1:mx*my),[mx,my]);
    ypsi(:,:,1) = reshape(y0(mx*my+1:2*mx*my),[mx,my]);
   
    for n = 1:mt-1
        y(:,n+1) = y(:,n);
        errest = -JF(y(:,n+1),ht1)\F(y(:,n+1),y(:,n),ht1);
        while norm(errest,inf)> tol
            y(:,n+1) = y(:,n+1)+errest;
            errest = -JF(y(:,n+1),ht1)\F(y(:,n+1),y(:,n),ht1);
        end
        y(:,n+1) = y(:,n+1)+errest;

        yphi(:,:,n+1) = reshape(y(1:mx*my,n+1),[mx,my]);
        ypsi(:,:,n+1) = reshape(y(mx*my+1:2*mx*my,n+1),[mx,my]);
       
    end
end



if k == 2 % explicit method 
    y(:,1) = y0;

for n = 1:mt-1
    y(:,n+1) = y(:,n)+ht1*(C1*y(:,n)+C2(y(:,n))+g(y(:,n)));
    yphi(:,:,n+1) = reshape(y(1:mx*my,n+1),[mx,my]);
    ypsi(:,:,n+1) = reshape(y(mx*my+1:2*mx*my,n+1),[mx,my]);
end

end

% time 
n1 = 1; % early 
n2 = mt/2; % middle 
n3 = mt; % late 

% initial time population A 
figure
contourf(X,Y,yphi(:,:,n1+1),100,'LineStyle','none')
shading flat
title(['Population A density at time t=' num2str(n1)])
xlabel('x','FontSize',12,'FontWeight','bold') 
ylabel('y','FontSize',12,'FontWeight','bold') 
axis equal
axis([0 100 0 50]) 
colormap(jet)
colorbar
set(gca,'FontSize',12,'FontWeight','bold')
caxis([0 0.3])

% initial time population B 
figure
contourf(X,Y,ypsi(:,:,n1+1),100,'LineStyle','none')
shading flat
title(['Population B density at time t=' num2str(n1)])
xlabel('x','FontSize',12,'FontWeight','bold') 
ylabel('y','FontSize',12,'FontWeight','bold') 
axis equal
axis([0 100 0 50]) 
colormap(jet)
colorbar
set(gca,'FontSize',12,'FontWeight','bold')
caxis([0 0.3])


% middle time population A 
figure
contourf(X,Y,yphi(:,:,n2+1),100,'LineStyle','none')
shading flat
title(['Population A density at time t=' num2str(2*n2)])
xlabel('x','FontSize',12,'FontWeight','bold') 
ylabel('y','FontSize',12,'FontWeight','bold') 
axis equal
axis([0 100 0 50]) 
colormap(jet)
colorbar
set(gca,'FontSize',12,'FontWeight','bold')
caxis([0 0.3])

% middle time population B 
figure
contourf(X,Y,ypsi(:,:,n2+1),100,'LineStyle','none')
shading flat
title(['Population B density at time t=' num2str(2*n2)])
xlabel('x','FontSize',12,'FontWeight','bold') 
ylabel('y','FontSize',12,'FontWeight','bold') 
axis equal
axis([0 100 0 50]) 
colormap(jet)
colorbar
set(gca,'FontSize',12,'FontWeight','bold')
caxis([0 0.3])


% final time population A 
figure
contourf(X,Y,yphi(:,:,n+1),100,'LineStyle','none')
shading flat
title(['Population A density at time t=' num2str(2*n3)])
xlabel('x','FontSize',12,'FontWeight','bold') 
ylabel('y','FontSize',12,'FontWeight','bold') 
axis equal
axis([0 100 0 50]) 
colormap(jet)
colorbar
set(gca,'FontSize',12,'FontWeight','bold')
caxis([0 0.3])

% final time population B 
figure
contourf(X,Y,ypsi(:,:,n+1),100,'LineStyle','none')
shading flat
title(['Population B density at time t=' num2str(2*n3)])
xlabel('x','FontSize',12,'FontWeight','bold') 
ylabel('y','FontSize',12,'FontWeight','bold') 
axis equal
axis([0 100 0 50])   
colormap(jet)
colorbar
set(gca,'FontSize',12,'FontWeight','bold')
caxis([0 0.3])
