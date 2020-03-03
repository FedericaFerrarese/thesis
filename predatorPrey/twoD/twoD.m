% TWO DIMENSIONAL PREDATORS PREYS MODEL: the stable equilibrium becomes
% unstable if we intruduce cross diffusion. This code solves the two
% dimensional problem for d=d_c where d_c is the bifurcation paramter. 
clc
clear all
close all

% transition rates 
b = 0.2;
a = 0.7;
m1 = 1; 
m2 = 0.05;
m = -1/(4*m2*m1);


% Jacobian matrix evaluated in the equilibrium (1-b,b)
J = [-b,-b;...
    a*(1-b),0];
% eigenvalues of J (the equilibrium is stable) 
eig(J)

fu = J(1,1);
fv = J(1,2);
gu = J(2,1);
gv = J(2,2);

% critical diffusion rate 
d = (m2*b^2*m-sqrt((m2*b^2*m)^2-b^2*m*(det(J)+m2^2*b^2*m)))/(b^2*m);

% diffusion matrix 
D = [m1,0;...
     d,m2];

q = m2*b-b*d;
% wave numbers 
k2 = -(m2*b-b*d)/(2*det(D));
k22 = (-q+sqrt(q^2-4*det(D)*det(J)))/(2*det(D));
k21 = (-q-sqrt(q^2-4*det(D)*det(J)))/(2*det(D));

% eigenvalues of the perturbed system (the equilibrium is unstable) 
eig(J-k2*D)

% determinant of J-k2*D 
h =@(k) k.^4*det(D)+k.^2*q+det(J);
k1=linspace(sqrt(k21)-0.2,sqrt(k22)+0.2,100);
plot(k1,h(k1),sqrt(k22),0,'+',sqrt(k21),0,'+',sqrt(k2),h(sqrt(k2)),'*')

% minimum of h in the critical wave number (negative) 
h(sqrt(k2)) 

% spatial discretization 
mx = 50;
my = 50;
ax = 0;
bx = 30;
ay = 0;
by = 30;
x = linspace(ax,bx,mx);
y1 = linspace(ay,by,my);
hx = (bx-ax)/(mx-1);
hy = (by-ay)/(my-1);
x = x';
y1 = y1';
[X,Y] = meshgrid(x,y1);

% time discretization 
mt = 1000;
at1 = 0;
bt1 = 50;
ht1 = (bt1-at1)/(mt-1);
t1 = linspace(at1,bt1,mt);


% finite differences matrices  
Cx = toeplitz(sparse([1,1],[1,2],[-2,1]/hx^2,1,mx));
Cy = toeplitz(sparse([1,1],[1,2],[-2,1]/hy^2,1,my));
idx = eye(mx);
idy = eye(my);


% periodic boundary condizions:
Cx(1,mx) = 1/hx^2;
Cx(mx,1) = 1/hx^2;
Cy(1,mx) = 1/hy^2;
Cy(my,1)=1/hy^2;

% matrix
C = kron(idy,Cx)+kron(Cy,idy);
C1 = [m1*C,sparse(zeros(mx*my));...
    d*C,m2*C];

% term without diffusion 
g = @(y) [y(1:mx*my).*(1-y(1:mx*my)-y(mx*my+1:2*mx*my));...
    a*y(mx*my+1:2*mx*my).*(y(1:mx*my)-b)];


% perturbed initial data 
y0phi = ones(mx*my,1)*b+(randn (mx*my,1))*0.01;
y0psi = ones(mx*my,1)*(1-b)+(randn (mx*my,1))*0.01;
y0phi(1) = y0phi(mx*my-1);
y0phi(mx*my) = y0phi(2); 
y0psi(1) = y0psi(mx*my-1);
y0psi(mx*my) = y0psi(2); 
y0 = [y0phi;y0psi];
yphi(:,:,1) = reshape(y0(1:mx*my),[mx,my]);
ypsi(:,:,1) = reshape(y0(mx*my+1:2*mx*my), [mx,my]);

% explicit method 
y(:,1)=y0;
for n = 1:mt-1
    y(:,n+1) = y(:,n)+ht1*(C1*y(:,n)+g(y(:,n)));
    
    yphi(:,:,n+1) = reshape(y(1:mx*my,n+1),[mx,my]);
    ypsi(:,:,n+1) = reshape(y(mx*my+1:2*mx*my,n+1),[mx,my]);
end

% preys population
figure
contourf(X,Y,yphi(:,:,n),100,'LineStyle','none')
shading flat
title('Prey population density')
xlabel('x','FontSize',12,'FontWeight','bold') 
ylabel('y','FontSize',12,'FontWeight','bold') 
axis equal
axis([ax bx ay by]) 
colormap(jet)
colorbar
set(gca,'FontSize',12,'FontWeight','bold')

% predators population
figure
contourf(X,Y,ypsi(:,:,n),100,'LineStyle','none')
shading flat
title('Predator population density')
xlabel('x','FontSize',12,'FontWeight','bold') 
ylabel('y','FontSize',12,'FontWeight','bold') 
axis equal
axis([ax bx ay by]) 
colormap(jet)
colorbar
set(gca,'FontSize',12,'FontWeight','bold')