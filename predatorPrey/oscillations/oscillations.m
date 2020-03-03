% PREDATOR PREY SIGNAL WITHOUT AND WITH NOISE
clc
clear all
close all

% Parameters 
mu = 0.5;
b = 0.1;
d1 = 0.1;
d2 = 0.;
p1 = 0.6;
p2 = 0.6;

% Scaled parameters
b = b*mu;
d1 = d1*mu;
d2 = d2*mu;
p1 = p1*mu;
p2 = p2*mu;

r=2*b-d2;
K=1-d2/(2*b);

% non trivial equilibrium 
f1 = (2*b*p1-b*d1-p1*d2)/(2*p1*(p1+p2+b));
f2 = d1/(2*p1);

% Jacobian matrix 
Jh=@(f1,f2)[2*p1*f2-d1, 2*p1*f1;...
    -2*(p1+p2+b)*f2, r*(1-(2*f2/K))-2*(p1+p2+b)*f1];

% Jacobian matrix evaluated in the equilibrium 
A = Jh(f1,f2);

% covariance noise matrix 
B = [2*d1*f1,-d1*f1;-d1*f1, 2*d1*(1+p2/p1)*f1+2*d2*f2];

% Spectrum
alfa = B(1,1)*A(2,2)^2+2*B(1,2)*A(1,2)*abs(A(2,2))+B(2,2)*A(1,2)^2;
beta = B(1,1);
Omega2 = A(1,2)*abs(A(2,1));
Gamma = abs(A(2,2));

P = @(o) (alfa+beta*o.^2)./((o.^2-Omega2).^2+Gamma^2*o.^2);

% Initial data 
y0 = [0.2;0.2];

% time interval
T = 10000;
t = 1:T;
ht = (t(end)-t(1))/(length(t)-1);

% signal with and without noise 
y(:,1) = y0;
y1(:,1) = y0;
for i = 1:length(t)-1
    r = randn(2,1);
    y(:,i+1) = y(:,i)+(A*y(:,i)+B*r); % signal with noise
    y1(:,i+1) = y1(:,i)+0.5*(A*y1(:,i)); % signal without noise
end

% signal with noise 
figure
plot(t,y(1,:),'b')
title('Signal with noise')
xlabel('t','FontSize',12,'FontWeight','bold')
ylabel('x_A','FontSize',12,'FontWeight','bold')
set(gca,'FontSize',12,'FontWeight','bold')

% signal without noise 
figure
plot(t,y1(1,:),'b')
title('Signal without noise')
xlabel('t','FontSize',12,'FontWeight','bold')
ylabel('x_A','FontSize',12,'FontWeight','bold')
set(gca,'FontSize',12,'FontWeight','bold')
