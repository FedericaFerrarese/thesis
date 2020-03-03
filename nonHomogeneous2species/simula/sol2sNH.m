% SIMULATED SOLUTION
% INPUT
% M1 = vector of initial random distribution of the population 
% M2 = vector of initial random distribution of the population 
% pn = vector of population A in each cell time i 
% pm = vector of population B in each cell time i 
% j = cell index 
% j1 = cell index (cell j comunicates with cells j-1,j,j+1)
% OUTPUT
% pn = vector of population A in each cell at time step i+1 
% pm = vector of population B in each cell at time step i+1
function [pn,pm] = sol2sNH(M1,M2,pn,pm,j,j1)
global N  q1 q2 m1 m2 b1 b2 d1 d2 c11 c12 c22 c21 mx

if j == 0 % cell 1 

    Mj = M1(j+1:j+N); % cell 1 
    Mjp1 = M1(j+N+1:j+2*N); % cell 2

    Mjc = M2(j+1:j+N); % cell 1 

    pjn = sum(Mj==1); % population A in cell 1
    pjp1n = sum(Mjp1==1); % population A in cell 2


    pjm = sum(Mj==2); % population B in cell 1 
    pjp1m = sum(Mjp1==2); % population B in cell 2

    pn(j1) = pjn; % population A in cell 1
    pn(j1+1) = pjp1n; % population A in cell 2

    pm(j1) = pjm; % population B in cell 1
    pm(j1+1) = pjp1m; % population B in cell 2


    Mq1j = Mj(1:round(N*q1));
    Mq1jc = Mjc(1:round(N*q1));
    Mq2j = Mj(round(N*q1)+1:round(N*q1+N*q2));
    Mq3j = Mj(round(N*q1+N*q2)+1:round(N*q1+N*q2+N*q1));


    Mq2jp1=Mjp1(round(N*q1)+1:round(N*q1+N*q2));


    % population A 
    % AjAj in AjEj
    p11 = sum(Mq1j==1 & Mq1jc==1);
    r1 = rand(1,p11)<c11;
    r = sum(r1);
    pn(j1) = pn(j1)-r;
    % AjEj in AjAj
    p2 = sum(Mq1j==1 & Mq1jc==0);
    r2 = rand(1,p2)<b1;
    r = sum(r2);
    pn(j1) = pn(j1)+r;
    % EjAj in AjAj
    p33 = sum(Mq1j==0 & Mq1jc==1);
    r3 = rand(1,p33)<b1;
    r = sum(r3);
    pn(j1) = pn(j1)+r;
    % Aj in Ej
    p44 = round(sum(Mq3j==1));
    r4 = rand(1,p44)<d1;
    r = sum(r4);
    pn(j1) = pn(j1)-r;

    % population B 
    % BjBj in BjEj
    p11 = sum(Mq1j==2 & Mq1jc==2);
    r1 = rand(1,p11)<c22;
    r = sum(r1);
    pm(j1) = pm(j1)-r;
    % BjEj in BjBj
    p2 = sum(Mq1j==2 & Mq1jc==0);
    r2 = rand(1,p2)<b2;
    r = sum(r2);
    pm(j1) = pm(j1)+r;
    % EjBj in BjBj
    p33 = sum(Mq1j==0 & Mq1jc==2);
    r3 = rand(1,p33)<b2;
    r = sum(r3);
    pm(j1) = pm(j1)+r;
    % Bj in Ej
    p44 = round(sum(Mq3j==2));
    r4 = rand(1,p44)<d2;
    r = sum(r4);
    pm(j1) = pm(j1)-r;

    % Interaction AB/BA
    % AjBj in AjEj
    p = sum(Mq1j==1 & Mq1jc==2);
    r2 = rand(1,round(p))<c21;
    r = sum(r2);
    r6 = rand(1,round(p))<c12;
    r1 = sum(r6);
    pm(j1) = pm(j1)-r;
    pn(j1) = pn(j1)-r1;

    % Migration
    % from cell 1 to cell 2 
    % EiAj in AiEj
    p33 = sum(Mq2j==1 & Mq2jp1==0);
    r3 = rand(1,round(p33))<m1;
    r = sum(r3);
    pn(j1) = pn(j1)-r;
    pn(j1+1) = pn(j1+1)+r;
    % EiAj in AiEj
    p33 = sum(Mq2j==2 & Mq2jp1==0);
    r3 = rand(1,round(p33))<m2;
    r = sum(r3);
    pm(j1) = pm(j1)-r;
    pm(j1+1) = pm(j1+1)+r;
    
end

if j>0 && j<N*mx-N % internal cells 
    
    Mj = M1(j+1:j+N); % cell j
    Mjp1 = M1(j+N+1:j+2*N); % cell j+1
    Mjm1 = M1(j-N+1:j); % cell j-1

    Mjc=M2(j+1:j+N); % cell j

    pjn = sum(Mj==1); % population A in cell j
    pjp1n = sum(Mjp1==1); % population A in cell j+1
    pjm1n = sum(Mjm1==1); % population A in cell j-1

    pjm = sum(Mj==2); % population B in cell j
    pjp1m = sum(Mjp1==2);% population B in cell j+1
    pjm1m = sum(Mjm1==2); % population B in cell j-1
    
   
    pn(j1) = pjn; % population A in cell j
    pn(j1+1) = pjp1n; % population A in cell j+1
    pn(j1-1) = pjm1n; % population A in cell j-1
    pm(j1) = pjm; % population B in cell j
    pm(j1+1) = pjp1m; % population B in cell j+1
    pm(j1-1) = pjm1m; % population B in cell j-1


    Mq1j = Mj(1:round(N*q1));
    Mq1jc = Mjc(1:round(N*q1));
    Mq2j = Mj(round(N*q1)+1:round(N*q1+N*q2));
    Mq3j = Mj(round(N*q1+N*q2)+1:round(N*q1+N*q2+N*q1));


    Mq2jp1 = Mjp1(round(N*q1)+1:round(N*q1+N*q2));
    Mq2jm1 = Mjm1(round(N*q1)+1:round(N*q1+N*q2));

    % Population A 
    % AjAj in AjEj
    p11 = sum(Mq1j==1 & Mq1jc==1);
    r1 = rand(1,p11)<c11;
    r = sum(r1);
    pn(j1) = pn(j1)-r;
    % AjEj in AjAj
    p2 = sum(Mq1j==1 & Mq1jc==0);
    r2 = rand(1,p2)<b1;
    r = sum(r2);
    pn(j1) = pn(j1)+r;
    % EjAj in AjAj
    p33 = sum(Mq1j==0 & Mq1jc==1);
    r3 = rand(1,p33)<b1;
    r = sum(r3);
    pn(j1) = pn(j1)+r;
    % Aj in Ej
    p44 = round(sum(Mq3j==1));
    r4 = rand(1,p44)<d1;
    r = sum(r4);
    pn(j1) = pn(j1)-r;

    % Population B 
    % BjBj in BjEj
    p11 = sum(Mq1j==2 & Mq1jc==2);
    r1 = rand(1,p11)<c22;
    r = sum(r1);
    pm(j1) = pm(j1)-r;
    % BjEj in BjBj
    p2 = sum(Mq1j==2 & Mq1jc==0);
    r2 = rand(1,p2)<b2;
    r = sum(r2);
    pm(j1) = pm(j1)+r;
    % EjBj in BjBj
    p33 = sum(Mq1j==0 & Mq1jc==2);
    r3 = rand(1,p33)<b2;
    r = sum(r3);
    pm(j1) = pm(j1)+r;
    % Bj in Ej
    p44 = round(sum(Mq3j==2));
    r4 = rand(1,p44)<d2;
    r = sum(r4);
    pm(j1) = pm(j1)-r;

    % Interaction AB/BA
    % AjBj in AjEj
    p = sum(Mq1j==1 & Mq1jc==2);
    r2 = rand(1,round(p))<c21;
    r = sum(r2);
    r6 = rand(1,round(p))<c12;
    r1 = sum(r6);
    pm(j1) = pm(j1)-r;
    pn(j1) = pn(j1)-r1;

    
    % Migration 
    % from cell j to cell j+1
    % EiAj in AiEj
    p33 = sum(Mq2j==1 & Mq2jp1==0);
    r3 = rand(1,round(p33))<m1;
    r = sum(r3);
    pn(j1) = pn(j1)-r;
    pn(j1+1) = pn(j1+1)+r;
    % EiAj in AiEj
    p33 = sum(Mq2j==2 & Mq2jp1==0);
    r3 = rand(1,round(p33))<m2;
    r = sum(r3);
    pm(j1) = pm(j1)-r;
    pm(j1+1) = pm(j1+1)+r;
    
    
    % from cell j to cell j-1
    % EiAj in AiEj
    p33 = sum(Mq2j==1 & Mq2jm1==0);
    r3 = rand(1,round(p33))<m1;
    r = sum(r3);
    pn(j1) = pn(j1)-r;
    pn(j1-1) = pn(j1-1)+r;
    % EiAj in AiEj
    p33 = sum(Mq2j==2 & Mq2jm1==0);
    r3 = rand(1,round(p33))<m2;
    r = sum(r3);
    pm(j1) = pm(j1)-r;
    pm(j1-1) = pm(j1-1)+r;
    
    % ridimensionalize the populations
    if pn(j1)<0
    pn(j1) = 0;
    end
    if pm(j1)<0
    pm(j1) = 0;
    end
end

if j == N*mx-N % cell N 

    Mj = M1(j+1:j+N); % cell N 
    Mjm1 = M1(j-N+1:j); % cell N-1

    Mjc = M2(j+1:j+N); % cell N 

    pjn = sum(Mj==1); % population A in cell N 
    pjm1n = sum(Mjm1==1); % population A in cell N-1

    pjm = sum(Mj==2); % population B in cell N 
    pjm1m = sum(Mjm1==2); % population B in cell N-1
    
    pn(j1) = pjn; % population A in cell N
    pn(j1-1) = pjm1n; % population A in cell N-1
    pm(j1) = pjm; % population B in cell N
    pm(j1-1) = pjm1m; % population B in cell N-1

   
    Mq1j = Mj(1:round(N*q1));
    Mq1jc = Mjc(1:round(N*q1));
    Mq2j = Mj(round(N*q1)+1:round(N*q1+N*q2));
    Mq3j = Mj(round(N*q1+N*q2)+1:round(N*q1+N*q2+N*q1));
    Mq2jm1 = Mjm1(round(N*q1)+1:round(N*q1+N*q2));

    % Population A 
    % AjAj in AjEj
    p11 = sum(Mq1j==1 & Mq1jc==1);
    r1 = rand(1,p11)<c11;
    r = sum(r1);
    pn(j1) = pn(j1)-r;
    % AjEj in AjAj
    p2 = sum(Mq1j==1 & Mq1jc==0);
    r2 = rand(1,p2)<b1;
    r = sum(r2);
    pn(j1) = pn(j1)+r;
    % EjAj in AjAj
    p33 = sum(Mq1j==0 & Mq1jc==1);
    r3 = rand(1,p33)<b1;
    r = sum(r3);
    pn(j1) = pn(j1)+r;
    % Aj in Ej
    p44 = round(sum(Mq3j==1));
    r4 = rand(1,p44)<d1;
    r = sum(r4);
    pn(j1) = pn(j1)-r;

    % Population B 
    % BjBj in BjEj
    p11 = sum(Mq1j==2 & Mq1jc==2);
    r1 = rand(1,p11)<c22;
    r = sum(r1);
    pm(j1) = pm(j1)-r;
    % BjEj in BjBj
    p2 = sum(Mq1j==2 & Mq1jc==0);
    r2 = rand(1,p2)<b2;
    r = sum(r2);
    pm(j1) = pm(j1)+r;
    % EjBj in BjBj
    p33 = sum(Mq1j==0 & Mq1jc==2);
    r3 = rand(1,p33)<b2;
    r = sum(r3);
    pm(j1) = pm(j1)+r;
    % Bj in Ej
    p44 = round(sum(Mq3j==2));
    r4 = rand(1,p44)<d2;
    r = sum(r4);
    pm(j1) = pm(j1)-r;

    % Interaction AB/BA
    % AjBj in AjEj
    p = sum(Mq1j==1 & Mq1jc==2);
    r2 = rand(1,round(p))<c21;
    r = sum(r2);
    r6 = rand(1,round(p))<c12;
    r1 = sum(r6);
    pm(j1) = pm(j1)-r;
    pn(j1) = pn(j1)-r1;


    % Migration
    % from cell N to cell N-1
    % EiAj in AiEj
    p33 = sum(Mq2j==1 & Mq2jm1==0);
    r3 = rand(1,round(p33))<m1;
    r = sum(r3);
    pn(j1) = pn(j1)-r;
    pn(j1-1) = pn(j1-1)+r;

    % EiAj in AiEj
    p33 = sum(Mq2j==2 & Mq2jm1==0);
    r3 = rand(1,round(p33))<m2;
    r = sum(r3);
    pm(j1) = pm(j1)-r;
    pm(j1-1) = pm(j1-1)+r;   

end
end