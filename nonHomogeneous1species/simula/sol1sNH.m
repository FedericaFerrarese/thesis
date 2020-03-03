% SIMULATED SOLUTION
% INPUT
% M1 = vector of initial random distribution of the population 
% M2 = vector of initial random distribution of the population 
% p = vector of population in each cell time i 
% j = cell index 
% j1 = cell index (cell j comunicates with cells j-1,j,j+1)
% OUTPUT
% pn = vector of population in each cell time i+1
function pn = sol1sNH(M1,M2,pn,j,j1)
global N mx m b c d q1 q2
if j > 0 && j < N*mx-N % internal cells 
        Mj = M1(j+1:j+N); % cell j
        Mjp1 = M1(j+N+1:j+2*N); % cell j+1
        Mjm1 = M1(j-N+1:j); % cell j-1

        Mjc = M2(j+1:j+N); % cell j

        pjn = sum(Mj==1); % population in cell j
        pjp1n = sum(Mjp1==1); % population in cell j+1
        pjm1n = sum(Mjm1==1); % population in cell j-1


        pn(j1) = pjn; % population in cell j
        pn(j1+1) = pjp1n; % population in cell j+1
        pn(j1-1) = pjm1n; % population in cell j-1


        Mq1j = Mj(1:round(N*q1));
        Mq1jc = Mjc(1:round(N*q1));
        Mq2j = Mj(round(N*q1)+1:round(N*q1+N*q2));
        Mq3j = Mj(round(N*q1+N*q2)+1:round(N*q1+N*q2+N*q1));


        Mq2jp1 = Mjp1(round(N*q1)+1:round(N*q1+N*q2));
        Mq2jm1 = Mjm1(round(N*q1)+1:round(N*q1+N*q2));


        % AjAj in AjEj
        p11 = sum(Mq1j==1 & Mq1jc==1);
        r1 = rand(1,p11)<c;
        r = sum(r1);
        pn(j1) = pn(j1)-r;
        
        % AjEj in AjAj
        p2 = sum(Mq1j==1 & Mq1jc==0);
        r2 = rand(1,p2)<b;
        r = sum(r2);
        pn(j1) = pn(j1)+r;
        
        % EjAj in AjAj
        p33 = sum(Mq1j==0 & Mq1jc==1);
        r3 = rand(1,p33)<b;
        r = sum(r3);
        pn(j1) = pn(j1)+r;
        % Aj in Ej
        p44 = round(sum(Mq3j==1));
        r4 = rand(1,p44)<d;
        r = sum(r4);
        pn(j1)=pn(j1)-r;



        % from cell j to j+1

        % EiAj in AiEj
        p33 = sum(Mq2j==1 & Mq2jp1==0);
        r3 = rand(1,p33)<m;
        r = sum(r3);
        pn(j1) = pn(j1)-r;
        pn(j1+1) = pn(j1+1)+r;

        % from cell j to j-1

        % EiAj in AiEj
        p33 = sum(Mq2j==1 & Mq2jm1==0);
        r3 = rand(1,p33)<m;
        r = sum(r3);
        pn(j1) = pn(j1)-r;
        pn(j1-1) = pn(j1-1)+r;
end
        
if j==0 % first cell
        j1 = 1;
        Mj = M1(j+1:j+N); % cell 1
        Mjp1 = M1(j+N+1:j+2*N); % cell 2


        Mjc = M2(j+1:j+N); % cell 1

        pjn = sum(Mj==1); % population in cell 1
        pjp1n = sum(Mjp1==1); % population in cell 2

        pn(j1) = pjn; % population in cell 1
        pn(j1+1) = pjp1n; % population in cell 2

        Mq1j = Mj(1:round(N*q1));
        Mq1jc = Mjc(1:round(N*q1));
        Mq2j = Mj(round(N*q1)+1:round(N*q1+N*q2));
        Mq3j = Mj(round(N*q1+N*q2)+1:round(N*q1+N*q2+N*q1));


        Mq2jp1 = Mjp1(round(N*q1)+1:round(N*q1+N*q2));



        % AjAj in AjEj
        p11 = sum(Mq1j==1 & Mq1jc==1);
        r1 = rand(1,p11)<c;
        r = sum(r1);
        pn(j1) = pn(j1)-r;
        % AjEj in AjAj
        p2 = sum(Mq1j==1 & Mq1jc==0);
        r2 = rand(1,p2)<b;
        r = sum(r2);
        pn(j1) = pn(j1)+r;
        % EjAj in AjAj
        p33 = sum(Mq1j==0 & Mq1jc==1);
        r3 = rand(1,p33)<b;
        r = sum(r3);
        pn(j1) = pn(j1)+r;
        % Aj in Ej
        p44 = sum(Mq3j==1);
        r4 = rand(1,p44)<d;
        r = sum(r4);
        pn(j1) = pn(j1)-r;


        % from cell 1 to 2
        % EiAj in AiEj
        p33 = sum(Mq2j==1 & Mq2jp1==0);
        r3 = rand(1,p33)<m;
        r = sum(r3);
        pn(j1) = pn(j1)-r;
        pn(j1+1) = pn(j1+1)+r;

end
if j == N*mx-N % last cell
    j1=mx;
    Mj=M1(j+1:j+N); % cell N

    Mjm1=M1(j-N+1:j); % cell N-1

    Mjc=M2(j+1:j+N); % cell N


    pjn=sum(Mj==1); % population in cell N
    pjm1n=sum(Mjm1==1); % population in cell N-1

    pn(j1)=pjn; % population in cell N
    pn(j1-1)=pjm1n; % population in cell N-1

    
    Mq1j=Mj(1:round(N*q1));
    Mq1jc=Mjc(1:round(N*q1));
    Mq2j=Mj(round(N*q1)+1:round(N*q1+N*q2));
    Mq3j=Mj(round(N*q1+N*q2)+1:round(N*q1+N*q2+N*q1));
    
    Mq2jm1=Mjm1(round(N*q1)+1:round(N*q1+N*q2));


    % AjAj in AjEj
    p11 = sum(Mq1j==1 & Mq1jc==1);
    r1 = rand(1,p11)<c;
    r = sum(r1);
    pn(j1) = pn(j1)-r;
    % AjEj in AjAj
    p2 = sum(Mq1j==1 & Mq1jc==0);
    r2 = rand(1,p2)<b;
    r = sum(r2);
    pn(j1) = pn(j1)+r;
    % EjAj in AjAj
    p33 = sum(Mq1j==0 & Mq1jc==1);
    r3 = rand(1,p33)<b;
    r = sum(r3);
    pn(j1) = pn(j1)+r;
    % Aj in Ej
    p44 = round(sum(Mq3j==1));
    r4 = rand(1,p44)<d;
    r = sum(r4);
    pn(j1) = pn(j1)-r;

    % from cell N to cell N-1

    % EiAj in AiEj
    p33 = sum(Mq2j==1 & Mq2jm1==0);
    r3 = rand(1,p33)<m;
    r = sum(r3);
    pn(j1) = pn(j1)-r;
    pn(j1-1) = pn(j1-1)+r;


end
end