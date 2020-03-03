% INITIAL DATA: PIECEWISE CONTINUOS AND GAUSSIAN
% INPUT
% k = selected initial data 
% k1 = initial concentration cell (Gaussian) 
% OUTPUT
% M1 = vector of initial random distribution of the population 
% M2 = vector of initial random distribution of the population 
% pn = vector of the initial population A in each cell 
% pm = vector of the initial population B in each cell
function [M1,M2,pn,pm] = initialData2sNH(k,k1)
global N mx n0 m0 

if k == 1 % piecewise continuos initial data   
    M1 = zeros(1,N*mx);
    M2 = zeros(1,N*mx);  
    j1 = 1;
    for j = 0:N:N*mx/2-N % population A is in the first half of cells
            I = randperm(N);
            J = randperm(N);
            p1n(j1) = n0;
            M = M1(j+1:j+N);
            Mt = M2(j+1:j+N);
            M(I(1:p1n(j1))) = 1;
            Mt(J(1:p1n(j1))) = 1;
            M1(j+1:j+N) = M;
            M2(j+1:j+N) = Mt;
            j1 = j1+1;
     end

     j1 = 1;
     for j = N*mx/2:N:mx*N-N % population B is in the second half of cells
            I = randperm(N);
            J = randperm(N);
            p1m(j1) = m0;
            M = M1(j+1:j+N);
            Mt = M2(j+1:j+N);
            M(I(1:p1m(j1))) = 2;
            Mt(J(1:p1m(j1))) = 2;
            M1(j+1:j+N) = M;
            M2(j+1:j+N) = Mt;
            j1 = j1+1;
     end
     
         pn = [p1n,zeros(1,mx/2)]; % population A in each cell
         pm = [zeros(1,mx/2),p1m]; % population B in each cell

end

if k == 2 % Gaussian initial data 
    I = randperm(N);
    J = randperm(N);
    p1n = zeros(1,mx);
    p1m = zeros(1,mx);
    
    jn = mx*N/2;
    M1 = zeros(1,N*mx);
    M2 = zeros(1,N*mx);  
    p1n(jn/N) = n0; % population A in cell jn


    M = M1(jn+1:jn+N);
    Mt = M2(jn+1:jn+N);
    M(I(1:p1n(jn/N))) = 1;
    Mt(J(1:p1n(jn/N))) = 1;
    M1(jn+1:jn+N) = M;
    M2(jn+1:jn+N) = Mt;

    j1 = jn/N-k1;

    for j = jn-k1*N:N:jn+k1*N-N % population B is in the cells nearest cell jn
            p1m(j1) = round(m0/(2*k1));
            I = randperm(N);
            J = randperm(N);
      
            M = M1(j+1:j+N);
            Mt = M2(j+1:j+N);
            M(I(1:p1m(j1))) = 2;
            Mt(J(1:p1m(j1))) = 2;
            M1(j+1:j+N) = M;
            M2(j+1:j+N) = Mt;
            j1 = j1+1;
            
    end
    pn = p1n; % population A in each cell
    pm = p1m; % population B in each cell
end
    
end