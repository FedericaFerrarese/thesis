% INITIAL DATA: PIECEWISE CONTINUOS AND GAUSSIAN
% INPUT
% k = selected initial data 
% OUTPUT
% M1 = vector of initial random distribution of the population 
% M2 = vector of initial random distribution of the population 
% p = vector of the initial population in each cell 
function [M1,M2,p] = initialData1sNH(k)

global N mx n0 

if k == 1 % piecewise continuous initial data 
        M1 = zeros(1,N*mx);
        M2 = zeros(1,N*mx);
        p = zeros(1,mx);
        j1 = 1;
        for j = 0:N:N*mx/2-N % population is in the first half of cells
            I = randperm(N);
            J = randperm(N);
            p(j1) = n0;
            M = M1(j+1:j+N);
            Mt = M2(j+1:j+N);
            M(I(1:p(j1))) = 1;
            Mt(J(1:p(j1))) = 1;
            M1(j+1:j+N) = M;
            M2(j+1:j+N) = Mt;
            j1 = j1+1;
        end
        pn = p; % population in each cell
end

if k == 2 % Gaussian initial data 
            M1 = zeros(1,N*mx);
            M2 = zeros(1,N*mx);
            p = zeros(1,mx);

            jn = N*mx/4; % cell in which population lives
            
            I = randperm(N);
            J = randperm(N);
         
            p1n = n0;
            
            Mn = M1(jn+1:jn+N);
            Mtn = M2(jn+1:jn+N);
            Mn(I(1:p1n)) = 1;
            Mtn(J(1:p1n)) = 1;
            
            
            M1(jn+1:jn+N) = Mn;
            M2(jn+1:jn+N) = Mtn;
         
            
            p = zeros(1,mx);
   
            p(jn/N) = p1n;
            
            pn = p; % population in each cell
end
    
end