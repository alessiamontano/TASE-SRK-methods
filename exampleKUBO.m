%% Main code: exampleKUBO.m
randn('state',100)
Funf = @funfKUBO; Fung = @fungKUBO;
Tmethod = {'EM','RKR1','RKP1','RKP15','TEM','TRKR1','TRKP1','TRKP15-1','TRKP15-2'}; % Select the methods
nTmethods = length(Tmethod);
jacup = 0; Jac = @jacKUBOfix;

% Initial conditions
global  a b
a = 10; b = 0.3;
tspan = [0 10];
X0 = [1; 0];
N = 2^15; % Number of grid intervals
M = 1000; % Number of realization

inR = 4; finR = 0; % We do the time integration using 2^(inR),...,2^(finR) time grid points

for r=1:M 
       % Compute a exact solution
       dW=sqrt(tspan(2)/N)*randn(1,N); W=cumsum(dW);
       XrefT=[cos(a*tspan(2)+b*W(end));sin(a*tspan(2)+b*W(end))];
     
for nm = 1:nTmethods % We apply all the selected TASE-SRK
    i = 1;
    for R = 2.^(inR:-1:finR) % For the simulations done
        [XTTRK,CPUtimeTRK] = TASESRK(N, R, tspan, X0, Funf, Fung, Jac, dW, dZ, Tmethod{nm} , jacup);
        errT_TRK(i,nm,r) = norm(XTTRK-XrefT); % Error
        i = i + 1;
    end    
end
end
errT_mean    = mean(errT_TRK,    3);                
pest = diff(-log10(errT_mean),1,1) ./ log10(2); % Compute the estimated order

% Print results
format short e
errT_Explicit = errT_mean(:,1:4)
errT_TASE = errT_mean(:,5:end)
format short
pest_TASE = pest(:,5:end)