%% Main code: exampleHEAT.m
randn('state',100)

Funf = @funfHEAT; Fung = @fungHEAT;
Tmethod = {'TEM'}; % Select the TASE-SRK methods
nTmethods = length(Tmethod);
jacup = 0; Jac = @jacHeatfix;

% Initial conditions
global X0 A B 
d=50;
tspan = [0 1];
X0 = ones(d,1);
e = ones(d,1);
delta=1/(d+1);
A = (1/(delta)^2) * spdiags([e -2*e e], [-1 0 1], d, d);
B = spdiags(ones(d,1), 0, d, d);

N = 2^15; % Number of grid intervals
M = 100; % Number of realization

inR = 9; finR = 0; % We do the time integration using 2^(inR),...,2^(finR) time grid points

for r=1:M 
       % Compute a exact solution
       [dW,dZ,W]=wiener(tspan(2),N);
       XrefT=expm((A-(1/2)*B)+B*W(end))*X0;
     
for nm = 1:nTmethods % We apply all the selected TASE-SRK
    i = 1;
    for R = 2.^(inR:-1:finR) % For the simulations done
        [XTTRK,CPUtimeTRK] = TASESRK(N, R, tspan, X0, Funf, Fung, Jac, dW, dZ, Tmethod{nm}, jacup);
        errT_TRK(i,nm,r) = norm(XTTRK-XrefT); % Error
        CPUtime_TRK(i,nm,r) = CPUtimeTRK; % CPU time
        i = i + 1;
    end
      
end
end
errT_mean    = mean(errT_TRK,    3);   
CPUtime_mean = mean(CPUtime_TRK, 3);   

% Print results
format short e
errT_Explicit = errT_mean (: ,1:4)
errT_TASE = errT_mean (: ,5:end)
format short
CPUtime_Explicit = CPUtime_mean (: ,1:4)
CPUtime_TASE = CPUtime_mean (: ,5:end)