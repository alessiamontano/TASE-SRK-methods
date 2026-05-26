function [XT, CPUtime] = TASESRK (N,R,tspan,X0,Funf,Fung,Jac,Method,jacup)

% Choice of the method and TASE operator coefficients
switch Method
    case 'EM'
        KernelFun = @(Xn,h,f,g,W,Z) Kernel_EM(Xn,h,Funf,Fung,Winc);
    case 'TEM'
        p = 1; alpha = 1/2;
        KernelFun = @(Xn,h,d,f,g,W,Z,Ll,Ul) Kernel_TEM(Xn,h,d,Funf,Fung,W,Z,Ll,Ul);
    case 'RKR1'
        KernelFun = @(Xn,h,f,g,W,Z) Kernel_RKR1(Xn,h,Funf,Fung,Winc);
    case 'TRKR1'
        p = 1; alpha = 2;
        KernelFun = @(Xn,h,d,f,g,W,Z,Ll,Ul) Kernel_TRKR1(Xn,h,d,Funf,Fung,W,Z,Ll,Ul);
    case 'RKP1'
        KernelFun = @(Xn,h,d,f,g,W,Z) Kernel_RKP1(Xn,h,Funf,Fung,Winc);
    case 'TRKP1'
        p = 1; alpha = 3/2; 
        KernelFun = @(Xn,h,d,f,g,W,Z,Ll,Ul) Kernel_TRKP1(Xn,h,d,Funf,Fung,W,Z,Ll,Ul);
    case 'RKP1.5'
        KernelFun = @(Xn,h,d,f,g,W,Z) Kernel_RKP15(Xn,h,Funf,Fung,Winc,Zinc);
    case 'TRKP1.5-1'
        p = 1; alpha = 51/100;
        KernelFun = @(Xn,h,d,f,g,W,Z,Ll,Ul) Kernel_TRKP15_1(Xn,h,d,Funf,Fung,Winc,Zinc,Ll,Ul);
    case 'TRKP1.5-2'
        p = 2; alpha = [16/20 16/10];
        KernelFun = @(Xn,h,d,f,g,W,Z,Ll,Ul) Kernel_TRKP15_2(Xn,h,d,Funf,Fung,Winc,Zinc,Ll,Ul);
 end
      
%Initialization 
t=linspace(tspan(1),tspan(2),N+1); %t: discrete time grid 
dt=(tspan(2)-tspan(1))/N; %dt: constant step-size Weiner process 
dW=sqrt(dt)*randn(1,N); % dW: elementary Wiener increments
dZ=sqrt(dt)*randn(1,N); % dZ: auxiliary Wiener increments indipendent of dW
h=R*dt; %h: constant time step-size of method
d=length(X0); %d: dimension of the problem
Id=eye(d); %Id: Identity matrix of order d 
X=X0;
L=N/R;

TASEmethods = {'TEM','TRKR1','TRKP1','TRKP1.5_1','TRKP1.5_2'};

if ismember(Method, TASEmethods)
   if (jacup==0)  % If we want a 'constant Jacobian'
       C=cputime;
       Jn = Jac();
       for l = 1:p % Compute, outside the loop, the LU factorizations of Id-alpha(l)*(h*Jn)
           [Ll(1:d,(l-1)*d+1:l*d),Ul(1:d,(l-1)*d+1:l*d)] = lu(Id-alpha(l)*(h*Jn));
       end
       for n=1:L
    
           Winc = sum(dW(R*(n-1)+1:R*n));
           Zinc = sum(dZ(R*(n-1)+1:R*n));

           X(:,n+1) = KernelFun(X(:,n),h,d,Funf,Fung,Winc,Zinc,Ll,Ul);
       end
       Cf = cputime;
   elseif (jacup==1)  % If we want exact Jacobian
        C=cputime;
       for n = 1:L
           Jn = Jac(t(n),X(:,n)); % Update Jn at each step
           for l = 1:p % Compute, at each step, the LU factorizations of Id-alpha(l)*(h*Jn)
               [Ll(1:d,(l-1)*d+1:l*d),Ul(1:d,(l-1)*d+1:l*d)] = lu(Id-alpha(l)*(h*Jn));
           end

           Winc = sum(dW(R*(n-1)+1:R*n));
           Zinc = sum(dZ(R*(n-1)+1:R*n));
    
           X(:,n+1) = KernelFun(X(:,n),h,d,Funf,Fun,Winc,Zinc,Ll,Ul);
      end
        Cf = cputime;
   end

   CPUtime = Cf - C;
   XT = X(:,end);

else
    C=cputime;
       for n=1:L
    
           Winc = sum(dW(R*(n-1)+1:R*n));
           Zinc = sum(dZ(R*(n-1)+1:R*n));

           X(:,n+1) = KernelFun(X(:,n),h,Funf,Fung,Winc,Zinc);
       end
       Cf = cputime;

       CPUtime = Cf - C;
       XT = X(:,end);
end
end