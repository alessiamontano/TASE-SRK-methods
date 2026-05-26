function Xnext = Kernel_TRKP15_2(Xn,h,d,Funf,Fung,Winc,Zinc,Ll,Ul)
        I11 = 1/2*(Winc^2-h); 
        I10 = h/2*(Winc+1/sqrt(3)*Zinc);
        I01 = Winc*h-I10;
        I111 = 1/6*(Winc^3-3*h*Winc);

        K110 = 2*(Ul(1:d,1:d)\(Ll(1:d,1:d)\Funf(Xn))); 
        K120 = -Ul(1:d,d+1:2*d)\(Ll(1:d,d+1:2*d)\Funf(Xn));
        K111 = 2*(Ul(1:d,1:d)\(Ll(1:d,1:d)\Fung(Xn))); 
        K121 = -Ul(1:d,d+1:2*d)\(Ll(1:d,d+1:2*d)\Fung(Xn));
        X1 = Xn + h*K10 + sqrt(h)*K11;
        X2 = Xn + h*K10 - sqrt(h)*K11;
        K210 = 2*(Ul(1:d,1:d)\(Ll(1:d,1:d)\Funf(X2))); 
        K220 = -Ul(1:d,d+1:2*d)\(Ll(1:d,d+1:2*d)\Funf(X2));
        K211 = 2*(Ul(1:d,1:d)\(Ll(1:d,1:d)\Fung(X2))); 
        K221 = -Ul(1:d,d+1:2*d)\(Ll(1:d,d+1:2*d)\Fung(X2));
        K310 = 2*(Ul(1:d,1:d)\(Ll(1:d,1:d)\Funf(X2))); 
        K320 = -Ul(1:d,d+1:2*d)\(Ll(1:d,d+1:2*d)\Funf(X2));
        K311 = 2*(Ul(1:d,1:d)\(Ll(1:d,1:d)\Fung(X2))); 
        K321 = -Ul(1:d,d+1:2*d)\(Ll(1:d,d+1:2*d)\Fung(X2));
        X3 = X1 + sqrt(h)*K21;
        X4 = X1 - sqrt(h)*K21;
        K411 = 2*(Ul(1:d,1:d)\(Ll(1:d,1:d)\Fung(X3))); 
        K421 = -Ul(1:d,d+1:2*d)\(Ll(1:d,d+1:2*d)\Fung(X3));
        K511 = 2*(Ul(1:d,1:d)\(Ll(1:d,1:d)\Fung(X4))); 
        K521 = -Ul(1:d,d+1:2*d)\(Ll(1:d,d+1:2*d)\Fung(X4));

        Xnest = Xn + h*(K110+K120) + Winc*(K111+K121) + I10/(2*sqrt(h))*(K210+K220-K310-K320) + I11/(2*sqrt(h))*(K211+K221-K311-K321) + h/4*(K210+K220-2*K110-2*K120+K310+K320) + I01/(2*h)*(K211+K221-2*K111-2*K121-K311+K321) + I111/(2*h)*(K411+K421-K511-K521-K211-K221+K311+K321);
end 