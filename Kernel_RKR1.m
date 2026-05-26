function Xnext = Kernel_RKR1(Xn,h,Funf,Fung,Winc)
        I11 = 1/2*(Winc^2-h);
        
        X10 = Xn; X11 = Xn; X30 = Xn;
        K10 = Funf(X10); K11 = Fung(X11);
        X20 = Xn + h*K10;
        X21 = Xn + h*K10 + I11/sqrt(h)*K10; 
        X31 = Xn + h*K10 - I11/sqrt(h)*K10;
        K20 = Funf(X20);
        K21 = Fung(X21); K31 = Fung(X31);
    
        Xnext = Xn + 1/2*h*(K10+K20) + Winc*K11 + sqrt(h)/2*(K21+K31);
end