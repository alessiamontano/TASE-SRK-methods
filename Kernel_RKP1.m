function Xnext = Kernel_RKP1(Xn,h,Funf,Fung,Winc)
        I11 = 1/2*(Winc^2-h);
        
        K10 = Funf(Xn); K11 = Fung(Xn);
        X = Xn + h*K10 + sqrt(h)*K11;
        K21 = Fung(X);

        Xnest = Xn + h*K10 + Winc*K11 + I11/sqrt(h)*(K21-K11);
end