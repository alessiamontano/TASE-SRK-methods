function Xnext = Kernel_TRKR1(Xn,h,Funf,Fung,Winc,Ll,Ul)
        I11 = 1/2*(Winc^2-h);
        
        X10 = Xn; X11 = Xn; X30 = Xn;
        K10 = Ul\(Ll\Funf(X10)); K11 = Ul\(Ll\Fung(X11));
        X20 = Xn + h*K10;
        X21 = Xn + h*K10 + I11/sqrt(h)*K10; 
        X31 = Xn + h*K10 - I11/sqrt(h)*K10;
        K20 = Ul\(Ll\Funf(X20));
        K21 = Ul\(Ll\Fung(X21)); K31 = Ul\(Ll\Fung(X31));
    
        Xnext = Xn + 1/2*h*(K10+K20) + Winc*K11 + sqrt(h)/2*(K21+K31);
end