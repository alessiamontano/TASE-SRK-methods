function Xnext = Kernel_TRKP1(Xn,h,Funf,Fung,Winc,Ll,Ul)
        I11 = 1/2*(Winc^2-h);
        
        K10 = Ul\(Ll\Funf(Xn)); K11 = Ul\(Ll\Fung(Xn));
        X = Xn + h*K10 + sqrt(h)*K11;
        K21 = Ul\(Ll\Fung(X));

        Xnest = Xn + h*K10 + Winc*K11 + I11/sqrt(h)*(K21-K11);
end