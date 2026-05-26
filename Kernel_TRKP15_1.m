function Xnext = Kernel_TRKP15_1(Xn,h,Funf,Fung,Winc,Zinc,Ll,Ul)
        I11 = 1/2*(Winc^2-h); 
        I10 = h/2*(Winc+1/sqrt(3)*Zinc);
        I01 = Winc*h-I10;
        I111 = 1/6*(Winc^3-3*h*Winc);

        K10 = Ul\(Ll\Funf(Xn)); K11 = Ul\(Ll\Fung(Xn));
        X1 = Xn + h*K10 + sqrt(h)*K11;
        X2 = Xn + h*K10 - sqrt(h)*K11;
        K21 = Ul\(Ll\Fung(X1));
        X3 = X1 + sqrt(h)*K21;
        X4 = X1 - sqrt(h)*K21;
        K20 = Ul\(Ll\Funf(X1)); K30 = Ul\(Ll\Funf(X2));
        K31 = Ul\(Ll\Fung(X2)); 
        K41 = Ul\(Ll\Fung(X3)); K51 = Ul\(Ll\Fung(X4));

        Xnest = Xn + h*K10 + Winc*K11 + I10/(2*sqrt(h))*(K20-K30) + I11/(2*sqrt(h))*(K21-K31) + h/4*(K20-2*K10+K30) + I01/(2*h)*(K21-2*K11-K31) + I111/(2*h)*(K41-K51-K21+K31);
end