function Xnext = Kernel_TEM(Xn,h,d,Funf,Fung,Winc,Zinc,Ll,Ul)
       K10 = Ul\(Ll\Funf(Xn));
       K11 = Ul\(Ll\Fung(Xn));
       
       Xnext = Xn+ h*K10 + Winc*K11;       
end