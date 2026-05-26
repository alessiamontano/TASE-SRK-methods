function Xnext = Kernel_EM(Xn, h, Funf, Fung, Winc)
      K10 = Funf(Xn);
      K11 = Fung(Xn);
      Xnext = Xn + h*K10 + Winc*K11;       
end