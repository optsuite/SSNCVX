function Ay = matvec_ClassicLasso_Amap1(y,par,opts)
tmp = -opts.ATmap(y)*par.DD.^2;
tmp = (1 - par.rr).*tmp;
Ay = par.mu * y + opts.Amap(tmp);