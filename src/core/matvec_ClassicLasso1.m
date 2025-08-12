function Ay = matvec_ClassicLasso1(y,par,AP)
tmp = AP*y;
Ay = par.mu * y + (AP'*tmp);