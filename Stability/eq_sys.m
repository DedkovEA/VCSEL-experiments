function out = eq_sys(in, g, alpha, gd, k, gp, ga, mu)
% in = (Qp, Qm, a=Sin(2psi))
Qp = in(1);
Qm = in(2);
s = in(3);

A = 1/(gd + (Qp+Qm)*(g+gd) + 4*g*Qp*Qm);
c = sqrt(s);

out(1) = k*(mu*A*(2*g*Qm + gd) - 1) - sqrt(Qm/Qp)*(ga*c + gp*s);
out(2) = k*(mu*A*(2*g*Qp + gd) - 1) - sqrt(Qp/Qm)*(ga*c - gp*s);
out(3) = 2*alpha*k*A*mu*g*(Qm-Qp) + (sqrt(Qp/Qm) - sqrt(Qm/Qp))*gp*c + ...
                                    (sqrt(Qp/Qm) + sqrt(Qm/Qp))*ga*s; 