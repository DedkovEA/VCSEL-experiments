function [p, D, x] = hurwitz(x0, g, alpha, gd, k, gp, ga, mu)
%HURWITZ Calculate coefficients of characteristic polynomial and Hurwits
%determinants
%   p - charpoly coeffs
%   D - Hurwits determinants up to D5
%   x - stationary solution
%   Note that in case at least one p is negative, this function returns
%   D=zeros(1,5).

%% Solving stationary problem
eqs = @(x)eq_sys(x, g, alpha, gd, k, gp, ga, mu);
options = optimoptions('fsolve','Display','off');
options.FunctionTolerance = 1e-14;
options.StepTolerance = 1e-14;

x = fsolve(eqs,x0, options);

Qp = x(1);
Qm = x(2);
s2p = x(3);
c2p = sqrt(1 - s2p*s2p);
A = 1/(gd + (Qp+Qm)*(g+gd) + 4*g*Qp*Qm);
N = A*mu*(g*(Qp + Qm) + gd);
n = A*mu*g*(Qm - Qp);
% w = sqrt(Qm/Qp)*( (alpha*ga - gp)*c2p + (ga + alpha*gp)*s2p );
% w2 = sqrt(Qp/Qm)*( (alpha*ga - gp)*c2p - (ga + alpha*gp)*s2p ); % the same


%% Defining matrix of linearize system
% It's better to use Qp Qm and s/c2p than N and n
M = zeros(6);

% 4x4 matrix for a_pm
M(1,3) = -ga*c2p - gp*s2p;      M(2,4) = M(1,3);
M(1,1) = -sqrt(Qm/Qp)*M(1,3);   M(2,2) = M(1,1);
M(1,4) = gp*c2p - ga*s2p;       M(2,3) = -M(1,4);
M(1,2) = -sqrt(Qm/Qp)*M(1,4);   M(2,1) = -M(1,2);
M(3,1) = -ga*c2p + gp*s2p;      M(4,2) = M(3,1);
M(3,3) = -sqrt(Qp/Qm)*M(3,1);   M(4,4) = M(3,3);
M(4,1) = -gp*c2p - ga*s2p;      M(3,2) = -M(4,1);
M(3,4) = sqrt(Qp/Qm)*M(4,1);    M(4,3) = -M(3,4);

% 2x2 matrix for N,n
M(5,5) = -g*(Qp + Qm) - g;
M(6,6) = -g*(Qp + Qm) - gd;
M(5,6) = g*(Qm - Qp);           M(6,5) = M(5,6);

%rest matrix
M(1,5) = k*sqrt(Qp);            M(1,6) = M(1,5);
M(2,5) = alpha*k*sqrt(Qp);      M(2,6) = M(2,5);
M(3,5) = k*sqrt(Qm);            M(3,6) = -M(3,5);
M(4,5) = alpha*k*sqrt(Qm);      M(4,6) = -M(4,5);
M(5,1) = -2*g*sqrt(Qp)*(N+n);   M(6,1) = M(5,1);
M(5,3) = -2*g*sqrt(Qm)*(N-n);   M(6,3) = -M(5,3);


%% Performing stability analysis
% Here p(i) corresponds to a_{i-1} 
p=fastCharPoly(M);
D = zeros(1,5);
if all(p >= 0) 
    D(1) = p(2);
    D(2) = p(2)*p(3) - p(1)*p(4);
    D(3) = p(4)*D(2) - p(2)*p(2)*p(5) + p(6)*p(1)*p(2);
    D(4) = p(5)*D(3) + p(1)*p(6)*(p(3)*p(4) - p(1)*p(6)) - ...
           p(2)*(p(3)*p(3)*p(6) - p(1)*p(5)*p(6) - p(7)*D(2));
    D(5) = p(6)*D(4) + p(7)*(p(1)*p(4)*p(4)*p(4) - ...
           p(2)*p(4) *(p(3)*p(4) + 2*p(1)*p(6)) + ...
           p(2)*p(2)*(p(4)*p(5) + p(3)*p(6)) - p(2)*p(2)*p(2)*p(7));
end

end

