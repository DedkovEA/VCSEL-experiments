alpha = 5;              % Linewidth enhancement factor
kappa = 300;           % Field decay rate
gamma = 1;            % Carrier decay rate
gamma_d = 50;         % Spin-flip relaxation rate
gamma_a = -0.1;         % Linear dichroism 
gamma_p = 13.0008;          % Linear birefringence

mu = 2;                  % Pump current

C_sp = 1*10^-5;         % Intensety of noise
N_th = 6.25e6;    % Carrier number at threshold
N_tr = 5.935e6;        % Carrier number at transparency
M = N_tr/(N_th - N_tr);

Dt = 1e-5;          % Time step
T = 10;             % Time for solving in ns

tic;
[means, stds] = FindPsiStd(T, Dt, gamma, kappa, alpha, gamma_d, gamma_p, gamma_a, mu, C_sp, N_th, N_tr);
toc;
means;
stds;
mean(stds)