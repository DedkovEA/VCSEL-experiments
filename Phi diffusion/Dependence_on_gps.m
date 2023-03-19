alpha = 5;              % Linewidth enhancement factor
kappa = 300;           % Field decay rate
gamma = 1;            % Carrier decay rate
gamma_d = 50;         % Spin-flip relaxation rate
gamma_a = -0.1;         % Linear dichroism 
gamma_p = 2*pi*10;          % Linear birefringence

mu = 2;                  % Pump current

C_sp = 1*10^-5;         % Intensety of noise
N_th = 6.25e6;    % Carrier number at threshold
N_tr = 5.935e6;        % Carrier number at transparency
M = N_tr/(N_th - N_tr);

tic;
gps = 2*pi*logspace(-3, 1.5, 20);
Dts = logspace(-5, -6, 3);
T = 10;
stds = zeros(length(Dts), length(gps));
for ndt = 1:length(Dts)
    curDt = Dts(ndt);
    parfor ngp = 1:length(gps)
        [means, stds] = FindPsiStd(T, curDt, gamma, kappa, alpha, gamma_d, gps(ngp), gamma_a, mu, C_sp, N_th, N_tr);
        disp("gp" + gps(ngp))
        disp(stds)
        mn = mean(means);
        mnsq = mean(stds.^2 + means.^2);
        stds(ndt, ngp) = sqrt(mnsq - mn^2);
    end
end
toc;