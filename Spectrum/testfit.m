alpha = 3;              % Linewidth enhancement factor
kappa = 80;           % Field decay rate
gamma = 1;            % Carrier decay rate
gamma_d = 1000;         % Spin-flip relaxation rate
gamma_a =  2.5 ;         % Linear dichroism    #1.34495
gamma_p = 2*pi*9;          % Linear birefringence
beta = 0%0.7  ;               % Angle between birefriginces

mu = 2;                  % Pump current

C_sp = 5*10^-4;         % Intensety of noise
N_th = 6.25e6;    % Carrier number at threshold
N_tr = 5.935e6;        % Carrier number at transparency
M = N_tr/(N_th - N_tr);

Dt = 1e-6;          % Time step for solving in ns
rnd_chunk_sz = ceil(1e6); % For each simulation process memory used will be 4 times larger


mu = (N_th*3.179 - N_tr)/(N_th - N_tr);
% ang = 10^(-4.424);

f0 = 194871;
ampl = 3.1927e-18;
% W = evalspec(freqs, f0, 1, sigma, gamma, kappa, alpha, gamma_d, gamma_p, beta, gamma_a, 3.179, C_sp, N_th, N_tr);
% half = ceil(length(W)/2);
% plot(freqs, log10(W(1:half)), freqs, log10(W(half+1:end)))


datax = table2array(X450mA);
datay = table2array(Y450mA);
fqs = datax(:,1) * 1e3;
Wx = 10.^datax(:,2);
Wy = 10.^datay(:,2);
plot(fqs, log10(Wy))
[fqs, Wx] = prepareCurveData(fqs, Wx);
[fqs, Wy] = prepareCurveData(fqs, Wy);
W = [Wx; Wy];
tic
ftype = fittype('evalspec(x, f0, ampl, sigma, gamma, kappa, alpha, gamma_d, gamma_p, beta, gamma_a, mutilde, C_sp, N_th, N_tr)');
answfit = fit(fqs, W, ftype, 'StartPoint', [f0, ampl, sigma, gamma, kappa, alpha, gamma_d, gamma_p, beta, gamma_a, mutilde, C_sp, N_th, N_tr])
toc