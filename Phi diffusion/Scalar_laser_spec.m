alpha = 5;              % Linewidth enhancement factor
kappa = 300;           % Field decay rate
gamma = 1;            % Carrier decay rate
Gamma = 1;         % Spin-flip relaxation rate

mu = 1.05;                  % Pump current

C_sp = 1*10^-5;         % Intensety of noise
R_sp = 1;

N_th = 6.25e6;    % Carrier number at threshold
N_tr = 5.935e6;        % Carrier number at transparency
M = N_tr/(N_th - N_tr);

Dt = 1e-5;          % Time step
T = 200;             % Time for solving in ns


%% Initializating variables---------------------------------------------

time = 0:Dt:(T-Dt);
L = length(time);

flag_noise = true; % Flag of noise
Q = zeros(1,L);
phi = zeros(1,L);
G = zeros(1,L);


Q(1) = Gamma/gamma/(2*kappa) * (mu - 1);

G(1) = 1;

%% ---------------------------------------------------------------------
%% Simulating of time series of the X-LP and Y-LP light intensities-----
AV = 100;
avspec = zeros(1, L);

tic;
for j=1:AV

    ksi_Q = wgn(1, L, 0);
    ksi_phi = wgn(1, L, 0);
    Fq = 0; Fp = 0;
    
    flag_noise = true;
    for i = 1:1:L-1
        if(flag_noise)
            Fq = sqrt(2*C_sp*R_sp*Q(i))*ksi_Q(i);
            Fp = sqrt(C_sp*R_sp/2/Q(i))*ksi_phi(i);
        end
        Q(i+1) = Q(i) + (2*kappa*(G(i)-1)*Q(i) + C_sp*R_sp)*Dt + Fq*sqrt(Dt);
        phi(i+1) = phi(i) + (alpha*kappa*(G(i)-1))*Dt + Fp*sqrt(Dt);
        G(i+1) = G(i) + (mu*gamma - gamma*G(i) - 2*kappa/Gamma*G(i)*Q(i))*Dt;
    end

    spec = abs(fftshift(fft(Q.*exp(1i * phi)))).^2;
    avspec = avspec + spec;
end
toc;

avspec = avspec ./ L;
freqs = 2*pi*(-L/2:L/2-1)/L/Dt;

wnd = 2000;
plot(freqs(floor(L/2-wnd):ceil(L/2+wnd)), abs(avspec(floor(L/2-wnd):ceil(L/2+wnd))));


% plot(time,Q,'k',"LineWidth",1);
% hold on;
% plot(time,phi,'b');
% plot(time,G,'r');
% hold off;