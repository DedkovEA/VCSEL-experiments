%close a19*ll; 
%clc; clear;

alpha = 5;              % Linewidth enhancement factor
kappa = 300;           % Field decay rate
gamma = 1;            % Carrier decay rate
gamma_d = 500;         % Spin-flip relaxation rate
gamma_a = -0.1;         % Linear dichroism 
gamma_p = 2*pi*32;          % Linear birefringence

mu = 2;                  % Pump current

C_sp = 1*10^-5;         % Intensety of noise
N_th = 6.25e6;    % Carrier number at threshold
N_tr = 5.935e6;        % Carrier number at transparency
M = N_tr/(N_th - N_tr);

Dt = 1e-6;          % Time step for solving in ns
rnd_chunk_sz = ceil(1e6); % For each simulation process memory used will be 4 times larger

%% Initializating variables---------------------------------------------

AV = 4;                     % Number for samples to obtain average
T = 100;                      % Window
offset = 0.2;                 % Window offset

tau = 2e-4;                   % Sampling time
tauDt = ceil(tau/Dt);
tau = Dt*tauDt;                 % Real sampling time
wndfreq = pi/tau;               % Spectra will be saved only in +-wmdfreq
wnd = floor(wndfreq*T/2/pi);

L = round(T/tau);                      % Points N in spectre
if wnd >= L/2
    wnd = floor(L/2) -1;
end

wndslice = floor(L/2)+1-wnd:floor(L/2)+1+wnd;
freqs = 2*pi*(-floor(L/2):ceil(L/2-1))/T; 
freqswnd = freqs(wndslice);


tic;

specsx = zeros(AV, length(freqswnd));
specsy = zeros(AV, length(freqswnd));

% Initialization
Qp_prev = 0;
Qm_prev = 0;
phi_prev = 0;
psi_prev = 0;
G_prev = 0;
d_prev = 0;

Q = (-gamma_a + kappa*(mu-1+2*M*C_sp) + sqrt(4*(2*C_sp-1)*kappa*mu*(gamma_a+kappa) + (gamma_a+kappa*(1+2*C_sp*M+mu))*(gamma_a+kappa*(1+2*C_sp*M+mu))) )/(4*(gamma_a+kappa));

G_prev = mu/(1 + 2*Q);
Qp_prev = (1+1e-4)*Q;
Qm_prev = (1-1e-4)*Q;

Lmat = [2*kappa*(G_prev-1), -8*Q*gamma_p, 4*kappa*(C_sp+Q);
        gamma_p/2/Q, 2*gamma_a, alpha*kappa;
        -G_prev*gamma, 0, -gamma_d-2*gamma*Q];
cp = charpoly(Lmat);
hurwitz = [cp(2), cp(4), cp(2)*cp(3)-cp(4)];
if sum(hurwitz <= 0) > 0
    psi_prev = pi/2;
    Lmat2 = [2*kappa*(G_prev-1), 8*Q*gamma_p, 4*kappa*(C_sp+Q);
        -gamma_p/2/Q, -2*gamma_a, alpha*kappa;
        -G_prev*gamma, 0, -gamma_d-2*gamma*Q];
    cp2 = charpoly(Lmat2);
    hurwitz2 = [cp2(2), cp2(4), cp2(2)*cp2(3)-cp2(4)];
    if sum(hurwitz <= 0) > 0
        disp("NOT STABLE")
        disp(gps(jj))
        disp(mus(kk))
        eval = false;
        failnum = failnum + 1;
    end
else
    psi_prev = 0;
end



ns = 1;

Ex = zeros(1,L); % X-LP mode
Ey = zeros(1,L); % Y-LP mode

ei = 1;

d_from_prev = 0;
while ns <= AV
    if d_from_prev >= tauDt
        d_from_prev = 0;
        Ep = sqrt(Qp).*exp(1i.*(phi + psi));
        Em = sqrt(Qm).*exp(1i.*(phi - psi));
        Ex(ei) = (Ep + Em)/sqrt(2);
        Ey(ei) = (Ep - Em)/sqrt(2)/1i;
        ei = ei+1;
        if ei >= L
            specx = abs(fftshift(fft(Ex))).^2;
            specy = abs(fftshift(fft(Ey))).^2;
            specsx(ns,:) = specx(wndslice);
            specsy(ns,:) = specy(wndslice);
            ns = ns + 1;
            ei = ceil(offset*L) + 1;
            Ex(1:ceil(offset*L)) = Ex(end+1-ceil(offset*L):end);
            Ey(1:ceil(offset*L)) = Ey(end+1-ceil(offset*L):end);
        end
    end

    Qp = Qp_prev + 2*( kappa*(G_prev+d_prev-1)*Qp_prev + kappa*C_sp*(G_prev+d_prev+M) - sqrt(Qp_prev*Qm_prev)*(gamma_a*cos(2*psi_prev) + gamma_p*sin(2*psi_prev) ) )*Dt;
    Qm = Qm_prev + 2*( kappa*(G_prev-d_prev-1)*Qm_prev + kappa*C_sp*(G_prev-d_prev+M) - sqrt(Qp_prev*Qm_prev)*(gamma_a*cos(2*psi_prev) - gamma_p*sin(2*psi_prev) ) )*Dt;
    phi = phi_prev + ( (G_prev-1)*alpha*kappa - 1/2/sqrt(Qm_prev*Qp_prev)*( (Qp_prev+Qm_prev)*gamma_p*cos(2*psi_prev) + (Qp_prev-Qm_prev)*gamma_a*sin(2*psi_prev) ))*Dt;
    psi = psi_prev + (alpha*kappa*d_prev + 1/2/sqrt(Qm_prev*Qp_prev)*( (Qp_prev-Qm_prev)*gamma_p*cos(2*psi_prev) + (Qp_prev+Qm_prev)*gamma_a*sin(2*psi_prev) ))*Dt;
    G = G_prev + gamma*Dt*( (mu-G_prev) - G_prev*(Qp_prev + Qm_prev) - d_prev*(Qp_prev - Qm_prev) );
    d = d_prev + Dt*(-gamma_d*d_prev -gamma*G_prev*(Qp_prev - Qm_prev) - gamma*d_prev*(Qp_prev + Qm_prev) );

    if abs(G) > 100
        disp("ALERT!")
        disp(gds(jj))
        disp(mus(kk))
        break
    end

    Qp_prev = Qp;
    Qm_prev = Qm;
    phi_prev = phi;
    psi_prev = psi;
    G_prev = G;
    d_prev = d;

    d_from_prev = d_from_prev + 1;
end
specx = mean(specsx);
specy = mean(specsy);

toc;

plot(freqswnd, log10(specx));
hold on;
plot(freqswnd, log10(specy));
hold off;

save("TestVCSELspecRK.mat",'freqswnd','specx','specy','mu','kappa','alpha','gamma','gamma_d','gamma_a','gamma_p','C_sp','N_th','N_tr')