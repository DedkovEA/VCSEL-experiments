alpha = 3;              % Linewidth enhancement factor
kappa = 80;           % Field decay rate
gamma = 1;            % Carrier decay rate
gamma_d = 1000;         % Spin-flip relaxation rate
gamma_a =  2.5 ;         % Linear dichroism    #1.34495
gamma_p = 2*pi*9;          % Linear birefringence
beta = 0%0.7  ;               % Angle between birefriginces

mu = 6;                  % Pump current

C_sp = 5*10^-4;         % Intensety of noise
N_th = 6.25e6;    % Carrier number at threshold
N_tr = 5.935e6;        % Carrier number at transparency
M = N_tr/(N_th - N_tr);

Dt = 1e-6;          % Time step for solving in ns
rnd_chunk_sz = ceil(1e6); % For each simulation process memory used will be 4 times larger


%mu = (N_th*3.179 - N_tr)/(N_th - N_tr);
% ang = 10^(-4.424);

%% Initializating variables---------------------------------------------

AV = 1;                     % Number for samples to obtain average
T = 10;                      % Window
offset = 0.2;                 % Window offset

tau = 1e-4;                   % Sampling time


% init
M = N_tr/(N_th - N_tr);

tauDt = round(tau/Dt);
tau = Dt*tauDt;                 % Real sampling time
wndfreq = 1/2/tau;               % Spectra will be saved only in +-wndfreq
wnd = round(wndfreq*T);

L = round(T/tau);                      % Points N in spectre
% if wnd >= L/2
%     wnd = floor(L/2) -1;
% end

wndslice = round(L/2)+1-wnd:round(L/2)+wnd;
freqs = (-floor(L/2):ceil(L/2-1))/T; 
freqswnd = freqs(wndslice);



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
Qp_prev = Q;
Qm_prev = Q;

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
    if sum(hurwitz2 <= 0) > 0
        disp("NOT STABLE")
    end
else
    psi_prev = 0;
end



ns = 1;

Ex = zeros(1,L); % X-LP mode
Ey = zeros(1,L); % Y-LP mode
QQpp = zeros(1,L);
QQmm = zeros(1,L);
pphi = zeros(1,L);
ppsi = zeros(1,L);
GG = zeros(1,L);
dd = zeros(1,L);


ei = 1;
ni = 1;
cntr = 1;
% ksi_plus = normrnd(0, 1, 1, rnd_chunk_sz);
% ksi_minus = normrnd(0, 1, 1, rnd_chunk_sz);
% ksi_phi = normrnd(0, 1, 1, rnd_chunk_sz);
% ksi_psi = normrnd(0, 1, 1, rnd_chunk_sz);
ksi_plus = SDErandom(1:rnd_chunk_sz, 1);
ksi_minus = SDErandom(1:rnd_chunk_sz, 2);
ksi_phi = SDErandom(1:rnd_chunk_sz, 3);
ksi_psi = SDErandom(1:rnd_chunk_sz, 4);

d_from_prev = 0;
while ns <= AV

    if ni > rnd_chunk_sz
        ksi_plus = SDErandom(cntr*rnd_chunk_sz+1:(cntr+1)*rnd_chunk_sz, 1);
        ksi_minus = SDErandom(cntr*rnd_chunk_sz+1:(cntr+1)*rnd_chunk_sz, 2);
        ksi_phi = SDErandom(cntr*rnd_chunk_sz+1:(cntr+1)*rnd_chunk_sz, 3);
        ksi_psi = SDErandom(cntr*rnd_chunk_sz+1:(cntr+1)*rnd_chunk_sz, 4);
        ni = 1;
        cntr = cntr + 1;
    end

    Fp = 2*sqrt(C_sp*kappa*Qp_prev*(G_prev+d_prev+M))*ksi_plus(ni);
    Fm = 2*sqrt(C_sp*kappa*Qm_prev*(G_prev-d_prev+M))*ksi_minus(ni);
    Fphi = 1/2*sqrt(C_sp*(G_prev+d_prev+M)/Qp_prev/2)*(ksi_phi(ni)+ksi_psi(ni)) + 1/2*sqrt(C_sp*(G_prev-d_prev+M)/Qm_prev/2)*(ksi_phi(ni)-ksi_psi(ni));
    Fpsi = 1/2*sqrt(C_sp*(G_prev+d_prev+M)/Qp_prev/2)*(ksi_phi(ni)+ksi_psi(ni)) - 1/2*sqrt(C_sp*(G_prev-d_prev+M)/Qm_prev/2)*(ksi_phi(ni)-ksi_psi(ni));

    Qp = Qp_prev + 2*( kappa*(G_prev+d_prev-1)*Qp_prev + kappa*C_sp*(G_prev+d_prev+M) - sqrt(Qp_prev*Qm_prev)*(gamma_a*cos(2*psi_prev-2*beta) + gamma_p*sin(2*psi_prev) ) )*Dt + Fp*sqrt(Dt);
    Qm = Qm_prev + 2*( kappa*(G_prev-d_prev-1)*Qm_prev + kappa*C_sp*(G_prev-d_prev+M) - sqrt(Qp_prev*Qm_prev)*(gamma_a*cos(2*psi_prev-2*beta) - gamma_p*sin(2*psi_prev) ) )*Dt + Fm*sqrt(Dt);
    phi = phi_prev + ( (G_prev-1)*alpha*kappa - 1/2/sqrt(Qm_prev*Qp_prev)*( (Qp_prev+Qm_prev)*gamma_p*cos(2*psi_prev) + (Qp_prev-Qm_prev)*gamma_a*sin(2*psi_prev-2*beta) ))*Dt + Fphi*sqrt(Dt);
    psi = psi_prev + (alpha*kappa*d_prev + 1/2/sqrt(Qm_prev*Qp_prev)*( (Qp_prev-Qm_prev)*gamma_p*cos(2*psi_prev) + (Qp_prev+Qm_prev)*gamma_a*sin(2*psi_prev-2*beta) ))*Dt + Fpsi*sqrt(Dt);
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

    ni = ni + 1;
    d_from_prev = d_from_prev + 1;

        if d_from_prev >= tauDt
        d_from_prev = 0;
        Ep = sqrt(Qp).*exp(1i.*(phi + psi));
        Em = sqrt(Qm).*exp(1i.*(phi - psi));
         Ex(ei) = (Ep + Em)/sqrt(2);
         Ey(ei) = (Ep - Em)/sqrt(2)/1i;
        % NEED Additional rotation on beta
        % Ex(ei) = cos(ang)*(Ep + Em)/sqrt(2) + sin(ang)*(Ep - Em)/sqrt(2)/1i;
        % Ey(ei) = -sin(ang)*(Ep + Em)/sqrt(2) + cos(ang)*(Ep - Em)/sqrt(2)/1i;

        QQpp(ei) = Qp;
        QQmm(ei) = Qm;
        ppsi(ei) = psi;
        pphi(ei) = phi;
        GG(ei) = G;
        dd(ei) = d;

        ei = ei+1;
        if ei > L
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
end
specx = mean(specsx, 1);
specy = mean(specsy, 1);

plot(freqs, log10(specx), freqs, log10(specy));