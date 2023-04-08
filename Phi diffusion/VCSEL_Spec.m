%close all; 
%clc; clear;

alpha = 5;              % Linewidth enhancement factor
kappa = 300;           % Field decay rate
gamma = 1;            % Carrier decay rate
gamma_d = 1500;         % Spin-flip relaxation rate
gamma_a = -0.1;         % Linear dichroism 
gamma_p = 2*pi*32;          % Linear birefringence

mu = 2;                  % Pump current

C_sp = 1*10^-5;         % Intensety of noise
N_th = 6.25e6;    % Carrier number at threshold
N_tr = 5.935e6;        % Carrier number at transparency
M = N_tr/(N_th - N_tr);

Dt = 1e-5;          % Time step for solving in ns
rnd_chunk_sz = ceil(1e6); % For each simulation process memory used will be 4 times larger


%% Initializating variables---------------------------------------------

AV = 20;                     % Number for samples to obtain average
T = 100;                      % Window
offset = 0.5;                 % Window offset

tau = 1e-3;                   % Sampling time
tauDt = ceil(tau/Dt);
tau = Dt*tauDt;                 % Real sampling time
wndfreq = 2000;               % Spectra will be saved only in +-wmdfreq
wnd = ceil(wndfreq*T/2/pi);

L = round(T/tau);                      % Points N in spectre

wndslice = floor(L/2)+1-wnd:floor(L/2)+1+wnd;
freqs = 2*pi*(-floor(L/2):ceil(L/2-1))/T; 
freqswnd = freqs(wndslice);


specsx = zeros(AV, length(freqswnd));
specsy = zeros(AV, length(freqswnd));

tic
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
else
    psi_prev = 0;
end

ns = 1;

Ex = zeros(1,L); % X-LP mode
Ey = zeros(1,L); % Y-LP mode

ei = 1;
ni = 1;
ksi_plus = normrnd(0, 1, 1, rnd_chunk_sz);
ksi_minus = normrnd(0, 1, 1, rnd_chunk_sz);
ksi_phi = normrnd(0, 1, 1, rnd_chunk_sz);
ksi_psi = normrnd(0, 1, 1, rnd_chunk_sz);

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

    if ni > rnd_chunk_sz
        ksi_plus = normrnd(0, 1, 1, rnd_chunk_sz);
        ksi_minus = normrnd(0, 1, 1, rnd_chunk_sz);
        ksi_phi = normrnd(0, 1, 1, rnd_chunk_sz);
        ksi_psi = normrnd(0, 1, 1, rnd_chunk_sz);
        ni = 1;
    end

    Fp = 2*sqrt(C_sp*kappa*Qp_prev*(G_prev+d_prev+M))*ksi_plus(ni);
    Fm = 2*sqrt(C_sp*kappa*Qm_prev*(G_prev-d_prev+M))*ksi_minus(ni);
    Fphi = 1/2*sqrt(C_sp*(G_prev+d_prev+M)/Qp_prev/2)*(ksi_phi(ni)+ksi_psi(ni)) + 1/2*sqrt(C_sp*(G_prev-d_prev+M)/Qm_prev/2)*(ksi_phi(ni)-ksi_psi(ni));
    Fpsi = 1/2*sqrt(C_sp*(G_prev+d_prev+M)/Qp_prev/2)*(ksi_phi(ni)+ksi_psi(ni)) - 1/2*sqrt(C_sp*(G_prev-d_prev+M)/Qm_prev/2)*(ksi_phi(ni)-ksi_psi(ni));

    Qp = Qp_prev + 2*( kappa*(G_prev+d_prev-1)*Qp_prev + kappa*C_sp*(G_prev+d_prev+M) - sqrt(Qp_prev*Qm_prev)*(gamma_a*cos(2*psi_prev) + gamma_p*sin(2*psi_prev) ) )*Dt + Fp*sqrt(Dt);
    Qm = Qm_prev + 2*( kappa*(G_prev-d_prev-1)*Qm_prev + kappa*C_sp*(G_prev-d_prev+M) - sqrt(Qp_prev*Qm_prev)*(gamma_a*cos(2*psi_prev) - gamma_p*sin(2*psi_prev) ) )*Dt + Fm*sqrt(Dt);
    phi = phi_prev + ( (G_prev-1)*alpha*kappa - 1/2/sqrt(Qm_prev*Qp_prev)*( (Qp_prev+Qm_prev)*gamma_p*cos(2*psi_prev) + (Qp_prev-Qm_prev)*gamma_a*sin(2*psi_prev) ))*Dt + Fphi*sqrt(Dt);
    psi = psi_prev + (alpha*kappa*d_prev + 1/2/sqrt(Qm_prev*Qp_prev)*( (Qp_prev-Qm_prev)*gamma_p*cos(2*psi_prev) + (Qp_prev+Qm_prev)*gamma_a*sin(2*psi_prev) ))*Dt + Fpsi*sqrt(Dt);
    G = G_prev + gamma*Dt*( (mu-G_prev) - G_prev*(Qp_prev + Qm_prev) - d_prev*(Qp_prev - Qm_prev) );
    d = d_prev + Dt*(-gamma_d*d_prev -gamma*G_prev*(Qp_prev - Qm_prev) - gamma*d_prev*(Qp_prev + Qm_prev) );
    
%     if Qp < 0 || Qm < 0
%         nd = 100;
%         ksi_plus_tmp = normrnd(0, 1, 1, nd);
%         ksi_minus_tmp = normrnd(0, 1, 1, nd);
%         ksi_phi_tmp = normrnd(0, 1, 1, nd);
%         ksi_psi_tmp = normrnd(0, 1, 1, nd);
%         for k = 1:nd
%             Fp = 2*sqrt(C_sp*kappa*Qp_prev*(G_prev+d_prev+M))*ksi_plus_tmp(k);
%             Fm = 2*sqrt(C_sp*kappa*Qm_prev*(G_prev-d_prev+M))*ksi_minus_tmp(k);
%             Fphi = 1/2*sqrt(C_sp*(G_prev+d_prev+M)/Qp_prev/2)*(ksi_phi_tmp(k)+ksi_psi_tmp(k)) + 1/2*sqrt(C_sp*(G_prev-d_prev+M)/Qm_prev/2)*(ksi_phi_tmp(k)-ksi_psi_tmp(k));
%             Fpsi = 1/2*sqrt(C_sp*(G_prev+d_prev+M)/Qp_prev/2)*(ksi_phi_tmp(k)+ksi_psi_tmp(k)) - 1/2*sqrt(C_sp*(G_prev-d_prev+M)/Qm_prev/2)*(ksi_phi_tmp(k)-ksi_psi_tmp(k));
%         
%             Qp = Qp_prev + 2*( kappa*(G_prev+d_prev-1)*Qp_prev + kappa*C_sp*(G_prev+d_prev+M) - sqrt(Qp_prev*Qm_prev)*(gamma_a*cos(2*psi_prev) + gamma_p*sin(2*psi_prev) ) )*Dt/nd + Fp*sqrt(Dt/nd);
%             Qm = Qm_prev + 2*( kappa*(G_prev-d_prev-1)*Qm_prev + kappa*C_sp*(G_prev-d_prev+M) - sqrt(Qp_prev*Qm_prev)*(gamma_a*cos(2*psi_prev) - gamma_p*sin(2*psi_prev) ) )*Dt/nd + Fm*sqrt(Dt/nd);
%             phi = phi_prev + ( (G_prev-1)*alpha*kappa - 1/2/sqrt(Qm_prev*Qp_prev)*( (Qp_prev+Qm_prev)*gamma_p*cos(2*psi_prev) + (Qp_prev-Qm_prev)*gamma_a*sin(2*psi_prev) ))*Dt/nd + Fphi*sqrt(Dt/nd);
%             psi = psi_prev + (alpha*kappa*d_prev + 1/2/sqrt(Qm_prev*Qp_prev)*( (Qp_prev-Qm_prev)*gamma_p*cos(2*psi_prev) + (Qp_prev+Qm_prev)*gamma_a*sin(2*psi_prev) ))*Dt/nd + Fpsi*sqrt(Dt/nd);
%             G = G_prev + gamma*Dt/nd*( (mu-G_prev) - G_prev*(Qp_prev + Qm_prev) - d_prev*(Qp_prev - Qm_prev) );
%             d = d_prev + Dt/nd*(-gamma_d*d_prev -gamma*G_prev*(Qp_prev - Qm_prev) - gamma*d_prev*(Qp_prev + Qm_prev) );
%             
%             Qp_prev = Qp;
%             Qm_prev = Qm;
%             phi_prev = phi;
%             psi_prev = psi;
%             G_prev = G;
%             d_prev = d;
%         end
%     end

    if abs(G) > 100
        disp("ALERT!")
        break
    end
    if imag(Qp) > 1e-17
        disp("ALERT!")
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
end
toc
  
avspecx = mean(specsx);
avspecy = mean(specsy);






% tic
% for j = 1:AV
%     
%     Qp = zeros(1,L); % X-LP mode
%     Qm = zeros(1,L); % Y-LP mode
%     psi = zeros(1,L);
%     phi = zeros(1,L);
%     G = zeros(1,L);  % Total carrier inversion
%     d = zeros(1,L);  % The difference between carrier inversions with opposite spins
%     
%     Q = 1/2*(kappa*mu/(gamma_a + kappa) - 1);
%     
%     Q = (-gamma_a + kappa*(mu-1+2*M*C_sp) + sqrt(4*(2*C_sp-1)*kappa*mu*(gamma_a+kappa) + (gamma_a+kappa*(1+2*C_sp*M+mu))*(gamma_a+kappa*(1+2*C_sp*M+mu))) )/(4*(gamma_a+kappa));
%     
%     G(1) = mu/(1 + 2*Q);
%     Qp(1) = Q;
%     Qm(1) = Q;
%     
%     
%     Lmat = [2*kappa*(G(1)-1), -8*Q*gamma_p, 4*kappa*(C_sp+Q);
%             gamma_p/2/Q, 2*gamma_a, alpha*kappa;
%             -G(1)*gamma, 0, -gamma_d-2*gamma*Q];
%     cp = charpoly(Lmat);
%     hurwitz = [cp(2), cp(4), cp(2)*cp(3)-cp(4)];
%     if sum(hurwitz <= 0) > 0
%         psi(1) = pi/2;
%     end
%     
%     %% ---------------------------------------------------------------------
%     %% Simulating of time series of the X-LP and Y-LP light intensities-----
%     
% %     ksi_plus = wgn(1, L, 0);
% %     ksi_minus = wgn(1, L, 0);
% %     ksi_phi = wgn(1, L, 0);
% %     ksi_psi = wgn(1, L, 0);
%     ksi_plus = normrnd(0, 1, 1,L);
%     ksi_minus = normrnd(0, 1, 1,L);
%     ksi_phi = normrnd(0, 1, 1,L);
%     ksi_psi = normrnd(0, 1, 1,L);
%     Fp = 0; Fm = 0; Fpsi = 0; Fphi = 0;
%     
%     flag_noise = true;
%     for i = 1:1:L-1
%         if(flag_noise)
%             Fp = 2*sqrt(C_sp*kappa*Qp(i)*(G(i)+d(i)+M))*ksi_plus(i);
%             Fm = 2*sqrt(C_sp*kappa*Qm(i)*(G(i)-d(i)+M))*ksi_minus(i);
%             % Fpsi = 1/2*sqrt(C_sp*( (G(i)+d(i)+M)/Qp(i) + (G(i)-d(i)+M)/Qm(i) ))*ksi_psi(i);
%             Fphi = 1/2*sqrt(C_sp*(G(i)+d(i)+M)/Qp(i)/2)*(ksi_phi(i)+ksi_psi(i)) + 1/2*sqrt(C_sp*(G(i)-d(i)+M)/Qm(i)/2)*(ksi_phi(i)-ksi_psi(i));
%             Fpsi = 1/2*sqrt(C_sp*(G(i)+d(i)+M)/Qp(i)/2)*(ksi_phi(i)+ksi_psi(i)) - 1/2*sqrt(C_sp*(G(i)-d(i)+M)/Qm(i)/2)*(ksi_phi(i)-ksi_psi(i));
%         end
%         Qp(i+1) = Qp(i) + 2*( kappa*(G(i)+d(i)-1)*Qp(i) + kappa*C_sp*(G(i)+d(i)+M) - sqrt(Qp(i)*Qm(i))*(gamma_a*cos(2*psi(i)) + gamma_p*sin(2*psi(i)) ) )*Dt + Fp*sqrt(Dt);
%         Qm(i+1) = Qm(i) + 2*( kappa*(G(i)-d(i)-1)*Qm(i) + kappa*C_sp*(G(i)-d(i)+M) - sqrt(Qp(i)*Qm(i))*(gamma_a*cos(2*psi(i)) - gamma_p*sin(2*psi(i)) ) )*Dt + Fm*sqrt(Dt);
%         phi(i+1) = phi(i) + ( (G(i)-1)*alpha*kappa - 1/2/sqrt(Qm(i)*Qp(i))*( (Qp(i)+Qm(i))*gamma_p*cos(2*psi(i)) + (Qp(i)-Qm(i))*gamma_a*sin(2*psi(i)) ))*Dt + Fphi*sqrt(Dt);
%         psi(i+1) = psi(i) + (alpha*kappa*d(i) + 1/2/sqrt(Qm(i)*Qp(i))*( (Qp(i)-Qm(i))*gamma_p*cos(2*psi(i)) + (Qp(i)+Qm(i))*gamma_a*sin(2*psi(i)) ))*Dt + Fpsi*sqrt(Dt);
%         G(i+1) = G(i) + gamma*Dt*( (mu-G(i)) - G(i)*(Qp(i) + Qm(i)) - d(i)*(Qp(i) - Qm(i)) );
%         d(i+1) = d(i) + Dt*(-gamma_d*d(i) -gamma*G(i)*(Qp(i) - Qm(i)) - gamma*d(i)*(Qp(i) + Qm(i)) );
%         if abs(G(i+1)) > 1000
%             disp("ALERT!")
%         end
%     end
%     
%     Ep = sqrt(Qp).*exp(1i.*(phi + psi));
%     Em = sqrt(Qm).*exp(1i.*(phi - psi));
%     Ex = (Ep + Em)/sqrt(2);
%     Ey = (Ep - Em)/sqrt(2)/1i;
%     plot
%     specx = abs(fftshift(fft(Ex))).^2;
%     specy = abs(fftshift(fft(Ey))).^2;
%     specsx(j,:) = specx(floor(L/2-wnd):floor(L/2-wnd)+ceil(2*wnd-1));
%     specsy(j,:) = specy(floor(L/2-wnd):floor(L/2-wnd)+ceil(2*wnd-1));
% end
% toc;
% avspecx = mean(specsx);
% avspecy = mean(specsy);



plot(freqswnd, log(abs(avspecx)));
hold on
plot(freqswnd, log(abs(avspecy)));
hold off


% save("VCSELspec1.mat",'specsx','specsy','mu','kappa','alpha','gamma','gamma_d','gamma_a','gamma_p','C_sp','N_th','N_tr')


% plot(freqs(floor(L/2-wnd):ceil(L/2+wnd)), abs(avspecx(floor(L/2-wnd):ceil(L/2+wnd))));
% hold on
% plot(freqs(floor(L/2-wnd):ceil(L/2+wnd)), abs(avspecy(floor(L/2-wnd):ceil(L/2+wnd))));
% hold off

%% Evaluating on different gps

% T = 200;
% Ngps = 12;
% gps = 2*pi*linspace(100,1000,Ngps);
% devs = zeros(1, Ngps);
% tic;
% parfor j = 1:Ngps
%     %Initializating variables---------------------------------------------
%     flag_noise = true; % Flag of noise
%     Ep = zeros(1,L); % X-LP mode
%     Em = zeros(1,L); % Y-LP mode
%     G = zeros(1,L);  % Total carrier inversion
%     d = zeros(1,L);  % The difference between carrier inversions with opposite spins
%     
%     G(1)=1;
%     gamma_p = gps(j);
%     ksi_plus = wgn(1, L, 0) + 1i*wgn(1, L, 0);
%     ksi_minus = wgn(1, L, 0) + 1i*wgn(1, L, 0);
%     Fp = 0; Fm = 0;
%     
%     for i = 1:1:L-1
%         if(flag_noise)
%             Fp = sqrt(C_sp*kappa*(G(i)+d(i)+M))*ksi_plus(i);
%             Fm =sqrt(C_sp*kappa*(G(i)-d(i)+M))*ksi_minus(i);
%         end
%         Ep(i+1) = Ep(i) + (kappa*(1+1i*alpha)*(G(i)+d(i)-1)*Ep(i) - (gamma_a + 1i*gamma_p)*Em(i))*Dt + Fp*sqrt(Dt);
%         Em(i+1) = Em(i) + (kappa*(1+1i*alpha)*(G(i)-d(i)-1)*Em(i) - (gamma_a + 1i*gamma_p)*Ep(i))*Dt + Fm*sqrt(Dt);
%         G(i+1) = G(i) + gamma*Dt*( (mu-G(i)) - G(i)*(conj(Ep(i))*Ep(i) + conj(Em(i))*Em(i)) - d(i)*(conj(Ep(i))*Ep(i) - conj(Em(i))*Em(i)) );
%         d(i+1) = d(i) + Dt*(-gamma_d*d(i) -gamma*G(i)*(conj(Ep(i))*Ep(i) - conj(Em(i))*Em(i)) - gamma*d(i)*(conj(Ep(i))*Ep(i) + conj(Em(i))*Em(i)) );
%     end
% 
%     psi = 0.5*angle(Em./Ep);
%     devs(j) = std(psi(ceil(0.02*L):end));
% end
% toc;
% 
% wolfStr = "{";
% for i = 1:Ngps
%     wolfStr = wolfStr + "{" + gps(i) + ", " + devs(i) + "}, ";
% end
% wolfStr = wolfStr + "}";
% disp(wolfStr)

%% Psi stats
% tic;
% Npsi = 200;
% psi_devs = zeros(1,Npsi);
% parfor j = 1:Npsi
%     %Initializating variables---------------------------------------------
%     flag_noise = true; % Flag of noise
%     Ep = zeros(1,L); % X-LP mode
%     Em = zeros(1,L); % Y-LP mode
%     G = zeros(1,L);  % Total carrier inversion
%     d = zeros(1,L);  % The difference between carrier inversions with opposite spins
%     
%     G(1)=1;
%     
%     %---------------------------------------------------------------------
%     %Simulating of time series of the X-LP and Y-LP light intensities-----
%     ksi_plus = wgn(1, L, 1) + 1i*wgn(1, L, 1);
%     ksi_minus = wgn(1, L, 1) + 1i*wgn(1, L, 1);
%     Fp = 0; Fm = 0;
%     
%     for i = 1:1:L-1
%         if(flag_noise)
%             Fp = sqrt(C_sp*kappa*(G(i)+d(i)+M))*ksi_plus(i);
%             Fm =sqrt(C_sp*kappa*(G(i)-d(i)+M))*ksi_minus(i);
%         end
%         Ep(i+1) = Ep(i) + (kappa*(1+1i*alpha)*(G(i)+d(i)-1)*Ep(i) - (gamma_a + 1i*gamma_p)*Em(i))*Dt + Fp*sqrt(Dt);
%         Em(i+1) = Em(i) + (kappa*(1+1i*alpha)*(G(i)-d(i)-1)*Em(i) - (gamma_a + 1i*gamma_p)*Ep(i))*Dt + Fm*sqrt(Dt);
%         G(i+1) = G(i) + gamma*Dt*( (mu-G(i)) - G(i)*(conj(Ep(i))*Ep(i) + conj(Em(i))*Em(i)) - d(i)*(conj(Ep(i))*Ep(i) - conj(Em(i))*Em(i)) );
%         d(i+1) = d(i) + Dt*(-gamma_d*d(i) -gamma*G(i)*(conj(Ep(i))*Ep(i) - conj(Em(i))*Em(i)) - gamma*d(i)*(conj(Ep(i))*Ep(i) + conj(Em(i))*Em(i)) );
%     end
%     psi = 0.5*angle(Em./Ep);
%     psi_devs(j) = std(psi(ceil(0.2*L):end));
% end
% 
% toc;
% 
% histogram(psi_devs, 20)


%% Plotting
% devPsi = std(psi(ceil(0.2*L):end))
% plot(time,Qp,'k',"LineWidth",1);
% hold on;
% plot(time,Qm,'r');
% plot(time,psi,'b');
% plot(time,G,'k');
% plot(time,d,'m');
% hold off;
% 
% tic;
% FindPsiStd(T, Dt, gamma, kappa, alpha, gamma_d, gamma_p, gamma_a, mu, C_sp, N_th, N_tr)
% toc;