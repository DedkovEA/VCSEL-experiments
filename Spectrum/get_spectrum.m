%close a19*ll; 
%clc; clear;

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

%% Initializating variables---------------------------------------------

AV = 5;                     % Number for samples to obtain average
T = 10;                      % Window
offset = 0.2;                 % Window offset

tau = 2e-4;                   % Sampling time
% tauDt = ceil(tau/Dt);
% tau = Dt*tauDt;                 % Real sampling time
% wndfreq = pi/tau;               % Spectra will be saved only in +-wmdfreq
% wnd = floor(wndfreq*T/2/pi);
% 
% L = round(T/tau);                      % Points N in spectre
% if wnd >= L/2
%     wnd = floor(L/2) -1;
% end
% 
% wndslice = floor(L/2)+1-wnd:floor(L/2)+1+wnd;
% freqs = (-floor(L/2):ceil(L/2-1))/T; 
% freqswnd = freqs(wndslice);
% 
% 
% tic;
% 
% specsx = zeros(AV, length(freqswnd));
% specsy = zeros(AV, length(freqswnd));
% 
% % Initialization
% Qp_prev = 0;
% Qm_prev = 0;
% phi_prev = 0;
% psi_prev = 0;
% G_prev = 0;
% d_prev = 0;
% 
% Q = (-gamma_a + kappa*(mu-1+2*M*C_sp) + sqrt(4*(2*C_sp-1)*kappa*mu*(gamma_a+kappa) + (gamma_a+kappa*(1+2*C_sp*M+mu))*(gamma_a+kappa*(1+2*C_sp*M+mu))) )/(4*(gamma_a+kappa));
% 
% G_prev = mu/(1 + 2*Q);
% Qp_prev = Q;
% Qm_prev = Q;
% 
% Lmat = [2*kappa*(G_prev-1), -8*Q*gamma_p, 4*kappa*(C_sp+Q);
%         gamma_p/2/Q, 2*gamma_a, alpha*kappa;
%         -G_prev*gamma, 0, -gamma_d-2*gamma*Q];
% cp = charpoly(Lmat);
% hurwitz = [cp(2), cp(4), cp(2)*cp(3)-cp(4)];
% if sum(hurwitz <= 0) > 0
%     psi_prev = pi/2;
%     Lmat2 = [2*kappa*(G_prev-1), 8*Q*gamma_p, 4*kappa*(C_sp+Q);
%         -gamma_p/2/Q, -2*gamma_a, alpha*kappa;
%         -G_prev*gamma, 0, -gamma_d-2*gamma*Q];
%     cp2 = charpoly(Lmat2);
%     hurwitz2 = [cp2(2), cp2(4), cp2(2)*cp2(3)-cp2(4)];
%     if sum(hurwitz2 <= 0) > 0
%         disp("NOT STABLE")
%     end
% else
%     psi_prev = 0;
% end
% 
% 
% 
% ns = 1;
% 
% Ex = zeros(1,L); % X-LP mode
% Ey = zeros(1,L); % Y-LP mode
% 
% ei = 1;
% ni = 1;
% ksi_plus = normrnd(0, 1, 1, rnd_chunk_sz);
% ksi_minus = normrnd(0, 1, 1, rnd_chunk_sz);
% ksi_phi = normrnd(0, 1, 1, rnd_chunk_sz);
% ksi_psi = normrnd(0, 1, 1, rnd_chunk_sz);
% 
% d_from_prev = 0;
% while ns <= AV
%     if d_from_prev >= tauDt
%         d_from_prev = 0;
%         Ep = sqrt(Qp).*exp(1i.*(phi + psi));
%         Em = sqrt(Qm).*exp(1i.*(phi - psi));
%          Ex(ei) = (Ep + Em)/sqrt(2);
%          Ey(ei) = (Ep - Em)/sqrt(2)/1i;
%         % NEED Additional rotation on beta
%         % Ex(ei) = cos(ang)*(Ep + Em)/sqrt(2) + sin(ang)*(Ep - Em)/sqrt(2)/1i;
%         % Ey(ei) = -sin(ang)*(Ep + Em)/sqrt(2) + cos(ang)*(Ep - Em)/sqrt(2)/1i;
%         ei = ei+1;
%         if ei >= L
%             specx = abs(fftshift(fft(Ex))).^2;
%             specy = abs(fftshift(fft(Ey))).^2;
%             specsx(ns,:) = specx(wndslice);
%             specsy(ns,:) = specy(wndslice);
%             ns = ns + 1;
%             ei = ceil(offset*L) + 1;
%             Ex(1:ceil(offset*L)) = Ex(end+1-ceil(offset*L):end);
%             Ey(1:ceil(offset*L)) = Ey(end+1-ceil(offset*L):end);
%         end
%     end
% 
%     if ni > rnd_chunk_sz
%         ksi_plus = normrnd(0, 1, 1, rnd_chunk_sz);
%         ksi_minus = normrnd(0, 1, 1, rnd_chunk_sz);
%         ksi_phi = normrnd(0, 1, 1, rnd_chunk_sz);
%         ksi_psi = normrnd(0, 1, 1, rnd_chunk_sz);
%         ni = 1;
%     end
% 
%     Fp = 2*sqrt(C_sp*kappa*Qp_prev*(G_prev+d_prev+M))*ksi_plus(ni);
%     Fm = 2*sqrt(C_sp*kappa*Qm_prev*(G_prev-d_prev+M))*ksi_minus(ni);
%     Fphi = 1/2*sqrt(C_sp*(G_prev+d_prev+M)/Qp_prev/2)*(ksi_phi(ni)+ksi_psi(ni)) + 1/2*sqrt(C_sp*(G_prev-d_prev+M)/Qm_prev/2)*(ksi_phi(ni)-ksi_psi(ni));
%     Fpsi = 1/2*sqrt(C_sp*(G_prev+d_prev+M)/Qp_prev/2)*(ksi_phi(ni)+ksi_psi(ni)) - 1/2*sqrt(C_sp*(G_prev-d_prev+M)/Qm_prev/2)*(ksi_phi(ni)-ksi_psi(ni));
% 
%     Qp = Qp_prev + 2*( kappa*(G_prev+d_prev-1)*Qp_prev + kappa*C_sp*(G_prev+d_prev+M) - sqrt(Qp_prev*Qm_prev)*(gamma_a*cos(2*psi_prev-2*beta) + gamma_p*sin(2*psi_prev) ) )*Dt + Fp*sqrt(Dt);
%     Qm = Qm_prev + 2*( kappa*(G_prev-d_prev-1)*Qm_prev + kappa*C_sp*(G_prev-d_prev+M) - sqrt(Qp_prev*Qm_prev)*(gamma_a*cos(2*psi_prev-2*beta) - gamma_p*sin(2*psi_prev) ) )*Dt + Fm*sqrt(Dt);
%     phi = phi_prev + ( (G_prev-1)*alpha*kappa - 1/2/sqrt(Qm_prev*Qp_prev)*( (Qp_prev+Qm_prev)*gamma_p*cos(2*psi_prev) + (Qp_prev-Qm_prev)*gamma_a*sin(2*psi_prev-2*beta) ))*Dt + Fphi*sqrt(Dt);
%     psi = psi_prev + (alpha*kappa*d_prev + 1/2/sqrt(Qm_prev*Qp_prev)*( (Qp_prev-Qm_prev)*gamma_p*cos(2*psi_prev) + (Qp_prev+Qm_prev)*gamma_a*sin(2*psi_prev-2*beta) ))*Dt + Fpsi*sqrt(Dt);
%     G = G_prev + gamma*Dt*( (mu-G_prev) - G_prev*(Qp_prev + Qm_prev) - d_prev*(Qp_prev - Qm_prev) );
%     d = d_prev + Dt*(-gamma_d*d_prev -gamma*G_prev*(Qp_prev - Qm_prev) - gamma*d_prev*(Qp_prev + Qm_prev) );
% 
%     if abs(G) > 100
%         disp("ALERT!")
%         disp(gds(jj))
%         disp(mus(kk))
%         break
%     end
% 
%     Qp_prev = Qp;
%     Qm_prev = Qm;
%     phi_prev = phi;
%     psi_prev = psi;
%     G_prev = G;
%     d_prev = d;
% 
%     ni = ni + 1;
%     d_from_prev = d_from_prev + 1;
% end
% specx = mean(specsx);
% specy = mean(specsy);
% 
% toc;    % AV = 100 -> 534.354746 s

% plot(freqswnd, log10(specx));
% hold on;
% plot(freqswnd, log10(specy));
% hold off;

[specx,specy,freqs] = getspec(T, Dt, tau, rnd_chunk_sz, AV, offset, gamma, kappa, alpha, gamma_d, gamma_p, beta, gamma_a, mu, C_sp, N_th, N_tr);

sigma=0.75;
gauss = 1/2/pi/sigma*exp(-freqswnd.^2/2/sigma^2);
plot(freqswnd, log10(conv(specx,gauss,"same")), freqswnd, log10(conv(specy,gauss,"same")));

%save("TestVCSELspecWithStochastic.mat",'freqswnd','specx','specy','mu','kappa','alpha','gamma','gamma_d','gamma_a','gamma_p','C_sp','N_th','N_tr')

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