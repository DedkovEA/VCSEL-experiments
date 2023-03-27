%close all; 
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

Dt = 1e-5;          % Time step
T = 200;             % Time for solving in ns


%% Initializating variables---------------------------------------------

time = 0:Dt:T;
L = length(time);

AV = 20;
wnd = 100000;
specsx = zeros(AV, 2*wnd);
specsy = zeros(AV, 2*wnd);

tic
for j = 1:AV
    
    flag_noise = true; % Flag of noise
    Qp = zeros(1,L); % X-LP mode
    Qm = zeros(1,L); % Y-LP mode
    psi = zeros(1,L);
    phi = zeros(1,L);
    G = zeros(1,L);  % Total carrier inversion
    d = zeros(1,L);  % The difference between carrier inversions with opposite spins
    
    Q = 1/2*(kappa*mu/(gamma_a + kappa) - 1);
    
    Q = (-gamma_a + kappa*(mu-1+2*M*C_sp) + sqrt(4*(2*C_sp-1)*kappa*mu*(gamma_a+kappa) + (gamma_a+kappa*(1+2*C_sp*M+mu))*(gamma_a+kappa*(1+2*C_sp*M+mu))) )/(4*(gamma_a+kappa));
    
    G(1) = mu/(1 + 2*Q);
    Qp(1) = Q;
    Qm(1) = Q;
    
    
    Lmat = [2*kappa*(G(1)-1), -8*Q*gamma_p, 4*kappa*(C_sp+Q);
            gamma_p/2/Q, 2*gamma_a, alpha*kappa;
            -G(1)*gamma, 0, -gamma_d-2*gamma*Q];
    cp = charpoly(Lmat);
    hurwitz = [cp(2), cp(4), cp(2)*cp(3)-cp(4)];
    if sum(hurwitz <= 0) > 0
        psi(1) = pi/2;
    end
    
    %% ---------------------------------------------------------------------
    %% Simulating of time series of the X-LP and Y-LP light intensities-----
    
    ksi_plus = wgn(1, L, 0);
    ksi_minus = wgn(1, L, 0);
    ksi_phi = wgn(1, L, 0);
    ksi_psi = wgn(1, L, 0);
    Fp = 0; Fm = 0; Fpsi = 0; Fphi = 0;
    
    flag_noise = true;
    for i = 1:1:L-1
        if(flag_noise)
            Fp = 2*sqrt(C_sp*kappa*Qp(i)*(G(i)+d(i)+M))*ksi_plus(i);
            Fm = 2*sqrt(C_sp*kappa*Qm(i)*(G(i)-d(i)+M))*ksi_minus(i);
            % Fpsi = 1/2*sqrt(C_sp*( (G(i)+d(i)+M)/Qp(i) + (G(i)-d(i)+M)/Qm(i) ))*ksi_psi(i);
            Fphi = 1/2*sqrt(C_sp*(G(i)+d(i)+M)/Qp(i)/2)*(ksi_phi(i)+ksi_psi(i)) + 1/2*sqrt(C_sp*(G(i)-d(i)+M)/Qm(i)/2)*(ksi_phi(i)-ksi_psi(i));
            Fpsi = 1/2*sqrt(C_sp*(G(i)+d(i)+M)/Qp(i)/2)*(ksi_phi(i)+ksi_psi(i)) - 1/2*sqrt(C_sp*(G(i)-d(i)+M)/Qm(i)/2)*(ksi_phi(i)-ksi_psi(i));
        end
        Qp(i+1) = Qp(i) + 2*( kappa*(G(i)+d(i)-1)*Qp(i) + kappa*C_sp*(G(i)+d(i)+M) - sqrt(Qp(i)*Qm(i))*(gamma_a*cos(2*psi(i)) + gamma_p*sin(2*psi(i)) ) )*Dt + Fp*sqrt(Dt);
        Qm(i+1) = Qm(i) + 2*( kappa*(G(i)-d(i)-1)*Qm(i) + kappa*C_sp*(G(i)-d(i)+M) - sqrt(Qp(i)*Qm(i))*(gamma_a*cos(2*psi(i)) - gamma_p*sin(2*psi(i)) ) )*Dt + Fm*sqrt(Dt);
        phi(i+1) = phi(i) + ( (G(i)-1)*alpha*kappa - 1/2/sqrt(Qm(i)*Qp(i))*( (Qp(i)+Qm(i))*gamma_p*cos(2*psi(i)) + (Qp(i)-Qm(i))*gamma_a*sin(2*psi(i)) ))*Dt + Fphi*sqrt(Dt);
        psi(i+1) = psi(i) + (alpha*kappa*d(i) + 1/2/sqrt(Qm(i)*Qp(i))*( (Qp(i)-Qm(i))*gamma_p*cos(2*psi(i)) + (Qp(i)+Qm(i))*gamma_a*sin(2*psi(i)) ))*Dt + Fpsi*sqrt(Dt);
        G(i+1) = G(i) + gamma*Dt*( (mu-G(i)) - G(i)*(Qp(i) + Qm(i)) - d(i)*(Qp(i) - Qm(i)) );
        d(i+1) = d(i) + Dt*(-gamma_d*d(i) -gamma*G(i)*(Qp(i) - Qm(i)) - gamma*d(i)*(Qp(i) + Qm(i)) );
    end
    
    Ep = sqrt(Qp).*exp(1i.*(phi + psi));
    Em = sqrt(Qm).*exp(1i.*(phi - psi));
    Ex = (Ep + Em)/sqrt(2);
    Ey = (Ep - Em)/sqrt(2)/1i;
    
    specx = abs(fftshift(fft(Ex))).^2;
    specy = abs(fftshift(fft(Ey))).^2;
    specsx(j,:) = specx(floor(L/2-wnd):floor(L/2-wnd)+2*wnd-1);
    specsy(j,:) = specy(floor(L/2-wnd):floor(L/2-wnd)+2*wnd-1);
end
toc;
avspecx = mean(specsx);
avspecy = mean(specsy);

freqs = 2*pi*(-L/2:L/2-1)/L/Dt; 
freqswnd = freqs(floor(L/2-wnd):floor(L/2-wnd)+2*wnd-1);

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