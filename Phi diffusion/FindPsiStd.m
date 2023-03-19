function [psi_means, psi_stds] = FindPsiStd(T, Dt, gamma, kappa, alpha, gamma_d, gamma_p, gamma_a, mu, C_sp, N_th, N_tr)
%FINDPSISTD Finds std(psi) in stationary process
%   Return means and stds of serial chunks of 1e7 points
    %% Initializating variables---------------------------------------------
    M = N_tr/(N_th - N_tr);
    
    Q = (-gamma_a + kappa*(mu-1+2*M*C_sp) + sqrt(4*(2*C_sp-1)*kappa*mu*(gamma_a+kappa) + (gamma_a+kappa*(1+2*C_sp*M+mu))*(gamma_a+kappa*(1+2*C_sp*M+mu))) )/(4*(gamma_a+kappa)); 
    G0 = mu/(1 + 2*Q);
    psi0 = 0;

    Lmat = [2*kappa*(G0-1), -8*Q*gamma_p, 4*kappa*(C_sp+Q);
            gamma_p/2/Q, 2*gamma_a, alpha*kappa;
            -G0*gamma, 0, -gamma_d-2*gamma*Q];
    cp = charpoly(Lmat);
    hurwitz = [cp(2), cp(4), cp(2)*cp(3)-cp(4)];
    if sum(hurwitz <= 0) > 0
        psi0 = pi/2;
    end
    
    %% ---------------------------------------------------------------------
    %% Simulating of time series of the X-LP and Y-LP light intensities-----
    
    %[~, ~, ~, psi, ~, ~] = EvolveSDE(T, Dt, Q, Q, psi0, G0, 0, gamma, kappa, alpha, gamma_d, gamma_p, gamma_a, mu, C_sp, N_th, N_tr);
    
%% TODO : split into chunks no more than 10 mln
    cycle_num = ceil(T/Dt/1e7);
    psi_means = zeros(1, cycle_num);
    psi_stds = zeros(1, cycle_num);

    if T/Dt < 1e7
        time = 0:Dt:T;
    else
        time = 0:Dt:1e7*Dt;
    end
    L = length(time);

    Qp0 = Q;
    Qm0 = Q;
    d0 = 0;
    for n = 1:cycle_num
        Qp = zeros(1,L); % Plus mode
        Qm = zeros(1,L); % Minus mode
        psi = zeros(1,L); % Angle
        G = zeros(1,L);  % Total carrier inversion
        d = zeros(1,L);  % The difference between carrier inversions with opposite spins
        Qp(1) = Qp0;
        Qm(1) = Qm0;
        psi(1) = psi0;
        G(1) = G0;
        d(1) = d0;
        
        ksi_plus = normrnd(0, 1, 1, L);
        ksi_minus = normrnd(0, 1, 1, L);
        ksi_psi = normrnd(0, 1, 1, L);
        
        for i = 1:1:L-1
            Fp = 2*sqrt(C_sp*kappa*Qp(i)*(G(i)+d(i)+M))*ksi_plus(i);
            Fm = 2*sqrt(C_sp*kappa*Qm(i)*(G(i)-d(i)+M))*ksi_minus(i);
            Fpsi = 1/2*sqrt(C_sp*( (G(i)+d(i)+M)/Qp(i) + (G(i)-d(i)+M)/Qm(i) ))*ksi_psi(i);
    
            Qp(i+1) = Qp(i) + 2*( kappa*(G(i)+d(i)-1)*Qp(i) + kappa*C_sp*(G(i)+d(i)+M) - sqrt(Qp(i)*Qm(i))*(gamma_a*cos(2*psi(i)) + gamma_p*sin(2*psi(i)) ) )*Dt + Fp*sqrt(Dt);
            Qm(i+1) = Qm(i) + 2*( kappa*(G(i)-d(i)-1)*Qm(i) + kappa*C_sp*(G(i)-d(i)+M) - sqrt(Qp(i)*Qm(i))*(gamma_a*cos(2*psi(i)) - gamma_p*sin(2*psi(i)) ) )*Dt + Fm*sqrt(Dt);
            psi(i+1) = psi(i) + (alpha*kappa*d(i) + 1/2/sqrt(Qm(i)*Qp(i))*( (Qp(i)-Qm(i))*gamma_p*cos(2*psi(i)) + (Qp(i)+Qm(i))*gamma_a*sin(2*psi(i)) ))*Dt + Fpsi*sqrt(Dt);
            G(i+1) = G(i) + gamma*Dt*( (mu-G(i)) - G(i)*(Qp(i) + Qm(i)) - d(i)*(Qp(i) - Qm(i)) );
            d(i+1) = d(i) + Dt*(-gamma_d*d(i) -gamma*G(i)*(Qp(i) - Qm(i)) - gamma*d(i)*(Qp(i) + Qm(i)) );
        end
        Qp0 = Qp(end);
        Qm0 = Qm(end);
        psi0 = psi(end);
        G0 = G(end);
        d0 = d(end);
        psi_stds(n) = std(psi);
        psi_means(n) = mean(psi);
    end
end

