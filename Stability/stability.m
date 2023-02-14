matrix = [0.677517, 0.232726, 0.029356, 0.850996, 0.301431, 0.848926;
          0.713054, 0.260037, 0.509707, 0.869024, 0.413981, 0.0924;
          0.731208, 0.748885, 0.619961, 0.664913, 0.344366, 0.166563;
          0.553577, 0.47379, 0.611191, 0.175111, 0.782086, 0.707329;
          0.933674, 0.942385, 0.752729, 0.856333, 0.632243, 0.0934586;
          0.0396751, 0.596296, 0.122535, 0.224435, 0.625694, 0.503896];


%% Parametres
g = 1.;
alpha = 5.;         % Henry factor
gd = 50.;           % d-decay rate
k = 1000.;          % photon decay rate
gp = 8.*2*pi;       % linear anisotropy
ga = -0.016*2*pi;   % birifrigence
mu = 2.;            % pump current


g = 1.;
alpha = 3.;         % Henry factor
gd = 50.;           % d-decay rate
k = 300.;          % photon decay rate
gp = 16.6681;       % linear anisotropy
ga = -0.01;   % birifrigence
mu = 1.9;            % pump current

% Initial values
qp0 = 0.1;      
qm0 = 0.9*2;
a0 = 0.;


% some setup
mus = linspace(1,4,100);
gps = logspace(1,2,100);

stab = zeros(100,100);
lin = zeros(100,100);
psi = zeros(100,100);
qp = zeros(100,100);
qm = zeros(100,100);

tic;
for i = 1:100
    for j = 1:100
        [p, D, x] = hurwitz([qp0, qm0, a0], g, alpha, gd, k, gps(j), ga, mus(i));
        if all(p >= 0)
            if all(D >= 0)
                stab(i,j) = 2;
            else
                stab(i,j) = 1;
            end
        else
            stab(i,j) = 0;
        end
        if x(1)-x(2) < 0.001*min(x(1:2))
            lin(i,j) = 1;
        end
        psi(i,j) = x(3);
        qp(i,j) = x(1);
        qm(i,j) = x(2);
    end
end
toc;



% [p, D, x] = hurwitz([qp0, qm0, a0], g, alpha, gd, k, gp, ga, mu);
% x
% 
% dmu = mu*1e-7;
% for n=1:10
%     if all([p, D] > 0)
%         pb = p; Db = D; xb=x;
%         [p, D, x] = hurwitz(xb, g, alpha, gd, k, gp, ga, mu+dmu);
%         d = ([p, D] - [pb, Db])/dmu;
%         musteps = [p, D] ./ d;
%         [~, ind] = min(abs(musteps));
%         mustep = musteps(ind);
%         disp(">");
%     else
%         pb = p; Db = D; xb=x;
%         [p, D, x] = hurwitz(xb, g, alpha, gd, k, gp, ga, mu+dmu);
%         tmpb = [pb, Db];
%         tmp = [p, D];
%         inds = find([pb, Db] < 0);
%         d = (tmp(inds) - tmpb(inds))./dmu;
%         musteps = tmp(inds) ./ d;
%         [~, ind] = max(abs(musteps));
%         mustep = musteps(ind);
%         disp(musteps);
%     end
%     
%     mu = mu + mustep;
%     disp([mu, mustep]);
%     dmu = mu*1e-7;
%     [p, D, x] = hurwitz(x, g, alpha, gd, k, gp, ga, mu); 
% end
% 
% 


% [p, D, x] = hurwitz([qp0, qm0, a0], g, alpha, gd, k, gp, ga, mu);
% x
% 
% if all(p > 0) 
%     if all(D > 0)
%         disp("System is stable");
%         fmt = format("shortE");
%         disp(D);
%         format(fmt)
%     else
%         disp("Some Hurwitz determinants is negative, so system is unstable:");
%         fmt = format("shortE");
%         disp(D);
%         format(fmt)
%     end
% else
%     disp("Some coefficients of charpoly is negative, so system is unstable:");
%     format short e;
%     disp(p);
%     format default;
% end

