Dt = 1e-3;          % Time step
T = 10000;             % Time for solving in ns

time = 0:Dt:(T-Dt);
L = length(time);
N = 10^-3;


ksi = wgn(1, L, 0);
x = zeros(1, L);
for i = 1:1:L-1
    x(i+1) = x(i) - x(i)*x(i)*x(i)*Dt + sqrt(N)*ksi(i)*sqrt(Dt);
end
    
plot(time, x)
std(x)

% 
% genspec = zeros(1,L);
% freqs = 2*pi*(-L/2:L/2-1)/L/Dt;
% AV = 1;
% 
% tic
% for j = 1:AV
%     ksi = wgn(1, L, 0);
%     x = zeros(1, L);
%     for i = 1:1:L-1
%         x(i+1) = x(i) - x(i)*x(i)*x(i)*Dt + sqrt(N)*ksi(i)*sqrt(Dt);
%     end
%     
%     % plot(time, x)
%     
%     [ac, lags] = xcorr(x,ceil(L/2),'unbiased');
%     spec = fftshift(fft(ac(1:L)))/L;
%     genspec = genspec + spec;
% end
% toc
% 
% plot(freqs(L/2-2000:L/2+2000), abs(genspec(L/2-2000:L/2+2000))/AV);