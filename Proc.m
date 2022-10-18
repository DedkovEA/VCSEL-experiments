% Parametres
find_Tspan = true;      % Set true in order to autiomatically find impulse period
find_Tstart = true;     % Set true in order to autiomatically find the beginning of 1st pulse
find_Tarea = true;      % Set true in order to autiomatically find dT for area measurement
find_Toff = 0;          % Set true and T_off will be automatically set to floor(Tarea/2)

Tspan = 32;             % Pulse period
Tstart = 28;             % Beginning of the first pulse area measurment
Tarea = 14;             % dT in area measurement
Toff = 5;               % time offset from Tstart to comparator trigger time

Np = 10;                % Number of check pulses to be plotted

% Load waveform
waveform = load("17.10.22\waveform_6,72mA_1023(36).dat");

V = waveform(:,2);
T = waveform(:,1);

% Finding period
if find_Tspan
    SItspan = 1000;  % Timespan with several indicies
    [r, l] = xcorr(V, V, SItspan, "normalized");
    [pksy, pksx] = findpeaks(r);
    Tspan = round(mean(pksx(2:end) - pksx(1:end-1)));
end

% Finding Tstart
% TODO

% Finding Tarea
% TODO

% Finding Toff
if find_Toff
    Toff = floor(Tarea / 2);
end


dt = T(2) - T(1);
N = floor( (length(T) - Tstart) / Tspan) - 1;
% Plotting several pulses (almost equally spaced in time) with area region and Toff
close all hidden;
fPulses = figure;
for k = linspace(0, N, Np)
    k = floor(k);
    plot((max(Tstart-5, 1) : Tstart + Tspan), ...
         V(max(Tstart-5, 1) + Tspan*k : Tstart + Tspan*(k+1)))
    hold on;
    xline((Tstart), "r", "Tstart");
    xline((Tarea+Tstart), "r", "Tarea");
    xline((Toff+Tstart), "g", "Toff");
end
hold off;

% Area calculations
% Possible TODO - convolution FFT speedup
S = zeros(N, 1);
for k = 0:N-1
    S(k+1) = sum(V(Tstart + Tspan*k : Tstart + Tarea + Tspan*k)) * dt;
end

% Binary sequence based on area
Sbin = zeros(N, 1);
Smedian = median(S);
Sbin(S > Smedian) = 1;

% Selfcorrelation plot
[Sr, Sl] = autocorr_func(S, "method", "matlab", "positiveonly", true);
[Sbinr, Sbinl] = autocorr_func(Sbin, "method", "matlab", "positiveonly", true);

fScorr = figure;
plot(Sl, Sr, "b");
hold on;
plot(Sbinl, Sbinr, "r");
hold off;