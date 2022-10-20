% Parametres
find_Tspan = true;      % Set true in order to autiomatically find impulse period
find_Tstart = true;     % Set true in order to autiomatically find the beginning of 1st pulse
find_Tarea = true;      % Set true in order to autiomatically find dT for area measurement
find_Toff = true;       % Set true and T_off will be automatically set to floor(Tarea/2)

Tspan = 32;             % Pulse period
Tstart = 28;            % Beginning of the first pulse area measurment
Tarea = 14;             % dT in area measurement
Toff = 5;               % time offset from Tstart to comparator trigger time

Np = 10;                % Number of check pulses to be plotted
level = 0.1;             % Level of comparison for tstart and area calculation

% Load waveform
waveform = load("17.10.22\waveform_6,23mA_1023(36).dat");

V = waveform(:,2);
T = waveform(:,1);

close all hidden;

% Finding period
if find_Tspan
    SItspan = 1000;  % Timespan with several indicies
    [r, l] = xcorr(V, V, SItspan, "normalized");
    [pksy, pksx] = findpeaks(r);
    Tspan = round(mean(pksx(2:end) - pksx(1:end-1)));
end

% Finding Tstart and Tarea by the highest pike
if find_Tstart || find_Tarea
    pulse = V(1 : Tspan*min(floor(length(T)/Tspan), 100));
    pulse = pulse - min(pulse);
    [amp, ind] = max(pulse);
    fronts = diff(pulse(Tspan*fix(ind/Tspan - 1) + 1 : Tspan * fix(ind/Tspan + 2)) > level*amp);
    start0 = rem(find(fronts(1:rem(ind, Tspan)+Tspan+1) > 0.5, 1, "last"), Tspan);
    dt0 = find(fronts(start0:end) < -0.5, 1, "first");
    %startinds = zeros(1, min(floor(length(T)/Tspan), 100) + 1);
    %dtinds = zeros(1, min(floor(length(T)/Tspan), 100) + 1);
    %wnd = floor(Tspan/8);
    %for k = 0:min(floor(length(T)/Tspan), 100)
    %    pulse = V(Tspan*k + startind0 - wnd : Tspan*(k+1) + startind0 + wnd);
    %    pulse = pulse - min(pulse);
    %    amp = max(pulse);
    %    fronts = diff(pulse > level*amp);
    %    startinds(k+1) = find(fronts > 0.5, 1, 'first') + 1 - wnd + startind0;
    %    %dtinds(k+1) = find(fronts(startinds(k+1):end) < -0.5, 1, 'first') + 1;
    %    plot(fronts)
    %    hold on;
    %end
    %hold off;
    if find_Tstart
        Tstart = round(start0);
    end
    if find_Tarea
        Tarea = round(dt0);
    end
end


% Finding Toff
if find_Toff
    Toff = floor(Tarea / 2);
end


dt = T(2) - T(1);
N = floor( (length(T) - Tstart) / Tspan) - 1;
% Plotting several pulses (almost equally spaced in time) with area region and Toff

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