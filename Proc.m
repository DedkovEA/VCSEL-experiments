% Parametres
find_Tspan = true;      % Set true in order to autiomatically find impulse period
find_Tstart = true;     % Set true in order to autiomatically find the beginning of 1st pulse
find_Tarea = false;      % Set true in order to autiomatically find dT for area measurement
find_Toff = true;       % Set true and T_off will be automatically set to floor(Tarea/2)

Tspan = 32;             % Pulse period
Tstart = 28;            % Beginning of the first pulse area measurment
Tarea = 32;             % dT in area measurement
Toff = 5;               % time offset from Tstart to comparator trigger time

Np = 10;                % Number of check pulses to be plotted
level = 0.1;             % Level of comparison for tstart and area calculation

% Load waveform

V = load("24112022_qrng\C2--C2_2_56mA_551_9mA_30303030_--00000.dat");

%waveform = load("sim_data.mat");
%wfm2 = load("sim_data_Pxx2.mat");

%V = waveform.Pxosc;
%V = wfm2.Pxx2;
%T = waveform(:,1);

%V = load("31.11\C3--Trace--3030303-c473-b2d95---00000.dat");


close all hidden;

% Finding period
if find_Tspan
    Tspan = find_period(V, 1000);
end

% Finding Tstart and Tarea by the highest pike
if find_Tstart || find_Tarea
    [start0, dt0] = find_arearegion(V, Tspan, level);
    if find_Tstart
        Tstart = start0;
    end
    if find_Tarea
        Tarea = dt0;
    end
end


% Finding Toff
if find_Toff
    Toff = floor(Tarea / 2);
end


try
    dt = T(2) - T(1);
catch excp
    fprintf("Something goes wrong wih dt, omitting it and dt=1 \n");
    dt = 1;
end

N = floor( (length(V) - Tstart) / Tspan) - 1;
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
    S(k+1) = sum(V(Tstart + Tspan*k : Tstart + Tarea + Tspan*k));
end

% Binary sequence based on area
Sbin = zeros(N, 1);
%Smedian = median(S);
%Smedian = 1.75e-11;
Smedian = (max(S)+min(S))/2;
Sbin(S > Smedian) = 1;

% Von Neuman extraction
newS = zeros(floor(N/2), 1);
ck = 1;
for i = 1:floor(N/2)
    if Sbin(2*i-1) ~= Sbin(2*i)
        newS(ck) = Sbin(2*i); % 10 -> 0, 01 -> 1
        ck = ck + 1;
    end
end
newS = newS(1:ck-1);

% Selfcorrelation plot
[Sr, Sl] = autocorr_func(S, "method", "matlab", "positiveonly", true);
[Sbinr, Sbinl] = autocorr_func(Sbin, "method", "matlab", "positiveonly", true);
[newSr, newSl] = autocorr_func(newS, "method", "matlab", "positiveonly", true);

fScorr = figure;
plot(Sl, Sr, "b");
hold on;
plot(Sbinl, Sbinr, "r");
plot(newSl, newSr, "g");
hold off;