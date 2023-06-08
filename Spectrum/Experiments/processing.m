file1 = load("20.04.23\Polarization_resolved_spectra_2,04mA_1.mat");
freqs = table2array(file1.S20230420_113702(:,1));
file3 = load("20.04.23\Polarization_resolved_spectra_2,04mA_3.mat");

file = load("20.04.23\current_corr_fix.mat");
currents = [1.26, 2.04, 3.01, 4.0, 4.99];
xspecs = [table2array(file.X_20230420_1_26_12_26(:,3)), ...
          table2array(file1.S20230420_113702(:,3)), ...
          table2array(file.X_20230420_3_01_12_29(:,3)), ... 
          table2array(file.X_20230420_4_00_12_34(:,3)), ...
          table2array(file.X_20230420_4_99_12_35(:,3))];

yspecs = [table2array(file.Y_20230420_1_26_12_25(:,3)), ...
          table2array(file3.S20230420_113901(:,3)), ...
          table2array(file.Y_20230420_3_01_12_28(:,3)), ... 
          table2array(file.Y_20230420_4_00_12_30(:,3)), ...
          table2array(file.Y_20230420_4_99_12_35_54(:,3))];


df = freqs(2) - freqs(1);
f0 = freqs(1);
f1 = 194.75;
f2 = 195.15;
i1 = ceil((f1 - f0)/df);
i2 = floor((f2-f0)/df);
wnd = i1:i2;

C = {'k','r','g','b','m'};
for i = 1:5
    plot(freqs(wnd), log(xspecs(wnd,i)), "Color", C{i}, "LineStyle","--")
    hold on;
    plot(freqs(wnd), log(yspecs(wnd,i)), "Color", C{i})
    hold on;
end
hold off;

ypikes = [195.038, 195.005, 194.956, 194.901, 194.834];
xpikes = [195.019, 194.986, 194.937, 194.882, 194.832];
fwmpikes = [194.97, 194.95, 194.902, 194.847, 194.789];

d1 = ypikes - xpikes; % X_Y pikes distance
d2 = xpikes - fwmpikes; % X_FWM pikes distance