function [start, dt] = find_arearegion(data, Tspan, level)
%FIND_AREAREGION Summary of this function goes here
%   Detailed explanation goes here

    arguments
        data (:,1) double 
        Tspan (1,1) double {mustBeInteger}
        level (1,1) double = 0.1 
    end

    pulse = data(1 : Tspan*min(floor(length(data)/Tspan), 100));
    pulse = pulse - min(pulse);
    [amp, ind] = max(pulse);
    fronts = diff(pulse(Tspan*fix(ind/Tspan - 1) + 1 : Tspan * fix(ind/Tspan + 2)) > level*amp);
    start = rem(find(fronts(1:rem(ind, Tspan)+Tspan+1) > 0.5, 1, "last"), Tspan);
    dt = find(fronts(start:end) < -0.5, 1, "first");

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
end

