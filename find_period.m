function [Tspan] = find_period(data, SItspan)
%FIND_PERIOD Summary of this function goes here
%   Detailed explanation goes here
    arguments
        data (:,1) double 
        SItspan (1,1) double {mustBeInteger} = round(length(data)/10) % Timespan with several indicies
    end
    
    [r, ~] = xcorr(data, data, SItspan, "normalized");
    [~, pksx] = findpeaks(r);
    Tspan = round(mean(pksx(2:end) - pksx(1:end-1)));
end

