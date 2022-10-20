function [r, lags] = autocorr_func(data, options)
    arguments
        data (:,1) double 
        options.maxlag (1,1) double {mustBePositive, mustBeInteger} = length(data)-1
        options.method string {mustBeMember(options.method,["matlab","manual"])} = "matlab"
        options.positiveonly logical = false
    end;
    m = mean(data);
    ndata = (data - m) / sqrt(mean(data.^2) - m^2);
    if options.method == "matlab"
        [r, lags] = xcorr(ndata, ndata, options.maxlag, 'biased');
    elseif options.method == "manual"
        N = length(ndata);
        if options.maxlag >= N
            options.maxlag = N-1;
        end
        r = zeros(2*options.maxlag+1, 1);
        lags = -options.maxlag:options.maxlag;
        for k = 0:options.maxlag
            r(options.maxlag+1+k) = sum(ndata(1:N-k) .* conj(ndata(1+k:N))) / N;
            r(options.maxlag+1-k) = conj(r(options.maxlag+1+k));
        end
    end

    if options.positiveonly
        r = r(options.maxlag+1:end);
        lags = lags(options.maxlag+1:end);
    end
end