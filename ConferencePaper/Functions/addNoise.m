function noisy_signal = addNoise(signal,SNR,varargin)
% Adds a gaussian noise to a signal proportional to a specific signal-to-noise ratio
% (in dB). Enter '[]' to disable the function

% Also add an extra variable: "distribution", which has options "uniform"
% or "normal", which changes the type of distribution applied
if ~isempty(varargin)
    if strcmpi('distribution',varargin{1})
        dist = varargin{2};
    else
        error("The two options for the 'distribution' argument ar 'normal' (default) or 'uniform'");
    end
else
    dist='normal'; %default to the normal distribution
end
%warning if SNR = 0
if SNR == 0
    warning('An SNR of 0 is not a noiseless signal! Set SNR to [] if you do not want any noise')
end

%check if SNR is given, if not, do nothing
if ~exist("SNR","var")
    SNR = [];
end

% If SNR is given as '[]', do nothing
if isempty(SNR)
    noisy_signal = signal;
else
    for i = 1:size(signal, 1)
        % Calculate the signal power
        P_signal = rms(signal(i,:)).^2;
        
        % Step 3: Determine the noise power
        P_noise = P_signal / (10^(SNR/20));
        
        % Step 4: Generate noise with the calculated noise power
        if strcmpi(dist, "normal")
            noise = sqrt(P_noise) * randn(1,size(signal,2));
        elseif strcmpi(dist, "uniform")
            noise = sqrt(P_noise) * rand(1,size(signal,2));
        end
        % Step 5: Add noise to the original signal
        noisy_signal(i,:) = signal(i,:) + noise;
    end
end
end