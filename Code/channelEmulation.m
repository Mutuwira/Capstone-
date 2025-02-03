function rx_signal = channelEmulation(tx_signal, SNR_dB, delay, h);
    % Add delay
    y_delayed = [zeros(1, delay), tx_signal]; % Create offset
    y_channel = conv(h, y_delayed); % Convolve with channel response
    signalPower = mean(abs(y_channel).^2); % Compute signal power
    noisePower = signalPower / (10^(SNR_dB/10));  % Compute noise power based on SNR
    
    % Generate AWGN noise
    noise = sqrt(noisePower/2) .* (randn(size(y_channel)) + 1i .* randn(size(y_channel)));
    rx_signal = y_channel + noise;     % Add noise to the signal
end

