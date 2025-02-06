function rx_signal = AWGN_with_CFO(tx_signal, SNR_dB, cfo, fs)

    % delay = randi([1, 100]);  % Random delay between 1 and 100 samples
    % tx_signal1 = [zeros(1, delay), tx_signal];  % Add delay to the signal

    % Apply CFO by introducing a phase rotation
    t = (0:length(tx_signal)-1) / fs;  % Time vector
    cfo_phase = exp(1i * 2 * pi * cfo * t);  % CFO-induced phase rotation
    tx_signal_with_cfo = tx_signal .* cfo_phase;  % Apply CFO to the signal

    SNR = 10^(SNR_dB/10);  % Convert SNR from dB to linear scale
    noise_power = (1/SNR)* var(tx_signal_with_cfo);

    % Generate AWGN noise
    noise = sqrt(noise_power/2) * (randn(size(tx_signal_with_cfo)) + 1i * randn(size(tx_signal_with_cfo)));

    rx_signal = tx_signal_with_cfo + noise;     % Add noise to the signal with CFO

end

% Next Step do the frequency offset with OFDM cps
% Plot MSE with SNR for frequency offset with (OFDM only, ZC only, OFDM + ZC)
