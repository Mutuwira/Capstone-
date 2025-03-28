function rx_signal = awgn_noise(tx_signal, SNR_dB)
    delay = randi([1, 100]);  % Random delay between 1 and 100 samples
    tx_signal1 = [zeros(1, delay), tx_signal];  
    SNR = 10^(SNR_dB/10);  % Convert SNR from dB to linear scale
    noise_power = (1/SNR) * var(tx_signal1);
    noise = sqrt(noise_power/2) * (randn(size(tx_signal1)) + 1i*randn(size(tx_signal1)));
    rx_signal = tx_signal1 + noise;  % Received signal with noise
end
