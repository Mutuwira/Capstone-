function [time_estimate, corr_output, cfo_estimate] = synchronization_with_CFO(rx_signal, tx_signal, cp_length, seq_length)
    % Cross-correlation for synchronization
    corr_output = abs(corr1(rx_signal, tx_signal));
    % Find the peak of cross-correlation to detect the start of OFDM symbols
    [~, peak_index] = max(corr_output);
    time_estimate = peak_index; % Adjust for lag in cross-correlation

    % Extract the CP portion from the received signal excluding the delay
    cp_portion = rx_signal(peak_index:peak_index+cp_length-1);  % cp part of the signal
    zc_portion = rx_signal(peak_index + seq_length: peak_index + seq_length + cp_length - 1); % last cp_length values of the ZC sequencee
    phase_diff = conj(cp_portion).*(zc_portion);
    cfo_estimate = angle(sum(phase_diff))/(2*pi*seq_length);
end