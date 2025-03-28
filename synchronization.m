function [time_estimate, corr_output] = synchronization(rx_signal, tx_signal)
    % Cross-correlation for synchronization
    corr_output = abs(corr1(rx_signal, tx_signal));
    % Find the peak of cross-correlation to detect the start of OFDM symbols
    [~, peak_index] = max(corr_output);
    time_estimate = peak_index; % Adjust for lag in cross-correlation
end