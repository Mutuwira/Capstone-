function [time_estimate, corr_output, cfo_estimate] = synchronization(rx_signal, preamble, CP, seq_length)
    corr_output = abs(corr1(rx_signal, preamble));
    [~, peak_index] = max(corr_output);
    time_estimate = peak_index;

    cp_portion = rx_signal(peak_index:peak_index+CP-1);
    zc_portion = rx_signal(peak_index + seq_length: peak_index + seq_length + CP - 1);
    phase_diff = conj(cp_portion) .* zc_portion;
    cfo_estimate = angle(sum(phase_diff)) / (2 * pi * seq_length);
    %cfo_estimate = (10000 / seq_length) * angle(sum(phase_diff)) / (2 * pi);


    function [c, lags] = corr1(rx_signal, tx_signal)
    len_tx = length(tx_signal);  % Get length of tx_signal
    c = [];     % Initialize empty arrays for correlation and lags
    lags = [];
    i = 0;        % Initialize index for lags
    
    while true 
        try       % extract a segment of rx_signal of length len_tx
            rx_signal_i = rx_signal(i + 1 : i + len_tx); 
        catch
           break; % Break loop if we reach the end of rx_signal

        end
        % Compute the correlation at this delay
        c(i + 1) = sum(conj(rx_signal_i) .* tx_signal); 
        lags(i + 1) = i;
        
        i = i + 1; % Increment the index for the next lag
    end
    
    % Normalize the correlation output
    c = c / sqrt(sum(abs(rx_signal(1:len_tx)).^2) * sum(abs(tx_signal).^2));
end
end
