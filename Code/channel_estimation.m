function H_est = channel_estimation(rx_signal, preamble, seq_length, N, noiseVariance, signalPower)
    % Determine the CP length (assumes preamble = [CP, clean_sequence])
    CP_length = length(preamble) - seq_length;
    
    % Check that the received signal is long enough to extract the full preamble
    if length(rx_signal) < length(preamble)
        error('Received signal is shorter than the preamble length.');
    end
    
    % Extract the preamble from the received signal.
    % (We assume that the preamble is the first length(preamble) samples.)
    rx_preamble = rx_signal(1:length(preamble));
    
    % Remove CP from both the received and transmitted preamble
    rx_preamble_noCP = rx_preamble(CP_length+1:end);
    tx_preamble_noCP = preamble(CP_length+1:end);
    
    % Compute the FFT (using FFT length = seq_length)
    RX_FFT = fft(rx_preamble_noCP, seq_length);
    TX_FFT = fft(tx_preamble_noCP, seq_length);
    %RX_FFT = fft(rx_preamble_noCP, N);
    %TX_FFT = fft(tx_preamble_noCP, N);
    
    % LS Channel Estimation: elementwise division in the frequency domain
    H_ls = RX_FFT ./ TX_FFT;
    
    % If seq_length does not equal N, interpolate to obtain N channel estimates.
    if seq_length ~= N
        % Create normalized frequency grids for the original and target FFT bins
        freqOrig = linspace(0, 1, seq_length);
        freqInterp = linspace(0, 1, N);
        H_est = interp1(freqOrig, H_ls, freqInterp, 'linear', 'extrap').';
    else
        H_est = H_ls;
    end
    
    % Optional MMSE correction: weight the LS estimate based on the estimated SNR.
    H_est = H_est .* (abs(H_est).^2 ./ (abs(H_est).^2 + (noiseVariance / signalPower)));
    
    % True Channel
    h = [1+1i*0.5 0.2+1i*0.6 0.1+1i*0.4 0.05+1i*0.01]; 
    H_true = fft(h,N);

    % Plot the channel in logarithmic scale
    H_est_dB= 10*log10(abs(H_est));
    H_true_dB = 10*log10(abs(H_true));
    plot(H_est_dB, 'r');
    hold on;
    plot(H_true_dB, 'b')
    hold off;
    title("Estimated Channel (in dB)")
    xlabel('Subcarrier Index');
    ylabel('Magnitude (dB)');
    legend('Estimated Channel','True Channel');
    grid on;

end
