function H_est = channel_estimation(rx_signal, preamble, seq_length, N, noiseVariance, signalPower)
% CHANNEL_ESTIMATION estimates the channel frequency response using a preamble.
%
%   H_est = channel_estimation(rx_signal, preamble, seq_length, N, noiseVariance, signalPower)
%
%   INPUTS:
%     rx_signal     : Received (CFO-corrected) signal that contains the preamble
%                     at the very beginning.
%     preamble      : Transmitted preamble signal (with cyclic prefix) used for estimation.
%     seq_length    : Length of the preamble sequence (without the cyclic prefix).
%     N             : Number of OFDM subcarriers (FFT size for data symbols).
%     noiseVariance : Estimated noise variance.
%     signalPower   : Average signal power.
%
%   OUTPUT:
%     H_est         : Estimated channel frequency response (N x 1 vector).
%
%   The function works as follows:
%     1. Determine the cyclic prefix (CP) length from the known preamble.
%     2. Extract the received preamble from rx_signal (assumed to start at sample 1).
%     3. Remove the CP from both the received preamble and the known (transmitted) preamble.
%     4. Compute the FFT of both (using FFT length = seq_length).
%     5. Compute the LS channel estimate: H_ls = FFT(received_preamble) ./ FFT(transmitted_preamble).
%     6. If seq_length and N differ, interpolate the estimate to yield N estimates.
%     7. Optionally, apply an MMSE-style correction.
%
%   NOTE:
%     Ensure that the first argument is the received signal—not the transmitted one.
%     (If needed, modify your main code accordingly.)

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
    
end
