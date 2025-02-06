
seq_length = 64;         % Length of ZC sequence (prime number)
root = 1;               % Root of ZC sequence
SNR_dB = -30:1:30;            % Signal-to-noise ratio (in dB)
cp_length = 16;          % Cyclic prefix length
CFO = randi([1 20]);     % Carrier Frequency Offset (in Hz)
fs = 10000;              % Sampling frequency (in Hz)

tx_signal = generate_signal(seq_length, root, cp_length);  % Generate transmitted signal (ZC sequence with CP and data)

t_MSE = zeros(1, length(SNR_dB));  % preallocation for MSE for time estimate
f_MSE = zeros(1, length(SNR_dB));  % preallocation for MSE for frequency estimate
num_t = 40; % Number of trials to calculate MSE for each SNR

for i = 1:length(SNR_dB)
    t_errors = zeros(1, num_t);  % preallocation for time estimate errors for each SNR
    f_errors = zeros(1, num_t);  % preallocation for frequency estimate errors for each SNR
    for trial = 1:num_t
        delay = randi([1 100]); % Random delay between 1 and 100 samples
        tx_signal1 = [zeros(1, delay), tx_signal];  % Add delay to the signal

        rx_signal = AWGN_with_CFO(tx_signal1, SNR_dB(i), CFO, fs);   % Generate recieved signal with CFO, delay, and noise
        
        % Call synchronization method to get the time and frequency estimate along
        % with the correlation output
        [time_estimate, corr_output, cfo_estimate] = synchronization_with_CFO(rx_signal, tx_signal, cp_length, seq_length);
        cfo_estimate_hz = cfo_estimate * fs;

        % Calculate errors
        t_errors(trial) = (time_estimate - (delay+1))^2;  % Squared time error
        f_errors(trial) = (cfo_estimate_hz - CFO)^2;   % Squared CFO error
    end
    % Compute mean square error for current SNR
    t_MSE(i) = mean(t_errors);
    f_MSE(i) = mean(f_errors);
end

% Plots for MSE
figure;
subplot(2, 1, 1);
plot(SNR_dB, t_MSE, 'k-*');
xlabel('SNR (dB)');
ylabel('MSE (Time Estimate)');
title('MSE vs SNR for Time Estimate');
grid on;

subplot(2, 1, 2);
plot(SNR_dB, f_MSE, 'b-o');
xlabel('SNR (dB)');
ylabel('MSE (CFO Estimate)');
title('MSE vs SNR for CFO Estimate');
grid on;



