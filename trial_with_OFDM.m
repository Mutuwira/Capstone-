clc
clear
close all

%% Parameters
% OFDM Parameters
N = 64;                % Number of OFDM subcarriers
M = 16;                % 16-QAM
numBits = N * log2(M); % Total number of bits
bitsPerSymbol = log2(M);
CP = N * 0.5;          % Cyclic Prefix Length

% Synchronization Parameters
seq_length = 64;       % Length of ZC sequence
root = 1;              % Root of ZC sequence
SNR_dB = 10;           % Signal-to-noise ratio (in dB)
cp_length = 16;        % Cyclic prefix length for ZC sequence
CFO = randi([1 20]);   % Carrier Frequency Offset (in Hz)
fs = 10000;            % Sampling frequency (in Hz)
delay = randi([1 100]);% Random delay (in samples)

%% Transmitter

% Generate random input data bits
dataBits = randi([0, 1], numBits, 1);

% QAM Mapping
numSymbols = numBits / bitsPerSymbol;
qamSymbols = zeros(numSymbols, 1);
for i = 1:numSymbols
    bits = dataBits((i-1)*bitsPerSymbol + (1:bitsPerSymbol))';
    grayIndex = bi2de(bits, 'left-msb');
    qamSymbols(i) = gray_to_qam(grayIndex, M);
end

% OFDM Modulation
Data_tx = reshape(qamSymbols, N, []);  % Reshape to N subcarriers x symbols
Data_tx = ifft_function(Data_tx, N);   % Perform IFFT
Data_tx = Data_tx * sqrt(N);           % Normalize
Data_tx = [Data_tx(end-CP+1:end, :); Data_tx];  % Add Cyclic Prefix
Data_tx = Data_tx(:);                 % Convert to serial stream


% Generate ZC Sequence for Synchronization
zc_signal = generate_signal(seq_length, root, cp_length);

% Combine Synchronization and Data Signals
tx_signal = [zc_signal, Data_tx.'];  % Concatenate ZC and OFDM data
% tx_signal = tx_signal(:);            % Ensure column vector

%% Channel: Add Delay, Noise, and CFO
tx_signal_delayed = [zeros(1, delay), tx_signal];  % Add delay
rx_signal = AWGN_with_CFO(tx_signal_delayed, SNR_dB, CFO, fs);  % Add noise and CFO

%% Receiver

% Synchronization: Detect ZC Sequence
[time_estimate, corr_output, cfo_estimate] = synchronization_with_CFO(rx_signal, tx_signal_delayed, cp_length, seq_length);
cfo_estimate_hz = cfo_estimate * fs;

% Extract OFDM Data after Synchronization
rx_signal_aligned = rx_signal(time_estimate + length(zc_signal):end);

% Demodulation (OFDM Receiver)
rx_signal_parallel = reshape(rx_signal_aligned, CP + N, []);  % Parallelize
rx_signal_no_cp = rx_signal_parallel(CP+1:end, :);           % Remove CP
Data_rx = fft(rx_signal_no_cp, N) / sqrt(N);                 % Perform FFT

% Plot Transmitted, Received Signal, and Correlation
my_plot(tx_signal, rx_signal, corr_output);

% Print Results
fprintf('Start of received signal detected at sample index: %d\n', time_estimate);
fprintf('Estimated CFO in Hz: %.2f Hz\n', cfo_estimate_hz);

%% Methods from QAM

% Function for IFFT calculation
function timeDomainSymbols = ifft_function(Data_tx, N_sub)
    [Nrows, Ncols] = size(Data_tx);

    timeDomainSymbols = zeros(N_sub, Ncols);
    for col = 1:Ncols
        X = Data_tx(:, col);
        N = length(X);
        x = zeros(N, 1);
        for n = 0:N-1
            for k = 0:N-1
                x(n+1) = x(n+1) + X(k+1) * exp(1i * 2 * pi * k * n / N);
            end
        end
        timeDomainSymbols(:, col) = x / N;
    end

end

% Function for FFT calculation
function frequencyDomainSymbols = fft_function(Data_rx, N_sub)
    [Nrows, Ncols] = size(Data_rx);

  frequencyDomainSymbols = zeros(N_sub, Ncols);
    for col = 1:Ncols
        x = Data_rx(:, col);
        N = length(x);
        X = zeros(N, 1);
        for k = 0:N-1
            for n = 0:N-1
                X(k+1) = X(k+1) + x(n+1) * exp(-1i * 2 * pi * k * n / N);
            end
        end
        frequencyDomainSymbols(:, col) = X;
    end
end

%Function for Gray to QAM Mapping (16-QAM)
function qamSymbol = gray_to_qam(grayIndex, M)
    if M == 16
        % Predefined QAM constellation points for 16-QAM
        constellation = [-3 -1 1 3]/ sqrt(10);
        realPart = constellation(mod(grayIndex, 4) + 1);
        imagPart = constellation(floor(grayIndex / 4) + 1);
        qamSymbol = (realPart + 1i * imagPart) ;
    else
        error('Mapping not implemented for this value of M');
    end
end

% Function for QAM to Gray Demapping (16-QAM)
function bits = qam_to_gray(qamSymbol, M)
    if M == 16
        constellation = [-3 -1 1 3];  % Original constellation points
        [~, realIdx] = min(abs(real(qamSymbol*sqrt(10)) - constellation));
        [~, imagIdx] = min(abs(imag(qamSymbol*sqrt(10)) - constellation));
        grayIndex = (imagIdx - 1) * 4 + (realIdx - 1);
        bits = de2bi(grayIndex, log2(M), 'left-msb')';
    else
        error('Demapping not implemented for this value of M');
    end
end