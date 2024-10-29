clc;
clear;
close all;

% Parameters
SNR_range = 0:5:30;  % SNR values for BER vs. SNR plot
QAM_Orders = [16];  % QAM orders: 4-QAM, 16-QAM, 64-QAM
Equalizer_Type = '1';  % MMSE equalizer

% Initialize storage for BER results
BER_normalized = zeros(length(SNR_range), length(QAM_Orders));
BER_non_normalized = zeros(length(SNR_range), length(QAM_Orders));

% Simulation Loop for Each QAM Order
for q = 1:length(QAM_Orders)
    qam_order = QAM_Orders(q);  % Current QAM order

    for s = 1:length(SNR_range)
        SNR = SNR_range(s);  % Current SNR value
        h = [1 0.2 0.1 0.05];
        % Run simulation for normalized and non-normalized cases
        BER_normalized(s, q) = simulate_ofdm(qam_order, Equalizer_Type, true, 64, 32, SNR, h, 0);
        BER_non_normalized(s, q) = simulate_ofdm(qam_order, Equalizer_Type, false, 64, 32, SNR, h, 0);
    end
end

% Plotting BER vs. SNR for All QAM Orders
figure;
for q = 1:length(QAM_Orders)
    qam_order = QAM_Orders(q);

    % Plot for normalized BER
    %subplot(2, 1, q);
    plot(SNR_range, BER_normalized(:, q), 'o-', 'LineWidth', 2, 'DisplayName', 'Normalized');
    hold on;
    
    % Plot for non-normalized BER
    %plot(SNR_range, BER_non_normalized(:, q), 'x--', 'LineWidth', 2, 'DisplayName', 'Non-Normalized');
    
    % Plot Formatting
    xlabel('SNR (dB)');
    ylabel('Bit Error Rate (BER)');
    title(sprintf('BER vs. SNR for %d-QAM', qam_order));
    %legend('show');
    grid on;
end

% Function to Simulate OFDM for Given Parameters
function BER = simulate_ofdm(qam_order, eq_type, normalize, N, CP, SNR, h, to)
    bitsPerSymbol = log2(qam_order);  % Bits per QAM symbol
    numBits = N * bitsPerSymbol;  % Total bits per OFDM symbol
    dataBits = randi([0, 1], numBits, 1);  % Generate random bit stream

    % QAM Mapping
    qamSymbols = gray_map(dataBits, qam_order);

    % OFDM Modulation with IFFT
    Data_tx = reshape(qamSymbols, N, []);
    Data_tx = ifft(Data_tx) * sqrt(N) * normalize;

    % Add Cyclic Prefix (CP)
    Data_tx = [Data_tx(end-CP+1:end, :); Data_tx];
    Data_tx = Data_tx(:);  % Convert to serial for transmission

    % Channel + Noise
    y = channelEmulation(Data_tx, 10^(SNR/10), to, h) .* sqrt(10^(-SNR/10));

    % Receiver: Remove CP and apply FFT
    extraSamples = length(h) - 1;  % Due to channel convolution
    Data_rx = reshape(y, N + CP + extraSamples, []);
    Data_rx = Data_rx(CP+1:CP+N, :);  % Remove CP
    Data_rx = fft(Data_rx) / (sqrt(N) * normalize);  % Apply FFT and normalize

    % Equalization
    h_hat = fft(h, N) / sqrt(N);  % Channel estimate
    signalPower = mean(abs(Data_tx(:)).^2);  % Signal power
    noiseVariance = signalPower / (10^(SNR / 10));  % Noise variance
    equalizedSymbols = equalize(Data_rx, h_hat, eq_type, noiseVariance, signalPower);

    % QAM Demapping and BER Calculation
    receivedBits = qam_demap(equalizedSymbols, qam_order);
    BER = sum(dataBits ~= receivedBits) / numBits;  % Compute BER
end

% QAM Mapping (Gray Encoding)
function qamSymbols = gray_map(bits, M)
    numSymbols = length(bits) / log2(M);
    qamSymbols = zeros(numSymbols, 1);

    for i = 1:numSymbols
        bitGroup = bits((i-1)*log2(M) + (1:log2(M)))';
        grayIndex = bi2de(bitGroup, 'left-msb');
        qamSymbols(i) = gray_to_qam(grayIndex, M);
    end
end

% QAM Demapping (Gray Decoding)
function bits = qam_demap(qamSymbols, M)
    bits = [];
    for sym = qamSymbols.'
        bits = [bits; qam_to_gray(sym, M)];
    end
end

% Gray to QAM Mapping Function
function qamSymbol = gray_to_qam(grayIndex, M)
    constellation = qammod(0:M-1, M, 'gray');  % Gray mapping
    constellation = constellation / sqrt(mean(abs(constellation).^2));  % Normalize
    qamSymbol = constellation(grayIndex + 1);  % Map to QAM symbol
end

% QAM to Gray Demapping Function
function bits = qam_to_gray(qamSymbol, M)
    constellation = qammod(0:M-1, M, 'gray');  % Gray mapping
    constellation = constellation / sqrt(mean(abs(constellation).^2));  % Normalize
    [~, idx] = min(abs(constellation - qamSymbol));  % Find closest symbol
    bits = de2bi(idx - 1, log2(M), 'left-msb')';  % Convert to bit vector
end

% Equalizer Function (ZF or MMSE)
function equalizedSymbols = equalize(Data_rx, h_hat, eq_type, noiseVariance, signalPower)
    switch eq_type
        case '1'  % Zero Forcing Equalizer
            equalizedSymbols = Data_rx ./ h_hat;
        case '2'  % MMSE Equalizer
            mmseFilter = conj(h_hat) ./ (abs(h_hat).^2 + noiseVariance / signalPower);
            equalizedSymbols = Data_rx .* mmseFilter;
        otherwise
            error('Invalid Equalizer Type');
    end

    % Normalize power after equalization
    equalizedSymbols = equalizedSymbols / sqrt(mean(abs(equalizedSymbols).^2));
end
