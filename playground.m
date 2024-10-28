clc;
clear;
close all;

%% Parameters
N = 64;  % Number of OFDM subcarriers
CP = N * 0.50;  % Cyclic prefix length
SNR = 0:5:40;  % SNR values for BER vs. SNR plot
h = [1 0.2 0.1 0.05];  % Channel impulse response
extraSamples = length(h) - 1;  % Extra samples from channel convolution

% QAM orders to test
QAM_Orders = [4, 16, 64];  % 4-QAM, 16-QAM, 64-QAM
Equalizers = {'1', '2'};  % Zero Forcing (1) and MMSE equalizers(2)

% Simulation across QAM Orders and Equalizers
for qam_order = QAM_Orders
    for eq_type = Equalizers
        for SNR_range = SNR
            eq_type = char(eq_type); 
            N = 64;  % Number of OFDM subcarriers
            CP = N * 0.50;  % Cyclic prefix length
            h = [1 0.2 0.1 0.05];  % Channel impulse response
            extraSamples = length(h) - 1;  % Extra samples from channel convolution
            to = 0;
    
            figure;
            BER_normalized = simulate_ofdm(qam_order, eq_type, true, N, CP, SNR_range, h, to);
            BER_non_normalized = simulate_ofdm(qam_order, eq_type, false, N, CP, SNR_range, h, to);
    
            semilogy(SNR_range, BER_normalized, 'o-', 'DisplayName', 'Normalized');
            hold on;
            semilogy(SNR_range, BER_non_normalized, 'x--', 'DisplayName', 'Non-Normalized');
            xlabel('SNR (dB)');
            ylabel('Bit Error Rate (BER)');
            title(sprintf('BER vs. SNR for %d-QAM with %s Equalizer', qam_order, eq_type));
            legend('show');
            grid on;
        end
    end
end

% Function to simulate OFDM for different parameters
function BER_results = simulate_ofdm(qam_order, eq_type, normalize, N, CP, SNR, h, to)
    rng(42); 
    bitsPerSymbol = log2(qam_order);  % Bits per QAM symbol
    numBits = N * bitsPerSymbol;  % Total bits per OFDM symbol
    dataBits = randi([0, 1], numBits, 1);  % Random bit stream

    % QAM Mapping
    qamSymbols = gray_map(dataBits, qam_order);

    % OFDM Modulation
    Data_tx = reshape(qamSymbols, N, []);
    Data_tx = ifft(Data_tx) * sqrt(N) * normalize;

    % Add Cyclic Prefix (CP)
    Data_tx_cp = [Data_tx(end-CP+1:end, :); Data_tx];
    Data_tx_cp = Data_tx_cp(:);  % Convert to serial for transmission

    % Channel + Noise
    y = channelEmulation(Data_tx, 10.^(SNR/10), to, h).*sqrt(10.^(-SNR/10));

    % Receiver: Remove CP and apply FFT
    extraSamples = length(h) - 1; % Due to convolution in channelEmulation
    Data_rx = reshape(y, N + extraSamples + CP , []);  % Serial to parallel conversion
    Data_rx = Data_rx(CP+1: CP + N, :); % Remove CP
    Data_rx = fft(Data_rx) / (sqrt(N) * normalize);  % FFT and normalize

    % Equalization
    h_hat = fft(h, N) / sqrt(N);  % Channel estimate
    signalPower = mean(abs(Data_tx(:)).^2);  % Signal power
    noiseVariance = signalPower / (10^(SNR / 10));  % Noise variance
    equalizedSymbols = equalize(Data_rx, h_hat,noiseVariance,signalPower );

    % QAM Demapping and BER calculation
    receivedBits = qam_demap(equalizedSymbols, qam_order);
    BER_results = sum(dataBits ~= receivedBits) / numBits;  % BER
end

% Helper Functions


% QAM Mapping (Gray Encoding)
function qamSymbols = gray_map(bits, M)
    numSymbols = length(bits) / log2(M);
    qamSymbols = zeros(numSymbols, 1);

    for i = 1:numSymbols
        bitGroup = bits((i-1) * log2(M) + (1:log2(M)))';
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
    qamSymbol = constellation(grayIndex + 1);  % Map gray index to QAM symbol
end

% QAM to Gray Demapping Function
function bits = qam_to_gray(qamSymbol, M)
    constellation = qammod(0:M-1, M, 'gray');  % Gray mapping
    constellation = constellation / sqrt(mean(abs(constellation).^2));  % Normalize
    [~, idx] = min(abs(constellation - qamSymbol));  % Closest symbol
    bits = de2bi(idx - 1, log2(M), 'left-msb')';  % Convert to bit vector
end

% Equalizer Function (ZF or MMSE)
function equalizedSymbols = equalize(Data_rx, H, eq_type)
    switch eq_type
        case '1'
            equalizedSymbols = Data_rx ./ h_hat;
        case '2'
            mmseFilter = conj(H) ./ (abs(H).^2 + noiseVariance / signalPower);
            equalizedSymbols = Data_rx .* mmseFilter;
        otherwise
            error('Invalid Equalizer Type');
    end

    % Normalize power after equalization
    equalizedSymbols = equalizedSymbols / sqrt(mean(abs(equalizedSymbols).^2));
end