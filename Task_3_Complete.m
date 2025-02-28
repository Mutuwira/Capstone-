clc
clear
close all
%% Transmitter

% Random seed for reproducibility
rng(42); 

% Parameters
N = 64;  % Number of OFDM subcarriers
M = 16;  % 16-QAM
numBits = N * log2(M);  % Total number of bits
bitsPerSymbol = log2(M);  % Bits per symbol (4 bits for 16-QAM)
CP = N * 0.25;


% Roberto' Parameters
SNR = 25; 
to = 0;
h = [1 0.2 0.1 0.05];

%h = [1];

% Input data
dataBits = randi([0, 1], numBits, 1);  % Random bit stream

% QAM Mapping (Gray Mapping for 16-QAM)
numSymbols = numBits / bitsPerSymbol;  % Number of QAM symbols
qamSymbols = zeros(numSymbols, 1);  % Initialize QAM symbols

for i = 1:numSymbols
    bits = dataBits((i-1)*bitsPerSymbol + (1:bitsPerSymbol))';  % Extract 4 bits
    grayIndex = bi2de(bits, 'left-msb');  
    qamSymbols(i) = gray_to_qam(grayIndex, M); 
end

% OFDM Modulation (TX): Serial to Parallel, IFFT
Data_tx = reshape(qamSymbols, N, []);  % Reshape to N subcarriers x symbols
Data_tx = ifft_function(Data_tx, N)* sqrt(N); % IFFT & Normalization

fprintf('Power before CP: %f\n', mean(abs(Data_tx(:)).^2));

% Add Cyclic Prefix (CP)
Data_tx = [Data_tx(end-CP+1:end, :); Data_tx];
fprintf('Power after CP: %f\n', mean(abs(Data_tx(:)).^2));

Data_tx = Data_tx(:);  % Convert to serial for transmission

% Pass through the channel (provided function)
y = channelEmulation(Data_tx, 10.^(SNR/10), to, h).*sqrt(10.^(-SNR/10));

%% Receiver

extraSamples = length(h) - 1; % Due to convolution in channelEmulation
Data_rx = reshape(y, N + extraSamples + CP , []);  % Serial to parallel conversion
Data_rx = Data_rx(CP+1: CP + N, :); % Remove CP
Data_rx = (fft_function(Data_rx, N))/ sqrt(N); % FFT & Normalization

fprintf('Received Power before EQ: %f\n', mean(abs(Data_rx(:)).^2));

% Equalization Using Known Channel (h_hat)
h_hat = fft(h, N) / sqrt(N); % estimation


signalPower = mean(abs(Data_tx(:)).^2);  % Signal power
noiseVariance = signalPower / (10^(SNR / 10));  % Noise variance

equalizedSymbols = mmse_equalization(Data_rx, h_hat, noiseVariance, signalPower);
equalizedSymbols = equalizedSymbols / sqrt(mean(abs(equalizedSymbols).^2)); % normalization
fprintf('Received Power after EQ: %f\n', mean(abs(Data_rx(:)).^2));

% QAM Demapping: Convert Equalized Symbols Back to Bits
receivedBits = zeros(numBits, 1);
for i = 1:numSymbols
    receivedBits((i-1)*bitsPerSymbol + (1:bitsPerSymbol)) = qam_to_gray(equalizedSymbols(i), M);
end

% Calculate Bit Error Rate (BER)
numErrors = sum(dataBits ~= receivedBits);  % Count bit errors
ber = numErrors / numBits;  % Compute BER
fprintf('Bit Error Rate (BER): %f\n', ber);

% Constellation Plots
scatterplot(qamSymbols);  % Original transmitted symbols
title('Transmitted QAM Constellation');

scatterplot(equalizedSymbols(:));  % Received symbols after equalization
title('Received Constellation after Equalization');

% Function for MMSE equalization
function equalizedSymbols = mmse_equalization(Data_rx, H, noiseVariance, signalPower)
    mmseFilter = conj(H) ./ (abs(H).^2 + noiseVariance / signalPower);
    equalizedSymbols = Data_rx .* mmseFilter;
end

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


% Visualizing BER vs SNR

SNR_range = 0:5:40; % Range of SNR values
BER = zeros(length(SNR_range), 1); % Store BER for each SNR

for i = 1:length(SNR_range)
    SNR = SNR_range(i); % Set current SNR
    y = channelEmulation(Data_tx, 10^(SNR/10), to, h); % Channel
    y = y .* sqrt(10^(-SNR/10)); % Scale by noise power

    % Receiver Processing (reuse your existing receiver code here)
    Data_rx = reshape(y, N + CP + extraSamples, []);
    Data_rx = Data_rx(CP + 1:CP + N, :);
    Data_rx = fft_function(Data_rx, N) / sqrt(N);

    % Equalization
    equalizedSymbols = mmse_equalization(Data_rx, h_hat, noiseVariance, signalPower);
    equalizedSymbols = equalizedSymbols / sqrt(mean(abs(equalizedSymbols).^2)); % normalization
    
    % QAM Demapping and BER Calculation
    receivedBits = zeros(numBits, 1);
    for j = 1:numSymbols
        receivedBits((j-1)*bitsPerSymbol + (1:bitsPerSymbol)) = qam_to_gray(equalizedSymbols(j), M);
    end
    numErrors = sum(dataBits ~= receivedBits);
    BER(i) = numErrors / numBits; % Store BER
end

% Plot BER vs. SNR
figure;
plot(SNR_range, BER, 'X - ', LineWidth=2);
yscale("log")
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
title('BER vs. SNR');
grid on;

