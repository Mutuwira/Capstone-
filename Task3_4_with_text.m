clc
clear
close all

%% Parameters

% OFDM Parameters
N = 64;                % Number of OFDM subcarriers
M = 16;                % 16-QAM
bitsPerSymbol = log2(M);
numSymbols = 1;        % Number of OFDM symbols

% Pilot Parameters
pilotSpacing = 4;
pilotIndices = 1:pilotSpacing:N;
dataIndices = setdiff(1:N, pilotIndices);
numDataSubcarriers = length(dataIndices);

% Synchronization Parameters
seq_length = 128;       % Length of ZC sequence
root = 1;              % Root of ZC sequence
SNR_dB =30 ;           % Signal-to-noise ratio (in dB)

CFO = randi([1 20]);   % Carrier Frequency Offset (in Hz)
fs = 10000;            % Sampling frequency (in Hz)
delay = randi([1 100]);% Random delay (in samples)

CP = 0.25*N;        % Cyclic prefix length
%CP = max(delay)+1; % if  we make the CP greater than the max delay spread of the channel

h = [1 ];  % Channel response
%delay = 0;
%h = 1;

% Convert "text" to bits
textMessage = ['Research Proposal, demonstrating anticipated research.'];
fprintf('Transmitted Text: %s \n', textMessage);
dataBits = reshape(de2bi(uint8(textMessage), 8, 'left-msb')', [], 1);

% Pad bits to match OFDM structure
numBits = length(dataBits);
numDataSymbols = ceil(numBits / bitsPerSymbol);
numDataSymbols = ceil(numDataSymbols / numDataSubcarriers) * numDataSubcarriers;
numBitsPadded = numDataSymbols * bitsPerSymbol;
dataBits = [dataBits; zeros(numBitsPadded - numBits, 1)];

%% Transmitter

% QAM Mapping
qamSymbols = zeros(numDataSymbols, 1);
for i = 1:numDataSymbols
    bits = dataBits((i-1)*bitsPerSymbol + (1:bitsPerSymbol))';
    grayIndex = bi2de(bits, 'left-msb');
    qamSymbols(i) = gray_to_qam(grayIndex, M);
end

% Map data symbols to data subcarriers
numOFDMsymbols = numDataSymbols / numDataSubcarriers;
Data_tx = zeros(N, numOFDMsymbols);
Data_tx(dataIndices, :) = reshape(qamSymbols, numDataSubcarriers, numOFDMsymbols);

% Insert pilots
Data_tx_with_pilots = insert_pilots(Data_tx, pilotSpacing);

% OFDM Modulation
Data_tx_ifft = ifft_function(Data_tx_with_pilots, N);
Data_tx_cp = [Data_tx_ifft(end-CP+1:end, :); Data_tx_ifft];

% Combine Synchronization and Data
zc_signal = generate_signal(seq_length, root, CP);
tx_signal = [zc_signal, Data_tx_cp(:).'];

%% Channel (AWGN and Multipath)

% Simulate channel effects
rx_signal = channelEmulation(tx_signal, SNR_dB, delay, h);

% Apply CFO
t = (0:length(rx_signal)-1) / fs;  % Time vector
cfo_phase = exp(1i * 2 * pi * CFO * t);  % CFO-induced phase rotation
rx_signal = rx_signal .* cfo_phase;

%% Receiver
% Synchronization: Detect ZC Sequence
[time_estimate, corr_output, cfo_estimate] = synchronization_with_CFO(rx_signal, tx_signal, CP, seq_length);
cfo_corrected_signal = rx_signal .* exp(-1i * 2 * pi * cfo_estimate * (0:length(rx_signal)-1));

% Extract OFDM Data
rx_signal_aligned = cfo_corrected_signal(time_estimate + length(zc_signal):end);
Data_rx_cp = reshape(rx_signal_aligned(1:numOFDMsymbols*(N+CP)), N+CP, numOFDMsymbols);
Data_rx = Data_rx_cp(CP+1:end, :);
Data_rx_fft = fft_function(Data_rx, N);

% Channel Estimation 
signalPower = mean(abs(Data_tx_with_pilots(:)).^2);  % Signal power
noiseVariance = signalPower / (10^(SNR_dB / 10));  % Noise variance
H_LS = ls_channel_estimation(Data_rx_fft, Data_tx_with_pilots, pilotSpacing, noiseVariance, signalPower);

% Equalization
equalizedSymbols = mmse_equalization(Data_rx_fft, H_LS, noiseVariance, signalPower);

% Remove pilots from equalized symbols
equalizedDataSymbols = equalizedSymbols(dataIndices, :);
equalizedDataSymbols_serial = equalizedDataSymbols(:);

% QAM Demapping
receivedBits = zeros(numBits, 1);
for i = 1:numDataSymbols
    bits = qam_to_gray(equalizedDataSymbols_serial(i), M);
    receivedBits((i-1)*bitsPerSymbol + (1:bitsPerSymbol)) = bits';
end

% Convert bits back to text
receivedText = char(bi2de(reshape(receivedBits(1:numBits), 8, []).', 'left-msb'))';
disp(['Received Text: ', receivedText]);

%% Results

% Plot Transmitted, Received Signal, and Correlation
my_plot(tx_signal, rx_signal, corr_output);

% Calculate Bit Error Rate (BER)
numErrors = sum(dataBits ~= receivedBits);  % Count bit errors
ber = numErrors / numBits;  % Compute BER
fprintf('Bit Error Rate (BER): %f\n', ber);

% Constellation Plots
scatterplot(qamSymbols);  % Original transmitted symbols
title('Transmitted QAM Constellation');

scatterplot(equalizedDataSymbols_serial(:));  % Received symbols after equalization
title('Received Constellation after Equalization');

% Visualizing BER vs SNR
SNR_range = 0:2:40;
BER = zeros(length(SNR_range), 1);

for i = 1:length(SNR_range)
    SNR = SNR_range(i); % Use current SNR
    
    % Simulate channel
    rx_signal = channelEmulation(tx_signal, SNR, delay, h);
    
    % Apply CFO
    t = (0:length(rx_signal)-1) / fs;
    cfo_phase = exp(1i * 2 * pi * CFO * t);
    rx_signal = rx_signal .* cfo_phase;

    % Synchronization
    [time_estimate, ~, cfo_estimate] = synchronization_with_CFO(rx_signal, tx_signal, CP, seq_length);
    cfo_corrected_signal = rx_signal .* exp(-1i * 2 * pi * cfo_estimate * (0:length(rx_signal)-1));

    % Extract OFDM Data
    rx_signal_aligned = cfo_corrected_signal(time_estimate + length(zc_signal):end);
    
    if length(rx_signal_aligned) < numOFDMsymbols * (N+CP)
        continue; % Skip this SNR point if reshaping fails
    end
    
    Data_rx_cp = reshape(rx_signal_aligned(1:numOFDMsymbols*(N+CP)), N+CP, numOFDMsymbols);
    Data_rx = Data_rx_cp(CP+1:end, :);
    Data_rx_fft = fft_function(Data_rx, N);
    
    % Channel Estimation
    signalPower = mean(abs(Data_tx_with_pilots(:)).^2);
    noiseVariance = signalPower / (10^(SNR / 10)); 
    H_LS = ls_channel_estimation(Data_rx_fft, Data_tx_with_pilots, pilotSpacing, noiseVariance, signalPower);
    
    % MMSE Equalization 
    equalizedSymbols = mmse_equalization(Data_rx_fft, H_LS, noiseVariance, signalPower);
    
    % Remove pilots
    equalizedDataSymbols = equalizedSymbols(dataIndices, :);
    equalizedDataSymbols_serial = equalizedDataSymbols(:);
    
    % QAM Demapping
    receivedBits = zeros(numBits, 1);
    for j = 1:numDataSymbols
        bits = qam_to_gray(equalizedDataSymbols_serial(j), M);
        receivedBits((j-1)*bitsPerSymbol + (1:bitsPerSymbol)) = bits';
    end
    
    % Compute BER
    numErrors = sum(dataBits ~= receivedBits);
    BER(i) = numErrors / numBits;
end

% Plot BER vs. SNR
figure;
plot(SNR_range, BER, 'X - ', 'LineWidth', 2);
set(gca, 'YScale', 'log') % Correct logarithmic scale
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
title('BER vs. SNR');
grid on;

%% Functions for OFDM

% LS estimation
function H_LS = ls_channel_estimation(Data_rx_fft, Data_tx_with_pilots, pilotSpacing, noiseVariance, signalPower )
    [N, numSymbols] = size(Data_rx_fft);
    pilotIndices = 1:pilotSpacing:N;


    H_LS = zeros(N, numSymbols);
    for sym = 1:numSymbols
        % Extract received pilot symbols
        Y_pilots = Data_rx_fft(pilotIndices, sym);
        % Extract transmitted pilot symbols
        X_pilots = Data_tx_with_pilots(pilotIndices, sym);
        % LS estimation for pilot subcarriers
        H_pilots = Y_pilots ./ X_pilots;
        % Interpolate channel estimates for data subcarriers
        H_LS(:, sym) = interp1(pilotIndices, H_pilots, (1:N)', 'linear', 'extrap');
    end

    % Apply MMSE correction
    H_LS = H_LS .* (abs(H_LS).^2 ./ (abs(H_LS).^2 + (noiseVariance / signalPower)));

end


% Function to Insert Pilots
function Data_tx_with_pilots = insert_pilots(Data_tx, pilotSpacing)
    pilotValue = 1 + 1i; % Known pilot value
    [N, numSymbols] = size(Data_tx);
    Data_tx_with_pilots = Data_tx;
    pilotIndices = 1:pilotSpacing:N;
    for sym = 1:numSymbols
        Data_tx_with_pilots(pilotIndices, sym) = pilotValue; % Insert pilots
    end
end

function y = channelEmulation(x, snr_dB, to, h)
    % Add delay
    y_delayed = [zeros(1, to), x]; % Create offset
    y_channel = conv(h, y_delayed); % Convolve with channel response
    signalPower = mean(abs(y_channel).^2); % Compute signal power
    noisePower = signalPower / (10^(snr_dB/10));  % Compute noise power based on SNR
    
    % Generate AWGN noise
    noise = sqrt(noisePower/2) .* (randn(size(y_channel)) + 1i .* randn(size(y_channel)));
    y = y_channel + noise;     % Add noise to the signal
end


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
        timeDomainSymbols(:, col) = x/N ;
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

function qamSymbol = gray_to_qam(grayIndex, M)
    if M == 16
        % Gray code mapping for 16-QAM
        mapping = [
            -3 + 3i, -3 + 1i, -3 - 3i, -3 - 1i, ...
            -1 + 3i, -1 + 1i, -1 - 3i, -1 - 1i, ...
             3 + 3i,  3 + 1i,  3 - 3i,  3 - 1i, ...
             1 + 3i,  1 + 1i,  1 - 3i,  1 - 1i
        ] / sqrt(10);
        qamSymbol = mapping(grayIndex + 1);
    else
        error('Mapping not implemented for this value of M');
    end
end

function bits = qam_to_gray(qamSymbol, M)
    if M == 16
        % Gray code mapping for 16-QAM
        mapping = [
            -3 + 3i, -3 + 1i, -3 - 3i, -3 - 1i, ...
            -1 + 3i, -1 + 1i, -1 - 3i, -1 - 1i, ...
             3 + 3i,  3 + 1i,  3 - 3i,  3 - 1i, ...
             1 + 3i,  1 + 1i,  1 - 3i,  1 - 1i
        ] / sqrt(10);
        [~, idx] = min(abs(qamSymbol - mapping));
        grayIndex = idx - 1;
        bits = de2bi(grayIndex, log2(M), 'left-msb')';
    else
        error('Demapping not implemented for this value of M');
    end
end

% Functions for synchronization

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

function tx_signal_plain = generate_signal(seq_length, root, cp_length) % Function to generate ZC sequence with CP and data
    % Generate Zadoff-Chu (ZC) sequence
    n = 0:seq_length-1;
    zc_sequence = exp(-1i*pi*root*n.*n/seq_length);  % ZC sequence formula
    % Add Cyclic Prefix (CP) to ZC sequence
    zc_with_cp = [zc_sequence(end-cp_length+1:end), zc_sequence];

    %data_val = randi([1, 100]); 
    %random_data = (randn(1, data_val) + 1i * randn(1, data_val)) / sqrt(2);  % complex data
    tx_signal_plain = zc_with_cp; % ZC with cp and random length zero data
end

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


function my_plot(tx_signal, rx_signal, corr_output)
    figure;
    subplot(3,1,1);
    plot(real(tx_signal));
    title('Transmitted ZC Signal');
    xlabel('Sample Index');
    ylabel('Amplitude');
    
    subplot(3,1,2);
    plot(real(rx_signal));
    title('Received ZC Signal');
    xlabel('Sample Index');
    ylabel('Amplitude');
    
    subplot(3,1,3);
    plot(corr_output);
    title('Cross-Correlation Output for Synchronization');
    xlabel('Lag');
    ylabel('Correlation Magnitude');
end



