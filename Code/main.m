%% Main code structure

%{
- main.m                     % Main script: Sets parameters and calls functions
- bit_padding.m              % To Match OFDM Structure
- qam_modulation.m           % QAM mapping
- qam_demodulation.m         % QAM demapping
- equalization.m             % MMSE equalization
- ifft_function.m            % Inverse Fast Fourier Transform 
- fft_function.m             % Fast Fourier Transform 
- generate_signal.m          % Function to generate ZC sequence with CP and data
- synchronization.m          % Synchronization using ZC sequence
- channel_estimation.m       % Channel estimation using preamble
- channelEmulation.m         % Channel simulation (AWGN, multipath)
- display_results.m          % Display results and plots
- decoder.m                  % Systematic Linear dencoder
- encoder.m                  % Systematic Linear encoder
%}

clc
clear
close all

%% Parameters

% OFDM Parameters
N = 64;                % Number of OFDM subcarriers (FFT Size)
M = 16;                % QAM order
bitsPerSymbol = log2(M);

% Linear Encoder Parameters (Code rate = k/n)
n = 8; % Codeword length
k = 4; % Message length
genMatrix = [eye(k), randi([0, 1], k, n - k)]; % Generator matrix


% Preamble and Synchronization Parameters
preamble_length = 128;  % Length of preamble sequence
seq_length = 128;       % Length of ZC sequence (used in preamble)
root = 1;              % Root of ZC sequence
SNR_dB = 30;           % Signal-to-noise ratio (in dB)
CFO = randi([1 20]);   % Carrier Frequency Offset (in Hz)
fs = 10000;            % Sampling frequency (in Hz)
delay = randi([1 100]);% Random delay (in samples)
CP = 0.25*N;           % Cyclic prefix length

% Channel Response
h = [1+1i*0.5 0.2+1i*0.6 0.1+1i*0.4 0.05+1i*0.01]; 
%h = 1;

% Input Transmitted
textMessage = ['Hello Kundai'];
fprintf('Transmitted Text: %s \n', textMessage);
dataBits = reshape(de2bi(uint8(textMessage), 8, 'left-msb')', [], 1);
numBits = length(dataBits);


%% Transmitter

% Linear encoding
encodedBits = encoder(dataBits, n, k, genMatrix);

% Bit Padding to match OFDM structure
[encodedBits] = bit_padding(encodedBits, bitsPerSymbol, N);

% QAM Modulation
% qamSymbols = qam_modulation(encodedBits, bitsPerSymbol, N, M);
qamSymbols = qammod(encodedBits,M,'InputType','bit','UnitAveragePower',true);

% OFDM Modulation
numofdmsymbols = length(qamSymbols)/N;
Data_tx = reshape(qamSymbols, N, numofdmsymbols);
Data_tx_ifft = ifft_function(Data_tx, N);
Data_tx_cp = [Data_tx_ifft(end-CP+1:end, :); Data_tx_ifft];

% Generate preamble, combine synchronization and Data
zc_signal = generate_signal(seq_length, root, CP);
tx_signal = [zc_signal, Data_tx_cp(:).'];

%% Channel (AWGN and Multipath)

rx_signal = channelEmulation(tx_signal, SNR_dB, delay, h);

% Apply CFO
t = (0:length(rx_signal)-1) / fs;  % Time vector
cfo_phase = exp(1i * 2 * pi * CFO * t);  % CFO-induced phase rotation
rx_signal = rx_signal .* cfo_phase;

%% Receiver

% Synchronization: Detect ZC Sequence
[time_estimate, corr_output, cfo_estimate] = synchronization(rx_signal, zc_signal, CP, seq_length);
cfo_corrected_signal = rx_signal .* exp(-1i * 2 * pi * cfo_estimate * (0:length(rx_signal)-1));

% Extracting OFDM Data
rx_signal_aligned = cfo_corrected_signal(time_estimate + length(zc_signal):end);
% Data_rx_cp = reshape(rx_signal_aligned(1:numDataSymbols*(N+CP)/N), N+CP, numDataSymbols/N);
Data_rx_cp = reshape(rx_signal_aligned(1:end-(length(h)-1)), N+CP, numofdmsymbols);
Data_rx = Data_rx_cp(CP+1:end, :);
Data_rx_fft = fft_function(Data_rx, N);

% Channel Estimation
signalPower = mean(abs(tx_signal(:)).^2);  % Signal power
noiseVariance = signalPower / (10^(SNR_dB / 10));  % Noise variance
preamble_rx = cfo_corrected_signal(time_estimate : time_estimate + length(zc_signal) - 1);
H_est = channel_estimation(preamble_rx, zc_signal, seq_length, N, noiseVariance, signalPower);

% MMSE Equalization
equalizedSymbols = equalization(Data_rx_fft, H_est, signalPower, noiseVariance);

% QAM Demodulation gives you receivedEncodedBits:
receivedEncodedBits = qam_demodulation(equalizedSymbols, bitsPerSymbol, N, M);

% Decoding using our linear systematic decoder:
decodedBits = decoder(receivedEncodedBits, n, k, genMatrix);

% Extract the original bits (if padding was added, use the first numBits bits).
receivedText = char(bi2de(reshape(decodedBits(1:numBits), 8, []).', 'left-msb'))';
disp(['Received Text: ', receivedText]);

% Extract only the originally transmitted bits from the decoded stream
decodedBits_truncated = decodedBits(1:numBits);

% Calculate Bit Error Rate (BER)
numErrors = sum(dataBits ~= decodedBits_truncated);  % Count bit errors
ber = numErrors / numBits;  % Compute BER
fprintf('Bit Error Rate (BER): %f\n', ber);

% Displaying results
display_results(tx_signal, rx_signal, corr_output, qamSymbols, equalizedSymbols, SNR_dB);
