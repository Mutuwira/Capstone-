clc
clear
close all
rng(12345); % Set random seed for reproducibility

%% Convert text to binary representation
[data, ~] = text2bin('This Capstone Sucks very bad. But I am still doing it. That is fair or unfair.');

%% SDR Transmitter Configuration
transmit = comm.SDRuTransmitter(...
              'Platform','B200', ... % Specify USRP hardware
              'SerialNum','329AEA7', ... % Define specific serial number
              'ChannelMapping',1, ... % Use first channel
              'Gain', 50, ... % Set transmission gain
              'CenterFrequency', 2.45e9, ... % Set frequency to 2.45 GHz
              'InterpolationFactor', 64); % Define interpolation factor

disp(transmit) % Display transmitter settings

%% Transmitter Parameters
N = 2; % QPSK modulation order
FFT_N = 64; % FFT length for OFDM
sym_len = 90; % Length of OFDM symbols
CP_size = 16; % Cyclic prefix length
PL = 8; % Number of payload symbols
pilot = 1; % Pilot signal

% Ensure data length matches payload length by padding if necessary
data = [data zeros(1, PL * sym_len - length(data))];  

tx = []; % Initialize transmission signal buffer

%% OFDM Symbol Generation
for i = 1:PL
    sym = []; % Initialize symbol buffer
    for j = 1:6:sym_len
        m = mapper(data((i - 1)*sym_len +j:(i - 1)*sym_len + j + 5));
        sym = [sym m];
        
        % Insert pilot symbols at specified positions
        if ismember(length(sym), [12 21 42 55])
            sym = [sym pilot];
        end
    end
    
    % Convert frequency domain symbols to time domain using IFFT
    time = ifft(sym, FFT_N);
    
    % Add cyclic prefix and append to transmission signal
    tx = [tx time(FFT_N - CP_size + 1:end) time];
end

%% Generate Preamble
preamble = 1 - 2*randi([0 1], FFT_N/2, 1)' + 1i - 1i*2*randi([0 1], FFT_N/2, 1)';
preamble_time = ifft(preamble, FFT_N/2); % Convert to time domain
preamble_time = [preamble_time preamble_time]; % Duplicate for robustness
preamble_cp = [preamble_time(FFT_N - CP_size + 1:end) preamble_time];

% Append preamble to transmission signal
tx = [preamble_cp/2 tx];
preamble = fft(preamble_time, FFT_N); % Convert back to frequency domain

%% Generate Additional Preamble for Synchronization
preamble_2 = 1 - 2*randi([0 1], FFT_N, 1)' + 1i - 1i*2*randi([0 1], FFT_N, 1)';
preamble_2([2:4:end]) = 0 + 1i*0;
preamble_2([3:4:end]) = 0 + 1i*0;
preamble_2([4:4:end]) = 0 + 1i*0;
preamble_time_2 = ifft(preamble_2, FFT_N);
preamble_cp_2 = [preamble_time_2(FFT_N - CP_size + 1:end) preamble_time_2];

tx = [preamble_cp_2*sqrt(2) tx]; % Add second preamble to transmission signal

%% Transmit Continuously
while (true)
    transmit(tx'); % Send the signal to the SDR
end

release(transmit); % Release the SDR transmitter

%% Function: QPSK Mapper
function m = mapper(sub_data)
    sig = 1 - 2*sub_data(3:6);
    s1 = sig(2) + 1i*sig(1);
    s2 = sig(4) + 1i*sig(3);
    
    % Map input bits to QPSK symbols with zero padding
    if sub_data(1:2) == [0 0]
        m = [s1 s2 0 0];
    elseif sub_data(1:2) == [0 1]
        m = [0 s1 s2 0];
    elseif sub_data(1:2) == [1 0]
        m = [0 0 s1 s2];
    elseif sub_data(1:2) == [1 1]
        m = [s1 0 0 s2];
    end
end

%% Function: Convert Text to Binary
function [binV, binS] = text2bin(text)
    binS = dec2bin(text,8); % Convert characters to 8-bit binary
    binS = binS';
    binS = binS(:)'; % Flatten to a single row
    binV = binS - 48; % Convert ASCII to numeric binary values
end
