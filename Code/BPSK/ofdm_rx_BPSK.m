clc
clear all
close all

rng(12345); % Set random seed for reproducibility

%% ========== 1) SDRu Receiver Setup & Data Acquisition ==========

% Create and configure the SDRu receiver object
receive = comm.SDRuReceiver( ...
    'Platform','B200', ...
    'SerialNum','329AEA8', ...
    'ChannelMapping',1, ...
    'Gain', 60, ...
    'CenterFrequency', 2.45e9, ...
    'DecimationFactor', 64, ...
    'SamplesPerFrame', 1000);
disp(receive)

% Acquire data from SDRu
rx = [];
try
    for i = 1:20
        [receivedSignal, len] = step(receive);
        % Check if valid samples are received
        if len > 0 && any(receivedSignal(:) ~= 0)
            rx = [rx; receivedSignal];
        end
    end
catch exception
    disp('error: ');
    disp(exception.message);
end

% Release SDRu receiver resource
release(receive);

% Convert to double, normalize, and store in row vector
rx = double(rx);
rx = rx' ./ max(rx); % Normalize amplitude


%% ========== 2) Parameter Definitions ==========

PL       = 8;             % Payload length (# of OFDM symbols)
N        = 2;             % Bits per sample (QPSK)
FFT_N    = 64;            % FFT size
CP_size  = FFT_N/4;       % Cyclic prefix length = 16
L        = FFT_N/2;       % Half of the FFT size (for Schmidl & Cox)
frame_len = (PL + 2) * (FFT_N + CP_size); 
xq       = 1:64;          % Used for interpolation indices

% Convert text to binary and pad to match desired payload length
[data, ~] = text2bin('This Capstone Sucks very bad. But I am still doing it. That is fair or unfair.');
data = [data zeros(1, PL*90 - length(data))];  


%% ========== 3) Generate Preamble ==========

% Create random half-FFT preamble in frequency domain
preamble = 1 - 2*randi([0 1], FFT_N/2, 1)' + ...
           1i - 1i*2*randi([0 1], FFT_N/2, 1)';
       
% Convert to time domain (half-FFT), then duplicate
preamble_time = ifft(preamble, FFT_N/2);
preamble_time = [preamble_time preamble_time]; 

% Add cyclic prefix for correlation-based sync
preamble_cp = [preamble_time(FFT_N - CP_size + 1:end), preamble_time];

% Store a frequency-domain version of the entire 64-length preamble
preamble = fft(preamble_time, FFT_N);


%% ========== 4) Coarse Time Synchronization (Schmidl & Cox) ==========

% Compute Schmidl & Cox metric: M(d) = |P(d)|^2 / (R(d))^2
%  where P(d) = sum_{m=0..L-1}[ conj(rx(d+m)) * rx(d+m+L) ]
%  and   R(d) = sum_{m=0..L-1}[ |rx(d+m+L)|^2 ]
P = zeros(length(rx), 1);
R = zeros(length(rx), 1);

for d = 1:(length(rx) - 2*L + 1)
    for m = 0:(L-1)
        P(d) = P(d) + conj(rx(d + m)) * rx(d + m + L);
        R(d) = R(d) + abs(rx(d + m + L))^2;
    end
end

M = abs(P).^2 ./ (R.^2);  % Schmidl & Cox metric

% Smooth or average M(d) with a rectangular pulse (length = CP_size) for peak detection
rectpulse = ones(CP_size, 1)';
a_pick = conv(rectpulse, M);

% Determine maximum from the convolved metric
[~, max_cf] = max(a_pick);
if (max_cf > 80) && (a_pick(max_cf - 80) >= 0.5 * a_pick(max_cf))
    max_cf = max_cf - 80;
end
disp("Max CFO index (coarse): " + max_cf);


%% ========== 5) Coarse Frequency Synchronization ==========

% Use the correlation angle at max_cf to correct CFO
df = angle(P(max_cf)) / pi;
rx = rx .* exp(-1i * 2 * pi * df * (1:length(rx)) / 64);


%% ========== 6) Fine Time Synchronization via Cross-Correlation ==========

[c, lags] = xcorr(rx, preamble_cp(1:40)); 
c = abs(c).^2;
padding = length(rx) - 40 + 1;
c_pad = c(padding:end);

[val_max, max_f] = max(c_pad);
th = val_max * 0.65;    % threshold value for peak detection

% Additional check for a secondary peak ~32 samples later
if (max_f + 32 <= length(c_pad)) && (0.75 * c_pad(max_f) <= c_pad(max_f + 32))
    max_f = max_f + 32;
end
disp("Max Fine Sync index: " + max_f);


%% ========== 7) Frame Extraction ==========

frames = [];
for i = 161 : (length(rx) - frame_len)
    % Condition: c_pad(i) >= threshold, must be local peak, and
    % the region after must remain below threshold to confirm
    if c_pad(i) >= th && max(c_pad(i + 30 : i + 35)) < th && ...
            c_pad(i) >= max(c_pad(i - 3 : i + 3))
        
        peak_val = c_pad(i);
        disp("Peak found at: " + i + " with value = " + peak_val);
        
        % Extract one frame, skipping the repeated preamble portion
        frame = rx(i - FFT_N + 1 : i + frame_len - 2*FFT_N - 2*CP_size);
        frames = [frames; frame];
    end
end


%% ========== 8) Fine Frequency Synchronization (Per-frame) ==========

for i = 1:size(frames,1)
    pre_detected = frames(i,1:FFT_N);
    
    % Similar approach to Schmidl & Cox for half-sym
    sum_val = 0;
    for j = 1:32
        sum_val = sum_val + conj(pre_detected(j)) * pre_detected(32 + j);
    end
    
    df_2 = angle(sum_val) / pi;
    frames(i,:) = frames(i,:) .* exp(-1i * 2 * pi * df_2 * (1:length(frames(i,:))) / 64);
end


%% ========== 9) Channel Estimation (Using Known Preamble) ==========

H   = zeros(size(frames,1), FFT_N);
H_q = zeros(size(frames,1), FFT_N/2);

for i = 1:size(frames,1)
    pre_detected = frames(i, 1:FFT_N);
    pre_y = fft(pre_detected, FFT_N);
    
    % Raw channel estimate in freq domain
    H(i,:) = pre_y ./ preamble;
    
    % Interpolate from half-size to full-size
    H_q(i,:) = H(i,1:2:end);
    H(i,:) = interp1(1:2:FFT_N, H_q(i,:), xq, 'spline','extrap');
    
    % Simple boundary fix: average edges
    H(i,end) = (H(i,end - 1) + H(i,1))/2;
end


%% ========== 10) Demodulation & Channel Equalization ==========

sym_all = [];
ber     = 0;
active_carriers = [];

for i = 1:size(frames,1)
    frame = frames(i, FFT_N + 1 : end);
    rx_data = [];
    sym     = [];
    
    for j = 1:PL
        % Extract one OFDM symbol (CP removed)
        sym_time = frame((j-1)*FFT_N + j*CP_size + 1 : j*FFT_N + j*CP_size);
        sym_freq = fft(sym_time, FFT_N);
        
        % First-level channel equalization
        sym_freq = sym_freq .* conj(H(i,:)) ./ abs(H(i,:)).^2;
        
        % Pilot-based phase compensation
        pilots = [sym_freq(13), sym_freq(22), sym_freq(43), sym_freq(56)];
        pilots = interp1([13, 22, 43, 56], pilots, xq, 'linear', 'extrap');
        pilots = smoothdata(pilots);
        
        % Adjust by pilot estimate
        sym_freq = sym_freq ./ pilots .* abs(pilots);
        
        % Remove pilot subcarriers
        sym_freq([13,22,43,56]) = [];
        
        % QPSK demapping on 60 subcarriers, 4 at a time
        for k = 1 : 4 : 60
            sub_block = sym_freq(k : k+3);
            [sub_data, active_carrier] = demapper(sub_block);
            rx_data         = [rx_data sub_data];
            active_carriers = [active_carriers active_carrier];
        end
        
        sym = [sym sym_freq];
    end
    
    % Convert bits -> text for quick check
    result = bin2text(rx_data);
    disp("Recovered Text (Frame " + i + "): " + result);
    
    % Compute BER relative to known data
    ber_frame = biterr(data, rx_data) / length(data);
    ber       = ber + ber_frame;
    
    % Collect all demodulated freq-domain symbols (for final constellation)
    sym_all = [sym_all sym];
end

ber = ber / size(frames,1);
disp("Average BER over all frames = " + ber);


%% ========== 11) Consolidated Plots with Colors, Legends, Limits ==========

%
% FIGURE 1: Received Signal & Synchronization Metrics
%
figure('Name','Signal & Sync Metrics','NumberTitle','off');

% --- Subplot (2,2,1): abs(rx) & real(rx)
subplot(2,2,1)
plot(abs(rx), 'b', 'DisplayName','|rx|'); 
hold on
plot(real(rx), 'r', 'DisplayName','Re\{rx\}');
hold off
legend('Location','best')
title('Received Signal (Magnitude & Real Part)')
xlabel('Sample Index')
ylabel('Amplitude')

% --- Subplot (2,2,2): M(d) & abs(rx)
subplot(2,2,2)
plot(abs(M), 'b', 'DisplayName','|M(d)|');
hold on
plot(abs(rx), 'r', 'DisplayName','|rx|');
hold off
legend('Location','best')
title('Schmidl & Cox Metric vs. |rx|')
xlabel('Sample Index')
ylabel('Magnitude')

% --- Subplot (2,2,3): c\_pad
subplot(2,2,3)
plot(c_pad, 'g', 'DisplayName','|xcorr|^2')
legend('Location','best')
title('Fine Time Sync Cross-Corr (c\_pad)')
xlabel('Sample Index (shifted)')
ylabel('Correlation')

% --- Subplot (2,2,4): a\_pick
subplot(2,2,4)
plot(a_pick, 'm', 'DisplayName','Rectified M(d)')
legend('Location','best')
title('Rectified Schmidl & Cox Metric (a\_pick)')
xlabel('Sample Index')
ylabel('Convolution Result')


%
% FIGURE 2: Channel Estimation & Constellation
%
figure('Name','Channel & Constellation','NumberTitle','off');

% --- Subplot (2,2,1): Channel Scatter (example: frame 2)
subplot(2,2,1)
if size(frames,1) >= 2
    scatter(real(H(2,:)), imag(H(2,:)), 'o', 'MarkerEdgeColor','b', ...
        'DisplayName','Frame 2 Channel');
    axis equal
    title('Estimated Channel Scatter (Frame 2)')
    xlabel('Real\{H\}')
    ylabel('Imag\{H\}')
    grid on
    legend('Location','best')
else
    text(0.5,0.5,'Not enough frames to show Frame 2','HorizontalAlignment','center')
end

% --- Subplot (2,2,2): Channel Magnitude (frame 2)
subplot(2,2,2)
if size(frames,1) >= 2
    plot(abs(H(2,:)), 'b','DisplayName','|H|')
    title('Magnitude of Channel Estimate (Frame 2)')
    xlabel('Subcarrier Index')
    ylabel('|H|')
    grid on
    legend('Location','best')
else
    text(0.5,0.5,'Not enough frames to show Frame 2','HorizontalAlignment','center')
end

% --- Subplot (2,2,3): Channel Phase (frame 2)
subplot(2,2,3)
if size(frames,1) >= 2
    plot(angle(H(2,:)), 'r','DisplayName','angle(H)')
    title('Phase of Channel Estimate (Frame 2)')
    xlabel('Subcarrier Index')
    ylabel('Phase(H)')
    grid on
    legend('Location','best')
else
    text(0.5,0.5,'Not enough frames to show Frame 2','HorizontalAlignment','center')
end

% --- Subplot (2,2,4): Constellation of All Recovered Symbols & Carriers
subplot(2,2,4)
hold on
% Plot recovered symbols in one color/marker
scatter(real(sym_all), imag(sym_all), 'o', ...
    'MarkerEdgeColor','blue','DisplayName','Recovered Symbols')
% Plot the "active carriers" in another color/marker
scatter(real(active_carriers), imag(active_carriers), 'x', ...
    'MarkerEdgeColor','red','DisplayName','Active Carriers')
hold off
legend('Location','best')
title('Combined Constellation')
xlabel('In-Phase')
ylabel('Quadrature')
grid on
xlim([-4 4])
ylim([-4 4])


%% ========== Helper Functions ==========

%%  A) QPSK Demapper
function [d, p] = demapper(sub_block)
    d = zeros(1,6);
    [p, x] = maxk(sub_block, 2);
    x = sort(x);
    i1 = x(1);
    i2 = x(2);

    % Bits from real/imag sign of top 2 picks
    d(6) = real(sub_block(i2)) < 0;
    d(5) = imag(sub_block(i2)) < 0;
    d(4) = real(sub_block(i1)) < 0;
    d(3) = imag(sub_block(i1)) < 0;

    % Decide first 2 bits by index pair
    if isequal([i1 i2], [1 2])
        d(1:2) = [0 0];
    elseif isequal([i1 i2], [2 3])
        d(1:2) = [0 1];
    elseif isequal([i1 i2], [3 4])
        d(1:2) = [1 0];
    elseif isequal([i1 i2], [1 4])
        d(1:2) = [1 1];
    end
end

%%  B) Convert Binary Vector to Text
function text = bin2text(binVS)
    btxt = reshape(binVS, [8, length(binVS)/8])';
    text = char(bin2dec(char(btxt + 48)))';
end

%%  C) Convert Text to Binary
function [binV, binS] = text2bin(text)
    binS = dec2bin(text,8);  
    binS = binS';           
    binS = binS(:)';        
    binV = binS - 48;       
end
