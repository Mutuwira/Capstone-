% Define the receiver object
rx = comm.SDRuReceiver( ...
    'Platform', 'B200', ...
    'SerialNum', '329AEA8', ... 
    'CenterFrequency', 915e6, ... % 915 MHz 
    'Gain', 30, ... % Adjust gain for better sensitivity
    'SampleRate', 1e6, ... % 1 MS/s
    'OutputDataType', 'double', ...
    'SamplesPerFrame', 10000); % Number of samples to receive per frame

% Open the USRP stream
rxData = [];
for i = 1:100
    data = rx(); % Receive samples
    if ~isempty(data)
        rxData = [rxData; data]; % Store data
    end
end

% Release USRP resource
release(rx);

% Plot received signal (raw noise)
figure;
plot(real(rxData));
title('Received Noise Signal');
xlabel('Sample Index');
ylabel('Amplitude');

% Plot Power Spectral Density (PSD) using FFT
Fs = 1e6; % Sample rate (modify based on your setting)
N = length(rxData);
frequencies = (-N/2:N/2-1)*(Fs/N);

figure;
plot(frequencies/1e6, 10*log10(abs(fftshift(fft(rxData))).^2))
title('Power Spectral Density of Received Noise');
xlabel('Frequency (MHz)');
ylabel('Power (dB)');
grid on;
