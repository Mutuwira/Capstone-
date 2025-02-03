% MMSE Equalization

function equalizedSymbols = equalization(Data_rx_fft, H_est, signalPower, noiseVariance)
    mmseFilter = conj(H_est) ./ (abs(H_est).^2 + noiseVariance / signalPower);
    equalizedSymbols = Data_rx_fft .* mmseFilter;
end
