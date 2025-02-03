function display_results(tx_signal, rx_signal, corr_output, qamSymbols, equalizedSymbols, SNR_dB)
    figure;
    subplot(3,1,1);
    plot(real(tx_signal));
    title('Transmitted Signal');
    xlabel('Sample Index');
    ylabel('Amplitude');

    subplot(3,1,2);
    plot(real(rx_signal));
    title('Received Signal');
    xlabel('Sample Index');
    ylabel('Amplitude');

    subplot(3,1,3);
    plot(corr_output);
    title('Cross-Correlation Output for Synchronization');
    xlabel('Lag');
    ylabel('Correlation Magnitude');

    scatterplot(qamSymbols);
    title('Transmitted QAM Constellation');

    scatterplot(equalizedSymbols(:));
    title(['Received QAM Constellation (SNR = ', num2str(SNR_dB), ' dB)']);
end
