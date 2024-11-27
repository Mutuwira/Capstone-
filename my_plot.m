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