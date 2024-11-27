%% Channel Estimation

    % Inputs:
        % Rx  -  Received Signal
        % xp  -  data signal (xp = xp' - CP) CP = cyclic prefix
    % Output:
        % H_hat - 

function H_hat = ChannelEstimation(rx,xp)
    % r = x * h + w
        % r, received signal
        % x, data
        % h, channel
        % w, additive white gaussian noise (awgn)
    

    % first, let us generate the awgn
    

    % we can approximate xp' * h to xp * h
    % which we can express in the frequency domain
    % Xp . H

    % transfer to frequency domain (fast fourier transform)
    Xp = fft(xp);
    Rx = fft(rx);
    
    H_hat = Rx / Xp;
    % or if we consider awgn
    % H_hat = (Rx-W) / Xp;

end


%% Mean Squared Error (MSE)

%Placeholder
actual = [3.0, -0.5, 2.0, 7.0];
predicted = [2.5, 0.0, 2.1, 7.8];
mse = calculateMSE(actual, predicted);
fprintf('Mean Squared Error (MSE): %.4f\n', mse);
%%%-------%%%

% Inputs:
    %   H     - Vector of true values of the channel
    %   H_hat - Vector of estimated values of the channel
    %
    % Output:
    %   MSE   - Mean Squared Error between the input vectors

function MSE = calculateMSE(H_hat, H)
    if length(H_hat) ~= length(H)
        error('Vectors "H_hat" and "H" must be of the same length.');
    end
    MSE = mean((H - H_hat).^2);
end

%% MSE vs. SNR

% This function plots the MSE against SNR.
function MSEvsSNR(SNR, MSE)
    if length(snrValues) ~= length(mseValues)
        error('Vectors "snrValues" and "mseValues" must be of the same length.');
    end
    
    % Plot the MSE vs SNR
    figure;
    plot(SNR, MSE, '--', 'LineWidth', 1.5);
    xlabel('SNR (dB)');
    ylabel('Mean Squared Error (MSE)');
    title('MSE vs SNR');
    grid on;
end