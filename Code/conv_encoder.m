function encoded_bits = conv_encoder(data_bits, K, generator)
% DYNAMIC_SYSTEMATIC_CONV_ENCODER  Rate-1/2 systematic convolutional encoder.
%
%   encoded_bits = dynamic_systematic_conv_encoder(data_bits, K, generator)
%
%   Inputs:
%     data_bits - Column vector of binary data (0s and 1s).
%     K         - Constraint length (must be >= 2).
%     generator - A vector of length K defining the generator polynomial taps.
%                 The first coefficient applies to the current bit.
%
%   Output:
%     encoded_bits - Column vector of encoded bits (2 bits per input bit).
%
%   Example:
%     For a simple code with K = 3, use:
%         generator = [1 1 1];  
%     which computes parity = mod(b[n] + b[n-1] + b[n-2], 2).

    % Convert inputs to double to avoid integer type conflicts.
    data_bits = double(data_bits);
    generator = double(generator);
    
    % Check input consistency
    if length(generator) ~= K
        error('Generator polynomial length must equal K.');
    end

    N = length(data_bits);
    encoded_bits = zeros(2 * N, 1);  % For each input bit, output 2 bits

    % Initialize the shift register (K-1 memory elements)
    reg = zeros(1, K - 1);

    for n = 1:N
        current_bit = data_bits(n);
        % Systematic output is just the input bit.
        encoded_bits(2 * n - 1) = current_bit;
        % Compute parity:
        %   parity = mod( generator(1)*current_bit + sum(generator(2:end).*reg), 2 )
        parity = mod(generator(1) * current_bit + sum(generator(2:end) .* reg), 2);
        encoded_bits(2 * n) = parity;
        % Update shift register: shift in current_bit
        reg = [current_bit, reg(1:end-1)];
    end
end
