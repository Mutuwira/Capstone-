function encoded_bits = systematic_conv_encoder(data_bits)
% systematic_conv_encoder  Rate-1/2 systematic convolutional encoder.
%   encoded_bits = systematic_conv_encoder(data_bits)
%
%   Input:
%       data_bits - Column vector of binary bits (0s and 1s)
%
%   Output:
%       encoded_bits - Column vector of encoded bits, where for each
%                      input bit two output bits are generated:
%                        [data_bit; parity_bit]

    N = length(data_bits);
    encoded_bits = zeros(2*N, 1);
    
    % Initialize shift register (assume zero initial state)
    % The register holds the previous two bits.
    reg = [0 0];
    
    for n = 1:N
        current_bit = data_bits(n);
        % Systematic output is the input bit.
        encoded_bits(2*n-1) = current_bit;
        % Parity bit is computed as mod(current_bit + reg(1) + reg(2), 2)
        parity = mod(current_bit + reg(1) + reg(2), 2);
        encoded_bits(2*n) = parity;
        % Update shift register: new bit enters and the older bits shift
        reg = [current_bit, reg(1)];
    end
end
