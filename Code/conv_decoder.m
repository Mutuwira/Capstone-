function decoded_bits = conv_decoder(encoded_bits, K, generator)
% DYNAMIC_SYSTEMATIC_CONV_DECODER  Viterbi decoder for a rate-1/2 systematic convolutional code.
%
%   decoded_bits = dynamic_systematic_conv_decoder(encoded_bits, K, generator)
%
%   Inputs:
%     encoded_bits - Column vector of hard-decision bits (0s and 1s).
%                    Its length should be even (2 bits per encoded symbol).
%     K            - Constraint length used in the encoder.
%     generator    - Generator vector of length K.
%
%   Output:
%     decoded_bits - Column vector of decoded bits.
%
%   This implementation uses the full trellis (with 2^(K-1) states) and
%   a Hamming distance metric for branch metrics.

    % Check that the length of encoded_bits is even.
    if mod(length(encoded_bits), 2) ~= 0
        error('Length of encoded_bits must be even.');
    end

    N = length(encoded_bits) / 2;  % Number of input bits
    numStates = 2^(K - 1);         % Number of states

    % Precompute trellis: for each state and each possible input (0 or 1)
    nextState = zeros(numStates, 2);    % nextState(s+1, input+1)
    expectedParity = zeros(numStates, 2); % expected parity bit for each transition

    % Loop over all states
    for s = 0:numStates - 1
        % Convert state to binary vector (of length K-1, MSB first)
        state_bits = de2bi(s, K - 1, 'left-msb');
        for input = 0:1
            % New state: shift the register: new state = [input, state_bits(1:end-1)]
            if K - 1 > 0
                newState_bits = [input, state_bits(1:end-1)];
            else
                newState_bits = [];
            end
            newState = bi2de(newState_bits, 'left-msb');
            nextState(s + 1, input + 1) = newState;
            % Expected parity: using generator polynomial
            % The systematic bit is the input.
            % Parity = mod(generator(1)*input + sum(generator(2:end) .* state_bits), 2)
            expectedParity(s + 1, input + 1) = mod(generator(1) * input + sum(generator(2:end) .* state_bits), 2);
        end
    end

    % Initialize path metrics
    pathMetric = inf(numStates, N + 1);
    pathMetric(1, 1) = 0;  % Assume initial state is 0

    % To store predecessor and decision for traceback
    predecessor = zeros(numStates, N);
    decision = zeros(numStates, N);

    % Viterbi algorithm (forward recursion)
    for t = 1:N
        r_sys = encoded_bits(2 * t - 1);  % Received systematic bit
        r_par = encoded_bits(2 * t);        % Received parity bit
        for s = 0:numStates - 1
            if isfinite(pathMetric(s + 1, t))
                for input = 0:1
                    ns = nextState(s + 1, input + 1);
                    exp_sys = input;  % Systematic bit expected is the input bit
                    exp_par = expectedParity(s + 1, input + 1);
                    % Use Hamming distance as branch metric
                    branchMetric = (exp_sys ~= r_sys) + (exp_par ~= r_par);
                    metric = pathMetric(s + 1, t) + branchMetric;
                    if metric < pathMetric(ns + 1, t + 1)
                        pathMetric(ns + 1, t + 1) = metric;
                        predecessor(ns + 1, t) = s;
                        decision(ns + 1, t) = input;
                    end
                end
            end
        end
    end

    % Traceback: choose the state with minimum metric at time N+1
    [~, final_state] = min(pathMetric(:, N + 1));
    final_state = final_state - 1;
    decoded_bits = zeros(N, 1);
    state = final_state;
    for t = N:-1:1
        input = decision(state + 1, t);
        decoded_bits(t) = input;
        state = predecessor(state + 1, t);
    end
end
