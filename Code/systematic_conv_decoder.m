function decoded_bits = systematic_conv_decoder(encoded_bits)
% systematic_conv_decoder  Viterbi decoder for the rate-1/2 systematic convolutional code.
%   decoded_bits = systematic_conv_decoder(encoded_bits)
%
%   Input:
%       encoded_bits - Column vector of hard decisions (0s and 1s),
%                      length should be an even number.
%
%   Output:
%       decoded_bits - Column vector of decoded bits.

    % Check that the input length is even.
    if mod(length(encoded_bits),2) ~= 0
        error('Length of encoded_bits must be even.');
    end

    N = length(encoded_bits) / 2;   % Number of input bits
    numStates = 4;                  % For constraint length 3
    
    % Precompute state transitions and outputs.
    % Represent state as an integer (0,1,2,3) corresponding to the two-bit register:
    % e.g., state 0 -> [0 0], state 1 -> [1 0], state 2 -> [0 1], state 3 -> [1 1].
    nextState = zeros(numStates,2);    % nextState(s+1, input+1)
    outParity = zeros(numStates,2);      % outParity(s+1, input+1)
    
    for s = 0:numStates-1
        % Extract register bits. (Using bitget: bit 1 is least-significant.)
        s1 = bitget(s,1);  % most recent (position 1)
        s2 = bitget(s,2);  % older (position 2)
        for input = 0:1
            % New state: shift the register; new state = [input, s1]
            newState = input + 2*s1;
            nextState(s+1, input+1) = newState;
            % Expected parity: mod(input + s1 + s2, 2)
            outParity(s+1, input+1) = mod(input + s1 + s2, 2);
        end
    end
    
    % Initialize path metrics (set to a large number for unreachable states)
    pathMetric = inf(numStates, N+1);
    pathMetric(1,1) = 0;  % Assume starting in state 0
    % To store decisions and predecessors:
    predecessor = zeros(numStates, N);
    inputDecision = zeros(numStates, N);
    
    % Viterbi algorithm (forward recursion)
    for t = 1:N
        r_sys = encoded_bits(2*t-1);  % Received systematic bit
        r_par = encoded_bits(2*t);      % Received parity bit
        for s = 0:numStates-1
            if isfinite(pathMetric(s+1, t))
                for input = 0:1
                    newState = nextState(s+1, input+1);
                    exp_sys = input;  % Expected systematic output equals the input.
                    exp_par = outParity(s+1, input+1);
                    % Hamming distance metric
                    branchMetric = (exp_sys ~= r_sys) + (exp_par ~= r_par);
                    metric = pathMetric(s+1, t) + branchMetric;
                    if metric < pathMetric(newState+1, t+1)
                        pathMetric(newState+1, t+1) = metric;
                        predecessor(newState+1, t) = s;
                        inputDecision(newState+1, t) = input;
                    end
                end
            end
        end
    end
    
    % Traceback: Choose the state with the minimum metric at time N+1.
    [~, final_state] = min(pathMetric(:, N+1));
    final_state = final_state - 1;
    decoded_bits = zeros(N, 1);
    state = final_state;
    for t = N:-1:1
        input = inputDecision(state+1, t);
        decoded_bits(t) = input;
        state = predecessor(state+1, t);
    end
end
