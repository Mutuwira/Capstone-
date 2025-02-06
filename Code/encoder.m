% Linear Encoder
function encodedBits = encoder(dataBits, code_rate, k)
    k = 4;
    % Encoded bits per block (including redundancy)
    n = round(k/code_rate);
    if n <= k
        error('n must be greater than k for valid encoding')
    end
    
    dataBits = dataBits(:);

    % Construct the Generator Matrix G dynamically
    I_k = eye(k); % Systematic identity matrix
    P = randi([0 1], k, n-k); % Random parity check bits (adjustable)
    G = [I_k, P]; % Final generator matrix

    % Ensure input is a multiple of k
    [rows, cols] = size(dataBits);
    if rows ~= k 
        error('Input data_bits should have %d rows for correct encoding.', k)
    end

    % Encode using matrix multiplication (mod 2)
    encoded_blocks = mod(G * double(dataBits), 2);

    % Convert to a column vector
    encodedBits = encoded_blocks(:);
end