function decodedBits = decoder(receivedEncodedBits, n, k )

    % Ensure receivedEncodedBits is a column vector
    if isrow(receivedEncodedBits)
        receivedEncodedBits = receivedEncodedBits';
    end

    % Ensure the input is a multiple of n
    numBlocks = ceil(length(receivedEncodedBits) / n);
    paddedLength = numBlocks * n;

    if length(receivedEncodedBits) < paddedLength
        padding = zeros(paddedLength - length(receivedEncodedBits), 1); % Zero-padding
        receivedEncodedBits = [receivedEncodedBits; padding];
    end

    % Reshape receivedEncodedBits into a matrix with n columns
    receivedEncodedBits = reshape(receivedEncodedBits, n, []).';

    genMatrix = [eye(k), randi([0, 1], k, n - k)]; % Generator matrix

    % Derive the parity-check matrix from the generator matrix
    paritySubMatrix = genMatrix(:, k+1:end); % Extract parity submatrix
    hMatrix = [paritySubMatrix'; eye(n - k)]; % Parity-check matrix [P' | I_(n-k)]
    hMatrix = double(hMatrix);
    hMatrix = hMatrix.';
  

    % Calculate the syndrome: S = receivedEncodedBits * H' mod 2
    syndrome = mod(receivedEncodedBits * hMatrix, 2);
    % Create a decoding table (syndrome table) for error correction
    syndromeTable = syndtable(hMatrix);

    % Find error locations using the syndrome table
    errorLocations = syndromeTable(bi2de(syndrome, 'left-msb') + 1, :);

    % Correct errors in the received code
    correctedCode = mod(receivedEncodedBits + errorLocations, 2);

    % Extract the message bits
     decodedBits = correctedCode(:, 1:k);

    % Reshape decodedBits back to a single column vector
    decodedBits = decodedBits.';
    decodedBits = decodedBits(:);
end
