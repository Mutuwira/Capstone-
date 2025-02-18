function decodedBits = decoder(receivedEncodedBits, n, k, genMatrix)

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
        
    % Derive the parity-check matrix from the generator matrix
    paritySubMatrix = genMatrix(:, k+1:end); % Extract parity submatrix
    hMatrix = [paritySubMatrix.' eye(n - k)]; % Parity-check matrix [P' | I_(n-k)]

    % Calculate the syndrome: S = receivedEncodedBits * H' mod 2
    syndrome = mod(receivedEncodedBits * hMatrix', 2);
    
    % Create a decoding table (syndrome table) for error correction
    syndromeTable = syndtable(hMatrix);

    % Find error locations using the syndrome table
    syndromeIndex = bi2de(syndrome, 'left-msb') + 1;

    % Look up error patterns
    errorLocations = syndromeTable(syndromeIndex, :);

    % Correct errors in the received code
    %correctedCode = mod(receivedEncodedBits + errorLocations, 2);
    correctedCode = mod(receivedEncodedBits + reshape(errorLocations, size(receivedEncodedBits)), 2);

    % Extract the message bits
     decodedBits = correctedCode(:, 1:k);

    % Reshape decodedBits back to a single column vector
    decodedBits = decodedBits.';
    decodedBits = decodedBits(:);
end
