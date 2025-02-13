function decodedBits = decoder(receivedEncodedBits, n, k, codeType)
    % Ensure receivedEncodedBits is a column vector for consistency
    if isrow(receivedEncodedBits)
        receivedEncodedBits = receivedEncodedBits';
    end

    % Reshape receivedEncodedBits into a matrix with n columns
    receivedEncodedBits = reshape(receivedEncodedBits, n, []).';

    switch lower(codeType)
        case 'hamming'
            % Hamming Decoding
            m = n - k;                 % Number of parity bits
            [h, genMatrix] = hammgen(m); % Generate parity-check and generator matrix
            syndromeTable = syndtable(h); % Generate syndrome decoding table

        case 'linear'
            % Linear Block Decoding
            genMatrix = [eye(k), randi([0, 1], k, n - k)]; % Example generator matrix
            h = gen2par(genMatrix);         % Generate parity-check matrix
            syndromeTable = syndtable(h);   % Generate syndrome decoding table

        case 'cyclic'
            % Cyclic Decoding
            genPoly = cyclpoly(n, k);       % Generate generator polynomial
            [h, genMatrix] = cyclgen(n, genPoly); % Generate parity-check and generator matrix
            syndromeTable = syndtable(h);   % Generate syndrome decoding table

        otherwise
            error('Unsupported codeType. Use "hamming", "linear", or "cyclic".');
    end

    % Calculate the syndrome
    syndrome = mod(receivedEncodedBits * h', 2);

    % Find error locations using the syndrome table
    errorLocations = syndromeTable(bi2de(syndrome, 'left-msb') + 1, :);

    % Correct errors in the received code
    correctedCode = mod(receivedEncodedBits + errorLocations, 2);

    % Extract the message bits
    if isequal(genMatrix(:, 1:k), eye(k))
        % Systematic generator matrix (message in first k columns)
        decodedBits = correctedCode(:, 1:k);
    elseif isequal(genMatrix(:, end-k+1:end), eye(k))
        % Systematic generator matrix (message in last k columns)
        decodedBits = correctedCode(:, end-k+1:end);
    else
        error('Unsupported generator matrix form. Ensure it is systematic.');
    end

    % Reshape decodedBits back to a single column vector
    decodedBits = decodedBits.';
    decodedBits = decodedBits(:);
end
