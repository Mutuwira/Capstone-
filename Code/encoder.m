function encodedBits = encoder(dataBits, n, k, genMatrix)
    % Ensure dataBits is a column vector for consistency
    if isrow(dataBits)
        dataBits = dataBits';
    end

    % Preprocessing the input into a K-column binary matrix
    % If the total number of bits in dataBits is not divisible by k, pad with zeros
    addedBits = mod(size(dataBits, 1), k);
    if addedBits ~= 0
        padding = k - addedBits;
        dataBits = [dataBits; zeros(padding, 1)];
    end

    dataBits = reshape(dataBits, k, []).'; % Reshape into matrix with k columns
    encodedBits = mod(double(dataBits) * double(genMatrix), 2);

    % Reshape encodedBits back to a single column vector
    encodedBits = encodedBits.';
    encodedBits = encodedBits(:);
end
