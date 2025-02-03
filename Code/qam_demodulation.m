function receivedBits = qam_demodulation(equalizedSymbols, bitsPerSymbol, numDataSymbols, M)
    receivedBits = zeros(numDataSymbols * bitsPerSymbol, 1);
    for i = 1:numDataSymbols
        bits = qam_to_gray(equalizedSymbols(i), M);
        receivedBits((i-1)*bitsPerSymbol + (1:bitsPerSymbol)) = bits';
    end

    function bits = qam_to_gray(qamSymbol, M)
    if M == 16
        % Gray code mapping for 16-QAM
        mapping = [
            -3 + 3i, -3 + 1i, -3 - 3i, -3 - 1i, ...
            -1 + 3i, -1 + 1i, -1 - 3i, -1 - 1i, ...
             3 + 3i,  3 + 1i,  3 - 3i,  3 - 1i, ...
             1 + 3i,  1 + 1i,  1 - 3i,  1 - 1i
        ] / sqrt(10);
        [~, idx] = min(abs(qamSymbol - mapping));
        grayIndex = idx - 1;
        bits = de2bi(grayIndex, log2(M), 'left-msb')';
    else
        error('Demapping not implemented for this value of M');
    end
end
end
