function qamSymbols = qam_modulation(dataBits, bitsPerSymbol, N, M)
    qamSymbols = zeros(N, 1);
    for i = 1:N
        bits = dataBits((i-1)*bitsPerSymbol + (1:bitsPerSymbol))';
        grayIndex = bi2de(bits, 'left-msb');
        qamSymbols(i) = gray_to_qam(grayIndex, M);
    end

    function qamSymbol = gray_to_qam(grayIndex, M)
    if M == 16
        % Gray code mapping for 16-QAM
        mapping = [
            -3 + 3i, -3 + 1i, -3 - 3i, -3 - 1i, ...
            -1 + 3i, -1 + 1i, -1 - 3i, -1 - 1i, ...
             3 + 3i,  3 + 1i,  3 - 3i,  3 - 1i, ...
             1 + 3i,  1 + 1i,  1 - 3i,  1 - 1i
        ] / sqrt(10);
        qamSymbol = mapping(grayIndex + 1);
    else
        error('Mapping not implemented for this value of M');
    end
end
end
