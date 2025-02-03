% Pad bits to match OFDM structure
function [dataBits, numDataSymbols] = bit_padding(dataBits, bitsPerSymbol, N)
    numBits = length(dataBits);
    numDataSymbols = ceil(numBits / bitsPerSymbol);
    numDataSymbols = ceil(numDataSymbols / N) * N;
    numBitsPadded = numDataSymbols * bitsPerSymbol;
    dataBits = [dataBits; zeros(numBitsPadded - numBits, 1)];
end
