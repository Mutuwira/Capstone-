%% source_decoder.m - Converts binary sequence back to text
function decodedText = source_decoder(decodedBits, numBits)
    decodedBits = decodedBits(1:numBits); % Extract original bits (handle padding if needed)
    decodedText = char(bi2de(reshape(decodedBits, 8, []).', 'left-msb'))';
end