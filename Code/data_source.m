%% data_source.m - Generates the source data
function [dataBits, numBits] = data_source(textMessage)
    fprintf('Transmitted Text: %s \n', textMessage);
    dataBits = reshape(de2bi(uint8(textMessage), 8, 'left-msb')', [], 1);
    numBits = length(dataBits);  
end