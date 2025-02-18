% function encodedBits = encoder(dataBits, n, k)
%     % Check if parameters are valid
%     if n <= k
%         error('n must be greater than k for valid encoding');
%     end
% 
%     dataBits = dataBits(:); % Ensure column vector
% 
%     % Select appropriate generator matrix based on codeType
%     switch codeType
%         case 'Hamming'
%             % Generate Hamming (n, k) Code Generator Matrix
%             G = generateHammingG(n, k);
% 
%         case 'BCH'
%             % Generate BCH (n, k) Code Generator Matrix
%             G = generateBCHG(n, k);
% 
%         case 'RandomLinear'
%             % Default: Systematic Linear Block Code
%             I_k = eye(k);
%             P = randi([0 1], k, n-k); % Random parity check bits
%             G = [I_k, P];
% 
%         otherwise
%             error('Unsupported encoding type. Choose Hamming, BCH, or RandomLinear.');
%     end
% 
%     % Convert G to double for matrix multiplication
%     G = double(G);
% 
%     % Ensure input is a multiple of k
%     numBlocks = ceil(length(dataBits) / k);
%     paddedLength = numBlocks * k;
% 
%     if length(dataBits) < paddedLength
%         padding = zeros(paddedLength - length(dataBits), 1); % Zero-padding
%         dataBits = [dataBits; padding];
%     end
% 
%     % Reshape into k-bit blocks
%     dataMatrix = reshape(dataBits, k, [])';
% 
%     % Convert dataMatrix to double for multiplication
%     dataMatrix = double(dataMatrix);
% 
%     % Encode using matrix multiplication (mod 2) OR optimized XOR
%     encoded_blocks = mod(dataMatrix * G, 2);
% 
%     % Convert to a column vector
%     encodedBits = encoded_blocks';
%     encodedBits = encodedBits(:);
% end
% 
% function G = generateHammingG(n, k)
%     % Hamming Code (7,4) Example
%     if n == 7 && k == 4
%         G = [1 0 0 0 1 1 0; 
%              0 1 0 0 1 0 1; 
%              0 0 1 0 0 1 1; 
%              0 0 0 1 1 1 1];
%     else
%         error('Only Hamming(7,4) is supported for now.');
%     end
% end
% 
% function G = generateBCHG(n, k)
%     % Placeholder for BCH code (use poly2trellis or MATLAB's bchgenpoly)
%     % BCH codes require a generator polynomial, this is an example.
%     [~, G] = bchgenpoly(n, k);
%     G = double(G);
% end

function encodedBits = encoder(dataBits, n, k, genMatrix)
    % Ensure dataBits is a column vector for consistency
    if isrow(dataBits)
        dataBits = dataBits';
    end

    % Preprocess the input into a K-column binary matrix
    % If the total number of bits in dataBits is not divisible by k, pad with zeros
    addedBits = mod(size(dataBits, 1), k);
    if addedBits ~= 0
        padding = k - addedBits;
        dataBits = [dataBits; zeros(padding, 1)];
    end

    % Reshape the input into a matrix with k columns
    dataBits = reshape(dataBits, k, []).';

    % Generating a simple generator matrix
    encodedBits = mod(double(dataBits) * double(genMatrix), 2);

    % Reshape encodedBits back to a single column vector
    encodedBits = encodedBits.';
    encodedBits = encodedBits(:);
end
