% Function to generate ZC sequence with CP and data

function tx_signal = generate_signal(seq_length, root, cp_length) 
    % Generate Zadoff-Chu (ZC) sequence
    n = 0:seq_length-1;
    zc_sequence = exp(-1i*pi*root*n.*n/seq_length);  % ZC sequence formula

    % Add Cyclic Prefix (CP) to ZC sequence
    zc_with_cp = [zc_sequence(end-cp_length+1:end), zc_sequence];
    tx_signal = zc_with_cp; % ZC with cp and random length zero data
end