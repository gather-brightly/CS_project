function bw_optimal = bw_finder(data_rate, bandwidth_list, k, W)
    
    % initialize the bandwidth-snr pair matrix
    bw_snr = zeros(length(bandwidth_list), 2);
    
    % loop throught the candidates of the bandwidth to find the optimal
    % bandwidth
    for bandwidth_idx = 1:length(bandwidth_list)
        bandwidth = bandwidth_list(bandwidth_idx);
        
        % generate X_temp from the given data_rate
        X_temp = X_generator(data_rate, bandwidth, 1, k);

        % generate snr_temp using the trained model parameter 'W'
        snr_temp = transpose(W) * X_temp;
        
        % store the predicted snr value into the 'bw_snr' matrix
        bw_snr(bandwidth_idx, 1) = bandwidth;
        bw_snr(bandwidth_idx, 2) = snr_temp;

    end

    % find the optimal bandwidth for the given data_rate
    % the optimal bandwidth maximizes the predicted snr
    snr_max = max(bw_snr(:, 2));
    snr_max_idx = find(bw_snr(:, 2) == snr_max);
    bw_optimal = bw_snr(snr_max_idx, 1);
end