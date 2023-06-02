function X = X_generator(rate_list, bandwidth_list, n_iter, k)

    % number of experiments
    n_experiment = length(rate_list) * length(bandwidth_list); 
    % initialize X (design matrix)
    X = zeros(2*k, n_iter*n_experiment);
    
    cnt = 1;

    % generate X through the grid of data rates and bandwidths
    for rate_idx = 1:length(rate_list)
        rate = rate_list(rate_idx);
        for bandwidth_idx = 1:length(bandwidth_list)
            bandwidth = bandwidth_list(bandwidth_idx);
            
            % in a single experiment condition, iterate 'n_iter' times
            for iter = 1:n_iter
                X_col_temp = zeros(2*k,1);
                for k_idx = 1:k
                    X_col_temp(k_idx) = rate^k_idx;
                    X_col_temp(k_idx + k) = bandwidth^k_idx;
                end
                X(:, cnt) = X_col_temp; % store 'X_col_temp' as a column in X
                cnt = cnt + 1;
            end
    
        end
    end

    % add ones at the first row of X to get bias weight
    X = [ones(1, n_iter*n_experiment); X];
end