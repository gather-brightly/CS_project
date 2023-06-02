%% compare the models with different 'k' parameter
%% basic parameters
addpath('./subfunctions/')

rate_list = [0.1 : 0.1 : 1]; % list of data rate (in Hz)
bandwidth_list = [0.1 : 0.1 : 1*2]; % list of bandwidth (in Hz)
k_list = [1:1:5]; % list of k (hyperparameter)

n_iter_test = 2; % (for test) number of iterations for each experiment

% parameters about FM modulation
A_c = 1; % amplitude of the carrier
f_c = 10; % carrier frequency
f_s = 20; % sampling frequency
k_f = 7.5; % frequency deviation constant


%% implementation

for k_idx = 1:length(k_list)
    k = k_list(k_idx);
    
    disp(['k = ', num2str(k), ';']);

%% generate test X
X_test = X_generator(rate_list, bandwidth_list, n_iter_test, k);
X_test = transpose(X_test);

%% generate test Y
N = length(rate_list) * length(bandwidth_list); % number of experiments

% initialize Y_train
Y_test_true = zeros(1, N*n_iter_test);

cnt = 1;
% generate Y through the grid of data rates and bandwidths
for rate_idx = 1:length(rate_list)
    rate = rate_list(rate_idx);
    for bandwidth_idx = 1:length(bandwidth_list)
        bandwidth = bandwidth_list(bandwidth_idx);
        
        % simulink parameter
        f_m = rate; % frequency of the single tone message
        % FM modulation via simulink
        simOut = sim('modulation_531.slx'); % Replace 'your_simulink_model' with the actual name of your Simulink model
        % extract the modulation output signal from the Simulink output
        x_c = simOut.yout{1}.Values.Data; % FM modulation output signal
        r_t_timevector = simOut.tout; % time vector of the FM modulation output signal 'x_c'

        for iter_test = 1:n_iter_test

            % AWGN (Additive White Gaussian Noise)
            r = awgn(x_c, 10); % signal after AWGN channel = received signal
            % standardize the received signal 'r'
            r_standardized = standardize(r);
           
            % predetection filter design
            w_c2 = 2*pi * (f_c + (bandwidth/2)); % upper cutoff frequency
            w_c1 = 2*pi * (f_c - (bandwidth/2)); % lower cutoff frequency
            prefilter_num = [w_c2 0]; % numerator of the transfer function of the predetection filter
            prefilter_den = [1 w_c2 w_c1*w_c2]; % denominator of the transfer function of the predetection filter
            prefilter = tf(prefilter_num, prefilter_den); % transfer function of the predetection filter
            % apply predetection filter to the received signal
            e_1 = lsim(prefilter, r_standardized, r_t_timevector);
            e_1_standardized = standardize(e_1);
    
            % demodulation
            y_D = fmdemod(e_1_standardized, f_c, f_s, k_f);
            % delete outliers of y_D (outlier threshold: 2 times the IQR)
            y_D = outlier_delete(y_D);
            % standardize y_D
            y_D = standardize(y_D);
    
            % SNR calculation
            % signal vector
            m = sin(2*pi*f_m*r_t_timevector);
            % calculate the noise vector
            noise = m - y_D;
            % calculate the power of the signal and noise vectors
            signalPower = sum(m.^2);
            noisePower = sum(noise.^2);
            % calculate the SNR
            snr = signalPower / noisePower;
            
            % put the obtained SNR into Y_test
            Y_test_true(1, cnt) = snr;

            cnt = cnt+1;
        end

    end
end

%% MSE that shows the model accuracy
% Compute predictions for test data
W = load(['./data/W_k-', num2str(k), '.mat']).W;
Y_test_pred = X_test * W;

% Calculate RMSE
Y_test_true = transpose(Y_test_true);
mse = mean((Y_test_pred - Y_test_true).^2);
rmse = sqrt(mse);

% Display RMSE
fprintf('Under k=%d / Root Mean Squared Error (RMSE) on test data: %.4f\n', k, rmse);

end