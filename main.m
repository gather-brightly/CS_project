%% 0. basic parameters
addpath('./subfunctions/')

n_iter_train = 10; % (for train) number of iterations for each experiment
n_iter_test = 2; % (for test) number of iterations for each experiment
k = 1; % highest order of power applied to each GLM input

rate_list = [0.1 : 0.1 : 1]; % list of data rate (in Hz)
bandwidth_list = [0.1 : 0.1 : 1*2]; % list of bandwidth (in Hz)

% parameters about FM modulation
A_c = 1; % amplitude of the carrier
f_c = 10; % carrier frequency
f_s = 20; % sampling frequency
k_f = 7.5; % frequency deviation constant


%% 1. implementation - generate input of the GLM (X)

X_train = X_generator(rate_list, bandwidth_list, n_iter_train, k);

X_filename = ['X_train_k-', num2str(k), '.mat'];
save(['./data/', X_filename], 'X_train');

%% 2. implementation - generate ground truth of the GLM (Y), via simulink

N = length(rate_list) * length(bandwidth_list); % number of experiments
% initialize Y_train
Y_train = zeros(1, N*n_iter_train);

cnt = 1;
% generate Y through the grid of data rates and bandwidths
for rate_idx = 1:length(rate_list)
    rate = rate_list(rate_idx);
    for bandwidth_idx = 1:length(bandwidth_list)
        bandwidth = bandwidth_list(bandwidth_idx);
        
        % display the current grid of data rate and bandwidth
        disp(['rate: ', num2str(rate),' / bandwidth: ', num2str(bandwidth)]);
        
        % simulink parameter
        f_m = rate; % frequency of the single tone message
        % FM modulation via simulink
        simOut = sim('modulation_531.slx'); % Replace 'your_simulink_model' with the actual name of your Simulink model
        % extract the modulation output signal from the Simulink output
        x_c = simOut.yout{1}.Values.Data; % FM modulation output signal
        r_t_timevector = simOut.tout; % time vector of the FM modulation output signal 'x_c'

        for iter_train = 1:n_iter_train

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
            
            % put the obtained SNR into Y_train
            Y_train(1, cnt) = snr;

            cnt = cnt+1;
        end

    end
end

Y_filename = ['Y_train_k-', num2str(k), '.mat'];
save(['./data/', Y_filename], 'Y_train');

%% 3. implementation - GLM

% basic parameters to run the encoding model
X_train = transpose(X_train); % feature matrix (size: num_instances x num_features)
Y_train = transpose(Y_train); % ground truth vector (size: num_instances x 1)
num_folds = 6; % num_folds for cross-validation
lambda_candidates = [0.1:0.1:0.9]; % candidates of regularization parameters

% run the encoding model with obtained X_train and Y_train
[W, Lambda, Rmat] = linear_regression_cv_vw(X_train, Y_train, lambda_candidates, num_folds);

% save the result of the encoding model
W_filename = ['W_k-', num2str(k), '.mat'];
save(['./data/' W_filename], 'W', 'Lambda', 'Rmat');

%% 4. get desired bandwidth

% loop through the candidates of data rates to get the optimum BW of the predetection filter
for rate_idx = 1:length(rate_list)
    rate = rate_list(rate_idx);

    bw_optimal = bw_finder(rate, bandwidth_list, k, W);

    disp(['optimal BW: ', num2str(bw_optimal), ...
      ' / under: data_rate=', num2str(rate), ' & k=', num2str(k)]);
end