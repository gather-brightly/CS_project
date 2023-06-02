function desired_signal = outlier_delete(signal)
    % compute the median and interquartile range (IQR)
    median_val = median(signal);
    iqr_val = iqr(signal);
    % define a threshold for outlier detection (2 times the IQR)
    threshold = 2 * iqr_val;
    % find the indices of outliers
    outlier_indices = find(abs(signal - median_val) > threshold);
    % set outlier values to zero
    signal(outlier_indices) = 0;
    desired_signal = signal;
end