function standardized_signal = standardize(signal)
    % calculate minimum and maximum values
    min_val = min(signal);
    max_val = max(signal);
    % shift the range to start from 0
    shifted_signal = signal - min_val;
    % normalize by dividing by the range
    range = max_val - min_val;
    normalized_signal = shifted_signal / range;
    % scale the signal by multiplying by 2
    scaled_signal = normalized_signal * 2;
    % shift the range to -1 to 1
    standardized_signal = scaled_signal - 1;
end