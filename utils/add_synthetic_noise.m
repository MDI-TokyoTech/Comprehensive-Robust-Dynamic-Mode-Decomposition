function [X_noisy, X_clean_gt] = add_synthetic_noise(X, sigma, outlier_ratio)
    % Add synthetic Gaussian and Salt-and-Pepper noise to the data
    fprintf('Adding synthetic noise...\n');

    % 1. Normalize to [-0.5, 0.5]
    min_val = min(X(:));
    max_val = max(X(:));
    X_clean_gt = X / max(abs(min_val), abs(max_val)) * 0.5;

    % 2. Add Gaussian noise
    noise_g = sigma * randn(size(X_clean_gt));
    X_noisy = X_clean_gt + noise_g;

    % 3. Add Salt-and-Pepper noise (Outliers)
    % Values set to 0.5 or -0.5 based on outlier_ratio
    Sp = 0.5 * ones(size(X_noisy));
    
    % Use imnoise if available, otherwise manual implementation could be used
    if exist('imnoise', 'file')
        Sp = imnoise(Sp, 'salt & pepper', outlier_ratio);
    else
        warning('Image Processing Toolbox not found. Using simple random noise for outliers.');
        mask = rand(size(X_noisy)) < outlier_ratio;
        Sp(mask) = (rand(sum(mask(:)), 1) > 0.5); % simple binary noise
    end
    
    X_noisy(Sp == 1) = 0.5;
    X_noisy(Sp == 0) = -0.5;

    fprintf('  - Normalized range: [-0.5, 0.5]\n');
    fprintf('  - Gaussian Noise: sigma = %.2f\n', sigma);
    fprintf('  - Outliers: %.1f%% (Values set to +/- 0.5)\n', outlier_ratio * 100);
end