function visualize_preprocessing(crdmd, X_clean_gt, X_noisy, X_clean, S_outliers, Lambda, h, v, sigma, outlier_ratio)
    % Visualize preprocessing results and eigenvalues
    fprintf('Visualizing results...\n');
    [~, M] = size(X_noisy);

    % Compute Ground Truth DMD for comparison
    fprintf('Computing Ground Truth DMD...\n');
    [~, Lambda_gt, ~] = crdmd.dmd(X_clean_gt);

    % Plot Eigenvalues
    figure('Name', 'DMD Eigenvalues Comparison', 'Position', [1200, 100, 600, 600]);
    theta = linspace(0, 2*pi, 100);
    plot(cos(theta), sin(theta), 'k--'); hold on;
    plot(real(Lambda_gt), imag(Lambda_gt), 'go', 'MarkerFaceColor', 'g', 'DisplayName', 'Ground Truth');
    plot(real(Lambda), imag(Lambda), 'rx', 'LineWidth', 2, 'DisplayName', 'CR-DMD');
    xlabel('Real(\lambda)');
    ylabel('Imag(\lambda)');
    title('Eigenvalues Comparison');
    axis equal;
    legend('Unit Circle', 'Ground Truth', 'CR-DMD', 'Location', 'best');
    grid on;

    % Select a snapshot to visualize
    snapshot_idx = round(M / 2);

    x_gt_snapshot = reshape(X_clean_gt(:, snapshot_idx), v, h);
    x_noisy_snapshot = reshape(X_noisy(:, snapshot_idx), v, h);
    x_clean_snapshot = reshape(X_clean(:, snapshot_idx), v, h);
    s_snapshot = reshape(S_outliers(:, snapshot_idx), v, h);

    figure('Name', 'CR-DMD Preprocessing Result', 'Position', [100, 100, 1500, 400]);

    % Helper to plot cylinder flow
    c_scale = [-0.5, 0.5]; 

    % Define custom colormap
    colors_recon = [0 0 1; 1 1 1; 1 0 0]; 
    map_recon = interp1(linspace(-1, 1, 3), colors_recon, linspace(-1, 1, 256));

    subplot(1, 4, 1);
    imagesc(real(x_gt_snapshot));
    hold on; draw_cylinder(); hold off;
    clim(c_scale);
    colormap(map_recon); 
    title('Ground Truth (Normalized)');
    axis off equal; 
    colorbar;

    subplot(1, 4, 2);
    imagesc(real(x_noisy_snapshot));
    hold on; draw_cylinder(); hold off;
    clim(c_scale);
    title(sprintf('Noisy Input (sigma=%.2f, sp=%.0f%%)', sigma, outlier_ratio*100));
    axis off equal;
    colorbar;

    subplot(1, 4, 3);
    imagesc(real(x_clean_snapshot));
    hold on; draw_cylinder(); hold off;
    clim(c_scale);
    title('Cleaned Data');
    axis off equal;
    colorbar;

    subplot(1, 4, 4);
    imagesc(real(s_snapshot));
    hold on; draw_cylinder(); hold off;
    clim(c_scale);
    title('Separated Outliers');
    axis off equal;
    colorbar;

    sgtitle('CR-DMD Preprocessing: Cylinder Wake with Synthetic Noise');
end