function visualize_dimensional_reduction(X_clean_gt, xi, S_outliers_dr, Phi, Lambda, h, v)
    % Visualize dimensional reduction results
    [~, M] = size(X_clean_gt);
    snapshot_idx = round(M / 2);
    
    % Count non-zero xi
    non_zero_xi = sum(abs(xi) > 0);
    fprintf('Number of non-zero xi: %d / %d\n', non_zero_xi, length(xi));

    % Reconstruct snapshots for visualization
    C_dr = Lambda .^ (0:M-1);
    X_recon_dr = real(Phi * diag(xi) * C_dr);
    S_outliers_dr = real(S_outliers_dr);

    x_gt_snapshot = reshape(X_clean_gt(:, snapshot_idx), v, h);
    x_recon_dr_snapshot = reshape(X_recon_dr(:, snapshot_idx), v, h);
    s_dr_snapshot = reshape(S_outliers_dr(:, snapshot_idx), v, h);

    % Visualize Dimensional Reduction
    figure('Name', 'CR-DMD Dimensional Reduction Result', 'Position', [100, 600, 1200, 400]);
    c_scale = [-0.5, 0.5];

    % Define custom colormap
    colors_recon = [0 0 1; 1 1 1; 1 0 0]; 
    map_recon = interp1(linspace(-1, 1, 3), colors_recon, linspace(-1, 1, 256));

    subplot(1, 3, 1);
    imagesc(real(x_gt_snapshot));
    hold on; draw_cylinder(); hold off;
    clim(c_scale);
    colormap(map_recon); 
    title('Ground Truth');
    axis off equal;
    colorbar;

    subplot(1, 3, 2);
    imagesc(x_recon_dr_snapshot);
    hold on; draw_cylinder(); hold off;
    clim(c_scale);
    title(sprintf('Low-dim Rep(Non-zero xi: %d)', non_zero_xi));
    axis off equal;
    colorbar;

    subplot(1, 3, 3);
    imagesc(s_dr_snapshot);
    hold on; draw_cylinder(); hold off;
    clim(c_scale);
    title('Separated Outliers');
    axis off equal;
    colorbar;

    sgtitle('Dimensional Reduction Results');
end