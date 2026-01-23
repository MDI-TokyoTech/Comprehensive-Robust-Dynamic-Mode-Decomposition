function z = diff_st(x, w, h, v)
    % diff_st : Compute Spatial-Temporal Differences
    % Inputs: x - Input data matrix (m x n)
    %         w - Weight between spatial and temporal differences
    %         h - number of rows in the spatial dimension
    %         v - number of columns in the spatial dimension
    % Outputs: z - Output difference matrix (3*m x n)

    [~, M] = size(x);
    X = reshape(x, [v, h, M]); % Reshape x into 3D array

    % Compute vertical differences Dv
    dv = diff(X, 1, 1); 
    dv(end+1, :, :) = 0; % Neumann boundary

    % Compute horizontal differences Dh
    dh = diff(X, 1, 2);
    dh(:, end+1, :) = 0; % Neumann boundary

    % Compute temporal differences Dt
    dt = diff(X, 1, 3);
    dt(:, :, end+1) = 0; % Neumann boundary

    % Combine differences with weights
    z = [w * reshape(dv, [], M); ...
         w * reshape(dh, [], M); ...
         (1-w) * reshape(dt, [], M)];
end