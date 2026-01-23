function px = prox_l12(x, gamma)
    % prox_l12 : Proximal operator for the l12 norm
    % Inputs: x - Input matrix (3*m x n)
    %         gamma - Step size parameter
    % Outputs: px - Output matrix after applying the proximal operator (3*m x n)

    N = size(x, 1) / 3; % Number of pixels (m)
    Dv = x(1:N, :);         % Vertical differences
    Dh = x(N+1:2*N, :);     % Horizontal differences
    Dt = x(2*N+1:end, :);   % Temporal differences

    % Compute the shrinkage for each pixel
    shrinkage = max(1 - gamma ./ sqrt(Dv.^2 + Dh.^2 + Dt.^2 + 1e-10), 0);

    % Apply the shrinkage
    px = [Dv .* shrinkage; ...
          Dh .* shrinkage; ...
          Dt .* shrinkage];
end