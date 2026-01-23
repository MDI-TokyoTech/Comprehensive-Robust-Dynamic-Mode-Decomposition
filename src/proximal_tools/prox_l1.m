function px = prox_l1(x, gamma)
    % prox_l1 : Proximal operator for the l1 norm
    % Inputs: x - Input matrix
    %         gamma - Step size parameter
    % Outputs: px - Output matrix after applying the proximal operator

    px = sign(x) .* max(abs(x) - gamma, 0);
end