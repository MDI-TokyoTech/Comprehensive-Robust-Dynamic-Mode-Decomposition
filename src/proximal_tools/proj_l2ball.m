function px = proj_l2ball(x, c, epsilon)
    % proj_l2ball : Project onto the l2 ball of radius epsilon centered at c
    % Inputs: x - Input matrix
    %         c - Center of the l2 ball
    %         epsilon - Radius of the l2 ball
    % Outputs: px - Output matrix after projection

    dist = x - c;
    norm_dist = norm(dist, 'fro');
    
    if norm_dist <= epsilon
        px = x; % Already inside the ball
    else
        px = c + (epsilon * dist) / norm_dist; % Project onto the boundary
    end
end 