function px = proj_l1ball(X, eta)
    % proj_l1ball : Fast projection onto the l1 ball of radius eta
    % Inputs: X - Input matrix
    %         eta - Radius of the l1 ball
    % Outputs: px - Output matrix after projection

    if sum(abs(X(:))) <= eta
        px = X; % Already inside the l1 ball
        return;
    end

    % The trick used here is that we compute the threshold like above, and if it is negative, 
    % then y is inside the ball, so there is nothing to do and the threshold is set to zero.
    x = X(:);
    x = max(abs(x)-max(max((cumsum(sort(abs(x),1,'descend'),1)-eta)./(1:size(x,1))'),0),0).*sign(x);
    px = reshape(x,size(X));
end