function x = adj_diff_st(z, w, h, v, M)
    % adj_diff_st : Compute Adjoint of Spatial-Temporal Differences
    % Inputs: z - Input difference matrix (3*m x n)
    %         w - Weight between spatial and temporal differences
    %         h - number of rows in the spatial dimension
    %         v - number of columns in the spatial dimension
    %         M - number of time frames
    % Outputs: x - Output data matrix (m x n)

    num_pix = h * v;
    zv = reshape(z(1:num_pix, :), [v, h, M]);
    zh = reshape(z(num_pix+1:2*num_pix, :), [v, h, M]);
    zt = reshape(z(2*num_pix+1:end, :), [v, h, M]);
    
    % Adjoint of vertical differences
    dvT = zv; 
    dvT(1,:,:) = -zv(1,:,:);
    dvT(2:end-1,:,:) = zv(1:end-2,:,:) - zv(2:end-1,:,:);
    dvT(end,:,:) = zv(end-1,:,:);
    
    % Adjoint of horizontal differences
    dhT = zh;
    dhT(:,1,:) = -zh(:,1,:);
    dhT(:,2:end-1,:) = zh(:,1:end-2,:) - zh(:,2:end-1,:);
    dhT(:,end,:) = zh(:,end-1,:);

    % Adjoint of temporal differences
    dtT = zt;
    dtT(:,:,1) = -zt(:,:,1);
    dtT(:,:,2:end-1) = zt(:,:,1:end-2) - zt(:,:,2:end-1);
    dtT(:,:,end) = zt(:,:,end-1);

    % Combine adjoint differences with weights
    x = w * (dvT + dhT) + (1 - w) * dtT;
    x = reshape(x, [num_pix, M]);
end