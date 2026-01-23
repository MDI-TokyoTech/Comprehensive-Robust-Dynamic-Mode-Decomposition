function c = redblue(m)
    %REDBLUE    Shades of red and blue color map
    %   REDBLUE(M), is an M-by-3 matrix that defines a colormap.
    %   The colors begin with bright blue, range through white, and end
    %   with bright red.
    %   REDBLUE, by itself, is the same length as the current figure's
    %   colormap. If no figure exists, MATLAB creates one.
    %
    %   For example, to reset the colormap of the current figure:
    %
    %             colormap(redblue)
    %
    %   See also HSV, GRAY, HOT, BONE, COPPER, PINK, FLAG, 
    %   COLORMAP, RGBPLOT.
    
    if nargin < 1, m = size(get(gcf,'colormap'),1); end
    
    if (mod(m,2) == 0)
        % From [0 0 1] to [1 1 1], then [1 1 1] to [1 0 0];
        m1 = m*0.5;
        r = (0:m1-1)'/max(m1-1,1);
        g = r;
        b = r*0 + 1;
        
        % Concatenate
        c = [ [r; ones(m1,1)], [g; flipud(g)], [b; zeros(m1,1)] ];
    else
        m1 = ceil(m*0.5);
        r = (0:m1-1)'/max(m1-1,1);
        g = r;
        b = r*0 + 1;
        
        % Concatenate
        % Middle is white. r ends at 1.
        
        r_bot = r; 
        g_bot = g; 
        b_bot = b;
        
        r_top = ones(m1-1, 1);
        g_top = flipud(g(1:end-1));
        b_top = zeros(m1-1, 1);
        
        c = [ [r_bot; r_top], [g_bot; g_top], [b_bot; b_top] ];
    end
end
