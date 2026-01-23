function draw_cylinder()
    % シリンダーの描写
    t_cyl = (1:100)/100'*2*pi;
    y_cyl = 99 + 25*cos(t_cyl); 
    x_cyl = 49 + 25*sin(t_cyl);
    
    fill(x_cyl, y_cyl, [.3 .3 .3]); 
    plot(x_cyl, y_cyl, 'k', 'LineWidth', 1.2);
end