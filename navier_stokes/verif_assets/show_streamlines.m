function show_streamlines(x_sampled, y_sampled, u, v, sample_factor, xstart, ystart)
    fig = gcf; 
    ax = fig.Children; 
    if numel(ax) ~= 2 % checking if were displaying only velocity or velocity and pressure 
        subplot(2,1,1); 
    end 
    sample_u = flipud(u(1:sample_factor:end, 1:sample_factor:end));
    sample_v = flipud(v(1:sample_factor:end, 1:sample_factor:end));
    hold on; 
    streamlines = streamline(x_sampled, y_sampled, sample_u, sample_v, xstart, ystart); 
    set(streamlines, 'color', 'black'); 
    hold off; 
end 