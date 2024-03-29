function show_streamlines(x_sampled, y_sampled, u, v, sample_factor, xstart, ystart)
    sample_u = flipud(u(1:sample_factor:end, 1:sample_factor:end));
    sample_v = flipud(v(1:sample_factor:end, 1:sample_factor:end));
    streamlines = streamline(x_sampled, y_sampled, sample_u, sample_v, xstart, ystart); 
    set(streamlines, 'color', 'black'); 
end 