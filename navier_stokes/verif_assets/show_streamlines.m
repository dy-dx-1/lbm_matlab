function show_streamlines(x_sampled, y_sampled, u, v, sample_factor)
    fig = gcf; 
    ax = fig.Children; 
    if numel(ax) ~= 2 % checking if were displaying only velocity or velocity and pressure 
        subplot(2,1,1); 
    end 
    sample_u = u(1:sample_factor:end, 1:sample_factor:end);
    sample_v = v(1:sample_factor:end, 1:sample_factor:end);
    hold on; 
    ss = streamslice(x_sampled, y_sampled, sample_u, sample_v);
    set(ss, 'Color', 'k'); 
    hold off; 
end 