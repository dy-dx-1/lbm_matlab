function show_density(rho)
    %AVERAGE DENSITY VISUALISATION
    average_density = mean(rho, 'all');
    fileID = fopen('C:/Users/Nicolas/Downloads/average_density.txt', 'a');
    fprintf(fileID, '%f\n', average_density); 
    fclose(fileID); 
    disp(average_density); % displaying average density to check for conservation