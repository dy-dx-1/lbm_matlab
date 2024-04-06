function show_aero_coeffs(f, b_cyl_indices, boundary_links, dh, dt, calc_coeff, rho, pressure_calc_coeff, p_inf, p_divider)
    [cd, cl, cp] = aero_coeffs(f, b_cyl_indices, boundary_links, dh, dt, calc_coeff, rho, pressure_calc_coeff, p_inf, p_divider); 
    % SINCE OPENING IN APPEND MODE MAKE SURE OLD FILES ARE DELETED ELSE YOU
    % WILL JUST APPEND TO OLD DATA!!!! 
    cp = cp.'; % transposing to write [array] at each newline
    fileID1 = fopen('C:/Users/Nicolas/Desktop/PI3_calcs/results_coeffs/CD.txt', 'a');
    fileID2 = fopen('C:/Users/Nicolas/Desktop/PI3_calcs/results_coeffs/CL.txt', 'a');
    fileID3 = fopen('C:/Users/Nicolas/Desktop/PI3_calcs/results_coeffs/CP.txt', 'a');
    fprintf(fileID1, '%f\n', cd); 
    fprintf(fileID2, '%f\n', cl); 
    fprintf(fileID3, '%f ', cp); 
    fprintf(fileID3, '\n'); 
    fclose(fileID1); 
    fclose(fileID2);
    fclose(fileID3);
    %disp([cd, cl]); % not displaying cp cause array too big 
  