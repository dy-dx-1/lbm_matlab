function show_aero_coeffs(f, b_cyl_indices, boundary_links, dh, dt, calc_coeff)
    [cd, cl] = aero_coeffs(f, b_cyl_indices, boundary_links, dh, dt, calc_coeff); 
    fileID1 = fopen('C:/Users/Nicolas/Downloads/CD.txt', 'a');
    fileID2 = fopen('C:/Users/Nicolas/Downloads/CL.txt', 'a');
    fprintf(fileID1, '%f\n', cd); 
    fprintf(fileID2, '%f\n', cl); 
    fclose(fileID1); 
    fclose(fileID2);
    disp([cd, cl]); 