function [dta, dtb, S, S_prime] = Surface_distances_3D(ROI1, ROI2, spatial)

    conn1 = [0, 0 ,0; 0, 1, 0; 0, 0, 0];
    conn2 = [0, 1, 0; 1, 1, 1 ;0, 1, 0];
    conn3 = [0, 0, 0; 0, 1, 0; 0, 0, 0];
    
    Z = cat(3, conn1, conn2, conn3);
    
    border1 = bwperim(ROI1, Z);
    border2 = bwperim(ROI2, Z);

    v_x = spatial.PixelSpacings(1);
    v_y = spatial.PixelSpacings(2);
    v_z = spatial.PatientPositions(2, 3) - spatial.PatientPositions(1, 3);

    dta = bwdistsc1(border1, [v_x, v_y, v_z]);
    dtb = bwdistsc1(border2, [v_x, v_y, v_z]);
    
    S = border1;
    S_prime = border2;
    
    %dta = dta(S_prime ~= 0);
    %dtb = dtb(S ~= 0);

    %dta(S_prime == 0) = 0;
    %dtb(S == 0) = 0;
    
    %dta = reshape(dta, [1, length(dta)]);
    %dtb = reshape(dtb, [1, length(dtb)]);
    
    %sds = abs(dta-dtb);
    
end