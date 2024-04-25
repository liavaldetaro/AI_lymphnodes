dir_list = [{'DAH01'}, {'DAH02'}, {'DAH03'}, {'DAH04'}, {'DAH05'}, ...
    {'DAH07'}, {'DAH08'}, {'DAH09'}, {'DAH10'}, {'DAH11'}, {'DAH12'}, {'DAH13'}, ...
    {'DAH14'}, {'DAH15'}, {'DAH16'}, {'DAH17'}, {'DAH18'}, {'DAH19'}, {'DAH20'}];

MDS_per_patient = [];
S_per_patient =[];
dice_median = [];


nameToFind = [ {'Level Ia' }, {'Level Ib L' }, {'Level Ib R' }, ...
        {'Level II L' }, {'Level II R' }, {'Level III L'}, ...
        {'Level III R'}, {'Level IV L' }, {'Level IV R' }, ...
        {'Level V L'  }, {'Level V R'  }];


for loop=2:length(nameToFind)

    loop

    for i = 1:1%length(dir_list)
        
        dir_name = append('\\prapprflstg01\Research\LBV\Level_data_MIM\', dir_list(i));
        [CTimage, spatial, dim, rtContours, info] = read_dicom(dir_name);
        [ROI_masks] = ROI_to_mask(rtContours, CTimage, spatial);
            
        c = 1;
        m = 1;    
        
        j = loop; % LEVEL WE ARE CALCULATING
          
        % First we check which contour has the lowest DICE (cheaper calculation) to all others
        for k=1:length(rtContours)
            for l=k+1:length(rtContours)
                
                [x, y, z] = find_edges(ROI_masks{1}.Mask{j});
                ROI1 = ROI_masks{k}.Mask{j};
                ROI1 = ROI1(x(1):x(end), y(1):y(end), z(1):z(end));
            
                ROI2 = ROI_masks{l}.Mask{j};
                ROI2 = ROI2(x(1):x(end), y(1):y(end), z(1):z(end));
        
                dice_median(m, 1) = dice(ROI1, ROI2);
                dice_median(m, 2) = k;
                dice_median(m, 3) = l;
                m = m + 1;

            end 
        end

   end
        
   % find the contour with lowest DICE to all others:
        
   sorted =  sortrows(dice_median, 1);
   best = mode(reshape(sorted(1:5, 2:3), 1, []));
        
   MDS_per_patient = [];
   S_per_patient =[];
        
   for i = 1:length(dir_list)
        
       dir_name = append('\\', dir_list(i));
       [CTimage, spatial, dim, rtContours, info] = read_dicom(dir_name);
       [ROI_masks] = ROI_to_mask(rtContours, CTimage, spatial);
            
       c = 1;
       MDS_mean = [];
       S_mean = [];
        
       for j = loop:loop%length(nameToFind) % number of ROIs
                
           [x, y, z] = find_edges(ROI_masks{best}.Mask{j});
           ROI1 = ROI_masks{best}.Mask{j};
           ROI1 = ROI1(x(1):x(end), y(1):y(end), z(1):z(end));
        
           for k=1:length(rtContours) % number of doctors drawing
               if k == best
                   continue
               end
        
               ROI2 = ROI_masks{k}.Mask{j};
               ROI2 = ROI2(x(1):x(end), y(1):y(end), z(1):z(end));
        
               [dta, dtb, S, S_prime] = Surface_distances_3D(ROI1, ROI2, spatial);
        
               MDS = dtb;
        
               if isempty(MDS_mean)
                   MDS_mean = MDS;
               else
                   MDS_mean = MDS_mean + MDS;
               end
                    
               c = c + 1;
           end
       end
        
       MDS_per_patient{i} =  MDS_mean ./ c;
       S_per_patient{i} = S;
        
   end
        
   R_MDS_per_patient = [];
   R_MEAN = MDS_per_patient{1};
        
        
   [optimizer,metric] = imregconfig("multimodal")
   optimizer.InitialRadius = 0.001;
   optimizer.Epsilon = 1.5e-4;
   optimizer.GrowthFactor = 1.01;
   optimizer.MaximumIterations = 400;
        
   for i = 2:length(dir_list)
        
       fixed = int8(S_per_patient{1});
       moving = int8(S_per_patient{i});
           
       tform = imregtform(moving, fixed, "affine", optimizer, metric);
       
       R_MDS_per_patient{i} = imwarp(MDS_per_patient{i}, tform, "OutputView", imref3d(size(fixed)));
       R_contour_per_patient{i} = imwarp(S_per_patient{i}, tform, "OutputView", imref3d(size(fixed)));
        
   end
        
        
   R_contour_per_patient{1} = S_per_patient{1};
   R_MDS_per_patient{1} = MDS_per_patient{1};
        
        
   v_x = spatial.PixelSpacings(1);
   v_y = spatial.PixelSpacings(2);
   v_z = spatial.PatientPositions(2, 3) - spatial.PatientPositions(1, 3);
        
   A = [v_x 0 0 0; 0 v_y 0 0; 0 0 v_z 0; 0 0 0 1];
   tform = affinetform3d(A);
        
        
   MDS_median = R_MDS_per_patient{1};
        
   for i=2:length(dir_list)
       temp = R_MDS_per_patient{i};
       temp(isnan(temp)) = 0;
       temp(isinf(temp)) = 0;
       MDS_median = MDS_median + temp;
   end
       
   MDS_median = MDS_median ./19;
        
   SURF = isosurface(uint8(R_contour_per_patient{1}), 0);
   DATA_SURF = interp3(MDS_median,SURF.vertices(:,1),SURF.vertices(:,2),SURF.vertices(:,3));
        

   file_name = string(append('\\', nameToFind(loop), '.mat'));
   save(file_name, 'SURF', 'MDS_median');

   disp('FINISHED A RUN')
        
end
