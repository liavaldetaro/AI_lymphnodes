%% 3D Surface plots 

dir_list = [{'DAH01'}, {'DAH02'}, {'DAH03'}, {'DAH04'}, {'DAH05'},{'DAH06'}, ...
    {'DAH07'}, {'DAH08'}, {'DAH09'}, {'DAH10'}, {'DAH11'}, {'DAH12'}, {'DAH13'}, ...
    {'DAH14'}, {'DAH15'}, {'DAH16'}, {'DAH17'}, {'DAH18'}, {'DAH19'}, {'DAH20'}];

nameToFind = [ {'Level Ia' }, {'Level Ib L' }, {'Level Ib R' }, ...
        {'Level II L' }, {'Level II R' }, {'Level III L'}, ...
        {'Level III R'}, {'Level IV L' }, {'Level IV R' }, ...
        {'Level V L'  }, {'Level V R'  }];



for i = 1:length(dir_list)

    dir_name = append('\\', dir_list(i))
    [CTimage, spatial, dim, rtContours, info] = read_dicom(dir_name);
    [ROI_masks] = ROI_to_mask(rtContours, CTimage, spatial);
    
    c = 1;
    MDS_mean = [];
    S_mean = [];
    for j = 10:10%length(nameToFind) % number of ROIs
        
        [x, y, z] = find_edges(ROI_masks{1}.Mask{j});

        for k=2:length(rtContours) % number of doctors drawing

            ROI1 = ROI_masks{1}.Mask{j};
            ROI1 = ROI1(x(1):x(end), y(1):y(end), z(1):z(end));

            ROI2 = ROI_masks{k}.Mask{j};
            ROI2 = ROI2(x(1):x(end), y(1):y(end), z(1):z(end));

            [dta, dtb, S, S_prime] = Surface_distances_3D(ROI1, ROI2, spatial);

            MDS = dta;
            MDS(int8(S_prime) == 0) = 0;
            MDS_mean{c} = MDS;
            S_mean{1} = S;
            S_mean{c+1} = S_prime;
            c = c + 1;
        end
        
        mean = MDS_mean{1};
        Smean = S_mean{1};
        for s = 2:c-1
            mean = mean + MDS_mean{s};
            Smean = Smean + S_mean{s};
        end

        MDS_per_patient{i} = mean./(c-1);
        contour_per_patient{i} = Smean./(c-1);
    end
end


R_MDS_per_patient = [];
R_MEAN = MDS_per_patient{1};


[optimizer,metric] = imregconfig("multimodal")
optimizer.InitialRadius = 0.009;
optimizer.Epsilon = 1.5e-9;
optimizer.GrowthFactor = 1.01;
optimizer.MaximumIterations = 300;

for i = 2:length(dir_list)

    fixed = int8(contour_per_patient{1});
    moving = int8(contour_per_patient{i});
   
    tform = imregtform(moving, fixed, "affine", optimizer, metric);

    R_MDS_per_patient{i} = imwarp(MDS_per_patient{i}, tform, "OutputView", imref3d(size(fixed)));
    R_contour_per_patient{i} = imwarp(contour_per_patient{i}, tform, "OutputView", imref3d(size(fixed)));

    R_MEAN = R_MEAN + R_MDS_per_patient{i};

end

R_MEAN = R_MEAN./20;


%%

