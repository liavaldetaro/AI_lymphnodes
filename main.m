%% Main 

% load the dicom files

dir_list = [{'DAH01'}, {'DAH02'}, {'DAH03'}, {'DAH04'}, {'DAH05'},  {'DAH07'}, ...
     {'DAH07'}, {'DAH08'}, {'DAH09'}, {'DAH10'}, {'DAH11'}, {'DAH12'}, {'DAH13'}, ...
    {'DAH14'}, {'DAH15'}, {'DAH16'}, {'DAH17'}, {'DAH18'}, {'DAH19'}, {'DAH20'}];


nameToFind = [ {'Level Ia' }, {'Level Ib L' }, {'Level Ib R' }, ...
        {'Level II L' }, {'Level II R' }, {'Level III L'}, ...
        {'Level III R'}, {'Level IV L' }, {'Level IV R' }, ...
        {'Level V L'  }, {'Level V R'  }];

%all_contours = [];
%all_masks = [];
for i = 1:length(dir_list)

    dir_name = append('\\prapprflstg01\Research\LBV\Level_data_MIM\', dir_list(i))
    [CTimage, spatial, dim, rtContours, info] = read_dicom(dir_name);
    [ROI_masks] = ROI_to_mask(rtContours, CTimage, spatial);
    
    %all_contours{i} = rtContours;
    %all_masks{i} = ROI_masks;

    %save(char(append('\\Rsyd.net\homedir$\0044\ten9ti\My Documents\Uddannelse\dahanca levels\Data_analysis\', append(dir_list(i), '.mat'))), 'ROI_masks', '-v7.3')
    %save(char(append('\\Rsyd.net\homedir$\0044\ten9ti\My Documents\Uddannelse\dahanca levels\Data_analysis\', append(dir_list(i), 'spatial.mat'))), 'spatial')

    for j = 1:length(nameToFind) % number of ROIs
       
        c = 1;

        dice_median = [];
        jaccard_median = [];
        bfscore_median = [];

        for k=1:length(rtContours) % number of doctors drawing
            
            for l=k+1:length(rtContours)

                dice_median(c) = dice(ROI_masks{k}.Mask{j}, ROI_masks{l}.Mask{j});
                jaccard_median(c) = jaccard(ROI_masks{k}.Mask{j}, ROI_masks{l}.Mask{j});
                bfscore_median(c) = bfscore(ROI_masks{k}.Mask{j}, ROI_masks{l}.Mask{j});
                [MSD, RMS, HD, HD95] = Surface_distances(ROI_masks{k}.Mask{j}, ROI_masks{l}.Mask{j}, spatial);
                MSD_median(c) = MSD;
                RMS_median(c) = RMS;
                HD_median(c) = HD;
                HD95_median(c) = HD95;
                % median og interq. range instead
                
                c = c + 1;
            end

        end

        similarity.patient(i,j) = dir_list(i);
        similarity.struct(i,j) = nameToFind(j);

        similarity.dice(i,j) = median(nonzeros(dice_median));
        similarity.jaccard(i,j) = median(nonzeros(jaccard_median));
        similarity.bfscore(i,j) = median(nonzeros(bfscore_median));

        similarity.MSD(i,j) = median(nonzeros(MSD_median));
        similarity.RMS(i,j) = median(nonzeros(RMS_median));
        similarity.HD(i,j) = median(nonzeros(HD_median));
        similarity.HD95(i,j) = median(nonzeros(HD95_median));

        similarity.dice_iqr(i,j) = iqr(nonzeros(dice_median));
        similarity.jaccard_iqr(i,j) = iqr(nonzeros(jaccard_median));
        similarity.bfscore_iqr(i,j) = iqr(nonzeros(bfscore_median)); 

        similarity.MSD_iqr(i,j) = iqr(nonzeros(MSD_median));
        similarity.RMS_iqr(i,j) = iqr(nonzeros(RMS_median));
        similarity.HD_iqr(i,j) = iqr(nonzeros(HD_median));
        similarity.HD95_iqr(i,j) = iqr(nonzeros(HD95_median));

    end


    clear CTimage rtContours ROI_masks

    %clear ans dim dir_name i info rtContours spatial CTimage
    %save(string(append('\\prapprflstg01\Research\LBV\Masks\', dir_list(i))))

end