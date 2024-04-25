function [ROI_masks] = ROI_to_mask(rtContours, CTimage, spatial)

    set(0,'DefaultFigureVisible','off');
    
    nameToFind = [ {'Level Ia' }, {'Level Ib L' }, {'Level Ib R' }, ...
        {'Level II L' }, {'Level II R' }, {'Level III L'}, ...
        {'Level III R'}, {'Level IV L' }, {'Level IV R' }, ...
        {'Level V L'  }, {'Level V R'  }];
    
    for i=1:length(rtContours)
        for j=1:length(nameToFind)
            
            contourIndex = int32(find(strcmp(rtContours{i}.ROIs.Name, ...
                nameToFind(j))));
      
            rtMask = createMask(rtContours{i},contourIndex, spatial);
           
            ROI_masks{i}.Name{j} = nameToFind(j);
            ROI_masks{i}.Mask{j} = rtMask;

        end
    end
    set(0,'DefaultFigureVisible','on');
end