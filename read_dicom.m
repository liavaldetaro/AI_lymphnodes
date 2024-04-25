% Function input: folder containing your DICOM CT and rt contours
% outputs CTimage and rtContours
function [CTimage, spatial, dim, rtContours, info] = read_dicom(dataFolder)

    % List all DICOM files in the folder
    dicomFiles = dir(fullfile(string(dataFolder), '*.dcm'));
    
    % Initialize an empty cell array to store DICOM images
    dicomImages = {};
    
    % Initialize an empty cell array to store RT contours
    rtContours = {};
    
    c1 = 1;
    c2 = 1;
    % Loop through each DICOM file and read it
    for i = 1:numel(dicomFiles)
        path = fullfile(dicomFiles(i).folder, dicomFiles(i).name);
        if startsWith(dicomFiles(i).name, 'CT')
            % It's a CT image file
            dicomVpath{c2} = path;
            c2 = c2 + 1;
        else
            % It's an RT structure file
            rtContours{c1} = dicomContours(dicominfo(path,'UseVRHeuristic', false)); % Extract contour data
            info{c1} = dicominfo(path);
            c1 = c1 + 1;
        end
    end
    
    [V,spatial,dim] = dicomreadVolume(dicomVpath);
    CTimage = squeeze(V);
    clear c1 c2
end
