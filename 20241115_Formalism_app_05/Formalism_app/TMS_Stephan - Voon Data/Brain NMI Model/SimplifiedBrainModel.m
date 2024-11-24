%This code attempts to identify the vector direction as outlined by S.Goetz


%% 
% Load NIfTI file for the Brain model data
[fileName, filePath] = uigetfile('*.nii', 'Select a NIfTI file');
if fileName ~= 0
    fullFilePath = fullfile(filePath, fileName);
    disp(['Selected file: ', fullFilePath]);

    % Load the file using the path
    hdr = nifti_load(fullFilePath);
else
    disp('File selection canceled.');
end
%% Load NIfTI file for the amydala ROI (Voon Data)
[fileName_mask, filePath_mask] = uigetfile('*.nii', 'Select a NIfTI file');
if fileName_mask ~= 0
    filePath_mask = fullfile(filePath_mask, fileName_mask);
    disp(['Selected file: ', filePath_mask]);

    % Load the file using the path
    hdr_mask = nifti_load(filePath_mask);
else
    disp('File selection canceled.');
end

%% Plotting the Voxel Data - using locations at voxel indices
%define the relevant brain region using the mask
brainmodel = hdr.vol;
mask = hdr_mask.vol;
brainmodel(mask == 0) = 0;
volData = brainmodel;
volshow(volData);
% Convert matrix in to grayscale image & Normalize intensity values to range [0, 1]
% volData = mat2gray(volData); 
% volData(volData < 10) = 0; % Set low-intensity values to transparent

% test = volData(:);  
% disp(max(test));
% volshow(mask);
%% 
