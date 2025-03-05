%% 05/03/2025 - Combine targets together. As described in Goetz meeting on 04/03/2025

%% Load Amygdala DLPFC target
[fileName, filePath] = uigetfile('*.nii', 'Select a NIfTI file');
if fileName ~= 0
    fullFilePath = fullfile(filePath, fileName);
    disp(['Selected file: ', fullFilePath]);

    amyg_targ = nifti_load(fullFilePath);
    
else
    disp('File selection canceled.');
end
clear fileName filePath fullFilePath;
%% Load Subcallosal FOX target
[fileName, filePath] = uigetfile('*.nii', 'Select a NIfTI file');
if fileName ~= 0
    fullFilePath = fullfile(filePath, fileName);
    disp(['Selected file: ', fullFilePath]);

    FOX_targ = nifti_load(fullFilePath);

else
    disp('File selection canceled.');
end

clear fileName filePath fullFilePath;
%% Load VMPFC target
[fileName, filePath] = uigetfile('*.nii', 'Select a NIfTI file');
if fileName ~= 0
    fullFilePath = fullfile(filePath, fileName);
    disp(['Selected file: ', fullFilePath]);

    VMPFC_targ = nifti_load(fullFilePath);
    
else
    disp('File selection canceled.');
end
clear fileName filePath fullFilePath;
%% Extract the data from each of the NIFTI files
data_amyg_tar = amyg_targ.vol;
data_FOX_tar = FOX_targ.vol;
data_VMPFC_tar = VMPFC_targ.vol;

%% Add the targets up --> Combined Target
data_com_targ = data_amyg_tar + data_FOX_tar + data_VMPFC_tar;

%% Freq. Plot (to check)

voxelData = data_com_targ(:);
voxelDataNoZero = voxelData(voxelData ~= 0);

[uniqueValues, ~, indices] = unique(voxelDataNoZero);

% Calculate the frequency of each unique value
frequencies = histcounts(voxelDataNoZero, [uniqueValues; uniqueValues(end) + 1]);

% Create the frequency plot
figure;
bar(uniqueValues, frequencies, 'FaceColor', [0.7, 0.7, 0.7]);

% Customize the plot
xlabel('Voxel Intensity');
ylabel('Frequency');
title('Frequency of Voxel Intensities in hdr.vol');
grid on;
    % Optionally, display how many times a specific value appears, e.g., for value '4'

disp(sum(voxelDataNoZero));  
%%
com_targ = amyg_targ;
com_targ.vol = data_com_targ;

%% Save Combined Target
nifti_save(com_targ, 'combined_target_111.nii'); % saves the combined target (with 111 weightings) as NIFTI file