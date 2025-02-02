%Adapted code on 05/01/2025
% The code uses the gradient function to find the 

%% Load NIfTI file for the Brain model data - use the 2mm Brain model for comparison with Voon Data
[fileName, filePath] = uigetfile('*.nii', 'Select a NIfTI file');
if fileName ~= 0
    fullFilePath = fullfile(filePath, fileName);
    disp(['Selected file: ', fullFilePath]);

    % Load the file using the path
    hdr = nifti_load(fullFilePath);
else
    disp('File selection canceled.');
end
%% Load NIfTI file for the DLPFC ROI MASK (Voon Data)
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
% volshow(brainmodel);
%% 
% brainmodel(mask == 0) = brainmodel(mask == 0) * 10; %for visualisation - makes the NOT ROI darker
brainmodel(mask == 1) = brainmodel(mask == 1) * 3;
volData = brainmodel;
viewer = volshow(volData,RenderingStyle="VolumeRendering");
viewer.Parent.BackgroundColor = [1 1 1];


%% Generate a Binary Mask
% This is redundant is the DLPFC mask is used.
binary_mask = mask ~= 0; % should be the same as mask if a mask is used 
ROI = binary_mask;
%volshow(ROI);
% masked_brain = brainmodel .* binary_mask; 
%volshow(masked_brain); 

%% 
brain_mask = brainmodel > 0.001;
volshow(brain_mask);

% Check if all the ROI mask elements are within the brain mask
is_contained = all(brain_mask(binary_mask > 0) == 1);
disp(is_contained);

% Combine the masks for visualization: brain mask is green, ROI mask is red
overlay = brain_mask;
overlay(binary_mask > 0) = 100;  % Set ROI voxels to 2 for differentiation

% Display the result using a volume viewer (or any visualization method)
volshow(overlay);  % You can use volshow, or any other visualization you prefer

% voxel indices where the ROI mask is non-zero
[roi_x, roi_y, roi_z] = ind2sub(size(binary_mask), find(binary_mask > 0));

% Check if these indices are inside the brain mask
for i = 1:length(roi_x)
    if brain_mask(roi_x(i), roi_y(i), roi_z(i)) == 0
        disp(['ROI voxel at (' num2str(roi_x(i)) ', ' num2str(roi_y(i)) ', ' num2str(roi_z(i)) ') is outside the brain mask.']);
    end
end

%% 
%% 
boundary_mask = bwperim(brain_mask,26); %returns the boundary voxels of brain_mask
% 26 indicates the 3D connectivity
volshow(boundary_mask);

% Logical check: boundary_mask is entirely within brain_mask
is_within = all(brain_mask(boundary_mask));
if is_within
    disp('The boundary_mask is entirely within the brain_mask.');
else
    disp('The boundary_mask is NOT entirely within the brain_mask.');
end

%% 
%% 
% Assumes `mask` is the binary ROI mask (1 for ROI, 0 elsewhere)
% `origin` is the brain origin (e.g., center of the brain volume)
origin = [size(brain_mask, 1)/2, size(brain_mask, 2)/2, size(brain_mask, 3)/2]; % Brain center
[roi_x, roi_y, roi_z] = ind2sub(size(mask), find(mask > 0.1)); % ROI voxel coordinates

roi_centroid = mean([roi_x, roi_y, roi_z], 1);% ROI centroid
% %% 
% 
% % Parameters for the cone
% cone_height = 5;
% max_cone_radius = 2;
% [roi_x, roi_y, roi_z] = ind2sub(size(mask), find(mask));
% roi_centroid = mean([roi_x, roi_y, roi_z], 1);
% 
% [X, Y, Z] = ndgrid(-max_cone_radius:max_cone_radius);
% sphere_mask = (X.^2 + Y.^2 + Z.^2) <= max_cone_radius^2;
% sphere_indices = find(sphere_mask);
% [sx, sy, sz] = ind2sub(size(sphere_mask), sphere_indices);
% sphere_offsets = [sx - (max_cone_radius + 1), sy - (max_cone_radius + 1), sz - (max_cone_radius + 1)];
% 
% extended_mask = mask;
% for i = 1:length(roi_x)
%     radial_vector = ([roi_x(i), roi_y(i), roi_z(i)] - roi_centroid) / norm([roi_x(i), roi_y(i), roi_z(i)] - roi_centroid);
%     for step = 1:cone_height
%         new_pos = round([roi_x(i), roi_y(i), roi_z(i)] + step * radial_vector);
% 
%         valid_offsets = bsxfun(@plus, sphere_offsets, new_pos);
%         valid_indices = sub2ind(size(mask), valid_offsets(:,1), valid_offsets(:,2), valid_offsets(:,3));
%         valid_indices = valid_indices(valid_offsets(:,1) > 0 & valid_offsets(:,1) <= size(mask,1) & ...
%                                       valid_offsets(:,2) > 0 & valid_offsets(:,2) <= size(mask,2) & ...
%                                       valid_offsets(:,3) > 0 & valid_offsets(:,3) <= size(mask,3));
%         extended_mask(valid_indices) = 1;
% 
%     end
% end
% % volshow(extended_mask, 'Colormap', 'red', 'Lighting', 'Gouraud');
% 
% 
% % Visualize the extended ROI (cone shape)
% % volshow(extended_mask);

%% 
%% A vectorised approach to the cone-extension method
% Avoids nested loops
cone_height = 2;
max_cone_radius = 1;
[roi_x, roi_y, roi_z] = ind2sub(size(mask), find(mask));
roi_centroid = mean([roi_x, roi_y, roi_z], 1);

[X, Y, Z] = ndgrid(-max_cone_radius:max_cone_radius);
sphere_mask = (X.^2 + Y.^2 + Z.^2) <= max_cone_radius^2;
sphere_indices = find(sphere_mask);
[sx, sy, sz] = ind2sub(size(sphere_mask), sphere_indices);
sphere_offsets = [sx - (max_cone_radius + 1), sy - (max_cone_radius + 1), sz - (max_cone_radius + 1)];

extended_mask = mask;
for i = 1:length(roi_x)
    radial_vector = ([roi_x(i), roi_y(i), roi_z(i)] - roi_centroid) / norm([roi_x(i), roi_y(i), roi_z(i)] - roi_centroid);
    for step = 1:cone_height
        new_pos = round([roi_x(i), roi_y(i), roi_z(i)] + step * radial_vector);
        
        valid_offsets = bsxfun(@plus, sphere_offsets, new_pos);
        valid_indices = sub2ind(size(mask), valid_offsets(:,1), valid_offsets(:,2), valid_offsets(:,3));
        valid_indices = valid_indices(valid_offsets(:,1) > 0 & valid_offsets(:,1) <= size(mask,1) & ...
                                      valid_offsets(:,2) > 0 & valid_offsets(:,2) <= size(mask,2) & ...
                                      valid_offsets(:,3) > 0 & valid_offsets(:,3) <= size(mask,3));
        extended_mask(valid_indices) = 1;

    end
end
volshow(extended_mask);



%% 
% Assuming `boundary_mask` is the binary mask for the boundary voxels
% and `brain_mask` is the original brain binary mask.

% Step 1: Compute the gradient of the brain mask

[grad_x, grad_y, grad_z] = gradient((brain_mask));

% Step 2: Extract gradients at the boundary voxels
[x, y, z] = ind2sub(size(boundary_mask), find(boundary_mask));
normal_x = grad_x(boundary_mask);
normal_y = grad_y(boundary_mask);
normal_z = grad_z(boundary_mask);

% Step 3: Normalize the normal vectors
norms = sqrt(normal_x.^2 + normal_y.^2 + normal_z.^2);
normal_x = normal_x ./ norms;
normal_y = normal_y ./ norms;
normal_z = normal_z ./ norms;


figure;
quiver3(x, y, z, normal_x, normal_y, normal_z, 0.5, 'r'); % Adjust scaling
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Normal Vectors at Boundary Voxels (Using Gradient)');
axis equal;
grid on;
view(3); % 3D view
camlight; lighting gouraud; % Improve lighting

%% Find the components related to the mask
[grad_x, grad_y, grad_z] = gradient(double(brain_mask));

[x, y, z] = ind2sub(size(brain_mask), find(extended_mask));
normal_x = grad_x(extended_mask == 1);
normal_y = grad_y(extended_mask ==1);
normal_z = grad_z(extended_mask ==1 );


norms = sqrt(normal_x.^2 + normal_y.^2 + normal_z.^2);
normal_x = normal_x ./ norms;
normal_y = normal_y ./ norms;
normal_z = normal_z ./ norms;

figure;
quiver3(x, y, z, normal_x, normal_y, normal_z, 2, 'r'); % Adjust scaling
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Normal Vectors at Boundary Voxels (Using Gradient)');
axis equal;
grid on;
view(3); % 3D view
camlight; lighting gouraud; % Improve lighting


%% 
filtered_x = normal_x(~isnan(normal_x));
filtered_y = normal_y(~isnan(normal_y));
filtered_z = normal_z(~isnan(normal_z));

sum_norm_x = sum(filtered_x);
sum_norm_y = sum(filtered_y);
sum_norm_z = sum(filtered_z);

filt_norms = sqrt(sum_norm_x.^2 + sum_norm_y.^2 + sum_norm_z.^2);

av_norm_x = sum_norm_x ./ filt_norms;
av_norm_y = sum_norm_y ./ filt_norms;
av_norm_z = sum_norm_z ./ filt_norms;

real_coord = hdr.sform * [av_norm_x,av_norm_y,av_norm_z,1]'
%disp([av_norm_x,av_norm_y,av_norm_z]);
disp(real_coord);

real_vec = [real_coord(1,:),real_coord(2,:),real_coord(3,:)];
norm = vecnorm(real_vec)
real_vec = real_vec ./ norm;
disp(real_vec);