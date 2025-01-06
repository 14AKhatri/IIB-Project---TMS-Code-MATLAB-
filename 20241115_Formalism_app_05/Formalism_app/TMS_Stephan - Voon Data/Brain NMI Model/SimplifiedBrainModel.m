%This code attempts to identify the vector direction as outlined by S.Goetz
% on 21/11/2024

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
brainmodel(mask == 0) = brainmodel(mask == 0) * 0.1; %for visualisation - makes the NOT ROI darker
volData = brainmodel;
volshow(volData);

%% Generate a Binary Mask
% This is redundant is the DLPFC mask is used.
binary_mask = mask ~= 0; % should be the same as mask if a mask is used 
ROI = binary_mask;
%volshow(ROI);
% masked_brain = brainmodel .* mask; % same as mask
%volshow(masked_brain); 
%% 
% % Define the structuring element (sphere for radial extension)
% % The radius of the sphere determines how far the ROI will extend.
% radius = 6; % Radius in voxels for radial extension
% se = strel('sphere', radius);
% 
% % Perform dilation on the ROI
% extended_ROI = imdilate(ROI > 0, se);
% 
% % Overlay the extended ROI on the brain model
% brain_overlay = brainmodel; 
% brain_overlay(extended_ROI > 0) = max(brain_overlay(:)); % Highlight extended ROI
% 
% % Use volshow to visualize the result
% volshow(brain_overlay);

%% 
% Assumes `mask` is the binary ROI mask (1 for ROI, 0 elsewhere)
% `origin` is the brain origin (e.g., center of the brain volume)
origin = [size(mask, 1)/2, size(mask, 2)/2, size(mask, 3)/2]; % Brain center
[roi_x, roi_y, roi_z] = ind2sub(size(mask), find(mask > 0)); % ROI voxel coordinates

% Calculate the ROI centroid
roi_centroid = mean([roi_x, roi_y, roi_z], 1);

% Parameters for the cone
cone_height = 10; % Height of the cone (in voxels)
max_cone_radius = 5; % Maximum radius of the cone at its farthest extent

% Create a new mask for the extended ROI
extended_mask = mask;

% Loop through each voxel in the ROI
for i = 1:length(roi_x)
    % Current voxel position
    voxel_pos = [roi_x(i), roi_y(i), roi_z(i)];
    
    % Compute radial direction vector from ROI centroid
    radial_vector = voxel_pos - roi_centroid;
    radial_vector = radial_vector / norm(radial_vector); % Normalize to unit vector
    
    % Step outward to extend in a cone shape
    for step = 1:cone_height
        % Compute the current cone radius
        current_radius = (step / cone_height) * max_cone_radius;
        
        % Compute the new voxel position along the radial vector
        new_pos = voxel_pos + step * radial_vector;
        
        % Check all voxels within the current radius of the cone slice
        [x_sphere, y_sphere, z_sphere] = ndgrid(-current_radius:current_radius, ...
                                                -current_radius:current_radius, ...
                                                -current_radius:current_radius);
        sphere_voxels = [x_sphere(:), y_sphere(:), z_sphere(:)];
        sphere_voxels = sphere_voxels(vecnorm(sphere_voxels, 2, 2) <= current_radius, :);
        
        % Translate sphere voxels to new position and update mask
        for j = 1:size(sphere_voxels, 1)
            voxel = round(new_pos + sphere_voxels(j, :));
            if all(voxel > 0) && voxel(1) <= size(mask, 1) && ...
               voxel(2) <= size(mask, 2) && voxel(3) <= size(mask, 3)
                extended_mask(voxel(1), voxel(2), voxel(3)) = 1;
            end
        end
    end
end

% Visualize the extended ROI (cone shape)
volshow(extended_mask);






%% 

% brainmodel = hdr.vol;  % Brain volume data
% % tbrain = brainmodel > 0.1; %makes all the voxels 1
% % tbrain = ROI>0;
% tbrain = (brainmodel > 0.1) | (ROI > 0);
% 
% volshow(tbrain);
volshow(brainmodel);
[faces, vertices] = isosurface(brainmodel);

%% 
% Plot the isosurface
figure;
k = patch('Faces', faces, 'Vertices', vertices);
set(k, 'FaceColor', 'red', 'EdgeColor', 'none', 'FaceAlpha', 0.5);  % Adjust transparency here
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Binary Field Visualization (Isosurface)');
axis equal;
grid on;
view(3); %3D view
camlight; lighting gouraud;  % Improve lighting for better 3D effect
% 

%% 

function normals = computef_normals(faces, vertices)
% Get vertex coordinates for each face
v1 = vertices(faces(:,1), :);
v2 = vertices(faces(:,2), :);
v3 = vertices(faces(:,3), :);
normals = cross(v2 - v1, v3 - v1, 2);
magnitudes = sqrt(sum(normals.^2, 2));
normals = normals ./ magnitudes; %returns unitvector


end

normals = computef_normals(faces, vertices); % Normals will have the same size as `vertices`
%% 
face_centroids = (vertices(faces(:,1), :) + ...
                  vertices(faces(:,2), :) + ...
                  vertices(faces(:,3), :)) / 3;

%% Flip the vectors to always be pointing inwards
for i = 1:size(normals, 1)
    vector_to_origin = -face_centroids(i, :);  %vector facing inwards
    
    if dot(normals(i, :), vector_to_origin) < 0 %'0' corresponds to 90
        normals(i, :) = -normals(i, :);
    end
end
%% Sampled quiver plot of vectors

sampling_ratio = 0.05; 
num_faces = size(face_centroids, 1);
selected_indices = randperm(num_faces, round(num_faces * sampling_ratio));

sampled_centroids = face_centroids(selected_indices, :);
sampled_normals = normals(selected_indices, :);

quiver3(sampled_centroids(:,1), sampled_centroids(:,2), sampled_centroids(:,3), ...
         sampled_normals(:,1), sampled_normals(:,2), sampled_normals(:,3), 0.9, 'black');


%% Only want the ones in the ROI
% 
% mask = ROI;  % Mask for ROI
% roi_normals = []; 
% roi_face_centroids = [];
% 
% % Iterating through each face
% for i = 1:size(faces, 1)
%     face_vertices = vertices(faces(i,:), :); %real coordinates of vertices
% 
%     % Convert the vertex coordinates to voxel indices (rounded to nearest voxel)
%     voxel_indices = round(face_vertices);
%     voxel_indices = max(min(voxel_indices, size(mask)), 1);  % Ensure valid indices
%     % Convert the face vertex coordinates (in mm) to voxel indices, accounting for 2mm voxel size
% 
% 
% 
% 
%     % Check if the vertices of the current face lie within the ROI (mask)
%     mask_values = mask(sub2ind(size(mask), voxel_indices(:,1), voxel_indices(:,2), voxel_indices(:,3)));
% 
%     % If all vertices in the face are inside the mask, add the normal vector
%     if all(mask_values == 1)
%         roi_normals = [roi_normals; normals(i, :)];
%         face_centroid = mean(face_vertices, 1);  
%         roi_face_centroids = [roi_face_centroids; face_centroid];
%     end
% end
%% 
mask = ROI;  % Mask for ROI
roi_normals = []; 
roi_face_centroids = [];

% Define the neighborhood range (in voxels)
neighborhood = 3;

% Get the size of the mask for boundary checks
mask_size = size(mask);

% Iterating through each face
for i = 1:size(faces, 1)
    face_vertices = vertices(faces(i, :), :); % Real coordinates of vertices
    
    % Convert the vertex coordinates to voxel indices (rounded to nearest voxel)
    voxel_indices = round(face_vertices);
    voxel_indices = max(min(voxel_indices, mask_size), 1);  % Ensure valid indices

    % Check if any vertex is within two voxels of the ROI
    in_roi = false;  % Initialize flag for proximity to ROI
    for j = 1:size(voxel_indices, 1)
        % Extract the current voxel index
        vx = voxel_indices(j, 1);
        vy = voxel_indices(j, 2);
        vz = voxel_indices(j, 3);
        
        % Define the range for checking neighborhood
        x_range = max(vx - neighborhood, 1):min(vx + neighborhood, mask_size(1));
        y_range = max(vy - neighborhood, 1):min(vy + neighborhood, mask_size(2));
        z_range = max(vz - neighborhood, 1):min(vz + neighborhood, mask_size(3));
        
        % Check if any voxel in the neighborhood is part of the ROI
        if any(mask(x_range, y_range, z_range), 'all')
            in_roi = true;
            break;  % Exit loop early if proximity is confirmed
        end
    end
    
    % If the face is near the ROI, add its normal vector and centroid
    if in_roi
        roi_normals = [roi_normals; normals(i, :)];
        face_centroid = mean(face_vertices, 1);  
        roi_face_centroids = [roi_face_centroids; face_centroid];
    end
end




%% find the average of the vectors in the ROI

avg_normal = mean(roi_normals, 1);

figure;
k = patch('Faces', faces, 'Vertices', vertices);
set(k, 'FaceColor', 'red', 'EdgeColor', 'none', 'FaceAlpha', 0.1);  % Adjust transparency here
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Plot of Vectors normal to ROI');
axis equal;
grid on;
view(3); %3D view
camlight; lighting gouraud;  % Improve lighting for better 3D effect
% 
hold on;
quiver3(roi_face_centroids(:,1), roi_face_centroids(:,2), roi_face_centroids(:,3), ...
         roi_normals(:,1), roi_normals(:,2), roi_normals(:,3), 1, 'black'); 
hold off;
