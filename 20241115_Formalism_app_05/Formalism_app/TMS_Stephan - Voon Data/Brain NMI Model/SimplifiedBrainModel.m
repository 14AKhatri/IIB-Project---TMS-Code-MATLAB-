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
brainmodel(mask == 0) = brainmodel(mask == 0) * 0.1; %for visualisation purposes
volData = brainmodel;
volshow(volData);

%% 
binary_mask = mask ~= 0; %this should be the same as mask if a mask is used
ROI = binary_mask;
volshow(ROI);
masked_brain = brainmodel .* mask;
volshow(masked_brain);
%% 
%% 

% Define the structuring element (sphere for radial extension)
% The radius of the sphere determines how far the ROI will extend.
radius = 5; % Radius in voxels for radial extension
se = strel('sphere', radius);

% Perform dilation on the ROI
extended_ROI = imdilate(ROI > 0, se);

% Overlay the extended ROI on the brain model
brain_overlay = brainmodel; 
brain_overlay(extended_ROI > 0) = max(brain_overlay(:)); % Highlight extended ROI

% Use volshow to visualize the result
volshow(brain_overlay);




%% 
% Input: ROI (binary mask where ROI = 1, background = 0)

% Compute the distance transform from both inside and outside the ROI
distance_out = bwdist(ROI == 0); % Distance of background voxels to ROI boundary
distance_in = bwdist(ROI == 1);  % Distance of ROI voxels to background

% Specify the radial distance for inward and outward extension (in voxels)
radial_distance = 10; % Adjust based on desired distance

% Create the extended ROI
extended_ROI = (distance_out <= radial_distance) | (distance_in <= radial_distance);

% Visualization
volshow(extended_ROI);




%% 

brainmodel = hdr.vol;  % Brain volume data
% tbrain = brainmodel > 0.1; %makes all the voxels 1
% tbrain = ROI>0;
tbrain = (brainmodel > 0.1) | (ROI > 0);

volshow(tbrain);

[faces, vertices] = isosurface(tbrain);

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

mask = extended_ROI;  % Mask for ROI
roi_normals = []; 
roi_face_centroids = [];

% Iterating through each face
for i = 1:size(faces, 1)
    face_vertices = vertices(faces(i,:), :); %real coordinates of vertices
    
    % Convert the vertex coordinates to voxel indices (rounded to nearest voxel)
    voxel_indices = round(face_vertices);
    voxel_indices = max(min(voxel_indices, size(mask)), 1);  % Ensure valid indices
    
    % Check if the vertices of the current face lie within the ROI (mask)
    mask_values = mask(sub2ind(size(mask), voxel_indices(:,1), voxel_indices(:,2), voxel_indices(:,3)));
    
    % If all vertices in the face are inside the mask, add the normal vector
    if all(mask_values == 1)
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
title('Binary Field Visualization (Isosurface)');
axis equal;
grid on;
view(3); %3D view
camlight; lighting gouraud;  % Improve lighting for better 3D effect
% 
hold on;
quiver3(roi_face_centroids(:,1), roi_face_centroids(:,2), roi_face_centroids(:,3), ...
         roi_normals(:,1), roi_normals(:,2), roi_normals(:,3), 1, 'black'); 
hold off;
