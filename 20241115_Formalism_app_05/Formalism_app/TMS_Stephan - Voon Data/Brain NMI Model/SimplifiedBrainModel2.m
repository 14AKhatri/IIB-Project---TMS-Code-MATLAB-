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
% volshow(brainmodel);
%% 
brainmodel(mask == 0) = brainmodel(mask == 0) * 0.1; %for visualisation - makes the NOT ROI darker
volData = brainmodel;
volshow(volData);

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
% Get the boundary voxels of the binary mask - finds exact voxel locations
% of boundaries
boundary_mask = bwperim(brain_mask,26); %returns the boundary voxels of brain_mask
% 26 indicates the 3D connectivity
volshow(boundary_mask);

%% 
% Extract external surface using isosurface
[faces, vertices] = isosurface(brain_mask, 0.5);
function normals = compute_normals(faces, vertices)
    v1 = vertices(faces(:, 1), :);
    v2 = vertices(faces(:, 2), :);
    v3 = vertices(faces(:, 3), :);

    % cross product
    normals = cross(v2 - v1, v3 - v1, 2);
    % magnitudes = vecnorm(normals, 2, 2);
    magnitudes = sqrt(sum(normals.^2, 2));
    normals = normals ./ magnitudes;
end

% Compute normals
normals = compute_normals(faces, vertices);
%% 
[faces, vertices] = isosurface(brain_mask);
normals = compute_normals(faces, vertices);
accumulated_normals = zeros(size(brain_mask)); % Accumulate normals in a 3D grid corresponding to vertices
normal_count = zeros(size(brain_mask)); % Count how many normals are assigned to each voxel

i=1;
v1 = faces(i, 1);v2 = faces(i, 2);v3 = faces(i, 3);

vertex1 = round(vertices(v1, :));
vertex2 = round(vertices(v2, :));
vertex3 = round(vertices(v3, :));
vertices_list = {vertex1, vertex2, vertex3};

vertex = vertices_list{1}; %gets the vector
idx = sub2ind(size(brain_mask), vertex(1), vertex(2), vertex(3));

accumulated_normals(vertex, 1) = accumulated_normals(vertex, 1) + normals(i, 1);
accumulated_normals(vertex, 2) = accumulated_normals(vertex, 2) + normals(i, 2);
accumulated_normals(vertex, 3) = accumulated_normals(vertex, 3) + normals(i, 3);

%% 


accumulated_normals = zeros(size(brain_mask)); % Accumulate normals in a 3D grid corresponding to vertices
normal_count = zeros(size(brain_mask)); % Count how many normals are assigned to each voxel

for i = 1:size(faces, 1)
    % Get the current face's vertices indexes
    v1 = faces(i, 1);v2 = faces(i, 2);v3 = faces(i, 3);
    
    % Get the vertex coordinates 
    vertex1 = round(vertices(v1, :));
    vertex2 = round(vertices(v2, :));
    vertex3 = round(vertices(v3, :));

    vertices_list = {vertex1, vertex2, vertex3}; 
    for j = 1:3
        vertex = vertices_list{j};
        
        % Check if the vertex is within the bounds of the mask
        if all(vertex >= 1) && all(vertex <= size(brain_mask))
            % Convert the voxel coordinates to linear indices
            idx = sub2ind(size(brain_mask), vertex(1), vertex(2), vertex(3));
            
            % Accumulate the normal
            accumulated_normals(vertex, 1) = accumulated_normals(vertex, 1) + normals(i, 1);
            accumulated_normals(vertex, 2) = accumulated_normals(vertex, 2) + normals(i, 2);
            accumulated_normals(vertex, 3) = accumulated_normals(vertex, 3) + normals(i, 3);
            
            % Count the number of normals assigned to the voxel
            normal_count(vertex) = normal_count(vertex) + 1;
        end
    end
end

% Normalize the accumulated normals by dividing by the normal count
valid_voxels = normal_count ~= 0; % Only normalize voxels that have received normals
accumulated_normals(valid_voxels, :) = accumulated_normals(valid_voxels, :) ./ ...
                                         sqrt(sum(accumulated_normals(valid_voxels, :).^2, 2));

% Extract the filtered normals and corresponding voxel locations
valid_voxels_linear = find(valid_voxels);  % Find linear indices where normal_count is non-zero
[filtered_voxels_x, filtered_voxels_y, filtered_voxels_z] = ind2sub(size(brain_mask), valid_voxels_linear);

% Check the dimensions of filtered_voxels_x, filtered_voxels_y, filtered_voxels_z
disp(size(filtered_voxels_x));  % Should be the same as the number of valid voxels

% Extract the corresponding filtered normals
filtered_normals = accumulated_normals([filtered_voxels_x, filtered_voxels_y, filtered_voxels_z]);
%% 

% Plot the normals at valid voxel locations
figure;
quiver3(filtered_voxels_x, filtered_voxels_y, filtered_voxels_z, ...
        filtered_normals(:, 1), filtered_normals(:, 2), filtered_normals(:, 3), ...
        0.5, 'Color', 'k', 'LineWidth', 1.5);

% Set labels and title
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Filtered Normals at Voxel Locations on Brain Surface');

% Adjust visualization
axis equal;
grid on;
view(3);  % Set 3D view

%% 
%% 


%% 

% Calculate centroids of triangular faces
face_centroids = (vertices(faces(:, 1), :) + ...
                  vertices(faces(:, 2), :) + ...
                  vertices(faces(:, 3), :)) / 3;

% Convert centroids to voxel indices
voxel_indices = round(face_centroids);

% Ensure indices are within bounds
% voxel_indices = max(min(voxel_indices, size(brain_mask)), 1);

% Initialize a 4D array to store normals (X, Y, Z, normal components)
normal_field = zeros([size(brain_mask), 3]);

% Assign normals to corresponding voxels
for i = 1:size(faces, 1)
    x = voxel_indices(i, 1);
    y = voxel_indices(i, 2);
    z = voxel_indices(i, 3);

    % Store the normal vector at the voxel location
    normal_field(x, y, z, :) = normals(i, :);
end
%% 
% Sample a subset of the surface points
sampling_ratio = 0.5; % Adjust for visualization
num_faces = size(face_centroids, 1);
selected_indices = randperm(num_faces, round(num_faces * sampling_ratio));

sampled_centroids = face_centroids(selected_indices, :);
sampled_normals = normals(selected_indices, :);

% Plot the vectors
figure;
scatter3(sampled_centroids(:, 1), sampled_centroids(:, 2), sampled_centroids(:, 3), 'red');
hold on;
quiver3(sampled_centroids(:, 1), sampled_centroids(:, 2), sampled_centroids(:, 3), ...
        sampled_normals(:, 1), sampled_normals(:, 2), sampled_normals(:, 3), 0.9, 'black');
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Surface Normals');
grid on;
axis equal;

%%
% Flatten the 3D binary mask into a 1D column vector
filtered_normals = normal_field;
filtered_normals(mask ~= 0) = 0;
volshow(filtered_normals);

% % vector normals stored in normals; voxel locations in face_centroids
% % Find indices of valid voxels within the binary mask
% valid_voxel_indices = find(mask);
% 
% % Create a logical array to track valid normals
% valid_normals = false(size(voxel_indices, 1), 1);
% 
% % Check each voxel index against the mask
% for i = 1:size(voxel_indices, 1)
%     voxel_linear_index = sub2ind(size(mask), ...
%                                  voxel_indices(i, 1), ...
%                                  voxel_indices(i, 2), ...
%                                  voxel_indices(i, 3));
%     if ismember(voxel_linear_index, valid_voxel_indices)
%         valid_normals(i) = true;
%     end
% end
% 
% % Filter normals and centroids
% filtered_normals = normals(valid_normals, :);
% filtered_centroids = face_centroids(valid_normals, :);
%% 
% Visualize filtered normals
figure;
scatter3(filtered_centroids(:, 1), filtered_centroids(:, 2), filtered_centroids(:, 3), 'red');
hold on;
quiver3(filtered_centroids(:, 1), filtered_centroids(:, 2), filtered_centroids(:, 3), ...
        filtered_normals(:, 1), filtered_normals(:, 2), filtered_normals(:, 3), 0.9, 'black');
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Filtered Surface Normals (Binary Mask)');
grid on;
axis equal;



%% 
% Assume `brain_mask` or `mask` is the brain region mask

% Get the linear indices of the non-zero elements of the mask
[rows, cols, pages] = ind2sub(size(mask), find(mask ~= 0));

% Extract the normals at those positions
filtered_normals = normal_field(mask ~= 0);

[faces, vertices] = isosurface(brain_mask);

% Create the surface plot
figure;
patch('Faces', faces, 'Vertices', vertices, 'FaceColor', 'red', 'EdgeColor', 'none', 'FaceAlpha', 0.5);
hold on;

% Overlay filtered normals (quiver plot)
quiver3(rows, cols, pages, ...
        filtered_normals(:, 4), filtered_normals(:, 5), filtered_normals(:, 6), ...
        scale_factor, 'Color', 'k', 'LineWidth', 1.5);

% Set axis labels and title
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Filtered Normals with Brain Surface');

axis equal;
grid on;
view(3);  % Set 3D view
camlight; lighting gouraud;  % Improve lighting for better 3D effect

%% 
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
avg_normal = mean(filtered_normals, 1);
% avg_normal = mean(roi_normals, 1);

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
% quiver3(roi_face_centroids(:,1), roi_face_centroids(:,2), roi_face_centroids(:,3), ...
%          roi_normals(:,1), roi_normals(:,2), roi_normals(:,3), 1, 'black'); 
quiver3(filtered_centroids(:,1), filtered_centroids(:,2), filtered_centroids(:,3), ...
          filtered_normals(:,1), filtered_normals(:,2), filtered_normals(:,3), 1, 'black'); 
hold off;
