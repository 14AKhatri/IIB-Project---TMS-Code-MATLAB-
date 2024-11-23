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
%% Load NIfTI file for the Brain model data
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
brainmodel(~mask) = 0;
volData = brainmodel;

% Convert matrix in to grayscale image & Normalize intensity values to range [0, 1]
% volData = mat2gray(volData); 
% volData(volData < 10) = 0; % Set low-intensity values to transparent

% test = volData(:);  
% disp(max(test));
volshow(volData);
%% 
[dimX, dimY, dimZ] = size(volData);
[x, y, z] = ndgrid(1:dimX, 1:dimY, 1:dimZ); %grid of indices

voxelIndices = [x(:), y(:), z(:), ones(numel(x), 1)]';  % Homogeneous coordinates
realCoords = hdr.sform * voxelIndices;  % The real coordinates are stored in columns
% Reshape to original 3D array format
realX = reshape(realCoords(1, :), dimX, dimY, dimZ); realY = reshape(realCoords(2, :), dimX, dimY, dimZ); realZ = reshape(realCoords(3, :), dimX, dimY, dimZ);

%% 
% finding the gradient for the 3D space -> make into unit vectors
[grad_x, grad_y, grad_z] = gradient((volData));
grad_magnitude = sqrt(grad_x.^2 + grad_y.^2 + grad_z.^2);
uvec_x = grad_x ./ grad_magnitude; uvec_y = grad_y ./ grad_magnitude; uvec_z = grad_z ./ grad_magnitude;

% Replace NaNs (from dividing by zero magnitude)
uvec_x(isnan(uvec_x)) = 0; uvec_y(isnan(uvec_y)) = 0; uvec_z(isnan(uvec_y)) = 0;


% % Display filtered gradients and their corresponding real-world coordinates
% figure;
% quiver3(realX, realY, realZ, grad_x, grad_y, grad_z);
% title('Gradients of Non-Zero Voxel Data');
% xlabel('Real X');
% ylabel('Real Y');
% zlabel('Real Z');
% grid on;

%% Using Gradient Coherence for Sulci and Gyri Detection
%Compute the six unique components of the structure tensor and smooth them (e.g., using a Gaussian filter)
T11 = imgaussfilt3(grad_x.^2, 1);  % <g_x^2>
T22 = imgaussfilt3(grad_y.^2, 1);  % <g_y^2>
T33 = imgaussfilt3(grad_z.^2, 1);  % <g_z^2>
T12 = imgaussfilt3(grad_x .* grad_y, 1);  % <g_x * g_y>
T13 = imgaussfilt3(grad_x .* grad_z, 1);  % <g_x * g_z>
T23 = imgaussfilt3(grad_y .* grad_z, 1);  % <g_y * g_z>

% Preallocate coherence map
coherence = zeros(dimX, dimY, dimZ);

% Loop through each voxel to compute coherence
for x = 2:dimX-1
    for y = 2:dimY-1
        for z = 2:dimZ-1
            % Extract structure tensor components at this voxel
            T = [
                T11(x, y, z), T12(x, y, z), T13(x, y, z);
                T12(x, y, z), T22(x, y, z), T23(x, y, z);
                T13(x, y, z), T23(x, y, z), T33(x, y, z)
            ];

            % Compute eigenvalues of the structure tensor
            eigenValues = eig(T);

            % Sort eigenvalues (ensure lambda1 >= lambda2 >= lambda3)
            eigenValues = sort(eigenValues, 'descend');
            lambda1 = eigenValues(1);
            lambda2 = eigenValues(2);
            lambda3 = eigenValues(3);

            % Compute coherence measure
            coherence(x, y, z) = (lambda1 - lambda2) / (lambda1 + lambda2 + lambda3 + eps);
        end
    end
end

volshow(coherence, 'Colormap', parula(256));

%% 
% Assuming you've already computed 'coherence' and 'gradient_magnitude'

% Define thresholds for classification
threshold_coherence = 0.90;  % Adjust based on your data
threshold_magnitude = 0.5;  % Adjust based on your data

%Define the gyri and sulci
gyri_mask = (coherence > threshold_coherence);% & (gradient_magnitude > threshold_magnitude);
sulci_mask = (coherence < threshold_coherence);% & (gradient_magnitude < threshold_magnitude);
sulci_mask(~mask) = 0;

% Visualize the results
figure;
subplot(1,3,1);
imshow(volData(:,:,round(size(volData,3)/2)), []); % Original brain slice
title('Original Brain Slice');

subplot(1,3,2);
imshow(gyri_mask(:,:,round(size(gyri_mask,3)/2)), []); % Gyri Mask
title('Gyri (Ridges)');

subplot(1,3,3);
imshow(sulci_mask(:,:,round(size(sulci_mask,3)/2)), []); % Sulci Mask
title('Sulci (Valleys)');
%% volshow plots of the gyri/sulci
% volshow(coherence, 'Colormap', parula(256));  

% Visualize the Gyri Mask
volshow(gyri_mask, 'Colormap', 'white');  

% Visualize the Sulci Mask
% volshow(sulci_mask, 'Colormap', 'yellow'); 
%% 
% [gyri_surface, gyri_faces] = isosurface(gyri_mask);

% p = patch(isosurface(gyri_mask,0.1))
% figure;
% set(p, 'FaceColor', 'red', 'EdgeColor', 'none');  % Remove triangle edges
% xlabel('X');
% ylabel('Y');
% zlabel('Z');
% title('Binary Field Visualization (Isosurface)');
% axis equal;
% grid on;
% view(3);
% camlight; lighting gouraud;  % Improve lighting for better 3D effect

%% 
[gyri_faceindex, gyri_vertices] = isosurface(gyri_mask); 
%gyri_faceindex (M X 3) returns indices of vertices of each face
%gyri_vertices (N X 3) real coordinates of each vertex  

%% 

function normals = computef_normals(faces, vertices)
% Get vertex coordinates for each face
v1 = vertices(faces(:,1), :);
v2 = vertices(faces(:,2), :);
v3 = vertices(faces(:,3), :);
normals = cross(v2 - v1, v3 - v1, 2);
magnitudes = sqrt(sum(normals.^2, 2));
normals = normals ./ magnitudes;
end

normals = computef_normals(gyri_faceindex, gyri_vertices);
%% 
% Calculate face centroids
face_centroids = (gyri_vertices(gyri_faceindex(:,1), :) + ...
                  gyri_vertices(gyri_faceindex(:,2), :) + ...
                  gyri_vertices(gyri_faceindex(:,3), :)) / 3;

% Plot the isosurface
figure;
% k = patch('Faces', gyri_faceindex, 'Vertices', gyri_vertices);
% set(k, 'FaceColor', 'red', 'EdgeColor', 'none', 'FaceAlpha', 0.1);  % Adjust transparency here
% xlabel('X');
% ylabel('Y');
% zlabel('Z');
% title('Binary Field Visualization (Isosurface)');
% axis equal;
% grid on;
% view(3); %3D view
% camlight; lighting gouraud;  % Improve lighting for better 3D effect
% % 
% hold on;
% % quiver3(face_centroids(:,1), face_centroids(:,2), face_centroids(:,3), ...
% %          normals(:,1), normals(:,2), normals(:,3), 0.1, 'black'); % Adjust the scale factor (0.1) as needed
% % % hold off;


% Sample a subset of faces (adjust the sampling ratio as needed)
sampling_ratio = 0.05; % Adjust this value to control the sampling density
num_faces = size(face_centroids, 1);
selected_indices = randperm(num_faces, round(num_faces * sampling_ratio));

% Select the corresponding face centroids and normal vectors
sampled_centroids = face_centroids(selected_indices, :);
sampled_normals = normals(selected_indices, :);

% Plot the sampled quiver plot
quiver3(sampled_centroids(:,1), sampled_centroids(:,2), sampled_centroids(:,3), ...
         sampled_normals(:,1), sampled_normals(:,2), sampled_normals(:,3), 0.9, 'black');
% hold off;

%% 


%% 

%% 

%% 

normals = compute_normals(gyri_surface, gyri_faces);


function normals = compute_normals(vertices, faces)
    % Preallocate the normals matrix
    normals = zeros(size(faces, 1), 3);

    % Loop through each face of the surface
    for i = 1:size(faces, 1)
        % Get the vertices for the current face
        v1 = vertices(faces(i, 1), :);
        v2 = vertices(faces(i, 2), :);
        v3 = vertices(faces(i, 3), :);

        % Compute two edge vectors of the triangle
        edge1 = v2 - v1;
        edge2 = v3 - v1;

        % Compute the cross product to get the normal vector
        normals(i, :) = cross(edge1, edge2);
    end
    % Normalize the normal vectors (make them unit vectors)
    normals = normals ./ vecnorm(normals, 2, 2);  % Normalize each normal vector
end





%% 



% Plot the gyri surface
figure;
patch('Faces', gyri_faces, 'Vertices', gyri_surface, 'FaceColor', 'g', 'EdgeColor', 'none');
hold on;

% Add the normal vectors as arrows
quiver3(gyri_surface(:, 1), gyri_surface(:, 2), gyri_surface(:, 3), normals(:, 1), normals(:, 2), normals(:, 3), 0.5, 'r');
axis equal;
title('Gyri Surface with Normal Vectors');






%%

% % Smooth gradients to reduce noise
% sigma = 1; % Gaussian kernel size
% grad_x_smooth = imgaussfilt3(grad_x, sigma);
% grad_y_smooth = imgaussfilt3(grad_y, sigma);
% grad_z_smooth = imgaussfilt3(grad_z, sigma);
% 
% % Compute structure tensor components
% J_xx = grad_x_smooth .* grad_x_smooth;
% J_yy = grad_y_smooth .* grad_y_smooth;
% J_zz = grad_z_smooth .* grad_z_smooth;
% J_xy = grad_x_smooth .* grad_y_smooth;
% J_xz = grad_x_smooth .* grad_z_smooth;
% J_yz = grad_y_smooth .* grad_z_smooth;
% 
% % Smooth tensor components for coherence calculation
% J_xx = imgaussfilt3(J_xx, sigma);
% J_yy = imgaussfilt3(J_yy, sigma);
% J_zz = imgaussfilt3(J_zz, sigma);
% J_xy = imgaussfilt3(J_xy, sigma);
% J_xz = imgaussfilt3(J_xz, sigma);
% J_yz = imgaussfilt3(J_yz, sigma);
% 
% coherence = J_xx + J_yy + J_zz; % Sum of diagonal elements approximates the largest eigenvalue
% 
% % Thresholds for gyri and sulci
% threshold_high = 0.0780*0.7; % Empirical value for sulci
% threshold_low = 0.0780*0.2;  % Empirical value for gyri
% 
% % Generate masks
% sulci_mask = coherence > threshold_high;
% gyri_mask = coherence < threshold_low;
% 
% volshow(volData); % Original volume
% hold on;
% % Visualize gyri and sulci masks with different colors
% volshow(sulci_mask, 'Color', [1, 0, 0], 'AlphaData', 0.5); % Sulci in red
% volshow(gyri_mask, 'Color', [0, 0, 1], 'AlphaData', 0.5); % Gyri in blue


%% 

coherence_values = coherence(:);

% Remove NaN or Inf values (if any)
coherence_values = coherence_values(~isnan(coherence_values) & ~isinf(coherence_values));

% Optional: Remove zeros if they're irrelevant to your analysis
% coherence_values = coherence_values(coherence_values ~= 0);

% Plot Kernel Density Estimate (KDE)
disp(max(coherence_values))



