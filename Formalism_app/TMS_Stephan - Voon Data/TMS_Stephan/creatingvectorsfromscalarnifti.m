% Load NIfTI file
[fileName, filePath] = uigetfile('*.nii', 'Select a NIfTI file');
if fileName ~= 0
    fullFilePath = fullfile(filePath, fileName);
    disp(['Selected file: ', fullFilePath]);

    % Load the file using the path
    hdr = nifti_load(fullFilePath);
else
    disp('File selection canceled.');
end

%% Plotting using the real coordinates
% Ensure hdr.vol is properly loaded
if isfield(hdr, 'vol') && ~isempty(hdr.vol)
    % Get the size of the volume
    [dimX, dimY, dimZ] = size(hdr.vol);
    
    % Create a grid of voxel indices
    [x, y, z] = ndgrid(1:dimX, 1:dimY, 1:dimZ);

    % Flatten the voxel data and coordinates into 1D arrays
    voxelData = hdr.vol(:);
    x = x(:);
    y = y(:);
    z = z(:);

    % Apply the affine transformation (hdr.sform) to get real-world coordinates
    % Voxel coordinates are in [x, y, z, 1] homogeneous coordinates.
    % The sform matrix is a 4x4 matrix
    realCoords = hdr.sform * [x, y, z, ones(length(x), 1)]';  % Apply the sform transformation
    realX = realCoords(1, :);  % Extract real-world x coordinates
    realY = realCoords(2, :);  % Extract real-world y coordinates
    realZ = realCoords(3, :);  % Extract real-world z coordinates

    % Separate the zero and non-zero voxel data
    nonZeroIdx = voxelData ~= 0;  % Indices of non-zero voxels
    zeroIdx = voxelData == 0;     % Indices of zero voxels

    % Non-zero voxel data and coordinates in real-world space
    voxelDataNonZero = voxelData(nonZeroIdx);
    realXNonZero = realX(nonZeroIdx);
    realYNonZero = realY(nonZeroIdx);
    realZNonZero = realZ(nonZeroIdx);

    % Zero voxel data and coordinates in real-world space
    voxelDataZero = voxelData(zeroIdx);
    realXZero = realX(zeroIdx);
    realYZero = realY(zeroIdx);
    realZZero = realZ(zeroIdx);

    % Create a 3D scatter plot
    figure;

    % Plot non-zero voxel data (color by voxel value)
    scatter3(realXNonZero, realYNonZero, realZNonZero, 5, voxelDataNonZero, 'filled'); 
    hold on;  % Keep the plot open to add zero voxels

    % Plot zero voxel data with transparency (alpha = 0)
    % scatter3(realXZero, realYZero, realZZero, 10, 'r', 'filled', 'MarkerFaceAlpha', 0);  % Transparent red for zeros

    % Customize the plot
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    title('3D Plot of Voxel Values (Zero and Non-Zero)');
    colorbar;  % Show color bar for non-zero voxel values
    axis equal;
    grid on;

else
    error('hdr.vol is either missing or empty. Please check the hdr structure.');
end

%%
%% Plot (quiver) the gradient values using real coordinates
if isfield(hdr, 'vol') && ~isempty(hdr.vol)

    [dimX, dimY, dimZ] = size(hdr.vol);
    [x, y, z] = ndgrid(1:dimX, 1:dimY, 1:dimZ);
    voxelData = hdr.vol(:);
    x = x(:);
    y = y(:);
    z = z(:);

    realCoords = hdr.sform * [x, y, z, ones(length(x), 1)]';  % Apply the sform transformation
    realX = realCoords(1, :);  % Extract real-world x coordinates
    realY = realCoords(2, :);  % Extract real-world y coordinates
    realZ = realCoords(3, :);  % Extract real-world z coordinates

    realX = realX(:);  % Extract real-world x coordinates
    realY = realY(:);  % Extract real-world y coordinates
    realZ = realZ(:);  % Extract real-world z coordinates

    % Compute gradient of the voxel data (scalar field)
    [gradX, gradY, gradZ] = gradient(hdr.vol);  % Call the function to compute gradien
    gradX = gradX(:);
    gradY = gradY(:);
    gradZ = gradZ(:);

     % Normalize the gradient vectors
    % magnitude = sqrt(gradX.^2 + gradY.^2 + gradZ.^2);
    % gradX = gradX ./ magnitude;
    % gradY = gradY ./ magnitude;
    % gradZ = gradZ ./ magnitude;

    % Only plot non-zero voxel data
    nonZeroIdx = voxelData ~= 0;
    realXNonZero = realX(nonZeroIdx);
    realYNonZero = realY(nonZeroIdx);
    realZNonZero = realZ(nonZeroIdx);
    voxelDataNonZero = voxelData(nonZeroIdx);

    % % Create a 3D scatter plot with reduced marker size
    figure;
    % scatter3(realXNonZero, realYNonZero, realZNonZero, 1, voxelDataNonZero, 'filled');
 

    % % Add gradient vectors as TMS targets
    % quiver3(realX, realY, realZ, gradX, gradY, gradZ, 0.5, 'Color', 'k', 'LineWidth', 1);
    
    step = 5;
    realXStep = realX(1:step:end);
    realYStep = realY(1:step:end);
    realZStep = realZ(1:step:end);
    gradXStep = gradX(1:step:end);
    gradYStep = gradY(1:step:end);
    gradZStep = gradZ(1:step:end);
    % Add the reduced set of gradient vectors (every 5th point)
    quiver3(realXStep, realYStep, realZStep, gradXStep, gradYStep, gradZStep, 4, 'Color', 'k', 'LineWidth', 1);
    
    % Customize the plot
    colormap(parula);
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    title('3D Plot of Voxel Values and TMS Target Vectors');
    colorbar;  % Show color bar for voxel values
    axis equal;
    grid on;
else
    error('hdr.vol is either missing or empty. Please check the hdr structure.');
end
%% 
% Plot gradient magnitudes against voxel positions
gradMagnitude = sqrt(gradX(:).^2 + gradY(:).^2 + gradZ(:).^2);

% Filter out zero gradient magnitudes
nonZeroGrad = gradMagnitude(gradMagnitude > 0);

% Plot a histogram of the non-zero gradient magnitudes
figure;
histogram(nonZeroGrad, 50); % Adjust the number of bins as needed
xlabel('Gradient Magnitude');
ylabel('Frequency');
title('Distribution of Non-Zero Gradient Magnitudes');
grid on;


%%
% Ensure hdr.vol is loaded
if isfield(hdr, 'vol') && ~isempty(hdr.vol)

    % Compute gradients
    [gradX, gradY, gradZ] = gradient(hdr.vol);

    % Get volume dimensions
    [dimX, dimY, dimZ] = size(hdr.vol);

    % Initialize a correlation map
    correlationMap = zeros(dimX, dimY, dimZ);

    % Loop through voxels (avoiding edges for simplicity)
    for x = 2:dimX-1
        for y = 2:dimY-1
            for z = 2:dimZ-1
                % Current voxel gradient
                currentGradient = [gradX(x, y, z), gradY(x, y, z), gradZ(x, y, z)];

                % Adjacent voxel gradients
                neighbors = [
                    gradX(x+1, y, z), gradY(x+1, y, z), gradZ(x+1, y, z);
                    gradX(x-1, y, z), gradY(x-1, y, z), gradZ(x-1, y, z);
                    gradX(x, y+1, z), gradY(x, y+1, z), gradZ(x, y+1, z);
                    gradX(x, y-1, z), gradY(x, y-1, z), gradZ(x, y-1, z);
                    gradX(x, y, z+1), gradY(x, y, z+1), gradZ(x, y, z+1);
                    gradX(x, y, z-1), gradY(x, y, z-1), gradZ(x, y, z-1);
                ];

                % Compute correlation with each neighbor
                correlations = zeros(size(neighbors, 1), 1);
                for n = 1:size(neighbors, 1)
                    neighborGradient = neighbors(n, :);

                    % Avoid division by zero
                    if norm(currentGradient) > 0 && norm(neighborGradient) > 0
                        correlations(n) = dot(currentGradient, neighborGradient) / ...
                                          (norm(currentGradient) * norm(neighborGradient));
                    else
                        correlations(n) = 0;  % Assign zero correlation if magnitude is zero
                    end
                end

                % Store average correlation
                correlationMap(x, y, z) = mean(correlations);
            end
        end
    end

    % Visualization
    % figure;
    % slice = round(dimZ / 2);  % Select a middle slice for visualization
    % imagesc(correlationMap(:, :, slice));  % Display correlation for the middle slice
    % colorbar;
    % title('Correlation of Gradient Vectors with Adjacent Voxels');
    % xlabel('X');
    % ylabel('Y');

    % 3D Visualization using volshow (requires Image Processing Add-On)
    %volshow(correlationMap, 'Colormap', jet(256), 'Alphamap', linspace(0, 1, 256));
    volshow(correlationMap, 'Colormap', hot(256),'Alphamap', linspace(0, 1, 256));  % Example for "hot" colormap

else
    error('hdr.vol is either missing or empty. Please check the hdr structure.');
end

%% Gradient Coherence

%Compute the six unique components of the structure tensor and smooth them (e.g., using a Gaussian filter)
T11 = imgaussfilt3(gradX.^2, 1);  % <g_x^2>
T22 = imgaussfilt3(gradY.^2, 1);  % <g_y^2>
T33 = imgaussfilt3(gradZ.^2, 1);  % <g_z^2>
T12 = imgaussfilt3(gradX .* gradY, 1);  % <g_x * g_y>
T13 = imgaussfilt3(gradX .* gradZ, 1);  % <g_x * g_z>
T23 = imgaussfilt3(gradY .* gradZ, 1);  % <g_y * g_z>

% Preallocate coherence map
[dimX, dimY, dimZ] = size(hdr.vol);
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

%test


