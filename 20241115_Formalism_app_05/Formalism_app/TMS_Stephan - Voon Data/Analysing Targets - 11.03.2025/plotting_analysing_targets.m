%% Following S.Goetz meetin on 10.03.2025 - Plot & Analyse the FC peak data to compare

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

%%

% Get the size of the volume
[dimX, dimY, dimZ] = size(hdr.vol);
    
[x, y, z] = ndgrid(1:dimX, 1:dimY, 1:dimZ); % Grid of voxel indices

% Flatten the voxel data and coordinates 
voxelData = hdr.vol(:);
x = x(:);
y = y(:);
z = z(:);

% Apply (hdr.sform) to get real-world coordinates
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

%%
figure;

% Plot non-zero voxel data (color by voxel value)
scatter3(realXNonZero, realYNonZero, realZNonZero, 5, voxelDataNonZero, 'filled'); 
hold on;

% Plot zero voxel data with transparency (alpha = 0)
% scatter3(realXZero, realYZero, realZZero, 10, 'r', 'filled', 'MarkerFaceAlpha', 0);  % Transparent red for zeros

xlabel('X');
ylabel('Y');
zlabel('Z');
% title('3D Plot of Voxel Values (Zero and Non-Zero)');
title('Amygdala DLPFC target (mm)');
title('FOX DLPFC target (mm)');
title('VMPFC DLPFC target (mm)');
title('Combined (111) DLPFC target (mm)');
% title('3D Plot of Voxel Values using Real Coordinates (mm)');
colorbar;  % Show color bar for non-zero voxel values
axis equal;
grid on;
%%
figure;

% Plot non-zero voxel data (color by voxel value)
scatter3(realXNonZero, realYNonZero, realZNonZero, 5, voxelDataNonZero, 'filled'); 
hold on;

% Plot zero voxel data with transparency (alpha = 0)
% scatter3(realXZero, realYZero, realZZero, 10, 'r', 'filled', 'MarkerFaceAlpha', 0);  % Transparent red for zeros

xlabel('X');
ylabel('Y');
zlabel('Z');
% title('3D Plot of Voxel Values (Zero and Non-Zero)');
title('3D Plot of Voxel Values using Real Coordinates (mm)');

%Adjust the colour range
minValue = min(voxelDataNonZero);
maxValue = max(voxelDataNonZero);


caxis([minValue, maxValue * 0.15]);  % Control the data range mapped to colours
colormap(jet);  % Change to jet colormap for a more vibrant color scheme

colorbar;  % Show color bar for non-zero voxel values
axis equal;
grid on;

%% Extract the View Matrix
% ax = gca;  % Get current axes handle
% 
% % Extract rotation matrix from View (3x3)
% R = ax.View;  
% 
% % Extract camera parameters
% camPos = ax.CameraPosition(:);   % Convert to column vector (3x1)
% camTgt = ax.CameraTarget(:);     % Convert to column vector (3x1)
% 
% % Compute view direction (negative Z-axis in camera space)
% zAxis = (camPos - camTgt);
% zAxis = zAxis / norm(zAxis);   % Normalize
% 
% % Compute right (X) axis (perpendicular to up and z-axis)
% xAxis = cross(ax.CameraUpVector, zAxis);
% xAxis = xAxis / norm(xAxis);   % Normalize
% 
% % Compute corrected up (Y) axis
% yAxis = cross(zAxis, xAxis);
% 
% % Construct the 3x3 rotation matrix (ensuring proper orientation)
% R_corrected = [xAxis(:), yAxis(:), zAxis(:)];
% 
% % Compute the translation vector
% T = -R_corrected' * camPos;
% 
% % Construct the full 4×4 view transformation matrix
% fullViewMatrix = [R_corrected', T; 0 0 0 1];
% 
% disp('Full 4×4 View Matrix:');
% disp(fullViewMatrix);
%   -0.6406   -0.7679         0   10.6047
%     0.3909   -0.3261    0.8607   -0.3735
%    -0.6610    0.5513    0.5091 -584.3728
%          0         0         0    1.0000
%%
camPos = [-3.793153716956582e+02;3.302138370409575e+02;2.977966156117114e+02];
camTgt = [-29;38;28];
camUp = [0,0,1];
%%
figure;

% Plot non-zero voxel data (color by voxel value)
scatter3(realXNonZero, realYNonZero, realZNonZero, 5, voxelDataNonZero, 'filled'); 
hold on;

% Plot zero voxel data with transparency (alpha = 0)
% scatter3(realXZero, realYZero, realZZero, 10, 'r', 'filled', 'MarkerFaceAlpha', 0);  % Transparent red for zeros

xlabel('X');
ylabel('Y');
zlabel('Z');
% title('3D Plot of Voxel Values (Zero and Non-Zero)');
% title('3D Plot of Voxel Values using Real Coordinates (mm)');
title('Amygdala DLPFC target (mm)');
title('FOX DLPFC target (mm)');
title('VMPFC DLPFC target (mm)');
title('Combined (111) DLPFC target (mm)');
colorbar;  % Show color bar for non-zero voxel values
axis equal;
grid on;
% 
caxis([minValue*100, maxValue ]);  % Control the data range mapped to colours
colormap(flipud(jet   )    );

axNew = gca;  % Get current axes handle for the new figure
set(axNew, 'CameraPosition', camPos);
set(axNew, 'CameraTarget', camTgt);
set(axNew, 'CameraUpVector', camUp);

















%% Access the image data, and convert to real coordinate system
vol_data = hdr.vol;

% Voxel indices
[i, j, k] = ndgrid(0:hdr.dim(2)-1, 0:hdr.dim(3)-1, 0:hdr.dim(4)-1);

% Apply the affine transformation
xx = hdr.sform(1,1) * i + hdr.sform(1,2) * j + hdr.sform(1,3) * k + hdr.sform(1,4);  %%first row of sform
yy = hdr.sform(2,1) * i + hdr.sform(2,2) * j + hdr.sform(2,3) * k + hdr.sform(2,4);  %%second row of sform
zz = hdr.sform(3,1) * i + hdr.sform(3,2) * j + hdr.sform(3,3) * k + hdr.sform(3,4);  %%third row of sform
% 
% xx = hdr.qform(1,1) * i + hdr.qform(1,2) * j + hdr.qform(1,3) * k + hdr.qform(1,4);  %%first row of sform
% yy = hdr.qform(2,1) * i + hdr.qform(2,2) * j + hdr.qform(2,3) * k + hdr.qform(2,4);  %%second row of sform
% zz = hdr.qform(3,1) * i + hdr.qform(3,2) * j + hdr.qform(3,3) * k + hdr.qform(3,4);  %%third row of sform


%% frequency plot of values in hdr.vol

% Ensure hdr.vol is properly loaded
if isfield(hdr, 'vol') && ~isempty(hdr.vol)
    % Flatten the volume data into a 1D array
    voxelData = hdr.vol(:);
    voxelDataNoZero = voxelData(voxelData ~= 0);

    % Get the unique voxel values and their corresponding frequencies
    %[uniqueValues, ~, indices] = unique(voxelData);
    [uniqueValues, ~, indices] = unique(voxelDataNoZero);

    % Calculate the frequency of each unique value
    %frequencies = histcounts(voxelData, [uniqueValues; uniqueValues(end) + 1]);
    frequencies = histcounts(voxelDataNoZero, [uniqueValues; uniqueValues(end) + 1]);

    % Create the frequency plot
    figure;
    bar(uniqueValues, frequencies, 'FaceColor', [0.7, 0.7, 0.7]);

    % Customize the plot
    xlabel('FC peaks per Voxel');
    ylabel('Number of Voxels');
    title('Amygdala FC Peak Distribution');
    title('FOX FC Peak Distribution');
    title('VMPFC FC Peak Distribution');
    title('Combined (111) FC Peak Distribution');
    grid on;

    % Optionally, display how many times a specific value appears, e.g., for value '4'
    valueToCheck = 4;
    countOfValue = sum(voxelData == valueToCheck);
    disp(['Value ', num2str(valueToCheck), ' appears ', num2str(countOfValue), ' times in hdr.vol.']);
else
    error('hdr.vol is either missing or empty. Please check the hdr structure.');
end

disp(sum(voxelDataNoZero));

%% Plot surface where hdr.vol != 0 , using voxel indices as coordinates

% Plot isosurface where hdr.vol equals 1
figure;
p = patch(isosurface(hdr.vol, 0.1));  % Threshold is set just below 1

% isosphere function extracts 3D surface from volumetric data; the function
% finds a f=surface where the value is 1
% isosurface returns a vertices of surface where value is 1, and faces
% connecting these vertices.

%patch creates a 3D surface plot
    
set(p, 'FaceColor', 'red', 'EdgeColor', 'none');  % Remove triangle edges
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Binary Field Visualization (Isosurface)');
axis equal;
grid on;
view(3);
camlight; lighting gouraud;  % Improve lighting for better 3D effect


%% 3D plot using the voxel indices



[dimX, dimY, dimZ] = size(hdr.vol);

% Create a grid of voxel coordinates
[x, y, z] = ndgrid(1:dimX, 1:dimY, 1:dimZ);

% Flatten the voxel data and coordinates into 1D arrays
voxelData = hdr.vol(:);
x = x(:);
y = y(:);
z = z(:);

% Separate the zero and non-zero voxel data
nonZeroIdx = voxelData ~= 0;  % Indices of non-zero voxels
zeroIdx = voxelData == 0;     % Indices of zero voxels

% Non-zero voxel data and coordinates
voxelDataNonZero = voxelData(nonZeroIdx);
xNonZero = x(nonZeroIdx);
yNonZero = y(nonZeroIdx);
zNonZero = z(nonZeroIdx);

% Zero voxel data and coordinates
voxelDataZero = voxelData(zeroIdx);
xZero = x(zeroIdx);
yZero = y(zeroIdx);
zZero = z(zeroIdx);

% Create a 3D scatter plot
figure;

% Plot non-zero voxel data (color by voxel value)
scatter3(xNonZero, yNonZero, zNonZero, 5, voxelDataNonZero, 'filled'); 
hold on;  % Keep the plot open to add zero voxels

% Plot zero voxel data with transparency (alpha = 0)
scatter3(xZero, yZero, zZero, 10, 'r', 'filled', 'MarkerFaceAlpha', 0);  % Transparent red for zeros

% Customize the plot
xlabel('X');
ylabel('Y');
zlabel('Z');
title('3D Plot of Voxel Values (Zero and Non-Zero)');
colorbar;  % Show color bar for non-zero voxel values
axis equal;
grid on;



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
    % voxelDataNonZero = voxelData(nonZeroIdx);
    % realXNonZero = realX(nonZeroIdx);
    % realYNonZero = realY(nonZeroIdx);
    % realZNonZero = realZ(nonZeroIdx);

    realXNonZero = xx(nonZeroIdx);
    realYNonZero = yy(nonZeroIdx);
    realZNonZero = zz(nonZeroIdx);
    % Zero voxel data and coordinates in real-world space
    voxelDataZero = voxelData(zeroIdx);
    realXZero = xx(zeroIdx);
    realYZero = yy(zeroIdx);
    realZZero = zz(zeroIdx);

    % Create a 3D scatter plot
    figure;

    % Plot non-zero voxel data (color by voxel value)
    scatter3(realXNonZero, realYNonZero, realZNonZero, 5, voxelDataNonZero, 'filled'); 
    hold on;  % Keep the plot open to add zero voxels

    % Plot zero voxel data with transparency (alpha = 0)
    % scatter3(realXZero, realYZero, realZZero, 10, 'r', 'filled', 'MarkerFaceAlpha', 0);  % Transparent red for zeros

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

%%



%%
