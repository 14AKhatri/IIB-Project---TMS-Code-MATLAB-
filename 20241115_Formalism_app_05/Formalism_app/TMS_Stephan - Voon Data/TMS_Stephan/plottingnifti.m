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
% Access the image data, and convert to real coordinate system
vol_data = hdr.vol;

voxel_sizes = hdr.pixdim; %displays the voxel dimensions in the x,y,z directions
num2str(voxel_sizes(2:4))
disp(['Data type: ', num2str(hdr.datatype)]); %%data type 2 corresponds to uint8


%matrix, shifting coordinate system in the z-direction - i think this is a
%special case
% [xx,yy,zz]=meshgrid(linspace(-((hdr.dim(2)-1)*hdr.sform(1,1)+hdr.sform(1,4)),(hdr.dim(2)-1)*hdr.sform(1,1)+hdr.sform(1,4),hdr.dim(2)),...
%                                 linspace(-((hdr.dim(3)-1)*hdr.sform(2,2)+hdr.sform(2,4)),(hdr.dim(3)-1)*hdr.sform(2,2)+hdr.sform(2,4),hdr.dim(3)),...
%                                 linspace((hdr.dim(4)-1)*hdr.sform(3,3)+hdr.sform(3,4),hdr.sform(3,4),hdr.dim(4)));


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
%%
% Define which slices you want to plot
% These can be changed based on the size of your volume
slice_idz = round(linspace(1, size(vol_data, 3), 6));  % indices for slices in axial (Z direction)
slice_idy = round(linspace(1, size(vol_data, 2), 5));  % Coronal slices (Y direction)
slice_idx = round(linspace(1, size(vol_data, 1), 5));  % Sagittal slices (X direction)

slice_idz = [54,55,56,57,58,59,60]
% Create a figure to hold the slices
figure;

% Plot axial slices (Z direction)
subplot(2, 3, 1);
imshow(squeeze(vol_data(:, :, slice_idz(1))), [min(vol_data(:)), max(vol_data(:))]);
title('Axial Slice 1');

subplot(2, 3, 2);
imshow(squeeze(vol_data(:, :, slice_idz(2))), [min(vol_data(:)), max(vol_data(:))]);
title('Axial Slice 2');

subplot(2, 3, 3);
imshow(squeeze(vol_data(:, :, slice_idz(3))), [min(vol_data(:)), max(vol_data(:))]);
title('Axial Slice 3');

subplot(2, 3, 4);
imshow(squeeze(vol_data(:, :, slice_idz(4))), [min(vol_data(:)), max(vol_data(:))]);
title('Axial Slice 4');

subplot(2, 3, 5);
imshow(squeeze(vol_data(:, :, slice_idz(5))), [min(vol_data(:)), max(vol_data(:))]);
title('Axial Slice 5');

subplot(2, 3, 6);
imshow(squeeze(vol_data(:, :, slice_idz(6))), [min(vol_data(:)), max(vol_data(:))]);
title('Axial Slice 6');

% % Plot coronal slices (Y direction)
% subplot(2, 3, 4);
% imshow(squeeze(vol_data(:, slice_idy(1), :)), []);
% title('Coronal Slice 1');
% 
% subplot(2, 3, 5);
% imshow(squeeze(vol_data(:, slice_idy(2), :)), []);
% title('Coronal Slice 2');
% 
% % Plot sagittal slices (X direction)
% subplot(2, 3, 6);
% imshow(squeeze(vol_data(slice_idxz(1), :, :)), []);
% title('Sagittal Slice 1');

%% 

%  % Flatten the 3D volume data and voxel coordinates for scatter3
% %[x, y, z] = ndgrid(1:size(hdr.vol, 1), 1:size(hdr.vol, 2), 1:size(hdr.vol, 3));
% values = hdr.vol(:);  % Flatten the volume data to a vector
% %x = x(:); y = y(:); z = z(:);  % Flatten coordinates to vectors
% 
% max(values)
% min(values)
% 
% % Filter out zero or low-intensity values (optional)
% threshold = 0.1 %0.1 * max(values);  % Set a threshold as needed
% valid = values > threshold;  % Logical indexing for valid data points
% x = xx(valid); y = yy(valid); z = zz(valid);
% values = values(valid);
% 
% % Normalize the values for color mapping (0 to 1 range)
% normalizedValues = (values - min(values)) / (max(values) - min(values));
% 
% % Plot using scatter3 with color mapped to the magnitude of each point
% figure;
% scatter3(x, y, z, 10, values, 'filled');  % Adjust size (20) as needed
% colormap hot;  % Use the 'hot' colormap as specified
% colorbar;  % Display color scale
% xlabel('X');
% ylabel('Y');
% zlabel('Z');
% title('Scalar Field Visualization (Dot Plot with Color Mapping)');
% axis tight;
% %grid on;


%%
% Ensure hdr.vol is properly loaded

% Flatten the volume data into a 1D array
voxelData = hdr.vol(:);
len = length(voxelData);
b = linspace(0, length(voxelData),length(voxelData));
scatter(b,voxelData);


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
    xlabel('Voxel Intensity');
    ylabel('Frequency');
    title('Frequency of Voxel Intensities in hdr.vol');
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
