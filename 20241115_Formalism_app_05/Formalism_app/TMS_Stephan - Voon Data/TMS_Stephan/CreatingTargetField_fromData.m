% 18/01/2024
% Creates target field suitable for the comatose applications 

%% Load the NIFTI file - the data file (not mask) associated with DLPFC
[fileName, filePath] = uigetfile('*.nii', 'Select a NIfTI file');
if fileName ~= 0
    fullFilePath = fullfile(filePath, fileName);
    disp(['Selected file: ', fullFilePath]);

    % Load the file using the path
    hdr = nifti_load(fullFilePath);
else
    disp('File selection canceled.');
end

%% Translate the voxel locations to centre at origin - original NIFTI file only in the 1st quadrant
dims = size(hdr.vol);
middle_voxel = floor(dims / 2);
middle_real = hdr.sform * [middle_voxel, 1]'; % finds equivalent real coordinates of centre point

%%
[x, y, z] = ndgrid(1:dims(1), 1:dims(2), 1:dims(3)); % Grid of voxel indices (not translation applied to NIFTI voxel grid
voxel_coords = [x(:), y(:), z(:), ones(numel(x), 1)]; 
real_coords = (hdr.sform * voxel_coords')'; % real coordinates
  
adjusted_coords = real_coords(:, 1:3) - middle_real(1:3)'; % real coordinates translated 
%% Find coordinates of non-zero voxel

voxelData = hdr.vol(:);
x = adjusted_coords(:,1);
y = adjusted_coords(:,2);
z = adjusted_coords(:,3);

nonZeroIdx = voxelData ~= 0;  % Indices of non-zero voxels
zeroIdx = voxelData == 0;     % Indices of zero voxels

voxelDataNonZero = voxelData(nonZeroIdx);
realXNonZero = x(nonZeroIdx);
realYNonZero = y(nonZeroIdx);
realZNonZero = z(nonZeroIdx);
% store real coordinates of non-zero voxels in coords
coord(:,1)=reshape(realXNonZero,[],1);
coord(:,2)=reshape(realYNonZero,[],1);
coord(:,3)=reshape(realZNonZero,[],1);

noc = vecnorm(coord,2,2);
    %% 

% % Flatten the voxel data and coordinates into 1D arrays
% voxelData = hdr.vol(:);
% x = x(:);
% y = y(:);
% z = z(:);
% 
% % Apply the affine transformation (hdr.sform) to get real-world coordinates
% % Voxel coordinates are in [x, y, z, 1] homogeneous coordinates.
% % The sform matrix is a 4x4 matrix
% realCoords = hdr.sform * [x, y, z, ones(length(x), 1)]';  % Apply the sform transformation
% realX = realCoords(1, :);  % Extract real-world x coordinates
% realY = realCoords(2, :);  % Extract real-world y coordinates
% realZ = realCoords(3, :);  % Extract real-world z coordinates
% 
% % Separate the zero and non-zero voxel data
% nonZeroIdx = voxelData ~= 0;  % Indices of non-zero voxels
% zeroIdx = voxelData == 0;     % Indices of zero voxels
% 
% % Non-zero voxel data and coordinates in real-world space
% voxelDataNonZero = voxelData(nonZeroIdx);
% realXNonZero = realX(nonZeroIdx);
% realYNonZero = realY(nonZeroIdx);
% realZNonZero = realZ(nonZeroIdx);
% 
% % Zero voxel data and coordinates in real-world space
% voxelDataZero = voxelData(zeroIdx);
% realXZero = realX(zeroIdx);
% realYZero = realY(zeroIdx);
% realZZero = realZ(zeroIdx);
% 
% coord(:,1)=reshape(realXNonZero,[],1);
% coord(:,2)=reshape(realYNonZero,[],1);
% coord(:,3)=reshape(realZNonZero,[],1);


%% Define the coordinates as per comatose application

ROI = 130; % ROI and Res as per comatose application
Res = 2;
gridpoints = (((ROI*2)/Res)+1); %no. of grid points in one direction (diameter / resolution)

[xx,yy,zz]= meshgrid(linspace(-ROI,ROI,gridpoints),linspace(-ROI,ROI,gridpoints),linspace(-ROI,ROI,gridpoints));

sup(:,1)=reshape(xx,[],1);
sup(:,2)=reshape(yy,[],1);
sup(:,3)=reshape(zz,[],1);

no = vecnorm(sup,2,2);
sup(no>ROI,:) =  []; %remove those coordinate points outside of ROI

           
%% Check if each row in coord is in sup (each target point is in sup)
    %Used to check the entire MNI 2mm brain model fits in spherical region

isInSup = ismember(coord, sup, 'rows'); %check each row of coord in sup

allInSup = all(isInSup); % Verify if all rows of coord are in sup - returns true or false

% Display result
if allInSup
    disp('All coordinates in coord are present in sup.');
else
    disp('Not all coordinates in coord are present in sup.');
end

missingCoords = coord(~isInSup, :);
disp('Coordinates in coord but not in sup:');
disp(missingCoords);

%%
% mask = false(size(sup, 1), 1);

% Use `ismember` to check which rows of `sup` exist in `coord`
mask = ismember(sup, coord, 'rows'); % generate a for sup - identify coordinates of sup which are in coord
% mask = logical(mask);
% mask = double(mask);

% sup = (unique(sup, 'rows', 'stable'));
%% 
vol = zeros(size(sup)); % initialise vol grid - for data
% vector = [1;2;3];
vector = [0.5285;-0.7362;-0.4226]; % target normal vector of surface

data = voxelDataNonZero * vector'; %multiply the scalar Voon data from normal vector
vol(mask,:) = data; % put in vectors in vol
%% PLOT quiver plot with sphere
% sup is an N x 3 matrix of coordinates
% vol is an N x 3 matrix of vector components

% real coordinates
x = sup(:, 1); 
y = sup(:, 2); 
z = sup(:, 3); 

% vector data
u = vol(:, 1); 
v = vol(:, 2); 
w = vol(:, 3); 

figure;
quiver3(x, y, z, u, v, w, 'AutoScale', 'on', 'AutoScaleFactor', 10, 'LineWidth', 1);
hold on;

radius = 130; % Sphere radius
[sphereX, sphereY, sphereZ] = sphere(50); % Generate sphere coordinates
surf(radius * sphereX, radius * sphereY, radius * sphereZ, ...
    'FaceAlpha', 0.1, 'EdgeColor', 'none', 'FaceColor', 'blue'); % Adjust colour and transparency


% Customise the plot
title('DLPFC Target Vector Field');
xlabel('X');
ylabel('Y');
zlabel('Z');
grid on;
axis equal;

% legend('Vectors', 'Points');

%% 
% Assuming:
% sup is an N x 3 matrix of coordinates
% vol is an N x 3 matrix of vector components

% Extract coordinates
x = sup(:, 1); % X-coordinates
y = sup(:, 2); % Y-coordinates
z = sup(:, 3); % Z-coordinates

% Extract vector components
u = vol(:, 1); % X-component of vectors
v = vol(:, 2); % Y-component of vectors
w = vol(:, 3); % Z-component of vectors

% Calculate vector magnitudes
magnitudes = sqrt(u.^2 + v.^2 + w.^2);

% Normalise the magnitudes for colour mapping
minMag = min(magnitudes);
maxMag = max(magnitudes);
normalisedMag = (magnitudes - minMag) / (maxMag - minMag);

% Create a colormap
cmap = parula(256); % 'parula' is MATLAB's default colormap
colours = interp1(linspace(0, 1, size(cmap, 1)), cmap, normalisedMag);

% Plot the vectors with colours
figure;
hold on;
for i = 1:length(x)
    % Plot each vector with a different colour
    quiver3(x(i), y(i), z(i), u(i), v(i), w(i), 'AutoScale', 'on', ...
        'AutoScaleFactor', 1, 'LineWidth', 1, 'Color', colours(i, :));
end

% Add a colourbar to indicate magnitudes
colormap(cmap);
caxis([minMag, maxMag]);
colorbar;
title('3D Vector Field with Magnitude-Based Colouring');
xlabel('X');
ylabel('Y');
zlabel('Z');
grid on;
axis equal;
