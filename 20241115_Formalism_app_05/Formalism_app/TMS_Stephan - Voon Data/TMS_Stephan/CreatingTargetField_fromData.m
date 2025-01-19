% 18/01/2024
% This code aims to create the target field suitable for the comatose
% application - MASK funciton not working

%% Load the NIFTI file - the data file associated with DLPFC
[fileName, filePath] = uigetfile('*.nii', 'Select a NIfTI file');
if fileName ~= 0
    fullFilePath = fullfile(filePath, fileName);
    disp(['Selected file: ', fullFilePath]);

    % Load the file using the path
    hdr = nifti_load(fullFilePath);
else
    disp('File selection canceled.');
end

%% Translate the voxels to centre at origin

dims = size(hdr.vol);
middle_voxel = floor(dims / 2);
middle_real = hdr.sform * [middle_voxel, 1]'; % Homogeneous coordinates

%%
    % Get the size of the volume and create a grid of voxel indices
    % [dimX, dimY, dimZ] = size(hdr.vol);
[x, y, z] = ndgrid(1:dims(1), 1:dims(2), 1:dims(3));
voxel_coords = [x(:), y(:), z(:), ones(numel(x), 1)]; 
real_coords = (hdr.sform * voxel_coords')';
  
adjusted_coords = real_coords(:, 1:3) - middle_real(1:3)';
%%
%Flatten data & coord
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

ROI = 130;
Res = 2;
gridpoints = (((ROI*2)/Res)+1); %no. of grid points in one direction (diameter / resolution)

[xx,yy,zz]= meshgrid(linspace(-ROI,ROI,gridpoints),linspace(-ROI,ROI,gridpoints),linspace(-ROI,ROI,gridpoints));

sup(:,1)=reshape(xx,[],1);
sup(:,2)=reshape(yy,[],1);
sup(:,3)=reshape(zz,[],1);

no = vecnorm(sup,2,2);
sup(no>ROI,:) =  []; %remove those points outside of ROI

           
%% Check if each row in coord is in sup

isInSup = ismember(coord, sup, 'rows'); %check each row of coord in sup

% Verify if all rows of coord are in sup
allInSup = all(isInSup);

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

% Use `ismember` to check which rows in `sup` exist in `coord`
mask = ismember(sup, coord, 'rows');
% mask = logical(mask);
% mask = double(mask);

sup = (unique(sup, 'rows', 'stable'));
%% 

vol = zeros(size(sup));
vector = [1;2;3];

data = voxelDataNonZero * vector';

%% 

vol(mask,:) = data;
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

% Plot the vectors
figure;
quiver3(x, y, z, u, v, w, 'AutoScale', 'on', 'AutoScaleFactor', 1, 'LineWidth', 1);
hold on;

% Optionally, plot the points in sup
% scatter3(x, y, z, 10, 'filled', 'MarkerFaceAlpha', 0.5);

% Customise the plot
title('3D Vector Field');
xlabel('X');
ylabel('Y');
zlabel('Z');
grid on;
axis equal;

legend('Vectors', 'Points');
