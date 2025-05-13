% % 09/05/2025 - Post S.Goetz meeting on 08.05.2025, need to plot the coil ontop of target

%% Load the NIFTI file - the full head model
[fileName, filePath] = uigetfile('*.nii', 'Select a NIfTI file');
if fileName ~= 0
    fullFilePath = fullfile(filePath, fileName);
    disp(['Selected file: ', fullFilePath]);

    % Load the file using the path
    hdr = nifti_load(fullFilePath);
else
    disp('File selection canceled.');
end
%% Load the NIFTI file - DLPFC target the full head model
[fileName, filePath] = uigetfile('*.nii', 'Select a NIfTI file');
if fileName ~= 0
    fullFilePath = fullfile(filePath, fileName);
    disp(['Selected file: ', fullFilePath]);

    % Load the file using the path
    targ = nifti_load(fullFilePath);
else
    disp('File selection canceled.');
end

%% remove spurious point from ONLY FULL HEAD MODEL
hdr.vol(hdr.vol < 1800 ) = 0; 
% volshow(hdr.vol);
%% Plotting in scalar
dims = size(hdr.vol);
[x, y, z] = ndgrid(1:dims(1), 1:dims(2), 1:dims(3)); % Grid of voxel indices (not translation applied to NIFTI voxel grid
voxel_coords = [x(:), y(:), z(:), ones(numel(x), 1)]; 
real_coords = (hdr.sform * voxel_coords'); % real coordinates - each column has vector
clear x y z voxel_coords;
%%
voxelData = hdr.vol(:);
x = real_coords(1,:);
y = real_coords(2,:);
z = real_coords(3,:);

nonZeroIdx = voxelData ~= 0;  % Indices of non-zero voxels
zeroIdx = voxelData == 0;     % Indices of zero voxels

brain_model_data = voxelData(nonZeroIdx);

realXNonZero = x(nonZeroIdx);
realYNonZero = y(nonZeroIdx);
realZNonZero = z(nonZeroIdx);
% store real coordinates of non-zero voxels in coords
head_coord(:,1)=reshape(realXNonZero,[],1);
head_coord(:,2)=reshape(realYNonZero,[],1);
head_coord(:,3)=reshape(realZNonZero,[],1);
clear x y z voxelData realXNonZero realYNonZero realZNonZero;
%%
figure
s = scatter3(head_coord(:,1), head_coord(:,2), head_coord(:,3), ...
             10, brain_model_data, 'filled');  

% head model transparency
s.MarkerFaceAlpha = 0.051;  
hold on;
%%
voxelData = targ.vol(:);
x = real_coords(1,:);
y = real_coords(2,:);
z = real_coords(3,:);

nonZeroIdx = voxelData ~= 0;  % Indices of non-zero voxels
zeroIdx = voxelData == 0;     % Indices of zero voxels

targ_data = voxelData(nonZeroIdx);

realXNonZero = x(nonZeroIdx);
realYNonZero = y(nonZeroIdx);
realZNonZero = z(nonZeroIdx);
% store real coordinates of non-zero voxels in coords
targ_coord(:,1)=reshape(realXNonZero,[],1);
targ_coord(:,2)=reshape(realYNonZero,[],1);
targ_coord(:,3)=reshape(realZNonZero,[],1);
clear x y z voxelData realXNonZero realYNonZero realZNonZero;
%%
scatter3(targ_coord(:,1), targ_coord(:,2), targ_coord(:,3), ...
             10, targ_data, 'filled');  % 'filled' helps show color better
view(-142.3015,   16.1799); 
hold off;
%%
clear targ rgbVolume nonZeroIdx hdr fileName filePath fullFilePath base zeroIdx s
%% Open and Plot Coil
[file, path] = uigetfile({'*.mat', 'MAT-files (*.mat)'}, 'Select a MAT File');

if isequal(file, 0)
    disp('User canceled file selection.');
else
    fullFilePath = fullfile(path, file);
    disp(['Selected file: ', fullFilePath]);

    % Load the MAT file
    loadedData = load(fullFilePath); 
    disp('MAT file loaded successfully.');
end
clear file path fullFilePath
%%
app = loadedData.app_data_coil;
%%
app.roi_radius = 100;
app.resolution = 2;
app.PlanesizeEditField.Value = 300;

%% Recreate the current distribution (assumes 10 modes)
syms x y l m ps %symbolic variables - indicates that they are variables; m & l are mode numbers of sinusoidal functions; ps is planar surface size
% Y(x,y,l,m,ps) = (sin(l*(x+ps/2)*pi/ps)).*(sin(m*(y+ps/2)*pi/ps)); % basis function for the current distribution;
% Y_grad(x,y,l,m,ps) = [Y(x-0.5,y,l,m,ps)-Y(x+0.5,y,l,m,ps), Y(x,y-0.5,l,m,ps)-Y(x,y+0.5,l,m,ps), 0]; % gradient of basis function - vector field

Y(x,y,l,m,ps) = (sin(l*(x+ps/2)*pi/ps)) .* (sin(m*(y+ps/2)*pi/ps));  % basis function for the current distribution

f_cd(x,y) = x + y;
f_cd(x,y) = 0;

for n = 1:10
    for m = 1:10
        % Compute mode contribution using symbolic gradient
        % f_cd = f_cd + app.CSUMM(app.lookupinv(m, n), 10) * Y(x,y,m,n,app.PlanesizeEditField.Value);      
        f_cd = f_cd + app.opticoeff(app.lookupinv(m, n), 1) * Y(x,y,n,m,app.PlanesizeEditField.Value);      
        % f_cd = f_cd + app.opticoeff(app.lookupinv(m, n), 1) * Y(x,y,m,n,app.PlanesizeEditField.Value); 
    end
end
f_cd_numeric = matlabFunction(f_cd,'Vars',[x,y]); %convert to Numerical Function - meant to be quicker; converting to double is also a numerical
%%
clear f_cd
%%
x_range = linspace(-app.PlanesizeEditField.Value/2, app.PlanesizeEditField.Value/2, 100); % X-axis
y_range = linspace(-app.PlanesizeEditField.Value/2, app.PlanesizeEditField.Value/2, 100); % Y-axis
[X, Y] = meshgrid(x_range, y_range); 

Sr = f_cd_numeric(X, (Y)); %use same variables as paper

Nw = 30; % No. of windings

Sr_min = min(Sr(:));
Sr_max = max(Sr(:));

% levels in magnitude of Sr (stream funciton) for each
% contour line
Lk = Sr_min + ((1:Nw) - (1/2)) * (Sr_max - Sr_min) / Nw; % produces a vector of levels

C = contourc(x_range, (y_range), (Sr), Lk);
                
% Iterate through  + find points in each contour line
contours = {}; % array for paths for each Lk step
index = 1;
while index < size(C, 2)
    level = C(1, index); % Current contour level
    num_points = C(2, index); % Number of points in this contour
    contour_points = C(:, index+1:index+num_points); % Extract contour points
    contours{end+1} = struct('level', level, 'points', contour_points);
    index = index + num_points + 1; % Move to the next contour
end

figure;
hold on;
imagesc(x_range, y_range, (Sr)); % Plot the stream function as a heatmap
set(gca, 'YDir', 'reverse');

colormap(jet);
colorbar;
for i = 1:length(contours)
    contour_points = contours{i}.points;
    plot(contour_points(1, :), contour_points(2, :), 'b', 'LineWidth', 0.8);
end
%title('Stream Function & Equally Spaced Coil Windings using Optimised Coefficients');
title('Post-Optimisation Current Distribution');
% title('Pre-Optimisation Current Distribution');
xlabel('X (mm)');
ylabel('Y (mm)');
% [~, h] = contour(x_range, y_range, Sr, Lk);
% clabel(C, h);
axis equal;
hold off;

%%
% % %%
% % test_contour = contours{15}.points;
% % test_contour_3D = [test_contour; zeros(1, size(test_contour, 2))]
% % transformed_points = (contour_points_3D' * app.my_rot)'; 
% % % plot3(test_contour_3D(1, :), test_contour_3D(2, :), test_contour_3D(3, :), plot_colour, 'LineWidth', 0.8)
% % 
% % plot3(transformed_points(1, :), transformed_points(2, :), transformed_points(3, :), plot_colour, 'LineWidth', 0.8)

%%
%%

%% (Find & Plot) the transformed_points
figure
s = scatter3(head_coord(:,1), head_coord(:,2), head_coord(:,3), ...
             10, brain_model_data, 'filled');  

s.MarkerFaceAlpha = 0.051;  
hold on;

scatter3(targ_coord(:,1), targ_coord(:,2), targ_coord(:,3), ...
             10, targ_data, 'filled');  % 'filled' helps show color better
view(-142.3015,   16.1799); 
hold on;

for i = 1:length(contours)
    contour_points = contours{i}.points; % Contour points in 2D
    % contour_points = contour_points([2, 1], :);
    contour_points_3D = [contour_points; zeros(1, size(contour_points, 2))]; % Convert to 3D points
    % transformed_points = (app.my_rot * contour_points_3D); 
    % transformed_points = (contour_points_3D' * app.my_rot'); 

    z_lift = 100;  % Desired height
    % contour_points_3D = [contour_points; z_lift * ones(1, size(contour_points, 2))];
    coilpath = (contour_points_3D' * (app.my_rot)) + app.my_coillift;

    % coilpath = (app.my_rot * contour_points_3D);% + app.my_coillift;
    % coilpath = (contour_points_3D) + app.my_coillift';
    coilpath = coilpath';
    
    if contours{i}.level < 0
        plot_colour = 'g'; % Yellow for negative levels
    else
        plot_colour = 'b'; % Blue for others
    end
    % plot(contour_points_3D(1, :), contour_points_3D(2, :), plot_colour, 'LineWidth', 0.8);
   plot3(coilpath(1, :), coilpath(2, :), coilpath(3, :), plot_colour, 'LineWidth', 1)
end

title('Stream Function & Equally Spaced Coil Windings using Optimised Coefficients');
xlabel('X');
ylabel('Y');
view(-142.3015,   16.1799); 
axis equal;
hold off;

%% % Analyse Coil A-Field
coil_field = zeros(size(app.A_sup,1),3); % initialise field from coil
Mu0=            1E-7;           % Magnetic field constant/4*pi
preconst=       Mu0;
sups=size(app.A_sup,1);
supi=app.A_sup; %coordinates

%%
function isClockwise = checkClockwise(points)
    x = points(1, :); % x coords
    y = points(2, :); % y coords
    area = sum((x(2:end) - x(1:end-1)) .* (y(2:end) + y(1:end-1)));
    isClockwise = area < 0; % Clockwise if area is negative
end
%% Flip order of the contour points to ensure all the contours run clockwise
figure;
hold on;

for i = 1:length(contours)
    contour_points = contours{i}.points;
    isClockwise = checkClockwise(contour_points);

    if isClockwise
        fprintf('Contour %d (Level %.2f) is CLOCKWISE\n', i, contours{i}.level); %Display for each contour
    else
        % fprintf('Contour %d (Level %.2f) is ANTICLOCKWISE\n', i, contours{i}.level);
        contours{i}.points = fliplr(contours{i}.points);
    end
    plot(contours{i}.points(1, :), contours{i}.points(2, :), 'g', 'LineWidth', 0.8);
end

title('Stream Function & Equally Spaced Coil Windings using Optimised Coefficients');
xlabel('X');
ylabel('Y');
axis equal;
% All the contours (in contours) are clockwise now

%% Generate the coil field
for i = 1:length(contours)
    contour_points = (contours{i}.points);
    contour_level = (contours{i}.level);
    
    contour_points = contour_points([2, 1], :);
    contour_points_3D = [contour_points; zeros(1, size(contour_points, 2))]; % Convert to 3D points
    
    coilpath = (contour_points_3D' * (app.my_rot)) + app.my_coillift;
    % coilpath = coilpath';
    
    coilpath(end+1,:) = coilpath(1,:); %extend by one so 'coilvec' is same dimensions
    coilvec = diff(coilpath);
    coilpath(end,:) = []; %remove the additional coilpath entry

    parfor i=1:sups 
        normdist = sqrt((coilpath(:,1)-supi(i,1)).*(coilpath(:,1)-supi(i,1))...
                       +(coilpath(:,2)-supi(i,2)).*(coilpath(:,2)-supi(i,2))...
                       +(coilpath(:,3)-supi(i,3)).*(coilpath(:,3)-supi(i,3))); 
        % 'normdist' should have same dimensions as the 'coilpath'
        dAdtcoil(i,:) = preconst*sum((1./normdist).*coilvec(:,:),1); % same size as supi
    end

    if contour_level >= 0
        coil_field=coil_field+dAdtcoil/size(coilpath,3); %divide by the number of strands
    else
        coil_field=coil_field-dAdtcoil/size(coilpath,3); %subtract the field if negative level
    end

    clear dAdtcoil;
end

%%
%% Field Plotting Function
function plotVectorField(coord_array, data_array, Title)
    figure;            
    currentColormap = colormap(jet);
    n = 20; % Sampling step

    % Create quiver3 plot
    q = quiver3(coord_array(:,1:n:end), coord_array(:,2:n:end), coord_array(:,3:n:end), ...
                data_array(:,1:n:end), data_array(:,2:n:end), data_array(:,3:n:end), 2);

    % Compute magnitudes of the vectors
    mags = sqrt(sum(cat(2, q.UData(:), q.VData(:), ...
                         reshape(q.WData, numel(q.UData), [])).^2, 2));

    % Map magnitudes to colormap indices
    [~, ~, ind] = histcounts(mags, size(currentColormap, 1));

    % Convert indices to RGB color mapping
    cmap = uint8(ind2rgb(ind(:), currentColormap) * 255);
    cmap(:,:,4) = 255; % Set alpha channel
    cmap = permute(repmat(cmap, [1 3 1]), [2 1 3]);

    % Apply colors to arrow heads
    set(q.Head, ...
        'ColorBinding', 'interpolated', ...
        'ColorData', reshape(cmap(1:3,:,:), [], 4).');

    % Apply colors to arrow tails
    set(q.Tail, ...
        'ColorBinding', 'interpolated', ...
        'ColorData', reshape(cmap(1:2,:,:), [], 4).'); 

    % Add labels and title
    title(Title);
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    axis equal;
end
%%
plotVectorField(app.A_sup, app.A_Aim, 'Target Field'); % Target Field 
view(-142.3015,   16.1799); 
%%
hold on;
plotVectorField(app.A_sup, coil_field.*10, 'From Closed Loop Coils'); % Coil Field
view(-142.3015,   16.1799); 
hold off;
%%
field = @(myccsum) sum(repmat(reshape(myccsum(:), 1, 1, length(myccsum)),size(app.A_modes, 1), size(app.A_modes, 2), 1) .* app.A_modes, 3);
optim_field = field(app.opticoeff);


plotVectorField(app.A_sup, optim_field, 'Optimised Field');
view(-142.3015,   16.1799);
%%
%%
function corr = calculate_corr(field, target_field)
        corr = sum(abs(target_field .* field), 'all');
end
 
function results = energy_in_field(field)
    results = sum(abs(vecnorm(field,2,2)),"all");
end

function average_coord = central_direction(field)
        %weighted average of the coordinates
        field_binmap = any(field ~= 0, 2); % returns a binary map of all non-zero vectors

        scaled_supp = field_binmap .* (app.A_sup); % get the coordinates Nx3
        scaled_vec = field_binmap .* vecnorm(field,2,2); % get the magnitudes of the vectors Nx1

        total_field_mag = sum(scaled_vec);
        weighted_sup = (scaled_supp.*scaled_vec);
        average_coord = sum(scaled_supp.*scaled_vec,1)/total_field_mag;
end

function magcorr = energy_in_target(field, target_field)
        %Returns the total field magnitude in target region
         target_field_binmap = any(target_field ~= 0, 2);
         magcorr = sum((target_field_binmap .* abs(vecnorm(field,2,2))), 'all');
end

%%
optim_field_correlation_target = calculate_corr(optim_field,app.A_Aim);
fprintf('Correlation of Target with Optimised Field: %.2f J\n', optim_field_correlation_target);

coil_field_correlation_target = calculate_corr(coil_field,app.A_Aim);
fprintf('Correlation of Target with Coil Field: %.2f J\n', coil_field_correlation_target*1000000);