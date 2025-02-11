%% Open closed loop coils and find the corresponding A-field
% load('closed_loop_data.mat');

% loadedData = load('closed_loop_data.mat');
loadedData = load('test_save.mat');
% app = loadedData.app;
% Assign back to app properties
% app.planarflag = planarflag;
% app.Spinner.Value = spinnerValue;
% app.PlanesizeEditField.Value = plansize;
% app.opticoeff = opticoeff;
% app.lookupinv = lookupinv;
% app.my_rot = my_rot;
% app.my_coillift = my_coillift;
% app.roi = roi_radius;
% app.res = resolutions;
%% (Find & Plot) the transformed_points
figure;
hold on;
for i = 1:length(contours)
    contour_points = contours{i}.points;
    contour_points_3D = [contour_points; zeros(1, size(contour_points, 2))];
    transformed_points = (app.my_rot * contour_points_3D);
    if contours{i}.level < 0
        plot_colour = 'g'; % Yellow for negative levels
    else
        plot_colour = 'bl'; % Blue for others
    end
    plot(transformed_points(1, :), transformed_points(2, :), plot_colour, 'LineWidth', 0.8);
end

title('Stream Function & Equally Spaced Coil Windings using Optimised Coefficients');
xlabel('X');
ylabel('Y');
axis equal;
hold off;

%%
gridpoints = ((app.roi*2)/app.res)+1;

[xx,yy,zz]= meshgrid(linspace(-app.roi,app.roi,gridpoints),linspace(-app.roi,app.roi,gridpoints),linspace(-app.roi,app.roi,gridpoints));

app.sup(:,1)=reshape(xx,[],1);
app.sup(:,2)=reshape(yy,[],1);
app.sup(:,3)=reshape(zz,[],1);

no = vecnorm(app.sup,2,2);
app.sup(no>app.roi,:) =  [];
app.vol=zeros(size(app.sup,1),3);

Mu0=            1E-7;           % Magnetic field constant/4*pi
preconst=       Mu0;

sups=size(app.sup,1);
supi=app.sup;

%% Function to check if the contour lines are clockwise or not
function isClockwise = checkClockwise(points)
    x = points(1, :); % x coords
    y = points(2, :); % y coords
    area = sum((x(2:end) - x(1:end-1)) .* (y(2:end) + y(1:end-1)));
    isClockwise = area < 0; % Clockwise if area is negative
end

%% How to access the contours
% first_contour = contours(50);
% testcontour_points = first_contour{1}.points;
% testcontour_level = first_contour{1}.level;

%% Flip order of the contour points to ensure all the contours run clockwise
figure;
hold on;

for i = 1:length(contours)
    contour_points = contours{i}.points;
    contour_points_3D = [contour_points; zeros(1, size(contour_points, 2))];
    transformed_points = (app.my_rot * contour_points_3D);
    transformed_2D = transformed_points(1:2, :);
    isClockwise = checkClockwise(transformed_2D);
    % Display result in console
    if isClockwise
        fprintf('Contour %d (Level %.2f) is CLOCKWISE\n', i, contours{i}.level);
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
%%
for i = 1:length(contours)
    contour_points = (contours{i}.points);
    contour_level = (contours{i}.level);
    contour_points_3D = [contour_points; zeros(1, size(contour_points, 2))];
    transformed_points = (contour_points_3D' * app.my_rot);
    coilpath2 = transformed_points;
    coilpath = coilpath2 + app.my_coillift;
    coilpath(end+1,:) = coilpath(1,:);

    coilvec = diff(coilpath);
    coilpath(end,:) = [];
    parfor i=1:sups 
        normdist = sqrt((coilpath(:,1)-supi(i,1)).*(coilpath(:,1)-supi(i,1))...
                       +(coilpath(:,2)-supi(i,2)).*(coilpath(:,2)-supi(i,2))...
                       +(coilpath(:,3)-supi(i,3)).*(coilpath(:,3)-supi(i,3))); 
        % 'normdist' should have same dimensions as the 'coilpath'
        dAdtcoil(i,:) = preconst*sum((1./normdist).*coilvec(:,:),1); % same size as supi
    end

    
    if contour_level >= 0
        app.vol=app.vol+dAdtcoil/size(coilpath,3); %divide by the number of strands
    else
        app.vol=app.vol-dAdtcoil/size(coilpath,3);
    end

    clear dAdtcoil;
end
%% Plot app.vol - same as process_coil.mlapp
figure;            
currentColormap = colormap(hot);
n=20;
q=quiver3(app.sup(:,1:n:end),app.sup(:,2:n:end),app.sup(:,3:n:end),app.vol(:,1:n:end),app.vol(:,2:n:end),app.vol(:,3:n:end),2);

mags = sqrt(sum(cat(2, q.UData(:), q.VData(:), ...
                        reshape(q.WData, numel(q.UData), [])).^2, 2));

[~, ~, ind] = histcounts(mags, size(currentColormap, 1));

%// Now map this to a colormap to get RGB
cmap = uint8(ind2rgb(ind(:), currentColormap) * 255);
cmap(:,:,4) = 255;
cmap = permute(repmat(cmap, [1 3 1]), [2 1 3]);

set(q.Head, ...
    'ColorBinding', 'interpolated', ...
    'ColorData', reshape(cmap(1:3,:,:), [], 4).');  % Apply colors to arrow heads

set(q.Tail, ...
    'ColorBinding', 'interpolated', ...
    'ColorData', reshape(cmap(1:2,:,:), [], 4).'); 

title('3D Vector Field with Colored Magnitudes');
xlabel('X');
ylabel('Y');
zlabel('Z');
axis equal;

% the 'app.vol' equivalent is generated again from the closed loop contour coils




