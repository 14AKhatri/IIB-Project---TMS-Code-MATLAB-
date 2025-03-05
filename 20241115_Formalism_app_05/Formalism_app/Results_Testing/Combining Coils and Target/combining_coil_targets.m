%% 23/02/2025
% Code to produce weighted coils from different targets

%% Load the Data Files from 10 modes, 200mm plane
% amyg_data = matfile('12.02.2025_run.mat');
% FOX_data = matfile('22.02.2025_FoxACCsub_10M_200Plane_coil.mat');
% VMPFC_data = matfile('22.02.2025_VMPFC_10M_200Plane_Coil_coil.mat');

%% Get opticoeff & lookupinv table from Amygdala data
app_data_amyg = load('12.02.2025_run.mat', 'app_data');
opticoeff_amyg = app_data_amyg.app_data.opticoeff;

coord_amy = app_data_amyg.app_data.A_sup;
lookupinv_amy = app_data_amyg.app_data.lookupinv;

clear amyg_data app_data_amyg;
%%

% Get opticoeff from FOX data
app_data_FOX = load('12.02.2025_run.mat','app_data')

opticoeff_FOX = app_data_FOX.app_data.opticoeff;
lookupinv_FOX = app_data_FOX.app_data.lookupinv;
clear FOX_data app_data_FOX;

%%
% Get opticoeff from VMPFC_data
app_data_VMPFC = load('22.02.2025_VMPFC_10M_200Plane_Coil_coil.mat')

opticoeff_VMPFC = app_data_VMPFC.app_data_coil.opticoeff;
% lookupinv_VMPFC = app_data_VMPFC.app_data_coil.lookupinv;
clear VMPFC_data app_data_VMPFC;

%%
function f_cd_numeric = generate_CD(opticoeff,lookupinv)
    syms x y l m ps %symbolic variables - indicates that they are variables; m & l are mode numbers of sinusoidal functions; ps is planar surface size
        Y(x,y,l,m,ps) = (sin(l*(x+ps/2)*pi/ps)).*(sin(m*(y+ps/2)*pi/ps)); % basis function for the current distribution;
        Y_grad(x,y,l,m,ps) = [Y(x-0.5,y,l,m,ps)-Y(x+0.5,y,l,m,ps), Y(x,y-0.5,l,m,ps)-Y(x,y+0.5,l,m,ps), 0]; % gradient of basis function - vector field
        
        % symbolic functions are initialized
        f_cd(x,y) = x + y;
        f_cd(x,y) = 0;
        f_grad(x,y) = x+y; %gradient component
        f_grad(x,y) = 0;
    
        for m = 1:10
            for n = 1:10
                % Compute mode contribution using symbolic gradient
                % f_cd = f_cd + app.CSUMM(app.lookupinv(m, n), app.Spinner.Value) * Y(x,y,m,n,app.PlanesizeEditField.Value);      
                p = lookupinv(m, n);
                f_cd = f_cd + opticoeff(lookupinv(m, n), 1) * Y(x,y,m,n,200);          
            end
        end
        f_cd_numeric = matlabFunction(f_cd, 'Vars', [x, y]); %convert to Numerical Function - meant to be quicker; converting to double is also a numerical
end


function [Sr_min, Sr_max] = min_max_CD(f_cd_numeric)
    x_range = linspace(-200/2, 200/2, 100); % X-axis
    y_range = linspace(-200/2, 200/2, 100); % Y-axis
    [X, Y] = meshgrid(x_range, y_range); 
    Sr = f_cd_numeric(X, Y); %use same variables as paper
    
    Sr_min = min(Sr(:));
    Sr_max = max(Sr(:));
end


function [x_range,y_range,contours,Sr] = generate_contours(Sr_min,Sr_max,f_cd_numeric)
    x_range = linspace(-200/2, 200/2, 100); % X-axis
    y_range = linspace(-200/2, 200/2, 100); % Y-axis
    [X, Y] = meshgrid(x_range, y_range);
    Nw = 30; % No. of windings
    Lk = Sr_min + ((1:Nw) - (1/2)) * (Sr_max - Sr_min) / Nw; % produces a vector of levels
    
    Sr = f_cd_numeric(X, Y);
    C = contourc(x_range, y_range, Sr, Lk);
    
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
end

    %% Generate the current distributions for all 
    % f_cd_amyg = generate_CD(opticoeff_amyg, lookupinv_amy);
    % f_cd_FOX = generate_CD(opticoeff_FOX, lookupinv_amy);
    % f_cd_VMPFC = generate_CD(opticoeff_VMPFC, lookupinv_amy);

    %% Weighted Coeff.
    a = -1;
    b = 1;
    c = -1;
    % 

    a1 = a/(0.212332559588605);
    b1 = b/(0.16425044);
    c1 = c/(0.118689998);
    app.opticoeff = (a.*(opticoeff_amyg)) + (b.*(opticoeff_FOX)) + (c.*(opticoeff_VMPFC));
    % app.opticoeff = (a1.*(opticoeff_amyg)) + (b1.*(opticoeff_FOX)) + (c1.*(opticoeff_VMPFC));
    app.lookupinv = lookupinv_amy;
    %%
    % f_cd = generate_CD(app.opticoeff,lookupinv_amy);
    f_cd = generate_CD(app.opticoeff,lookupinv_FOX);
    [min,max] = min_max_CD(f_cd)
%%

%%
% [amy_min,amy_max] = min_max_CD(f_cd_amyg)
% [FOX_min,FOX_max] = min_max_CD(f_cd_FOX)
% [VMPFC_min,VMPFC_max] = min_max_CD(f_cd_vmpfc)


%%
[x_range,y_range,contours,Sr] = generate_contours(min,max,f_cd);
%%
figure;
hold on;
imagesc(x_range, y_range, Sr); % Plot the stream function as a heatmap
colormap(jet);
colorbar;
for i = 1:length(contours)
    contour_points = contours{i}.points;
    plot(contour_points(1, :), contour_points(2, :), 'b', 'LineWidth', 0.8);
end
% title('Stream Function & Equally Spaced Coil Windings using Optimised Coefficients');
title(sprintf('a = %.2f, b = %.2f, c = %.2f', a, b, c));
xlabel('X');
ylabel('Y');
axis equal;
hold off;
