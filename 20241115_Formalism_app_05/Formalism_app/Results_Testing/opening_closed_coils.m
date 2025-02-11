% Load the Closed Loop Coils 
load('closed_loop_data.mat');

% Assign back to app properties
app.planarflag = planarflag;
app.Spinner.Value = spinnerValue;
app.PlanesizeEditField.Value = plansize;
app.opticoeff = opticoeff;
app.lookupinv = lookupinv;
app.my_rot = my_rot;
app.my_coillift = my_coillift;
%%
syms x y l m ps %symbolic variables - indicates that they are variables; m & l are mode numbers of sinusoidal functions; ps is planar surface size
Y(x,y,l,m,ps) = (sin(l*(x+ps/2)*pi/ps)).*(sin(m*(y+ps/2)*pi/ps)); % basis function for the current distribution;
Y_grad(x,y,l,m,ps) = [Y(x-0.5,y,l,m,ps)-Y(x+0.5,y,l,m,ps), Y(x,y-0.5,l,m,ps)-Y(x,y+0.5,l,m,ps), 0]; % gradient of basis function - vector field

% symbolic functions are initialized
f_cd(x,y) = x + y;
f_cd(x,y) = 0;
f_grad(x,y) = x+y; %gradient component
f_grad(x,y) = 0;

for m = 1:app.Spinner.Value
    for n = 1:app.Spinner.Value
        % Compute mode contribution using symbolic gradient
        % f_cd = f_cd + app.CSUMM(app.lookupinv(m, n), app.Spinner.Value) * Y(x,y,m,n,app.PlanesizeEditField.Value);      
        p = app.lookupinv(m, n);
        f_cd = f_cd + app.opticoeff(app.lookupinv(m, n), 1) * Y(x,y,m,n,app.PlanesizeEditField.Value);          
    end
end
f_cd_numeric = matlabFunction(f_cd, 'Vars', [x, y]); %convert to Numerical Function - meant to be quicker; converting to double is also a numerical

%%
% Find approx. min/max values of the stream function
x_range = linspace(-app.PlanesizeEditField.Value/2, app.PlanesizeEditField.Value/2, 100); % X-axis
y_range = linspace(-app.PlanesizeEditField.Value/2, app.PlanesizeEditField.Value/2, 100); % Y-axis
[X, Y] = meshgrid(x_range, y_range); 
Sr = f_cd_numeric(X, Y); %use same variables as paper
Nw = 30; % No. of windings
Sr_min = min(Sr(:));
Sr_max = max(Sr(:));

%%
figure;
hold on;
imagesc(x_range, y_range, Sr); % Plot the stream function as a heatmap
colormap(jet);
colorbar;

for i = 1:length(contours)
    contour_points = contours{i}.points;
   
    if contours{i}.level < 0
        plot_colour = 'w'; % Yellow for negative levels
    else
        plot_colour = 'bl'; % Blue for others
    end

    plot(contour_points(1, :), contour_points(2, :), plot_colour, 'LineWidth', 0.8);
 
end
title('Stream Function & Equally Spaced Coil Windings using Optimised Coefficients');
xlabel('X');
ylabel('Y');
axis equal;
hold off;

%% Plot the rotated coils
figure;
hold on;

for i = 1:length(contours)
    contour_points = contours{i}.points;
    contour_points_3D = [contour_points; zeros(1, size(contour_points, 2))];
    transformed_points = (app.my_rot * contour_points_3D);
    % + app.my_coillift;
    
    if contours{i}.level < 0
        plot_colour = 'g'; % Yellow for negative levels
    else
        plot_colour = 'bl'; % Blue for others
    end

    % plot(contour_points(1, :), contour_points(2, :), plot_colour, 'LineWidth', 0.8);
    plot(transformed_points(1, :), transformed_points(2, :), plot_colour, 'LineWidth', 0.8);
end
title('Stream Function & Equally Spaced Coil Windings using Optimised Coefficients');
xlabel('X');
ylabel('Y');
axis equal;
hold off;
