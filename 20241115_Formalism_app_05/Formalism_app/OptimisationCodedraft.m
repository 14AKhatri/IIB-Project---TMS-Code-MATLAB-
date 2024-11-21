%% Script to draft different optimisation functions

%existing optimisation code - unchanges

myresid = @(myccsum) sum((app.A_Aim - sum( repmat( reshape( myccsum(:),1,1,length(myccsum) ), size(app.A_modes,1), size(app.A_modes,2), 1).*app.A_modes, 3)).^2, "all");
% "@(myccsum) -defines anonymous function- indicated that myccsum is input
% variable

myoptions = optimset('TolX', app.StepToleranceEditField.Value, 'TolFun', app.FunctionToleranceEditField.Value, 'maxfunevals', app.MaxFunEvalsEditField.Value, 'maxiter', app.MaxIterationsEditField.Value, 'disp', 'iter');

myinitvalue = app.CSUMM(:,app.maxresSlider.Value);
%initial guess for the coeff

[app.opticoeff, myresidval, myexitflag, myoutput] = fmincon(myresid, myinitvalue, [], [], [], [], -Inf, Inf, [], myoptions);
%%constrained optimisation function
%3rd 6th arguments - empty arrays indicate no linear inequality or
    %equality constraints on the coefficients in the optimisation
%7th & 8th arguments set lower and upper bounds of myccsum

%9th arguments indicates no nonlinear constraints
%myoptions - options for optimization process
    % Tolx - step tolerance - smallest change in coeff that fmincon will
        % attempt
    % TolFun - function tolerance - min. imporvement in objective function
        % before fmincon stops
    % maxfunevals - maximum no. of iterations
    %disp level set to 'iter' - display information at each iteration

    %% 


%Define the logical mask again
r_min = eval(app.R_minEditField.Value);% Minimum radius
r_max = eval(app.R_maxEditField.Value);  % Maximum radius
theta_min = eval(app.theta_minEditField.Value);  % Minimum azimuth (rad) - [-pi, pi]
theta_max = eval(app.theta_maxEditField.Value);   % Maximum azimuth (rad)
phi_min = eval(app.phi_minEditField.Value);    % Minimum elevation (rad)-[-pi/2,pi/2]
phi_max = eval(app.phi_maxEditField.Value);     % Maximum elevation (rad)
gridpoints = (((app.roi*2)/app.res)+1); %no. of grid points in one direction (diameter / resolution)
[xx,yy,zz]= meshgrid(linspace(-app.roi,app.roi,gridpoints),linspace(-app.roi,app.roi,gridpoints),linspace(-app.roi,app.roi,gridpoints));


app.sup(:,1)=reshape(xx,[],1);
app.sup(:,2)=reshape(yy,[],1);
app.sup(:,3)=reshape(zz,[],1);

no = vecnorm(app.sup,2,2);
app.sup(no>app.roi,:) =  []; %remove those points outside of ROI


[azimuth, elevation, radius] = cart2sph(app.sup(:, 1), app.sup(:, 2), app.sup(:, 3));

%Filter points within the specified spherical region
region_mask = (radius >= r_min) & (radius <= r_max) & (azimuth >= theta_min) & (azimuth <= theta_max) & (elevation >= phi_min) & (elevation <= phi_max);
region_mask = zeros(length(radius),1);  % Create a logical array of the same size as radius

for i = 1:length(radius);
    if ((radius(i) >= r_min) & (radius(i) <= r_max) & (azimuth(i) >= theta_min) & (azimuth(i) <= theta_max) & (elevation(i)>= phi_min) & (elevation(i) <= phi_max));
  Check if the radius meets the condition
        region_mask(i) = 1;  % Set the corresponding mask value to true
    end
end

app.vol = zeros(size(app.sup));  %app.sup is the coordinates
app.vol(logical(region_mask), :) = 1;


app.target = app.vol;





myresid = @(myccsum) sum((app.target - sum( repmat( reshape( myccsum(:),1,1,length(myccsum) ), size(app.A_modes,1), size(app.A_modes,2), 1).*app.A_modes, 3)).^2, "all");

%a proxy for the overall energy of approximated field
%squared sum of field with itself
approx_A_field = sum(repmat(reshape(app.myccsum(:), 1, 1, length(app.myccsum)), size(app.A_modes, 1), size(app.A_modes, 2), 1) .* app.A_modes, 3);
field_energy = sum(approx_A_field.^2, 'all');

f = myresid / field_energy


%% function to compare the app.opticoeff & app.ccsum

%function to calculate squared normalised element wise differences in the
%mode coeff.

%ensure same size of the two
if ~isequal(size(app.opticoeff), size(app.CCSUM(:,app.maxresSlider.Value))) %%may need to change the maxresslider.value
    error('opticoeff and CCSUM must have the same dimensions.');
end

coeff_diff = app.opticoeff - app.CCSUM;
squared_diff_sum = sum(coeff_diff.^2, 'all');