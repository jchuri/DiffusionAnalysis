function [diff_x, diff_y, G_x, G_z, diff_xFit, diff_yFit] = CalculateDiffusionCoefs(times, x_var, y_var, PXtoCM)
%% Description
% This function calculates the diffusion coefficients for the atomic cloud,
% both in terms of conventional units and in the unitless factor form that
% Grynberg often used (the "G-factors", we'll call them). It also fits the
% diffusion coefficients with a best fit line that will be plotted later.
% THIS MUST BE SAVED IN THE SAME DIRECTORY AS "DiffusionAnalysis.m"!
%% INPUTS
    % times = List of expansion times
    % x_var, y_var = Variance along x and y, obtained from "TwoDGaussianFitting.m"
%% OUTPUTS
    % diff_x, diff_y = Diffusion coefficients along the z and x directions,
    % respectively. Kinda weird, I know, but that's what we've determined
    % the correct axes to be.
    % G_x, G_z = Grynberg's unitless factors corresponding to the diffusion
    % along the x and z directions
    % diff_xFit, diff_yFit = Terms obtained after evaluating the linear fit
    % at each expansion time (I think??)
%%
% We're going to find it for x and y because it might be different
% First start by pulling out variables for a linear fit
% Not quite sure how polyfit works, but it seems to, so...
linFitCoef_x = polyfit(times, x_var, 1);
linFitCoef_y = polyfit(times, y_var, 1);
% This fit is to pixels, so convert it into cm
conversionFactor = PXtoCM^2;
diff_x = linFitCoef_x(1) * conversionFactor * 0.5;
diff_y = linFitCoef_y(1) * conversionFactor * 0.5;
% Convert the diffusion constants into the unitless factors that Grynberg
% often used in his diffusion graphs (let's call them the "G-Factors")
h = 6.63e-34;
RbMass = 1.4e-25;
G_x = 2*pi*RbMass*(diff_y/10000)/h;
G_z = 2*pi*RbMass*(diff_x/10000)/h;
% Print out the diffusion constants and G-factors
disp(['Z Diffusion = ', num2str(diff_x)]);
disp(['X Diffusion = ', num2str(diff_y)]);
disp(['Gz = ', num2str(G_z)]);
disp(['Gx = ', num2str(G_x)]);
% Once you have those variables for your line, plug them into a function to
% get out a set of actualy fitted lines
diff_xFit = polyval(linFitCoef_x, times);
diff_yFit = polyval(linFitCoef_y, times);
end