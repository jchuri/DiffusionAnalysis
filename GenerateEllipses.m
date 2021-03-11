function [x_points, y_points] = GenerateEllipses(x_var, y_var, x_cen, y_cen)
%% Description
% This function generates the data points used to plot the atomic clouds in
% DiffusionAnalysis.m. This takes the variances along x and y, uses them to
% calculate the full-width half-max (FWHM) along those directions, and
% generates lists of points for the x and y axes based on those values.
% THIS MUST BE SAVED IN THE SAME DIRECTORY AS "DiffusionAnalysis.m"!
%% INPUTS
    % x_var, y_var = Variance along x and y, which act as the "radii" for
    % the ellipse
    % x_cen, y_cen = Coordinates of the center of the ellipse
        % Note: These values are outputs of "TwoDGaussianFitting.m"
%% OUTPUTS
    % x_points, y_points = Arrays of values that, when combined, define the
    % coordinates of each point used to draw the ellipse
%%
% Full-width half-max of a Gaussian is based on standard deviation
% Specifically, about 2.355 standard deviations away
FWHM_X = 2.355*sqrt(x_var);
FWHM_Y = 2.355*sqrt(y_var);
% Generate the points to draw the ellipse representing the FWHM of the
% cloud. Here's one way I found to do this, based on this webpage:
% https://www.mathworks.com/matlabcentral/answers/129273-to-plot-an-ellipse
theta = linspace(0, 2*pi, 1000);
x_points = FWHM_X*cos(theta) + x_cen;
y_points = FWHM_Y*sin(theta) + y_cen;
end