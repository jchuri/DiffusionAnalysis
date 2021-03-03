function [correctedPic, x_var, y_var, x_cen, y_cen] = TwoDGaussianFitting(parsedTXTFiles)
%% Description
% This function takes the list of parsed .txt files, compiles them into an
% average picture of the atomic cloud, then fits that image with a 2D
% Gaussian. From that fitting, the variance along x and y and the x-y
% coordinates of the cloud's center of mass is extracted and returned. THIS
% MUST BE SAVED IN THE SAME DIRECTORY AS "DiffusionAnalysis.m"!
%% INPUTS
    % parsedTXTFiles = Array of data from the .txt files evaluated by
    % "parseTXTFiles.m"
%% OUPUTS
    % correctedPic = Averaged and rescaled picture of the atomic cloud
    % x_var, y_var = Variance along x and y, respectively, as determined by
    % the 2D Gaussian fitting
    % x_cen, y_cen = Coordinates of the cloud's center of mass, as
    % determined by the 2D Gaussian fitting
%%
image = cat(3,parsedTXTFiles{1:end});
avgPic = mean(image, 3);
% Then pull away whatever the background picks up
% This is the "corrected picture" of the averaged clouds, taking away
% whatever light we saw from the background
% Pull off the first column, it's just zeros
% CorPic = AvgPic(2:end,1:end) - AvgBack{i}(2:end,1:end);
correctedPic = avgPic(2:end,1:end);
% Then scale it so most intense point is 1
correctedPic = correctedPic./max(max((correctedPic)));
% Fit the corrected picture with a 2D Gaussian
% This is started with the assumption that the center of the cloud is
% probably right around the middle. The starting variance is just kind of a
% guess. If the fitting isn't working well, maybe try changing these.
% The other option for "Robust" is "LAR"
% LAR works best if there aren't many outliers in the data; Bisquare allows
% for some occasional outliers
[xo,yo,zo] = prepareSurfaceData(1:125,1:125,correctedPic);
ft = fit([xo,yo],zo,'exp(-(((x-xcen)^2)/vx + ((y-ycen)^2)/vy))','StartPoint',[100,100,62,62],'Robust','Bisquare');
% Order currently vx, vy, xcen, ycen
% Ease of computation, both vx and vy are missing a 1/2 in the fitting
% Extract vx and vy from the fitting, save them to their respective variance lists
x_var = ft.vx./2;
y_var = ft.vy./2;
% Extract xcen and ycen, save them to their respective lists. These
% aren't used outside of the loop, but it may be worthwhile to save
% them in case we want to look at them later.
x_cen = ft.xcen;
y_cen = ft.ycen;
end