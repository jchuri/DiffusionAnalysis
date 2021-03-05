%% Intro and Overview
% This is the 2021 version of "DiffAnalNoBack120619.m". The people who named
% it were obviously very mature.
% Originally written by Ethan Clements and Anthony Rapp, slightly refined by
% Alex Staron, and most recently overhauled by Ian Dilayrd and Jordan Churi.

% While the original code worked just fine, there was definitely room for
% improvement and cleaning-up. This version is meant to do just that while
% maintaining readability and making maintainence easier in the future.

% This script takes a full set of diffusion data and produces average images
% of the MOT at each stage of expansion, as well as a graph of the variance
% along x and z with respect to expansion time. It also calculates the
% diffusion constants along those directions and converts them into
% "G-factors" that Grynberg sometimes reported in his papers on diffusion.

%% Change Log (as of March 5, 2021):
    % Revised some variable names and changed some spacings
    % Preallocated some arrays that hold important data (times, variance_x/y, etc.)
    % Broke apart the operations within the big FOR loop over each timing folder
        % -"parseTXTFile.m" checks to see if each .txt file is acceptable to use
        % -"TwoDGaussianFitting.m" compiles the results of "parseTXTFile.m",
        % creates the averaged and corrected image, fits it with a 2D
        % Gaussian, then extracts the variance along x and y and the
        % coordinates of the cloud's center of mass
        % -"GenerateEllipses.m" takes the variance along x and y, uses them
        % to calculate the FWHM for each direction, then generates lists of
        % coordinates to plot each cloud
%% DATA REQUIREMENTS
% All data folders should be named in the scheme '**ms' for data and '**bg'
% for backgrounds. These folders go into a folder whose name is the date
% that the measurements were taken. (Previous format was MonthDayYear
% (example: 073020), but any format should be OK)

% Assuming that only occasionally you have bad data, it should get averaged
% out by the other 'good' data sets. In a best case scenario, you totally
% missed the cloud and just got the background, in which case it doesn't do
% anything bad at all, since all of these values will get scaled anyways.
% The cloud itself doesn't need to be centered (in fact, it can be
% partially off-screen) since the script fits everything.
%% Selecting the Data Folder
% Manually select the folder containing the folders for each diffusion
% timing. Start in the folder this script is saved to.
clear;
initDir = 'C:\Users\staronal\Desktop\Raw Data\Diffusion'; %IF YOU RUN THIS SCRIPT ON ANOTHER COMPUTER, CHANGE THIS PATH BEFOREHAND!
dataDir = uigetdir('Where are all the data folders?');
cd(dataDir);
% Add initDir as an extra file path.
% This is necessary in order to properly call the custom functions used to
% analyze the data (like parseTXTFile, TwoDGaussianFitting, etc.)
% There's probably a better way to do this whole analysis that would avoid
% this step...
addpath(genpath(initDir));
%% Load the Data and Background (if there is one)
% Create a list of all the data folders and background folders (if there
% are any). Nothing too complicated here.
dataFolderList = dir('*ms');
if exist('bg','dir')
    backgroundFolderList = dir('*bg*');
end
%% Threshold and Preallocation of Arrays
% Set the threshold for the minimum intensity. The maximum intensity within
% a given .txt file must be greater than this threshold to be usable.
% You can pick the 'all' option to sum over every pixel
% The 'one' option just finds the brightest pixel and compares to that
%thresholdall = 5e+7;
thresholdOne = 2500;
% Preallocate arrays to hold all the important data from each expansion
% time (the times themselves, the variances along x and y, averaged images,
% etc.). Unless otherwise specified, each expansion time gets its own
% column in these arrays.
times = zeros(1, length(dataFolderList));
timesMS = zeros(1, length(dataFolderList));
x_var = zeros(1, length(dataFolderList));
y_var = zeros(1, length(dataFolderList));
PrettyPic = cell(1, length(dataFolderList));
x_cen = zeros(1, length(dataFolderList));
y_cen = zeros(1, length(dataFolderList));
ellipse_X = cell(1, length(dataFolderList));
ellipse_Y = cell(1, length(dataFolderList));
%% Extracting and Compiling the Raw Data
% Go through each folder in the list, read through and evaluate the .txt
% files within, and store the results in an array.
% Create a cell array to hold the data from the good files for each timing.
% Each timing gets its own row, and each good file gets its own column.
% MATLAB can automatically extend cell arrays to fit the number of good
% files.
goodFiles = cell(length(dataFolderList), 1);
for i = 1:length(dataFolderList);
    % Extract the timing (in ms) from the folder name. It's saved in a cell.
    timeFromFolderName = regexp(dataFolderList(i).name,'\d+\.?\d*|-\d+\.?\d*|\.?\d+|-\.?\d+','match');
    % Convert the string witin the cell into a number. This number will be
    % used when plotting stuff later.
    times(i) = str2num(timeFromFolderName{1}) * 0.001;
    timesMS(i) = str2num(timeFromFolderName{1});
    % Go into the timing folder and create a list of all the .txt files inside
    cd(dataFolderList(i).name);
    fileList = dir('*.txt');
    % Create a counter to both track the number of good files and to act as
    % an index to place the contents of those files in the array of good files.
    counter = 1;
    % Go through the .txt files and check them for validity. The maximum
    % intensity should be above thresholdOne and the minimum intensity
    % should be positive or 0. See "parseTXTFile.m" for more details.
    for j = 1:length(fileList)
        goodFiles{i, counter} = parseTXTFile(fileList(j).name, thresholdOne, j, dataFolderList(i).name);
        if ~isempty(goodFiles{i, counter})
            counter = counter + 1;
        end
    end
    % Move back up in directory so we can go into the next folder
    cd('..')
end
fprintf('Data extraction complete\n')
%% Generating Composite Images and 2D Gaussian Fitting
% For each timing (each row in goodFiles), create an average image of the
% cloud, fit that image with a 2D Gaussian, and extract the necessary
% information from that fitting.
for i = 1:length(dataFolderList);
    % Compile the data in the good .txt files to create an average and
    % intensity-corrected image of the cloud. Fit that with a 2D Gaussian,
    % then extract the variance along x and y and the coordinates of the
    % cloud's center of mass from the fitting. See "TwoDGaussianFitting.m"
    % for more details.
    [PrettyPic{i}, x_var(i), y_var(i), x_cen(i), y_cen(i)] = TwoDGaussianFitting(goodFiles(i,1:length(goodFiles)));
    % Use the extracted variances to calculate the full-width half-max
    % (FWHM). See "GenerateEllipses.m" for more details.
    [ellipse_X{1,i}, ellipse_Y{1,i}] = GenerateEllipses(x_var(i), y_var(i), x_cen(i), y_cen(i));
end
%% Plotting the Images
% Go through the lists of points and centers for each cloud and plot them
% on one graph. For now, label each center with its corresponding expansion
% time.
figure;
hold on
for i = 1:length(dataFolderList);
    colorNum = rand(1,3); % Make ellipse border, center, and label the same color for easy identification
    plot(ellipse_X{1,i}, ellipse_Y{1,i}, 'Color', colorNum, 'Linewidth', 2);
    plot(x_cen(i), y_cen(i), '*', 'Color', colorNum);
    text(x_cen(i)+12, y_cen(i), [num2str(timesMS(i)) 'ms'], 'HorizontalAlignment', 'center', 'Color', colorNum);
end
xlabel('Pixels (-Z)');
ylabel('Pixels (-X)');
title('Full-Width Half-Max of Atomic Cloud');
hold off
%% Calculating Diffusion Constants and G-Factors
% Now we want to figure out that diffusion constant
% We're going to find it for x and y because it might be different
% First start by pulling out variables for a linear fit
linfitcoeffx = polyfit(times,x_var,1);
linfitcoeffy = polyfit(times,y_var,1);
% This fit is to pixels, so convert it into cm. The conversion factor here
% is subject to change and comes from comparing the size of the cloud in
% the image to the size of the cloud on the CRT TV. Will fill in with more
% details when available.
% Ethan got a value of 5.2e-3, while Schiavoni used 4e-3 in their thesis.
% Other values include 3.6e-3 and 2.85e-3 (the current value)
PXtoCM = 2.85e-3;
ConversionFactor = PXtoCM^2;
diffx = linfitcoeffx(1) * ConversionFactor * 0.5;
diffy = linfitcoeffy(1) * ConversionFactor * 0.5;

disp(['Z Diffusion = ',num2str(diffx)]);
disp(['X Diffusion = ',num2str(diffy)]);

% Adding in a feature that calculates the unitless factor that Grynberg uses
% in most of his diffusion graphs.
% Convert the diffusion constants into the unitless factors that Grynberg
% often used in his diffusion graphs (let's call them the "G-Factors")
h = 6.63e-34;
RbMass = 1.4e-25;
Gx = 2*pi*RbMass*(diffy/10000)/h;
Gz = 2*pi*RbMass*(diffx/10000)/h;

disp(['Gx = ',num2str(Gx)]);
disp(['Gz = ',num2str(Gz)]);

% Once you have those variables for your line, plug them into a function to
% get out a set of actualy fitted lines
vxfit = polyval(linfitcoeffx,times);
vyfit = polyval(linfitcoeffy,times);

% Make a new figure separate from the FWHM plot
figure;
% Plotting the actual data points with *
plot(timesMS,x_var.*ConversionFactor,'*',timesMS,y_var.*ConversionFactor,'*');
legend show
title('Variance vs. Expansion Time');
xlabel('Expansion Time (ms)');
ylabel('Variance \sigma^{2} (cm^{2})');

hold on
% Plot the lienar fit overtop them
plot(timesMS,vxfit.*ConversionFactor,timesMS,vyfit.*ConversionFactor);

% Finally save all of the data
% The first is a set of the variables that we solved for

%save('Solved.mat','vx','vy','xcen','ycen','times');
%save('Diffusions.mat','diffx','diffy');

% The second is a list of the arrays if you wanted to look at the 3D plots
% of each timing separately
save('Plots.mat','PrettyPic');
hold off
%Reset the directory for the next analysis
cd(initDir)
