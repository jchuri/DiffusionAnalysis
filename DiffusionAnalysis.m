%% Intro and Change Log
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

% Change log (as of March 3, 2021):
    % Revised some variable names and changed some spacings
    % Preallocated some arrays that hold important data (times, variance_x/y, etc.)
    % Broke apart the operations within the big FOR loop over each timing folder
        % -"parseTXTFile.m" checks to see if each .txt file is acceptable to use
        % -"TwoDGaussianFitting.m" compiles the results of "parseTXTFile.m",
        % creates the averaged and corrected image, fits it with a 2D
        % Gaussian, then extracts the variance along x and y and the
        % coordinates of the cloud's center of mass

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
%% GUI for finding the right folder
% Start by figuring out where all of our folders are 
clear;
initDir = 'C:\Users\staronal\Desktop\Raw Data\Diffusion'; %IF YOU RUN THIS SCRIPT ON ANOTHER COMPUTER, CHANGE THIS PATH BEFOREHAND!
dataDir = uigetdir('Where are all the data folders?');
cd(dataDir);
%% Load the Data and Background (if there is one)
% Use this section of the code if dealing with one, constant background.
% Looks for the folder named 'Background' and pulls the fitted file from
% that
% In case it doesn't find it, this just stops it from crashing by giving it
% a bunch of zeros to work from.
% if ~exist('Background','dir')
%     disp('No background folder! Using zeros')
%     background = zeros(126,125);
% else
%     load('Background/AverageBackground.mat')
% end
folderList = dir('*ms');
% backgroundlist = dir('*bg');
%% Load in Multiple Backgrounds
% Use this section of the code if dealing with a changing background.
% Looks for folders with 'bg' at the end of their name and builds
% background images from the data they contain. Stores these background
% images in an array AvgBack with index i. Index 'i' should always match up
% with the index 'i' used in the next 'folderlist' for loop--that is, your 20ms
% background will always be subtracted from your 20ms data. AvgBack{i} is
% subtracted from AvgPic in the next for loop.

% for i = 1:length(backgroundlist)
%     cd(backgroundlist(i).name);
%     backgroundfilelist = dir('*.txt');
%     
%     for j = 1:length(backgroundfilelist)
%         backcell{j} = load(backgroundfilelist(j).name); 
%     end
%     
%     backimage = cat(3,backcell{1:end});
%     AvgBack{i} = mean(backimage,3);
%     cd('..')
% end
%%
% This is threshold for picking out good images
% You can pick the 'all' option to sum over every pixel
% The 'one' option just finds the brightest pixel and compares to that
%thresholdall = 5e+7;
thresholdOne = 2500;
%Preallocate arrays to hold the timing and variance info, averaged
%pictures, and centers of mass for each expansion time
times = zeros(1, length(folderList));
timesMS = zeros(1, length(folderList));
x_var = zeros(1, length(folderList));
y_var = zeros(1, length(folderList));
PrettyPic = cell(1, length(folderList));
x_cen = zeros(1, length(folderList));
y_cen = zeros(1, length(folderList));
% Add a path to the Diffusion folder, which was set as initDir.
% This is necessary in order to properly call the custom functions used to
% analyze the data (like parseTXTFile, etc.)
% There's probably a better way to do this whole analysis that would avoid
% this step...
addpath(genpath(initDir));

% Now we're going to go through the folders found in the directory
figure;
for i = 1:length(folderList);  
    % Extract the timing (in ms) from the folder name. It's saved in a cell.
    timeFromFolderName = regexp(folderList(i).name,'\d+\.?\d*|-\d+\.?\d*|\.?\d+|-\.?\d+','match');
    % Convert the string witin the cell into a number. This number will be
    % used when plotting stuff later.
    times(i) = str2num(timeFromFolderName{1}) * 0.001;
    timesMS(i) = str2num(timeFromFolderName{1});
    % Go into the timing folder and create a list of all the .txt files inside
    cd(folderList(i).name);
    fileList = dir('*.txt');
    % Create a cell array to hold the good data cells and a counter to track how many there are
    goodFiles = cell(1, length(fileList));
    counter = 1;
    % Go through the .txt files and check them for validity. The maximum
    % intensity should be above thresholdOne and the minimum intensity
    % should be positive or 0. See "parseTXTFile.m" for more details.
    for j = 1:length(fileList)
        if parseTXTFile(fileList(j).name, thresholdOne, j, folderList(i).name) == true
            goodFiles{counter} = load(fileList(j).name);
            counter = counter + 1;
        end
    end
    % Compile the data in the good .txt files to create an average and
    % intensity-corrected image of the cloud. Fit that with a 2D Gaussian,
    % then extract the variance along x and y and the coordinates of the
    % cloud's center of mass from the fitting. See "TwoDGaussianFitting.m"
    % for more details.
    [PrettyPic{i}, x_var(i), y_var(i), x_cen(i), y_cen(i)] = TwoDGaussianFitting(goodFiles);
    % Use the extracted variances to calculate the full-width half-max
    % (FWHM) along x and y
    % Full-width half-max of a Gaussian is based on standard deviation
    % Specifically, about 2.355 standard deviations away
    FWHM_X = 2.355*sqrt(x_var(i));
    FWHM_Y = 2.355*sqrt(y_var(i));
    % This is an X generated for plotting out the different FWHM
    X = linspace(-FWHM_X + x_cen(i), FWHM_X + x_cen(i), 1000);
    % There's a better way to do this, but it creates the top and bottom
    % halves of the ellipse of the FWHM...
    for j = 1:length(X)
        tophalf(j) = (FWHM_Y/FWHM_X)*sqrt(FWHM_X^2 - (X(j) - x_cen(i))^2) + y_cen(i);
        bothalf(j) =-(FWHM_Y/FWHM_X)*sqrt(FWHM_X^2 - (X(j) - x_cen(i))^2) + y_cen(i);
    end
    
    %disp(class(tophalf))
    %TopList(i) = tophalf
    %BottomList(i) = bothalf
    
    % ...And then plots them
    legendname = num2str(timesMS(i));
    v = timesMS(i)/max(timesMS);
    w = 1-v;

    %stores graphs for each expansion time in an array to be plotted later
    Xarray{i} = X;
    Toparray{i} = tophalf;
    Bottomarray{i} = bothalf;
    legendarray{i} = legendname;
    
    hold on
    Top = plot(X,tophalf,'Color',rand(1,3),'Linewidth',2,'DisplayName',legendname);
    Bottom = plot(X,bothalf,'Color',rand(1,3),'Linewidth',2,'DisplayName',legendname);
    CenterofMass = plot(x_cen(i),y_cen(i),'*','Color',rand(1,3));
    xlabel('Pixels (-Z)');
    ylabel('Pixels (-X)');
    title('Full-width half-max of cloud');
    
    %disp(class(legendname))
    legend show
    
    % Move back up in directory so we can go into the next folder
    cd('..')
end
hold off

% figure;
% Xmat = cell2mat(Xarray);
% Topmat = cell2mat(Toparray);
% Bottommat = cell2mat(Bottomarray);
% legendmat = cell2mat(legendarray);
% disp(legendmat)
% 
% plot(Xmat,Topmat)
% legend show
% hold on
% plot(Xmat,Bottommat)

PXtoCM = 2.85e-3;%3.6e-3%5.2e-3;%4e-3 gives X diffusion to match Schiavoni. Ethan's value is 5.2e-3.
ConversionFactor = PXtoCM*PXtoCM;

% Now we want to figure out that diffusion constant
% We're going to find it for x and y because it might be different
% First start by pulling out variables for a linear fit
linfitcoeffx = polyfit(times,x_var,1);
linfitcoeffy = polyfit(times,y_var,1);
% And this fit is to pixels, so convert it into cm
diffx = linfitcoeffx(1) * ConversionFactor * 0.5;
diffy = linfitcoeffy(1) * ConversionFactor * 0.5;

disp(['Z Diffusion = ',num2str(diffx)]);
disp(['X Diffusion = ',num2str(diffy)]);

%Adding in a feature that calculates the unitless factor (which I'm calling
%G) that Grynberg uses in most of his diffusion graphs.
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
