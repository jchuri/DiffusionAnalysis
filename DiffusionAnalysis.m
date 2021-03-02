function DiffAnal = DiffAnal()

%% Intro and news
% A mature naming scheme for mature people

% This is the 2019 iteration of the old code "FWHMMulAdv.m"
% That code was started by Ethan and put together by Anthony
% Neither of which knew what they were doing

% New code is meant to be as comprehensive as possible, provide easy
% analysis

% Assuming that only occasionally you have bad data, it should get averaged
% out by the other 'good' data sets. In a best case scenario, you totally
% missed the cloud and just got the background, in which case it doesn't do
% anything bad at all, since all of these values will get scaled anyways

% Another big improvement is that the cloud doesn't have to be centered
% anymore! This used to be a VERY big issue. Since it works by fitting now,
% you can have off-center or even partially off-screen clouds and it should
% still fit pretty well

%% DATA REQUIREMENTS
% All data folders should be named in the scheme '**ms' for data and '**bg'
% for backgrounds.
% This program will find all folders with that name and pull them to use in
% data. Other things are allowed to be in the folder now

%% GUI for finding the right folder
% Start by figuring out where all of our folders are 

startdir = uigetdir('Where are all the data folders?');
cd(startdir);

%% Load in the background
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

folderlist = dir('*ms');
% backgroundlist = dir('*bg');
%%

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
thresholdone = 2500;

figure;
% Now we're going to go through the folders found in the directory
% Start at 3 because the first two are '.' and '..'
for i = 1:length(folderlist);  
    % Take the name of the folder, look at it and pull out a number
    % That should be the timing in ms, and it'll save to temp as a string
    % inside of a cell
    temp = regexp(folderlist(i).name,'\d+\.?\d*|-\d+\.?\d*|\.?\d+|-\.?\d+','match');
    % Call that cell, convert the string inside to a number so we can plot
    % it later
    times(i) = str2num(temp{1}) * 0.001;
    timesMS(i) = str2num(temp{1});
   
    % Now we're actually going to go into the folder...
    cd(folderlist(i).name);
    % ...And look at all the filse inside
    Filelist = dir('*.txt');
    % One by one, load them all into a cell
    j2 = 1;
    for j = 1:length(Filelist)
        % We're going to check how much light each imaged actually picked
        % up
        % If it's below the threshold, the image is just going to be
        % removed
        tempcell{j2} = load(Filelist(j).name);
        % This sums up every single pixel intensity to one value
        intensitysum = sum(sum(tempcell{j2}));
        % And then this will find the maximum intensity of any one bin AND
        % it's index
        [intensitymax,index] = max(tempcell{j2}(:));
        [intensitymin,index] = min(tempcell{j2}(:));
        % It turns out that index is worthless on its own for a 2D list
        % Instead, convert it to an x and a y coordinate
        [xguess,yguess] = ind2sub(size(tempcell{j2}),index);     
        
        % Are we above threshold? 
        %if intensitysum < thresholdall
        if intensitymax < thresholdone
            % If not, toss it
            tempcell{j2} = [];
            disp(['Tossed file number ',num2str(j),' from folder ',folderlist(i).name]);
        else
            % If yes, keep it, update j2
            j2 = j2 + 1;
        end
        
        % Are there negative numbers?
        if intensitymin < 0
            % If so, toss it
            tempcell{j2} = [];
            disp(['Just tossed file number ',num2str(j),' from folder ',folderlist(i).name]);
        else
            % If not, keep it, update j2
            j2 = j2 + 1;
        end
    end
    % Then take all of the cells and convert them into a big
    % three-dimensional array
    image = cat(3,tempcell{1:end});
    % Now we have everything in one, smash it all down to an average
    % picture
    AvgPic = mean(image,3);
    
    % Then pull away whatever the background picks up
    % Pull off the first column, it's just zeros
    % This is the "corrected picture" of the averaged clouds, taking away
    % whatever light we saw from the background
    
%     CorPic = AvgPic(2:end,1:end) - AvgBack{i}(2:end,1:end);
    CorPic = AvgPic(2:end,1:end);
    
    % Then scale it so most intense point is 1
    CorPic = CorPic./max(max((CorPic)));
    % PrettyPic is just a cell of the scaled and cleaned up array. We won't
    % plot it in this program because it would just be a lot of plots, but
    % it's saved so you can look at them later
    PrettyPic{i} = CorPic;
    % Now we want to fit a two-dimensional Gaussian to this data
    % This essentially turns xo, yo, and zo into these really big lists
    % that the fitting program will pull from in the next step
    [xo,yo,zo] = prepareSurfaceData(1:125,1:125,CorPic);
    % Then we actually fit them to a basic two-dimensional Gaussian
    % This is started with the assumption that the center of the cloud is
    % probably right around the middle. The starting variance is just kind
    % of a guess. If the fitting isn't working well, maybe try changing
    % these.
    % The other option for "Robust" is "LAR"
    % As I understand, use LAR if you don't have many outliers in your
    % data. Use Bisquare in case you might have some occasional "out there"
    % values.
    ft = fit([xo,yo],zo,'exp(-(((x-xcen)^2)/vx + ((y-ycen)^2)/vy))','StartPoint',[100,100,62,62],'Robust','Bisquare');
    % Order currently vx,vy,xcen,ycen
    % Ease of computation, both vx and vy are missing a 1/2 in the fitting
    % Then vx and vy are lists of the variance in x and y
    vx(i) = ft.vx./2;
    vy(i) = ft.vy./2;
    % While xcen and ycen are lists of the center of mass
    xcen(i)=ft.xcen;
    ycen(i)=ft.ycen;
    % Move back up in directory so we can go into the next folder
    cd('..')
    
    % Full-width half-max of a Gaussian is based on standard deviation
    % Specifically, about 2.355 standard deviations away
    FWHMX = 2.355*sqrt(vx(i));
    FWHMY = 2.355*sqrt(vy(i));
    % This is an X generated for plotting out the different FWHM
    X = linspace(-FWHMX+xcen(i),FWHMX+xcen(i),1000);
    % There's a better way to do this, but it creates the top and bottom
    % halves of the ellipse of the FWHM...
    for j = 1:length(X)
        tophalf(j) = (FWHMY/FWHMX)*sqrt(FWHMX^2-(X(j)-xcen(i))^2) + ycen(i);
        bothalf(j) =-(FWHMY/FWHMX)*sqrt(FWHMX^2-(X(j)-xcen(i))^2) + ycen(i);
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
    
    Top = plot(X,tophalf,'Color',rand(1,3),'Linewidth',2,'DisplayName',legendname);
    hold on
    Bottom = plot(X,bothalf,'Color',rand(1,3),'Linewidth',2,'DisplayName',legendname);
    hold on
    CenterofMass = plot(xcen(i),ycen(i),'*','Color',rand(1,3));
    hold on
    xlabel('Pixels (-Z)');
    ylabel('Pixels (-X)');
    title('Full-width half-max of cloud');
    
    %disp(class(legendname))
    legend show
    hold on;
end

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

PXtoCM = 2.85e-3;%3.6e-3%5.2e-3;%2.85e-3;%4e-3 gives X diffusion to match Schiavoni. Ethan's value is 5.2e-3.
ConversionFactor = PXtoCM*PXtoCM;

% Now we want to figure out that diffusion constant
% We're going to find it for x and y because it might be different
% First start by pulling out variables for a linear fit
linfitcoeffx = polyfit(times,vx,1);
linfitcoeffy = polyfit(times,vy,1);
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

% Once you have those variables for your line plug them into a function to
% get out a set of actualy fitted lines
vxfit = polyval(linfitcoeffx,times);
vyfit = polyval(linfitcoeffy,times);

% Make a new figure separate from the FWHM plot
figure;
% Plotting the actual data points with *
plot(timesMS,vx.*ConversionFactor,'*',timesMS,vy.*ConversionFactor,'*');
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

