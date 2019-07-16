% The program uses the cropped images to calculate the frequency of
% rotation of beads. If not done already, run the GetCrop script
% to perform the cropping of the images to regions containing the rotating
% cell to be analyzed.

% The code is modified from the code written by Siyu He, who employed
% wavelet transform of the xy location of the cell centroid to measure the
% rotation frequency. That method worked well for high frequency, but not
% for low frequency. I kept the skeleton of that code, but changed the core
% to use the methods adapted from Lele, Hosu, and Berg 2013, who fitted
% ellipses on the cell body and made a fit on the cloud of dots
% representing the cell. Instead, I locate the center of rotation and find
% the vector pointing from the rotation center to the cell centroid. By
% measuring the angle between that vector and the x-axis, I am able to
% compute the frequency of rotation.


clc
close all
clearvars

data_path = uigetdir(pwd,'Select folder containig the files to analyze'); %select the folder wanted
cd(data_path); %set the path
Files=dir('*.avi'); %Collect a list of all the files to be analyzed in the current directory
NFiles=length(Files); %Measure how long the list of files is

for bb = 1:NFiles
%%
% Read the cropped region of the tif files in the folder, and proceed with
% mesurement of the centroid location, rotation center, angle, and
% frequency.

% Load the tif stack
filename= [Files(bb).name]; % get the name and path of the TIF file

 if exist(strcat(filename(1:end-4),'_centroid.mat'),'file')==0   % If the TIF stack was analyzed previously, it is not necessary to do it again
    
    disp('Reading Images...')   % display progress
    vid = VideoReader(filename);       % start a video object
    FrameRate = vid.FrameRate;
    NFrames = floor(vid.Duration*vid.FrameRate);
    T = NaN(1,NFrames);   % Initialize array for experiment time
    F = NaN(1,NFrames);   % Initialize array for calculated frequency
    
    %% Measure the location of the cenroid of the cell in each frame
    ImCenter = zeros(NFrames,2);    % Make a 2-column array of length NFrames
    % Loop through all frames and get the centriod of the cell in each
    % frame
    kk = 0;
    while hasFrame(vid)
        kk = kk+1;
        clc
        disp(['Processing, ' filename ', progress = ' strcat(num2str(100*vid.CurrentTime/vid.Duration, '%.f'),'%')]) % display progress
        Im = gpuArray(rgb2gray(readFrame(vid)));
        ImNorm = mat2gray(Im); % convert the RGB image to grayscale
        ImSmooth = imgaussfilt(ImNorm,2);   % smooth the grayscale image
        thresh = multithresh(gather(ImSmooth),6);   % find a high threshold to binarize the image
        imLabel=gpuArray(imbinarize(gather(ImSmooth),thresh(end))); % binarize the gray image
        %       imshow(ImMask)
        % imLabel = logical(ImMask);  % Convert to 0s and 1s: Jinming ? I
        % donot think this is useful
        L = bwlabel(imLabel);       % Label the binary image
        props = regionprops(imLabel,'Area','centroid'); % Extract the properties: Area and Centroid of each blob
        areas = [props.Area]; % unpack the areas
        centroids = cat(1, props.Centroid); % unpack the centroids
        % Pick the cell by using both the area and the centroid.
        judge_center=areas'.*(1./((centroids(:,1)-size(Im,1)/2).^2+(centroids(:,2)-size(Im,2)/2).^2));
        [~,max_area] = max(judge_center);   % get the position of cell
        ImCenter(kk,:) = centroids(max_area,:);
        T(kk) = vid.CurrentTime;
    end     % end of centroid measurement
    save(strcat(filename(1:end-4),'_centroid'),'ImCenter','T','FrameRate'); % save the centroid measurement in a file
else    % if the centroid measurement already exists, load directly
    disp([filename ' centroid measurement already exists. Loading from file..'])
    CENTROID=load(strcat(filename(1:end-4),'_centroid.mat'));
    ImCenter=CENTROID.ImCenter;
    T = CENTROID.T;
    FrameRate = CENTROID.FrameRate;
end

%% Frequency measurement 
filt_order = 100;
X = ImCenter(:,1);
Y = ImCenter(:,2);
WindowSize = 2*FrameRate;

XCenter = movmean(ImCenter(:,1),floor(FrameRate/2));
YCenter = movmean(ImCenter(:,2),floor(FrameRate/2));

% Measure the frequecy using the GetFreq function
Omega = GetFreq(X,Y,FrameRate,XCenter,YCenter);

% Store the result in appropriate vectors
F=medfilt1([NaN Omega'],filt_order,'omitnan','truncate');
plot(T,F,'-k');

% Make a time-freq plot for the whole data and save it
% ylim([min([0 min(F)]) 20]) % max expected freq = ~ 20. change this if needed
xlabel('Time (s)')
ylabel('Rotation Frequency (Hz)')
xlim([0 max(T)])
set(gca,'Fontsize',24)

set(gcf,'PaperUnits','centimeters')
xSize = 10; ySize = 5;
xLeft = 1; yTop = 1;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[100 200 xSize*100 ySize*100])
pause(0.5)
saveas(gcf,strcat(filename(1:end-4),'_Plot'),'png'); % save the plot as png

save(strcat(filename(1:end-4),'_Freq'),'T','F'); % save the frequency measurement in a file
close(gcf)
end
profile on
profile viewer
profsave
