%% visualize ID'ed particles to determine trimming parameters


% goal: measure cell size from a phase image (150x),
%       overlay measurement onto phase image,
%       visualize and repeat to determine desired values
%       (e.g. values that identify single cells)


% important terminology:
%
%       1. strain = the fixed cell sample used for imaging/labeling
%                   e.g. stationary G6 and exponential phase G6 are
%                   different fixed cell samples that may need different
%                   segmentation parameters
%
%       2. sample = the labeling condition in that experiment
%                   e.g. unlabeled and labeled G6 are different conditions
%
%     * i realize this is potentially confusing, but wish to remain consistent with
%       the terminology used in enterMetaData.m and metadata.mat


% strategy:
%
%   For each strain of interest, user is prompted to:
%       1. choose one test image to identify single cells based on width
%             - have image name ready to choose its index from displayed list
%       2. user is prompted to enter guesses of min and mix width
%       3. user then checks particle overlaps onto phase image to determine
%           if input widths are accurate
%       4. if not accurate, user can adjust guesses and go again
%       5. if acceptable, user can then choose to save data
%       6. this repeats for all unique strain samples (fixed cultures)
%
%       7. when finished, edit the "last updated" and "commit" fields
%          immediately below with date and a note about which experiment
%          was added to segdata.mat in this run 
%
%           example commit: "add 2021-06-15 experiment to segdata"


% last updated: kisa, 2021 June 28
% commit: Add 2021-06-25 experiment to segdata 

% ok, let's go!

%% Part 0. intitialize user specific file paths & experiment of interest

%  Note to user: do these items before running whos_a_cell.m
%
%       1. comment out "path2meta" line(s) not appropriate to your system
%       2. comment out "prepath" line(s) not appropriate to your system
%       3. input metadata index of experiment of interest

clc
clear

% 1. define path to metadata.mat and segdata.mat files
%path2meta = '/Users/jen/summero-strains'; % jen's Mac OS
path2meta = 'C:/Users/Kisa Naqvi/Documents/TropiniLab/summero-strains-master'; % kisa's PC


% 2. define prefix for path to image data
%prepath = '/Users/jen/Documents/TropiniLab/Molecular_tools/HiPR_fish/'; % jen's MacOS (non Kisa data)
%prepath = '/Users/jen/Documents/TropiniLab/Data/Kisa/';                %jen's MacOS (Kisa data)
prepath = 'C:/Users/Kisa Naqvi/Documents/TropiniLab/Data/';            % kisa's PC


% 3. load stored data and define experiment of interest
cd(path2meta)
load('metadata.mat')
load('segdata.mat')
tempdata = segdata;
index = 7; % index of experiment in metadata


%% Part 1. collect and store segmentation parameters for a new experiment

% 1. initialize experiment meta data
exp_metadata = metadata{index};
magn = exp_metadata.magnification;
exp_strains = exp_metadata.strains;
exp_samples = exp_metadata.samples;
experiment = exp_metadata.date;
exp_stage = exp_metadata.growthStage;


% 2. initialize image meta data
px_size = 11/magn; % 11 um pixels with experiment specific magnification
prefix = 'img_';
suffix = '_position000_time000000000_z000.tif';
name_phase = strcat(prefix,'channel000',suffix);
clear prefix suffix

% 3. access image data
data_folder = strcat(prepath,experiment);
cd(data_folder)


% 4. define sample image(s) per unique strain in experiment
%
%    note: this image should contain both single cells and clumps so that
%          we can determine how well the segmentation parameters classify
%          particles into the categories of "noise", "single cells" and "clumps"

% 4a. determine how many unique strains in experiment
expstrains_string = cellfun(@unique,exp_strains); % convert cell array to string
strain_groups = findgroups(expstrains_string); % classify strings by text
unique_groups = unique(strain_groups); % list unique classes
unique_strains = {};
for ug = 1:length(unique_groups) % loop through unique classes and store strain info
    
    group = unique_groups(ug);
    group_idx = find(strain_groups==group,1,'first');
    unique_strains{ug} = exp_strains{group_idx};
end
clear expstrains_string group group_idx ug

% 4b. for each unique strain, input index chosen image in list of image names
for strain = 1:length(unique_strains)
    
    % 4b,i. make a directory of images with current strain
    current_strain = unique_strains{strain};
    samples_with_current_strain = exp_samples(strain_groups==strain);
    
    sDirectory = [];
    for ss = 1:length(samples_with_current_strain)
        current_sample = samples_with_current_strain{ss};
        current_dir = dir(strcat(current_sample,'_*'));
        sDirectory = [sDirectory;current_dir];   
    end
    disp(strcat('Image names of current strain sample:',current_strain, ', ',exp_stage))
    names = {sDirectory.name} % display all names for images of current strain
    
    % 4b,ii. prompt user for image to analyse
    disp(strcat('Current strain sample = ',current_strain, ', ',exp_stage))
    prompt = 'Enter one image index in "names" as a double: ';
    img2test = input(prompt);
    strain_img = names{img2test}; % chosen image
    testset{strain} = strain_img;
end
clear strain_groups ss prompt
clear current_strain current_sample strain_img



% 5. loop through each strain in experiment to find segmentation parameters
%    and store parameters in segdata.mat
uniqcounter = 1;
while uniqcounter <= length(testset)
    
    unisample = testset{uniqcounter};
    unistr = unique_strains{uniqcounter};
    
    % 5a. display phase image & mask
    cd(data_folder)
    cd(strcat(unisample,'/Default'))
    img_phase = imread(name_phase); % read phase image
    figure(1) % display phase image
    imshow(img_phase, 'DisplayRange',[1000 3000]); %lowering right # increases num sat'd pxls
    title(strcat(unistr,'-phase'))
    
    phase_smoothed = imgaussfilt(img_phase,0.8);
    bw = edge(phase_smoothed,'sobel');
    se = strel('disk',1); % structuring element = disk of radius x pixels
    bw_dil = imdilate(bw, se); % dilate
    bw_fill = imfill(bw_dil, 'holes'); % fill
    bw_final = imerode(imerode(bw_fill,se),se); % erode twice
    figure(2)
    imshow(bw_final)
    title(strcat(unistr,'-mask'))
    figure(3) % display phase image
    imshow(img_phase, 'DisplayRange',[1000 3000]); %lowering right # increases num sat'd pxls
    title(strcat(unistr,'-phase-overlay'))
    
    % 5b. measure particle parameters
    cc = bwconncomp(bw_final);
    stats = regionprops(cc,'Area','Centroid','MajorAxisLength','MinorAxisLength','Eccentricity','Orientation');
    majorAxes = extractfield(stats,'MajorAxisLength')'.*px_size;
    minorAxes = extractfield(stats,'MinorAxisLength')'.*px_size;
    angles = extractfield(stats,'Orientation')';
    clear bw bw_dil bw_fill se phase_smoothed
    
    % 5c. prompt user to input min and max width thresholds for single cells
    disp(strcat('Current strain sample = ',unistr,', ',exp_stage))
    prompt_min = 'Enter minimum width (um) to define a single cell as a double: ';
    minWidth = input(prompt_min);
    prompt_max = 'Enter maximum width (um) to define a single cell as a double: ';
    maxWidth = input(prompt_max);
    clear prompt_min prompt_max
    
    % 5d. overlay particle outlines with blue indicating particles within width thresholds
    for pp = 1:length(stats)
        
        centroid = stats(pp).Centroid.*px_size;
        [x_rotated, y_rotated] = drawEllipse2(pp, majorAxes, minorAxes, centroid(1), centroid(2), angles, px_size);
        lineVal = 1;
        
        if minorAxes(pp) < minWidth  % noise threshold
            color = rgb('HotPink');
            color_text = rgb('DarkMagenta');
        elseif minorAxes(pp) > maxWidth % clump threshold
            color = rgb('DeepPink');
            color_text = rgb('DarkMagenta');
        else
            color = rgb('DarkTurquoise'); % single cells
            color_text = rgb('MidnightBlue');
        end
        
        figure(3)
        hold on
        plot(x_rotated,y_rotated,'Color',color,'lineWidth',lineVal)
        text((centroid(1)+1)/px_size, (centroid(2)+1)/px_size, num2str(pp),'Color',color_text,'FontSize',10)
        
    end
    
    
    % 5d. prompt user to accept or edit thresholds
    %     i) if accept, store thresholds and proceed to next unique strain
    %    ii) if edit, repeat with current strain until satisfied
    disp(strcat('Current strain sample = ',unistr,', ',exp_stage))
    prompt_accept = 'Satisfied? Input "Yes" (stores seg data) or "No" (repeats viz) as a string: ';
    isSatified = input(prompt_accept);
    if strcmp(isSatified, "Yes") == 1
        disp('Saving min and max widths to segdata.mat')
        
        % saving data to segdata.mat
        % a. determine which row to add new data 
        uniqcounter = uniqcounter + 1;
        newrow = length(tempdata) + 1;
        newdata(1).sample_strain = unistr;
        newdata(1).sample_stage = exp_stage;
        newdata(1).tested_img = unisample;
        newdata(1).test_experiment = experiment;
        newdata(1).minWidth = minWidth;
        newdata(1).maxWidth = maxWidth;
        
        prompt_sampledate = 'Input date sample was prepared/collected as a string: ';
        sample_date = input(prompt_sampledate);
        newdata(1).sample_date = sample_date;
        
        % b. store data into temporary structure
        tempdata{newrow,1} = newdata;
        
    elseif strcmp(isSatified, "No") == 1
        disp('Thanks for being thorough! Please try again, adjusting min and max widths')
    else
        disp('Error: must input either "Yes" or "No", please restart :P')
    end
    close(figure(1),figure(2),figure(3))
    
end

tempdata
segdata

%% Part 2. if happy, save data!
cd(path2meta)
segdata = tempdata;
save('segdata.mat','segdata')


%% 381 - Lactobacillus plantarum ATCC BAA-793 (Firmicute)

% clc
% clear

% 0. initialize experiment meta data
% px_size = 11/150; % 11 um pixels with 150x magnification
% experiment = '2021-03-31';
% data_folder = strcat('/Users/jen/Documents/TropiniLab/Molecular_tools/HiPR_fish/',experiment);
% cd(data_folder)
sample = 's1';

% prefix = 'img_';
% suffix = '_position000_time000000000_z000.tif';
% name_phase = strcat(prefix,'channel000',suffix);
% clear prefix suffix
    

% % 0. manually choose one sample image (contains single cells and clumps)
% cd(data_folder)
% sDirectory = dir(strcat(sample,'_*'));
% names = {sDirectory.name};
% img_381 = names{2}; % chosen image


% 1. load one 150x phase image
% cd(data_folder)
% cd(strcat(img_381,'/Default'))
% img_phase = imread(name_phase); % read phase image
% figure(1) % display phase image
% imshow(img_phase, 'DisplayRange',[500 2500]); %lowering right # increases num sat'd pxls
% title('381 Lactobacillus plantarum (150x), phase')
% 
% figure(3)
% imshow(img_phase, 'DisplayRange',[500 2500]);
% title('381 Lactobacillus plantarum (150x), phase')

% 2. create mask using parameters in current data matrix
%    see segmentIntensity_GFP.m for full comments / annotation
% phase_smoothed = imgaussfilt(img_phase,0.8);
% bw = edge(phase_smoothed,'sobel');
% se = strel('disk',1); % structuring element = disk of radius x pixels
% bw_dil = imdilate(bw, se); % dilate
% bw_fill = imfill(bw_dil, 'holes'); % fill
% bw_final = imerode(imerode(bw_fill,se),se); % erode twice
% figure(2)
% imshow(bw_final)
% title('381 Lactobacillus plantarum (150x), mask')


% 3. measure particle features
% cc = bwconncomp(bw_final);
% stats = regionprops(cc,'Area','Centroid','MajorAxisLength','MinorAxisLength','Eccentricity','Orientation');
% clear bw bw_dil bw_fill se phase_smoothed
        

% 4. overlay outlines of all tracked particles
% majorAxes = extractfield(stats,'MajorAxisLength')'.*px_size;
% minorAxes = extractfield(stats,'MinorAxisLength')'.*px_size;
% angles = extractfield(stats,'Orientation')';

% 5. vary min and max measurements to determine...
%      i) below what threshold is background noise?
%     ii) within what thresholds is a single-cell?
%    iii) above what threshold is a clump of cells?
% for pp = 1:length(stats)
%     
%     centroid = stats(pp).Centroid.*px_size;
%     [x_rotated, y_rotated] = drawEllipse2(pp, majorAxes, minorAxes, centroid(1), centroid(2), angles, px_size);
%     lineVal = 1;
%     
%     if minorAxes(pp) < 1  % noise threshold
%         color = rgb('HotPink');
%         color_text = rgb('DarkMagenta');
%     elseif minorAxes(pp) > 1.6 % clump threshold
%         color = rgb('DeepPink');
%         color_text = rgb('DarkMagenta');
%     else
%         color = rgb('DarkTurquoise'); % single cells
%         color_text = rgb('MidnightBlue');
%     end
% 
%     figure(3)
%     hold on
%     plot(x_rotated,y_rotated,'Color',color,'lineWidth',lineVal)
%     text((centroid(1)+1)/px_size, (centroid(2)+1)/px_size, num2str(pp),'Color',color_text,'FontSize',10)
%     
% end

%% 505 - Enterobacter cancerogenus ATCC 35316 (Proteobacteria)

clc
clear

% 0. initialize experiment meta data
px_size = 11/150; % 11 um pixels with 150x magnification
experiment = '2021-03-31';
data_folder = strcat('/Users/jen/Documents/TropiniLab/Molecular_tools/HiPR_fish/',experiment);
cd(data_folder)
sample = 's5'; % fixed, no encoding probe

prefix = 'img_';
suffix = '_position000_time000000000_z000.tif';
name_phase = strcat(prefix,'channel000',suffix);
clear prefix suffix
    

% 0. manually choose one sample image (contains single cells and clumps)
cd(data_folder)
sDirectory = dir(strcat(sample,'_*'));
names = {sDirectory.name};
% for stk = 1:length(names)
%     
%     cd(data_folder)
%     current_stack = names{stk};
%     cd(strcat(current_stack,'/Default'))
%     
%     img_phase = imread(name_phase);
%   
%     figure(stk)
%     imshow(img_phase, 'DisplayRange',[500 2500]); %lowering right # increases num sat'd pxls
% end
img_505 = names{9}; % chosen image


% 1. load one 150x phase image
cd(data_folder)
cd(strcat(img_505,'/Default'))
img_phase = imread(name_phase); % read phase image
figure(1) % display phase image
imshow(img_phase, 'DisplayRange',[500 2500]); %lowering right # increases num sat'd pxls
title('505 Enterobacter cancerogenus (150x), phase')

figure(3) % this one will get overlaid with particle outlines
imshow(img_phase, 'DisplayRange',[500 2500]);
title('505 Enterobacter cancerogenus (150x), phase')


% 2. create mask using parameters in current data matrix
%    see segmentIntensity_GFP.m for full comments / annotation
phase_smoothed = imgaussfilt(img_phase,0.8);
bw = edge(phase_smoothed,'sobel');
se = strel('disk',1); % structuring element = disk of radius x pixels
bw_dil = imdilate(bw, se); % dilate
bw_fill = imfill(bw_dil, 'holes'); % fill
bw_final = imerode(imerode(bw_fill,se),se); % erode twice
figure(2)
imshow(bw_final)
title('505 Enterobacter cancerogenus (150x), mask')


% 3. measure particle features
cc = bwconncomp(bw_final);
stats = regionprops(cc,'Area','Centroid','MajorAxisLength','MinorAxisLength','Eccentricity','Orientation');
clear bw bw_dil bw_fill se phase_smoothed
        

% 4. overlay outlines of all tracked particles
majorAxes = extractfield(stats,'MajorAxisLength')'.*px_size;
minorAxes = extractfield(stats,'MinorAxisLength')'.*px_size;
angles = extractfield(stats,'Orientation')';

% 5. vary min and max measurements to determine...
%      i) below what threshold is background noise?
%     ii) within what thresholds is a single-cell?
%    iii) above what threshold is a clump of cells?
for pp = 1:length(stats)
    
    centroid = stats(pp).Centroid.*px_size;
    [x_rotated, y_rotated] = drawEllipse2(pp, majorAxes, minorAxes, centroid(1), centroid(2), angles, px_size);
    lineVal = 1;
    
    if minorAxes(pp) < 0.7  % noise threshold
        color = rgb('HotPink');
        color_text = rgb('DarkMagenta');
    elseif minorAxes(pp) > 1 % clump threshold
        color = rgb('DeepPink');
        color_text = rgb('DarkMagenta');
    else
        color = rgb('DarkTurquoise'); % single cells
        color_text = rgb('MidnightBlue');
    end

    figure(3)
    hold on
    plot(x_rotated,y_rotated,'Color',color,'lineWidth',lineVal)
    text((centroid(1)+1)/px_size, (centroid(2)+1)/px_size, num2str(pp),'Color',color_text,'FontSize',10)
    
end

%% 488 - Bacteroides fragilis NCTC 9343

clc
clear

% 0. initialize experiment meta data
px_size = 11/150; % 11 um pixels with 150x magnification
experiment = '2021-03-31';
data_folder = strcat('/Users/jen/Documents/TropiniLab/Molecular_tools/HiPR_fish/',experiment);
cd(data_folder)
sample = 's12'; % fixed, no encoding probe

prefix = 'img_';
suffix = '_position000_time000000000_z000.tif';
name_phase = strcat(prefix,'channel000',suffix);
clear prefix suffix
    

% 0. manually choose one sample image (contains single cells and clumps)
cd(data_folder)
sDirectory = dir(strcat(sample,'_*'));
names = {sDirectory.name};
% for stk = 1:length(names)
%     
%     cd(data_folder)
%     current_stack = names{stk};
%     cd(strcat(current_stack,'/Default'))
%     
%     img_phase = imread(name_phase);
%   
%     figure(stk)
%     imshow(img_phase, 'DisplayRange',[500 2500]); %lowering right # increases num sat'd pxls
% end
img_488 = names{10}; % chosen image


% 1. load one 150x phase image
cd(data_folder)
cd(strcat(img_488,'/Default'))
img_phase = imread(name_phase); % read phase image
figure(1) % display phase image
imshow(img_phase, 'DisplayRange',[500 2500]); %lowering right # increases num sat'd pxls
title('488 Bacteroides fragilis (150x), phase')

figure(3) % this one will get overlaid with particle outlines
imshow(img_phase, 'DisplayRange',[500 2500]);
title('488 Bacteroides fragilis (150x), phase')


% 2. create mask using parameters in current data matrix
%    see segmentIntensity_GFP.m for full comments / annotation
phase_smoothed = imgaussfilt(img_phase,0.8);
bw = edge(phase_smoothed,'sobel');
se = strel('disk',1); % structuring element = disk of radius x pixels
bw_dil = imdilate(bw, se); % dilate
bw_fill = imfill(bw_dil, 'holes'); % fill
bw_final = imerode(imerode(bw_fill,se),se); % erode twice
figure(2)
imshow(bw_final)
title('488 Bacteroides fragilis (150x), mask')


% 3. measure particle features
cc = bwconncomp(bw_final);
stats = regionprops(cc,'Area','Centroid','MajorAxisLength','MinorAxisLength','Eccentricity','Orientation');
clear bw bw_dil bw_fill se phase_smoothed
        

% 4. overlay outlines of all tracked particles
majorAxes = extractfield(stats,'MajorAxisLength')'.*px_size;
minorAxes = extractfield(stats,'MinorAxisLength')'.*px_size;
angles = extractfield(stats,'Orientation')';

% 5. vary min and max measurements to determine...
%      i) below what threshold is background noise?
%     ii) within what thresholds is a single-cell?
%    iii) above what threshold is a clump of cells?
for pp = 1:length(stats)
    
    centroid = stats(pp).Centroid.*px_size;
    [x_rotated, y_rotated] = drawEllipse2(pp, majorAxes, minorAxes, centroid(1), centroid(2), angles, px_size);
    lineVal = 1;
    
    if minorAxes(pp) < 0.6  % noise threshold
        color = rgb('HotPink');
        color_text = rgb('DarkMagenta');
    elseif minorAxes(pp) > 1 % clump threshold
        color = rgb('DeepPink');
        color_text = rgb('DarkMagenta');
    else
        color = rgb('DarkTurquoise'); % single cells
        color_text = rgb('MidnightBlue');
    end

    figure(3)
    hold on
    plot(x_rotated,y_rotated,'Color',color,'lineWidth',lineVal)
    text((centroid(1)+1)/px_size, (centroid(2)+1)/px_size, num2str(pp),'Color',color_text,'FontSize',10)
    
end

%% BW25113 - E. coli (exponential phase)

clc
clear

% 0. initialize experiment meta data
px_size = 11/150; % 11 um pixels with 150x magnification
experiment = '2021-04-12';
data_folder = strcat('/Users/jen/Documents/TropiniLab/Molecular_tools/HiPR_fish/',experiment);
cd(data_folder)
sample = '1_no'; % fixed, no encoding probe

prefix = 'img_';
suffix = '_position000_time000000000_z000.tif';
name_phase = strcat(prefix,'channel000',suffix);
clear prefix suffix
    

% 0. manually choose one sample image (contains single cells and clumps)
cd(data_folder)
sDirectory = dir(strcat(sample,'_*'));
names = {sDirectory.name};
img_bw25113 = names{6}; % chosen image


% 1. load one 150x phase image
cd(data_folder)
cd(strcat(img_bw25113,'/Default'))
img_phase = imread(name_phase); % read phase image
figure(1) % display phase image
imshow(img_phase, 'DisplayRange',[500 2500]); %lowering right # increases num sat'd pxls
title('BW25113 E. coli exponential (150x), phase')

figure(3) % this one will get overlaid with particle outlines
imshow(img_phase, 'DisplayRange',[500 2500]);
title('BW25113 E. coli exponential (150x), phase')


% 2. create mask using parameters in current data matrix
%    see segmentIntensity_GFP.m for full comments / annotation
phase_smoothed = imgaussfilt(img_phase,0.8);
bw = edge(phase_smoothed,'sobel');
se = strel('disk',1); % structuring element = disk of radius x pixels
bw_dil = imdilate(bw, se); % dilate
bw_fill = imfill(bw_dil, 'holes'); % fill
bw_final = imerode(imerode(bw_fill,se),se); % erode twice
figure(2)
imshow(bw_final)
title('BW25113 E. coli exponential (150x), mask')


% 3. measure particle features
cc = bwconncomp(bw_final);
stats = regionprops(cc,'Area','Centroid','MajorAxisLength','MinorAxisLength','Eccentricity','Orientation');
clear bw bw_dil bw_fill se phase_smoothed
        

% 4. overlay outlines of all tracked particles
majorAxes = extractfield(stats,'MajorAxisLength')'.*px_size;
minorAxes = extractfield(stats,'MinorAxisLength')'.*px_size;
angles = extractfield(stats,'Orientation')';

% 5. vary min and max measurements to determine...
%      i) below what threshold is background noise?
%     ii) within what thresholds is a single-cell?
%    iii) above what threshold is a clump of cells?
for pp = 1:length(stats)
    
    centroid = stats(pp).Centroid.*px_size;
    [x_rotated, y_rotated] = drawEllipse2(pp, majorAxes, minorAxes, centroid(1), centroid(2), angles, px_size);
    lineVal = 1;
    
    if minorAxes(pp) < 1.3  % noise threshold
        color = rgb('HotPink');
        color_text = rgb('DarkMagenta');
    elseif minorAxes(pp) > 2.4 % clump threshold
        color = rgb('DeepPink');
        color_text = rgb('DarkMagenta');
    else
        color = rgb('DarkTurquoise'); % single cells
        color_text = rgb('MidnightBlue');
    end

    figure(3)
    hold on
    plot(x_rotated,y_rotated,'Color',color,'lineWidth',lineVal)
    text((centroid(1)+1)/px_size, (centroid(2)+1)/px_size, num2str(pp),'Color',color_text,'FontSize',10)
    
end

%% BW25113 - E. coli (stationary phase)

clc
clear

% 0. initialize experiment meta data
px_size = 11/150; % 11 um pixels with 150x magnification
experiment = '2021-04-12';
data_folder = strcat('/Users/jen/Documents/TropiniLab/Molecular_tools/HiPR_fish/',experiment);
cd(data_folder)
sample = '4_TP001'; % fixed, no encoding probe

prefix = 'img_';
suffix = '_position000_time000000000_z000.tif';
name_phase = strcat(prefix,'channel000',suffix);
clear prefix suffix
    

% 0. manually choose one sample image (contains single cells and clumps)
cd(data_folder)
sDirectory = dir(strcat(sample,'_*'));
names = {sDirectory.name};
img_bw25113 = names{2}; % chosen image


% 1. load one 150x phase image
cd(data_folder)
cd(strcat(img_bw25113,'/Default'))
img_phase = imread(name_phase); % read phase image
figure(1) % display phase image
imshow(img_phase, 'DisplayRange',[500 2500]); %lowering right # increases num sat'd pxls
title('BW25113 E. coli exponential (150x), phase')

figure(3) % this one will get overlaid with particle outlines
imshow(img_phase, 'DisplayRange',[500 2500]);
title('BW25113 E. coli exponential (150x), phase')


% 2. create mask using parameters in current data matrix
%    see segmentIntensity_GFP.m for full comments / annotation
phase_smoothed = imgaussfilt(img_phase,0.8);
bw = edge(phase_smoothed,'sobel');
se = strel('disk',1); % structuring element = disk of radius x pixels
bw_dil = imdilate(bw, se); % dilate
bw_fill = imfill(bw_dil, 'holes'); % fill
bw_final = imerode(imerode(bw_fill,se),se); % erode twice
figure(2)
imshow(bw_final)
title('BW25113 E. coli exponential (150x), mask')


% 3. measure particle features
cc = bwconncomp(bw_final);
stats = regionprops(cc,'Area','Centroid','MajorAxisLength','MinorAxisLength','Eccentricity','Orientation');
clear bw bw_dil bw_fill se phase_smoothed
        

% 4. overlay outlines of all tracked particles
majorAxes = extractfield(stats,'MajorAxisLength')'.*px_size;
minorAxes = extractfield(stats,'MinorAxisLength')'.*px_size;
angles = extractfield(stats,'Orientation')';

% 5. vary min and max measurements to determine...
%      i) below what threshold is background noise?
%     ii) within what thresholds is a single-cell?
%    iii) above what threshold is a clump of cells?
for pp = 1:length(stats)
    
    centroid = stats(pp).Centroid.*px_size;
    [x_rotated, y_rotated] = drawEllipse2(pp, majorAxes, minorAxes, centroid(1), centroid(2), angles, px_size);
    lineVal = 1;
    
    if minorAxes(pp) < 1  % noise threshold
        color = rgb('HotPink');
        color_text = rgb('DarkMagenta');
    elseif minorAxes(pp) > 1.9 % clump threshold
        color = rgb('DeepPink');
        color_text = rgb('DarkMagenta');
    else
        color = rgb('DarkTurquoise'); % single cells
        color_text = rgb('MidnightBlue');
    end

    figure(3)
    hold on
    plot(x_rotated,y_rotated,'Color',color,'lineWidth',lineVal)
    text((centroid(1)+1)/px_size, (centroid(2)+1)/px_size, num2str(pp),'Color',color_text,'FontSize',10)
    
end

%% 381 - L. plantarum ATCC BAA-793 (exponential phase)

clc
clear

% 0. initialize experiment meta data
px_size = 11/100; % 11 um pixels with 100x magnification
experiment = '2021-04-20';
data_folder = strcat('/Users/jen/Documents/TropiniLab/Molecular_tools/HiPR_fish/',experiment);
cd(data_folder)
sample = '1i'; % fixed, no encoding probe

prefix = 'img_';
suffix = '_position000_time000000000_z000.tif';
name_phase = strcat(prefix,'channel000',suffix);
clear prefix suffix
    

% 0. manually choose one sample image (contains single cells and clumps)
cd(data_folder)
sDirectory = dir(strcat(sample,'_*'));
names = {sDirectory.name};
img_381_exp = names{1}; % chosen image


% 1. load one 150x phase image
cd(data_folder)
cd(strcat(img_381_exp,'/Default'))
img_phase = imread(name_phase); % read phase image
figure(1) % display phase image
imshow(img_phase, 'DisplayRange',[2000 5000]); %lowering right # increases num sat'd pxls
title('381 L. plantarum exponential (100x), phase')

figure(3) % this one will get overlaid with particle outlines
imshow(img_phase, 'DisplayRange',[2000 5000]);
title('381 L. plantarum exponential (100x), phase')


% 2. create mask using parameters in current data matrix
%    see segmentIntensity_GFP.m for full comments / annotation
phase_smoothed = imgaussfilt(img_phase,0.8);
bw = edge(phase_smoothed,'sobel');
se = strel('disk',1); % structuring element = disk of radius x pixels
bw_dil = imdilate(bw, se); % dilate
bw_fill = imfill(bw_dil, 'holes'); % fill
bw_final = imerode(imerode(bw_fill,se),se); % erode twice
figure(2)
imshow(bw_final)
title('381 L. plantarum exponential (100x), mask')


% 3. measure particle features
cc = bwconncomp(bw_final);
stats = regionprops(cc,'Area','Centroid','MajorAxisLength','MinorAxisLength','Eccentricity','Orientation');
clear bw bw_dil bw_fill se phase_smoothed
        

% 4. overlay outlines of all tracked particles
majorAxes = extractfield(stats,'MajorAxisLength')'.*px_size;
minorAxes = extractfield(stats,'MinorAxisLength')'.*px_size;
angles = extractfield(stats,'Orientation')';

% 5. vary min and max measurements to determine...
%      i) below what threshold is background noise?
%     ii) within what thresholds is a single-cell?
%    iii) above what threshold is a clump of cells?
for pp = 1:length(stats)
    
    centroid = stats(pp).Centroid.*px_size;
    [x_rotated, y_rotated] = drawEllipse2(pp, majorAxes, minorAxes, centroid(1), centroid(2), angles, px_size);
    lineVal = 1;
    
    if minorAxes(pp) < 1  % noise threshold
        color = rgb('HotPink');
        color_text = rgb('DarkMagenta');
    elseif minorAxes(pp) > 1.8 % clump threshold
        color = rgb('DeepPink');
        color_text = rgb('DarkMagenta');
    else
        color = rgb('DarkTurquoise'); % single cells
        color_text = rgb('MidnightBlue');
    end

    figure(3)
    hold on
    plot(x_rotated,y_rotated,'Color',color,'lineWidth',lineVal)
    text((centroid(1)+1)/px_size, (centroid(2)+1)/px_size, num2str(pp),'Color',color_text,'FontSize',10)
    
end

%% H07 - S24-7 G6 (exponential phase)

clc
clear

% 0. initialize experiment meta data
px_size = 11/150; % 11 um pixels with 150x magnification
experiment = '2021-06-04';
data_folder = strcat('/Users/jen/Documents/TropiniLab/Data/Kisa/',experiment);
cd(data_folder)
sample = 'G2'; % fixed, no encoding probe

prefix = 'img_';
suffix = '_position000_time000000000_z000.tif';
name_phase = strcat(prefix,'channel000',suffix);
clear prefix suffix
    

% 0. manually choose one sample image (contains single cells and clumps)
cd(data_folder)
sDirectory = dir(strcat(sample,'_*'));
names = {sDirectory.name};
img_h07_exp = names{4}; % chosen image


% 1. load one 150x phase image
cd(data_folder)
cd(strcat(img_h07_exp,'/Default'))
img_phase = imread(name_phase); % read phase image
figure(1) % display phase image
imshow(img_phase, 'DisplayRange',[1000 3000]); %lowering right # increases num sat'd pxls
title('h07 S24-7 G6 (150x), phase')


% 2. create mask using parameters in current data matrix
%    see segmentIntensity_GFP.m for full comments / annotation
phase_smoothed = imgaussfilt(img_phase,0.8);
bw = edge(phase_smoothed,'sobel');
se = strel('disk',1); % structuring element = disk of radius x pixels
bw_dil = imdilate(bw, se); % dilate
bw_fill = imfill(bw_dil, 'holes'); % fill
bw_final = imerode(imerode(bw_fill,se),se); % erode twice
figure(2)
imshow(bw_final)
title('h07 S24-7 G6 (150x), mask')


% 3. measure particle features
cc = bwconncomp(bw_final);
stats = regionprops(cc,'Area','Centroid','MajorAxisLength','MinorAxisLength','Eccentricity','Orientation');
clear bw bw_dil bw_fill se phase_smoothed
        

% 4. overlay outlines of all tracked particles
majorAxes = extractfield(stats,'MajorAxisLength')'.*px_size;
minorAxes = extractfield(stats,'MinorAxisLength')'.*px_size;
angles = extractfield(stats,'Orientation')';

% 5. vary min and max measurements to determine...
%      i) below what threshold is background noise?
%     ii) within what thresholds is a single-cell?
%    iii) above what threshold is a clump of cells?

figure(3) % this one will get overlaid with particle outlines
imshow(img_phase, 'DisplayRange',[1000 3000]);
title('h07 S24-7 G6 (150x), phase')

for pp = 1:length(stats)
    
    centroid = stats(pp).Centroid.*px_size;
    [x_rotated, y_rotated] = drawEllipse2(pp, majorAxes, minorAxes, centroid(1), centroid(2), angles, px_size);
    lineVal = 1;
    
    if minorAxes(pp) < 0.6  % noise threshold
        color = rgb('HotPink');
        color_text = rgb('DarkMagenta');
    elseif minorAxes(pp) > 0.9 % clump threshold
        color = rgb('DeepPink');
        color_text = rgb('DarkMagenta');
    else
        color = rgb('DarkTurquoise'); % single cells
        color_text = rgb('MidnightBlue');
    end

    figure(3)
    hold on
    plot(x_rotated,y_rotated,'Color',color,'lineWidth',lineVal)
    text((centroid(1)+1)/px_size, (centroid(2)+1)/px_size, num2str(pp),'Color',color_text,'FontSize',10)
    
end

%% H01 - S24-7 NM-65 (exponential phase)


clc
clear

% 0. initialize experiment meta data
px_size = 11/150; % 11 um pixels with 150x magnification
experiment = '2021-06-04';
data_folder = strcat('C:/Users/Kisa Naqvi/Documents/TropiniLab/Data/',experiment);
cd(data_folder)
sample = 'N1'; % fixed with NM65 encoding probe

prefix = 'img_';
suffix = '_position000_time000000000_z000.tif';
name_phase = strcat(prefix,'channel000',suffix);
clear prefix suffix
    

% 0. manually choose one sample image (contains single cells and clumps)
cd(data_folder)
sDirectory = dir(strcat(sample,'_*'));
names = {sDirectory.name};
img_h01_exp = names{2}; % chosen image


% 1. load one 150x phase image
cd(data_folder)
cd(strcat(img_h01_exp,'/Default'))
img_phase = imread(name_phase); % read phase image
figure(1) % display phase image
imshow(img_phase, 'DisplayRange',[1000 3000]); %lowering right # increases num sat'd pxls
title('h01 S24-7 NM65 (150x), phase')


% 2. create mask using parameters in current data matrix
%    see segmentIntensity_GFP.m for full comments / annotation
phase_smoothed = imgaussfilt(img_phase,0.8);
bw = edge(phase_smoothed,'sobel');
se = strel('disk',1); % structuring element = disk of radius x pixels
bw_dil = imdilate(bw, se); % dilate
bw_fill = imfill(bw_dil, 'holes'); % fill
bw_final = imerode(imerode(bw_fill,se),se); % erode twice
figure(2)
imshow(bw_final)
title('h01 S24-7 NM65 (150x), mask')


% 3. measure particle features
cc = bwconncomp(bw_final);
stats = regionprops(cc,'Area','Centroid','MajorAxisLength','MinorAxisLength','Eccentricity','Orientation');
clear bw bw_dil bw_fill se phase_smoothed
        

% 4. overlay outlines of all tracked particles
majorAxes = extractfield(stats,'MajorAxisLength')'.*px_size;
minorAxes = extractfield(stats,'MinorAxisLength')'.*px_size;
angles = extractfield(stats,'Orientation')';

% 5. vary min and max measurements to determine...
%      i) below what threshold is background noise?
%     ii) within what thresholds is a single-cell?
%    iii) above what threshold is a clump of cells?

figure(3) % this one will get overlaid with particle outlines
imshow(img_phase, 'DisplayRange',[1000 3000]);
title('h01 S24-7 NM65 (150x), phase')

for pp = 1:length(stats)
    
    centroid = stats(pp).Centroid.*px_size;
    [x_rotated, y_rotated] = drawEllipse2(pp, majorAxes, minorAxes, centroid(1), centroid(2), angles, px_size);
    lineVal = 1;
    
    if minorAxes(pp) < 0.45  % noise threshold
        color = rgb('HotPink');
        color_text = rgb('DarkMagenta');
    elseif minorAxes(pp) > 0.9 % clump threshold
        color = rgb('DeepPink');
        color_text = rgb('DarkMagenta');
    else
        color = rgb('DarkTurquoise'); % single cells
        color_text = rgb('MidnightBlue');
    end

    figure(3)
    hold on
    plot(x_rotated,y_rotated,'Color',color,'lineWidth',lineVal)
    text((centroid(1)+1)/px_size, (centroid(2)+1)/px_size, num2str(pp),'Color',color_text,'FontSize',10)
    
end

%% G6 (H07) - S1 2021-06-15 (stationary phase)


clc
clear

% 0. initialize experiment meta data
px_size = 11/100; % 11 um pixels with 100x magnification
experiment = '2021-06-15';

%data_folder = strcat('C:/Users/Kisa Naqvi/Documents/TropiniLab/Data/',experiment);
data_folder = strcat('/Users/jen/Documents/TropiniLab/Data/Kisa/',experiment);

cd(data_folder)
sample = 'S1'; % fixed with NM65 encoding probe

prefix = 'img_';
suffix = '_position000_time000000000_z000.tif';
name_phase = strcat(prefix,'channel000',suffix);
clear prefix suffix
    

% 0. manually choose one sample image (contains single cells and clumps)
cd(data_folder)
sDirectory = dir(strcat(sample,'_*'));
names = {sDirectory.name};
img_h01_exp = names{10}; % chosen image


% 1. load one 100x phase image
cd(data_folder)
cd(strcat(img_h01_exp,'/Default'))
img_phase = imread(name_phase); % read phase image
figure(1) % display phase image
imshow(img_phase, 'DisplayRange',[1000 3000]); %lowering right # increases num sat'd pxls
title('G6 S1 (100x), phase')


% 2. create mask using parameters in current data matrix
%    see segmentIntensity_GFP.m for full comments / annotation
phase_smoothed = imgaussfilt(img_phase,0.8);
bw = edge(phase_smoothed,'sobel');
se = strel('disk',1); % structuring element = disk of radius x pixels
bw_dil = imdilate(bw, se); % dilate
bw_fill = imfill(bw_dil, 'holes'); % fill
bw_final = imerode(imerode(bw_fill,se),se); % erode twice
figure(2)
imshow(bw_final)
title('G6 S1 (100x), mask')


% 3. measure particle features
cc = bwconncomp(bw_final);
stats = regionprops(cc,'Area','Centroid','MajorAxisLength','MinorAxisLength','Eccentricity','Orientation');
clear bw bw_dil bw_fill se phase_smoothed
        

% 4. overlay outlines of all tracked particles
majorAxes = extractfield(stats,'MajorAxisLength')'.*px_size;
minorAxes = extractfield(stats,'MinorAxisLength')'.*px_size;
angles = extractfield(stats,'Orientation')';

% 5. vary min and max measurements to determine...
%      i) below what threshold is background noise?
%     ii) within what thresholds is a single-cell?
%    iii) above what threshold is a clump of cells?

figure(3) % this one will get overlaid with particle outlines
imshow(img_phase, 'DisplayRange',[1000 3000]);
title('G6 S1 (100x), phase')

for pp = 1:length(stats)
    
    centroid = stats(pp).Centroid.*px_size;
    [x_rotated, y_rotated] = drawEllipse2(pp, majorAxes, minorAxes, centroid(1), centroid(2), angles, px_size);
    lineVal = 1;
    
    if minorAxes(pp) < 0.6  % noise threshold
        color = rgb('HotPink');
        color_text = rgb('DarkMagenta');
    elseif minorAxes(pp) > 0.9 % clump threshold
        color = rgb('DeepPink');
        color_text = rgb('DarkMagenta');
    else
        color = rgb('DarkTurquoise'); % single cells
        color_text = rgb('MidnightBlue');
    end

    figure(3)
    hold on
    plot(x_rotated,y_rotated,'Color',color,'lineWidth',lineVal)
    text((centroid(1)+1)/px_size, (centroid(2)+1)/px_size, num2str(pp),'Color',color_text,'FontSize',10)
    
end

%% G6 (H07) - S2 2021-06-15 (stationary phase)


clc
clear

% 0. initialize experiment meta data
px_size = 11/100; % 11 um pixels with 100x magnification
experiment = '2021-06-15';

%data_folder = strcat('C:/Users/Kisa Naqvi/Documents/TropiniLab/Data/',experiment);
data_folder = strcat('/Users/jen/Documents/TropiniLab/Data/Kisa/',experiment);

cd(data_folder)
sample = 'S2'; % fixed with NM65 encoding probe

prefix = 'img_';
suffix = '_position000_time000000000_z000.tif';
name_phase = strcat(prefix,'channel000',suffix);
name_gfp = strcat(prefix,'channel001',suffix);

clear prefix suffix
    

% 0. manually choose one sample image (contains single cells and clumps)
cd(data_folder)
sDirectory = dir(strcat(sample,'_*'));
names = {sDirectory.name};

img_h01_exp = names{8}; % chosen image



% 1. load one 100x phase image
cd(data_folder)
cd(strcat(img_h01_exp,'/Default'))
img_phase = imread(name_phase); % read phase image

img_gfp = imread(name_gfp); 
figure(1) % display phase image
imshow(img_gfp, 'DisplayRange',[500 1000]); %lowering right # increases num sat'd pxls
title('G6 S2 (100x), phase')


% 2. create mask using parameters in current data matrix
%    see segmentIntensity_GFP.m for full comments / annotation
phase_smoothed = imgaussfilt(img_phase,0.8);
bw = edge(phase_smoothed,'sobel');
se = strel('disk',1); % structuring element = disk of radius x pixels
bw_dil = imdilate(bw, se); % dilate
bw_fill = imfill(bw_dil, 'holes'); % fill
bw_final = imerode(imerode(bw_fill,se),se); % erode twice
figure(2)
imshow(bw_final)
title('G6 S2 (100x), mask')


% 3. measure particle features
cc = bwconncomp(bw_final);
stats = regionprops(cc,'Area','Centroid','MajorAxisLength','MinorAxisLength','Eccentricity','Orientation');
clear bw bw_dil bw_fill se phase_smoothed
        

% 4. overlay outlines of all tracked particles
majorAxes = extractfield(stats,'MajorAxisLength')'.*px_size;
minorAxes = extractfield(stats,'MinorAxisLength')'.*px_size;
angles = extractfield(stats,'Orientation')';

% 5. vary min and max measurements to determine...
%      i) below what threshold is background noise?
%     ii) within what thresholds is a single-cell?
%    iii) above what threshold is a clump of cells?

figure(3) % this one will get overlaid with particle outlines

imshow(img_gfp, 'DisplayRange',[500 1000]);
title('G6 S2 (100x), phase')

for pp = 1:length(stats)
    
    centroid = stats(pp).Centroid.*px_size;
    [x_rotated, y_rotated] = drawEllipse2(pp, majorAxes, minorAxes, centroid(1), centroid(2), angles, px_size);
    lineVal = 1;
    
    if minorAxes(pp) < 0.6  % noise threshold
        color = rgb('HotPink');
        color_text = rgb('DarkMagenta');
    elseif minorAxes(pp) > 0.9 % clump threshold
        color = rgb('DeepPink');
        color_text = rgb('DarkMagenta');
    else
        color = rgb('DarkTurquoise'); % single cells
        color_text = rgb('MidnightBlue');
    end

    figure(3)
    hold on
    plot(x_rotated,y_rotated,'Color',color,'lineWidth',lineVal)
    text((centroid(1)+1)/px_size, (centroid(2)+1)/px_size, num2str(pp),'Color',color_text,'FontSize',10)
    
end

%% H03 + H06 - S3 2021-06-15 (stationary phase)


clc
clear

% 0. initialize experiment meta data
px_size = 11/100; % 11 um pixels with 100x magnification
experiment = '2021-06-15';

%data_folder = strcat('C:/Users/Kisa Naqvi/Documents/TropiniLab/Data/',experiment);
data_folder = strcat('/Users/jen/Documents/TropiniLab/Data/Kisa/',experiment);
cd(data_folder)
sample = 'S3'; 


prefix = 'img_';
suffix = '_position000_time000000000_z000.tif';
name_phase = strcat(prefix,'channel000',suffix);
clear prefix suffix
    

% 0. manually choose one sample image (contains single cells and clumps)
cd(data_folder)
sDirectory = dir(strcat(sample,'_*'));
names = {sDirectory.name};
img_h01_exp = names{4}; % chosen image


% 1. load one 100x phase image
cd(data_folder)
cd(strcat(img_h01_exp,'/Default'))
img_phase = imread(name_phase); % read phase image
figure(1) % display phase image
imshow(img_phase, 'DisplayRange',[1000 3000]); %lowering right # increases num sat'd pxls
title('H03 + H06 (100x), phase')


% 2. create mask using parameters in current data matrix
%    see segmentIntensity_GFP.m for full comments / annotation
phase_smoothed = imgaussfilt(img_phase,0.8);
bw = edge(phase_smoothed,'sobel');
se = strel('disk',1); % structuring element = disk of radius x pixels
bw_dil = imdilate(bw, se); % dilate
bw_fill = imfill(bw_dil, 'holes'); % fill
bw_final = imerode(imerode(bw_fill,se),se); % erode twice
figure(2)
imshow(bw_final)
title('H03 + H06  (100x), mask')


% 3. measure particle features
cc = bwconncomp(bw_final);
stats = regionprops(cc,'Area','Centroid','MajorAxisLength','MinorAxisLength','Eccentricity','Orientation');
clear bw bw_dil bw_fill se phase_smoothed
        

% 4. overlay outlines of all tracked particles
majorAxes = extractfield(stats,'MajorAxisLength')'.*px_size;
minorAxes = extractfield(stats,'MinorAxisLength')'.*px_size;
angles = extractfield(stats,'Orientation')';

% 5. vary min and max measurements to determine...
%      i) below what threshold is background noise?
%     ii) within what thresholds is a single-cell?
%    iii) above what threshold is a clump of cells?

figure(3) % this one will get overlaid with particle outlines
imshow(img_phase, 'DisplayRange',[1000 3000]);
title('H03 + H06  (100x), phase')

for pp = 1:length(stats)
    
    centroid = stats(pp).Centroid.*px_size;
    [x_rotated, y_rotated] = drawEllipse2(pp, majorAxes, minorAxes, centroid(1), centroid(2), angles, px_size);
    lineVal = 1;
    
    if minorAxes(pp) < 0.7  % noise threshold
        color = rgb('HotPink');
        color_text = rgb('DarkMagenta');
    elseif minorAxes(pp) > 1.6 % clump threshold
        color = rgb('DeepPink');
        color_text = rgb('DarkMagenta');
    else
        color = rgb('DarkTurquoise'); % single cells
        color_text = rgb('MidnightBlue');
    end

    figure(3)
    hold on
    plot(x_rotated,y_rotated,'Color',color,'lineWidth',lineVal)
    text((centroid(1)+1)/px_size, (centroid(2)+1)/px_size, num2str(pp),'Color',color_text,'FontSize',10)
    
end

