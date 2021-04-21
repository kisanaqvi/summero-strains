%% visualize ID'ed particles to determine trimming parameters


% goal: measure cell size from a phase image (150x),
%       overlay measurement onto phase image,
%       visualize and repeat to determine desired values
%       (e.g. values that identify single cells)

% strategy:
%
%   For each strain of interest
%
%       0. intialize data (experiment meta data)
%       1. load one 150x phase image
%       2. create mask using parameters in current data matrix
%       3. measure particle features
%       4. overlay outlines of all tracked particles
%       5. vary min and max measurements to determine...
%            i) below what threshold is background noise?
%           ii) within what thresholds is a single-cell?
%          iii) above what threshold is a clump of cells?

% ok, let's go!

% last updated: jen, 2021 April 20
% commit: L. plantarum in exponential phase


%% 381 - Lactobacillus plantarum ATCC BAA-793 (Firmicute)

clc
clear

% 0. initialize experiment meta data
px_size = 11/150; % 11 um pixels with 150x magnification
experiment = '2021-03-31';
data_folder = strcat('/Users/jen/Documents/TropiniLab/Molecular_tools/HiPR_fish/',experiment);
cd(data_folder)
sample = 's1';

prefix = 'img_';
suffix = '_position000_time000000000_z000.tif';
name_phase = strcat(prefix,'channel000',suffix);
clear prefix suffix
    

% 0. manually choose one sample image (contains single cells and clumps)
cd(data_folder)
sDirectory = dir(strcat(sample,'_*'));
names = {sDirectory.name};
img_381 = names{2}; % chosen image


% 1. load one 150x phase image
cd(data_folder)
cd(strcat(img_381,'/Default'))
img_phase = imread(name_phase); % read phase image
figure(1) % display phase image
imshow(img_phase, 'DisplayRange',[500 2500]); %lowering right # increases num sat'd pxls
title('381 Lactobacillus plantarum (150x), phase')

figure(3)
imshow(img_phase, 'DisplayRange',[500 2500]);
title('381 Lactobacillus plantarum (150x), phase')

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
title('381 Lactobacillus plantarum (150x), mask')


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
img_bw25113 = names{10}; % chosen image


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
