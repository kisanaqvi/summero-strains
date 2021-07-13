%% segment phase and fluorescent image and build data matrix

% goal: extract cell size and HiPR-FISH signal from phase and GFP images

% strategy:
%
% Part ONE: measurements from raw images
%
%   0. initialize experiment data
%   0. define image name of each channel
%   0. for each sample, build directory and loop through stacks
%   1. for each stack per sample, make a mask to isolate pixels representing cells
%   2. quantify fluorescence intensity inside mask 
%   3. quantify fluorescence intenstiy outside mask
%
% Part TWO:
%
%   4. trim measured particules by size
%   5. save new data matrix
%
% Part THREE:
%
%   6. plot channel by channel comparison


% ok, let's go!

% last updated: jen, 2021 June 16
% commit: 2021-06-04 analysis, testing collab merge with kisa


%% Part ONE: measurements from raw images 

clc
clear
cd('/Users/jen/Documents/TropiniLab/Data/such-hipr/sourcedata')
load('metadata.mat')

% 0. initialize experiment data
index = 4; % 2021-06-04
date = metadata{index}.date;
magnification = metadata{index}.magnification;
samples = metadata{index}.samples;

data_folder = strcat('/Users/jen/Documents/TropiniLab/Data/Kisa/',date);
cd(data_folder)
px_size = 11/magnification; % 11 um pixels with 150x magnification



% 0. define image name of each channel
prefix = 'img_';
suffix = '_position000_time000000000_z000.tif';
name_phase = strcat(prefix,'channel000',suffix);
name_gfp = strcat(prefix,'channel001',suffix);
clear prefix suffix



% 0. for each sample, build directory and loop through stacks
for ss = 1:length(samples)
    
    cd(data_folder)
    sDirectory = dir(strcat(samples{ss},'_*'));
    names = {sDirectory.name};
    
    
    %   1. for each stack per sample, make a mask to isolate pixels representing cells
    for stk = 1:length(names)
        
        cd(data_folder)
        current_stack = names{stk};
        cd(strcat(current_stack,'/Default'))
        
        %   1a. read phase, gfp images
        img_phase = imread(name_phase);
        img_gfp = imread(name_gfp);
        
        
        %   1b. make mask from phase image
        figure(1)
        imshow(img_phase, 'DisplayRange',[1000 3000]); %lowering right # increases num sat'd pxls
        
        %   i. gaussian smoothing of phase image
        phase_smoothed = imgaussfilt(img_phase,0.8);
        %figure(2)
        %imshow(phase_smoothed, 'DisplayRange',[500 2500]);
        
        %   ii. edge detection of cells in phase image
        bw = edge(phase_smoothed,'sobel');
        %figure(3)
        %imshow(bw);%, 'DisplayRange',[2000 6000]);
        
        %   iii. clean up edge detection
        se = strel('disk',1); % structuring element; disk-shaped with radius of 3 pixels
        bw_dil = imdilate(bw, se); % dilate
        %figure(4)
        %imshow(bw_dil);
        
        bw_fill = imfill(bw_dil, 'holes'); % fill
        %figure(5)
        %imshow(bw_fill);
        
        bw_final = imerode(imerode(bw_fill,se),se); % erode 2x to smooth; this is your final mask
        figure(6)
        imshow(bw_final)
        
        %   iv. segment cells
        cc = bwconncomp(bw_final);
        
        %   v. gather properties for each identified particle
        stats = regionprops(cc,'Area','Centroid','MajorAxisLength','MinorAxisLength','Eccentricity','Orientation');
        clear bw bw_dil bw_fill se phase_smoothed
        
        
        
        %   2. quantify fluorescence intensity inside mask (in GFP channel)
        
        %   2a. overlay mask with fluorescence image by dot-product
        masked_gfp = bw_final .* double(img_gfp); % convert uint16 image to double class

        
        %   2b. compute the mean fluorescence intensity for each particle in mask
        for pp = 1:cc.NumObjects
            
            pixel_id_of_cell = []; % initialize
            pixel_id_of_cell = cc.PixelIdxList{pp}; % pixel index where the cell of interest is
            
            cell_intensity_gfp(pp) = mean(masked_gfp(pixel_id_of_cell)); % compute the mean intensity of pixels in the cell
            
        end
        clear pp pixel_id_of_cell
        

        
        %   3. quantify fluorescence intenstiy outside mask (all channels)
        
        %   3a. background mask is inverse of cell mask
        bg_mask = imcomplement(bw_final);
        %figure(10)
        %imshow(bg_mask)
        
        %   3b. overlay background mask with fluorescence image by dot-product
        bg_gfp = bg_mask .* double(img_gfp); % convert uint16 image to double class

        
        %   3c. compute the mean background intensity for each channel
        bg_gfp_mean = mean(mean(bg_gfp));
        clear bg_mask
        
        
       
        % 4. store mean fluorescent intensity of each cell and background in data structure
        for particle = 1:cc.NumObjects
            
            % cell intensity
            stats(particle).gfp_cell = cell_intensity_gfp(particle);
            
            % mean background intensity
            stats(particle).gfp_bg = bg_gfp_mean;
           
        end
        clear cell_intensity_gfp img_phase img_gfp  bg_gfp 
        clear masked_gfp  bw_final particle cc
        
        
        % 5. store particle measurements into cell per stack
        dm{stk,ss} = stats;
        clear bg_gfp_mean
        
    end
    clear stats
    
end
clear name_gfp name_phase names 
clear ss stk current_stack

cd('/Users/jen/Documents/TropiniLab/Data/Kisa')
save(strcat('dm-segmentIntensity-',date,'.mat'),'dm')

%% Part TWO: trim measured data and create data structure


clear
clc
cd('/Users/jen/Documents/TropiniLab/Data/such-hipr/sourcedata')
%cd('/Users/jen/Documents/TropiniLab/Data/Kisa') % move metadata to this path
load('metadata.mat')

% 0. initialize experiment data
index = 4; % 2021-06-04
date = metadata{index}.date;
cd('/Users/jen/Documents/TropiniLab/Data/Kisa')
load(strcat('dm-segmentIntensity-',date,'.mat'))

samples = metadata{index}.samples;
magnification = metadata{1}.magnification;
px_size = 11/magnification; % 11 um pixels with 150x magnification


% 1. concatenate data from same sample
for col = 1:length(samples)
    
    current_sample = dm(:,col);
    num_stacks = sum(~cellfun(@isempty,current_sample));
    
    combined_sample_data = [];
    for sstack = 1:num_stacks
        stk_data = current_sample{sstack,1};
        combined_sample_data = [combined_sample_data; stk_data];
    end
    clear sstack
    combined_particles{1,col} = combined_sample_data;
    
end
clear col combined_sample_data num_stacks current_sample stk_data



% 2. convert measurements to microns based on imaging specifications     
for sample = 1:length(samples)
    
    sample_particles = combined_particles{1,sample};
    if isempty(sample_particles) == 1
        continue
    end
    
    
    % 2a. convert x,y coordinate of particle centroid
    x_position = [];
    y_position = [];
    for ii = 1:length(sample_particles)
        
        centroid = sample_particles(ii).Centroid.*px_size;
        x_position = [x_position; centroid(1)];
        y_position = [y_position; centroid(2)];
        
    end
    parameter_unit.X = x_position;
    parameter_unit.Y = y_position;
    clear particle x_position y_position centroid ii
    
    
    % 2b. convert area
    area = extractfield(sample_particles,'Area')';
    parameter_unit.A = area.*px_size^2;
    clear area
    
    
    % 2c. major & minor axis
    majorAxis = extractfield(sample_particles,'MajorAxisLength')';
    minorAxis = extractfield(sample_particles,'MinorAxisLength')';
    parameter_unit.MajAx = majorAxis.*px_size;
    parameter_unit.MinAx = minorAxis.*px_size;
    clear majorAxis minorAxis
    
    
    % 2d. cell intensity
    GFP_cell = extractfield(sample_particles,'gfp_cell')';
    parameter_unit.gfp_cell = GFP_cell;
    clear GFP_cell
    
    
    % 2e. background intensity of corresponding image
    GFP_bg = extractfield(sample_particles,'gfp_bg')';
    parameter_unit.gfp_bg = GFP_bg;
    clear GFP_bg
    
    
    % 2f. eccentricity and angle
    ecc = extractfield(sample_particles,'Eccentricity')';
    angle = extractfield(sample_particles,'Orientation')';
    parameter_unit.Ecc = ecc;
    parameter_unit.Angle = angle;
    clear ecc angle
    
    
    % 3. trim particles by width
    %    values are set as recorded in whos_a_cell.m
    
    % 3b. trim by width
    TrimField = 'MinAx';  % choose relevant characteristic to restrict, run several times to apply for several fields
    if sample < 4
        LowerBound = 0.6;     % lower bound for exponential (see whos_a_cell.m)
        UpperBound = 0.9;
    elseif sample >= 4
        LowerBound = 0.6;       % lower bound for stationary (see whos_a_cell.m)
        UpperBound = 0.9;
    end
    p_trim = ParticleTrim_glycogen(parameter_unit,TrimField,LowerBound,UpperBound);
    

    % 4. store final data 
    converted_data{1,sample} = p_trim;
    
end
clear sample sample_particles p_trim UpperBound LowerBound
    
%% Part THREE: visualize measured data

% 1. isolate single cells from clumps
% 2. plot absolute intensities of background, single cells, clumps
% 3. plot single cell and clump intensities normalized by background

% a counter for each
counter_G1 = 0; counter_G2 = 0; counter_G3 = 0;
ct_G1 = 0; ct_G2 = 0; ct_G3 = 0;

% for each sample
for smpl = 1:3%length(samples)
    
    % 0. gather cell width and intensity data
    sample_data = converted_data{1,smpl};
    cell_width = sample_data.MinAx;
    cell_gfp = sample_data.gfp_cell; % mean cell fluorescence 
    bg_gfp = sample_data.gfp_bg; % mean bg fluorescence
    

    % 1. isolate single cells from clumps
    if smpl < 5
        clumpThresh = 1.6; % min width of clumps for 381 (see whos_a_cell.m)
    elseif smpl < 9
        clumpThresh = 1;   % min width of clumps for 505 (see whos_a_cell.m)
    else
        clumpThresh = 1;   % min width of clumps for 488 (see whos_a_cell.m)
    end
    
    single_gfp = cell_gfp(cell_width <= clumpThresh);
    single_bg = bg_gfp(cell_width <= clumpThresh);
    clump_gfp = cell_gfp(cell_width > clumpThresh);
    clump_bg = bg_gfp(cell_width > clumpThresh);
    
    
    % 2. box plots of absolute intensities of background, single cells, clumps
      
     n_single = length(single_bg);
     n_clump = length(clump_bg);
    
    % group subplots by strain
    if smpl < 5
        counter_G1 = counter_G1 + 1;
        
        figure(7)
        %subplot(1,length(samples),counter_G1)
        subplot(1,3,counter_G1)
        x = [single_bg; single_gfp; clump_bg; clump_gfp];
        g = [zeros(length(single_bg), 1); ones(length(single_gfp), 1); 2*ones(length(clump_bg), 1); 3*ones(length(clump_gfp), 1)];
        boxplot(x,g)
        set(gca,'xticklabel',{'BG','1x'})
        title(strcat(samples{smpl},', n =',num2str(n_single)))
        ylim([200 800])
        
    elseif smpl < 9
        counter_G2 = counter_G2 + 1;
        
        figure(8)
        subplot(1,length(samples),counter_G2)
        x = [single_bg; single_gfp; clump_bg; clump_gfp];
        g = [zeros(length(single_bg), 1); ones(length(single_gfp), 1); 2*ones(length(clump_bg), 1); 3*ones(length(clump_gfp), 1)];
        boxplot(x,g)
        set(gca,'xticklabel',{'BG','1x'})
        title(strcat(samples{smpl},', n =',num2str(n_single)))
        ylim([200 1600])
        
    else
        counter_G3 = counter_G3 + 1;
        
        figure(9)
        subplot(1,4,counter_G3)
        x = [single_bg; single_gfp; clump_bg; clump_gfp];
        g = [zeros(length(single_bg), 1); ones(length(single_gfp), 1); 2*ones(length(clump_bg), 1); 3*ones(length(clump_gfp), 1)];
        boxplot(x,g)
        set(gca,'xticklabel',{'BG','1x'})
        title(strcat(samples{smpl},', n =',num2str(n_single)))
        ylim([200 1600])
        
    end
    
    
    % 3. plot single cell and clump intensities normalized by background
    
    % cell fluorescence normalized by bg fluorescence
    norm_single = single_gfp./single_bg;
    norm_clump = clump_gfp./clump_bg;
    norm_n = [length(norm_single); length(norm_clump)];

    
    if smpl < 5
        ct_G1 = ct_G1 + 1;
        
        figure(17)
        %subplot(1,4,ct_G1)
        subplot(1,3,counter_G1)
        xx = [norm_single; norm_clump];
        gg = [zeros(length(norm_single), 1); ones(length(norm_clump), 1)];
        boxplot(xx,gg)
        set(gca,'xticklabel',{'1x'})
        title(strcat(samples{smpl},', n =',num2str(norm_n(1))))
        ylim([0.8 3])
        
    elseif smpl < 9
        ct_G2 = ct_G2 + 1;
        
        figure(18)
        subplot(1,4,ct_G2)
        xx = [norm_single; norm_clump];
        gg = [zeros(length(norm_single), 1); ones(length(norm_clump), 1)];
        boxplot(xx,gg)
        set(gca,'xticklabel',{'1x'})
        title(strcat(samples{smpl},', n =',num2str(norm_n(1))))
        ylim([0.8 5])
        
    else
        ct_G3 = ct_G3 + 1;
        
        figure(19)
        subplot(1,4,ct_G3)
        xx = [norm_single; norm_clump];
        gg = [zeros(length(norm_single), 1); ones(length(norm_clump), 1)];
        boxplot(xx,gg)
        set(gca,'xticklabel',{'1x'})
        title(strcat(samples{smpl},', n =',num2str(norm_n(1))))
        ylim([0.8 5])
        
    end
    
    
end

%% Part FOUR. save boxplots

cd('/Users/jen/Documents/TropiniLab/Data/HiPR_fish')

figure(7)
saveas(gcf,strcat(date,'-absolute-381'),'epsc')
close(gcf)

figure(8)
saveas(gcf,strcat(date,'-absolute-505'),'epsc')
close(gcf)

figure(9)
saveas(gcf,strcat(date,'-absolute-488'),'epsc')
close(gcf)

figure(17)
saveas(gcf,strcat(date,'-norm-381'),'epsc')
close(gcf)

figure(18)
saveas(gcf,strcat(date,'-norm-505'),'epsc')
close(gcf)

figure(19)
saveas(gcf,strcat(date,'-norm-488'),'epsc')
close(gcf)