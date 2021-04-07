%% segment phase and fluorescent image and huild data matrix

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

% last updated: jen, 2021 April 6
% commit: first commit, use whos_a_cell to quantify intensity in single cells and clumps in 2021-03-31 experiment


%% Part ONE: measurements from raw images 

clc
clear


% 0. initialize experiment data
px_size = 11/150; % 11 um pixels with 150x magnification

experiment = '2021-03-31';
data_folder = strcat('/Users/jen/Documents/TropiniLab/Molecular_tools/HiPR_fish/',experiment);
cd(data_folder)

%sample_directory = dir('s*');
%samples = {samples.name};
samples = {'s1', 's2', 's3', 's4', 's5', 's6', 's7', 's8', 's9', 's10', 's11', 's12'};



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
        
        %   1a. read phase, gfp, mcherry and dapi images
        img_phase = imread(name_phase);
        img_gfp = imread(name_gfp);
        
        
        %   1b. make mask from phase image
        figure(1)
        imshow(img_phase, 'DisplayRange',[500 2500]); %lowering right # increases num sat'd pxls
        
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
        
        
        
        %   2. quantify fluorescence intensity inside mask (all channels)
        
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


%% Part TWO: trim measured data and create data structure

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
    
    % 2a. convert x,y coordinate data
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
    if sample < 5
        LowerBound = 1;     % lower bound for 381 (see whos_a_cell.m)
    elseif sample < 9
        LowerBound = 0.7;   % lower bound for 505 (see whos_a_cell.m)
    else
        LowerBound = 0.6;   % lower bound for 488 (see whos_a_cell.m)
    end
    UpperBound = 176;     % whole image length
    p_trim = ParticleTrim_glycogen(parameter_unit,TrimField,LowerBound,UpperBound);
    

    % 4. store final data 
    converted_data{1,sample} = p_trim;
    
end
clear sample sample_particles p_trim UpperBound LowerBound
     
    
%% Part THREE: visualize measured data

% 1. isolate single cells from clumps
% 2. plot absolute intensities of background, single cells, clumps
% 3. plot single cell and clump intensities normalized by background


counter_381 = 0; counter_505 = 0; counter_488 = 0;
ct_381 = 0; ct_505 = 0; ct_488 = 0;

% for each sample
for smpl = 1:length(samples)
    
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
    
    
    % 2. plot absolute intensities of background, single cells, clumps
    
    abs_single = [mean(single_bg); mean(single_gfp)];
    abs_clump = [mean(clump_bg); mean(clump_gfp)];
    
    std_single = [std(single_bg); std(single_gfp)];
    std_clump = [std(clump_bg); std(clump_gfp)];
    
    n_single = length(single_bg);
    n_clump = length(clump_bg);
    
    sem_single = std_single./sqrt(n_single);
    sem_clump = std_clump./sqrt(n_clump);
    
    % group subplots by strain
    if smpl < 5
        counter_381 = counter_381 + 1;
        
        figure(7)
        subplot(1,4,counter_381)
        bar([abs_single; abs_clump])
        hold on
        errorbar([abs_single; abs_clump], [std_single; std_clump],'.')
        set(gca,'xticklabel',{'BG','1x','BG', 'Clump'})
        title(strcat(samples{smpl},', n =',num2str(n_single),' and n =',num2str(n_clump)))
        ylim([0 650])
        
    elseif smpl < 9
        counter_505 = counter_505 + 1;
        
        figure(8)
        subplot(1,4,counter_505)
        bar([abs_single; abs_clump])
        hold on
        errorbar([abs_single; abs_clump], [std_single; std_clump],'.')
        set(gca,'xticklabel',{'BG','1x','BG', 'Clump'})
        title(strcat(samples{smpl},', n =',num2str(n_single),' and n =',num2str(n_clump)))
        ylim([0 650])
        
    else
        counter_488 = counter_488 + 1;
        
        figure(9)
        subplot(1,4,counter_488)
        bar([abs_single; abs_clump])
        hold on
        errorbar([abs_single; abs_clump], [std_single; std_clump],'.')
        set(gca,'xticklabel',{'BG','1x','BG', 'Clump'})
        title(strcat(samples{smpl},', n =',num2str(n_single),' and n =',num2str(n_clump)))
        ylim([0 650])
        
    end
    
    
    % 3. plot single cell and clump intensities normalized by background
    
    % cell fluorescence normalized by bg fluorescence
    norm_single = single_gfp./single_bg;
    norm_clump = clump_gfp./clump_bg;
    
    % mean of normalized values
    norm_means = [mean(norm_single); mean(norm_clump)];
    
    % error of normalized values
    norm_std = [std(norm_single); std(norm_clump)];
    norm_n = [length(norm_single); length(norm_clump)];
    norm_sem = norm_std./sqrt(norm_n);
    
    if smpl < 5
        ct_381 = ct_381 + 1;
        
        figure(17)
        subplot(1,4,ct_381)
        bar(norm_means)
        hold on
        errorbar(norm_means, norm_std,'.')
        set(gca,'xticklabel',{'1x','Clump'})
        title(strcat(samples{smpl},', n =',num2str(norm_n(1)),' and n =',num2str(norm_n(2))))
        ylim([0 2.5])
        
    elseif smpl < 9
        ct_505 = ct_505 + 1;
        
        figure(18)
        subplot(1,4,ct_505)
        bar(norm_means)
        hold on
        errorbar(norm_means, norm_std,'.')
        set(gca,'xticklabel',{'1x','Clump'})
        title(strcat(samples{smpl},', n =',num2str(norm_n(1)),' and n =',num2str(norm_n(2))))
        ylim([0 2.5])
        
    else
        ct_488 = ct_488 + 1;
        
        figure(19)
        subplot(1,4,ct_488)
        bar(norm_means)
        hold on
        errorbar(norm_means, norm_std,'.')
        set(gca,'xticklabel',{'1x','Clump'})
        title(strcat(samples{smpl},', n =',num2str(norm_n(1)),' and n =',num2str(norm_n(2))))
        ylim([0 2.5])
        
    end
    
    
end





