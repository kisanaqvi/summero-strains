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

% last updated: Kisa, 2021 June 23
% commit: 2021-06-23 analysis with stronger signal but weak for novel probe


%% Part ONE: measurements from raw images 

clc
clear
%cd('/Users/jen/summero-strains')
cd('C:/Users/Kisa Naqvi/Documents/TropiniLab/summero-strains-master')
load('metadata.mat')

% 0. initialize experiment data
index = 6; % 2021-06-22
date = metadata{index}.date;
magnification = metadata{index}.magnification;
samples = metadata{index}.samples;

%data_folder = strcat('/Users/jen/Documents/TropiniLab/Data/Kisa/',date);
data_folder = strcat('C:/Users/Kisa Naqvi/Documents/TropiniLab/Data/',date);
cd(data_folder)
px_size = 11/magnification; % 11 um pixels before magnification



% 0. define image name of each channel
prefix = 'img_';
suffix = '_position000_time000000000_z000.tif';
name_phase = strcat(prefix,'channel000',suffix);
name_gfp = strcat(prefix,'channel001',suffix);
name_mcherry = strcat(prefix,'channel002',suffix);
name_dapi = strcat(prefix,'channel003',suffix);
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
        
        %   1a. read phase, gfp, mCherry, DAPI images
        img_phase = imread(name_phase);
        img_gfp = imread(name_gfp);
        img_mcherry = imread(name_mcherry);
        img_dapi = imread(name_dapi);
        
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
        mask_name = strcat('mask_',current_stack);
        cd(data_folder)
        saveas(gcf,mask_name)
        close(gcf)
        
        %   iv. segment cells
        cc = bwconncomp(bw_final);
        
        %   v. gather properties for each identified particle
        stats = regionprops(cc,'Area','Centroid','MajorAxisLength','MinorAxisLength','Eccentricity','Orientation');
        clear bw bw_dil bw_fill se phase_smoothed
        
        
        
        %   2. quantify fluorescence intensity inside mask (in GFP channel)
        
        %   2a. overlay mask with fluorescence image by dot-product
        masked_gfp = bw_final .* double(img_gfp); % convert uint16 image to double class
        masked_mcherry = bw_final .* double(img_mcherry);
        masked_dapi = bw_final .* double(img_dapi); 
        
        %   2b. compute the mean fluorescence intensity for each particle in mask
        for pp = 1:cc.NumObjects
            
            pixel_id_of_cell = []; % initialize
            pixel_id_of_cell = cc.PixelIdxList{pp}; % pixel index where the cell of interest is
            
            cell_intensity_gfp(pp) = mean(masked_gfp(pixel_id_of_cell)); % compute the mean intensity of pixels in the cell
            cell_intensity_mcherry(pp) = mean(masked_mcherry(pixel_id_of_cell));
            cell_intensity_dapi(pp) = mean(masked_dapi(pixel_id_of_cell));
            
        end
        clear pp pixel_id_of_cell
        

        
        %   3. quantify fluorescence intenstiy outside mask (all channels)
        
        %   3a. background mask is inverse of cell mask
        bg_mask = imcomplement(bw_final);
        %figure(10)
        %imshow(bg_mask)
        
        %   3b. overlay background mask with fluorescence image by dot-product
        bg_gfp = bg_mask .* double(img_gfp); % convert uint16 image to double class
        bg_mcherry = bg_mask .* double(img_mcherry);
        bg_dapi = bg_mask .* double(img_dapi);
        
        %   3c. compute the mean background intensity for each channel
        bg_gfp_mean = mean(mean(bg_gfp));
        bg_mcherry_mean = mean(mean(bg_mcherry));
        bg_dapi_mean = mean(mean(bg_dapi));
        clear bg_mask
        
        
       
        % 4. store mean fluorescent intensity of each cell and background in data structure
        for particle = 1:cc.NumObjects
            
            % cell intensity
            stats(particle).gfp_cell = cell_intensity_gfp(particle);
            stats(particle).mcherry_cell = cell_intensity_mcherry(particle);
            stats(particle).dapi_cell = cell_intensity_dapi(particle);
            
            % mean background intensity
            stats(particle).gfp_bg = bg_gfp_mean;
            stats(particle).mcherry_bg = bg_mcherry_mean;
            stats(particle).dapi_bg = bg_dapi_mean;
        end
        clear cell_intensity_gfp cell_intensity_mcherry cell_intensity_dapi
        clear masked_gfp masked_mcherry masked_dapi bw_final particle cc
        clear img_phase img_gfp img_mcherry img_dapi bg_gfp bg_mcherry bg_dapi
        
        
        % 5. store particle measurements into cell per stack
        dm{stk,ss} = stats;
        clear bg_gfp_mean bg_mcherry_mean bg_dapi_mean
        
    end
    clear stats
    
end
clear name_gfp name_phase names name_dapi name_mcherry
clear ss stk current_stack

%cd('/Users/jen/Documents/TropiniLab/Data/Kisa')
cd('C:/Users/Kisa Naqvi/Documents/TropiniLab/Data')
save(strcat('dm-segmentIntensity-',date,'.mat'),'dm')

%% Part TWO: trim measured data and create data structure


clear
clc
%cd('/Users/jen/summero-strains')
cd('C:/Users/Kisa Naqvi/Documents/TropiniLab/summero-strains-master')
load('metadata.mat')
%cd('/Users/jen/Documents/TropiniLab/Data/Kisa') % move metadata to this path



% 0. initialize experiment data
index = 6; % 2021-06-22
date = metadata{index}.date;
%cd('/Users/jen/Documents/TropiniLab/Data/Kisa')
cd('C:/Users/Kisa Naqvi/Documents/TropiniLab/Data')
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
    mCherry_cell = extractfield(sample_particles,'mcherry_cell')';
    DAPI_cell = extractfield(sample_particles,'dapi_cell')';
    parameter_unit.gfp_cell = GFP_cell;
    parameter_unit.mcherry_cell = mCherry_cell;
    parameter_unit.dapi_cell = DAPI_cell;
    clear GFP_cell mCherry_cell DAPI_cell
    clear GFP_cell
    
    
    % 2e. background intensity of corresponding image
    GFP_bg = extractfield(sample_particles,'gfp_bg')';
    mCherry_bg = extractfield(sample_particles,'mcherry_bg')';
    DAPI_bg = extractfield(sample_particles,'dapi_bg')';
    parameter_unit.gfp_bg = GFP_bg;
    parameter_unit.mcherry_bg = mCherry_bg;
    parameter_unit.dapi_bg = DAPI_bg;
    clear GFP_bg mCherry_bg DAPI_bg
    
    
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
    if sample < 3
        LowerBound = 0.5;       % bounds for stationary G6 (see whos_a_cell.m)
        UpperBound = 1.2;
    elseif sample == 3
        LowerBound = 0.7;       % bounds for stationary H03 + H06 (see whos_a_cell.m)
        UpperBound = 1.6;
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


% for each sample
for smpl = 1:length(samples)
    
    % 0. gather cell width and intensity data
    sample_data = converted_data{1,smpl};
    cell_width = sample_data.MinAx;
    
    % mean cell fluorescence 
    cell_gfp = sample_data.gfp_cell;
    cell_mcherry = sample_data.mcherry_cell;
    cell_dapi = sample_data.dapi_cell;
    
    
    % mean bg fluorescence
    bg_gfp = sample_data.gfp_bg;
    bg_mcherry = sample_data.mcherry_bg;
    bg_dapi = sample_data.dapi_bg;
    

    % 1. isolate single cells from clumps
    if smpl < 3
        clumpThresh = 1.2; % min width of clumps for S1 & S2 (see whos_a_cell.m)
    elseif smpl == 3 
        clumpThresh = 1.6;
    end
    
    single_gfp = cell_gfp(cell_width <= clumpThresh);
    single_mcherry = cell_mcherry(cell_width <= clumpThresh);
    single_dapi = cell_dapi(cell_width <= clumpThresh);
    single_bg_gfp = bg_gfp(cell_width <= clumpThresh);
    single_bg_mcherry = bg_mcherry(cell_width <= clumpThresh);
    single_bg_dapi = bg_dapi(cell_width <= clumpThresh);

    
    % 2. box plots of absolute intensities of background vs single cells
     n_single = length(single_bg_gfp);
     
     figure(10+smpl)
     subplot(1,length(samples),1)
     x = [single_bg_gfp; single_gfp]; %clump_bg; clump_gfp];
     g = [zeros(length(single_bg_gfp), 1); ones(length(single_gfp), 1)];%; 2*ones(length(clump_bg), 1); 3*ones(length(clump_gfp), 1)];
     boxplot(x,g)
     set(gca,'xticklabel',{'BG','1x'})%'BG', 'Clump'})
     title(strcat('GFP:',samples{smpl},', n =',num2str(n_single)))%' and n =',num2str(n_clump)))
     ylim([200 3000])
     
     subplot(1,length(samples),2)
     x = [single_bg_mcherry; single_mcherry]; %clump_bg; clump_gfp];
     g = [zeros(length(single_bg_mcherry), 1); ones(length(single_mcherry), 1)];% 2*ones(length(clump_bg), 1); 3*ones(length(clump_gfp), 1)];
     boxplot(x,g)
     set(gca,'xticklabel',{'BG','1x'})%,'BG', 'Clump'})
     title(strcat('mCherry:',samples{smpl},', n =',num2str(n_single)))%,' and n =',num2str(n_clump)))
     ylim([200 3000])
     
     subplot(1,length(samples),3)
     x = [single_bg_dapi; single_dapi]; %clump_bg; clump_gfp];
     g = [zeros(length(single_bg_dapi), 1); ones(length(single_dapi), 1)];% 2*ones(length(clump_bg), 1); 3*ones(length(clump_gfp), 1)];
     boxplot(x,g)
     set(gca,'xticklabel',{'BG','1x'})%'BG', 'Clump'})
     title(strcat('DAPI:',samples{smpl},', n =',num2str(n_single)))%,' and n =',num2str(n_clump)))
     ylim([200 3000])
     
     
     % 3. plot single cell and clump intensities normalized by background
    
    % cell fluorescence normalized by bg fluorescence
    norm_single_gfp = single_gfp./single_bg_gfp;
    norm_single_mcherry = single_mcherry./single_bg_mcherry;
    norm_single_dapi = single_dapi./single_bg_dapi;
      
    %figure(20+smpl)
    figure(20)
    
    subplot(1,length(samples),smpl)
    xx = [norm_single_gfp; norm_single_mcherry; norm_single_dapi];
    gg = [zeros(length(norm_single_gfp), 1); ones(length(norm_single_mcherry), 1); 2*ones(length(norm_single_dapi), 1)];
    boxplot(xx,gg)
    set(gca,'xticklabel',{'novel','Carolinas','dapi'})
    title(strcat(samples{smpl},', n =',num2str(length(norm_single_gfp))))
    ylim([0.8 5])
        
    
     
end

