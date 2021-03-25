%% calculate cell vs bg signal intensity

% goal: determine extent of HiPR-FISH signal from cell mask vs background

% strategy:
%
% Part ONE: measurements from raw images
%
%   0. initialize experiment data
%   0. define image name of each channel
%   0. for each sample, build directory and loop through stacks
%   1. for each stack per sample, make a mask to isolate pixels representing cells
%   2. quantify fluorescence intensity inside mask (all channels)
%   3. quantify fluorescence intenstiy outside mask (all channels)
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

% last updated: jen, 2021 March 24
% commit: first commit, normalized fluorescence in 2021-03-23 experiment

%% Part ONE: measurements from raw images 

clc
clear


% 0. initialize experiment data
px_size = 11/100; % 11 um pixels with 100x magnification

experiment = '2021-03-23';
samples = {'bt404', 'cb', 'cc', 'cf455', 'cf538', 'entero6', 'ful7'};

data_folder = strcat('/Users/jen/Documents/TropiniLab/Molecular_tools/HiPR_fish/',experiment);
cd(data_folder)


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
    sDirectory = dir(strcat(samples{ss},'_100x*'));
    names = {sDirectory.name};
    
    
    %   1. for each stack per sample, make a mask to isolate pixels representing cells
    for stk = 1:length(names)
        
        cd(data_folder)
        current_stack = names{stk};
        cd(strcat(current_stack,'/Default'))
        
        %   1a. read phase, gfp, mcherry and dapi images
        img_phase = imread(name_phase);
        img_gfp = imread(name_gfp);
        img_mcherry = imread(name_mcherry);
        img_dapi = imread(name_dapi);
        
        %   1b. make mask from phase image
        %figure(1)
        %imshow(img_phase, 'DisplayRange',[2000 6000]); %lowering right # increases num sat'd pxls
 
        %   i. gaussian smoothing of phase image
        phase_smoothed = imgaussfilt(img_phase,0.8);
        %figure(2)
        %imshow(phase_smoothed, 'DisplayRange',[2000 6000]);
        
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
        
        bw_final = imerode(imerode(imerode(bw_fill,se),se),se); % erode 3x to smooth; this is your final mask
        %figure(6)
        %imshow(bw_final)
 
        %   iv. segment cells
        cc = bwconncomp(bw_final);
        
        %   v. gather properties for each identified particle
        stats = regionprops(cc,'Area','Centroid','MajorAxisLength','MinorAxisLength','Eccentricity','Orientation');
        clear bw bw_dil bw_fill se phase_smoothed
        
        
        
        %   2. quantify fluorescence intensity inside mask (all channels)
        
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
clear name_dapi name_gfp name_mcherry name_phase names 
clear ss stk current_stack


%% Part TWO: trim measured data and save final data structure


% 3. concatenate data from same sample
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



% 4. convert measurements to microns based on imaging specifications     
for sample = 1:length(samples)
    
    sample_particles = combined_particles{1,sample};
    
    % 4a. convert x,y coordinate data
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
    
    
    % 4b. convert area
    area = extractfield(sample_particles,'Area')';
    parameter_unit.A = area.*px_size^2;
    clear area
    
    
    % 4c. major & minor axis
    majorAxis = extractfield(sample_particles,'MajorAxisLength')';
    minorAxis = extractfield(sample_particles,'MinorAxisLength')';
    parameter_unit.MajAx = majorAxis.*px_size;
    parameter_unit.MinAx = minorAxis.*px_size;
    clear majorAxis minorAxis
    
    
    % 4d. cell intensity
    GFP_cell = extractfield(sample_particles,'gfp_cell')';
    mCherry_cell = extractfield(sample_particles,'mcherry_cell')';
    DAPI_cell = extractfield(sample_particles,'dapi_cell')';
    parameter_unit.gfp_cell = GFP_cell;
    parameter_unit.mcherry_cell = mCherry_cell;
    parameter_unit.dapi_cell = DAPI_cell;
    clear GFP_cell mCherry_cell DAPI_cell
    
    
    % 4e. background intensity of corresponding image
    GFP_bg = extractfield(sample_particles,'gfp_bg')';
    mCherry_bg = extractfield(sample_particles,'mcherry_bg')';
    DAPI_bg = extractfield(sample_particles,'dapi_bg')';
    parameter_unit.gfp_bg = GFP_bg;
    parameter_unit.mcherry_bg = mCherry_bg;
    parameter_unit.dapi_bg = DAPI_bg;
    clear GFP_bg mCherry_bg DAPI_bg
    
    
    % 4f. eccentricity and angle
    ecc = extractfield(sample_particles,'Eccentricity')';
    angle = extractfield(sample_particles,'Orientation')';
    parameter_unit.Ecc = ecc;
    parameter_unit.Angle = angle;
    clear ecc angle
    
    
    % 5. trim particles by (1) area and (2) width
    
    % 5a. trim by area
    TrimField = 'A';    % choose relevant characteristic to restrict, run several times to apply for several fields
    LowerBound = 0.8;   % lower bound for restricted field, or -Inf
    UpperBound = 176;   % whole image length
    p_trim1 = ParticleTrim_glycogen(parameter_unit,TrimField,LowerBound,UpperBound);
    
    
    % 5b. trim by width
    TrimField = 'MinAx';  % choose relevant characteristic to restrict, run several times to apply for several fields
    LowerBound = 0.7;     % lower bound for restricted field, or -Inf
    UpperBound = 176;     % whole image length
    p_trim2 = ParticleTrim_glycogen(p_trim1,TrimField,LowerBound,UpperBound);
    

    
    % 6. store final data 
    converted_data{1,sample} = p_trim2;
    
end
clear sample sample_particles
    
    
%% Part THREE: visualize measured data


% 4. plot mean intensity in cells as a fraction of mean background
for smpl = 1:length(samples)
    
    sample_data = converted_data{1,smpl};
    
    % mean cell fluorescence 
    cell_gfp = sample_data.gfp_cell;
    cell_mcherry = sample_data.mcherry_cell;
    cell_dapi = sample_data.dapi_cell;
    
    % mean bg fluorescence
    bg_gfp = sample_data.gfp_bg;
    bg_mcherry = sample_data.mcherry_bg;
    bg_dapi = sample_data.dapi_bg;
    
    % cell fluorescence normalized by bg fluorescence
    norm_gfp = cell_gfp./bg_gfp;
    norm_mcherry = cell_mcherry./bg_mcherry;
    norm_dapi = cell_dapi./bg_dapi;
    
    % mean of normalized values
    gfp = mean(norm_gfp);
    mch = mean(norm_mcherry);
    dap = mean(norm_dapi);
    
    % std of normalized values
    gfp_std = std(norm_gfp);
    mch_std = std(norm_mcherry);
    dap_std = std(norm_dapi);
    
    % standard error of normalized values
    gfp_n = length(norm_gfp);
    mch_n = length(norm_mcherry);
    dap_n = length(norm_dapi);
    
    gfp_sem = gfp/sqrt(gfp_n);
    mch_sem = mch/sqrt(mch_n);
    dap_sem = dap/sqrt(dap_n);
 
    % plot
    figure(smpl+10)
    bar([gfp, mch, dap])
    hold on
    errorbar([gfp, mch, dap], [gfp_sem, mch_sem, dap_sem],'.')
    set(gca,'xticklabel',{'CFB286','Eub338','LGC354A'})
    title(samples{smpl})
    ylim([1 2.5])
    ylabel('Normalized intensity')
    
    
end





