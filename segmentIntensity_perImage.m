%% segmentIntensity_perImage


% goal: boxplots of single-cell intensity divided by background intensity
%       per image

%  - uses data matrix saved from from Part One of segmentIntensity_3channels
%  - makes one figure per condition, with subplots being different images
%    per condition


% ok, let's go!

% last updated: Kisa, 2021 June 23
% commit: plotted data from 2021-06-22 experiment


%% cell vs background intensity per image


clear
clc
%cd('/Users/jen/summero-strains')
cd('C:/Users/Kisa Naqvi/Documents/TropiniLab/summero-strains-master')
load('metadata.mat')



% 0. initialize experiment data
index = 6; % 2021-06-22
date = metadata{index}.date;
%cd('/Users/jen/Documents/TropiniLab/Data/Kisa')
cd('C:/Users/Kisa Naqvi/Documents/TropiniLab/Data')
load(strcat('dm-segmentIntensity-',date,'.mat'))

samples = metadata{index}.samples;
magnification = metadata{1}.magnification;
px_size = 11/magnification; % 11 um pixels with 150x magnification


% 1. convert measurements to microns based on imaging specifications

% plot boxplot per image per sample
for sample = 1:length(samples)
    
    sample_images = dm(:,sample);
    
    for image = 1:length(sample_images)
        
        img_particles = sample_images{image,1};
        
        if isempty(img_particles) == 1
            continue
        end
        
        
        % 1a. convert x,y coordinate of particle centroid
        x_position = [];
        y_position = [];
        for ii = 1:length(img_particles)
            
            centroid = img_particles(ii).Centroid.*px_size;
            x_position = [x_position; centroid(1)];
            y_position = [y_position; centroid(2)];
            
        end
        parameter_unit.X = x_position;
        parameter_unit.Y = y_position;
        clear particle x_position y_position centroid ii
        
        
        % 1b. convert area
        area = extractfield(img_particles,'Area')';
        parameter_unit.A = area.*px_size^2;
        clear area
        
        
        % 1c. major & minor axis
        majorAxis = extractfield(img_particles,'MajorAxisLength')';
        minorAxis = extractfield(img_particles,'MinorAxisLength')';
        parameter_unit.MajAx = majorAxis.*px_size;
        parameter_unit.MinAx = minorAxis.*px_size;
        clear majorAxis minorAxis
        
        
        % 1d. cell intensity
        GFP_cell = extractfield(img_particles,'gfp_cell')';
        mCherry_cell = extractfield(img_particles,'mcherry_cell')';
        DAPI_cell = extractfield(img_particles,'dapi_cell')';
        parameter_unit.gfp_cell = GFP_cell;
        parameter_unit.mcherry_cell = mCherry_cell;
        parameter_unit.dapi_cell = DAPI_cell;
        clear GFP_cell mCherry_cell DAPI_cell
        clear GFP_cell
        
        
        % 1e. background intensity of corresponding image
        GFP_bg = extractfield(img_particles,'gfp_bg')';
        mCherry_bg = extractfield(img_particles,'mcherry_bg')';
        DAPI_bg = extractfield(img_particles,'dapi_bg')';
        parameter_unit.gfp_bg = GFP_bg;
        parameter_unit.mcherry_bg = mCherry_bg;
        parameter_unit.dapi_bg = DAPI_bg;
        clear GFP_bg mCherry_bg DAPI_bg
        
        
        % 1f. eccentricity and angle
        ecc = extractfield(img_particles,'Eccentricity')';
        angle = extractfield(img_particles,'Orientation')';
        parameter_unit.Ecc = ecc;
        parameter_unit.Angle = angle;
        clear ecc angle
        
        
        % 2. trim particles by width
        %    values are set as recorded in whos_a_cell.m
        
        % 2b. trim by width
        TrimField = 'MinAx';  % choose relevant characteristic to restrict, run several times to apply for several fields
        if image < 3
            LowerBound = 0.6;       % bounds for stationary G6 (see whos_a_cell.m)
            UpperBound = 0.9;
        elseif image == 3
            LowerBound = 0.7;       % bounds for stationary H03 + H06 (see whos_a_cell.m)
            UpperBound = 1.6;
        end
        p_trim = ParticleTrim_glycogen(parameter_unit,TrimField,LowerBound,UpperBound);
        
        
 
        
        % 3. plotting!
        
        % 3a. gather cell width and intensity data
        %sample_data = converted_data{1,sample};
        cell_width = p_trim.MinAx;
        
        % mean cell fluorescence
        cell_gfp = p_trim.gfp_cell;
        cell_mcherry = p_trim.mcherry_cell;
        cell_dapi = p_trim.dapi_cell;
        
        
        % mean bg fluorescence
        bg_gfp = p_trim.gfp_bg;
        bg_mcherry = p_trim.mcherry_bg;
        bg_dapi = p_trim.dapi_bg;
        
        
        % 3b. isolate single cells from clumps
        if sample < 3
            clumpThresh = 0.9; % min width of clumps for S1 & S2 (see whos_a_cell.m)
        elseif sample == 3
            clumpThresh = 1.6;
        end
        
        single_gfp = cell_gfp(cell_width <= clumpThresh);
        single_mcherry = cell_mcherry(cell_width <= clumpThresh);
        single_dapi = cell_dapi(cell_width <= clumpThresh);
        single_bg_gfp = bg_gfp(cell_width <= clumpThresh);
        single_bg_mcherry = bg_mcherry(cell_width <= clumpThresh);
        single_bg_dapi = bg_dapi(cell_width <= clumpThresh);
        
        norm_single_gfp = single_gfp./single_bg_gfp;
        norm_single_mcherry = single_mcherry./single_bg_mcherry;
        norm_single_dapi = single_dapi./single_bg_dapi;
        
        figure(sample)
        subplot(1,length(sample_images),image)
        xx = [norm_single_gfp; norm_single_mcherry];%; norm_single_dapi];
        gg = [zeros(length(norm_single_gfp), 1); ones(length(norm_single_mcherry), 1)];% 2*ones(length(norm_single_dapi), 1)];
        boxplot(xx,gg)
        set(gca,'xticklabel',{'new','old'})%,'dapi'})
        title(strcat(samples{sample},', n =',num2str(length(norm_single_gfp))))
        ylim([0.8 8])
        
    end
    
end





