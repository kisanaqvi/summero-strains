%% FISH n flow distributions

% goal: visualize distributions of event fluorescence

% strategy:
%
%       1. load data tables
%       2. plot histograms
%       3. save


% okie, let's go!

% last updated: jen, 2021 April 20
% commit: histograms for 2021-04-20 flow experiment

%% 

% 0. initialize experiment data
clc
clear
cd('/Users/jen/Documents/TropiniLab/Data/such-hipr/sourcedata')
load('metadata.mat')

index = 3; % 2021-04-20
date = metadata{index}.date;
samples = metadata{index}.samples; % not yet streamlined with flow data nomenclature

data_folder = strcat('/Users/jen/Documents/TropiniLab/Molecular_tools/HiPR_fish/',date);
cd(data_folder)
cd('Exp_20210420_2') % Experiment data from cytoflex


% for each sample of interest
soi = {'s1','s1_ii','s2','s2_ii','s3','s3_ii','s4','s4_ii_diluted'};
for ss = 1:4
    
    % 1. load data tables
    table = readtable(strcat('export_',soi{ss},'.csv'));
    fitc_a = table2array(table(:,6));
    fitc_h = table2array(table(:,7));
    
    % 2. plot histograms
    figure(1)
    histogram(fitc_a)
    %hold on
    %histogram(fitc_h)
    legend('fitc-a','fitc-h')
    title(soi{ss})
    ylabel('counts')
    xlabel('value')
    
    % 3. save
    figure(1)
    saveas(gcf,strcat(date,'-',soi{ss},'-histogram'),'epsc')
    close(gcf)
    
end
