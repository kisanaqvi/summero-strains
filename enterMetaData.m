%% enterMetaData.m

% goal: prompt user for inputs and store to experiment specific structures
%       stored data:
%               1. experiment date
%               2. magnification
%               3. sample names
%               4. sample type (pure culture, mixed culture, fecal sample)
%               5. number of channels
%               6. growth stage ('lag','exponential','stationary','mixed','fecal')
%
%       each row is a different experimental "replicate"
%       each column represents a type of image data that requires the same analysis
%
%            column 1 = pure cultures
%            column 2 = mixed samples
%            column 3 = fecal samples
%


% last updated: 2021 April 20
% commit message: edit location of metadata.mat & add "mixed" as growth
%                 stage option


% OK let's go!


%% 0. initialize dimensions of current data structure

clc
clear
cd('/Users/jen/Documents/TropiniLab/Data/such-hipr/sourcedata')
load('metadata.mat')


%% 1. enter meta data for a new experiment

% 1. prompt user for sample type
prompt = 'Enter experiment type as a string (pure/mixed/fecal): ';
sampleType = input(prompt);
newdata(1).sampleType = sampleType;


% 2. prompt user for number of channels data
prompt = 'Enter number of fluorescent channels (i.e. exclude phase): ';
numFluors = input(prompt);
newdata(1).numFluors = numFluors;



% 3. determine column of experiment to add
if strcmp(sampleType,'pure') == 1 
    column = 1;                       
elseif strcmp(sampleType,'mixed') == 1
    column = 2;
elseif strcmp(sampleType,'fecal') == 1
    column = 3;
else
    disp('Error: unknown sampleType!');
end


% 4. prompt user for experiment date
prompt = 'Enter experiment date as a string: ';
date = input(prompt);
newdata(1).date = date;


% 5. prompt user for growth stage
prompt = 'Enter growth stage as a string (lag,exponential,stationary,mixed,fecal): ';
growthStage = input(prompt);
newdata(1).growthStage = growthStage;

%% still to finish!

% 6. prompt user for sample names
prompt = 'Enter sample names as a cell array of strings ({name1,name2,...}): ';
samples = input(prompt);
newdata(1).samples = samples;

%%

% 7. prompt user for strain info
prompt = 'Enter strains as a cell array of strings ({strain1,strain2,...}): ';
strains = input(prompt);
newdata(1).strains = strains;


% 8. prompt user for sample magnification
prompt = 'Enter sample magnification as an integer (100, 150,...): ';
magnification = input(prompt);
newdata(1).magnification = magnification;


% 9. assign data structure to new (experiment-specific cell)
priorEntries = size(metadata);
if column > priorEntries % new column, row = 1
    newrow = 1;
else % column already exists!
    columndata = metadata(:,column);
    priorRows = sum(~cellfun(@isempty,columndata));
    newrow = priorRows + 1;
end
metadata{newrow,column} = newdata;


%% 10. save storedMetaData
save('metadata.mat','metadata')

