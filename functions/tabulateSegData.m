% tabulate segmentation data

function segTable = tabulateSegData(segdata)

% 1. gather field names
fnames = fieldnames(segdata{1});

% 2. loop through each index in segdata and build matrix, with each column
%    representing one field and each row representing one index
segTable = cell(length(segdata),length(fnames));
for idx = 1:length(segdata)
    
    idxdata = segdata{idx};
    segTable{idx,1} = idxdata.sample_strain;
    segTable{idx,2} = idxdata.sample_stage;
    segTable{idx,3} = idxdata.sample_date;
    segTable{idx,4} = idxdata.test_experiment;
    segTable{idx,5} = idxdata.tested_img;
    segTable{idx,6} = idxdata.minWidth;
    segTable{idx,7} = idxdata.maxWidth;
    
end

end