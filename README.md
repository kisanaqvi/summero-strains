# such-hipr

A collection of scripts analyzing HiPR-FISH imaging data

A work in progress by Jen Nguyen

## determining single cells from junk and clumps
1. whos_a_cell.m: threshold determination (by cell width) to distinguish single cells from other particles
	- a different section of code for each strain, enables strain specific image processing
	- these values are important to inform intepretations of signal intensity
	- for example, clumps have more autofluorescence than single cells

## measuring signal intensity
2. intensity_cell_v_bg: plots mean per particle fluor int normalized by background int
	- does not yet account for single vs clumps of cells
3. segmentIntensity_GFP: uses whos_a_cell to quantify fluorescence intensity in single cells and clumps
