# summero-strains

A collection of scripts analyzing FISH imaging data

A work in progress by Jen Nguyen

## Summary of Scripts 
**1. enterMetaData.m:** creates and/or builds from data structure called metadata.mat to store experiment conditions 
- experiment date
- magnification
- sample names
- sample type 
	- "pure": fixed cultures of known strains
	- "mixed": undefined microbial community
	- "fecal": from fecal sample
- number of channels
- growth stage ('lag','exponential','stationary','mixed','fecal')

**2. whos_a_cell.m:** threshold determination (by cell width) to distinguish single cells from other particles. User is prompted to enter and validate threshold values for each strain. Paramaters are stored in data structure called segdata.mat
- These values are important to inform intepretations of signal intensity (for example, clumps have more autofluorescence than single cells)

**3. intensity_cell_v_bg.m:** plots mean per particle fluorescense intensity normalized by background intensity 
- Does not yet account for single vs clumps of cells

**4. segmentIntensity_GFP.m:** uses threshold paramaters in segdata.mat to quantify fluorescence intensity in single cells and clumps

**5. segmentIntensity_3channels.m:** uses threshold pramatares in segdata.mat to quantify and plot fluorescene intensity in single cells normalized by background intensity for experiments with three fluorecent channels

**6. segmentIntensity_GFP.m:** uses threshold paramaters in segdata.mat to quantify and plot fluorescence intensity in single cells sepearately for each image
- used for troubleshooting purposes or to find outliers in the imaging data 

## Workflow 

1. For each new experiment, use  **enterMetaData.m** to determine and store experiment conditions 
2. For each new experiment, use **whos_a_cell.m** to determine and store threshold values to distinguish single cells 
3. Use the appropriate segmentIntensity script to calculate and plot fluorescene intensity and plot data 
