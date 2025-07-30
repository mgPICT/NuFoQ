# NuFoQ
Yeast foci analysis from fluorescence microscopy images.

Created by Mickaël Garnier (mickael.garnier@curie.fr) on June 15, 2021.

Following the work done by the PICT@Pasteur and the UMR 3664 Taddei lab.

# Script purpose
The NuFoQ matlab script performs quantitative analysis of foci in different yeast strains and/or conditions. A few parameters allow the users to target specific spots inside the nuclei and adapt to their image quality and specificities. A second step can be used to filter out spots and nuclei that are badly segmented or should not be  included in the following analysis. Finally NuFoQ can generate descriptive graphs to represent the data.

# Details
Operating System: tested on Debian 11 and 12, Windows 7 and 10, Max OS 10.15.

MATLAB: Requires version 2020b or newer.

# Installation
Simply put all the files in a directory with the utilities folder next to them.

# Start the script
Prepare your image data:
  - have all the images inside a single folder (another folder for deconvolved images if needed),
  - the image names in both folders should be the same,
  - the yeast strain and any other relevant condition information should appear in the image names.

Start the NuFoQ script by running the NuFoQ.m script for the editor from its folder. The user interface should appear shortly after.

<img width="514" height="440" alt="maingui" src="https://github.com/user-attachments/assets/ab0ef0f6-5b39-478b-8f51-c1cb354efcc6" />

# Set parameters
Choose wether you want to analyse a single file or a complete folder.

The option *Use deconvolution images* will have you select two folders (or images):
  - first for quantification, usually the raw data
  - then for segmentation support, usually deconvolved images

Check the *Save segmented images* if you want to be able to supervise the segmentation outcomes.

The *Set parameters* button opens a new dialog box with the following entries to fill:
  - _XY resolution (µm)_ height and width of a pixel in micrometers,
  - *Z resolution (µm)* depth of a voxel in micrometers,
  - *Foci size (pixels)* approximate dimension of the targeted foci, it can be a MATLAB range like 1:2:5,
  - *Z-slices to skep* number of Z slices fomr the top of the stack to exclude from the analysis,
  - *Signal name* indicator of the channel to quantify, the script will search these characters in the files to perform the analysis,
  - *Low signal* flag indicating very bad signal to noise ratio and will remove some denoising steps,
  - *Number of std above the mean for spots* filter to remove peaks with an intensity lower than the average of their nucleus plus a certain amount of the nucleus intensity standard deviation. This value can be any rational number.

<img width="255" height="452" alt="setparams" src="https://github.com/user-attachments/assets/061c8010-e3c4-4b25-9f7e-b380dc311463" />


# Check and supervise segmentation outputs
<img width="2442" height="1402" alt="check" src="https://github.com/user-attachments/assets/06a87fca-557e-4dfe-baaf-0c5591aafc21" />

# Output files
- Analysis:
  - *filename*.nucs: all quantification of nuclei,
  - *filename*.spots: all quantification of foci, includes their respective nuclei measurements,
  - *filename*_segNuc.tif: segmentation mask of nuclei with labels (optional),
  - *filename*_segSpots.tif: segmentation mask of foci with labels (optional).

- Supervision:
  - *filename*.nucs: updates the *To keep* flag,
  - *filename*.spots: updates the *To keep* flag.

- Graphs
  - *all_measurements.txt*: exhaustive table including all foci quantifications with image names and unique labels for both foci and nuclei,
  - several boxplots, histograms and probability density functions.

