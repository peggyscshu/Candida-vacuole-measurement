# Candida-vacuole-measurement                          
[![DOI](https://zenodo.org/badge/951821520.svg)](https://doi.org/10.5281/zenodo.15309072) v1.1.3
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15347737.svg)](https://doi.org/10.5281/zenodo.15347737) v1.1.0  

![image](https://github.com/user-attachments/assets/29302f53-1753-4ff3-a76c-433c661d8310)


# Introduction
This is an Image J[1] macro-based analysis tool for identifying Candida-containing vacuoles (CCVs) and quantifying LC3 intensity at vacuole boundaries. 
This tool is designed for the automated detection and quantification of LC3 signal surrounding CCVs in fluorescence microscopy images. It processes multi-channel images, including DAPI, GFP-labeled Candida (strain OG-1), and LC3, and outputs a summary table of the cell number, vacuole number and the LC3 intensity. 

# Examples
Candida-containing vacuoles (CCVs) images, captured by Mr. Shun-Chih Liu, are the courtesy from Dr. Li-Chung Hsu, Institute of Molecular Medicine, National Taiwan University. 
The test image is composed by four channels in the following order—DAPI, GFP-labeled Candida, Actin and LC3. Please use the same order to arrange the input image to avoid the quantification from the wrong channel.
Image 13.lsm is used to test the script version 1.1.0, while Image10.lsm is used for version 1.1.3.

# Version
•  __20250429_Vacuole counting and measure intensity_WithMask.ijm__: 
The intermediate files for all processed images will be saved in the same folder for further data inspection. Nine additional intermediate files will be generated from one raw image. This version secures all the data but increases the requirement for data storage.  

•  __20250430_Vacuole counting and measure intensity_NoMask.ijm__: 
  Only the intermediate files for the last processed image will be saved in the same folder for a quick inspection. 

•  __20250506_Vacuole counting and measure intensity_ver1.1.0.ijm__:   v1.1.0: DOI  [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15347737.svg)](https://doi.org/10.5281/zenodo.15347737)

  Main changes in
   1. Process file format with .tif, .lsm, .czi
   2. Store all intermediate files in another folder created under the parent folder.
   3. Report vacuole no. directly if no LC3 holes inside cells are detected.
   4. Elevate the minimum size requirement for fungal DNA.

•  __20250603_Vacuole counting and measure intensity_ver1.1.3.ijm__: v1.1.3: DOI [![DOI](https://zenodo.org/badge/951821520.svg)](https://doi.org/10.5281/zenodo.15309072)

  Main change in 
  1. Change the channel order from "DAPI, Candida, actin, LC3" to "DAPI, actin, Candia and LC3."
  2. Use the actin channel to identify cell contour and fiilter with the area size.
  3. Candida was defined by the fungal DNA surrounded by Candida GFP. Because fungal DNA were not surrounded by Candida GFP, these regions were defined by the DAPI staining channel.
  4. The cell number was not counted by the DAPI channel to avoid the double nuclei over counting.
  5. Count Candida cell wall inside of cells.
  6. The nucleus number is counted by the DAPI channel.
  7. The cell area is measured from the actin channel.




# Description
1.	__Image Preprocessing__:
Image channels are split and renamed for downstream analysis.
2.	__Cell Boundary Detection__:
The LC3 image in the fourth channel is used to detect the cell boundary. Gaussian blur, thresholding, and convex hull operations are applied to the LC3 channel to approximate cell contours. A binary mask of the cell area is then generated based on these outlines.
3.  __Cell Counting__:
Gaussian blur and Otsu[2] thresholding are applied to the DAPI channel to detect nuclei. The resulting image is multiplied by the cell mask to isolate nuclei within cells, allowing for accurate cell counting.
4.	__LC3 Hole Extraction__:
CLIJ[3]-based tubeness[4] filtering is applied to the LC3 channel to enhance and identify vacuole-like structures.
5.	__Candida Localization__:
The DAPI and Candida channels are combined to localize Candida and nuclei. Nuclear regions are then removed to generate a Candida location mask.
6.	__Candida-Containing Vacuole (CCV) Identification__:
LC3 holes are screened for overlap with Candida localization to identify Candida-containing vacuoles.
7.	__LC3 Intensity Quantification in CCVs__:
The LC3 signal surrounding each identified vacuole is measured to assess LC3 recruitment.
8.	__Result Export__:
 Depending on the selected edition, the results may include intermediate masks and ROIs. In all cases, an Excel file is generated containing summary measurements, including LC3 intensity per vacuole, cell count, and the number of LC3-positive vacuoles.

# Instructions
1.	Install ImageJ/Fiji with CLIJ2 plugin. 
2.	Place all your multi-channel TIFF images in a folder.
3.	According to your needs, choose a version of the macro to download. 
4.	Open the macro in Fiji.
5.	Run the script; It will prompt you to select your folder.
6.	Depending on the version, intermediate files including masks and ROIs will be saved under the same directory or under the parental folder. The measurements for each image will be collected in a summary table.

# Tutorial  
[![YouTube](https://img.youtube.com/vi/GqjaSe0SBtk/0.jpg)](https://youtu.be/GqjaSe0SBtk)

# Workflow Design (v. 1.1.0)
![image](https://github.com/user-attachments/assets/2cfde3d2-89ca-4bc1-8371-514c30491f21)

# Workflow Design (v. 1.1.3)
![image](https://github.com/user-attachments/assets/e41abf09-202c-4f89-965c-4cd20e9ffe96)


# Acknowledgements
We thank Mr.Shun-Chih Liu (劉舜治) and Dr. Li-Chung Hsu(徐立中) for providing the demonstration images used during the development of this workflow.

# Reference
1.	Schindelin, J., Arganda-Carreras, I., Frise, E., et al. (2012). Fiji: an open-source platform for biological-image analysis. Nature Methods, 9(7), 676–682. https://doi.org/10.1038/nmeth.2019
2.	Otsu, N. (1979). A threshold selection method from gray-level histograms. IEEE Transactions on Systems, Man, and Cybernetics, 9(1), 62–66. https://doi.org/10.1109/TSMC.1979.4310076
3.	Haase, R., Royer, L. A., Steinbach, P., et al. (2020). CLIJ: GPU-accelerated image processing for everyone. Nature Methods, 17, 5–6. https://doi.org/10.1038/s41592-019-0650-1
4.	Sato, Y., Nakajima, S., Shiraga, N., Atsumi, H., Yoshida, S., Koller, T., Gerig, G., & Kikinis, R. (1998). Three-dimensional multi-scale line filter for segmentation and visualization of curvilinear structures in medical images. Medical Image Analysis, 2(2), 143–168.
5.	Phansalkar, N., More, S., Sabale, A., & Joshi, M. (2011). Adaptive local thresholding for detection of nuclei in diversity stained cytology images. In 2011 International Conference on Communications and Signal Processing (pp. 218–220). IEEE. https://doi.org/10.1109/ICCSP.2011.5739305



