## About this project
This project contains an R script and associated data files to generate all figures and Tables 3 and 4 presented in Nishiaki et al. (2021, Quaternary International). All codes were written by Yutaka Kobayashi, the corresponding author of the paper.

## Main R script
Run "PaleoAsiaDB2021QI.R" with R Studio.

## Associated files
The main script refers to a dataset file (PADBMode200801_excluding_Africa.csv) as well as shape files included in four separate directories. The shape files were obtained from the website of Natural Earth; I remark that those shape files are all in the public domain and free to disseminate. The dataset file and the directories of the shape files must be located in the same directory as this source file when it is executed.
 
## Figures
To plot each figure, one needs to uncomment a specific part of the main script. One can find it searching for the word "Figure". Note that one may need to adjust the size of each figure manually because the figures are tailored for publication, not for displaying on a monitor.
 
## Tables
Tables 3 and 4 are automatically output as TeX files to the same directory as the main script.

## Programming environment
The program has been developed mainly using R ver 3.6.3 and RStudio ver 1.3.1093. It should work with similar versions of R and RStudio.

## Packages
The following packages are required:
 
+ ggplot2:
H. Wickham. ggplot2: Elegent Graphics for Data Analysis. Springer-Verlag New York, 2016.

+ geosphere:
Robert J. Hijmans (2019). geosphere: Spherical Trigonometry. R package version 1.5-10. https://CRAN.R-project.org/package=geosphere

+ philentropy:
Robert J. Hijmans (2019). geosphere: Spherical Trigonometry. R package version 1.5-10. https://CRAN.R-project.org/package=geosphere

+ ggrepel:
Drost HG.  Philentropy: Information Theory and Distance Quantification with R. Journal of Open Source Software (2018). doi:10.21105/joss.00765

+ sf:
Pebesma, E., 2018. Simple Features for R: Standardized Support for Spatial Vector Data. The R Journal 10 (1), 439-446, https://doi.org/10.32614/RJ-2018-009

+ xtable:
David B. Dahl, David Scott, Charles Roosen, Arni Magnusson and Jonathan Swinton (2019). xtable: Export Tables to LaTeX or HTML. R package version 1.8-4. https://CRAN.R-project.org/package=xtable

+ gridExtra:
Baptiste Auguie (2017). gridExtra: Miscellaneous Functions for "Grid" Graphics. R package version 2.3. https://CRAN.R-project.org/package=gridExtra

+ scales:
Hadley Wickham and Dana Seidel (2019). scales: Scale Functions for Visualization. R package version 1.1.0. https://CRAN.R-project.org/package=scales

## Contact info
If any questions as to this project, please send an email to kobayashi.yutaka@kochi-tech.ac.jp

<div style="text-align: right;">
Yutaka Kobayashi, Feb. 14, 2021.
</div>
