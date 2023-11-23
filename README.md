# hCoCena - Horizontal integration and analysis of transcriptomics datasets [[paper](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btac589/6677225)]

hCoCena is an R-package that allows you to integrate and jointly analyse multiple transcriptomic datasets or simply analyse a single dataset if you don't want to do any data integration! hCoCena uses network representations of the data as the basis for integration. You can find more details of how that works in our [paper](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btac589/6677225) . Below, you will find some info on how to install the package and tips for using it. 

## Installation
To install hcocena from this repo, run the codeline provided in the `install_hcocena.R` script.
To install versioned dependencies, use the script `install_versioned_dependecies.R`.

## Usage
**hCoCena is divided into 2 parts:** 

**1.** the main analysis that comprises the mandatory steps to process and integrate the data and

**2.** the satellite functions that offer you a plethora of analysis options in a pick & mix kind of fashion. 

The figure below illustrates this: the main analysis is at the center, while the satellite functions can be found in the orbits around it. 
A step-by-step walkthrough of the main analysis steps can be found in the `hcocena_main.Rmd`, the satellite functions are in the `hcocena_saltellite.Rmd`. 

hCoCena was written with user-friendliness and customizability in mind. We are doing our best to provide you with plenty of supplementary information that make the usage of the tool easy for you. You can also always extend the tool's functionalities with your on custom scripts and functions to adapt the analysis to your needs! For more details on hCoCena's object structure and where to find the outputs of different analysis steps for customization, please refer to the overview in the [Wiki](https://github.com/MarieOestreich/hCoCena/wiki/Structure-of-the-hcobject) and the extensive function documentations you can access from within R Studio.


![hCoCenaFig1](https://user-images.githubusercontent.com/50077786/158609782-2048c06e-0420-4c3f-8680-5d99f91d6905.jpg)
*Marie Oestreich, Lisa Holsten, Shobhit Agrawal, Kilian Dahm, Philipp Koch, Han Jin, Matthias Becker, Thomas Ulas, hCoCena: horizontal integration and analysis of transcriptomics datasets, Bioinformatics, Volume 38, Issue 20, 15 October 2022, Pages 4727–4734, https://doi.org/10.1093/bioinformatics/btac589*

## Showcase
To rerun the showcase example from our original publication, please refer to the branch of version [1.0.1](https://github.com/MarieOestreich/hCoCena/tree/v-1.0.1).

## Wiki
For loads of additional information regarding the [satellite functions](https://github.com/MarieOestreich/hCoCena/wiki/Satellite-Functions), [community detection](https://github.com/MarieOestreich/hCoCena/wiki/Background-Info-on-the-Community-Detection-Algorithms) algorithms etc. please check out our carefully curated Wiki pages!
