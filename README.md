
Start here!

Author: **Dr Miguel Lurgi** 
Computational Ecology Lab. 
Swansea University, UK. 

Centre for Biodiversity Theory and Modelling. 
Theoretical and Experimental Ecology Station, CNRS, France. 

Date Created: 19-12-2019

Copyright (c) Miguel Lurgi, 2019-2021. 
email: miguel.lurgi@swansea.ac.uk

This repository contains the data and source code used to produce the results
presented in our paper entitled:

Lurgi, M., Galiana, N., Broitman, B. R., Kéfi, S., Wieters, E. A., and Navarrete, S. A. (2020)
**Geographical variation of multiplex ecological networks in marine intertidal communities**. 
***Ecology*** 101(11):e03165. https://doi.org/10.1002/ecy.3165

In this file I describe the contents of the archive (i.e. files contained here) and provide
a brief description of what each data file contains and what each source file does when executed.

DATA FILES:
Data files contained in this archive contain empirical data collected by researchers at the 
Coastal Marine Research Station of the Pontifical Catholic University of Chile over the last few
years. These data contain information on the presence and abundance of different species in the rocky shore
coasts of Chile, and also on different environmental factors. Some data comes from previous publications
by others and are included here for completeness and to facilitate the execution of the source code.

1.- 'chilean-rocky-shores-samples.csv' : This is the main database files that contain the records
of presence and abundance (plus some extra metadata) of each species in the marine intertidal rocky
shores at each location along the latitudinal spatial extent of the central coast of Chile

2.- 'sites-with-coordinates.csv' : This table contains the lat-lon coordinates and codes of each of the
sampling sites used in this study

3.- 'chilean_metadata.csv' : This file is the original metadata file of the multiplex network of ecological
interactions from the paper Kéfi et al. (2015) Network structure beyond food webs: Mapping non-trophic and 
trophic interactions on Chilean rocky shores, Ecology.

4.- 'chilean_TI.txt' : This file contains the original trophic interaction network from the paper 
Kéfi et al. (2015) Network structure beyond food webs: Mapping non-trophic and trophic interactions 
on Chilean rocky shores, Ecology.

5.- 'chilean_NTIneg.txt' : This file contains the original non-trophic negative interaction network from the paper 
Kéfi et al. (2015) Network structure beyond food webs: Mapping non-trophic and trophic interactions 
on Chilean rocky shores, Ecology.

6.- 'chilean_NTIpos.txt' : This file contains the original non-trophic positive interaction network from the paper 
Kéfi et al. (2015) Network structure beyond food webs: Mapping non-trophic and trophic interactions 
on Chilean rocky shores, Ecology.

7.- 'species-dict.csv' : This file contains a dictionary that links the species scientific names as they appear
in the multiplex network and the codes used to identify the species in the database of species occurrences contained
in 'chilean-rocky-shores-samples.csv'

8.- 'avhrr_environmental.csv' : This data file contains the values of three of the environmental variables for each
of the study sites considered in this study: (1) long term sea surface temperature (LT SST), (2) fraction of the
yearly variability in SST not explained by seasonality (Fr Annual), and (3) climatology (Clim).

9.- 'Data_Fig5b_Tapia_et_al_2009.csv' : This contains the data for presented in Figure 5b of the paper Tapia et al.
(2009) Thermal indices of upwelling effects on inner-shelf habitats, Progress in Oceanography; and used here to 
obtain upwelling values using cooling rate (from that figure) as a guidance. See Methods section in the manuscript 
for details.

10.- 'data-for-chilean-analysis.rda' : This R data file contains 4 data structures needed to run the structural
equation models. Three data frames contain the measured properties of each layer of the multiplex interaction 
network at each spatial location and the fourth one contains the environmental variables used as predictors
of these network properties.

In addition to these data files, an extra data file that must be downloaded from the NOAA website: 
https://coastwatch.noaa.gov/ is needed to run the script 'get-environmental-variables.r'. See details within the
script. The file, which contains satellite image information on sea surface temperature (SST) is too big to be 
included here.


SOURCE CODE FILES:
There are three main executable scripts contained here:

1.- 'building-networks.r' : Run this script to obtain each local network at each one of the study sites based on the 
species occurrences and the regional network of species interactions. See Methods in the manuscript for details on how
this procedure works. The code in this script contains the implementation of the network inference based on empirical
data as explained in the paper.

2.- 'get-environmental-variables.r' : In this script the routines to obtain the different environmental variables used
in this study are implemented. It contains the algorithm used to calculate upwelling index from satellite imagery as
explained in the methods.

3.- 'SEMs.r' : This script contains all the statistical analyses performed as part of the structural equation models,
including their fits.

4.- 'rda-procrustes.r' : This script presents the redundancy analyses (RDA) performed over the qualitative and
quantitative versions of the networks constrained by the environmental predictors and their comparisons using
procrustes and protest analyses. This tests were performed to analyse whether qualitative and quantitative
versions of the networks are differentially affected by environmental conditions.

In addition to the three scripts mentioned above, to other source code files containing helper functions are needed to 
support some calculations: auxiliary-functions.r, which contains modified functions for miscellaneous purposes, and 
utils.r, which contains source code to calculate different properties of the networks of ecological interactions.

I hope you enjoy the code!

If you have any questions / run into any issues, please contact me at miguel.lurgi@swansea.ac.uk






