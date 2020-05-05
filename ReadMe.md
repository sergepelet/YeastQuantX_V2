#  Read Me YeastQuant Version 2

YeastQuant is an image analysis platform designed to extract single cell measurements from time-lapse movies.

Experiments are documented in a FileMaker database and the image analysis is performed with Matlab.

The initial description of the platform has been published in *Integrative Biology*: 
Pelet, S., Dechant, R., Lee, S.S., van Drogen, F., and Peter, M. (2012). An integrated image analysis platform to quantify signal transduction in single cells. **Integr. Biol.** *4*, 1274. [Integrative Biology] (https://doi.org/10.1039/c2ib20139a), [ETHZ repository] ()

Since this first published version, multiple changes have been made:
 - The general flow of the image analysis has been changed to facilite the parallel processing of the analysis on a cluster.
 - The general structure of the output have been simplified.
 - A new segmentation method based on an image of fluorescent nuclei and two bright-field images has been implemented.
 - New secondary objects can been generated. 
These modification are explained in YeastQuantX_info.pdf file.