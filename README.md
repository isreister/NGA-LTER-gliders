# NGA-LTER-gliders
This repo holds all the code that was used to process gliders deployed in the northern gulf of alaska during the timespan of 2021 to 2024. This code is what PhD candidate Isaac Reister used to standardize variable names, identify upcasts and downcasts, and split upcasts and downcasts up, then complete various quality control measures. 

List number corresponds to what script should be run first, second, and so on.
	File name										Main task
1) Data ProcessingGlider_AllGliderVariablesMaker2	  --------------- Standardizes the glider variable names and applies some ranged QC.
2) Data ProcessingGlider_Profilesplitterv4		  --------------- splits up the glider into upcasts and downcasts.
3) DataProcessingGlider_QCglider_and_verticalgridding	  --------------- runs a rolling QC window. Produces a vertically gridded product of upcasts and downcasts. 
4) DataProcessingGlider_more_QCglider_metrics_finalgrid   --------------- runs a few more QC things, makes a final gridded product, and saves out prelim plots of the gridded data

Functions that are run in these scripts are in the glider function folder. As we make adjustments to the code, this will also serve as a good place to put additional functions.

              
