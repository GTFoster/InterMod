if (!"remotes"%in%installed.packages()){install.packages("remotes")}
remotes::install_github("cran/rgeos")
remotes::install_github("cran/rgdal")
remotes::install_github("andrefaa/ENMTML")  #Installation of ENTML

ToucanGBIF <- ENMTML::ENMTML(pred_dir="data/wc2.1_country/",
                             proj_dir = NULL,
                             result_dir = NULL,
                             occ_file = "data/occ_Dat/Toucan_gbif.txt",
                             sp="species",
                             x="decimalLongitude",
                             y="decimalLatitude",
                             min_occ=20,
                             thin_occ=NULL,
                             eval_occ=NULL,
                             colin_var = c(method='PCA'),
                             imp_var = TRUE,
                             sp_accessible_area = NULL,
                             pseudoabs_method=c(method='GEO_CONST', width='100'), #100 km buffer is what cleber does
                             pres_abs_ratio = 0.5,
                             part=c(method= 'KFOLD', folds='2'), #Kfolds valifation, where K=5
                             save_part = FALSE,
                             save_final = TRUE,
                             algorithm="MXD", #Using Maxent with Default Params
                             thr="MAX_TSS",#Threshold by maximizing TSS
                             msdm = NULL,
                             ensemble = NULL,
                             extrapolation = FALSE,
                             cores = 14)