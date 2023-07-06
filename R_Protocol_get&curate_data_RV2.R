################################################################################################################################
#         BIODIVERSITY ONLINE DATABASES: AN APPLIED R PROTOCOL TO GET AND CURATESPATIAL AND CLIMATIC DATA
################################################################################################################################

#####################################################################################
# Authors:  Marina Coca-de-la-Iglesia, Nagore G. Medina, Virginia Valcárcel
# Date: 10/11/2022
# GitHub

# PROTOCOL USED TO ELABORATE FOLLOWING PUBLICATIONS:
# PAPER pre-print: Marina Coca de la Iglesia, Nagore G. Medina, Jun Wen, Virginia Valcárcel (2022). Tropical-temperate dichotomy falls apart in the Asian Palmate Group of Araliaceae. bioRxiv 2021.10.20.465102; doi: https://doi.org/10.1101/2021.10.20.465102
# PAPER: Coca-de-la-Iglesia, M., Medina, N. G., Wen, J., and Valcárcel, V. 2022. Evaluation of tropical–temperate transitions: an example of climatic characterization in the Asian Palmate group of Araliaceae. American Journal of Botany 109( 9): 1488– 1507. https://doi.org/10.1002/ajb2.16059
# DATABASE: Coca de la Iglesia, Marina, Medina, Nagore G., Wen, Jun, & Valcárcel, Virginia. (2021). Spatial and climatic worldwide database of the Asian Palmate Group of Araliaceae (1.0.0) [Data set]. Zenodo. https://doi.org/10.5281/zenodo.5578149

# PROTOCOL RUN ENTIRELY IN THE FOLLOWING OPERATING SYSTEMS AND VERSIONS OF R AND RStudio:
# MacOS Monterey v12.6.7. (Intel) / (Silicon M1) / R v4.3.1 / Rstudio 2023.06.0+421
# Windows® 10 education / R v4.1.2 / Rstudio 2022.12.00
# Kubuntu 22.04 / R v4.3.0 


###############################################
#               STEPS PROCEDURE
###############################################
#### A.  BUILT A CHECKLIST OF TAXA NATIVE RANGE####
  # Create a txt file with the names of all taxa and the country codes according to TDWG standard and save it with the name "Natural_Distribution_Checklist_TDWG.txt"

#### B.	CREATE AN ACCOUNT IN GBIF DATABASE ####
  # https://www.gbif.org/

#### C.  INITIAL PREPARATION ####
  # 2. Set the paths for the input data and the output of results 

  # run getwd() to check the working directory
  # To change the working directory uncomment the following line
  # and replace path with the path of the main folder in your computer.
  # setwd("C:/Users/User/Documents/") 
  
  if(!dir.exists("input")){dir.create("input")}
  path.input <- paste0(getwd(), "/input/")
  
  # path.output is the folder with the results of R code
  
  # Otherwise you can create the output forlder in the working directory (works with any OS)
  if(!dir.exists("output")){dir.create("output")}
  path.output <- paste0(getwd(), "/output/")
  
  
  # 3a. Install R packages (run only the first time, after that comment these lines)
  install.packages("BIEN")          # Functions used: BIEN_occurrence_genus()
  install.packages("countrycode")   # Functions used: countrycode()
  install.packages("data.table")    # Functions used: data.table()
  install.packages("dplyr")         # Functions used: sapply(), filter(), select()
  install.packages("plyr")          # Functions used: ldply()
  install.packages("sf")            # Functions used: st_read()
  install.packages("raster")        # Functions used: crs(), getData(),crop()
  install.packages("readr")         # Functions used: readLines()
  install.packages("rgbif")         # Functions used: occ_search()
  install.packages("rgdal")         # Functions used: readOGR()
  install.packages("spocc")         # Functions used: occ()
  install.packages("spThin")        # Functions used: thin()
  install.packages("remotes")       # Functions used: install_github(). Alternatively you can install library "devtools" 
  install.packages("tidyr")			    # Functions used: separate()
  remotes::install_github("SEEG-Oxford/seegSDM")   # Functions used: nearestLand(). Alternativerly you can use devtoosl::install_github("SEEG-Oxford/seegSDM")
  remotes::install_github("barnabywalker/kewr")    # Functions used: search_ipni(). Alternativerly you can use devtoosl::install_github("barnabywalker/kewr")
  
  # 3b. Load R packages
  library(BIEN)          # Functions used: BIEN_occurrence_genus()
  library(countrycode)   # Functions used: countrycode()
  library(data.table)    # Functions used: data.table()
  library(dplyr)         # Functions used: sapply(), filter(), select()
  library(plyr)          # Functions used: ldply()
  library(sf)            # Functions used: st_read()
  library(raster)        # Functions used: crs(), getData(),crop()
  library(readr)         # Functions used: readLines()
  library(rgbif)         # Functions used: occ_search()
  library(rgdal)         # Functions used: readOGR()
  library(spocc)         # Functions used: occ()
  library(spThin)        # Functions used: thin()
  library(remotes)       # Functions used: install_github(). Alternatively you can call library devtools
  library(seegSDM)       # Functions used: nearestLand()
  library(tidyr)         # Functions used: separate()
  library(kewr)          # Functions used: search_ipni()
  
  # 4. Read names list or create names list. For this protocol, we used the "Alternative 2".
  # Notice that names of taxa must match those included in the "Natural_Distribution_Checklist_TDWG.txt" file
  # and must be arranged in the same order as in the "Natural_Distribution_Checklist_TDWG.txt" file.
  # "Alternative 1":
  # taxa.names <- c("Brassaiopsis","Chengiopanax","Dendropanax","Eleutherococcus","Fatsia","Gamblea",
  #               "Heteropanax","Kalopanax","Macropanax","Merrilliopanax","Metapanax","Oplopanax","Oreopanax",
  #               "Sinopanax","Tetrapanax","Trevesia")

  # "Alternative 2":
  # read the file from your computer
  # taxa.names <-  sort(readLines(paste0(path.input, "namelist.csv"))) 
  
  # "Alternative 3":
  # read the file from the github page
  taxa.names <- sort(readLines("https://raw.githubusercontent.com/NiDEvA/R-protocols/main/namelist.csv"))
    
  # 5. Read the checklist of taxa native range built in step A. 
  
  # Load the checklist in R as vector
  # "Alternative 1":
  # read the file from your computer
  # checklist <-  readLines(paste0(path.input,"Natural_Distribution_Checklist_TDWG.txt"))
  # Convert the vector into a list. 
  # checklist <- strsplit(checklist,";") # The first character of item is the name of genus, then there are the country codes
  
  # "Alternative 2":
  # As an example, we provide the code to retrieve the native ranges of the genera from POWO
  # using the library kewr
  
  # CAUTION: Please note that only a few families are listed in POWO, so it is often necessary to create the checklist manually 
  # by referring to national checklists see step A
  
  # Retrieve taxonomic information from IPNI for each taxon name in the list
  # The function search_ipni() is used with the filter set to "genera"
  taxa_ipni <- lapply(taxa.names, search_ipni, filters = "genera")
  
  # Extract the taxonomic IDs from the IPNI search results
  taxa_id <- lapply(taxa_ipni, tidy)

  # Select the taxonomic ID codes that are in POWO
  # CAUTION: The database of distributions in POWO may not have complete data. Missing data is common.
  # If a taxon is missing, the code will return NA. In such cases, the distributions
  # need to be checked manually
  taxa_id <- lapply(taxa_id, dplyr::filter, inPowo == TRUE)
  
  # check if there is any missing taxa
  length(taxa_id) == length(taxa.names)
  
  taxa_id <- lapply(taxa_id, dplyr::select, "id") # Extract the taxonomic IDs from the retrieved data
  
  # Retrieve additional information from Plants of the World Online (POWO) for each taxonomic ID
  # The function lookup_powo() is used with the parameter distribution set to TRUE
  record <- lapply(taxa_id, function(id) {
    print(id) # Print the taxonomic ID for reference
    
    # Try to retrieve information from POWO using the taxonomic ID
    # If an error occurs during the retrieval, assign NA (missing value) to the record
    tryCatch(
      lookup_powo(id, distribution = TRUE),
      error = function(e) NA
    )
  })
  
  # Tidy up the retrieved records
  tidied <- lapply(record, tidy)
  
  # Create a list of genus and distribution data for native species
  list <- lapply(tidied, function(t) {
    t %>%
      dplyr::select(genus, distribution) %>%
      tidyr::unnest(cols = distribution) %>% 
      tidyr::unnest(cols = natives) %>%
      dplyr::filter(establishment == "Native")
  })
  
  # Create a checklist by combining unique genus names and TDWG codes
  checklist <- lapply(list, function(l) {
    checkl <- c(unique(l$genus), l$tdwgCode)
  })
  
    
#### D. DOWNLOAD DATA ####
  # 1. Records from GBIF (https://www.gbif.org/)
	# a. Indicate the username, email and password created in step B-2 (see step descriptions in the main manuscript)
    user <- "XX" # your gbif.org username 
    pwd <- "XX" # your gbif.org password
    email <- "XX" # your email 
      
    # b. Search the taxon keys of each taxon
    # using the library rgbif
    taxon.keys <- sapply(taxa.names,function(x) name_backbone(x, rank="genus")$usageKey) # in this example, the sample units are "genus"
    # If your sample unit is not genus, then you need to replace rank="genus" 
    # by rank="species"or rank="family"

    # c. Prepare a download request and download data
	  # i. Create a path to the folder that will contain the files downloaded from GBIF 
      path.download.gbif <- paste0(path.output,"download_GBIF")
      
      # Create the folder named "download_GBIF" inside "input" folder
      dir.create(path.download.gbif)
    
      # Select "download_GBIF" as working directory to proceed with download
      setwd(path.download.gbif)
    
      # Specify the number of taxa contained in your list automatically
      n.spp <- length(taxa.names)
      
      # ii. Prepare and download the occurrences
      # Request download
      gbif.download <- occ_download(pred_in("taxonKey", taxon.keys), # Download only records with the taxon.keys recorded for the taxa
                                    pred("hasCoordinate", TRUE), # Download only records with coordinates
                                    pred_not(pred("establishmentMeans","INTRODUCED")), # Discard INTRODUCED establishment means of a taxon
                                    pred_not(pred("establishmentMeans","INVASIVE")), # Discard INVASIVE establishment means of a taxon
                                    pred_not(pred("establishmentMeans","MANAGED")), # Discard MANAGED establishment means of a taxon
                                    pred_not(pred("establishmentMeans","NATURALISED")), # Discard NATURALISED establishment means of a taxon
                                    user = user, pwd = pwd, email = email, format = "SIMPLE_CSV")
        
      
      # The download may take some time. You can use this function to check the status of the download 
      # Once the download is finished "status: Succeeded. Download is done" the function provides you with the "status: Succeeded. Download is done" message and a GBIF reference of the download
      occ_download_wait(gbif.download) 
      # Do not worry if a bit after retrieving the "status: running" message, you recovered the following error:
      # Error in curl::curl_fetch_memory(x$url$url, handle = x$url$handle) : 
      # HTTP/2 stream 7 was not closed cleanly before end of the underlying stream
      # the download is still running. You just have to wait (dowload can take time) and run the function in line 198 again. When the download is finished
      # you will recover the "status: succeeded" message
      
      # The object gbif.download.key allows you to keep the GBIF reference of the donwload
      gbif.download.key <- as.character(gbif.download)
      
      # Download as many zip files to the "download_GBIF" folder as the number of taxa specified in L98  
      occ_get <- occ_download_get(key=gbif.download.key,path = path.download.gbif, overwrite = TRUE) 
	    # once the download is finished check in "download_GBIF" folder that there is a zip file named "gbif.download.key". See Note 4 in the manuscript
      list.files()
      
      # Import to R all zip files downloaded from GBIF
      raw.GBIF.dataset <- occ_download_import(occ_get)
      
  # 2. Records from BIEN (https://bien.nceas.ucsb.edu/bien/)
  raw.BIEN.dataset <- list()
      
  # Download BIEN occurrences for each element of taxa.names. This step may take some time
  # using the library BIEN
  for (i in 1:n.spp ){
    # Name of taxa
    x <-taxa.names[i]
      
    # Occurrences of each genus in data frame (remember to deselect cultivated occurrences by setting "cultivated = FALSE")
    df<- BIEN_occurrence_genus(genus=x, cultivated = FALSE,all.taxonomy = TRUE,collection.info = TRUE, 
                               # If sample unit is "species" then use function "BIEN_occurrence_species" instead of "BIEN_occurrence_genus"
                               # and "species" as the first argument instead of "genus".See Note 5 in the manuscript
                               observation.type = TRUE,political.boundaries = TRUE, natives.only = TRUE)
      
    # Introduce the data frame downloaded from BIEN as a list element in R
    raw.BIEN.dataset[[i]]<-df 
    }
      
  # Reformat the data from list to data.frame type. 
  raw.BIEN.dataset <- do.call(rbind,raw.BIEN.dataset) 
  # If there are no records in BIEN this object will be appear as "0 Obs". 
  # Then, see Note 6 in the manuscript. Only the object "raw.GBIF.dataset" will be further used in the R scrip.
      
  # 3. Save workspace with downloaded data in output folder.
  save.image(paste0(path.output,"1_Workspace_Download.RData"))
  # load(paste0(path.output,"1_Workspace_Download.RData")) # This functions allows you to load the workspace with all R objects done so far. If needed, uncomment it and run it.
    
  # Select "output" as working directory
  setwd(path.output)
    
#### E. UNIFY AND SIMPLIFY DOWNLOADED DATA ####
  # load(paste0(path.output,"1_Workspace_Download.RData")) # This line is commented because it is only neccessary if there is a problem in later steps. 
  # If there is a problem later and you need to startover again with the original data downloaded and merged, remove "#" to load the previous workspace. 
 
  # 1. Simplify raw GBIF list
  simple.GBIF.dataset <- raw.GBIF.dataset
      
    # a. Create new columns for further merging simple.GBIF.dataset and simple.BIEN.dataset
      # i. Country names column from 2-letter ISO 3166-1 to country name
      # using the library countrycode
      simple.GBIF.dataset$countryName<- countrycode(simple.GBIF.dataset$countryCode, origin = "iso2c", destination = "country.name") 
	    # Do not worry if the following warning appears 'Some values were not matched unambiguously:  , ZZ'. 
	    # See Note 7 in the manuscript. The records with ZZ and empty values in ISO code will have NAs in the “countryName” column.
        
      # ii. Add column with data origin
      simple.GBIF.dataset$dataOrigin <- "GBIF"
    
    # b.Select useful columns to simplify the dataset. 
    # Columns with useful information
    useful.col.GBIF <- c("gbifID","dataOrigin","basisOfRecord","genus","species","scientificName","decimalLongitude","decimalLatitude","elevation",
                         "countryName","countryCode","locality","eventDate","institutionCode","collectionCode","catalogNumber")
        
    # Simplify dataframe by keeping only useful columns
    simple.GBIF.dataset <- dplyr::select(simple.GBIF.dataset,all_of(useful.col.GBIF)) 
    
    # c. Export data.frame as csv file
    write.csv(simple.GBIF.dataset, paste0(path.output,"1_simple_GBIF_dataset.csv"), row.names = FALSE)
      
      
  # 2. Simplify raw BIEN dataset
    # a. Remove "data_collected". This columns is duplicated because "collection.info" is TRUE in download.
    simple.BIEN.dataset <- raw.BIEN.dataset[,-28] 
	  # If your sample unit is species, then replace raw.BIEN.dataset[,-28] by raw.BIEN.dataset[,-27]. See Note 8 in the manuscript.
      
    # b. Remove records with GBIF data source to avoid replicates
    simple.BIEN.dataset <- dplyr::filter(simple.BIEN.dataset, datasource !="GBIF")
      
    # c. Create new columns for further merging between simple.GBIF.dataset and simple.BIEN.dataset
      # i. Country codes according 2-letter ISO 3166-1
      simple.BIEN.dataset$country_ISOcode<- countrycode(simple.BIEN.dataset$country, origin = "country.name", destination = "iso2c")
        
      # ii. Elevation
      # Create an empty column, this information is not available in the BIEN database, but it is include in GBIF and we want to keep it.
      simple.BIEN.dataset$elevation <- NA
        
      # iii. ID Origin
      simple.BIEN.dataset$ID_Originin <- NA
        
      # iv. dataOrigin
      simple.BIEN.dataset$dataOrigin <- "BIEN"
        
      # v. If your sample unit is species uncomment the following line (209) and run it. 
	    # This function creates a column to include the name of the Genus in a separate column. # See Note 9 in the manuscript
      # simple.BIEN.dataset <- simple.BIEN.dataset %>% separate(scrubbed_species_binomial, c("scrubbed_genus","temp_spp"),sep= " ", remove=FALSE)
        
        
    # d. Select useful columns from all available variables for our dataset
    useful.col.BIEN <- c("ID_Originin","dataOrigin","observation_type","scrubbed_genus","scrubbed_species_binomial","verbatim_scientific_name","longitude","latitude","elevation",
                         "country","country_ISOcode","locality","date_collected","datasource","collection_code","catalog_number")
      
    simple.BIEN.dataset <- dplyr::select(simple.BIEN.dataset,all_of(useful.col.BIEN)) 
      
    # e. Export data.frame as csv file
    write.csv(simple.BIEN.dataset, paste0(path.output,"1_simple_BIEN_dataset.csv"), row.names = FALSE)
        
  # 3. Rename the columns that will be used for merging between simple.GBIF.dataset and simple.BIEN.dataset. 
  # This step is necessary even if there was no record downloaded from BIEN
  # Rename the columns of GBIF and BIEN
  colnames.matched <- c("ID_Originin","Data_Origin","Basis_of_Record","Genus","Spp","Scientific_name","Longitude","Latitude","Elevation",
                        "Country_Name","Country_ISOcode","Locality","Date","Institution_code","Collection_code","Catalog_number")
  colnames(simple.GBIF.dataset) <- colnames.matched # If there are no records in "raw.BIEN.dataset", then run only this line, skip line 227 and go to line 231. See Note 10 in the manuscript
  colnames(simple.BIEN.dataset) <- colnames.matched
      
  # 4. Merge matched.GBIF.dataset and matched.BIEN.dataset. 
  # This step is necessary even if there was no record downloaded from BIEN
  merged.dataset <- rbind(simple.GBIF.dataset,simple.BIEN.dataset) # If there are no occurrences in BIEN, replace with: merged.dataset <- simple.GBIF.dataset. See Note 10 in the manuscript
      
  # 5. Save merged.dataset as csv file. This step is necessary even if there was no record downloaded from BIEN
  write.csv(merged.dataset, paste0(path.output,"2_merged_dataset.csv"), row.names = FALSE)
      
#### F. Add occurrences from other sources #### 
# (This step is an example, it has not been performed in this protocol) but it is provided
# because it might be possible that the records gathered in BIEN and GBIF do not accurately 
# reflect the distribution range of your case study
# and you may need to complete the merge.dataset with records comming from other sources. 
# If you want to run this part of the code uncomment lines 244-253
  
  # 1. For online databases we use spocc R package that allows downloading records from multiple online databases aside from (including GBIF) 
  # using the library spocc
  # raw.spocc.dataset <- list()
  # for (i in 1:length(taxa.names) ){
  #   # Name of taxa
  #   y <-taxa.names[i]
  #   # Occurrences of each genus in data frame (remember to deselect cultivated occurrences by setting "cultivated = FALSE")
  #   df<- occ(query = y, from = c("gbif", "bison", "inat","ebird", "ecoengine", "vertnet"), limit = 2000,has_coords = TRUE)
  #   # Introduce the data frame as list element
  #   raw.spocc.dataset[[i]]<-df 
  # }
  # 2. Add records manually obtained from herbaria and literature using Microsoft Excel and save the file as "merged_dataset_Version2.csv"). 
  # The separator for latitude and longitude columns must be a period. See Note 11 in the manuscript
      
#### G. CHECK THE DATASET ####
  # 1. If new occurrences have been manually added in excel (F-2) then you need to uncomment line 259 to import the csv created in step F-2. 
  # merged.dataset.2 <- read.csv(paste0(path.input,"merged_dataset_Version2.csv"),sep=";",dec=".")
  # merged.dataset.2 <- read.csv(paste0(path.input,"merged_dataset_Version2.csv"),sep=",",dec=".") # if you have different separators (not comma nor semicolons), then indicate the type of separator in sep=""
  
  # 2. Visualize the dataset (merged.dataset or merged.dataset.2 if you run line 259) 
  # This is done to check that the dataset meets the criteria needed to run the following part of the protocol
  
  # Visualize
  View(merged.dataset) # If new occurrences have been manually added in excel (F-2), then replace merged.dataset by merged.dataset.2. See Note 12 in the manuscript
  
  # Number of rows and columns
  dim(merged.dataset) # If new occurrences have been manually added in excel (F-2), then replace merged.dataset by merged.dataset.2. See Note 12 in the manuscript
  # If different than 16 the dataset is not OK and you need to solve the issue to make sure that you only have the 16 columns in E-2-d (see manuscript)
      
  # Structure
  str(merged.dataset) # If new occurrences have been manually added in excel (F-2), then replace merged.dataset by merged.dataset.2. See Note 12 in the manuscript   
  # The object must be a data.frame, columns "Latitude", "Longitude" and "Elevation" must be numeric, and the remaining columns must be characters.
  # If any of these requirements is not meet, the dataset is not OK and you need to solve whatever the issue detected before you go to step H
  
#### H. DATA CLEANING ####
  # 1. Remove the records that lack coordinates (common in BIEN)
  filtered.data <- merged.dataset[!(is.na(merged.dataset$Longitude) & is.na(merged.dataset$Latitude)),] # If new occurrences have been manually added in excel (F-2), then replace merged.dataset by merged.dataset.2. See Note 12 in the manuscript
    
  # 2. Remove the records with invalid coordinates (Lon or Lat = 0) (common in GBIF)
  filtered.data <- filtered.data[!(filtered.data$Longitude==0 & filtered.data$Latitude==0),]
    
  # 3. Remove the records with coordinates < 2 decimal places. 
  # For our geographical scale (world) and the ultimate gold (climatic characterization of sample units), two decimals are enough precision.
  # If your geographical scale is different or you need more precission, you can modified it as needed.
  filtered.data <- filtered.data[grep("\\.[0-9]", filtered.data$Longitude), ] 
  
  # If you want to keep only records with three decimals, then replace "\\.[0-9][0-9]" by "\\.[0-9][0-9][0-9]";
  # Alternatively, if you one decimal is enough to your purpose, then replace replace "\\.[0-9][0-9]" by "\\.[0-9]"
  
  filtered.data <- filtered.data[grep("\\.[0-9]", filtered.data$Latitude), ] 
  # If you want to keep only records with three decimals, then replace "//.[0-9][0-9]" by "\\.[0-9][0-9][0-9]";
  # Alternatively, if you one decimal is enough to your purpose, then replace replace "\\.[0-9][0-9]" by "\\.[0-9]"
    
  # 4. Remove replicated records 
  filtered.data <- filtered.data[!duplicated(filtered.data[,c("Genus","Spp","Longitude","Latitude","Date","Basis_of_Record","Elevation","Catalog_number","Country_ISOcode")]),] 
    
  # 5. Remove records outside natural distribution of the genera 
    # a. Add level-3 TDWG code to merged.dataset for future comparation with the checklist created in Step-A-5
      # ii. Load the world shapefile and select the projection       
     
      # using the library sf
      # using the library raster
  
      shape.level3TDWG <- st_read("/vsicurl/https://github.com/tdwg/wgsrpd/raw/master/level3/level3.shp")
      shape.level3TDWG <- as(shape.level3TDWG[,"LEVEL3_COD"], "Spatial")
    
      # iii. Extract level-3 TDWG code for the occurrences
      # Convert your dataset (filtered.data) to spatial object and select the projection
      col_coord <- c("Longitude","Latitude") 
      filtered.data.spatial <- filtered.data
      coordinates(filtered.data.spatial)<-col_coord 
      projection(filtered.data.spatial) <- crs("+proj=longlat +datum=WGS84")
     
      # Extract all values of shape.level3TDWG for each record in your dataset (filtered.data)
      values.shape <- raster::extract(shape.level3TDWG,filtered.data.spatial) # If a warning may appear, it does not affect the result. See Note 13 in the manuscript
        
      # Replace the values extracted in line 319 to the level-3 TDWG codes and convert the spatial object into a dataframe
      Country_TDWGcode <- as.character(values.shape$LEVEL3_COD)
      filtered.data.TDWG<-cbind.data.frame(filtered.data.spatial,Country_TDWGcode) 
        
      # iv. Remove records outside level-3 TDWG shapefile (occurrence that fall in the sea)
      filtered.data.TDWG <- filtered.data.TDWG[complete.cases(filtered.data.TDWG$Country_TDWGcode), ]
        
      # Remove optional column created by default when running "extract". See Note 14 in the manuscript
      filtered.data.TDWG$optional <- NULL
        
    # b. Divide filtered.data.TDWG in as many dataframes as genera (or the sample unit you are using) and store them in a list
    # Vector with taxa names downloaded
    taxa.downloaded <- sort(unique(filtered.data.TDWG$Genus)) # Replace f=filtered.data.TDWG$Genus with filtered.data.TDWG$Spp if sample units is species
   
    # Number of taxa downloaded
    n.taxon.downloaded <- length(taxa.downloaded)
    
    # Divide dataset
    split.taxon <- split(filtered.data.TDWG, f=filtered.data.TDWG$Genus) # If your sample units is species, then replace f=filtered.data.TDWG$Genus by f=filtered.data.TDWG$Spp. 
    
    # Create a new list to add filtered taxa
    filtered.list <- list()
   
    # Add taxa names to checklist created in A-5
    names(checklist) <- taxa.names 
    
    # Extract the distribution of the sample units downloaded.
    # This is an essential step because your original checklist may include more sample units (genera in this example) than the sample units that are included in your dataset (split.taxon)
    # This happens when you did not recover records for all the sample units you have in your case study
    checklist.filtered <- checklist[taxa.downloaded] 
    
    # Check that the names of the taxa names
    View(checklist.filtered) # The names in "Name" column of checklist.filtered must be the same as those in the "value"
    
    
    # c. Filter your dataset (filtered.data.TDWG) to retain only natural records. 
    for (x in 1:n.taxon.downloaded) { # Each sample unit
      taxa<-taxa.downloaded[x]
      natdist <- checklist.filtered[[x]][2:length(checklist.filtered[[x]])] # The first element's vector is the name of the sample unit
      filtered.taxa <- split.taxon[[x]] # Records of each sample unit
      taxa.nat <- list() # This list will include only the natural records
      
      for (i in 1:length(natdist)){ # Each TDWG countrycode
        taxa.nat[[i]] <- subset(filtered.taxa,filtered.taxa$Country_TDWGcode==natdist[i])}
      
      filtered.taxa <- do.call(rbind,taxa.nat)
      filtered.list[[x]]<-filtered.taxa
    }
    
    filtered.dataset.WCSP <- do.call(rbind,filtered.list) # Convert the filtered list into a data.frame
    
    # d. Export filtered.dataset.WCSP as csv file
    # csv
    write.csv(filtered.dataset.WCSP, paste0(path.output,"3_Cleaning_dataset_WCSP.csv"), row.names = FALSE)
    
    # Workspace with cleaned data
    save.image(paste0(path.output,"2_Workspace_Cleaning.RData"))

#### I. DISTRIBUTION MAPS ####
  # 1. Create global distribution map with all sample units in PDF format
  # Create and open PDF
  pdf(file=paste0(path.output,"4_global_distribution_map.pdf"), width = 8,height = 8)
  
  # Create base map 
  plot(shape.level3TDWG, col="grey40", bg="white", lwd=0.25, border=1,legend = FALSE)
  
  # Add ocurrences to base map
  points(filtered.dataset.WCSP$Longitude, filtered.dataset.WCSP$Latitude, pch = 21, cex= 0.8, col= "black", bg="red") 
  
  # Close PDF
  dev.off()
   
  # 2. Crear a distribution map for each sample unit in PDF format
  # This function creates a list from the original data.frame with as many elements as the number of sample units you have
  split.maps <- split(filtered.dataset.WCSP,filtered.dataset.WCSP$Genus) # If your sample unit is species, then replaced $Genus by $Spp  
  
  # Create and open PDF
  pdf(file=paste0(path.output,"4_distribution_maps.pdf"), width = 8,height = 8)

  # Built distribution maps for each sample unit
    for (t in 1:n.taxon.downloaded){ # Each sample unit
      title.name <- taxa.downloaded[t] # Title of the map, it corresponds with the name of the sample unit
      taxon <- split.maps[[t]] # Records of the sample unit
      coordinates(taxon)<-col_coord # Convert to spatial object
      projection(taxon) <- crs("+proj=longlat +datum=WGS84") # Select WGS84 projection
      ext <- extent(taxon)+10 # Expand the map to the extent of each sample unit
      plot(crop(shape.level3TDWG,extent(ext)), col="grey40", bg="white", lwd=0.25, border=1,legend = FALSE, main=title.name) # Base map
      points(taxon$Longitude, taxon$Latitude, pch = 21, cex= 0.8, col= "black", bg="red") # Add occurrences
    }
	
  # Close PDF
  dev.off()

#### J. DATA THINNING #### 
  # If the case study is not spatially biased (that is, there is no sample unit or region in the world that is significantly
  # oversampled, see the exported maps from step I),then go directly to step K. 
  
  # 1. Remove occurrences ramdomly with a buffer of 10 kilometres for all sample units
  # This function creates a list from the original data.frame with as many elements as the number of sample units you have
  split.thin <- split(filtered.dataset.WCSP,filtered.dataset.WCSP$Genus) # If your sample unit is species, then replaced $Genus by $Spp  
    
  # Path to export thin files
  path.thin <- paste0(path.output,"thin/")
  
  # Create the folder
  dir.create(path.thin)
    
    # a. Thinning process for all sample units. Note 
  # CAUTION: Please be aware that this step consumes a significant amount of memory and may fail if the available memory is low.
  # Additionally, it can take a substantial amount of time to complete. Ensure that the number of CSV files obtained matches the number of sample units (refer to Note 15 in the manuscript).
  
  # In this thinning we apply a buffer of 10 km (thin.par=10), but you may need to modify this buffer and adapt it to your geographical scale
  # using the library spThin
    lapply(1:n.taxon.downloaded,function(x){
      df.taxa<-split.thin[[x]]
      thin(df.taxa,lat.col = "Latitude",long.col = "Longitude",thin.par = 10,reps = 1,#replace 10 in thin.par=10 by the kms you want to apply as buffer
           write.files = TRUE,out.dir = path.thin,spec.col = "Genus",out.base = paste0("5_thinned_data_",taxa.names[x]), max.files = 1)
      })
    
    # b. Thinning process for one sample unit. There are two cases in which you may need to perform this step 
    # Case 1. The sampling bias in your dataset does not affect the majority of sample units (thus, you did not run step J-1-a) 
    # but you need to do the thinning in the sample unit for which you have detected a spatial bias
    # Case 2. You did perform step J-1-a, but there is still one of the sample units that suffer from spatial bias
    # In any of the two cases, uncomment lines 444-449 and run them to apply the thinning to the given sample unit.
    # If you are in case 1, you need to adapt the buffer to your geographical scale (see manuscript step-J).
    # If you are in case 2, you need to increase the buffer in relation to the one you used in step J-1-a (see manuscript step-J).
    # 
	# using the library spThin
    # x=1 # Replace "1" by the number of the position of the taxa in n.taxon.downloaded
    # df.taxa<-split.thin[[x]]
    # thin(df.taxa,lat.col = "Latitude",long.col = "Longitude",thin.par = 10,reps = 1,#replace 10 in thin.par=10 by the kms you want to apply as buffer
    #        write.files = TRUE,out.dir = path.thin,spec.col = "Genus",out.base = paste0("4_thinned_data_",taxa.names[x]), max.files = 1) 
    # If you run this step because you are in case 2 (that is, after running step G-1-a there is still one sample unit with spatial bias), 
	# then it is highly important that you go to "thin" folder look for the csv created in step K-1-a for the sample unit for which you repeating the thinning and delete it before go to step 2. 
	# Note that the number of csv files in "thin" folder must be equal to the number of sample units in taxa.names. See Note 16 in manuscript.
    
  # 2. Import data after thinning process
    # a. Import csv files of all taxa thinned (that is, if you run step-J-1-a only or followed by running step-J-1-b (Case 2))
    files.thin <- list.files(path.thin,pattern="*.csv") 
    setwd(path.thin)
    thin.occ <- sapply(files.thin, read.csv, simplify=FALSE) %>% 
      bind_rows(.id = "id")
	  # If a warning appears, it does not affect the results, check the object "thin.occ" to confirm that it contains as many all the information 
	  # (three columns: "Genus" or "Spp","Longitude" and "Latitude"). See Note 17 in the manuscript
    
	  # b. Import csv file for one taxon thinned (that is, if you run step-J-1-b only (Case 1))
    # thin.taxon.file <- list.files(path.thin,pattern="*.csv")
    # setwd(path.thin)
    # thin.taxon <- read.csv(thin.taxon.file, header = TRUE, sep=",")
    # 
    
  # 3. Join columns before thinning. After the thinning the files only have three columns
  # You need to add all the remaining 13 columns of filtered.dataset.WCSP
    
    # a. Join columns of all taxa thinned. Perform this step if you run step J-1-a or J-1-a and J-1-b
    # Vector with thin.occ colnames
    nms1 <- colnames(thin.occ)
      
    # Convert columns into character format
    thin.occ[nms1] <- lapply(thin.occ[nms1], as.character) 
      
    # Vector with filtered.data.WCSP colnames
    nms2 <- colnames(filtered.dataset.WCSP)
    filtered.dataset.WCSP[nms2] <- lapply(filtered.dataset.WCSP[nms2], as.character) 
      
    # Join datasets
    # using the library data.table
    dt1 <- data.table(thin.occ)
    dt2 <- data.table(filtered.dataset.WCSP)
    setkey(dt1,Genus,Longitude,Latitude) # If your sample unit is species, then replaced $Genus by $Spp 
    setkey(dt2,Genus,Longitude,Latitude) # If your sample unit is species, then replaced $Genus by $Spp 
    joined.dataset<- dt2[dt1] 
	
    # Check that the number of rows of "thin.occ" and "joined.dataset" is the same
    # If the number of rows do not match, it is possible that joining the two tables may generate duplicates. Uncomment next line 491. See Note 18 in manuscript.
    # joined.dataset<-joined.dataset[!duplicated(joined.dataset[,c("Genus","Longitude","Latitude")]),] # If your sample unit is species, then replaced $Genus by $Spp 
    # this remove possible duplicates generated when merging the "joined.dataset" and "thin.occ"
    # once run line 493, that the number of rows of "thin.occ" and "joined.dataset" are the same
    # Check that "joined.dataset" contains the columns of "filtered.dataset.WCSP"
    
    # b. Join columns of one taxon thinned and the rest of taxa
    # # Vector with thin.taxon colnames
    # nms1 <- colnames(thin.taxon)
    # 
    # # Convert columns into character format
    # thin.taxon[nms1] <- lapply(thin.taxon[nms1], as.character) 
    # 
    # # Vector with filtered.data.WCSP colnames
    # nms2 <- colnames(filtered.dataset.WCSP)
    # filtered.dataset.WCSP[nms2] <- lapply(filtered.dataset.WCSP[nms2], as.character) 
    # 
    # # Join datasets
    # using the library data.table
    # dt1 <- data.table(thin.taxon)
    # dt2 <- data.table(filtered.dataset.WCSP)
    # setkey(dt1,Genus,Longitude,Latitude) # If your sample unit is species, then replaced $Genus by $Spp 
    # setkey(dt2,Genus,Longitude,Latitude) # If your sample unit is species, then replaced $Genus by $Spp 
    # joined.dataset.1taxon<- dt2[dt1] 
	  #
    # Check that the number of rows of "thin.taxon" and "joined.dataset.1taxon" is the same 
    # If the number of rows do not match, uncomment next line 519. See Note 18 in the manuscript
    # joined.dataset.1taxon<-joined.dataset.1taxon[!duplicated(joined.dataset.1taxon[,c("Genus","Longitude","Latitude")]),] # If your sample unit is species, then replaced $Genus by $Spp 
    # this remove possible duplicates generated when merging the "joined.dataset" and "thin.occ"
    # once run line 519, that the number of rows of "thin.occ" and "joined.dataset" are the same
    # Check that "joined.dataset" contains all columns in "filtered.dataset.WCSP"
  
    # # Add joined.dataset.1taxon to the rest of taxa
    # joined.dataset <- rbind(joined.dataset.1taxon, filtered.dataset.WCSP)
      
  # 4. Export dataset as csv file and workspace
  # csv file
  write.csv(joined.dataset, paste0(path.output,"6_Joined_dataset.csv"), row.names = FALSE)
    
  # Save worckspace
  save.image(paste0(path.output,"3_Workspace_Thinning.RData"))
    
#### K.  LOAD BIOCLIMATIC VARIABLES FROM WORLDCLIM ####
  # a. Alternative 1: Download Bioclimatic variables layers directly from Worldclim webpage.
    # a. Download variables in: https://www.worldclim.org/data/worldclim21.html. Unzip the file in "input" folder and rename "Bioclimatic _variables_WC2'
    # b. Rename bioclimatic variables 1 to 9 (see manuscript K-a)
    # c. Import climatic variables from Alternative 1 of 30 seconds resolution
    # path.clim <- "C:/Users/User/Documents/input/Bioclimatic _variables_WC2/"
	  # If using iOS the path will be written as this (uncomment line 541):
    # path.clim <- "/Users/User/Documents/input/Bioclimatic _variables_WC2/"
    # 	  files.clim <- list.files(path.clim,pattern="*.tif", full.names = TRUE) # Check that there are 19 files
    #     using the library raster
    #     bioclim <- stack(files.clim) # this makes a list of all the .tif files, then calls raster() on them, then stacks them together
    #     names(bioclim) <- gsub("wc2.1_30s_", "", names(bioclim)) # Remove "wc2.1_30s_" from bioclimatic names
    #     names(bioclim) <- gsub("_", "", names(bioclim)) # Remove "_" from bioclimatic names
    #     names(bioclim) <- gsub(".tif", "", names(bioclim))# Remove ".tif" from bioclimatic names
    #     projection(bioclim) <- crs("+proj=longlat +datum=WGS84") # Select WGS84 projection
    
  # b. Alternative 2: Download variables in R for 10, 5 or 2.5 minutes.
  # If using Alternative 1. Skip lines 637 and 638 of this code and got to STEP L.
  bioclim <- getData("worldclim",var="bio",res=10)
  projection(bioclim) <- crs("+proj=longlat +datum=WGS84")
    
#### L.  SPATIAL CORRECTION FOR TERRESTRIAL ORGANISMS####
  #  This step is to identify if there are occurrences in the "joined.dataset" that fall in the sea
  #  Set a distance buffer to the coastal limit of a given template. See Note 19 in the manuscript
  #  Remove all occurrences outside the buffer
  #  Recalculate coordinates for occurrences within the buffer so that the new coordinates fall in the nearest climatic cell of the template
   
  # 1. Visualize the distribution of the joined.dataset occurrences
  # Uncomment the following line ONLY if you have not performed step J
  # joined.dataset <- filtered.dataset.WCSP
      
  # Color palette for variable 1 of the bioclimatic variables
  colfunc <- colorRampPalette(c("black", "grey"))
    
  # Plot of bioclimatic variable 1
  # We assume that all bioclimatic variables have the same geographical boundaries, so we only use variable 1 to make the comparison.
  pdf(file=paste0(path.output,"7_total_points.pdf"), width = 8,height = 8)
  plot(bioclim[[1]],col=colfunc(10),legend = FALSE)
    
  # Plot all the occurrences of the joined.dataset
  points(joined.dataset$Longitude, joined.dataset$Latitude, pch = 21, cex= 0.8, col= "black", bg="blue") 
  dev.off()
    
  # 2. Convert joined.dataset into spatial object
  # New columns for corrected coordinates
  joined.dataset$LongitudeCorrected <- as.numeric(joined.dataset$Longitude) #to avoid formatting issues
  joined.dataset$LatitudeCorrected <- as.numeric(joined.dataset$Latitude)#to avoid formatting issues
  joined.dataset$dist <- rep(0, dim(joined.dataset)[1]) # then enter in this column the distance between the nearest climatic cell and the occurrences outside the terrestrial boundary of the bioclimatic layer used as template.
  joined.dataset <- as.data.frame(joined.dataset)
    
  # Convert  to spatial object
  joined.dataset$Longitude <- as.numeric(joined.dataset$Longitude) #to avoid formatting issues
  joined.dataset$Latitude <- as.numeric(joined.dataset$Latitude) #to avoid formatting issues
  joined.dataset.spatial <- as.data.frame(joined.dataset)
  coordinates(joined.dataset.spatial)<-col_coord # convert to spatial object
  projection(joined.dataset.spatial) <- crs("+proj=longlat +datum=WGS84") # Select WGS84 projection
    
  # 3. Check if all joined.dataset.spatial points are within the boundaries of the bioclimatic variable layers
  # Extract climatic data for all occurrences of joined.dataset.spatial, the NAs correspond to the occurrences outside the boundaries of the layer.
  vals <- raster::extract(bioclim[[1]], joined.dataset.spatial,na.rm=TRUE) 
    
  # Logical vector that if TRUE there is a NA
  outside_mask <- is.na(vals)
    
  # Extract coordinates of occurrences with NA values (occurrences outside limits)
  outside_pts <- joined.dataset.spatial[outside_mask,,drop = FALSE] 
    
  # When occurrences with NA values = 0, success= TRUE
  success <- sum(is.na(vals)) == 0 
    
  # Distance from which we correct the coordinates (1000 m)
  dist <- 1000 
    
  # Plot occurrences outside layer limits (NA values) and save
  # If you want to extract these maps, uncomment lines 610 and 613
  # If occurrences outside limits are not obtained, go to step M-1. See Note 20 in the manuscript
  pdf(file=paste0(path.output,"7_points_outside.pdf"), width = 8,height = 8) 
  plot(bioclim[[1]],col=colfunc(10),legend = FALSE,  main="Ocurrences outside limits")
  points(outside_pts, pch = 21, col= "black", bg="red", cex = 1)
  dev.off()
    
  # 4. Coordinates correction
  # This loop calculates the climatic cell closest to the point and corrects the coordinates, using "dist" as reference. 

  while (!success) { # When success=TRUE, the loop stops. 
      
    print(paste(Sys.time(),"Starting spatial correction for distance:",dist))
    # Coordinates of all occurrences of joined.dataset
    newcoord <- coordinates(joined.dataset[,c("LongitudeCorrected","LatitudeCorrected")])
      
    # Extract climatic data of variable 1 for all occurrences 
    newvals <- raster::extract(bioclim[[1]], newcoord)
      
    # Extract records with NA values (occurrences outside limits)
    newoutside_mask <- is.na(newvals)
    newoutside_pts <- newcoord[newoutside_mask,,drop = FALSE]
      
    # Calculation of the nearest neighbor with maximum distance="dist"
    land <- seegSDM::nearestLand(newoutside_pts, bioclim[[1]], dist) 
      
    # Plot of corrected and uncorrected occurrences by distance="dist"
    plot(crop(bioclim[[1]], extent(newcoord)),col=colfunc(10), main=dist,legend = FALSE) # Plot with variable 1 of worldclim as template
    points(newoutside_pts, col = "red",pch = 20, cex = 1) # Outside points (red)
    points(land, col = "green",pch = 20, cex = 1) # Corrected points (green)
      
    # Save the coordinates of the corrected points in the new columns
    joined.dataset[newoutside_mask,"LongitudeCorrected"] <- land[,1]
    joined.dataset[newoutside_mask,"LatitudeCorrected"] <- land[,2]
      
    # Complete the NA values of the new coordinates columns with previous coordinates
    joined.dataset$LongitudeCorrected[is.na(joined.dataset$LongitudeCorrected)] <- joined.dataset$Longitude[is.na(joined.dataset$LongitudeCorrected)]
    joined.dataset$LatitudeCorrected[is.na(joined.dataset$LatitudeCorrected)] <- joined.dataset$Latitude[is.na(joined.dataset$LatitudeCorrected)]
      
    # Indicate the distance of the corrected points from the climatic cells   
    joined.dataset[outside_mask,"dist"] <- dist
      
    # Check if there are more points with NA values 
    success <- sum(is.na(land[, 1])) == 0
      
    # Distance limit to stop correcting coordinates (5000 m)
    success <- ifelse(dist > 5000, TRUE, success) # If you need to modify the limit of 5km, replace "5000" by the new distance. See Note 19 in tha manuscript
    print(paste(Sys.time(),"Spatial correction for distance",dist,"finished"))
    dist <- dist + 1000 # Next distance to correct the points
    # The plots will be updated with the points that are corrected and those that are not.
    }
  
#### M. EXTRACT CLIMATIC DATA ####
  # 1. Convert joined.dataset to spatial object. 
  col.coord.corrected <- c("LongitudeCorrected","LatitudeCorrected") 
  # Uncomment line 664  if Note 20 applies to your case 
  # col.coord.corrected <- c("Longitude","Latitude") 
  coordinates(joined.dataset)<-col.coord.corrected
  projection(joined.dataset) <- crs("+proj=longlat +datum=WGS84") # Select WGS84 projection
    
  # 2. Extract climatic values
  bioclim.values <- raster::extract(bioclim,joined.dataset, method="bilinear") 
  # If your dataset is large (high number of occurrences and/or high number of sample units, this step may take time
  # then consider replace method="bilinear" buy method="simple". See Note 21 in the manuscript
  
  # 3. Join the climatic data with the occurrences data
  climatic.data<-cbind.data.frame(joined.dataset,bioclim.values) 
    
  # Remove optional column
  climatic.data$optional <- NULL
    
#### N. VISUALIZE AND EXPORT THE FINAL DATASET ####
  # Visualize
  View(climatic.data)
    
  # Export as csv file
  write.csv(climatic.data, paste0(path.output,"8_Final_dataset.csv"), row.names = FALSE)
    
  # Save final workspace
   save.image(paste0(path.output,"4_Workspace_Final_Data.RData"))
   
