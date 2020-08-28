# IMECsentinel
MSc IMEC project on Sentinel2 data

## Project Background
This project was completed in 2020, as the final project for MSc International Marine Environmental Consultancy at Newcastle University.

The purpose of the project was to identify whether satellite derived bathymetry would be able to support the UK legal requirements to monitor habitat 1110, offshore seabanks covered by water at all times. It was focussed on a specific area of the North Sea.

## Set-up
The project is run solely using the R language, run through RStudio. All required packages are detailed within the code. RStudio version 1.3 using R version 3.6

The project was completed with the support of JNCC, who provided access to all required satellite and in-situ data. 
The satellite data used is from ESA Sentinel-2 mission, satellite tiles 31UCU, 31UCV, 31UDU, 31UDV, which have been pre-processed using the ACOLITE and ARCSI methods. 
In-situ data was kindly provided by Cefas for the purposes of this data. 

The expected directory structure is:

```
Parent directory  > SentinelData > Acolite > 'Year' > 'Date' > 'Tile' > 'filename'
                  > SentinelData > Arcsi > 'Year' > 'Date' > 'Tile' > 'filename'
                  > CefasData > Cefas_acoustic_2016 > 'filename'
                  > CefasData > ShapeFile > 'sand_ID'
                  > CefasData > Statistics > 
                  > DataSets
                  > DataSets > results.csv
                  > Graphs
```
'Year', 'Date', 'Tile' are all derived from the filenames, to allow for organisation. 'filename' is not important for the satellite data, as the code will read based on the year, date, tile data, however will need to be updated in the cefas data if they do not match.

For the sand_ID, this refers to a simple ability to recognise each unique sandbank, identified arbitrarily in the cefas_mapping at the beginning of analysisfunctions.R .

The results.csv needs to be generated prior to be populated, the order of results is on line 279 in analysisfunctions.R .

The code is designed to be self sufficient once the package code is run, at the beginning of each script. 
The working directories need to be updated to be suitable for the work environment. 
Once the working directories, folder paths and relevant files are created and accurately located, the code can be run.


## Usage
Basic order of running required:
```
1. Analysisfunctions.R - this is the main code base which creates the csv files required for the other scripts
2. 31UCV.R / Turbidity.R - these both work of the products of analysisfunctions
```
For analysisfunctions, running the final function, main(), will run the complete file for all variations. If a smaller set of data is wanted, add a hash tag to the beggining of each unrequired line  prior to running the file will stop it from being run. 

For 31UCV and turbidity, these will require the relevant code to be highlighted and then either 'run' or use 'source'. If numerical statistics are required from either of these scripts, this will have to be retrieved at intermediate stages, so only highlight code to where you want the results, or run each line individually.

## Future Ideas
This project has generated a large volume of potentials for development. These may be completed within this project, or branched by JNCC or other interested parties. All known links will be added.

## Acknowledgements
Many thanks to Helen Lillis and Dr Paula Lightfoot for proposing the project and Dr Lightfoot for
providing access to the satellite data required and her ongoing support throughout the project. Also,
to Cefas for supplying the in-situ data required and to the NERC Earth Observation Data
Acquisition and Analysis Service (NEODAAS) for supplying data for this study for the processing of
Sentinel-2 data. Academic support was gratefully received from Dr Clare Fitzsimmons and Dr
Rachel Gaulton, Newcastle University. The waterquality package used originates from Johansen, R. A., Nowosad, J., Reif, M., and Emery, E. (Nov. 2018). waterquality: Satellite Derived Water Quality
Detection Algorithms. R package. Version v0.2.2

I would also like to thank many family members, for providing
emotional and technical support throughout the slightly weird experience of a Masters project in
lock down. Without them, this project would have been a lot harder
