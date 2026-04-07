# Introduction 

The following python code correlates geologic data with areas of known nickel mineralization. This correlation is used to create predictor maps designed to represent the critical components of magmatic nickel mineralizaing systems. These predictor maps are combined to produce a final prospectivity map for Western Australia. This code also includes functions to validate the results of the prospectivity map. 

This approach may be useful to mineral explorers in identifying areas of potential mineral prospectivity, either through the supplementary data linked in this code, or integrating their own datasets.

# Getting Started

### Installation process

This project was carried out using ArcGIS Pro 2.4.0.1 with the spatial analyst extension. All python code used in this project was written and run using ArcGIS Pro's built in python notebooks. This code was collated and published on this GitHub page as four .py files:

1. **Data_Preperation.py** takes data from [GSWA's Data and Software Centre](https://dasc.dmirs.wa.gov.au/) and applies various filters and transformations to extract mappable proxies for a magmatic Ni-Cu-PGE sulfide mineralizing system.
2. **Data_Correlation_Predictor_Maps.py** calculates the spatial correlation between mappable proxies and areas of mineralization, and produces predictor maps based on that correlation.
3. **Map_Processing.py** combines the predictor maps into a mineral prospectivity map, and test the validity of the map aginst known mineral deposits
4. **tools.py** contains functions for the above files

###	Software dependencies
- ArcGIS Pro 2.4.0.1 with the *Spatial analyst* extension
- ArcPy 3.1
- Numpy 1.20.1
- Pandas 1.4.4

# Maintenance and Support

This codebase was developed for research purposes and is not actively maintained. If you encounter issues or have suggestions for improvement, feel free to submit them to [GitHub Issues](#link) or to [nathan.bowman@dmpe.wa.gov.au](#nathan.bowman@dmpe.wa.gov.au ). Feedback is appreciated and may inform future work, but any changes will be addressed in a new project, not this one.

# License

[MIT](https://choosealicense.com/licenses/mit/)