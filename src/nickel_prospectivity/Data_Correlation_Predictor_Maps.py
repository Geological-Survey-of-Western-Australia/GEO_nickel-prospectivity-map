# Import libraries
import arcpy
import arcpy.sa
import numpy as np
import pandas as pd

from tools import SpatialCorrelation, AssociationCorrelation, BuffMap, AsscRasters

arcpy.CheckOutExtension("Spatial")

# Set up folder pathways
root_dir = "C:/project_folder"  # Your project folder here
prjGDB = f"{root_dir}/Project.gdb"  # The .gdb for your project

# Set file pathways
grav1km = f"{prjGDB}/Grav_1km"  # Download from DASC Geophysical Surveys Geophysical compilations: https://dasc.dmirs.wa.gov.au/
geologyMerged = f"{prjGDB}/Geology_MERGED"  # Download from DASC Minerals Systems Atlas Primary Data Layers: https://dasc.dmirs.wa.gov.au/

# Set up workspace, all layers must be in prjGDB geodatabase
arcpy.env.workspace = prjGDB
arcpy.env.overwriteOutput = True
arcpy.env.addOutputsToMap = False
arcpy.env.snapRaster = ""  # Set snap raster here
arcpy.env.cellAlignment = ""  # Set prefered cell size here
arcpy.env.cellSize = "MINOF"
arcpy.env.mask = ""  # If using a mask set it here
arcpy.env.extent = ""  # Set extent of work area here
arcpy.env.outputCoordinateSystem = arcpy.SpatialReference("GCS_GDA_1994")

# Correlating Source layers with Ni mineralization
stateLayer = f"{prjGDB}/WA_Grid"
minLayer = f"{prjGDB}/Ni_mineralized_DHRC_grd"
featureLayer = [
    "um_poly",
    "um_fobs_grd",
    "um_chem_grd",
    "cumm_fobs_grd",
    "cumm_chem_grd",
    "Ticumm_chem_grd",
    "thltc_assc_grd",
    "flsc_assc_grd",
    "lip_poly",
]
exp = 0.00115

SpatialCorrelation(stateLayer, minLayer, featureLayer, exp)

# The following code creates Source predictor maps.
# Buffered distances and scores are determined by the outputs above
layer = "um_poly"
buffers = [0.01, 0.04]  # add your buffer ranges here
score = [2, 1]
BuffMap(layer, buffers, score)

layer = "um_fobs_grd"
buffers = [0.01, 0.03, 0.11]  # add your buffer ranges here
score = [3, 2, 1]  # add your realative prospectivity values here
BuffMap(layer, buffers, score)

layer = "um_chem_grd"
buffers = [0.01, 0.05, 0.08, 0.25]  # add your buffer ranges here
score = [4, 3, 2, 1]  # add your realative prospectivity values here
BuffMap(layer, buffers, score)

layer = "cumm_fobs_grd"
buffers = [0.01, 0.04, 0.09, 0.14, 0.24]  # add your buffer ranges here
score = [5, 4, 3, 2, 1]  # add your realative prospectivity values here
BuffMap(layer, buffers, score)

layer = "cumm_chem_grd"
buffers = [0.01, 0.03, 0.09, 0.14, 0.4]  # add your buffer ranges here
score = [5, 4, 3, 2, 1]  # add your realative prospectivity values here
BuffMap(layer, buffers, score)

layer = "Ticumm_chem_grd"
buffers = [0.04, 0.09, 0.13]  # add your buffer ranges here
score = [3, 2, 1]  # add your realative prospectivity values here
BuffMap(layer, buffers, score)

layer = "flsc_assc_grd"
buffers = [0.01, 0.04, 0.09]  # add your buffer ranges here
score = [3, 2, 1]  # add your realative prospectivity values here
BuffMap(layer, buffers, score)

layer = "thltc_assc_grd"
buffers = [0.01, 0.04, 0.08, 0.11, 0.26]  # add your buffer ranges here
score = [5, 4, 3, 2, 1]  # add your realative prospectivity values here
BuffMap(layer, buffers, score)

layer = "lip_poly"
buffers = [0.01, 0.08, 0.14]  # add your buffer ranges here
score = [3, 1, 2]  # add your realative prospectivity values here
BuffMap(layer, buffers, score)

# Correlating Pathway layers with Ni mineralization
stateLayer = f"{prjGDB}/WA_Grid"
minLayer = f"{prjGDB}/Ni_mineralized_DHRC_grd"
featureLayer = ["faults_shears", "faults_intersection", "FOBS_mantle_rock_grd"]
exp = 0.00115

SpatialCorrelation(stateLayer, minLayer, featureLayer, exp)

# The following code creates Pathways predictor maps.
# Buffered distances and scores are determined by the outputs above
layer = "faults_shears"
buffers = [0.01, 0.03, 0.04]  # add your buffer ranges here
score = [3, 2, 1]  # add your realative prospectivity values here
BuffMap(layer, buffers, score)

layer = "faults_intersection"  # Point file created using the Intersect tool in ArcGIS
buffers = [0.02, 0.04, 0.06, 0.09]  # add your buffer ranges here
score = [4, 3, 2, 1]  # add your realative prospectivity values here
BuffMap(layer, buffers, score)

layer = "FOBS_mantle_rock_grd"
buffers = [0.01, 0.03, 0.09]  # add your buffer ranges here
score = [3, 2, 1]  # add your realative prospectivity values here
BuffMap(layer, buffers, score)


# Determine preferential host rock
# Geology polygons are avalible as Geology_MERGED from DASC Minerals Systems Atlas Primary Data Layers: https://dasc.dmirs.wa.gov.au/
targetFeatures = f"{prjGDB}/Ni_mineralized_DHRC_grd"
outFeatureClass = f"{prjGDB}/HostRock"

hostRock = arcpy.analysis.SpatialJoin(targetFeatures, geologyMerged, outFeatureClass)

arr = arcpy.da.TableToNumPyArray(hostRock, "ROCKTYPE1")
name, count = np.unique(arr, return_counts=True)
df = pd.DataFrame([name, count]).transpose()
df.to_excel(f"{root_dir}/HostRock.xlsx")

targetFeatures = f"{prjGDB}/WA_Grid"
outFeatureClass = f"{prjGDB}/AllRock"

allRock = arcpy.analysis.SpatialJoin(targetFeatures, geologyMerged, outFeatureClass)

arr = arcpy.da.TableToNumPyArray("AllRock", "ROCKTYPE1")
name, count = np.unique(arr, return_counts=True)
df = pd.DataFrame([name, count]).transpose()
df.to_excel(f"{root_dir}/AllRock.xlsx")

hostLithQry = (
    "(ROCKTYPE1 LIKE 'igneous ultramafic volcanic' OR "
    "ROCKTYPE1 LIKE 'igneous ultramafic intrusive' OR "
    "ROCKTYPE1 LIKE 'meta-igneous ultramafic volcanic or undivided' OR "
    "ROCKTYPE1 LIKE 'metasedimentary other chemical' OR "
    "ROCKTYPE1 LIKE 'meta-igneous felsic volcanic' OR "
    "ROCKTYPE1 LIKE 'meta-igneous ultramafic intrusive' OR "
    "ROCKTYPE1 LIKE 'igneous felsic volcanic' OR "
    "ROCKTYPE1 LIKE 'igneous mafic intrusive' OR "
    "ROCKTYPE1 LIKE 'igneous mafic volcanic' OR "
    "ROCKTYPE1 LIKE 'metasedimentary carbonate' OR "
    "ROCKTYPE1 LIKE 'meta-igneous mafic volcanic' OR "
    "ROCKTYPE1 LIKE 'meta-igneous mafic intrusive' OR "
    "ROCKTYPE1 LIKE 'meta-igneous mafic' OR "
    "ROCKTYPE1 LIKE 'sedimentary other chemical or biochemical')"
)

hostLithology = arcpy.management.SelectLayerByAttribute(
    geologyMerged, "NEW_SELECTION", hostLithQry
)
arcpy.management.CopyFeatures(hostLithology, "HostLithology")

# Determine cell count
stateLayer = f"{prjGDB}/WA_grid"
minLayer = f"{prjGDB}/Ni_mineralized_DHRC_grd"
layers = [
    "HostLithology",
    "Ni_depleted_poly",
    "Ni_enriched_poly",
    "PGE_enriched_poly",
    "PGE_depleted_poly",
    "Cu_depleted_poly",
    "S_enriched_poly",
    "CuZr_enriched_poly",
    "CrustalContamination_poly",
]

AssociationCorrelation(stateLayer, minLayer, layers)

# Reclassify Host Lithology
# Reclassify can't convert between text and numeric values
# remap values are obtained from the attribute table of 'HostLithology_rst'
# each lithology in ROCKTYPE1 is assinged a unique "value".

polygon = arcpy.conversion.PolygonToRaster(
    "HostLithology", "ROCKTYPE1", "HostLithology_rst", "MAXIMUM_COMBINED_AREA"
)
remap = arcpy.sa.RemapValue(
    [
        [1, 1],
        [2, 1],
        [3, 1],
        [4, 3],
        [5, 1],
        [6, 3],
        [7, 1],
        [8, 3],
        [9, 1],
        [10, 1],
        [11, 2],
        [12, 1],
        [13, 1],
        [14, 1],
        ["NODATA", 0],
    ]
)

rstRcl = arcpy.sa.Reclassify(polygon, "Value", remap, "NODATA")
rstRcl.save("HostLithology_rst_rcl")


# Reclassify gravity. Relative prospectivity values based on correlation between Ni mineralization sites and gravity anomaly

gravRange = arcpy.sa.RemapRange(
    [
        [-1500, -715, 0],
        [-715, -583, 1],
        [-583, -395, 2],
        [-395, 25, 1],
        [25, 900, 0],
        ["NODATA", 0],
    ]
)
gravRcl = arcpy.sa.Reclassify("Grav_1km", "Value", gravRange, "NODATA")
gravRcl.save("Grav_1km_rcl")


# Reclassify polygon layers

rasters = [
    "Ni_enriched_poly",
    "Ni_depleted_poly",
    "CuZr_enriched_poly",
    "Cu_depleted_poly",
    "PGE_enriched_poly",
    "PGE_depleted_poly",
    "S_enriched_poly",
    "CrustalContamination_poly",
]

AsscRasters(rasters)
