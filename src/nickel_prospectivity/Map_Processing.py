# import libraries
import arcpy
import arcpy.sa
import numpy as np

from tools import RasterProcessing, RasterCorrelation, RasterROC, RasterTornado

arcpy.CheckOutExtension("Spatial")

# Set up folder pathways
root_dir = "C:/project_folder"  # Your project folder here
prjGDB = f"{root_dir}/Project.gdb"  # The .gdb for your project

# Set up workspace, all layers must be in prjGDB geodatabase
arcpy.env.workspace = prjGDB
arcpy.env.overwriteOutput = True
arcpy.env.addOutputsToMap = False
arcpy.env.snapRaster = ""  # If using a snap raster set it here.
arcpy.env.cellAlignment = ""  # Set prefered cell size here.
arcpy.env.cellSize = "MINOF"
arcpy.env.mask = ""  # If using a mask set it here
arcpy.env.extent = ""  # Set extent here
arcpy.env.outputCoordinateSystem = arcpy.SpatialReference("GCS_GDA_1994")


# Source
umPoly = f"{prjGDB}/um_poly_bff_rst_rcl"
umChem = f"{prjGDB}/um_chem_grd_bff_rst_rcl"
umFobs = f"{prjGDB}/um_fobs_grd_bff_rst_rcl"
cummFobs = f"{prjGDB}/cumm_fobs_grd_bff_rst_rcl"
cummChem = f"{prjGDB}/cumm_chem_grd_bff_rst_rcl"
cummTi = f"{prjGDB}/Ticumm_chem_grd_bff_rst_rcl"
flscAssc = f"{prjGDB}/flsc_assc_grd_bff_rst_rcl"
thltcAssc = f"{prjGDB}/thltc_assc_grd_bff_rst_rcl"
lipPoly = f"{prjGDB}/lip_poly_bff_rst_rcl"

# Pathway
faultShear = f"{prjGDB}/faults_shears_bff_rst_rcl"
faultShearInter = f"{prjGDB}/faults_intersection_bff_rst_rcl"
mntlRck = f"{prjGDB}/FOBS_mantle_rock_grd_bff_rst_rcl"

# Trap
hostLith = f"{prjGDB}/HostLithology_rst_rcl"
grav = f"{prjGDB}/Grav_1km_rcl"
niEnrch = f"{prjGDB}/Ni_enriched_poly_rst_rcl"
niDplt = f"{prjGDB}/Ni_depleted_poly_rst_rcl"
cuEnrch = f"{prjGDB}/CuZr_enriched_poly_rst_rcl"
cuDplt = f"{prjGDB}/Cu_depleted_poly_rst_rcl"
pgeEnrch = f"{prjGDB}/PGE_enriched_poly_rst_rcl"
pgeDplt = f"{prjGDB}/PGE_depleted_poly_rst_rcl"
crstCntmntn = f"{prjGDB}/CrustalContamination_poly_rst_rcl"
sulfEnrch = f"{prjGDB}/S_enriched_poly_rst_rcl"

# Weight assigned to each predictor map are based on its importance to Ni mineralization and confidence in the data
layers = [
    [umPoly, 0.49],
    [umChem, 0.58],
    [umFobs, 0.47],
    [cummFobs, 0.49],
    [cummChem, 0.59],
    [cummTi, 0.6],
    [flscAssc, 0.28],
    [thltcAssc, 0.29],
    [lipPoly, 0.4],
    [faultShear, 0.47],
    [faultShearInter, 0.47],
    [mntlRck, 0.31],
    [hostLith, 0.54],
    [grav, 0.5],
    [niEnrch, 0.45],
    [niDplt, 0.49],
    [cuEnrch, 0.37],
    [cuDplt, 0.42],
    [pgeEnrch, 0.46],
    [pgeDplt, 0.42],
    [crstCntmntn, 0.44],
    [sulfEnrch, 0.46],
]

for lyr, wt in layers:
    RasterProcessing(lyr, wt)

# Apply fuzzy OR function to each predictor map in each targetting element
# Apply fuzzy SUM function to each targetting element in each critical process
umLayers = [
    "um_chem_grd_bff_rst_rcl_fzz_con_wt",
    "um_fobs_grd_bff_rst_rcl_fzz_con_wt",
    "um_poly_bff_rst_rcl_fzz_con_wt",
]
fzzOvrly = arcpy.sa.FuzzyOverlay(umLayers, "OR")
fzzOvrly.save("UltramaficRock_FuzzOr")

cumulates = [
    "cumm_chem_grd_bff_rst_rcl_fzz_con_wt",
    "cumm_fobs_grd_bff_rst_rcl_fzz_con_wt",
    "Ticumm_chem_grd_bff_rst_rcl_fzz_con_wt",
]
fzzOvrly = arcpy.sa.FuzzyOverlay(cumulates, "OR")
fzzOvrly.save("Cumulate_FuzzOr")

bimodal = [
    "flsc_assc_grd_bff_rst_rcl_fzz_con_wt",
    "thltc_assc_grd_bff_rst_rcl_fzz_con_wt",
]
fzzOvrly = arcpy.sa.FuzzyOverlay(bimodal, "OR")
fzzOvrly.save("Bimodal_FuzzOr")

layers = [
    "UltramaficRock_FuzzOr",
    "Cumulate_FuzzOr",
    "Bimodal_FuzzOr",
    "lip_poly_bff_rst_rcl_fzz_con_wt",
]
fzzOvrly = arcpy.sa.FuzzyOverlay(layers, "SUM")
fzzOvrly.save("Ni_source_rst")

layers = [
    "faults_shears_bff_rst_rcl_fzz_con_wt",
    "faults_intersection_bff_rst_rcl_fzz_con_wt",
]

fzzOvrly = arcpy.sa.FuzzyOverlay(layers, "OR")
fzzOvrly.save("Faulting_FuzzOr")

layers = ["Faulting_FuzzOr", "FOBS_mantle_rock_grd_bff_rst_rcl_fzz_con_wt"]
fzzOvrly = arcpy.sa.FuzzyOverlay(layers, "SUM")
fzzOvrly.save("Ni_pathway_rst")

sRich = ["CuZr_enriched_poly_rst_rcl_fzz_con_wt", "S_enriched_poly_rst_rcl_fzz_con_wt"]

fzzOvrly = arcpy.sa.FuzzyOverlay(sRich, "OR")
fzzOvrly.save("S_Enriched_FuzzOr")

deplt = [
    "Ni_depleted_poly_rst_rcl_fzz_con_wt",
    "Cu_depleted_poly_rst_rcl_fzz_con_wt",
    "PGE_depleted_poly_rst_rcl_fzz_con_wt",
]

fzzOvrly = arcpy.sa.FuzzyOverlay(deplt, "OR")
fzzOvrly.save("depleted_FuzzOr")

enrch = ["Ni_enriched_poly_rst_rcl_fzz_con_wt", "PGE_enriched_poly_rst_rcl_fzz_con_wt"]

fzzOvrly = arcpy.sa.FuzzyOverlay(enrch, "OR")
fzzOvrly.save("enriched_FuzzOr")

layers = [
    "Grav_1km_rcl_fzz_con_wt",
    "HostLithology_rst_rcl_fzz_con_wt",
    "CrustalContamination_poly_rst_rcl_fzz_con_wt",
    "S_Enriched_FuzzOr",
    "depleted_FuzzOr",
    "enriched_FuzzOr",
]

fzzOvrly = arcpy.sa.FuzzyOverlay(layers, "SUM")
fzzOvrly.save("Ni_trap_rst")

layers = ["Ni_source_rst", "Ni_pathway_rst", "Ni_trap_rst"]

prospMap = arcpy.sa.FuzzyOverlay(layers, "GAMMA", 0.95)
prospMap.save("Ni_prospectivity_rst")

minResult = arcpy.GetRasterProperties_management(prospMap, "MINIMUM")
minValue = float(minResult.getOutput(0))
maxResult = arcpy.GetRasterProperties_management(prospMap, "MAXIMUM")
maxValue = float(maxResult.getOutput(0))

norm = (arcpy.sa.Raster(prospMap) - minValue) / (maxValue - minValue)
norm.save("Ni_prospectivity_rst_norm")

# Intersects the location of Ni mines and barren locations with the propsecitivty map and records the value
stateLayer = f"{prjGDB}/Ni_barren_DHRC_pt"
minLayer = f"{prjGDB}/Ni_MinesDeposits"
rasterLayer = "Ni_prospectivity_rst_norm"

RasterROC(stateLayer, minLayer, rasterLayer)

# Intersects the location of Ni mines and points of state grid with the propsecitivty map and records the value
stateLayer = f"{prjGDB}/WA_Points"
minLayer = f"{prjGDB}/Ni_MinesDeposits"
rasterLayer = "Ni_prospectivity_rst_norm"
rasterMin = 0
rasterMax = 1

RasterCorrelation(stateLayer, minLayer, rasterLayer, rasterMin, rasterMax)

# Monte carlo simulation to test the effect of varying weights on the predictor maps
layer = [
    f"{prjGDB}/um_poly_bff_rst_rcl_fzz_con",
    f"{prjGDB}/um_chem_grd_bff_rst_rcl_fzz_con",
    f"{prjGDB}/um_fobs_grd_bff_rst_rcl_fzz_con",
    f"{prjGDB}/cumm_fobs_grd_bff_rst_rcl_fzz_con",
    f"{prjGDB}/cumm_chem_grd_bff_rst_rcl_fzz_con",
    f"{prjGDB}/Ticumm_chem_grd_bff_rst_rcl_fzz_con",
    f"{prjGDB}/flsc_assc_grd_bff_rst_rcl_fzz_con",
    f"{prjGDB}/thltc_assc_grd_bff_rst_rcl_fzz_con",
    f"{prjGDB}/lip_poly_bff_rst_rcl_fzz_con",
    f"{prjGDB}/faults_shears_bff_rst_rcl_fzz_con",
    f"{prjGDB}/faults_intersection_bff_rst_rcl_fzz_con",
    f"{prjGDB}/FOBS_mantle_rock_grd_bff_rst_rcl_fzz_con",
    f"{prjGDB}/CrustalContamination_poly_rst_rcl_fzz_con",
    f"{prjGDB}/HostLithology_rst_rcl_fzz_con",
    f"{prjGDB}/Grav_1km_rcl_fzz_con",
    f"{prjGDB}/Ni_depleted_poly_rst_rcl_fzz_con",
    f"{prjGDB}/Cu_depleted_poly_rst_rcl_fzz_con",
    f"{prjGDB}/PGE_depleted_poly_rst_rcl_fzz_con",
    f"{prjGDB}/Ni_enriched_poly_rst_rcl_fzz_con",
    f"{prjGDB}/PGE_enriched_poly_rst_rcl_fzz_con",
    f"{prjGDB}/S_enriched_poly_rst_rcl_fzz_con",
    f"{prjGDB}/CuZr_enriched_poly_rst_rcl_fzz_con",
]

for lyr in layer:
    arcpy.management.Resample(lyr, lyr + "_rs", 0.03, "BILINEAR")

sc_um1_arr = arcpy.RasterToNumPyArray(
    f"{prjGDB}/um_poly_bff_rst_rcl_fzz_con_rs", nodata_to_value=0
)
sc_um2_arr = arcpy.RasterToNumPyArray(
    f"{prjGDB}/um_chem_grd_bff_rst_rcl_fzz_con_rs", nodata_to_value=0
)
sc_um3_arr = arcpy.RasterToNumPyArray(
    f"{prjGDB}/um_fobs_grd_bff_rst_rcl_fzz_con_rs", nodata_to_value=0
)
sc_cm1_arr = arcpy.RasterToNumPyArray(
    f"{prjGDB}/cumm_fobs_grd_bff_rst_rcl_fzz_con_rs", nodata_to_value=0
)
sc_cm2_arr = arcpy.RasterToNumPyArray(
    f"{prjGDB}/cumm_chem_grd_bff_rst_rcl_fzz_con_rs", nodata_to_value=0
)
sc_cm3_arr = arcpy.RasterToNumPyArray(
    f"{prjGDB}/Ticumm_chem_grd_bff_rst_rcl_fzz_con_rs", nodata_to_value=0
)
sc_flsc_arr = arcpy.RasterToNumPyArray(
    f"{prjGDB}/flsc_assc_grd_bff_rst_rcl_fzz_con_rs", nodata_to_value=0
)
sc_thltc_arr = arcpy.RasterToNumPyArray(
    f"{prjGDB}/thltc_assc_grd_bff_rst_rcl_fzz_con_rs", nodata_to_value=0
)
sc_lip_arr = arcpy.RasterToNumPyArray(
    f"{prjGDB}/lip_poly_bff_rst_rcl_fzz_con_rs", nodata_to_value=0
)
pth1_arr = arcpy.RasterToNumPyArray(
    f"{prjGDB}/faults_shears_bff_rst_rcl_fzz_con_rs", nodata_to_value=0
)
pth2_arr = arcpy.RasterToNumPyArray(
    f"{prjGDB}/faults_intersection_bff_rst_rcl_fzz_con_rs", nodata_to_value=0
)
pth3_arr = arcpy.RasterToNumPyArray(
    f"{prjGDB}/FOBS_mantle_rock_grd_bff_rst_rcl_fzz_con_rs", nodata_to_value=0
)
tp_cc_arr = arcpy.RasterToNumPyArray(
    f"{prjGDB}/CrustalContamination_poly_rst_rcl_fzz_con_rs", nodata_to_value=0
)
tp_hl_arr = arcpy.RasterToNumPyArray(
    f"{prjGDB}/HostLithology_rst_rcl_fzz_con_rs", nodata_to_value=0
)
tp_gr_arr = arcpy.RasterToNumPyArray(
    f"{prjGDB}/Grav_1km_rcl_fzz_con_rs", nodata_to_value=0
)
tp_d1_arr = arcpy.RasterToNumPyArray(
    f"{prjGDB}/Ni_depleted_poly_rst_rcl_fzz_con_rs", nodata_to_value=0
)
tp_d2_arr = arcpy.RasterToNumPyArray(
    f"{prjGDB}/Cu_depleted_poly_rst_rcl_fzz_con_rs", nodata_to_value=0
)
tp_d3_arr = arcpy.RasterToNumPyArray(
    f"{prjGDB}/PGE_depleted_poly_rst_rcl_fzz_con_rs", nodata_to_value=0
)
tp_e1_arr = arcpy.RasterToNumPyArray(
    f"{prjGDB}/Ni_enriched_poly_rst_rcl_fzz_con_rs", nodata_to_value=0
)
tp_e2_arr = arcpy.RasterToNumPyArray(
    f"{prjGDB}/PGE_enriched_poly_rst_rcl_fzz_con_rs", nodata_to_value=0
)
tp_s1_arr = arcpy.RasterToNumPyArray(
    f"{prjGDB}/S_enriched_poly_rst_rcl_fzz_con_rs", nodata_to_value=0
)
tp_s2_arr = arcpy.RasterToNumPyArray(
    f"{prjGDB}/CuZr_enriched_poly_rst_rcl_fzz_con_rs", nodata_to_value=0
)

values = []
sourceValues = []
pathwayValues = []
trapValues = []
sc_um1Values = []
sc_um2Values = []
sc_um3Values = []
sc_cm1Values = []
sc_cm2Values = []
sc_cm3Values = []
sc_flscValues = []
sc_thltcValues = []
sc_lipValues = []
pth1Values = []
pth2Values = []
pth3Values = []
tp_ccValues = []
tp_hlValues = []
tp_grValues = []
tp_d1Values = []
tp_d2Values = []
tp_d3Values = []
tp_e1Values = []
tp_e2Values = []
tp_s1Values = []
tp_s2Values = []

# This generates 500 maps which randomoly sample weights from a triangualr distribtuion.
# Distribution defined by the min, max, and avg of weights given by geologist as a part of this project
for _ in range(500):
    sc_um1_wt = sc_um1_arr * float(np.random.triangular(0.22, 0.49, 0.65, 1))
    sc_um2_wt = sc_um2_arr * float(np.random.triangular(0.45, 0.58, 0.73, 1))
    sc_um3_wt = sc_um3_arr * float(np.random.triangular(0.39, 0.47, 0.57, 1))
    sc_um = np.maximum.reduce([sc_um1_wt, sc_um2_wt, sc_um3_wt])
    sc_cm1_wt = sc_um1_arr * float(np.random.triangular(0.46, 0.49, 0.50, 1))
    sc_cm2_wt = sc_um2_arr * float(np.random.triangular(0.58, 0.59, 0.61, 1))
    sc_cm3_wt = sc_um3_arr * float(np.random.triangular(0.58, 0.60, 0.65, 1))
    sc_cm = np.maximum.reduce([sc_cm1_wt, sc_cm2_wt, sc_cm3_wt])
    sc_flsc_wt = sc_flsc_arr * float(np.random.triangular(0.15, 0.28, 0.39, 1))
    sc_thltc_wt = sc_thltc_arr * float(np.random.triangular(0.19, 0.29, 0.34, 1))
    sc_bv = np.maximum.reduce([sc_flsc_wt, sc_thltc_wt])
    sc_lip_wt = sc_lip_arr * float(np.random.triangular(0.25, 0.4, 0.51, 1))
    source = 1 - ((1 - sc_um) * (1 - sc_cm) * (1 - sc_bv) * (1 - sc_lip_wt))

    pth1_wt = pth1_arr * float(np.random.triangular(0.45, 0.47, 0.51, 1))
    pth2_wt = pth2_arr * float(np.random.triangular(0.39, 0.47, 0.58, 1))
    pth_fl = np.maximum.reduce([pth1_wt, pth2_wt])
    pth3_wt = pth3_arr * float(np.random.triangular(0.19, 0.31, 0.45, 1))
    pathway = 1 - ((1 - pth_fl) * (1 - pth3_wt))

    tp_s1_wt = tp_s1_arr * float(np.random.triangular(0.23, 0.46, 0.77, 1))
    tp_s2_wt = tp_s2_arr * float(np.random.triangular(0.19, 0.37, 0.50, 1))
    tp_s = np.maximum.reduce([tp_s1_wt, tp_s2_wt])
    tp_d1_wt = tp_d1_arr * float(np.random.triangular(0.18, 0.49, 0.73, 1))
    tp_d2_wt = tp_d2_arr * float(np.random.triangular(0.18, 0.42, 0.56, 1))
    tp_d3_wt = tp_d3_arr * float(np.random.triangular(0.18, 0.42, 0.56, 1))
    tp_d = np.maximum.reduce([tp_d1_wt, tp_d2_wt, tp_d3_wt])
    tp_e1_wt = tp_e1_arr * float(np.random.triangular(0.18, 0.45, 0.61, 1))
    tp_e2_wt = tp_e2_arr * float(np.random.triangular(0.18, 0.46, 0.65, 1))
    tp_e = np.maximum.reduce([tp_e1_wt, tp_e2_wt])
    tp_cc_wt = tp_cc_arr * float(np.random.triangular(0.34, 0.44, 0.50, 1))
    tp_hl_wt = tp_hl_arr * float(np.random.triangular(0.45, 0.54, 0.73, 1))
    tp_gr_wt = tp_gr_arr * float(np.random.triangular(0.48, 0.50, 0.51, 1))
    trap = 1 - (
        (1 - tp_s)
        * (1 - tp_d)
        * (1 - tp_e)
        * (1 - tp_cc_wt)
        * (1 - tp_hl_wt)
        * (1 - tp_gr_wt)
    )

    fzzsum = 1 - ((1 - source) * (1 - pathway) * (1 - trap))
    fzzprd = source * pathway * trap
    gmm = 0.95
    fzzgmm = (fzzsum**gmm) * (fzzprd ** (1 - gmm))

    values.append(fzzgmm)
    sourceValues.append(source)
    pathwayValues.append(pathway)
    trapValues.append(trap)
    sc_um1Values.append(sc_um1_wt)
    sc_um2Values.append(sc_um2_wt)
    sc_um3Values.append(sc_um3_wt)
    sc_cm1Values.append(sc_cm1_wt)
    sc_cm2Values.append(sc_cm2_wt)
    sc_cm3Values.append(sc_cm3_wt)
    sc_flscValues.append(sc_flsc_wt)
    sc_thltcValues.append(sc_thltc_wt)
    sc_lipValues.append(sc_lip_wt)
    pth1Values.append(pth1_wt)
    pth2Values.append(pth2_wt)
    pth3Values.append(pth3_wt)
    tp_ccValues.append(tp_cc_wt)
    tp_hlValues.append(tp_hl_wt)
    tp_grValues.append(tp_gr_wt)
    tp_d1Values.append(tp_d1_wt)
    tp_d2Values.append(tp_d2_wt)
    tp_d3Values.append(tp_d3_wt)
    tp_e1Values.append(tp_e1_wt)
    tp_e2Values.append(tp_e2_wt)
    tp_s1Values.append(tp_s1_wt)
    tp_s2Values.append(tp_s2_wt)

arr_Values = np.array(values)
arr_SourceValues = np.array(sourceValues)
arr_PathwayValue = np.array(pathwayValues)
arr_TrapValues = np.array(trapValues)
arr_sc_um1Values = np.array(sc_um1Values)
arr_sc_um2Values = np.array(sc_um2Values)
arr_sc_um3Values = np.array(sc_um3Values)
arr_sc_cm1Values = np.array(sc_cm1Values)
arr_sc_cm2Values = np.array(sc_cm2Values)
arr_sc_cm3Values = np.array(sc_cm3Values)
arr_sc_flscValues = np.array(sc_flscValues)
arr_sc_thltcValues = np.array(sc_thltcValues)
arr_sc_lipValues = np.array(sc_lipValues)
arr_pth1Values = np.array(pth1Values)
arr_pth2Values = np.array(pth2Values)
arr_pth3Values = np.array(pth3Values)
arr_tp_ccValues = np.array(tp_ccValues)
arr_tp_hlValues = np.array(tp_hlValues)
arr_tp_grValues = np.array(tp_grValues)
arr_tp_d1Values = np.array(tp_d1Values)
arr_tp_d2Values = np.array(tp_d2Values)
arr_tp_d3Values = np.array(tp_d3Values)
arr_tp_e1Values = np.array(tp_e1Values)
arr_tp_e2Values = np.array(tp_e2Values)
arr_tp_s1Values = np.array(tp_s1Values)
arr_tp_s2Values = np.array(tp_s2Values)

inRas = arcpy.Raster(
    f"{prjGDB}/um_poly_bff_rst_rcl_fzz_con_rs"
)  # chose a raster with the extents that match the desired output rasters
lowerLeft = arcpy.Point(inRas.extent.XMin, inRas.extent.YMin)

arrays = [
    [arr_Values, "Map"],
    [arr_SourceValues, "Source"],
    [arr_PathwayValue, "Pathway"],
    [arr_TrapValues, "Trap"],
    [arr_sc_um1Values, "sc_um1"],
    [arr_sc_um2Values, "sc_um2"],
    [arr_sc_um3Values, "sc_um3"],
    [arr_sc_cm1Values, "sc_cm1"],
    [arr_sc_cm2Values, "sc_cm2"],
    [arr_sc_cm3Values, "sc_cm3"],
    [arr_sc_flscValues, "sc_flsc"],
    [arr_sc_thltcValues, "sc_thlt"],
    [arr_sc_lipValues, "sc_lip"],
    [arr_pth1Values, "pth1"],
    [arr_pth2Values, "pth2"],
    [arr_pth3Values, "pth3"],
    [arr_tp_ccValues, "tp_cc"],
    [arr_tp_hlValues, "tp_hl"],
    [arr_tp_grValues, "tp_gr"],
    [arr_tp_d1Values, "tp_d1"],
    [arr_tp_d2Values, "tp_d2"],
    [arr_tp_d3Values, "tp_d3"],
    [arr_tp_e1Values, "tp_e1"],
    [arr_tp_e2Values, "tp_e2"],
    [arr_tp_s1Values, "tp_s1"],
    [arr_tp_s2Values, "tp_s2"],
]

for arr, name in arrays:
    p10 = np.percentile(arr, 10, axis=0)
    p10Raster = arcpy.NumPyArrayToRaster(p10, lowerLeft, 0.03)
    p10Raster.save(f"{prjGDB}/{name}_p10_map")

    med = np.median(arr, axis=0)
    medRaster = arcpy.NumPyArrayToRaster(med, lowerLeft, 0.03)
    medRaster.save(f"{prjGDB}/{name}_p50_map")

    p90 = np.percentile(arr, 90, axis=0)
    p90Raster = arcpy.NumPyArrayToRaster(p90, lowerLeft, 0.03)
    p90Raster.save(f"{prjGDB}/{name}_p90_map")


minLayer = f"{prjGDB}/Ni_MinesDeposits"
maps = [
    "Map_p10_map",
    "Source_p10_map",
    "Pathway_p10_map",
    "Trap_p10_map",
    "sc_um1_p10_map",
    "sc_um2_p10_map",
    "sc_um3_p10_map",
    "sc_cm1_p10_map",
    "sc_cm2_p10_map",
    "sc_cm3_p10_map",
    "sc_flsc_p10_map",
    "sc_thlt_p10_map",
    "sc_lip_p10_map",
    "pth1_p10_map",
    "pth2_p10_map",
    "pth3_p10_map",
    "tp_cc_p10_map",
    "tp_hl_p10_map",
    "tp_gr_p10_map",
    "tp_d1_p10_map",
    "tp_d2_p10_map",
    "tp_d3_p10_map",
    "tp_e1_p10_map",
    "tp_e2_p10_map",
    "tp_s1_p10_map",
    "tp_s2_p10_map",
    "Map_p50_map",
    "Source_p50_map",
    "Pathway_p50_map",
    "Trap_p50_map",
    "sc_um1_p50_map",
    "sc_um2_p50_map",
    "sc_um3_p50_map",
    "sc_cm1_p50_map",
    "sc_cm2_p50_map",
    "sc_cm3_p50_map",
    "sc_flsc_p50_map",
    "sc_thlt_p50_map",
    "sc_lip_p50_map",
    "pth1_p50_map",
    "pth2_p50_map",
    "pth3_p50_map",
    "tp_cc_p50_map",
    "tp_hl_p50_map",
    "tp_gr_p50_map",
    "tp_d1_p50_map",
    "tp_d2_p50_map",
    "tp_d3_p50_map",
    "tp_e1_p50_map",
    "tp_e2_p50_map",
    "tp_s1_p50_map",
    "tp_s2_p50_map",
    "Map_p90_map",
    "Source_p90_map",
    "Pathway_p90_map",
    "Trap_p90_map",
    "sc_um1_p90_map",
    "sc_um2_p90_map",
    "sc_um3_p90_map",
    "sc_cm1_p90_map",
    "sc_cm2_p90_map",
    "sc_cm3_p90_map",
    "sc_flsc_p90_map",
    "sc_thlt_p90_map",
    "sc_lip_p90_map",
    "pth1_p90_map",
    "pth2_p90_map",
    "pth3_p90_map",
    "tp_cc_p90_map",
    "tp_hl_p90_map",
    "tp_gr_p90_map",
    "tp_d1_p90_map",
    "tp_d2_p90_map",
    "tp_d3_p90_map",
    "tp_e1_p90_map",
    "tp_e2_p90_map",
    "tp_s1_p90_map",
    "tp_s2_p90_map",
]
RasterTornado(minLayer, maps)
