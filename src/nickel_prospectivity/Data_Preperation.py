# Import libraries
import arcpy
import arcpy.sa

from tools import AnnMonteCarlo, AsscMonteCarlo

# Set up folder pathways
root_dir = "C:/project_folder"  # Your project folder here
prjGDB = f"{root_dir}/Project.gdb"  # The .gdb for your project

# Set up workspace, all layers must be in prjGDB geodatabase
arcpy.env.workspace = prjGDB
arcpy.env.overwriteOutput = True
arcpy.env.addOutputsToMap = False
arcpy.env.outputCoordinateSystem = arcpy.SpatialReference("GCS_GDA_1994")

# Set file pathways
stateGrid = (
    f"{prjGDB}/WA_Grid"  # a 1km grid was created using the fishnet tool in ArcGIS
)
statePoints = f"{prjGDB}/WA_Points"  # a 1km point grid was created using the fishnet tool in ArcGIS
mindex = f"{prjGDB}/MINEDEX_DERIVED"  # Download from DASC Minerals Systems Atlas Primary Data Layers: https://dasc.dmirs.wa.gov.au/
dh = f"{prjGDB}/ni_dh_nsmax_point"  # Download from DASC Minerals Systems Atlas Ni_Cu_PGE: https://dasc.dmirs.wa.gov.au/
wachem = f"{prjGDB}/Geochemistry_RECALC"  # Download from DASC Minerals Systems Atlas Primary Data Layers: https://dasc.dmirs.wa.gov.au/
geologyMerged = f"{prjGDB}/Geology_MERGED"  # Download from DASC Minerals Systems Atlas Primary Data Layers: https://dasc.dmirs.wa.gov.au/
fobs = f"{prjGDB}/WAROX_Lithology"  # Download from DASC Geology WAROX Lithology: https://dasc.dmirs.wa.gov.au/
dykes = f"{prjGDB}/interpdykel_2_5m"  # Download from DASC Minerals Systems Atlas Ni_Cu_PGE: https://dasc.dmirs.wa.gov.au/
lineStrct = f"{prjGDB}/Geology_linear_structures_500K"  # Download from DASC Minerals Systems Atlas Primary Data Layers: https://dasc.dmirs.wa.gov.au/

# Create standard layers

# The following code filters a list of Ni mines and deposits from Mindex derived.
# This can be downloaded from DASC Minerals Systems Atlas Primary Data Layers: https://dasc.dmirs.wa.gov.au/
niMine = """
(SITE_TYPE = 'Mine' Or SITE_TYPE = 'Deposit') 
And Ni_SITES IS NOT NULL
"""
niDep = arcpy.management.SelectLayerByAttribute(mindex, "NEW_SELECTION", niMine)
arcpy.management.CopyFeatures(niDep, f"{prjGDB}/Ni_MinesDeposits")


# Creates a grid representing the spatial distribution of Ni mineralization across the state
# First dh are filtered for RC or DD holes that are > 3000 ppm Ni.
# The prefiltered data is avalible from DASC Minerals Systems Atlas Magmatic Nickel: https://dasc.dmirs.wa.gov.au/
# These drillholes are intersected with a 1km2 grid and points created using the fishnet tool in ArcGIS.
# You will need to create your own grid and points to replicate this step
dhMin = """
Ni_ppm >= 3000 AND 
(HOLETYPE = 'DD' OR 
HOLETYPE = 'RC')
"""
minDH = arcpy.management.SelectLayerByAttribute(dh, "NEW_SELECTION", dhMin)
arcpy.management.CopyFeatures(minDH, f"{prjGDB}/Ni_mineralized_DHRC")
minDHGrd = arcpy.management.SelectLayerByLocation(
    stateGrid,
    "INTERSECT",
    f"{prjGDB}/Ni_mineralized_DHRC",
    selection_type="NEW_SELECTION",
)
arcpy.management.CopyFeatures(minDHGrd, f"{prjGDB}/Ni_mineralized_DHRC_grd")
minDHPt = arcpy.management.SelectLayerByLocation(
    statePoints,
    "INTERSECT",
    f"{prjGDB}/Ni_mineralized_DHRC_grd",
    selection_type="NEW_SELECTION",
)
arcpy.management.CopyFeatures(minDHPt, f"{prjGDB}/Ni_mineralized_DHRC_pt")

# Creates a grid representing the spatial distribution of barren Ni areas across the state
# First dh are filtered for RC or DD holes that are < 60 ppm Ni (average crustal value).
# These drillholes are intersected with a 1km2 grid and points created using the fishnet tool in ArcGIS.
# You will need to create your own grid and points to replicate this step
# Any barren cell with < 2 Ni mineralized drillholes are removed from the sample
# Any barren cell that intersects a Ni mineralized cell, or a cell with a Ni  mine or deposit is removed from the sample
# Any barren cell that < 2km from a Ni mineralized cell, or a cell with a Ni  mine or deposit is removed from the sample
dhBrrn = """
Ni_ppm <= 60 AND 
(HOLETYPE = 'DD' OR 
HOLETYPE = 'RC')
"""
brrnDH = arcpy.management.SelectLayerByAttribute(dh, "NEW_SELECTION", dhBrrn)
arcpy.management.CopyFeatures(brrnDH, f"{prjGDB}\\Ni_barren_DHRC")
arcpy.analysis.SpatialJoin(
    f"{prjGDB}\\WA_Grid", f"{prjGDB}\\Ni_barren_DHRC", f"{prjGDB}\\Ni_Barren_cells"
)
dhBrrn2 = """Join_Count >= 2"""
brrnDH2 = arcpy.management.SelectLayerByAttribute(
    f"{prjGDB}\\Ni_Barren_cells", "NEW_SELECTION", dhBrrn2
)
brrnDH2km = arcpy.management.SelectLayerByLocation(
    brrnDH2,
    "WITHIN_A_DISTANCE",
    f"{prjGDB}\\Ni_mineralized_DHRC_grd",
    0.02,
    "REMOVE_FROM_SELECTION",
)
brrnDH2kmMine = arcpy.management.SelectLayerByLocation(
    brrnDH2km,
    "WITHIN_A_DISTANCE",
    f"{prjGDB}\\Ni_MinesDeposits",
    0.02,
    "REMOVE_FROM_SELECTION",
)
arcpy.management.CopyFeatures(brrnDH2km, f"{prjGDB}\\Ni_barren_DHRC_grd")
brrnDHPt = arcpy.management.SelectLayerByLocation(
    f"{prjGDB}\\WA_Points",
    "INTERSECT",
    f"{prjGDB}\\Ni_barren_DHRC_grd",
    selection_type="NEW_SELECTION",
)
arcpy.management.CopyFeatures(brrnDHPt, f"{prjGDB}\\Ni_barren_DHRC_pt")

# The following code filters WACHEM data for mafic lithologies which have SiO2 between 0 and 57 wt% and MgO greater than 3 wt%
# # This ensures major elements are included which allows for Mg number to be calculated, and that truly mafic rocks are included
wachemQuery = """
(LOWER (LITHNAME) LIKE LOWER ('%chromitit%') OR
LOWER (LITHNAME) LIKE LOWER ('%dunit%') OR
LOWER (LITHNAME) LIKE LOWER ('%harzburgit%') OR
LOWER (LITHNAME) LIKE LOWER ('%hornblendit') OR
LOWER (LITHNAME) LIKE LOWER ('%komatiit%') OR
LOWER (LITHNAME) LIKE LOWER ('%lherzolit%') OR
LOWER (LITHNAME) LIKE LOWER ('%magnetitit%') OR
LOWER (LITHNAME) LIKE LOWER ('%peridotit%') OR
LOWER (LITHNAME) LIKE LOWER ('%pyroxenit%') OR
LOWER (LITHNAME) LIKE LOWER ('%serpentin%') OR
LOWER (LITHNAME) LIKE LOWER ('%tremolit%') OR
LOWER (LITHNAME) LIKE LOWER ('%mafic%') OR
LOWER (LITHNAME) LIKE LOWER ('%websterit%') OR
LOWER (LITHNAME) LIKE LOWER ('%wehrlit%') OR
LOWER (LITHNAME) LIKE LOWER ('%amphibolite%') OR
LOWER (LITHNAME) LIKE LOWER ('%anorthosite%') OR
LOWER (LITHNAME) LIKE LOWER ('%gabbro%') OR
LOWER (LITHNAME) LIKE LOWER ('%dolerite%') OR
LOWER (LITHNAME) LIKE LOWER ('%norite%') OR
LOWER (LITHNAME) LIKE LOWER ('%basalt%')) AND
(SiO2_pct > 0 AND SiO2_pct <= 57 AND MgO_pct >= 3)  
"""
wachemFilter = arcpy.management.SelectLayerByAttribute(
    wachem, "NEW_SELECTION", wachemQuery
)
arcpy.management.CopyFeatures(wachemFilter, f"{prjGDB}/WACHEM_filter")

# Claculating Mg number. First Fe2O3T is converted to FeO (Fe2O3T*0.89998)
# FeO and MgO are converted to Fe and Mg molar fraction (FeO/71.84 and MgO/40.31). Mg number is a ratio of Mg/Mg+Fe expressed as a percentage
arcpy.management.AddFields(
    f"{prjGDB}/WACHEM_filter",
    [["FeO", "DOUBLE"], ["FeMol", "DOUBLE"], ["MgMol", "DOUBLE"], ["MgNum", "DOUBLE"]],
)
fields = ["FeO", "FeMol", "MgMol", "MgNum"]
expressions = [
    "!Fe2O3T_pct!*0.89998",
    "!FeO!/71.84",
    "!MgO_pct!/40.31",
    "(!MgMol!/(!MgMol!+!FeMol!)*100)",
]

for field, expression in zip(fields, expressions):
    arcpy.management.CalculateField(
        f"{prjGDB}/WACHEM_filter", field, expression, "PYTHON3"
    )
arcpy.conversion.TableToExcel(
    f"{prjGDB}/WACHEM_filter", f"{root_dir}/WachemFilter.xlsx"
)

# Create layers to test spatial correlation with Ni mineralization

# The following code filters geology polygons for ultramafic lithologies
# Geology polygons are avalible as Geology_MERGED from DASC Minerals Systems Atlas Primary Data Layers: https://dasc.dmirs.wa.gov.au/
umQuery = """
ROCKTYPE1 LIKE '%igneous ultramafic%' OR 
ROCKTYPE2 LIKE '%igneous ultramafic%'
"""
umPoly = arcpy.management.SelectLayerByAttribute(
    geologyMerged, "NEW_SELECTION", umQuery
)
arcpy.management.CopyFeatures(umPoly, f"{prjGDB}/um_poly")

# The following code filters field observations for ultramafic litholgies and minerals
umQuery = """
(LOWER (OBSNOTES) LIKE LOWER ('%olivine%') OR
LOWER (OBSNOTES) LIKE LOWER ('%chromite%') OR
LOWER (OBSNOTES) LIKE LOWER ('%pyroxene%') OR
LOWER (OBSNOTES) LIKE LOWER ('%chromitit%') OR
LOWER (OBSNOTES) LIKE LOWER ('%dunit%') OR
LOWER (OBSNOTES) LIKE LOWER ('%harzburgit%') OR
LOWER (OBSNOTES) LIKE LOWER ('%hornblendit') OR
LOWER (OBSNOTES) LIKE LOWER ('%komatiit%') OR
LOWER (OBSNOTES) LIKE LOWER ('%lherzolit%') OR
LOWER (OBSNOTES) LIKE LOWER ('%magnetitit%') OR
LOWER (OBSNOTES) LIKE LOWER ('%peridotit%') OR
LOWER (OBSNOTES) LIKE LOWER ('%pyroxenit%') OR
LOWER (OBSNOTES) LIKE LOWER ('%serpentin%') OR
LOWER (OBSNOTES) LIKE LOWER ('%tremolit%') OR
LOWER (OBSNOTES) LIKE LOWER ('%ultramafic%') OR
LOWER (OBSNOTES) LIKE LOWER ('%websterit%') OR
LOWER (OBSNOTES) LIKE LOWER ('%wehrlit%')) AND
(LOWER (LITH_DESC) LIKE LOWER ('%olivine%') OR
LOWER (LITH_DESC) LIKE LOWER ('%chromite%') OR
LOWER (LITH_DESC) LIKE LOWER ('%pyroxene%') OR
LOWER (LITH_DESC) LIKE LOWER ('%chromitit%') OR
LOWER (LITH_DESC) LIKE LOWER ('%dunit%') OR
LOWER (LITH_DESC) LIKE LOWER ('%harzburgit%') OR
LOWER (LITH_DESC) LIKE LOWER ('%hornblendit') OR
LOWER (LITH_DESC) LIKE LOWER ('%komatiit%') OR
LOWER (LITH_DESC) LIKE LOWER ('%lherzolit%') OR
LOWER (LITH_DESC) LIKE LOWER ('%magnetitit%') OR
LOWER (LITH_DESC) LIKE LOWER ('%peridotit%') OR
LOWER (LITH_DESC) LIKE LOWER ('%pyroxenit%') OR
LOWER (LITH_DESC) LIKE LOWER ('%serpentin%') OR
LOWER (LITH_DESC) LIKE LOWER ('%tremolit%') OR
LOWER (LITH_DESC) LIKE LOWER ('%ultramafic%') OR
LOWER (LITH_DESC) LIKE LOWER ('%websterit%') OR
LOWER (LITH_DESC) LIKE LOWER ('%wehrlit%')) AND
(LOWER (ROCKTYPE) LIKE LOWER ('%olivine%') OR
LOWER (ROCKTYPE) LIKE LOWER ('%chromite%') OR
LOWER (ROCKTYPE) LIKE LOWER ('%pyroxene%') OR
LOWER (ROCKTYPE) LIKE LOWER ('%chromitit%') OR
LOWER (ROCKTYPE) LIKE LOWER ('%dunit%') OR
LOWER (ROCKTYPE) LIKE LOWER ('%harzburgit%') OR
LOWER (ROCKTYPE) LIKE LOWER ('%hornblendit') OR
LOWER (ROCKTYPE) LIKE LOWER ('%komatiit%') OR
LOWER (ROCKTYPE) LIKE LOWER ('%lherzolit%') OR
LOWER (ROCKTYPE) LIKE LOWER ('%magnetitit%') OR
LOWER (ROCKTYPE) LIKE LOWER ('%peridotit%') OR
LOWER (ROCKTYPE) LIKE LOWER ('%pyroxenit%') OR
LOWER (ROCKTYPE) LIKE LOWER ('%serpentin%') OR
LOWER (ROCKTYPE) LIKE LOWER ('%tremolit%') OR
LOWER (ROCKTYPE) LIKE LOWER ('%ultramafic%') OR
LOWER (ROCKTYPE) LIKE LOWER ('%websterit%') OR
LOWER (ROCKTYPE) LIKE LOWER ('%wehrlit%') OR
LOWER (ROCKNAME) LIKE LOWER ('%olivine%') OR
LOWER (ROCKNAME) LIKE LOWER ('%chromite%') OR
LOWER (ROCKNAME) LIKE LOWER ('%pyroxene%') OR
LOWER (ROCKNAME) LIKE LOWER ('%chromitit%') OR
LOWER (ROCKNAME) LIKE LOWER ('%dunit%') OR
LOWER (ROCKNAME) LIKE LOWER ('%harzburgit%') OR
LOWER (ROCKNAME) LIKE LOWER ('%hornblendit') OR
LOWER (ROCKNAME) LIKE LOWER ('%komatiit%') OR
LOWER (ROCKNAME) LIKE LOWER ('%lherzolit%') OR
LOWER (ROCKNAME) LIKE LOWER ('%magnetitit%') OR
LOWER (ROCKNAME) LIKE LOWER ('%peridotit%') OR
LOWER (ROCKNAME) LIKE LOWER ('%pyroxenit%') OR
LOWER (ROCKNAME) LIKE LOWER ('%serpentin%') OR
LOWER (ROCKNAME) LIKE LOWER ('%tremolit%') OR
LOWER (ROCKNAME) LIKE LOWER ('%ultramafic%') OR
LOWER (ROCKNAME) LIKE LOWER ('%websterit%') OR
LOWER (ROCKNAME) LIKE LOWER ('%wehrlit%'))
"""
umFobs = arcpy.management.SelectLayerByAttribute(fobs, "NEW_SELECTION", umQuery)
arcpy.management.CopyFeatures(umFobs, "um_fobs")
umFobsGrd = arcpy.management.SelectLayerByLocation(
    stateGrid, "INTERSECT", "um_fobs", selection_type="NEW_SELECTION"
)
arcpy.management.CopyFeatures(umFobsGrd, f"{prjGDB}/um_fobs_grd")


umQuery = """
MgO_pct >= 18
"""
umChem = arcpy.management.SelectLayerByAttribute(
    "WACHEM_filter", "NEW_SELECTION", umQuery
)
arcpy.management.CopyFeatures(umChem, f"{prjGDB}/um_chem")
umChemGrd = arcpy.management.SelectLayerByLocation(
    stateGrid, "INTERSECT", "um_chem", selection_type="NEW_SELECTION"
)
arcpy.management.CopyFeatures(umChemGrd, f"{prjGDB}/um_chem_grd")

# Cumulates
cummQuery = """
(LOWER (OBSNOTES) LIKE LOWER ('%cumulate%') OR
LOWER (LITH_DESC) LIKE LOWER ('%cumulate%')) AND NOT
(LOWER(OBSNOTES) LIKE '%accumulated%' OR
LOWER(OBSNOTES) LIKE '%accumulates%')
"""
cummFobs = arcpy.management.SelectLayerByAttribute(fobs, "NEW_SELECTION", cummQuery)
arcpy.management.CopyFeatures(cummFobs, f"{prjGDB}/cumm_fobs")
cummFobsGrd = arcpy.management.SelectLayerByLocation(
    stateGrid, "INTERSECT", "cumm_fobs", selection_type="NEW_SELECTION"
)
arcpy.management.CopyFeatures(cummFobsGrd, f"{prjGDB}/cumm_fobs_grd")

cummQuery = """
MgO_pct >= 32
"""
cummChem = arcpy.management.SelectLayerByAttribute(
    "WACHEM_filter", "NEW_SELECTION", cummQuery
)
arcpy.management.CopyFeatures(cummChem, f"{prjGDB}/cumm_chem")
cummChemGrd = arcpy.management.SelectLayerByLocation(
    stateGrid, "INTERSECT", "cumm_chem", selection_type="NEW_SELECTION"
)
arcpy.management.CopyFeatures(cummChemGrd, f"{prjGDB}/cumm_chem_grd")

tiCummQuery = """
SAMPLETYPE = 'OUTC' And 
((Al2O3_pct / TiO2_pct >= 0) 
And (Al2O3_pct / TiO2_pct <= 50))
And MgNum >= 65
"""
tiCummChem = arcpy.management.SelectLayerByAttribute(
    "WACHEM_filter", "NEW_SELECTION", tiCummQuery
)
arcpy.management.CopyFeatures(tiCummChem, f"{prjGDB}/Ticumm_chem")
tiCummChemGrd = arcpy.management.SelectLayerByLocation(
    stateGrid, "INTERSECT", "Ticumm_chem", selection_type="NEW_SELECTION"
)
arcpy.management.CopyFeatures(tiCummChemGrd, f"{prjGDB}/Ticumm_chem_grd")


# Bimodal volcanism
thltcQuery = """
LOWER (DESCRIPTN) LIKE LOWER ('%tholeiit%') OR
LOWER (ROCKTYPE1) LIKE LOWER ('%tholeiit%') OR
LOWER (LITHNAME1) LIKE LOWER ('%tholeiit%') OR
LOWER (ROCKTYPE2) LIKE LOWER ('%tholeiit%') OR 
LOWER (LITHNAME2) LIKE LOWER ('%tholeiit%')
"""
thltcPoly = arcpy.management.SelectLayerByAttribute(
    geologyMerged, "NEW_SELECTION", thltcQuery
)

thltcQuery = """
LOWER (OBSNOTES) LIKE LOWER ('%tholeiit%') OR
LOWER (LITH_DESC) LIKE LOWER ('%tholeiit%') OR
LOWER (ROCKNAME) LIKE LOWER ('%tholeiit%') OR
LOWER (ROCKTYPE) LIKE LOWER ('%tholeiit%')
"""
thltcFobs = arcpy.management.SelectLayerByAttribute(fobs, "NEW_SELECTION", thltcQuery)

thltcPolyGrd = arcpy.management.SelectLayerByLocation(
    stateGrid, "INTERSECT", thltcPoly, selection_type="NEW_SELECTION"
)
thltcFobsGrd = arcpy.management.SelectLayerByLocation(
    thltcPolyGrd, "INTERSECT", thltcFobs, selection_type="ADD_TO_SELECTION"
)
thltAssc = arcpy.management.SelectLayerByLocation(
    thltcFobsGrd,
    "WITHIN_A_DISTANCE",
    "um_poly",
    selection_type="SUBSET_SELECTION",
    search_distance="1 Kilometers",
)
arcpy.management.CopyFeatures(thltAssc, f"{prjGDB}/thltc_assc_grd")

flscQuery = """
LOWER (DESCRIPTN) LIKE LOWER ('%felsic volcanic%') OR
LOWER (ROCKTYPE1) LIKE LOWER ('%felsic volcanic%') OR
LOWER (LITHNAME1) LIKE LOWER ('%felsic volcanic%') OR
LOWER (ROCKTYPE2) LIKE LOWER ('%felsic volcanic%') OR 
LOWER (LITHNAME2) LIKE LOWER ('%felsic volcanic%')
"""
flscPoly = arcpy.management.SelectLayerByAttribute(
    geologyMerged, "NEW_SELECTION", flscQuery
)

flscQuery = """
LOWER (OBSNOTES) LIKE LOWER ('%felsic volcanic%') OR
LOWER (LITH_DESC) LIKE LOWER ('%felsic volcanic%') OR
LOWER (ROCKNAME) LIKE LOWER ('%felsic volcanic%') OR
LOWER (ROCKTYPE) LIKE LOWER ('%felsic volcanic%')
"""
flscFobs = arcpy.management.SelectLayerByAttribute(fobs, "NEW_SELECTION", flscQuery)

flscPolyGrd = arcpy.management.SelectLayerByLocation(
    stateGrid, "INTERSECT", flscPoly, selection_type="NEW_SELECTION"
)
flscFobsGrd = arcpy.management.SelectLayerByLocation(
    flscPolyGrd, "INTERSECT", flscFobs, selection_type="ADD_TO_SELECTION"
)
flscAssc = arcpy.management.SelectLayerByLocation(
    flscFobsGrd,
    "WITHIN_A_DISTANCE",
    "um_poly",
    selection_type="SUBSET_SELECTION",
    search_distance="1 Kilometers",
)
arcpy.management.CopyFeatures(flscAssc, f"{prjGDB}/flsc_assc_grd")

# Selecting polygons from geology_MERGED thought to represent Large Igneous Provinces
# Specific codes cover the 2.8 and 2.7 Ga layered intrusions of the Meeline, Boodanoo, Warriedar and Yalgowra Suite LIP
lipQuery = """
WATECTUNIT LIKE LOWER ('%Large igneous province%') OR 
SSUITE_TEC LIKE LOWER ('%Large igneous province%') OR 
SUITE_TEC LIKE LOWER ('%Large igneous province%') OR
CODE LIKE '%A-KGak-uk' OR
CODE LIKE '%A-WAa-uk' OR
CODE LIKE '%A-_li-xbk-uk' OR
CODE LIKE '%A-KEc-uk' OR
CODE LIKE '%A-KGC-xbb-uk' OR
CODE LIKE '%A-_dk-uk' OR
CODE LIKE '%A-uk-YEG' OR
CODE LIKE '%A-KEe-uk' OR
CODE LIKE '%A-KGhh-uk' OR
CODE LIKE '%A-KGhh-uko' OR
CODE LIKE '%A-KGH-xbb-uk' OR
CODE LIKE '%A-KGhl-uk' OR
CODE LIKE '%A-KGhw-uk' OR
CODE LIKE '%A-SCh-xcx-uk' OR
CODE LIKE '%A-KG-uk' OR
CODE LIKE '%A-KG-uko' OR
CODE LIKE '%A-KGkm-uk' OR
CODE LIKE '%A-SSk-uk' OR
CODE LIKE '%A-BKl-xbk-uk' OR
CODE LIKE '%A-WAm-uk' OR
CODE LIKE '%A-KGmc-uk' OR
CODE LIKE '%A-KGmk-uk' OR
CODE LIKE '%A-_rt-uk' OR
CODE LIKE '%A-KGws-uk' OR
CODE LIKE '%A-KGks-xup-uk' OR
CODE LIKE '%A-_sf-uk' OR
CODE LIKE '%A-_sf-xbb-uk' OR
CODE LIKE '%A-KGkt-uk' OR
CODE LIKE '%A-KGwo-uk' OR
CODE LIKE '%A-uk-YYO' OR
CODE LIKE '%A-xcx-uk-YYO' OR
CODE LIKE '%A-xmbs-uk-YYO' OR
CODE LIKE '%A-xog-uk-YYO' OR
CODE LIKE '%A-KGak-xu-bk' OR
CODE LIKE '%A-KGak-uk' OR
CODE LIKE '%A-KGak-bk' OR
CODE LIKE '%A-_dk-uk' OR
CODE LIKE '%A-_dk-must' OR
CODE LIKE '%A-_dk-xuk-must' OR
CODE LIKE '%A-KGkm-uk' OR
CODE LIKE '%A-KGkm-xshh-cc' OR
CODE LIKE '%A-KGkm-mus' OR
CODE LIKE '%A-KGmc-xu-bk' OR
CODE LIKE '%A-KGmc-bk' OR
CODE LIKE '%A-KGmc-must' OR
CODE LIKE '%A-KGmc-uk' OR
CODE LIKE '%A-KGmc-up' OR
CODE LIKE '%A-_rt-xuk-bk' OR
CODE LIKE '%A-_rt-bk' OR
CODE LIKE '%A-_rt-uk' OR
CODE LIKE '%A-_rt-mus' OR
CODE LIKE '%A-_rt-xuk-mus' OR
CODE LIKE '%A-KGws-bk' OR
CODE LIKE '%A-KGws-musr' OR
CODE LIKE '%A-KGws-uk' OR
CODE LIKE '%A-KGws-mus' OR
CODE LIKE '%A-ANa%' OR
CODE LIKE '%A-ANb%' OR
CODE LIKE '%A-ANu%' OR
CODE LIKE '%A-ANl%' OR
CODE LIKE '%A-ANw%' OR
CODE LIKE '%A-ANy%' OR
CODE LIKE '%A-ANg%' OR
CODE LIKE '%A-ANo%' OR
CODE LIKE '%A-ANn%' OR
CODE LIKE '%A-SDw%' OR
CODE LIKE '%A-SDo%' OR
CODE LIKE '%A-SDf%' OR
CODE LIKE '%A-SDk%' OR
CODE LIKE '%A-ANf%' OR
CODE LIKE '%A-SDY%' OR
CODE LIKE '%FM-xmog%'
"""
lipPoly = arcpy.management.SelectLayerByAttribute(
    geologyMerged, "NEW_SELECTION", lipQuery
)
arcpy.management.CopyFeatures(lipPoly, f"{prjGDB}/lip")

# Adding in mapped dykes
dykeBuffer = arcpy.analysis.MultipleRingBuffer(
    dykes, "dykes_poly", [50], "Meters", "Dist", "ALL", "FULL", "GEODESIC"
)
arcpy.management.Merge(["dykes_poly", "lip"], f"{prjGDB}/lip_poly")

# Field observations of mantle derived rock
mantleRckQuery = """
(LOWER (OBSNOTES) LIKE LOWER ('%carbonatite%') Or
LOWER (OBSNOTES) LIKE LOWER ('%pyroxenite%') Or
LOWER (OBSNOTES) LIKE LOWER ('%kimberlite%') Or
LOWER (OBSNOTES) LIKE LOWER ('%lamphrophyre%') Or
LOWER (OBSNOTES) LIKE LOWER ('%lamphrophyr%') Or
LOWER (OBSNOTES) LIKE LOWER ('%sanukitoid%') Or
LOWER (OBSNOTES) LIKE LOWER ('%peridotite%') Or
LOWER (OBSNOTES) LIKE LOWER ('%dunite%') Or
LOWER (OBSNOTES) LIKE LOWER ('%lherzolite%') Or
LOWER (OBSNOTES) LIKE LOWER ('%harzburgite%')) Or
(LOWER (LITH_DESC) LIKE LOWER ('%carbonatite%') Or
LOWER (LITH_DESC) LIKE LOWER ('%pyroxenite%') Or
LOWER (LITH_DESC) LIKE LOWER ('%kimberlite%') Or
LOWER (LITH_DESC) LIKE LOWER ('%lamphrophyre%') Or
LOWER (LITH_DESC) LIKE LOWER ('%lamphrophyr%') Or
LOWER (LITH_DESC) LIKE LOWER ('%sanukitoid%') Or
LOWER (LITH_DESC) LIKE LOWER ('%peridotite%') Or
LOWER (LITH_DESC) LIKE LOWER ('%dunite%') Or
LOWER (LITH_DESC) LIKE LOWER ('%lherzolite%') Or
LOWER (LITH_DESC) LIKE LOWER ('%harzburgite%'))
"""
mantleRckFilter = arcpy.management.SelectLayerByAttribute(
    fobs, "NEW_SELECTION", mantleRckQuery
)
arcpy.management.CopyFeatures(mantleRckFilter, f"{prjGDB}/FOBS_mantle_rock")
mantleRckGrd = arcpy.management.SelectLayerByLocation(
    stateGrid, "INTERSECT", "FOBS_mantle_rock", selection_type="NEW_SELECTION"
)
arcpy.management.CopyFeatures(mantleRckGrd, f"{prjGDB}/FOBS_mantle_rock_grd")

# The following code filters 500k line structure for lines recorded as faults or shears
# A point layer is then generated representing fault intersections
fltsShrsQuery = """
FEATURE LIKE '%Fault or shear zone%'
"""
fltsShrsFilter = arcpy.management.SelectLayerByAttribute(
    lineStrct, "NEW_SELECTION", fltsShrsQuery
)
arcpy.management.CopyFeatures(fltsShrsFilter, f"{prjGDB}/faults_shears")
arcpy.analysis.Intersect(
    "faults_shears", f"{prjGDB}/faults_intersection", output_type="POINT"
)

# The following code filters for mafic/ultramafic lithologies

mum = """
ROCKTYPE1 LIKE '%mafic%' Or ROCKTYPE1 LIKE '%ultramafic%'
"""
mumPoly = arcpy.management.SelectLayerByAttribute(geologyMerged, "NEW_SELECTION", mum)
arcpy.management.CopyFeatures(mumPoly, f"{prjGDB}/mum_poly")

# The following code filters for mafic/ultramafic with 70% more or 70% less Ni than predicted by Mg number
# Ni prediction is based on an exponetial line of best fit between Mg # and Ni ppm
mgNum = """
MgNum > 0 AND Ni_ppm > 0
"""
mgPts = arcpy.management.SelectLayerByAttribute("WACHEM_filter", "NEW_SELECTION", mgNum)
arcpy.management.CopyFeatures(mgPts, f"{prjGDB}/Ni_Mgnum_pts")

mg65 = """
MgNum >= 65
"""
niVal = arcpy.management.SelectLayerByAttribute("Ni_Mgnum_pts", "NEW_SELECTION", mg65)
ni_pred_alpha = 0.223
ni_pred_beta = 0.1027
niPred = f"ni_pred_alpha * math.exp(ni_pred_beta * !MgNum!)"
arcpy.management.CalculateField(niVal, "NiPred", niPred, "PYTHON3", field_type="DOUBLE")

mg65 = """
(MgNum < 65 AND MgNum > 0)
"""
niVal = arcpy.management.SelectLayerByAttribute("Ni_Mgnum_pts", "NEW_SELECTION", mg65)
ni_pred_alpha = 15.44
ni_pred_beta = 0.0362
niPred = f"ni_pred_alpha * math.exp(ni_pred_beta * !MgNum!)"
arcpy.management.CalculateField(niVal, "NiPred", niPred, "PYTHON3", field_type="DOUBLE")

diff = "((!Ni_ppm! - !NiPred!)/!NiPred!)*100"
arcpy.management.SelectLayerByAttribute("Ni_Mgnum_pts", "CLEAR_SELECTION")
arcpy.management.CalculateField(
    "Ni_Mgnum_pts", "NiDiff", diff, "PYTHON3", field_type="DOUBLE"
)

niDpltd = """NiDiff <=-70"""
arcpy.management.SelectLayerByAttribute("Ni_Mgnum_pts", "CLEAR_SELECTION")
niDpltdPts = arcpy.management.SelectLayerByAttribute(
    "Ni_Mgnum_pts", "NEW_SELECTION", niDpltd
)
arcpy.management.CopyFeatures(niDpltdPts, f"{prjGDB}/Ni_Mgnum_depleted_pts")

niEnrchd = """NiDiff >=70"""
arcpy.management.SelectLayerByAttribute("Ni_Mgnum_pts", "CLEAR_SELECTION")
niEnrchdPts = arcpy.management.SelectLayerByAttribute(
    "Ni_Mgnum_pts", "NEW_SELECTION", niEnrchd
)
arcpy.management.CopyFeatures(niEnrchdPts, "Ni_Mgnum_enriched_pts")

niDpltdPoly = arcpy.management.SelectLayerByLocation(
    "mum_poly", "INTERSECT", "Ni_Mgnum_depleted_pts", selection_type="NEW_SELECTION"
)
arcpy.management.CopyFeatures(niDpltdPoly, "Ni_depleted_poly")

niEnrchdPoly = arcpy.management.SelectLayerByLocation(
    "mum_poly", "INTERSECT", "Ni_Mgnum_enriched_pts", selection_type="NEW_SELECTION"
)
arcpy.management.CopyFeatures(niEnrchdPoly, f"{prjGDB}/Ni_enriched_poly")


# The follwoing code filters for Mafic/ultramafic samples with Pt and Pd > 17ppm and < 4ppm (top and bottom 10th percentile of Pt and Pd).
Pt = '("Pt_ppb" > 2)'
positivePt = arcpy.management.SelectLayerByAttribute(
    "WACHEM_filter", "NEW_SELECTION", Pt
)
arcpy.conversion.TableToExcel(positivePt, f"{root_dir}/PtValues.xlsx")

Pd = '("Pd_ppb" > 2)'
positivePd = arcpy.management.SelectLayerByAttribute(
    "WACHEM_filter", "NEW_SELECTION", Pd
)
arcpy.conversion.TableToExcel(positivePd, f"{root_dir}/PdValues.xlsx")

pge90 = '("Pt_ppb" > 17 OR "Pd_ppb" > 17)'

pgeValues = arcpy.management.SelectLayerByAttribute(
    "WACHEM_filter", "NEW_SELECTION", pge90
)
arcpy.management.CopyFeatures(pgeValues, f"{prjGDB}/PGE_enriched_pts")

pgeEnrchd = arcpy.management.SelectLayerByLocation(
    "mum_poly", "INTERSECT", "PGE_enriched_pts", selection_type="NEW_SELECTION"
)
arcpy.management.CopyFeatures(pgeEnrchd, "PGE_enriched_poly")

pge10 = '(("Pt_ppb" > 2 AND "Pt_ppb" < 4) AND ("Pd_ppb" > 2 AND "Pd_ppb" < 4))'

pgeValues = arcpy.management.SelectLayerByAttribute(
    "WACHEM_filter", "NEW_SELECTION", pge10
)
arcpy.management.CopyFeatures(pgeValues, "PGE_depleted_pts")

pgeDpltd = arcpy.management.SelectLayerByLocation(
    "mum_poly", "INTERSECT", "PGE_depleted_pts", selection_type="NEW_SELECTION"
)
arcpy.management.CopyFeatures(pgeDpltd, f"{prjGDB}/PGE_depleted_poly")


# The follwoing code filters for Mafic/ultramafic samples with Cu < 14ppm and (bottom 10th percentile of Cu).
Cu = '("Cu_ppm" > 0)'

positiveCu = arcpy.management.SelectLayerByAttribute(
    "WACHEM_filter", "NEW_SELECTION", Cu
)
arcpy.conversion.TableToExcel(positiveCu, f"{root_dir}/CuValues.xlsx")


Cu10 = '("Cu_ppm" > 0 AND "Cu_ppm" < 14)'

cuValues = arcpy.management.SelectLayerByAttribute(
    "WACHEM_filter", "NEW_SELECTION", Cu10
)
arcpy.management.CopyFeatures(cuValues, f"{prjGDB}/Cu_depleted_pts")

cuDpltd = arcpy.management.SelectLayerByLocation(
    "mum_poly", "INTERSECT", "Cu_depleted_pts", selection_type="NEW_SELECTION"
)
arcpy.management.CopyFeatures(cuDpltd, f"{prjGDB}/Cu_depleted_poly")


# The follwoing code filters for Mafic/ultramafic samples with S > 0.2 wt% and (top 10th percentile of S).
Spct = '("S_pct" > 0)'

positiveS = arcpy.management.SelectLayerByAttribute(
    "WACHEM_filter", "NEW_SELECTION", Spct
)
arcpy.conversion.TableToExcel(positiveS, f"{root_dir}/SValues.xlsx")

S90 = '("S_pct" > 0.2)'

SValues = arcpy.management.SelectLayerByAttribute("WACHEM_filter", "NEW_SELECTION", S90)
arcpy.management.CopyFeatures(cuValues, f"{prjGDB}/S_enriched_pts")

sEnrchd = arcpy.management.SelectLayerByLocation(
    "mum_poly", "INTERSECT", "S_enriched_pts", selection_type="NEW_SELECTION"
)
arcpy.management.CopyFeatures(sEnrchd, f"{prjGDB}/S_enriched_poly")


# The follwoing code filters for Mafic/ultramafic samples with ratio Cu/Zr > 3.56 (exceeding Primative Mantle ratio).

cuZr = '("Cu_ppm" > 0 AND "Zr_ppm" > 0)'

positiveCuZr = arcpy.management.SelectLayerByAttribute(
    "WACHEM_filter", "NEW_SELECTION", cuZr
)
arcpy.conversion.TableToExcel(positiveCuZr, f"{root_dir}/CuZrValues.xlsx")


cuZr90 = '("Cu_ppm" > 0 AND "Zr_ppm" > 0 AND "Cu_ppm" / "Zr_ppm" > 3.56)'

cuZrValues = arcpy.management.SelectLayerByAttribute(
    "WACHEM_filter", "NEW_SELECTION", cuZr90
)
arcpy.management.CopyFeatures(cuValues, f"{prjGDB}/CuZr_enriched_pts")

cuZrEnrchd = arcpy.management.SelectLayerByLocation(
    "mum_poly",
    "INTERSECT",
    f"{prjGDB}/CuZr_enriched_pts",
    selection_type="NEW_SELECTION",
)
arcpy.management.CopyFeatures(cuZrEnrchd, "CuZr_enriched_poly")


# The follwoing code filters for Mafic/ultramafic samples with ratio Th/Yb > 0.17 (exceeding Primative Mantle ratio).
crstCntmntnQry = '("Th_ppm" / "Yb_ppm" > 0.17)'

crstCntmntn = arcpy.management.SelectLayerByAttribute(
    "WACHEM_filter", "NEW_SELECTION", crstCntmntnQry
)
arcpy.management.CopyFeatures(crstCntmntn, "CrustalContamination_pts")

crstCntmntn = arcpy.management.SelectLayerByLocation(
    "mum_poly", "INTERSECT", "CrustalContamination_pts", selection_type="NEW_SELECTION"
)
arcpy.management.CopyFeatures(crstCntmntn, "CrustalContamination_poly")

# The following code creates random points across map area
# Points used to test statistical significance of spatial correlations
# Set up a blank .gdb to hold random point layers
gdbName = "monetcarlo.gdb"
arcpy.management.CreateFileGDB(root_dir, gdbName)
arcpy.env.workspace = f"{root_dir}/monetcarlo.gdb"

# Create 500 layers of 2607 random points across the state
extent = ""  # Set extent of work area here
for i in range(500):
    arcpy.management.CreateRandomPoints(
        out_path=arcpy.env.workspace,
        out_name=f"testPoint_{i}",
        constraining_feature_class=extent,
        number_of_points_or_field=2607,
    )


# Create pathways to input files created in the Correlation step
minLayer = f"{prjGDB}/Ni_mineralized_DHRC_pt"

nearLayer = [
    f"{prjGDB}/um_poly",
    f"{prjGDB}/um_fobs_grd",
    f"{prjGDB}/um_chem_grd",
    f"{prjGDB}/cumm_chem_grd",
    f"{prjGDB}/cumm_fobs_grd",
    f"{prjGDB}/Ticumm_chem_grd",
    f"{prjGDB}/thltc_assc_grd",
    f"{prjGDB}/flsc_assc_grd",
    f"{prjGDB}/lip_poly",
    f"{prjGDB}/FOBS_mantle_rock_grd",
    f"{prjGDB}/faults_shears",
    f"{prjGDB}/faults_intersection",
]

AnnMonteCarlo(nearLayer, minLayer)

minLayer = f"{prjGDB}/Ni_mineralized_DHRC_pt"

nearLayer = [
    f"{prjGDB}/Ni_enriched_poly",
    f"{prjGDB}/Ni_depleted_poly",
    f"{prjGDB}/PGE_enriched_poly",
    f"{prjGDB}/PGE_depleted_poly",
    f"{prjGDB}/CuZr_enriched_poly",
    f"{prjGDB}/Cu_depleted_poly",
    f"{prjGDB}/S_enriched_poly",
    f"{prjGDB}/CrustalContamination_poly",
    f"{prjGDB}/HostLithology",
]

AsscMonteCarlo(nearLayer, minLayer)

arcpy.env.workspace = prjGDB  # change workspace back to project .gdb
