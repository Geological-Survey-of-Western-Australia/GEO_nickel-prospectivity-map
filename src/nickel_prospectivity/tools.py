# Import libraries
import arcpy
import arcpy.sa
import numpy as np
import pandas as pd

# Set up folder pathways
root_dir = "C:/project_folder"  # Your project folder here
arcpy.CheckOutExtension("Spatial")


def SpatialCorrelation(stateLayer, minLayer, featureLayer, exp):
    """Calculates spatial correlation between mineralized cells and input data.
    Counts the number of mineralized cells and total cells around input data in 1km intervals.
    For each interval the number expected mineralized cell is calculated, and compared to the observed number of mineralized cells.

    Args:
        stateLayer (vector): A fishnet grid covering the area of interest
        minLayer (vector): Grid cells that represent the distribution of mineralization in the area of interest
        featureLayer (list, vector): Vector files representing geologic data of interest
        exp (int): expected number of mineralized cells across area of interest (exp = mineralized/total cells)

    Returns:
        excel: For each file in featureLayer there is a row of count of total cells per interval, cumulative count of total cells per interval,
        count of mineralized cells per interval, cumulative count of mineralized cells per interval, a calculation of (Observed - Expected)/Expected
        mineralized cells.
    """
    np.seterr(divide="ignore", invalid="ignore")
    count = []
    labels = []

    for layer in featureLayer:
        calNN = arcpy.analysis.GenerateNearTable(
            in_features=stateLayer, near_features=layer, out_table="CalNN"
        )
        calNNArray = arcpy.da.TableToNumPyArray(calNN, "NEAR_DIST").astype(float)
        calHist = np.histogram(calNNArray, bins=2000, range=(0, 20))
        calCount, calBins = calHist
        count.append(calCount)
        cmlt = np.cumsum(calCount)
        cfd = cmlt / np.max(cmlt)
        count.append(cfd)
        minNN = arcpy.analysis.GenerateNearTable(
            in_features=minLayer, near_features=layer, out_table="InputNN"
        )
        minNNArray = arcpy.da.TableToNumPyArray(minNN, "NEAR_DIST").astype(float)
        minHist = np.histogram(minNNArray, bins=2000, range=(0, 20))
        minCount, minBins = minHist
        count.append(minCount)
        cmlt = np.cumsum(minCount)
        cfd = cmlt / np.max(cmlt)
        count.append(cfd)
        expected = calCount * exp
        ove = (minCount - expected) / expected
        count.append(ove)

    for name in featureLayer:
        labels.append(name + " Tot count")
        labels.append(name + " Total cmlt")
        labels.append(name + " Min count")
        labels.append(name + " Min cmlt")
        labels.append(name + " (O - E)/E")

    result = pd.DataFrame(count).transpose()
    result.columns = labels
    dist = pd.DataFrame(calBins).head(-1)
    result["Dist"] = dist
    result = result.set_index("Dist")
    result.to_excel(f"{root_dir}/Spatial_Stats.xlsx")


def AssociationCorrelation(stateLayer, minLayer, featureLayer):
    """Counts the number of mineralized cells and total cells within input data.

    Args:
        stateLayer (vector): A fishnet grid covering the area of interest
        minLayer (vector): Grid cells that represent the distribution of mineralization in the area of interest
        featureLayer (list, vector): Vector files representing geologic data of interest

    Returns:
        excel: For each file in featureLayer a count of total cells and mineralized cells which fall within input data.
    """
    np.seterr(divide="ignore", invalid="ignore")
    statecount = []
    mincount = []
    labels = []

    for layer_names in featureLayer:
        minGrid = arcpy.management.SelectLayerByLocation(
            minLayer, "INTERSECT", layer_names, selection_type="NEW_SELECTION"
        )
        stateGrid = arcpy.management.SelectLayerByLocation(
            stateLayer, "INTERSECT", layer_names, selection_type="NEW_SELECTION"
        )
        minCount = arcpy.management.GetCount(minGrid)
        stateCount = arcpy.management.GetCount(stateGrid)
        statecount.append(stateCount)
        mincount.append(minCount)
        labels.append(layer_names)

        Result = pd.DataFrame([labels, statecount, mincount]).transpose()
        Result.columns = ["Greenstone", "Total cells", "Min cells"]
        Result.to_excel(f"{root_dir}/Association_Stats.xlsx")


def AnnMonteCarlo(nearLayer, minLayer):
    """
    Calculates the average nearest neighbour (ANN) distance between input data and the random points in each input layer.
    Then the ANN distance between input data and mineralized points is calculated.

    Args:
        nearLayer (vector data, list): List of input geologic data.
        minLayer (points): Points representing mineralization.

    Returns:
        excel: Each column is for each dataset in nearLayer, each row is the ANN from the points in the random layers. The last row is the ANN from the mineralization points.

    """
    points = arcpy.ListFeatureClasses()
    monteCarlo = pd.DataFrame()
    nn = []

    for layer in nearLayer:
        ann = []
        layerNN = arcpy.analysis.GenerateNearTable(
            in_features=minLayer, near_features=layer, out_table="layerNN"
        )
        layerArray = arcpy.da.TableToNumPyArray(layerNN, "NEAR_DIST").astype(float)
        avg = np.average(layerArray)
        nn.append(avg)

        for i in points:
            monte = arcpy.analysis.GenerateNearTable(
                in_features=i, near_features=layer, out_table="monteNN"
            )
            monteArray = arcpy.da.TableToNumPyArray(monte, "NEAR_DIST").astype(float)
            avg = np.average(monteArray)
            ann.append(avg)

        monteCarlo[layer + "_ann"] = ann
    monteCarlo.loc[len(monteCarlo)] = nn
    monteCarlo.to_excel(f"{root_dir}/ANN_MonteCarlo.xlsx")


def AsscMonteCarlo(nearLayer, minLayer):
    """
    Calculates the average count of random points falling within each input layer. Then the count of mineralized points within each input layer is calculated.

    Args:
        nearLayer (vector data, list): List of input geologic data.
        minLayer (points): Points representing mineralization.

    Returns:
        Excel file: Each column is for each dataset in nearLayer, each row is the count of points in the random layers. The last row is the count of mineralized points.

    """
    points = arcpy.ListFeatureClasses()
    monteCarlo = pd.DataFrame()
    count = []

    for layer in nearLayer:
        avg = []
        layerSelect = arcpy.management.SelectLayerByLocation(
            in_layer=minLayer,
            overlap_type="INTERSECT",
            select_features=layer,
            selection_type="NEW_SELECTION",
        )
        layerCount = arcpy.management.GetCount(layerSelect)
        count.append(layerCount)

        for i in points:
            monteSelect = arcpy.management.SelectLayerByLocation(
                in_layer=i,
                overlap_type="INTERSECT",
                select_features=layer,
                selection_type="NEW_SELECTION",
            )
            monteCount = arcpy.management.GetCount(monteSelect)
            avg.append(monteCount)

        monteCarlo[layer + "_count"] = avg
    monteCarlo.loc[len(monteCarlo)] = count
    monteCarlo.to_excel(f"{root_dir}/Assc_MonteCarlo.xlsx")


def BuffMap(layer, buffers, score):
    """
    Creates multiple ring buffers around the input layers. Polygons are then rasterised and reclassified based on their realtive prospectivity score
    Areas outside the buffers are given a score of 0

    Args:
        layer (list): The input layers to have buffered rings.
        buffers (list, int): List of int representing the distances of buffers, in decimal degrees (0.01 dd is ~1km)
        score (list, int): List of int representing the realtive prospectivity score of each buffer. Each buffer should have an associated score.

    Returns:
          polygon: polygon of ring buffers (_bff suffix)
          raster: raster of ring buffers (rst_bff suffix)
          raster: raster of ring buffers reclassified with realtive prospectivity scores (rst_bff suffix)
    """
    remap = [list(i) for i in zip(buffers, score)]
    remap.append(["NODATA", 0])
    rescore = arcpy.sa.RemapValue(remap)
    bff = arcpy.analysis.MultipleRingBuffer(
        layer,
        layer + "_bff",
        buffers,
        "DecimalDegrees",
        "Dist",
        "ALL",
        "FULL",
        "GEODESIC",
    )
    rst = arcpy.conversion.PolygonToRaster(
        bff, "dist", layer + "_bff_rst", "CELL_CENTER"
    )
    rcl = arcpy.sa.Reclassify(rst, "Value", rescore, "NODATA")
    rcl.save(layer + "_bff_rst_rcl")


def AsscRasters(layers):
    """
    Creates rasters from input polygons. Cells in the presence of a polygons are given a prospectivity score of 1
    Cells outside the polygons are given a score of 0

    Args:
        layer (list): The input layers to have buffered rings.

    Returns:
          raster: raster of polygons (_rst suffix)
          raster: raster of polygons reclassified with a prospectivity score of 1 (_rst_rcl suffix)
    """
    remap = arcpy.sa.RemapValue([[1, 1], ["NODATA", 0]])
    for layer in layers:
        arcpy.management.AddField(layer, "RasterCode")
        arcpy.management.CalculateField(layer, "RasterCode", "1")
        rst = arcpy.conversion.PolygonToRaster(
            layer, "RasterCode", layer + "_rst", "MAXIMUM_COMBINED_AREA"
        )
        rcl = arcpy.sa.Reclassify(rst, "Value", remap, "NODATA")
        rcl.save(layer + "_rst_rcl")


def RasterProcessing(layer, weight):
    """
    Produces an intermediate predictor map. A linear fuzzy membership function is applied to the precitor maps, normalising the scores to range between 0 and 1.
    Conditional replacement is used to remove all 0 values with an arbitraily small value (0.001)
    Predictor maps are then weighted based on their importance to Ni mineralization, and confidence in the input data.

    Args:
        layer (list): The input predictor maps.
        weight (int): A weight between 0 and 1 based on the imprtance and confidence of the input data

    Returns:
        raster: Fuzzy membership raster normalized between 0 and 1 (_fzz suffix)
        raster: Raster with all 0 values replaed with 0.001 (_fzz_con suffix)
        raster: Raster multipled by weight (_fzz_con_wt suffix)
    """
    layerName = arcpy.Describe(layer).name
    fuzzMem = arcpy.sa.FuzzyMembership(layer, arcpy.sa.FuzzyLinear())
    fuzzMem.save(layerName + "_fzz")
    con = arcpy.sa.Con(fuzzMem, 0.001, fuzzMem, "VALUE=0")
    con.save(layerName + "_fzz_con")
    wt = arcpy.sa.Raster(con) * weight
    wt.save(layerName + "_fzz_con_wt")


def RasterCorrelation(stateLayer, minLayer, rasterLayer, rasterMin, rasterMax):
    """
    Samples the raster value intersected by a point layer. Raster values are divided into 100 bins, and the number of points in each bin is collated.

    Args:
        stateLayer (points): A statewide grid of points
        minLayer (points): A point layer representing Ni mines and deposits.
        rasterLayer (raster): Ni prospectivity raster
        rasterMin (int): Minimum value of the raster.
        rasterMax (int): Maximum value of the raster.

    Returns:
        excel sheet: Binned counts for points in stateLayer and minLayer
    """
    count = []
    stateValues = arcpy.sa.ExtractValuesToPoints(
        stateLayer, rasterLayer, rasterLayer + "_samples"
    )
    stateArray = arcpy.da.TableToNumPyArray(stateValues, "RASTERVALU").astype(float)
    stateHist = np.histogram(stateArray, bins=100, range=(rasterMin, rasterMax))
    stateCount, stateBins = stateHist
    count.append(stateCount)
    minValues = arcpy.sa.ExtractValuesToPoints(
        minLayer, rasterLayer, rasterLayer + "_samples"
    )
    minArray = arcpy.da.TableToNumPyArray(minValues, "RASTERVALU").astype(float)
    minHist = np.histogram(minArray, bins=100, range=(rasterMin, rasterMax))
    minCount, minBins = minHist
    count.append(minCount)
    result = pd.DataFrame(count).transpose()
    result.columns = ["Tot count", " Min count"]
    value = pd.DataFrame(stateBins).head(-1)
    result["Raster Value"] = value
    result = result.set_index("Raster Value")
    result.to_excel(f"{root_dir}/Raster_Correlation.xlsx")


def RasterROC(BarrenLayer, MinLayer, RasterLayer):
    """
    Samples the raster value intersected by a point layer. These points are exported to excel to create an ROC curve.

    Args:
        BarrenLayer (points): Points representing areas barren of Ni
        minLayer (points): A point layer representing Ni mines and deposits.
        rasterLayer (raster): Ni prospectivity raster

    Returns:
        excel sheet: raster values for each point in BarrenLayer and MinLayer
    """
    count = []
    barrenValues = arcpy.sa.ExtractValuesToPoints(
        BarrenLayer, RasterLayer, RasterLayer + "_samples"
    )
    barrenArray = arcpy.da.TableToNumPyArray(barrenValues, "RASTERVALU").astype(float)
    count.append(barrenArray)
    minValues = arcpy.sa.ExtractValuesToPoints(
        MinLayer, RasterLayer, RasterLayer + "_samples"
    )
    minArray = arcpy.da.TableToNumPyArray(minValues, "RASTERVALU").astype(float)
    count.append(minArray)
    result = pd.DataFrame(count).transpose()
    result.columns = ["Barren count", " Min count"]
    result.to_excel(f"{root_dir}/Raster_ROC.xlsx")


def RasterTornado(MinLayer, RasterLayer):
    """
    Samples a list of raster value intersected by a point layer. These points are exported to excel to create a tornado plot.

    Args:
        minLayer (points): A point layer representing Ni mines and deposits.
        rasterLayer (list; raster): List of rasters to be sampled

    Returns:
        excel sheet: raster values of each point for each raster in the list
    """
    count = []
    labels = []

    for raster in RasterLayer:
        minValues = arcpy.sa.ExtractValuesToPoints(
            MinLayer, raster, raster + "_samples"
        )
        minArray = arcpy.da.TableToNumPyArray(minValues, "RASTERVALU").astype(float)
        count.append(minArray)
    for name in RasterLayer:
        labels.append(name)

    result = pd.DataFrame(count).transpose()
    result.columns = labels
    result.to_excel(f"{root_dir}/Raster_tornado.xlsx")
