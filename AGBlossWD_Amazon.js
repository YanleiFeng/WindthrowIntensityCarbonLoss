
///////////////// ESA_CCI_2020 /////////////////////////////
// Organize windthrow biomass and reference biomass
// Sample Code
// Organized on May 7, 2025, by Yanlei Feng

// This script prepares data for windthrow biomass and reference biomass,
// and extracts associated environmental variables for regression analysis.
// Variables include: CAPE, DTFE, DTMR, FH, HAND, MAP, MAT, MeanLat, MeanLon,
// NumMCS, RefAviBiomass, etc.

// --------------------------------------------------------------------------------------
//  Updates:
//  - Mar 27: boundary to 1000m; NPV wt detection to Hansen disturbance detection
//  - Mar 30: forest mask to 50%; filter wt 2015-2020; handle 2008 case
//  - Jul  9: add NPV calculation for windthrow intensity comparison
//  - Jul 10: updated wt mask combining Hansen and custom NPV; choose best imagery
//  - Jul 28: added David's datasets; moved 2014 code into 2010-2015 block
//  - Aug  7: added more variables for download
//  - Oct 16: added DtMR, DtFE, HAND, FH, soil fertility variables
// --------------------------------------------------------------------------------------

// ======================================
// Section 1: Map Setup and Biomass Layers
// ======================================

// Add Amazon boundary layer (toggle off by default)
Map.addLayer(amz_boundary, {}, "amz_boundary", false);

// Import and display ESA CCI biomass maps (2019, 2020)
var biomassMap_2020 = ee.Image("projects/wind-project-309500/assets/ESA_CCI/biomass_agb_2020");
Map.addLayer(biomassMap_2020, {min:150, max:250, palette:["yellow","green","cyan","blue"]}, "biomassMap_2020", false);

var biomassMap_2019 = ee.Image("projects/wind-project-309500/assets/ESA_CCI/biomass_agb_2019");
Map.addLayer(biomassMap_2019, {min:150, max:250, palette:["yellow","green","cyan","blue"]}, "biomassMap_2019", false);

// ======================================
// Section 2: Windthrow Polygons by Year
// ======================================

// Load windthrow polygons and create subsets by year
var polygons = ee.FeatureCollection('projects/wind-project-309500/assets/carbon_windthrow/polygon_with_years');
var polygons_et2005 = polygons.filter(ee.Filter.lte("year", 2005));
var polygons_et2010 = polygons.filter(ee.Filter.lte("year", 2010)).filter(ee.Filter.gt("year", 2005));
var polygons_et2015 = polygons.filter(ee.Filter.lte("year", 2015)).filter(ee.Filter.gt("year", 2010));
var polygons_et2020 = polygons.filter(ee.Filter.lte("year", 2020)).filter(ee.Filter.gt("year", 2015));

// Add polygon layers (all and by period)
Map.addLayer(polygons, {palette:'black'}, "polygon_with_years", false);
Map.addLayer(polygons_et2005, {color:'blue'}, "polygon_<=2005", false);
Map.addLayer(polygons_et2010, {color:'green'}, "polygon_2006-2010", false);
Map.addLayer(polygons_et2015, {color:'yellow'}, "polygon_2011-2015", false);
Map.addLayer(polygons_et2020, {color:'red'}, "polygon_2016-2020", false);

// ======================================
// Section 3: Cloud Masking Function
// ======================================

// maskL8Clouds(): Masks clouds, shadows, cirrus, and dilated clouds in Landsat 8 QA_PIXEL band
function maskL8Clouds(image) {
  var qa = image.select('QA_PIXEL');
  var cloudBitMask      = 1 << 3;  // high confidence cloud
  var cloudshadowBitMask= 1 << 4;  // high confidence cloud shadow
  var dilatedCloudMask  = 1 << 1;  // cloud dilation
  var cirrusMask        = 1 << 15; // high confidence cirrus
  var clearMask         = 1 << 6;  // clear

  // create individual masks
  var cloud      = qa.bitwiseAnd(cloudBitMask).neq(0);
  var shadow     = qa.bitwiseAnd(cloudshadowBitMask).neq(0);
  var dilated    = qa.bitwiseAnd(dilatedCloudMask).neq(0);
  var cirrus     = qa.bitwiseAnd(cirrusMask).neq(0);
  var notClear   = qa.bitwiseAnd(clearMask).eq(0);

  // combine and invert mask to keep clear pixels
  var mask = cloud.or(shadow).or(dilated).or(cirrus).or(notClear);
  return image.updateMask(mask.not());
}

// ======================================
// Section 4: Spectral Mixture Analysis
// ======================================

// smaima(): Performs spectral unmixing using predefined endmembers (NPV, GV, Shade)
var smaima = function(image) {
  var gv =   [0.1164, 0.0904, 0.0726, 0.0459, 0.3635, 0.1489, 0.0476];
  var npv=   [0.1263, 0.1059, 0.0916, 0.0992, 0.2374, 0.2782, 0.1451];
  var shd=   [0.1218, 0.0981, 0.0772, 0.0649, 0.0470, 0.0102, 0.0040];
  return image.unmix([npv, gv, shd]);
};

// selectBands(): Chooses Landsat bands 1-7 for analysis
function selectBands(image) {
  return image.select(["B1","B2","B3","B4","B5","B6","B7"]);
}

// ======================================
// Section 5: Kernel Functions for Buffers
// ======================================

// kernel fuction, to make a buffer with user-defined radius around raster pixels
// use 70m to complement that 30m resolution for Landsat 30m
var kernel_70 = function(image){
  var kernelimage = image.focal_max({radius:70, units:'meters'});
  return kernelimage;
}

// used for excluding 100m from the windthrow
// kernel fuction, to make a buffer with user-defined radius around raster pixels
var kernel_1070 = function(image){
  var kernelimage = image.focal_max({radius:1070, units:'meters'});
  return kernelimage;
}


// ======================================
// Section 6: Forest Masks (Hansen Data)
// ======================================

// Load Hansen forest change dataset for initial forest mask
var hansen = ee.Image("UMD/hansen/global_forest_change_2022_v1_10").clip(amz_boundary);
var hansen_fc2000 = hansen.select("treecover2000").gte(50);        // forest in 2000
var hansen_loss   = hansen.select("loss");
var hansen_gain   = hansen.select("gain");
// Compute forest at 2022: original forest minus loss plus gain
var hansen_fc2022 = hansen_fc2000.subtract(hansen_loss).add(hansen_gain).gte(1);
Map.addLayer(hansen_fc2022, {}, "hansen_forest_2022", false);

// Forest mask at 2015 for reference mask
var gfc2014      = ee.Image("UMD/hansen/global_forest_change_2015");
var loss2015     = gfc2014.select("loss");
var gain2015     = gfc2014.select("gain");
var treeCover2000= gfc2014.select("treecover2000");
var hansen_fc2015= treeCover2000.gte(50).subtract(loss2015).add(gain2015).gte(1);
Map.addLayer(hansen_fc2015, {}, "hansen_forest_2015", false);

// ======================================
// Section 7: Define Windthrow Mask
// ======================================

// Invert 2022 forest mask to get windthrow areas, clip to boundary
var only_wt_mask   = hansen_fc2022.not().clip(amz_boundary);
var fm_amz         = hansen_fc2015.clip(amz_boundary);
var full_wt        = kernel_70(only_wt_mask.updateMask(fm_amz));
var wt_and_surround= kernel_1070(only_wt_mask);

// Add windthrow layers (toggle off)
Map.addLayer(full_wt, {}, "full_wt", false);
Map.addLayer(wt_and_surround, {}, "wt_and_surround", false);

// ======================================
// Section 8: Reference Mask Function
// ======================================

// calculateReferenceMask(): Identifies reference pixels based on elevation and soil percentiles
function calculateReferenceMask(region, wt_mask, fm_mask) {
  var elev = ee.Image("USGS/SRTMGL1_003").select("elevation").clip(amz_boundary).clip(region);
  var soil= ee.Image("projects/soilgrids-isric/ocd_mean").select("ocd_0-5cm_mean").clip(amz_boundary).clip(region);
  // Create windthrow buffer
  var wt_buf = kernel_70(wt_mask.updateMask(fm_mask));
  // Compute 5th and 95th elevation percentiles within wt area
  var elevPct = elev.mask(wt_buf).reduceRegion({
    reducer: ee.Reducer.percentile([5,95]),
    geometry: region, scale:30, bestEffort:true
  });
  // Compute soil percentiles
  var soilPct = soil.mask(wt_buf).reduceRegion({
    reducer: ee.Reducer.percentile([5,95]),
    geometry: region, scale:30, bestEffort:true
  });
  // Define lower/upper bounds
  var elevLow = elev.gte(elevPct.get('elevation_p5'));
  var elevHigh= elev.lte(elevPct.get('elevation_p95'));
  var soilLow = soil.gte(soilPct.get('ocd_0-5cm_mean_p5'));
  var soilHigh= soil.lte(soilPct.get('ocd_0-5cm_mean_p95'));
  // Combine and exclude wt surround
  var combined= elevLow.and(elevHigh).and(soilLow).and(soilHigh);
  var surround= kernel_1070(wt_mask);
  return combined.updateMask(surround.not()).updateMask(fm_mask).clip(region);
}

// ======================================
// Section 9: Biomass Extraction Functions
// ======================================

// getBiomass: Retrieves ESA CCI biomass for masked region
function getBiomass(mask, region) {
  return biomassMap_2020.updateMask(mask).clip(region);
}

// reduceBiomass: Computes mean biomass value within region
function reduceBiomass(layer, region) {
  var meanDict = layer.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: region, scale:100
  });
  return meanDict.get('b1');
}

// ======================================
// Section 10: NPV Computation Functions
// ======================================

// NPV_wt: Calculates NPV and windthrow intensity
// Function to calculate NPV (Non-Photosynthetic Vegetation) of the windthrow area for a given region and year
function NPV_wt(region, wt_mask, year, img, cm) {
  // Select the Landsat 8 image collection for the specified year
  
  var landsatCollection = ee.Image(img)
  var cloudMaskedlayer = maskL8Clouds(landsatCollection).clip(region)
  // Apply cloud mask to the image collection
  
  // var cloudMaskedlayer = landsatCollection.map(maskL8Clouds).mean().clip(region)
  Map.addLayer(cloudMaskedlayer, {'bands': 'B6,B5,B4', 'gain': '800, 500, 800'}, "cloudMaskedlayer",false)
  
  // Select only the bands of interest
  var selectedBandslayer = selectBands(cloudMaskedlayer)
  
  // Apply SMA (Spectral Mixture Analysis) to the selected bands
  var smaresults = smaima(selectedBandslayer);
  
  // Select NPV and GV (Green Vegetation) fractions
  var npv_1 = ee.Image(smaresults.select(0));
  var gv_1  = ee.Image(smaresults.select(1));
  var nrmlzd_npv_1 =npv_1.divide(npv_1.add(gv_1));
  
  // Add the normalized NPV layer to the map (optional)
  Map.addLayer(nrmlzd_npv_1, {min: 0, max: 1, palette: ['white', 'brown']}, 'Normalized NPV',false);
  
  ////////////////////////////////////////////////////////////////////////////////
  // In addition  to Hansen's map, we build additional NPV mask to get the wt
  // We union both results
  var meanValue = nrmlzd_npv_1.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: region,
    scale: 30, // Specify the scale depending on the resolution you want
    maxPixels: 1e9 // Adjust as necessary to manage computation limits
  });
  
  var stdValue = nrmlzd_npv_1.reduceRegion({
    reducer: ee.Reducer.stdDev(),
    geometry: region,
    scale: 30, // Specify the scale depending on the resolution you want
    maxPixels: 1e9 // Adjust as necessary to manage computation limits
  });
  
  var meanNPVvalue = ee.Number(meanValue.get('band_0'))
  var stdNPVvalue = ee.Number(stdValue.get('band_0')).multiply(2)
  var NPVthreshold = meanNPVvalue.add(stdNPVvalue) //95%
  var extremeNPVmask = nrmlzd_npv_1.gt(NPVthreshold) // Get the wt mask from larger than 95% of NPV value
  
  if (cm === "Y"){
    var wtmask = extremeNPVmask.or(wt_mask); // combined
  } else if (cm === "N") {
    var wtmask = wt_mask // hansen
  }
  // var wtmask = extremeNPVmask.or(wt_mask)
  // Map.addLayer(wtmask, {}, "combined wt mask")
  // Map.addLayer(wt_mask, {}, "single hansen mask")
  //////////////////////////////////////////////////////////////////////////
  
  // Apply the mean reduction over the specified region
  var meanValue = nrmlzd_npv_1.mask(wtmask).reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: region,
    scale: 30, // Specify the scale depending on the resolution you want
    maxPixels: 1e9 // Adjust as necessary to manage computation limits
  });

  var meanNPVValue = meanValue.get('band_0');
  // Return the mean value as a dictionary
  // Return both mask and value as a dictionary
  return {
    'wtmask': wtmask,
    'meanNPVValue': meanNPVValue
  };
}

// NPV_ref: Calculates mean NPV in reference area
// Function to calculate NPV (Non-Photosynthetic Vegetation) of the windthrow area for a given region and year
function NPV_ref(region, ref_mask, year, img) {

  var landsatCollection = ee.Image(img)
  var cloudMaskedlayer = maskL8Clouds(landsatCollection).clip(region)
  // Apply cloud mask to the image collection
  
  // var cloudMaskedlayer = landsatCollection.map(maskL8Clouds).mean().clip(region)
  Map.addLayer(cloudMaskedlayer, {'bands': 'B6,B5,B4', 'gain': '700, 500, 800'}, "cloudMaskedlayer", false)
  
  // Select only the bands of interest
  var selectedBandslayer = selectBands(cloudMaskedlayer)
  
  // Apply SMA (Spectral Mixture Analysis) to the selected bands
  var smaresults = smaima(selectedBandslayer);
  
  // Select NPV and GV (Green Vegetation) fractions
  var npv_1 = ee.Image(smaresults.select(0));
  var gv_1  = ee.Image(smaresults.select(1));
  var nrmlzd_npv_1 =npv_1.divide(npv_1.add(gv_1));
  
  // Add the normalized NPV layer to the map (optional)
  Map.addLayer(nrmlzd_npv_1, {min: 0, max: 1, palette: ['white', 'brown']}, 'Normalized NPV',false);
  

  // Apply the mean reduction over the specified region
  var meanValue = nrmlzd_npv_1.mask(ref_mask).reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: region,
    scale: 30, // Specify the scale depending on the resolution you want
    maxPixels: 1e9 // Adjust as necessary to manage computation limits
  });

  var meanNPVValue = meanValue.get('band_0');
  // Return the mean value as a dictionary
  return meanNPVValue;
}

// ======================================
// Section 11: Environmental Variables
// ======================================

// Weather
var MAP = ee.ImageCollection("NASA/GPM_L3/IMERG_MONTHLY_V07").filterDate("2000-01-01","2020-01-01").select(['precipitation']).mean().clip(amz).float()
var TP = ee.Image('users/ylfeng/deepconvections/TotalPrecip_1990_2020').clip(amz).float()
var NumMCS = ee.Image('projects/wind-project-309500/assets/MCSs_num_2001_2010_every5yrs_average').clip(amz).float()
var CAPE = ee.Image('users/ylfeng/deepconvections/CAPE_1990_2019_afternoon_mean').clip(amz).float()
var WS = ee.Image('users/ylfeng/deepconvections/ERA5-vertical-wind-shear').clip(amz).float()
var WTdensity = ee.Image('users/ylfeng/Heatmap2dg_wt_density_RNJ').clip(amz).float()
var MAT = ee.ImageCollection("MODIS/061/MOD11A1").filterDate("2000-01-01","2020-01-01").select(['LST_Day_1km']).mean().clip(amz).float()
Map.addLayer(WTdensity, {}, "WTdensity", false)


//Topographic
// Elevation
var elev = ee.Image("USGS/SRTMGL1_003").clip(amz).float()
// Aspect
var aspect = ee.Terrain.aspect(elevation).clip(amz).float()
// Slope
var slope = ee.Terrain.slope(elevation).clip(amz).float()

// Forest characteristics
// Forest Type
var foresttype = ee.Image("ESA/GLOBCOVER_L4_200901_200912_V2_3").select('landcover').rename("foresttype").clip(amz).float()
// Topographic Diversity
var topodiv = ee.Image("CSP/ERGo/1_0/Global/SRTM_topoDiversity").clip(amz).float().rename("topodiversity")
// Landforms
var landform = ee.Image("CSP/ERGo/1_0/Global/SRTM_landforms").rename("landform").clip(amz).float()
// TC only have data till 2014, so I use 2014 same reason data 
var tc = ee.ImageCollection("Oxford/MAP/TCW_5km_Monthly").filterDate("2014-09-01","2014-09-30").select('Mean').mosaic().clip(amz).rename("wetness").float()

// Soil grids
var soc = ee.Image("projects/soilgrids-isric/soc_mean").select("soc_0-5cm_mean").clip(amz).float();
// Map.addLayer(soc, {}, "soc", false)
var nitrogen = ee.Image("projects/soilgrids-isric/nitrogen_mean").select("nitrogen_0-5cm_mean").clip(amz).float();

// Additional variables added
var DtMR = ee.Image("projects/wind-project-309500/assets/Amazon_dist_to_main_river").rename("DtMR").clip(amz).float()
var DtFE = ee.Image("projects/wind-project-309500/assets/Amazon_dis_forest_edge").rename("DtFE").clip(amz).float()
// Map.addLayer(DtMR, {}, "DtMR")
var hand30_100 = ee.ImageCollection("users/gena/global-hand/hand-100").mosaic().rename("HAND").clip(amz).float()
var FCH_GEDI = ee.ImageCollection("users/potapovpeter/GEDI_V27").mosaic().rename("FH").clip(amz).float()
var Soil_cation = ee.Image("projects/wind-project-309500/assets/KrigeSoil_fernRcropNEW").rename("SoilFert").clip(amz).float()

// ======================================
// Section 12: Variable Extraction Function
// ======================================

// Function to get the average values of various environmental variables over a masked and clipped region
function getVariables(mask, region) {
  // Apply mask and clip
  var maskedMAP = MAP.clip(region);
  var maskedTP = TP.clip(region);
  var maskedCAPE = CAPE.clip(region);
  var maskedWS = WS.clip(region);
  var maskedWTdensity = WTdensity.clip(region);
  var maskedMAT = MAT.clip(region);
  var maskedNumMCS = NumMCS.clip(region);
  var maskedelev = elev.updateMask(mask).clip(region);
  var maskedaspect = aspect.updateMask(mask).clip(region);
  var maskedslope = slope.updateMask(mask).clip(region);
  var maskedForestType = foresttype.updateMask(mask).clip(region);
  var maskedSOC = soc.updateMask(mask).clip(region);
  var maskedNitrogen = nitrogen.updateMask(mask).clip(region);
  var maskedDtR = DtMR.updateMask(mask).clip(region);
  var maskedDtE = DtFE.updateMask(mask).clip(region);
  var maskedHAND = hand30_100.updateMask(mask).clip(region);
  var maskedFH = FCH_GEDI.updateMask(mask).clip(region);
  var maskedSF = Soil_cation.updateMask(mask).clip(region);
  
  
  // Reduce each masked and clipped layer to get the mean
  var meanValues = {
    meanMAP: maskedMAP.reduceRegion({
      reducer: ee.Reducer.mean(),
      geometry: region,
      scale: 11133
    }).get('precipitation'),
    meanTP: maskedTP.reduceRegion({
      reducer: ee.Reducer.mean(),
      geometry: region,
      scale: 27000
    }).get('b1'),
    meanCAPE: maskedCAPE.reduceRegion({
      reducer: ee.Reducer.mean(),
      geometry: region,
      scale: 27000
    }).get('b1'),
    meanWTdensity: maskedWTdensity.reduceRegion({
      reducer: ee.Reducer.mean(),
      geometry: region,
      scale: 1100
    }).get('b1'),
    meanMAT: maskedMAT.reduceRegion({
      reducer: ee.Reducer.mean(),
      geometry: region,
      scale: 1000
    }).get('LST_Day_1km'),
     meanNumMCS: maskedNumMCS.reduceRegion({
      reducer: ee.Reducer.mean(),
      geometry: region,
      scale: 11000
    }).get('b1'),
    meanElev: maskedelev.reduceRegion({
      reducer: ee.Reducer.mean(),
      geometry: region,
      scale: 30
    }).get('elevation'),
    meanAspect: maskedaspect.reduceRegion({
      reducer: ee.Reducer.mean(),
      geometry: region,
      scale: 30
    }).get('aspect'),
    meanSlope: maskedslope.reduceRegion({
      reducer: ee.Reducer.mean(),
      geometry: region,
      scale: 30
    }).get('slope'),
    meanForestType: maskedForestType.reduceRegion({
      reducer: ee.Reducer.mean(),
      geometry: region,
      scale: 100
    }).get('foresttype'),
    meanSOC: maskedSOC.reduceRegion({
      reducer: ee.Reducer.mean(),
      geometry: region,
      scale: 100
    }).get('soc_0-5cm_mean'),
    meanNitrogen: maskedNitrogen.reduceRegion({
      reducer: ee.Reducer.mean(),
      geometry: region,
      scale: 100
    }).get('nitrogen_0-5cm_mean'),
    meanDtR: maskedDtR.reduceRegion({
      reducer: ee.Reducer.mean(),
      geometry: region,
      bestEffort: true, 
      scale: 100
    }).get('DtMR'),
    meanDtE: maskedDtE.reduceRegion({
      reducer: ee.Reducer.mean(),
      geometry: region,
      bestEffort: true, 
      scale: 100
    }).get('DtFE'),
    meanHAND: maskedHAND.reduceRegion({
      reducer: ee.Reducer.mean(),
      geometry: region,
      scale: 30
    }).get('HAND'),
    meanFH: maskedFH.reduceRegion({
      reducer: ee.Reducer.mean(),
      geometry: region,
      scale: 30
    }).get('FH'),
    meanSF: maskedSF.reduceRegion({
      reducer: ee.Reducer.mean(),
      geometry: region,
      scale: 30
    }).get('SoilFert'),
  };

  return meanValues;
}


// ======================================
// Section 13: Main Loop to Assemble Features
// ======================================

// Updated list of regions for more results
var regions = [
  {geometry: ee.Geometry.Rectangle([-62.36832, -6.86386, -62.2659, -6.7712]), year: 2017, img: 'LANDSAT/LC08/C02/T1_TOA/LC08_232065_20170721', cm: "Y"},
  //.. more cases are added for the analysis. Here we just show one windthrow case for the simplicity of sample code
];

// Define a function to attach region information to raster pixels
function assignRegion(region, regionalization_) {
  // Filter features that intersect with feature1
  var intersectingFeature = regionalization_.filterBounds(region);

  // Extract properties of the intersecting feature
  var properties = intersectingFeature.first().get("Region");

  return properties
}

var results = [];
regions.forEach(function(item) {
  var region = item.geometry;
  var year = item.year;
  var img = item.img
  var cm = item.cm
  // Get the combined windthrow mask, I calculate it using 2 std from mean NPV from mean
  // Add this wtmask with Hansen's mask, union them, updated the mask to be a combined wt mask
  // Get the NPV of windthrow
  
  var NPV_wt_results = NPV_wt(region, only_wt_mask, year, img, cm)
  var wt_com_mask = NPV_wt_results.wtmask // this is decided by last step if combined or separate
  
  var wt_NPV = NPV_wt_results.meanNPVValue
  
  // Map.addLayer(wt_com_mask, {}, "wt_com_mask")
  // Map.addLayer(only_wt_mask, {}, "only_wt_mask")
  
  // Calculate the reference mask for the current region
  var refMask = calculateReferenceMask(region, wt_com_mask, fm_amz);
 
  // Get the NPV of reference parcel
  var ref = refMask
  var ref_NPV = NPV_ref(region, ref, year, img)
  
  //////////////////////////////////////////////////////
  // Note that the windthrow mask used for biomass and NPV should be different because they have
  // different resolution. We apply 70m kernel to the wtmask to get the biomass
  // It doesn't matter that much since I am calculating the mean biomass
  
   
  var full_wt = kernel_70(wt_com_mask.updateMask(fm_amz))
  // Get the biomass for the reference area
  var refBiomass = getBiomass(refMask, region);

  Map.addLayer(refBiomass, {min: 150, max: 300, palette:["yellow", "green", "blue"]}, "refBiomass")
  // Get the biomass for the windthrow area
  var wtBiomass = getBiomass(full_wt, region);
  Map.addLayer(wtBiomass, {min: 150, max: 300, palette:["yellow", "green", "blue"]}, "wtBiomass")
  
  // Reduce biomass to get mean values
  var reducedRefBiomass = reduceBiomass(refBiomass, region);
  var reducedWtBiomass = reduceBiomass(wtBiomass, region);
  print(reducedRefBiomass, "reducedRefBiomass")
  var refBiomass_constant = ee.Image.constant(reducedRefBiomass).mask(full_wt).clip(region)
  var DiffBiomass = refBiomass_constant.subtract(wtBiomass)
  // Map.addLayer(wtBiomass, {}, "wtBiomass") // You can print for checking
  // Map.addLayer(DiffBiomass, {}, "DiffBiomass")
  ///////////////Add evaluation from Avitabile biomass
  var refAviBiomass = getAviBiomass(refMask, region);
  var reducedAviRefBiomass = reduceBiomass(refAviBiomass, region);
  
  ///////////////////////////////////////////////////////
  // Variables
  var reducedVariables = getVariables(full_wt, region)

  
 // updated the windthrow mask to update the windthrow size
  var wt_size = wt_com_mask.clip(region)
                .reduceRegion({
                  reducer: ee.Reducer.sum(),
                  geometry: region,
                  scale: 30, // Adjust this scale to the resolution of your image, in meters
                  maxPixels: 1e9 // Adjust if necessary to accommodate large areas
                });
         
  // Get the centroid of the region
  var centroid = region.centroid();
  var centroidCoords = centroid.coordinates();
  var meanLon = centroidCoords.get(0);
  var meanLat = centroidCoords.get(1);
   
  // add region from Robinson's regionalization results
  var RegionProp = assignRegion(region, regionalization);
  
  // Construct an ee.Feature with year included
  var feature = ee.Feature(null, {
    'Year': year,
    'MeanLat': meanLat,
    'MeanLon': meanLon,
    'RefBiomass': ee.Number(reducedRefBiomass),
    'WtBiomass': ee.Number(reducedWtBiomass),
    'RefAviBiomass': ee.Number(reducedAviRefBiomass),
    'WtNPV': ee.Number(wt_NPV),
    'RefNPV': ee.Number(ref_NPV),
    'wtSize':ee.Number(ee.List(wt_size.values()).get(0)).multiply(900),
    'MAP':ee.Number(reducedVariables.meanMAP),
    'TP':ee.Number(reducedVariables.meanTP),
    'CAPE':ee.Number(reducedVariables.meanCAPE),
    'WS':ee.Number(reducedVariables.meanWS),
    'WTdensity':ee.Number(reducedVariables.meanWTdensity),
    'MAT':ee.Number(reducedVariables.meanMAT),
    'NumMCS':ee.Number(reducedVariables.meanNumMCS),
    'foresttype':ee.Number(reducedVariables.meanForestType),
    'TopoDiv':ee.Number(reducedVariables.meanTopoDiversity),
    'landform':ee.Number(reducedVariables.meanLandform),
    'TCwetness':ee.Number(reducedVariables.meanTC),
    'SOC':ee.Number(reducedVariables.meanSOC),
    'SN':ee.Number(reducedVariables.meanNitrogen),
    'elevation':ee.Number(reducedVariables.meanElev),
    'aspect':ee.Number(reducedVariables.meanAspect),
    'slope':ee.Number(reducedVariables.meanSlope),
    'regionalization':ee.Number(RegionProp),
    'DTMR':ee.Number(reducedVariables.meanDtR),
    'DTFE':ee.Number(reducedVariables.meanDtE),
    'HAND':ee.Number(reducedVariables.meanHAND),
    'FH':ee.Number(reducedVariables.meanFH),
    'SF':ee.Number(reducedVariables.meanSF),
  });
  
  // Add the feature to the results list
  results.push(feature); // Now directly adding ee.Feature objects to results
});


// Create FeatureCollection and print
var featureCollection = ee.FeatureCollection(results);
print(featureCollection);

// ======================================
// Section 14: Export Results to Drive
// ======================================

Export.table.toDrive({
  collection: featureCollection,
  description: 'windthrow_biomass_NPV_size_updatedNPV_Img_updatedWtMask_2015_2019',
  folder: 'GEEexport',
  fileFormat: 'CSV'
});
