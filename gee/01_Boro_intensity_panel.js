// ============================================================
// STEP 1 — BANGLADESH STUDY DISTRICTS
// ============================================================
var admin2 = ee.FeatureCollection('FAO/GAUL/2025/level2');
var bangladesh = admin2.filter(ee.Filter.eq('GAUL0_NAME', 'Bangladesh'));

var studyDistricts = bangladesh.filter(
  ee.Filter.inList('GAUL2_NAME', ['Rajshahi', 'Rangpur', 'Mymensingh', 'Sunamganj'])
);

Map.centerObject(studyDistricts, 7);
Map.addLayer(studyDistricts.style({color: 'FF0000', fillColor: '00000000', width: 3}), {}, 'Study Districts');


// ============================================================
// STEP 2 — CREATE STUDY GEOMETRY
// ============================================================
var studyGeometry = studyDistricts.geometry();
print('Study geometry area (m2):', studyGeometry.area());


// ============================================================
// STEP 3 — MULTI-SENSOR LANDSAT COLLECTION SETUP (2001–2025)
// ============================================================
var l5 = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2');
var l7 = ee.ImageCollection('LANDSAT/LE07/C02/T1_L2');
var l8 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2');
var l9 = ee.ImageCollection('LANDSAT/LC09/C02/T1_L2');

function harmonizeL5L7(image) {
  var qa = image.select('QA_PIXEL');
  var cloudShadow = qa.bitwiseAnd(1 << 4).eq(0);
  var cloud = qa.bitwiseAnd(1 << 3).eq(0);
  
  var optical = image.select(['SR_B1', 'SR_B3', 'SR_B4', 'SR_B5']) // Blue, Red, NIR, SWIR1 for L5/7
                     .multiply(0.0000275).add(-0.2)
                     .rename(['Blue', 'Red', 'NIR', 'SWIR1']);
                     
  var thermal = image.select('ST_B6')
                     .multiply(0.00341802).add(149.0)
                     .rename('LST');

  return image.addBands(optical, null, true)
              .addBands(thermal, null, true)
              .updateMask(cloud.and(cloudShadow))
              .select(['Blue', 'Red', 'NIR', 'SWIR1', 'LST'])
              .copyProperties(image, ['system:time_start']);
}

function harmonizeL8L9(image) {
  var qa = image.select('QA_PIXEL');
  var cloudShadow = qa.bitwiseAnd(1 << 4).eq(0);
  var cloud = qa.bitwiseAnd(1 << 3).eq(0);
  
  var optical = image.select(['SR_B2', 'SR_B4', 'SR_B5', 'SR_B6']) // Blue, Red, NIR, SWIR1 for L8/9
                     .multiply(0.0000275).add(-0.2)
                     .rename(['Blue', 'Red', 'NIR', 'SWIR1']);
                     
  var thermal = image.select('ST_B10')
                     .multiply(0.00341802).add(149.0)
                     .rename('LST');

  return image.addBands(optical, null, true)
              .addBands(thermal, null, true)
              .updateMask(cloud.and(cloudShadow))
              .select(['Blue', 'Red', 'NIR', 'SWIR1', 'LST'])
              .copyProperties(image, ['system:time_start']);
}

var landsatCol = l5.map(harmonizeL5L7)
  .merge(l7.map(harmonizeL5L7))
  .merge(l8.map(harmonizeL8L9))
  .merge(l9.map(harmonizeL8L9));


// ============================================================
// STEP 4 — ADD 30m INDICES (NDVI, EVI, LSWI, NIRv)
// ============================================================
function addLandsatIndices(image) {
  var blue = image.select('Blue');
  var red = image.select('Red');
  var nir = image.select('NIR');
  var swir1 = image.select('SWIR1');

  var ndvi = nir.subtract(red).divide(nir.add(red)).rename('NDVI');
  var evi = nir.subtract(red).multiply(2.5)
              .divide(nir.add(red.multiply(6)).subtract(blue.multiply(7.5)).add(1))
              .rename('EVI');
  var lswi = nir.subtract(swir1).divide(nir.add(swir1)).rename('LSWI');
  var nirv = nir.multiply(ndvi).rename('NIRv');

  return image.addBands([ndvi, evi, lswi, nirv]);
}


// ============================================================
// STEP 5 — MONTHLY LANDSAT COMPOSITING FUNCTION
// ============================================================
function monthlyLandsatComposite(startDate, endDate) {
  var collection = landsatCol
    .filterBounds(studyGeometry)
    .filterDate(startDate, endDate)
    .map(addLandsatIndices);

  var composite = collection.median().select(['NDVI', 'EVI', 'LSWI', 'NIRv', 'LST']);
  var monthNumber = ee.Number(startDate.get('month'));
  var monthBand = composite.select('NDVI').multiply(0).add(monthNumber).rename('month');

  return composite.addBands(monthBand)
                  .select(['NDVI', 'EVI', 'LSWI', 'NIRv', 'LST', 'month'])
                  .set('system:time_start', startDate.millis())
                  .set('month', monthNumber);
}


// ============================================================
// STEP 6 — BORO SEASON MONTHLY COLLECTION (Nov → May)
// ============================================================
function boroMonthlyCollection(boroYear) {
  boroYear = ee.Number(boroYear);
  var start = ee.Date.fromYMD(boroYear.subtract(1), 11, 1);
  var months = ee.List.sequence(0, 6);

  var images = months.map(function(offset) {
    var mStart = start.advance(ee.Number(offset), 'month');
    var mEnd = mStart.advance(1, 'month');
    return monthlyLandsatComposite(mStart, mEnd).set('boro_year', boroYear);
  });

  return ee.ImageCollection.fromImages(images);
}


// ============================================================
// STEP 7 — 30m BORO PHENOLOGY METRICS
// ============================================================
function boroPhenology(boroYear) {
  var col = boroMonthlyCollection(boroYear);
  
  var maxNDVI = col.select('NDVI').max().rename('maxNDVI');
  var minNDVI = col.select('NDVI').min().rename('minNDVI');
  var amplitude = maxNDVI.subtract(minNDVI).rename('NDVI_amplitude');
  var integratedNDVI = col.select('NDVI').sum().rename('integratedNDVI');
  var maxNIRv = col.select('NIRv').max().rename('maxNIRv');
  var meanLST = col.select('LST').mean().rename('meanLST');
  
  var peakImage = col.qualityMosaic('NDVI');
  var peakMonth = peakImage.select('month').rename('peakMonth');

  return ee.Image.cat([maxNDVI, minNDVI, amplitude, integratedNDVI, maxNIRv, meanLST, peakMonth])
                 .set('boro_year', boroYear);}
                 
                 
// ============================================================
// STEP 8 — ADD LSWI PEAK TO BORO PHENOLOGY
// ============================================================
function boroPhenology2(boroYear) {

  var col = boroMonthlyCollection(boroYear);

  var maxNDVI = col.select('NDVI').max().rename('maxNDVI');
  var minNDVI = col.select('NDVI').min().rename('minNDVI');

  var amplitude = maxNDVI
    .subtract(minNDVI)
    .rename('NDVI_amplitude');

  var integratedNDVI = col.select('NDVI')
    .sum()
    .rename('integratedNDVI');

  var maxNIRv = col.select('NIRv')
    .max()
    .rename('maxNIRv');

  var maxLSWI = col.select('LSWI')
    .max()
    .rename('maxLSWI');

  var meanLST = col.select('LST')
    .mean()
    .rename('meanLST');

  var peakImage = col.qualityMosaic('NDVI');

  var peakMonth = peakImage.select('month')
    .rename('peakMonth');

  return ee.Image.cat([
    maxNDVI,
    minNDVI,
    amplitude,
    integratedNDVI,
    maxNIRv,
    maxLSWI,
    meanLST,
    peakMonth
  ]).set('boro_year', boroYear);
}


// ============================================================
// STEP 9 — BORO CLASSIFICATION
// ============================================================
function boroMask(boroYear) {

  var p = boroPhenology2(boroYear);

  var boro = p.select('maxNDVI').gt(0.55)
    .and(p.select('NDVI_amplitude').gt(0.25))
    .and(p.select('maxLSWI').gt(0.05))

    // Main Boro peak generally occurs during Feb–Apr
    .and(p.select('peakMonth').gte(2))
    .and(p.select('peakMonth').lte(4));

  return boro
    .rename('Boro')
    .set('boro_year', boroYear);
}


// ============================================================
// STEP 10 — TEST BORO MAP
// ============================================================
var boro2025 = boroMask(2025);

Map.addLayer(
  boro2025.selfMask(),
  {palette: ['00AA00']},
  'Boro 2025'
);

print('Boro 2025:', boro2025);


// ============================================================
// STEP 11 — CREATE ANNUAL BORO MAPS
// ============================================================
var years = ee.List.sequence(2001, 2025);

var boroIC = ee.ImageCollection.fromImages(
  years.map(function(y) {
    return boroMask(ee.Number(y))
      .set('system:time_start',
        ee.Date.fromYMD(ee.Number(y), 4, 1).millis());
  })
);

print('Annual Boro collection:', boroIC);


// ============================================================
// STEP 12 — HISTORICAL AGRICULTURAL LAND MASK
// ============================================================
// We use Landsat phenology over the full study period to identify
// land that repeatedly shows agricultural/cultivated behavior.
//
// This is a stable denominator for comparing Boro intensity
// through time. It is NOT a yearly cropland map.

// Create annual agricultural proxy
function agriculturalMask(year) {

  year = ee.Number(year);

  // Full agricultural year
  var start = ee.Date.fromYMD(year, 1, 1);
  var end   = ee.Date.fromYMD(year, 12, 31);

  var col = landsatCol
    .filterBounds(studyGeometry)
    .filterDate(start, end)
    .map(addLandsatIndices);

  var maxNDVI = col.select('NDVI').max();
  var minNDVI = col.select('NDVI').min();

  var amplitude = maxNDVI.subtract(minNDVI);

  // Agricultural land generally has repeated vegetation signal
  var ag = maxNDVI.gt(0.55)
    .and(amplitude.gt(0.20));

  return ag.rename('Agriculture')
    .set('year', year);
}


// ============================================================
// STEP 13 — CREATE STABLE AGRICULTURAL MASK
// ============================================================
// Pixel is considered agricultural if classified as agricultural
// in at least 40% of available years.

var agIC = ee.ImageCollection.fromImages(
  years.map(function(y) {
    return agriculturalMask(ee.Number(y));
  })
);

var agriculturalFrequency = agIC.mean();

var agriculturalMaskStable = agriculturalFrequency
  .gte(0.40)
  .rename('Agriculture');

Map.addLayer(
  agriculturalMaskStable.selfMask(),
  {palette: ['yellow']},
  'Stable Agricultural Land'
);


// ============================================================
// STEP 14 — CALCULATE BORO INTENSITY
// ============================================================
// Boro intensity (%) = Boro area / Agricultural area × 100

var pixelAreaHa = ee.Image.pixelArea()
  .divide(10000);

var annualResults = ee.FeatureCollection(
  years.map(function(y) {

    y = ee.Number(y);

    var boro = boroIC
      .filter(ee.Filter.eq('boro_year', y))
      .first();

    // Boro area
    var boroArea = boro
      .multiply(pixelAreaHa)
      .rename('BoroArea');

    // Agricultural area
    var agArea = agriculturalMaskStable
      .multiply(pixelAreaHa)
      .rename('AgArea');

    // Combine
    var image = boroArea
      .addBands(agArea);

    return image.reduceRegions({
      collection: studyDistricts,
      reducer: ee.Reducer.sum(),
      scale: 30,
      tileScale: 4
    }).map(function(f) {

      var boroHa = ee.Number(f.get('BoroArea'));
      var agHa   = ee.Number(f.get('AgArea'));

      var intensity = boroHa
        .divide(agHa)
        .multiply(100);

      return f.set({
        'year': y,
        'district': f.get('GAUL2_NAME'),
        'Boro_area_ha': boroHa,
        'Agricultural_area_ha': agHa,
        'Boro_intensity_pct': intensity
      });
    });
  })
).flatten();


// ============================================================
// STEP 15 — VIEW RESULTS
// ============================================================
print(
  'Annual Boro intensity:',
  annualResults
);


// ============================================================
// STEP 16 — EXPORT
// ============================================================
Export.table.toDrive({
  collection: annualResults,
  description: 'Boro_Intensity_Landsat_2001_2025',
  folder: 'Boro_Intensification',
  fileNamePrefix: 'Boro_Intensity_District_2001_2025',
  fileFormat: 'CSV'
});

// ============================================================
// STEP 17 — CREATE PIXEL-LEVEL BORO INTENSITY RASTER (2001–2025)
// ============================================================
// Calculates how frequently each 30m pixel was cultivated in Boro (0 to 100%)

var boroFrequencyRaster = boroIC.mean()
  .multiply(100)
  .rename('Boro_Intensity_Pct')
  .clip(studyGeometry);

// Visualize intensity on the map (0% = transparent, 100% = dense green)
Map.addLayer(
  boroFrequencyRaster.selfMask(),
  {min: 10, max: 100, palette: ['#ffffcc', '#a1dab4', '#41b6c4', '#225ea8']},
  'Boro Multi-Year Cultivation Intensity (%)'
);


// ============================================================
// STEP 18 — EXPORT INTENSITY RASTER TIF TO GOOGLE DRIVE
// ============================================================
Export.image.toDrive({
  image: boroFrequencyRaster.toByte(),
  description: 'Boro_MultiYear_Intensity_Raster_2001_2025',
  folder: 'Nature_sustain_boro',
  region: studyGeometry,
  scale: 30, // 30 meters resolution
  crs: 'EPSG:4326',
  maxPixels: 1e10
});

// ============================================================
// STEP 8 — 30m BORO RICE CANDIDATE MASK & AREA (2024 Benchmark)
// ============================================================
var phenology2024 = boroPhenology(2024);

var boroCandidate2024 = phenology2024.select('maxNDVI').gte(0.55)
  .and(phenology2024.select('NDVI_amplitude').gte(0.20))
  .and(phenology2024.select('peakMonth').gte(3))
  .and(phenology2024.select('peakMonth').lte(5))
  .rename('BoroCandidate')
  .clip(studyGeometry);

Map.addLayer(boroCandidate2024.selfMask(), {min: 0, max: 1, palette: ['green']}, 'Landsat 30m Boro Mask 2024');


// ============================================================
// ACCURACY ASSESSMENT — STRATIFIED RANDOM SAMPLING
// ============================================================

// Generate stratified random points based on the 2024 Boro candidate classification mask
var validationPoints = boroCandidate2024.stratifiedSample({
  numPoints: 250,              // Number of points per class
  classBand: 'BoroCandidate',   // Band representing classification (0 or 1)
  region: studyGeometry,
  scale: 30,
  geometries: true,
  seed: 42                     // Ensures reproducibility
});

// Export validation points to Google Drive as a CSV
Export.table.toDrive({
  collection: validationPoints,
  description: 'Boro_Classification_Validation_Points_2024',
  folder: 'Nature_sustain_boro',
  fileNamePrefix: 'Boro_Validation_Points_2024',
  fileFormat: 'CSV'
});

print('Validation points generated successfully. Total sample count:', validationPoints.size());