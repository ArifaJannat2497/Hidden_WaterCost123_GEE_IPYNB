/***************************************************************************
 * PRODUCTIVITY–WATER TRAJECTORIES OF BORO RICE, 2001–05 vs 2021–25
 * Study districts: Mymensingh, Sunamganj, Rangpur, Rajshahi (Bangladesh)
 *
 * Standalone Google Earth Engine (JavaScript) script. It produces:
 *   A. ET-based trajectory classes (ERA5-Land ET)            -> main Fig. 5
 *   B. LST-based trajectory classes (Landsat LST, 30 m)      -> Supp. Fig. 7
 *   C. Class areas per district with correct pixel areas    -> Supp. Tables 10–11
 *   D. Diagnostics: ET-change layer, valid-month counts,
 *      productivity change on the Boro footprint vs all pixels
 *   E. GeoTIFFs in UTM 46N (true 30-m pixels) for ArcGIS
 *
 * Class codes (same for ET and LST versions):
 *   1 = productivity up,   water use down
 *   2 = productivity up,   water use up
 *   3 = productivity down, water use up
 *   4 = productivity down, water use down
 * Water use: ET version -> dET > 0 means water use up.
 *            LST version -> dLST < 0 means water use up (cooler surface =
 *            more evaporation), so LST is sign-flipped before classifying.
 *            NOTE: this makes LST classes 2 and 4 swap relative to the old
 *            "Prod/LST" legend (old Type 4 Prod-/LST- = new class 3).
 ***************************************************************************/

// ============================================================
// 0. SETTINGS
// ============================================================
var BASE_START = 2001, BASE_END = 2005;
var REC_START  = 2021, REC_END  = 2025;
var FOOTPRINT_YEAR = 2024;
// 'BI'  = same rule as the Boro-intensity paper (NDVI>0.55, amp>0.25, LSWI>0.05, peak Feb–Apr)
// 'OLD' = rule used for the earlier trajectory maps (NDVI>=0.55, amp>=0.20, peak Mar–May)
var FOOTPRINT_RULE = 'BI';
// Productivity metric: 'integratedNDVI' (sum of monthly NDVI; used so far)
//                      'meanNDVI' (mean of monthly NDVI; robust to missing months)
var PROD_METRIC = 'integratedNDVI';
var FOLDER = 'Nature_sustain_boro';
var UTM = 'EPSG:32646';

// ============================================================
// 1. STUDY DISTRICTS
// ============================================================
var districts = ee.FeatureCollection('FAO/GAUL/2025/level2')
  .filter(ee.Filter.eq('GAUL0_NAME', 'Bangladesh'))
  .filter(ee.Filter.inList('GAUL2_NAME', ['Rajshahi', 'Rangpur', 'Mymensingh', 'Sunamganj']));
var region = districts.geometry();
Map.centerObject(districts, 7);
Map.addLayer(districts.style({color: 'black', fillColor: '00000000', width: 1}), {}, 'Districts');

// ============================================================
// 2. LANDSAT 5/7/8/9 COLLECTION 2 LEVEL-2 (harmonised band names)
// ============================================================
function prepL57(img) {
  var qa = img.select('QA_PIXEL');
  var clear = qa.bitwiseAnd(1 << 3).eq(0).and(qa.bitwiseAnd(1 << 4).eq(0));
  var opt = img.select(['SR_B1', 'SR_B3', 'SR_B4', 'SR_B5']).multiply(0.0000275).add(-0.2)
               .rename(['Blue', 'Red', 'NIR', 'SWIR1']);
  var lst = img.select('ST_B6').multiply(0.00341802).add(149.0).rename('LST');
  return opt.addBands(lst).updateMask(clear).copyProperties(img, ['system:time_start']);
}
function prepL89(img) {
  var qa = img.select('QA_PIXEL');
  var clear = qa.bitwiseAnd(1 << 3).eq(0).and(qa.bitwiseAnd(1 << 4).eq(0));
  var opt = img.select(['SR_B2', 'SR_B4', 'SR_B5', 'SR_B6']).multiply(0.0000275).add(-0.2)
               .rename(['Blue', 'Red', 'NIR', 'SWIR1']);
  var lst = img.select('ST_B10').multiply(0.00341802).add(149.0).rename('LST');
  return opt.addBands(lst).updateMask(clear).copyProperties(img, ['system:time_start']);
}
function addIdx(img) {
  var ndvi = img.normalizedDifference(['NIR', 'Red']).rename('NDVI');
  var lswi = img.normalizedDifference(['NIR', 'SWIR1']).rename('LSWI');
  var nirv = img.select('NIR').multiply(ndvi).rename('NIRv');
  return img.addBands([ndvi, lswi, nirv]);
}
var landsat = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2').map(prepL57)
  .merge(ee.ImageCollection('LANDSAT/LE07/C02/T1_L2').map(prepL57))
  .merge(ee.ImageCollection('LANDSAT/LC08/C02/T1_L2').map(prepL89))
  .merge(ee.ImageCollection('LANDSAT/LC09/C02/T1_L2').map(prepL89))
  .filterBounds(region)
  .map(addIdx);

// ============================================================
// 3. MONTHLY COMPOSITES FOR ONE BORO SEASON (Nov of y-1 to May of y)
// ============================================================
function boroMonths(y) {
  y = ee.Number(y);
  var start = ee.Date.fromYMD(y.subtract(1), 11, 1);
  return ee.ImageCollection.fromImages(ee.List.sequence(0, 6).map(function (k) {
    var s = start.advance(ee.Number(k), 'month');
    var comp = landsat.filterDate(s, s.advance(1, 'month'))
      .select(['NDVI', 'LSWI', 'NIRv', 'LST']).median();
    var month = ee.Image.constant(s.get('month')).toFloat().rename('month');
    // guarantee the bands exist even when a month has no scene
    var empty = ee.Image.constant([0, 0, 0, 0]).rename(['NDVI', 'LSWI', 'NIRv', 'LST']).toFloat().updateMask(0);
    comp = ee.Image(ee.Algorithms.If(comp.bandNames().size().gt(0), comp.toFloat(), empty));
    return comp.addBands(month.updateMask(comp.select('NDVI').mask()))
               .set('system:time_start', s.millis());
  }));
}

// ============================================================
// 4. SEASONAL METRICS PER PIXEL
// ============================================================
function seasonMetrics(y) {
  var col = boroMonths(y);
  var ndvi = col.select('NDVI');
  var maxNDVI = ndvi.max().rename('maxNDVI');
  var amp = maxNDVI.subtract(ndvi.min()).rename('amp');
  return ee.Image.cat([
    maxNDVI, amp,
    ndvi.sum().rename('integratedNDVI'),
    ndvi.mean().rename('meanNDVI'),
    ndvi.count().rename('nMonths'),
    col.select('LSWI').max().rename('maxLSWI'),
    col.select('LST').mean().rename('meanLST'),
    col.qualityMosaic('NDVI').select('month').rename('peakMonth')
  ]).set('year', y);
}

// ============================================================
// 5. 2024 BORO FOOTPRINT
// ============================================================
var fp = seasonMetrics(FOOTPRINT_YEAR);
var footprint = (FOOTPRINT_RULE === 'BI')
  ? fp.select('maxNDVI').gt(0.55).and(fp.select('amp').gt(0.25))
      .and(fp.select('maxLSWI').gt(0.05))
      .and(fp.select('peakMonth').gte(2)).and(fp.select('peakMonth').lte(4))
  : fp.select('maxNDVI').gte(0.55).and(fp.select('amp').gte(0.20))
      .and(fp.select('peakMonth').gte(3)).and(fp.select('peakMonth').lte(5));
footprint = footprint.selfMask().rename('footprint').clip(region);
Map.addLayer(footprint, {palette: ['00AA00']}, 'Boro footprint ' + FOOTPRINT_YEAR, false);

// ============================================================
// 6. PRODUCTIVITY AND LST CHANGE (recent mean minus baseline mean)
// ============================================================
function windowMean(y0, y1) {
  return ee.ImageCollection(ee.List.sequence(y0, y1).map(seasonMetrics)).mean();
}
var baseL = windowMean(BASE_START, BASE_END);
var recL  = windowMean(REC_START, REC_END);
var dProd = recL.select(PROD_METRIC).subtract(baseL.select(PROD_METRIC)).rename('dProd');
var dLST  = recL.select('meanLST').subtract(baseL.select('meanLST')).rename('dLST');

// ============================================================
// 7. ET CHANGE (ERA5-Land, Boro season, mm)
// ============================================================
var era5 = ee.ImageCollection('ECMWF/ERA5_LAND/MONTHLY_AGGR');
function boroET(y) {
  y = ee.Number(y);
  var s = ee.Date.fromYMD(y.subtract(1), 11, 1), e = ee.Date.fromYMD(y, 6, 1);
  return era5.filterDate(s, e).select('total_evaporation_sum')
    .map(function (i) { return i.multiply(-1000).max(0); })   // m (negative) -> mm
    .sum().rename('ET').set('year', y);
}
var etBase = ee.ImageCollection(ee.List.sequence(BASE_START, BASE_END).map(boroET)).mean();
var etRec  = ee.ImageCollection(ee.List.sequence(REC_START, REC_END).map(boroET)).mean();
var dET = etRec.subtract(etBase).rename('dET');   // native ~9-km grid

// DIAGNOSTIC: this layer must look blocky (~9 km). Every class map built from it
// must follow the same blocks: inside a block, only {1,4} (ET down) or {2,3} (ET up).
Map.addLayer(dET.clip(region), {min: -40, max: 40, palette: ['2166ac', 'f7f7f7', 'b2182b']}, 'dET (mm), ERA5-Land');

// ============================================================
// 8. CLASSIFY: 1 P+W-, 2 P+W+, 3 P-W+, 4 P-W-
// ============================================================
function classify(dP, waterUp) {          // waterUp: 1 where water use increased
  return ee.Image(0)
    .where(dP.gt(0).and(waterUp.eq(0)), 1)
    .where(dP.gt(0).and(waterUp.eq(1)), 2)
    .where(dP.lte(0).and(waterUp.eq(1)), 3)
    .where(dP.lte(0).and(waterUp.eq(0)), 4)
    .updateMask(footprint).clip(region).toByte();
}
// dET stays on its native ERA5 grid; GEE samples it at 30 m by nearest neighbour,
// so every 30-m pixel inside an ERA5 cell gets that cell's value.
var classET  = classify(dProd, dET.gt(0)).rename('class_ET');
var classLST = classify(dProd, dLST.lt(0)).rename('class_LST');   // cooler = more water use

var pal = ['009E73', 'E69F00', 'CC79A7', '0072B2'];   // Okabe–Ito, colour-blind safe
Map.addLayer(classET,  {min: 1, max: 4, palette: pal}, 'A. ET-based trajectories');
Map.addLayer(classLST, {min: 1, max: 4, palette: pal}, 'B. LST-based trajectories', false);

// ============================================================
// 9. CLASS AREAS PER DISTRICT (ha, from ee.Image.pixelArea)
// ============================================================
function classAreas(classImg, label) {
  var img = ee.Image.pixelArea().divide(1e4).rename('ha').addBands(classImg);
  return districts.map(function (d) {
    var groups = ee.List(img.reduceRegion({
      reducer: ee.Reducer.sum().group({groupField: 1, groupName: 'class'}),
      geometry: d.geometry(), scale: 30, crs: UTM, maxPixels: 1e10, tileScale: 8
    }).get('groups'));
    var props = ee.Dictionary(groups.iterate(function (g, acc) {
      g = ee.Dictionary(g);
      return ee.Dictionary(acc).set(
        ee.String('class').cat(ee.Number(g.get('class')).format('%d')).cat('_ha'), g.get('sum'));
    }, ee.Dictionary({})));
    var district_ha = ee.Image.pixelArea().divide(1e4).reduceRegion({
      reducer: ee.Reducer.sum(), geometry: d.geometry(), scale: 30, crs: UTM,
      maxPixels: 1e10, tileScale: 8}).get('area');
    return ee.Feature(null, props).set({
      District: d.get('GAUL2_NAME'), Version: label, district_ha: district_ha});
  });
}
var areasET  = classAreas(classET, 'ET');
var areasLST = classAreas(classLST, 'LST');
print('Class areas, ET version (ha)', areasET);
print('Class areas, LST version (ha)', areasLST);

// ============================================================
// 10. DIAGNOSTICS PER DISTRICT
// ============================================================
var diagImg = ee.Image.cat([
  dProd, dLST,
  baseL.select('nMonths').rename('nMonths_base'),
  recL.select('nMonths').rename('nMonths_recent'),
  dProd.gt(0).rename('share_prod_up')]);
var diag = districts.map(function (d) {
  var onFp = diagImg.updateMask(footprint).reduceRegion({
    reducer: ee.Reducer.mean(), geometry: d.geometry(), scale: 30, crs: UTM,
    maxPixels: 1e10, tileScale: 8});
  var allPx = diagImg.reduceRegion({
    reducer: ee.Reducer.mean(), geometry: d.geometry(), scale: 30, crs: UTM,
    maxPixels: 1e10, tileScale: 8});
  var et = dET.reduceRegion({
    reducer: ee.Reducer.minMax().combine(ee.Reducer.mean(), '', true),
    geometry: d.geometry(), scale: 11132});
  var fpHa = ee.Image.pixelArea().divide(1e4).updateMask(footprint).reduceRegion({
    reducer: ee.Reducer.sum(), geometry: d.geometry(), scale: 30, crs: UTM,
    maxPixels: 1e10, tileScale: 8}).get('area');
  return ee.Feature(null, {
    District: d.get('GAUL2_NAME'), footprint_ha: fpHa,
    fp_dProd: onFp.get('dProd'), all_dProd: allPx.get('dProd'),
    fp_share_prod_up: onFp.get('share_prod_up'), all_share_prod_up: allPx.get('share_prod_up'),
    fp_nMonths_base: onFp.get('nMonths_base'), fp_nMonths_recent: onFp.get('nMonths_recent'),
    fp_dLST: onFp.get('dLST'),
    dET_min: et.get('dET_min'), dET_max: et.get('dET_max'), dET_mean: et.get('dET_mean')});
});
print('Diagnostics', diag);

// ============================================================
// 11. EXPORTS
// ============================================================
Export.table.toDrive({collection: areasET.merge(areasLST), description: 'Trajectory_class_areas_ET_LST',
  folder: FOLDER, fileFormat: 'CSV'});
Export.table.toDrive({collection: diag, description: 'Trajectory_diagnostics',
  folder: FOLDER, fileFormat: 'CSV'});
Export.image.toDrive({image: classET, description: 'ET_trajectory_classes_30m_UTM46N',
  folder: FOLDER, region: region, scale: 30, crs: UTM, maxPixels: 1e13});
Export.image.toDrive({image: classLST, description: 'LST_trajectory_classes_30m_UTM46N',
  folder: FOLDER, region: region, scale: 30, crs: UTM, maxPixels: 1e13});
Export.image.toDrive({image: dET.clip(region).toFloat(), description: 'dET_ERA5Land_mm',
  folder: FOLDER, region: region, scale: 11132, crs: 'EPSG:4326', maxPixels: 1e10});

/* HOW TO READ THE OUTPUTS
 * - 'A. ET-based trajectories' should look blocky (~9-km cells): within a block
 *   only classes 1+4 or 2+3 occur. If your ArcGIS map is not blocky, it was not
 *   made from ET_trajectory_classes_30m_UTM46N.
 * - dET_min/dET_max per district: if both have the same sign, ET changed in one
 *   direction across the whole district and the ET classes there reduce to a
 *   productivity-up/down split.
 * - fp_dProd < 0 but all_dProd > 0: the district-mean rise comes from non-Boro land.
 * - fp_nMonths_recent differs clearly from fp_nMonths_base: integrated NDVI is
 *   biased by data availability -> set PROD_METRIC = 'meanNDVI' and re-run.
 * Send me: Trajectory_class_areas_ET_LST.csv and Trajectory_diagnostics.csv.
 */
