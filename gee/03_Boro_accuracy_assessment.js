/***************************************************************************
 * ACCURACY ASSESSMENT OF THE BORO RICE MAP (Supplementary Table 2)
 * Study districts: Mymensingh, Sunamganj, Rangpur, Rajshahi (Bangladesh)
 *
 * Good-practice design (Olofsson et al. 2014, RSE 148:42–57; Stehman 2014,
 * IJRS 35:4923–4939):
 *   - Map classes (= strata within each district):
 *       1 = Boro rice (same rule as the Boro-intensity panel)
 *       2 = Other agriculture (agricultural mask, not Boro)
 *       3 = Non-agriculture (everything else: settlement, trees, water, sand)
 *   - Stratified random sample by district x map class (12 strata)
 *   - Reference labels by BLIND visual interpretation of Sentinel-2
 *     (10 m) and Sentinel-1 (radar) time series plus high-resolution
 *     basemap imagery - i.e., data independent of the Landsat rule
 *   - Area-weighted error matrix, user's / producer's accuracy, F1,
 *     overall accuracy and area-adjusted class areas with 95% CIs
 *   - Kappa is NOT reported (Pontius & Millones 2011; Olofsson et al. 2014)
 *
 * RUN IN THREE STEPS - change MODE and re-run each time:
 *   MODE = 'SAMPLE' -> exports (i) sample points and (ii) stratum areas
 *                      (Tasks tab -> Run both; they go to Drive + Assets)
 *   MODE = 'LABEL'  -> labelling app. Label every point, then click
 *                      "Export labels" (Tasks tab -> Run)
 *   MODE = 'ASSESS' -> error matrix + accuracy table (Console) and CSVs
 *                      for Supplementary Table 2
 ***************************************************************************/

// ============================================================
// 0. SETTINGS
// ============================================================
var MODE = 'SAMPLE';                    // 'SAMPLE' | 'LABEL' | 'ASSESS'
var YEAR = 2024;                        // Boro season Nov (YEAR-1) – May (YEAR); >= 2019 for Sentinel-2 L2A
var SEED = 2024;
var FOLDER = 'Nature_sustain_boro';
var ASSET_ROOT = 'projects/YOUR_PROJECT/assets/';   // <-- edit
var UTM = 'EPSG:32646';

// Sample size per district and map class (total 4 x 225 = 900 points).
// Boro gets most points because its accuracy matters most; 50 points is
// the minimum per class recommended by Olofsson et al. (2014).
var N_PER_CLASS = {1: 100, 2: 75, 3: 50};

// Optional: agricultural mask exported from the Boro-intensity script
// (agriculturalMaskStable, 1 = agriculture). Leave null to recompute here.
var AG_MASK_ASSET = null;

// Asset names written in step 1 and 2 (read in steps 2 and 3)
var POINTS_ASSET = ASSET_ROOT + 'Boro_accuracy_points_' + YEAR;
var AREAS_ASSET  = ASSET_ROOT + 'Boro_accuracy_strata_areas_' + YEAR;
var LABELS_ASSET = ASSET_ROOT + 'Boro_accuracy_labels_' + YEAR;

// To resume labelling: paste the JSON printed by "Show labels" here.
var RESUME_LABELS = {};

// ============================================================
// 1. STUDY DISTRICTS
// ============================================================
var DCODE = {Mymensingh: 1, Sunamganj: 2, Rangpur: 3, Rajshahi: 4};
var DNAME = {1: 'Mymensingh', 2: 'Sunamganj', 3: 'Rangpur', 4: 'Rajshahi'};
var CNAME = {1: 'Boro rice', 2: 'Other agriculture', 3: 'Non-agriculture'};
var districts = ee.FeatureCollection('FAO/GAUL/2025/level2')
  .filter(ee.Filter.eq('GAUL0_NAME', 'Bangladesh'))
  .filter(ee.Filter.inList('GAUL2_NAME', Object.keys(DCODE)))
  .map(function (f) { return f.set('dcode', ee.Dictionary(DCODE).get(f.get('GAUL2_NAME'))); });
var region = districts.geometry();

// ============================================================
// 2. LANDSAT 5/7/8/9 C2 L2 (identical to the trajectory script)
// ============================================================
function prepL57(img) {
  var qa = img.select('QA_PIXEL');
  var clear = qa.bitwiseAnd(1 << 3).eq(0).and(qa.bitwiseAnd(1 << 4).eq(0));
  var opt = img.select(['SR_B1', 'SR_B3', 'SR_B4', 'SR_B5']).multiply(0.0000275).add(-0.2)
               .rename(['Blue', 'Red', 'NIR', 'SWIR1']);
  return opt.updateMask(clear).copyProperties(img, ['system:time_start']);
}
function prepL89(img) {
  var qa = img.select('QA_PIXEL');
  var clear = qa.bitwiseAnd(1 << 3).eq(0).and(qa.bitwiseAnd(1 << 4).eq(0));
  var opt = img.select(['SR_B2', 'SR_B4', 'SR_B5', 'SR_B6']).multiply(0.0000275).add(-0.2)
               .rename(['Blue', 'Red', 'NIR', 'SWIR1']);
  return opt.updateMask(clear).copyProperties(img, ['system:time_start']);
}
function addIdx(img) {
  return img.addBands([img.normalizedDifference(['NIR', 'Red']).rename('NDVI'),
                       img.normalizedDifference(['NIR', 'SWIR1']).rename('LSWI')]);
}
var landsat = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2').map(prepL57)
  .merge(ee.ImageCollection('LANDSAT/LE07/C02/T1_L2').map(prepL57))
  .merge(ee.ImageCollection('LANDSAT/LC08/C02/T1_L2').map(prepL89))
  .merge(ee.ImageCollection('LANDSAT/LC09/C02/T1_L2').map(prepL89))
  .filterBounds(region).map(addIdx).select(['NDVI', 'LSWI']);

function monthly(start, nMonths) {
  return ee.ImageCollection.fromImages(ee.List.sequence(0, nMonths - 1).map(function (k) {
    var s = ee.Date(start).advance(ee.Number(k), 'month');
    var comp = landsat.filterDate(s, s.advance(1, 'month')).median();
    var empty = ee.Image.constant([0, 0]).rename(['NDVI', 'LSWI']).toFloat().updateMask(0);
    comp = ee.Image(ee.Algorithms.If(comp.bandNames().size().gt(0), comp.toFloat(), empty));
    var month = ee.Image.constant(s.get('month')).toFloat().rename('month')
                  .updateMask(comp.select('NDVI').mask());
    return comp.addBands(month).set('system:time_start', s.millis());
  }));
}

// ============================================================
// 3. MAP CLASSES FOR THE ASSESSED SEASON
// ============================================================
// 3a. Boro (same rule as the Boro-intensity panel and Fig. 5 footprint)
var col = monthly(ee.Date.fromYMD(YEAR - 1, 11, 1), 7);
var maxNDVI = col.select('NDVI').max();
var amp = maxNDVI.subtract(col.select('NDVI').min());
var maxLSWI = col.select('LSWI').max();
var peak = col.qualityMosaic('NDVI').select('month');
var boro = maxNDVI.gt(0.55).and(amp.gt(0.25)).and(maxLSWI.gt(0.05))
  .and(peak.gte(2)).and(peak.lte(4)).unmask(0);

// 3b. Agricultural land — identical to Steps 12–13 of the Boro-intensity script:
//     scene-level calendar-year max NDVI > 0.55 and amplitude > 0.20,
//     in >= 40% of years with data (2001–2025)
var agMask;
if (AG_MASK_ASSET) {
  agMask = ee.Image(AG_MASK_ASSET).gt(0).unmask(0);
} else {
  var agIC = ee.ImageCollection.fromImages(ee.List.sequence(2001, 2025).map(function (y) {
    y = ee.Number(y);
    var c = landsat.filterDate(ee.Date.fromYMD(y, 1, 1), ee.Date.fromYMD(y, 12, 31)).select('NDVI');
    var mx = c.max(), mn = c.min();
    return mx.gt(0.55).and(mx.subtract(mn).gt(0.20)).rename('Agriculture');
  }));
  agMask = agIC.mean().gte(0.40).unmask(0);
}

// 3c. Map class and stratum (district code x 10 + map class)
var mapClass = ee.Image(3).where(agMask, 2).where(boro, 1).rename('mapClass').toInt();
var dImg = districts.reduceToImage(['dcode'], ee.Reducer.first()).rename('dcode').toInt();
var stratum = dImg.multiply(10).add(mapClass).rename('stratum').toInt()
  .updateMask(dImg.mask()).clip(region);

// ============================================================
// 4. STEP 1 — SAMPLE POINTS AND STRATUM AREAS
// ============================================================
if (MODE === 'SAMPLE') {
  Map.centerObject(districts, 7);
  Map.addLayer(mapClass.clip(region), {min: 1, max: 3, palette: ['00A600', 'FFD37F', 'BDBDBD']},
               'Map class ' + YEAR);
  Map.addLayer(districts.style({color: 'black', fillColor: '00000000', width: 1}), {}, 'Districts');

  var classValues = [], classPoints = [];
  [1, 2, 3, 4].forEach(function (d) {
    [1, 2, 3].forEach(function (c) { classValues.push(d * 10 + c); classPoints.push(N_PER_CLASS[c]); });
  });

  var pts = stratum.addBands(mapClass).addBands(dImg).stratifiedSample({
    numPoints: 0, classBand: 'stratum', region: region,
    scale: 30, projection: ee.Projection(UTM).atScale(30),
    seed: SEED, classValues: classValues, classPoints: classPoints,
    geometries: true, tileScale: 16
  });

  // Shuffle so that interpreters meet points in random order (blind labelling)
  pts = pts.randomColumn('rnd', SEED).sort('rnd');
  var lst = pts.toList(5000);
  pts = ee.FeatureCollection(ee.List.sequence(0, lst.size().subtract(1)).map(function (i) {
    var f = ee.Feature(lst.get(i));
    var xy = f.geometry().coordinates();
    return f.set({pid: ee.Number(i).add(1), lon: xy.get(0), lat: xy.get(1),
                  district: ee.Dictionary(DNAME).get(ee.Number(f.get('dcode')).format('%d')),
                  ref: -1});
  })).select(['pid', 'stratum', 'dcode', 'district', 'mapClass', 'lon', 'lat', 'ref']);

  Export.table.toAsset({collection: pts, description: 'Boro_accuracy_points_' + YEAR,
                        assetId: POINTS_ASSET});
  Export.table.toDrive({collection: pts, description: 'Boro_accuracy_points_' + YEAR + '_csv',
                        folder: FOLDER, fileFormat: 'CSV'});

  // Mapped area of each stratum (true pixel areas in UTM 46N)
  var areaDict = ee.Image.pixelArea().divide(1e4).rename('ha').addBands(stratum)
    .reduceRegion({reducer: ee.Reducer.sum().group({groupField: 1, groupName: 'stratum'}),
                   geometry: region, scale: 30, crs: UTM, maxPixels: 1e13, tileScale: 16});
  var areas = ee.FeatureCollection(ee.List(areaDict.get('groups')).map(function (g) {
    g = ee.Dictionary(g);
    return ee.Feature(null, {stratum: g.get('stratum'), area_ha: g.get('sum')});
  }));
  Export.table.toAsset({collection: areas, description: 'Boro_accuracy_strata_areas_' + YEAR,
                        assetId: AREAS_ASSET});
  Export.table.toDrive({collection: areas, description: 'Boro_accuracy_strata_areas_' + YEAR + '_csv',
                        folder: FOLDER, fileFormat: 'CSV'});
  print('Step 1: run the four export tasks in the Tasks tab, then set MODE = "LABEL".');
}

// ============================================================
// 5. STEP 2 — LABELLING APP (blind to the map class)
// ============================================================
/* Interpretation key (label the 30-m cell outlined in red, by majority cover
 * during the Boro season Dec (YEAR-1) – Apr (YEAR)):
 *  1 Boro rice      : field flooded/puddled at transplanting (Sentinel-1 VV
 *                     drops, typically < -17 dB, Dec–Feb), then NDVI rises to a
 *                     peak (> ~0.6) in Feb–Apr and falls at harvest (Apr–May);
 *                     bright red, uniform paddy plots in the false-colour images.
 *  2 Other agric.   : cropland not under Boro in this season – mustard, wheat,
 *                     maize, potato, vegetables (NDVI peak Dec–Feb, no flooding),
 *                     or cropland left fallow.
 *  3 Non-agriculture: settlements and homestead trees (high NDVI all season),
 *                     roads, rivers, beels/haors under water all season, sand
 *                     bars, forest, brick kilns.
 *  9 Unsure         : revisit later (excluded from the estimates if unresolved).
 */
if (MODE === 'LABEL') {
  var points = ee.FeatureCollection(POINTS_ASSET);
  var start = ee.Date.fromYMD(YEAR - 1, 11, 1), end = ee.Date.fromYMD(YEAR, 6, 1);

  // Sentinel-2 L2A, cloud-masked with SCL
  var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
    .filterBounds(region).filterDate(start, end)
    .map(function (i) {
      var scl = i.select('SCL');
      var ok = scl.eq(4).or(scl.eq(5)).or(scl.eq(6)).or(scl.eq(7)).or(scl.eq(2));
      return i.updateMask(ok).divide(10000)
        .addBands(i.normalizedDifference(['B8', 'B4']).rename('NDVI'))
        .copyProperties(i, ['system:time_start']);
    });
  var s1 = ee.ImageCollection('COPERNICUS/S1_GRD')
    .filterBounds(region).filterDate(start.advance(-1, 'month'), end)
    .filter(ee.Filter.eq('instrumentMode', 'IW'))
    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
    .select('VV');

  Map.setOptions('SATELLITE');
  var months = [[YEAR - 1, 12], [YEAR, 1], [YEAR, 2], [YEAR, 3], [YEAR, 4]];
  months.forEach(function (ym, k) {
    var s = ee.Date.fromYMD(ym[0], ym[1], 1);
    var comp = s2.filterDate(s, s.advance(1, 'month')).median();
    Map.addLayer(comp, {bands: ['B8', 'B4', 'B3'], min: 0.02, max: 0.40},
                 'S2 false colour ' + ym[0] + '-' + (ym[1] < 10 ? '0' : '') + ym[1], k === 3);
  });

  var labels = RESUME_LABELS;
  var feats = [];
  var idx = 0;
  var cellLayer = null;

  var panel = ui.Panel({style: {width: '420px', padding: '8px'}});
  ui.root.insert(0, panel);
  var title = ui.Label('Boro reference labelling – season ' + (YEAR - 1) + '/' + YEAR,
                       {fontWeight: 'bold', fontSize: '15px'});
  var info = ui.Label('Loading points…');
  var prog = ui.Label('');
  var chartBox = ui.Panel();
  var jump = ui.Textbox({placeholder: 'Go to point id…', onChange: function (v) {
    var k = Number(v) - 1; if (k >= 0 && k < feats.length) { idx = k; show(); }
  }});

  function setLabel(code) {
    labels[String(feats[idx].properties.pid)] = code;
    // move to next unlabelled point
    for (var k = 1; k <= feats.length; k++) {
      var j = (idx + k) % feats.length;
      if (!(String(feats[j].properties.pid) in labels)) { idx = j; show(); return; }
    }
    show();
  }
  var btns = ui.Panel([
    ui.Button('1 Boro rice', function () { setLabel(1); }),
    ui.Button('2 Other agriculture', function () { setLabel(2); }),
    ui.Button('3 Non-agriculture', function () { setLabel(3); }),
    ui.Button('9 Unsure', function () { setLabel(9); })
  ], ui.Panel.Layout.flow('horizontal', true));
  var nav = ui.Panel([
    ui.Button('◀ Back', function () { idx = Math.max(0, idx - 1); show(); }),
    ui.Button('Next ▶', function () { idx = Math.min(feats.length - 1, idx + 1); show(); }),
    ui.Button('Show labels', function () { print('RESUME_LABELS =', JSON.stringify(labels)); }),
    ui.Button('Export labels', exportLabels)
  ], ui.Panel.Layout.flow('horizontal', true));
  panel.add(title).add(info).add(prog).add(btns).add(nav).add(jump).add(chartBox);

  function show() {
    var f = feats[idx], pid = String(f.properties.pid);
    var pt = ee.Geometry.Point(f.geometry.coordinates);
    var cell = pt.transform(UTM, 0.01).buffer(15, 0.01, UTM).bounds(0.01, UTM);
    if (cellLayer) Map.layers().remove(cellLayer);
    cellLayer = ui.Map.Layer(ee.FeatureCollection([ee.Feature(cell)])
      .style({color: 'red', fillColor: '00000000', width: 2}), {}, 'Cell');
    Map.layers().add(cellLayer);
    Map.centerObject(pt, 17);
    info.setValue('Point ' + pid + ' of ' + feats.length + '   (' +
                  f.geometry.coordinates[1].toFixed(5) + ', ' + f.geometry.coordinates[0].toFixed(5) + ')' +
                  (pid in labels ? '   current label: ' + labels[pid] : ''));
    prog.setValue('Labelled: ' + Object.keys(labels).length + ' / ' + feats.length);
    chartBox.clear();
    chartBox.add(ui.Chart.image.series({imageCollection: s2.select('NDVI'), region: cell,
        reducer: ee.Reducer.mean(), scale: 10})
      .setOptions({title: 'Sentinel-2 NDVI', vAxis: {viewWindow: {min: -0.2, max: 1}},
                   lineWidth: 1, pointSize: 3, legend: {position: 'none'}}));
    chartBox.add(ui.Chart.image.series({imageCollection: s1, region: cell.buffer(10),
        reducer: ee.Reducer.mean(), scale: 10})
      .setOptions({title: 'Sentinel-1 VV (dB): flooding shows as a dip',
                   lineWidth: 1, pointSize: 3, legend: {position: 'none'}}));
  }

  function exportLabels() {
    var done = feats.filter(function (f) { return String(f.properties.pid) in labels; })
      .map(function (f) {
        var p = f.properties;
        return ee.Feature(ee.Geometry.Point(f.geometry.coordinates),
          {pid: p.pid, stratum: p.stratum, dcode: p.dcode, district: p.district,
           mapClass: p.mapClass, ref: labels[String(p.pid)]});
      });
    var fc = ee.FeatureCollection(done);
    Export.table.toAsset({collection: fc, description: 'Boro_accuracy_labels_' + YEAR,
                          assetId: LABELS_ASSET});
    Export.table.toDrive({collection: fc, description: 'Boro_accuracy_labels_' + YEAR + '_csv',
                          folder: FOLDER, fileFormat: 'CSV'});
    print('Export tasks created for ' + done.length + ' labelled points (Tasks tab -> Run).');
  }

  // Map class is deliberately NOT loaded, so the interpreter stays blind.
  points.select(['pid', 'stratum', 'dcode', 'district', 'mapClass']).sort('pid')
    .evaluate(function (fc) {
      feats = fc.features;
      for (var k = 0; k < feats.length; k++) {
        if (!(String(feats[k].properties.pid) in labels)) { idx = k; break; }
      }
      show();
    });
}

// ============================================================
// 6. STEP 3 — ERROR MATRIX AND ACCURACY (Stehman 2014 estimators)
// ============================================================
// Generic stratified estimator; strata (district x map class) may differ
// from the reported classes. Verified against Olofsson et al. (2014),
// Table 8 example: OA 0.947, UA(1) 0.88 ± 0.07, PA(1) 0.75 ± 0.21,
// area(1) 235,086 ± 68,418 pixels.
function estimate(samples, strata, classes) {
  var Ntot = 0, hs = Object.keys(strata);
  hs.forEach(function (h) { Ntot += strata[h]; });
  function stats(fy, fx) {
    var Y = 0, X = 0, parts = [];
    hs.forEach(function (h) {
      var s = samples.filter(function (u) { return String(u.h) === String(h); });
      var n = s.length; if (n === 0) return;
      var W = strata[h] / Ntot;
      var y = s.map(fy), x = fx ? s.map(fx) : s.map(function () { return 1; });
      var my = y.reduce(function (a, b) { return a + b; }, 0) / n;
      var mx = x.reduce(function (a, b) { return a + b; }, 0) / n;
      var syy = 0, sxx = 0, sxy = 0;
      for (var k = 0; k < n; k++) {
        syy += (y[k] - my) * (y[k] - my); sxx += (x[k] - mx) * (x[k] - mx);
        sxy += (y[k] - my) * (x[k] - mx);
      }
      var d = Math.max(n - 1, 1);
      parts.push({W: W, n: n, syy: syy / d, sxx: sxx / d, sxy: sxy / d});
      Y += W * my; X += W * mx;
    });
    var R = fx ? Y / X : Y, v = 0;
    parts.forEach(function (p) {
      v += fx ? p.W * p.W * (p.syy + R * R * p.sxx - 2 * R * p.sxy) / p.n
              : p.W * p.W * p.syy / p.n;
    });
    if (fx) v = v / (X * X);
    return {est: R, ci: 1.96 * Math.sqrt(v)};
  }
  var out = {Ntot: Ntot, n: samples.length};
  out.OA = stats(function (u) { return u.m === u.r ? 1 : 0; });
  out.matrix = classes.map(function (i) {
    return classes.map(function (j) {
      return stats(function (u) { return (u.m === i && u.r === j) ? 1 : 0; }).est;
    });
  });
  out.counts = classes.map(function (i) {
    return classes.map(function (j) {
      return samples.filter(function (u) { return u.m === i && u.r === j; }).length;
    });
  });
  out.cls = classes.map(function (c) {
    var ua = stats(function (u) { return (u.m === c && u.r === c) ? 1 : 0; },
                   function (u) { return u.m === c ? 1 : 0; });
    var pa = stats(function (u) { return (u.m === c && u.r === c) ? 1 : 0; },
                   function (u) { return u.r === c ? 1 : 0; });
    var ar = stats(function (u) { return u.r === c ? 1 : 0; });
    var mp = stats(function (u) { return u.m === c ? 1 : 0; });
    return {c: c, UA: ua, PA: pa, F1: 2 * ua.est * pa.est / (ua.est + pa.est),
            mapped: mp.est * Ntot, adj: ar.est * Ntot, adjCI: ar.ci * Ntot};
  });
  // Boro intensity (Boro / agricultural land): map vs area-adjusted
  out.BI_map = out.cls[0].mapped / (out.cls[0].mapped + out.cls[1].mapped);
  out.BI_ref = stats(function (u) { return u.r === 1 ? 1 : 0; },
                     function (u) { return (u.r === 1 || u.r === 2) ? 1 : 0; });
  return out;
}

if (MODE === 'ASSESS') {
  var lab = ee.FeatureCollection(LABELS_ASSET);
  var ar = ee.FeatureCollection(AREAS_ASSET);
  ee.Dictionary({lab: lab.toList(5000).map(function (f) { return ee.Feature(f).toDictionary(); }),
                 ar: ar.toList(100).map(function (f) { return ee.Feature(f).toDictionary(); })})
  .evaluate(function (d) {
    var areas = {};
    d.ar.forEach(function (a) { areas[String(a.stratum)] = a.area_ha; });
    var all = d.lab.map(function (p) {
      return {h: String(p.stratum), m: Number(p.mapClass), r: Number(p.ref), d: Number(p.dcode)};
    });
    var sm = all.filter(function (u) { return u.r >= 1 && u.r <= 3; });
    print('Labelled points: ' + all.length + '; unsure/unlabelled excluded: ' + (all.length - sm.length));

    var rows = [];
    function pct(x) { return (100 * x).toFixed(1); }
    function report(name, samples, strata) {
      var o = estimate(samples, strata, [1, 2, 3]);
      print('==== ' + name + ' (n = ' + o.n + ') ====');
      print('Error matrix, sample counts (rows = map, columns = reference: Boro, Other agric., Non-agric.)', o.counts);
      print('Error matrix, estimated area proportions', o.matrix.map(function (r) {
        return r.map(function (v) { return Number(v.toFixed(4)); }); }));
      print('Overall accuracy: ' + pct(o.OA.est) + '% ± ' + pct(o.OA.ci));
      o.cls.forEach(function (c) {
        print(CNAME[c.c] + ': UA ' + pct(c.UA.est) + ' ± ' + pct(c.UA.ci) +
              '; PA ' + pct(c.PA.est) + ' ± ' + pct(c.PA.ci) + '; F1 ' + c.F1.toFixed(2) +
              '; mapped ' + Math.round(c.mapped) + ' ha; adjusted ' + Math.round(c.adj) +
              ' ± ' + Math.round(c.adjCI) + ' ha');
        rows.push(ee.Feature(null, {unit: name, cls: CNAME[c.c],
          UA: c.UA.est, UA_ci: c.UA.ci, PA: c.PA.est, PA_ci: c.PA.ci, F1: c.F1,
          mapped_ha: c.mapped, adjusted_ha: c.adj, adjusted_ci_ha: c.adjCI,
          OA: o.OA.est, OA_ci: o.OA.ci, n: o.n}));
      });
      print('Boro intensity: map ' + pct(o.BI_map) + '%, area-adjusted ' + pct(o.BI_ref.est) +
            '% ± ' + pct(o.BI_ref.ci));
      rows.push(ee.Feature(null, {unit: name, cls: 'Boro intensity', BI_map: o.BI_map,
                                  BI_adjusted: o.BI_ref.est, BI_adjusted_ci: o.BI_ref.ci, n: o.n}));
      var mrow = [];
      o.matrix.forEach(function (r, i) { r.forEach(function (v, j) {
        mrow.push(ee.Feature(null, {unit: name, map: CNAME[i + 1], reference: CNAME[j + 1],
                                    proportion: v, count: o.counts[i][j]})); }); });
      return mrow;
    }
    var mats = report('All four districts', sm, areas);
    [1, 2, 3, 4].forEach(function (dc) {
      var st = {};
      Object.keys(areas).forEach(function (h) { if (Math.floor(Number(h) / 10) === dc) st[h] = areas[h]; });
      mats = mats.concat(report(DNAME[dc], sm.filter(function (u) { return u.d === dc; }), st));
    });
    // Two-class version (Boro vs non-Boro), as used by the Boro-intensity panel
    var sm2 = sm.map(function (u) { return {h: u.h, m: u.m === 1 ? 1 : 2, r: u.r === 1 ? 1 : 2, d: u.d}; });
    var o2 = estimate(sm2, areas, [1, 2]);
    print('Two-class (Boro vs not Boro): OA ' + pct(o2.OA.est) + '% ± ' + pct(o2.OA.ci) +
          '; Boro UA ' + pct(o2.cls[0].UA.est) + ' ± ' + pct(o2.cls[0].UA.ci) +
          '; Boro PA ' + pct(o2.cls[0].PA.est) + ' ± ' + pct(o2.cls[0].PA.ci));

    Export.table.toDrive({collection: ee.FeatureCollection(rows),
      description: 'Boro_accuracy_summary_' + YEAR, folder: FOLDER, fileFormat: 'CSV'});
    Export.table.toDrive({collection: ee.FeatureCollection(mats),
      description: 'Boro_accuracy_error_matrix_' + YEAR, folder: FOLDER, fileFormat: 'CSV'});
    print('Run the two export tasks; send both CSVs to fill Supplementary Table 2.');
  });
}
