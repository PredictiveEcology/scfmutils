# scfmutils 2.0.7

- scfmDriver `calibrateFireRegimePolys()` gracefully exits if polygons have no fires;

# scfmutils 2.0.6

- fix attempt to use fire regime attributes when not using FRT/FRU polygons;

# scfmutils 2.0.5

- add Erni et al. (2020) fire regime polygon attributes to FRT and FRU polygons;

# scfmutils 2.0.4

- add `withr` to Suggests (used with testing);
- fix plotting to handle cases with no escaped fires;
- more robust deslivering;

# scfmutils 2.0.3

- require `reproducible >= 2.1.0`;
- updates for recent changes to `reproducible`: pass `to` arguments to `reproducible::prepInputs()` and related functions, instead of passing `studyArea` or `rasterToMatch`;
- new option `scfmutils.driver.plot.scam` (default `TRUE`) to allow user to disable creation of scam plots in `calibrateFireRegimePolys()`;

# scfmutils 2.0.2

- warn instead of error if spread model calibration fails in `calibrateFireRegimePolys()`;
- improve documentation;

# scfmutils 2.0.1

- add `fpCompare` to Imports, for improved floating point number comparisons;
- bump minimum versions requirements for dependency packages `SpaDES.core`, `SpaDES.tools`, and `reproducible`;
- fix issue with per-polygon flammable cell calculations in  `.makeLandscapeAttr()` that was introduced by switching to `sf` for fire regime polygons;
- modified `comparePredictions_fireDistribution()` to use arguments `fireRegimePoints`, `burnSummary`, and `size` instead of `dt`;
- improve CRS checks in `checkForIssues()`;
- documentation improvements;

# scfmutils 2.0.0

- fully implement use of `sf` polygons object instead of lists of polygon attributes;
- bump minimum versions of `Require`, `reproducible`, `SpaDES.core`, `SpaDES.tools`;
- new ggplot wrappers for plotting scfm rasters;
- documentation improvements;

# scfmutils v1.0.0

**Transitional version: introduces several breaking changes**

- transitioning to using `sf` polygons object instead of lists of polygon attributes;
- remove arguments `cellSize` and `neighbours` from `calibrateFireRegimePolys()`;
- `calcZonalRegimePars()` now returns a `data.table` instead of a list;
- fix issue with `prepInputsFireRegimePolys()` when url supplied by user;
- warn instead of error if scam plots can't be created;

# scfmutils 0.0.13

- drop support for R 4.1;
- remove `LandR` dependency;
- improved documentation;
- improved messaging;

# scfmutils 0.0.12

- add package `bcdata` to Suggests to allow use of BEC zones etc. `prepInputsfireRegimePolys()`;

# scfmutils 0.0.11

- fix for 0 neighbour instance in `makeLandscapeAttr()`;

# scfmutils 0.0.10

- add package `googledrive` to Suggests;
- increase the minimum `SpaDES.core` and `SpaDES.tools` versions following `terra` migration and other changes;
- update `prepInputsFireRegimePolys()` example to use `terra`;
- documentation improvements;

# scfmutils 0.0.9
- transition to using `terra` instead of `raster` package;
- remove `fasterize` dependency;
- use R native pipe instead of `magrittr`'s;
- drop R 4.0 support (reproducible requires R >= 4.1)
- add plotting functions for escapes and historical fire distributions;
- improved diagnostic plot calculations using `studyAreaReporting`;
- internal package restructuring;

# scfmutils 0.0.8

- use `%>%` for pipe (R < v4.2);
- add Erni et al. (2020) FRTs + FRUs to `prepInputsFireRegimePolys()`;

# scfmutils 0.0.7

- corrections calculating stats for diagnostics;

# scfmutils 0.0.6

- `prepInputsFireRegimePolys()` from BEC NDTs: allow either col name;
- ensure `plotPath` exists for scam plots;
- fix bug in `comparePredictions_summaryDT()`;
- fix bug in `deSliver()` and clarify documentation;
- pass 'optimizer' argument to `scam()` in `calibrateFireRegimePolys()`;
- additional messaging for driver calibration;
- save scam plots per poly during driver calibration;

# scfmutils 0.0.5

- add `scfmRegime` functions;
- add functions `checkForIssues()` and `deSliver()` for working with fire regime polygons;
- add BECNDT to `prepInputsFireRegimePolys()`;

# scfmutils 0.0.4

- increase minimum `SpaDES.core` version requirement;
- rename `PolyId` to `PolyID`;
- fix bug in `executeDesign()`;
- improved messaging;

# scfmutils 0.0.3

- use `theme_bw()` instead of `theme_minimal()` in diagnostic plots;

# scfmutils 0.0.2

- export `fireRegimePolyTypes()`;
- export `comparePredictions_annualIgnitions()`;

# scfmutils 0.0.1

- requires R >= 4.0;
- add `times` argument to `comparePredictions_fireReturnInterval()`;
- pass only what's used from `simList` to `comparePredictions_summaryDT()`;
- add note re: output crs using `prepInputsFireRegimePolys()`;
- pass `startTime` to `.executeDesignInternal()`;
- fix `.makeLandscapeAttr()`;
- add functions from `scfmLandCoverInit`;
- fix use of `inherits`;
- add `prepInputsFireRegimePolys()`;
- use `pkgdown` for website;
