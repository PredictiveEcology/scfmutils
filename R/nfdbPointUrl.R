## The National Fire Database point archive (all fires, shapefile) on the CFS server.
##
## The file was `NFDB_point.zip`; CFS renamed it `NFDB_point_shp.zip` (with `_txt` and `large_fires`
## variants beside it), and the old name returns HTTP 404, so getFirePoints_NFDB_scfm() could not
## download a release newer than the copy already on disk.
nfdbPointUrl <- function() {
  "https://cwfis.cfs.nrcan.gc.ca/downloads/nfdb/fire_pnt/current_version/NFDB_point_shp.zip"
}
