% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/code.R
\name{calibrationGlabTop2}
\alias{calibrationGlabTop2}
\title{calibration of glabTop2 model}
\usage{
calibrationGlabTop2(g.thickness, cs)
}
\arguments{
\item{g.thickness}{[object of class thickness] see \link[glabTop2]{thicknessGlabTop2} for details}

\item{cs}{cross-sections of the glacier, whose thickness along the central line of the glacier is to be used for calibration  (Formal class SpatialPolygons). Must have same projection as the DEM}
}
\value{
Calibrated Shape Factor for the glacier.
}
\description{
calibration of shape factor 'f' parameter for glabTop2 model (based on the methodology proposed by Ramsankaran et al., 2018)
}
\examples{
#creating a sample DEM
proj<-"+proj=utm +zone=44N +datum=WGS84 +units=m +no_defs "
values<-round(runif(2793,4000,6000),digits = 2)
r.extent<-raster::extent(308004.8, 336707.3, 3399079, 3423596)
r<-raster::raster(ext=r.extent,resolution=500, crs=proj,vals=values)

#creating a sample glacier boundary
shp <- as(raster::extent(309004.8, 334707.3, 3400579, 3420596), "SpatialPolygons")
sp::proj4string(shp) <- "+proj=utm +zone=44N +datum=WGS84 +units=m +no_defs"

#estimating glacier thickness
thickness<-thicknessGlabTop2(DEM=r,glacier.outline=shp,f.range=c(0.65,0.66),n=1,ncores=2)

#cross-section
lines.shp <- as(raster::extent(309004.8, 334707.3, 3410000, 3410000), "SpatialPolygons")
sp::proj4string(lines.shp) <- "+proj=utm +zone=44N +datum=WGS84 +units=m +no_defs"
lines = as(lines.shp, "SpatialLines")

#calibration of model
calibrationGlabTop2(g.thickness=thickness,cs=lines)

}
\references{
RAAJ Ramsankaran, Ankur Pandit & Mohd Farooq Azam (2018) Spatially distributed ice-thickness modelling for Chhota Shigri Glacier in western Himalayas, India, International Journal of Remote Sensing, 39:10, 3320-3343, DOI: 10.1080/01431161.2018.1441563
}
\author{
Aniruddha
}
