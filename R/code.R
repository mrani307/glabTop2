
#' glabTop2 model for glacier thickness estimation
#' @description  Estimation of glacier thickness, based on the Glacier Bed Topography (GlabTop2) model (Frey et al., 2014)
#' @export
#'
#' @param DEM glacier DEM (Formal class RasterLayer).Must be projected
#' @param glacier.outline glacier boundary shapefile (Formal class SpatialPolygons).Must have same projection as the DEM
#' @param g gravitational acceleration in ms-2 [default 9.8]
#' @param rho density of ice in kgm-3[default 900]
#' @param h.min minimum elevation difference(in m) in the buffer formed around randomly selected cell, for calculation of mean surface slope[default 50]
#' @param rc fraction of randomly selected cells from glacier inner cells to be used for calculation of thickness[default 0.3]
#' @param f.range range of values of shape factors 'f', for which thickness is to be estimated [default seq(0.6,0.9,0.01)]
#' @param n number of thickness estimation for a particular value of shape factor 'f' [default n=3]
#' @param idw.P power of the IDW(inverse distant weighted) Interpolation [default 3]
#' @param hga value assigned to glacier adjacent cells [default 15]
#' @param ncores number of cores of your CPU, to be used for running the code [default: total number of codes of CPU-3]
#'
#' @importFrom raster raster extent res crs mask as.matrix terrain sampleRandom
#' @importFrom sp coordinates gridded spTransform
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach getDoParRegistered registerDoSEQ foreach
#' @importFrom data.table fsetdiff as.data.table
#' @importFrom stats na.omit
#'
#' @returns  Returns a list of two objects -
#'           (i)'averageThickness'  :  (Formal Class RasterLayer) average of the thickness estimated, using the selected values of shape factor 'f'
#'           (ii)'thickness.shapeFactors'  : (List of Formal class RasterLayer)thickness estimated for each
#' @examples #creating a sample glacier DEM
#' proj<-"+proj=utm +zone=44N +datum=WGS84 +units=m +no_defs "
#' values<-round(runif(2793,4000,6000),digits = 2)
#' r.extent<-raster::extent(308004.8, 336707.3, 3399079, 3423596)
#' r<-raster::raster(ext=r.extent,resolution=500, crs=proj,vals=values)
#' #creating a sample glacier boundary
#' shp <- as(raster::extent(309004.8, 334707.3, 3400579, 3420596), "SpatialPolygons")
#' sp::proj4string(shp) <- "+proj=utm +zone=44N +datum=WGS84 +units=m +no_defs"
#'
#' #estimating glacier thickness
#' thickness<-thicknessGlabTop2(DEM=r,glacier.outline=shp,f.range=c(0.65,0.66),n=1,ncores=2)
#'
#' @author Aniruddha
#' @references RAAJ Ramsankaran, Ankur Pandit & Mohd Farooq Azam (2018) Spatially distributed ice-thickness modelling for Chhota Shigri Glacier in western Himalayas, India, International Journal of Remote Sensing, 39:10, 3320-3343, DOI: 10.1080/01431161.2018.1441563
#' @references Frey, H., Machguth, H., Huss, M., Huggel, C., Bajracharya, S., Bolch, T., Kulkarni, A., Linsbauer, A., Salzmann, N., and Stoffel, M.: Estimating the volume of glaciers in the Himalayan–Karakoram region using different methods, The Cryosphere, 8, 2313–2333, https://doi.org/10.5194/tc-8-2313-2014, 2014.

thicknessGlabTop2<-function(DEM,glacier.outline,g,rho,h.min,rc,f.range,n,idw.P,hga,ncores)
{
  ###########################################################################################################
  #setting default values 
  ###########################################################################################################
  
    if(missing(g)) g<-9.8
    if(missing(rho)) rho<-900 #density of ice in kgm-3
    if(missing(h.min)) h.min<-50 #maximum elevation difference in the buffer
    if(missing(f.range)) f.range<-seq(0.6,0.9,0.01)  #range of f for which it thickness is to be estimated
    if(missing(n)) n<-3        #number of iteration for each 'f'
    if(missing(hga)) hga<-15
    if(missing(rc)) rc<-0.3
    if(missing(idw.P)) idw.P<-3
    {
      cores<-parallel::detectCores()
      if(missing(ncores)) n.cores<-cores-3
      else
      {
        if(ncores>(cores-1)) n.cores<-cores-3
        else  n.cores<-ncores
      }
    }
  }
  
  #########################################################################################################
  #functions
  #########################################################################################################
  
    #making nxn buffer
    buffer<-function(r,c,size,data.mat){
      if(r-size<=0 | c-size<=0)
        return(F)
      else{
        mat<-matrix(nrow = (2*size+1), ncol = (2*size+1))
        count.i<-count.j<-1
        for( r1 in (r-size):(r+size)){
          count.j<-1
          for(c1 in (c-size):(c+size)){
            mat[count.i,count.j]<-data.mat[r1,c1]
            count.j<-count.j+1
          }
          count.i<- count.i+1
        }
        return(mat)
      }
    }
    #make raster from matrix data
    makeRaster<-function(ras=g.DEM,ras.data){
      r<-raster::raster(ext=raster::extent(ras),resolution= raster::res(ras),
                crs=as.character(raster::crs(ras)),vals=ras.data)
      return(r)
    }
    #make raster from coordinates data
    makeRaster2<-function(r){
      sp::coordinates(r)<- ~ x+y
      sp::gridded(r)<- TRUE
      raster::crs(r)<-raster::crs(g.DEM)
      r<-raster::raster(r)
      return(r)
    }
    #identify the boundary cells
    checkBoundaryCell<-function(x){
      r<-as.numeric(g.DEM.rowcol[x,1])
      c<-as.numeric(g.DEM.rowcol[x,2])
      if(r==1 | c==1 | r==nrow(g.DEM.data) | c==ncol(g.DEM.data))
        return(NA)
      else{
        mat<-buffer(r=r,c=c,size=1,data.mat = g.DEM.data)
        if(is.na(g.DEM.data[r,c])==F & all(is.na(mat))==F & any(is.na(mat))==T) return(0) else return(NA)
      }
    }
    #identify each type of cells
    checkAllCell<-function(x){
      r<-as.numeric(g.DEM.rowcol[x,1])
      c<-as.numeric(g.DEM.rowcol[x,2])
      if(r==1 | c==1 | r==nrow(g.DEM.data) | c==ncol(g.DEM.data))
        return(NA)
      else{
        mat<-buffer(r=r,c=c,size=1, data.mat = g.DEM.data.boundarycells)
        if(is.na(g.DEM.data.boundarycells[r,c])==F & g.DEM.data.boundarycells[r,c]==0) #glacier boundary cells marked as zero
          return(0)
        else if(is.na(g.DEM.data[r,c])==F) #glacier inner cells given their DEM heights
          return(g.DEM.data[r,c])
        else if(all(is.na(mat))==F & is.na(g.DEM.data[r,c]==T)){
          if(any(mat==0)==T)
            return(15)
        }
        else
          return(NA)
      }
    }
    #idw interpolation
    idw<-function(x,c){
      w<-((as.numeric(empty.inner.cells[x,3])-known.data[,1])^2
          +(as.numeric(empty.inner.cells[x,4])-known.data[,2])^2)^(-idw.P/2)
      return(sum(w*known.data[,3])/sum(w))
    }
  }

######################################################################################################
#model inputs preperation
######################################################################################################
  
    #checking if the inouts have same projection, else correct
      if(as.character(raster::crs(DEM))=="+proj=longlat +datum=WGS84 +no_defs")
        stop("DEM must be projected..")
      if(as.character(raster::crs(DEM))!=as.character(raster::crs(glacier.outline)))
        stop("glacier.outline must have same projection as the DEM..")

    #glacier DEM
      g.DEM<-raster::mask(DEM,glacier.outline)
      g.DEM.data<-raster::as.matrix(g.DEM)

    #glacier slope
      glacier.slope<-raster::terrain(DEM, opt="slope", unit="degrees", neighbours=8)
      glacier.slope<-raster::mask(glacier.slope,glacier.outline)
      g.slope.data<-raster::as.matrix(glacier.slope)
  



#######################################################################################################
#datas and parameters
#######################################################################################################
{
  del.h<-(max(g.DEM.data, na.rm= T)-min(g.DEM.data, na.rm= T))/1000 #in kms
  if(del.h<=1.6)
    b.sh.stress<-0.5+159.8*del.h-43.5*del.h^2 # in kPa
  else
    b.sh.stress<-150
}

#######################################################################################################
#ONE TIME CODE
#######################################################################################################

  #row col data for all the masked DEM data
    g.DEM.rowcol<-raster::sampleRandom(g.DEM, size = length((as.vector(g.DEM.data))), na.rm=F, rowcol=T)

  #identifying the glacier boundary cells
    g.DEM.data.boundarycells<-matrix(sapply(c(1:(nrow(g.DEM.data)*ncol(g.DEM.data))), FUN = checkBoundaryCell), nrow = nrow(g.DEM.data), ncol = ncol(g.DEM.data), byrow = T)

  #identifying all types of cells
    g.DEM.data.allcells<-matrix(sapply(c(1:(nrow(g.DEM.data)*ncol(g.DEM.data))), FUN = checkAllCell), nrow = nrow(g.DEM.data), ncol = ncol(g.DEM.data), byrow = T)
    g.DEM.allcells<-makeRaster(ras.data = g.DEM.data.allcells)
    g.DEM.allcells[g.DEM.allcells==15]<-hga

  #each type of cells
    g.DEM.data.cells<-list()
    g.DEM.cells.rowcol<-list()
    #1:boundary,2:adjacent,3:inner,4:non-glacier
    for(i in 1:4){
      store<-g.DEM.data.allcells
      if(i==1)store<-g.DEM.data.boundarycells
      if(i==2)store[store!=15]<-NA
      if(i==3)store[store==0 | store==15]<-NA
      if(i==4){
        store[is.na(store)==T]<--6
        store[store!=(-6)]<-NA
      }
      g.DEM.data.cells[[i]]<-store
      store.raster<-makeRaster(ras.data = g.DEM.data.cells[[i]])
      g.DEM.cells.rowcol[[i]]<-raster::sampleRandom(store.raster, size = length(as.vector(store.raster)), na.rm=T,rowcol=T, xy=T)
    }

  #slope for all cells
    #register parallel code
      my.cluster <- parallel::makeCluster(n.cores,type = "PSOCK")
      print(my.cluster)
      doParallel::registerDoParallel(cl = my.cluster)
      if(foreach::getDoParRegistered()==F)
        stop("Restart R Studio")

  message("Calculating the mean surface slope of the randomly selected cells.....")
  s.r.c<-foreach::foreach(i=1:nrow( g.DEM.cells.rowcol[[3]]), .combine = 'c')  %dopar% {
    three=0
    size=0
    while(three<=h.min){
      size=size+1
      a<-buffer(r=g.DEM.cells.rowcol[[3]][i,1],c=g.DEM.cells.rowcol[[3]][i,2], size = size, data.mat = g.DEM.data)
      if(is.logical(a)==T){
        size<-NA
        break
      }
      else
        three<-max(a, na.rm = T)-min(a,na.rm = T)
    }
    if(is.na(size)==F)
      size<-mean(buffer(r=g.DEM.cells.rowcol[[3]][i,1],c=g.DEM.cells.rowcol[[3]][i,2], size = size, data.mat = g.slope.data),na.rm=T)
    size
  }
  message("Mean surface slope of the randomly selected cells.....DONE..")

    #end parallel computing
      foreach::registerDoSEQ()
      parallel::stopCluster(cl = my.cluster)

  data<-cbind(g.DEM.cells.rowcol[[3]],s.r.c)
  data.na<-na.omit(data)



#######################################################################################################
# glacial lake thickness estimations - iterations
#######################################################################################################


  #register parallel cores
    my.cluster <- parallel::makeCluster(n.cores,type = "PSOCK")
    print(my.cluster)
    doParallel::registerDoParallel(cl = my.cluster)
    if(foreach::getDoParRegistered()==F)
      stop("Restart R Studio")

 ######################################################################################################
 #thickness algorithms
 ######################################################################################################
  
  message("Calculating Glacier Thickness.....")
  sum<-0
  count<-1
  thickness.f<-list()
  
  #thickness estimation for a various shape-factor'f' values
  for(f in f.range){
    
    r.data<-list()
    sum.r.data<-0
    
    #number of iteration of thickness estimations for each shape factor values
    for(i in 1:n){

      print(paste0(f,"-",i,"     :: ",Sys.time()))

      #random cell determination
      rand.cells<-data[sample(nrow(data.na), as.integer(nrow(g.DEM.cells.rowcol[[3]])*rc)), ]
      empty.inner.cells<-data.table::fsetdiff(data.table::as.data.table(data), data.table::as.data.table(rand.cells))

      #thickness estimation of random cells
      thickess_random<-b.sh.stress/(f*rho*g*sin(rand.cells[,6]*pi/180))*1000
      rand.cells<-cbind(rand.cells,thickess_random)
      colnames( rand.cells)[7]<-"thickness"

      known.data<-na.omit(rbind(rand.cells[,c(3,4,7)],g.DEM.cells.rowcol[[2]][,c(3:5)]))

      #idw interpolation
      interpolated.heights<-foreach::foreach(j=1:nrow(empty.inner.cells), .combine = 'c') %dopar%{
        idw(x=j,c=i)
      }

      empty.inner.cells<-cbind(empty.inner.cells,interpolated.heights)
      colnames( empty.inner.cells)[7]<-"thickness"

      #making raster data from the thickness values
      emp<-g.DEM.cells.rowcol[[4]][,c(3,4,5)]
      emp[,3]<-NA
      colnames(emp)[3]<-"thickness"
      raster.data<-rbind(rand.cells[,c(3,4,7)],empty.inner.cells[,c(3,4,7)],emp)
      
      #thickness raster for particular iteration for a value of shape-factor
      r.data[[i]]<-makeRaster2(as.data.frame(raster.data))
       
      sum.r.data<-sum.r.data+r.data[[i]]
    }
    
    #thickness raster for a partuicular value of shape-factor (mean of all iterations)
    r.mean<-sum.r.data/n
    thickness.f[[as.character(f)]]<-r.mean
    count<-count+1
    sum<-sum+r.mean

  }
 

  #stop parallel
    foreach::registerDoSEQ()
    parallel::stopCluster(cl = my.cluster)
 
  #average thickness of glacier - mean of thickness calculated for all the shape-factors
  avg.thickness<-makeRaster(ras.data = as.matrix(sum/length(f.range)))

  message("Modelling Over.....")
  return(list(averageThickness=avg.thickness,thickness.shapeFactors=thickness.f))


#########################################################################################################
#MODELLING OVER
#########################################################################################################
}


#'calibration of glabTop2 model
#' @description calibration of shape factor 'f' parameter for glabTop2 model (based on the methodology proposed by Ramsankaran et al., 2018)
#' @export
#'
#' @param g.thickness [object of class thickness] see \link[glabTop2]{thicknessGlabTop2} for details
#' @param cs cross-sections of the glacier, whose thickness along the central line of the glacier is to be used for calibration  (Formal class SpatialPolygons). Must have same projection as the DEM
#'
#' @importFrom raster crs extract
#' @importFrom rgeos gLength
#' @importFrom maptools SpatialLinesMidPoints
#' @returns Calibrated Shape Factor for the glacier.
#'
#' @examples #creating a sample DEM
#' proj<-"+proj=utm +zone=44N +datum=WGS84 +units=m +no_defs "
#' values<-round(runif(2793,4000,6000),digits = 2)
#' r.extent<-raster::extent(308004.8, 336707.3, 3399079, 3423596)
#' r<-raster::raster(ext=r.extent,resolution=500, crs=proj,vals=values)
#'
#' #creating a sample glacier boundary
#' shp <- as(raster::extent(309004.8, 334707.3, 3400579, 3420596), "SpatialPolygons")
#' sp::proj4string(shp) <- "+proj=utm +zone=44N +datum=WGS84 +units=m +no_defs"
#'
#' #estimating glacier thickness
#' thickness<-thicknessGlabTop2(DEM=r,glacier.outline=shp,f.range=c(0.65,0.66),n=1,ncores=2)
#'
#' #cross-section
#' lines.shp <- as(raster::extent(309004.8, 334707.3, 3410000, 3410000), "SpatialPolygons")
#' sp::proj4string(lines.shp) <- "+proj=utm +zone=44N +datum=WGS84 +units=m +no_defs"
#' lines = as(lines.shp, "SpatialLines")
#'
#' #calibration of model
#' calibrationGlabTop2(g.thickness=thickness,cs=lines)
#'
#' @author Aniruddha
#' @references RAAJ Ramsankaran, Ankur Pandit & Mohd Farooq Azam (2018) Spatially distributed ice-thickness modelling for Chhota Shigri Glacier in western Himalayas, India, International Journal of Remote Sensing, 39:10, 3320-3343, DOI: 10.1080/01431161.2018.1441563
calibrationGlabTop2<-function(g.thickness,cs)
{
  thickness<-g.thickness[[1]]
  if(as.character(raster::crs(thickness))!=as.character(raster::crs(cs)))
    stop("cs must have same projection as the DEM..")
  else
  {
    cs.length<-rgeos::gLength(cs, byid=T)
    cs.mid<-maptools::SpatialLinesMidPoints(cs)
    cs.thickness<-raster::extract(thickness,cs.mid)

    f<-2/pi*atan(cs.length/(2*cs.thickness))
    calib.f<-round(mean(f),digits = 2)
    print(paste0("Calibrated Shape factor is ",calib.f))
  }
}


