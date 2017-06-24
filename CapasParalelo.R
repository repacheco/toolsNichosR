library(raster)
library(rgdal)
library(maptools)
library(parallel)

rm(list = ls())

LeeArchivoShape<- function(Archivo,directorio,datum)
{
  ame.map<- readShapeSpatial(paste(directorio,Archivo,sep=''))
  proj4string(ame.map)<- datum
  return(ame.map)
}

LeeRecortaCapa<- function(especie.shape,directorio,datum)
{
  capa.raster<- raster(directorio)
  proj4string(capa.raster)<- datum
  extespecie<- as(extent(especie.shape),'SpatialPolygons')
  capa.especie<- crop(capa.raster, extespecie, byid=TRUE)
  capa.especie<- mask(capa.especie,especie.shape)
  return(capa.especie)
}

GuardaCapaRaster<- function(capa.especie,Archivo,directorio)
{
     writeRaster(capa.especie,paste(directorio,Archivo,sep=''),format='raster',overwrite=T)
     writeRaster(capa.especie,paste(directorio,Archivo,sep=''),format='ascii',overwrite=T)
}

GuardaAsc<- function(capa.asc,Archivo,directorio)
{
   writePolyShape(capa.asc, paste(directorio,Archivo,sep=''))
}


findArchivosShape<- function(k)
{
      arch<- carpeta[k]
      ext<- substr(arch,nchar(arch)-2,nchar(arch))
      if (ext=='shp')
      {
        Nom<- substr(arch,1,nchar(arch)-4)
        especie.shape<- LeeArchivoShape(arch,paste(dirshpres,dircarpetas[i],'/',sep=''),datum)
        GuardaAsc(especie.shape,Nom,dirguardar)
        for (j in 1:19)
        {
          capa.especie<- LeeRecortaCapa(especie.shape,paste(dircapas,j,'/','bio',j,'.asc',sep=''),datum)
          GuardaCapaRaster(capa.especie,paste(Nom,'_',j,sep=''),dirguardar)
          plot(capa.especie)
        }
      }
}

findCompPrin<- function(k)
{
  arch<- carpeta[k]
  ext<- substr(arch,nchar(arch)-2,nchar(arch))
  if (ext=='shp')
  {
    Nom<- substr(arch,1,nchar(arch)-4)
    especie.shape<- LeeArchivoShape(arch,paste(dirshpres,dircarpetas[i],'/',sep=''),datum)
    GuardaAsc(especie.shape,Nom,dirguardar)
    for (j in 1:4)
    {
      capa.especie<- LeeRecortaCapa(especie.shape,paste(dircapas,j,'.asc',sep=''),datum)
      GuardaCapaRaster(capa.especie,paste(Nom,'_',j,sep=''),dirguardar)
      #plot(capa.especie)
    }
  }
}

#ciclo para cada especie
datum <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")  # geographical, datum WGS84
dirshpres<- '~/Desktop/Proyecto/Especies/ResultM/'
dirsave<- '~/Desktop/Proyecto/Especies/ResultCapas/'
dircapas<- '~/Desktop/Proyecto/Especies/Componentes/Comp'

dircarpetas<- list.files(dirshpres)
numdircar<- length(dircarpetas)

#Parte paralela
num.cores<- detectCores()-1;
cl<- makeCluster(num.cores)
clusterEvalQ(cl, library(raster))
clusterEvalQ(cl, library(maptools))
clusterEvalQ(cl, library(rgdal))

for (i in 1:numdircar)
{
   carpeta<- list.files(paste(dirshpres,dircarpetas[i],sep=''))
   numcar<- length(carpeta)
   dirguardar<- paste(dirsave,dircarpetas[i],'/', sep='')
   if (dir.exists(dirguardar)==F)
      dir.create(dirguardar)
   
   clusterExport(cl,c('dirshpres','dircarpetas','dircapas','datum','dirguardar','i','carpeta','LeeArchivoShape','GuardaAsc','GuardaCapaRaster','LeeRecortaCapa'))
   parLapply(cl,1:numcar,findCompPrin)
   
}
stopCluster(cl)