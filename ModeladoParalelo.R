library(raster)
library(parallel)
library(dismo)

LeeArchivoEspecie<- function(Archivo,directorio,datum)
{
  coordesp<- read.csv(paste(directorio,'/',Archivo,sep=''),header=T)
  coordinates(coordesp) <- c("Longitud", "Latitud")  # set spatial coordinates
  proj4string(coordesp) <- datum  # define projection system of our data
  return(coordesp)
}

ObtenerStackCapasBio<- function(Archivo, directorio)
{
   for (k in 1:4)
   {
      raster.capa<- raster(paste(directorio,'/',Archivo,'_',k,'.grd',sep=''))
      if (k ==1)
          stkcapas<- stack(raster.capa)
      else
         stkcapas<- stack(stkcapas,raster.capa)
   }
  return(stkcapas)
}

GuardaAsc<- function(capa.asc,Archivo,directorio)
{
  writeRaster(capa.asc,paste(directorio,Archivo,sep=''),format='ascii',overwrite=T)
}


genModelado<- function(j)
{
    ocurrencias<- LeeArchivoEspecie(archivo[j],paste(dircsv,carpeta,sep=''),datum)
   # ocurrenciaslval<- LeeArchivoEspecie(archivo[j],paste(dircsvv,carpeta,sep=''),datum)
    dirbio<- paste(dirbiorecor,carpeta,sep='')
    arch<- substr(archivo[j],1,nchar(archivo[j])-4)
    if (file.exists(paste(dirbio,'/',arch,'_',1,'.grd',sep='')) && length(ocurrencias)>10)
    {
      predictores<- ObtenerStackCapasBio(arch,dirbio)
      
      #Eliminamos los valores NA que son los puntos que no estan en las capas
      val <- extract(predictores, ocurrencias)
      indNA = which(is.na(val)==T)
      if (length(indNA) > 0)
           ocurrencias <- ocurrencias[-indNA, ]		   
    
      dyn.load('/Library/Java/JavaVirtualMachines/jdk1.8.0_131.jdk/Contents/Home/jre/lib/server/libjvm.dylib')
      #Modelo Maxend
      ### Maxent with noclamping and no extrapolation. 
      mxent <- maxent(predictores, ocurrencias, args=c("nowarnings",
                                              "noprefixes",
                                              "responsecurves",
                                              "outputformat=cumulative",
                                              "randomseed",
                                              "nowarnings", 
                                              "nowriteclampgrid",
                                              "replicates=5",
                                              "randomtestpoints=50",
                                              "replicatetype=bootstrap",
                                              "nodoclamp",
                                              "noextrapolate"))
      
      #Prediccion
      pmx <- predict(mxent, predictores, args=c("outputformat=cumulative")) 
      ModProm = calc(pmx, mean)
      ModMedian=calc(pmx,median)
      ModStd = calc(pmx, sd)
      GuardaAsc(ModProm,paste(arch,'_maxent_prom',sep=''),dirguardar)
      GuardaAsc(ModMedian,paste(arch,'_maxent_median',sep=''),dirguardar)
      GuardaAsc(ModStd,paste(arch,'_maxent_std',sep=''),dirguardar)
    }
}


dirbiorecor<- '~/Desktop/Proyecto/Especies/ResultCapas/';
dircsv<- '~/Desktop/Proyecto/Especies/PuntosSA/';
dirsave<- '~/Desktop/Proyecto/Especies/Modelado/'
datum <- CRS('+proj=longlat +ellps=WGS84 +datum=WGS84')  # geographical, datum WGS84
carpcsv<- list.files(dircsv)


#Parte paralela
num.cores<- detectCores()-1;
cl<- makeCluster(num.cores)
clusterEvalQ(cl, library(raster))
clusterEvalQ(cl, library(dismo))


#lee datos de ocurrencias length(carpcsv)
for (i in 1:length(carpcsv))
{
   carpeta<- carpcsv[i]
   archivo<- list.files(paste(dircsv,carpeta,sep=''))
   dirguardar<- paste(dirsave,carpeta,'/', sep='')
   if (dir.exists(dirguardar)==F)
     dir.create(dirguardar)
   
   clusterExport(cl,c('dircsv','dirbiorecor','datum','carpeta','archivo','dirguardar','LeeArchivoEspecie','ObtenerStackCapasBio','GuardaAsc'))
   parLapply(cl,1:length(archivo),genModelado)
}
stopCluster(cl)
