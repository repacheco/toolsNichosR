library(maptools)
library(rgdal)
library(spdep)
library(parallel)
library(sp)
library(raster)
library(mgcv)
library(data.table)

LeeArchivoEspecie<- function(directorio,Archivo,datum)
{
  coordesp<- read.csv(paste(directorio,'/',Archivo,sep=''),header=T)
  coordinates(coordesp) <- c("Longitud", "Latitud")  # set spatial coordinates
  proj4string(coordesp) <- datum  # define projection system of our data
  return(coordesp)
}


LeeArchivoAsc<- function(Directorio,Archivo)
{
  capa.raster<- raster(paste(Directorio,Archivo,sep=''))
  return(capa.raster)
}

DistanciaCentroide<- function(mapa.binario,stack.componentes,metodo='Euclideana')
{
  #cortamos al mapa adecuado
  mapa.binario<- reclassify(mapa.binario,c(-Inf,0,NA))
  stack.cortados<- crop(stack.componentes,mapa.binario@extent)
  
  # Crear una máscara de las coberturas con base en el mapa binario
  msk = mask(stack.cortados, mapa.binario)
  # Calcula la media de cada cobertura 
  media = cellStats(msk,mean)
  
  if (metodo=='Euclideana')
  {
    # Calcula la diferencia entre cada pixel y la media. Luego la eleva al cuadrado
    dif = (msk - media)^2
    # Genera la raíz cuadrada de la sumatoria de las diferencias. 
    capa.dc = sqrt(sum(dif))
  }
  else
  {
    matriz.covarianza=layerStats(msk,'cov',na.rm=T)
    comp1<- msk$Comp1@data@values
    comp1<- comp1[!is.na(comp1)]
    
    comp2<- msk$Comp2@data@values
    comp2<- comp2[!is.na(comp2)]
    
    comp3<- msk$Comp3@data@values
    comp3<- comp3[!is.na(comp3)]
    
    comp4<- msk$Comp4@data@values
    comp4<- comp4[!is.na(comp4)]
    
    matriz<- cbind(comp1,comp2,comp3,comp4)
    
    valores.mahalanobis<- mahalanobis(matriz,matriz.covarianza$mean,matriz.covarianza$covariance)
    msk$Comp1@data@values[!is.na(msk$Comp1@data@values)]<- valores.mahalanobis
    capa.dc<- msk$Comp1
    
  }
  return(capa.dc)
}

GeneraStackComponentes<- function(dirComponentes)
{
  for (l in 1:4)
  {
    comp<- raster(paste(dirComponentes,l,'.asc',sep=''))
    if (l==1)
      stack.componentes<- stack(comp)
    else
      stack.componentes<- stack(stack.componentes,comp)
  }
  return(stack.componentes)
}

MapaBinario<- function (mapa.maxent,ocurrencias,Threshold)
{
  favorabilidad<- extract(mapa.maxent,ocurrencias,na.rm=T)
  favorabilidad.ordenada<- sort(favorabilidad)
  threshold.favorabilidad<- favorabilidad.ordenada[round((length(favorabilidad)*Threshold))]
  mapa.binario<- mapa.maxent>threshold.favorabilidad
  return(mapa.binario)
}

GuardaAsc<- function(capa,Archivo,directorio)
{
  writeRaster(capa,paste(directorio,Archivo,sep=''),format='ascii',overwrite=T)
}

genAnalisisFinales<- function(j)
{
  arch<- carpeta[j]
  NomArch<- substr(arch,1,nchar(arch)-4)
  dirarchmodelo<- paste(NomArch,'_maxent_median.asc', sep='')
  
  if (file.exists(paste(dirModelado,dircarpetas[i],'/',dirarchmodelo,sep=''))==T)
  {    
    ocurrencias<- LeeArchivoEspecie(directorioCA,arch,datum)
    ocurrenciasSA<- LeeArchivoEspecie(directorioSA,arch,datum)
    mapa.maxent<- LeeArchivoAsc(paste(dirModelado,dircarpetas[i],'/',sep=''),dirarchmodelo)
    
    #Creamos el mapa binario 
    mapa.binario<- MapaBinario(mapa.maxent,ocurrencias,0.05)
    GuardaAsc(mapa.binario,paste(NomArch,'_binario',sep=''),dirguardar)
    
    #Hallamos la distancia centroide
    distancia.centroide<- DistanciaCentroide(mapa.binario,stack.componentes,'Mahalanobis')
    GuardaAsc(distancia.centroide,paste(NomArch,'_centroide',sep=''),dirguardar)
    
    #Correlacion 
    favorabilidad<- extract(mapa.maxent,ocurrencias,na.rm=T)
    pruebacor<- cor.test(favorabilidad,ocurrencias$Abundancias)
    correlacion<- pruebacor$estimate
    ajustelineal<- lm(ocurrencias$Abundancias~favorabilidad)
    media<- mean(ocurrencias$Abundancias,na.rm=T)
    r2<- sum((ajustelineal$fitted.values-media)^2)/sum((ocurrencias$Abundancias-media)^2)
    #summary(ajustelineal)
    
    favorabilidad.dc<- extract(distancia.centroide,ocurrencias,na.rm=T)
    pruebacor.dc<- cor.test(favorabilidad.dc,ocurrencias$Abundancias)
    correlacion.dc<- pruebacor.dc$estimate
    ajustelineal.dc<- lm(ocurrencias$Abundancias~favorabilidad.dc)
    media.dc<- mean(ocurrencias$Abundancias,na.rm=T)
    r2.dc<- sum((ajustelineal.dc$fitted.values-media.dc)^2)/sum((ocurrencias$Abundancias-media.dc)^2)
    #summary(ajustelineal.dc)
    
    #Grafica de los puntos
    png(paste(dirguardar,NomArch,'_favorabilidad.png',sep=''))
    plot(favorabilidad,ocurrencias$Abundancias,type='p')
    abline(ajustelineal,col='red')
    dev.off()
    
    png(paste(dirguardar,NomArch,'_favorabilidadCentroide.png',sep=''))
    plot(favorabilidad.dc,ocurrencias$Abundancias)
    abline(ajustelineal.dc,col='red')
    dev.off()
    
    #ROC Parcial
    roc<- iPartialROC(paste(directorioCA,'/',arch,sep=''), paste(dirModelado,dircarpetas[i],'/',dirarchmodelo,sep=''), 0.9, 50, 100, paste(dirguardar,NomArch,"_TestRoc.txt",sep=''))
    mean.ratio<- mean(roc$AUC_ratio)
    ratiom1<- length(which(roc$AUC_ratio<1))
    pvalor<- ratiom1/length(roc$AUC_ratio)
    
    #GAMS
    gam.fav<- gam(ocurrencias$Abundancias~favorabilidad, family = poisson(link = "log"))
    summary(gam.fav)
    devianza.fav<- ((gam.fav$null.deviance-gam.fav$deviance)/gam.fav$null.deviance)
    gam.r2<- sum((gam.fav$fitted.values-media.dc)^2)/sum((ocurrencias$Abundancias-media.dc)^2)
    png(paste(dirguardar,NomArch,'_GAM.png',sep=''))
    plot(favorabilidad,ocurrencias$Abundancias)
    favorabilidad<- favorabilidad[!is.na(favorabilidad)]
    lines(favorabilidad,gam.fav$fitted.values,col='red')
    dev.off()
    
    gam.dc<- gam(ocurrencias$Abundancias~favorabilidad.dc, family = poisson(link = "log"))
    summary(gam.dc)
    devianza.dc<- ((gam.dc$null.deviance-gam.dc$deviance)/gam.dc$null.deviance)
    gam.r2.dc<- sum((gam.dc$fitted.values-media.dc)^2)/sum((ocurrencias$Abundancias-media.dc)^2)
    png(paste(dirguardar,NomArch,'_GAMDC.png',sep=''))
    plot(favorabilidad.dc,ocurrencias$Abundancias)
    favorabilidad.dc<- favorabilidad.dc[!is.na(favorabilidad.dc)]
    lines(favorabilidad.dc,gam.dc$fitted.values,col='red')
    dev.off()
    
    #Tabla
    total.abundancias<- length(ocurrencias$Abundancias)
    total.modelado<- length(ocurrenciasSA$X)
    total=total.abundancias+total.modelado
    
    newrow<- data.frame(Especie=NomArch,Total=total,DatosAbundancias=total.abundancias,DatosModelado=total.modelado,Mean_ratio=mean.ratio,Pvalor=pvalor,Ratiomenor1=ratiom1,Correlacion=correlacion,R2.lineal=r2,CorrelacionDC=correlacion.dc,R2.lineal.DC=r2.dc,R2.GAM=gam.r2,R2.GAM.DC=gam.r2.dc,Devianza=devianza.fav,Devianza.DC=devianza.dc)
  }
  else
     newrow<- data.frame(Especie='',Total=0,DatosAbundancias=0,DatosModelado=0,Mean_ratio=0,Pvalor=0,Ratiomenor1=0,Correlacion=0,R2.lineal=0,CorrelacionDC=0,R2.lineal.DC=0,R2.GAM=0,R2.GAM.DC=0,Devianza=0,Devianza.DC=0)
  return(newrow)
}



dirpuntosCA<- '~/Desktop/Proyecto/Especies/PuntosCA/'
dirpuntosSA<- '~/Desktop/Proyecto/Especies/PuntosSA/'
dirsave<- '~/Desktop/Proyecto/Especies/AnalisisFinales/'
dirModelado<- '~/Desktop/Proyecto/Especies/ModeladoCumulative/'
dirComponentes<- '~/Desktop/Proyecto/Especies/Componentes/Comp'
datum <- CRS('+proj=longlat +ellps=WGS84 +datum=WGS84')  # geographical, datum WGS84

dircarpetas<- list.files(dirpuntosCA)
numdircar<- length(dircarpetas)
stack.componentes<- GeneraStackComponentes(dirComponentes)
source('~/Desktop/Proyecto/Especies/PartialROC.r')
tabla<- data.frame(Especie='',Total=0,DatosAbundancias=0,DatosModelado=0,Mean_ratio=0,Pvalor=0,Ratiomenor1=0,Correlacion=0,R2.lineal=0,CorrelacionDC=0,R2.lineal.DC=0,R2.GAM=0,R2.GAM.DC=0,Devianza=0,Devianza.DC=0)

#Parte paralela
num.cores<- detectCores()-1;
cl<- makeCluster(num.cores)
clusterEvalQ(cl, library(raster))
clusterEvalQ(cl, library(maptools))
clusterEvalQ(cl, library(rgdal))
clusterEvalQ(cl, library(spdep))
clusterEvalQ(cl, library(sp))
clusterEvalQ(cl, library(mgcv))
clusterEvalQ(cl, source('~/Desktop/Proyecto/Especies/PartialROC.r'))


for (i in 1:numdircar)
{
  carpeta<- list.files(paste(dirpuntosCA,dircarpetas[i], sep=''))
  dirguardar<- paste(dirsave,dircarpetas[i],'/', sep='')
  directorioCA<- paste(dirpuntosCA,dircarpetas[i],sep='')
  directorioSA<- paste(dirpuntosSA,dircarpetas[i],sep='')
  if (dir.exists(dirguardar)==F)
    dir.create(dirguardar)
  
  clusterExport(cl,c('i','dircarpetas','directorioSA','directorioCA','datum','carpeta','dirModelado','dirComponentes','stack.componentes','dirguardar','LeeArchivoEspecie','LeeArchivoAsc','DistanciaCentroide','MapaBinario','GuardaAsc'))
  newrow<- parLapply(cl,1:length(carpeta),genAnalisisFinales)
  nw<- rbindlist(newrow)
  tabla<- rbind(tabla,nw)
}

indices<- which(tabla$Especie=='')
tabla<- tabla[-indices,]
write.csv(tabla,paste(dirsave,'Tabla de Resultados.csv',sep=''))
stopCluster(cl)


