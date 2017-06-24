library(maptools)
library(rgdal)
library(spdep)
library(sp)
library(raster)

LeeArchivoGeneral<- function(directorio,Archivo,datum)
{
  ame.map<- shapefile(paste(directorio,Archivo,sep=''))
  proj4string(ame.map)<- datum
  return(ame.map)
}

LeeCoordenadasEspaciales<- function(directorio,Archivo,datum)
{
  coordesp<- read.csv(paste(directorio,'/',Archivo,sep=''),header=T)
  coordinates(coordesp) <- c("Longitud", "Latitud")  # set spatial coordinates
  proj4string(coordesp) <- datum  # define projection system of our data
  return(coordesp)
}

#inicia el tiempo de procesamiento
ptm <- proc.time();

#Datos Generales
setwd('~/Desktop/Proyecto/Especies/')
datum <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")  # geographical, datum WGS84
dirshpgral<- '~/Desktop/Proyecto/Especies/wwfreg/'
ame.map.Original<- LeeArchivoGeneral(dirshpgral,"wwf_NAme",datum)

#poligonos adyacentes a los poligonos de los puntos
pt.ady.tot<- poly2nb(ame.map.Original)
pt.card<- card(pt.ady.tot)

#leer el directorio donde estan los archivos
dirgral<- '~/Desktop/Proyecto/Especies/PuntosSA/';
Dirsave<- '~/Desktop/Proyecto/Especies/ResultMs/';
dircarpetas<- list.files(dirgral)
numdircar<- length(dircarpetas)

Errores<- data.frame(Especie='',Puntos=0)

for (i in 1:numdircar)
{
  carpeta<- list.files(paste(dirgral,dircarpetas[i],sep=''))
  numcar<- length(carpeta)
  dirguardar<- paste(Dirsave,dircarpetas[i],'/', sep='')
  directorio<- paste(dirgral,dircarpetas[i],sep='')
  
  if (dir.exists(dirguardar)==F)
    dir.create(dirguardar)
  
  for (j in 1:numcar)
  {
    arch<- carpeta[j]
    coordesp<- LeeCoordenadasEspaciales(directorio,arch,datum)
    NomArc<- substr(arch,1,nchar(arch)-4)
    
    #determina los puntos sobre los poligonos
    ame.map<- ame.map.Original
    pt.in.poly <- unique(over(coordesp,ame.map))
    
    #determina los poligonos donde estan los puntos 
    idpuntos<- pt.in.poly$OBJECTID
    if (length(idpuntos)>10)
    {
      if (length(idpuntos)==1 && is.na(min(idpuntos))==T)
        Errores<- rbind(Errores,data.frame(Especie=NomArc,Puntos=length(idpuntos)))
      else
      {
        total<- ame.map@data$OBJECTID
        idtotal<- 1:length(total)
        npuntos<- length(idpuntos)
        t<- max(pt.card)*npuntos
        vecid<- rep(0,t)
        p<- 1
        for (k in 1:npuntos)
        {
          
          id<- which(total==idpuntos[k])
          if (length(id)>0 && id!=0)
          {
            #Encuentra a los vecinos
            vecid[p]<- id
            id.ady<- pt.ady.tot[[id]]
            p<-p+1
            #agrega a los vecinos al vector
            if (min(id.ady)!= 0)
            {
              for (l in 1:length(id.ady))
              {
                vecid[p]<- id.ady[l]
                p<- p+1
              }
            }
            rm(id.ady)
          }
        }
        
        #Elimina a los poligonos que no son vecinos
        vecid<- unique(vecid)
        vecid<- vecid[vecid != 0]
        del<- idtotal[-vecid]
        ame.map<- ame.map[-del,]
        
        plot(ame.map)
        points(coordesp,col='steelblue', pch=20)
        
        #guarda el shapefile con los poligonos adecuados
        writePolyShape(ame.map, paste(dirguardar,NomArc,sep=''))
      }
    }
    else
      Errores<- rbind(Errores,data.frame(Especie=NomArc,Puntos=length(idpuntos)))
  }
}

#Guardamos la tabla de errores
Errores<- Errores[-1,]
write.csv(Errores,paste(Dirsave,'Errores.csv',sep=''))

#Calcula el tiempo de procesamiento 
(proc.time() - ptm)
