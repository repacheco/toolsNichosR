library(maptools)
library(rgdal)
library(spdep)
library(parallel)
library(sp)


rm(list = ls())

LeeArchivoGeneral<- function(directorio,Archivo,datum)
{
  ame.map<- readShapeSpatial(paste(directorio,Archivo,sep=''))
  proj4string(ame.map)<- datum
  return(ame.map)
}


mkworker<- function (directorio,carpeta,datum,ame.map,dirguardar,pt.ady.tot,pt.card)
{
  
  force(carpeta)
  force(directorio)
  force(datum)
  force(ame.map)
  force(dirguardar)
  force(pt.ady.tot)
  force(pt.card)
  
    EncuentraPoligonos<- function (x,directorio,carpeta,datum,ame.map,pt.ady.tot,pt.card)
    {
        arch<- carpeta[x]
        coordesp<- read.csv(paste(directorio,'/',arch,sep=''),header=T)
        coordinates(coordesp) <- c("Longitud", "Latitud")  # set spatial coordinates
        proj4string(coordesp) <- datum  # define projection system of our data
        NomArc<- substr(arch,1,nchar(arch)-4)
      
        
        #determina los puntos sobre los poligonos
        pt.in.poly <- unique(over(coordesp,ame.map))
        
        #determina los poligonos donde estan los puntos 
        idpuntos<- pt.in.poly$OBJECTID
        if (length(idpuntos)>10)
        {
            if (length(idpuntos)==1 && is.na(min(idpuntos))==T)
               write.table(NomArc,'Errores.txt',append = T)
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
                          vecid[p]<- id
                          id.ady<- pt.ady.tot[[id]]
                          p<-p+1
                          #agrega a los vecinos
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
              
            vecid<- unique(vecid)
            vecid<- vecid[vecid != 0]
            del<- idtotal[-vecid]
            ame.map<- ame.map[-del,]
            
            plot(ame.map)
            plot(coordesp,col='steelblue', pch=20, add=T)
            
            #guarda el shapefile con los poligonos adecuados
            writePolyShape(ame.map, paste(dirguardar,NomArc,sep=''))
           }
        }
        else
          write.table(NomArc,'Errores.txt',append = T)
        #return(list(mapa=ame.map,puntos=coordesp))
    }
    
    worker<- function(x)
    {
       EncuentraPoligonos(x,directorio,carpeta,datum,ame.map,pt.ady.tot,pt.card)
    }
    return(worker)
}


#inicia el tiempo de procesamiento
ptm <- proc.time();

#Datos Generales
setwd('~/Desktop/Proyecto/Especies/')
datum <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")  # geographical, datum WGS84
dirshpgral<- '~/Desktop/Proyecto/Especies/wwfreg/'
ame.map<- LeeArchivoGeneral(dirshpgral,"wwf_NAme",datum)

#poligonos adyacentes a los poligonos de los puntos
pt.ady.tot<- poly2nb(ame.map)
pt.card<- card(pt.ady.tot)

#leer el directorio donde estan los archivos
dirgral<- '~/Desktop/Proyecto/Especies/PuntosSA/';
Dirsave<- '~/Desktop/Proyecto/Especies/ResultM/';
dircarpetas<- list.files(dirgral)
numdircar<- length(dircarpetas)

num.cores <- detectCores() - 1
cl <- makeCluster(num.cores)
clusterEvalQ(cl, library(sp))
clusterEvalQ(cl, library(maptools))


for (i in 1:numdircar)
{
     carpeta<- list.files(paste(dirgral,dircarpetas[i],sep=''))
     numcar<- length(carpeta)
     dirguardar<- paste(Dirsave,dircarpetas[i],'/', sep='')
     directorio<- paste(dirgral,dircarpetas[i],sep='')
     
     if (dir.exists(dirguardar)==F)
          dir.create(dirguardar)
     
     parLapply(cl,1:numcar,mkworker(directorio,carpeta,datum,ame.map,dirguardar,pt.ady.tot,pt.card))
}

stopCluster(cl)

#Calcula el tiempo de procesamiento 
(proc.time() - ptm)



