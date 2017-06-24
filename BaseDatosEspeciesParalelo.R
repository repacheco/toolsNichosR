library (XLConnect)    #Leer archivos de excel
library(parallel)
  
  genarchivos<- function (x)
  {
    e<- tipo.especie[x];
    NomEsp<- substr(e,1,nchar(e)-4)
    Abundancias<- read.csv(paste(Diresp,e,sep=''))
    
    #Declaramos el dataframe para guardar los elementos que querramos
    MatAbun<- data.frame(Especie="",Ruta=0,Latitud=0,Longitud=0,Abundancias=0)
    MatAbun$Especie<- as.character(MatAbun$Especie)
    
    #Ciclo para las especies
    for (i in 1:length(Rutas$countrynum)){
      #Encontramos sus coordenadas en el archivo de rutas
      Latitud<- Rutas$Lati[i]
      Longitud<- Rutas$Longi[i]
      Rut<- Rutas$Route[i]
      
      #Encontramos el promedio de las abundancias
      cond1<- Abundancias$countrynum==Rutas$countrynum[i] & Abundancias$statenum==Rutas$statenum[i] & Abundancias$Route==Rutas$Route[i]
      SelAbun<- Abundancias[cond1,]
      ValAou<- unique(SelAbun$Aou)
      for (l in 1:length(ValAou)){
        cond2<- SelAbun$Aou==ValAou[l]
        tam<- length(SelAbun[cond2,]$SpeciesTotal)                 
        if (tam!=0){
          
          #Encontramos su abundancia
          PromAbun<- mean(SelAbun[cond2,]$SpeciesTotal)
          
          #Encontramos su especie
          cond3<- Especies$AOU==ValAou[l]
          Specie<- Especies[cond3,]$Species[1]
          genero<- Especies[cond3,]$Genus[1]
          
          #realizamos un tipo de formato para la especie
          gen<- substr(genero,1,3)
          esp<- substr(Specie,1,3)
          especie<- paste(gen,esp,sep="_");
          
          #Rellenamos el dataframe
          MatAbun<- rbind(MatAbun,c(especie,Rut,Latitud,Longitud,PromAbun))
        }
      }          
    }
    
    #Eliminamos el primer elemento
    MatAbun<- MatAbun[-1,];
    
    #Guardamos valores segun su especie y el 50% de los datos aleatorios
    Especies.unicas<- unique(MatAbun$Especie)
    totesp<- length(Especies.unicas)
    
    #Creamos el directorio
    if (Tipo==1)
    {
      dir.create(paste(Dirsave,"PuntosSA/",NomEsp,sep=""))
      dir.create(paste(Dirsave,"PuntosCA/",NomEsp,sep=""))
    }
    else
    {
      dir.create(paste(Dirsave,'Puntos/',NomEsp,sep=''))
    }
    
    for (i in 1:totesp){
      #Obtenemos cada especie y encontramos el 50% de manera aleatoria
      MatEspecie<- MatAbun[MatAbun$Especie==Especies.unicas[i],]
      
      if (Tipo==1)
      {
          tamesp<- length(MatEspecie$Especie)
          total<- 1:tamesp
          aleatorizacion<- sample(1:tamesp,ceiling(tamesp/2))
          Muestra1<- MatEspecie[aleatorizacion,1:5]
          Muestra2<- MatEspecie[total[-aleatorizacion],] 
          
          #Checamos la correlacion espacial por medio del maximo de la ruta
          mu<- as.integer(unique(Muestra1$Rut))
          n.rutas<- length(mu)
          Muestra1.rutas.max<- data.frame(Especie="",Latitud=0,Longitud=0)
          for (k in 1:n.rutas)
          {
            condr<- subset(Muestra1, Muestra1$Rut==mu[k])
            rutas.maxim<- condr[condr$Abundancias==max(condr$Abundancias,na.rm=T),c(1,3,4)]
            Muestra1.rutas.max<- rbind(Muestra1.rutas.max,rutas.maxim)
          }
          Muestra1.rutas.max<- Muestra1.rutas.max[-1,]
          
          #Checamos la correlacion espacial por medio del maximo de la ruta
          mu<- as.integer(unique(Muestra2$Rut))
          n.rutas<- length(mu)
          Muestra2.rutas.max<- data.frame(Especie="",Latitud=0,Longitud=0,Abundancias=0)
          for (k in 1:n.rutas)
          {
            condr<- subset(Muestra2, Muestra2$Rut==mu[k])
            mx<- max(condr$Abundancias,na.rm = T)
            rutas.maxim<- condr[condr$Abundancias==mx,c(1,3,4,5)]
            Muestra2.rutas.max<- rbind(Muestra2.rutas.max,rutas.maxim)
          }
          Muestra2.rutas.max<- Muestra2.rutas.max[-1,]
      }
      else
      {
        if (Tipo==2)
        {
           Muestra<- MatEspecie
           mu<- as.integer(unique(Muestra$Rut))
           n.rutas<- length(mu)
           Muestra.rutas.max<- data.frame(Especie="",Latitud=0,Longitud=0,Abundancias=0)
           for (k in 1:n.rutas)
           {
             condr<- subset(Muestra, Muestra$Rut==mu[k])
             mx<- max(condr$Abundancias,na.rm = T)
             rutas.maxim<- condr[condr$Abundancias==mx,c(1,3,4,5)]
             Muestra.rutas.max<- rbind(Muestra.rutas.max,rutas.maxim)
           }
           Muestra.rutas.max<- Muestra.rutas.max[-1,]
        }
        else
            Muestra<- MatEspecie[,-2]
      }
      
      #Verifica si no tiene un caracter no valido
      if (length(grep("/",Especies.unicas[i]))==1)
        Especies.unicas[i]<- sub("/","-",Especies.unicas[i])
      
      if (length(grep(".",Especies.unicas[i]))==1)
         Especies.unicas[i]<- sub("\\.","",Especies.unicas[i])
      
      if (length(grep('"',Especies.unicas[i]))==1)
          Especies.unicas[i]<- sub('"',"",Especies.unicas[i])
      
      #Guardamos los resultados en archivos .csv
      if (Tipo==1)
      {
          if (nrow(Muestra1.rutas.max)>=10 && nrow(Muestra2.rutas.max)>=10)
          {
             write.csv(Muestra1.rutas.max,paste(Dirsave,'PuntosSA/',NomEsp,'/',Especies.unicas[i],".csv",sep=""))
             write.csv(Muestra2.rutas.max,paste(Dirsave,'PuntosCA/',NomEsp,'/',Especies.unicas[i],".csv",sep=""))
          }
          else
          {
            if (nrow(Muestra1.rutas.max)<10)
                   write.table(paste(NomEsp,'_',Especies.unicas[i],sep=''),paste(Dirsave,'PocosPuntosSA.txt',sep=''),append = T)
            
            if (nrow(Muestra2.rutas.max)<10) 
                  write.table(paste(NomEsp,'_',Especies.unicas[i],sep=''),paste(Dirsave,'PocosPuntosCA.txt',sep=''),append = T)
          }
      }
      else
      {
        if (Tipo==2)
        {
          if (nrow(Muestra.rutas.max)>=10)
                write.csv(Muestra.rutas.max,paste(Dirsave,'Puntos/',NomEsp,"/",Especies.unicas[i],".csv",sep=""))
          else
                write.table(paste(NomEsp,'_',Especies.unicas[i],sep=''),paste(Dirsave,'PocosPuntos.txt',sep=''),append = T) 
        }
        else
        {
          if (nrow(Muestra.rutas.max)>=10)
              write.csv(Muestra,paste(Dirsave,'Puntos/',NomEsp,"/",Especies.unicas[i],".csv",sep=""))
          else
             write.table(paste(NomEsp,'_',Especies.unicas[i],sep=''),paste(Dirsave,'PocosPuntos.txt',sep=''),append = T) 
        }
      }
    }  
 }



#####################################

#inicia el tiempo de procesamiento
ptm <- proc.time();

Rutas<- read.csv("~/Desktop/Proyecto/Especies/routes.csv")
#Lectura de los datos procedentes de excel
Especies<- readWorksheetFromFile("~/Desktop/Proyecto/Especies/SpeciesList1.xlsx",1)

#Buscara de un directorio todos los archivos con las especies
Diresp<- '~/Desktop/Proyecto/Especies/restoespecies/'
Dirsave<- '~/Desktop/Proyecto/Especies/'
#Tipo 1 Realiza el 50% de los datos para validacion y 50% para prueba tomando los maximos por ruta
#Tipo 2 Toma el total de los datos y los maximos por ruta
#Tipo 3 Toma el total
Tipo<- 1
tipo.especie<- list.files(Diresp)

#Parte parelela
num.cores <- detectCores() - 1
cl <- makeCluster(num.cores)
clusterExport(cl,c('tipo.especie','Rutas','Especies','Diresp','Dirsave','Tipo'))
parLapply(cl,1:length(tipo.especie),genarchivos)
stopCluster(cl)

#Calcula el tiempo de procesamiento 
(proc.time() - ptm)

############################################################################
