#---------------------------LABORATORIO 44-A------------------------
#----------ACADÉMICO: CARLA CAROLINA PÉREZ HERNÁNDEZ----------------
#---------------ALUMNA: ANA GRISEL SANJUAN MERIDA-------------------

#------------------A tale of two heatmap functions------------------
#https://divingintogeneticsandgenomics.rbind.io/post/a-tale-of-two-heatmap-functions/

#Se instala la paquetería
install.packages("gplots")
install.packages("heatmaps")
install.packages("pheatmaps")

#Se cargan las librerías
library(stats)
library(gplots)

#Medida de la expresión génica de 4 genes (h1, h2, l1 y l2) en 8 puntos de tiempo
h1 <- c(10,20,10,20,10,20,10,20)
h2 <- c(20,10,20,10,20,10,20,10)
l1 <- c(1,3,1,3,1,3,1,3)
l2 <- c(3,1,3,1,3,1,3,1)

#Generamos matriz denominada mat para enlazar los genes
mat <- rbind(h1,h2,l1,l2)

#Se corren los siguientes comandos para generar el plot
par(mfrow =c(1,1), mar=c(4,4,1,1))
plot(1:8,rep(0,8), ylim=c(0,35), pch="", xlab="Time", ylab="Gene Expression")

for (i in 1:nrow(mat)) {
  lines(1:8,mat[i,], lwd=3, col=i)
}

legend(1,35,rownames(mat), 1:4, cex=0.7)

#Para calcular la distancia
dist(mat)

#Se usa el método predeterminado para el enlace: completo
plot(hclust(dist(mat)))

#Headmap predeterminado, para obtener explicitamente parámetros
#Escala en los renglones con escala de color verde y roja
heatmap(mat, Colv=NA, col=greenred(10), scale = "row")

#Heatmap con la escala desactivada
heatmap(mat, Colv = NA, col=greenred(10), scale = "none")

#Escalando los genes antes de introducirlos al heatmap
mat.scaled<- t(scale(t(mat), center=TRUE, scale = TRUE))
mat.scaled

#Cambio e la distancia entre genes
dist(mat.scaled)

#Obtener el plot
#Ahora h1 y l1 están agrupados juntos; l2 y h2 están agrupados juntos
plot(hclust(dist(mat.scaled)))

#Generando el heatmap
heatmap(mat.scaled, Colv = NA, col=greenred(10), scale = "none")

#Si aún no se escalan los datos pero se desea que l1 y h1 se agrupen juntos,
#a igual que l2 y h2, se puede usar la medida de distancia diferentes
#Correlación entre genes para asignar valores de 1 y -1
cor(t(mat))

#Correlación para definir la distancia
1- cor(t(mat))

#Generamos el dendrograma
hc <- hclust(as.dist(1-cor(t(mat))))
plot(hc)

#Heatmap sin escala
heatmap(mat, Colv = NA, Rowv=as.dendrogram(hc), col=greenred(10), scale = "none")

#Heatmap con la escala en los renglones
heatmap(mat, Colv = NA, Rowv=as.dendrogram(hc), col=greenred(10), scale = "row")

#Valores predeterminados de heatmap.2 con ninguna escala
heatmap.2(mat, trace = "none", Colv= NA, dendrogram = "row", scale = "none")

#Funciones de heatmap en R: primero se agrupa y luego usa el argumento de escala
#(si está configurado) para representar los datos
heatmap.2(mat, trace = "none", Colv= NA, dendrogram = "row", scale = "row")

#Datos escalados explícitamente primero y uso de la distancia euclidiana
heatmap.2(t(scale(t(mat), center=TRUE, scale=TRUE)), trace = "none", Colv= NA, dendrogram = "row", scale = "none")

#Usando 1- cor(x) como distancia y no escalar
heatmap.2(mat, trace = "none", 
          Colv= NA, dendrogram = "row",
          scale = "none",
          hclust=function(x) hclust(x, method='complete'), distfun=function(x) as.dist(1-cor(t(x))))

#Usando 1- cor(x) como distancia y no escalar
#Usar la escala en la función heatmap.2 para representar los colores
heatmap.2(mat, trace = "none", 
          Colv= NA, dendrogram = "row",
          scale = "row",
          hclust=function(x) hclust(x, method='complete'), distfun=function(x) as.dist(1-cor(t(x))))

#Escala y uso de 1- cor(x) como distancia
heatmap.2(t(scale(t(mat), center=TRUE, scale=TRUE)), trace = "none", 
          Colv= NA, dendrogram = "row",
          hclust=function(x) hclust(x, method='complete'), distfun=function(x) as.dist(1-cor(t(x))))

#-----------------------------FIN DE LABORATORIO 44-A----------------------------------