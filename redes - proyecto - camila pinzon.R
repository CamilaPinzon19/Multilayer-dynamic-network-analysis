
#                                            PROYECTO ANALISIS DE REDES
#
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# ::::::::::::::::::::::::::::::::::::::::::::::::
####            0.1 Carga de librerias           ####
# ::::::::::::::::::::::::::::::::::::::::::::::::

library(igraph)
library(tensorTS)
library(fable)
library(tsibble)
library(feasts)
library(dplyr)
library(timetk)
library(zoo)
library(expm)
library(ggplot2)
library(psych)

# ::::::::::::::::::::::::::::::::::::::::::::::::
####            0.2 Funciones auxiliares           ####
# ::::::::::::::::::::::::::::::::::::::::::::::::

# - - - - - - - Función para visualizar las redes - - - - - - - - - - - - - 

mplot.ts <- function(xx,
                     titulos_filas = dimnames(xx)[[3]],
                     titulos_columnas = dimnames(xx)[[2]],
                     xlab = "Tiempo") {
  
  d <- dim(xx)
  
  op <- par(
    mfrow = c(d[3], d[2]),
    mar = c(1, 2.5, 1.5, 0.5),
    oma = c(4, 3, 1, 1)
  )
  on.exit(par(op))
  
  for(j in seq_len(d[3])) {
    for(i in seq_len(d[2])) {
      
      plot(xx[, i, j],
           type = "l",
           lwd = 1,
           xlab = "",
           ylab = "",
           xaxt = "n")
      
      ## Eje x solo en la última fila
      if(j == d[3])
        axis(1)
      
      ## Eje y solo en la primera columna
      if(i == 1)
        axis(2)
      
      ## Títulos de las columnas
      if(j == 1)
        title(main = titulos_columnas[i], cex.main = 1.2)
      
      ## Títulos de las filas
      if(i == 1)
        mtext(titulos_filas[j], side = 2, line = 2, cex = 0.8, font = 2)
    }
  }
  
  mtext(xlab, side = 1, outer = TRUE, line = 2)
}

# - - - - - - - Grafica los ACF de los datos tipo matriz - - - - - - - - - - - - - 

mplot.acf <- function(xx,
                      titulos_filas = dimnames(xx)[[2]],
                      titulos_columnas = dimnames(xx)[[3]]) {
  
  if (isS4(xx)) xx <- xx@data
  d <- dim(xx)
  
  op <- par(mfrow = c(d[2], d[3]), mar = c(0.5, 0.5, 1.5, 0.5), oma = c(4, 3, 1, 1))
  on.exit(par(op))
  
  for(i in seq_len(d[2])) {
    for(j in seq_len(d[3])) {
      
      acf(xx[, i, j],
          ylim = c(0, 1),
          xaxt = if(i == d[2]) "s" else "n",
          yaxt = if(j == 1) "s" else "n",
          xlab = "",
          ylab = "",
          main = "")
      
      ## Títulos de las columnas
      if(i == 1)
        title(main = titulos_columnas[j], cex.main = 1.2, font.main = 2)
      
      ## Títulos de las filas
      if(j == 1)
        mtext(titulos_filas[i], side = 2, line = 2, cex = 0.8, font = 2)
    }
  }
  
  ## Etiquetas comunes
  mtext("Rezago", side = 1, outer = TRUE, line = 2)
}

# - - - - - - - Función para construir la tabla de significancia - - - - - - - - - - - - - 

tabla_significancia <- function(coef_mat, sd_mat) {
  
  t_mat <- coef_mat / sd_mat
  p_mat <- 2 * (1 - pnorm(abs(t_mat)))  # aproximación normal (Wald)
  
  df <- expand.grid(Fila = rownames(coef_mat), Columna = colnames(coef_mat))
  df$Coef <- as.vector(coef_mat)
  df$SE   <- as.vector(sd_mat)
  df$t    <- as.vector(t_mat)
  df$p    <- as.vector(p_mat)
  
  df$Sig <- cut(df$p,
                breaks = c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf),
                labels = c("***", "**", "*", ".", ""))
  df <- df[order(df$Fila, df$Columna), ]
  df[, c("Fila", "Columna", "Coef", "SE", "t", "p", "Sig")]
}

# - - - - - - - Función impulso-respuesta - - - - - - - - - - - - - 
mar_oirf <- function(A1, A2, Sigma, i, j, H = 10) {
  d1 <- nrow(A1); d2 <- nrow(A2)
  idx     <- d1 * (j - 1) + i        # posición de (i,j) en vec(Et), columna-mayor
  col_sig <- Sigma[, idx]
  sd_ij   <- sqrt(Sigma[idx, idx])   # para normalizar a choque de 1 sd
  
  irf <- array(NA, dim = c(H + 1, d1, d2))
  for (k in 0:H) {
    Ak <- if (k == 0) diag(d1) else A1 %^% k
    Bk <- if (k == 0) diag(d2) else A2 %^% k
    Fk_vec      <- kronecker(Bk, Ak) %*% col_sig / sd_ij
    irf[k+1,,]  <- matrix(Fk_vec, nrow = d1, ncol = d2)  # column-major, coincide con vec()
  }
  irf
}

# - - - - - - - Graficar las curvas impulso-respuesta generadas - - - - - - - - - - - - - 
plot_oirf_panel <- function(irf, capa_labels = NULL, caract_labels = NULL,
                            ref_lines       = c(-0.1, 0.1),
                            titulo_shock    = "",
                            same_scale      = TRUE) {
  
  H1 <- dim(irf)[1]; d1 <- dim(irf)[2]; d2 <- dim(irf)[3]
  H  <- H1 - 1
  
  yr <- if (same_scale) range(irf, ref_lines, na.rm = TRUE) else NULL
  
  op <- par(mfrow = c(d1, d2), mar = c(0.3, 0.3, 0.3, 0.5),
            oma = c(3.5, 3, 3, 1), xaxs = "i")
  
  for (i in 1:d1) {
    for (j in 1:d2) {
      y_range <- if (same_scale) yr else range(irf[, i, j], ref_lines)
      
      # 1) lienzo vacío (sin dibujar la curva todavía)
      plot(0:H, irf[, i, j], type = "n",
           ylim = y_range, xlab = "", ylab = "",
           xaxt = if (i == d1) "s" else "n",
           yaxt = if (j == 1)  "s" else "n",
           cex.axis = 0.7)
      
      # 2) líneas de referencia -- quedan al fondo
      abline(h = ref_lines, lty = 3, col = "gray40")
      abline(h = 0, col = "gray80")
      
      # 3) la curva se dibuja al final -- queda por encima
      lines(0:H, irf[, i, j], lwd = 1)
      
      box()
      
      if (j == 1) {
        text(par("usr")[1] + 0.05 * diff(par("usr")[1:2]),
             par("usr")[3] + 0.08 * diff(par("usr")[3:4]),
             capa_labels[i], adj = 0, cex = 1, font = 2)
      }
      if (i == 1) {
        text(par("usr")[2] - 0.05 * diff(par("usr")[1:2]),
             par("usr")[4] - 0.12 * diff(par("usr")[3:4]),
             caract_labels[j], adj = 1, cex = 1, font = 2)
      }
    }
  }
  
  mtext("Rezago (k)", side = 1, outer = TRUE, line = 2, cex = 0.9)
  mtext(titulo_shock, side = 3, outer = TRUE, line = 1, cex = 1.05, font = 2)
  
  par(op)
}

# - - - - - - - Gráficos para el diagnostico de los residuales - - - - - - - - - - - - - 
mplot.diag <- function(xx,
                       type = c("acf", "pacf", "qq"),
                       titulos_filas = dimnames(xx)[[2]],
                       titulos_columnas = dimnames(xx)[[3]]) {
  
  type <- match.arg(type)
  if (isS4(xx)) xx <- xx@data
  d <- dim(xx)
  
  op <- par(mfrow = c(d[2], d[3]), mar = c(0.5, 0.5, 1.5, 0.5), oma = c(4, 3, 1, 1))
  on.exit(par(op))
  
  for (i in seq_len(d[2])) {
    for (j in seq_len(d[3])) {
      
      serie <- xx[, i, j]
      
      if (type == "acf") {
        acf(serie, ylim = c(0, 1),
            xaxt = if (i == d[2]) "s" else "n",
            yaxt = if (j == 1) "s" else "n",
            xlab = "", ylab = "", main = "")
        
      } else if (type == "pacf") {
        pacf(serie, ylim = c(-1, 1),
             xaxt = if (i == d[2]) "s" else "n",
             yaxt = if (j == 1) "s" else "n",
             xlab = "", ylab = "", main = "")
        
      } else if (type == "qq") {
        qqnorm(serie,
               xaxt = if (i == d[2]) "s" else "n",
               yaxt = if (j == 1) "s" else "n",
               xlab = "", ylab = "", main = "", pch = 19, cex = 0.5)
        qqline(serie, col = "red", lwd = 1.5)
      }
      
      ## Títulos de las columnas
      if (i == 1)
        title(main = titulos_columnas[j], cex.main = 1.2, font.main = 2)
      
      ## Títulos de las filas
      if (j == 1)
        mtext(titulos_filas[i], side = 2, line = 2, cex = 0.8, font = 2)
    }
  }
  
  ## Etiquetas comunes según el tipo
  etiqueta_x <- switch(type,
                       acf  = "Rezago",
                       pacf = "Rezago",
                       qq   = "Cuantiles teóricos")
  mtext(etiqueta_x, side = 1, outer = TRUE, line = 2)
}

# ::::::::::::::::::::::::::::::::::::::::::::::::
####            0.3 Preparar los datos           ####
# ::::::::::::::::::::::::::::::::::::::::::::::::

load("icews_tensor_proyecto.RData")

n_tiempos <- dim(tensor_array)[4]
n_capas <- dim(tensor_array)[3]
n_paises <- dim(tensor_array)[1]

nombres_capas <- c("Material-", "Material+", "Verbal-", "Verbal+")
nombres_paises <- dimnames(tensor_array)$Source

# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
####            1. Elegir los 25 países más centrales           ####
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Arreglo para guardar las centralidades
page_cent <- array(NA,
                  dim = c(n_paises, n_capas, n_tiempos),
                  dimnames = list(nombres_paises, paste0("Layer_",1:n_capas), paste0("Time_",1:n_tiempos)))

for(m in 1:n_capas){
  for(t in 1:n_tiempos){
    A <- tensor_array[,,m,t]
    g <- graph_from_adjacency_matrix(A, mode = "directed", weighted = TRUE, diag = FALSE)
    page_cent[,m,t] <- page_rank(g, directed = TRUE, weights = E(g)$weight)$vector
  }
}

# Promedio de centralidades a través del tiempo
cent_mean <- apply(page_cent, c(1,2), mean, na.rm = TRUE)

# Obtener el top de los 25 países más centrales
top25_layers <- vector("list", n_capas)
for(m in 1:n_capas){
  ord <- order(cent_mean[,m], decreasing = TRUE)
  top25_layers[[m]] <- data.frame(Country = rownames(cent_mean)[ord[1:25]],
                                  MeanEigenCentrality = cent_mean[ord[1:25], m])
}
top25_names <- lapply(top25_layers, \(x) x$Country)
# Los países más centrales para todas las capas 
Reduce(intersect, top25_names) 
# Son 19 países

# - - - - - - - Top 25 global - - - - - - - - - - - - - 
overall_mean <- rowMeans(cent_mean)
ord <- order(overall_mean, decreasing = TRUE)
top25_global <- data.frame(Country = names(overall_mean)[ord[1:25]],
                           MeanEigenCentrality = overall_mean[ord[1:25]])
top25_global # Son los 19 paises de arriba más otros 6

# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
####            2. Obtener las matrices de estadísticas           ####
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

paises_final <- top25_global$Country

# Indices de los países a usar
indices <- match(paises_final, nombres_paises)

# Estadísticos a calcular
estadisticos <- c("densidad",
                  "asortatividad",
                  "transitividad",
                  "reciprocidad",
                  "dist_prom",
#                  "diametro",
                  "fuerza_media"
#                 ,"grado_medio"
)
n_estad <- length(estadisticos)

stats_array <- array(NA, 
                     dim = c(n_capas, length(estadisticos), n_tiempos),
                     dimnames = list(paste0("Layer_",1:n_capas), estadisticos, paste0("Time_",1:n_tiempos) ) )

for(m in 1:n_capas){
  for(t in 1:n_tiempos){
    A <- tensor_array[indices, indices, m, t]
    g <- graph_from_adjacency_matrix(A, mode = "directed", weighted = TRUE, diag = FALSE)
    comp <- igraph::components(as_undirected(g))
    
    stats_array[m,"densidad",t] <- edge_density(g)
    stats_array[m,"asortatividad",t] <- assortativity_degree(g, directed = TRUE)
    stats_array[m,"transitividad",t] <- transitivity(g, type="global")
    stats_array[m,"reciprocidad",t] <- reciprocity(g)
    stats_array[m,"dist_prom",t] <- mean_distance(g, directed = TRUE, unconnected = TRUE)
#    stats_array[m,"diametro",t] <- diameter(g, directed = TRUE, unconnected = TRUE)
    stats_array[m,"fuerza_media",t] <- mean(strength(g, mode = "all"))
    #stats_array[m,"grado_medio",t] <- mean(degree(g, mode = "all"))
  }
  
}
stats_array <- aperm(stats_array, c(3, 1, 2))

# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
####        3. Evaluar estacionariedad de las matrices         ####
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

mplot.ts(stats_array,
         titulos_columnas = nombres_capas,
         titulos_filas = c("Densidad", "Asortatividad", "Transitividad", "Reciprocidad", "Distancia media", "Fuerza media"),
         xlab = "Tiempo")
mplot.acf(stats_array, 
          titulos_filas = nombres_capas,
          titulos_columnas = c("Densidad", "Asortatividad", "Transitividad", "Reciprocidad", "Distancia media", "Fuerza media")) 
# Los estadísticos presentan tendencia a través del tiempo

# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#####          3.1  Estimar la tendencia y eliminarla           ####
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Se prueban 3 métodos para estimar y eliminar la tendencia

tend_decompose <- array(NA, dim = dim(stats_array))
tend_STL_fable <- array(NA, dim = dim(stats_array))
tend_LOESS <- array(NA, dim = dim(stats_array))

sin_tend_decompose <- array(NA, dim = dim(stats_array))
sin_tend_STL_fable <- array(NA, dim = dim(stats_array))
sin_tend_LOESS <- array(NA, dim = dim(stats_array))

# Obtener las tendencias con cada método

for(i in 1:n_capas){
  for(j in 1:n_estad){
    serie <- ts(stats_array[,i,j], frequency = 12)
    # decompose
    fit <- stats::decompose(serie)$trend
    tend_decompose[,i,j] <- fit
    # STL fable
    df <- data.frame(
      fecha = yearmonth(seq.Date(from = as.Date("2014-01-01"), by = "month", length.out = 132)),
      value = as.numeric(serie) )
    tsb <- as_tsibble(df, index = fecha)
    fit_fable <- tsb %>%
      model(STL(value ~ trend() + season(window="periodic")) ) %>%
      components()
    tend_STL_fable[,i,j] <- fit_fable$trend
    # LOESS
    lo <- loess(value ~ as.numeric(fecha), data = df, span = 0.75)
    tend_LOESS[,i,j] <- predict(lo)
  }
}

# Series sin tendencia para cada método
sin_tend_decompose <- stats_array - tend_decompose
sin_tend_STL_fable <- stats_array - tend_STL_fable
sin_tend_LOESS <- stats_array - tend_LOESS

# - - - - - - - Gráficos de las series sin tendencia - - - - - - - - - - - - - 
par(mfrow = c(n_capas, n_estad), mar = c(2,2,2,1))

for(i in 1:n_capas){
  for(j in 1:n_estad){
    ts.plot(
      ts(sin_tend_decompose[,i,j]),
      ts(sin_tend_STL_fable[,i,j]),
      ts(sin_tend_LOESS[,i,j]),
      col = 1:3, lwd = 1.5, main = paste("Serie", i, j), xlab = "", ylab = "")
  }
}
legend("topright", legend = c("Sin decompose", "Sin STL-fable", "Sin LOESS"),
       col = 1:3, lty = 1, cex = 0.5)
# De momento me quedo con el STL para quitar la tendencia

# - - - - - - - Gráfica de series finales sin tendencia - - - - - - - - - - - - - 
mplot.ts(sin_tend_STL_fable, 
         titulos_filas = c("Densidad", "Asortatividad", "Transitividad", "Reciprocidad", "Distancia media", "Fuerza media"), 
         titulos_columnas = nombres_capas)
mplot.acf(xx = sin_tend_STL_fable)

# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#####                 3.2  Estandarización previa              ####
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Para que el modelo no esté donimado por variables altamente dispersas
d1 <- n_capas # capas
d2 <- n_estad  # características

datos_std <- sin_tend_STL_fable
escalas <- numeric(d2)

for (j in 1:d2) {
  valores <- as.vector(sin_tend_STL_fable[, , j])  # toda capa y tiempo para esa característica
  escalas[j] <- sd(valores)
  datos_std[, , j] <- sin_tend_STL_fable[, , j] / escalas[j]
}
names(escalas) <- dimnames(sin_tend_STL_fable)[[3]]  
print(round(escalas, 4))

mplot.acf(datos_std, 
          titulos_filas = nombres_capas,
          titulos_columnas = c("Densidad", "Asortatividad", "Transitividad", "Reciprocidad", "Distancia media", "Fuerza media")) 

# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#####                     3.3 Justificar MAR(1)                 ####
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Ajustar modelos con varios rezagos
aa <- tenAR.est(datos_std, R = 1, P = 1, method = "LSE")
bb <- tenAR.est(datos_std, R = 1, P = 2, method = "LSE")
cc <- tenAR.est(datos_std, R = 1, P = 3, method = "LSE")
dd <- tenAR.est(datos_std, R = 1, P = 4, method = "LSE")
ee <- tenAR.est(datos_std, R = 1, P = 5, method = "LSE")
ff <- tenAR.est(datos_std, R = 1, P = 6, method = "LSE")

# Armar tabla comparativa de BIC
modelos <- list(P1 = aa, P2 = bb, P3 = cc, P4 = dd, P5 = ee, P6 = ff)

tabla_bic <- data.frame(P   = 1:6,
                        BIC = sapply(modelos, function(m) m$BIC) )
tabla_bic
# El mejor modelo (menor BIC) es el de p=1 rezagos

# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
####                     4. Estimación del MAR                 ####
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Reajusta el modelo con los datos estandarizados
modelo_MAR_std <- aa
mplot.acf(modelo_MAR_std$res)

# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#####          4.1  Significancia de los coeficientes         ####
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Matriz 1 (Capas con Capas)
mat1 <- modelo_MAR_std$A[[1]][[1]][[1]]
sd1  <- modelo_MAR_std$sd[[1]][[1]][[1]]
rownames(mat1) <- rownames(sd1) <- nombres_capas
colnames(mat1) <- colnames(sd1) <- nombres_capas

tabla_significancia(mat1, sd1)

# Matriz 2 (variables de red)
mat2 <- modelo_MAR_std$A[[1]][[1]][[2]]
sd2  <- modelo_MAR_std$sd[[1]][[1]][[2]]
rownames(mat2) <- rownames(sd2) <- c("densidad","asortatividad","transitividad", "reciprocidad","dist_prom","fuerza_media")
colnames(mat2) <- colnames(sd2) <- c("densidad","asortatividad","transitividad", "reciprocidad","dist_prom","fuerza_media")

tabla_significancia(mat2, sd2)

# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#####                 4.2 Análisis impulso-respuesta             ####
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

caract_labels <- c("Densidad", "Asortatividad", "Transitividad", "Reciprocidad", "Dist promedio", "Fuerza media")

# Corroborar condicion de causalidad
rho_A1 <- max(Mod(eigen(modelo_MAR_std$A[[1]][[1]][[1]])$values))
rho_A2 <- max(Mod(eigen(modelo_MAR_std$A[[1]][[1]][[2]])$values))
cat("rho(A1) =", round(rho_A1, 3), " | rho(A2) =", round(rho_A2, 3),
    " | producto =", round(rho_A1 * rho_A2, 3), "\n")

# Choque de 1 sd en (capa 1, característica 1)
irf_calc <- mar_oirf(modelo_MAR_std$A[[1]][[1]][[1]], modelo_MAR_std$A[[1]][[1]][[2]], modelo_MAR_std$Sig, i = 1, j = 6, H = 6)

plot_oirf_panel(irf_calc,
                capa_labels   = nombres_capas,
                caract_labels = caract_labels,
                ref_lines     = c(-0.1, 0.1)#,
                #titulo_shock  = "s1-oIRF: choque de 1 sd en (Capa 1, Característica 1)"
)

# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#####                 4.3 Diagnostico de residuales             ####
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

res <- modelo_MAR_std$res
Tt  <- dim(res)[1]

mplot.diag(res, type = "acf",
           titulos_filas = nombres_capas,
           titulos_columnas = c("Densidad", "Asortatividad", "Transitividad", "Reciprocidad", "Distancia media", "Fuerza media")) 
mplot.diag(res, type = "pacf",
           titulos_filas = nombres_capas,
           titulos_columnas = c("Densidad", "Asortatividad", "Transitividad", "Reciprocidad", "Distancia media", "Fuerza media")) 

mplot.diag(res, type = "qq",
           titulos_filas = nombres_capas,
           titulos_columnas = c("Densidad", "Asortatividad", "Transitividad", "Reciprocidad", "Distancia media", "Fuerza media")) 



# Residuales
res_mat <- t(apply(res, 1, as.vector))  # T x N
colnames(res_mat) <- 1:ncol(res_mat)

lb_pvalores <- sapply(1:ncol(res_mat), function(j) {
  Box.test(res_mat[, j], lag = 10, type = "Ljung-Box")$p.value
})

# Corrección FDR (Benjamini-Hochberg) por las 24 pruebas simultáneas
lb_ajustado <- p.adjust(lb_pvalores, method = "BH")

# Para cada serie, contar cuántos rezagos (de los primeros 20) cruzan la banda de 95%
banda <- 1.96 / sqrt(nrow(res_mat))  # aprox banda de significancia del ACF

conteo_significativos <- sapply(1:ncol(res_mat), function(j) {
  acf_vals <- acf(res_mat[, j], lag.max = 20, plot = FALSE)$acf[-1]  # excluye rezago 0
  sum(abs(acf_vals) > banda)
})

tabla_diagnostico <- data.frame(
  Prueba = c("Ljung-Box univariado (24 series, corrección FDR)",
             "Conteo de cruces ACF vs. esperado bajo ruido blanco",
             "Normalidad multivariada (Mardia, asimetría)",
             "Normalidad multivariada (Mardia, curtosis)"),
  Resultado = c(
    paste0(sum(lb_ajustado < 0.05), " de 24 series con autocorrelación tras ajuste FDR"),
    paste0(sum(conteo_significativos), " observados vs. ", round(24*20*0.05,1), " esperados bajo H0"),
    paste0("Estad. = ", round(mardia_test$skew, 2), ", p = ", round(mardia_test$p.skew, 4)),
    paste0("Estad. = ", round(mardia_test$kurtosis, 2), ", p = ", round(mardia_test$p.kurt, 4))
  )
)
print(tabla_diagnostico, row.names = FALSE)
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
####                     5. Simulación                 ####
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

d1 <- 4
coupling_true <- matrix(0, d1, d1)
coupling_true[1, 1] <- 0.6
coupling_true[3, 1] <- 0.5
coupling_true[3, 3] <- 0.3
coupling_true[4, 3] <- 0.4
coupling_true[4, 4] <- 0.5

# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Generador con DOS fuentes de variación independientes:
#  (a) 'coupling' controla la densidad general (como antes)
#  (b) 'triangle_bias' controla el sesgo hacia cerrar triángulos, con su PROPIA dinámica temporal, independiente de (a)
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
simulate_temporal_network_v2 <- function(n_nodes, Tn, coupling,
                                         noise_sd = 0.05,
                                         triangle_coupling = NULL,  # matriz d1 x d1, propia
                                         burn = 50, threshold = 0.15,
                                         directed = TRUE) {
  
  d1 <- nrow(coupling)
  if (is.null(triangle_coupling)) triangle_coupling <- diag(0.4, d1)  # memoria propia moderada
  
  Ttot <- Tn + burn
  W    <- array(0, dim = c(Ttot, d1, n_nodes, n_nodes))
  Bias <- matrix(0.3, Ttot, d1)   # sesgo de triangulación por capa y tiempo
  
  for (l in 1:d1) W[1, l, , ] <- matrix(runif(n_nodes^2, 0, 0.3), n_nodes, n_nodes)
  
  for (t in 2:Ttot) {
    # (b) el sesgo de triangulación tiene su propia recursión, con SU PROPIO ruido
    for (l in 1:d1) {
      Bias[t, l] <- pmin(pmax(
        sum(triangle_coupling[l, ] * (Bias[t - 1, ])) + rnorm(1, 0, 0.05),
        0), 1)
    }
    
    for (l in 1:d1) {
      base <- matrix(0, n_nodes, n_nodes)
      for (lp in 1:d1) {
        if (coupling[l, lp] != 0) base <- base + coupling[l, lp] * W[t - 1, lp, , ]
      }
      <- matrix(rnorm(n_nodes^2, 0, noise_sd), n_nodes, n_nodes)
      Wt_raw <- pmin(pmax(base + noise, 0), 1)
      if (!directed) Wt_raw[lower.tri(Wt_raw)] <- t(Wt_raw)[lower.tri(Wt_raw)]
      diag(Wt_raw) <- 0
      
      Wt <- Wt_raw * (Wt_raw > threshold)
      
      # --- boost de triangulación, impulsado por Bias[t,l], independiente de la densidad ---
      adj_bin <- (Wt_raw > threshold) * 1
      common_neighbors <- adj_bin %*% t(adj_bin)
      diag(common_neighbors) <- 0
      
      # Fuerza directamente una fracción de los pares con vecinos comunes a conectarse
      n_boost <- sum(common_neighbors > 0 & adj_bin == 0)
      if (n_boost > 0) {
        candidatos <- which(common_neighbors > 0 & adj_bin == 0, arr.ind = FALSE)
        n_forzar <- round(length(candidatos) * Bias[t, l] * 0.5)
        if (n_forzar > 0) {
          elegidos <- sample(candidatos, min(n_forzar, length(candidatos)))
          Wt_raw[elegidos] <- runif(length(elegidos), threshold + 0.05, 0.6)
        }
      }
      Wt <- Wt_raw
      diag(Wt) <- 0
      
      # --- reciprocidad parcial, independiente de lo anterior ---
      mask_recip <- matrix(runif(n_nodes^2) < 0.4, n_nodes, n_nodes)
      Wt[mask_recip] <- pmax(Wt[mask_recip], t(Wt)[mask_recip])
      diag(Wt) <- 0
      
      W[t, l, , ] <- Wt
    }
  }
  W[(burn + 1):Ttot, , , ]
}

# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#                   Características
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
compute_network_features_v2 <- function(W) {
  Tn <- dim(W)[1]; d1 <- dim(W)[2]; n <- dim(W)[3]
  feat_names <- c("densidad", "asortatividad", "transitividad",
                  "reciprocidad", "dist_prom", "fuerza_media")
  d2 <- length(feat_names)
  X <- array(NA, dim = c(Tn, d1, d2), dimnames = list(NULL, NULL, feat_names))
  
  for (t in 1:Tn) {
    for (l in 1:d1) {
      g <- graph_from_adjacency_matrix(W[t, l, , ], mode = "directed",
                                       weighted = TRUE, diag = FALSE)
      X[t, l, "densidad"]      <- edge_density(g)
      X[t, l, "asortatividad"] <- tryCatch(assortativity_degree(g, directed = TRUE), error = function(e) 0)
      X[t, l, "transitividad"] <- tryCatch(transitivity(g, type = "global"), error = function(e) 0)
      X[t, l, "reciprocidad"]  <- tryCatch(reciprocity(g), error = function(e) 0)
      X[t, l, "dist_prom"]     <- tryCatch(mean_distance(g, directed = TRUE, unconnected = TRUE, weights = NA),
                                           error = function(e) 0)
      X[t, l, "fuerza_media"]  <- mean(strength(g, mode = "all", weights = E(g)$weight))
    }
  }
  X[is.na(X)] <- 0
  X
}

# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#         Una réplica: simula, estandariza, ajusta MAR
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
run_one_replicate_v2 <- function(n_nodes = 25, Tn = 132,
                                 coupling = coupling_true,
                                 noise_sd = 0.05,
                                 triangle_coupling = NULL,
                                 threshold = 0.15) {
  
  W_sim <- simulate_temporal_network_v2(n_nodes, Tn, coupling,
                                        noise_sd = noise_sd,
                                        triangle_coupling = triangle_coupling,
                                        threshold = threshold)
  
  X_sim <- compute_network_features_v2(W_sim)
  
  # --- estandariza cada característica ---
  X_std <- X_sim
  for (j in 1:dim(X_sim)[3]) {
    s <- sd(as.vector(X_sim[, , j]))
    if (s > 0) X_std[, , j] <- X_sim[, , j] / s
  }
  
  fit <- tryCatch(matAR.RR.est(X_std, method = "LSE"), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  
  se  <- get_se_full(fit, d1 = dim(coupling)[1], d2 = dim(X_sim)[3])
  sig <- sig_table(fit$A1, se$se_A1)
  
  verdad_capas <- coupling != 0
  detectado    <- sig$p < 0.05
  
  data.frame(
    sensibilidad  = mean(detectado[verdad_capas]),
    especificidad = mean(!detectado[!verdad_capas])
  )
}

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#####                     5.1 Monte Carlo                    ####
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
set.seed(123)
B <- 60
resultados_red <- do.call(rbind, lapply(1:B, function(b) run_one_replicate_v2()))

cat("Sensibilidad promedio (recupera acoplamientos reales):",
    round(mean(resultados_red$sensibilidad, na.rm = TRUE), 3), "\n")
cat("Especificidad promedio (no inventa acoplamientos falsos):",
    round(mean(resultados_red$especificidad, na.rm = TRUE), 3), "\n")

boxplot(resultados_red, main = "Detección de acoplamiento entre capas\n(generador: red dinámica real)")

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#####    5.2 Escenarios de simulación                           ####
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

run_scenario_red <- function(scenario_name, B = 40, coupling, n_nodes = 25, Tn = 132,
                             noise_sd = 0.05, triangle_coupling = NULL, threshold = 0.15) {
  
  resultados <- vector("list", B)
  for (b in 1:B) {
    set.seed(1000 + b)
    resultados[[b]] <- tryCatch({
      W_sim <- simulate_temporal_network_v2(n_nodes, Tn, coupling, noise_sd,
                                            triangle_coupling, threshold = threshold)
      X_sim <- compute_network_features_v2(W_sim)
      
      X_std <- X_sim
      for (j in 1:dim(X_sim)[3]) {
        s <- sd(as.vector(X_sim[, , j]))
        if (s > 0) X_std[, , j] <- X_sim[, , j] / s
      }
      fit <- matAR.RR.est(X_std, method = "LSE")
      se  <- get_se_full(fit, d1 = nrow(coupling), d2 = dim(X_sim)[3])
      sig <- sig_table(fit$A1, se$se_A1)
      
      verdad <- coupling != 0
      detectado <- sig$p < 0.05
      
      data.frame(sensibilidad = mean(detectado[verdad]),
                 especificidad = mean(!detectado[!verdad]))
    }, error = function(e) NULL)
  }
  out <- do.call(rbind, resultados)
  if (is.null(out)) return(NULL)
  out$scenario <- scenario_name
  out
}

# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#             Escenarios de simulación
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# --- Escenario 1: caso base (el que ya calibraste) ---
esc1 <- run_scenario_red("Base", coupling = coupling_true)

# --- Escenario 2: SIN ningún acoplamiento real
coupling_null <- matrix(0, d1, d1)
esc2 <- run_scenario_red("Control negativo (sin acoplamiento)", coupling = coupling_null)

# --- Escenario 3: acoplamiento débil
coupling_debil <- coupling_true
coupling_debil[coupling_debil != 0] <- coupling_debil[coupling_debil != 0] * 0.5
esc3 <- run_scenario_red("Acoplamiento débil (mitad de fuerza)", coupling = coupling_debil)

# --- Escenario 4: más ruido en la dinámica de aristas ---
esc4 <- run_scenario_red("Ruido alto", coupling = coupling_true, noise_sd = 0.15)

# --- Escenario 5: red más pequeña
esc5 <- run_scenario_red("Red pequeña (n=12)", coupling = coupling_true, n_nodes = 12)

# --- Escenario 6: red más grande
esc6 <- run_scenario_red("Red grande (n=50)", coupling = coupling_true, n_nodes = 50)

# --- Escenario 7: T más corto
esc7 <- run_scenario_red("T corto (T=60)", coupling = coupling_true, Tn = 60)

# --- Escenario 8: el sesgo de triangulación (2da fuente de variación)
tri_coupling_acoplado <- matrix(0, d1, d1)
diag(tri_coupling_acoplado) <- 0.4
tri_coupling_acoplado[3, 1] <- 0.3   # el sesgo de triangulación también fluye 1->3
esc8 <- run_scenario_red("Acoplamiento también en triangulación",
                         coupling = coupling_true,
                         triangle_coupling = tri_coupling_acoplado)

# --- Escenario 9: umbral de arista más exigente -- 
esc9 <- run_scenario_red("Umbral más alto (redes más dispersas)",
                         coupling = coupling_true, threshold = 0.25)

resultados_todos <- rbind(esc1, esc2, esc3, esc4, esc5, esc6, esc7, esc8, esc9)

# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#####         5.3 Gráficos de sensibilidad y especificidad.    ####
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Crear identificador numérico para cada escenario
niveles <- c(esc1$scenario[1], esc2$scenario[1], esc3$scenario[1], esc4$scenario[1],
             esc5$scenario[1], esc6$scenario[1], esc7$scenario[1], esc8$scenario[1], esc9$scenario[1])
resultados_todos$escenario_id <- factor(resultados_todos$scenario,
                                        levels = niveles,
                                        labels = seq_along(niveles) )

# Especificidad
ggplot(resultados_todos,
       aes(x = escenario_id, y = especificidad)) +
  geom_boxplot(fill = "steelblue", alpha = 0.8, outlier.size = 1, staplewidth = 0.4) +
  geom_hline(yintercept = 0.95, linetype = 2, colour = "red") +
  labs(x = "Escenario",
       y = "Especificidad") +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0),
    plot.title = element_text(hjust = 0.5)
  )

# Sensibilidad
ggplot(resultados_todos,
       aes(x = escenario_id, y = sensibilidad)) +
  geom_boxplot(fill = "grey70", alpha = 0.8, outlier.size = 1, staplewidth = 0.4) +
  geom_hline(yintercept = 0.95, linetype = 2, colour = "red") +
  labs(x = "Escenario",
       y = "Sensibilidad") +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0),
    plot.title = element_text(hjust = 0.5)
  )
