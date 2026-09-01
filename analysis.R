# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#               MULTILAYER DYNAMIC NETWORK ANALYSIS - MAIN SCRIPT
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
# Description: 
# Main code for structural analysis of the multi-layered dynamic network of 
# interactions between countries (ICEWS). 
# This file contains the selection of the 25 most central countries, the fitting 
# of the Matrix Autoregressive (MAR) model, and the final simulation.

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
####                           0.1 Library loading                          ####
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

packages <- c('igraph','tensorTS','fable','tsibble','feasts', 'dplyr', 'timetk', 
              'zoo', 'expm', 'ggplot2', 'psych')
new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
invisible(lapply(packages, library, character.only = TRUE))
invisible(gc())

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
####                         0.2 Auxiliary functions                        ####
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# - - - - - - - - - - Network visualization function - - - - - - - - - - - - - -

mplot.ts <- function(xx,
                     titulos_filas = dimnames(xx)[[3]],
                     titulos_columnas = dimnames(xx)[[2]],
                     xlab = "Time") {
  d <- dim(xx)
  
  # Arrange one plot for each combination of rows and columns
  op <- par(mfrow = c(d[3], d[2]),
            mar = c(1, 2.5, 1.5, 0.5),
            oma = c(4, 3, 1, 1))
  on.exit(par(op))
  
  for(j in seq_len(d[3])) {
    for(i in seq_len(d[2])) {
      
      # Plot the time series corresponding to the (i, j) entry
      plot(xx[, i, j], type = "l", lwd = 1, xlab = "", ylab = "", xaxt = "n")
      
      # Display the x-axis only in the bottom row of panels.
      if(j == d[3])
        axis(1)
      
      # Display the y-axis only in the first column of panels
      if(i == 1)
        axis(2)
      
      # Add column labels to the top row of panels
      if(j == 1)
        title(main = titulos_columnas[i], cex.main = 1.2)
      
      # Add row labels to the first column of panels
      if(i == 1)
        mtext(titulos_filas[j], side = 2, line = 2, cex = 0.8, font = 2)
    }
  }
  # Add a common x-axis label to the entire figure
  mtext(xlab, side = 1, outer = TRUE, line = 2)
}


# - - - - - - - - - - - Plot the ACF of matrix-valued data - - - - - - - - - - -

mplot.acf <- function(xx,
                      titulos_filas = dimnames(xx)[[2]],
                      titulos_columnas = dimnames(xx)[[3]]) {
  
  # Extract the underlying data when xx is an S4 object
  if (isS4(xx)) xx <- xx@data
  d <- dim(xx)
  
  # Arrange the ACF plots in a matrix corresponding to the data dimensions
  op <- par(mfrow = c(d[2], d[3]),
            mar = c(0.5, 0.5, 1.5, 0.5),
            oma = c(4, 3, 1, 1))
  on.exit(par(op))
  
  for(i in seq_len(d[2])) {
    for(j in seq_len(d[3])) {
      
      # Compute and plot the autocorrelation function for each series
      acf(xx[, i, j],
          ylim = c(0, 1),
          xaxt = if(i == d[2]) "s" else "n",
          yaxt = if(j == 1) "s" else "n",
          xlab = "",
          ylab = "",
          main = "")
      
      # Add column labels to the top row of panels
      if(i == 1)
        title(main = titulos_columnas[j], cex.main = 1.2, font.main = 2)
      
      # Add row labels to the first column of panels
      if(j == 1)
        mtext(titulos_filas[i], side = 2, line = 2, cex = 0.8, font = 2)
    }
  }
  # Add a common x-axis label to the entire figure
  mtext("Lag", side = 1, outer = TRUE, line = 2)
}


# - - - - - - - Construct the coefficient significance table - - - - - - - - - -

tabla_significancia <- function(coef_mat, sd_mat) {
  
  # Compute the standardized coefficient statistics
  t_mat <- coef_mat / sd_mat
  
  # Compute two-sided p-values under the standard normal approximation
  p_mat <- 2 * (1 - pnorm(abs(t_mat)))
  
  # Convert the coefficient matrices into a long-format data frame
  df <- expand.grid(Row = rownames(coef_mat),
                    Column = colnames(coef_mat))
  
  df$Coef <- as.vector(coef_mat)
  df$SE   <- as.vector(sd_mat)
  df$t    <- as.vector(t_mat)
  df$p    <- as.vector(p_mat)
  
  # Assign significance codes according to conventional p-value thresholds
  df$Sig <- cut(df$p,
                breaks = c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf),
                labels = c("***", "**", "*", ".", ""))
  
  # Order the table by row and column names
  df <- df[order(df$Row, df$Column), ]
  
  # Return the final significance table.
  df[, c("Row", "Column", "Coef", "SE", "t", "p", "Sig")]
}


# - - - - - - - - - - - - Impulse-response function - - - - - - - - - - - - - -

mar_oirf <- function(A1, A2, Sigma, i, j, H = 10) {
  
  d1 <- nrow(A1)
  d2 <- nrow(A2)
  
  # Determine the position of the (i, j) element in vec(E_t) using column-major 
  # ordering
  idx <- d1 * (j - 1) + i
  
  # Extract the covariance vector associated with the selected shock
  col_sig <- Sigma[, idx]
  
  # Normalize the shock to have a magnitude of one standard deviation
  sd_ij <- sqrt(Sigma[idx, idx])
  
  # Initialize the impulse-response array
  irf <- array(NA, dim = c(H + 1, d1, d2))
  
  for (k in 0:H) {
    
    # Compute the k-th powers of the autoregressive coefficient matrices
    # At k = 0, the identity matrix is used
    Ak <- if (k == 0) diag(d1) else A1 %^% k
    Bk <- if (k == 0) diag(d2) else A2 %^% k
    
    # Propagate the selected shock through the MAR system
    # The Kronecker product follows the column-major vectorization convention
    Fk_vec <- kronecker(Bk, Ak) %*% col_sig / sd_ij
    
    # Reshape the vectorized response back into matrix form
    irf[k + 1, , ] <- matrix(Fk_vec, nrow = d1, ncol = d2)
  }
  # Return the impulse-response array.
  irf
}


# - - - - - - - - - - Plot the impulse-response curves - - - - - - - - - - - - -

plot_oirf_panel <- function(irf,
                            capa_labels = NULL,
                            caract_labels = NULL,
                            ref_lines = c(-0.1, 0.1),
                            titulo_shock = "",
                            same_scale = TRUE) {
  
  H1 <- dim(irf)[1]
  d1 <- dim(irf)[2]
  d2 <- dim(irf)[3]
  H <- H1 - 1
  
  # Use a common y-axis range across all panels when requested
  yr <- if (same_scale) range(irf, ref_lines, na.rm = TRUE) else NULL
  
  # Arrange one panel for each matrix element
  op <- par(mfrow = c(d1, d2),
            mar = c(0.3, 0.3, 0.3, 0.5),
            oma = c(3.5, 3, 3, 1),
            xaxs = "i")
  
  for (i in 1:d1) {
    for (j in 1:d2) {
      
      # Determine the y-axis range for the current panel
      y_range <- if (same_scale) {
        yr
      } else {
        range(irf[, i, j], ref_lines)
      }
      
      # Create an empty plot to establish the axes and limits
      plot(0:H,
           irf[, i, j],
           type = "n",
           ylim = y_range,
           xlab = "",
           ylab = "",
           xaxt = if (i == d1) "s" else "n",
           yaxt = if (j == 1) "s" else "n",
           cex.axis = 0.7)
      
      # Add reference lines before plotting the impulse-response curve
      abline(h = ref_lines, lty = 3, col = "gray40")
      abline(h = 0, col = "gray80")
      
      # Plot the impulse-response curve on top of the reference lines
      lines(0:H, irf[, i, j], lwd = 1)
      
      # Add a border around the panel
      box()
      
      # Add the row label to the first column of panels
      if (j == 1) {
        text(par("usr")[1] + 0.05 * diff(par("usr")[1:2]),
             par("usr")[3] + 0.08 * diff(par("usr")[3:4]),
             capa_labels[i],
             adj = 0,
             cex = 1,
             font = 2)
      }
      
      # Add the column label to the top row of panels
      if (i == 1) {
        text(par("usr")[2] - 0.05 * diff(par("usr")[1:2]),
             par("usr")[4] - 0.12 * diff(par("usr")[3:4]),
             caract_labels[j],
             adj = 1,
             cex = 1,
             font = 2)
      }
    }
  }
  
  # Add common axis and figure labels
  mtext("Lag (k)", side = 1, outer = TRUE, line = 2, cex = 0.9)
  mtext(titulo_shock, side = 3, outer = TRUE, line = 1, cex = 1.05, font = 2)
  
  # Restore the previous graphical parameters
  par(op)
}


# - - - - - - - - - - - Residual diagnostic plots - - - - - - - - - - - - - - -

mplot.diag <- function(xx,
                       type = c("acf", "pacf", "qq"),
                       titulos_filas = dimnames(xx)[[2]],
                       titulos_columnas = dimnames(xx)[[3]]) {
  
  type <- match.arg(type)
  
  # Extract the underlying data when xx is an S4 object
  if (isS4(xx)) xx <- xx@data
  d <- dim(xx)
  
  # Arrange one diagnostic plot for each combination of rows and columns
  op <- par(mfrow = c(d[2], d[3]),
            mar = c(0.5, 0.5, 1.5, 0.5),
            oma = c(4, 3, 1, 1))
  on.exit(par(op))
  
  for (i in seq_len(d[2])) {
    for (j in seq_len(d[3])) {
      
      # Extract the time series corresponding to the current matrix element
      serie <- xx[, i, j]
      
      if (type == "acf") {
        
        # Plot the autocorrelation function.
        acf(serie,
            ylim = c(0, 1),
            xaxt = if (i == d[2]) "s" else "n",
            yaxt = if (j == 1) "s" else "n",
            xlab = "",
            ylab = "",
            main = "")
        
      } else if (type == "pacf") {
        
        # Plot the partial autocorrelation function
        pacf(serie,
             ylim = c(-1, 1),
             xaxt = if (i == d[2]) "s" else "n",
             yaxt = if (j == 1) "s" else "n",
             xlab = "",
             ylab = "",
             main = "")
        
      } else if (type == "qq") {
        
        # Create a normal Q-Q plot to assess the distribution of residuals
        qqnorm(serie,
               xaxt = if (i == d[2]) "s" else "n",
               yaxt = if (j == 1) "s" else "n",
               xlab = "",
               ylab = "",
               main = "",
               pch = 19,
               cex = 0.5)
        
        # Add a reference line to the Q-Q plot
        qqline(serie, col = "red", lwd = 1.5)
      }
      
      # Add column labels to the top row of panels
      if (i == 1)
        title(main = titulos_columnas[j], cex.main = 1.2, font.main = 2)
      
      # Add row labels to the first column of panels
      if (j == 1)
        mtext(titulos_filas[i], side = 2, line = 2, cex = 0.8, font = 2)
    }
  }
  # Set the common x-axis label according to the diagnostic plot type
  etiqueta_x <- switch(type,
                       acf  = "Lag",
                       pacf = "Lag",
                       qq   = "Theoretical Quantiles")
  
  # Add the common x-axis label to the entire figure
  mtext(etiqueta_x, side = 1, outer = TRUE, line = 2)
}

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
####                          0.3 Data preparation                          ####
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Load the data
load("icews_tensor_proyecto.RData")

# Number of times, layers and countries
n_tiempos <- dim(tensor_array)[4]
n_capas <- dim(tensor_array)[3]
n_paises <- dim(tensor_array)[1]

nombres_capas <- c("Material-", "Material+", "Verbal-", "Verbal+")
nombres_paises <- dimnames(tensor_array)$Source

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
####                 1. SELECT THE 25 MOST CENTRAL COUNTRIES                ####
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Array to store centrality measures
page_cent <- array(NA,
                   dim = c(n_paises, n_capas, n_tiempos),
                   dimnames = list(nombres_paises, paste0("Layer_",1:n_capas), 
                                   paste0("Time_",1:n_tiempos)))

# Calculate Page Rank centrality for every graph
for(m in 1:n_capas){
  for(t in 1:n_tiempos){
    A <- tensor_array[,,m,t]
    g <- graph_from_adjacency_matrix(A, mode = "directed", weighted = TRUE, diag = FALSE)
    page_cent[,m,t] <- page_rank(g, directed = TRUE, weights = E(g)$weight)$vector
  }
}

# Average centrality measures across time
cent_mean <- apply(page_cent, c(1,2), mean, na.rm = TRUE)

# Identify the top 25 most central countries
top25_layers <- vector("list", n_capas)
for(m in 1:n_capas){
  ord <- order(cent_mean[,m], decreasing = TRUE)
  top25_layers[[m]] <- data.frame(Country = rownames(cent_mean)[ord[1:25]],
                                  MeanEigenCentrality = cent_mean[ord[1:25], m])
}
top25_names <- lapply(top25_layers, \(x) x$Country)

# Countries that are among the top 25 in all layers 
Reduce(intersect, top25_names) 
# There are 19 countries

# - - - - - - - - - - - - - - - - Global top 25 - - - - - - - - - - - - - - - - 
overall_mean <- rowMeans(cent_mean)
ord <- order(overall_mean, decreasing = TRUE)
top25_global <- data.frame(Country = names(overall_mean)[ord[1:25]],
                           MeanEigenCentrality = overall_mean[ord[1:25]])
top25_global 
# Includes the 19 countries above plus 6 additional countries

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
####                          2. STATISTICS MATRICES                        ####
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

paises_final <- top25_global$Country

# Indices of the selected countries
indices <- match(paises_final, nombres_paises)

# Network statistics to compute
estadisticos <- c("densidad",
                  "asortatividad",
                  "transitividad",
                  "reciprocidad",
                  "dist_prom",
                  "fuerza_media")
n_estad <- length(estadisticos)

# Array for the statistics
stats_array <- array(NA, 
                     dim = c(n_capas, length(estadisticos), n_tiempos),
                     dimnames = list(paste0("Layer_",1:n_capas), estadisticos, 
                                     paste0("Time_",1:n_tiempos) ) )

for(m in 1:n_capas){
  for(t in 1:n_tiempos){
    A <- tensor_array[indices, indices, m, t]
    g <- graph_from_adjacency_matrix(A, mode = "directed", weighted = TRUE, 
                                     diag = FALSE)
    comp <- igraph::components(as_undirected(g))
    
    stats_array[m,"densidad",t] <- edge_density(g)
    stats_array[m,"asortatividad",t] <- assortativity_degree(g, directed = TRUE)
    stats_array[m,"transitividad",t] <- transitivity(g, type="global")
    stats_array[m,"reciprocidad",t] <- reciprocity(g)
    stats_array[m,"dist_prom",t] <- mean_distance(g, directed = TRUE, unconnected = TRUE)
    stats_array[m,"fuerza_media",t] <- mean(strength(g, mode = "all"))
  }
  
}
# Re-organize the array to the appropiate dimensions
stats_array <- aperm(stats_array, c(3, 1, 2))

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
####                3. EVALUATE STATIONARITY OF THE MATRICES                ####
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

mplot.ts(stats_array,
         titulos_columnas = nombres_capas,
         titulos_filas = c("Density", "Assortativity", "Transitivity", "Reciprocity", 
                           "Mean distance", "Mean strength"),
         xlab = "Time")

mplot.acf(stats_array, 
          titulos_filas = nombres_capas,
          titulos_columnas = c("Density", "Assortativity", "Transitivity", "Reciprocity", 
                               "Mean distance", "Mean strength")) 
# The network statistics exhibit trends over time

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#####          3.1  Estimate and remove the trend           ####
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Three methods are tested to estimate and remove the trend
#   1. Time series decomposition
#   2. STL
#   3. LOESS

# Empty arrays for trend
tend_decompose <- array(NA, dim = dim(stats_array))
tend_STL_fable <- array(NA, dim = dim(stats_array))
tend_LOESS <- array(NA, dim = dim(stats_array))

# Empty arrays for detrended data
sin_tend_decompose <- array(NA, dim = dim(stats_array))
sin_tend_STL_fable <- array(NA, dim = dim(stats_array))
sin_tend_LOESS <- array(NA, dim = dim(stats_array))

# Estimate the trend using each method

for(i in 1:n_capas){
  for(j in 1:n_estad){
    serie <- ts(stats_array[,i,j], frequency = 12)
    
    # Decompose
    fit <- stats::decompose(serie)$trend
    tend_decompose[,i,j] <- fit
    
    # STL using fable
    df <- data.frame(fecha = yearmonth(seq.Date(from = as.Date("2014-01-01"), 
                                                by = "month", length.out = 132)),
                     value = as.numeric(serie) )
    tsb <- as_tsibble(df, index = fecha)
    fit_fable <- tsb %>%
      model(STL(value ~ trend() + season(window = "periodic")) ) %>%
      components()
    tend_STL_fable[,i,j] <- fit_fable$trend
    
    # LOESS
    lo <- loess(value ~ as.numeric(fecha), data = df, span = 0.75)
    tend_LOESS[,i,j] <- predict(lo)
  }
}

# Detrended series obtained with each method
sin_tend_decompose <- stats_array - tend_decompose
sin_tend_STL_fable <- stats_array - tend_STL_fable
sin_tend_LOESS <- stats_array - tend_LOESS

# - - - - - - - - - - - - - - Detrended series plots - - - - - - - - - - - - - - 
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
legend("topright", legend = c("Decompose detrended", "STL-fable detrended", 
                              "LOESS detrended"),
       col = 1:3, lty = 1, cex = 0.5)
# For the analysis, STL is used to remove the trend

# - - - - - - - - - - - - - Final detrended series plots - - - - - - - - - - - -

mplot.ts(sin_tend_STL_fable, 
         titulos_filas = c("Density", "Assortativity", "Transitivity", "Reciprocity", 
                           "Mean distance", "Mean strength"), 
         titulos_columnas = nombres_capas)

mplot.acf(xx = sin_tend_STL_fable)

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#####                        3.2  PRIOR STANDARDIZATION                     ####
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Prevent the model from being dominated by highly variable features

d1 <- n_capas # layers
d2 <- n_estad  # features

datos_std <- sin_tend_STL_fable
escalas <- numeric(d2)

for (j in 1:d2) {
  # all layers and time points for the given feature
  valores <- as.vector(sin_tend_STL_fable[, , j])  
  escalas[j] <- sd(valores)
  datos_std[, , j] <- sin_tend_STL_fable[, , j] / escalas[j]
}
names(escalas) <- dimnames(sin_tend_STL_fable)[[3]]  
print(round(escalas, 4))

# ACF of the final data
mplot.acf(datos_std, 
          titulos_filas = nombres_capas,
          titulos_columnas = c("Density", "Assortativity", "Transitivity", 
                               "Reciprocity", "Mean distance", "Mean strength")) 

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#####                         3.3 How many lags to use?                     ####
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Fit models with different lag orders
aa <- tenAR.est(datos_std, R = 1, P = 1, method = "LSE")
bb <- tenAR.est(datos_std, R = 1, P = 2, method = "LSE")
cc <- tenAR.est(datos_std, R = 1, P = 3, method = "LSE")
dd <- tenAR.est(datos_std, R = 1, P = 4, method = "LSE")
ee <- tenAR.est(datos_std, R = 1, P = 5, method = "LSE")
ff <- tenAR.est(datos_std, R = 1, P = 6, method = "LSE")

# Construct the BIC comparison table
modelos <- list(P1 = aa, P2 = bb, P3 = cc, P4 = dd, P5 = ee, P6 = ff)

tabla_bic <- data.frame(P   = 1:6,
                        BIC = sapply(modelos, function(m) m$BIC) )
tabla_bic
# The best model, according to the lowest BIC, has p = 1 lag

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
####                              4. MAR ESTIMATION                         ####
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Refit the model using the standardized data
modelo_MAR_std <- aa
mplot.acf(modelo_MAR_std$res)

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#####                      4.1  Coefficient significance                    ####
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Matrix 1 (Layer-to-layer effects)
mat1 <- modelo_MAR_std$A[[1]][[1]][[1]]
sd1  <- modelo_MAR_std$sd[[1]][[1]][[1]]
rownames(mat1) <- rownames(sd1) <- nombres_capas
colnames(mat1) <- colnames(sd1) <- nombres_capas

# Significance table
tabla_significancia(mat1, sd1)

# Matrix 2 (Network features)
mat2 <- modelo_MAR_std$A[[1]][[1]][[2]]
sd2  <- modelo_MAR_std$sd[[1]][[1]][[2]]
rownames(mat2) <- rownames(sd2) <- c("Density", "Assortativity", "Transitivity", 
                                     "Reciprocity", "Mean distance", "Mean strength")
colnames(mat2) <- colnames(sd2) <- c("Density", "Assortativity", "Transitivity", 
                                     "Reciprocity", "Mean distance", "Mean strength")

# Significance table
tabla_significancia(mat2, sd2)

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#####                       4.2 Impulse-response analysis                   ####
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

caract_labels <- c("Density", "Assortativity", "Transitivity", "Reciprocity", 
                   "Mean distance", "Mean strength")

# Check the stability condition
rho_A1 <- max(Mod(eigen(modelo_MAR_std$A[[1]][[1]][[1]])$values))
rho_A2 <- max(Mod(eigen(modelo_MAR_std$A[[1]][[1]][[2]])$values))
cat("rho(A1) =", round(rho_A1, 3), " | rho(A2) =", round(rho_A2, 3),
    " | product =", round(rho_A1 * rho_A2, 3), "\n")

# One-standard-deviation shock in (layer 1, feature 1)
irf_calc <- mar_oirf(modelo_MAR_std$A[[1]][[1]][[1]], 
                     modelo_MAR_std$A[[1]][[1]][[2]], 
                     modelo_MAR_std$Sig, 
                     i = 1, 
                     j = 6, 
                     H = 6)

plot_oirf_panel(irf_calc,
                capa_labels   = nombres_capas,
                caract_labels = caract_labels,
                ref_lines     = c(-0.1, 0.1)#,
                #titulo_shock  = "s1-oIRF: shock of 1 sd in (Layer 1, Feauture 1)"
)

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#####                          4.3 Residual diagnostics                     ####
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Residuals in tensor form
res <- modelo_MAR_std$res
Tt  <- dim(res)[1]

# Diagnostics plots: ACF, PACF and QQ-plot

mplot.diag(res, type = "acf",
           titulos_filas = nombres_capas,
           titulos_columnas = c("Density", "Assortativity", "Transitivity", "Reciprocity", 
                                "Mean distance", "Mean strength")) 
mplot.diag(res, type = "pacf",
           titulos_filas = nombres_capas,
           titulos_columnas = c("Density", "Assortativity", "Transitivity", "Reciprocity", 
                                "Mean distance", "Mean strength")) 

mplot.diag(res, type = "qq",
           titulos_filas = nombres_capas,
           titulos_columnas = c("Density", "Assortativity", "Transitivity", "Reciprocity", 
                                "Mean distance", "Mean strength")) 



# Residuals in matrix form
res_mat <- t(apply(res, 1, as.vector))  # T x N
colnames(res_mat) <- 1:ncol(res_mat)

# Ljung-Box test
lb_pvalores <- sapply(1:ncol(res_mat), function(j) {
  Box.test(res_mat[, j], lag = 10, type = "Ljung-Box")$p.value
})

# FDR correction (Benjamini-Hochberg) for the 24 simultaneous tests
lb_ajustado <- p.adjust(lb_pvalores, method = "BH")

# For each series, count how many of the first 20 lags exceed the 95% significance band
banda <- 1.96 / sqrt(nrow(res_mat))  # approximate ACF significance band

conteo_significativos <- sapply(1:ncol(res_mat), function(j) {
  acf_vals <- acf(res_mat[, j], lag.max = 20, plot = FALSE)$acf[-1]  # exclude lag 0
  sum(abs(acf_vals) > banda)
})

# Summary table with diagnostics
tabla_diagnostico <- data.frame(
  Prueba = c("Univariate Ljung-Box (24 series, FDR correction)",
             "Count of ACF crossings vs. white-noise expectation",
             "Multivariate normality (Mardia, skewness)",
             "Multivariate normality (Mardia, kurtosis)"),
  Resultado = c(
    paste0(sum(lb_ajustado < 0.05), " of 24 series with autocorrelation after FDR adjustment"),
    paste0(sum(conteo_significativos), " observed vs. ", round(24*20*0.05,1), " expected under H0"),
    paste0("Estad. = ", round(mardia_test$skew, 2), ", p = ", round(mardia_test$p.skew, 4)),
    paste0("Estad. = ", round(mardia_test$kurtosis, 2), ", p = ", round(mardia_test$p.kurt, 4))
  )
)
print(tabla_diagnostico, row.names = FALSE)

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
####                                5. SIMULATION                           ####
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Coupling matrix
d1 <- 4
coupling_true <- matrix(0, d1, d1)
coupling_true[1, 1] <- 0.6
coupling_true[3, 1] <- 0.5
coupling_true[3, 3] <- 0.3
coupling_true[4, 3] <- 0.4
coupling_true[4, 4] <- 0.5

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Generator with TWO independent sources of variation:
#  (a) 'coupling' controls the overall density (as before)
#  (b) 'triangle_bias' controls the tendency to close triangles, with its OWN  
#       temporal dynamics, independent of (a)
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

simulate_temporal_network_v2 <- function(n_nodes, Tn, coupling,
                                         noise_sd = 0.05,
                                         triangle_coupling = NULL, # d1 x d1 matrix
                                         burn = 50, threshold = 0.15,
                                         directed = TRUE) {
  
  d1 <- nrow(coupling)
  
  # Moderate autoregressive memory
  if (is.null(triangle_coupling)) triangle_coupling <- diag(0.4, d1)  
  
  Ttot <- Tn + burn
  W    <- array(0, dim = c(Ttot, d1, n_nodes, n_nodes))
  
  # Triangle-closing bias by layer and time
  Bias <- matrix(0.3, Ttot, d1)   
  
  for (l in 1:d1) W[1, l, , ] <- matrix(runif(n_nodes^2, 0, 0.3), n_nodes, n_nodes)
  
  for (t in 2:Ttot) {
    # (b) the triangle-closing bias follows its own recursion, with its OWN noise
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
      noise <- matrix(rnorm(n_nodes^2, 0, noise_sd), n_nodes, n_nodes)
      Wt_raw <- pmin(pmax(base + noise, 0), 1)
      if (!directed) Wt_raw[lower.tri(Wt_raw)] <- t(Wt_raw)[lower.tri(Wt_raw)]
      diag(Wt_raw) <- 0
      
      Wt <- Wt_raw * (Wt_raw > threshold)
      
      # Triangle-closing boost driven by Bias[t,l], independently of density ---
      adj_bin <- (Wt_raw > threshold) * 1
      common_neighbors <- adj_bin %*% t(adj_bin)
      diag(common_neighbors) <- 0
      
      # Force a fraction of node pairs with common neighbors to become connected
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
      
      # ------- Partial reciprocity, independent of the mechanisms above -------
      mask_recip <- matrix(runif(n_nodes^2) < 0.4, n_nodes, n_nodes)
      Wt[mask_recip] <- pmax(Wt[mask_recip], t(Wt)[mask_recip])
      diag(Wt) <- 0
      
      W[t, l, , ] <- Wt
    }
  }
  W[(burn + 1):Ttot, , , ]
}

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#                                 Features
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

compute_network_features_v2 <- function(W) {
  Tn <- dim(W)[1]; d1 <- dim(W)[2]; n <- dim(W)[3]
  feat_names <- c("densidad", "asortatividad", "transitividad", "reciprocidad", 
                  "dist_prom", "fuerza_media")
  d2 <- length(feat_names)
  X <- array(NA, dim = c(Tn, d1, d2), dimnames = list(NULL, NULL, feat_names))
  
  for (t in 1:Tn) {
    for (l in 1:d1) {
      g <- graph_from_adjacency_matrix(W[t, l, , ], mode = "directed",
                                       weighted = TRUE, diag = FALSE)
      X[t, l, "densidad"]      <- edge_density(g)
      X[t, l, "asortatividad"] <- tryCatch(assortativity_degree(g, directed = TRUE), 
                                           error = function(e) 0)
      X[t, l, "transitividad"] <- tryCatch(transitivity(g, type = "global"), 
                                           error = function(e) 0)
      X[t, l, "reciprocidad"]  <- tryCatch(reciprocity(g), error = function(e) 0)
      X[t, l, "dist_prom"]     <- tryCatch(mean_distance(g, directed = TRUE, 
                                                         unconnected = TRUE, weights = NA),
                                           error = function(e) 0)
      X[t, l, "fuerza_media"]  <- mean(strength(g, mode = "all", weights = E(g)$weight))
    }
  }
  X[is.na(X)] <- 0
  X
}

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#         One replicate: simulate, standardize, and fit the MAR model
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

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
  
  # --- Standardize each feature ---
  X_std <- X_sim
  for (j in 1:dim(X_sim)[3]) {
    s <- sd(as.vector(X_sim[, , j]))
    if (s > 0) X_std[, , j] <- X_sim[, , j] / s
  }
  
  # Fit MAR(1) model
  fit <- tryCatch(tenAR.est(X_std, R = 1, P = 1, method = "LSE"), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  
  # A1 and its SE: A[[lag=1]][[term=1]][[mode=1]]
  A1_hat <- fit$A[[1]][[1]][[1]]
  se_A1  <- fit$sd[[1]][[1]][[1]]
  
  # p-value
  p_A1 <- 2 * (1 - pnorm(abs(A1_hat / se_A1)))
  
  verdad_capas <- coupling != 0
  detectado    <- p_A1 < 0.05
  
  data.frame(
    sensibilidad  = mean(detectado[verdad_capas]),
    especificidad = mean(!detectado[!verdad_capas])
  )
}

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#####                             5.1 Monte Carlo                           ####
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

set.seed(123)
B <- 60
resultados_red <- do.call(rbind, lapply(1:B, function(b) run_one_replicate_v2()))

cat("Mean sensibility (recovers true couplings):",
    round(mean(resultados_red$sensibilidad, na.rm = TRUE), 3), "\n")
cat("Mean specificity (does not invent false couplings):",
    round(mean(resultados_red$especificidad, na.rm = TRUE), 3), "\n")

boxplot(resultados_red, main = "Detection of layer coupling\n(generator: dynamic network)")

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#####                       5.2 Simulation scenarios                        ####
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

run_scenario_red <- function(scenario_name, B = 40, coupling, n_nodes = 25, 
                             Tn = 132, noise_sd = 0.05, triangle_coupling = NULL, 
                             threshold = 0.15) {
  
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
      fit <- tenAR.est(X_std, R = 1, P = 1, method = "LSE")
      
      A1_hat <- fit$A[[1]][[1]][[1]]
      se_A1  <- fit$sd[[1]][[1]][[1]]
      p_A1   <- 2 * (1 - pnorm(abs(A1_hat / se_A1)))
      
      verdad    <- coupling != 0
      detectado <- p_A1 < 0.05
      
      data.frame(sensibilidad = mean(detectado[verdad]),
                 especificidad = mean(!detectado[!verdad]))
    }, error = function(e) NULL)
  }
  out <- do.call(rbind, resultados)
  if (is.null(out)) return(NULL)
  out$scenario <- scenario_name
  out
}

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#             Simulation scenarios
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# --- Scenario 1: baseline case ---
esc1 <- run_scenario_red("Baseline", coupling = coupling_true)

# --- Scenario 2: no true coupling ---
coupling_null <- matrix(0, d1, d1)
esc2 <- run_scenario_red("Negative control (no coupling)", coupling = coupling_null)

# --- Scenario 3: weak coupling ---
coupling_debil <- coupling_true
coupling_debil[coupling_debil != 0] <- coupling_debil[coupling_debil != 0] * 0.5
esc3 <- run_scenario_red("Weak coupling (half strength)", coupling = coupling_debil)

# --- Scenario 4: higher noise in edge dynamics ---
esc4 <- run_scenario_red("High noise", coupling = coupling_true, noise_sd = 0.15)

# --- Scenario 5: smaller network ---
esc5 <- run_scenario_red("Small network (n=12)", coupling = coupling_true, n_nodes = 12)

# --- Scenario 6: larger network ---
esc6 <- run_scenario_red("Large network (n=50)", coupling = coupling_true, n_nodes = 50)

# --- Scenario 7: shorter time series ---
esc7 <- run_scenario_red("Short T (T=60)", coupling = coupling_true, Tn = 60)

# --- Scenario 8: triangle-closing bias (second source of variation) ---
tri_coupling_acoplado <- matrix(0, d1, d1)
diag(tri_coupling_acoplado) <- 0.4
tri_coupling_acoplado[3, 1] <- 0.3
esc8 <- run_scenario_red("Coupling also in triangle closure",
                         coupling = coupling_true,
                         triangle_coupling = tri_coupling_acoplado)

# --- Scenario 9: more stringent edge threshold --- 
esc9 <- run_scenario_red("Higher threshold (sparser networks)",
                         coupling = coupling_true, threshold = 0.25)

resultados_todos <- rbind(esc1, esc2, esc3, esc4, esc5, esc6, esc7, esc8, esc9)

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#####                 5.3 Sensitivity and specificity plots                 ####
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Create a numeric identifier for each scenario
niveles <- c(esc1$scenario[1], esc2$scenario[1], esc3$scenario[1], esc4$scenario[1],
             esc5$scenario[1], esc6$scenario[1], esc7$scenario[1], esc8$scenario[1], 
             esc9$scenario[1])
resultados_todos$escenario_id <- factor(resultados_todos$scenario,
                                        levels = niveles,
                                        labels = seq_along(niveles) )

# Specificity
ggplot(resultados_todos,
       aes(x = escenario_id, y = especificidad)) +
  geom_boxplot(fill = "steelblue", alpha = 0.8, outlier.size = 1, staplewidth = 0.4) +
  geom_hline(yintercept = 0.95, linetype = 2, colour = "red") +
  labs(x = "Scenario",
       y = "Specificity") +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0),
    plot.title = element_text(hjust = 0.5)
  )

# Sensitivity
ggplot(resultados_todos,
       aes(x = escenario_id, y = sensibilidad)) +
  geom_boxplot(fill = "grey70", alpha = 0.8, outlier.size = 1, staplewidth = 0.4) +
  geom_hline(yintercept = 0.95, linetype = 2, colour = "red") +
  labs(x = "Scenario",
       y = "Sensitivity") +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0),
    plot.title = element_text(hjust = 0.5)
  )
