
###########################################################################################################
# CODE PART 1: formatting, spatial smoothing, collection of all variables needed for subsequent modeling: #
###########################################################################################################


# list of sample names
list.of.names <- paste0("sample_",
                        LETTERS[1:8])



# matrix where each variable of interest gets assigned corresponding transformation to make more gaussian
# (transformations were chosen after inspecting distributions of said variable in each region)
var_n_tr <- rbind(
  cbind(variable = "red3", transform ="inv.sqrt"),
  cbind(variable = "red7var", transform ="sqrt"),
  
  cbind(variable = "green3", transform =""),
  cbind(variable = "green7var", transform ="sqrt"),
  
  cbind(variable = "blue3", transform ="inv.sqrt"),
  cbind(variable = "blue7var", transform ="sqrt"),
  
  cbind(variable = "area", transform ="sqrt"),
  cbind(variable = "perimeter", transform =""),
  cbind(variable = "orientation", transform =""),
  cbind(variable = "red", transform =""),
  cbind(variable = "blue", transform ="sqrt"),
  cbind(variable = "green", transform ="sqrt"),
  cbind(variable = "variance_red", transform ="log"),
  cbind(variable = "variance_blue", transform ="log"),
  cbind(variable = "variance_green", transform ="log"),
  cbind(variable = "intensity", transform ="log"))




# data frame where all samples get grouped and all cells get binary classifications
big.df <- do.call(rbind, lapply(list.of.names, function(sample.n){
  cat(sample.n, "\n")
  
  name1 <- paste0("https://raw.githubusercontent.com/jaricl/roitumordetect/refs/heads/main/", sample.n, "_part1.txt")
  name2 <- paste0("https://raw.githubusercontent.com/jaricl/roitumordetect/refs/heads/main/", sample.n, "_part2.txt")
  
  part1 <- read.table(file = name1)
  part2 <- read.table(file = name2)
  
  data <- rbind(part1, part2)
  
  class_c <- rep(NA, length(data$class))
  class_c[data$class != "Tumor"] = 0
  class_c[data$class == "Tumor"] = 1
  
  data$class = class_c; rm(class_c)
  return(data %>% mutate(sample = sample.n))
}))

# Changing the naming convention for interpretability (1-7 signifying degrees of neighborhood size)

colnames(big.df) <- gsub("_1_76µ_mean", "1", colnames(big.df))
colnames(big.df) <- gsub("_3_52µm_mean", "2", colnames(big.df))
colnames(big.df) <- gsub("_5_28µm_mean", "3", colnames(big.df))
colnames(big.df) <- gsub("7_04µm_mean", "4", colnames(big.df))
colnames(big.df) <- gsub("_8_80µm_mean", "5", colnames(big.df))
colnames(big.df) <- gsub("_10_56µm_mean", "6", colnames(big.df))
colnames(big.df) <- gsub("_12_32µm_mean", "7", colnames(big.df))
colnames(big.df) <- gsub("_var_12_32µm", "7var", colnames(big.df))
big.df <- big.df %>% rename(X = location_x, Y = location_y)



# Function used to filter out cells that are likely artefacts:
# white, black, gray, or any color too close to the diagonal gray line in the RGB cube
dist.to.gray.line <- function(x){ # matrix with 3 columns as the rgb dimensions
  x1 <- matrix(rep(0, 3*nrow(x)), ncol = 3); x2 <- matrix(rep(255, 3*nrow(x)), ncol = 3)
  x01 <- x - x1
  x02 <- x - x2
  te1 <- cbind(x01[,2]*x02[,3]-x01[,3]*x02[,2],
               x01[,3]*x02[,1]-x01[,1]*x02[,3],
               x01[,1]*x02[,2]-x01[,2]*x02[,1])
  no1 <- x2 - x1
  te <- lapply(1:nrow(te1), function(j) sqrt(sum(te1[j,]^2))) %>% unlist
  no <- lapply(1:nrow(no1), function(j) sqrt(sum(no1[j,]^2))) %>% unlist
  return(te/no)
}


# Inference -----
full.pattern.inference <- function(sample, variable, nx = 25, transform){
  # Read in sample:
  st <- Sys.time()
  p.copy <- big.df[big.df$sample == sample, ]
  
  # make hull around sample
  hull <- geojson_sf(paste0("https://raw.githubusercontent.com/JariCL/ROITumorDetect/refs/heads/main/GEOJSON/",
                            sample, ".geojson"))$geometry[[1]] %>%
    as.matrix %>% as.data.frame %>% rename(X = 1, Y = 2)
  
  yrange <- range(hull$Y)[2]-range(hull$Y)[1]
  xrange <- range(hull$X)[2]-range(hull$X)[1]
  
  min.Y <- min(hull$Y)
  min.X <- min(hull$X)
  
  
  # standardize to have y-range be 0 to 1 (keeping dimensionality)
  
  hull$Y <- hull$Y - min.Y; hull$X <- hull$X - min.X; hull$Y <- hull$Y / yrange; hull$X <- hull$X / yrange
  
  p.copy$Y <- p.copy$Y - min.Y
  p.copy$X <- p.copy$X - min.X
  
  p.copy$Y <- p.copy$Y / yrange
  p.copy$X <- p.copy$X / yrange
  
  # filter artefacts:
  p.copy <- p.copy[!(p.copy$red == 0 & p.copy$green == 0 & p.copy$blue == 0),]
  p.copy <- p.copy %>%
    mutate(distances = dist.to.gray.line(p.copy[,c("red", "green", "blue")] %>% as.matrix)) %>%
    filter(distances > 50) # filter enige vorm van grijswaarden met bepaalde speling 
  
  
  # coordinates of cells:
  coords.big <- p.copy[,c("X","Y")]
  
  if (variable != "intensity"){
    
    #variable transformation:
    if (transform == "log"){
      Z.big <- log(p.copy[, variable] %>% as.numeric() + 0.000001)
    } else if (transform == "inv.log"){
      Z.big <- log(max(p.copy[, variable] %>% as.numeric) - p.copy[, variable] %>% as.numeric + 0.000001)
      
    } else if (transform == "sqrt"){
      Z.big <- sqrt(p.copy[, variable] %>% as.numeric())
      
    } else if (transform == "inv.sqrt"){
      Z.big <- sqrt(max(p.copy[, variable] %>% as.numeric()) - p.copy[, variable]%>% as.numeric())
      
    } else {
      Z.big <- p.copy[, variable]%>% as.numeric()
    }
    
    
    
    Z.big <- as.numeric(Z.big)
    
    # hull into sf object:
    hull <- Polygon(hull)
    hull = Polygons(list(hull),1)
    hull = SpatialPolygons(list(hull))
    hull <- st_as_sf(hull, coords = c("X", "Y"))
    
    # tessellation:
    tess <- st_make_grid(hull, n = nx)
    tess <- st_as_sf(tess)
    
    # making points (cells) into sf:
    points <- st_as_sf(x = coords.big, coords = c("X", "Y"))
    
    # cell that per cell gives index of grid cell it is positioned in:
    which.cell <- as.integer(st_intersects(points, tess))
    
    # centroids of grid cell:
    coords <- st_coordinates(st_centroid(tess))
    
    # gives per grid cell the median of the variable of interest in all cells found there:
    Z <- sapply(1:length(st_geometry(tess)), function(x){
      median(Z.big[which(which.cell == x)])
    })
    
    
    # Filter out the grid cells where value is NA or not applicable:
    coords <- coords[!is.na(Z) & is.finite(Z),]
    tess <- tess[!is.na(Z) & is.finite(Z),]
    Z <- Z[!is.na(Z) & is.finite(Z)]
    
  } else {
    
    # Similar steps as previous chunk, just for log intensity (log of counts/area) instead of specific variable per cell
    
    hull <- Polygon(hull)
    hull = Polygons(list(hull),1)
    hull = SpatialPolygons(list(hull))
    hull <- st_as_sf(hull, coords = c("X", "Y"))
    
    
    tess <- st_make_grid(hull, n = nx)
    tess <- st_as_sf(tess)
    points <- st_as_sf(x = coords.big, coords = c("X", "Y"))
    
    tess <- tess
    which.cell <- as.integer(st_intersects(points, tess))
    coords <- st_coordinates(st_centroid(tess))
    
    # takes log of (counts divided by area)
    Z <- log(sapply(1:length(st_geometry(tess)), function(x){
      a <- which.cell == x
      a <- a[!is.na(a)]
      sum(a)
    }) / st_area(tess))
    
    coords <- coords[!is.na(Z) & is.finite(Z),]
    tess <- tess[!is.na(Z) & is.finite(Z),]
    Z <- Z[!is.na(Z) & is.finite(Z)]
    
  }
  
  
  # Non-informative X covariate for subsequent regression (done such that all spatial information is captured by GP term):
  X=matrix(1,nr=nrow(coords))
  
  
  # Do n chains of MCMC of 10000 iterations, 1/n of last 5000 iterations kept per chain to form 5000 final iterations (using 'combnthin')
  # (ensures proper mixing):
  nchains = 3
  report = 100 # reports completion per 100 iterations
  mc_sp <- combnthin(coords, Z, nchains = nchains, report = report)
  
  
  # MC_SP files stored were saved as follows
  # (not stored in github repository due to excess number of files and size,
  # but are available upon reasonable request).
  # (Any relevant files derived from said files áre stored in the repository)
  # save(mc_sp, file = paste0("mc_sp_", sample, "_", variable, "_nx", nx, ".Rda"))
  # load(file = paste0("mc_sp_", sample, "_", variable, "_nx", nx, ".Rda"))
  
  
  # Inspect posterior distributions of all relevant parameters:  
  begin <- 1; end <- 5000
  print(mcmc_combo(matrix(c(mc_sp$parameters$post_phis, mc_sp$parameters$post_sigma2,
                            mc_sp$parameters$post_tau2, mc_sp$parameters$post_beta),
                          ncol = 4, dimnames = list(1:5000,c("phi", "sigma2", "tau2", "beta")))))
  
  
  # Make more detailed tessellation to do predictive inference of GP and first derivative of GP values 
  grid <- st_make_grid(hull, n = 100)
  grid <- st_as_sf(grid)
  grid.points <- st_centroid(grid)
  grid.points <- st_coordinates(grid.points)
  
  
  
  nburn = 0; niter = 5000
  samples <- (nburn+1):niter
  nbatch <- (length(samples))/2
  
  # Distance matrix for original grid centroids:
  Delta <- as.matrix(dist(coords))
  # Matrix with small error on diagonal to add (ensures ability to take inverse)
  d.factor <- 1e-10*diag(nrow(Delta))
  
  if(is.null(samples)) samples <- (nburn+1):niter
  
  # Distance from any detailed grid cell (column) to any original grid cell (row)
  dist.s0 <- sapply(1:nrow(grid.points),function(y) apply(coords,1,function(x) sqrt(sum((x-grid.points[y,])^2)) ))
  
  # List where each element corresponds to detailed grid cell. In each element, gives a data.frame that per detailed grid cell
  # shows x and y directional difference to original grid cells:
  delta.s0 <- sapply(1:nrow(grid.points),function(y) t(apply(coords,1,function(x) x-grid.points[y,])), simplify = F) # diff grid as coords
  
  # renaming MCMC posterior info:
  chain <- mc_sp
  
  # set covariance types
  cov.type = "matern2"
  
  # Split the posterior distribution into 10 subdistributions
  # (Rather than using each value of the posterior distributions of relevant parameters,
  # 10 subdistributions are made and a summary of said subdistribution is instead used in the following predictive inference)
  ngroups = 10
  grouplist <- split((nburn+1):niter, sample(1:ngroups, size = (niter-nburn), replace = T))
  
  
  # Function that inputs the posterior and all of the spatial information regarding the original grid on which inference was done
  # as well as the more detailed grid to do predictive inference on to perform said predictive inference
  st1 <- Sys.time()
  j <- lapply(1:length(grouplist), f1, grouplist, cov.type = cov.type, dist.s0 = dist.s0,
              delta.s0 = delta.s0, Delta = Delta, chain = chain, grid.points = grid.points, Z = Z, d.factor = d.factor)
  j <- lapply(j, function(x) mutate(x, grad = sqrt(V1^2 + V2^2))) # --> length of gradient vector 
  # output from j:
  # columns in order:: predictions for: gradient direction x (V1), gradient direction y (V2),
  # curvature direction xx (V3), curvature direction xy (V4), curvature direction yy (V5),
  # value of GP or "surface" (V6), value of gradient vector length (grad).
  
  
  
  
  
  # median of gradient strength across the iterations
  grad.est <-
    lapply(7, FUN = function(type) do.call(rbind, lapply(j, function(x) x[, type]))) %>%
    as.data.frame %>% apply(2,median) %>% as.vector
  
  # variance of gradient strength across the iterations
  grad.var <-
    lapply(7, FUN = function(type) do.call(rbind, lapply(j, function(x) x[, type]))) %>%
    as.data.frame %>% apply(2,var) %>% as.vector
  
  # median of GP surface across the iterations
  obs.surf <-
    lapply(6, FUN = function(type) do.call(rbind, lapply(j, function(x) x[, type]))) %>%
    as.data.frame %>% apply(2,median) %>% as.vector
  
  # variance of GP surface across the iterations
  obs.surf.var <-
    lapply(6, FUN = function(type) do.call(rbind, lapply(j, function(x) x[, type]))) %>%
    as.data.frame %>% apply(2,var) %>% as.vector
  
  cov.type <- "matern2"
  
  
  # grad.est, grad.var, obs.surf, and obs.surf.var were saved with this code
  # (not stored in github repository due to excess number of files,
  # but are available upon reasonable request).
  # (Any relevant files derived from said files áre stored in the repository)
  
  # save(grad.est,
  #      file = paste0("grad.est.", sample, '.variable.', variable, '.aggreg.', nx,
  #                    '.cov.type.', cov.type, '.grid.', 100, '.groups.', ngroups, '.Rda'))
  # 
  # save(grad.var,
  #      file = paste0("grad.var.", sample, '.variable.', variable, '.aggreg.', nx,
  #                    '.cov.type.', cov.type, '.grid.', 100, '.groups.', ngroups, '.Rda'))
  # save(obs.surf,
  #      file = paste0("obs.surf.", sample, '.variable.', variable, '.aggreg.', nx,
  #                    '.cov.type.', cov.type, '.grid.', 100,  '.groups.', ngroups,'.Rda'))
  # save(obs.surf.var,
  #      file = paste0("obs.surf.var.", sample, '.variable.', variable, '.aggreg.', nx,
  #                    '.cov.type.', cov.type, '.grid.', 100, '.groups.', ngroups, '.Rda'))
  en <- Sys.time()
  passed <- difftime(as.POSIXct(en), as.POSIXct(st), unit="secs")
  cat("\n --- PASSED: , ", passed, " seconds \n")
}





# Apply the inference function on all variables within all samples 
lapply(list.of.names, function(sample1){
  lapply(1:length(var_n_tr[,1]), function(variable1){
    cat("Sample: ", sample1, "---", var_n_tr[variable1,1], "---", var_n_tr[variable1,2], "\n")
    full.pattern.inference(sample = sample1, nx = 25, variable = var_n_tr[variable1,1], transform = var_n_tr[variable1,2])
  })
})

# Do the same for nx = 10 and 20
lapply(c(10,20), function(nx1){
  lapply(list.of.names, function(sample1){
    lapply(1:length(var_n_tr[,1]), function(variable1){
      cat("Sample: ", sample1, "---", var_n_tr[variable1,1], "---", var_n_tr[variable1,2], "\n")
      full.pattern.inference(sample = sample1, nx = nx1, variable = var_n_tr[variable1,1], transform = var_n_tr[variable1,2])
    })
  })
})







# Following bits of code simply aggregate the real observed values into the detailed grid: -----

make.obs.real <- function(sample){
  st <- Sys.time()
  p.copy <- big.df[big.df$sample == sample, ]
  
  hull <- geojson_sf(paste0("https://raw.githubusercontent.com/JariCL/ROITumorDetect/refs/heads/main/GEOJSON/",
                            sample, ".geojson"))$geometry[[1]] %>%
    as.matrix %>% as.data.frame %>% rename(X = 1, Y = 2)
  
  yrange <- range(hull$Y)[2]-range(hull$Y)[1]
  xrange <- range(hull$X)[2]-range(hull$X)[1]
  
  min.Y <- min(hull$Y)
  min.X <- min(hull$X)
  
  hull$Y <- hull$Y - min.Y; hull$X <- hull$X - min.X; hull$Y <- hull$Y / yrange; hull$X <- hull$X / yrange
  
  p.copy$Y <- p.copy$Y - min.Y
  p.copy$X <- p.copy$X - min.X
  
  p.copy$Y <- p.copy$Y / yrange
  p.copy$X <- p.copy$X / yrange
  
  p.copy <- p.copy[!(p.copy$red == 0 & p.copy$green == 0 & p.copy$blue == 0),]
  p.copy <- p.copy %>% filter(area <= 1000)
  
  p.copy <- p.copy %>%
    mutate(distances = dist.to.gray.line(p.copy[,c("red", "green", "blue")] %>% as.matrix)) %>%
    filter(distances > 50) # filter enige vorm van grijswaarden met bepaalde speling 
  
  classes <- p.copy$class
  
  coords.big <- p.copy[,c("X","Y")]
  
  hull <- Polygon(hull)
  hull = Polygons(list(hull),1)
  hull = SpatialPolygons(list(hull))
  hull <- st_as_sf(hull, coords = c("X", "Y"))
  
  tess <- st_make_grid(hull, n = 100)
  tess <- st_as_sf(tess)
  overlap <- st_intersection(tess, hull)
  points <- st_as_sf(x = coords.big, coords = c("X", "Y"))
  
  grid <- overlap
  
  which.cell <- as.integer(st_intersects(points, overlap))
  
  
  transf.var <- lapply(1:length(var_n_tr[,1]), function(i){
    cat(i)
    variable <- var_n_tr[i, 1]
    transform <- var_n_tr[i, 2]
    
    Z <- as.numeric(p.copy[, variable])
    
    if (transform == "log"){
      Z.big <- log(Z + 0.00000001)
    } else if (transform == "inv.log"){
      Z.big <- log(max(Z) - Z + 0.00000001)
      
    } else if (transform == "sqrt"){
      Z.big <- sqrt(Z)
      
    } else if (transform == "inv.sqrt"){
      Z.big <- sqrt(max(Z) - Z)
      
    } else {
      Z.big <- Z
    }
    Z.big <- as.numeric(Z.big)
    return(Z.big)
  })
  
  
  Z <- do.call(rbind, lapply(1:length(st_geometry(overlap)), function(x){
    cat(x, "-")
    index <- which(which.cell == x)
    
    c(lapply(1:length(transf.var), function(var){
      mean(transf.var[[var]][index])
    }) %>% unlist)
    
  }))
  
  grid <- 100
  lapply(1:length(transf.var), function(var){
    obs_real <- Z[,var]
    variable1 <- var_n_tr[var,1]
    # this stores the obs_real objects:
    # (not stored in github repository due to excess number of files,
    # but are available upon reasonable request).
    # (Any relevant files derived from said files áre stored in the repository)
    # save(obs_real, file = paste0("obs_real_", sample, '.variable.', variable1, '.grid.', grid, '.Rda'))
  })
  
  en <- Sys.time()
  passed <- difftime(as.POSIXct(en), as.POSIXct(st), unit="secs")
  cat("\n --- PASSED: , ", passed, " seconds \n")
}
lapply(list.of.names, make.obs.real)







# Following bits of code define the response variable into the detailed grid: -----
define.grids <- function(sample){
  st <- Sys.time()
  p.copy <- big.df[big.df$sample == sample, ]
  
  hull <- geojson_sf(paste0("https://raw.githubusercontent.com/JariCL/ROITumorDetect/refs/heads/main/GEOJSON/",
                            sample, ".geojson"))$geometry[[1]] %>%
    as.matrix %>% as.data.frame %>% rename(X = 1, Y = 2)
  
  yrange <- range(hull$Y)[2]-range(hull$Y)[1]
  xrange <- range(hull$X)[2]-range(hull$X)[1]
  
  min.Y <- min(hull$Y)
  min.X <- min(hull$X)
  
  hull$Y <- hull$Y - min.Y; hull$X <- hull$X - min.X; hull$Y <- hull$Y / yrange; hull$X <- hull$X / yrange
  
  p.copy$Y <- p.copy$Y - min.Y
  p.copy$X <- p.copy$X - min.X
  
  p.copy$Y <- p.copy$Y / yrange
  p.copy$X <- p.copy$X / yrange
  
  p.copy <- p.copy[!(p.copy$red == 0 & p.copy$green == 0 & p.copy$blue == 0),]
  p.copy <- p.copy %>% filter(area <= 1000)
  
  p.copy <- p.copy %>%
    mutate(distances = dist.to.gray.line(p.copy[,c("red", "green", "blue")] %>% as.matrix)) %>%
    filter(distances > 50) # filter enige vorm van grijswaarden met bepaalde speling 
  
  classes <- p.copy$class
  
  coords.big <- p.copy[,c("X","Y")]
  
  hull <- Polygon(hull)
  hull = Polygons(list(hull),1)
  hull = SpatialPolygons(list(hull))
  hull <- st_as_sf(hull, coords = c("X", "Y"))
  
  tess <- st_make_grid(hull, n = 100)
  tess <- st_as_sf(tess)
  overlap <- st_intersection(tess, hull)
  points <- st_as_sf(x = coords.big, coords = c("X", "Y"))
  
  grid <- overlap
  # This saves the detailed grid sf object, as seen in the github repository
  # save(grid, file = paste0(sample, "_grid_regionpred.Rda"))
  
  en <- Sys.time()
  passed <- difftime(as.POSIXct(en), as.POSIXct(st), unit="secs")
  cat("\n --- PASSED: , ", passed, " seconds \n")
}
lapply(list.of.names, function(sample1){
  define.grids(sample1)
})


# data frame where all samples get grouped and all cells get binary classifications
big.df <- do.call(rbind, lapply(list.of.names, function(sample.n){
  cat(sample.n, "\n")
  
  name1 <- paste0("https://raw.githubusercontent.com/jaricl/roitumordetect/refs/heads/main/", sample.n, "_part1.txt")
  name2 <- paste0("https://raw.githubusercontent.com/jaricl/roitumordetect/refs/heads/main/", sample.n, "_part2.txt")
  
  part1 <- read.table(file = name1)
  part2 <- read.table(file = name2)
  
  data <- rbind(part1, part2)
  
  return(data %>% mutate(sample = sample.n))
}))

# Changing the naming convention for interpretability (1-7 signifying degrees of neighborhood size)

colnames(big.df) <- gsub("_1_76µ_mean", "1", colnames(big.df))
colnames(big.df) <- gsub("_3_52µm_mean", "2", colnames(big.df))
colnames(big.df) <- gsub("_5_28µm_mean", "3", colnames(big.df))
colnames(big.df) <- gsub("7_04µm_mean", "4", colnames(big.df))
colnames(big.df) <- gsub("_8_80µm_mean", "5", colnames(big.df))
colnames(big.df) <- gsub("_10_56µm_mean", "6", colnames(big.df))
colnames(big.df) <- gsub("_12_32µm_mean", "7", colnames(big.df))
colnames(big.df) <- gsub("_var_12_32µm", "7var", colnames(big.df))
big.df <- big.df %>% rename(X = location_x, Y = location_y)

class_c <- big.df$class
class_c[big.df$class %in% c("unregistered")] <- "Unregistered"
class_c[big.df$class %in% c("Necrosis", "Necrose")] <- "Necrosis"
class_c[big.df$class %in% c("Plasma Cells", "Plasmacells", "Plasma cells")] <- "Plasma"
class_c[big.df$class %in% c("Bronchial Epithelium", "Alveolar Epithelium", "Alveolar epithelium")] <- "Epithelium"
class_c[big.df$class %in% c("Lymphocytes", "Lymfocytes")] <- "Lymphocytes"
class_c[big.df$class %in% c("Healthy (Alveoli)", "Healthy (Bronchus)")] <- "Healthy"

big.df$class <- class_c; rm(class_c)


colnames(big.df) <- gsub("_1_76µ_mean", "1", colnames(big.df))
colnames(big.df) <- gsub("_3_52µm_mean", "2", colnames(big.df))
colnames(big.df) <- gsub("_5_28µm_mean", "3", colnames(big.df))
colnames(big.df) <- gsub("7_04µm_mean", "4", colnames(big.df))
colnames(big.df) <- gsub("_8_80µm_mean", "5", colnames(big.df))
colnames(big.df) <- gsub("_10_56µm_mean", "6", colnames(big.df))
colnames(big.df) <- gsub("_12_32µm_mean", "7", colnames(big.df))
colnames(big.df) <- gsub("_var_12_32µm", "7var", colnames(big.df))
big.df <- big.df %>% rename(X = location_x, Y = location_y)

un.classes <- big.df$class %>% unique

makecellcounts <- function(sample){
  p.copy <- big.df[big.df$sample == sample, ]
  
  hull <- geojson_sf(paste0("https://raw.githubusercontent.com/JariCL/ROITumorDetect/refs/heads/main/GEOJSON/",
                            sample, ".geojson"))$geometry[[1]] %>%
    as.matrix %>% as.data.frame %>% rename(X = 1, Y = 2)
  
  yrange <- range(hull$Y)[2]-range(hull$Y)[1]
  xrange <- range(hull$X)[2]-range(hull$X)[1]
  
  min.Y <- min(hull$Y)
  min.X <- min(hull$X)
  
  hull$Y <- hull$Y - min.Y; hull$X <- hull$X - min.X; hull$Y <- hull$Y / yrange; hull$X <- hull$X / yrange
  
  p.copy$Y <- p.copy$Y - min.Y
  p.copy$X <- p.copy$X - min.X
  
  p.copy$Y <- p.copy$Y / yrange
  p.copy$X <- p.copy$X / yrange
  
  p.copy <- p.copy[!(p.copy$red == 0 & p.copy$green == 0 & p.copy$blue == 0),]
  
  p.copy <- p.copy %>%
    mutate(distances = dist.to.gray.line(p.copy[,c("red", "green", "blue")] %>% as.matrix)) %>%
    filter(distances > 50) # filter enige vorm van grijswaarden met bepaalde speling 
  
  hull <- Polygon(hull)
  hull = Polygons(list(hull),1)
  hull = SpatialPolygons(list(hull))
  hull <- st_as_sf(hull, coords = c("X", "Y"))
  
  coords.big <- p.copy[,c("X","Y")]
  
  
  tess <- st_make_grid(hull, n = 100)
  tess <- st_as_sf(tess)
  points <- st_as_sf(x = coords.big, coords = c("X", "Y"))
  
  overlap <- tess
  
  which.cell <- as.integer(st_intersects(points, overlap))
  
  cellcounts <- do.call(rbind, lapply(1:length(st_geometry(overlap)), function(x){
    cat(x, "--")
    index <- which(which.cell == x)
    
    return(unlist(lapply(un.classes, function(class){
      sum(p.copy$class[index] == class)
    })))
  })) %>% as.data.frame
  colnames(cellcounts) <- un.classes
  
  # Saving the counts per grid cell of all of the cell types given a 100x100 tessellation, as seen in the github repository
  # save(cellcounts, file = paste0("cellcounts100.", sample,".Rda"))
  
}

lapply(list.of.names, makecellcounts)







# Following bits of code define all cell type counts (with extra split of alveoli/bronchus) into the detailed grid:-----
# (nx 100, but also nx = 10 --> this is necessary for appendix MH analysis, which is explicitly done on rougher grid)
big.df <- do.call(rbind, lapply(list.of.names, function(sample.n){
  cat(sample.n, "\n")
  
  name1 <- paste0("https://raw.githubusercontent.com/jaricl/roitumordetect/refs/heads/main/", sample.n, "_part1.txt")
  name2 <- paste0("https://raw.githubusercontent.com/jaricl/roitumordetect/refs/heads/main/", sample.n, "_part2.txt")
  
  part1 <- read.table(file = name1)
  part2 <- read.table(file = name2)
  
  data <- rbind(part1, part2)
  
  return(data %>% mutate(sample = sample.n))
}))

class_c <- big.df$class
class_c[big.df$class %in% c("unregistered")] <- "Unregistered"
class_c[big.df$class %in% c("Necrosis", "Necrose")] <- "Necrosis"
class_c[big.df$class %in% c("Plasma Cells", "Plasmacells", "Plasma cells")] <- "Plasma"
class_c[big.df$class %in% c("Alveolar Epithelium", "Alveolar epithelium")] <- "Epithelium (Alveoli)"
class_c[big.df$class %in% c("Bronchial Epithelium")] <- "Epithelium (Bronchus)"
class_c[big.df$class %in% c("Lymphocytes", "Lymfocytes")] <- "Lymphocytes"

big.df$class <- class_c; rm(class_c)


colnames(big.df) <- gsub("_1_76µ_mean", "1", colnames(big.df))
colnames(big.df) <- gsub("_3_52µm_mean", "2", colnames(big.df))
colnames(big.df) <- gsub("_5_28µm_mean", "3", colnames(big.df))
colnames(big.df) <- gsub("7_04µm_mean", "4", colnames(big.df))
colnames(big.df) <- gsub("_8_80µm_mean", "5", colnames(big.df))
colnames(big.df) <- gsub("_10_56µm_mean", "6", colnames(big.df))
colnames(big.df) <- gsub("_12_32µm_mean", "7", colnames(big.df))
colnames(big.df) <- gsub("_var_12_32µm", "7var", colnames(big.df))
big.df <- big.df %>% rename(X = location_x, Y = location_y)

un.classes <- big.df$class %>% unique

makecellcounts_with_extrasplit <- function(sample, nx = 100){
  p.copy <- big.df[big.df$sample == sample, ]
  
  hull <- geojson_sf(paste0("https://raw.githubusercontent.com/JariCL/ROITumorDetect/refs/heads/main/GEOJSON/",
                            sample, ".geojson"))$geometry[[1]] %>%
    as.matrix %>% as.data.frame %>% rename(X = 1, Y = 2)
  
  yrange <- range(hull$Y)[2]-range(hull$Y)[1]
  xrange <- range(hull$X)[2]-range(hull$X)[1]
  
  min.Y <- min(hull$Y)
  min.X <- min(hull$X)
  
  hull$Y <- hull$Y - min.Y; hull$X <- hull$X - min.X; hull$Y <- hull$Y / yrange; hull$X <- hull$X / yrange
  
  p.copy$Y <- p.copy$Y - min.Y
  p.copy$X <- p.copy$X - min.X
  
  p.copy$Y <- p.copy$Y / yrange
  p.copy$X <- p.copy$X / yrange
  
  p.copy <- p.copy[!(p.copy$red == 0 & p.copy$green == 0 & p.copy$blue == 0),]
  
  p.copy <- p.copy %>%
    mutate(distances = dist.to.gray.line(p.copy[,c("red", "green", "blue")] %>% as.matrix)) %>%
    filter(distances > 50)
  hull <- Polygon(hull)
  hull = Polygons(list(hull),1)
  hull = SpatialPolygons(list(hull))
  hull <- st_as_sf(hull, coords = c("X", "Y"))
  
  coords.big <- p.copy[,c("X","Y")]
  
  tess <- st_make_grid(hull, n = nx)
  tess <- st_as_sf(tess)
  points <- st_as_sf(x = coords.big, coords = c("X", "Y"))
  
  overlap <- tess
  
  which.cell <- as.integer(st_intersects(points, overlap))
  
  cellcounts <- do.call(rbind, lapply(1:length(st_geometry(overlap)), function(x){
    cat(x, "--")
    index <- which(which.cell == x)
    
    return(unlist(lapply(un.classes, function(class){
      sum(p.copy$class[index] == class)
    })))
  })) %>% as.data.frame
  colnames(cellcounts) <- un.classes
  
  # Saving the counts per grid cell of all of the cell types given a 100x100 tessellation
  # specifically differentiating alveolar and bronchial types, as seen in the github repository
  # save(cellcounts, file = paste0("cellcounts_extrasplit_", nx, ".", sample,".Rda"))
}

lapply(c(10,100), function(nx1){
  lapply(list.of.names, function(sample1){
    makecellcounts_with_extrasplit(sample1, nx1)
  })
})











# Save per nx type (10, 20 and 25) all relevant variables in all samples in one data.frame
# .txt file ('complete_variables_dataframe_nx_*.txt', as seen in the repository):


make.df.bin <- function(p.name = "", nx = 20){
  load(file = paste0(p.name, "_grid_regionpred.Rda"))
  grid.points <- st_centroid(grid)
  grid.points <- st_coordinates(grid.points)
  
  aggr.df <- data.frame(X = grid.points[, 1], Y = grid.points[, 2])
  cov.type <- "matern2"
  grid <- 100; ngroups = 10
  
  list1 <- lapply(var_n_tr[,1], function(variable1){
    # cat(paste0(variable1, '_' ,c('obs_surf', 'obs_surf_var')))
    load(file = paste0("obs.surf.", p.name, '.variable.', variable1, '.aggreg.', nx,
                       '.cov.type.', cov.type, '.grid.', grid, '.groups.', ngroups, '.Rda'))
    load(file = paste0("obs.surf.var.", p.name, '.variable.', variable1, '.aggreg.', nx,
                       '.cov.type.', cov.type, '.grid.', grid, '.groups.', ngroups, '.Rda'))
    df <- cbind(V1 = obs.surf, V2 = obs.surf.var) %>% as.data.frame
    colnames(df) <- paste0(variable1, '_' ,c('obs_surf', 'obs_surf_var'))
    return(df)
  })
  surf.df <- do.call(cbind, list1)
  aggr.df <- cbind(aggr.df, surf.df)
  
  
  
  list1 <- lapply(var_n_tr[,1], function(variable1){
    load(file = paste0("grad.est.", p.name, '.variable.', variable1, '.aggreg.', nx,
                       '.cov.type.', cov.type, '.grid.', grid, '.groups.', ngroups, '.Rda'))
    load(file = paste0("grad.var.", p.name, '.variable.', variable1, '.aggreg.', nx,
                       '.cov.type.', cov.type, '.grid.', grid, '.groups.', ngroups, '.Rda'))
    df <- cbind(V1 = grad.est, V2 = grad.var) %>% as.data.frame
    colnames(df) <- paste0(variable1, '_' ,c('grad_est', 'grad_var'))
    return(df)
  })
  
  grad.df <- do.call(cbind, list1)
  aggr.df <- cbind(aggr.df, grad.df)
  
  list1 <- lapply(var_n_tr[1:(nrow(var_n_tr)),1], function(variable1){
    load(file = paste0("obs_real_", p.name, '.variable.', variable1, '.grid.', grid, '.Rda'))
    df <- data.frame(V1 = obs_real)
    colnames(df) <- paste0(variable1, '_' ,'obs_real')
    return(df)
  })
  
  grad.df <- do.call(cbind, list1)
  aggr.df <- cbind(aggr.df, grad.df)
  
  load(file = paste0("cellcounts100.", p.name, ".Rda"))
  nt <- cellcounts[,c(1:6,8:11)] %>% apply(1, sum) %>% as.vector
  t <- cellcounts[,7]
  aggr.df <- cbind(aggr.df, t = t, nt = nt) %>% as.data.frame
  
  df <- aggr.df; rm(aggr.df)
  
  df <- na.omit(df)
  
  if (standardize){
    df <- df %>% mutate(across(.cols = everything() & !nt & !t, ~ standardize.clmn(.x)))
  }

  cat(p.name, " done\n")
  df$sample = p.name
  return(df)
}



var_n_tr <- rbind(
  cbind(variable = "red3", transform ="inv.sqrt"),
  cbind(variable = "red7var", transform ="sqrt"),
  
  cbind(variable = "green3", transform =""),
  cbind(variable = "green7var", transform ="sqrt"),
  
  cbind(variable = "blue3", transform ="inv.sqrt"),
  cbind(variable = "blue7var", transform ="sqrt"),
  
  cbind(variable = "area", transform ="sqrt"),
  cbind(variable = "perimeter", transform =""),
  cbind(variable = "orientation", transform =""),
  cbind(variable = "red", transform =""),
  cbind(variable = "blue", transform ="sqrt"),
  cbind(variable = "green", transform ="sqrt"),
  cbind(variable = "variance_red", transform ="log"),
  cbind(variable = "variance_blue", transform ="log"),
  cbind(variable = "variance_green", transform ="log"),
  cbind(variable = "intensity", transform ="log"))



lapply(c(10,20,25), function(nx1){
  cat("NX: ", nx1, "\n")
  temp <- do.call(rbind, lapply(list.of.names, function(pattern.name){
    cat(pattern.name, " - ")
    df <- make.df.bin(pattern.name, nx = nx1)
    colnames(df) <- sub("obs_real", "rm", colnames(df))
    colnames(df) <- sub("obs_surf_var", "sv", colnames(df))
    colnames(df) <- sub("obs_surf", "sm", colnames(df))
    colnames(df) <- sub("grad_est", "gm", colnames(df))
    colnames(df) <- sub("grad_var", "gv", colnames(df))
    df
  }))
  
  # write.table(temp, file = paste0("complete_variables_dataframe_nx_", nx1, ".txt"))
})


