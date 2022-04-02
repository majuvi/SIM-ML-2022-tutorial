library(data.table)
library(magrittr)
library(igraph)
library(RANN)
library(sf)
library(mgcv)
library(stringr)
library(mvtnorm)
library(zoo)

laad_geodata_wfs <- function(region="buurt", year=2020) {
  fn.cache <- sprintf("%s_%d.rds", region, year)
  # Read a local copy if fetched before
  if (file.exists(fn.cache)) {
    geo_data <- readRDS(fn.cache)
  } else {
    print("Info: Possible regions are gemeente, ggdregio, buurt, wijk)")
    # Set WFS web address. You can look this up in QGIS through the PDOK plugin
    # Or http://www.nationaalgeoregister.nl/geonetwork/srv/dut/catalog.search#/metadata/effe1ab0-073d-437c-af13-df5c5e07d6cd
    baseurl <- "https://geodata.nationaalgeoregister.nl/cbsgebiedsindelingen/wfs?"
    wfsrequest <- sprintf("request=GetFeature&outputFormat=json&typeName=cbs_%s_%s_gegeneraliseerd", region, year)
    geo_data <- subset(st_read(str_c(baseurl, wfsrequest)), select=c(statcode, statnaam))
    geo_data$statcode <- as.factor(geo_data$statcode)
    saveRDS(geo_data, fn.cache)
  }
  return(geo_data)
}

DataGeometry <- function(ref_date, which="bu") {
  ref_year <- substr(ref_date, start = 1, stop = 4)
  if (which == "bu") {
    data.geometry <- st_read(sprintf("Data/Geodata/cbs_buurt/cbs_buurt_%s.geojson", ref_year))
    data.geometry <- data.geometry[data.geometry$statnaam != "",]
    colnames(data.geometry) <- c("bu_code", "bu_name", "geometry")
  } else if (which == "wk") {
    data.geometry <- st_read(sprintf("Data/Geodata/cbs_wijk/cbs_wijk_%s.geojson", ref_year))
    data.geometry <- data.geometry[data.geometry$statnaam != "",]
    colnames(data.geometry) <- c("wk_code", "wk_name", "geometry")
  } else {
    data.geometry <- st_read(sprintf("Data/Geodata/cbs_gemeente/cbs_gemeente_%s.geojson", ref_year))
    data.geometry <- data.geometry[data.geometry$statnaam != "",]
    colnames(data.geometry) <- c("gm_code", "gm_name", "geometry")
  }
  return(data.geometry)
}

# Create a list ggd: bu_codes
unique_labels <- function(x) levels(droplevels(unique(x)))
ggd_bucodes  <- function(df) with(df, tapply(bu_code, gg_code, unique_labels))

# THIS SCRIPT CONTAINS ALL FUNCTIONS THAT ARE NEEDED FOR THE ORIGINAL SMAP MODEL

sf2nb <- function(x, sparse = TRUE) {

  # Get centroid coordinates
  x.coords <- x %>% st_centroid %>% st_coordinates

  # Create rook-type adjacency matrix
  x.adj <- st_relate(x, x, pattern = "F***1****", sparse = FALSE)

  # Get connected components
  x.comp <- x.adj %>% graph.adjacency %>% components

  # While the number of subgraphs is > 1, connect subgraphs by connecting closest neighbours
  # Result: spatial neighbours list without islands
  while(x.comp$no > 1) {

    # Split coordinates by subgraph
    x.coords.split <- data.frame(x.coords) %>% split(f = x.comp$membership)

    # Distance matrix between all subgraphs
    dist.subgraph <- matrix(Inf, x.comp$no, x.comp$no)
    for (i in 1:(x.comp$no - 1)) {
      for (j in (i + 1):x.comp$no) {
        # Get distances between all points in x.coords.split[[j]] and nearest point in x.coords.split[[i]]
        # Use nn2 function from RANN package for fast nearest neighbour search
        nn.list <- nn2(
          data  = x.coords.split[[i]],
          query = x.coords.split[[j]],
          k = 1)
        # Return nearest distance between x.coords.split[[i]] and x.coords.split[[j]]
        dist.subgraph[i, j] <- with(nn.list, nn.dists[which.min(nn.dists)])
      }
    }

    # Which two subgraphs are the closest to eachother and should be connected?
    index1 <- which(dist.subgraph == min(dist.subgraph), arr.ind = TRUE)

    # Which nodes between the two subgraphs should be connected?
    nn.list1 <- nn2(
      data  = x.coords.split[[index1[1]]],
      query = x.coords.split[[index1[2]]],
      k = 1)
    nn.list2 <- nn2(
      data   = x.coords.split[[index1[2]]],
      query = x.coords.split[[index1[1]]],
      k = 1)
    index2 <- c(
      with(nn.list1, nn.idx[which.min(nn.dists)]),
      with(nn.list2, nn.idx[which.min(nn.dists)]))

    # Get index number of THE nodes within each subgraph
    x.comp$node.index <- x %>% nrow %>% integer
    for (i in seq_len(x.comp$no)) {
      x.comp <- within(x.comp, node.index[membership == i] <- csize[i] %>% seq_len)
    }

    # These two nodes are to be connected
    add <- with(x.comp, c(
      which(membership == index1[1] & node.index == index2[1]),
      which(membership == index1[2] & node.index == index2[2])))

    # Make the connection. Should be symmetric
    x.adj[add[1], add[2]] <- x.adj[add[2], add[1]] <- TRUE

    # Update connect subgraphs
    x.comp <- x.adj %>% graph.adjacency %>% components

  }

  # Return adjacency list or matrix
  if (sparse) {
    # sparse = TRUE -> list
    return(x.adj %>% apply(MARGIN = 1, FUN = which))
  } else {
    # sparse = FALSE -> matrix
    return(x.adj)
  }

}

# Calculate the neighbourhood of each region
bu_neighbours <- function(bu_codes.sf, bu_codes, buffer_dist=10000) {
  buffer = subset(bu_codes.sf, bu_code %in% bu_codes) %>% st_union %>% st_buffer(dist = buffer_dist, nQuadSegs = 6)
  bu_codes.sf.subset <- bu_codes.sf[st_intersects(buffer, bu_codes.sf %>% st_centroid) %>% unlist,]
  bu_codes.nb.subset <- bu_codes.sf.subset %>% sf2nb #!! Where does this function come from?
  names(bu_codes.nb.subset) <- bu_codes.sf.subset$bu_code %>% droplevels %>% levels
  return(bu_codes.nb.subset)
}


smapmodel <- function(ggd.bu_codes, bu_codes.sf) {

  # ggd.bu_codes             is a list ggd: bu_codes
  # ggd.bu_codes.neighbours  is a list of lists ggd: bu_codes (+buffer): neighours

  ggd.bu_codes.neighbours = list()
  for (ggd in names(ggd.bu_codes)) {
    print(ggd)
    # Given geometry and specified bu_codes, find a buffer and neighbours around the given bu_codes
    ggd.bu_codes.neighbours[[ggd]] <- bu_neighbours(bu_codes.sf, ggd.bu_codes[[ggd]])
  }

  # Construct the spatial information required by the SMAP model without a trained model
  model <- list(ggd.bu_codes = ggd.bu_codes,
                ggd.bu_codes.neighbours = ggd.bu_codes.neighbours,
                ggd.models = NULL)
  attr(model, "class") <- "smapmodel"
  return(model)
}

#fit <- function(x, ...) UseMethod("fit")

fit.smapmodel <- function(model, train, formula, family=binomial(), verbose=F) {

  ggd.bu_codes.neighbours <- model$ggd.bu_codes.neighbours
  ggd.bu_codes            <- model$ggd.bu_codes
  ggd.models              <- list()

  # Loop over ggd regions and fit a separate model to each
  for (ggd in names(ggd.bu_codes)) {
    if (verbose)
      print(ggd)

    bu_code.neighbours <- ggd.bu_codes.neighbours[[ggd]]
    bu_codes           <- ggd.bu_codes[[ggd]]
    bu_codes.buffer    <- names(bu_code.neighbours)

    ggd.train <- train[train$bu_code %in% bu_codes.buffer, ]
    ggd.train$bu_code <- factor(ggd.train$bu_code, levels=bu_codes.buffer)

    formula.base <- deparse(formula)
    formula.smap <- '~ . + s(bu_code, bs = "mrf", k = round(length(bu_code.neighbours)/5), xt = list(nb = bu_code.neighbours))'
    formula.new <- update(as.formula(formula.base), formula.smap)
    ggd.model <- try(bam(
      formula = formula.new,
      data = ggd.train,
      family = family,
      drop.unused.levels = FALSE,
      discrete = TRUE,
      select = TRUE))

    if (any(class(ggd.model) == "try-error")) {
      print(sprintf("Error fitting the model to %s, replacing by nullmodel.", ggd))
      ggd.model <- nullmodel(ggd.train)
    }

    ggd.models[[ggd]] <- ggd.model
  }

  model$ggd.models <- ggd.models
  return(model)
}

#predict <- function(x, ...) UseMethod("predict")

predict.smapmodel <- function(model, test, verbose=F) {

  ggd.bu_codes.neighbours <- model$ggd.bu_codes.neighbours
  ggd.bu_codes            <- model$ggd.bu_codes
  ggd.models              <- model$ggd.models

  if (is.null(ggd.models)) stop("Call train() on the data set first!")

  predictions <- vector(mode="double", length=nrow(test))
  for (ggd in names(ggd.bu_codes)) {
    if (verbose)
      print(ggd)

    bu_code.neighbours <- ggd.bu_codes.neighbours[[ggd]]
    bu_codes           <- ggd.bu_codes[[ggd]]
    bu_codes.buffer    <- names(bu_code.neighbours)

    ggd.test <- test[test$bu_code %in% bu_codes, ]
    ggd.test$bu_code <- factor(ggd.test$bu_code, levels=bu_codes.buffer)

    ggd.pred <- predict(ggd.models[[ggd]], newdata = ggd.test, type = "response")
    predictions[test$bu_code %in% bu_codes] <- ggd.pred
  }
  return(predictions)
}


terms.smapmodel <- function(model, test, verbose=F) {

  ggd.bu_codes.neighbours <- model$ggd.bu_codes.neighbours
  ggd.bu_codes            <- model$ggd.bu_codes
  ggd.models              <- model$ggd.models

  if (is.null(ggd.models)) stop("Call train() on the data set first!")

  terms <- list()
  for (ggd in names(ggd.bu_codes)) {
    if (verbose)
      print(ggd)

    bu_code.neighbours <- ggd.bu_codes.neighbours[[ggd]]
    bu_codes           <- ggd.bu_codes[[ggd]]
    bu_codes.buffer    <- names(bu_code.neighbours)

    ggd.test <- test[test$bu_code %in% bu_codes, ]
    ggd.test$bu_code <- factor(ggd.test$bu_code, levels=bu_codes.buffer)

    ggd.pred <- predict(ggd.models[[ggd]], newdata = ggd.test, type = "terms")
    terms[[ggd]] <- cbind(ggd.test, ggd.pred)
  }
  terms <- rbindlist(terms)
  return(terms)
}



plot.smap <- function(bu_code.geometry, fn.model = "drinker_interpretation_smap.txt", ncol=5) {
  library(gridExtra)

  # smapmodel: terms
  term.values <- read.csv(fn.model)
  term.values$column <- trimws(term.values$column)
  term.values$value <- trimws(term.values$value)
  term.values$by <- trimws(term.values$by)
  # Subplots are saved here
  p <- list()
  p.theme <- theme(plot.title=element_text(size=10),
                   axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_text(angle=90),
                   legend.title=element_text(size=8), legend.text=element_text(size=8), legend.position="right")
  # Plot all interaction terms
  for (col in c('s(age, by = sex,  bs = "ps", k = 10)', 's(age, by = ethnicity,  bs = "ps", k = 10)',
                's(age, by = marital_status, bs = "ps", k = 10)', 's(age, by = education, bs = "ps", k = 10)',
                's(sex, ethnicity,  bs = "re")', 's(sex, marital_status,  bs = "re")', 's(sex, education,  bs = "re")')) {
    term.subset <- subset(term.values, column == col)
    if (col %in% c('s(age, by = sex,  bs = "ps", k = 10)', 's(age, by = ethnicity,  bs = "ps", k = 10)',
                   's(age, by = marital_status, bs = "ps", k = 10)', 's(age, by = education, bs = "ps", k = 10)'))
      term.subset$value <- strtoi(term.subset$value)
    p[[col]] <- ggplot(term.subset, aes(x=value, y=shap, fill=by, color=by)) + geom_point() +
      labs(title=col, y="Term", x="Feature value") + p.theme
  }
  # Plot all feature terms
  for (col in c("hhtype" , "hhsize", "hhincomesource", "hhhomeownership", "hhincome", "hhassets", "oad")) {
    term.subset <- subset(term.values, column == sprintf("s(%s)", col))
    if (col %in% c("hhsize", "hhincome", "hhassets", "oad"))
      term.subset$value <- strtoi(term.subset$value)
    p[[col]] <- ggplot(term.subset, aes(x=value, y=shap)) + geom_point(color="blue") +
      labs(title=col, y="Term", x="Feature value") + p.theme
  }
  # smapmodel: bu_code markov random field
  term.subset <- subset(term.values, column == "s(bu_code)")
  bu_code.subset <- merge(x = bu_code.geometry, y = term.subset, by.x = "bu_code", by.y="value")
  bu_code.subset$quantile <- with(bu_code.subset, cut(shap, quantile(shap, na.rm=T, probs=seq(0,1,length.out=11))))
  p[["bu_code"]] <- ggplot(bu_code.subset) + geom_sf(aes(fill=quantile), color=NA) +
    theme(plot.title=element_text(size=10), legend.position="none") +
    scale_fill_brewer(palette="RdYlGn", direction=-1) + labs(title="bu_code")

  # plot all
  return(grid.arrange(grobs=p, ncol=ncol, gp=gpar(fontsize=16)))
}

ModelFeatures <- function(data.population) {

  # Get only bu_codes with N > 10
  bu_code.n   <- aggregate(gm_code ~ bu_code, data=data.population, FUN=length)
  bu_code.n10 <- bu_code.n[bu_code.n$gm_code >= 10, "bu_code"]
  data.population <- subset(data.population, subset=(bu_code %in% bu_code.n10))
  data.population <- droplevels(data.population)

  # One hot encode the columns
  levels.list <- list()
  contrasts.list <- list()
  cols <- colnames(data.population)
  for (col in cols[3:length(cols)]) {
    if (is.factor(data.population[[col]])) {
      contrasts.list[[col]] <- contrasts(data.population[[col]], contrasts=F)
      levels.list[[col]] <- levels(data.population[[col]])
    }
  }

  # Calculate the mean and covariance of data in each bu_code
  formula    <-  ~ age + sex + ethnicity + marital_status + education +
    hhtype + hhsize + hhhomeownership + hhincomesource + hhincome + hhassets + oad - 1
  bu_code.data <- split(data.population, data.population$bu_code)
  bu_code.mean <- list()
  bu_code.cov  <- list()
  bu_code.n  <- list()
  for (bu_code in names(bu_code.data)) {
    # Dummy encoded model matrix
    X = model.matrix(formula, data=bu_code.data[[bu_code]], contrasts.arg = contrasts.list)
    # stats
    bu_code.mean[[bu_code]] <- colMeans(X)#round(,4)
    bu_code.cov[[bu_code]]  <- as.data.frame(cov(X))
    bu_code.n[[bu_code]]    <- nrow(X)
  }
  data.bu_code <- NULL

  # Each factorial column got assigned to what in the model matrix
  col.index <- attr(X, "assign")
  col.label <- attr(terms(formula), "term.labels")
  col.map   <- col.label[col.index]

  bu_code.n <- data.frame(bu_code=names(bu_code.n), N=unlist(bu_code.n))
  bu_code.mean <- t(data.frame(bu_code.mean))
  bu_code.mean <- cbind(bu_code.n, bu_code.mean)

  bu_code.cov  <- rbindlist(bu_code.cov)
  cov.dims <- rep(length(col.map), nrow(bu_code.mean))
  bu_code.covn <- data.frame(bu_code = rep(bu_code.mean$bu_code, cov.dims),
                             N       = rep(bu_code.mean$N, cov.dims))
  bu_code.cov <- cbind(bu_code.covn, bu_code.cov)

  return(list(bu_code.mean=bu_code.mean, bu_code.cov=bu_code.cov))
}

ModelTerms <- function(data.population, bu_code.geometry) {
  # Train on observed labels
  train  <- data.population[!is.na(data.population$y),]

  # SMAP is fitted separately to each ggd region using spatial information
  ggd.bu_codes <- ggd_bucodes(data.population)
  model <- smapmodel(ggd.bu_codes, bu_code.geometry)
  formula = y ~
    s(age, by = sex,  bs = "ps", k = 10) +
    s(age, by = ethnicity,  bs = "ps", k = 10) +
    s(age, by = marital_status, bs = "ps", k = 10) +
    s(age, by = education, bs = "ps", k = 10) +
    s(sex, ethnicity,  bs = "re") +
    s(sex, marital_status, bs = "re") +
    s(sex, education, bs = "re") +
    s(hhtype, bs = "re") +
    s(hhsize, bs = "ps", k = 5) +
    s(hhincomesource, bs = "re") +
    s(hhhomeownership, bs = "re") +
    s(hhincome, bs = "ps", k = 10) +
    s(hhassets, bs = "ps", k = 10) +
    s(oad, bs = "ps", k = 10)
  model <- fit.smapmodel(model, train, formula)

  #  CALCULATE TERMS FOR SMAPMODEL ("DRINKER" HEALTH INDICATOR)
  predi <- terms.smapmodel(model, train)

  # Calculate the mean effect over different ggd models
  shap.values <- list()
  # age x sex
  term                <- 's(age, by = sex,  bs = "ps", k = 10)'
  predi$shap          <- rowSums(predi[,c('s(age):sexman', 's(age):sexwoman')])
  shap.df             <- predi[,.(shap=mean(shap), N=.N), by=.(sex,age)]
  colnames(shap.df)   <- c("by", "value", "shap", "n")
  shap.values[[term]] <- shap.df
  # age x ethnicity
  term                <- 's(age, by = ethnicity,  bs = "ps", k = 10)'
  predi$shap          <- rowSums(predi[,c("s(age):ethnicitynetherlands", "s(age):ethnicitymarokko", "s(age):ethnicityturkey",
                                          "s(age):ethnicitysuriname", "s(age):ethnicityantille", "s(age):ethnicityother_nonwestern",
                                          "s(age):ethnicityother_western")])
  shap.df             <- predi[,.(shap=mean(shap), N=.N), by=.(ethnicity,age)]
  colnames(shap.df)   <- c("by", "value", "shap", "n")
  shap.values[[term]] <- shap.df
  # age x marital_status
  term                <- 's(age, by = marital_status, bs = "ps", k = 10)'
  predi$shap          <- rowSums(predi[,c("s(age):marital_statusmarried", "s(age):marital_statussingle",
                                          "s(age):marital_statusdivorced", "s(age):marital_statuswidowed")])
  shap.df             <- predi[,.(shap=mean(shap), N=.N), by=.(marital_status,age)]
  colnames(shap.df)   <- c("by", "value", "shap", "n")
  shap.values[[term]] <- shap.df
  # age x education
  term                <- 's(age, by = education, bs = "ps", k = 10)'
  predi$shap          <- rowSums(predi[,c("s(age):educationbasis", "s(age):educationvmbo_bk", "s(age):educationvmbo_gt",
                                          "s(age):educationmbo_23", "s(age):educationmbo_4", "s(age):educationhavo_vwo",
                                          "s(age):educationhbo_wo_bc", "s(age):educationhbo_wo_ma")])
  shap.df             <- predi[,.(shap=mean(shap), N=.N), by=.(education,age)]
  colnames(shap.df)   <- c("by", "value", "shap", "n")
  shap.values[[term]] <- shap.df
  # sex x ethnicity
  term                <- 's(sex, ethnicity,  bs = "re")'
  predi$shap          <- predi[["s(sex,ethnicity)"]]
  shap.df             <- predi[,.(shap=mean(shap), N=.N), by=.(sex,ethnicity)]
  colnames(shap.df)   <- c("by", "value", "shap", "n")
  shap.values[[term]] <- shap.df
  # sex x marital_status
  term                <- 's(sex, marital_status,  bs = "re")'
  predi$shap          <- predi[["s(sex,marital_status)"]]
  shap.df             <- predi[,.(shap=mean(shap), N=.N), by=.(sex,marital_status)]
  colnames(shap.df)   <- c("by", "value", "shap", "n")
  shap.values[[term]] <- shap.df
  # sex x education
  term                <- 's(sex, education,  bs = "re")'
  predi$shap          <- predi[["s(sex,education)"]]
  shap.df             <- predi[,.(shap=mean(shap), N=.N), by=.(sex,education)]
  colnames(shap.df)   <- c("by", "value", "shap", "n")
  shap.values[[term]] <- shap.df
  # non-interactions
  for (col in c("hhtype" , "hhsize", "hhincomesource", "hhhomeownership", "hhincome", "hhassets", "oad", "bu_code")) {
    term <- sprintf("s(%s)", col)
    predi$shap          <- predi[[term]]
    shap.df             <- predi[,.(shap=mean(shap), N=.N), by=.(get(col))]
    colnames(shap.df)   <- c("value", "shap", "n")
    shap.values[[term]] <- shap.df
  }
  # Save
  shap.values <- rbindlist(shap.values, idcol = "column", fill=T)
  shap.values$shap <- round(shap.values$shap, 4)
  shap.values$value <- as.character(shap.values$value)
  shap.values[(shap.values$n < 10), c("shap", "n")] <- NA
  return(shap.values)
}

DataGenerate <- function(bu_code.geometry,
                         fn.mean  = "amsterdam_bu_mean.txt",
                         fn.cov   = "amsterdam_bu_cov.txt",
                         fn.model = "drinker_interpretation_smap.txt",
                         p.answer = 0.04, seed = 42) {
  set.seed(seed)

  # Feature distribution
  print("Reading statistics")
  bu_code.mean <- read.csv(fn.mean)
  bu_code.cov <- read.csv(fn.cov)

  col.map <-  c('bu_code',
                'age','sex','sex','ethnicity','ethnicity','ethnicity','ethnicity','ethnicity','ethnicity','ethnicity',
                'marital_status','marital_status','marital_status','marital_status',
                'education','education','education','education','education','education','education','education',
                'hhtype','hhtype','hhtype','hhtype','hhtype','hhtype','hhtype','hhsize',
                'hhhomeownership','hhhomeownership','hhhomeownership',
                'hhincomesource','hhincomesource','hhincomesource','hhincomesource','hhincomesource',
                'hhincomesource','hhincomesource','hhincomesource','hhincomesource','hhincomesource',
                'hhincome','hhassets','oad')

  # Generate data in each bu_code
  print("Generating model matrix")
  data.bu_code <- list()
  bu_codes <- bu_code.mean$bu_code
  for (bu_code in bu_codes) {
    temp.n    <- bu_code.mean[bu_code.mean$bu_code == bu_code,2]
    temp.mean <- as.matrix(bu_code.mean[bu_code.mean$bu_code == bu_code,3:ncol(bu_code.mean)])
    temp.cov  <- as.matrix(bu_code.cov[bu_code.cov$bu_code == bu_code,3:ncol(bu_code.cov)])
    data.bu_code[[bu_code]] <- as.data.frame(rmvnorm(temp.n, temp.mean, temp.cov))
  }
  data.bu_code <- rbindlist(data.bu_code, idcol="bu_code")
  colnames(data.bu_code) <- col.map

  print("Transforming to data frame")
  # Transform columns to valid numerical
  clip.between <- function(y, val.min, val.max) {
    x <- copy(y)
    x[x < val.min] <- val.min
    x[x > val.max] <- val.max
    return(x)
  }
  data          <- data.bu_code[,c("bu_code")]
  data$age      <- round(clip.between(data.bu_code$age, 18, 105),0)
  data$hhsize   <- round(clip.between(data.bu_code$hhsize, 1, 10),0)
  data$hhincome <- round(clip.between(data.bu_code$hhincome, 1, 100),0)
  data$hhassets <- round(clip.between(data.bu_code$hhassets, 1, 100),0)
  data$oad      <- round(clip.between(data.bu_code$oad, 5, 100),0)
  # Transform columns to categorical
  col.concat <- colnames(bu_code.mean)[3:ncol(bu_code.mean)]
  col        <- col.map[2:length(col.map)]
  label <- substr(col.concat, nchar(col)+1, nchar(col)+nchar(col.concat))
  levels.list <- split(label, col)
  for (col in c("sex", "ethnicity", "marital_status", "education", "hhtype", "hhhomeownership", "hhincomesource")) {
    col.sampled <- apply(subset(data.bu_code, select= col.map == col), 1, which.max)
    data[[col]] <- levels.list[[col]][col.sampled]
  }
  data <- as.data.frame(lapply(data, as.factor))

  # SMAP model specification (X => y)
  # =================================
  print("Sampling predicted labels")
  ffill <- function(values, keys, keys.fill=NULL) {
    x <- setNames(values, keys)
    if (!is.null(keys.fill))
      x <- setNames(x[keys.fill], keys.fill)
    x <- na.locf(na.locf(x[order(as.numeric(names(x)))], na.rm=F), fromLast=T)
    return(x)
  }

  # smapmodel: terms
  shap.values <- read.csv(fn.model)
  shap.values$column <- trimws(shap.values$column)
  shap.values$value <- trimws(shap.values$value)
  shap.values$by <- trimws(shap.values$by)
  # Columns with interactions
  ages <- as.character(18:105)
  shap.subset <- subset(shap.values, column == 's(age, by = sex,  bs = "ps", k = 10)')
  shap.sex.age <- lapply(split(shap.subset, shap.subset$by), function(df) ffill(df$shap, df$value, keys.fill=ages))
  shap.subset <- subset(shap.values, column == 's(age, by = ethnicity,  bs = "ps", k = 10)')
  shap.ethnicity.age <- lapply(split(shap.subset, shap.subset$by), function(df) ffill(df$shap, df$value, keys.fill=ages))
  shap.subset <- subset(shap.values, column == 's(age, by = marital_status, bs = "ps", k = 10)')
  shap.marital_status.age <- lapply(split(shap.subset, shap.subset$by), function(df) ffill(df$shap, df$value, keys.fill=ages))
  shap.subset <- subset(shap.values, column == 's(age, by = education, bs = "ps", k = 10)')
  shap.education.age <- lapply(split(shap.subset, shap.subset$by), function(df) ffill(df$shap, df$value, keys.fill=ages))
  shap.subset <- subset(shap.values, column == 's(sex, ethnicity,  bs = "re")')
  shap.sex.ethnicity <- lapply(split(shap.subset, shap.subset$by), function(df) setNames(df$shap, df$value))
  shap.subset <- subset(shap.values, column == 's(sex, marital_status,  bs = "re")')
  shap.sex.marital_status <- lapply(split(shap.subset, shap.subset$by), function(df) setNames(df$shap, df$value))
  shap.subset <- subset(shap.values, column == 's(sex, education,  bs = "re")')
  shap.sex.education <- lapply(split(shap.subset, shap.subset$by), function(df) setNames(df$shap, df$value))
  # Columns without interactions
  shap.subset <- subset(shap.values, column == 's(hhsize)')
  shap.hhsize <- ffill(shap.subset$shap, shap.subset$value)
  shap.subset <- subset(shap.values, column == 's(hhincome)')
  shap.hhincome <- ffill(shap.subset$shap, shap.subset$value)
  shap.subset <- subset(shap.values, column == 's(hhassets)')
  shap.hhassets <- ffill(shap.subset$shap, shap.subset$value)
  shap.subset <- subset(shap.values, column == 's(oad)')
  shap.oad <- ffill(shap.subset$shap, shap.subset$value)
  shap.subset <- subset(shap.values, column == 's(hhtype)')
  shap.hhtype <- setNames(shap.subset$shap, shap.subset$value)
  shap.subset <- subset(shap.values, column == 's(hhincomesource)')
  shap.hhincomesource <- setNames(shap.subset$shap, shap.subset$value)
  shap.subset <- subset(shap.values, column == 's(hhhomeownership)')
  shap.hhhomeownership <- setNames(shap.subset$shap, shap.subset$value)
  shap.subset <- subset(shap.values, column == 's(bu_code)')
  # Spatial effect, since many don't have 10 observations in training data, replace effect by 0
  shap.bu_code <- ffill(shap.subset$shap, shap.subset$value, keys.fill=unique(data$bu_code))
  shap.bu_code[is.na(shap.bu_code)] <- 0.0

  # Calculate predicted value for every observation
  shap.calculate <- function(row) {
    val <- shap.sex.age[[row[["sex"]]]][[row[["age"]]]] +
      shap.ethnicity.age[[row[["ethnicity"]]]][[row[["age"]]]] +
      shap.marital_status.age[[row[["marital_status"]]]][[row[["age"]]]] +
      shap.education.age[[row[["education"]]]][[row[["age"]]]] +
      shap.sex.ethnicity[[row[["sex"]]]][[row[["ethnicity"]]]] +
      shap.sex.marital_status[[row[["sex"]]]][[row[["marital_status"]]]] +
      shap.sex.education[[row[["sex"]]]][[row[["education"]]]] +
      shap.hhtype[[row[["hhtype"]]]] +
      shap.hhsize[[row[["hhsize"]]]] +
      shap.hhincomesource[[row[["hhincomesource"]]]] +
      shap.hhhomeownership[[row[["hhhomeownership"]]]] +
      shap.hhincome[[row[["hhincome"]]]] +
      shap.hhassets[[row[["hhassets"]]]] +
      shap.oad[[row[["oad"]]]] +
      shap.bu_code[[row[["bu_code"]]]]
    return(val)
  }
  # Calculate outcome from the simulated model
  log.odds    <- apply(data, 1, shap.calculate)
  probability <- 1/(1+exp(-log.odds))

  # to numerical
  data$age <- as.numeric(as.character(data$age))
  data$hhsize <- as.numeric(as.character(data$hhsize))
  data$hhincome <- as.numeric(as.character(data$hhincome))
  data$hhassets <- as.numeric(as.character(data$hhassets))
  data$oad <- as.numeric(as.character(data$oad))
  # coordinates
  bu_code.coords   <- cbind(bu_code=bu_code.geometry$bu_code,
                            as.data.frame(st_coordinates(st_centroid(bu_code.geometry$geometry))))
  data = merge(data, bu_code.coords)

  # This is the 'true' answers
  data$p <- probability
  data$y <- rbinom(length(probability), 1, probability)

  ## This is the observed 'survey' answers
  #data$y.survey <- NA
  #answers <- (rbinom(nrow(data), 1, p.answer) == 1)
  #data[answers, "y.survey"] <- data[answers, "y"]

  setDT(data)
  return(data)
}

