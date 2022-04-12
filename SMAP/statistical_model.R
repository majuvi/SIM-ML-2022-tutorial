library(igraph)
library(RANN)
library(sf)
library(mgcv)
library(gridExtra)

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

plot.smap <- function(bu_code.geometry, fn.model = "drinker_interpretation_smap.txt", ncol=5) {
  
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
