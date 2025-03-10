library(treesurgeon)
library(paleotree)
library(rphenoscate)
library(castor)

cols1 <- c("#CC79A7", "#0072B2", "#009E73", "#D55E00")
cols2 <- c("#999999", "#CC79A7", "#0072B2", "#009E73", "#D55E00")
cols3 <- c("#CC79A7", "#a3d3ee", "#00a3ff", "#00466e", "#009E73", "#ffdf75", "#ffc300", "#9c7800")
cols4 <- c("#999999", "#CC79A7", "#a3d3ee", "#00a3ff", "#00466e", "#009E73", "#ffdf75", "#ffc300", "#9c7800")


# set wd
#wd <- "/home/jk12521/Dropbox/Work/Elenora_sponges"
setwd(wd)

setwd("trees/")
mol_trees_names <- list.files()
mol_trees <- list()
for(i in 1:length(mol_trees_names)){
  mol_trees[[i]] <- ladderize(read.nexus(mol_trees_names[[i]]))
}
mol_trees <- .compressTipLabel(mol_trees)
class(mol_trees) <- "multiPhylo"

mol_con <- mol_trees[[1]]
mol_con$edge.length <- (mol_trees[[1]]$edge.length + mol_trees[[2]]$edge.length + mol_trees[[3]]$edge.length + mol_trees[[4]]$edge.length)/4
setwd(wd)

# load fossil trees
trees <- read.nexus("trees_fossils.nex")

# get ages
ages <- read.csv("ages.csv", row.names = 1)

#get tip states
tip_states <- read_nexdat("Sponge_tip_states.nex")
tip_states <- t(as.data.frame(tip_states))

# Loop through models

modelA <- list() # list to store model A rate variations
x <- 1

for(i in c("ER", "ARD")){
    calc <- get_transition_index_matrix(Nstates = 2, rate_model =i)$index_matrix
    rownames(calc) <- c("c*", "c")
    colnames(calc) <- c("c*", "c")
    for(j in  c("ER", "ARD")){
        silica <- get_transition_index_matrix(Nstates = 2, rate_model =j)$index_matrix
        silica[which(calc > 0)] <- silica[which(calc > 0)] + max(calc) 
        rownames(silica) <- c("s*", "s")
        colnames(silica) <- c("s*", "s")
        modelA[[x]] <- amaSMM(calc, silica, diag.as = 0)
        x <- x+1
    }
}
    
modelB <- list() # list to store model A rate variations
x <- 1

for(i in c("ER", "ARD")){
   spic <- get_transition_index_matrix(Nstates = 2, rate_model =i)$index_matrix 
   rownames(spic) <- c("A", "P")
   colnames(spic) <- c("A", "P")
   for(j in modelA){
        model_temp <- j
        model_temp[which(model_temp > 0)] <- model_temp[which(model_temp > 0)] + max(spic)
        modelB[[x]] <- amaED(spic, model_temp, phi = c(1, 1, 1, 0), diag.as = 0)
        x <- x+1
   }
}


modelC <- list() # list to store model A rate variations
x <- 1

for(i in c("ER", "ARD")){
    calc <- get_transition_index_matrix(Nstates = 2, rate_model =i)$index_matrix
    rownames(calc) <- c("c*", "c")
    colnames(calc) <- c("c*", "c")
    for(j in  c("ER", "ARD")){
        silica <- get_transition_index_matrix(Nstates = 2, rate_model =j)$index_matrix
        silica[which(calc > 0)] <- silica[which(calc > 0)] + max(calc) 
        rownames(silica) <- c("s*", "s")
        colnames(silica) <- c("s*", "s")
        for(k in  c("ER", "ARD")){
            pw <- get_transition_index_matrix(Nstates = 3, rate_model =k)$index_matrix
            pw[which(silica > 0)] <- pw[which(silica > 0)] + max(silica) 
            rownames(pw) <- c("1", "2", "3")
            colnames(pw) <- c("1", "2", "3")
            Spw <- amaED(silica, pw, diag.as = 0)
            modelC[[x]] <- amaSMM(calc, Spw, diag.as = 0)
            x <- x+1
        }
    }
}


modelD <- list() # list to store model A rate variations
x <- 1

for(i in c("ER", "ARD")){
    spic <- get_transition_index_matrix(Nstates = 2, rate_model =i)$index_matrix 
    rownames(spic) <- c("A", "P")
    colnames(spic) <- c("A", "P")
    for(j in modelC){
        model_temp <- j
        model_temp[which(model_temp > 0)] <- model_temp[which(model_temp > 0)] + max(spic)
        modelD[[x]] <- amaED(spic, model_temp, phi = c(1, 1, 1, 1, 1, 0, 0, 0), diag.as = 0)
        x <- x+1
    }
}

######################################## Model A run #################################################

setwd(wd)
res_modelA <- list()
df_res_modelA <- data.frame()

x <- 1
for(i in c("tree_fossils", "tree_fossils_Wang", "tree_no_fossils")){
    tp <- get_tip_priors(tip_states)
    colnames(tp[[1]]) <- colnames(modelA[[1]])
    setwd("trees/")
    mol_trees_names <- list.files()
    for(j in mol_trees_names){
        mol_tree_j <- read.nexus(j)
        root_age <- max(node.depth.edgelength(mol_tree_j))
        if(i == "tree_no_fossils"){
            timetree <- ladderize(mol_tree_j)
            timetree$edge.length <- timetree$edge.length * 100
        } else {
            if(i == "tree_fossils"){
                tree <- ladderize(trees[[1]])
            }
            if(i == "tree_fossils_Wang"){
                tree <- ladderize(trees[[2]])
            }
            node_ages <- rep(NA, Nnode(tree))
            for(k in 1:Nnode(mol_tree_j)){
                tips <- mol_tree_j$tip.label[unlist(Descendants(mol_tree_j, node = k + Ntip(mol_tree_j), type = "tips"))]
                node_x <- mrca.phylo(tree, node = tips)
                node_ages[[node_x - Ntip(tree)]] <- (root_age - (node.depth.edgelength(mol_tree_j)[[k + Ntip(mol_tree_j)]])) * 100
            }
            timetree <- timePaleoPhy(tree, timeData = ages, type = "equal", node.mins = node_ages)
        }
        timetree_scaled <- timetree
        timetree_scaled$edge.length <- timetree_scaled$edge.length / max(node.depth.edgelength(timetree))
        for(z in 1:length(modelA)){
            if(i == "tree_no_fossils"){
              fossils <- rownames(ages)[which(ages[,1] >0)]
              tp <- lapply(tp, function(x) x[!(row.names(x) %in% fossils),])
            }
            tp <- lapply(tp, function(x) x[match(timetree$tip.label, rownames(x)),]) #rearrange tip prior order. 
            fit_SMM <- fitMk(timetree_scaled, x = tp[[1]], model = modelA[[z]], pi="fitzjohn")
            res_modelA[[x]] <- fit_SMM
            newrow <- data.frame(tree = i, mol_ages = j, k = max(modelA[[z]]), logL = fit_SMM$logLik)
            df_res_modelA <- rbind(df_res_modelA, newrow)
            x <- x + 1
        }
    }
    setwd(wd)
}



######################################## Model B run #################################################

setwd(wd)
res_modelB <- list()
df_res_modelB <- data.frame()

x <- 1
for(i in c("tree_fossils", "tree_fossils_Wang", "tree_no_fossils")){
  tp <- get_tip_priors(tip_states)
  colnames(tp[[2]]) <- colnames(modelB[[1]])
  setwd("trees/")
  mol_trees_names <- list.files()
  for(j in mol_trees_names){
    mol_tree_j <- read.nexus(j)
    root_age <- max(node.depth.edgelength(mol_tree_j))
    if(i == "tree_no_fossils"){
      timetree <- ladderize(mol_tree_j)
      timetree$edge.length <- timetree$edge.length * 100
    } else {
      if(i == "tree_fossils"){
        tree <- ladderize(trees[[1]])
      }
      if(i == "tree_fossils_Wang"){
        tree <- ladderize(trees[[2]])
      }
      node_ages <- rep(NA, Nnode(tree))
      for(k in 1:Nnode(mol_tree_j)){
        tips <- mol_tree_j$tip.label[unlist(Descendants(mol_tree_j, node = k + Ntip(mol_tree_j), type = "tips"))]
        node_x <- mrca.phylo(tree, node = tips)
        node_ages[[node_x - Ntip(tree)]] <- (root_age - (node.depth.edgelength(mol_tree_j)[[k + Ntip(mol_tree_j)]])) * 100
      }
      timetree <- timePaleoPhy(tree, timeData = ages, type = "equal", node.mins = node_ages)
    }
    timetree_scaled <- timetree
    timetree_scaled$edge.length <- timetree_scaled$edge.length / max(node.depth.edgelength(timetree))
    for(z in 1:length(modelC)){
      if(i == "tree_no_fossils"){
        fossils <- rownames(ages)[which(ages[,1] >0)]
        tp <- lapply(tp, function(x) x[!(row.names(x) %in% fossils),])
      }
      tp <- lapply(tp, function(x) x[match(timetree$tip.label, rownames(x)),]) #rearrange tip prior order. 
      fit_SMM <- fitMk(timetree_scaled, x = tp[[2]], model = modelB[[z]], pi="fitzjohn")
      res_modelB[[x]] <- fit_SMM
      newrow <- data.frame(tree = i, mol_ages = j, k = max(modelB[[z]]), logL = fit_SMM$logLik)
      df_res_modelB <- rbind(df_res_modelB, newrow)
      x <- x + 1
    }
  }
  setwd(wd)
}


######################################## Model C run #################################################

setwd(wd)
res_modelC <- list()
df_res_modelC <- data.frame()

x <- 1
for(i in c("tree_fossils", "tree_fossils_Wang", "tree_no_fossils")){
    tp <- get_tip_priors(tip_states)
    colnames(tp[[5]]) <- colnames(modelC[[1]])
    setwd("trees/")
    mol_trees_names <- list.files()
    for(j in mol_trees_names){
        mol_tree_j <- read.nexus(j)
        root_age <- max(node.depth.edgelength(mol_tree_j))
        if(i == "tree_no_fossils"){
            timetree <- ladderize(mol_tree_j)
            timetree$edge.length <- timetree$edge.length * 100
        } else {
            if(i == "tree_fossils"){
                tree <- ladderize(trees[[1]])
            }
            if(i == "tree_fossils_Wang"){
                tree <- ladderize(trees[[2]])
            }
            node_ages <- rep(NA, Nnode(tree))
            for(k in 1:Nnode(mol_tree_j)){
                tips <- mol_tree_j$tip.label[unlist(Descendants(mol_tree_j, node = k + Ntip(mol_tree_j), type = "tips"))]
                node_x <- mrca.phylo(tree, node = tips)
                node_ages[[node_x - Ntip(tree)]] <- (root_age - (node.depth.edgelength(mol_tree_j)[[k + Ntip(mol_tree_j)]])) * 100
            }
            timetree <- timePaleoPhy(tree, timeData = ages, type = "equal", node.mins = node_ages)
        }
        timetree_scaled <- timetree
        timetree_scaled$edge.length <- timetree_scaled$edge.length / max(node.depth.edgelength(timetree))
        for(z in 1:length(modelC)){
            if(i == "tree_no_fossils"){
              fossils <- rownames(ages)[which(ages[,1] >0)]
              tp <- lapply(tp, function(x) x[!(row.names(x) %in% fossils),])
            }
            tp <- lapply(tp, function(x) x[match(timetree$tip.label, rownames(x)),]) #rearrange tip prior order. 
            fit_SMM <- fitMk(timetree_scaled, x = tp[[5]], model = modelC[[z]], pi="fitzjohn")
            res_modelC[[x]] <- fit_SMM
            newrow <- data.frame(tree = i, mol_ages = j, k = max(modelC[[z]]), logL = fit_SMM$logLik)
            df_res_modelC <- rbind(df_res_modelC, newrow)
            x <- x + 1
        }
    }
    setwd(wd)
}

######################################## Model D run #################################################

setwd(wd)
res_modelD <- list()
df_res_modelD <- data.frame()

x <- 1
for(i in c("tree_fossils", "tree_fossils_Wang", "tree_no_fossils")){
    tp <- get_tip_priors(tip_states)
    colnames(tp[[6]]) <- colnames(modelD[[1]])
    setwd("trees/")
    mol_trees_names <- list.files()
    for(j in mol_trees_names){
        mol_tree_j <- read.nexus(j)
        root_age <- max(node.depth.edgelength(mol_tree_j))
        if(i == "tree_no_fossils"){
            timetree <- ladderize(mol_tree_j)
            timetree$edge.length <- timetree$edge.length * 100
        } else {
            if(i == "tree_fossils"){
                tree <- ladderize(trees[[1]])
            }
            if(i == "tree_fossils_Wang"){
                tree <- ladderize(trees[[2]])
            }
            node_ages <- rep(NA, Nnode(tree))
            for(k in 1:Nnode(mol_tree_j)){
                tips <- mol_tree_j$tip.label[unlist(Descendants(mol_tree_j, node = k + Ntip(mol_tree_j), type = "tips"))]
                node_x <- mrca.phylo(tree, node = tips)
                node_ages[[node_x - Ntip(tree)]] <- (root_age - (node.depth.edgelength(mol_tree_j)[[k + Ntip(mol_tree_j)]])) * 100
            }
            timetree <- timePaleoPhy(tree, timeData = ages, type = "equal", node.mins = node_ages)
        }
        timetree_scaled <- timetree
        timetree_scaled$edge.length <- timetree_scaled$edge.length / max(node.depth.edgelength(timetree))
        for(z in 1:length(modelD)){
            if(i == "tree_no_fossils"){
              fossils <- rownames(ages)[which(ages[,1] >0)]
              tp <- lapply(tp, function(x) x[!(row.names(x) %in% fossils),])
            }
            tp <- lapply(tp, function(x) x[match(timetree$tip.label, rownames(x)),]) #rearrange tip prior order. 
            fit_SMM <- fitMk(timetree_scaled, x = tp[[6]], model = modelD[[z]], pi="fitzjohn")
            res_modelD[[x]] <- fit_SMM
            newrow <- data.frame(tree = i, mol_ages = j, k = max(modelD[[z]]), logL = fit_SMM$logLik)
            df_res_modelD <- rbind(df_res_modelD, newrow)
            x <- x + 1
        }
    }
    setwd(wd)
}

################## model averaging and plotting ################## 

#ModelA

setwd(wd)
dir.create("ModelA")
dir.create("ModelB")
dir.create("ModelC")
dir.create("ModelD")

library(strap)
x <- 1
for(i in c("tree_fossils", "tree_fossils_Wang", "tree_no_fossils")){
  setwd(wd)
  tp <- get_tip_priors(tip_states)
  colnames(tp[[1]]) <- colnames(modelA[[1]])
  if(i == "tree_no_fossils"){
    timetree <- mol_con
    timetree$edge.length <- timetree$edge.length * 100
  } else {
    if(i == "tree_fossils"){
      tree <- ladderize(trees[[1]])
    }
    if(i == "tree_fossils_Wang"){
      tree <- ladderize(trees[[2]])
    }
    node_ages <- rep(NA, Nnode(tree))
    for(k in 1:Nnode(mol_con)){
      tips <- mol_con$tip.label[unlist(Descendants(mol_con, node = k + Ntip(mol_con), type = "tips"))]
      node_x <- mrca.phylo(tree, node = tips)
      node_ages[[node_x - Ntip(tree)]] <- (root_age - (node.depth.edgelength(mol_con)[[k + Ntip(mol_con)]])) * 100
    }
    timetree <- timePaleoPhy(tree, timeData = ages, type = "equal", node.mins = node_ages)
  }
  LogL <- df_res_modelA[x:(x+15),4]
  k <- df_res_modelA[x:(x+15),3]
  BIC_wt <- comp_models(loglik = LogL, k = k, n = Ntip(timetree), method = "BIC")
  
  if(i == "tree_no_fossils"){
    timetree$root.time <- max(node.depth.edgelength(timetree))
    fossils <- rownames(ages)[which(ages[,1] >0)]
    tp <- lapply(tp, function(x) x[!(row.names(x) %in% fossils),])
  }
  tp <- lapply(tp, function(a) {
    idx <- match(rownames(a), timetree$tip.label)
    # Filter out NA values and reorder
    a[!is.na(idx), ][order(idx[!is.na(idx)]), ]
  })
  mdlA_avg <- matrix(0, Nnode(timetree), 4)
  q <- 1
  for(j in x:(x+15)){
    ancMk <- ancr(res_modelA[[j]])
    mdlA_avg <- mdlA_avg + (ancMk$ace * BIC_wt[q, 6])
    q <- q+1
  }
  
  cols <- cols1
  names(cols) <- colnames(mdlA_avg)
  
  setwd("ModelA/")
  
  pdf(paste(i, "_modelA.pdf", sep = ""), width = 9, height = 11)
  par(mar = c(1.5,1.5,1.5,1.5))
  geoscalePhylo(timetree, units = "Period", tick.scale = "Period", boxes ="Period", cex.tip = 0.5, label.offset = 10, width = 1.8, quat.rm = T)
  par(fg = "transparent")
  nodelabels(node=1:Nnode(timetree)+Ntip(timetree),pie= mdlA_avg ,piecol=cols,cex=0.4)
  tiplabels(pie=tp[[1]],piecol=cols,cex=0.3)
  par(fg = "black")
  legend("topright", names(cols),  pch = 16, col = cols, pt.cex = 1.2,  cex = 0.8)
  title(paste("\n", i, "_modelA.pdf", sep = ""))
  dev.off()
    
  N <- Ntip(timetree)
  Bilateria <- mrca.phylo(timetree, node = c("Lottia_gigantea", "Homo_sapiens")) - N
  Homoscleromorpha <- mrca.phylo(timetree, node = c("Corticium_candelabrum", "Oscarella_pearsei")) - N
  Calcarea <- mrca.phylo(timetree, node = c("Sycon_coactum", "Leucetta_giribeti")) - N
  Hexactinellida <- mrca.phylo(timetree, node = c("Rossella_fibulata", "Vazella_pourtalesi")) - N
  Demospongiae <- mrca.phylo(timetree, node = c("Aplysina_aerophoba", "Iophon_unicorne")) - N
  Homosc_Calc <- mrca.phylo(timetree, node = c("Corticium_candelabrum", "Leucetta_giribeti")) - N
  Hex_Dem <- mrca.phylo(timetree, node = c("Rossella_fibulata", "Iophon_unicorne")) - N
  Porifera <- mrca.phylo(timetree, node = c("Corticium_candelabrum", "Iophon_unicorne")) - N
  
  rownames(mdlA_avg)[c(Bilateria, Homoscleromorpha, Calcarea, Hexactinellida, Demospongiae, Homosc_Calc, Hex_Dem, Porifera)] <- c("Bilateria", "Homoscleromorpha", "Calcarea", "Hexactinellida", "Demospongiae", "Homosc_Calc", "Hex_Dem", "Porifera")
  
  write.csv(round(mdlA_avg, 3), file = paste(i, "_mdlA_avg.csv", sep = ""))

  
  x <- x + 16
}


#ModelB
    
x <- 1
for(i in c("tree_fossils", "tree_fossils_Wang", "tree_no_fossils")){
  setwd(wd)
  tp <- get_tip_priors(tip_states)
  colnames(tp[[2]]) <- colnames(modelB[[1]])
  if(i == "tree_no_fossils"){
    timetree <- mol_con
    timetree$edge.length <- timetree$edge.length * 100
  } else {
    if(i == "tree_fossils"){
      tree <- ladderize(trees[[1]])
    }
    if(i == "tree_fossils_Wang"){
      tree <- ladderize(trees[[2]])
    }
    node_ages <- rep(NA, Nnode(tree))
    for(k in 1:Nnode(mol_con)){
      tips <- mol_con$tip.label[unlist(Descendants(mol_con, node = k + Ntip(mol_con), type = "tips"))]
      node_x <- mrca.phylo(tree, node = tips)
      node_ages[[node_x - Ntip(tree)]] <- (root_age - (node.depth.edgelength(mol_con)[[k + Ntip(mol_con)]])) * 100
    }
    timetree <- timePaleoPhy(tree, timeData = ages, type = "equal", node.mins = node_ages)
  }
  LogL <- df_res_modelB[x:(x+31),4]
  k <- df_res_modelB[x:(x+31),3]
  BIC_wt <- comp_models(loglik = LogL, k = k, n = Ntip(timetree), method = "BIC")
  
  if(i == "tree_no_fossils"){
    timetree$root.time <- max(node.depth.edgelength(timetree))
    fossils <- rownames(ages)[which(ages[,1] >0)]
    tp <- lapply(tp, function(x) x[!(row.names(x) %in% fossils),])
  }
  tp <- lapply(tp, function(a) {
    idx <- match(rownames(a), timetree$tip.label)
    # Filter out NA values and reorder
    a[!is.na(idx), ][order(idx[!is.na(idx)]), ]
  })
  mdlB_avg <- matrix(0, Nnode(timetree), 5)
  q <- 1
  for(j in x:(x+31)){
    ancMk <- ancr(res_modelB[[j]])
    mdlB_avg <- mdlB_avg + (ancMk$ace * BIC_wt[q, 6])
    q <- q+1
  }
  
   cols <- cols2 
   names(cols) <- colnames(mdlB_avg)
  
  setwd("ModelB/")
  
  pdf(paste(i, "_modelB.pdf", sep = ""), width = 9, height = 11)
  par(mar = c(1.5,1.5,1.5,1.5))
  geoscalePhylo(timetree, units = "Period", tick.scale = "Period", boxes ="Period", cex.tip = 0.5, label.offset = 10, width = 1.8, quat.rm = T)
  par(fg = "transparent")
  nodelabels(node=1:Nnode(timetree)+Ntip(timetree),pie= mdlB_avg ,piecol=cols,cex=0.4)
  tiplabels(pie=tp[[2]],piecol=cols,cex=0.3)
  par(fg = "black")
  legend("topright", names(cols),  pch = 16, col = cols, pt.cex = 1.2,  cex = 0.8)
  title(paste("\n", i, "_modelB.pdf", sep = ""))
  dev.off()
  
  N <- Ntip(timetree)
  Bilateria <- mrca.phylo(timetree, node = c("Lottia_gigantea", "Homo_sapiens")) - N
  Homoscleromorpha <- mrca.phylo(timetree, node = c("Corticium_candelabrum", "Oscarella_pearsei")) - N
  Calcarea <- mrca.phylo(timetree, node = c("Sycon_coactum", "Leucetta_giribeti")) - N
  Hexactinellida <- mrca.phylo(timetree, node = c("Rossella_fibulata", "Vazella_pourtalesi")) - N
  Demospongiae <- mrca.phylo(timetree, node = c("Aplysina_aerophoba", "Iophon_unicorne")) - N
  Homosc_Calc <- mrca.phylo(timetree, node = c("Corticium_candelabrum", "Leucetta_giribeti")) - N
  Hex_Dem <- mrca.phylo(timetree, node = c("Rossella_fibulata", "Iophon_unicorne")) - N
  Porifera <- mrca.phylo(timetree, node = c("Corticium_candelabrum", "Iophon_unicorne")) - N
  
  rownames(mdlB_avg)[c(Bilateria, Homoscleromorpha, Calcarea, Hexactinellida, Demospongiae, Homosc_Calc, Hex_Dem, Porifera)] <- c("Bilateria", "Homoscleromorpha", "Calcarea", "Hexactinellida", "Demospongiae", "Homosc_Calc", "Hex_Dem", "Porifera")
  
  write.csv(round(mdlB_avg, 3), file = paste(i, "_mdlB_avg.csv", sep = ""))
  
  x <- x + 32
}
setwd(wd)

#ModelC

x <- 1
for(i in c("tree_fossils", "tree_fossils_Wang", "tree_no_fossils")){
  setwd(wd)
  tp <- get_tip_priors(tip_states)
  colnames(tp[[5]]) <- colnames(modelC[[1]])
  if(i == "tree_no_fossils"){
    timetree <- mol_con
    timetree$edge.length <- timetree$edge.length * 100
  } else {
    if(i == "tree_fossils"){
      tree <- ladderize(trees[[1]])
    }
    if(i == "tree_fossils_Wang"){
      tree <- ladderize(trees[[2]])
    }
    node_ages <- rep(NA, Nnode(tree))
    for(k in 1:Nnode(mol_con)){
      tips <- mol_con$tip.label[unlist(Descendants(mol_con, node = k + Ntip(mol_con), type = "tips"))]
      node_x <- mrca.phylo(tree, node = tips)
      node_ages[[node_x - Ntip(tree)]] <- (root_age - (node.depth.edgelength(mol_con)[[k + Ntip(mol_con)]])) * 100
    }
    timetree <- timePaleoPhy(tree, timeData = ages, type = "equal", node.mins = node_ages)
  }
  LogL <- df_res_modelC[x:(x+31),4]
  k <- df_res_modelC[x:(x+31),3]
  BIC_wt <- comp_models(loglik = LogL, k = k, n = Ntip(timetree), method = "BIC")
  
  if(i == "tree_no_fossils"){
    timetree$root.time <- max(node.depth.edgelength(timetree))
    fossils <- rownames(ages)[which(ages[,1] >0)]
    tp <- lapply(tp, function(x) x[!(row.names(x) %in% fossils),])
  }
  tp <- lapply(tp, function(a) {
    idx <- match(rownames(a), timetree$tip.label)
    # Filter out NA values and reorder
    a[!is.na(idx), ][order(idx[!is.na(idx)]), ]
  })
  mdlC_avg <- matrix(0, Nnode(timetree), 8)
  q <- 1
  for(j in x:(x+31)){
    ancMk <- ancr(res_modelC[[j]])
    mdlC_avg <- mdlC_avg + (ancMk$ace * BIC_wt[q, 6])
    q <- q+1
  }
  
  cols <- cols3
  names(cols) <- colnames(mdlC_avg)
  
  setwd("ModelC/")
  
  pdf(paste(i, "_modelC.pdf", sep = ""), width = 9, height = 11)
  par(mar = c(1.5,1.5,1.5,1.5))
  geoscalePhylo(timetree, units = "Period", tick.scale = "Period", boxes ="Period", cex.tip = 0.5, label.offset = 10, width = 1.8, quat.rm = T)
  par(fg = "transparent")
  nodelabels(node=1:Nnode(timetree)+Ntip(timetree),pie= mdlC_avg ,piecol=cols,cex=0.4)
  tiplabels(pie=tp[[5]],piecol=cols,cex=0.3)
  par(fg = "black")
  legend("topright", names(cols),  pch = 16, col = cols, pt.cex = 1.2,  cex = 0.8)
  title(paste("\n", i, "_modelC.pdf", sep = ""))
  dev.off()
  
  N <- Ntip(timetree)
  Bilateria <- mrca.phylo(timetree, node = c("Lottia_gigantea", "Homo_sapiens")) - N
  Homoscleromorpha <- mrca.phylo(timetree, node = c("Corticium_candelabrum", "Oscarella_pearsei")) - N
  Calcarea <- mrca.phylo(timetree, node = c("Sycon_coactum", "Leucetta_giribeti")) - N
  Hexactinellida <- mrca.phylo(timetree, node = c("Rossella_fibulata", "Vazella_pourtalesi")) - N
  Demospongiae <- mrca.phylo(timetree, node = c("Aplysina_aerophoba", "Iophon_unicorne")) - N
  Homosc_Calc <- mrca.phylo(timetree, node = c("Corticium_candelabrum", "Leucetta_giribeti")) - N
  Hex_Dem <- mrca.phylo(timetree, node = c("Rossella_fibulata", "Iophon_unicorne")) - N
  Porifera <- mrca.phylo(timetree, node = c("Corticium_candelabrum", "Iophon_unicorne")) - N
  
  rownames(mdlC_avg)[c(Bilateria, Homoscleromorpha, Calcarea, Hexactinellida, Demospongiae, Homosc_Calc, Hex_Dem, Porifera)] <- c("Bilateria", "Homoscleromorpha", "Calcarea", "Hexactinellida", "Demospongiae", "Homosc_Calc", "Hex_Dem", "Porifera")
 
  write.csv(round(mdlC_avg, 3), file = paste(i, "_mdlC_avg.csv", sep = ""))
  
  x <- x + 32
}


#ModelD

x <- 1
for(i in c("tree_fossils", "tree_fossils_Wang", "tree_no_fossils")){
  setwd(wd)
  tp <- get_tip_priors(tip_states)
  colnames(tp[[6]]) <- colnames(modelD[[1]])
  if(i == "tree_no_fossils"){
    timetree <- mol_con
    timetree$edge.length <- timetree$edge.length * 100
  } else {
    if(i == "tree_fossils"){
      tree <- ladderize(trees[[1]])
    }
    if(i == "tree_fossils_Wang"){
      tree <- ladderize(trees[[2]])
    }
    node_ages <- rep(NA, Nnode(tree))
    for(k in 1:Nnode(mol_con)){
      tips <- mol_con$tip.label[unlist(Descendants(mol_con, node = k + Ntip(mol_con), type = "tips"))]
      node_x <- mrca.phylo(tree, node = tips)
      node_ages[[node_x - Ntip(tree)]] <- (root_age - (node.depth.edgelength(mol_con)[[k + Ntip(mol_con)]])) * 100
    }
    timetree <- timePaleoPhy(tree, timeData = ages, type = "equal", node.mins = node_ages)
  }
  LogL <- df_res_modelD[x:(x+63),4]
  k <- df_res_modelD[x:(x+63),3]
  BIC_wt <- comp_models(loglik = LogL, k = k, n = Ntip(timetree), method = "BIC")
  
  if(i == "tree_no_fossils"){
    timetree$root.time <- max(node.depth.edgelength(timetree))
    fossils <- rownames(ages)[which(ages[,1] >0)]
    tp <- lapply(tp, function(x) x[!(row.names(x) %in% fossils),])
  }
  tp <- lapply(tp, function(a) {
    idx <- match(rownames(a), timetree$tip.label)
    # Filter out NA values and reorder
    a[!is.na(idx), ][order(idx[!is.na(idx)]), ]
  })
  mdlD_avg <- matrix(0, Nnode(timetree), 9)
  q <- 1
  for(j in x:(x+63)){
    ancMk <- ancr(res_modelD[[j]])
    mdlD_avg <- mdlD_avg + (ancMk$ace * BIC_wt[q, 6])
    q <- q+1
  }
  
  cols <- cols4
  names(cols) <- colnames(mdlD_avg)
  
  setwd("ModelD/")
  
  pdf(paste(i, "_modelD.pdf", sep = ""), width = 9, height = 11)
  par(mar = c(1.5,1.5,1.5,1.5))
  geoscalePhylo(timetree, units = "Period", tick.scale = "Period", boxes ="Period", cex.tip = 0.5, label.offset = 10, width = 1.8, quat.rm = T)
  par(fg = "transparent")
  nodelabels(node=1:Nnode(timetree)+Ntip(timetree),pie= mdlD_avg ,piecol=cols,cex=0.4)
  tiplabels(pie=tp[[6]],piecol=cols,cex=0.3)
  par(fg = "black")
  legend("topright", names(cols),  pch = 16, col = cols, pt.cex = 1.2,  cex = 0.8)
  title(paste("\n", i, "_modelD.pdf", sep = ""))
  dev.off()
  
  N <- Ntip(timetree)
  Bilateria <- mrca.phylo(timetree, node = c("Lottia_gigantea", "Homo_sapiens")) - N
  Homoscleromorpha <- mrca.phylo(timetree, node = c("Corticium_candelabrum", "Oscarella_pearsei")) - N
  Calcarea <- mrca.phylo(timetree, node = c("Sycon_coactum", "Leucetta_giribeti")) - N
  Hexactinellida <- mrca.phylo(timetree, node = c("Rossella_fibulata", "Vazella_pourtalesi")) - N
  Demospongiae <- mrca.phylo(timetree, node = c("Aplysina_aerophoba", "Iophon_unicorne")) - N
  Homosc_Calc <- mrca.phylo(timetree, node = c("Corticium_candelabrum", "Leucetta_giribeti")) - N
  Hex_Dem <- mrca.phylo(timetree, node = c("Rossella_fibulata", "Iophon_unicorne")) - N
  Porifera <- mrca.phylo(timetree, node = c("Corticium_candelabrum", "Iophon_unicorne")) - N
  
  rownames(mdlD_avg)[c(Bilateria, Homoscleromorpha, Calcarea, Hexactinellida, Demospongiae, Homosc_Calc, Hex_Dem, Porifera)] <- c("Bilateria", "Homoscleromorpha", "Calcarea", "Hexactinellida", "Demospongiae", "Homosc_Calc", "Hex_Dem", "Porifera")
 
  write.csv(round(mdlD_avg, 3), file = paste(i, "_mdlD_avg.csv", sep = ""))
  
  x <- x + 64
}


###multi figure plot

cols <- c(spicules_absent = "#999999", mineralisation_absent = "#CC79A7", 
  silica = "#0072B2", calcite = "#009E73", bimineralic = "#D55E00"
)
state_names <- c("spicules_absent", "mineralisation_absent", "silica", "calcite", 
"bimineralic")

par(mfrow=c(3,4))
reduced_tips <- c("Homo_sapiens", "Oscarella_pearsei", "Corticium_candelabrum", "Sycon_ciliatum", "Vauxia", "Eiffelia" , "Helicolocellus", "Cyathophycus", "Minitaspongia", "Protospongia", "Diagoniella",  "Rossella_fibulata", "Ircinia_fasciculata", "Vaceletia_sp", "Haliclona_penicillata", "Geodia_atlantica", "Clathrina_coriacea", "Chondrilla_caribensis", "Iophon_unicorne", "Vazella_pourtalesi")

for(i in c("tree_fossils", "tree_fossils_Wang", "tree_no_fossils")){
	for(j in c("ModelA", "ModelB", "ModelC", "ModelD")){
		setwd(wd)
 	 	tp <- get_tip_priors(tip_states)
  		colnames(tp[[1]]) <- colnames(modelA[[1]])
  		colnames(tp[[2]]) <- colnames(modelB[[1]])
  		colnames(tp[[5]]) <- colnames(modelC[[1]])
  		colnames(tp[[6]]) <- colnames(modelD[[1]])
  		names(tp)[1:2] <- c("ModelA", "ModelB")
  		names(tp)[5:6] <- c("ModelC", "ModelD")
  		if(i == "tree_no_fossils"){
    		timetree <- mol_con
    		timetree$edge.length <- timetree$edge.length * 100
  		} else {
    	if(i == "tree_fossils"){
      		tree <- ladderize(trees[[1]])
    	}
    	if(i == "tree_fossils_Wang"){
      		tree <- ladderize(trees[[2]])
    	}
    	node_ages <- rep(NA, Nnode(tree))
    	for(k in 1:Nnode(mol_con)){
      		tips <- mol_con$tip.label[unlist(Descendants(mol_con, node = k + Ntip(mol_con), type = "tips"))]
      		node_x <- mrca.phylo(tree, node = tips)
      		node_ages[[node_x - Ntip(tree)]] <- (root_age - (node.depth.edgelength(mol_con)[[k + Ntip(mol_con)]])) * 100
    	}
    		timetree <- timePaleoPhy(tree, timeData = ages, type = "equal", node.mins = node_ages)
  		}
  		
  		if(i == "tree_no_fossils"){
    		timetree$root.time <- max(node.depth.edgelength(timetree))
    		fossils <- rownames(ages)[which(ages[,1] >0)]
    		tp <- lapply(tp, function(x) x[!(row.names(x) %in% fossils),])
  		}
  		tp <- lapply(tp, function(a) {
    		idx <- match(rownames(a), timetree$tip.label)
    		# Filter out NA values and reorder
    		a[!is.na(idx), ][order(idx[!is.na(idx)]), ]
  		})
		setwd(paste(j, "/", sep = ""))
		
		t <- keep.tip(timetree, tip = intersect(reduced_tips, timetree$tip.label))
		anc_states <- read.csv(paste(i, "_", gsub("Model", "mdl", j), "_avg.csv", sep = ""), row.names = 1)
		anc_states <- as.matrix(anc_states)
		
		tp_x <- tp[[j]][intersect(reduced_tips, timetree$tip.label),]
		
		#get anc states for matching nodes
		anc_x <- matrix(0, Nnode(t), ncol(tp_x))
		colnames(anc_x) <- colnames(anc_states)
		for(k in 1:Nnode(t)){
			des <- Descendants(t, node = k+Ntip(t))[[1]]
			anc_node <- mrca.phylo(timetree, node = t$tip.label[des])
			anc_x[k,] <- anc_states[anc_node - Ntip(timetree),]
		}
		
		if(j == "ModelA"){
			tp_x <- cbind(0, tp_x)
			anc_x <- cbind(0, anc_x)
			colnames(tp_x) <- state_names
			colnames(anc_x) <- state_names
		}
		if(j == "ModelB"){
			colnames(tp_x) <- state_names
			colnames(anc_x) <- state_names
		}
		if(j == "ModelC"){
			tp_x <- cbind(0, tp_x)
			anc_x <- cbind(0, anc_x)
			tp_x <- cbind(tp_x[,1:2], rowSums(tp_x[,3:5]), tp_x[,6], rowSums(tp_x[,7:9]))
			anc_x <- cbind(anc_x[,1:2], rowSums(anc_x[,3:5]), anc_x[,6], rowSums(anc_x[,7:9]))
			colnames(tp_x) <- state_names
			colnames(anc_x) <- state_names
			tp_x[which(tp_x > 1)] <- 1
		}
		if(j == "ModelD"){
			tp_x <- cbind(tp_x[,1:2], rowSums(tp_x[,3:5]), tp_x[,6], rowSums(tp_x[,7:9]))
			anc_x <- cbind(anc_x[,1:2], rowSums(anc_x[,3:5]), anc_x[,6], rowSums(anc_x[,7:9]))
			colnames(tp_x) <- state_names
			colnames(anc_x) <- state_names
			tp_x[which(tp_x > 1)] <- 1
		}
        tp_x <- tp_x[t$tip.label,]
        setwd(wd)
        plot(t, no.margin = T, label.offset = 50, cex = 0.1, tip.color = "white")
        cladelabels(t,node= mrca.phylo(t, node = c("Chondrilla_caribensis", "Geodia_atlantica")),"Demo", offset = 1)
        cladelabels(t,node= mrca.phylo(t, node = c("Rossella_fibulata", "Vazella_pourtalesi")),"Hex", offset = 1)
        cladelabels(t,node= mrca.phylo(t, node = c("Sycon_ciliatum", "Clathrina_coriacea")),"Calc", offset = 1)
        cladelabels(t,node= mrca.phylo(t, node = c("Oscarella_pearsei",     "Corticium_candelabrum")),"Hom", offset = 1)
        par(fg = "transparent")
        nodelabels(node=1:Nnode(t)+Ntip(t),pie= anc_x ,piecol=cols,cex=1)
        tiplabels(pie=tp_x,piecol=cols,cex=1)
        par(fg = "black")
    }
}

  
