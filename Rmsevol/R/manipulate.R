################################################################################
## MSE1 RESISTRACK MODEL
##  A bacterial hospital ecosystem simulation framework
## (c)2018-2024 Jean-Philippe Rasigade, <jean-philippe.rasigade@univ-lyon1.fr>
##  University of Lyon, Centre International de Recherche en Infectiologie
################################################################################
# Author: JP Rasigade, N Lenuzza
################################################################################
# Additional Helper Functions to Visualize a Model, Run a Simulation, and
# Analyze Result
################################################################################



#'@title Parse results in some directory
#'@export
mse_parse <- function(datapath) {
  # datapath <- "bin/data_out2"
  filenames <- dir(datapath)
  filenames <- filenames[grepl("(.csv)$", filenames)]
  stopifnot(length(filenames) > 0)
  objectnames <- gsub("(.csv)$", "", filenames)

  is.vertex <- grepl("^(Vertex)", objectnames)
  is.edge <- grepl("^(Edge)", objectnames)
  vtxnames <- objectnames[is.vertex]
  edgnames <- objectnames[is.edge]

  vtx <- list()
  edg <- list()

  for(vtxname in vtxnames) {
    vtx[[ vtxname ]] <- fread( sprintf("%s/%s.csv", datapath, vtxname) , sep = "\t")
    # Replace IDs with object names
    setnames(vtx[[ vtxname ]], "id", vtxname)
  }

  for(edgname in edgnames) {
    edg[[ edgname ]] <- fread( sprintf("%s/%s.csv", datapath, edgname) , sep = "\t")
    # Replace IDs with object names
    splitnames <- strsplit(edgname, "_")[[1]]
    setnames( edg[[ edgname ]], c("source", "target"), splitnames[2:3] )
  }

  return(list(
    vtx = vtx, edg = edg
  ))
}

#'@title Merge multiplicity across two compatible edges
#'@export
mergemult <- function(l, r) {
  # Column names
  left  <- names(l)[1]
  mid   <- names(l)[2]
  right <- names(r)[2]
  stopifnot(mid == names(r)[1])

  tmp <- merge(l, r, by = c(mid, "step"), suffixes = c(".x", ".y"), allow.cartesian = T)
  tmp <- tmp[, multiplicity := multiplicity.x * multiplicity.y][ , .(multiplicity = sum(multiplicity)), by = c(left, right, "step")]
  tmp <- tmp[order(tmp[["step"]], tmp[[left]], tmp[[right]])]
  return( tmp  )
}

#'@title  Start counter from an existing edge
count_from <- function(left, right, mgraph) {
  edge <- sprintf("Edge_%s_%s_Multiplicity", left, right)
  return(mgraph$edg[[edge]])
}

#'@title Append multiplicity by name
#'@export
count_in <- function(l, right, mgraph) {
  # Column names
  left  <- names(l)[1]
  mid   <- names(l)[2]
  right_edge <- sprintf("Edge_%s_%s_Multiplicity", mid, right)

  return(
    mergemult(l, mgraph$edg[[right_edge]])
  )
}

##----------------------------------------------------------------------------##
## Plotting an input model as an interactive graph
##----------------------------------------------------------------------------##
## Note: color_map must be provided as hexadecimal (for automatic addition
## of transparency of 'archetype' nodes)


#'@title mse1.plot_graph
#'@export
mse1.plot_graph<- function(mseprm,
                           vertices_type = c("Patch", "Cell", "Chromosome",
                                             "Plasmid", "Gene"),
                           add_archetype = TRUE,
                           edges_type = c("inclusion", "plasmid_range",
                                          "diffusion"),
                           color_map = c("Patch" = "#FF7043",
                                         "Cell" = "#FFEE58",
                                         "Chromosome" = "#66BB6A",
                                         "Plasmid" = "#42A5F5",
                                         "Gene" = "#EF5350"),
                           layout = c("hr", "fr", "free")[1]){


  ## Check arguments
  ##---------------------------------
  if(!all(vertices_type %in% c("Patch", "Cell", "Chromosome", "Plasmid", "Gene"))){
    stop("Unknown vetices_type. Must be 'Patch', 'Cell', 'Chromosome', 'Plasmid' and/or 'Gene'")
  }


  if("plasmid_range" %in% edges_type){
    if(!add_archetype){
      edges_type <- edges_type[which(edges_type!= "plasmid_range")]
      message("Please include archetype to display Plasmid Range Edges")
    }

    if(!all( c("Plasmid", "Chromosome") %in% vertices_type)){
      edges_type <- edges_type[which(edges_type!= "plasmid_range")]
      message("Please include both 'Plasmid' and 'Chromosome' to display Plasmid Range Edges")
    }
  }


  if("diffusion" %in% edges_type){
    if(!("Patch"%in% vertices_type)){
      edges_type <- edges_type[which(edges_type!= "diffusion")]
      message("Please include 'Patch' to diffusion Edges.")
    }
  }




  ## Vertices
  ##--------------------------------

  nodeDt <-  data.table()

  ## Extraction of selected vertices

  for(type in vertices_type){


    if(add_archetype){

      ## archetypes (with labels for plot, aggregation of parameters)
      type_archetypeDt <- copy(mseprm[[paste("Vertex_", type, "Archetype", sep = "")]])
      colnames(type_archetypeDt)[1] <- "typed_id"
      type_archetypeDt[, label:= do.call(paste,
                                         c(Map(function(col, x) paste0(col, "=", x),
                                               colnames(.SD), .SD), sep = "\n"))]
      # .SDcols = -1]

      type_archetypeDt[, type:= type]
      type_archetypeDt[, category:= "archetype"]
      type_archetypeDt[, fullname := paste0(type, "Archetype_", typed_id, sep = "")]

      type_archetypeDt <- type_archetypeDt[, c("typed_id", "type", "category", "label", "fullname")]

      nodeDt <- rbind(nodeDt, type_archetypeDt)
    }


    ## agents (with labels for plot, = typed_id)
    type_vertexDt <- copy(mseprm[[paste("Vertex_", type, sep = "")]])
    colnames(type_vertexDt)[1] <- "typed_id"
    type_vertexDt[, type:= type]
    type_vertexDt[, category:= "agent"]
    type_vertexDt[, label := as.character(typed_id)]
    type_vertexDt[, fullname := paste0(type, "_", typed_id, sep = "")]
    nodeDt <- rbind(nodeDt, type_vertexDt)

  }


  ## Format : id & label
  nodeDt[is.na(label), label:= ""]
  nodeDt$id <- 0:(nrow(nodeDt)-1)

  ## Format: adding "level"
  level_map <- c("Patch" = 0,
                 "Cell" = 1,
                 "Chromosome" = 2,
                 "Plasmid" = 2,
                 "Gene" = 3)
  nodeDt[, level := level_map[type]]
  nodeDt[category == "archetype", level:= level+0.5]
  nodeDt[category == "archetype" & type == "Patch", level := -0.5]


  ## Format: adding color
  nodeDt[category=="agent", color := color_map[type]]
  nodeDt[category=="archetype", color := paste0(color_map[type], "80", sep = "")]


  ## Format: adding shape
  shape_map <- c("archetype" = "box",
                 "agent" = "circle")

  nodeDt[, shape := shape_map[category]]





  ## Edges
  ##--------------------------------

  edgelist <- data.frame(matrix(nrow = 0, ncol = 4))
  colnames(edgelist) <- c("containee", "container", "multiplicity", "edgetype")


  ##  extracting inclusion edges

  selected_type <- c()

  if("inclusion" %in% edges_type){

    all_types <- c("Cell_Patch", "Chromosome_Cell","Plasmid_Cell",
                   "Gene_Chromosome", "Gene_Plasmid")

    selected_type <- all_types[sapply(all_types,
                                      function(str){
                                        return(all(strsplit(str, split = "_" )[[1]] %in%
                                                     vertices_type))
                                      })]
  }

  if(add_archetype){
    selected_type <- c(selected_type,
                       paste(vertices_type, "Archetype_",
                             vertices_type, sep = "" ))
  }


  if(length(selected_type)){

    for(type in selected_type){

      containee <- strsplit(type, split = "_" )[[1]][1]
      container <- strsplit(type, split = "_" )[[1]][2]

      if(nrow(mseprm[[paste("Edge_", type, "_Multiplicity", sep = "")]])){
        edgelist<- rbind.data.frame(edgelist,
                                    cbind.data.frame(containee = paste(containee,
                                                                       mseprm[[paste("Edge_", type, "_Multiplicity", sep = "")]]$source,
                                                                       sep = "_"),
                                                     container = paste(container,
                                                                       mseprm[[paste("Edge_", type, "_Multiplicity", sep = "")]]$target,
                                                                       sep = "_"),
                                                     multiplicity = mseprm[[paste("Edge_", type, "_Multiplicity", sep = "")]]$multiplicity,
                                                     edgetype = "mult"))
      }
    }

  }







  ##  extracting PlasmidRange

  if("plasmid_range" %in% edges_type){
    if(nrow(mseprm[["Edge_PlasmidArchetype_ChromosomeArchetype_PlasmidRange"]])){
      edgelist<- rbind.data.frame(edgelist,
                                  cbind.data.frame(
                                    containee = paste("PlasmidArchetype",
                                                      mseprm[["Edge_PlasmidArchetype_ChromosomeArchetype_PlasmidRange"]]$source,
                                                      sep = "_"),
                                    container = paste("ChromosomeArchetype",
                                                      mseprm[["Edge_PlasmidArchetype_ChromosomeArchetype_PlasmidRange"]]$target,
                                                      sep = "_"),
                                    multiplicity = 1,
                                    edgetype = "plrange"))
    }
  }

  ## If needed, Diffusion
  if("diffusion" %in% edges_type){
    if(nrow(mseprm[["Edge_Patch_Patch_Diffusion"]])){
      edgelist<- rbind.data.frame(edgelist,
                                  cbind.data.frame(
                                    containee = paste("Patch",
                                                      mseprm[["Edge_Patch_Patch_Diffusion"]]$target,
                                                      sep = "_"),
                                    container = paste("Patch",
                                                      mseprm[["Edge_Patch_Patch_Diffusion"]]$source,
                                                      sep = "_"),
                                    multiplicity = mseprm[["Edge_Patch_Patch_Diffusion"]]$diffusion,
                                    edgetype = "diffusion"))
    }
  }




  ## Formatting

  edgelist$from <- sapply(edgelist$container,
                          function(i){nodeDt[fullname==i, ]$id})

  edgelist$to <- sapply(edgelist$containee,
                        function(i){nodeDt[fullname==i, ]$id})

  edgelist$label <- as.character(edgelist$multiplicity)


  edgelist$category <- ifelse(grepl("Archetype", edgelist$containee),
                              "archetype", "inclusion")
  edgelist$category[which(edgelist$edgetype == "plrange")] <- "plasmidrange"

  edgelist$label[which(edgelist$category == "archetype")] <- ""
  edgelist$label[which(edgelist$category == "plasmidrange")] <- ""


  edgelist$arch_type <- ifelse(edgelist$category == "archetype",
                               sub("Archetype.*", "", edgelist$containee), NA)

  edgelist$color <- ifelse(edgelist$category == "archetype",
                           color_map[edgelist$arch_type ],
                           "grey")

  edgelist$dashes <- ifelse(edgelist$category %in% c("inclusion", "diffusion"), FALSE, TRUE)


  edgelist$arrows  <- rep('to', nrow(edgelist))

  ## Add smoothness to diffusion edges to avoid overlap of bidirectional edges

  edgelist$smooth.enable <- rep(FALSE, nrow(edgelist))
  edgelist$smooth.type <- rep("curvedCW", nrow(edgelist))
  edgelist$smooth.roundness<- rep(0, nrow(edgelist))

  edgelist$smooth.enable[which( edgelist$edgetype == "diffusion")] <- TRUE
  edgelist$smooth.roundness[which( edgelist$edgetype == "diffusion")] <- 0.2





  ## Plotting graph
  ##-------------------------------------------------

  graph <- visNetwork(node = as.data.frame(nodeDt),
                      edges = edgelist)

  if (layout == "fr") {
    graph <- graph %>%
      visIgraphLayout(layout = "layout_with_fr",
                      smooth = TRUE)
  }

  if(layout == "hr") {
    graph <- graph %>%
      visHierarchicalLayout(direction = "UD", sortMethod = "directed")
  }

  print(graph)

  return(invisible(graph))

}





##----------------------------------------------------------------------------##
## Quick check for common errors frequently encountered during the model
## building process (uncomplete)
##----------------------------------------------------------------------------##

#'@title mse1.diagnose_model
#'@export
mse1.diagnose_model <- function(mseprm){

  mse1.type_agents <- c("Gene", "Plasmid", "Chromosome", "Cell", "Patch")
  mse1.type_edges <- c(
    paste(mse1.type_agents, "Archetype_", mse1.type_agents,
          "_Multiplicity", sep = ""),
    c("Gene_Plasmid_Multiplicity",
      "Gene_Chromosome_Multiplicity",
      "Plasmid_Cell_Multiplicity",
      "Chromosome_Cell_Multiplicity",
      "Cell_Patch_Multiplicity"),
    "PlasmidArchetype_ChromosomeArchetype_PlasmidRange",
    "Patch_Patch_Diffusion")


  ## Check that there is no duplicated "typed_id"

  cat("\n Check for duplicated ID : ")
  dupL <- FALSE

  for(typ in c(mse1.type_agents, paste(mse1.type_agents, "Archetype", sep = ""))){


    typed_id <-  mseprm[[paste0("Vertex_", typ, sep = "")]]$id
    repeated_id<- NULL

    if(length(typed_id)){
      count_table <- table(typed_id)
      repeated_id<- as.integer(names(count_table[count_table > 1]))
    }

    if(length(repeated_id)){
      dupL <- TRUE
      cat("\n   * ", typ, " ")
      message("[ duplication:", paste0(repeated_id, collapse =", "), "]", "\n")
    }
  }

  if(dupL){ message("!! Duplicated typed_id !!")}else{cat(" OK.")}



  ## Check if there exists duplicated assignation
  cat("\n Check for duplicated edges  :")
  edgeDupL <- FALSE
  for(typ in mse1.type_edges){
    typed_id <-  mseprm[[paste0("Edge_", typ, sep = "")]]

    if(nrow(typed_id)){
      n_typ <- unique(typed_id[, c("source", "target")])
      if(nrow(n_typ) != nrow(typed_id)){
        cat("\n   * ", typ, " ")
        message("[Duplicated edges !!!]")
        edgeDupL <- TRUE
      }
    }
  }

  if(edgeDupL){ message("!! Duplicated edges !!")}else{cat(" OK.")}



  ## Check if all "assignation" occurs in existing nodes

  cat("\n Check for unknown source(s)/target(s)  :")
  edgeValidityL <- FALSE

  for(typ in mse1.type_edges){
    source_typ <- strsplit(typ, split = "_" )[[1]][1]
    target_typ <- strsplit(typ, split = "_" )[[1]][2]

    typed_id <-  mseprm[[paste0("Edge_", typ, sep = "")]]

    if(nrow(typed_id)){

      missing_sources <- typed_id$source[!(typed_id$source %in% mseprm[[paste0("Vertex_", source_typ, sep = "")]]$id)]
      missing_targets <- typed_id$target[!(typed_id$target %in% mseprm[[paste0("Vertex_", target_typ, sep = "")]]$id)]

      if(length(missing_sources)){
        cat("\n   * ", typ, " ")
        message("[missing source(s):",
                paste0(as.character(missing_sources), collapse =", "),
                "]")
        edgeValidityL <- TRUE

      }

      if(length(missing_targets)){
        cat("\n   * ", typ, " ")
        message("[missing target(s):",
                paste0(as.character(missing_targets), collapse =", "),
                "]\n")
        edgeValidityL <- TRUE

      }
    }
  }
  if(edgeValidityL){ message("!! Edges with undeclared vertices !!")}else{cat(" OK.")}



  ## Check that plasmids with non-zero transfer rate has a defined PlasmidRange
  cat("\n Check for plasmidRange  : ")
  plasrangeL <- FALSE

  conjugative_plasmid <- mseprm[["Vertex_PlasmidArchetype"]][transfer>0, ]$id

  if(length(conjugative_plasmid)){
    missing_plasmidrange <-
      conjugative_plasmid[which(
        !(conjugative_plasmid %in%
            mseprm[["Edge_PlasmidArchetype_ChromosomeArchetype_PlasmidRange"]]$source)
      )]

    if(length(missing_plasmidrange)){
      cat("\n   * Conjugative plasmid : ")
      message("[missing plasmid range :",
              paste0(as.character(missing_plasmidrange), collapse =", "),
              "]\n")

      plasrangeL  <- TRUE
    }
  }
  if(plasrangeL){ message("!! Possibly missing plasmid ranges declaration !!")
  }else{cat(" OK.")}


}




##----------------------------------------------------------------------------##
## Extract one initial configuration from a simulation output
##----------------------------------------------------------------------------##
## For now, the inputs and outputs of the mse1 simulator have different formats.
## To allow restarting a simulation from any saved time step of an existing
## simulation, the function "mse1.extract_prm_from_mseres" extracts the model's
## state at that time and converts it into the input format
##----------------------------------------------------------------------------##

#'@title mse1.extract_prm_from_mseres
#'@export
mse1.extract_prm_from_mseres <- function(mseres, stepidx) {

  vtxnames <- names(mseres$vtx)
  edgnames <- names(mseres$edg)

  prmout <- emptymodel()

  for(vtxname in vtxnames) {
    prmout[[ vtxname ]] <- copy(mseres$vtx[[ vtxname ]])
    setnames(prmout[[ vtxname ]],  vtxname, "id")
    prmout[[ vtxname ]] <- prmout[[ vtxname ]][step == stepidx]
    prmout[[ vtxname ]][ , step:= NULL]
  }

  for(edgname in edgnames) {
    prmout[[ edgname ]] <- copy(mseres$edg[[ edgname ]])
    setnames( prmout[[ edgname ]],  1:2, c("source", "target"))
    prmout[[ edgname ]] <- prmout[[ edgname ]][step == stepidx]
    prmout[[ edgname ]][ , step:= NULL]
  }


  return(prmout)
}



##----------------------------------------------------------------------------##
## Add zeros to every time step
##----------------------------------------------------------------------------##
## Except for a few exceptions, mse1 only stores elements with a strictly
## positive multiplicity. This function fills in the data table with zero
## entries at time steps where they are not detected.
## The input data table must contain at least three columns, with two of them
## named 'multiplicity' and 'step' respectively. The function ensures that there
## exist one row (i.e. one multiplicity) per step per unique combination of
## columns except step and multiplicity.
##----------------------------------------------------------------------------##

#'@title mse1.add_zero_steps
#'@export
mse1.add_zero_steps <- function(dt){


  dt_wide <- dcast(
    dt,
    ... ~ step,
    value.var = "multiplicity",
    fill = 0
  )

  dt_with_zeros <- melt(
    dt_wide,
    id.vars = setdiff(names(dt), c("step", "multiplicity")),
    variable.name = "step",
    value.name = "multiplicity"
  )

  dt_with_zeros$step <- as.numeric(dt_with_zeros$step)

  return(dt_with_zeros)

}




##----------------------------------------------------------------------------##
## Describe cells' contents
##----------------------------------------------------------------------------##
## In simulations involving plasmids, certain cells may completely vanish from
## the ecosystem to reemerge later through horizontal transfers. Since the
## simulator does not retain memory of configurations that disappear, a
## reemerging configuration will be treated as new and assigned a new index.
## This can result in disruptions to temporal tracking.
##
## The function "mse1.describe_cells" describes each cell in terms of its
## archetype, plasmids, and chromosomes (= direct containees, sufficient to
## ensure identity, as plasmids and chromosomes are currently fixed in mse1)
## and provides the time window during which the configuration is detected.
## Each identical configuration is assigned an index, corresponding to the
## first one encountered in time (named CellPrototype)."
##----------------------------------------------------------------------------##

#'@title mse1.describe_cells
#'@export
mse1.describe_cells <-  function(mseres){

  ## Extracting unique indexes of cells throughout the simulation
  cell_description_dt <-
    mseres$vtx$Vertex_Cell[, .(first = min(step),
                               last = max(step)),
                           by = list(Vertex_Cell)]
  colnames(cell_description_dt)[1] <- "Cell"



  ## Adding their archetype
  cell_description_dt <-
    merge(cell_description_dt,
          unique(mseres$edg$Edge_CellArchetype_Cell_Multiplicity[, c("CellArchetype", "Cell")]),
          by = "Cell", all.x = TRUE)



  ## Assign content
  for(type in c("Chromosome", "Plasmid")){

    content_dt <- mseres$edg[[paste("Edge_", type, "_Cell_Multiplicity", sep = "")]][, 1:3]

    if(nrow(content_dt)){

      content_dt <- unique(content_dt)
      colnames(content_dt)[1] <- "type"

      content_dt <- dcast(content_dt, Cell ~ type, value.var = "multiplicity", fill = 0)
      colnames(content_dt)[2:ncol(content_dt)] <- paste(type, colnames(content_dt)[2:ncol(content_dt)],sep = "_")
    }


    cell_description_dt <- merge(cell_description_dt, content_dt, by = "Cell", all = TRUE)
  }
  cell_description_dt[is.na(cell_description_dt)] <- 0


  ## Assign a unique prototype
  var_cols <- setdiff(names(cell_description_dt), c("Cell", "first", "last", "CellArchetype"))
  cell_description_dt[, CellPrototype := min(Cell), by = c("CellArchetype", var_cols)]



  return(cell_description_dt)
}




##----------------------------------------------------------------------------##
## Clean up
##----------------------------------------------------------------------------##
## During a simulation, mse1 creates two folders (bin/xxx_in and bin/xxx_out,
## where 'xxx' is a user-provided name) and stores several .csv files. To avoid
## unnecessary storage of these files, it can be helpful to automatically delete
## them once the run is finished.
##----------------------------------------------------------------------------##


#'@title mse1.cleanup_csv
#'@export
mse1.cleanup_csv <-  function(csvres_folder){
  unlink(paste("bin/",csvres_folder, "_in", sep = ""), recursive = TRUE)
  unlink(paste("bin/",csvres_folder, "_out", sep = ""), recursive = TRUE)
}





##----------------------------------------------------------------------------##
## Functions for fluctuating models (i.e. models whose simulations is regularly
## interrupting and restarted in R)
##----------------------------------------------------------------------------##

#'@title mse1.shift_step
#'@export
mse1.shift_step <- function(mseres, shift_in_step = 0){

  for(tab_name in names(mseres$vtx)){
    mseres$vtx[[tab_name]][, step:= step + shift_in_step]
  }

  for(tab_name in names(mseres$edg)){
    mseres$edg[[tab_name]][, step:= step + shift_in_step]
  }

  return(mseres)

}


#'@title mse1.suppress_step
#'@export
mse1.suppress_step <- function(mseres, step_to_suppress){

  out <- list(vtx = list(), edg = list())

  for(tab_name in names(mseres$vtx)){
    out$vtx[[tab_name]]<- mseres$vtx[[tab_name]][step!= step_to_suppress, ]
  }

  for(tab_name in names(mseres$edg)){
    out$edg[[tab_name]]<- mseres$edg[[tab_name]][step!= step_to_suppress, ]
  }

  return(out)

}



#'@title mse1.merge_res
#'@export
mse1.merge_res <- function(first, second){

  out <- list(vtx = list(), edg =list())

  for(tab_name in names(first$vtx)){
    out$vtx[[tab_name]]<- rbind(first$vtx[[tab_name]],
                                second$vtx[[tab_name]])
  }

  for(tab_name in names(first$edg)){
    out$edg[[tab_name]]<- rbind(first$edg[[tab_name]],
                                second$edg[[tab_name]])
  }

  return(out)

}



