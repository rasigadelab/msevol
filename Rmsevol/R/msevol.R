################################################################################
## MSE1 RESISTRACK MODEL
##  A bacterial hospital ecosystem simulation framework
## (c)2018-2024 Jean-Philippe Rasigade, <jean-philippe.rasigade@univ-lyon1.fr>
##  University of Lyon, Centre International de Recherche en Infectiologie
################################################################################
# BUILD CONFIGURATION AND RUN A SIMULATION
################################################################################

# UPDATE 2024-09-19 - Suppressing transposons
# UPDATE 2022-02-08 to deal with  R-32 integer limitations (NL)
# UPDATE 2022-01-xx with random seed (JPR)
# UPDATE 2021-12-01 with transposon (JRP)
# Match Universe.h Cfg structure for consistency


#'@title Initialize an empty model
#'@export
emptymodel <- function() {
  prm <- list()

  ## multedge <- data.table(source = integer(), target = integer(), multiplicity = integer())
  multedge <- data.table(source = integer(), target = integer(), multiplicity = numeric())

  # VERTEX TYPES

  prm[["Vertex_GeneArchetype"]] <- data.table(id = integer(), susceptibility.0 = numeric(), susceptibility.1 = numeric(), fitness = numeric())
  prm[["Vertex_Gene"]] <- data.table(id = integer())

  prm[["Vertex_PlasmidArchetype"]] <- data.table(id = integer(), loss = numeric(), transfer = numeric(), max_count = integer(), fitness = numeric())
  prm[["Vertex_Plasmid"]] <- data.table(id = integer())

  prm[["Vertex_ChromosomeArchetype"]] <- data.table(id = integer(), fitness = numeric(), survival = numeric())
  prm[["Vertex_Chromosome"]] <- data.table(id = integer())

  prm[["Vertex_CellArchetype"]] <- data.table(id = integer())
  prm[["Vertex_Cell"]] <- data.table(id = integer())

  prm[["Vertex_PatchArchetype"]] <- data.table(id = integer(), capacity = numeric(), pressure.0 = numeric(), pressure.1 = numeric())
  prm[["Vertex_Patch"]] <- data.table(id = integer())

  # EDGE TYPES - ARCHETYPES

  prm[["Edge_GeneArchetype_Gene_Multiplicity"]] <- data.table(multedge)
  prm[["Edge_ChromosomeArchetype_Chromosome_Multiplicity"]] <- data.table(multedge)
  prm[["Edge_PlasmidArchetype_Plasmid_Multiplicity"]] <- data.table(multedge)
  prm[["Edge_CellArchetype_Cell_Multiplicity"]] <- data.table(multedge)
  prm[["Edge_PatchArchetype_Patch_Multiplicity"]] <- data.table(multedge)

  # EDGE TYPES - CONTAINERS

  prm[["Edge_Gene_Plasmid_Multiplicity"]] <- data.table(multedge)
  prm[["Edge_Gene_Chromosome_Multiplicity"]] <- data.table(multedge)
  prm[["Edge_Plasmid_Cell_Multiplicity"]] <- data.table(multedge)
  prm[["Edge_Chromosome_Cell_Multiplicity"]] <- data.table(multedge)
  prm[["Edge_Cell_Patch_Multiplicity"]] <- data.table(multedge)

  # EDGE TYPES - SPECIFICS

  prm[["Edge_Patch_Patch_Diffusion"]] <- data.table(source = integer(), target = integer(), diffusion = numeric())
  prm[["Edge_PlasmidArchetype_ChromosomeArchetype_PlasmidRange"]] <- data.table(source = integer(), target = integer())

  return(prm)
}

#'@title Generic container
#'@export
container <- function(prm, name, id) {
  vertex <- sprintf("Vertex_%s", name)
  prm[[vertex]] <- rbind(
    prm[[vertex]],
    data.table(id = as.integer(id))
  )
  return(prm)
}

#'@title Generic assignation
#'@export
assign <- function(prm, from, source, to, target, mult = 1) {
  edge <- sprintf("Edge_%s_%s_Multiplicity", from, to)
  prm[[edge]] <- rbind(
    prm[[edge]],
    ## data.table(source = as.integer(source), target = as.integer(target), multiplicity = as.integer(mult))
    data.table(source = as.integer(source), target = as.integer(target), multiplicity = mult)
  )
  return(prm)
}

###############################################
# ARCHETYPES

#'@title  Archetype: gene
#'@export
geneArchetype <- function(prm, id, susc, fitness) {
  names(susc) <- paste("susceptibility", 1:length(susc) - 1, sep = ".")
  prm[["Vertex_GeneArchetype"]] <- rbind(
    prm[["Vertex_GeneArchetype"]],
    data.table(id = as.integer(id), t(susc), fitness = fitness)
  )
  return(prm)
}

#'@title  gene
#'@export
gene <- function(prm, id, type) {
  prm <- prm %>% container("Gene", id) %>% assign("GeneArchetype", type, "Gene", id)
  return(prm)
}

# archGene(emptymodel(), 0, c(0.1, 0.2), 0.5)$Vertex_GeneArchetype

#'@title  Archetype: plasmid
#'@export
plasmidArchetype <- function(prm, id, loss, transfer, max_count, fitness) {
  prm[["Vertex_PlasmidArchetype"]] <- rbind(
    prm[["Vertex_PlasmidArchetype"]],
    data.table(id = as.integer(id), loss = loss, transfer = transfer, max_count = as.integer(max_count), fitness = fitness)
  )
  return(prm)
}

#'@title  plasmid
#'@export
plasmid <- function(prm, id, type) {
  prm <- prm %>% container("Plasmid", id) %>% assign("PlasmidArchetype", type, "Plasmid", id)
  return(prm)
}

# archPlasmid(emptymodel(), 0, 0, 0, 1)$Vertex_PlasmidArchetype


#'@title Archetype: chromosome
#'@export
chromosomeArchetype <- function(prm, id, fitness, survival) {
  prm[["Vertex_ChromosomeArchetype"]] <- rbind(
    prm[["Vertex_ChromosomeArchetype"]],
    data.table(id = as.integer(id), fitness = fitness, survival = survival)
  )
  return(prm)
}

#'@title chromosome
#'@export
chromosome <- function(prm, id, type) {
  prm <- prm %>% container("Chromosome", id) %>% assign("ChromosomeArchetype", type, "Chromosome", id)
  return(prm)
}

# archChromosome(emptymodel(), 0, 0.2, 0.95)$Vertex_ChromosomeArchetype


#'@title Archetype: cell
#'@export
cellArchetype <- function(prm, id) {
  prm[["Vertex_CellArchetype"]] <- rbind(
    prm[["Vertex_CellArchetype"]], data.table(id = as.integer(id))
  )
  return(prm)
}

#'@title cell
#'@export
cell <- function(prm, id, type) {
  prm <- prm %>% container("Cell", id) %>% assign("CellArchetype", type, "Cell", id)
  return(prm)
}

#'@title Archetype: patch
#'@export
patchArchetype <- function(prm, id, capacity, pressure) {
  names(pressure) <- paste("pressure", 1:length(pressure) - 1, sep = ".")
  prm[["Vertex_PatchArchetype"]] <- rbind(
    prm[["Vertex_PatchArchetype"]],
    data.table(id = as.integer(id), capacity = capacity, t(pressure))
  )
  return(prm)
}

#'@title  patch
#'@export
patch <- function(prm, id, type) {
  prm <- prm %>% container("Patch", id) %>% assign("PatchArchetype", type, "Patch", id)
  return(prm)
}

###############################################
# SPECIFICS

#'@title  diffusion edge
#'@export
diffusion <- function(prm, source, target, diffusion, symmetric = T) {
  prm[["Edge_Patch_Patch_Diffusion"]] <- rbind(
    prm[["Edge_Patch_Patch_Diffusion"]],
    data.table(source = as.integer(source), target = as.integer(target), diffusion = diffusion)
  )
  if(symmetric) {
    prm[["Edge_Patch_Patch_Diffusion"]] <- rbind(
      prm[["Edge_Patch_Patch_Diffusion"]],
      data.table(source = as.integer(target), target = as.integer(source), diffusion = diffusion)
    )
  }
  return(prm)
}

#'@title  plasmidRange edge
#'@export
plasmidRange <- function(prm, source, target) {
  prm[["Edge_PlasmidArchetype_ChromosomeArchetype_PlasmidRange"]] <- rbind(
    prm[["Edge_PlasmidArchetype_ChromosomeArchetype_PlasmidRange"]],
    data.table(source = as.integer(source), target = as.integer(target))
  )
  return(prm)
}


###############################################
# INTEGRATED LAUNCHER & PARSER

#'@title  Run & parse
#'@export
mse_run <- function(prm, path, n_steps = 100, n_writes = 1, s1 = 1, s2 = 2) {
  cat("Writing configuration...")
  indir <- sprintf("bin/%s_in", path)
  outdir <- sprintf("bin/%s_out", path)
  # unlink("some_directory", recursive = TRUE)
  dir.create(indir, showWarnings = F)
  dir.create(outdir, showWarnings = F)

  # Clean directories if needed
  filenames <- dir(indir)
  filenames <- filenames[grepl("(.csv)$", filenames)]
  unlink(filenames)

  filenames <- dir(outdir)
  filenames <- filenames[grepl("(.csv)$", filenames)]
  unlink(filenames)

  for(f in names(prm)) {
    fname <- sprintf("%s/%s.csv", indir, f)

    #fwrite(prm[[f]], file = fname, sep = "\t")

    ## if needed, adjust 'scipen' to ensure multiplicity is written as an integer
    if('multiplicity' %in% colnames(prm[[f]]) & nrow(prm[[f]])>0){
      fwrite(prm[[f]], file = fname, sep = "\t",
             scipen = nchar(format(max(prm[[f]]$multiplicity), scientific = FALSE)))
    } else {
      fwrite(prm[[f]], file = fname, sep = "\t")
    }

  }

  cat(" done.\nLaunching simulator:\n")
  cmd <- sprintf("bin/msevolcli.exe RUN %i %i %s %s %i %i", n_steps, n_writes, indir, outdir, s1, s2)
  system(cmd, invisible = F)

  cat("Parsing results...")
  res <- mse_parse(outdir)
  cat(" done.\n")
  return(res)
}
