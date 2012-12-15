GraphClust <-
function(mutation.data, position.data, insertion.type = "cheapest_insertion",
                      alpha = 0.05, MultComp = "Bonferroni", fix.start.pos = "Y",
                      Include.Culled = "Y", Include.Full = "Y"){
  

  ####################################
  #Gets the missing positions mtatrix#
  ####################################
  
  #If positions (x_i, x_i+n_i) are missing, skipped.position.left.vector will keep track of the x_i, and skipped.position.right.vector will keep track of the x_i+n_i
  LHS <- numeric()
  RHS <- numeric()
  
  missing.positions<-setdiff(c(1:dim(mutation.data)[2]), position.data[,2])
  
  LHS.found <- FALSE
  for(i in 1:dim(mutation.data)[2]){
    if(is.element(i, missing.positions) && LHS.found != TRUE){ #This block finds the left hand side
      LHS <- append(LHS, i)
      LHS.found <- TRUE
    }else if(!(is.element(i, missing.positions)) && LHS.found == TRUE){ #This block finds the right hand side
      RHS<- append(RHS, i-1)
      LHS.found  <- FALSE
    }else if(i == dim(mutation.data)[2] && is.element(i, missing.positions)){ #This block handles the case if the final aa is not included
      RHS <- append(RHS, i)
      LHS.found <- FALSE
    }
    
  }
  
  missing.values.matrix <- cbind(LHS, RHS)

  
  ########################################
  #Finds the TSP path and makes the graph#
  ########################################
  path.results <- Find.TSP.Path(position.data, mutation.data, insertion.type = insertion.type, fix.start.pos = fix.start.pos)
  
  protein.graph <- graph.empty(n = length(path.results$candidate.path))
  V(protein.graph)$name <- position.data$Can.Count
  V(protein.graph)$label <- V(protein.graph)$name
  
  for(i in 1:(length(path.results$candidate.path)-1)){
    current.vertex <- which(V(protein.graph)$label == path.results$candidate.path[i])
    destination.vertex <- which(V(protein.graph)$label == path.results$candidate.path[i+1])
    protein.graph <- add.edges(protein.graph, c(current.vertex, destination.vertex))
  }
  
  ###############################################################
  #Gets the remapped mutation matrix and then finds the clusters#
  ###############################################################
  remapped.mutation.matrix <- matrix(data= 0, nrow = dim(mutation.data)[1], ncol = dim(path.results$dist.matrix)[1])
  for(i in 1 :length(path.results$candidate.path)){
    remapped.mutation.matrix[,i] <- mutation.data[,path.results$candidate.path[i]]    
  }
  
  message(cat("Calculating Remapped Clusters."), appendLF = FALSE)
  if(inherits(try(cluster.results <- nmc(remapped.mutation.matrix, alpha = alpha, multtest = MultComp)),"try-error")==TRUE){
    unmapped.cluster <- "ERROR! Error While running clustering on the full linear protein! Verify that your mutation data matrix is not all 0's!"
  }else{
    if(!is.null(cluster.results)){
      if(is.null(dim(cluster.results))){ #This is the case when there is only one result
        column.names <- names(cluster.results)
        cluster.results<- matrix(cluster.results,nrow =1) #puts it into matrix format
        colnames(cluster.results) <- column.names
      }
    
      unmapped.cluster <- unmap.clusters(cluster.results,
                        path.results$candidate.path,
                        mutation.data,
                        position.data)
    }else{
      unmapped.cluster <- NULL
    }
  }
  
  ####################################
  #This calculates the Culled results#
  ####################################
  if(Include.Culled == "Y") {
    message(cat("Calculating Culled Clusters."), appendLF = FALSE)
    mutation.data.culled <- mutation.data[,position.data[,2]]
    if(inherits(try(culled.results<- rbind(nmc(mutation.data.culled, alpha=alpha,multtest = MultComp))),"try-error")==TRUE){
      culled.results <- "ERROR! Error While running clustering on the full linear protein! Verify that your mutation data matrix is not all 0's!"
    }else{
      path <- as.numeric(mapply(function(x){substring(x,2)}, names(mutation.data.culled)))
      if(!is.null(culled.results)){
        if(is.null(dim(culled.results))){
          column.names <- names(culled.results)
          culled.results <- matrix(culled.results, nrow = 1)
          colnames(culled.results) <- column.names
        }
        culled.results <- unmap.clusters(culled.results,path, mutation.data,position.data)
        
      }else{
        culled.results <- NULL
        
      }
    }
   }else{
    culled.results <- "None Requested!"
   }
  
  ##################################
  #This calculates the Full results#
  ##################################
  if(Include.Full == "Y") {
    message(cat("Calculating Full Clusters."), appendLF = FALSE)
    if(inherits(try(full.results <- rbind(nmc(mutation.data, alpha=alpha,multtest = MultComp))),"try-error")==TRUE){
      full.results <- "ERROR! Error While running clustering on the culled linear protein! Verify that your mutation data matrix is not all 0's in the remaining positions!"
    }else{
      if(!is.null(full.results)){
        if(is.null(dim(full.results))){
          column.names <- names(full.results)
          full.results <- matrix(full.results, nrow = 1)
          colnames(full.results)<- column.names
        }
      }else{
        full.results <-NULL
      }
    }
  }else{
    full.results <- "None Requested!"
  }
  
  return(list(Remapped = unmapped.cluster, OriginalCulled = culled.results, Original = full.results,
              candidate.path = path.results$candidate.path,
              path.distance = path.results$candidate.path.distance,
              linear.path.distance = path.results$linear.path.distance,
              protein.graph = protein.graph,
              missing.positions = missing.values.matrix))
}
