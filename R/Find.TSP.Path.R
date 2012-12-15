Find.TSP.Path <-
function(PositionList, mutation.matrix, insertion.type = "cheapest_insertion", fix.start.pos = "Y"){

  dist.matrix <- as.matrix(dist(PositionList[,c('XCoord','YCoord','ZCoord')],diag = T, method = "euclidean"))
  
  linear.path.distance <- sum(dist.matrix[col(dist.matrix)==row(dist.matrix)+1])
  
  dimnames(dist.matrix) <- list(PositionList$Can.Count,PositionList$Can.Count)
  mutation.positions <- as.vector(which(colSums(mutation.matrix)>0))
  positions.available <- PositionList$Can.Count
  #graph.vertices <- intersect(mutation.positions, positions.available)
  
  #start.positions <- which(rownames(dist.matrix)%in% graph.vertices[1:length(graph.vertices)])
  #end.positions <- which(rownames(dist.matrix) %in% graph.vertices[1:length(graph.vertices)])
  #dist.mut.matrix <- dist.matrix[start.positions, end.positions]
  
  start.pos <- 1
  end.pos <- dim(dist.matrix)[1]
  atsp<- ATSP(dist.matrix[-c(start.pos,end.pos),-c(start.pos,end.pos)])
  atsp<-insert_dummy(atsp,label = "TempCombinedAcid")
  combined <- which(labels(atsp)=="TempCombinedAcid")
  atsp[combined, ] <- c(dist.matrix[-c(start.pos, end.pos), start.pos], 0)
  atsp[, combined] <- c(dist.matrix[end.pos, -c(start.pos, end.pos)], 0)
  if(fix.start.pos == "Y"){
    tour <- solve_TSP(atsp, method = insertion.type, control = list(start=combined))
  }
  else{
    tour <- solve_TSP(atsp, method = insertion.type)
  }
  tour
  
  path_labels <- c(rownames(dist.matrix)[1], labels(cut_tour(tour,combined)), rownames(dist.matrix)[length(rownames(dist.matrix))])
  path_ids <- match(path_labels, rownames(dist.matrix))
  
  return(list(candidate.path = as.numeric(path_labels), candidate.path.distance  = attributes(tour)$tour_length,
              dist.matrix = dist.matrix, linear.path.distance = linear.path.distance))
  
}
