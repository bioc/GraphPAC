test_unmap.clusters  <- function() {  
  Residue <- c("MET","THR","GLU","TYR")
  Can.Count<- c(1,2,3,4)
  SideChain <- c("A","A","A","A")
  XCoord <- c(62.935,63.155,65.289,64.899)
  YCoord <- c(97.579,95.525,96.895,96.220)
  ZCoord <- c(30.223,27.079,24.308,20.615)
  
  path.data <- c(1,2,3,4)
  mutation.data <-as.data.frame(matrix(data = c(1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0),ncol = 4))
  names(mutation.data)<- c("V1","V2","V3","V4")
  position.data <- data.frame(Residue = Residue, Can.Count = Can.Count, SideChain = SideChain, XCoord = XCoord, YCoord = YCoord, ZCoord = ZCoord)
  
    
  cluster.results <- as.data.frame(matrix(data = c(3,1,1,1,1,3,3,1,3,24,12,12,0.002006723,0.007267573,0.007267573),nrow=3))
  names(cluster.results)<- c("cluster_size", "start", "end", "number", "p_value")
  
  checkEquals(dim(GraphPAC:::unmap.clusters(cluster.results, path.data, mutation.data, position.data))[1], 3)
  checkEquals(dim(GraphPAC:::unmap.clusters(cluster.results, path.data, mutation.data, position.data))[2], 5)
  checkTrue( sum(GraphPAC:::unmap.clusters(cluster.results, path.data, mutation.data, position.data)$start == c(1,1,3)) == 3)
  checkTrue( sum(GraphPAC:::unmap.clusters(cluster.results, path.data, mutation.data, position.data)$end == c(3,1,3)) == 3)
  checkTrue( sum(GraphPAC:::unmap.clusters(cluster.results, path.data, mutation.data, position.data)$start == c(1,1,1)) != 3)
  checkTrue( sum(GraphPAC:::unmap.clusters(cluster.results, path.data, mutation.data, position.data)$end == c(1,1,1)) != 3)           
}