test_GraphClust  <- function() {
  Residue <- c("MET","THR","GLU","TYR")
  Can.Count<- c(1,2,3,4)
  SideChain <- c("A","A","A","A")
  XCoord <- c(62.935,63.155,65.289,64.899)
  YCoord <- c(97.579,95.525,96.895,96.220)
  ZCoord <- c(30.223,27.079,24.308,20.615)
  
  mutation.matrix <-as.data.frame(matrix(data = c(1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0),ncol = 4))
  colnames(mutation.matrix)<- c("V1","V2","V3","V4")
  PositionList <- data.frame(Residue = Residue, Can.Count = Can.Count, SideChain = SideChain, XCoord = XCoord, YCoord = YCoord, ZCoord = ZCoord)
  
  checkEquals(length(GraphClust(mutation.matrix, PositionList)), 8)
  checkEquals(length(GraphClust(mutation.matrix, PositionList)$candidate.path), dim(mutation.matrix)[2])
  checkEquals(dim(GraphClust(mutation.matrix, PositionList)$Remapped)[1], 3)
  checkEquals(dim(GraphClust(mutation.matrix, PositionList)$OriginalCulled)[1], 3)
  checkEquals(dim(GraphClust(mutation.matrix, PositionList)$Original)[1], 3)
}