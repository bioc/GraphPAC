test_Find.TSP.Path  <- function() {
  Residue <- c("MET","THR","GLU","TYR")
  Can.Count<- c(1,2,3,4)
  SideChain <- c("A","A","A","A")
  XCoord <- c(62.935,63.155,65.289,64.899)
  YCoord <- c(97.579,95.525,96.895,96.220)
  ZCoord <- c(30.223,27.079,24.308,20.615)
  
  mutation.matrix <- matrix(data = c(1,0,0,1,0,0,1,1,1,0,1,1,0,1,1,0,1,1,0,1),nrow = 5, ncol = 4)
  PositionList <- data.frame(Residue = Residue, Can.Count = Can.Count, SideChain = SideChain, XCoord = XCoord, YCoord = YCoord, ZCoord = ZCoord)
  
  checkEquals(length(Find.TSP.Path(PositionList,mutation.matrix)), 4)
  checkEquals(length(Find.TSP.Path(PositionList,mutation.matrix)$candidate.path), dim(mutation.matrix)[2])
}