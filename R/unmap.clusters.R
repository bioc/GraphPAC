unmap.clusters <-
function (results, path, mutation.data, position.data ){
  mutation.data.culled<-mutation.data[,position.data[,2]]
  colnames(mutation.data.culled)<-mapply(function(x){substring(x,2)}, colnames(mutation.data.culled))
  for(i in 1:dim(results)[1]){
    results[i,2] <- path[results[i,2]]
    results[i,3] <- path[results[i,3]]
    
    if(results[i,2] > results[i,3] ){
      temp <- results[i,2]
      results[i,2]<- results[i,3]
      results[i,3]<- temp
    }
    
    
    start.sum <- which(colnames(mutation.data.culled) == results[i,2])
    end.sum <- which(colnames(mutation.data.culled ) == results[i,3])
    results[i,4] <- sum(mutation.data.culled[,start.sum:end.sum])
    results[i,1] <- results[i,3]-results[i,2]+1
  }
  return(results)
}
