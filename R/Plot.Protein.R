Plot.Protein <-
function(graph, path, vertex.size=5, color.palette="heat"){
  
  if(color.palette == "heat"){
    pal <- rev(heat.colors(vcount(graph)))
  }else if(color.palette == "gray"){
    pal <- rev(gray(seq(1/2, 1, length=vcount(graph))))
  }else if(color.palette == "topo"){
    pal <- rev(topo.colors(vcount(graph)))
  }else if(color.palette == "cm"){
    pal <- rev(cm.colors(vcount(graph)))
  }
  
  
  dist <- c()
  for(i in 1:vcount(graph)){
    current.pos <- as.numeric(V(graph)[i]$label)
    dist <- append(dist, which(path == (current.pos)))
  }  
  #start.vertex = which(V(graph)$label == 1)
  #dist<- shortest.paths(graph, start.vertex) +1
  
  tkplot(graph, layout = layout.circle, vertex.color = pal[dist], vertex.size =7) 
}
