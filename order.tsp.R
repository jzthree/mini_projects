require(TSP)
order.tsp<-function(d){
  # random: solve heatmap ordering as minimum Hamiltonian path 
  # problem which is solved as a special case of TSP. 
  # Take distance as input
  
  d= rbind(0,cbind(0,as.matrix(d)))
  concorde_path('/Genomics/ogtr04/jzthree/bin/')
  solved=solve_TSP(TSP(d),method = 'concorde')
  order.pre=as.numeric(solved)
  order = order.pre[2:length(order.pre)]-1
  return(order)
}

