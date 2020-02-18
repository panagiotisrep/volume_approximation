library("volesti")


for (d in c(25,50,75,100,125,150,175,200)) {
  str=sprintf("cube_%d.ine",d)
  P = volesti::file_to_polytope(path = str)
  start.time <- Sys.time()
  points <- sample_points(P, N=200000, distribution="uniform", random_walk = "RDHR", walk_length = 1)

  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(time.taken)
}


for (d in c(25,50,100,125,150)) {
  str=sprintf("cube_%d.ine",d)
  P = volesti::file_to_polytope(path = str)
  vol=0
  start.time <- Sys.time()
  for (i in 1:1) {
    vol=vol+volume(P,algo = 'CB',random_walk = 'CDHR')
  }
  vol = vol / 1
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(time.taken)
  print(vol)
}
