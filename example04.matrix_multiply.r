# create vector of 100M elems
vec <- 1:100000000

# create a matrix from the vector
mat <- matrix(vec,10000)

#start timer
start.time <- Sys.time()

# do matrix mult using R operator
mymult <- mat %*% mat

# end timer and print elapsed time
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

# check size
dim(mymult)
object.size(mymult)

# exit the r program 
# important when running as script)
q()
