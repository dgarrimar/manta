
## Generate response variables
set.seed(1234)
n <- 100
biomarkers <- cbind(rnorm(n)*sqrt(2), 
      runif(n)*pi, 
      rbeta(n, shape1 = runif(1), shape2 = runif(1)),
      rchisq(n, df = round(runif(1)*100)),
      rf(n, df1 = round(runif(1)*10), df2 = round(runif(1)*10) ))
colnames(biomarkers) <- paste0("biomarker", 1:ncol(biomarkers))

## Generate explanatory variables
patients <- data.frame("age" = runif(100)*10,
                       "gender" = sample(c("male", "female"), size = 100, replace = TRUE),
                       "status" = sample(c("healthy", "mild", "severe"), size = 100, replace = TRUE))
patients$status <- as.ordered(patients$status)

## Introduce NA's randomly in both datasets
k <- ceiling(0.02*n)
for (i in 1:k){
  a <- sample(1:n, size = 1)
  b <- sample(1:ncol(biomarkers), size = 1)
  biomarkers[a, b] <- NA
}

for (i in 1:k){
  a <- sample(1:n, size = 1)
  b <- sample(1:ncol(patients), size = 1)
  patients[a, b] <- NA
}

write.table(biomarkers, "biomarkers.tsv", sep = "\t", quote = FALSE)
write.table(patients, "patients.tsv", sep = "\t", quote = FALSE)

devtools::use_data(biomarkers, overwrite = TRUE)
devtools::use_data(patients, overwrite = TRUE)
