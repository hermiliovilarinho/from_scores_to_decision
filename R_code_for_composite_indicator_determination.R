## Determination of Composite Indicator 
## Based on: BoD Input Oriented + ARI Lower Bound weight restrictions
## Weights Formulation


#Data Reading

#install.packages("readxl")
library(readxl)
data <- read_excel("Data.xlsx", sheet="PT",range = "A1:AM491",col_names=TRUE)


#Reading Firm ID
rownames <- as.data.frame(read_excel("Data.xlsx", sheet="PT",range = "A2:A491",col_names=FALSE))

y<-data.frame(cbind(data$`S&O`,data$`SI`,data$`SO`,data$`SP`,data$`DDS`,data$`HR`))

# Select specific columns from the 'data' object to create 'exog_data'
exog_data <- data[, c("id", "Dimension", "Sector", "Age_Categ", "Op_Inc_per_Worker", "Result_per_Worker")]

exog_var<- data.frame(exog_data)
exog_var<- data.frame(x1=factor(exog_var$'Dimension'),
                      x3=factor(exog_var$'Sector'),x4=factor(exog_var$'Age_Categ'),
                      x5=as.numeric(exog_var$Op_Inc_per_Worker), x6=as.numeric(exog_var$Result_per_Worker))

dat <- data.frame(y)

#Replace null values by a very small value (0.001)
dat[dat == 0] <- 0.01


#averages for the indicators
avg<-colMeans(dat)# calculating averages to use in weight restrictions


n <- length(y[,1]) # of units
s <- length(y[1,])  # of indicators


###################################################################

#NON ROBUST UNCONDITIONAL COMPOSITE INDICATOR

##################################################################


# To store the results
optimum<-matrix(nrow=length(y[,1]),ncol=2)   # to store the results
optimum <- as.data.frame(optimum)
rownames(optimum) <- rownames[,1]
colnames(optimum) <- c("h","CI")

weights<-matrix(nrow=length(y[,1]),ncol=(s))   # to store the results
weights <- as.data.frame(weights)
rownames(weights) <- rownames[,1]
colnames(weights) <- c("u_S&O","u_SI","u_SO","u_SP","u_DDS","u_HR")

#create column labels for the coefficient matrices
col_label <- as.data.frame(1:(s))
colnames(col_label)=NULL
col_label <- t(col_label)

#Setting weights' restrictions
lb <-0.05

lower_bound <- matrix(nrow=s,ncol=s)   # to create the matrix as additional constraints
lower_bound <- as.data.frame(lower_bound)

data1=-lb*avg
data2=(1-lb)*avg # data on the diagonal

for (j in 1:s)
{
  lower_bound[j,] <- c(matrix(data1,ncol=s,nrow=1))
  for (k in 1:s)
  {
    if (j==k)  {
      lower_bound[j,j] <- data2[j]
    }
  }
}


colnames(lower_bound)<- col_label

# Calculating the efficiency and CI for each DMU i
i<-1
for (i in 1:n)
{
  yk <- dat[i,]   # observed desirable indicator
  
  # Setting the coefficients of decision variables
  obj <- yk
  
  # Constraint Matrix : 
  
  mat1 <-  cbind(dat)
  colnames(mat1)<- col_label
  
  
  mat <- rbind(mat1,lower_bound)
  #mat<- mat1
  
  dir <- c(
    matrix(data="<=",nrow=n,ncol =1),
    matrix(data=">=",nrow=s,ncol =1)
  )
  
  rhs <- c(
    matrix(data=1,nrow=n,ncol=1),
    matrix(data=0,nrow=s,ncol=1)
  )
  
  types <-  c(matrix(data="C",ncol=s,nrow=1))
  
  #install.packages("Rglpk")
  library(Rglpk)
  #install.packages("lpSolve")
  library(lpSolve)
  
  ans <- Rglpk_solve_LP(obj, mat, dir, rhs, types, max=TRUE)
  
  optimum[i,1] <- ans$optimum
  optimum[i,2] <- ans$optimum
  
  
  
  wgt <- matrix(ans$solution)
  wgt <- t(wgt)
  t(wgt)
  weights[i,1:s]<-wgt[,c(1:s)]
}


optimum <- cbind(Company=rownames(optimum),optimum)
weights <- cbind(Company=rownames(weights),weights)

##################################################################

# ROBUST UNCONDICIONAL COMPOSITE INDICATOR

#################################################################
linhas <- vector("list", length = length(seq(50, 450, by = 20))) # para calculo do valor de m
a <- 1

B = 2000 #number of bootstrap replicates


m = 200 #number of units included each time

# To store the results
optimum_rob<-matrix(nrow=length(y[,1]),ncol=2)   # to store the results
optimum_rob <- as.data.frame(optimum_rob)
rownames(optimum_rob) <- rownames[,1]
colnames(optimum_rob) <-  c("h_Rob","CI_Rob")


optimum_b <- matrix(NA, ncol = 2, nrow = B)
weights_b<- matrix(NA,ncol=s,nrow=B)


weights_rob<-matrix(nrow=length(y[,1]),ncol=(s))   # to store the results
weights_rob <- as.data.frame(weights_rob)
rownames(weights_rob) <- rownames[,1]
colnames(weights_rob) <- c("u_S&O","u_SI","u_SO","u_SP","u_DDS","u_HR")


for (i in 1:n) {
  
  print(i)
  
  yk <- dat[i,]   # observed desirable indicator
  
  # Setting the coefficients of decision variables
  
  obj <- yk
  
  
  for (bt in 1:B) {
    
    robust_sample <- sample(n, m, replace = TRUE)  #pick m random units out of the sample of n_rob units
    
    yy <- dat[robust_sample,]
    nn <- nrow(yy)
    ones <- matrix(data =1,nrow = nn, ncol=1)
    mat1 <- cbind(yy)
    colnames(mat1)<- col_label
    
    mat <- rbind(mat1,lower_bound)
    
    dir <- c(
      matrix(data="<=",nrow=nn,ncol =1),
      matrix(data=">=",nrow=s,ncol =1))
    
    rhs <- c(
      matrix(data=1,nrow=nn,ncol=1),
      matrix(data=0,nrow=s,ncol=1)
    )
    
    types <-  c(matrix(data="C",ncol=s,nrow=1))
    
    ans <- Rglpk_solve_LP(obj, mat, dir, rhs, types, max=TRUE)
    
    optimum_b[bt,1] <- ans$optimum
    optimum_b[bt,2] <-ans$optimum
    
    weights_b [bt,]<-ans$solution[1:s]
    
  }
  
  optimum_rob[i,1] <- mean(optimum_b[,1])
  
  optimum_rob[i,2] <- mean(optimum_b[,2])
  
  weights_rob[i,]<-apply(weights_b,2,mean)
  
}

library(beepr)
beep()

optimum_rob <- cbind(Company=rownames(optimum_rob),optimum_rob)
weights_rob <- cbind(Company=rownames(weights_rob),weights_rob)


##################################################################
# ROBUST CONDITIONAL COMPOSITE INDICATOR — BALANCED VERSION
# Better differentiation across contextual variables
##################################################################

# --- Preprocessing (less aggressive) ---
# Scale only numeric features, keep factors as-is for np
library(dplyr)
library(recipes)

rec <- recipe(~ ., data = exog_var) %>%
  step_center(all_numeric_predictors()) %>%
  step_scale(all_numeric_predictors())

prep_rec <- prep(rec, training = exog_var)
exog_scaled <- as.data.frame(bake(prep_rec, new_data = exog_var))

# --- KDE using mixed data: np handles it well ---
library(np)

# Estimate KDE bandwidth with cross-validation (adaptive)
bw <- npudensbw(dat = exog_scaled, ckertype = "epanechnikov", bwmethod = "cv.ml")

# Storage
optimum_rob_cond <- matrix(nrow = n, ncol = 2)
optimum_rob_cond <- as.data.frame(optimum_rob_cond)
rownames(optimum_rob_cond) <- rownames[, 1]
colnames(optimum_rob_cond) <- c("h_Rob_Cond", "CI_Rob_Cond")

weights_rob_cond <- matrix(nrow = n, ncol = s)
weights_rob_cond <- as.data.frame(weights_rob_cond)
rownames(weights_rob_cond) <- rownames[, 1]
colnames(weights_rob_cond) <- colnames(dat)

# --- Main Loop ---
for (i in 1:n) {
  print(i)
  
  # KDE conditional density for unit i
  tdat_i <- exog_scaled[rep(i, 1), , drop = FALSE]
  
  # Conditional KDE
  fhat <- npudens(
    bws = bw,
    cykertype = "epanechnikov",
    cxkertype = "epanechnikov",
    tdat = tdat_i,
    edat = exog_scaled
  )
  z_full <- fhat$dens
  z_full[is.na(z_full)] <- 0
  
  # Step 5: Define sampling pool
  if (sum(z_full > 0) >= 30) {
    other_indexes <- setdiff(1:n, i)
    z_conditional <- z_full[other_indexes]
    conditional_sample_pool <- other_indexes
  } else {
    z_conditional <- z_full
    conditional_sample_pool <- 1:n
  }
  
  # Step 6: Normalize again (after possible NA or zeros)
  if (sum(z_conditional, na.rm = TRUE) == 0) {
    z_conditional <- rep(1, length(z_conditional)) / length(z_conditional)
  } else {
    z_conditional <- z_conditional / sum(z_conditional)
  }
  
  # Step 7: Bootstrap resampling + optimization
  yk <- dat[i, ]
  obj <- yk
  optimum_b <- matrix(NA, nrow = B, ncol = 2)
  weights_b <- matrix(NA, nrow = B, ncol = s)
  
  for (bt in 1:B) {
    attempts <- 0
    repeat {
      conditional_sample <- sample(conditional_sample_pool, m, replace = TRUE, prob = z_conditional)
      if (length(unique(conditional_sample)) >= m * 0.5 || attempts >= 10) break
      attempts <- attempts + 1
    }
    
    yy <- dat[conditional_sample, ]
    mat1 <- cbind(yy)
    colnames(mat1)<- col_label
    mat <- rbind(mat1, lower_bound)
    
    dir <- c(rep("<=", nrow(yy)), rep(">=", s))
    rhs <- c(rep(1, nrow(yy)), rep(0, s))
    types <- rep("C", s)
    
    ans <- Rglpk_solve_LP(obj, mat, dir, rhs, types, max = TRUE)
    optimum_b[bt, ] <- ans$optimum
    weights_b[bt, ] <- ans$solution[1:s]
  }
  
  # Step 8: Store averages
  optimum_rob_cond[i, ] <- colMeans(optimum_b, na.rm = TRUE)
  weights_rob_cond[i, ] <- colMeans(weights_b, na.rm = TRUE)  
  # Optional warning
  if (optimum_rob_cond[i, 2] > 2) {
    cat("⚠️ High CI_Rob_Cond for unit", i, ":", optimum_rob_cond[i, 2], "\n")
  }
}

library(beepr)
beep()

# Final output
optimum_rob_cond <- cbind(Company = rownames(optimum_rob_cond), optimum_rob_cond)
weights_rob_cond <- cbind(Company = rownames(weights_rob_cond), weights_rob_cond)

ckertype = "epanechnikov", oxkertype = "liracine")

model_single <- npreg(bws = bw_single, gradients = TRUE)

for (qname in names(quartiles)) {
  xq_val <- quartiles[[qname]]
  
  png_filename <- file.path(output_dir, paste0("Plot_", qname, "_", var_name, ".png"))
  
  png(filename = png_filename, width = 1200, height = 900, res = 150)
  
  par(mfrow = c(1, 1), 
      mar = c(4, 4, 3, 1), 
      oma = c(0, 0, 0, 0), 
      font.lab = 2, cex.lab = 1.3,
      font.main = 1, cex.main = 1.2,
      font.axis = 3, cex.axis = 1)
  
  if (qname == "Median") {
    npplot(bws = bw_single, xq = xq_val,
           main = paste(qname, "-", var_name),
           ylab = "Score ratio",
           plot.errors.method = "bootstrap",
           plot.errors.center = "bias-corrected",
           plot.errors.type = "quantiles")
  } else {
    npplot(bws = bw_single, xq = xq_val,
           main = paste(qname, "-", var_name),
           ylab = "Score ratio")
  }
  
  dev.off()
}
}


###################################################################################

#Export results to Excel

###################################################################################


# Add results to the data table
data <- cbind(data,
              optimum[[3]],
              optimum_rob[[3]],
              optimum_rob_cond[[3]])

colnames(data)[(ncol(data)-2):ncol(data)] <- c("CI_std", "CI_rob", "CI_rob_cond") #giving names to the columns

cor_ci<-cor(data[, c("CI_std", "CI_rob", "CI_rob_cond")], method = "spearman")

ranks <- data.frame(
  rank_CI_std      = rank(-data$CI_std,      ties.method = "min"),
  rank_CI_rob      = rank(-data$CI_rob,      ties.method = "min"),
  rank_CI_rob_cond = rank(-data$CI_rob_cond, ties.method = "min"),
  rank_Class_Ger   = rank(-data[["Classif"]], ties.method = "min")
)

cor_ranks <- cor(ranks, method = "spearman")  # ou spearman/kendall, se desejar
print(round(cor_ranks, 3))


#install.packages("openxlsx")
library(openxlsx)
#l <- list("CI" = optimum, "Weights" = weights,"CI_Rob" = optimum_rob, "Weights_Rob"=weights_rob)
l <- list("Summary"= data,"CI" = optimum, "Weights" = weights,"CI_Rob" = optimum_rob, 
          "Weights_Rob"=weights_rob,"CI_Rob_Cond" = optimum_rob_cond, "Weights_Rob_Cond"=weights_rob_cond,
          "Betas"=summary(betas2), "SE"=summary(se2),"Ranks"= ranks,
          "Corr_CI" = cor_ci, "Corr_Rank" = cor_ranks)

write.xlsx(l, file = "CI_Results.xlsx")

