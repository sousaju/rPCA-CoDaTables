#############################
# Packages and data loading #
#############################
require(compositions)
require(robCompositions)
require(robustbase)
data <- read.csv2("unemployment.csv",sep =",")   #unemployment data set 
factor1lev <- rownames(data) <- c("15-24","25-39","40-54","55+")   #age groups
factor2lev <- c("M","F")  #gender

n <- length(data[1,])/2   #number of tables
I2 <- 4   #number of rows
J2 <- 2   #number of columns

# create separate table for each country
Q2 <- as.list(rep(0,n))
for (i in 1:n)
{
  q  <- matrix(0, I2, J2)
  for (j in 1:I2)
  {
    for (k in 1:J2)
    {q[j,k] <- data[j,J2*(i-1)+k]}
    Q2[[i]] <- q} 
}

#######################################
# Independence and interaction tables #
#######################################
# projection row orthogonal
row <- as.list(rep(0,n))
for (i in 1:n)
{
  r <- matrix(0, I2, J2)
  for (j in 1:I2)
  {
    for (k in 1:J2)
    {r[j,k] <- prod(Q2[[i]][j,])^(1/J2)}
  } 
  row[[i]] <- r
}

# projection col orthogonal
col <- as.list(rep(0,n))
for (i in 1:n)
{
  c <- matrix(0, I2, J2)
  for (j in 1:I2)
  {
    for (k in 1:J2)
    {c[j,k] <- prod(Q2[[i]][,k])^(1/I2)}
  } 
  col[[i]] <- c
}

# independence table
ind <- as.list(rep(0,n))
for (i in 1:n)
{
  d <- matrix(0, I2, J2)
  for (j in 1:I2)
  {
    for (k in 1:J2)
    {d[j,k] <- row[[i]][j,k]*col[[i]][j,k]}
  } 
  ind[[i]] <- d
}

# closure operation
ind_prop <- as.list(rep(0,n))
for (i in 1:n)
{
  dp <- matrix(0, I2, J2)
  for (j in 1:I2)
  {
    for (k in 1:J2)
    {dp[j,k] <- ind[[i]][j,k]/sum(ind[[i]])}
  } 
  ind_prop[[i]] <- dp
}

# interaction table
int <- as.list(rep(0,n))
for (i in 1:n)
{
  t <- matrix(0, I2, J2)
  for (j in 1:I2)
  {
    for (k in 1:J2)
    {t[j,k] <- Q2[[i]][j,k]/ind_prop[[i]][j,k]}
  } 
  int[[i]] <- t
}

# closure operation
int_prop <- as.list(rep(0,n))
for (i in 1:n)
{
  tp <- matrix(0, I2, J2)
  for (j in 1:I2)
  {
    for (k in 1:J2)
    {tp[j,k] <- int[[i]][j,k]/sum(int[[i]])}
  } 
  int_prop[[i]] <- tp
}

# since J<I, transpose the tables
if(J2<I2)
{
  for(i in 1:n)
  {
    ind_prop[[i]] <- t(ind_prop[[i]])
    int_prop[[i]] <- t(int_prop[[i]])
  }
  I <- J2
  J <- I2
  pom=factor1lev
  factor1lev=factor2lev
  factor2lev=pom
  rm(pom)
} else {
  I <- I2
  J <- J2
}

# vectorize the tables
ind_vec <- matrix(0,n,I*J)
for(i in 1:n)
{
  for(j in 1:I)
    for(k in 1:J)
      ind_vec[i,J*(j-1)+k] <- ind_prop[[i]][j,k]
}

int_vec <- matrix(0,n,I*J)
for(i in 1:n)
{
  for(j in 1:I)
    for(k in 1:J)
      int_vec[i,J*(j-1)+k] <- int_prop[[i]][j,k]
}


# rownames and colnames
colnames(ind_vec) <- c(1:(I*J))
i <- 1
  for (j in 1:length(factor1lev))
    for (k in 1:length(factor2lev)) 
      {
      colnames(ind_vec)[i] <- paste(factor1lev[j],factor2lev[k],sep="")
      i <- i+1
    }
colnames(int_vec) <- colnames(ind_vec)

rownames(ind_vec) = c(1:n)
k=1
for(i in seq(1,2*n,2))
{
  rownames(ind_vec)[k]=colnames(data)[i]
  k=k+1
}
rownames(int_vec) <- rownames(ind_vec) 

# pivot coordinates
z_int1 <- sqrt(3/8)*log((int_vec[,6]*int_vec[,7]*int_vec[,8])^(1/3)*int_vec[,1]/((int_vec[,2]*int_vec[,3]*int_vec[,4])^(1/3)*int_vec[,5]))
z_int2 <- sqrt(1/3)*log((int_vec[,7]*int_vec[,8])^(1/2)*int_vec[,2]/((int_vec[,3]*int_vec[,4])^(1/2)*int_vec[,6]))
z_int3 <- (1/2)*log(int_vec[,3]*int_vec[,8]/(int_vec[,4]*int_vec[,7]))
z_int <- cbind(z_int1,z_int2,z_int3)

z_r1 <- sqrt(2)*1/4*log(ind_vec[,1]*ind_vec[,2]*ind_vec[,3]*ind_vec[,4]/(ind_vec[,5]*ind_vec[,6]*ind_vec[,7]*ind_vec[,8]))
z_c1 <- sqrt(3/2)*1/2*log(ind_vec[,1]*ind_vec[,5]/((ind_vec[,2]*ind_vec[,3]*ind_vec[,4]*ind_vec[,6]*ind_vec[,7]*ind_vec[,8])^(1/3)))
z_c2 <- sqrt(4/3)*1/2*log(ind_vec[,2]*ind_vec[,6]/((ind_vec[,3]*ind_vec[,4]*ind_vec[,7]*ind_vec[,8])^(1/2)))
z_c3 <- 1/2*log(ind_vec[,3]*ind_vec[,7]/(ind_vec[,4]*ind_vec[,8]))
z_ind <- cbind(z_r1,z_c1,z_c2,z_c3)

z <- cbind(z_ind,z_int)

# transformation matrices 
V_int <- cbind(c(sqrt(3/8),-sqrt(3/8)*1/3,-sqrt(3/8)*1/3,-sqrt(3/8)*1/3,-sqrt(3/8),sqrt(3/8)*1/3,sqrt(3/8)*1/3,sqrt(3/8)*1/3),c(0,sqrt(1/3),-sqrt(1/3)*1/2,-sqrt(1/3)*1/2,0,-sqrt(1/3),sqrt(1/3)*1/2,sqrt(1/3)*1/2),c(0,0,1/2,-1/2,0,0,-1/2,1/2))
V_ind <- cbind(c(1/4*sqrt(2),1/4*sqrt(2),1/4*sqrt(2),1/4*sqrt(2),-1/4*sqrt(2),-1/4*sqrt(2),-1/4*sqrt(2),-1/4*sqrt(2)),c(1/2*sqrt(3/2),-1/6*sqrt(3/2),-1/6*sqrt(3/2),-1/6*sqrt(3/2),1/2*sqrt(3/2),-1/6*sqrt(3/2),-1/6*sqrt(3/2),-1/6*sqrt(3/2)))
V_ind <- cbind(V_ind,c(0,1/2*sqrt(4/3),-1/4*sqrt(4/3),-1/4*sqrt(4/3),0,1/2*sqrt(4/3),-1/4*sqrt(4/3),-1/4*sqrt(4/3)),c(0,0,1/2,-1/2,0,0,1/2,-1/2))

V <- cbind(V_ind,V_int)


###########
# Biplots #
###########
pdf("Biplots_Unemployment2.pdf",height=5,width=9)
par(mfrow=c(1,2))

# biplots for the entire compositional tables
# robust biplot
detMCD=c(1:1000)
for (i in 1:1000)
{
  set.seed(i)
  detMCD[i]=covMcd(z)$crit
}
set.seed(which.min(exp(detMCD))) #set.seed(1)

pcarob=princomp(z,cov=covMcd(z))
#summary(pcarob)
scorob=pcarob$scores
rownames(scorob)=rownames(ind_vec)
loarob=V%*%pcarob$loadings
rownames(loarob)=colnames(ind_vec)

biplot(x=scorob,y=loarob,xlab="PC 1 (52.29%)",ylab="PC 2 (36.83%)",cex=0.5)
title("Robust biplot for compositional table")

# classical biplot
pcacla=princomp(z)
#summary(pcacla)
scocla=pcacla$scores
rownames(scocla)=rownames(ind_vec)
loacla=V%*%pcacla$loadings
rownames(loacla)=colnames(ind_vec)

biplot(x=scocla,y=loacla,xlab="PC 1 (53.41%)", ylab="PC 2 (24.72%)",cex=0.5)
title("Classical biplot for compositional table")


# biplots for the independence tables
# robust biplot
detMCD=c(1:1000)
for (i in 1:1000)
{
  set.seed(i)
  detMCD[i]=covMcd(z_ind)$crit
}
set.seed(which.min(exp(detMCD))) #set.seed(4)

pcarob=princomp(z_ind,cov=covMcd(z_ind))
#summary(pcarob)
scorob=pcarob$scores
rownames(scorob)=rownames(ind_vec)
loarob=V_ind%*%pcarob$loadings
rownames(loarob)=colnames(ind_vec)

biplot(x=scorob,y=loarob,xlab="PC 1 (51.15%)",ylab="PC 2 (39.12%)",cex=0.5)
title("Robust biplot for independence table")

# classical biplot
pcacla=princomp(z_ind)
#summary(pcacla)
scocla=pcacla$scores
rownames(scocla)=rownames(ind_vec)
loacla=V_ind%*%pcacla$loadings
rownames(loacla)=colnames(ind_vec)

biplot(x=scocla,y=loacla,xlab="PC 1 (57.32%)", ylab="PC 2 (26.64%)",cex=0.5)
title("Classical biplot for independence table")

# biplots for the interaction tables
# robust biplot 
detMCD=c(1:1000)
for (i in 1:1000)
{
  set.seed(i)
  detMCD[i]=covMcd(z_int)$crit
}
set.seed(which.min(exp(detMCD))) #set.seed(1)

pcarob=princomp(z_int,cov=covMcd(z_int))
#summary(pcarob)
scorob=pcarob$scores
rownames(scorob)=rownames(int_vec)
loarob=V_int%*%pcarob$loadings
rownames(loarob)=colnames(int_vec)

biplot(x=scorob,y=loarob,xlab="PC 1 (87.22%)",ylab="PC 2 (8.43%)",cex=0.5)
title("Robust biplot for interaction table")

# classical biplot
pcacla=princomp(z_int)
#summary(pcacla)
scocla=pcacla$scores
rownames(scocla)=rownames(int_vec)
loacla=V_int%*%pcacla$loadings
rownames(loacla)=colnames(int_vec)

biplot(x=scocla,y=loacla,xlab="PC 1 (77.31%)", ylab="PC 2 (17.05%)",cex=0.5)
title("Classical biplot for interaction table")

dev.off()


