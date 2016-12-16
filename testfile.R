# some stuff to mess
n <- 100
m <- 1000
X <- matrix(sample(c(0, 1), size = n * m, replace = TRUE), nrow = n)
p <- colMeans(X)

print('you are a fool')

##############################
n <- 100
m <- 1000
X <- matrix(sample(c(0, 1), size = n * m, replace = TRUE), nrow = n)
p <- colMeans(X)


microbenchmark::microbenchmark(times = 1L, 
u <- LD_mult(X = X, matr = comb_all(1:m),
              is_phased = TRUE, any_na = TRUE, check = FALSE)
)




library(LDtools)
pos <- 1:20
comb_sliding(pos, start = 0, width = 2, advance = 1)



n <- 100L
m <- 5L
x <- matrix(sample(c(0, 1), size = n * m, replace = TRUE), nrow = n)
p <- colMeans(x)
LDtools:::center_matrix(x, p)
pos = 1:m
.LD_list(X = x, p = p, list = comb_all(pos), F, T)


microbenchmark::microbenchmark(times = 10L, 
LDtools:::center_matrix(x, p),
LDtools:::center_matrix2(x, p)
)





n1 <- 100
n2 <- 1000
pos1 <- sort(runif(n1))
pos2 <- sort(runif(n2))
all(dplyr::near(
  apply(abs(outer(X = pos1, Y = pos2, FUN = `-`)), 1L, which.min),
  .comb_nearest(1:n1, 1:n2, pos1, pos2)[[1L]][, 2]
))


pairs_collapse <- function(x) {
  stopifnot((n <- length(x)) %% 2 == 0)
  ix <- seq.int(1L, n, by = 2L)
  x[ix] + x[ix + 1L]
}
# Simulate correlation phase genotypes
n <- 10000L
o <- n %/% 2L
x <- sample(c(0, 1), size = n, replace = TRUE)
y <- x
ix <- sample.int(n, size = o)
y[ix] <- sample(y[ix])
px <- mean(x)
py <- mean(y)
Dmin = max(-px * py, -(1-px) * (1-py));
Dmax = min(px * (1-py), (1-px) * py);

get_LD_pair(x-px, y-py, px, py, FALSE, is_phased=T) %>% unlist()
D <- crossprod(x-px, y-py) / (n - 1)
D / sqrt(px * (1-px) * py * (1-py)) # r def
D / if (D > 0) Dmax else Dmin # Dprime
cor(x,y)

# Simulate related unphased genotypes by collapsing pairs
xu <- pairs_collapse(x)
yu <- pairs_collapse(y)
get_LD_pair(x-2*px, y-2*py, px, py, FALSE, is_phased=F)


############################################
library(genetics)
library('LDtools')
x <- genotype( c('C/A', 'C/A', 'C/C', 'C/A', 'C/C', 'C/A', 'C/A', 'C/A',
                  'C/A', 'C/C', 'C/A', 'A/A', 'C/A', 'A/A', 'C/A', 'C/C',
                  'C/A', 'C/A', 'C/A', 'A/A') )

xt <- drop(transform_geno(matrix(x), sep = '/')$geno)

y <- genotype( c('T/A', 'T/A', 'T/T', 'T/A', 'T/T', 'T/A', 'T/A', 'T/A',
                  'T/A', 'T/T', 'T/A', 'T/T', 'T/A', 'T/A', 'T/A', 'T/T',
                  'T/A', 'T/A', 'T/A', 'T/T') )

yt <- drop(transform_geno(matrix(y), sep = '/')$geno)

# Compute LD on a single pair
xt
yt[5] <- NA
.LD(xt, yt, px = mean(xt)/2, py = mean(yt)/2,
    is_phased = FALSE, any_na = FALSE)

microbenchmark::microbenchmark(times = 1000L, 
genetics::LD(x, y), 
.LD(xt, yt, px = mean(xt)/2, py = mean(yt)/2, is_phased = FALSE, any_na = FALSE),
LDtools::LD(xt, yt, is_phased = FALSE, any_na = FALSE, check = TRUE)
)

xt <- rep(0, length(xt))
.LD(xt, yt, px = mean(xt)/2, py = mean(yt)/2, is_phased = FALSE, any_na = FALSE)

x <- c(0,0,0,0)
y <- c(1,1,1,0)
.LD(x-mean(x),y-mean(y), mean(x), mean(y), F, T)

######################################################################33
set.seed(123L)
nucleotides <- c("A", "C", "G", "T")
m <- 10L
n <- 1000L
prob_na <- 0.2

geno <- gen_geno(nucleotides, n, m, prob_na)

new_geno <- data.frame(lapply(as.data.frame(geno), genetics::as.genotype))
tmp <- transform_geno(geno, sep  = '/')
microbenchmark::microbenchmark(times = 1L,
LD_mat <- genetics::LD(new_geno)[["R^2"]],
out <- .LD_list(tmp$geno, p = tmp$p, pos = 1:m,
                list = list(comb_all(n = m)[[1]] - 1),
                any_na = TRUE, is_phased = FALSE)
)
        
a <- gdata::lowerTriangle(t(LD_mat))
b <- out[,"r2"]

cor(a, b)
plot(a, b)

u <- gdata::upperTriangle(cor(tmp$geno, use = 'pairwise.complete.ob'))

cor(u,b)

data("wheat", package = "BGLR")
X <- do.call(cbind, replicate(n = 1L, wheat.X, simplify = FALSE))

n <- ncol(X)
pos <- sort(runif(n, 0, 100))
p <- colMeans(X)
X <- scale(X, center = TRUE, scale = FALSE)
sdv <- apply(X, 2L, sd)
keep <- sdv > sqrt(.Machine$double.eps)
X <- X[, keep]
pos <- pos[keep]
p <- p[keep]

microbenchmark::microbenchmark(times = 1,
lst <- comb_dist_wind(pos, min_dist = 0, max_dist = 100),
u <- .LD_list(X, p, pos, lst, F, T),
print(nrow(.LD_window(X, p, pos, 0, 100, is_phased = T, any_na = F))),
gdata::upperTriangle(cor(X)^2, byrow = TRUE)
)

lst <- comb_dist_adj(ncol(X))


tp <- t(combn(x = 1:500, m = 2, simplify = TRUE))
microbenchmark::microbenchmark(times=1,
u <- .LD_list(X, p, pos, list(tp,tp,tp,tp,tp,tp), F, T)
)
dim(u)

microbenchmark::microbenchmark(times=1,
LD1a <- (crossprod(scale(X, center = p, scale = sqrt(p * (1-p))))/nrow(X))^2,
LD1b <- cor(X)^2,
LD2 <- pair_LD_phased2(X, p, pos, min_dist = 0.0, max_dist = 100.0)
# LD2 <- test(X)
)

cor(LD2$r2, gdata::upperTriangle(cor(X)^2, byrow = TRUE))


sum((LD2 - LD1b)^2)

Xs2 <- test(X)
plot(colMeans(Xs2))
hist(matrixStats::colSds(Xs2))

plot(colMeans(Xs1))
hist(matrixStats::colSds(Xs1))

str(LD2)

sd <- apply(X, 2L, sd)




microbenchmark::microbenchmark(times=1,
uu <- pair_LD_phased(X, p, pos, 0,  5, phased = T, has_NA = T),
co <-  gdata::upperTriangle(cor(X)^2, byrow = TRUE),
u <- cor(X),
uu <- test(X)
)
dim(uu)
hist(uu$dist)

microbenchmark::microbenchmark(times=1,
# tr <- pair_LD(X, p, pos, 0.0, 1.0, phased = TRUE, hasNA = TRUE),
# u <- pair_LD(X, p, pos, 0.0, 1.0, phased = TRUE, hasNA = FALSE),
co <-  gdata::upperTriangle(cor(X)^2, byrow = TRUE),
colMeans(X),
unique(c(X))
)








cor(u$r2, co)

Xorg <- Xorg[-nrow(Xorg), ]
s <- seq.int(1L, nrow(Xorg) - 1L, by = 2L)
Xup <- Xorg[s,] + Xorg[s + 1,]
microbenchmark::microbenchmark(times=1,
uup <- pair_LD(Xup, p, pos, 0.0, 1.0, phased = FALSE, hasNA = TRUE)
)

cor(cbind(u$r2, uup$r2, co))

rr <-sapply(u, nrow)

x <- list(a = 5, b =7)
rbind_frames(x)


u <- 5.0
pair_LD(u)

library(genetics)
library(data.table)
library("Rcpp")


# g1 <- genotype( c('T/A',    NA, 'T/T',    NA, 'T/A',    NA, 'T/T', 'T/A',
#                   'T/T', 'T/T', 'T/A', 'A/A', 'T/T', 'T/A', 'T/A', 'T/T',
#                   NA, 'T/A', 'T/A',   NA) )
# 
# g2 <- genotype( c('C/A', 'C/A', 'C/C', 'C/A', 'C/C', 'C/A', 'C/A', 'C/A',
#                   'C/A', 'C/C', 'C/A', 'A/A', 'C/A', 'A/A', 'C/A', 'C/C',
#                   'C/A', 'C/A', 'C/A', 'A/A') )
# 
# g3 <- genotype( c('T/A', 'T/A', 'T/T', 'T/A', 'T/T', 'T/A', 'T/A', 'T/A',
#                   'T/A', 'T/T', 'T/A', 'T/T', 'T/A', 'T/A', 'T/A', 'T/T',
#                   'T/A', 'T/A', 'T/A', 'T/T') )

g1 <- c('T/A',    NA, 'T/T',    NA, 'T/A',    NA, 'T/T', 'T/A',
        'T/T', 'T/T', 'T/A', 'A/A', 'T/T', 'T/A', 'T/A', 'T/T',
        NA, 'T/A', 'T/A',   NA) 

g2 <-  c('C/A', 'C/A', 'C/C', 'C/A', 'C/C', 'C/A', 'C/A', 'C/A',
         'C/A', 'C/C', 'C/A', 'A/A', 'C/A', 'A/A', 'C/A', 'C/C',
         'C/A', 'C/A', 'C/A', 'A/A') 

g3 <-  c('T/A', 'T/A', 'T/T', 'T/A', 'T/T', 'T/A', 'T/A', 'T/A',
         'T/A', 'T/T', 'T/A', 'T/T', 'T/A', 'T/A', 'T/A', 'T/T',
         'T/A', 'T/A', 'T/A', 'T/T') 



founder <- data.frame(g1, g2, g3, stringsAsFactors = FALSE)
founder <- founder[rep(seq_len(nrow(founder)), 10),]
X <- founder[,sample.int(ncol(founder), 100, rep=TRUE)]
dim(X)


generateGenotypeMatrix <- function (nucleotides, n.loci, n.geno, prob.NA) {  
  genotype.matrix <- matrix(NA_character_, nrow = n.geno, ncol = n.loci)
  for (j in seq_len(n.loci)) {
    nuc <- sample(x = nucleotides, size = 2L, replace = FALSE)
    p <- rbeta(n = 1L, shape1 = 4, shape2 = 4)
    a <- sample(x = nuc, size = n.geno, replace = TRUE, prob = c(p, 1.0-p))
    b <- sample(x = nuc, size = n.geno, replace = TRUE, prob = c(p, 1.0-p))
    ab <- paste(a, b, sep = "/")
    genotype.matrix[,j] <- ab
  }
  genotype.matrix[runif(n = n.loci * n.geno) < prob.NA] <- NA_character_
  return(genotype.matrix)
}


transformGenotypes <- function (genotype.matrix, sep = "/") {
  
  transGeno <- function(geno, sep) {
    
    splt.list <- strsplit(x = geno, split = sep, fixed = TRUE)
    freq <- table(unlist(splt.list), useNA = "always")
    stopifnot(length(na.omit(names(freq))) == 2L)
    freq[is.na(names(freq))] <- 2L * freq[is.na(names(freq))]
    stopifnot(sum(freq) == 2L * length(geno))
    major.allele <- names(which.max(freq[!is.na(names(freq))]))
    p.major <- freq[major.allele] / sum(freq[!is.na(names(freq))])
    allele.count <- sapply(splt.list, function (x) sum(x == major.allele))
       
    return(list(freq = freq,
           allele.count = allele.count,
           major.allele = major.allele,
           p.major = p.major))
  }
   
  M <- ncol(genotype.matrix)
  genotypes.list <- vector("list", M)
  for (j in seq_len(M)) {
    genotypes.list[[j]] <- transGeno(geno = genotype.matrix[, j, drop = TRUE], 
                                     sep = sep)
  }
  
  return(genotypes.list)
}





set.seed(123L)
nucleotides <- c("A", "C", "G", "T")
n.loci <- 3L
n.geno <- 6L
prob.NA <- 0.2

(genotype.matrix <- generateGenotypeMatrix(nucleotides, n.loci, n.geno, prob.NA))

library("genetics")
library("gdata")
new.geno <- data.frame(lapply(as.data.frame(genotype.matrix), as.genotype))
LD.mat <- LD(new.geno)[["R^2"]]
a <- lowerTriangle(t(LD.mat))

out <- unphasedLD(genotype.matrix, sep = "/")
b <- out[,"r2"]

cor(a, b)



tmp <- transformGenotypes(genotype.matrix = genotype.matrix, sep = "/")
t.geno <- tmp$new.geno
major.freq <- tmp$major.freq

u <- LDwithMLE(cMajMat = t.geno, majFreqs = major.freq)
head(u)









g.a <- xt
g.b <- yt
g.b[5] <- NA
pA <- mean(g.a, na.rm=T) / 2
pB <- mean(g.b, na.rm=T)/ 2


loglik <- function(pAB,...)
  {
    (2*n3x3[1,1]+n3x3[1,2]+n3x3[2,1])*log(pAB) +
      (2*n3x3[1,3]+n3x3[1,2]+n3x3[2,3])*log(pA-pAB) +
      (2*n3x3[3,1]+n3x3[2,1]+n3x3[3,2])*log(pB-pAB) +
      (2*n3x3[3,3]+n3x3[3,2]+n3x3[2,3])*log(1-pA-pB+pAB) + 
      n3x3[2,2]*log(pAB*(1-pA-pB+pAB) + (pA-pAB)*(pB-pAB))
  }
curve(loglik, from = 0, to = 2)

.optimum_pxy(n3x3 = LDtools:::get_counts(g.a, g.b),
             px = pA, py = pB,
             lower = 0,
             upper = 10.0)
Dmin = max(-pA * pB, -(1-pA) * (1-pB))
Dmax = min(pA * (1-pB), (1-pA) * pB)

.LD(x = g.a, y = g.b, px = pA, py = pB, any_na = FALSE, is_phased = FALSE)
LDtools::LD(x = g.a, y = g.b, px = pA, py = pB, any_na = FALSE, is_phased = FALSE)

LDGenotypes <- function (g.a, g.b, pA, pB) {
  
  pa <- 1.0 - pA
  pb <- 1.0 - pB
  
  Dmin <- max(-pA * pB, -pa * pb)
  pmin <- pA * pB + Dmin
  
  Dmax <- min(pA * pb, pB * pa)
  pmax <- pA * pB + Dmax
    
  counts <- table(g.a, g.b)
  
  n3x3 <- matrix(0, nrow=3, ncol=3)
  colnames(n3x3) <- rownames(n3x3) <- 0:2
  
  # ensure the matrix is 3x3, with highest frequency values in upper left
  for(i in rownames(counts))
    for(j in colnames(counts))
      n3x3[3-as.numeric(i),3-as.numeric(j)] <- counts[i,j]
  
  
  loglik <- function(pAB,...)
  {
    (2*n3x3[1,1]+n3x3[1,2]+n3x3[2,1])*log(pAB) +
      (2*n3x3[1,3]+n3x3[1,2]+n3x3[2,3])*log(pA-pAB) +
      (2*n3x3[3,1]+n3x3[2,1]+n3x3[3,2])*log(pB-pAB) +
      (2*n3x3[3,3]+n3x3[3,2]+n3x3[2,3])*log(1-pA-pB+pAB) + 
      n3x3[2,2]*log(pAB*(1-pA-pB+pAB) + (pA-pAB)*(pB-pAB))
  }
  
  solution <- optimize(
    loglik,
    lower=pmin+.Machine$double.eps,
    upper=pmax-.Machine$double.eps,
    maximum=TRUE
  )
  pAB <- solution$maximum
  
  estD <- pAB - pA*pB
  if (estD>0)  
    estDp <- estD / Dmax
  else
    estDp <- estD / Dmin
  
  n <-  sum(n3x3)
  
  corr <- estD / sqrt( pA * pB * pa * pb )
  
  dchi <- (2*n*estD^2)/(pA * pa * pB* pb)
  dpval <- 1 - pchisq(dchi,1)
  
  retval <- list(
    call=match.call(),
    "D"=estD,
    "D'"=estDp,
    "r" = corr,
    "R^2" = corr^2,
    "n"=n,
    "X^2"=dchi,
    "P-value"=dpval
  )
  
  return(retval)
}

LDGenotypes(g.a, g.b)


  
cppFunction(code = '
  int test(const Rcpp::List & t_genotypes, 
                  const int i) {

  int n_geno = t_genotypes.size();





  Rcpp::List geno = t_genotypes[i];

  return n_geno;
  }
')


test(t.genotypes, 0L)




geno <- genotype.matrix[,1L]



cppFunction(code = '
  Rcpp::String test(const Rcpp::CharacterVector locus) {
            
  //  int n_geno = locus.size();
  //  Rcpp:CharacterMatrix geno_mat(n_geno, 2);

  //  for (int i = 0; i < n_geno; i++) {
  //    if (Rcpp::is_na(locus[i])) {
  //      geno_mat(i,0) = 
    //  }
   // }

  //  std::string gt = Rcpp::as<std::string>(locus[27]);
  //  int str_len = gt.length();
  Rcpp::String gt = locus[0];

  Rcpp:Rcout << gt.size() << std::endl;
    return gt;
  }
')

test(geno)






cppFunction(code = '
  double sacha(Rcpp::List L, int i) {
  double sum = 0;
     // for (int i=0; i<L.size(); i++) {
        Rcpp::NumericMatrix M = L[i];
        double topleft = M(0,0);
        Rcpp::Rcout << "Element is " << topleft << std::endl;
    //  }
      return 0;    
    }          
')

L <- list(matrix(rnorm(9),3), matrix(1:9,3), matrix(sqrt(1:4),2))
sacha(L, 2)












            IntegerMatrix getCounts(const NumericVector & x,
            const NumericVector & y){
            int N = x.size();
            IntegerMatrix counts(3, 3);
            std::fill(counts.begin(), counts.end(), 0);
            
            for(int i = 0; i != N; i++){
            if(!IntegerVector::is_na(x(i)) & !IntegerVector::is_na(y(i))){
            counts(2-x(i), 2-y(i)) += 1;
            }
            }
            return(counts);
            }
            ')



 
  





microbenchmark(z={
  
  
  # A(1), C(2), G(3), T(4)
  Desoxys <- c("A", "C", "G", "T")
  coding <- seq_len(length(Desoxys))
  genoTable <- outer(Desoxys, Desoxys, FUN=paste, sep="/")
  # codingTable <- as.integer(outer(coding, coding, FUN=paste, sep=""))
  codingTable <- outer(coding, coding, FUN=function(x,y){
    as.integer(paste0(x,y))
  })
  
  Z <- matrix(NA_integer_, nrow = nrow(X), ncol = ncol(X))
  for(j in seq_len(ncol(X))){
    for(i in seq_len(nrow(X))){
      Z[i,j] <- codingTable[match(x = X[i,j], table = genoTable)]
    }
  }
  
  
  
  cMajMat <- sapply(seq_len(ncol(Z)), FUN = function(j){
    u <- sapply(X = Z[,j,drop=TRUE], FUN = function(x){
      if(is.na(x)){
        return(c(NA_integer_, NA_integer_))
      } else {
        return(c(x%/%10, x%%10))
      }
    })
    ## get major allele plus frequncy
    ## get allele count of major allele including NA
    apply(u == as.integer(names(which.max(table(u)))), 2, sum)
  })
  
  
  majFreqs <- apply(cMajMat, 2, function(x){
    sum(x, na.rm = TRUE) / (2*sum(!is.na(x)))
  })
  
},times=1)



getCounts <- function(x, y){
  N <- length(x)
  counts <- matrix(0, nrow = 3, ncol = 3, dimnames = list(2:0,2:0))
  for(i in seq_len(N)){
    a <- x[i]
    b <- y[i]
    if(all(!is.na(c(a,b)))){
      counts[3-a, 3-b] <- counts[3-a, 3-b] + 1
    }
  }
  return(counts)
}


cppFunction(code = '
            
            IntegerMatrix getCounts(const NumericVector & x,
            const NumericVector & y){
            int N = x.size();
            IntegerMatrix counts(3, 3);
            std::fill(counts.begin(), counts.end(), 0);
            
            for(int i = 0; i != N; i++){
            if(!IntegerVector::is_na(x(i)) & !IntegerVector::is_na(y(i))){
            counts(2-x(i), 2-y(i)) += 1;
            }
            }
            return(counts);
            }
            ')


getCounts(x,y)





rotate = function(mat) t(mat[nrow(mat):1,,drop=FALSE])
n3x3 <- rotate(rotate(counts))


microbenchmark({
  for(i in 1:(ncol(X)-1)){
    for(j in (i + 1):ncol(X)){
      LD.genotype(genotype(X[,i]), genotype(X[,j]))
    }
    cat(sprintf("%5d\n",i))
  }
},times=1)


microbenchmark({
  LDwithMLE(cMajMat, majFreqs)
},times=1)





LD.genotype(genotype(X[,1]), genotype(X[,2]))

LD.genotype(genotype(g1), genotype(g2))



# Compute LD on a single pair

LD(g1,g2)

# Compute LD table for all 3 genotypes

data <- makeGenotypes(data.frame(g1,g2,g3))
LD(data)

fix(genotype)


LD(g1,g2)
LD.genotype(g2,g3)






LD.genotype <- function(g1,g2,...)
{
  prop.A <- summary(g1)$allele.freq[,2]
  prop.B <- summary(g2)$allele.freq[,2]
  ## obtain major allele by comparing frequencies
  major.A <- names(prop.A)[which.max(prop.A)]
  major.B <- names(prop.B)[which.max(prop.B)]
  ## extract allele frequncies of major allele
  pA <- max(prop.A, na.rm=TRUE)
  pB <- max(prop.B, na.rm=TRUE)
  ## compute minor allele frequencies
  pa <- 1-pA
  pb <- 1-pB
  ## compute stistics for calculating Dprime
  Dmin <- max(-pA*pB, -pa*pb)
  pmin <- pA*pB + Dmin;
  
  Dmax <- min(pA*pb, pB*pa);
  pmax <- pA*pB + Dmax;
  
  counts <- table(
    allele.count(g1, major.A),
    allele.count(g2, major.B)
  )
  
  n3x3 <- matrix(0, nrow=3, ncol=3)
  colnames(n3x3) <- rownames(n3x3) <- 0:2
  
  # ensure the matrix is 3x3, with highest frequency values in upper left
  for(i in rownames(counts))
    for(j in colnames(counts))
      n3x3[3-as.numeric(i),3-as.numeric(j)] <- counts[i,j]
  
  
  loglik <- function(pAB,...)
  {
    (2*n3x3[1,1]+n3x3[1,2]+n3x3[2,1])*log(pAB) +
      (2*n3x3[1,3]+n3x3[1,2]+n3x3[2,3])*log(pA-pAB) +
      (2*n3x3[3,1]+n3x3[2,1]+n3x3[3,2])*log(pB-pAB) +
      (2*n3x3[3,3]+n3x3[3,2]+n3x3[2,3])*log(1-pA-pB+pAB) + 
      n3x3[2,2]*log(pAB*(1-pA-pB+pAB) + (pA-pAB)*(pB-pAB))
  }
  
  # SAS code uses:
  #
  #s <- seq(pmin+0.0001,pmax-0.0001,by=0.0001)
  #lldmx <- loglik(s)
  #maxi <- which.max(lldmx)
  #pAB <- s[maxi]
  
  # but this should be faster:
  solution <- optimize(
    loglik,
    lower=pmin+.Machine$double.eps,
    upper=pmax-.Machine$double.eps,
    maximum=TRUE
  )
  pAB <- solution$maximum
  
  estD <- pAB - pA*pB
  if (estD>0)  
    estDp <- estD / Dmax
  else
    estDp <- estD / Dmin
  
  n <-  sum(n3x3)
  
  corr <- estD / sqrt( pA * pB * pa * pb )
  
  dchi <- (2*n*estD^2)/(pA * pa * pB* pb)
  dpval <- 1 - pchisq(dchi,1)
  
  
  evalCpp(pchisq(3.0, 1, 0, 1, 0))
  pchisq(3.0, 1, 0, 1, 0)
  
  
  "runit_pchisq" = list(
    signature( x = "numeric" ),
    '
    NumericVector xx(x) ;
    return List::create(_["lowerNoLog"] = pchisq( xx, 6.0 ),
    _["lowerLog"] = pchisq( xx, 6.0, true, true ),
    _["upperNoLog"] = pchisq( xx, 6.0, false ),
    _["upperLog"] = pchisq( xx, 6.0, false, true ));
    ')
  
  
  cppFunction(code = '
              double val(double xx){
              return pchisq(NumericVector::create(xx), 1.0, true, false )[0];
              }    
              ')
  
  
  val(3)
  
  evalCpp(pnorm(0.5, 0, 1, 1, 0))
  pnorm(q = 0.5, 0, 1, lower.tail = 1, log.p = 0)
  
  retval <- list(
    call=match.call(),
    "D"=estD,
    "D'"=estDp,
    "r" = corr,
    "R^2" = corr^2,
    "n"=n,
    "X^2"=dchi,
    "P-value"=dpval
  )
  
  class(retval) <- "LD"
  
  retval
  
}
