

L <- spam::read.MM('~/git_repos/RandPedigreePCA/datasets/pedLInv.mtx')
# K <- spam::read.MM('~/temp/pedLInv.mtx') # never finished
rppca(L, method = "randSVD")
rppca(L)
rppca(L, method=123)

ll <- Matrix::readMM('~/git_repos/RandPedigreePCA/datasets/pedLInv.mtx')
kk <- Matrix::readMM('~/temp/pedLInv.mtx')
typeof(kk)
class(kk)
kk
mode(kk)

ll
L
class(L)



pc <- rppca(ll)
t0 <- Sys.time()
svdL <- randSVD(L, 10, 3, 15)
Sys.time()- t0

t0 <- Sys.time()
svdll <- randSVD(ll, 10, 3, 15)
Sys.time()- t0



t0 <- Sys.time()
kkk <- as(kk, "dgCMatrix")
kkkk <-  as.spam.dgCMatrix(kkk)
svdkk <- randSVD(kkkk, 10, 3, 15)
Sys.time()- t0

plot(pc$T[,1:2])
MatrixClass(ll)
Matrix::MatrixClass(kk)
Matrix::MatrixClass(kkk)
Matrix::MatrixClass(ll)
Matrix::MatrixClass(L)
getClassdef(L)
getClassDef(L)




Linv <- importLinv("~/temp/pedLInv.mtx")

# smaller simulation
linv <- as.spam.dgCMatrix(as(readMM("~/git_repos/RandPedigreePCA/datasets/pedLInv.mtx"), "CsparseMatrix"))
pc <- rppca(linv)
plot(pc$scores[,1:2])
# try  naive approach
ainv <- readMM("~/git_repos/RandPedigreePCA/datasets/pedAInv.mtx")
ainv[1:10, 1:10]
aa <- solve(ainv)
llt <- t(linv) %*% linv
llti <- solve(llt)

# aa imported and llti are essentially identical:
all(zapsmall(llti - aa , digits = 10)==0)
all(zapsmall(llti - aa , digits = 13)==0)

hist(as.vector(aa))

# image(llti - aa) # takes a little time

hist(as.vector(llti - aa))
pcaa <- prcomp(aa, scale. = T)
pcaaNoscale <- prcomp(aa, scale. = F)
pcllti <- prcomp(llti, scale. = T)
pcainvPrincomp <- princomp(aa)

plot(pcaa$x[,1:2]) # somehow warped, 1 and 2 swapped
plot(pcllti$x[,1:2])  # somehow warped, 1 and 2 swapped
plot(pcainvPrincomp$scores[,1:2])  # not warped, 1 and 2 swapped

plot(pcaaNoscale$x[,1:2]) # not warped, 1 and 2 swapped
cor(pcaaNoscale$x[,1], pcaa$x[,1])
cor(pcllti$x[,1], pcaa$x[,1])
cor(pcaaNoscale$x[,1], pcainvPrincomp$scores[,1])
cor(pcaaNoscale$x[,1], pc$scores[,1])
cor(pcaaNoscale$x[,1], pc$scores[,2])
ccors <- cor(pcaaNoscale$x[,1:10], pc$scores)
ComplexHeatmap::Heatmap(as.matrix(ccors))
cor(pcainv$x[,1], pc$scores[,1])
cor(pcainv$x[,1], pc$scores[,2])
cor(pcainv$x[,2], pc$scores[,1])


pc <- rppca(Linv)
pc$T
smoothScatter(pc$T[,1:2])


svdObj <- randSVD(Linv, 10, 3, 15)
u <- svdObj$u
d <- svdObj$d
dim(u)
dim(d)
length(d)
TT <- u %*% diag(d)
smoothScatter(TT[,1:2])
?rppca
?randSVD
