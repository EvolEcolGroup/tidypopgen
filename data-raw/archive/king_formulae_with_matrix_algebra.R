#X <- show_genotypes(test_gt)
# code ignoring missing data
X0 <- X==0
X0[is.na(X0)]<-0
X1 <-X==1
X1[is.na(X1)]<-0
X2 <-X==2
X2[is.na(X2)]<-0
king_num <- (X1 %*% t(X1) - 2* ((X0) %*% t(X2) + (X2) %*% t(X0)) )

N_Aa_i <- matrix(rep(rowSums(X1), nrow(X)), nrow = nrow(X), byrow = T)
N_Aa_j <- matrix(rep(rowSums(X1), nrow(X)), nrow = nrow(X), byrow = F)
king.r.ignore.na <- king_num/(2* pmin(N_Aa_i,N_Aa_j))+0.5-0.25*(N_Aa_i+N_Aa_j)/pmin(N_Aa_i,N_Aa_j)

### code accounting for missing data




################################################################################
test_gt <- tidypopgen::gen_tibble("./data-raw/datasets/families.bed")
X_mat <- show_genotypes(test_gt)

X_mat0 <- X_mat==0
X_mat0[is.na(X_mat0)]<-0
X_mat1 <-X_mat==1
X_mat1[is.na(X_mat1)]<-0
X_mat2 <-X_mat==2
X_mat2[is.na(X_mat2)]<-0
king_num <- (X_mat1 %*% t(X_mat1) - 2* ((X_mat0) %*% t(X_mat2) + (X_mat2) %*% t(X_mat0)) )
X_mat_valid <- !is.na(X_mat)
N_mat_Aa_i <- X_mat1 %*% t(X_mat_valid)
N_mat_Aa_j <- t(N_mat_Aa_i)
king.r <- king_num/(2* pmin(N_mat_Aa_i,N_mat_Aa_j))+0.5-0.25*(N_mat_Aa_i+N_mat_Aa_j)/pmin(N_mat_Aa_i,N_mat_Aa_j)

X <- .gt_get_bigsnp(test_gt)
snp_king(X)

