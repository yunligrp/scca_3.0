scca_compute <-
function(X, Y, penalty="HL", olamxy, CR, nc=1, tuning="CV.alt",K=5, seed=NULL, center=TRUE, scale=FALSE){

    if (penalty=="SOFT" && min(lamx)>=1) {stop("Range of lamx for SOFT should be (0,1)")}
    if (penalty=="SOFT" && min(lamy)>=1) {stop("Range of lamy for SOFT should be (0,1)")}
 
    X<-scale(X,center=center,scale=scale);Y<-scale(Y,center=center,scale=scale)
 
    if (is.null(seed)) {seed<-Sys.time()}             
    if (nrow(X)!=nrow(Y)) {stop("X and Y should have same number of rows")}     

    svd.X<-svd(X) 
    
    n<-dim(X)[1];p<-dim(X)[2]; q<-dim(Y)[2]
    
    U<-matrix(0,n,nc)    
    V<-matrix(0,n,nc)    
    A<-matrix(0,p,nc)    
    B<-matrix(0,q,nc)    
    # CR<-rep(0,nc)      
    
    X.new<-X; Y.new<-Y
    
    for (cnt in 1:nc){
        olamx = olamxy[cnt, 1]
        olamy = olamxy[cnt, 2]
        CR1 = CR[cnt]
        if (penalty!="SOFT") {sscca.rst<-NIPALS.sparse(X.new,Y.new,olamx, olamy, penalty=penalty)}
        if (penalty=="SOFT") {sscca.rst<-NIPALS.soft(X.new,Y.new,olamx, olamy)}
              
        X.new<-X.new- sscca.rst$u1 %*%t(sscca.rst$a1)
        Y.new<-Y.new- sscca.rst$v1 %*%t(sscca.rst$b1)

        ytx.new<-t(Y.new)%*%X.new
        cat("Computing Component number ",cnt,"\n")
        
        U[,cnt]<-sscca.rst$u1; V[,cnt]<-sscca.rst$v1
        A[,cnt]<-sscca.rst$a1; B[,cnt]<-sscca.rst$b1
        CR[cnt]<-CR1
    }

return(list(U=U, V=V, A=A, B=B))
}