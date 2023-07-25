activity_plot = function(point.size, p.val, p.val.th = 0.05, shape = "s", size = c(1,2), alpha.min = .2,
                         title=NA, pw.names=NA, h.line=FALSE, v.line=FALSE){

  Max = 20
  Min = -log(.05)
  PSn = (point.size - Min) / (Max - Min)
  PSn[PSn<0]=0
  
  col_es = c() 
  col_p = c()
  size_p = c()
  N = nrow(PSn)
  M = ncol(PSn)
  for (i in seq(N,1,-1)) 
    for (j in 1:M){
      
      alpha = ifelse(p.val[i,j] < p.val.th , .8, alpha.min)
      
      # col_es = c(col_es, rgb(ESn[i,j],0,1-ESn[i,j], alpha = alpha))
      col_es = c(col_es, rgb(1,0,0, alpha = alpha))
      col_p = "red" #c(col_p, ifelse(p.val[i,j] > 10, rgb(.4,.4,.4,alpha=alpha),rgb(.8,.8,.8,alpha=alpha)))
      size_p = c(size_p, PSn[i,j]*(size[2]-size[1]) + size[1])
    }
  
  
  par(mar=c(1, 20 ,1 ,3))
  plot (rep(1:M,N),matrix(t(matrix(rep(1:N,M),N,M)),M*N,1), main=title,axes=T,xlab = "",ylab = "",xaxt="n",yaxt="n",
        cex = size_p, pch = ifelse(shape == "s", 22, 21), xlim=c(.5,M+ .5), ylim=c(.5,N+.5), lwd = 0, col = NA , bg = col_es)
  axis(2, at=length(pw.names):1, labels=pw.names, 
       col.axis="black", las=2,lwd = .2,lty=1,cex.axis=.5)
  if (h.line)
    abline(h = length(pw.names):1, lty = 1, lwd = .1, col=alpha("black",.2))
  if (v.line)
    abline(v = 1:4, lty = 1, lwd = .1, col=alpha("black",.2))
  
  
}

# activity_plot(matrix(runif(18,0,1),6,3), matrix(runif(18,.01,.1),6,3), alpha.min = 0, shape="c", size=2)
