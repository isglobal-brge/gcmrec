
*       * Functions to estimate General Class of Models for 
*       *                     Recurrent Event Data (Peña and Hollander, 2003)
*       *
*       *  Author: Juan R González
*       *  Date: July 2003
*       *  Note: Based on R functions written by Peña and Slate, April/July 2003
*       *  Parameters estimation procedures (Newton-Raphson with modifications 
*       *  to the Hessian) and Brent's method written by Juan R Gonzalez
*       *      
*       *    
*       *  Functions:
*       *       To be supplied
*       *
*       *
*       *
*       *  Input parameters:
*       *       To be supplied
*       *
*       *
*       *  Returned parameters:
*       *       To be supplied

*       * 
*       *  Work arrays:
*       *      arrays with OK are arrays for the subject i (not dinamic arrays)
*       *

*       *
*       *  Calls functions:  
*       *       To be supplied

*       *  Caution: 
*       *  the data must be ordered by id and event-censored
*       *      MAXIMUM NUMBER OF INDIVIDUAL RECURRENCES SET EQUAL TO 200
*       *      
*       /


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C       Jacknife procedures
C


c   Case without frailties  
       
       SUBROUTINE Jacknife(s,n,nvar,k,nk,tau,caltimes,gaptimes,
     *    censored,intercepts,slopes,lastperrep,perrepind,
     *    effagebegin,effage,cov,alphaSeed,betaSeed,Z,offset,rhoFunc,
     *    ns,tol,maxiter,estiEndJack)
       
       implicit none
       
       integer n,nvar,rhoFunc,searchBoth
       integer k(n),nk,i,j,ns(n),maxiter,method,kkiter
       double precision s,tau(n),caltimes(nk),gaptimes(nk),
     *    censored(n),offset(n),
     *    intercepts(nk),slopes(nk),lastperrep(nk),perrepind(nk),
     *    effagebegin(nk),effage(nk),cov(nvar,nk),
     *    alphaSeed,betaSeed(nvar),Z(n),tol,loglik,estiEnd(nvar+1),
     *    info(nvar+1,nvar+1) 
       
       integer kOK(n-1),nkOK,start(n),stop(n)
       double precision tauOK(n-1),caltimesOK(nk),gaptimesOK(nk),
     *  censoredOK(n-1),interceptsOK(nk),slopesOK(nk),lastperrepOK(nk),
     *  perrepindOK(nk),effagebeginOK(nk),effageOK(nk),covOK(nvar,nk),
     *  estiEndJack(n,nvar+1),offsetOK(n-1)
       
       
       stop(1)=k(1)
       start(1)=1
       do i=2,n
         stop(i)=stop(i-1)+k(i)
         start(i)=start(i-1)+k(i-1)
       end do 
       
       
       do i=1,n  
       
c  Leave-one-out
         call extract(n,k,i,kOK)
         nkOk=nk-k(i)
         call dextract(n,tau,i,tauOK)
         call dextract(n,offset,i,offsetOK)
         call extract2 (nk,caltimes,start(i),stop(i),caltimesOK)
         call extract2 (nk,gaptimes,start(i),stop(i),gaptimesOK)
         call dextract(n,censored,i,censoredOK)
         call extract2 (nk,intercepts,start(i),stop(i),interceptsOK)
         call extract2 (nk,slopes,start(i),stop(i),slopesOK)
         
         call extract2 (nk,lastperrep,start(i),stop(i),lastperrepOK)
         call extract2 (nk,perrepind,start(i),stop(i),perrepindOK)
         call extract2 (nk,effagebegin,start(i),stop(i),effagebeginOK)
         call extract2 (nk,effage,start(i),stop(i),effageOK)
         call extract3 (nvar,nk,cov,start(i),stop(i),covOK)
       
c Estimation 
       
         call newtraph(s,n-1,nvar,kOK,nkOK,tauOK,caltimesOK,gaptimesOK,
     *  censoredOK,interceptsOK,slopesOK,lastperrepOK,perrepindOK,
     *  effagebeginOK,effageOK,covOK,alphaSeed,betaSeed,Z,offsetOK,
     *  rhoFunc,ns,tol,maxiter,loglik,estiEnd,info,searchBoth,
     *  kkiter)
                  
         do j=1,nvar+1
          estiEndJack(i,j)=estiEnd(j)
         end do
          
       end do
       
       return
       
       END SUBROUTINE



c   Case with frailties  
       
       SUBROUTINE Jacknife2(s,n,nvar,k,nk,tau,caltimes,gaptimes,
     *    censored,intercepts,slopes,lastperrep,perrepind,
     *    effagebegin,effage,ndiseff,diseff,cov,alphaSeed,betaSeed,
     *    xi,Z,offset,rhoFunc,tol,maxiter,maxXi,estiEndJack)
       
       implicit none
       
       integer n,nvar,rhoFunc,control,ndiseff,maxXi
       integer k(n),nk,i,j,ns(n),maxiter,method,kkiter
       double precision s,tau(n),caltimes(nk),gaptimes(nk),
     *    censored(n),offset(n),
     *    intercepts(nk),slopes(nk),lastperrep(nk),perrepind(nk),
     *    effagebegin(nk),effage(nk),cov(nvar,nk),
     *    alphaSeed,betaSeed(nvar),Z(n),tol,loglik,
     *    estiEnd(nvar+2+n),info(nvar+1,nvar+1),xi,loglikEnd,
     *    diseff(ndiseff) 
       
       integer kOK(n-1),nkOK,start(n),stop(n)
       double precision tauOK(n-1),caltimesOK(nk),gaptimesOK(nk),
     *  censoredOK(n-1),interceptsOK(nk),slopesOK(nk),lastperrepOK(nk),
     *  perrepindOK(nk),effagebeginOK(nk),effageOK(nk),covOK(nvar,nk),
     *  estiEndJack(n,nvar+2),offsetOK(n-1)
       
       
       stop(1)=k(1)
       start(1)=1
       do i=2,n
         stop(i)=stop(i-1)+k(i)
         start(i)=start(i-1)+k(i-1)
       end do 
       
       
       do i=1,n  
       
c  Leave-one-out
         call extract(n,k,i,kOK)
         nkOk=nk-k(i)
         call dextract(n,tau,i,tauOK)
         call dextract(n,offset,i,offsetOK)
         call extract2 (nk,caltimes,start(i),stop(i),caltimesOK)
         call extract2 (nk,gaptimes,start(i),stop(i),gaptimesOK)
         call dextract(n,censored,i,censoredOK)
         call extract2 (nk,intercepts,start(i),stop(i),interceptsOK)
         call extract2 (nk,slopes,start(i),stop(i),slopesOK)
         
         call extract2 (nk,lastperrep,start(i),stop(i),lastperrepOK)
         call extract2 (nk,perrepind,start(i),stop(i),perrepindOK)
         call extract2 (nk,effagebegin,start(i),stop(i),effagebeginOK)
         call extract2 (nk,effage,start(i),stop(i),effageOK)
         call extract3 (nvar,nk,cov,start(i),stop(i),covOK)
       
c Estimation 
       
         call EstimWithFrailty(s,n-1,nvar,kOK,nkOK,tauOK,caltimesOK,
     *  gaptimesOK,censoredOK,interceptsOK,slopesOK,lastperrepOK,
     *  perrepindOK,effagebeginOK,effageOK,ndiseff,diseff,covOK,
     *  alphaSeed,betaSeed,xi,Z,offset,rhoFunc,tol,maxiter,maxXi,
     *  estiEnd,control,loglikEnd)
                  
         do j=1,nvar+1
          estiEndJack(i,j)=estiEnd(j)
         end do
         estiEndJack(i,nvar+2)=estiEnd(nvar+2) 
          
       end do
       
       return
       
       END SUBROUTINE
       


       
        SUBROUTINE extract(n,x,pos,xOK)
          implicit none
          integer n,pos,i,j
          integer x(n),xOk(n-1)
             
          i=1
          do j=1,n
            if (j.ne.pos) then
             xOk(i)=x(j)
             i=i+1
            end if       
          end do
          return  
       END SUBROUTINE
       
       SUBROUTINE dextract(n,x,pos,xOK)
          implicit none
          integer n,pos,i,j
          double precision x(n),xOk(n-1)
             
          i=1
          do j=1,n
            if (j.ne.pos) then
             xOk(i)=x(j)
             i=i+1
            end if      
          end do
          return  
       END SUBROUTINE
     
       
       SUBROUTINE extract2(n,x,start,stop,xOK)
          implicit none
          integer n,start,stop,i,j
          double precision x(n),xOk(n-(stop-start+1))
             
          i=1
          do j=1,n
            if ((j.lt.start).or.(j.gt.stop)) then
             xOk(i)=x(j)
             i=i+1
            end if       
          end do
          return  
       END SUBROUTINE
       
       SUBROUTINE extract3(nn,n,x,start,stop,xOK)
          implicit none
          integer n,nn,start,stop,i,j,k
          double precision x(nn,n),xOk(nn,n-(stop-start+1))
             
          do k=1,nn
           i=1
           do j=1,n
            if ((j.lt.start).or.(j.gt.stop)) then
             xOk(k,i)=x(k,j)
             i=i+1
            end if       
           end do
          end do
          return  
       END SUBROUTINE


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       
       
       
       
       
       
       SUBROUTINE EstimWithFrailty(s,n,nvar,k,nk,tau,caltimes,gaptimes,
     *  censored,intercepts,slopes,lastperrep,perrepind,effagebegin,
     *  effage,ndiseff,diseff,cov,alphaSeed,betaSeed,xi,Z,offset,
     *  rhoFunc,tol,maxiter,maxXi,estimEnd,control,loglikEnd)

       implicit none
       
       integer n,nvar,rhoFunc,ndiseff,method,searchxi
       integer k(n),nk,i,j,kiter,search,maxiter,searchNR,kiterNR
       double precision s,tau(n),caltimes(nk),gaptimes(nk),censored(n),
     *     intercepts(nk),slopes(nk),lastperrep(nk),perrepind(nk),
     *     effagebegin(nk),effage(nk),cov(nvar,nk),
     *     alphaSeed,betaSeed(nvar),Z(n),xi,diseff(ndiseff)

       double precision lambdaOld(ndiseff),survOld(ndiseff),
     *     lambdaNew(ndiseff),survNew(ndiseff),ZOld(n),xiOld,alphaOld,
     *     betaOld(nvar),AA(n),distAll,tol,ZNew(n),loglik,
     *     estiNew(nvar+1),
     *     info(nvar+1,nvar+1),alphaNew,betaNew(nvar),xiNew,norm,
     *     GammaLogLikOptim,distAlpha,distBeta,distXi,distZ,distLamb,
     *     deltalambdafuncOld(ndiseff),deltalambdaNew(ndiseff),BB(n)

       double precision estimEnd(nvar+2+n),HUGE,xiratOld,xiratNew,offset(n)
       real*8 loglikEnd,GammaLogLik,loglikMarg
       integer control(2),maxXi 
       
       integer KK(n),ns(n)  
       external GammaLogLikOptim,norm

       HUGE=1.0d60
       
c       Step 0: Initialization
       do i=1,n
          ZOld(i)=Z(i) 
       end do
       xiOld=xi
       xiratOld=xiOld/(1+xiOld)
       alphaOld=alphaSeed
       do i=1,nvar
          betaOld(i)=betaSeed(i)
       end do

       

       distAll=10.d0
       kiter=0
       search=0

c       Initialize the hazard estimate (same without frailties)
       
       call EstLambSurv(s,n,nvar,k,nk,tau,caltimes,gaptimes,
     *  censored,intercepts,slopes,lastperrep,perrepind,effagebegin,
     *  effage,ndiseff,diseff,cov,alphaOld,betaOld,ZOld,offset,
     *  rhoFunc,lambdaOld,deltalambdafuncOld,survOld)
       
       
       call CompAK(s,n,nvar,k,nk,tau,caltimes,gaptimes,
     *  censored,intercepts,slopes,lastperrep,perrepind,
     *  effagebegin,effage,ndiseff,diseff,cov,alphaOld,
     *  betaOld,deltalambdafuncOld,ZOld,offset,rhoFunc,KK,AA,BB)


c       Start of EM Iteration
       
      do while ((distAll.gt.tol).and.(kiter.lt.maxiter))

          kiter=kiter+1

c       Step-1 (E-Step): Given alphahat,betahat,xihat,Lambhat:
c       obtain new Zhat's
          call FrailtyValues(n,xiOld,KK,AA,ZNew)

          
c       Step-2 (M-Step#1): Given alphahat,betahat,xihat,Lambhat:
c       obtain new Lambhat
          
         call EstLambSurv(s,n,nvar,k,nk,tau,caltimes,gaptimes,
     *    censored,intercepts,slopes,lastperrep,perrepind,effagebegin,
     *    effage,ndiseff,diseff,cov,alphaOld,betaOld,ZNew,offset,
     *    rhoFunc,lambdaNew,deltalambdaNew,survNew)

c       Step-3 (M-Step#2): Given Lambhat,xihat,Zhat:
c       obtain new alphahat,betahat

c       $Newton-Raphson$                         

                if (rhoFunc.eq.2) then

              call newtraph(s,n,nvar,k,nk,tau,caltimes,gaptimes,
     *         censored,intercepts,slopes,lastperrep,perrepind,
     *         effagebegin,effage,cov,alphaOld,betaOld,ZNew,offset,
     *         rhoFunc,ns,tol,maxiter,loglik,estiNew,info,searchNR,
     *         kiterNR)

                     alphaNew=estiNew(1)
              do i=1,nvar
                 betaNew(i)=estiNew(i+1)
              end do   
            
             else

           call nrBeta(s,n,nvar,k,nk,tau,caltimes,gaptimes,
     *      censored,intercepts,slopes,lastperrep,perrepind,
     *      effagebegin,effage,cov,alphaOld,betaOld,ZNew,offset,rhoFunc,
     *      ns,tol,maxiter,loglik,estiNew,info,searchNR,kiterNR)

                     alphaNew=1.0d0
              do i=1,nvar
                 betaNew(i)=estiNew(i)
              end do   

                end if
 
          
          
c       Step-4 (M-Step#3): Given Lambhat,Zhat,alphahat,betahat:
c       obtain new xihat
          call CompAK(s,n,nvar,k,nk,tau,caltimes,gaptimes,
     *      censored,intercepts,slopes,lastperrep,perrepind,
     *      effagebegin,effage,ndiseff,diseff,cov,alphaNew,
     *      betaNew,deltalambdaNew,ZNew,offset,rhoFunc,KK,AA,BB)
          

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCC    Changes here:
CCCC    MaxWrtXi now accecpts new args maxiter and searchxi.
CCCC    Purpose to limit looping in MaxWrtXi.  If searchxi is 1 upon
CCCC    return, then MaxWrtXi stopped because of iteration limit.
CCCC    This is printed directly to R process using the subroutine intpr.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
           if (maxXi.eq.1) then           
            call MaxWrtXi(xiOld,n,KK,AA,BB,tol,xiNew,loglikMarg,
     *            maxiter,searchxi)
          else 
            xiNew=GammaLogLikOptim(xiOld,n,KK,AA,BB)
          end if

c          if (searchxi.eq.0) then
c            call intpr("MaxWrtXi failed, searchxi is:", -1, searchxi, 1)
c          end if
          if (xiNew.eq.HUGE) then
             xiratNew=1
          else 
             xiratNew=xiNew/(1+xiNew)
          end if


          
c       Step-5: Checking for Convergence
          distAlpha=norm(1,alphaOld,alphaNew)
          distBeta=norm(nvar,betaOld,betaNew)
c       distance based on xi/(1+xi)
          distXi=norm(1,xiratOld,xiratNew)
         distLamb=norm(ndiseff,lambdaOld,lambdaNew)    
c       no comparing Z's
c       distAll=max(max(distBeta,distAlpha),distXi) 
        distAll=max(max(max(distLamb,distXi),
     *             distBeta),distAlpha) 
          


c       Update the vectors
          alphaOld=alphaNew
          do i=1,nvar
             betaOld(i)=betaNew(i)
          end do

          xiOld=xiNew
          xiratOld=xiratNew

          do i=1,n
             ZOld(i)=ZNew(i)
          end do          
          do i=1,ndiseff
             lambdaOld(i)=lambdaNew(i)
          end do 
        

       end do 
       
       if (distAll.le.tol) then
          search=1
       end if
       

       control(1)=search
       control(2)=kiter
       estimEnd(1)=alphaOld
       do i=2,nvar+1
          estimEnd(i)=betaOld(i-1)
       end do
       estimEnd(nvar+2)=xiOld
       do i=nvar+3,nvar+2+n
          estimEnd(i)=ZOld(i-(nvar+2))
       end do
       loglikEnd=loglik
       
       
       return  
       END SUBROUTINE 
       
       
       
       
       SUBROUTINE MaxWrtXi(xiOld,n,K,A,B,tol,xiNew,G0,maxiter,search)
       implicit none

       integer n,i
       double precision xiOld,A(n),B(n),tol,dist,etaOld,
     *    xiratOld,G0,G1,G2,HUGE,etaNew,xiNew,xiratNew,tol1
       integer K(n),KK(n), niter, maxiter,search
       parameter (HUGE=1.0d60) 

c       Maximizing the log-likelihood for Xi

       do i=1,n
          KK(i)=K(i)-1
       end do

C       tol1=0.001
       tol1 = tol
       niter = 0
       search = 0
       
       dist=100.d0
       etaOld=dlog(xiOld)
       xiratOld=xiOld/(1+xiOld)

       do while ((dist.gt.tol1).and.(niter.lt.maxiter))
          niter = niter + 1
          call LogLikXi(xiOld,n,KK,A,B,G0,G1,G2)        
          etaNew=etaOld-(G1/(G1+xiOld*G2))
          xiNew=dexp(etaNew)
          if (xiNew.gt.HUGE) then
             xiNew=HUGE
             search = 1
             return
          else
             xiratNew=xiNew/(1+xiNew)
             dist=abs(xiratOld-xiratNew)
             xiOld=xiNew
             etaOld=etaNew
             xiratOld=xiratNew
          end if 

       end do
       if (dist.le.tol1) then
          search=1
       end if

       return
       END SUBROUTINE

       

       SUBROUTINE LogLikXi(xi,n,K,A,B,G0,G1,G2)

       implicit none

c       Marginal Log likelihood for Xi        

       integer i,j,n
       double precision xi,A(n),B(n),G0,G1,G2,q1,q2,q3
       integer K(n)

       
       G0=0.d0
       G1=0.d0
       G2=0.d0
       
       do i=1,n
          q1=0.d0
          q2=0.d0
          q3=0.d0
          if (K(i).gt.0) then
             do j=1,K(i)
              q1=q1+dlog(xi+dble(j-1))
              q2=q2+(1/(xi+dble(j-1)))
              q3=q3+(1/((xi+dble(j-1))**2))
             end do
          end if
          G0=G0+q1+(xi*dlog(xi))-((xi+K(i))*dlog(xi+A(i)))+B(i)
          G1=G1+q2+(dlog(xi)+1.d0)-dlog(xi+A(i))
     *       -((xi+dble(K(i)))/(xi+A(i)))
           G2=G2-q3+(1/xi)-(1/(xi+A(i)))-
     *       ((A(i)-dble(K(i)))/((xi+A(i))**2))
       end do
       
       return 
       END SUBROUTINE
      


       REAL*8 FUNCTION LogLikXi2(xi,n,K,A,B)

       implicit none

c       Marginal Log likelihood for Xi        

       integer i,j,n
       double precision xi,A(n),B(n),G0,q1
       integer K(n)

       
       G0=0.d0
       
       do i=1,n
          q1=0.d0
          if (K(i).gt.0) then
             do j=1,K(i)
              q1=q1+dlog(xi+dble(j-1))
             end do
          end if
          G0=G0+q1+(xi*dlog(xi))-((xi+K(i))*dlog(xi+A(i)))+B(i)
       end do
       
       LogLikXi2=-G0

       return 
       END FUNCTION
       
       
       DOUBLE PRECISION FUNCTION GammaLogLikOptim(seedXi,n,K,A,B)

       implicit none
       
       integer n,i
       double precision seedXi,A(n),B(n),HUGE
       integer K(n),KK(n)
       
       real*8 xmin,LogLikXi2
       external LogLikXi2

       real*8 ax,bx,cx,fa,fb,fc,temp

       parameter (HUGE=1.d030)

       do i=1,n
          KK(i)=K(i)-1
       end do
       
       ax=0.d-10
       bx=700.d0
       
       call mnbrak(ax,bx,cx,fa,fb,fc,LogLikXi2,n,KK,A,B)

       call brent(ax,bx,cx,LogLikXi2,1.d-6,xmin,n,KK,A,B)
       

       if (xmin.gt.HUGE) then
          GammaLogLikOptim=HUGE
          return 
       else
          GammaLogLikOptim=xmin
          return
       end if
       
       
       END FUNCTION
       
       
       
       
       
       SUBROUTINE FrailtyValues(n,xi,K,A,Z)
       
       implicit none
       
c       Estimates of the frailties values
c       We need to adjust K to reflect that 0 is included in count  
       integer n,i
       double precision xi,A(n),Z(n),HUGE
       integer K(n),KK(n)
       parameter (HUGE=1.0d60)

       do i=1,n
          KK(i)=K(i)-1
          Z(i)=0.d0
       end do
       
       if (xi.eq.HUGE) then
          do i=1,n
             Z(i)=1.d0 
          end do
          
       else

          do i=1,n
             Z(i)=(xi+dble(KK(i)))/(xi+A(i))
          end do
          
       end if

       return          
       END SUBROUTINE




       SUBROUTINE CompAK(s,n,nvar,k,nk,tau,caltimes,gaptimes,
     *    censored,intercepts,slopes,lastperrep,perrepind,
     *    effagebegin,effage,ndiseff,diseff,cov,alpha,
     *    beta,deltalambdafunc,Z,offset,rhoFunc,KK,AA,BB)
       implicit none

       integer n,nvar,rhoFunc,ndiseff,saca
       integer k(n),nk,i,j,r,t,u
       double precision s,tau(n),caltimes(nk),gaptimes(nk),
     *    censored(n),offset(n),
     *    intercepts(nk),slopes(nk),lastperrep(nk),perrepind(nk),
     *    effagebegin(nk),effage(nk),cov(nvar,nk),alpha,
     *    beta(nvar),Z(n),xi,diseff(ndiseff),deltalambdafunc(ndiseff)
       
       double precision S0,GrS0A1,Gr2S0A1,GrS0Be(nvar),
     *    Gr2S0A1Be(nvar),
     .       Gr2S0Be(nvar,nvar),Ysubj(n)
       
       double precision AA(n),BB(n),rho,psi,
     *    covariate(nvar),effageOK(200)
       integer KK(n),pos


c       Computes Ki's and Ai's to obtain new Z-values
       
       do i=1,n
          AA(i)=0.d0
       end do

       call nsm(s,n,nk,caltimes,k,KK)

       do i=1,ndiseff
          
          call AtRisk(s,KK,diseff(i),n,nvar,nk,k,tau,caltimes,
     *     gaptimes,censored,intercepts,slopes,lastperrep,
     *     perrepind,effagebegin,effage,cov,alpha,beta,Z,offset,
     *     rhoFunc,Ysubj,S0,GrS0A1,GrS0Be,Gr2S0A1,Gr2S0A1Be,Gr2S0Be)
          
          do j=1,n
             AA(j)=AA(j)+(Ysubj(j)*deltalambdafunc(i))
          end do

        end do                                                                
       
       do i=1,n
          BB(i)=0.d0
       end do
       

       pos=1

       do i=1,n
          do r=1,k(i)
             do t=1,nvar
              covariate(t)=cov(t,pos)           
             end do
             effageOK(r)=effage(pos)
             pos=pos+1
          end do

          if (KK(i).gt.1) then
              


c ---- Changed Feb'09 after CRAN requirement

c             do j=2,KK(i)
c               BB(i)=BB(i)+dlog(rho(j-2,alpha,rhoFunc))+
c     *          dlog(psi(nvar,covariate,beta,offset(i)))
c               do u=1,ndiseff
c                if (effageOK(j).eq.diseff(u)) then
c                 BB(i)=BB(i)+dlog(deltalambdafunc(u))
c                 goto 3000
c                end if
c               end do 

c 3000       end do 
           
             j=2
             saca=0 
             do while ((j.le.KK(i)).and.(saca.eq.0))
               BB(i)=BB(i)+dlog(rho(j-2,alpha,rhoFunc))+
     *          dlog(psi(nvar,covariate,beta,offset(i)))
               u=1 
               do while ((u.le.ndiseff).and.(saca.eq.0))
                if (effageOK(j).eq.diseff(u)) then
                 BB(i)=BB(i)+dlog(deltalambdafunc(u))
                 saca=1
                end if
               u=u+1
               end do 
             j=j+1 
             end do 
 
        end if
          
       end do       

       return      
       END SUBROUTINE 




       SUBROUTINE newtraph(s,n,nvar,k,nk,tau,caltimes,gaptimes,
     *    censored,intercepts,slopes,lastperrep,perrepind,
     *    effagebegin,effage,cov,alphaSeed,betaSeed,Z,offset,rhoFunc,
     *    ns,tol,maxiter,loglik,estiEnd,info,searchBoth,
     *    kkiter)
       
       implicit none
       
       integer n,nvar,rhoFunc
       integer k(n),nk,i,j
       double precision s,tau(n),caltimes(nk),gaptimes(nk),
     *    censored(n),offset(n),
     *    intercepts(nk),slopes(nk),lastperrep(nk),perrepind(nk),
     *    effagebegin(nk),effage(nk),cov(nvar,nk),
     *    alphaSeed,betaSeed(nvar),Z(n),
     *    loglik,score(nvar+1),info(nvar+1,nvar+1),
     *    scorealpha,infoalpha,scorebeta(nvar),
     *    infoalphabeta(nvar),infobeta(nvar,nvar)
       
       
       double precision estiOld(nvar+1),estiNew(nvar+1),
     *    DecDirec(nvar+1),
     *    distance,tol,alphaOld,betaOld(nvar),alpha,beta(nvar),
     *    estiNewAlpha,estiNewBeta(nvar),tol1,estiEnd(nvar+1)
       
       integer ns(n),kkiter,maxiter,
     *    search,searchAlpha,searchBeta,searchBoth,
     *    kkitera, kkiterb, kkiterab


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C    call realpr("loglik", -1, loglik, 1)
C    call intpr("kkiter", -1, kkiter, 1)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

       tol1=100.0*tol

c       search:
c       121: problems with alpha step
c       122: problems with beta step
c       0: no convergence both 
c       1: OK
        
       
c       initialization 
       distance=1000.d0
       kkiter=1

       alphaOld=alphaSeed
       do i=1,nvar
          betaOld(i)=betaSeed(i)
       end do
       search=0   

       
c       Newton-Raphson procedure alpha and beta 
       do while ((distance.gt.tol1).and.(kkiter.lt.maxiter)) 
          
          kkiter=kkiter+1
c       Maximum Alpha 
          
          estiOld(1)=alphaOld
          do i=2,nvar+1
             estiOld(i)=betaOld(i-1)
          end do
          
          call newtraphAlpha(s,n,nvar,k,nk,tau,caltimes,gaptimes,
     *       censored,intercepts,slopes,lastperrep,perrepind,
     *       effagebegin,effage,cov,alphaOld,betaOld,Z,offset,rhoFunc,
     *       ns,tol1,maxiter,estiNewAlpha,searchAlpha,kkitera)
          
          alphaOld=estiNewAlpha
          
          if (searchAlpha.ne.1) then
             search=121  
          end if         

c       Maximum Beta
          call newtraphBeta(s,n,nvar,k,nk,tau,caltimes,gaptimes,
     *      censored,intercepts,slopes,lastperrep,perrepind,
     *      effagebegin,effage,cov,alphaOld,betaOld,Z,offset,rhoFunc,
     *      ns,tol1,maxiter,estiNewBeta,searchBeta,kkiterb)
          
          
          do i=1,nvar
             betaOld(i)=estiNewBeta(i)
          end do
          
          if (searchBeta.ne.1) then
             search=122  
          end if         


          estiNew(1)=alphaOld
          do i=2,nvar+1
             estiNew(i)=betaOld(i-1)
          end do


c       stopping rule    
          distance=(estiNew(1)-estiOld(1))**2
          do i=2,nvar+1
             distance=distance+(estiNew(i)-estiOld(i))**2
          end do
          distance=distance**0.5

       end do 

       if (distance.le.tol1) then
          search=1
       end if

       alpha=estiNew(1)
       do i=1,nvar
          beta(i)=estiNew(i+1)
       end do
       

c       Maximum both
       if (search.eq.1) then
          call newtraphBoth(s,n,nvar,k,nk,tau,caltimes,gaptimes,
     *       censored,intercepts,slopes,lastperrep,perrepind,
     *       effagebegin,effage,cov,alpha,beta,Z,offset,rhoFunc,
     *       ns,tol,maxiter,loglik,estiEnd,
     *       info,searchBoth,kkiterab)

       end if
       
       return
       END SUBROUTINE 


       
       SUBROUTINE newtraphBoth(s,n,nvar,k,nk,tau,caltimes,gaptimes,
     *   censored,intercepts,slopes,lastperrep,perrepind,
     *   effagebegin,effage,cov,alphaSeed,betaSeed,Z,offset,rhoFunc,
     *   ns,tol,maxiter,loglik,estiNew,info,search,kiter)
       
       implicit none
       
       integer n,nvar,rhoFunc
       integer k(n),nk,i,j
       double precision s,tau(n),caltimes(nk),
     *    gaptimes(nk),censored(n),offset(n),
     *   intercepts(nk),slopes(nk),lastperrep(nk),perrepind(nk),
     *   effagebegin(nk),effage(nk),cov(nvar,nk),
     *   alphaSeed,betaSeed(nvar),Z(n),
     *   loglik,score(nvar+1),info(nvar+1,nvar+1),
     *   scorealpha,infoalpha,scorebeta(nvar),
     *   infoalphabeta(nvar),infobeta(nvar,nvar)
       
       
       double precision estiOld(nvar+1),estiNew(nvar+1),
     *    DecDirec(nvar+1),
     *   distance,tol,alphaOld,betaOld(nvar),
     *   dinc,wa(nvar+1,nvar+1)
       
       integer ns(n),kiter,maxiter,search
       real det
       
       
c       initialization 
       distance=1000.d0
       kiter=1
       alphaOld=alphaSeed
       do i=1,nvar
          betaOld(i)=betaSeed(i)
       end do
       search=0   


       
c       Newton-Raphson procedure 
       do while ((distance.gt.tol).and.(kiter.lt.maxiter)) 
          
          kiter=kiter+1
          call scorefunc(s,n,nvar,k,nk,tau,caltimes,gaptimes,
     *        censored,intercepts,slopes,lastperrep,perrepind,
     *        effagebegin,effage,cov,alphaOld,betaOld,Z,offset,rhoFunc,
     *        ns,loglik,score,info,scorealpha,scorebeta,
     *        infoalpha,infobeta)
          
          
          estiOld(1)=alphaOld
          do i=2,nvar+1
             estiOld(i)=betaOld(i-1)
          end do
          

          call MATINV(nvar+1,info,det)  
          
          
c       computes the descent direction 
          do i=1,nvar+1
             DecDirec(i)=0.d0
          end do
          do i=1,nvar+1
             do j=1,nvar+1
              DecDirec(i)=DecDirec(i)+info(i,j)*score(j)
             end do
          end do   
          
          do i=1,nvar+1
             estiNew(i)=estiOld(i)+DecDirec(i)
          end do
          
c       stopping rule    
          distance=(estiNew(1)-estiOld(1))**2
          do i=2,nvar+1
             distance=distance+(estiNew(i)-estiOld(i))**2
          end do
          distance=distance**0.5

          alphaOld=estiNew(1)
          do i=1,nvar
             betaOld(i)=estiNew(i+1)
          end do
          
          
          if (distance.le.tol) then
             search=1
          end if

       end do 
       
       return
       END SUBROUTINE 

       
       
       SUBROUTINE newtraphAlpha(s,n,nvar,k,nk,tau,caltimes,gaptimes,
     *   censored,intercepts,slopes,lastperrep,perrepind,
     *   effagebegin,effage,cov,alphaSeed,
     *   betaSeed,Z,offset,rhoFunc,
     *   ns,tol,maxiter,estiNew,search,kiter)
       
       implicit none
       
       integer n,nvar,rhoFunc
       integer k(n),nk,i,j
       double precision s,tau(n),caltimes(nk),
     *    gaptimes(nk),censored(n),offset(n),
     *   intercepts(nk),slopes(nk),lastperrep(nk),perrepind(nk),
     *   effagebegin(nk),effage(nk),cov(nvar,nk),
     *   alphaSeed,betaSeed(nvar),Z(n),
     *   loglik,score(nvar+1),info(nvar+1,nvar+1),
     *   scorealpha,infoalpha,scorebeta(nvar),
     *   infoalphabeta(nvar),infobeta(nvar,nvar)
       
       
       double precision estiOld,estiNew,distance,tol,alphaOld
       
       integer ns(n),kiter,maxiter,search
       real det
       
       
       
c       initialization 
       distance=1000.d0
       kiter=1
       alphaOld=alphaSeed
       search=0   


       
c       Newton-Raphson procedure 
       do while ((distance.gt.tol).and.(kiter.lt.maxiter)) 
          
          kiter=kiter+1
          call scorefunc(s,n,nvar,k,nk,tau,caltimes,gaptimes,
     *        censored,intercepts,slopes,lastperrep,perrepind,
     *        effagebegin,effage,cov,alphaOld,betaSeed,Z,offset,rhoFunc,
     *        ns,loglik,score,info,scorealpha,scorebeta,
     *        infoalpha,infobeta)
          
          
          estiOld=alphaOld
          
          estiNew=estiOld+scorealpha/infoalpha
          
          
c       stopping rule    
          distance=abs(estiNew-estiOld)

          alphaOld=estiNew
          
          
          if (distance.le.tol) then
             search=1
          end if

       end do 
       
       return
       END SUBROUTINE 


       SUBROUTINE newtraphBeta(s,n,nvar,k,nk,tau,caltimes,gaptimes,
     *    censored,intercepts,slopes,lastperrep,perrepind,
     *    effagebegin,effage,cov,alphaSeed,
     *    betaSeed,Z,offset,rhoFunc,
     *    ns,tol,maxiter,estiNew,search,kiter)
       
       implicit none
       
       integer n,nvar,rhoFunc
       integer k(n),nk,i,j
       double precision s,tau(n),caltimes(nk),
     *    gaptimes(nk),censored(n),offset(n),
     *    intercepts(nk),slopes(nk),lastperrep(nk),perrepind(nk),
     *    effagebegin(nk),effage(nk),cov(nvar,nk),
     *    alphaSeed,betaSeed(nvar),Z(n),
     *    loglik,score(nvar+1),info(nvar+1,nvar+1),
     *    scorealpha,infoalpha,scorebeta(nvar),
     *    infoalphabeta(nvar),infobeta(nvar,nvar)
       
       
       double precision estiOld(nvar),estiNew(nvar),DecDirec(nvar),
     *    distance,tol,betaOld(nvar)
       
       integer ns(n),kiter,maxiter,search
       real det
       
       
c       initialization 
       distance=1000.d0
       kiter=1
       do i=1,nvar
          betaOld(i)=betaSeed(i)
       end do
       search=0   

        
       
c       Newton-Raphson procedure 
       do while ((distance.gt.tol).and.(kiter.lt.maxiter)) 
          
          kiter=kiter+1
          call scorefunc(s,n,nvar,k,nk,tau,caltimes,gaptimes,
     *      censored,intercepts,slopes,lastperrep,perrepind,
     *      effagebegin,effage,cov,alphaSeed,betaOld,Z,offset,rhoFunc,
     *      ns,loglik,score,info,scorealpha,scorebeta,
     *      infoalpha,infobeta)
          
          
          do i=1,nvar
             estiOld(i)=betaOld(i)
          end do
          
          
          call MATINV(nvar,infobeta,det)
          
c       computes the descent direction 
          do i=1,nvar
             DecDirec(i)=0.d0
          end do
          do i=1,nvar
             do j=1,nvar
              DecDirec(i)=DecDirec(i)+infobeta(i,j)*scorebeta(j)
             end do
          end do   
          

          do i=1,nvar
             estiNew(i)=estiOld(i)+DecDirec(i)
          end do
          
c       stopping rule    
          distance=(estiNew(1)-estiOld(1))**2
          do i=2,nvar
             distance=distance+(estiNew(i)-estiOld(i))**2
          end do
          distance=distance**0.5

          
          do i=1,nvar
             betaOld(i)=estiNew(i)
          end do
       end do 
       
       if (distance.le.tol) then
          search=1
       end if
       
       return
       END SUBROUTINE 



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C             subroutines for rho=1 and for no covariates (no implemented this last one)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


cccc ----    rho=1

       SUBROUTINE nrBeta(s,n,nvar,k,nk,tau,caltimes,gaptimes,
     *    censored,intercepts,slopes,lastperrep,perrepind,
     *    effagebegin,effage,cov,alphaSeed,
     *    betaSeed,Z,offset,rhoFunc,
     *    ns,tol,maxiter,loglik,estiNew,infobeta,search,kiter)

       
       implicit none
       
       integer n,nvar,rhoFunc
       integer k(n),nk,i,j
       double precision s,tau(n),caltimes(nk),
     *    gaptimes(nk),censored(n),offset(n),
     *    intercepts(nk),slopes(nk),lastperrep(nk),perrepind(nk),
     *    effagebegin(nk),effage(nk),cov(nvar,nk),
     *    alphaSeed,betaSeed(nvar),Z(n),
     *    loglik,score(nvar+1),info(nvar+1,nvar+1),
     *    scorealpha,infoalpha,scorebeta(nvar),
     *    infoalphabeta(nvar),infobeta(nvar,nvar)
       
       
       double precision estiOld(nvar),estiNew(nvar),DecDirec(nvar),
     *    distance,tol,betaOld(nvar)
       
       integer ns(n),kiter,maxiter,search
       real det
       
       
c       initialization 
       distance=1000.d0
       kiter=1
       do i=1,nvar
          betaOld(i)=betaSeed(i)
       end do
       search=0   

        
       
c       Newton-Raphson procedure 
       do while ((distance.gt.tol).and.(kiter.lt.maxiter)) 
          
          kiter=kiter+1
          call scorefunc(s,n,nvar,k,nk,tau,caltimes,gaptimes,
     *      censored,intercepts,slopes,lastperrep,perrepind,
     *      effagebegin,effage,cov,alphaSeed,betaOld,Z,offset,rhoFunc,
     *      ns,loglik,score,info,scorealpha,scorebeta,
     *      infoalpha,infobeta)
          
          
                 do i=1,nvar
             estiOld(i)=betaOld(i)
          end do
          
          
          call MATINV(nvar,infobeta,det)
          
c       computes the descent direction 
          do i=1,nvar
             DecDirec(i)=0.d0
          end do
          do i=1,nvar
             do j=1,nvar
              DecDirec(i)=DecDirec(i)+infobeta(i,j)*scorebeta(j)
             end do
          end do   
          

          do i=1,nvar
             estiNew(i)=estiOld(i)+DecDirec(i)
          end do
          
c       stopping rule    
          distance=(estiNew(1)-estiOld(1))**2
          do i=2,nvar
             distance=distance+(estiNew(i)-estiOld(i))**2
          end do
          distance=distance**0.5

          
          do i=1,nvar
             betaOld(i)=estiNew(i)
          end do
       end do 
       
       if (distance.le.tol) then
          search=1
       end if
       
       return
       END SUBROUTINE 



ccc --------  No covariates

       SUBROUTINE nrAlpha(s,n,nvar,k,nk,tau,caltimes,gaptimes,
     *   censored,intercepts,slopes,lastperrep,perrepind,
     *   effagebegin,effage,cov,alphaSeed,
     *   betaSeed,Z,offset,rhoFunc,
     *   ns,tol,maxiter,loglik,estiNew,infoalpha,search,kiter)
       
       implicit none
       
       integer n,nvar,rhoFunc
       integer k(n),nk,i,j
       double precision s,tau(n),caltimes(nk),
     *    gaptimes(nk),censored(n),offset(n),
     *   intercepts(nk),slopes(nk),lastperrep(nk),perrepind(nk),
     *   effagebegin(nk),effage(nk),cov(nvar,nk),
     *   alphaSeed,betaSeed(nvar),Z(n),
     *   loglik,score(nvar+1),info(nvar+1,nvar+1),
     *   scorealpha,infoalpha,scorebeta(nvar),
     *   infoalphabeta(nvar),infobeta(nvar,nvar)
       
       
       double precision estiOld,estiNew,distance,tol,alphaOld
       
       integer ns(n),kiter,maxiter,search
       real det
       
       
       
c       initialization 
       distance=1000.d0
       kiter=1
       alphaOld=alphaSeed
       search=0   


       
c       Newton-Raphson procedure 
       do while ((distance.gt.tol).and.(kiter.lt.maxiter)) 
          
          kiter=kiter+1
          call scorefunc(s,n,nvar,k,nk,tau,caltimes,gaptimes,
     *        censored,intercepts,slopes,lastperrep,perrepind,
     *        effagebegin,effage,cov,alphaOld,betaSeed,Z,offset,rhoFunc,
     *        ns,loglik,score,info,scorealpha,scorebeta,
     *        infoalpha,infobeta)
          
          
          estiOld=alphaOld
          
          estiNew=estiOld+scorealpha/infoalpha
          
          
c       stopping rule    
          distance=abs(estiNew-estiOld)

          alphaOld=estiNew
          
          
          if (distance.le.tol) then
             search=1
          end if

       end do 
       
       return
       END SUBROUTINE 





CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


        SUBROUTINE EstLambSurv(s,n,nvar,k,nk,tau,caltimes,gaptimes,
     .  censored,intercepts,slopes,lastperrep,perrepind,effagebegin,
     .  effage,ndiseff,diseff,cov,alpha,beta,Z,offset,rhoFunc,
     .  lambdafunc,deltalambdafunc,survfunc)
      
        implicit none

       integer n,nvar,rhoFunc
       integer k(n),nk,i,j,ndiseff,ns(n)
       double precision s,tau(n),caltimes(nk),
     *   gaptimes(nk),censored(n),offset(n),
     *    intercepts(nk),slopes(nk),lastperrep(nk),perrepind(nk),
     *    effagebegin(nk),effage(nk),cov(nvar,nk),
     *    alpha,beta(nvar),Z(n),ONE,ZERO,diseff(ndiseff),
     *    DeltaNst(ndiseff),S0st(ndiseff) 
      
       double precision GrS0A1,Gr2S0A1,GrS0Be(nvar),Gr2S0A1Be(nvar),
     .       Gr2S0Be(nvar,nvar),Ysubj(n)
       
       double precision lambdafunc(ndiseff),survfunc(ndiseff),dellamb,
     .          deltalambdafunc(ndiseff)       

        ONE=1.d0
       ZERO=0.d0


c     This function need that alpha and beta are estimated

c   Initialization  
       do i=1,ndiseff
        lambdafunc(i)=ZERO 
        survfunc(i)=ONE
c   Delta N(s,t)
        DeltaNst(i)=ZERO
c   S0(s,t)
         S0st(i)=ZERO
      end do
      
       do i=2,ndiseff
        
         do j=1,nk
           if ((caltimes(j).le.s).and.(effage(j).eq.diseff(i))) then
             DeltaNst(i)=DeltaNst(i)+1  
           end if
          end do

          call nsm(s,n,nk,caltimes,k,ns)
        
      
         call AtRisk(s,ns,diseff(i),n,nvar,nk,k,tau,caltimes,
     .     gaptimes,censored,intercepts,slopes,lastperrep,perrepind,
     .     effagebegin,effage,cov,alpha,beta,Z,offset,rhoFunc,
     .     Ysubj,S0st(i),GrS0A1,GrS0Be,Gr2S0A1,Gr2S0A1Be,Gr2S0Be)
        
       
         if (S0st(i).gt.ZERO) then
           dellamb=min(DeltaNst(i)/S0st(i),ONE)
           lambdafunc(i)=lambdafunc(i-1)+dellamb
           survfunc(i)=survfunc(i-1)*(1-dellamb)
          else
           lambdafunc(i)=lambdafunc(i-1)
           survfunc(i)=survfunc(i-1)
         end if           
      
       end do
      
       deltalambdafunc(1)=0.d0
       do i=2,ndiseff
        deltalambdafunc(i)=lambdafunc(i)-lambdafunc(i-1)
       end do
         
       return
       END SUBROUTINE 



       SUBROUTINE scorefunc(s,n,nvar,k,nk,tau,caltimes,gaptimes,
     .         censored,intercepts,slopes,lastperrep,perrepind,
     .         effagebegin,effage,cov,alpha,beta,Z,offset,rhoFunc,
     .         ns,loglik,score,info,scorealpha,scorebeta,
     .         infoalpha,infobeta)
       

        implicit none
       
       integer n,nvar,rhoFunc
        integer k(n),nk,i,j,pos,r,t
       double precision s,tau(n),caltimes(nk),
     .       gaptimes(nk),censored(n),
     .       intercepts(nk),slopes(nk),lastperrep(nk),perrepind(nk),
     .       effagebegin(nk),effage(nk),cov(nvar,nk),
     .       alpha,beta(nvar),Z(n),w
       
       double precision loglik,scorealpha,infoalpha,
     .       scorebeta(nvar),infoalphabeta(nvar),infobeta(nvar,nvar),
     .       effageOK(200),ZERO
      
       integer ns(n)

       double precision S0,GrS0A1,Gr2S0A1,GrS0Be(nvar),
     .     Gr2S0A1Be(nvar),offset(n),
     .     Gr2S0Be(nvar,nvar),Zsubj,Ysubj(n),GrYsubjA1,Gr2YsubjA1,
     .     GrYsubjBe(nvar),Gr2YsubjA1Be(nvar),Gr2YsubjBe(nvar,nvar),
     .     covariate(nvar),ealpha,ebeta(nvar),ebeta2(nvar,nvar),
     .     score(nvar+1),info(nvar+1,nvar+1)
      
       real det
               
        double precision rho,psi
     
       ZERO=0.d0 
      
c     initialization

        loglik=ZERO
       scorealpha=ZERO
       infoalpha=ZERO
       do i=1,nvar
        scorebeta(i)=ZERO
         infoalphabeta(i)=ZERO
        ebeta(i)=ZERO
       end do
      
       do i=1,nvar
        do j=1,nvar
         infobeta(i,j)=ZERO
          ebeta2(i,j)=ZERO
        end do
       end do

       ealpha=ZERO

       
      
       call nsm(s,n,nk,caltimes,k,ns)
       
      
       pos=1
       do i=1,n
  
         do r=1,200
           effageOK(r)=0.d0
         end do
         
          do r=1,k(i)
           effageOK(r)=effage(pos)
           do t=1,nvar
           covariate(t)=cov(t,pos)           
          end do
          pos=pos+1
         end do
      
        

          if (ns(i).gt.1) then
           do j=2,ns(i)
            w=effageOK(j)
              call AtRisk(s,ns,w,n,nvar,nk,k,tau,caltimes,gaptimes,
     .        censored,intercepts,slopes,lastperrep,perrepind,
     .        effagebegin,effage,cov,alpha,beta,Z,offset,rhoFunc,
     .        Ysubj,S0,GrS0A1,GrS0Be,Gr2S0A1,Gr2S0A1Be,Gr2S0Be)
     
              loglik=loglik+dlog(rho(j-2,alpha,rhoFunc))+
     .               dlog(psi(nvar,covariate,beta,offset(i)))
     .              -dlog(S0)
     

              scorealpha=scorealpha+((dble(j-2)/alpha)-(GrS0A1/S0))

             ealpha=GrS0A1/S0
             infoalpha=infoalpha+(dble(j-2)/alpha**2)+
     .                          (Gr2S0A1/S0)-ealpha**2

                  do t=1,nvar
                  scorebeta(t)=scorebeta(t)+
     *            covariate(t)-(GrS0Be(t)/S0) 
                  end do
   
                  do t=1,nvar
                   ebeta(t)=GrS0Be(t)/S0 
                  end do 
  
                  do t=1,nvar 
                   infoalphabeta(t)=infoalphabeta(t)+(Gr2S0A1Be(t)/S0)-
     .                            ealpha*ebeta(t) 
                 end do
             
                 do t=1,nvar 
                  do r=1,nvar
                    ebeta2(t,r)=ebeta(t)*ebeta(r)
                  end do
                 end do

                  do t=1,nvar
                   do r=1,nvar
                    infobeta(t,r)=infobeta(t,r)+(Gr2S0Be(t,r)/S0)-
     *                            ebeta2(t,r)
                   end do
                  end do

            end do   
         end if

       
          end do
           



       score(1)=scorealpha
       info(1,1)=infoalpha
       
        do i=2,nvar+1
         score(i)=scorebeta(i-1)
        info(1,i)=infoalphabeta(i-1)
        info(i,1)=infoalphabeta(i-1)
       end do
       

       do i=2,nvar+1
        do j=2,nvar+1 
         info(i,j)=infobeta(i-1,j-1)
        end do
        end do 
       
       return
      END SUBROUTINE 




       
       SUBROUTINE nsm(s,n,nk,caltimes,k,nm)
       implicit none

       double precision s
       integer n,nk,pos,i,j
       double precision caltimes(nk)
        integer k(n)
       integer nm(n),nism

       double precision caltimesOK(200)


       do i=1,200
        caltimesOK(i)=0.d0
        end do

       pos=1 
       do i=1,n
         do j=1,k(i)
         caltimesOK(j)=caltimes(pos)
          pos=pos+1
        end do
         nm(i)=nism(s,caltimesOK,k(i))
       end do 
      
       return 
       END SUBROUTINE
       


       INTEGER FUNCTION nism(s,caltimes,k)
        implicit none

       double precision s
       integer k,i
       double precision caltimes(k)

c     routine that compute n_i^\dagger(s-) 
       if (s.gt.caltimes(k)) then
        nism=k
       
c     modification made by Juan (no break)
        else 
        nism=0
         i=1
         do while (caltimes(i).lt.s) 
           nism=nism+1
           i=i+1
         end do
       
       end if
       
       return
       END FUNCTION

      
       DOUBLE PRECISION FUNCTION norm(n,vec1,vec2)
        implicit none
      
       integer n,i
       double precision vec1(n),vec2(n),distance      

        distance=(vec1(1)-vec2(1))**2
        if (n.gt.1) then
         do i=2,n
           distance=distance+(vec1(i)-vec2(i))**2
         end do
         end if
         norm=distance**0.5
      
       return
       END FUNCTION 

       
       
       SUBROUTINE AtRisk(s,ns,w,n,nvar,nk,k,tau,caltimes,
     .    gaptimes,censored,intercepts,slopes,lastperrep,
     .    perrepind,effagebegin,effage,cov,alpha,beta,Z,offset,rhoFunc,
     .    Ysubj,S0,GrS0A1,GrS0Be,Gr2S0A1,Gr2S0A1Be,Gr2S0Be)

       implicit none
       
       integer n,nvar,nk,kk,i,j,t,pos,rhoFunc
       integer k(n),ns(n),nsi
       double precision w,s,tau(n),caltimes(nk),gaptimes(nk),
     .     censored(n),intercepts(nk),slopes(nk),lastperrep(nk),
     .     perrepind(nk),effagebegin(nk),effage(nk),cov(nvar,nk)

       double precision caltimesOK(200),gaptimesOK(200),
     .      censoredOK,interceptsOK(200),slopesOK(200),
     .      lastperrepOK(200),perrepindOK(200),effagebeginOK(200),
     .      effageOK(200),covOK(nvar,200) 
       
       double precision S0,GrS0A1,Gr2S0A1,GrS0Be(nvar),
     .    Gr2S0A1Be(nvar),offset(n),
     .    Gr2S0Be(nvar,nvar),alpha,beta(nvar),Z(n),Zsubj,ZERO,
     .    Ysubj(n),GrYsubjA1,Gr2YsubjA1,GrYsubjBe(nvar),
     .    Gr2YsubjA1Be(nvar),Gr2YsubjBe(nvar,nvar)      
       
       


c     Equations (25) (26) (27) S0 and gradients

c     initialization
       ZERO=0.d0
c     initialitation 
        S0=ZERO
        GrS0A1=ZERO
        Gr2S0A1=ZERO
        do i=1,nvar
         GrS0Be(i)=ZERO
         Gr2S0A1Be(i)=ZERO
        do j=1,nvar
         Gr2S0Be(i,j)=ZERO
        end do
        end do
        do i=1,n
          Ysubj(i)=ZERO
        end do

        pos=1
        do i=1,n
c     we select the data for each subject from the array
         censoredOK=censored(i)
        do j=1,k(i)
           caltimesOK(j)=caltimes(pos)
          gaptimesOK(j)=gaptimes(pos)
           interceptsOK(j)=intercepts(pos)
           slopesOK(j)=slopes(pos)
          slopesOK(j)=slopes(pos)
          lastperrepOK(j)=lastperrep(pos)
          perrepindOK(j)=perrepind(pos)
          effagebeginOK(j)=effagebegin(pos)
          effageOK(j)=effage(pos)
          do t=1,nvar
           covOK(t,j)=cov(t,pos)
          end do
          pos=pos+1
         end do
         
        nsi=ns(i)
         
         call AtRiskSubj(s,nsi,w,nvar,k(i),tau(i),caltimesOK,
     .     gaptimesOK,censoredOK,interceptsOK,slopesOK,lastperrepOK,
     .     perrepindOK,effagebeginOK,effageOK,covOK,alpha,beta,rhoFunc,
     .     offset(i),Ysubj(i),GrYsubjA1,GrYsubjBe,Gr2YsubjA1,
     .     Gr2YsubjA1Be,Gr2YsubjBe)
         
      

        Zsubj=Z(i)

         S0=S0+Zsubj*Ysubj(i)
       GrS0A1=GrS0A1+Zsubj*GrYsubjA1
       Gr2S0A1=Gr2S0A1+Zsubj*Gr2YsubjA1

       do j=1,nvar
         GrS0Be(j)=GrS0Be(j)+Zsubj*GrYsubjBe(j)
        Gr2S0A1Be(j)=Gr2S0A1Be(j)+Zsubj*Gr2YsubjA1Be(j)
        do t=1,nvar
          Gr2S0Be(j,t)=Gr2S0Be(j,t)+Zsubj*Gr2YsubjBe(j,t)
        end do 
       end do

            
       end do
       
       return
       END SUBROUTINE




      SUBROUTINE AtRiskSubj(s,nsi,w,nvar,k,tau,caltimes,
     .    gaptimes,censored,intercepts,slopes,lastperrep,
     .    perrepind,effagebegin,effage,cov,alpha,beta,rhoFunc,offset,
     .    Ysubj,GrYsubjA1,GrYsubjBe,Gr2YsubjA1,
     .    Gr2YsubjA1Be,Gr2YsubjBe)

      implicit none

      integer nsi,nvar,k,kk,i,j,t,rhoFunc
      double precision w,s,tau,caltimes(k),gaptimes(k),
     .   censored,intercepts(k),slopes(k),lastperrep(k),
     .   perrepind(k),effagebegin(k),effage(k),cov(nvar,k)
      
      double precision Q(nsi),R,GrQA1(nsi),GrQBe(nvar,nsi),GrRA1,
     . GrRBe(nvar),Gr2QA1(nsi),Gr2QA1Be(nvar,nsi),
     . Gr2QBe(nsi,nvar,nvar),Gr2RA1,Gr2RA1Be(nvar),Gr2RBe(nvar,nvar),
     . alpha,beta(nvar),ONE,ZERO
     
      double precision covariate(nvar),effageats,vect1(nsi),Ysubj,
     .    GrYsubjA1,Gr2YsubjA1,GrYsubjBe(nvar),Gr2YsubjA1Be(nvar),
     .    Gr2YsubjBe(nvar,nvar)
      
      double precision rho,psi,offset

      ONE=1.d0
      ZERO=0.d0

c     if nsi is NULL is not implemented  
c        take care with changes in nsi and array sizes !!!

      kk=nsi 


c     Initialization 
       R=ZERO
       GrRA1=ZERO
       Gr2RA1=ZERO
       do i =1,kk
        Q(i)=ZERO
c      Gradient of Q wrt alpha
        GrQA1(i)=ZERO
c      second derivative wrt alpha       
        Gr2QA1(i)=ZERO
        do j=1,nvar
c      Gradient of Q wrt beta
         GrQBe(j,i)=ZERO
         Gr2QA1Be(j,i)=ZERO
        end do
       end do
   
       do i=1,nvar
         GrRBe(i)=ZERO
         Gr2RA1Be(i)=ZERO
         do j=1,nvar
          Gr2RBe(j,i)=ZERO
         do t=1,kk
          Gr2QBe(t,i,j)=ZERO
         end do 
        end do
       end do

         
             
c     Equation (23) Q_ij
      if (kk.gt.1) then
        do j=2,kk
         if ((w.gt.effagebegin(j-1)).and.(w.le.effage(j))) then
          
              do i=1,nvar
               covariate(i)=cov(i,j-1)           
                end do

           Q(j)=(rho(j-2,alpha,rhoFunc)*psi(nvar,covariate,beta,offset))
     .                           /slopes(j-1)

c      change in the equations (26) (27) version June 6,2003
c      if rho is equal to $1$ or equal to $\alpha^j$
         if (rhoFunc.eq.1) then
c          $\rho(j;\alpha)=1$          
                 GrQA1(j)=Q(j)
          else
c          $\rho(j;\alpha)=\alpha^j$                     
                 GrQA1(j)=(dble(j-2)/alpha)*Q(j)
         end if  
          
                  do i=1,nvar
                  GrQBe(i,j)=covariate(i)*Q(j)
                  if (rhoFunc.eq.1) then
c             $\rho(j;\alpha)=1$                         
                   Gr2QA1Be(i,j)=covariate(i)*Q(j)
                  else
c             $\rho(j;\alpha)=\alpha^j$                     
                   Gr2QA1Be(i,j)=(dble(j-2)/alpha)*covariate(i)*Q(j)
                  end if
                  do t=1,nvar
             Gr2QBe(j,i,t)=Q(j)*covariate(i)*covariate(t) 
             end do
                end do
                if (rhoFunc.eq.1) then
c          $\rho(j;\alpha)=1$      
                Gr2QA1(j)=Q(j)
                else
c          $\rho(j;\alpha)=\alpha^j$                     
                Gr2QA1(j)=(dble((j-2)*(j-3))/alpha**2)*Q(j)
           end if 
          end if
         end do
      end if

       
       
      
c       Equation (24) R_i
       effageats=intercepts(kk)+slopes(kk)*min(s,tau)
     .     -caltimes(INT(lastperrep(kk)))


       if ((w.gt.effagebegin(kk)).and.(w.le.effageats)) then
          do i=1,nvar
             covariate(i)=cov(i,kk)           
          end do

          R=(rho(kk-1,alpha,rhoFunc)*
     *       psi(nvar,covariate,beta,offset))/slopes(kk)
          
          if (rhoFunc.eq.1) then
c       $\rho(j;\alpha)=1$  
             GrRA1=R
          else
c       $\rho(j;\alpha)=\alpha^j$          
             GrRA1=(dble(kk-1)/alpha)*R
          end if
          
          do i=1,nvar
             GrRBe(i)=covariate(i)*R
             if (rhoFunc.eq.1) then
c       $\rho(j;\alpha)=1$          
               Gr2RA1Be(i)=covariate(i)*R
             else
c       $\rho(j;\alpha)=\alpha^j$                  
               Gr2RA1Be(i)=(dble(kk-1)/alpha)*covariate(i)*R
             end if
             
             do j=1,nvar
                Gr2RBe(i,j)=(R*covariate(i))*covariate(j)
             end do
          end do

          if (rhoFunc.eq.1) then
c       $\rho(j;\alpha)=1$    
             Gr2RA1=R
          else 
c       $\rho(j;\alpha)=\alpha^j$             
             Gr2RA1=(dble((kk-1)*(kk-2))/alpha**2)*R
          end if   
       end if
       

c       OUTPUT 
c       Initialization 
       do i=1,kk
          vect1(i)=ONE
       end do
       Ysubj=ZERO
       GrYsubjA1=ZERO
       Gr2YsubjA1=ZERO
       do i=1,nvar
          GrYsubjBe(i)=ZERO
          Gr2YsubjA1Be(i)=ZERO
          do j=1,nvar
             Gr2YsubjBe(i,j)=ZERO
          end do
       end do

c       computes
       do i=1,kk
          Ysubj=Ysubj+Q(i)*vect1(i)
          GrYsubjA1=GrYsubjA1+GrQA1(i)*vect1(i)
          Gr2YsubjA1=Gr2YsubjA1+Gr2QA1(i)*vect1(i)
       end do
       
        Ysubj=Ysubj+R


        GrYsubjA1=GrYsubjA1+GrRA1
       Gr2YsubjA1=Gr2YsubjA1+Gr2RA1
       
       
       do i=1,nvar
          do j=1,kk
             GrYsubjBe(i)=GrYsubjBe(i)+GrQBe(i,j)*vect1(j)
             Gr2YsubjA1Be(i)=Gr2YsubjA1Be(i)+Gr2QA1Be(i,j)*vect1(j)
          end do
       end do  
       do i=1,nvar
          GrYsubjBe(i)=GrYsubjBe(i)+GrRBe(i)
          Gr2YsubjA1Be(i)=Gr2YsubjA1Be(i)+Gr2RA1Be(i)
       end do  
       
       if (kk.gt.1) then
          do j=2,kk
             do i=1,nvar
               do t=1,nvar
                 Gr2YsubjBe(i,t)=Gr2YsubjBe(i,t)+Gr2QBe(j,i,t)
               end do
             end do   
          end do   
       end if
       
       do i=1,nvar
          do j=1,nvar
             Gr2YsubjBe(i,j)=Gr2YsubjBe(i,j)+Gr2RBe(i,j)
          end do
       end do

       return 
       END SUBROUTINE

       
      
       DOUBLE PRECISION FUNCTION rho(k,alpha,rhoFunc)
       implicit none
       integer k,rhoFunc
       double precision alpha
      
       if (rhoFunc.eq.1) then
        rho=1.0 
       end if
       if (rhoFunc.eq.2) then
         rho=alpha**k
       end if 
       
       return
       END FUNCTION


       
       DOUBLE PRECISION FUNCTION psi(nvar,cov,beta,offset)
       implicit none
       
       integer nvar,i
       double precision cov(nvar),beta(nvar),offset
       
        psi=0.d0
       do i=1,nvar
          psi=psi+cov(i)*beta(i)
       end do

       psi=dexp(psi+offset)
       
       return
       END FUNCTION  


***********************************************************************
**
**       Initial functions to maximize Log-likelihood for Xi
** 
**
***********************************************************************
      
       real*8 FUNCTION GammaLogLik(xi,n,A,K)

       implicit none
       
        integer i,j,n
       double precision xi,A(n),G0,q1
        integer K(n),KK(n)
       real gamma
       external gamma
      
        do i=1,n
         KK(i)=K(i)-1
       end do
       G0=0.d0
      
       do i=1,n
         q1=0.d0
        if (KK(i).gt.0) then
          do j=1,KK(i)
           q1=q1+dlog(xi+dble(j)-1.d0)
         end do
        end if
         G0=G0+q1+xi*dlog(xi)-(xi+KK(i))*dlog(xi+A(i))  
       end do
      
       
       GammaLogLik=-G0
       
      
       return
       END FUNCTION





************************************************************************
**
**      Search descent direction procedures
**   
**
************************************************************************


************************************************************************
**
**  Invert a symmetric matrix and calculate its determinant.            
**                                                                      
**                                                                      
**  To call:      CALL MATINV(m,array,det)                   
**                                                                      
**                                                                      
**  m       : dimension of matrix                                       
**  array   : input matrix which is replaced by its inverse            
**  det     : determinant of input matrix                              
**  w1, w2  : work vectors of dimension m
**                                                                      
**                                                                      
**  Reference: Philip B Bevington, "Data Reduction and Error Analysis   
**             for the Physical Sciences", McGraw-Hill, New York, 1969, 
**             pp. 300-303.                                             
**                                                                      
**                                                                      
**                                                                     
************************************************************************

      SUBROUTINE matinv(m,array,det)
        double precision    array(m,m), ik(m), jk(m)
c
   10   det = 1.0
   11   do 100 k = 1, m
c       find largest element array(i,j) in rest of matrix.
        amax = 0.0
   21      do 30 i = k, m
              do 30 j = k, m
   23            if (abs(amax)-abs(array(i,j))) 24,24,30
   24            amax = array(i,j)
                 ik(k) = i
                 jk(k) = j
   30      continue
c          interchange rows and columns to put amax in array(k,k).
   31      if (amax) 41,32,41
   32      det = 0.0
           goto 140
   41      i = ik(k)
           if (i-k) 21,51,43
   43      do 50 j = 1, m
              save = array(k,j)
              array(k,j) = array(i,j)
   50      array(i,j) = -save
   51      j = jk(k)
           if (j-k) 21,61,53
   53      do 60 i = 1, m
              save = array(i,k)
              array(i,k) = array(i,j)
   60      array(i,j) = -save
c          accumulate elements of inverse matrix.
   61      do 70 i = 1, m
              if (i-k) 63,70,63
   63         array(i,k) = -array(i,k)/amax
   70      continue
   71      do 80 i = 1, m
              do 80 j = 1, m
                 if (i-k) 74,80,74
   74            if (j-k) 75,80,75
   75            array(i,j) = array(i,j) + array(i,k)*array(k,j)
   80      continue
   81      do 90 j = 1, m
              if (j-k) 83,90,83
   83         array(k,j) = array(k,j)/amax
   90      continue
           array(k,k) = 1.0/amax        
  100   det = det * amax
c       restore ordering of matrix.
  101   do 130 l = 1, m
           k = m - l + 1
           j = ik(k)
           if (j-k) 111,111,105
  105      do 110 i = 1, m
              save = array(i,k)
              array(i,k) = -array(i,j)
  110      array(i,j) = save
  111      i = jk(k)
           if (i-k) 130,130,113
  113      do 120 j = 1, m
              save = array(k,j)
              array(k,j) = -array(i,j)
  120      array(i,j) = save
  130   continue
  140   return
      
        END SUBROUTINE



      subroutine dschol( n,a,ratm,frac,wa,dinc,info)
      integer           info, n
      double precision  a(n,*), ratm, dinc, frac, wa(n,*)
*
*  n       (input) INTEGER The order of the matrix A.  
*
*  A       (input/output) DOUBLE PRECISION array, dimension (N,N)
*          On entry, the symmetric matrix A.  On output, the diagonal
*          elements of A may have been increased to make A pd.  Only the
*          upper triangle is needed on input.
*
*  ratm   min ratio allowed of smallest to largest diag elements of 
*          U (the Choleski factor)
*  frac   fraction of largest element to add to diag when not pd
*  wa     (output) if INFO = 0, wa is the factor U  from the Choleski
*          factorization of the output matrix A, A = U'U (the lower
*          triangle is completed with 0's)
*  dinc   (output), the amount added to the diag of A
*  INFO   (output) INTEGER
*         <= 0: wa=Choleski factor of output A
*          -2 : all diagonal elements of A were <=0, so A was replaced by
*               an identity matrix
*          > 0: factorization could not be performed on the final A matrix.
*
*  =====================================================================
*
      double precision   one, zero
      parameter          ( one = 1.0d+0, zero = 0.0d+0 )
      integer            i, j
      double precision   dma,dm2,dm3,dmi,doff
      double precision   ddot
      external           ddot
      intrinsic          max, min, sqrt
      dinc=zero
*     Quick return if possible
      if( n.le.0 ) return
      if (n.eq.1) then
         if (a(1,1).gt.zero) then 
            wa(1,1)=sqrt(a(1,1))
            return
         else
            dinc=1-a(1,1)
            wa(1,1)=1
            a(1,1)=1
            info=-2
            return
         endif
      endif
*     determine max and min diag and max abs off-diag elements
      dma=a(1,1)
      dmi=a(1,1)
      doff=zero
      do 5 i=1,n
         dma=max(dma,a(i,i))
         dmi=min(dmi,a(i,i))
         do 6 j=1,i-1
            doff=max(doff,abs(a(j,i)))
 6       continue 
 5    continue 
      if (dma.le.zero) then
         info=-2
         do 7 i=1,n
            do 8 j=1,n
               wa(j,i)=zero
               a(j,i)=zero
 8          continue 
            wa(i,i)=one
            a(i,i)=one
 7       continue 
         return
      endif
c make sure dma > doff, and dmi >= dma*ratm*ratm
      dinc=max(dma*ratm**2-dmi,doff-dma)/(1-ratm**2)
      if (dinc.gt.zero) then
         do 9 i=1,n
            a(i,i)=a(i,i)+dinc
 9       continue 
         dma=dma+dinc
      else
         dinc=zero
      endif
c dm3=base amount to add to diagonal if not pd
      dm3=dma*frac
c in ochol, diagonal elements of factorization required to be > sqrt(dm2)
c assuming largest diag element is approx sqrt(dma), and smallest
c should be > largest*ratm, need dm2=dma*ratm*ratm
 988  dm2=dma*ratm*ratm
c since # rows = n ...
      do 35 i=1,n*n
         wa(i,1)=a(i,1)
 35   continue 
      call ochol(wa,n,dm2,info)
      if (info.gt.0) then
c not pd -- double dm3 and add it to diagonal of A
c adjust dma and dinc accordingly
         dm3=dm3*2
         dma=dma+dm3
         dinc=dinc+dm3
         do 50 i=1,n
            a(i,i)=a(i,i)+dm3
 50      continue 
         go to 988
      endif
      return
      end

c ordinary Choleski decomposition -- return info=j if jth diagonal element 
c of factorization (prior to sqrt) is < damin
      subroutine ochol(a,n,damin,info)
      double precision a(n,n),ajj,ddot,damin
      integer n,info,j,i
      info=0
      do 10 j = 1, n
*  Update jth column
         do 15 i=1,j-1
            a(i,j)=(a(i,j)-ddot(i-1,a(1,i),1,a(1,j),1))/a(i,i)
            a(j,i)=0
 15      continue 
*  Compute U(J,J) and test for non-positive-definiteness.
         ajj = a(j,j)-ddot(j-1,a(1,j),1,a(1,j),1)
         if(ajj.le.damin) then
            info=j
            a(j,j)=ajj
            return
         endif
         a(j,j)=sqrt(ajj)
 10   continue
      return
      end


      
      double precision function ddot(n,dx,incx,dy,incy)
c
c     forms the dot product of two vectors.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double precision dx(*),dy(*),dtemp
      integer i,incx,incy,ix,iy,m,mp1,n
c
      ddot = 0.0d0
      dtemp = 0.0d0
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dtemp = dtemp + dx(ix)*dy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      ddot = dtemp
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dtemp + dx(i)*dy(i)
   30 continue
      if( n .lt. 5 ) go to 60
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) +
     * dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
   50 continue
   60 ddot = dtemp
      return
      end





************************************************************************
**
**      Compute Gamma function
**   
**
************************************************************************





CGAMMA
C   IMSL ROUTINE NAME   - MGAMA=GAMMA                                   MGAA0010
C                                                                       MGAA0020
C-----------------------------------------------------------------------MGAA0030
C                                                                       MGAA0040
C   COMPUTER            - CDCFT5/SINGLE                                 MGAA0050
C                                                                       MGAA0060
C   LATEST REVISION     - JUNE 1, 1982                                  MGAA0070
C                                                                       MGAA0080
C   PURPOSE             - EVALUATE THE GAMMA FUNCTION                   MGAA0090
C                                                                       MGAA0100
C   USAGE               - RESULT = GAMMA(X)                             MGAA0110
C                                                                       MGAA0120
C   ARGUMENTS    X      - INPUT ARGUMENT.                               MGAA0130
C                           GAMMA IS SET TO MACHINE INFINITY, WITH THE  MGAA0140
C                           PROPER SIGN, WHENEVER                       MGAA0150
C                             X IS ZERO,                                MGAA0160
C                             X IS A NEGATIVE INTEGER,                  MGAA0170
C                             ABS(X) .LE. XMIN, OR                      MGAA0180
C                             ABS(X) .GE. XMAX.                         MGAA0190
C                             XMIN IS OF THE ORDER OF 10**(-39) AND     MGAA0200
C                             XMAX IS AT LEAST 34.8. THE EXACT VALUES   MGAA0210
C                             OF XMIN AND XMAX MAY ALLOW LARGER RANGES  MGAA0220
C                             FOR X ON SOME COMPUTERS.                  MGAA0230
C                             SEE THE PROGRAMMING NOTES IN THE MANUAL   MGAA0240
C                             FOR THE EXACT VALUES.                     MGAA0250
C                GAMMA  - OUTPUT SINGLE PRECISION VALUE OF THE GAMMA    MGAA0260
C                           FUNCTION.                                   MGAA0270
C                                                                       MGAA0280
C   PRECISION/HARDWARE  - SINGLE/ALL                                    MGAA0290
C                         NOTE - GAMMA MAY NOT BE SUPPLIED BY IMSL IF   MGAA0300
C                           IT RESIDES IN THE MATHEMATICAL SUBPROGRAM   MGAA0310
C                           LIBRARY SUPPLIED BY THE MANUFACTURER.       MGAA0320
C                                                                       MGAA0330
C   REQD. IMSL ROUTINES - UERTST,UGETIO                                 MGAA0340
C                                                                       MGAA0350
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND           MGAA0360
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL      MGAA0370
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP  MGAA0380
C                                                                       MGAA0390
C   REMARKS      AN ERROR MESSAGE PRINTED BY UERTST FROM GAMMA SHOULD   MGAA0400
C                BE INTERPRETED AS FOLLOWS                              MGAA0410
C                IER    - ERROR INDICATOR                               MGAA0420
C                         TERMINAL ERROR                                MGAA0430
C                           IER = 129 INDICATES THAT THE ABSOLUTE VALUE MGAA0440
C                             OF THE INPUT ARGUMENT X WAS SPECIFIED     MGAA0450
C                             GREATER THAN OR EQUAL TO XMAX. GAMMA      MGAA0460
C                             IS SET TO MACHINE INFINITY.               MGAA0470
C                           IER = 130 INDICATES THAT THE INPUT ARGUMENT MGAA0480
C                             X WAS SPECIFIED AS ZERO OR A NEGATIVE     MGAA0490
C                             INTEGER OR THAT THE ABSOLUTE VALUE OF THE MGAA0500
C                             INPUT ARGUMENT WAS LESS THAN OR EQUAL TO  MGAA0510
C                             XMIN. GAMMA IS SET TO MACHINE INFINITY    MGAA0520
C                             WITH THE PROPER SIGN FOR THE GAMMA        MGAA0530
C                             FUNCTION. IF X IS ZERO OR AN EVEN         MGAA0540
C                             NEGATIVE INTEGER, GAMMA HAS A NEGATIVE    MGAA0550
C                             SIGN. OTHERWISE IT HAS A POSITIVE SIGN.   MGAA0560
C                                                                       MGAA0570
C   COPYRIGHT           - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.       MGAA0580
C                                                                       MGAA0590
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN MGAA0600
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,    MGAA0610
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.        MGAA0620
C                                                                       MGAA0630
C-----------------------------------------------------------------------MGAA0640
C                                                                       MGAA0650
      DOUBLE PRECISION FUNCTION GAMMA (X)                               MGAA0660
C                                  SPECIFICATIONS FOR ARGUMENTS         MGAA0670
      DOUBLE PRECISION   X                                              MGAA0680
C                                  SPECIFICATIONS FOR LOCAL VARIABLES   MGAA0690
      REAL               P,Q,P4,BIG1,PI,XMIN,XINF,XSIGN,Y,T,R,A,TOP,    MGAA0700
     *                   DEN,B                                          MGAA0710
      INTEGER            I,IEND,IEND1,IEND2,IER,J                       MGAA0720
      LOGICAL            MFLAG                                          MGAA0730
      DIMENSION          P(7),Q(6),P4(5)                                MGAA0740
C                                  COEFFICIENTS FOR MINIMAX             MGAA0750
C                                  APPROXIMATION TO GAMMA(X),           MGAA0760
C                                  2.0 .LE. X .LE. 3.0                  MGAA0770
      DATA               P(1)/3.4109112397125E+01/,                     MGAA0780
     1                   P(2)/ -4.8341273405598E+01/,                   MGAA0790
     2                   P(3)/4.3005887829594E+02/,                     MGAA0800
     3                   P(4)/-5.5688734338586E+01/,                    MGAA0810
     4                   P(5)/2.0585220673644E+03/,                     MGAA0820
     5                   P(6)/7.7192407739800E-01/,                     MGAA0830
     6                   P(7)/-3.1721064346240E+00/                     MGAA0840
      DATA               Q(1)/2.4455138217658E+02/,                     MGAA0850
     1                   Q(2)/-1.0174768492818E+03/,                    MGAA0860
     2                   Q(3)/1.1615998272754E+03/,                     MGAA0870
     3                   Q(4)/2.0512896777972E+03/,                     MGAA0880
     4                   Q(5)/6.8080353498091E-01/,                     MGAA0890
     5                   Q(6)/-2.5386729086746E+01/                     MGAA0900
C                                  COEFFICIENTS FOR MINIMAX             MGAA0910
C                                  APPROXIMATION TO LN(GAMMA(X)),       MGAA0920
C                                  12.0 .LE. X                          MGAA0930
      DATA               P4(1)/9.1893853320467E-01/,                    MGAA0940
     1                   P4(2)/8.3333333333267E-02/,                    MGAA0950
     2                   P4(3)/-2.7777776505151E-03/,                   MGAA0960
     3                   P4(4)/7.9358449768E-04/,                       MGAA0970
     4                   P4(5)/-5.82399983E-04/                         MGAA0980
      DATA               IEND/7/,IEND1/6/,IEND2/5/                      MGAA0990
      DATA               XINF/.12650140831E+30/                         MGAA1000
      DATA               PI/3.1415926535898/                            MGAA1010
C                                  GAMMA(XMIN) .APPROX. XINF            MGAA1020
C                                  GAMMA(BIG1) .APPROX. XINF            MGAA1030
      DATA               XMIN/0.0/                                      MGAA1040
      DATA               BIG1/177.803/                                  MGAA1050
C                                  FIRST EXECUTABLE STATEMENT           MGAA1060
      IER = 0                                                           MGAA1070
      MFLAG = .FALSE.                                                   MGAA1080
      T = X                                                             MGAA1090
      IF (ABS(T).GT.XMIN) GO TO 5                                       MGAA1100
      IER = 130                                                         MGAA1110
      GAMMA = XINF                                                      MGAA1120
      IF (T.LE.0.0) GAMMA = -XINF                                       MGAA1130
      GO TO 9000                                                        MGAA1140
    5 IF (ABS(T).LT.BIG1) GO TO 10                                      MGAA1150
      IER = 129                                                         MGAA1160
      GAMMA = XINF                                                      MGAA1170
      GO TO 9000                                                        MGAA1180
   10 IF (T.GT.0.0) GO TO 25                                            MGAA1190
C                                  ARGUMENT IS NEGATIVE                 MGAA1200
      MFLAG = .TRUE.                                                    MGAA1210
      T = -T                                                            MGAA1220
      R = AINT(T)                                                       MGAA1230
      XSIGN = 1.0                                                       MGAA1240
      IF (AMOD(R,2.0) .EQ. 0.0) XSIGN = -1.                             MGAA1250
      R = T-R                                                           MGAA1260
      IF (R.NE.0.0) GO TO 20                                            MGAA1270
      IER = 130                                                         MGAA1280
      GAMMA = XINF                                                      MGAA1290
      IF (XSIGN.EQ.-1.0) GAMMA = -XINF                                  MGAA1300
      GO TO 9000                                                        MGAA1310
C                                  ARGUMENT IS NOT A NEGATIVE INTEGER   MGAA1320
   20 R = PI/SIN(R*PI)*XSIGN                                            MGAA1330
      T = T+1.0                                                         MGAA1340
C                                  EVALUATE APPROXIMATION FOR GAMMA(T)  MGAA1350
C                                    T .GT. XMIN                        MGAA1360
   25 IF (T.GT.12.0) GO TO 60                                           MGAA1370
      I = T                                                             MGAA1380
      A = 1.0                                                           MGAA1390
      IF (I.GT.2) GO TO 40                                              MGAA1400
      I = I+1                                                           MGAA1410
      GO TO (30,35,50),I                                                MGAA1420
C                                  0.0 .LT. T .LT. 1.0                  MGAA1430
   30 A = A/(T*(T+1.0))                                                 MGAA1440
      T = T+2.0                                                         MGAA1450
      GO TO 50                                                          MGAA1460
C                                  1.0 .LE. T .LT. 2.0                  MGAA1470
   35 A = A/T                                                           MGAA1480
      T = T+1.0                                                         MGAA1490
      GO TO 50                                                          MGAA1500
C                                  3.0 .LE. T .LE. 12.0                 MGAA1510
   40 DO 45 J=3,I                                                       MGAA1520
         T = T-1.0                                                      MGAA1530
         A = A*T                                                        MGAA1540
   45 CONTINUE                                                          MGAA1550
C                                  2.0 .LE. T .LE. 3.0                  MGAA1560
   50 TOP = P(IEND1)*T+P(IEND)                                          MGAA1570
      DEN = T+Q(IEND1)                                                  MGAA1580
      DO 55 J=1,IEND2                                                   MGAA1590
         TOP = TOP*T+P(J)                                               MGAA1600
         DEN = DEN*T+Q(J)                                               MGAA1610
   55 CONTINUE                                                          MGAA1620
      Y = (TOP/DEN)*A                                                   MGAA1630
      IF (MFLAG) Y = R/Y                                                MGAA1640
      GAMMA = Y                                                         MGAA1650
      GO TO 9005                                                        MGAA1660
C                                  T .GT. 12.0                          MGAA1670
   60 TOP = ALOG(T)                                                     MGAA1680
      TOP = (T-1.5)*(TOP-1.0)+TOP-1.5                                   MGAA1690
      T = 1.0/T                                                         MGAA1700
      B = T*T                                                           MGAA1710
      Y = (((P4(5)*B+P4(4))*B+P4(3))*B+P4(2))*T+P4(1)+TOP               MGAA1720
      Y = EXP(Y)                                                        MGAA1730
      IF (MFLAG) Y = R/Y                                                MGAA1740
      GAMMA = Y                                                         MGAA1750
      GO TO 9005                                                        MGAA1760
 9000 CONTINUE                                                          MGAA1770
c      CALL UERTST(-IER,' MGAMA')                                        MGAA1780
c      CALL UERTST(IER,'GAMMA ')                                         MGAA1790
 9005 RETURN                                                            MGAA1800
      END                                                               MGAA1810



************************************************************************
**
**      Golden section search procedure
**   
**      Routines modified to use GammaLogLik
************************************************************************

      SUBROUTINE mnbrak(ax,bx,cx,fa,fb,fc,func,nn,KK,AA,BB)

       integer nn
       integer KK(nn)
       double precision AA(nn),BB(nn)
            
       REAL*8 ax,bx,cx,fa,fb,fc,func,GOLD,GLIMIT,TINY
      EXTERNAL func
      PARAMETER (GOLD=1.618034, GLIMIT=100., TINY=1.e-20)
c Given a function func, and given distinct initial points ax and bx, this routine searches
c in the downhill direction (dened by the function as evaluated at the initial points) and
c returns new points ax, bx, cx that bracket a minimum of the function. Also returned are
c the function values at the three points, fa, fb, and fc.
c Parameters: GOLD is the default ratio by which successive intervals are magnied; GLIMIT
c is the maximum magnication allowed for a parabolic-t step.
      REAL*8 dum,fu,q,r,u,ulim
      fa=func(ax,nn,KK,AA,BB)
      fb=func(bx,nn,KK,AA,BB)
      if(fb.gt.fa)then 
      dum=ax
      ax=bx
      bx=dum
      dum=fb
      fb=fa
      fa=dum
      endif
      cx=bx+GOLD*(bx-ax) 
      fc=func(cx,nn,KK,AA,BB)
    1 if(fb.ge.fc)then 
      r=(bx-ax)*(fb-fc) 
      q=(bx-cx)*(fb-fa)
      u=bx-((bx-cx)*q-(bx-ax)*r)/(2.*sign(max(abs(q-r),TINY),q-r))
      ulim=bx+GLIMIT*(cx-bx) 
      if((bx-u)*(u-cx).gt.0.)then 
      fu=func(u,nn,KK,AA,BB)
      if(fu.lt.fc)then 
      ax=bx
      fa=fb
      bx=u
      fb=fu
      return
      else if(fu.gt.fb)then 
      cx=u
      fc=fu
      return
      endif
      u=cx+GOLD*(cx-bx) 
      fu=func(u,nn,KK,AA,BB)
      else if((cx-u)*(u-ulim).gt.0.)then 
      fu=func(u,nn,KK,AA,BB)
      if(fu.lt.fc)then
      bx=cx
      cx=u
      u=cx+GOLD*(cx-bx)
      fb=fc
      fc=fu
      fu=func(u,nn,KK,AA,BB)
      endif
      else if((u-ulim)*(ulim-cx).ge.0.)then 
      u=ulim
      fu=func(u,nn,KK,AA,BB)
      else 
      u=cx+GOLD*(cx-bx)
      fu=func(u,nn,KK,AA,BB)
      endif
      ax=bx 
      bx=cx
      cx=u
      fa=fb
      fb=fc
      fc=fu
      goto 1
      endif
      return
      END
 

      SUBROUTINE brent(ax,bx,cx,f,tol,xmin,nn,KK,AA,BB)

c
c  Function modified to use GammaLogLik
c       
       integer nn
       integer KK(nn)
       double precision AA(nn),BB(nn)
       
       INTEGER ITMAX
       REAL*8 brentFUN,ax,bx,cx,tol,xmin,f,CGOLD,ZEPS
       EXTERNAL f
       PARAMETER (ITMAX=100,CGOLD=.3819660,ZEPS=1.0e-10)
c       Given a function f, and given a bracketing triplet of abscissas ax, bx, cx (such that bx is
c       between ax and cx, and f(bx) is less than both f(ax) and f(cx)), this routine isolates
c       the minimum to a fractional precision of about tol using Brent's method. The abscissa of
c       the minimum is returned as xmin, and the minimum function value is returned as brent,
c       the returned function value.
c       Parameters: Maximum allowed number of iterations; golden ratio; and a small number that
c       protects against trying to achieve fractional accuracy for a minimum that happens to be
c       exactly zero.
       INTEGER iter
       REAL*8 a,b,d,e,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm
       a=min(ax,cx) 
       b=max(ax,cx)
       v=bx 
       w=v
       x=v
       e=0. 
       fx=f(x,nn,KK,AA,BB)
       fv=fx
       fw=fx
       do iter=1,ITMAX 
       xm=0.5*(a+b)
       tol1=tol*abs(x)+ZEPS
       tol2=2.*tol1
       if(abs(x-xm).le.(tol2-.5*(b-a))) goto 3 
       if(abs(e).gt.tol1) then 
       r=(x-w)*(fx-fv)
       q=(x-v)*(fx-fw)
       p=(x-v)*q-(x-w)*r
       q=2.*(q-r)
       if(q.gt.0.) p=-p
       q=abs(q)
       etemp=e
       e=d
       if(abs(p).ge.abs(.5*q*etemp).or.p.le.q*(a-x).or.
     * p.ge.q*(b-x)) goto 1
c       The above conditions determine the acceptability of the parabolic t. Here it is o.k.:
c       d=p/q Take the parabolic step.
       d=p/q
       u=x+d
       if(u-a.lt.tol2 .or. b-u.lt.tol2) d=sign(tol1,xm-x)
       goto 2 
       endif
    1 if(x.ge.xm) then 
       e=a-x
       else
       e=b-x
       endif
       d=CGOLD*e 
    2 if(abs(d).ge.tol1) then 
       u=x+d 
       else 
       u=x+sign(tol1,d)
       endif
       fu=f(u,nn,KK,AA,BB) 
       if(fu.le.fx) then 
       if(u.ge.x) then
       a=x
       else
       b=x
       endif
       v=w
       fv=fw
       w=x
       fw=fx
       x=u
       fx=fu
       else
       if(u.lt.x) then
       a=u
       else
       b=u
       endif
       if(fu.le.fw .or. w.eq.x) then
       v=w
       fv=fw
       w=u
       fw=fu
       else if(fu.le.fv .or. v.eq.x .or. v.eq.w) then
       v=u
       fv=fu
       endif
       endif 
       enddo 
c       pause 'brent exceed maximum iterations'
    3 xmin=x 
       brentFUN=fx
       return
       END SUBROUTINE


      FUNCTION golden(ax,bx,cx,f,tol,xmin,nn,KK,AA,BB)
      
       integer nn
       integer KK(nn)
       double precision AA(nn),BB(nn)
       
       
       REAL*8 golden,ax,bx,cx,tol,xmin,f,R,C
      EXTERNAL f
      PARAMETER (R=.61803399,C=1.-R)
c Given a function f, and given a bracketing triplet of abscissas ax, bx, cx (such that bx is
c between ax and cx, and f(bx) is less than both f(ax) and f(cx)), this routine performs
c a golden section search for the minimum, isolating it to a fractional precision of about
c tol. The abscissa of the minimum is returned as xmin, and the minimum function value
c is returned as golden, the returned function value.
c Parameters: The golden ratios.
      REAL*8 f1,f2,x0,x1,x2,x3
      x0=ax
      x3=cx
      if(abs(cx-bx).gt.abs(bx-ax))then 
      x1=bx
      x2=bx+C*(cx-bx) 
      else
      x2=bx
      x1=bx-C*(bx-ax)
      endif
      f1=f(x1,nn,KK,AA,BB) 
      f2=f(x2,nn,KK,AA,BB)
    1 if(abs(x3-x0).gt.tol*(abs(x1)+abs(x2)))then 
      if(f2.lt.f1)then 
      x0=x1 
      x1=x2
      x2=R*x1+C*x3
      f1=f2
      f2=f(x2,nn,KK,AA,BB)
      else 
      x3=x2
      x2=x1
      x1=R*x2+C*x0
      f2=f1
      f1=f(x1,nn,KK,AA,BB) 
      endif
      goto 1 
      endif
      if(f1.lt.f2)then 
      golden=f1
      xmin=x1
      else
      golden=f2
      xmin=x2
      endif
      return
      END

   
