CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C    Procedures for penalized frailty
C
C                lambda(t)=lambda_0(t)exp(betaX+Zw)
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC	

	
 


      SUBROUTINE likpenal(sIn,nIn,nvar,kIn,nk,tau,caltimes,gaptimes,
     .  censored,intercepts,slopes,lastperrep,perrepind,
     .  effagebegin,effage,cov,alpha,beta,Z,rhoFunc,ns,xi,loglikpen,
     .  loglik,tol,X,hessian)


      implicit none
	
      integer nIn,nvar,rhoFunc
      integer kIn(nIn),nk,i,j,pos,r,t
      double precision sIn,tau(nIn),caltimes(nk),gaptimes(nk),
     .  censored(nIn),intercepts(nk),slopes(nk),lastperrep(nk),
     .  perrepind(nk),effagebegin(nk),effage(nk),cov(nvar,nk),
     .  covariate(nvar),alpha,beta(nvar),Z(nIn),Ysubj

      integer ns(nIn)    	     
      double precision S0,loglik,offset(nIn),pen,xi,sumZpen,loglikpen

      INTEGER IPRINT,LP,MAXFUN,MODE,NFUN
      COMMON /VA13BD/ IPRINT,LP,MAXFUN,MODE,NFUN
      
      EXTERNAL FUNC
      INTEGER NP,tamH
      parameter(tamH=100000)
      DOUBLE PRECISION X(nvar+1),SCORE(nvar+1),SCALE(nvar+1),ACC,
     .     hessian(tamH),tol

c Note: tol es ACC  y hessian es W en routinas BFGS


      NP=nvar+1

      SCALE(1)=0.1
      do i=2,NP
       SCALE(i)=1.0d0
      end do 

      
      X(1)=alpha
      do i=2,NP
       X(i)=beta(i-1) 
      end do    

      do i=1,tamH
        hessian(i)=0.d0
      end do 


      do i=1,nIn
        offset(i)=Z(i)
      end do 

      IPRINT=0

     
      call VA13AD(FUNC,NP,X,loglik,SCORE,SCALE,tol,hessian,
     .  sIn,nIn,nvar,kIn,nk,tau,caltimes,gaptimes,
     .  censored,intercepts,slopes,lastperrep,perrepind,
     .  effagebegin,effage,cov,offset,rhoFunc,
     .  ns,xi,loglikpen)

 
      return
      END SUBROUTINE 


	

   
      SUBROUTINE FUNC(NPAR,B,F0,SCORE,s,n,nvar,k,nk,tau,caltimes,
     . gaptimes,censored,intercepts,slopes,lastperrep,perrepind,
     . effagebegin,effage,cov,offset,rhoFunc,ns,xi,loglikpen)


      implicit none
	
      integer n,nvar,rhoFunc
      integer k(n),nk,i,j,pos,r,t
      double precision s,tau(n),caltimes(nk),gaptimes(nk),censored(n),
     .       intercepts(nk),slopes(nk),lastperrep(nk),perrepind(nk),
     .       effagebegin(nk),effage(nk),cov(nvar,nk),covariate(nvar),
     .       alpha,beta(nvar),Z(n),w,Ysubj
	
      double precision  effageOK(200),ZERO,rho,psi
      
      integer ns(n),NPAR
            	     
      double precision S0,loglik,offset(n),pen,xi,loglikpen,
     .     B(NPAR),F0,SCORE(NPAR),scorealpha,scorebeta(nvar),GrS0A1,
     .     GrS0Be(nvar)
         
     
      ZERO=0.d0 
      
c     initialization
  	
      loglik=ZERO
      loglikpen=ZERO
      pen=ZERO  
	
      scorealpha=ZERO	
      do i=1,nvar
	scorebeta(i)=ZERO
      end do      


      alpha=B(1)
      do i=1,NPAR-1
       beta(i)=B(i+1)
      end do  


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

	     call AtRisk2(s,ns,w,n,nvar,nk,k,tau,caltimes,gaptimes,
     .        censored,intercepts,slopes,lastperrep,perrepind,
     .        effagebegin,effage,cov,alpha,beta,offset,rhoFunc,
     .        S0,GrS0A1,GrS0Be)	     



              scorealpha=scorealpha+((dble(j-2)/alpha)-(GrS0A1/S0))
              do t=1,nvar
	        scorebeta(t)=scorebeta(t)+covariate(t)-(GrS0Be(t)/S0) 
   	      end do

              loglik=loglik+dlog(rho(j-2,alpha,rhoFunc))+
     .		            dlog(psi(nvar,covariate,beta,offset(i)))
     .                      -dlog(S0)
	      

            end do   
	  end if
            
        end do

        do i=1,n
         pen=pen+(Z(i)-dexp(Z(i)))
        end do
        pen=pen*(1/xi) 
	    
        loglikpen=loglik+pen  

        F0=-loglik

        score(1) = -scorealpha

        do i=1,NPAR-1
         score(i+1) = -scorebeta(i)
        end do 


	return

        END SUBROUTINE 





       SUBROUTINE AtRisk2(s,ns,w,n,nvar,nk,k,tau,caltimes,gaptimes,
     .  censored,intercepts,slopes,lastperrep,perrepind,
     .  effagebegin,effage,cov,alpha,beta,offset,rhoFunc,S0,GrS0A1,
     .  GrS0Be)

	implicit none
	
       integer n,nvar,nk,kk,i,j,t,pos,rhoFunc
       integer k(n),ns(n),nsi
       double precision w,s,tau(n),caltimes(nk),gaptimes(nk),
     .  censored(n),intercepts(nk),slopes(nk),lastperrep(nk),
     .  perrepind(nk),effagebegin(nk),effage(nk),cov(nvar,nk)

      double precision caltimesOK(200),gaptimesOK(200),
     .  censoredOK,interceptsOK(200),slopesOK(200),
     .  lastperrepOK(200),perrepindOK(200),effagebeginOK(200),
     .  effageOK(200),covOK(nvar,200),GrS0A1,GrS0Be(nvar),
     .  GrYsubjA1,GrYsubjBe(nvar) 
	
      double precision S0,alpha,beta(nvar),ZERO,Ysubj,offset(n)



c     initialization	
	ZERO=0.d0
c     initialitation 
        S0=ZERO
	Ysubj=ZERO
 	GrS0A1=ZERO
 	do i=1,nvar
	  GrS0Be(i)=ZERO
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

        call AtRiskSubj2(n,s,nsi,w,nvar,k(i),tau(i),caltimesOK,
     .   gaptimesOK,censoredOK,interceptsOK,slopesOK,lastperrepOK,
     .	 perrepindOK,effagebeginOK,effageOK,covOK,alpha,beta,offset(i),
     .   rhoFunc,Ysubj,GrYsubjA1,GrYsubjBe)
         
         
          S0=S0+Ysubj
          GrS0A1=GrS0A1+GrYsubjA1
          do j=1,nvar
            GrS0Be(j)=GrS0Be(j)+GrYsubjBe(j)
          end do

	end do
	
	return
	END SUBROUTINE






      SUBROUTINE AtRiskSubj2(n,s,nsi,w,nvar,k,tau,caltimes,gaptimes,
     .  censored,intercepts,slopes,lastperrep,perrepind,effagebegin,
     .  effage,cov,alpha,beta,offset,rhoFunc,Ysubj,GrYsubjA1,GrYsubjBe)


      implicit none
	
      integer n,nsi,nvar,k,kk,i,j,t,rhoFunc
      double precision w,s,tau,caltimes(k),gaptimes(k),
     .	      censored,intercepts(k),slopes(k),lastperrep(k),
     .	      perrepind(k),effagebegin(k),effage(k),cov(nvar,k)
      
      double precision Q(nsi),R,alpha,beta(nvar),ONE,ZERO
     
      double precision covariate(nvar),effageats,Ysubj,
     .   rho,psi,offset,GrYsubjA1,GrYsubjBe(nvar),GrQA1(nsi),
     .   GrQBe(nvar,nsi),GrRA1,GrRBe(nvar)    
  		
 
      ONE=1.d0
      ZERO=0.d0




c     if nsi is NULL is not implemented  
c        take care with changes in nsi and array sizes !!!

      kk=nsi 


c     Initialization 
	R=ZERO
        GrRA1=ZERO
	do i =1,kk
	  Q(i)=ZERO
          GrQA1(i)=ZERO 
          do j=1,nvar
	    GrQBe(j,i)=ZERO
          end do
	end do
        do i=1,nvar
         GrRBe(i)=ZERO
        end do  
   
        
            
c     Equation (23) Q_ij
      if (kk.gt.1) then
        do j=2,kk
         if ((w.gt.effagebegin(j-1)).and.(w.le.effage(j))) then
          
	    do i=1,nvar
              covariate(i)=cov(i,j-1)           
            end do


	    Q(j)=(rho(j-2,alpha,rhoFunc)*psi(nvar,covariate,beta,offset))
     .                         		/slopes(j-1)

            GrQA1(j)=(dble(j-2)/alpha)*Q(j)
            do i=1,nvar
 	       GrQBe(i,j)=covariate(i)*Q(j)
            end do
	   	   
	   end if
	  end do
        end if

	
	
      
c       Equation (24) R_i
	effageats=intercepts(kk)+slopes(kk)*min(s,tau)
     .	    -caltimes(lastperrep(kk))	
	

      if ((w.gt.effagebegin(kk)).and.(w.le.effageats)) then

	   do i=1,nvar
	      covariate(i)=cov(i,kk)           
	   end do

 
           R=(rho(kk-1,alpha,rhoFunc)*
     .	      psi(nvar,covariate,beta,offset))/slopes(kk)

           GrRA1=(dble(kk-1)/alpha)*R
	   do i=1,nvar
	      GrRBe(i)=covariate(i)*R
           end do


       end if
	

c       OUTPUT 

c       Initialization 
	Ysubj=ZERO
        GrYsubjA1=ZERO
        do i=1,nvar
	   GrYsubjBe(i)=ZERO
        end do


c       computes
	do i=1,kk
	   Ysubj=Ysubj+Q(i)
           GrYsubjA1=GrYsubjA1+GrQA1(i)
	end do

        do i=1,nvar
	   do j=1,kk
	      GrYsubjBe(i)=GrYsubjBe(i)+GrQBe(i,j)
	   end do
	end do  
	
        Ysubj=Ysubj+R
        GrYsubjA1=GrYsubjA1+GrRA1
        do i=1,nvar
	   GrYsubjBe(i)=GrYsubjBe(i)+GrRBe(i)
	end do	  

        

	return  	
	END SUBROUTINE





CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc
C
C   BFGS routines modified 
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc




      SUBROUTINE VA13AD(FUNC,N,X,F,G,SCALE,ACC,W,
     .  sIn,nIn,nvar,kIn,nk,tau,caltimes,gaptimes,
     .  censored,intercepts,slopes,lastperrep,perrepind,
     .  effagebegin,effage,cov,offset,rhoFunc,
     .  ns,xi,loglikpen)

      IMPLICIT integer*4 (i-n), REAL*8(A-H,O-Z)
      COMMON /VA13BD/ IPRINT,LP,MAXFUN,MODE,NFUN
      integer tamH
      parameter(tamH=100000)
      DIMENSION X(N),G(N),SCALE(N),W(tamH)
      EXTERNAL FUNC  
c ......... Esto es el COMMON en chapuza
         integer nIn,nvar,nk,rhoFunc 
         integer kIn(nIn),ns(nIn)
         double precision sIn,tau(nIn),caltimes(nk),gaptimes(nk),
     1   censored(nIn),intercepts(nk),slopes(nk),
     2   lastperrep(nk),perrepind(nk),effagebegin(nk),
     3   effage(nk),cov(nvar,nk),offset(nIn),xi,loglikpen

      ND=1+(N*(N+1))/2
      NW=ND+N
      NXA=NW+N
      NGA=NXA+N
      NXB=NGA+N
      NGB=NXB+N

      CALL VA13CD(FUNC,N,X,F,G,SCALE,ACC,W,W(ND),W(NW),
     1  W(NXA),W(NGA),W(NXB),W(NGB),
     .  sIn,nIn,nvar,kIn,nk,tau,caltimes,gaptimes,
     .  censored,intercepts,slopes,lastperrep,perrepind,
     .  effagebegin,effage,cov,offset,rhoFunc,
     .  ns,xi,loglikpen)
      RETURN
      END
      BLOCK DATA
      COMMON /VA13BD/ IPRINT,LP,MAXFUN,MODE,NFUN
      INTEGER*4 IPRINT/0/,LP/6/,MAXFUN/0/,MODE/1/
      END



      SUBROUTINE VA13CD (FUNC,N,X,F,G,SCALE,ACC,H,D,W,XA,GA,XB,GB,
     .  sIn,nIn,nvar,kIn,nk,tau,caltimes,gaptimes,
     .  censored,intercepts,slopes,lastperrep,perrepind,
     .  effagebegin,effage,cov,offset,rhoFunc,
     .  ns,xi,loglikpen)

      IMPLICIT integer*4 (i-n), REAL*8 (A-H,O-Z)
      COMMON /VA13BD/ IPRINT,LP,MAXFUN,MODE,NFUN

      integer tamH
      parameter(tamH=100000)
 
      DIMENSION X(N),G(N),SCALE(N),H(tamH),D(tamH),W(tamH),
     1XA(tamH),GA(tamH),XB(tamH),GB(tamH)
      EXTERNAL FUNC
c ......... Esto es el COMMON en chapuza
         integer nIn,nvar,nk,rhoFunc 
         integer kIn(nIn),ns(nIn)
         double precision sIn,tau(nIn),caltimes(nk),gaptimes(nk),
     1   censored(nIn),intercepts(nk),slopes(nk),
     2   lastperrep(nk),perrepind(nk),effagebegin(nk),
     3   effage(nk),cov(nvar,nk),offset(nIn),xi,loglikpen

C     BEGIN THE PRINTING FROM THE SUBROUTINE

      IF (IPRINT.EQ.0) GO TO 10
      WRITE (LP,1000)
 1000 FORMAT ('1ENTRY TO VA13AD')
C     CALCULATE THE INITIAL FUNCTION VALUE
   10 CALL FUNC (N,X,F,G,sIn,nIn,nvar,kIn,nk,tau,caltimes,
     . gaptimes,censored,intercepts,slopes,lastperrep,perrepind,
     .  effagebegin,effage,cov,offset,rhoFunc,ns,xi,loglikpen)

      NFUN=1
      ITR=0
      NP=N+1
C     SET THE HESSIAN TO A DIAGONAL MATRIX DEPENDING ON SCALE(.)
      IF (MODE.GE.2) GO TO 60
   20 C=0.D0
      DO 30 I=1,N
   30 C=DMAX1(C,DABS(G(I)*SCALE(I)))
      IF (C.LE.0.D0) C=1.D0
      K=(N*NP)/2
      DO 40 I=1,K
   40 H(I)=0.D0
      K=1
      DO 50 I=1,N
      H(K)=0.01D0*C/SCALE(I)**2
   50 K=K+NP-I
      GO TO 100
C     FACTORIZE THE GIVEN HESSIAN MATRIX
   60 IF (MODE.GE.3) GO TO 80
      CALL MC11BD (H,N,K)
      IF (K.GE.N) GO TO 100
   70 WRITE (LP,1010)
 1010 FORMAT (/5X,'BECAUSE THE HESSIAN GIVEN TO VA13AD IS NOT POS DEF,'/
     15X,'IT HAS BEEN REPLACED BY A POSITIVE DIAGONAL MATRIX')
      GO TO 20
C     CHECK THAT THE GIVEN DIAGONAL ELEMENTS ARE POSITIVE
   80 K=1
      DO 90 I=1,N
      IF (H(K).LE.0D0) GO TO 70
   90 K=K+NP-I
C     SET SOME VARIABLES FOR THE FIRST ITERATION
  100 DFF=0.D0
      IPRA=IABS(IPRINT)
      IP=IABS(IPRA-1)
  110 FA=F
      ISFV=1
      DO 120 I=1,N
      XA(I)=X(I)
  120 GA(I)=G(I)
C     BEGIN THE ITERATION BY GIVING THE REQUIRED PRINTING
  130 IP=IP+1
      IF (IP.NE.IPRA) GO TO 140
      IP=0
      WRITE (LP,1020) ITR,NFUN
 1020 FORMAT (/5X,'ITERATION =',I5,5X,'FUNCTIONS =',I5)
      WRITE (LP,1030) FA
 1030 FORMAT (5X,'F =',D24.16)
      IF (IPRINT.LE.0) GO TO 140
      WRITE (LP,1040) (XA(I),I=1,N)
 1040 FORMAT (5X,'X(.) =',(5D24.16))
      WRITE (LP,1050) (GA(I),I=1,N)
 1050 FORMAT (5X,'G(.) =',(5D24.16))
  140 ITR=ITR+1
C     CALCULATE THE SEARCH DIRECTION OF THE ITERATION
      DO 150 I=1,N
  150 D(I)=-GA(I)
      CALL MC11ED (H,N,D,W,N)
C     CALCULATE A LOWER BOUND ON THE STEP-LENGTH
C     AND THE INITIAL DIRECTIONAL DERIVATIVE
      C=0.D0
      DGA=0.D0
      DO 160 I=1,N
      C=DMAX1(C,DABS(D(I)/SCALE(I)))
  160 DGA=DGA+GA(I)*D(I)
C     TEST IF THE SEARCH DIRECTION IS DOWNHILL
      IF (DGA.GE.0.D0) GO TO 240
C     SET THE INITIAL STEP-LENGTH OF THE LINE SEARCH
      STMIN=0.D0
      STEPBD=0.D0
      STEPLB=ACC/C
      FMIN=FA
      GMIN=DGA
      STEP=1.D0
      IF (DFF.LE.0D0) STEP=DMIN1(STEP,1D0/C)
      IF (DFF.GT.0D0) STEP=DMIN1(STEP,(DFF+DFF)/(-DGA))
  170 C=STMIN+STEP
C     TEST WHETHER FUNC HAS BEEN CALLED MAXFUN TIMES
      IF (NFUN.EQ.MAXFUN) GO TO 250
      NFUN=NFUN+1
C     CALCULATE ANOTHER FUNCTION VALUE AND GRADIENT
      DO 180 I=1,N
  180 XB(I)=XA(I)+C*D(I)

       

      CALL FUNC (N,XB,FB,GB,sIn,nIn,nvar,kIn,nk,tau,caltimes,
     . gaptimes,censored,intercepts,slopes,lastperrep,perrepind,
     .  effagebegin,effage,cov,offset,rhoFunc,ns,xi,loglikpen)
C     STORE THIS FUNCTION VALUE IF IT IS THE SMALLEST SO FAR

c      write(*,*) XB
c      write(*,*) FB
c      write(*,*) GB

      ISFV=MIN0(2,ISFV)
      IF (FB.GT.F) GO TO 220
      IF (FB.LT.F) GO TO 200
      GL1=0.D0
      GL2=0.D0
      DO 190 I=1,N
      GL1=GL1+(SCALE(I)*G(I))**2
  190 GL2=GL2+(SCALE(I)*GB(I))**2
      IF (GL2.GE.GL1) GO TO 220
  200 ISFV=3
      F=FB
      DO 210 I=1,N
      X(I)=XB(I)
  210 G(I)=GB(I)
C     CALCULATE THE DIRECTIONAL DERIVATIVE AT THE NEW POINT
  220 DGB=0.D0
      DO 230 I=1,N
  230 DGB=DGB+GB(I)*D(I)
C     BRANCH IF WE HAVE FOUND A NEW LOWER BOUND ON THE STEP-LENGTH
      IF (FB-FA.LE.0.1D0*C*DGA) GO TO 280
C     FINISH THE ITERATION IF THE CURRENT STEP IS STEPLB
      IF (STEP.GT.STEPLB) GO TO 270
  240 IF (ISFV.GE.2) GO TO 110
C     AT THIS STAGE THE WHOLE CALCULATION IS COMPLETE
  250 IF (IPRINT.EQ.0) GO TO 260
      WRITE (LP,1070)
 1070 FORMAT (/5X,'THE RESULTS FROM VA13AD ARE AS FOLLOWS')
      WRITE (LP,1020) ITR,NFUN
      WRITE (LP,1030) F
      WRITE (LP,1040) (X(I),I=1,N)
      WRITE (LP,1050) (G(I),I=1,N)
  260 RETURN
C     CALCULATE A NEW STEP-LENGTH BY CUBIC INTERPOLATION
  270 STEPBD=STEP
      C=GMIN+DGB-3.D0*(FB-FMIN)/STEP
      CC=DSQRT(C*C-GMIN*DGB)
      C=(C-GMIN+CC)/(DGB-GMIN+CC+CC)
      STEP=STEP*DMAX1(0.1D0,C)
      GO TO 170
C     SET THE NEW BOUNDS ON THE STEP-LENGTH
  280 STEPBD=STEPBD-STEP
      STMIN=C
      FMIN=FB
      GMIN=DGB
C     CALCULATE A NEW STEP-LENGTH BY EXTRAPOLATION
      STEP=9.D0*STMIN
      IF (STEPBD.GT.0.D0) STEP=0.5D0*STEPBD
      C=DGA+3.D0*DGB-4.D0*(FB-FA)/STMIN
      IF (C.GT.0D0) STEP=DMIN1(STEP,STMIN*DMAX1(1D0,-DGB/C))
      IF (DGB.LT.0.7D0*DGA) GO TO 170
C     TEST FOR CONVERGENCE OF THE ITERATIONS
      ISFV=4-ISFV
      IF (STMIN+STEP.LE.STEPLB) GO TO 240
C     REVISE THE SECOND DERIVATIVE MATRIX
      IR=-N
      DO 290 I=1,N
      XA(I)=XB(I)
      XB(I)=GA(I)
      D(I)=GB(I)-GA(I)
  290 GA(I)=GB(I)
      CALL MC11AD(H,N,XB,1.D0/DGA,W,IR,1,0.D0)
      IR=-IR
      CALL MC11AD(H,N,D,1D0/(STMIN*(DGB-DGA)),D,IR,0,0D0)
C     BRANCH IF THE RANK OF THE NEW MATRIX IS DEFICIENT
      IF (IR.LT.N) GO TO 250
C     BEGIN ANOTHER ITERATION
      DFF=FA-FB
      FA=FB
      GO TO 130
      END


      SUBROUTINE MC11AD(A,N,Z,SIG,W,IR,MK,EPS)
C  STANDARD FORTRAN 66 (A VERIFIED PFORT SUBROUTINE)
      implicit integer*4 (i-n)
      DOUBLE PRECISION A(1),AL,B,EPS,GM,R,SIG,TI,TIM,V,W(1),Y,Z(1)
C   UPDATE FACTORS GIVEN IN A BY   SIG*Z*ZTRANSPOSE
      IF(N.GT.1)GOTO1
      A(1)=A(1)+SIG *Z(1)**2
      IR=1
      IF(A(1).GT.0.0D0)RETURN
      A(1)=0.D0
      IR=0
      RETURN
    1 CONTINUE
      NP=N+1
      IF(SIG.GT.0.0D0)GOTO40
      IF(SIG.EQ.0.0D0.OR.IR.EQ.0)RETURN
      TI=1.0D0/SIG
      IJ=1
      IF(MK.EQ.0)GOTO10
      DO 7 I=1,N
      IF(A(IJ).NE.0.0D0)TI=TI+W(I)**2/A(IJ)
    7 IJ=IJ+NP-I
      GOTO20
   10 CONTINUE
      DO 11 I=1,N
   11 W(I)=Z(I)
      DO 15 I=1,N
      IP=I+1
      V=W(I)
      IF(A(IJ).GT.0.0D0)GOTO12
      W(I)=0.D0
      IJ=IJ+NP-I
      GOTO15
   12 CONTINUE
      TI=TI+V**2/A(IJ)
      IF(I.EQ.N)GOTO14
      DO 13 J=IP,N
      IJ=IJ+1
   13 W(J)=W(J)-V*A(IJ)
   14 IJ=IJ+1
   15 CONTINUE
   20 CONTINUE
      IF(IR.LE.0)GOTO21
      IF(TI.GT.0.0D0)GOTO22
      IF(MK-1)40,40,23
   21 TI=0.D0
      IR=-IR-1
      GOTO23
   22 TI=EPS/SIG
      IF(EPS.EQ.0.0D0)IR=IR-1
   23 CONTINUE
      MM=1
      TIM=TI
      DO 30 I=1,N
      J=NP-I
      IJ=IJ-I
      IF(A(IJ).NE.0.0D0)TIM=TI-W(J)**2/A(IJ)
      W(J)=TI
   30 TI=TIM
      GOTO41
   40 CONTINUE
      MM=0
      TIM=1.0D0/SIG
   41 CONTINUE
      IJ=1
      DO 66 I=1,N
      IP=I+1
      V=Z(I)
      IF(A(IJ).GT.0.0D0)GOTO53
      IF(IR.GT.0.OR.SIG.LT.0.0D0.OR.V.EQ.0.0D0)GOTO52
      IR=1-IR
      A(IJ)=V**2/TIM
      IF(I.EQ.N)RETURN
      DO 51 J=IP,N
      IJ=IJ+1
   51 A(IJ)=Z(J)/V
      RETURN
   52 CONTINUE
      TI=TIM
      IJ=IJ+NP-I
      GOTO66
   53 CONTINUE
      AL=V/A(IJ)
      IF(MM)54,54,55
   54 TI=TIM+V*AL
      GOTO56
   55 TI=W(I)
   56 CONTINUE
      R=TI/TIM
      A(IJ)=A(IJ)*R
      IF(R.EQ.0.0D0)GOTO70
      IF(I.EQ.N)GOTO70
      B=AL/TI
      IF(R.GT.4.0D0)GOTO62
      DO 61 J=IP,N
      IJ=IJ+1
      Z(J)=Z(J)-V*A(IJ)
   61 A(IJ)=A(IJ)+B*Z(J)
      GOTO64
   62 GM=TIM/TI
      DO 63 J=IP,N
      IJ=IJ+1
      Y=A(IJ)
      A(IJ)=B*Z(J)+Y*GM
   63 Z(J)=Z(J)-V*Y
   64 CONTINUE
      TIM=TI
      IJ=IJ+1
   66 CONTINUE
   70 CONTINUE
      IF(IR.LT.0)IR=-IR
      RETURN
      END


      SUBROUTINE MC11BD(A,N,IR)
C   FACTORIZE A MATRIX GIVEN IN A
      implicit integer*4 (i-n)
      DOUBLE PRECISION A(1),AA,V
      IR=N
      IF(N.GT.1)GOTO100
      IF(A(1).GT.0.0D0)RETURN
      A(1)=0.D0
      IR=0
      RETURN
  100 CONTINUE
      NP=N+1
      II=1
      DO 104 I=2,N
      AA=A(II)
      NI=II+NP-I
      IF(AA.GT.0.0D0)GOTO101
      A(II)=0.D0
      IR=IR-1
      II=NI+1
      GOTO104
  101 CONTINUE
      IP=II+1
      II=NI+1
      JK=II
      DO 103 IJ=IP,NI
      V=A(IJ)/AA
      DO 102 IK=IJ,NI
      A(JK)=A(JK)-A(IK)*V
  102 JK=JK+1
  103 A(IJ)=V
  104 CONTINUE
      IF(A(II).GT.0.0D0)RETURN
      A(II)=0.D0
      IR=IR-1
      RETURN
      END


      SUBROUTINE MC11CD(A,N)
C   MULTIPLY OUT THE FACTORS GIVEN IN A
      implicit integer*4 (i-n)
      DOUBLE PRECISION A(1),AA,V
      IF(N.EQ.1)RETURN
      NP=N+1
      II=N*NP/2
      DO 202 NIP=2,N
      JK=II
      NI=II-1
      II=II-NIP
      AA=A(II)
      IP=II+1
      IF(AA.GT.0.0D0)GOTO203
      DO 204 IJ=IP,NI
  204 A(IJ)=0.D0
      GOTO202
  203 CONTINUE
      DO 201 IJ=IP,NI
      V=A(IJ)*AA
      DO 200 IK=IJ,NI
      A(JK)=A(JK)+A(IK)*V
  200 JK=JK+1
  201 A(IJ)=V
  202 CONTINUE
      RETURN
      END


      SUBROUTINE MC11DD(A,N,Z,W)
      implicit integer*4 (i-n)
      DOUBLE PRECISION A(1),W(1),Z(1),Y
C   MULTIPLY A VECTOR Z BY THE FACTORS GIVEN IN A
      IF(N.GT.1)GOTO300
      Z(1)=Z(1)*A(1)
      W(1)=Z(1)
      RETURN
  300 CONTINUE
      NP=N+1
      II=1
      N1=N-1
      DO 303 I=1,N1
      Y=Z(I)
      IF(A(II).EQ.0.0D0)GOTO302
      IJ=II
      IP=I+1
      DO 301 J=IP,N
      IJ=IJ+1
  301 Y=Y+Z(J)*A(IJ)
  302 Z(I)=Y*A(II)
      W(I)=Z(I)
  303 II=II+NP-I
      Z(N)=Z(N)*A(II)
      W(N)=Z(N)
      DO 311 K=1,N1
      I=N-K
      II=II-NP+I
      IF(Z(I).EQ.0.0D0)GOTO311
      IP=I+1
      IJ=II
      Y=Z(I)
      DO 310 J=IP,N
      IJ=IJ+1
  310 Z(J)=Z(J)+A(IJ)*Z(I)
  311 CONTINUE
      RETURN
      END


      SUBROUTINE MC11ED(A,N,Z,W,IR)
      implicit integer*4 (i-n)
      DOUBLE PRECISION A(1),V,W(1),Z(1)
C   MULTIPLY A VECTOR Z BY THE INVERSE OF THE FACTORS GIVEN IN A
      IF(IR.LT.N)RETURN
      W(1)=Z(1)
      IF(N.GT.1)GOTO400
      Z(1)=Z(1)/A(1)
      RETURN
  400 CONTINUE
      DO 402 I=2,N
      IJ=I
      I1=I-1
      V=Z(I)
      DO 401 J=1,I1
      V=V-A(IJ)*Z(J)
  401 IJ=IJ+N-J
      W(I)=V
  402 Z(I)=V
      Z(N)=Z(N)/A(IJ)
      NP=N+1
      DO 411 NIP=2,N
      I=NP-NIP
      II=IJ-NIP
      V=Z(I)/A(II)
      IP=I+1
      IJ=II
      DO 410 J=IP,N
      II=II+1
  410 V=V-A(II)*Z(J)
  411 Z(I)=V
      RETURN
      END


      SUBROUTINE MC11FD(A,N,IR)
      implicit integer*4 (i-n)
      DOUBLE PRECISION A(1),AA,V
C   COMPUTE THE INVERSE MATRIX FROM FACTORS GIVEN IN A
      IF(IR.LT.N)RETURN
      A(1)=1.0D0/A(1)
      IF(N.EQ.1)RETURN
      NP=N+1
      N1=N-1
      II=2
      DO 511 I=2,N
      A(II)=-A(II)
      IJ=II+1
      IF(I.EQ.N)GOTO502
      DO 501 J=I,N1
      IK=II
      JK=IJ
      V=A(IJ)
      DO 500 K=I,J
      JK=JK+NP-K
      V=V+A(IK)*A(JK)
  500 IK=IK+1
      A(IJ)=-V
  501 IJ=IJ+1
  502 CONTINUE
      A(IJ)=1.0D0/A(IJ)
      II=IJ+1
      AA=A(IJ)
      IJ=I
      IP=I+1
      NI=N-I
      DO 511 J=2,I
      V=A(IJ)*AA
      IK=IJ
      K=IJ-IP+J
      I1=IJ-1
      NIP=NI+IJ
      DO 510 JK=K,I1
      A(JK)=A(JK)+V*A(IK)
  510 IK=IK+NIP-JK
      A(IJ)=V
  511 IJ=IJ+NP-J
      RETURN
      END




CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc


	
