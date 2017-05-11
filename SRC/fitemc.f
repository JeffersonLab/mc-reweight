
      REAL*8 FUNCTION FITEMC(X,A,GOODFIT)                                

!---------------------------------------------------------------------  
! Fit to EMC effect.  Steve Rock 8/3/94                                 
! Funciton returns value of sigma(A)/sigma(d) for isoscaler nucleus     
!  with A/2 protons=neutrons.                                           
! A= atomic number                                                      
! x = Bjorken x.                                                        
!                                                                       
! Fit of sigma(A)/sigma(d) to form C*A**alpha where A is atomic number  
! First data at each x was fit to form C*A**alpha.  The point A=2(d)    
!  was included with a value of 1 and an error of about 2%.             
! For x>=.125 Javier Gomez fit of 7/93 to E139 data was used.           
! For .09 >=x>=.0085 NMC data from Amaudruz et al Z. Phys C. 51,387(91) 
!  Steve did the fit for alpha and C to the He/d. C/d and Ca/d NMC data.
! Alpha(x) was fit to a 9 term polynomial a0 +a1*x +a2*x**2 + a3*x**3 ..
! C(x) was fit to a 3 term polynomial in natural logs as                
!  Ln(C) = c0 + c1*Ln(x) + c2*[Ln(x)]**2.                               
!-----------------------------------------------------------------------
                                                                        
                                                                        
      IMPLICIT NONE                                                     
      INTEGER I                                                         
      REAL*8 ALPHA, C,LN_C,X,A,A_OVER_D
      LOGICAL GOODFIT
                                                                        
!Chisq=         19.   for 30 points                                     
!Term    Coeficient     Error                                           
*
* N.I. modified to have arrays indices starting from 1 (not 0)
      REAL*8 ALPHA_COEF(2,9)                                
!Chisq=         22.    for 30 points                                   
!Term    Coeficient     Error                                          
      REAL*8 C_COEF(2,3)         ! Value and error for 6 term fit to 
                              
*
* N.I. moved things around for f2c compat.
*
      DATA C_COEF /        ! Value and error for 6 term fit to 
     >  1.69029097D-02,    4.137D-03,                                   
     >  1.80889367D-02,    5.808D-03,                                   
     >  5.04268396D-03,    1.406D-03   /    
      DATA ALPHA_COEF    /                                     
     > -6.98871401D-02,    6.965E-03,                                   
     >  2.18888887D+00,    3.792E-01,                                   
     > -2.46673765D+01,    6.302E+00,                                   
     >  1.45290967D+02,    4.763E+01,                                   
     > -4.97236711D+02,    1.920E+02,                                   
     >  1.01312929D+03,    4.401E+02,                                   
     > -1.20839250D+03,    5.753E+02,                                   
     >  7.75766802D+02,    3.991E+02,                                   
     > -2.05872410D+02,    1.140E+02 /     
*
                                                                        
      IF( (X.GT.0.88).OR.(X.LT. 0.0085) ) THEN   !Out of range of fit   
       X =.0085
       GOODFIT=.FALSE.
      ELSE
       GOODFIT=.TRUE.
      ENDIF                                                            
                                                                        
      LN_C = C_COEF(1,0+1)                
      DO I =1,2                                                         
       LN_C = LN_C + C_COEF(1,I+1) * (DLOG(X))**I                         
      ENDDO                                                             
                                                                        
      C = EXP(LN_C)                                                     
                                                                        
      ALPHA = ALPHA_COEF(1,0+1)                                           
      DO I=1,8                                                          
       ALPHA = ALPHA + ALPHA_COEF(1,I+1) * X**I                           
      ENDDO                                                             
                                                                        
      FITEMC  =  C *A**ALPHA                                            
      RETURN                                                            
      END    
