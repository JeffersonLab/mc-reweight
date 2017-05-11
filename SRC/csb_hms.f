**         program poz
CCCC   Returns positron cross section in ub/sr/GeV/A  CCCC

	SUBROUTINE csb_hms(PYCH,EP,YG,TARG,CRSPOZ)
        IMPLICIT NONE

        integer targ        	
        REAL*4 PYCH,EP,YG,crspoz
	REAL*4 a1,a2,a3,a4
	REAL*4 b1,b2,b3,b4
	REAL*4 p1,p2  
	  

c        pych = 2.3
c        ep = 0.8
c        yg = 30.0
c        targ = 15

c        write(6,*) "CSB TEST  ",pych,ep,yg,targ

*********************************************
*******PARAMETERS FOR CARBON ****************
*********************************************

	if(targ.eq.4) then
	
**** PARAMETSR FOR A of THETA DEPENDENCE ***
	
	 a1 =  9.6620
	 a2 = -0.34965
	 a3 =  0.67087E-02
	 a4 = -0.34478E-04

**** PARAMETSR FOR B of THETA DEPENDENCE ***

	 b1 = -7.3074
	 b2 =  0.14937
	 b3 = -0.37033E-02
	 b4 =  0.
 
        endif

*********************************************
*******PARAMETERS FOR ALUMINIUM *************
*********************************************

	if(targ.eq.17) then
	
**** PARAMETSR FOR A of THETA DEPENDENCE ***
	
	 a1 =  12.203
	 a2 = -0.51778
	 a3 =  0.10980E-01
	 a4 = -0.71606E-04

**** PARAMETSR FOR B of THETA DEPENDENCE ***
		    
	 b1 = -12.730
	 b2 =  0.60047
	 b3 = -0.15836E-01
	 b4 =  0.10480E-03
	         
        endif

*********************************************
*******PARAMETERS FOR DEUTERIUM *************
*********************************************

	if(targ.eq.15.or.targ.eq.11) then
	
**** PARAMETSR FOR A of THETA DEPENDENCE ***
	
	 a1 =  4.8491
	 a2 = -0.33319E-01
	 a3 =  0.55768E-03
	 a4 =  0.43711E-05

**** PARAMETSR FOR B of THETA DEPENDENCE ***
		    
	 b1 = -1.1248
	 b2 = -0.20751
	 b3 =  0.27740E-02
	 b4 = -0.40153E-04
	         
        endif

*********************************************
********** PARAMETERS FOR FERUM *************
*********************************************
 
	if(targ.eq.3) then
	
**** PARAMETSR FOR A of THETA DEPENDENCE ***
	
	 a1 =  9.7137
	 a2 = -0.29383
	 a3 =  0.34126E-02
	 a4 = -0.28555E-05

**** PARAMETSR FOR B of THETA DEPENDENCE ***
 		    
	 b1 = -7.1223
	 b2 =  0.28352E-02
	 b3 =  0.28599E-02
	 b4 = -0.56969E-04
	         
        endif


*########################################################################
*############### HMS 2.3 Beam ENERGY ############  CARBON ###############
*########################################################################


	if(pych.gt.2.and.pych.lt.3) then

	  p1=a1+(a2*yg)+(a3*yg*yg)+(a4*yg*yg*yg)
	  p2=b1+(b2*yg)+(b3*yg*yg)+(b4*yg*yg*yg)

	  crspoz=exp(p1+p2*ep)	 
	  crspoz=crspoz/1000
          if(ep.gt.2.0) crspoz=0.   

	endif   
	   

c        write(6,*) "CSB TEST  ",pych,ep,yg,targ,crspoz


       END 
