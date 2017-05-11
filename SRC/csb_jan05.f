**         program csb_jan05
CCCC   Returns positron cross section in nb/sr/GeV/A  CCCC
	SUBROUTINE csb_jan05(beam,p,thave,target,poscs)
        IMPLICIT NONE
        
        INTEGER target	
        REAL*4 beam,p,thave,poscs,radcon
	REAL*4 a1,a2,a3,a4
	REAL*4 b1,b2,b3,b4
	REAL*4 p1,p2  
	 thave=thave*180.0/3.141592653

*********************************************
*******PARAMETERS FOR CARBON ****************
*********************************************

	if(target.eq.4) then
****Beam 1.2 GeV*****
          if(beam.lt.2) then 	

**** PARAMETERS FOR A of THETA DEPENDENCE ***

	 a1 = -1.6558 
	 a2 = -0.1646
	 a3 =  0
	 a4 =  0	
	 
**** PARAMETERS FOR B of THETA DEPENDENCE ***

	 b1 = 5.4913
	 b2 = 0.1365
	 b3 = 0
	 b4 = 0
	
	endif
	
	
****Beam 2.3 GeV*****      
           if(beam.gt.2.and.beam.lt.3) then
	
**** PARAMETERS FOR A of THETA DEPENDENCE ***

	 a1 = -8.1192
	 a2 =  0.0265
	 a3 = -0.0035
	 a4 =  0.0000

**** PARAMETERS FOR B of THETA DEPENDENCE ***

	 b1 =  6.9610
	 b2 = -0.0895
	 b3 =  0.0023
	 b4 =  0.0000
         endif
	 
****Beam 3.4 GeV*****      
           if(beam.gt.3.and.beam.lt.4) then	 
**** PARAMETERS FOR A of THETA DEPENDENCE ***
	
	 a1 =  26.9224
	 a2 = -1.8704
	 a3 =  0.0182
	 a4 =  0.0000	 	
	
**** PARAMETERS FOR B of THETA DEPENDENCE ***
	
	 b1 = -9.6744 
	 b2 =  0.7035
	 b3 = -0.0074
	 b4 =  0.0000		 	  
	 	
          endif

****Beam 4.6 GeV*****      
           if(beam.gt.4) then	 
**** PARAMETERS FOR A of THETA DEPENDENCE ***
	
	 a1 = -13.6523
	 a2 =  0.0000
	 a3 =  0.0000
	 a4 =  0.0000	 	
	
**** PARAMETERS FOR B of THETA DEPENDENCE ***
	
	 b1 =  3.8695
	 b2 =  0.0000
	 b3 =  0.0000 
	 b4 =  0.0000		 	  
	 	
          endif

	 endif

*********************************************
*******PARAMETERS FOR ALUMINUM **************
*********************************************

	if(target.eq.17) then
****Beam 1.2 GeV*****
          if(beam.lt.2) then 	

**** PARAMETERS FOR A of THETA DEPENDENCE ***

	 a1 = -1.7224
	 a2 = -0.1623
	 a3 =  0
	 a4 =  0	
	 
**** PARAMETERS FOR B of THETA DEPENDENCE ***

	 b1 = 6.6154
	 b2 = 0.1154
	 b3 = 0
	 b4 = 0
	
	endif
	
	
****Beam 2.3 GeV*****      
           if(beam.gt.2.and.beam.lt.3) then
	
**** PARAMETERS FOR A of THETA DEPENDENCE ***

	 a1 = -21.8214
	 a2 =  1.1190
	 a3 = -0.03116
	 a4 =  0.0002264

**** PARAMETERS FOR B of THETA DEPENDENCE ***

	 b1 =  15.3137
	 b2 = -0.7450
	 b3 =  0.01905
	 b4 = -0.0001376
         endif
	 
****Beam 3.4 GeV*****      
           if(beam.gt.3.and.beam.lt.4) then	 
**** PARAMETERS FOR A of THETA DEPENDENCE ***
	
	 a1 =  6.8910
	 a2 = -0.6439
	 a3 =  0.0000
	 a4 =  0.0000	 	
	
**** PARAMETERS FOR B of THETA DEPENDENCE ***
	
	 b1 = -1.2567 
	 b2 =  0.1994
	 b3 =  0.0000
	 b4 =  0.0000		 	  
	 	
          endif

****Beam 4.6 GeV*****      
           if(beam.gt.4) then	 
**** PARAMETERS FOR A of THETA DEPENDENCE ***
	
	 a1 = -13.7204
	 a2 =  0.0000
	 a3 =  0.0000
	 a4 =  0.0000	 	
	
**** PARAMETERS FOR B of THETA DEPENDENCE ***
	
	 b1 =  3.9694
	 b2 =  0.0000
	 b3 =  0.0000 
	 b4 =  0.0000		 	  
	 	
          endif

	 endif

*********************************************
*******PARAMETERS FOR IRON ******************
*********************************************

	if(target.eq.3) then
****Beam 1.2 GeV*****
          if(beam.lt.2) then 	

**** PARAMETERS FOR A of THETA DEPENDENCE ***

	 a1 = -1.2280 
	 a2 = -0.1755
	 a3 =  0
	 a4 =  0	
	 
**** PARAMETERS FOR B of THETA DEPENDENCE ***

	 b1 = 4.6206
	 b2 = 0.1555
	 b3 = 0
	 b4 = 0
	
	endif
	
	
****Beam 2.3 GeV*****      
           if(beam.gt.2.and.beam.lt.3) then
	
**** PARAMETERS FOR A of THETA DEPENDENCE ***

	 a1 = -6.4556
	 a2 = -0.0998
	 a3 = -0.0020
	 a4 =  0.0000

**** PARAMETERS FOR B of THETA DEPENDENCE ***

	 b1 =  5.7891
	 b2 = -0.0158
	 b3 =  0.0016
	 b4 =  0.0000
         endif
	 
****Beam 3.4 GeV*****      
           if(beam.gt.3.and.beam.lt.4) then	 
**** PARAMETERS FOR A of THETA DEPENDENCE ***
	
	 a1 =  1.8171
	 a2 = -0.4893
	 a3 =  0.0000
	 a4 =  0.0000	 	
	
**** PARAMETERS FOR B of THETA DEPENDENCE ***
	
	 b1 =  0.4638 
	 b2 =  0.1448
	 b3 =  0.0000
	 b4 =  0.0000		 	  
	 	
          endif

****Beam 4.6 GeV*****      
           if(beam.gt.4) then	 
**** PARAMETERS FOR A of THETA DEPENDENCE ***
	
	 a1 = -13.3629
	 a2 =  0.0000
	 a3 =  0.0000
	 a4 =  0.0000	 	
	
**** PARAMETERS FOR B of THETA DEPENDENCE ***
	
	 b1 =  3.7369
	 b2 =  0.0000
	 b3 =  0.0000 
	 b4 =  0.0000		 	  
	 	
          endif

	 endif

*********************************************
*******PARAMETERS FOR DEUTERIUM *************
*********************************************

	if(target.eq.15) then
****Beam 1.2 GeV*****
          if(beam.lt.2) then 	

**** PARAMETERS FOR A of THETA DEPENDENCE ***

	 a1 =  2.6517
	 a2 = -0.2514
	 a3 =  0
	 a4 =  0	
	 
**** PARAMETERS FOR B of THETA DEPENDENCE ***

	 b1 = -0.3717
	 b2 =  0.2546
	 b3 =  0
	 b4 =  0
	
	endif
	
	
****Beam 2.3 GeV*****      
           if(beam.gt.2.and.beam.lt.3) then
	
**** PARAMETERS FOR A of THETA DEPENDENCE ***

	 a1 = -0.2654
	 a2 = -0.3321
	 a3 =  0.0000
	 a4 =  0.0000

**** PARAMETERS FOR B of THETA DEPENDENCE ***

	 b1 =  2.0729
	 b2 =  0.1433
	 b3 =  0.0000
	 b4 =  0.0000
         endif
	 
****Beam 3.4 GeV*****      
           if(beam.gt.3.and.beam.lt.4) then	 
**** PARAMETERS FOR A of THETA DEPENDENCE ***
	
	 a1 =  3.6954
	 a2 = -0.5599
	 a3 =  0.0000
	 a4 =  0.0000	 	
	
**** PARAMETERS FOR B of THETA DEPENDENCE ***
	
	 b1 =  0.1669 
	 b2 =  0.1634
	 b3 =  0.0000
	 b4 =  0.0000		 	  
	 	
          endif

****Beam 4.6 GeV*****      
           if(beam.gt.4) then	 
**** PARAMETERS FOR A of THETA DEPENDENCE ***
	
	 a1 = -14.6394
	 a2 =  0.0000
	 a3 =  0.0000
	 a4 =  0.0000	 	
	
**** PARAMETERS FOR B of THETA DEPENDENCE ***
	
	 b1 =  4.0406
	 b2 =  0.0000
	 b3 =  0.0000 
	 b4 =  0.0000		 	  
	 	
          endif

	 endif

*********************************************
*******PARAMETERS FOR HYDROGEN **************
*********************************************

	if(target.eq.11) then
****Beam 1.2 GeV*****
          if(beam.lt.2) then 	

**** PARAMETERS FOR A of THETA DEPENDENCE ***

	 a1 = -6.3602 
	 a2 =  0.2518
	 a3 = -0.0060
	 a4 =  0	
	 
**** PARAMETERS FOR B of THETA DEPENDENCE ***

	 b1 =  12.4827
	 b2 = -0.4320 
	 b3 =  0.0082
	 b4 =  0
	
	endif
	
	
****Beam 2.3 GeV*****      
           if(beam.gt.2.and.beam.lt.3) then
	
**** PARAMETERS FOR A of THETA DEPENDENCE ***

	 a1 =  16.8560
	 a2 = -0.5786
	 a3 =  0.0000
	 a4 =  0.0000

**** PARAMETERS FOR B of THETA DEPENDENCE ***

	 b1 = -6.8712
	 b2 =  0.2732
	 b3 =  0.0000
	 b4 =  0.0000
         endif
	 
****Beam 3.4 GeV*****      
           if(beam.gt.3.and.beam.lt.4) then	 
**** PARAMETERS FOR A of THETA DEPENDENCE ***
	
	 a1 =  5.5181
	 a2 = -0.6242
	 a3 =  0.0000
	 a4 =  0.0000	 	
	
**** PARAMETERS FOR B of THETA DEPENDENCE ***
	
	 b1 = -0.4809 
	 b2 =  0.1899
	 b3 =  0.0000
	 b4 =  0.0000		 	  
	 	
          endif

****Beam 4.6 GeV*****      
           if(beam.gt.4) then	 
**** PARAMETERS FOR A of THETA DEPENDENCE ***
	
	 a1 =  0
	 a2 =  0.0000
	 a3 =  0.0000
	 a4 =  0.0000	 	
	
**** PARAMETERS FOR B of THETA DEPENDENCE ***
	
	 b1 =  0
	 b2 =  0.0000
	 b3 =  0.0000 
	 b4 =  0.0000		 	  
	 	
          endif

	 endif
*********************************************
*********************************************

	  p1=a1+(a2*thave)+(a3*thave*thave)+(a4*thave*thave*thave)
	  p2=b1+(b2*thave)+(b3*thave*thave)+(b4*thave*thave*thave)

	  poscs=(exp(p1))*(exp(p2*(beam-p))-1.)	 
c	  poscs=poscs/1000

cc Set Quasi-elastic contribution to ZERO.
	  
c	if(beam.lt.2) then	  
c          if(p.gt.1.0) poscs=0.
c	endif   	  
	  
c	if(beam.gt.2.and.beam.lt.3) then	  
c          if(p.gt.2.0) poscs=0.
c	endif   
	   
c	if(beam.gt.3.and.beam.lt.4) then	  
c          if(p.gt.3.0) poscs=0.
c	endif   
	
       END 
