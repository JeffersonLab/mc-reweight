      PROGRAM CALLMOD

      implicit none

      INCLUDE 'logicals.cmn'
      INCLUDE 'flags.cmn'
      INCLUDE 'rad.cmn'
        

**************************************************************************************
*             Calls model cross section for use in extract_cs.kumac                  *  
**************************************************************************************  


      real*4 e,ep,th,kin(5),sig,w2,q2,thc,t1,trad
      integer target
      
      open(unit=25,file='../output/data/cskin.dat',status='old')
      read(25,*) kin
      e = kin(1)
      ep = kin(2)
      th = kin(3)
      thc = th*3.14159/180.
      q2 = 4.*e*ep*sin(thc/2.)*sin(thc/2.)
      w2 = 0.9382727*0.9382727+2.*0.9382727*(e-ep)-q2      
      target = int(kin(4))
      
      th = th*3.141592654/180.

      call model_new(e,ep,th,target,sig)

CCCC   Now interpolate model from RC table  CCCC

c      trad = 3.14159/180.0*thc

      call rc_mod(.true.,thc,thc,ep,1,t1,sig)
      
      open(unit=26,file='../output/data/csmod.dat',status='new')
      write(26,*) W2,Q2,sig
c      write(6,*) e,ep,th,target,sig
    
      end







