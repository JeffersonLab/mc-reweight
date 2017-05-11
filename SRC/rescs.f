CCC  Version 050521  -  Author:  M.E. Christy                        CCC
CCC  Subroutine to get Transvese and Longitudinal eP cross sections  CCC 
CCC  from fits to L/T cross sections.  The subroutine resmod.f is    CCC
CCC  required.  Units are in ub/Sr/Gev.                              CCC


      SUBROUTINE rescs(W2,Q2,sigt,sigl)

      IMPLICIT NONE

      real*8 w2,q2,fn,fnerr,xval1(50),xvall(50),w2inc,temp(4)
      real*8 mp,mp2,pi,alpha,xb,sigt,sigl
      integer i,npts,sf
 
      mp = .9382727
      mp2 = mp*mp
      pi = 3.141593
      alpha = 1./137.036

      open(unit=15,file='f1parms.dat',status='old')
      do i=1,50
        read(15,*) temp(1)  ! par #
        read(15,*) xval1(i) ! starting value
        read(15,*) temp(2)   ! initial step (0 means fixed parm)
        read(15,*) temp(3) ! low limit
        read(15,*) temp(4) ! high limit 
      enddo
      close(15)

      open(unit=16,file='flparms.dat',status='old')
      do i=1,50
        read(16,*) temp(1)  ! par #
        read(16,*) xvall(i) ! starting value
        read(16,*) temp(2)   ! initial step (0 means fixed parm)
        read(16,*) temp(3) ! low limit
        read(16,*) temp(4) ! high limit 
      enddo
      close(16)

      call resmod(1,w2,q2,xval1,sigt)
      call resmod(2,w2,q2,xvall,sigl)

    
      end

       




