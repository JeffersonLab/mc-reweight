      SUBROUTINE MODEL_NEW(epass,eppass,thpass,tar,SIG4)

      implicit none
       
      INCLUDE 'flags.cmn'

      real*8 E,EP,TH
      real*8 WW,QQ,x,xtemp,nu,a(20),z(20)
      real*8 W1,W2,MP,MP2,mott,q4,sig,f2
      real*8 sin2,cos2,tan2,alpha,alpha2,dr,qqtemp
      real*8 f1,fl
      real*4 epass,eppass,thpass,sig4
      integer tar
      sig = 0.

      e = epass
      ep = eppass
      th = thpass      

CCCCCC            Target to Amu conversion table                    CCCCCCC

      a(2) = 63.                   ! Cu  !!!  Fix, should be gold
      a(3) = 55.                   ! Fe
      a(4) = 12.                   ! C12
      a(5) = 12.                   ! C12
      a(11) = 1.                   ! H
      a(15) = 2.                   ! D
      a(17) = 27.                  ! Al  

      z(2) = 30.
      z(3) = 27.
      z(4) = 6.
      z(5) = 6.
      z(11) = 1.
      z(15) = 1.
      z(17) = 13.

      mp = .9382727
      mp2 = mp*mp
      alpha = 1/137.03599
      alpha2 = alpha*alpha  
      sin2 = sin(th/2.)*sin(th/2.)
      cos2 = cos(th/2.)*cos(th/2.)
      tan2 = tan(th/2.)*tan(th/2.)
      qq = 4.*e*ep*sin2
      q4 = qq*qq
      nu = e-ep
      ww = mp2 + 2.*mp*nu -qq
      x = qq/2./mp/nu

      mott = (4*alpha2*ep*ep*cos2/q4)*0.38938e6     !  In nb/GeV/sr  !   

      qqtemp = qq

      call f1f209(z(tar),a(tar),qq,ww,f1,f2)
c      f2 = 0.25
c      f1 = 0.25/x 
     
      sig = mott*(f2/nu + 2.*mp*f1*tan2)

      sig4 = sig


        
      end














