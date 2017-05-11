      subroutine samp_eloss(Z,A,mass,thick,ex,pin,emean,de)       

CCC   samples eloss for charged particles passing through materials.     CCC
CCC   thick is in g/cm^2, ex = excitation potential for material ,       CCC 
CCC   pin = initial particle momentum in MeV, de = energy loss in MeV    CCC
CCC   The average energy loss (emean) and the particular energy loss     CCC
CCC   (de) is returned.                                                  CCC

      implicit none

      real*4 A,Z,mass,thick,ex,pin,de,ran3,vavrnd,gausin,ranlan
      real*4 rkappa,beta2
      real*4 x,f,ftest,f1,f2
      real*4 vviden,xl,xu,p,etot,lambda,ein
      real*4 mp,me,q,ex2,r,beta,gamma,gamma2
      real*4 del_ba,w,wmax,h1,h2,h3,Emean
      integer seed,mode

      de = 0.

      seed = 124571

      mp = 938.2727
      me = 0.511

      p = pin
      Etot = sqrt(p*p+mass*mass)     !!!  Calculate total energy  !!!
      ein = Etot-mass

      q = 1.                         !!!  Charge of incoming particle  !!!   
      ex2 = ex*ex

      r = (Z/A)

      beta = p/Etot
      beta2 = beta*beta
      gamma = Etot/mass
      gamma2 = gamma*gamma

      del_ba = (0.1534*r*thick*q*q/beta2)

      w = 1.+(2.*gamma*me/mass)+(me/mass)**2
      wmax = (2.*me*beta2*gamma2)/w

      h1 = (2.*me*beta2*gamma2*wmax)/ex2
      h2 = log(h1)
      h3 = h2-(2.*beta2)

      Emean = (del_ba*h3)       !!!  Mean energy loss      !!!
      rkappa = (del_ba/wmax)  

      mode = 1
      
      f = 0.

      x = ran3(seed)            !!!  random lambda value   !!!
      if (rkappa.GE.0.001.AND.rkappa.LE.12.) then  
        call VAVSET(rkappa,beta2,1)
        f = vavrnd(x)           !!!  sample Vavilov        !!!  
        de = (del_ba*((f/rkappa)+0.422784+beta2))+Emean
      elseif(rkappa.GT.12.) then
        f = gausin(x)           !!!  sample Gaussian       !!! 
        f1 = sqrt(1.-beta2/2.)
        f2 = del_ba*f1/sqrt(rkappa)
        de = (f*f2)+Emean
      else
        f = ranlan(x)           !!!  sample Landau         !!!
        de = del_ba*(f+0.422784+beta2 + log(del_ba/wmax))+Emean          
      endif

      if(de.GT.100) de = 0.     !!!  energy loss cutoff    !!! 
                                !!!  to keep things stable !!!
      return
      end
