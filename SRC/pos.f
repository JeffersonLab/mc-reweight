       subroutine pos(e,hsp,thetac,posw)
       
ccccc  give charge symetric ratios for E94-110 at E=3.12 GeV  ccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  
       implicit none

       real*4 e,hsp,thetac,tt,posw,post(10),t(10)
       real*4 tlow,thigh,dt
       integer i,j,nangles,index
       logical extraplow,extraphigh

c       write(6,*) e,hsp,thetac
       
c       read(5,*) e,hsp,thetac
       j = 1
       do i=1,10
         t(i) = 0.
         post(i) = 0.
       enddo    

       if(abs(e-3.1183).lt.0.4) then

         post(1)=1.0-1.0/exp(2.581+1.1076*hsp
     &     +0.3212*hsp**2-0.0923*hsp**3+0.1478*hsp**4)


         post(2)=1.0-1.0/exp(-0.1792+4.3226*hsp
     &     +2.7167*hsp**2-4.1914*hsp**3+1.4280*hsp**4)
         post(3)=1.0-1.0/exp(-6.4848+11.2818*hsp
     &     +6.7417*hsp**2-10.8566*hsp**3+3.2703*hsp**4)
         post(4)=1.0-1.0/exp(-8.6611+19.3912*hsp
     &     +9.6906*hsp**2-30.5748*hsp**3+13.8177*hsp**4)
         post(5)=1.0-1.0/exp(-8.1676+25.7517*hsp
     &     +1.6765*hsp**2-44.3315*hsp**3+29.9229*hsp**4)
         post(6)=1.0-1.0/exp(-13.6448+50.9522*hsp
     &     +20.1701*hsp**2-235.312*hsp**3+232.311*hsp**4)

         nangles = 6
         t(1) = 23.0
         t(2) = 33.0
         t(3) = 41.0
         t(4) = 50.0
         t(5) = 62.0
         t(6) = 78.0
       endif

       extraplow = .false.
       extraphigh = .false.
    
       if(thetac.LT.t(1)) extraplow = .true.
       if(thetac.GT.t(nangles)) extraphigh = .true.
      
       if(.not.extraplow.AND..not.extraphigh) then
         do i=1,nangles
           if(thetac.GT.t(i)) then
             j = i        
           endif
         enddo
       endif 

       tlow = t(j)
       thigh = t(j+1)
       dt = thigh - tlow

c       write(6,*) j

       posw =  post(j)*(thigh-thetac)/dt 
     &           +post(j+1)*(thetac-tlow)/dt
 
c       write(6,*) e,hsp,thetac,extraplow,extraphigh,tlow,thigh
    

c       write(6,*) e,hsp,thetac,posw


       if(abs(e-3.5).GT.0.4) posw = 1.0
    
       return  
       end



