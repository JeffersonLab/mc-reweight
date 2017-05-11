      subroutine rc_mod(firstr,theta,thspect,xin,tarid,rce,sig)
      implicit none
      INCLUDE 'rad.cmn'
      character*80 infile
      integer*4 i,j,k
      real*4 xin,rce,theta,thspect,thcentdeg,sig
      real*8 x(6),rc(6),m(6)
      real*8 radtab_temp(5),thetadeg,thcent,thetalow,thetahigh
      real*8 thetatab,diffxL1,diffxL2,diffxH1,diffxH2
      integer*4 tarnum,tarid,tar(20),eof,tdiff,tdiff_min,tdiff_max
      real*8 mp,mp2,radcon,thetarad
      real*8 xtab,xtab_next,xtab_pre
      logical firstr,endof,extrap1,extrap2,extrap_x_hi,extrap_x_lo

      infile = 'rad_corr.dat'

      if(firstr) open(unit=34,file=infile,status='old')    

      radcon = 180./3.141593
      thetadeg = theta*radcon       !!! convert rad to deg !!!
      thcentdeg = thspect*radcon    !!! convert rad to deg !!!

      do i=1,6
        rc(i) = 0.
        m(i) = 0.
      enddo

CCCCCC              read in radcor table              CCCCCC

c      write(6,*)"here",tarid,tar(tarid),firstr

      if (firstr) then 
       i = 1
       eentries = 0 
       endof = .false.
       do while(.not.endof)
         read(34,*,END=1001) radtab_temp
         do j=1,5
           exttab(i,j) = radtab_temp(j)
         enddo 
         eentries = eentries + 1
         i = i + 1 
       enddo
       write(6,*) "Nentries in radcor table is:  ",eentries
      endif

 1001 endof = .true.

      close(34) 

c      write(6,*) "Nentries in radcor table is:  ",eentries

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

    
CCCCCC       Calculate radiative correction and model by doing    CCCCCC
CCCCCC       linear interpolation in theta and xin                CCCCCC

        tdiff_min = 10.0
        tdiff_max = 10.0

        do j=1,eentries
          thetatab = exttab(j,3)
          tdiff = abs(thetatab - thetadeg)
          if(thetatab.LT.thetadeg) then
            if(tdiff.LT.tdiff_min) then 
              tdiff_min = tdiff 
              thetalow = thetatab
            endif
          else
            if(tdiff.LT.tdiff_max) then 
              tdiff_max = tdiff
              thetahigh = thetatab
            endif 
          endif
        enddo

c        thetalow = int(thetadeg)-tdiff_min     !!! find integer angle below !!! 
c        thetahigh = int(thetadeg)+tdiff_max    !!! find integer angle above !!!



CCCCCC     do search for rcs to interpolate in theta and xin.     CCCCCC 
CCCCCC     thetahigh is the integer theta above the               CCCCCC
CCCCCC     central theta.                                         CCCCCC
 

        extrap_x_lo = .true.
        extrap_x_hi = .true.
        diffxh1 = 1000.
        diffxL1  = 1000.
        diffxh2 = 1000.
        diffxL2  = 1000.

        do j=1,eentries
          thetatab = exttab(j,3)

          if(int(thetatab).EQ.thetalow) then
c            extrap_th_lo = .false.
            if(exttab(j,2).LT.xin) then
              if(abs(exttab(j,2)-xin).LE.diffxL1) then
                diffxL1 = abs(exttab(j,2)-xin)
                x(1) = exttab(j,2)
                m(1) = exttab(j,4)
                rc(1) = exttab(j,5)
              endif
            elseif(exttab(j,2).GT.xin) then
              if(abs(exttab(j,2)-xin).LE.diffxH1) then       
                diffxH1 = abs(exttab(j,2)-xin)       
                x(2) = exttab(j,2)
                m(2) = exttab(j,4)
                rc(2) = exttab(j,5)
              endif
            endif 
          endif

          if(int(thetatab).EQ.thetahigh) then
c            extrap_th_hi = .false.
            if(exttab(j,2).LT.xin) then
              if(abs(exttab(j,2)-xin).LE.diffxL2) then
                diffxL2 = abs(exttab(j,2)-xin)
                x(3) = exttab(j,2)
                m(3) = exttab(j,4)
                rc(3) = exttab(j,5)
              endif
            elseif(exttab(j,2).GT.xin) then
              if(abs(exttab(j,2)-xin).LE.diffxH2) then       
                diffxH2 = abs(exttab(j,2)-xin)       
                x(4) = exttab(j,2)
                m(4) = exttab(j,4)
                rc(4) = exttab(j,5)
              endif
            endif 
          endif
     
        enddo 

        m(5) = (m(2)*(xin-x(1))+m(1)*(x(2)-xin))/(x(2)-x(1))
        rc(5) = (rc(2)*(xin-x(1))+rc(1)*(x(2)-xin))/(x(2)-x(1))
        m(6) = (m(4)*(xin-x(3))+m(3)*(x(4)-xin))/(x(4)-x(3))
        rc(6) =(rc(4)*(xin-x(3))+rc(3)*(x(4)-xin))/(x(4)-x(3))


c        write(6,*)"here:  ",m(3),m(4),m(6)


        sig = m(6)*(thetadeg-thetalow)+m(5)*(thetahigh-thetadeg)
        sig = sig/(thetahigh-thetalow)      
        rce = rc(6)*(thetadeg-thetalow)+rc(5)*(thetahigh-thetadeg)
        rce = rce/(thetahigh-thetalow)



CCCCCC                          End search                            CCCCCC

   
CCCCCC             Now do interpolation in theta                      CCCCCC


c        sig = (m(3)-m(6))*(thetadeg-thetahigh)/
c     &        (thetalow-thetahigh)+rc(6)
c        rce = (rcr(3)-rcr(6))*(thetadeg-thetahigh)/
c     &        (thetalow-thetahigh)+rcr(6)

c        if(rce.LE.0) rce = 0.000001
c        if(sig.LE.0) sig = 0.000001

c        write(6,*) thetadeg,xin,sig,rce

 8000 format(a80) 

      return

      end





















