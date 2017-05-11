      subroutine rcint(firstr,theta,thspect,xin,tarid,rci)
      implicit none
      INCLUDE 'rad.cmn'
      character*80 infile,text
      integer*4 i,j,endof
      real*4 xin,rci,theta,thspect,thcentdeg
      real*8 rc(15),rc_cent
      real*8 radtab_temp(9),thetadeg,thcent,thetalow,thetahigh
      real*8 thetatab
      integer*4 tarnum,tarid,tar(20)
      real*8 mp,mp2,radcon,thetarad
      real*8 xtab,xtab_next
      logical firstr

      infile = 'radi.dat'
      if(firstr) open(unit=35,file=infile,status='old')    
c      read(35,8000) text

      radcon = 180./3.141593
      thetadeg = theta*radcon
      thcentdeg = thspect*radcon     

      do i=1,15
       rc(i) = 0.
      enddo


CCCCCC       The following are the column #'s for a given target    CCCCCC

      tar(1) = 9              !  Au           !
      tar(2) = 8              !  Cu2%         !
      tar(3) = 8              !  Cu3%         !
      tar(4) = 6              !  C1%          !
      tar(5) = 6              !  C3%          ! 
      tar(11) = 4             !  H2           !
      tar(15) = 5             !  D2           !
      tar(17) = 7             !  Al           !

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  
      tarnum = tar(tarid)

CCCCCC              read in internal radcor table              CCCCCC

      if (firstr) then
       i = 1
       ientries = 0 
       endof = 1
       dowhile(endof.NE.-1)
        read(35,*,END=1001) radtab_temp
        do j=1,9
         inttab(i,j) = radtab_temp(j)
        enddo 
        ientries = ientries + 1 
        i = i + 1
       enddo
       write(6,*) "Nentries in Internal correction table is:  ",ientries
      endif
      close(35) 


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


CCCCCC         Calculate internal correction                      CCCCCC
CCCCCC         by doing interpolation in theta and xin            CCCCCC

       thetalow = int(thetadeg)
       thetahigh = thetalow+1

CCCCCC     do search for rcs and interpolate in theta and xin.      CCCCCC 
CCCCCC     Currently the expansion is done about thetahigh,        CCCCCC
CCCCCC     Where thetahigh is the integer theta above the          CCCCCC
CCCCCC     central theta.                                          CCCCCC
 

       do j=1,ientries
        thetatab = inttab(j,3)
        xtab = inttab(j,2)
        xtab_next = inttab(j+1,2)       
        if(thetatab.eq.thetalow) then
         if(xin.GE.xtab_next.AND.xin.LE.xtab) then
          rc(1) = inttab(j,tarnum)
          rc(2) = inttab(j+1,tarnum) 
          rc(3) = ((xin-xtab_next)*rc(1)+(xtab-xin)*rc(2))/
     &          (xtab-xtab_next)                      
         endif
        endif
        if(thetatab.eq.thetahigh) then
         if(xin.GE.xtab_next.AND.xin.LE.xtab) then
          rc(4) = inttab(j,tarnum)
          rc(5) = inttab(j+1,tarnum)
          rc(6) = ((xin-xtab_next)*rc(4)+(xtab-xin)*rc(5))/
     &          (xtab-xtab_next)
         endif
        endif  
       enddo


CCCCCC                          End search                            CCCCCC


CCCCC             Now do interpolation in theta                      CCCCCC

       rci = ((thetadeg-thetahigh)*rc(3) +
     &       (thetalow-thetadeg)*rc(6))/(thetalow-thetahigh)

c       write(6,*) xin,thetadeg,rc(3),rc(6),rci


       rci = 1./(1.+rci/100.)
       

       if(rci.LE.0) rci = 1.

 8000 format(a80) 
 1001 endof = -1
      return

      end





















