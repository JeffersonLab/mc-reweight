      subroutine rcext(firstr,theta,thspect,xin,tarid,rce)
      implicit none
      INCLUDE 'rad.cmn'
      character*80 infile
      integer*4 i,j,endof
      real*4 xin,rce,theta,thspect,thcentdeg
      real*8 rc(15),rc_cent
      real*8 radtab_temp(12),thetadeg,thcent,thetalow,thetahigh
      real*8 thetatab
      integer*4 tarnum,tarid,tar(20)
      real*8 mp,mp2,radcon,thetarad
      real*8 xtab,xtab_next
      logical firstr

      infile = 'rade.dat'
      if(firstr) open(unit=34,file=infile,status='old')    

      radcon = 180./3.141593
      thetadeg = theta*radcon
      thcentdeg = thspect*radcon
   
      do i=1,15
       rc(i) = 0.
      enddo


CCCCCC       The following are the column #'s for a given target    CCCCCC

      tar(1) = 12             !  Au           !
      tar(2) = 10             !  Cu2%         !
      tar(3) = 11             !  Cu3%         !
      tar(4) = 6              !  C1%          !
      tar(5) = 7              !  C3%          !
      tar(11) = 4             !  H2           !
      tar(15) = 5             !  D2           !
      tar(17) = 8             !  Al           !
      tar(20) = 9             !  cryo endcaps !

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  
      tarnum = tar(tarid)

CCCCCC              read in external radcor table              CCCCCC

      if (firstr) then 
       i = 1
       eentries = 0 
       endof = 1
       dowhile(endof.NE.-1)
        read(34,*,END=1001) radtab_temp
        do j=1,12
         exttab(i,j) = radtab_temp(j)
        enddo  
        eentries = eentries + 1
        i = i + 1 
       enddo
       write(6,*) "Nentries in External correction table is:  ",eentries
      endif
      close(34) 

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


CCCCCC         Calculate external correction                      CCCCCC
CCCCCC         by doing interpolation in theta and xin            CCCCCC

       thetalow = int(thetadeg)
       thetahigh = thetalow+1

CCCCCC     do search for rcs to interpolate in theta and xin.     CCCCCC 
CCCCCC     thetahigh is the integer theta above the               CCCCCC
CCCCCC     central theta.                                         CCCCCC
 
       do j=1,eentries
        thetatab = exttab(j,3)
        xtab = exttab(j,2)
        xtab_next = exttab(j+1,2)
        if(thetatab.eq.thetalow) then
         if(xin.GE.xtab_next.AND.xin.LE.xtab) then
          rc(1) = exttab(j,tarnum)
          rc(2) = exttab(j+1,tarnum)
          rc(3) = ((xin-xtab_next)*rc(1)+(xtab-xin)*rc(2))/
     &          (xtab-xtab_next)

         endif
        endif
        if(thetatab.eq.thetahigh) then
         if(xin.GE.xtab_next.AND.xin.LE.xtab) then
          rc(4) = exttab(j,tarnum)
          rc(5) = exttab(j+1,tarnum) 
          rc(6) = ((xin-xtab_next)*rc(4)+(xtab-xin)*rc(5))/
     &          (xtab-xtab_next)                      
         endif
        endif   

       enddo


CCCCCC                          End search                            CCCCCC

   
CCCCCC             Now do interpolation in theta                      CCCCCC

       rce = ((thetadeg-thetahigh)*rc(3) + 
     &       (thetalow-thetadeg)*rc(6))/(thetalow-thetahigh)  
     
c      write(6,*) xin,thetadeg,rc(3),rc(6),rce
       

       if(rce.LE.0) rce = 1.
 
 8000 format(a80) 
 1001 endof = -1
      return

      end





















