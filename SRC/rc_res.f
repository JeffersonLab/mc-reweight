      subroutine rc_res(firstr,theta,thspect,xin,tarid,rce)
      implicit none
      INCLUDE 'rad.cmn'
      character*80 infile
      integer*4 i,j
      real*4 xin,rce,theta,thspect,thcentdeg
      real*8 rc(15),rc_cent
      real*8 radtab_temp(6),thetadeg,thcent,thetalow,thetahigh
      real*8 thetatab
      integer*4 tarnum,tarid,tar(20)
      real*8 mp,mp2,radcon,thetarad
      real*8 xtab,xtab_next
      logical firstr,endof,extrap1,extrap2

      infile = 'rc94.dat'

      if(firstr) open(unit=34,file=infile,status='old')    

      extrap1 = .false.           !!!  initialize to false       !!!
      extrap2 = .false.           !!!  if true then extrapolate  !!!

      radcon = 180./3.141593
      thetadeg = theta*radcon
      thcentdeg = thspect*radcon

   
      do i=1,15
       rc(i) = 0.
      enddo

CCCCCC       The following are the column #'s for a given target    CCCCCC


      tar(1) = 5             !  H2, or solid targets  !
      tar(2) = 5             !  Al                    !
      tar(3) = 5             !  cryo endcaps          !

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  
      tarnum = tar(tarid)   

CCCCCC              read in radcor table              CCCCCC

c      write(6,*)"here",tarid,tar(tarid),firstr


      if (firstr) then 
       i = 1
       eentries = 0 
       endof = .false.
       dowhile(.not.endof)
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

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


CCCCCC         Calculate radiative correction by doing            CCCCCC
CCCCCC         linear interpolation in theta and xin              CCCCCC

       thetalow = int(thetadeg)           !!! find integer angle below !!! 
       thetahigh = thetalow+1             !!! find integer angle above !!!

CCCCCC     do search for rcs to interpolate in theta and xin.     CCCCCC 
CCCCCC     thetahigh is the integer theta above the               CCCCCC
CCCCCC     central theta.                                         CCCCCC
 
       do j=1,eentries
        thetatab = exttab(j,3)
        xtab = exttab(j,2)
        xtab_next = exttab(j+1,2)  
        if(thetatab.eq.thetalow) then 
         if(xin.GE.xtab) then
          if(exttab(j-1,3).NE.thetatab) then  !!!  extrapolate  !!!
           extrap1 = .true.          
           rc(1) = exttab(j,tarnum)
           rc(2) = exttab(j+1,tarnum)   
           rc(3) = (rc(2)-rc(1))/(xtab_next-xtab)*(xin-xtab)+rc(1)

c           write(6,*) thetalow,xin,xtab,xtab_next,rc(1),rc(2),rc(3)

          endif
         endif
         if(xin.LE.xtab.and.xin.GE.xtab_next) then !!!  interpolate  !!!
           rc(1) = exttab(j,tarnum)
           rc(2) = exttab(j+1,tarnum)
           rc(3) = ((xin-xtab_next)*rc(1)+(xtab-xin)*rc(2))/
     &          (xtab-xtab_next) 
 
c           write(6,*) xin,xtab,xtab_next,rc(1),rc(2),rc(3)

         endif
       
c         write(6,*) thetalow,xin,rc(1),rc(2),rc(3)    

        endif

        if(thetatab.eq.thetahigh) then

         if(xin.GE.xtab) then
          if(exttab(j-1,3).NE.thetatab) then  !!!  extrapolate  !!!
           extrap2 = .true.
           rc(4) = exttab(j,tarnum)
           rc(5) = exttab(j+1,tarnum)
           rc(6) = (rc(5)-rc(4))/(xtab_next-xtab)*(xin-xtab)+rc(4)

c           write(6,*) thetahigh,xin,xtab,xtab_next,rc(4),rc(5),rc(6)
          endif
         endif

         if(xin.LE.xtab.and.xin.GE.xtab_next) then !!!  interpolate  !!!
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
   
c       write(6,*) xin,thetalow,rc(1),rc(2),rc(3)  

c       write(6,*) xin,thetadeg,rc 
   
       if(rce.LE.0) rce = -1.
 
 8000 format(a80) 

      return

      end





















