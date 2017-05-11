      subroutine chsymfit(q2,theta,tarid,posweight)
      implicit none
      integer*4 i,tarid
      real*4 q2,theta,posweight
      real*8 jinvar,p1(20),p2(20),p3(20),p4(20)


CCCCCC       The following are the id #'s for a given target    CCCCCC
C
C      tar(1)            !  Au           !
C      tar(2)            !  Cu2%         !
C      tar(3)            !  Cu3%         !
C      tar(4)            !  C1%          !
C      tar(5)            !  C3%          ! 
C      tar(11)           !  H2           !
C      tar(15)           !  D2           !
C      tar(20)           !  cryo endcaps !
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      do i=1,20
       p1(i)=1.
       p2(i) = 1.
       p3(i) = 1.
       p4(i) = 1.
      enddo

      p1(2) = -0.2069e-04
      p2(2) = -309.9
      p3(2) = -0.3319e-1
      p4(2) = 32.38
  
      p1(3) = -0.1228e-3
      p2(3) = -331.3
      p3(3) = 1.499
      p4(3) = 41.51

      p1(4) = -0.1482e-5
      p2(4) = -309.9
      p3(4) = -0.9096
      p4(4) = 25.76
      
      p1(5) = 1
      p2(5) = 1
      p3(5) = 1
      p4(5) = 1

      p1(11) = -0.1406e-4
      p2(11) = -283.2
      p3(11) = -0.1612
      p4(11) = 28.22

      p1(15) = -0.1077e-4
      p2(15) = -295.1
      p3(15) = -0.4896e-1
      p4(15) = 28.14

c      write(6,*) tarid

      jinvar = q2/sin(theta)-theta/0.12+3.
      posweight = 1. + p1(tarid)*exp(-p2(tarid)/(jinvar*jinvar+
     &        p3(tarid)*jinvar+p4(tarid)))      
    
c      write(6,*) q2,theta,sin(theta),jinvar,posweight
      return

      end





















