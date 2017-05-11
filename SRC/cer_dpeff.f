1;3804;0cc       Program yptest

CCCCC     xfp in cm, ypfp and ypcor in mrad,   CCCCC

	Subroutine yp_optcor(xfp,ypfp,ypcor)
	implicit none

* Original version made on 01/07/06 by E. Christy
* to correct the MC_HMS Recon MEs (hms_cosey_refit_1.57.dat).  
* Interpolate in Xfp.

	real*4 xfp,ypfp,ypcor,coef(4,5),xfpb(5),xfplo,xfphi,dxfp
        real*4 ypcorlo,ypcorhi,slope
        integer i,j
        logical llo,hhi

        data eff/


        data xfpb/-38.,-20.,0.,20.,38./  

        llo = .false.
        hhi = .true. 
        xfphi = 1000.
        i = 1   

        if(abs(xfp).GT.50..OR.abs(ypfp).GT.50.) then
         ypcor = 0.0
         return
        endif  


c        write(6,*) "Enter xfp,ypfp"
c        read(5,*) xfp,ypfp

        if(xfp.LT.xfpb(1)) then 
          xfphi = xfpb(2)
          xfplo = xfpb(1)
          dxfp = xfp - xfpb(1)  
          i=2
          llo = .true.
          hhi = .false.
        elseif(xfp.GT.xfpb(5)) then
          xfphi = xfpb(5)
          xfplo = xfpb(4)
          dxfp = xfp - xfpb(5)
          i=5
          llo = .false.
          hhi = .true.
        else 
          dowhile(xfphi.GE.100)
            i = i+1
            if(xfpb(i).GT.xfp) then
              xfphi = xfpb(i)
              xfplo = xfpb(i-1)
              dxfp = xfp-xfplo
              hhi = .false.
            endif
          enddo
        endif

        ypcorlo = 0.0
        ypcorhi = 0.0

        do j=1,4
          ypcorlo = ypcorlo+coef(j,i-1)*ypfp**float(j-1)
          ypcorhi = ypcorhi+coef(j,i)*ypfp**float(j-1)
        enddo

        slope = (ypcorhi-ypcorlo)/(xfphi-xfplo)
        ypcor = ypcorlo + slope*dxfp
        if(hhi) ypcor = ypcorhi + slope*dxfp

c        write(6,*) xfp,ypfp,xfplo,xfphi,i,ypcorlo,ypcorhi,ypcor

	return
	end




