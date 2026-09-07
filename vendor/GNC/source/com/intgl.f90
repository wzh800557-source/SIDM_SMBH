!     IDID        REPORTS ON SUCCESSFULNESS UPON RETURN:
!                   IDID= 1  COMPUTATION SUCCESSFUL,
!                   IDID= 2  COMPUT. SUCCESSFUL (INTERRUPTED BY SOLOUT)
!                   IDID=-1  INPUT IS NOT CONSISTENT,
!                   IDID=-2  LARGER NMAX IS NEEDED,
!                   IDID=-3  STEP SIZE BECOMES TOO SMALL.
!                   IDID=-4  PROBLEM IS PROBABLY STIFF (INTERRUPTED).

module my_intgl
	implicit none
contains    
	subroutine my_integral_none(xs, xe, y, FCN,idid)
	!!  subroutine FCN(N,X,Y,F,RPAR,IPAR)
		IMPLICIT NONE
		integer ITOL,IDID
		integer n_y,ipar(100)
		real(8) y,yout(1)
		real(8) RTOL(1),ATOL(1),rpar(100)
		integer,parameter::LWORK=400,LIWORK=400
		real(8)	WORK(LWORK)
		integer IWORK(LIWORK)
		real(8) x_start,x_end,xs,xe
		integer IOUT
		real(8) hrlow,hrhigh
		external::dopri5,FCN,my_solout_empty

		RTOL=1d-10
		ATOL=1d-10
		ITOL=0;;WORK=0;IWORK=0
		!WORK(1)=1d-30;  
		!IWORK(3)=6;   ! print messages
		IWORK(1)=10000000;
		WORK(7)=0.
		IOUT=0
		x_start=xs
		x_end=xe
        idid=0
		if(x_start.eq.x_end) then
			y=0
			return
		end if
		yout(1)=y
		call dopri5(1,FCN,x_start,yout,x_end,RTOL,ATOL,ITOL,my_solout_empty,IOUT,WORK,LWORK,IWORK,LIWORK,RPAR,ipar,IDID)
		y=yout(1)
	end subroutine
	subroutine my_integral_acc(xs, xe, y,ATOL_in, RTOL_in, FCN,IDID)
	!!  subroutine FCN(N,X,Y,F,RPAR,IPAR)
		IMPLICIT NONE
		integer LWORK,ITOL,LIWORK,IDID
		integer n_y,ipar(100)
		real(8) y,yout(1)
		real(8) RTOL(1),ATOL(1),rpar(100),atol_in,rtol_in
		real(8)	WORK(400)
		real(8) x_start,x_end,xs,xe
		integer IOUT
		integer IWORK(150)
		real(8) hrlow,hrhigh
		external::dopri5,FCN,my_solout_empty

		ITOL=0;LWORK=400;WORK=0;IWORK=0;LIWORK=400
		!WORK(1)=1d-30;  
		!IWORK(3)=6;   ! print messages
		IWORK(1)=10000000;
		WORK(7)=0.
		IOUT=0

		x_start=xs
		x_end=xe
		rtol=rtol_in
		atol=atol_in
		yout(1)=y
		call dopri5(1,FCN,x_start,yout,x_end,RTOL,ATOL,ITOL,my_solout_empty,IOUT,WORK,LWORK,IWORK,LIWORK,RPAR,ipar,IDID)
		y=yout(1)
	end subroutine
	
end module

subroutine my_solout_empty(NR,XOLD,X,Y,N,CON,ICOMP,ND,RPAR,IPAR,IRTRN)
!     SOLOUT      NAME (EXTERNAL) OF SUBROUTINE PROVIDING THE
!                 NUMERICAL SOLUTION DURING INTEGRATION. 
!                 IF IOUT.GE.1, IT IS CALLED AFTER EVERY SUCCESSFUL STEP.
!                 SUPPLY A DUMMY SUBROUTINE IF IOUT=0. 
!                 IT MUST HAVE THE FORM
!                    SUBROUTINE SOLOUT (NR,XOLD,X,Y,N,CON,ICOMP,ND,
!                                       RPAR,IPAR,IRTRN)
!                    DIMENSION Y(N),CON(5*ND),ICOMP(ND)
!                   ....  
!                 SOLOUT FURNISHES THE SOLUTION "Y" AT THE NR-TH
!                    GRID-POINT "X" (THEREBY THE INITIAL VALUE IS
!                    THE FIRST GRID-POINT).
!                 "XOLD" IS THE PRECEEDING GRID-POINT.
!                 "IRTRN" SERVES TO INTERRUPT THE INTEGRATION. IF IRTRN
!                    IS SET <0, DOPRI5 WILL RETURN TO THE CALLING PROGRAM.
!                    IF THE NUMERICAL SOLUTION IS ALTERED IN SOLOUT,
!                    SET  IRTRN = 2
!           
!          -----  CONTINUOUS OUTPUT: -----
!                 DURING CALLS TO "SOLOUT", A CONTINUOUS SOLUTION
!                FOR THE INTERVAL [XOLD,X] IS AVAILABLE THROUGH
!                 THE FUNCTION
!                        >>>   CONTD5(I,S,CON,ICOMP,ND)   <<<
!                 WHICH PROVIDES AN APPROXIMATION TO THE I-TH
!                 COMPONENT OF THE SOLUTION AT THE POINT S. THE VALUE
!                 S SHOULD LIE IN THE INTERVAL [XOLD,X].
	implicit none
	INTEGER N, ND, NR
	real(8) X, XOLD
    real(8) Y(N),CON(5*ND),ICOMP(ND)
	integer IRTRN,IPAR(100)
	real(8) rpar(100)
	
end subroutine
!subroutine FCN(N,X,Y,F,RPAR,IPAR)
!integer N
!integer IPAR(100)
!real(8) rpar(100)
!real(8) X,Y(N),F(N)
!end subroutine
