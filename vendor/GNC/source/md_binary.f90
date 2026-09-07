module md_particle
!	use init_constant
!	use commonfunc
	use constant
	type::particle
        real(8) :: x(3)
        real(8) :: vx(3)
        real(8) :: M
        real(8) radius
 		integer id, obtype, obidx  
	end type
	!integer::PARTICLE_MPI_TYPE_
contains
	real(8) function get_distance(P1,P2)
		type(particle)::P1,P2
		get_distance=sqrt((P1%x(1)-P2%x(1))**2+(P1%x(2)-P2%x(2))**2+(P1%x(3)-P2%x(3))**2)
	end function
	real(8) function get_VelMag(P)
		type(particle)::P
		get_VelMag=sqrt(P%vx(1)**2+P%vx(2)**2+P%vx(3)**2)
	end function	
	real(8) function get_RMag(P)	
		type(particle)::P
		get_RMag=sqrt(P%x(1)**2+P%x(2)**2+P%x(3)**2)
	end function
end module

module md_binary
	use md_particle
	type::binary
		type(particle)::Ms
		type(particle)::Mm
		type(particle)::rd           !the relative position and velocity, from massive to light one
		real(8) E,l,k,miu,Mtot, Jc!, tl  !tl is the true longitude (tl=pe+f0)
		real(8) a_bin,e_bin, lum(3), f0
		real(8) Inc, Om, pe, t0, me ! Inc, Om, pe in unit of degree
		character*(100) bname
		integer an_in_mode
	!contains
	!	procedure::print=>print_binary
	!	procedure::get_period=>by_get_period
	!---
	!---be careful when e=0
	!---
	end type
	integer,parameter::an_in_mode_f0=1, an_in_mode_t0=2, an_in_mode_mean=3
	integer,parameter::nint_by=10, nreal_by=17+24, nstr_by=100
	!integer::BINARY_MPI_TYPE_
contains
	
	subroutine by_st2em(by,f_return)
        !!!------
		!by%t0 is not defined as the zero point is unkown
		!!!------
		IMPLICIT NONE
		type(binary) by
		integer f_return

		real(8) car(6), x,y,z,vx,vy,vz,a,e,Inc,Om, pe, tau ! ele(Ieles)
		real(8) gm, q, p,n,l,lm
		real(8) mco_kep
		f_return=0
		gm=by%Ms%m+by%Mm%m
		x=by%rd%x(1); vx=by%rd%vx(1)
		y=by%rd%x(2); vy=by%rd%vx(2)
		z=by%rd%x(3); vz=by%rd%vx(3)
!		print*, by%rd%x, by%rd%vx
		call mco_x2el(gm,x,y,z,vx,vy,vz, q,e,inc,p,n,l)
		if(e.gt.1d-6)then
			call get_true_anomaly(gm, x,y,z,vx,vy,vz,tau)
		else 
			tau=l
		end if
		om=n; pe=p-n;lm=l
		by%a_bin=q/(1-e); by%e_bin=e
		by%inc=inc; by%om=om; 
		by%pe=pe
		by%f0=tau; by%me=l

		if(isnan(tau))then
			print*, "st2em:nan!"
			call print_binary(by)
			f_return=1
!			stop
		end if
!		ele(8)=ele(4)+ele(5)+ele(7)      
		! t is the current time, n is the mean angular velocity 
!		ele(10)=mco_kep(e,l)*180d0/pi
	
	end subroutine
    subroutine get_center_particle_p(m1,m2,mc)
        type(particle)::m1,m2,mc
        real(8) mtot
        mtot=m1%m+m2%m
       	mc%x=(m1%x*m1%m+m2%x*m2%m)/mtot
		mc%vx=(m1%vx*m1%m+m2%vx*m2%m)/mtot
		mc%m=mtot
    end subroutine
	function get_center_particle(by)
		type(binary)::by
		type(particle) :: gcp, get_center_particle
		by%mtot=by%ms%m+by%mm%m
		gcp%x=(by%Ms%x*by%ms%m+by%Mm%x*by%mm%m)/by%Mtot
		gcp%vx=(by%Ms%vx*by%ms%m+by%Mm%vx*by%mm%m)/by%Mtot
		gcp%m=by%mtot
		get_center_particle=gcp
	end function

	subroutine by_move_to_mass_center(by) !!origin move to mass center
		type(binary) by
		type(particle) center_particle
 		center_particle=get_center_particle(by)
		by%Ms%x=by%Ms%x-center_particle%x
		by%Ms%vx=by%Ms%vx-center_particle%vx
		by%Mm%x=by%Mm%x-center_particle%x
		by%Mm%vx=by%Mm%vx-center_particle%vx
	end subroutine
    subroutine by_get_el(by,kind_in)
        implicit none
        type(binary)::by
        integer,optional::kind_in
        integer kind
        kind=0
        if(present(kind_in))kind=kind_in
        select case(kind)
        case(0)
            by%k=by%Ms%m*by%Mm%m
            by%Mtot=by%Ms%m+by%Mm%m
            by%miu=by%k/by%Mtot	
            by%e=-by%k/2/by%a_bin
            by%l=sqrt(by%miu*by%k*by%a_bin*(1-by%e_bin**2))
        case(1)
            !by%k=by%Mm%m
            by%Mtot=by%Ms%m+by%Mm%m
            !by%miu=by%k/by%Mtot	
            by%e=-by%mtot/2/by%a_bin
            by%l=sqrt(by%mtot*by%a_bin*(1-by%e_bin**2))
        end select
    end subroutine
	subroutine by_get_energy_lum(by)  !!get energy and l from position and velocity
		IMPLICIT NONE
		type(binary)::by
		type(particle)::bin_em,Ms,Mm
		real(8) l_m(3),get_l_m
		Ms=by%Ms;Mm=by%Mm	
		bin_em%x=Ms%x-Mm%x
		bin_em%vx=Ms%vx-Mm%vx
		by%e=-by%k/get_Rmag(bin_em)+by%miu*get_VelMag(bin_em)**2/2
		call by_get_lm(by)
		l_m=by%lum
		by%l=by%miu*sqrt(l_m(1)**2+l_m(2)**2+l_m(3)**2)
	end subroutine
	subroutine by_get_lm(by)
	IMPLICIT NONE
		type(binary)::by
		type(particle)::Ms,Mm
		real(8) get_l_m(3)
		Ms=by%Ms;Mm=by%Mm
		get_l_m(1)=sin(by%inc)*sin(by%om)
		get_l_m(2)=-sin(by%inc)*cos(by%om)
		get_l_m(3)=cos(by%inc)
		by%lum=get_l_m
	end subroutine

	subroutine by_get_Jc(by)
	IMPLICIT NONE
		type(binary)::by
		type(particle)::Ms,Mm
		real(8) get_l_m(3)
		by%Jc=((by%ms%m+by%mm%m)*by%a_bin*(1-by%e_bin**2))**0.5
	end subroutine

	real(8) function by_get_period(by)
		class(binary) ::by
		by_get_period=2*PI*by%a_bin*sqrt(by%a_bin/by%mtot)
	end function

	subroutine by_em2st(by)
		use constant
		implicit NONE
		type(binary)::by
		integer way_an_in
		real(8) gm, q,e,inc,p,n,l,x,y,z,u,v,w
		real(8) ecc_ano, r,theta,phi
		real(8) pos(3),vel(3), outvel(3)
		!if by%e_bin=0, we should avoid an_in_mode to be an_in_mode_f0

!  gm = grav const * (central + secondary mass)
!  q = perihelion distance
!  e = eccentricity
!  i = inclination                 )
!  p = longitude of perihelion !!! )   in
!  n = longitude of ascending node ) radians
!  l = mean anomaly       
!  x axis points to n=0
!  velocity in units of sqrt(gm/a) 

		gm=by%ms%m+by%mm%m;
		e=by%e_bin
		q=by%a_bin*(1-e)
		inc=by%inc
		n=by%om
		p=by%pe+n
		
		select case(by%an_in_mode)
			case (an_in_mode_t0)
				l=-sqrt(by%Mtot/by%a_bin**3)*by%t0
			case (an_in_mode_mean)
				l=by%me
				by%t0=-l/sqrt(by%Mtot/(by%a_bin**3))
			case (an_in_mode_f0)
				if(by%e_bin<1d0)then
					ecc_ano=atan(tan(by%f0/2d0)*(1-by%e_bin)**0.5/(1+by%e_bin)**0.5)*2
					if (by%f0.gt.pi-30d0/180d0*pi.and.by%f0.le.pi)then
						ecc_ano=acos((cos(by%f0)+by%e_bin)/(1+by%e_bin*cos(by%f0)))
					else if (by%f0.lt.pi+30d0/180d0*pi.and.by%f0.gt.pi)then
						ecc_ano=-acos((cos(by%f0)+by%e_bin)/(1+by%e_bin*cos(by%f0)))
					end if
					if(by%f0.eq.pi) ecc_ano= pi
					l=ecc_ano-by%e_bin*sin(ecc_ano)	
				else if(by%e_bin>1d0)then

					ecc_ano=atanh(tan(by%f0/2d0)*(by%e_bin-1)**0.5/(1+by%e_bin)**0.5)*2
					if (by%f0.gt.pi-30d0/180d0*pi.and.by%f0.lt.pi)then
						ecc_ano=acosh((cos(by%f0)+by%e_bin)/(1+by%e_bin*cos(by%f0)))
					else if (by%f0.lt.pi+30d0/180d0*pi.and.by%f0.gt.pi)then
						ecc_ano=-acosh((cos(by%f0)+by%e_bin)/(1+by%e_bin*cos(by%f0)))
					end if
					if(by%f0.eq.pi) ecc_ano=pi
					l=by%e_bin*sinh(ecc_ano)-ecc_ano	
				else
					ecc_ano=tan(by%f0/2d0)
					l=ecc_ano+ecc_ano**3/3d0
				end if
				by%t0=-l/sqrt(by%mtot/by%a_bin**3)
			case default
				print*, "ERROR: an_in_mode not defined!"
				stop
		end select
		by%me=l
		!print*, "l=",l
		call mco_el2x(gm,q,e,inc,p,n,l,x,y,z,u,v,w)

		if(by%an_in_mode.eq.an_in_mode_t0) then
			if(by%e_bin.gt.1d-6)then
				call get_true_anomaly(gm, x,y,z,u,v,w,by%f0)
			else
				by%f0=l
			end if
		end if
		!print*, x,y,z,u,v,w,by%f0,by%an_in_mode.eq.an_in_mode_t0
!		stop
		by%rd%x=(/x,y,z/);
		by%rd%vx=(/u,v,w/);
		!print*,"I, P, N, L=", inc, p ,n ,l 
		!print*, "x=",by%rd%x
		!stop
		by%rd%m=by%ms%m*by%mm%m/gm
!		by%mm%x=0;	by%mm%vx=0
		!by%rd%x=by%ms%x; by%rd%vx=by%ms%vx; 
		!print*, by%ms%x,by%ms%vx, by%mm%x,by%mm%vx
!		call by_move_to_mass_center(by)
		!print*, by%ms%x,by%ms%vx, by%mm%x,by%mm%vx

		if(isnan(by%f0).and.by%an_in_mode.eq.an_in_mode_f0)then
			print*, "em2st:nan!"
			call print_binary(by)
!			stop
		end if

	end subroutine

	subroutine by_split_from_rd(by)
		implicit none
		type(binary)::by
		
 		by%ms%x=by%rd%x*by%mm%m/(by%ms%m+by%mm%m)
		by%mm%x=-by%rd%x*by%ms%m/(by%ms%m+by%mm%m)
		by%ms%vx=by%rd%vx*by%mm%m/(by%ms%m+by%mm%m)
		by%mm%vx=-by%rd%vx*by%ms%m/(by%ms%m+by%mm%m)
		
	end subroutine
	subroutine by_get_rd(by)
		implicit none
		type(binary)::by
		by%rd%x=by%ms%x-by%mm%x
		by%rd%vx=by%ms%vx-by%mm%vx
	end subroutine
	real(8) function by_get_rdmag(by)
		implicit none
		type(binary)::by
		call vector_mag(by%rd%x, by_get_rdmag)
	end function
	subroutine vector_mag(v,mf)
		implicit none
        real(8) v(3),mf
		!if(v(1)**2+v(2)**2+v(3)**2.eq.0d0)then
		!	mf=0d0
		!	return
		!end if
        mf=sqrt(v(1)**2+v(2)**2+v(3)**2)
end subroutine
	subroutine print_binary(by,funit)
		implicit none
		class(binary)::by
        integer,optional::funit
        character*(10) str_ms_obtype, str_mm_obTYPe

        if(present(funit))then
            write(unit=funit,fmt="(2A20)") "NAME=", trim(adjustl(by%bname))
            write(unit=funit,fmt="(A20,I4)") "ANMODE=",by%an_in_mode
            write(unit=funit,fmt="(A20,1P20E20.8)") "MS,MM=",by%ms%m,by%mm%m
            call get_star_type(by%ms%obtype, str_ms_obtype)
            call get_star_type(by%mm%obtype, str_mm_obtype)
            write(unit=funit,fmt="(3A20)") "TYPE_MS,TYPE_MM=",str_ms_obtype,str_mm_obtype
            
            write(unit=funit,fmt="(A20,1P20E20.8)") "Sma,Ecc, Inc=",by%a_bin,by%e_bin,by%Inc
            write(unit=funit,fmt="(A20,1P20E20.8)") "E,L=",by%e,by%l
            if(by%a_bin<0d0)then
                write(*,fmt="(A20,1P20E20.8)") "Vinf(kms)=",sqrt(-(by%ms%m+by%mm%m)/by%a_bin)*29.79
            end if
            write(unit=funit,fmt="(A20,20F20.10)") "Ome,Pe,Me=",by%Om,by%pe,by%me
            write(unit=funit,fmt="(A20,20F20.10)") "f=",by%f0

            write(unit=funit,fmt="(A20,20F20.10)") "RD X=",by%rd%x
            write(unit=funit,fmt="(A20,20F20.10)") "RD MAG(X)=",sqrt(by%rd%x(1)**2+by%rd%x(2)**2+by%rd%x(3)**2)
            write(unit=funit,fmt="(A20,20F20.10)") "RD VX=",by%rd%vx

            write(unit=funit,fmt="(A20,20F20.10)") "MS X=",by%ms%x
            write(unit=funit,fmt="(A20,20F20.10)") "MS VX=",by%ms%vx

            write(unit=funit,fmt="(A20,20F20.10)") "MM X=",by%mm%x
            write(unit=funit,fmt="(A20,20F20.10)") "MM VX=",by%mm%vx
        else
            write(*,fmt="(2A20)") "NAME=", trim(adjustl(by%bname))
            write(*,fmt="(A20,I4)") "ANMODE=",by%an_in_mode
            write(*,fmt="(A20,1P20E20.8)") "MS,MM=",by%ms%m,by%mm%m
            call get_star_type(by%ms%obtype, str_ms_obtype)
            call get_star_type(by%mm%obtype, str_mm_obtype)
            write(*,fmt="(3A20)") "TYPE_MS,TYPE_MM=",str_ms_obtype,str_mm_obtype
            
            write(*,fmt="(A20,1P20E20.8)") "Sma,Ecc, Inc=",by%a_bin,by%e_bin,by%Inc
            write(*,fmt="(A20,1P20E20.8)") "E,L=",by%e,by%l
            if(by%a_bin<0d0)then
                write(*,fmt="(A20,1P20E20.8)") "Vinf(kms)=",sqrt(-(by%ms%m+by%mm%m)/by%a_bin)*29.79
            end if
            write(*,fmt="(A20,20F20.10)") "Ome,Pe,Me=",by%Om,by%pe,by%me
            write(*,fmt="(A20,20F20.10)") "f=",by%f0

            write(*,fmt="(A20,20F20.10)") "RD X=",by%rd%x
            write(*,fmt="(A20,20F20.10)") "RD MAG(X)=",sqrt(by%rd%x(1)**2+by%rd%x(2)**2+by%rd%x(3)**2)
            write(*,fmt="(A20,20F20.10)") "RD VX=",by%rd%vx

            write(*,fmt="(A20,20F20.10)") "MS X=",by%ms%x
            write(*,fmt="(A20,20F20.10)") "MS VX=",by%ms%vx

            write(*,fmt="(A20,20F20.10)") "MM X=",by%mm%x
            write(*,fmt="(A20,20F20.10)") "MM VX=",by%mm%vx
        end if
	end subroutine
end module

!!---------------------------system definition
!!--------------------------->>>>>>>>>>>>>>>
!module mydef
!	use binary_def
!	type system
!		type(binary)::by
!		type(binary)::mbhb
!		type(binary)::out_by
!		real(8) E
!	end type
!contains
!	subroutine init_sys(sys)
!		type(system) sys
!		sys%cm=out_by%Mm
!		call get_system_energy(sys)
!	end subroutine
!
!	subroutine get_system_energy(sys)
!		real(8) R12,R23,R13, R14, R24, R34
!		type(system)::sys
!		type(particle)::Ms,Mm,mbhbm, mbhbs
!		Ms=sys%by%Ms;Mm=sys%by%Mm;mbhbs=sys%mbhb%ms; mbhbm=sys%mbhb%mm
!		R12=get_distance(Ms,Mm)
!		R13=get_distance(Ms,mbhbs)
!		R14=get_distance(Ms,mbhbm)
!		R23=get_distance(Mm,mbhbs)
!		R24=get_distance(Mm,mbhbm)
!		R34=get_distance(mbhbs,mbhbm)
!		
!		sys%e=Ms%m*get_Velmag(Ms)**2/2+Mm%m*get_Velmag(Mm)**2/2 &
!				+mbhbs%m*get_Velmag(mbhbs)**2/2+mbhbm%m*get_Velmag(mbhbm)**2/2 &
!		sys%e=sys%e-Ms%m*Mm%m/R12-Ms%m*mbhbs%m/R13-Ms%m*mbhbm%m/R14 &
!			-Mm%m*mbhbs%m/R23-Mm%m*mbhbm%m/R24 -mbhbs%m*mbhbm%m/R34
!	end subroutine
!	subroutine com_sys(sys)
!		IMPLICIT NONE
!		type(system) sys
!		type(particle) cn_p
!
!		sys%by%ms=sys%by%ms+out_by%ms
!		sys%by%mm=sys%by%mm+out_by%ms
!		sys%mbhb%ms=sys%mbhb%ms+out_by%mm
!		sys%mbhb%mm=sys%mbhb%mm+out_by%mm
!
!	end subroutine
!
!	subroutine reflsh_sys(sys,in_by,cm)
!		IMPLICIT NONE
!		type(system) sys
!		type(binary),optional:: in_by,out_by
!		type(particle),optional:: cm
!		type(particle) cn_p
!
!		if(present(in_by).and.present(cm))then
!			sys%out_by%ms=get_center_particle(in_by)
!			sys%out_by%mm=cm
!			call reflsh_by(sys%out_by)
!			call get_system_energy(sys)
!		end if
!
!	end subroutine
!end module
!
!

module md_triple
	use constant
	use md_particle
	use md_binary
contains
	subroutine split_by3p(byin,byot, p1,p2,pc)
		implicit none
		type(particle)::p1,p2,pc
		type(binary)::byin,byot
		p1%x=byin%ms%x+byot%ms%x
		p2%x=byin%mm%x+byot%ms%x
		pc%x=byot%mm%x

		p1%vx=byin%ms%vx+byot%ms%vx
		p2%vx=byin%mm%vx+byot%ms%vx
		pc%vx=byot%mm%vx
	
		p1%m=byin%ms%m
		p2%m=byin%mm%m
		pc%m=byot%mm%m
	end subroutine
	real(8) function Get_triple_energy(p1,p2,p3)
		implicit none
		type(particle)::p1,p2,p3
		real(8) r12,r23,r31, v1,v2,v3
		r12=get_distance(p1,p2);
		r23=get_distance(p2,p3);
		r31=get_distance(p3,p1);
		v1=get_velmag(p1); v2=get_velmag(p2);v3=get_velmag(p3)
		get_triple_energy=0.5d0*p1%m*v1**2+0.5d0*p2%m*v2**2+0.5d0*p3%m*v3**2 &
						-p1%m*p2%m/r12-p1%m*p3%m/r31-p2%m*p3%m/r23
	end function

end module
