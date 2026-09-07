module md_coeff
	!use my_intgl
	!use model_basic
	use com_sts_type
	implicit none
	integer,parameter::Invns=101
	real(8),parameter::inf=-50.1d0, tiny_=1d-9
    real(8):: emin_factor, emax_factor
    real(8):: log10emin_factor, log10emax_factor
	type coeff_type
		real(8) e,j,ee,jj,ej, e_110,e_0,j_111,j_rest, m_avg
	end type
	!type(coeff_type)::coeff_DEJ

    integer,parameter::Jbin_type_lin=1, Jbin_type_log=2, Jbin_type_sqr=3

	type diffuse_coeffient_type
		type(s2d_type)::s2_de, s2_dee,s2_dj, s2_djj,s2_dej
		type(s2d_type)::s2_dRRJJ, s2_dRRJ
		type(s2d_type)::s2_de_110, s2_de_0, s2_dj_111, s2_dj_rest
        type(s1d_type)::s1_de, s1_dee
        real(8) emin, emax, jmin, jmax
		integer nbin!, jbin_type		
		contains
		procedure::init=>init_diffuse_coeffient_grid
		procedure::write_grid=>write_diffuse_coeffient_grid
		procedure::read_grid=>read_diffuse_coeffient_grid
	end type
	type(diffuse_coeffient_type),allocatable::df(:)
	type(diffuse_coeffient_type)::df_tot
	type Inv_type
		integer n
		real(8),allocatable:: xprime(:), I_nvs(:), y2(:)
	end type
	type(Inv_type)::Invs2, Invs3,Invs5, Invs6, Invs7, Invs8, Invs9, Invs10, Invs11
	type(Inv_type)::GBW76, eps1, eps2
	integer coeff_chattery
	private::init_diffuse_coeffient_grid
contains
	subroutine init_invs(inv, n)
		implicit none
		integer n
		type(Inv_type) inv
		if(allocated(inv%xprime))then
			deallocate(inv%xprime, inv%I_nvs, inv%y2)
		end if
		allocate(inv%xprime(n), inv%I_nvs(n), inv%y2(n))
		inv%n=n
	end subroutine
	subroutine conv_j_space(jin,jout,jbin_type)
		implicit none
		real(8) jin,jout
		integer jbin_type
		select case(jbin_type)
		case(jbin_type_log)
			jout=log10(jin)
		case(jbin_type_lin)
			jout=jin
		case(jbin_type_sqr)
			jout=jin**2
		case default
			print*, "define jbin type", jbin_type
			stop
		end select
	end subroutine
	subroutine init_diffuse_coeffient_grid(this, nbin, emin, emax, jmin,jmax, sts_type_dc)
		implicit none
		class(diffuse_coeffient_type)::this
		integer::sts_type_dc
		integer ::nbin
		real(8) emin, emax
		real(8) jmin, jmax
        !real(8) tmin, tmax
        integer jbin_type

		this%nbin=nbin
		this%emin=emin; this%emax=emax; this%jmin=jmin; this%jmax=jmax
        !this%jbin_type=jbin_type
		!call conv_j_space(jmin,tmin,jbin_type)
		!call conv_j_space(jmax,tmax,jbin_type)

		call this%s2_de_110%init( nbin, nbin, emin, emax, jmin,jmax, sts_type_dc)
		call this%s2_de_110%set_range()

		call this%s2_de_0%init( nbin, nbin, emin, emax, jmin,jmax, sts_type_dc)
		call this%s2_de_0%set_range()
        !print*, tmin,tmax
        !print*, this%s2_de_0%ycenter
        !stop
		call this%s2_dee%init(nbin, nbin, emin, emax, jmin,jmax, sts_type_dc)
		call this%s2_dee%set_range()

		call this%s2_dj_111%init( nbin, nbin, emin, emax, jmin,jmax, sts_type_dc)
		call this%s2_dj_111%set_range()

		call this%s2_dj_rest%init(nbin, nbin, emin, emax, jmin,jmax, sts_type_dc)
		call this%s2_dj_rest%set_range()

		call this%s2_djj%init(nbin, nbin, emin, emax, jmin,jmax, sts_type_dc)
		call this%s2_djj%set_range()

		call this%s2_dej%init(nbin, nbin, emin, emax, jmin,jmax, sts_type_dc)
		call this%s2_dej%set_range()

		call this%s2_drrjj%init(nbin, nbin, emin,emax, jmin,jmax, sts_type_dc)
		call this%s2_drrjj%set_range()

		call this%s2_drrj%init(nbin, nbin, emin,emax, jmin,jmax, sts_type_dc)
		call this%s2_drrj%set_range()
		
        call this%s1_de%init(emin,emax, nbin, sts_type_dc)
        this%s1_de%xb=this%s2_de_0%xcenter
!
        call this%s1_dee%init(emin,emax, nbin, sts_type_dc)
        this%s1_dee%xb=this%s2_dee%xcenter
	end subroutine
	subroutine write_diffuse_coeffient_grid(this, file_unit)
		implicit none
		class(diffuse_coeffient_type)::this
		integer file_unit
		write(unit=file_unit) this%nbin
		write(unit=file_unit) this%s2_de_110,this%s2_de_0,this%s2_dee
		write(unit=file_unit) this%s2_dj_111, this%s2_dj_rest, this%s2_djj
		write(unit=file_unit) this%s2_dej, this%s2_drrjj, this%s2_drrj
	end subroutine
	subroutine read_diffuse_coeffient_grid(this, file_unit)
		implicit none
		class(diffuse_coeffient_type)::this
		integer file_unit
		read(unit=file_unit) this%nbin
		read(unit=file_unit) this%s2_de_110,this%s2_de_0,this%s2_dee
		read(unit=file_unit) this%s2_dj_111, this%s2_dj_rest, this%s2_djj
		read(unit=file_unit) this%s2_dej, this%s2_drrjj, this%s2_drrj
	end subroutine
	
end module
