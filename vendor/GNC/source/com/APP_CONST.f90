module Astron_constant
	real(8),parameter::pc=206264.98d0      !AU
	real(8),parameter::rd_sun=4.65d-3      !AU
	real(8),parameter::AU_GS=1.4959787d13  !cm
	real(8),parameter::AU_SI=1.4959787d11  !m  
	real(8),parameter::m_sun_GS=1.98855d33    !g
	real(8),parameter::m_sun_SI=1.98855d30    !kg
	real(8),parameter::one_year=3.1556d7    !s
	real(8),parameter::pi4degree2=41252.96125 ! deg**2
end module
module unitless_value
	real(8),parameter::PI=3.141592653589793d0   
	real(8),parameter::TWO_PI=3.141592653589793d0*2d0   
	real(8),parameter::e_nature=2.718281828459045d0
!	real(8),parameter::moer_const=6.02214129e23 !avogadro constant
!	real(8),parameter::PI=acos(-1d0)
!	real(8),parameter::alpha_FC=7.297352570d-3     !Fine-structure constant
end module

module my_unit
	use unitless_value
	!G=1,M=1msun, Length=1AU
	!then the velocity_unit=29.79kms-1
	real(8),parameter::my_unit_vel_c=3d5/29.784d0
!    real(8),parameter::my_unit_of_time=58.123 !days
!	real(8),parameter::my_unit_of_energy=6.138d5 !solar_lumi also =3.53d32 W
contains
	real(8) function myu_conv_t2second(t)
		implicit none
		real(8) t
		myu_conv_t2second=t*365.2425d0/2/pi*86400      !(second)
	end function 
	 
end module
module constant 
use unitless_value     !Unitless values
use Astron_constant        
use my_unit

end module

