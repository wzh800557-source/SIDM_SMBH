module md_chain
    use md_chain_pointer
    type chain_type
		integer n 
        type(chain_pointer_type),pointer::head => null()
        contains
        procedure::init=>init_chain
        procedure::destory=>destory_chain
        procedure::get_length=>chains_get_length
        procedure::copy=>copy_a_chain
        procedure::output_bin=>output_chains_bin
        procedure::output_txt=>output_chains_txt
        procedure::input_bin=>input_chains_bin
        procedure::insert_after_chain
        procedure::output_screen=>chain_output_list_chain_type
	end type

    private::chains_get_length,copy_a_chain,chain_output_list_chain_type
    private::output_chains_bin, input_chains_bin,destory_chain,output_chains_txt
contains

subroutine init_chain(sps,n)
    implicit none
    class(chain_type)::sps
    integer n,i
    type(chain_pointer_type),pointer::p, pc
    if(associated(sps%head)) then
        call sps%destory()
    end if
    sps%n=n
    
    allocate(sps%head)
    call sps%head%init_head()
    
    p=>sps%head
    pc=>sps%head
    do i=1, n-1
        allocate(pc%next)
        p=>pc%next
        p%prev=>pc
        p%bg=>pc%bg
        p%idx=pc%idx+1
        !print*, "i=",i, p%idx
        pc=>pc%next
    end do
    call p%set_end()
end subroutine


subroutine attach_chain_end(this, c)
    implicit none
    type(chain_type),target:: this, c
    type(chain_pointer_type),pointer::phead, pend,p, chead, cend
    phead=>c%head
    chead=>this%head
    pend=>this%head%ed
    cend=>c%head%ed

    pend%next=>phead
    phead%prev=>pend
    p=>phead
    do while(associated(p))
        p%idx=p%prev%idx+1
        p=>p%next
    end do
    call set_item_chain_head_chain_type(chead)
    call cend%set_end()
end subroutine
subroutine attach_two_chains(this, c2)
    ! after this operation, c2 is messed up 
    ! and can not be used anymore
    implicit none
    type(chain_type) this, c2
    type(chain_pointer_type),pointer::p
    !logical::condition

    call attach_chain_end(this,c2)  
    this%n=this%n+c2%n
    
end subroutine

subroutine insert_after_chain(chain,length) 
    class(chain_type),target::chain
    type(chain_pointer_type),pointer::p, pc
    integer length,i

    pc=>chain%head%ed
    p=>chain%head%ed
    !print*, "create_chain_single:", chain%idx, chain%ob%id
    if(.not.associated(chain%head))then
        print*, "pc not associated"
        stop
    end if	        
    do i=1, length
        allocate(pc%next)
        p=>pc%next
        p%prev=>pc
        p%bg=>pc%bg
        p%idx=pc%idx+1
        pc=>pc%next
    end do
    call chain%head%ed%next%set_end()
    !print*, "p%idx=",p%idx, chain%ed%idx
    !call chain%output()
    !call chain%next%next%output()
    !read(*,*)
end subroutine

subroutine destory_chain(chain)
    implicit none
    class(chain_type),target::chain
    type(chain_pointer_type),pointer::p, pc
    if(associated(chain%head))then
        p=>chain%head%ed
        pc=>p
        do while(associated(pc%prev))
            p=>p%prev
            call destroy_attach_pointer_chain_type(pc)			
            pc=>p
        end do
        call destroy_attach_pointer_chain_type(chain%head)
        !print*, "finished destroy"
    end if
end subroutine

subroutine smmerge(sma,n,smam)
    !after this subroutine, sma(i) is modified to be connected with each other
    implicit none
    integer n,i,j, nsam
    type(chain_type)::sma(n), smam
    type(chain_pointer_type),pointer::pt
    allocate(smam%head)

       smam%head=>sma(1)%head
    do i=2, n
        call attach_two_chains(smam,sma(i))
    end do
    
end subroutine

subroutine chains_get_length(chain, n,type)
    implicit none
    class(chain_type),target::chain
    type(chain_pointer_type),pointer::p
    integer,optional::type
    integer n, typeI
    p=>chain%head
    n=0
    typeI=0
    if(present(type))then
        typeI=type
    end if
    !print*, "typeI=", typeI
    !do while(associated(p).and.allocated(p%ob))
    do while(associated(p))
        select case(typeI)
        case(0)
            n=n+1
        case(1)
            select type(ca=>p%ob)
            type is(particle_sample_type)
                n=n+1
            end select
        !case(2)
        !    select type(ca=>p%ob)
        !    type is(sample_type)
        !        n=n+1
        !    end select
        !    !print*, "n=", n
        end select
        p=>p%next
    end do
end subroutine

subroutine copy_a_chain(ch,ch_copy)
    implicit none
    class(chain_type),target::ch
    type(chain_type),target::ch_copy
    type(chain_pointer_type),pointer::p, p_copy
    integer length !, i

    call ch%get_length(length)
    call ch_copy%init(length)
    p=>ch%head
    p_copy=>ch_copy%head
    do while(associated(p))
        if(.not.allocated(p_copy%ob)) then
            allocate(p_copy%ob,source=p%ob)
        end if
        call p%copy(p_copy)
        p=>p%next
        p_copy=>p_copy%next
    end do
end subroutine

subroutine split_two_chains_by_types(sma,smbk, smby)
    implicit none
    integer ntot, nbk, nby,i,j, nsam
    type(chain_type)::sma,smbk, smby
    type(chain_pointer_type),pointer::pt
    integer type_flag, n
    call sma%get_length(ntot)
    call sma%get_length(nbk,1)
    call sma%get_length(nby,2)
    print*, "ntot, nbk, nby=", ntot, nbk, nby
    type_flag=1
    call chain_select_by_condition(sma,smbk, selection)
    print*, "bk length=", n
    type_flag=2
    call chain_select_by_condition(sma,smby, selection)
    print*, "by length=", n
contains
    logical function selection(pt)
        implicit none
        type(chain_pointer_type)::pt
        logical cond
        cond=.false.
        select case(type_flag)
        case(1)
            select type(ca=>pt%ob)
            type is(particle_sample_type)
                cond=.true.
            end select
        !case(2)
        !    select type(ca=>pt%ob)
        !    type is(sample_type)
        !        cond=.true.
        !    end select
        end select
        selection=cond
    end function

end subroutine

subroutine chain_select_by_condition(ch, ch_out, selection)
    implicit none	
    !type(particle_sample_type)::sps, sps_out
    type(chain_type),target::ch, ch_out
    type(chain_pointer_type),pointer::ps,psout
    logical,external::selection
    integer nsel, n

    ps=>ch%head
    nsel=0
    do while (associated(ps).and.allocated(ps%ob))
        if(selection(ps))then
            nsel=nsel+1
        end if
        ps=>ps%next
    end do
    call ch_out%init(nsel)
    psout=>ch_out%head
    ps=>ch%head
    do while(associated(ps).and.allocated(ps%ob))
        if(selection(ps))then
            !print*, "associated(psout)=", associated(psout)
            if(associated(psout))then
                if(.not.allocated(psout%ob)) allocate(psout%ob,source=ps%ob)
                call ps%copy(psout)
                !print*, "copied"
                psout=>psout%next
            end if
        end if
        ps=>ps%next
    end do
    !call psout%get_length(n)
   ! print*, "chain created, length=",n
end subroutine

subroutine output_chains_bin(sps, fl)
	!use com_main_gw
	implicit none
	character*(*) fl
	integer i,n
	class(chain_type)::sps
	type(chain_pointer_type),pointer::pt
    integer,parameter::flag_sg=1,flag_by=2

	open(unit=999,file=trim(adjustl(fl))//".bin", form='unformatted',access='stream')
	
    call sps%get_length(n)
    print*, "chain length=",n
    write(999) n
    !print*, "n=",n
    pt=>sps%head
    do while(associated(pt))
        select type (ca=>pt%ob)
        type is(particle_sample_type)
            write(999) flag_sg
            call ca%write_info(999)
        !type is(sample_type)
        !    write(999) flag_by
        !    call ca%write_info(999)
        end select        
        pt=>pt%next
    end do
	close(unit=999)
end subroutine

subroutine input_chains_bin(sps, fl)
	!use com_main_gw
	implicit none
	character*(*) fl
	integer i,n
	class(chain_type)::sps
	type(chain_pointer_type),pointer::pt
    integer flag
    integer,parameter::flag_sg=1,flag_by=2

	open(unit=999,file=trim(adjustl(fl))//".bin", form='unformatted',access='stream',status='old')
	read(unit=999) n
    print*, "chain length=", n
    call sps%init(n)
    pt=>sps%head
    do while(associated(pt))
        read(999) flag
        select case(flag)
        case(flag_sg)
            allocate(particle_sample_type::pt%ob)
            call pt%ob%read_info(999)
        !case(flag_by)
        !    allocate(sample_type::pt%ob)
        !    call pt%ob%read_info(999)
        end select        
        pt=>pt%next
    end do
	close(unit=999)
end subroutine

subroutine output_chains_txt(ch, fl)
    implicit none
    class(chain_type)::ch
    character*(*) fl
    character*(6) str_type
    integer i
    type(chain_pointer_type),pointer::pt
    open(unit=999,file=trim(adjustl(fl))//".txt")

    pt=>ch%head
    i=0
    write(unit=999,fmt="(2A6,33A15)") "i", "type", "idx", "id",&
         "m", "w_N", "ao","eo", "in_ms", "in_mm", "ain","ein"
    do while(associated(pt))
        i=i+1
        select type (ca=>pt%ob)
        !type is(sample_type)
        !    write(unit=999,fmt="(I6, A6,3I15,30E15.5)") i, trim(adjustl(binary_types(ca%obtype))), &
        !        ca%idx, ca%id,ca%N_gene, ca%m, ca%weight_N, &
        !        ca%byot%a_bin,ca%byot%e_bin, ca%byin%ms%m, ca%byin%mm%m, ca%byin%a_bin,ca%byin%e_bin
        type is(particle_sample_type)
            write(unit=999,fmt="(I6, A6,3I15,30E15.5)") i, trim(adjustl(star_type(ca%obtype))),&
                ca%idx, ca%id, ca%m, ca%weight_N, ca%byot%a_bin,ca%byot%e_bin
        end select
        pt=>pt%next
    end do
    close(999)
end subroutine

subroutine output_chains_txt_bdsample(ch,ebd, fl)
    implicit none
    class(chain_type)::ch
    character*(*) fl
    character*(6) str_type
    integer i
	real(8) ebd
    type(chain_pointer_type),pointer::pt
    open(unit=999,file=trim(adjustl(fl))//".txt")

    pt=>ch%head
    i=0
    write(unit=999,fmt="(2A6,33A15)") "i", "en", "jm", "en0", "jm0"
    do while(associated(pt))
        i=i+1
        select type (ca=>pt%ob)
        !type is(sample_type)
            
        type is(particle_sample_type)
			if(ca%en>ebd)then
            	write(unit=999,fmt="(I6, 30E15.5)") i, ca%en, ca%jm, ca%en0, ca%jm0
			end if
        end select
        pt=>pt%next
    end do
    close(999)
end subroutine

subroutine output_chains_txt_if(ch, fl)
    implicit none
    character*(*) fl
    integer i,n
    class(chain_type)::ch
    type(chain_pointer_type),pointer::pt

    open(unit=999,file=trim(adjustl(fl))//".txt")
    write(unit=999,fmt="(13A20, 3A10)") "t0", "e0", "tf", "ef"
    pt=>ch%head
    do while(associated(pt))
        write(unit=999,fmt="(13E20.10, 3I10)") pt%ob%create_time, pt%ob%en0, pt%ob%exit_time, pt%ob%en
        pt=>pt%next
    end do
    close(unit=999)
end subroutine

subroutine	chain_output_list_chain_type(list,flag_output_in, maxnum_in)
	IMPLICIT NONE
		class(chain_type),target:: list
		type(chain_pointer_type),pointer::p
		character*4 itmp
		character*20 fmts,fmts2
        integer flag, flag_output, i
        integer,optional::flag_output_in
        integer,optional::maxnum_in
        flag_output=0
        if(present(flag_output_in))flag_output=flag_output_in
		p=>list%head
        i=0
		do while(associated(p))
             
			if(associated(p%prev))then
				write(*,fmt="(A2)", advance='no') "<-"
			else
				write(*,fmt="(A6)", advance='no') "Null<-"
			end if
			if(p%idx.eq.0)then
				itmp="1"
			else
				write(unit=itmp, fmt="(I4)") int(log10(real(p%idx)))+1
			end if
			fmts="(I"//trim(adjustl(itmp))//",L1,I1,F10.1)"
            fmts2="(I"//trim(adjustl(itmp))//",L1,2I2,I10)"
            flag=0
            select type(ca=>p%ob)
            type is(particle_sample_type)
                flag=1
            !type is(sample_type)
            !    flag=2
            end select
            select case(flag_output)
            case(0)
                write(*,fmt="(L1)", advance='no') allocated(p%ob)
            case(1)
			    write(*,fmt=fmts, advance='no') P%idx, allocated(p%ob), flag, p%ob%en
            case(2)
			    write(*,fmt=fmts2, advance='no') P%idx, allocated(p%ob), flag,p%ob%exit_flag,&
                     p%ob%id
            end select
			if(associated(p%next))then
				write(*,fmt="(A2)", advance='no') "->"
			else
				write(*,fmt="(A6)", advance='no') "->Null"
			end if
			p=>p%next

            i=i+1
            if(present(maxnum_in))then
                if(i.ge.maxnum_in) exit 
            end if
		end do
		write(*,*) 	
		write(*,*) "--------------------"
	end subroutine

    subroutine chain_select(sps, sps_out, exitflag, obj_type)
        implicit none	
        type(chain_type)::sps, sps_out
        type(chain_pointer_type),pointer::ps,psout
        integer nsel,i, exitflag
        integer,optional:: obj_type
        integer objtype
        objtype=0
        if(present(obj_type)) objtype=obj_type
    
        call chain_select_by_condition(sps,sps_out, selection)
    
    contains
        logical function selection(pt)
            implicit none
            type(chain_pointer_type)::pt
            logical cond
            cond=.false.
            select case(objtype)
            case(0)
                cond=.true.
            case(1)
                select type(ca=>pt%ob)
                type is(particle_sample_type)
                   cond=.true.
                end select
            end select
            if(((pt%ob%exit_flag.eq.exitflag).or.exitflag.eq.-1)&
                .and.cond)then
                selection=.True.
            else
                selection=.False.
            end if
        end function
    end subroutine

end module