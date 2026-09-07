subroutine gethering_samples_single(fdir, nsnap, bks, bksar, ex)
	use com_main_gw
	implicit none
	integer nsnap, nvalid,i,n
	type(particle_samples_arr_type)::bksar
	type(chain_type)::bks
	type(chain_type)::bys
	character*(*) fdir
	type(particle_samples_arr_type),allocatable::smsa(:)
	type(chain_type),allocatable::sms(:)	
	character*(4) tmprid, tmpspid
	logical ex
	character*(9) str_
	character*(200) fl

	allocate(smsa(ctl%ntask_total))
	allocate(sms(ctl%ntask_total))
	
	write(unit=tmpspid,fmt="(I4)") nsnap

	nvalid=0
	do i=1, ctl%ntask_total
		write(unit=tmprid,fmt="(I4)") i
		str_=trim(adjustl(tmprid))//"_"//trim(adjustl(tmpspid))
		fl=trim(adjustl(fdir))//"/bin/single/samchn"//trim(adjustl(str_))
		inquire(file=trim(adjustl(fl))//".bin",exist=ex)
		if(ex)then
			nvalid=nvalid+1
			!call smsa(nvalid)%input_bin(trim(adjustl(fl)))
			!print*, "arr nvalid=",nvalid
			call sms(nvalid)%input_bin(trim(adjustl(fl)))	
            call all_chain_to_arr_single(sms(nvalid),smsa(nvalid))
		!	print*, "ch nvalid=",nvalid		
		else
			print*, trim(adjustl(fl))//" does not exist"
            return
		endif
	end do
	print*, "sms chain gathering finished"

	call smmerge_arr_single(smsa,nvalid,bksar)
	!print*, "gathering smsa"
	!print*, smsa(1)%sp(1:10)%weight_real
	!stop
	call smmerge(sms, nvalid,bks)
	!call dms%input_bin(trim(adjustl(fdir))//"/dms/dms_"//trim(adjustl(tmpspid)))
	!==================================
	!print*,"by reading skiped"
	if(nvalid.ge.1)then
		call bks%destory()
	end if
	!return
	nvalid=0
    
end subroutine
