module MPI_comu
	use model_basic
	!use md_sts
	use mpi
	implicit none
	!include 'mpif.h'
	integer rid, proc_id
	integer:: mpi_master_id=0
    integer,parameter::MPI_COMMAND_END=-1
	integer,parameter::MPI_COMMAND_comu_sent=4
	integer,parameter::MPI_COMMAND_comu_recv=5
	integer,parameter::mpi_command_own=6
contains
	subroutine init_mpi()
		implicit none
		integer ierr, getpid
		!print*, "start"
		call mpi_init(ierr)
		!print*, "init"
		call mpi_comm_rank(mpi_comm_world, rid,ierr)
		!print*, "rid=", rid
		!print*, "size",ctl%ntasks
		call mpi_comm_size(mpi_comm_world, ctl%ntasks, ierr)
		proc_id=getpid()
	! print*, "size", rid
		if(ctl%chattery>4)then 
			print*, "mpi initializtion finished"
		end if
        !print*, "declare particle"
		!call declare_particle_mpi_type()
        !print*, "declare binary"
		!call declare_binary_mpi_type()
		
	end subroutine
	subroutine stop_mpi()
		implicit none
		integer ierr
		!print*, "start stop",rid
		!call MPI_Type_free(PARTICLE_MPI_TYPE_, ierr)
		!print*, "particle free",rid
		!call MPI_Type_free(BINARY_MPI_TYPE_, ierr)
		!print*, "binary free",rid
		call mpi_finalize(ierr)
		print*, "mpi finalized",rid
	end subroutine
	subroutine collect_data_mpi(fxy, nbin, nbg, ned, nblocks,ntasks)
		!use com_main_gw
		implicit none
		!type(diffuse_coeffient_type)::dc
		integer nbin
		real(8) fxy(nbin, nbin)
		real(8),allocatable::reivbuffer(:)
		real(8),allocatable::sentbuffer(:)
		integer nd_tot,i, j, idbg, ided, sent_count, ierr
		integer nbg, ned,nblocks, ntasks

		nd_tot=nbin*nbin
		allocate(reivbuffer(nd_tot))
		allocate(sentbuffer(nd_tot))
		reivbuffer=0
		if(mod(nbin, ntasks).ne.0)then
			print*, "error in collect data mpi!: ntasks is not times of nbin"
			stop
		end if
		sent_count=nbin*nblocks
		!print*, "rid, sent_count, nd_tot=", rid, sent_count, nd_tot
		do i=1, ntasks
			do j=nbg, ned
				idbg=(i-1)*sent_count+(j-nbg)*nbin+1
				ided=idbg+nbin-1
				!print*, "idbg, ided=", idbg, ided
				sentbuffer(idbg:ided)=fxy(j,1:nbin)
			end do
		end do
		!print*, "rid,sentbuffer=", rid,  sentbuffer
		!print*
		call mpi_alltoall(sentbuffer,sent_count, MPI_DOUBLE, reivbuffer, sent_count, &
			MPI_DOUBLE, mpi_comm_world, ierr)
		!print*, "alltoall finished"
		!print*, "rid,reibuffer=", rid,  reivbuffer
		!print*
		
		do i=1, nbin
			!do j=1, dc%nbin
			idbg=(i-1)*nbin+1
			ided=i*nbin
			!print*, "rid,idbg, ided=", rid, idbg, ided
			fxy(i, 1:nbin)=reivbuffer(idbg:ided)
			!end do
		end do
		!call mpi_BARRIER(mpi_comm_world, ierr)
		!print*, "rid=", rid
		!!call dc%s2_de_110%print()

	end subroutine
	subroutine collect_to_root_sps_single(sps_send,sps, n)
		!use com_main_gw
		implicit none
		integer i,n
		type(particle_samples_arr_type)::sps(n)
		type(particle_samples_arr_type) sps_send
		!allocate(sps(ctl%ntasks))
		do i=1, ctl%ntasks
			if(i.ne.mpi_master_id+1)then
				!print*, "rid,i=",rid,i
				call send_particle_sample_arr_mpi(sps_send,sps(i), i-1, mpi_master_id)
				!print*, "collect_to_root_sps_single:i=",i
			!end if
			else
				sps(mpi_master_id+1)=sps_send
				!print*, "direct :i=",i
			end if
		end do
	end subroutine

	subroutine bcast_dms_barge(dm)
		!use com_main_gw
		implicit none
		type(diffuse_mspec)::dm
		integer i,ierr!, n, flag_log

		do i=1, dm%n
			call bcast_s1d_mpi(dm%mb(i)%all%barge)
            call bcast_s1d_mpi(dm%mb(i)%sbh%barge)
            call bcast_s1d_mpi(dm%mb(i)%star%barge)
            call bcast_s1d_mpi(dm%mb(i)%bbh%barge)       
		end do

		call bcast_s1d_mpi(dm%all%all%barge)
	end subroutine
	subroutine bcast_dms_fden(dm)
		!use com_main_gw
		implicit none
		type(diffuse_mspec)::dm
		integer i,ierr!, n, flag_log

		do i=1, dm%n
			call bcast_s1d_mpi(dm%mb(i)%all%fden)
            call bcast_s1d_mpi(dm%mb(i)%all%fden_simu)
		end do
		call bcast_s1d_mpi(dm%all%all%fden)
        call bcast_s1d_mpi(dm%all%all%fden_simu)

	end subroutine

	subroutine send_particle_sample_arr_mpi(sps_send, sps_recv, proc_id_source,proc_id_dest)
		!use com_main_gw
		implicit none
		type(particle_samples_arr_type),intent(in)::sps_send
		type(particle_samples_arr_type),intent(out)::sps_recv
		integer ierr, proc_id_source,proc_id_dest, i, nintarr_sp, nrealarr_sp, nintarr_by,nrealarr_by, nstrarr_by
		integer status(MPI_Status_size)
		integer,allocatable::intarr_sp(:,:), intarr_by(:,:)
		real(8),allocatable::realarr_sp(:,:), realarr_by(:,:)
		character*(nstr_by),allocatable:: strarr_by(:)
		!type(binary),allocatable::by(:), byini(:),bybf(:)
		

		if(rid.eq.proc_id_source)then
			!print*, "source:rid=",rid,proc_id_dest
			allocate(intarr_sp(nint_particle,sps_send%n))
			allocate(realarr_sp(nreal_particle,sps_send%n))
			!allocate(by(sps_send%n),byini(sps_send%n),bybf(sps_send%n))
			allocate (intarr_by(nint_by,sps_send%n))
			allocate(realarr_by(nreal_by,sps_send%n))
			allocate(strarr_by(sps_send%n))

			do i=1, sps_send%n
				call conv_sp_int_real_arrays(sps_send%sp(i), intarr_sp(1:nint_particle,i), realarr_sp(1:nreal_particle,i))
				!by(i)=sps_send%sp(i)%byot
				!byini(i)=sps_send%sp(i)%byot_ini
				!bybf(i)=sps_send%sp(i)%byot_bf
			end do
			call mpi_send(sps_send%n,1, MPI_INTEGER,proc_id_dest,0,MPI_COMM_WORLD,ierr)
			!print*, "sps_send:n=", sps_send%n
			call mpi_send(intarr_sp, nint_particle*sps_send%n, MPI_INTEGER, proc_id_dest, 0, MPI_COMM_WORLD, ierr)
			!!print*, "intarr sent:rid=",rid
			call mpi_send(realarr_sp, nreal_particle*sps_send%n, MPI_DOUBLE_PRECISION, proc_id_dest, 0, MPI_COMM_WORLD, ierr)
			!!print*, "realarr sent:rid=",rid,realarr(1:2,1)
            !!print*, "by sent:rid=",rid,by(1)%a_bin, proc_id_source, proc_id_dest
			do i=1, sps_send%n
				call conv_by_arrays(sps_send%sp(i)%byot,intarr_by(1:nint_by,i),realarr_by(1:nreal_by,i),strarr_by(i))
			end do
			call mpi_send(intarr_by, nint_by*sps_send%n, MPI_INTEGER, proc_id_dest, 0, MPI_COMM_WORLD, ierr)
			call mpi_send(realarr_by, nreal_by*sps_send%n, MPI_DOUBLE_PRECISION, proc_id_dest, 0, MPI_COMM_WORLD, ierr)
			call mpi_send(strarr_by, nstr_by*sps_send%n, MPI_CHARACTER, proc_id_dest, 0, MPI_COMM_WORLD, ierr)

			do i=1, sps_send%n
				call conv_by_arrays(sps_send%sp(i)%byot_ini,intarr_by(1:nint_by,i),realarr_by(1:nreal_by,i),strarr_by(i))
			end do
			call mpi_send(intarr_by, nint_by*sps_send%n, MPI_INTEGER, proc_id_dest, 0, MPI_COMM_WORLD, ierr)
			call mpi_send(realarr_by, nreal_by*sps_send%n, MPI_DOUBLE_PRECISION, proc_id_dest, 0, MPI_COMM_WORLD, ierr)
			call mpi_send(strarr_by, nstr_by*sps_send%n, MPI_CHARACTER, proc_id_dest, 0, MPI_COMM_WORLD, ierr)

			do i=1, sps_send%n
				call conv_by_arrays(sps_send%sp(i)%byot_bf,intarr_by(1:nint_by,i),realarr_by(1:nreal_by,i),strarr_by(i))
			end do
			call mpi_send(intarr_by, nint_by*sps_send%n, MPI_INTEGER, proc_id_dest, 0, MPI_COMM_WORLD, ierr)
			call mpi_send(realarr_by, nreal_by*sps_send%n, MPI_DOUBLE_PRECISION, proc_id_dest, 0, MPI_COMM_WORLD, ierr)
			call mpi_send(strarr_by, nstr_by*sps_send%n, MPI_CHARACTER, proc_id_dest, 0, MPI_COMM_WORLD, ierr)
			!do i=1, 2
			!	call sps_send%sp(i)%print("sps_send")
			!end do

		elseif(rid.eq.proc_id_dest)then
			!print*, "dest:rid=",rid,proc_id_source
			call mpi_recv(sps_recv%n,1, MPI_INTEGER,proc_id_source,0,MPI_COMM_WORLD,status,ierr)
			!print*, "sps_recv:n=", sps_recv%n, rid

			call sps_recv%init(sps_recv%n)
			allocate(intarr_sp(nint_particle,sps_recv%n))
			allocate(realarr_sp(nreal_particle,sps_recv%n))
			!allocate(by(sps_recv%n))
			!allocate(byini(sps_recv%n),bybf(sps_recv%n))
			allocate (intarr_by(nint_by,sps_recv%n))
			allocate(realarr_by(nreal_by,sps_recv%n))
			allocate(strarr_by(sps_recv%n))
			
			nintarr_sp=nint_particle*sps_recv%n; nrealarr_sp=nreal_particle*sps_recv%n
			nintarr_by=nint_by*sps_recv%n;    nrealarr_by=nreal_by*sps_recv%n
			nstrarr_by=nstr_by*sps_recv%n

			call mpi_recv(intarr_sp, nintarr_sp, MPI_INTEGER, proc_id_source, 0, mpi_comm_world,status, ierr)
			!!print*, "intarr recv:rid=",rid
			call mpi_recv(realarr_sp, nrealarr_sp, MPI_DOUBLE_PRECISION, proc_id_source, 0, mpi_comm_world, status,ierr)
		
			do i=1, sps_recv%n
				call conv_int_real_arrays_sp(sps_recv%sp(i), intarr_sp(1:nint_particle, i), realarr_sp(1:nreal_particle,i))
			end do
			call mpi_recv(intarr_by, nintarr_by, MPI_INTEGER, proc_id_source, 0, mpi_comm_world,status, ierr)
			call mpi_recv(realarr_by, nrealarr_by, MPI_DOUBLE_PRECISION, proc_id_source, 0, mpi_comm_world, status,ierr)
			call mpi_recv(strarr_by, nstrarr_by, MPI_CHARACTER, proc_id_source, 0, mpi_comm_world, status,ierr)
			do i=1, sps_recv%n
				call conv_arrays_by(sps_recv%sp(i)%byot, intarr_by(1:nint_by, i), realarr_by(1:nreal_by,i), strarr_by(i))
			end do

			call mpi_recv(intarr_by, nintarr_by, MPI_INTEGER, proc_id_source, 0, mpi_comm_world,status, ierr)
			call mpi_recv(realarr_by, nrealarr_by, MPI_DOUBLE_PRECISION, proc_id_source, 0, mpi_comm_world, status,ierr)
			call mpi_recv(strarr_by, nstrarr_by, MPI_CHARACTER, proc_id_source, 0, mpi_comm_world, status,ierr)
			do i=1, sps_recv%n
				call conv_arrays_by(sps_recv%sp(i)%byot_ini, intarr_by(1:nint_by, i), realarr_by(1:nreal_by,i), strarr_by(i))
			end do
			call mpi_recv(intarr_by, nintarr_by, MPI_INTEGER, proc_id_source, 0, mpi_comm_world,status, ierr)
			call mpi_recv(realarr_by, nrealarr_by, MPI_DOUBLE_PRECISION, proc_id_source, 0, mpi_comm_world, status,ierr)
			call mpi_recv(strarr_by, nstrarr_by, MPI_CHARACTER, proc_id_source, 0, mpi_comm_world, status,ierr)
			do i=1, sps_recv%n
				call conv_arrays_by(sps_recv%sp(i)%byot_bf, intarr_by(1:nint_by, i), realarr_by(1:nreal_by,i), strarr_by(i))
			end do
			!do i=1, 2
			!	call sps_recv%sp(i)%print("sps_recv")
			!end do
			!stop
		end if

	end subroutine
	

	subroutine conv_sp_int_real_arrays(sp, intarr, realarr)
		!use com_main_gw
		implicit none
		class(particle_sample_type)::sp
		integer nint, nreal
		integer intarr(nint_particle)
		real(8) realarr(nreal_particle)
		intarr(1:5)=(/sp%id,sp%obtype,sp%obidx,sp%state_flag_last,sp%exit_flag/)
		intarr(6:10)=(/ sp%within_jt, sp%rid, sp%idx, sp%length,sp%write_down_track/)

		realarr(1:5)=(/sp%create_time,sp%den,sp%djp,sp%djp0,sp%elp/)
		realarr(6:10)=(/sp%En,sp%en0,sp%exit_time,sp%Jm,sp%weight_clone/)
		realarr(11:14)=(/sp%weight_real,sp%m,sp%r_td,sp%jm0/)
		realarr(15:19)=(/sp%pd,sp%rp, sp%tgw, sp%weight_N,sp%simu_bgtime/)	
	end subroutine

	subroutine conv_int_real_arrays_sp(sp, intarr, realarr)
		!use com_main_gw
		implicit none
		class(particle_sample_type)::sp
		integer nint, nreal
		integer intarr(nint_particle)
		real(8) realarr(nreal_particle)
		sp%id=intarr(1); sp%obtype=intarr(2); sp%obidx=intarr(3)
		sp%state_flag_last=intarr(4);sp%exit_flag=intarr(5)
		!intarr(1:5)=(/sp%id,sp%obtype,sp%obidx,sp%state_flag_last,sp%exit_flag/)
		sp%within_jt=intarr(6);sp%rid=intarr(7);
		sp%idx=intarr(8)
		sp%length=intarr(9); sp%write_down_track=intarr(10)
		call track_init(sp, 0)
        
		sp%create_time=realarr(1); sp%den=realarr(2); sp%djp=realarr(3)
		sp%djp0=realarr(4);sp%elp=realarr(5)
		!realarr(1:5)=(/sp%create_time,sp%den,sp%djp,sp%djp0,sp%elp/)
		sp%en=realarr(6);sp%en0=realarr(7);sp%exit_time=realarr(8);
		sp%jm=realarr(9);sp%weight_clone=realarr(10)
		!realarr(6:10)=(/sp%En,sp%en0,sp%exit_time,sp%Jm,sp%weight/)
		!sp%weight_asym=realarr(11);
        sp%weight_real=realarr(11); sp%m=realarr(12)
		sp%r_td=realarr(13);sp%jm0=realarr(14)
		!realarr(11:15)=(/sp%weight0,sp%weight_real,sp%m,sp%r_td,sp%jm0/)
		sp%pd=realarr(15);sp%rp=realarr(16);sp%tgw=realarr(17)
		sp%weight_N=realarr(18)
        sp%simu_bgtime=realarr(19)
		!call sp%print("conv_int_real_arrays_sp")
		!realarr(16:18)=(/sp%pd,sp%rp, sp%tgw/)	
	end subroutine
	subroutine conv_arrays_by(by,intarr,realarr,strarr)
		implicit none
		type(binary)::by
		integer intarr(nint_by)
		real(8) realarr(nreal_by)
		character*(nstr_by) strarr
		by%ms%x=realarr(1:3)
		by%ms%vx=realarr(4:6)
		by%ms%m=realarr(7); by%ms%radius=realarr(8)
		by%ms%id=intarr(1); by%ms%obtype=intarr(2); by%ms%obidx=intarr(3)

		by%mm%x=realarr(9:11)
		by%mm%vx=realarr(12:14)
		by%mm%m=realarr(15); by%mm%radius=realarr(16)
		by%mm%id=intarr(4); by%mm%obtype=intarr(5); by%mm%obidx=intarr(6)
		
		by%rd%x=realarr(17:19)
		by%rd%vx=realarr(20:22)
		by%rd%m=realarr(23); by%rd%radius=realarr(24)
		by%rd%id=intarr(7); by%rd%obtype=intarr(8); by%rd%obidx=intarr(9)

		by%e=realarr(25); by%l=realarr(26); by%k=realarr(27)
		by%miu=realarr(28); by%mtot=realarr(29); by%Jc=realarr(30)
		by%a_bin=realarr(31); by%e_bin=realarr(32); by%lum=realarr(33:35)
		by%f0=realarr(36); by%inc=realarr(37); by%om=realarr(38)
		by%pe=realarr(39); by%t0=realarr(40); by%me=realarr(41)

		by%an_in_mode=intarr(10)
		by%bname=strarr
	end subroutine
	subroutine conv_by_arrays(by,intarr,realarr,strarr)
		implicit none
		type(binary)::by
		integer intarr(nint_by)
		real(8) realarr(nreal_by)
		character*(nstr_by) strarr
		realarr(1:3)=by%ms%x
		realarr(4:6)=by%ms%vx
		realarr(7:8)=(/by%ms%m, by%ms%radius/)
		intarr(1:3)=(/by%ms%id,by%ms%obtype,by%ms%obidx/)

		realarr(9:11)=by%mm%x
		realarr(12:14)=by%mm%vx
		realarr(15:16)=(/by%mm%m, by%mm%radius/)
		intarr(4:6)=(/by%mm%id,by%mm%obtype,by%mm%obidx/)

		realarr(17:19)=by%rd%x
		realarr(20:22)=by%rd%vx
		realarr(23:24)=(/by%rd%m, by%rd%radius/)
		intarr(7:9)=(/by%rd%id,by%rd%obtype,by%rd%obidx/)

		realarr(25:41)=(/by%e, by%l, by%k, by%miu, by%mtot, by%Jc, &
			by%a_bin, by%e_bin, by%lum(1:3), by%f0, by%inc, by%om, by%pe, &
			by%t0, by%me/)

		intarr(10)=by%an_in_mode
		strarr=by%bname
	end subroutine


	subroutine bcast_s1d_mpi(s1d)
		!use com_main_gw
		implicit none
		type(s1d_type)::s1d
		integer::intarr(s1d%type_int_size)
		integer nint, nreal, ierr
		real(8)::realarr(s1d%type_real_size)
		logical::logarr(s1d%type_log_size)

		!nint=1
		!nreal=
		!allocate(intarr(nint), realarr(nreal))
		if(rid.eq.mpi_master_id)then
			call conv_s1d_int_real_arrays(s1d, intarr, realarr,logarr)
		end if
		call mpi_bcast(intarr, s1d%type_int_size, MPI_INTEGER, mpi_master_id, mpi_comm_world,ierr)
		call mpi_bcast(realarr, s1d%type_real_size, MPI_DOUBLE_PRECISION, mpi_master_id, mpi_comm_world,ierr)
		call mpi_bcast(logarr, s1d%type_log_size, MPI_LOGICAL, mpi_master_id, mpi_comm_world,ierr)

		if(rid.ne.mpi_master_id)then
			call conv_int_real_arrays_s1d(s1d, intarr, realarr,logarr)
		end if
	end subroutine
	
	subroutine bcast_fc_mpi(fc)
		!use com_main_gw
		implicit none
		type(sts_fc_type)::fc
		integer::intarr(fc%type_int_size)
		integer nint, nreal, ierr
		real(8)::realarr(fc%type_real_size)
		logical::logarr(fc%type_log_size)

		!nint=
		!nreal=
		!allocate(intarr(nint), realarr(nreal))
		if(rid.eq.mpi_master_id)then
			call conv_fc_int_real_arrays(fc, intarr, realarr,logarr)
		end if
		call mpi_bcast(intarr, fc%type_int_size, MPI_INTEGER, mpi_master_id, mpi_comm_world,ierr)
		call mpi_bcast(realarr, fc%type_real_size, MPI_DOUBLE_PRECISION, mpi_master_id, mpi_comm_world,ierr)
		call mpi_bcast(logarr, fc%type_log_size, MPI_LOGICAL, mpi_master_id, mpi_comm_world,ierr)

		if(rid.ne.mpi_master_id)then
			call conv_int_real_arrays_fc(fc, intarr, realarr,logarr)
		end if
	end subroutine

	subroutine bcast_s2d_mpi(s2d)
		!use com_main_gw
		implicit none
		type(s2d_type)::s2d
		!integer nx, ny
		!real(8)::a2(nx,ny)
		!integer,allocatable::intarr(:)
		integer nreal, ierr
		!real(8),allocatable::realarr(:)

		nreal=s2d%nx*s2d%ny
		!print*, "rid, nreal=", rid, nreal
		call mpi_bcast(s2d%fxy, nreal, MPI_DOUBLE_PRECISION, mpi_master_id, mpi_comm_world,ierr)
		!print*, "bcast success"
	end subroutine
	subroutine bcast_dms_gxj(dm)
		implicit none
		type(diffuse_mspec)::dm
		integer i
		do i=1, dm%n
			call  bcast_s2d_mpi(dm%mb(i)%all%gxj)
		end do
	end subroutine
	subroutine get_dms(dm)
		!use com_main_gw
		implicit none
		type(diffuse_mspec)::dm
		integer ierr
        call mpi_BARRIER(mpi_comm_world, ierr)
        !call dms%input_bin("output/ecev/dms/dms_"//trim(adjustl(tmpssnapid))//"_"//trim(adjustl(tmpj)))
        !call dms%input_all_barge_bin("output/ecev/dms/barge_"//trim(adjustl(tmpssnapid))//"_"//trim(adjustl(tmpj)))		
        !print*, "rid input finished", rid
        !if(pid.ne.0)then
        !	call dms%mb(1)%barge%print("barge from file")
        !end if
        !print*, "bcast begin!"
        call bcast_dms_barge(dm)
        call bcast_dms_fden(dm)
        
        print*, "start get diffuse coefficients", rid	
		
!=================================		
        call dm_get_dc_mpi(dm)
        print*, "cal dms finished",rid
        call mpi_BARRIER(mpi_comm_world, ierr)

	end subroutine
	subroutine bcast_dms_asym_weights(dm)
		implicit none
		type(diffuse_mspec)::dm
		real(8),allocatable:: weights(:,:)
		integer ierr,i

		allocate(weights(6,dm%n))
		if(rid.eq.mpi_master_id)then
			do i=1, dm%n
				weights(1:6,i)=(/dms%mb(i)%all%asymp, dms%mb(i)%star%asymp, &
				dms%mb(i)%bstar%asymp,dms%mb(i)%sbh%asymp,dms%mb(i)%bbh%asymp,&
				dms%weight_asym/)
			end do
		end if
		call mpi_bcast(weights, 6*dm%n, MPI_DOUBLE_PRECISION, 0, mpi_comm_world,ierr)
		if(rid.ne.mpi_master_id)then
			do i=1, dm%n
				associate(mb=>dms%mb(i))
					mb%all%asymp=weights(1, i); mb%star%asymp=weights(2,i)
					mb%bstar%asymp=weights(3, i); mb%sbh%asymp=weights(4,i)
					mb%bbh%asymp=weights(5, i);
				end associate
			end do
            dms%weight_asym=weights(6,1)
		end if
	end subroutine

subroutine collection_data_single_bympi(smsa, n)
	!use com_main_gw
	implicit none 
	integer ierr, n, ntot, i
	type(particle_samples_arr_type)::smsa(n)

	!call bksams%sp(1)%get_length(ntot)
	!print*, "bksams%n=", ntot
	call update_arrays_single()
	call mpi_BARRIER(mpi_comm_world, ierr)
	!if(rid.ne.0)then
	!	print*, "collection: before collection"
	!	call get_memo_usage(proc_id)
	!end if
	if(rid.eq.mpi_master_id)then
			print*, "single:start_collect"
			!print*, "rid:bksams_arr_norm%n=",rid, bksams_arr_norm%n
			!read(*,*)
			call collect_to_root_sps_single(bksams_arr_norm, smsa,ctl%ntasks)
			!print*, bksams_arr%n
			print*, "single:end_collect"
	else
		!print*, "rid:bksams_arr_norm%n=",rid, bksams_arr_norm%n
		call send_particle_sample_arr_mpi(bksams_arr_norm,smsa(rid+1), rid, mpi_master_id)
	end if
	print*, "collection finished rid=",rid
end subroutine

end module
