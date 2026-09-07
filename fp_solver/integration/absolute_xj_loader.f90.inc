subroutine load_absolute_xj_samples(samples)
    use com_main_gw
    use, intrinsic :: ieee_arithmetic
    implicit none
    type(particle_samples_arr_type)::samples
    character(len=256)::fname
    character(len=4)::rankstr
    integer iu,ios,nfile,i
    real(8) x,j,importance,sum_importance,mean_importance
    logical have_file

    write(rankstr,'(I4.4)') rid
    fname='initial_xj_rank'//rankstr//'.dat'
    inquire(file=trim(fname),exist=have_file)
    if(.not.have_file)then
        if(rid.eq.0) print*,"ABSOLUTE_CLOSURE samples absent; native sampler retained"
        return
    end if
    open(newunit=iu,file=trim(fname),status='old',action='read',iostat=ios)
    if(ios.ne.0) error stop "ABSOLUTE_CLOSURE: cannot open xj sample file"
    read(iu,*,iostat=ios)nfile
    if(ios.ne.0.or.nfile.ne.samples%n)then
        print*,"ABSOLUTE_CLOSURE sample count",rid,nfile,samples%n
        error stop "ABSOLUTE_CLOSURE: xj sample count mismatch"
    end if
    sum_importance=0d0
    do i=1,samples%n
        read(iu,*,iostat=ios)x,j,importance
        if(ios.ne.0) error stop "ABSOLUTE_CLOSURE: truncated xj sample file"
        if(.not.ieee_is_finite(x).or.x.le.emin_factor.or.x.ge.emax_factor) &
            error stop "ABSOLUTE_CLOSURE: sample x outside energy grid"
        if(.not.ieee_is_finite(j).or.j.le.jmin_value.or.j.ge.jmax_value) &
            error stop "ABSOLUTE_CLOSURE: sample j outside angular-momentum grid"
        if(.not.ieee_is_finite(importance).or.importance.le.0d0) &
            error stop "ABSOLUTE_CLOSURE: invalid importance weight"
        samples%sp(i)%byot%a_bin=0.5d0*rh/x
        samples%sp(i)%byot%e_bin=sqrt(max(0d0,1d0-j*j))
        samples%sp(i)%weight_n=samples%sp(i)%weight_n*importance
        samples%sp(i)%weight_real=1d0
        call init_particle_sample_one_model_rnd(samples%sp(i),flag_ini_or)
        sum_importance=sum_importance+importance
    end do
    close(iu)
    mean_importance=sum_importance/dble(samples%n)
    ! The Python preflight verifies unit mean over the complete MPI ensemble.
    ! A rank-local unit-mean requirement would independently renormalize each
    ! shard and bias the physical DF by its finite-sample fluctuation.
    print*,"ABSOLUTE_CLOSURE_XJ_LOADED",rid,samples%n, &
        sum_importance,mean_importance
end subroutine
