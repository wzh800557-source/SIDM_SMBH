program cfuns
    use com_main_gw
    implicit none
    !type(cfuns_type)::cfs
    character*(200) str_num_bin,str_jmin, str_jmax, str_output_file
    integer num_bin
    real(8) logjmin, logjmax

    call getarg(1,str_num_bin)
    read(unit=str_num_bin, fmt="(I5)") num_bin

    call getarg(2,str_jmin)
    read(unit=str_jmin, fmt=*)  logjmin
    call getarg(3,str_jmax)
    read(unit=str_jmax, fmt=*)  logjmax
    call getarg(4, str_output_file)

    if(len(str_output_file).eq.0)then 
        str_output_file="../common_data/"
    end if
    call cfs%init(num_bin,num_bin, logjmin,logjmax)

    print*, "start cfuns ..., nx, ny=", cfs%nj, cfs%ns 
    call init_mpi()
    if(mod(num_bin,ctl%ntasks).ne.0)then
        print*, "error! mpi_threads should be integer times of num_bin!", num_bin, ctl%ntasks
        stop
    end if

    call cfs%get_cfs_s_mpi()
    if(rid.eq.0)then
        call cfs%output_bin(trim(adjustl(str_output_file)))
    end if
    call stop_mpi()
end