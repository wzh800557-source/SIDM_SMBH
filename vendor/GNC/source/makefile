USE=MPI

#--global

COMPILE=gfortran
#COMPILE=ifort
use_hdf5=1
hdf5dirlib?=$(HDF5_LIBDIR)
hdf5dirinc?=$(HDF5_INCLUDEDIR)
performance_analyze=0
maximum_compatible=1
FC= mpif90
#FC= mpiifort
outdir:=./
ifeq ($(COMPILE), ifort)
	FLAG=  -O3 -traceback -check bounds -check uninit \
	      -check pointers -fpp  
	#FLAG+=   -f90=ifort
	ifeq ($(maximum_compatible), 1)
		FLAG+=-align none
	endif
    ifeq ($(use_hdf5),1)
    	FLAG+= -DHDF5
    	FLAG+= -I$(hdf5dirinc) -L$(hdf5dirlib) -lhdf5_fortran -lhdf5_hl_fortran -lhdf5 -lhdf5_hl
    endif
    ifeq ($(performance_analyze), 1)
		FLAG+= -g -DNDEBUG
    endif
else
	FLAG= -O3 -fbacktrace -fcheck=all -Waliasing  -fbounds-check   -cpp \
	    -malign-double
    ifeq ($(use_hdf5),1)
    	FLAG+= -DHDF5
    	FLAG+= -I$(hdf5dirinc) -L$(hdf5dirlib) -lhdf5_fortran -lhdf5_hl_fortran -lhdf5 -lhdf5_hl
    endif
    ifeq ($(performance_analyze), 1)
		FLAG+= -pg --param=max-vartrack-size=60000000 -g
	endif
    ifeq ($(performance_analyze), 2)
		FLAG+= -fprofile-arcs -ftest-coverage
	endif
endif
VPATH=com

OUT=libgnc.a

OBJS=fun_math.o APP_CONST.o sts_type_com.o \
sts_type_s1d_com.o sts_type_s2d_com.o md_fc.o intgl.o \
 FILE_STRING.o str.o gaussian.o FUN_RANDOM.o linear_interpolate.o \
 sts_fc.o ingl_dopri5.o app_celestial.o md_s1d_hst.o sts_basic_type.o

ifeq ($(use_hdf5),1)
OBJS+= file_md_hdf5_file.o file_md_hdf5_group.o file_md_hdf5.o\
file_md_hdf5_table.o
endif

OBJS+= events.o md_binary.o md_particle_sample.o 
OBJS+= md_bk_particle.o 
OBJS+=chain_pointer_type.o chain_type.o   coeff.o  dms.o md_cfuns.o modules.o 
OBJS+= sidm_params.o md_main_gw.o  stellar_obj.o
OBJS+= io_txt.o ini_single.o Sigma_funcs.o ini_single_user.o
OBJS+= EJ_funcs_ba16.o 
OBJS+= ac_ec_evl_single.o 
OBJS+= sams_sub.o sams_sub_single.o clone.o 
OBJS+= sams_sts.o sams_sts_single.o MPI_comu.o
OBJS+= read_ini_single.o  io_bin.o pro_single.o flux_diag.o gen_ge.o 
OBJS+=  io_hdf5.o 
OBJS+= load_df_from_table.o
#OBJS+=analy_solu_2D.o 
#4_body.o

ifeq ($(USE),MPI)
OBJS+=RR_mpi.o
else
OBJS+=RR_omp.o 
endif

#---make target
all: $(OUT) 

ifeq ($(use_hdf5),1)
file_md_hdf5.o: file_md_hdf5_table.o file_md_hdf5_group.o file_md_hdf5_file.o

md_fc.o: sts_type_s1d_com.o md_s1d_hst.o file_md_hdf5.o
sts_type_s1d_com.o: sts_basic_type.o  file_md_hdf5.o
else
md_fc.o: sts_type_s1d_com.o md_s1d_hst.o
sts_type_s1d_com.o: sts_basic_type.o 
endif

md_s1d_hst.o:  sts_type_s1d_com.o

sts_type_s2d_com.o: sts_type_s1d_com.o md_s1d_hst.o md_fc.o
sts_type_com.o : sts_type_s1d_com.o sts_type_s2d_com.o 

md_particle_sample.o : md_binary.o

md_sample.o:md_particle_sample.o

coeff.o: sts_type_com.o

chain_pointer_type.o: md_particle_sample.o

md_bk_particle.o: md_particle_sample.o coeff.o chain_type.o

steller_obj.o: sts_type_com.o coeff.o chain_type.o md_bk_particle.o 

dms.o: coeff.o stellar_obj.o

modules.o: events.o dms.o chain_type.o

md_cfuns.o: MPI_comu.o

MPI_comu.o : modules.o  

mobse_subs.o: mobse.o

chain_type.o: chain_pointer_type.o


md_main_gw.o: modules.o MPI_comu.o md_cfuns.o   \
	coeff.o events.o  sts_type_com.o

sidm_params.o: APP_CONST.o
md_main_gw.o: sidm_params.o

io_txt.o ini_single.o ac_ec_evl.o Sigma_funcs.o EJ_funcs_ba16.o \
sams_sub.o sams_sub_single.o  clone.o   sams_sts.o \
 read_ini_single.o  io_bin.o pro_single.o io_hdf5.o sams_sts_single.o flux_diag.o gen_ge.o \
 RR_mpi.o ac_ec_evl_single.o  sigma_funcs_aniso.o  ini_single_user: md_main_gw.o 

$(OUT):$(OBJS)
	ar -rc $(OUT) $(OBJS)
#	cp $(OUT)  $(outdir)
#	cp $(MODS) $(outdir)


%.o:%.f90
#   $(FC) -c -fpic $<
	$(FC) -c $(FLAG) $<

%.o:%.for
#   $(FC) -c -fpic $<
	$(FC) -c $(FLAG) $<

%.o:%.f
#   $(FC) -c -fpic $<
	$(FC) -c $(FLAG) $<

.PHONY:clean
clean:
	rm *.o
	rm *.mod
	rm *.a
	rm *.gcda
	rm *.gcno
	rm com/*.mod
	

