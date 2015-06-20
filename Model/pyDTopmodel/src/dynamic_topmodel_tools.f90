INCLUDE 'mkl_dss.f90' ! Include the standard DSS "header file."
subroutine update(recharge,storage,qout,qin,recharge1,storage1,qout1,qin1,&
                  area,dx,dt,celerity,celerity1,wvalues,wcolumns,wrowindex,&
                  nthreads,maxntt,nhsu,nvalues)

 use mkl_dss
 implicit none
 integer :: nhsu,nvalues
 real*8,intent(in) :: dt
 real*4,intent(in),dimension(nhsu) :: dx,area,celerity,celerity1
 real*4,intent(inout),dimension(nhsu) :: recharge,storage,recharge1,storage1,qout,qin,qout1,qin1
 integer*4,intent(in),dimension(nvalues) :: wcolumns
 integer*4,intent(in),dimension(nhsu+1) :: wrowindex
 real*4,intent(in),dimension(nvalues) :: wvalues
 integer,intent(in) :: nthreads,maxntt
 real*4,dimension(nhsu) :: scarea
 real*4 :: dt_minimum,dtt
 integer :: ntt,itime
 type(MKL_DSS_HANDLE) :: dss_handle ! Allocate storage for the solver handle.
 integer :: opt,error,perm(nhsu)

 !Set the number of threads
 call mkl_set_num_threads(nthreads)

 !Fix the number of threads
 call mkl_set_dynamic(1)

 ! Initialize the solver (Move eventually to beginning and end of entire simulation?)
 opt = MKL_DSS_MSG_LVL_WARNING + MKL_DSS_TERM_LVL_ERROR + MKL_DSS_SINGLE_PRECISION + MKL_DSS_ZERO_BASED_INDEXING
 error = dss_create(dss_handle,opt)

 ! Define the non-zero structure of the matrix.
 error = DSS_DEFINE_STRUCTURE(dss_handle,MKL_DSS_NON_SYMMETRIC,wrowindex,nhsu,nhsu,wcolumns,nvalues)

 ! Reorder the matrix.
 error = DSS_REORDER(dss_handle,MKL_DSS_DEFAULTS, perm )

 !Define the specific catchment area
 scarea = area/dx

 !Determine the required time step
 dt_minimum = minval(dx/celerity)
 ntt = int(ceiling(dt/dt_minimum))
 if (ntt .gt. maxntt) ntt = maxntt
 if (ntt .le. 0) ntt = 1
 dtt = dt/ntt
          
 !Process the smaller time steps
 do itime = 1,ntt
 
   !Solve the kinematic wave for this time step
   call solve_kinematic_wave(nhsu,nvalues,storage,qout,qin,recharge,storage1,qout1,qin1,&
                             recharge1,scarea,dx,dtt,celerity,celerity1,&
                             wvalues,wcolumns,wrowindex,dss_handle)

 enddo

 ! Finalize the solver
 error = DSS_DELETE(dss_handle,opt)

end subroutine update

subroutine solve_kinematic_wave(nhsu,nvalues,storage,qout,qin,recharge,storage1,qout1,qin1,&
                                recharge1,scarea,dx,dtt,celerity,celerity1,&
                                wvalues,wcolumns,wrowindex,dss_handle)
 use mkl_dss
 use omp_lib
 implicit none
 integer,intent(in) :: nhsu,nvalues
 real*4,intent(in) :: dtt
 real*4,intent(inout),dimension(nhsu) :: storage,qout,qin,recharge,storage1,qout1,qin1,recharge1
 real*4,intent(in),dimension(nhsu) :: scarea,dx,celerity,celerity1
 integer*4,intent(in),dimension(nvalues) :: wcolumns
 integer*4,intent(in),dimension(nhsu+1) :: wrowindex
 real*4,intent(in),dimension(nvalues) :: wvalues
 type(MKL_DSS_HANDLE),intent(inout) :: dss_handle ! Allocate storage for the solver handle.
 real*4,dimension(nvalues) :: A
 real*4,dimension(nhsu) :: denominator,numerator1,numerator2,part1,part2
 real*4 :: w
 real*4 :: alpha
 integer*4 :: i,j
 integer :: opt,error,perm(nhsu)
 real*8 :: t0,t1

 !Define implicit scheme parameters
 w = 0.5 
 alpha = 1.0

 !Calculate as many constants as possible before going into the iterations
 denominator = (1 + dtt*w*celerity/dx)
 numerator1 = qout1 + dtt*w*celerity*recharge + dtt*(1.0 - w)*celerity1*((qin1 - qout1)/dx + recharge1)
 numerator2 = dtt*w*celerity/dx
 part1 = numerator1/denominator
 part2 = numerator2/denominator

 !Prepare the A matrix (I - A) where A = ((F*(I*P2/SCAREA)).T*(I*SCAREA)).T
 do i = 1,nhsu
  A(wrowindex(i)+1:wrowindex(i+1)) = wvalues(wrowindex(i)+1:wrowindex(i+1))*&
                                     scarea(wcolumns(wrowindex(i)+1:wrowindex(i+1))+1)*&
                                     part2(i)/scarea(i)
  !A(wrowindex(i)+1:wrowindex(i+1)) = wvalues(wrowindex(i)+1:wrowindex(i+1))*scarea(i)*&
  !                                   part2(wcolumns(wrowindex(i)+1:wrowindex(i+1))+1)/&
  !                                   scarea(wcolumns(wrowindex(i)+1:wrowindex(i+1))+1)
 enddo

 !I - A
 i = 0
 do j = 1,nvalues
  if ((wcolumns(j) .eq. i) .and. (j .ge. wrowindex(i+1) + 1))then
   A(j) = 1 - A(j)
   i = i + 1
  else
   A(j) = - A(j)
  endif
 enddo
   
 ! Factor the matrix.
 error = DSS_FACTOR_REAL( dss_handle, MKL_DSS_DEFAULTS, A)

 ! Solve the system of linear equations
 error = DSS_SOLVE_REAL( dss_handle, MKL_DSS_DEFAULTS, part1, 1, qout)

 !Set all negative fluxes to 0
 where (qout < 0.0) qout = 0.0

 !Calculate qin
 !call mkl_cspblas_scsrgemv('T',nhsu,wvalues,wrowindex,wcolumns,scarea*qout,qin)
 call mkl_cspblas_scsrgemv('N',nhsu,wvalues,wrowindex,wcolumns,scarea*qout,qin)
 qin = qin/scarea

 !Adjust the storages
 storage = storage + dtt*((qin - qout)/dx + recharge)

 !Set the next time step's info
 qout1 = qout
 qin1 = qin

end subroutine



