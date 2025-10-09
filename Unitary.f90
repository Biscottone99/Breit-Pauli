program evolve_with_unitary
implicit none
    !> integer dim of matrix
    integer  dim
    !> hamiltonian (dim, dim)
    doubleprecision, allocatable :: eigenvalues(:)
    !> reduce density matrix (dim,dim)
    doublecomplex, allocatable :: rho(:,:)
    !> integer number of properties
    integer  n_prop
    !> property tensor (n_prop, dim, dim)
    doublecomplex, allocatable :: props(:,:,:)
    !> name of the properties  tensor (n_prop)
    character(len=405), allocatable :: name_props(:)
    !>doubleprecision timestep
    doubleprecision Deltat
       !> integer number of timesteps
    integer timesteps
    !> dimension of the krylov space
    integer nkrylov

    doublecomplex, allocatable :: gammas(:,:), Hdiag(:,:)
    doubleprecision, allocatable :: density(:,:)

    integer i,j
    doubleprecision, parameter :: e_SI=1.6021766D-19, kB_SI=1.380649D-23,kB_au=3.166811563d0*(10**(-6)), KB_ev = 8.617333d-5, h_SI=6.62607015D-34, h_ev=4.135667696D-15,&
       hbar_SI=1.054571817D-34, hbar_ev=6.582119569D-16, epsilon0_SI= 8.8541878128D-12, me_SI=9.1093897D-31, c_SI = 2.99892D8, c_au = 137.d0, epsilon0_ev=55.26349406D6, &
       pi=4.D0*DATAN(1.D0), e_nepero=2.7182818284590d0
    doublecomplex, parameter ::  imagine=(0.d0, 1.d0), reality=(1.d0, 0.d0)

    character(len=400) name
    
    read(*,*) dim
    allocate(rho(dim,dim), eigenvalues(dim), Hdiag(dim,dim))
    read(*,'(A400)') name
    open(1,file=name, form='unformatted', access='stream')
    read(1) rho
    close(1)
      
    name = ''
    read(*,'(A400)') name
    open(1,file=name, form='unformatted', access='stream')
    read(1) eigenvalues
    close(1)


!!$    do i = 1, dim
!!$       write(*,*) eigenvalues(i)
!!$    enddo
    read(*,*) n_prop
    
    if(n_prop.gt.0)then
            allocate(props(n_prop, dim, dim), name_props(n_prop))
            do i=1, n_prop
                name = ''
                read(*,'(A400)') name
                name_props(i) = name
                open(1, file=name, form='unformatted',access='stream')
                read(1) props(i,:,:)
                close(1)
            end do
        else
            allocate(props(0, dim, dim), name_props(0))
        end if
        
    read(*,*) Deltat
    read(*,*) timesteps
    
    call EvolveWithUnitary(rho, eigenvalues ,dim, Deltat, timesteps, n_prop, props,'*')

contains
 subroutine EvolveWithUnitary(rho, H, dim, Deltat, timesteps,n_prop, props,name)
        !> integer dimension of all the operator
        integer, intent(in) :: dim
        !> doublecomplex(dim, dim) density matrix to be evolved
        doublecomplex, intent(inout) :: rho(dim,dim)
        !> doublecomplex(dim,dim) Hamiltonian
        doubleprecision, intent(in) :: H(dim)
        !> doubleprecision timesteps in femtoseconds
        doubleprecision, intent(in) :: Deltat
        !> integer number of timesteps
        integer, intent(in) :: timesteps
        !>integer additional numbero f properties i want to calculate
        integer, intent(in) :: n_prop
        !> doublecomplex (n_prop, dim, dim) matrix of properties
        doublecomplex, intent(in) :: props(n_prop, dim, dim)
        !> do you want to write? if "*" write on output, if name, write on filename
        character(len=*), optional :: name

        doublecomplex, allocatable :: propts(:,:)
        character(len=50) :: spettro
        doublecomplex property
        integer i,j

        do i=1, n_prop
            spettro= ''
            write(spettro,*) 'prop', i,'.dat'
            call StripSpaces(spettro)
            open(i, file=spettro)
        end do
        write(*,*) 'Start unitary evolution'
        if(present(name))then
            write(*,*) 'Save at: ', name
            if(name.ne."*")then
                open(1, file=name, form='unformatted')
                do i=0, timesteps-1
                    write(1) rho(1:,1:)
                    call Unitary_Evolution(rho, H, dim, Deltat)
                end do
                write(1) rho(1:,1:)
                close(1)
            else
                do i=0, timesteps-1
                    write(*,"(A, F12.4,A,F12.4)") 'Time: ', Deltat*i, ',trace= ', real(trace_complex(rho, dim))
                    write(*,*) ' Population:'
                    write(*,'(<dim>(F7.3, X))') (real(rho(j,j)), j=1, dim)
                    if(n_prop.gt.0)then
                        allocate(propts(dim, dim))
                        do j=1, n_prop
                            propts = props(j,:,:)
                            call npmatmul_complex(propts, propts, rho, dim, dim, dim)
                            property = trace_complex(propts, dim)
                            write(*,'(A, I2)') 'Property number ', j
                            write(*,'(A, e18.6E5, 2X, e18.6E5)') "Value ", real(property), aimag(property)
                            write(j,'(F12.4, x, e18.6E5, x, e18.6E5)') Deltat*i, real(property), aimag(property)
                        end do
                        deallocate(propts)
                    end if
                    call Unitary_Evolution(rho, H, dim, Deltat)
                end do
                write(*,"(A, F12.4)") 'Time: ', Deltat*(i+1)
                call write_matrix_complex(rho, dim, dim, 'e11.3')
                write(*,*)
            end if
        else
            do i=0, timesteps-1
                call Unitary_Evolution(rho, H, dim, Deltat)
            end do
        end if
        write(*,*) 'Evolved with Unitary Evolution!'
        do  i=1, n_prop
            close(i)
        end do
    end subroutine EvolveWithUnitary
!>this subroutine compute the matrix product between 3 doublecomplex matrix
    !> @note parallelized with openmp
    subroutine npmatmul_complex(prod,h1, h2, a, b, c)
        integer, intent(in) :: a, b, c
        integer i, j, k
        doublecomplex :: h1(a,b), h2(b,c)
        doublecomplex, intent(inout) ::prod(a,c)
        doublecomplex, allocatable :: prod1(:,:)
        allocate(prod1(a,c))
        prod1 = (0.d0, 0.d0)
        !$OMP parallel do default(none), &
        !$OMP private(i,j,k), &
        !$OMP shared(h1,h2,a,b,c),&
        !$OMP reduction(+:prod1)
        do i=1, a
            do j=1, c
                do k=1, b
                    prod1(i,j) = prod1(i,j) + h1(i, k)*h2(k,j)
                enddo
            enddo
        enddo
        !$OMP end parallel do

        !$OMP parallel do default(none), &
        !$OMP private(i,j), &
        !$OMP shared(prod1, prod,a,c)
        do i=1, a
            do j=1, c
                prod(i,j) = prod1(i,j)
            end do
        end do
        !$OMP end parallel do
    endsubroutine npmatmul_complex
    !> this function computes the trace of a complex operator
    complex*16 function trace_complex(mat, dim)
        integer i
        !> integer dimension of the operator
        integer, intent(in) :: dim
        !> doublecomplex(dim,dim) operator
        complex*16, intent(in) :: mat(dim, dim)

        trace_complex = (0.d0, 0.d0)
        !$OMP parallel do default(none),&
        !$OMP private(i),&
        !$OMP shared(dim, mat), &
        !$omp reduction(+:trace_complex)
        do i=1, dim
            trace_complex = trace_complex + mat(i,i)
        enddo
        !$omp end parallel do
    end function trace_complex


subroutine write_matrix_complex(mat, a1, a2, fin)
      integer i,j,a1,a2
      complex*16 mat(a1, a2)
      character(len=50) f
      character(len=*), optional :: fin
      if(present(fin))then
         write(f,*) '(A1,'//fin//',A1,'//fin//',A1)'
      else
         f='(A1,f5.2, A1, f5.2, A1)'
      end if
      do i=1, a1
         do j=1, a2
           write(*,f, advance='no') '(',real(mat(i,j)),',',aimag(mat(i,j)),')'
         enddo
         write(*,*) ''
      enddo
      write(*,*)
    end subroutine write_matrix_complex
  
subroutine StripSpaces(string)				!thank you stackoverflow!
      character(len=*) :: string
      integer :: stringLen
      integer :: last, actual
      
      stringLen = len (string)
      last = 1
      actual = 1
      
      do while (actual < stringLen)
         if (string(last:last) == ' ') then
            actual = actual + 1
            string(last:last) = string(actual:actual)
            string(actual:actual) = ' '
         else
            last = last + 1
            if (actual < last) &
                 actual = last
         endif
      end do
      
    end subroutine StripSpaces
    
    subroutine Unitary_Evolution(rho,H, dim,timestep)
        !> integer dimension of the operators
        integer, intent(in) :: dim
        !> complex operator or complex matrix of dimension(dim, dim) to be evolvedù
        doublecomplex, intent(inout) :: rho(dim,dim)
        !> doubleprecision timestep
        doubleprecision, intent(in) :: timestep
        !> doublecomplex hamiltonian of dimension (dim, dim)
        doubleprecision, intent(in) :: H(dim)

        
        integer i,j, i_time
        doubleprecision hbar_evfs, omega_ij

        hbar_evfs = hbar_ev*1.d15
        do i=1,dim
            do j=1,dim
                omega_ij = (H(i)-H(j))/hbar_evfs
                rho(i,j) = rho(i,j)*zexp(-imagine*omega_ij*timestep)
            end do
        end do
    end subroutine Unitary_Evolution


end program evolve_with_unitary
