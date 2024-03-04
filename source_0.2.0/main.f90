program main

use precision
use varsub
use constants
use co2cool

    implicit none
!    character(len=*), parameter :: infile = 'input.dat'
!    character(len=*), parameter :: outfile = 'output.dat'
!    character(len=*), parameter :: infile = 'input_xgrid_co2_3.dat'
!    character(len=*), parameter :: outfile = 'output_xgrid_co2_3_gx.dat'
!    integer, parameter :: nlev=71

!    character(len=*), parameter :: infile  = 'input_par_s2b_24_fr.dat' 
!    character(len=*), parameter :: outfile = 'output_par_s2b_24_fr.dat'
!    integer, parameter :: nlev=273

!    character(len=*), parameter :: infile  = 'input_par_s2b_24.dat' 
!    character(len=*), parameter :: outfile = 'output_par_s2b_24.dat'
!    integer, parameter :: nlev=141

    character(len=*), parameter :: infile  = 'input_par_s2b_24_fri.dat' 
    character(len=*), parameter :: outfile = 'output_par_s2b_24_fri.dat'
    integer, parameter :: nlev=273

    real(dp), dimension(nlev) :: pres, X, temp, co2vmr, ovmr, o2vmr, n2vmr
    real(dp), dimension(nlev) :: cr_new
    real(dp) :: surf_temp
    integer :: lev0, i

    ! Read input data
    !----------------
    open(unit=10, file=infile,action = 'read')
    do i = 1, 8 
       read(10,*)
    end do
    read(10,*) lev0, surf_temp
!    do i = nlev,1,-1  
    do i = 1, nlev
       read(10,*)pres(i), temp(i), co2vmr(i), ovmr(i), o2vmr(i), n2vmr(i)
    end do 
    cr_new=0.0_dp 
    X = log(p0/pres)  

    ! Set surf_temp if provided
    if (surf_temp > 0.0) then
        print *, 'Set surf_temp: ', surf_temp
    else
        if (pres(1) > pres(2)) surf_temp = temp(1)
        if (pres(2) > pres(1)) surf_temp = temp(nlev)
   end if

write(*,*) 'Surf_temp ',surf_temp,lev0

    ! Calculate new parameters
    print *, ' ####         start : ',seconds(),' s'
    do i = 1,1 !000 
    call NEW_PARAM(temp, pres, co2vmr, ovmr, o2vmr, n2vmr, lev0, &
                      surf_temp, cr_new)
    enddo
    print *, ' ####         end : ',seconds()/1000,' s'

    ! Write output
    open(unit=10, file=outfile, status='replace', action='write')
    do i = 1, size(pres)
        write(10, '(f10.4, 1x, E12.4, 1x, f12.4)') X(i), pres(i), cr_new(i)
    end do
    close(10)

    contains

end program main

