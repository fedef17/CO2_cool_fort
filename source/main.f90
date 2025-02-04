program main

use precision
use varsub
use constants
use co2cool

    implicit none
    character(len=*), parameter :: infile = 'input.dat'
    character(len=*), parameter :: outfile = 'output.dat'
    real(dp), dimension(:), allocatable :: pres, X, temp, co2vmr, ovmr, o2vmr, n2vmr
    real(dp), dimension(:), allocatable :: cr_new
    real(dp) :: surf_temp
    integer :: lev0, i, nlev

    ! Read input data
    !----------------
    open(unit=10, file=infile,action = 'read')
    do i = 1, 8 
       read(10,*)
    end do
    read(10,*) nlev, lev0, surf_temp
    do i = 1, 2
       read(10,*)
    end do
    allocate(pres(nlev), X(nlev), temp(nlev), &
        co2vmr(nlev), ovmr(nlev), o2vmr(nlev), n2vmr(nlev), cr_new(nlev))    
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

    ! Calculate new parameters
!    print *, ' ####         start : ',seconds(),' s'
!    do i = 1,1000 
    call CO2_NLTE_COOL(temp, pres, co2vmr, ovmr, o2vmr, n2vmr, lev0, &
                      surf_temp, cr_new)
!    enddo
!    print *, ' ####         end : ',seconds()/1000,' s'

    ! Write output
    open(unit=10, file=outfile, status='replace', action='write')
    write(10, '(a2, a8, 1x, a12, 1x, a20)') "##", "X", "Pres. (hPa)", "Heating rate (K/day)"
    do i = 1, size(pres)
        write(10, '(f10.4, 1x, E12.4, 1x, f12.4)') X(i), pres(i), cr_new(i)
    end do
    close(10)

    deallocate(pres, X, temp, co2vmr, ovmr, o2vmr, n2vmr, cr_new)    

    contains

end program main

