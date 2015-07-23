module bank_header

  implicit none
  
  integer, parameter :: BANK_REALLOC_SIZE = 3

!===============================================================================
! BANK is used for storing fission sites in eigenvalue calculations. Since all
! the state information of a neutron is not needed, this type allows sites to be
! stored with less memory
!===============================================================================

  type Bank
    ! The 'sequence' attribute is used here to ensure that the data listed
    ! appears in the given order. This is important for MPI purposes when bank
    ! sites are sent from one processor to another.
    sequence

    real(8)    :: wgt    ! weight of bank site
    real(8)    :: xyz(3) ! location of bank particle
    real(8)    :: uvw(3) ! diretional cosines
    real(8)    :: E      ! energy
  end type Bank
  
  type BankArray
    type(Bank), allocatable :: neutron(:)
    integer :: count = 0 ! The number of neutrons banked
    integer :: size  = 0 ! The number of neutrons that can be banked without reallocation
  end type BankArray
  
contains

  subroutine bank_neutron(this, wgt, xyz, uvw, E)
    type(BankArray), intent(inout) :: this
    real(8), intent(in)            :: wgt
    real(8), intent(in)            :: xyz(3)
    real(8), intent(in)            :: uvw(3)
    real(8), intent(in)            :: E
    
    type(Bank), allocatable :: temporary(:)
    
    if (this % count == this % size) then
      if (this % size == 0) then
        ! We need to create the bank array
        allocate(this % neutron(BANK_REALLOC_SIZE))
      else
        ! We need to expand the bank array
        call move_alloc(this % neutron, temporary)
        allocate(this % neutron(this % size + BANK_REALLOC_SIZE))
        this % neutron(1:this % size) = temporary
        deallocate(temporary)
      end if
      this % size = this % size + BANK_REALLOC_SIZE
    end if
    
    ! Copy data into bank
    this % count = this % count + 1
    this % neutron(this % count) % wgt = wgt
    this % neutron(this % count) % xyz = xyz
    this % neutron(this % count) % uvw = uvw
    this % neutron(this % count) % E = E
    
  end subroutine bank_neutron

end module bank_header
