! ================================================================================================================================ !
module wignerd__system
  !! Contains the definitions of stderr and procedures to stop the execution of the program while
  !! printing error messages.

  use wignerd__types,                only: rp
  use, intrinsic :: iso_fortran_env, only: error_unit

  implicit none

  private

  public die

  integer, parameter, public :: stderr = error_unit
    !! The file unit associated with standard error

  interface die
    !! Stop program execution with one or two messages
    module procedure die_1
    module procedure die_2
  end interface die

! ================================================================================================================================ !
contains
! ================================================================================================================================ !

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine die_1(message)
    !! Stop program execution with a message
    implicit none
    character(*), intent(in), optional :: message
    write(stderr,*)
    write(stderr,'("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")')
    write(stderr,'("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")')
    write(stderr,'("@@@@                                                           @@@@")')
    write(stderr,'("@@@@                          ERROR                            @@@@")')
    write(stderr,'("@@@@                                                           @@@@")')
    write(stderr,'("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")')
    write(stderr,'("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")')
    write(stderr,*)
    if(.not.present(message)) error stop
    write(stderr,'("STOP",X,"::",X,A)') message
    write(stderr,*)
    error stop
  end subroutine die_1
  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine die_2(message1,message2)
    !! Stop program execution with two messages
    implicit none
    character(*), intent(in) :: message1
    character(*), intent(in) :: message2
    write(stderr,*)
    write(stderr,'("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")')
    write(stderr,'("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")')
    write(stderr,'("@@@@                                                           @@@@")')
    write(stderr,'("@@@@                          ERROR                            @@@@")')
    write(stderr,'("@@@@                                                           @@@@")')
    write(stderr,'("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")')
    write(stderr,'("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")')
    write(stderr,*)
    write(stderr,'("STOP",X,"::",X,A,/,A)') message1, message2
    write(stderr,*)
    error stop
  end subroutine die_2

! ================================================================================================================================ !
end module wignerd__system
! ================================================================================================================================ !
