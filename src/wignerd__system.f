! WignerD: calculate the Wigner d and D matrices
! Copyright (C) 2025 Josh Forer <j.forer@posteo.net>
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
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
