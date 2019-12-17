! Copyright (c) 2013 Alberto Otero de la Roza <aoterodelaroza@ucmerced.edu>,
! Julia Conteras-Garcia <julia.contreras.garcia@gmail.com>,
! Erin R. Johnson <ejohnson29@ucmerced.edu>, and Weitao Yang
! <weitao.yang@duke.edu>
!
! nciplot is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

! Param: Parameters and enviroment variables

module param
   implicit none

   public
   integer, parameter :: mline = 2048 !< line length

   ! program limits
   integer, parameter :: atomic_zmax = 18

   ! paths
   character*(mline) :: nciplot_home
   character*(mline) :: nciplot_dat

   ! logical units
   integer, parameter :: stderr = 0 !< standard error lu
   integer, parameter :: stdin = 5 !< standard input lu
   integer, parameter :: stdout = 6 !< standard output lu
   integer :: uin, uout, uerr

   ! math
   real*8, parameter :: pi = 3.14159265358979323846d0 !< pi
   real*8, parameter :: const = 2.D0*(3.D0*PI**2)**(1.D0/3.D0)
   real*8, parameter :: eps = epsilon(1d0)

   ! physics
   real*8, parameter :: bohrtoa = 0.52917720859 !< bohr to angstrom conversion factor (nist2006)

   ! error types
   integer, parameter :: faterr = -1 !< fatal error flag
   integer, parameter :: warning = 1 !< warning flag
   integer, parameter :: noerr = 0   !< info flag
   integer :: nwarns = 0
   integer :: ncomms = 0

contains

   subroutine param_init()

      integer :: isenv, nn

      uin = stdin
      uout = stdout
      uerr = stderr
      call get_environment_variable("NCIPLOT_HOME", nciplot_home, status=isenv)
      if (isenv /= 0) then
         write (uout, '("(!) WARNING ")')
         write (uout, '("(!) The environment variable NCIPLOT_HOME is not set")')
         write (uout, '("(!) The library of atomic density grids is not available.")')
         write (uout, '("(!) WARNING ")')
         nciplot_home = "./"
      end if
      nn = len(trim(adjustl(nciplot_home)))
      nciplot_dat = trim(adjustl(nciplot_home))
      if (nciplot_home(nn:nn) /= "/") then
         nciplot_dat = trim(nciplot_dat)//"/"
      end if
      nciplot_dat = trim(nciplot_dat)//"dat/"

   end subroutine param_init

end module param
