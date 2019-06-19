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

! reader: Tools for reading .wfn and .xyz files


module reader
  use param, only: mline
  implicit none
  private
  
  public :: readfile
  public :: molecule, ifile_xyz, ifile_wfn, ifile_grd, ifile_wfx
  public :: grid1, readgrid,readxyz, reference, read_integers
  public :: read_reals1

  !> Radial grid type.
  type grid1
     logical :: init !< Is initialized?
     real*8 :: a !< Logarithmic grid parameter ri = a * exp(b * (i-1))
     real*8 :: b !< Logarithmic grid parameter ri = a * exp(b * (i-1))
     real*8 :: rmax !< Max. grid distance
     real*8 :: rmax2 !< Squared max. grid distance
     integer :: ngrid !< Number of nodes
     real*8, allocatable :: r(:) !< Node positions
     real*8, allocatable :: f(:) !< Grid values, f = 4*pi*r^2*rho
     real*8, allocatable :: fp(:) !< First derivative of f
     real*8, allocatable :: fpp(:) !< Second derivative of f 
  end type grid1

  !> Molecule type
  type molecule
     ! file name
     character(mline) :: name
     ! common for xyz and wfn
     integer :: ifile
     integer :: n
     real*8, allocatable :: x(:,:)
     integer, allocatable :: z(:)
     integer, allocatable :: q(:)
     ! fragment id
     integer, allocatable :: ifrag(:)
     ! only for wfn
     integer :: nmo, npri
     integer, allocatable :: icenter(:) 
     integer, allocatable :: itype(:)
     integer :: ntyp(35), maxntyp
     integer, allocatable :: intyp(:)
     real*8, allocatable :: e(:)
     real*8, allocatable :: occ(:)
     real*8, allocatable :: c(:,:) 
     real*8 :: mult  
     integer :: nelec
  end type molecule
  integer, parameter :: ifile_xyz = 1
  integer, parameter :: ifile_grd = 2
  integer, parameter :: ifile_wfn = 3
  integer, parameter :: ifile_wfx = 4


! Roberto 
!> Reference type
 type reference
      real*8 :: xinit(3),xmax(3) 
      integer :: nstep(3)
      real*8, allocatable, dimension(:,:,:) :: cgrad,crho 
      real*8, allocatable, dimension(:,:,:) :: celf,cxc,ctp
      real*8, allocatable, dimension(:,:,:,:) :: cheigs 
      real*8 :: rho,dimgrad
 endtype reference  
contains

  function readfile(file) result(m)
    use param
    use tools_io
    implicit none
    
    character*(mline), intent(in) :: file
    type(molecule) :: m
   
    character*(mline) :: ext

    ext = file(index(file,'.',.true.)+1:)
    call upper(ext)
    select case(adjustl(trim(ext)))
    case('XYZ')
       m = readxyz(file)
    case('WFN')
       m = readwfn(file) 
    case('WFX')
       m=readwfx(file)
    case default
       call error('readfile','Not recognized molecular format',faterr)
    end select
    
  end function readfile

  !> Read xyz file
  function readxyz(file) result(m)
    use param
    use tools_io
    implicit none

    character*(mline), intent(in) :: file
    type(molecule) :: m

    integer :: istat, i, lp
    integer, parameter :: lu = 10
    character(mline) :: name, line
    logical :: ok

    ! set type
    m%ifile = ifile_xyz
    m%name = file
    ! read number of atoms
    open(lu,file=file,status='old')
    read(lu,*) m%n

    ! allocate arrays
    allocate(m%x(3,m%n),stat=istat)
    if (istat /= 0) call error('readxyz','could not allocate memory for atomic positions',faterr)
    allocate(m%z(m%n),stat=istat)
    if (istat /= 0) call error('readxyz','could not allocate memory for atomic numbers',faterr)
    allocate(m%q(m%n),stat=istat)
    if (istat /= 0) call error('readxyz','could not allocate memory for atomic charges',faterr)
    allocate(m%ifrag(m%n),stat=istat)
    if (istat /= 0) call error('readxyz','could not allocate memory for fragment ids',faterr)

    ! read atomic coordinates and atomic numbers
    read(lu,*)
    do i = 1, m%n
       read (lu,'(a)') line
       line = trim(line) // char(0)
       lp = 1
       name = getword(name,line,lp)
       ok = isreal(m%x(1,i),line,lp)
       ok = ok .and. isreal(m%x(2,i),line,lp)
       ok = ok .and. isreal(m%x(3,i),line,lp)
       if (.not.ok) call error('readxyz','error reading xyz file',faterr)
       m%z(i) = zatguess(name)
       if (m%z(i) == -1) call error('readxyz','atom type not recognized: '//trim(name),faterr)
       ok = isinteger(m%q(i),line,lp)
       if (.not.ok) m%q(i) = 0
       m%x(:,i) = m%x(:,i) / bohrtoa
    end do
       m%nelec=sum(m%z)+sum(m%q)
    ! close
    close(lu)

  end function readxyz

  !> Read wfn file
  function readwfn(file) result(m)
    use param
    use tools_io
    implicit none

    character*(mline), intent(in) :: file
    type(molecule) :: m

    integer, parameter :: lu = 10

    character*4 :: orbtyp
    integer :: i, j, istat, imax, icount
    real*8 :: zreal

    ! set type
    m%ifile = ifile_wfn
    m%name = file

    ! read number of atoms, primitives, orbitals
    open(lu,file=file,status='old')
    read (lu,*)
    read (lu,101) orbtyp, m%nmo, m%npri, m%n

    ! atomic positions and numbers
    allocate(m%x(3,m%n),stat=istat)
    if (istat /= 0) call error('readwfn','could not allocate memory for atomic positions',faterr)
    allocate(m%z(m%n),stat=istat)
    if (istat /= 0) call error('readwfn','could not allocate memory for atomic numbers',faterr)
    allocate(m%q(m%n),stat=istat)
    if (istat /= 0) call error('readwfn','could not allocate memory for atomic charges',faterr)
    allocate(m%ifrag(m%n),stat=istat)
    if (istat /= 0) call error('readwfn','could not allocate memory for fragment ids',faterr)
    do i = 1, m%n
       read(lu,106) m%x(:,i), zreal
       m%z(i) = nint(zreal)
       m%q(i) = 0
    end do

    ! center assignments, types of primitives
    allocate(m%icenter(m%npri),stat=istat)
    if (istat /= 0) call error('readwfn','could not allocate memory for icenter',faterr)
    allocate(m%itype(m%npri),stat=istat)
    if (istat /= 0) call error('readwfn','could not allocate memory for intyp',faterr)
    read(lu,102) (m%icenter(i),i=1,m%npri)
    read(lu,102) (m%itype(i),i=1,m%npri)

    ! primitive exponents
    allocate(m%e(m%npri),stat=istat)
    if (istat /= 0) call error('readwfn','could not allocate memory for exponents',faterr)
    read(lu,103) (m%e(i),i=1,m%npri)
    
    ! occupations and orbital coefficients
    allocate(m%occ(m%nmo),stat=istat)
    if (istat /= 0) call error('readwfn','could not allocate memory for occupations',faterr)
    allocate(m%c(m%nmo,m%npri),stat=istat)
    if (istat /= 0) call error('readwfn','could not allocate memory for orbital coefficients',faterr)
    do i = 1, m%nmo
       read(lu,104) m%occ(i)
       read(lu,105) (m%c(i,j),j=1,m%npri)
    end do
    m%nelec=sum(m%occ)
    close(lu)

    ! order by primitive type
    imax = 0
    do i = 1, 35
       imax = max(count(m%itype == i),imax)
       if (count(m%itype == i) == 0) exit
       m%maxntyp = i
    enddo
    allocate(m%intyp(m%npri))
    m%ntyp = 0
    icount = 0
    do i = 1, m%maxntyp
       do j = 1, m%npri
          if (m%itype(j) == i) then
             icount = icount + 1
             m%ntyp(i) = m%ntyp(i) + 1
             m%intyp(icount) = j
          end if
       enddo
    enddo

101 format (4X,A4,10X,3(I5,15X))
102 format(20X,20I3)
103 format(10X,5E14.7)
104 format(35X,F12.8)
105 format(5(E16.8))
106 format(24X,3(F12.8),10X,F5.1)

  end function readwfn 

  !> Read wfx file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	WFX READER, ARIAS	!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function readwfx(file) result(m)
    use param
    use tools_io
    implicit none

    character*(mline), intent(in)   :: file
    type(molecule)                  :: m
    integer, parameter :: lu = 10
    integer            :: i,j, istat    !counter
    character*(mline)  :: line    !pivot to search
    real*8             :: zreal 
    integer :: imax,icount 

    !set type
     m%ifile = ifile_wfx
     m%name = file

     !open the file
     open(lu,file=file,status='old')
     !start the reader

     do 
        read(lu,'(a)') line
        if (line == '<Number of Nuclei>') then          !READ THE NUCLEI
                read(lu,*) m%n
                        if (m%n == 0) call error('readwfx','Number of Nuclei is equal to zero',faterr)  
                        !Roberto
                        allocate(m%z(m%n), stat=istat)  
                        allocate(m%q(m%n))  
                        do i=1,m%n 
                           m%z(i)=0  !Initialize atomic numbers
                           m%q(i)=0  !Initialize atomic charges
                        enddo 

        elseif (line == '<Number of Occupied Molecular Orbitals>') then !READ THE MOL. ORB
                        read(lu,*) m%nmo
                        if (m%nmo == 0) call error('readwfx','Number of Occupied MO is equal to zero',faterr) 

        elseif (line == '<Atomic Numbers>') then    !ARRAY OF ATOMIC NUMBERS
                        if (.not. allocated (m%z)) allocate(m%z(m%n), stat=istat)
                        if (istat /= 0) call error('readwfx','could not allocate memory for Atomic Numbers',faterr) 

                        do i=1, m%n
                           read(lu,*) m%z(i) 
                          ! Roberto Check
                           if (m%z(i) == 0) call error('readwfx','Any atomic number is equal to zero',faterr)
                        enddo 

        elseif (line == '<Nuclear Charges>') then        !ARRAY OF NUCLEAR CHARGES
                        if (.not. allocated (m%q)) allocate(m%q(m%n), stat=istat)
                        do i=1, m%n
                           m%q(i) = 0  
                        enddo
                        allocate(m%ifrag(m%n))       
                 !NEW LINE, FOR ERROR IN RUN
        elseif (line == '<Nuclear Cartesian Coordinates>') then   !ARRAY OF CARTESIAN COORD
                        allocate(m%x(3,m%n),stat=istat)
                        if (istat /= 0) call error('readwfx','could not allocate memory for Cartesian Coordinates',faterr)
                        do i=1, m%n
                           read(lu,*) m%x(:,i)
                        enddo

        elseif (line == '<Number of Primitives>') then
                        read(lu,*) m%npri 
                        if (m%npri == 0) call error('readwfx','Number of Primitives is equal to zero',faterr) 

        elseif (line == '<Primitive Centers>') then      !ARRAY OF PRIMITIVES
                        allocate(m%icenter(m%npri),stat=istat)
                        if (istat /= 0) call error('readwfx','could not allocate memory for Primitives Centres',faterr)
                        read(lu,*) (m%icenter(i),i=1,m%npri) 

        elseif (line == '<Primitive Types>') then     !TYPES OF PRIM
                       allocate(m%itype(m%npri), stat=istat)
                       if (istat /= 0) call error('readwfx','could not allocate memory for Primitives Types',faterr)
                       read(lu,*) (m%itype(i),i=1,m%npri)  

        elseif (line == '<Primitive Exponents>') then  !EXPONENTS
                        allocate(m%e(m%npri), stat=istat)
                        if (istat /= 0) call error('readwfx','could not allocate memory for Primitive Exponents',faterr)
                        read(lu,*) (m%e(i),i=1,m%npri) 

        elseif (line == '<Molecular Orbital Occupation Numbers>') then !OCC NUMBER OF EACH MO
                        allocate(m%occ(m%nmo), stat=istat)
                        if (istat /= 0) call error('readwfx','could not allocate memory for Occupation Numbers',faterr)
                        read(lu,*) (m%occ(i),i=1,m%nmo)      

        elseif (line == '<Molecular Orbital Primitive Coefficients>') then  !COEF OF EACH PRIM OF EACH MO
                        allocate(m%c(m%nmo,m%npri), stat=istat)
                        if (istat /= 0) call error('readwfx','could not allocate memory for Primitive Coefficients',faterr)
                        do i=1, m%nmo
                           read(lu,*)
                           read(lu,*)
                           read(lu,*)                            !MUST BE READ 3 LINES, CORRECT FORMAT
                           read(lu,*) (m%c(i,j),j=1,m%npri) 
                        enddo
        endif 

        if (line == '<Virial Ratio (-V/T)>') then       !END OF FILE
           exit
        endif

     enddo

    
    ! Roberto. Check m has all the information.  
     do i=1, m%n    
          if (m%z(i) == 0) call error('readwfx','Any atomic number is equal to zero',faterr)
     enddo  

     if (m%n == 0) call error('readwfx','Number of Nuclei is equal to zero',faterr) 
     if (m%nmo == 0) call error('readwfx','Number of Occupied MO is equal to zero',faterr)

    ! Roberto.  Order by primitive type.
    imax = 0
    do i = 1, 35
       imax = max(count(m%itype == i),imax)
       if (count(m%itype == i) == 0) exit
       m%maxntyp = i
    enddo
    allocate(m%intyp(m%npri))
    m%ntyp = 0
    icount = 0
    do i = 1, m%maxntyp
       do j = 1, m%npri
          if (m%itype(j) == i) then
             icount = icount + 1
             m%ntyp(i) = m%ntyp(i) + 1
             m%intyp(icount) = j
          end if
       enddo
    enddo
 
    close(lu)

end function readwfx 


  !> Read a density core file from the database. 
  subroutine readgrid(g,z)
    use tools_io
    use param

    type(grid1), intent(out) :: g !< Output radial grid
    integer, intent(in) :: z !< Number of electrons to be read

    integer :: i, j
    real*8 :: xmin, zz, dx, r, r1, r2, r3 ,r4, delta, delta2
    integer :: ngrid, ns, ic
    integer, allocatable :: occ(:)
    real*8, allocatable :: wfcin(:), rr(:,:)
    character*2, allocatable :: wfcl(:)
    character*(mline) :: econf
    integer :: nn
    logical :: exist
    character*(mline) :: file

    integer, parameter :: lu = 10
    real*8, parameter :: core_cutdens = 1d-12 !< Cutoff contribution for core radial grids
    ! radial grid derivation formulas
    integer, parameter :: noef(6,3) = reshape((/&
       0,  1,  2,  3,  4,  5,&
       -2, -1,  0,  1,  2,  3,&
       -5, -4, -3, -2, -1,  0/),shape(noef)) !< Node offsets for 6-point derivation formulas.
    real*8, parameter :: coef1(6,3) = reshape((/&
       -274,  600, -600,  400,  -150,  24,&
       6,  -60,  -40,  120,   -30,   4,&
       -24, 150, -400,  600,  -600, 274 /),shape(coef1)) !< Coefficients for first derivative.
    real*8, parameter :: coef2(6,3) = reshape((/&
       225, -770,  1070,  -780,   305,   -50,&
       -5,   80,  -150,    80,    -5,     0,&
       -50,  305,  -780 ,  1070,  -770,   225/),shape(coef2)) !< Coefficients for second derivative.
    real*8, parameter :: fac1=1d0/120d0 !< Prefactor for first derivative.
    real*8, parameter :: fac2=2d0/120d0 !< Prefactor for second derivative.

    ! the filename
    file = trim(nciplot_dat) // nameguess(z) // "_lda" // ".wfc"
    ! check that the file exists
    inquire(file=file,exist=exist)
    if (.not.exist) then
       write (uout,'("File: ",A)') trim(file)
       call error("readgrid","Atomic density file not found",warning)
       g%init = .false.
       return
    end if

    ! Read header and allocate arrays
    open (unit=lu,file=file,status='old')
    read (lu,*) nn
    
    allocate(wfcin(nn),wfcl(nn),occ(nn))
    occ = 0
    read (lu,*) (wfcl(i),i=1,nn)
    read (lu,*) (occ(i),i=1,nn)
    read (lu,*) xmin, zz, dx, ngrid

    if (sum(occ) /= z) then
       ns = 0
       do i = 1, nn
          if (ns + occ(i) > z) then
             occ(i) = z - ns
             occ(i+1:nn) = 0
             exit
          else
             ns = ns + occ(i)
          end if
       end do
    end if

    ! Read the grid and build the density
    allocate(g%r(ngrid),rr(ngrid,0:2))
    rr = 0d0
    do i = 1, ngrid
       read (lu,*) r, (wfcin(j),j=1,nn)
       g%r(i) = r
       rr(i,0) = dot_product(occ(1:nn),wfcin(1:nn)**2)
       if (rr(i,0)/(4d0*pi*r**2) < core_cutdens .and. i > 1) then
          ngrid = i
          call realloc(g%r,ngrid)
          exit
       end if
    end do

    ! fill grid info
    g%init = .true.
    g%a = exp(xmin) / zz
    g%b = dx
    g%ngrid = ngrid
    g%rmax = g%r(ngrid)
    g%rmax2 = g%r(ngrid)**2

    ! calculate derivatives
    allocate(g%f(ngrid),g%fp(ngrid),g%fpp(ngrid))
    do i = 1, ngrid
       if (i <= 2) then
          ic = 1
       else if (i >= ngrid-2) then
          ic = 3
       else
          ic = 2
       end if
       do j = 1, 6
          rr(i,1) = rr(i,1) + coef1(j,ic) * rr(i+noef(j,ic),0)
          rr(i,2) = rr(i,2) + coef2(j,ic) * rr(i+noef(j,ic),0)
       end do
       rr(i,1) = rr(i,1) * fac1
       rr(i,2) = rr(i,2) * fac2

       r = g%r(i)
       r1 = 1d0 / r
       r2 = r1 * r1
       r3 = r2 * r1
       r4 = r3 * r1
       delta=1.d0/g%b
       delta2=delta*delta
       
       g%f(i) = rr(i,0) * r2
       g%fp(i) = (rr(i,1) * delta - 2.d0 * rr(i,0)) * r3
       g%fpp(i) = (rr(i,2) * delta2 - 5.d0 * rr(i,1) * delta + 6.d0 * rr(i,0)) * r4
    end do
    g%f = g%f / (4d0*pi)
    g%fp = g%fp / (4d0*pi)
    g%fpp = g%fpp / (4d0*pi)

    ! close the density file
    close(lu)

    econf = ""
    do i = 1, nn
       write (econf,'(A,X,A2)') trim(econf), wfcl(i)
       if (occ(i) < 10) then
          write (econf,'(A,"(",I1,")")') trim(econf), occ(i)
       else
          write (econf,'(A,"(",I2,")")') trim(econf), occ(i)
       end if
    end do

    write (uout,'("+ Read density file : ", A)') trim(file)
    write (uout,'("  Log grid (r = a*e^(b*x)) with a = ",1p,E12.4,", b = ",E12.4)') g%a, g%b
    write (uout,'("  Num. grid points = ",I5, ", rmax (bohr) = ",F10.4)') g%ngrid, g%rmax
    write (uout,'("  Integrated charge = ",F20.8)') sum(g%f * g%r**3 * g%b * 4d0 * pi)
    write (uout,'("  El. conf. : ",A)') trim(econf)
    write (uout,*)

    ! cleanup
    deallocate(rr,wfcin,wfcl,occ)

  end subroutine readgrid 


 function read_integers(lu,n) result(x)
    use param
    use tools_io
    implicit none 

   integer, intent(in) :: lu, n
   integer :: x(n)

   integer :: kk, lp, idum
   character*(mline) :: line

   kk = 0
   lp = 1
   read(lu,'(A)',end=999) line
   do while(.true.)
      if (.not.isinteger(idum,line,lp)) then
         lp = 1
         read(lu,'(A)',end=999) line
         line = adjustl(line)
         if (line(1:2) == "</") exit
      else
         kk = kk + 1
         if (kk > n) call error("read_integers","exceeded size of the array",2)
         x(kk) = idum
      endif
   enddo
   
   return
999 call error("read_integers","unexpected end of file",2)

 endfunction read_integers


 function read_reals1(lu,n) result(x)
    use param
    use tools_io
    implicit none

   integer, intent(in) :: lu, n
   real*8 :: x(n)

   integer :: kk, lp
   real*8 :: rdum
   character*(mline) :: line

   kk = 0
   lp = 1
   read(lu,'(A)',end=999) line
   do while(.true.)
      if (.not.isreal(rdum,line,lp)) then
         lp = 1
         read(lu,'(A)',end=999) line
         line = adjustl(line)
         if (line(1:1) == "<") exit
      else
         kk = kk + 1
         if (kk > n) call error("read_reals1","exceeded size of the array",2)
         x(kk) = rdum
      endif
   enddo
   
   return
999 call error("read_reals1","unexpected end of file",2) 

 endfunction read_reals1 


end module reader
