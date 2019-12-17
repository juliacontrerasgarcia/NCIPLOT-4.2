! Julia Conteras-Garcia <julia.contreras.garcia@gmail.com>,
! Erin R. Johnson <ejohnson29@ucmerced.edu>,
! A. Otero-de-la-Roza <aoterodelaroza@ucmerced.edu>
! Weitao Yang <weitao.yang@duke.edu>,
! Roberto A. Boto <robalboto@gmail.com>,
! Chaoyou Quan  <quanchaoyu@gmail.com>,
! Ruben Laplaza <rlaplaza@lct.jussieu.fr>
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
!
! This is NCIPLOT Ver. 4.0
! which was curated by R.L. for Github distribution
! after development by R.A.B, C.Q. and others.

program nciplot
   use param
   use tools_io
   use tools_math
   use reader
   use props

   implicit none

   integer, parameter :: mfiles = 100 ! max number of geometry files

   character*(mline) :: argv(2), oname
   integer :: argc, nfiles, ifile, idx, istat, ntotal
   integer :: i, j, k, nn0, nnf, lp, n1
   character*(mline) :: filein, line, oline, word
   logical :: ok, ispromol
   real*8 :: rdum, deltag
   ! the molecular info
   type(molecule), allocatable :: m(:)
   ! logical units
   integer :: lugc, ludc, luvmd, ludat
   ! cubes
   real*8, allocatable, dimension(:, :, :) :: crho, cgrad
   ! ligand, intermolecular keyword
   logical :: ligand, inter, intra
   real*8 :: rthres
   integer :: udat0
   ! radius and cube keywords
   logical :: autor
   real*8 :: x(3), xinit(3), xmax(3), xinc(3)
   integer :: nstep(3)
   ! noutput keyword
   integer :: noutput
   ! cutoffs
   real*8 :: rhocut, dimcut
   ! cutplot
   real*8 :: rhoplot, isordg
   ! discarding rho parameter
   real*8 :: rhoparam, rhoparam2
   ! properties of rho
   real*8 :: rho, grad(3), dimgrad, grad2, hess(3, 3)
   integer, parameter :: mfrag = 100 ! max number of fragments
   real*8 :: rhom(mfrag)
   ! eispack
   real*8 :: wk1(3), wk2(3), heigs(3), hvecs(3, 3)
   ! Fragments
   integer :: nfrag
   logical :: autofrag
   ! chk file
   real*8 :: xinc_init(3)
   ! Modification
   logical ::  dointeg
   ! Atcube
   integer :: natommax, igroup, natom
   integer, allocatable, dimension(:, :) :: group
   real*8, allocatable, dimension(:, :) :: xinitat, xmaxat, nstepat

   ! Initializing multilevel grids
   integer :: ind_g, ng
   integer :: indx(3), nstep_coarse(3)
   real*8, allocatable :: fginc(:)
   real*8 :: xinc_coarse(3)
   real*8, allocatable, dimension(:, :, :) :: tmp_crho, tmp_cgrad
   logical :: flag, firstgrid
   logical, allocatable :: rmbox_coarse(:, :, :), tmp_rmbox(:, :, :)
   integer :: cr, c0, c1, c2, c3, c4, c5, c6
   integer :: i0, j0, k0
   integer :: lumesh, lusol
   logical, allocatable :: vert_use(:, :, :)
   real*8 :: sum_rhon_vol(9), sum_signrhon_vol(7)
   real*8 :: sum_rhon_area(9)
   real*8, allocatable ::  rho_n(:)
   real*8, allocatable :: crho_n(:, :, :, :), cheig(:, :, :)
   ! crho_n is a variable created for multilevel grids. It is equivalen to rhom (density per fragment)
   real*8, allocatable :: tmp_crho_n(:, :, :, :), tmp_cheigs(:, :, :)
   real*8 :: percent

   ! Range integration variables
   logical :: dorange, flag_range(2, 2, 2)
   integer :: nranges, l, l1(2), j1(2), k1(2)
   real*8, allocatable :: srhorange(:, :)
   real*8 :: upperbound, lowerbound

   !===============================================================================!
   ! System clock to measure run time.
   !===============================================================================!
   call system_clock(count_rate=cr)
   call system_clock(count=c0)

   call param_init()

   !===============================================================================!
   ! Reading positional parameters/arguments given to code. Input file to read.
   !===============================================================================!
   call getargs(argc, argv)
   if (argc >= 1) then
      open (uin, file=argv(1), status='old')
      if (argc == 2) then
         open (uout, file=argv(2), status='unknown')
      endif
   endif

   !===============================================================================!
   ! Internal clock starts! Header drops.
   !===============================================================================!
   call header()
   call tictac(' # Start')

   !===============================================================================!
   ! Check number of structure/wavefunction files, read them, check them, group them
   !===============================================================================!
   read (uin, *) nfiles   ! number of files to read
   if (nfiles > mfiles) call error('nciplot', 'too many files, increase mfiles', faterr)
   allocate (m(nfiles), stat=istat)
   if (istat /= 0) call error('nciplot', 'could not allocate memory for molecules', faterr)
   do ifile = 1, nfiles ! read files
      read (uin, '(a)') filein ! filein is read
      filein = trim(adjustl(filein))
      inquire (file=filein, exist=ok) ! check for existence
      if (.not. ok) &
         call error('nciplot', 'requested file does not exist: '//trim(filein), faterr)
      m(ifile) = readfile(filein)  ! reading wave function file (wfn or wfx) or structure (xyz file)
      if (ifile == 1) oname = filein(1:index(filein, '.', .true.) - 1)
      do i = 1, m(ifile)%n       ! assigning atoms to fragments, every atom in ifile to fragment ifile
         m(ifile)%ifrag(i) = ifile
      end do
   enddo
   nfrag = nfiles ! by default each file defines a fragment
   autofrag = .true.
   if (nfrag > mfrag) then ! too many fragments?
      call error('nciplot', 'too many fragments. Increase mfrag', faterr)
   end if

   !===============================================================================!
   ! Stop if user tries to mix computed and promolecular densities.
   !===============================================================================!
   if (any(m(:)%ifile == ifile_xyz) .and. any(m(:)%ifile == ifile_wfn)) then
      call error('nciplot', 'mixing xyz and wfn  not allowed', faterr)
   end if
   if (any(m(:)%ifile == ifile_xyz) .and. any(m(:)%ifile == ifile_wfx)) then
      call error('nciplot', 'mixing xyz and wfx not allowed', faterr)
   end if
   if (any(m(:)%ifile == ifile_wfn) .and. any(m(:)%ifile == ifile_wfx)) then
      call error('nciplot', 'mixing wfn and wfx is not advised', warning)
   end if

   !===============================================================================!
   ! Implicitely checking if the run mode is promolecular or not.
   !===============================================================================!
   ispromol = .not. (all(m(:)%ifile == ifile_wfn) .or. (all(m(:)%ifile == ifile_wfx)))

   ! read density grids (props.f90)
   call init_rhogrid(m, nfiles)

   ! by default, use density grids for heavier or charged atoms if promolecularity is on
   if (ispromol) then
      do i = 1, nfiles
         if (m(i)%ifile == ifile_xyz .and. (any(m(i)%z > atomic_zmax) .or. any(m(i)%q > 0))) then
            m(i)%ifile = ifile_grd
            do j = 1, m(i)%n
               if (m(i)%z(j) > atomic_zmax .and. .not. grd(iztype(m(i)%z(j), m(i)%q(j)))%init) then
                  call error('nciplot', 'Some atomic density grids for heavy atoms are needed but not initialized', faterr)
               end if
            end do
         end if
      end do
   end if

   !===============================================================================!
   ! Input files read and processed. Set defaults for running NCIPLOT now.
   !===============================================================================!
   rhocut = 0.5d0 ! density cutoff
   xinc = 0.1d0/bohrtoa   ! grid step
   if (any(m(:)%ifile == ifile_wfn) .or. any(m(:)%ifile == ifile_wfx)) then
      dimcut = 1.0d0  ! RDG cutoff
      isordg = 0.5d0  ! RDG isosurface
      rhoplot = 0.05d0 ! Density isosurface
   else
      dimcut = 1.0d0  ! RDG cutoff
      isordg = 0.3d0  ! RDG isosurface
      rhoplot = 0.07d0 ! Density isosurface
   end if
   rhoparam = 0.95d0  ! cutoffs for inter or intramolecularity
   rhoparam2 = 0.75d0 !
   noutput = 3        ! number of outputs
   udat0 = 1
   autor = .true.     ! build the cube automatically
   ligand = .false.   ! ligand keyword
   inter = .false.    ! intermolecular keyword
   rthres = 1.5d0     ! box limits around the molecule
   dointeg = .false.  ! integrating properties or not
   dorange = .false.  ! do not integrate range
   firstgrid = .true. ! flag for the first adaptive grid run
   if (.not. allocated(fginc)) then ! default setting CG2FG 3 4 2 1
      ng = 4
      allocate (fginc(ng))
      fginc = (/8, 4, 2, 1/)
   end if

   !===============================================================================!
   ! Estimating box around the molecule using xinit and xmax for the main system.
   !===============================================================================!
   xinit = m(1)%x(:, 1)
   xmax = m(1)%x(:, 1)
   do i = 1, nfiles
      do j = 1, m(i)%n
         xinit = min(xinit, m(i)%x(:, j))
         xmax = max(xmax, m(i)%x(:, j))
      enddo
   enddo
   ntotal = 0
   do i = 1, nfiles
      ntotal = ntotal + m(i)%n    ! compute the total number of atoms
   enddo

   !===============================================================================!
   ! Read optional keywords.
   ! Accepted keywords in this version:
   ! - RTHRES
   ! - LIGAND
   ! - RADIUS
   ! - INTERMOLECULAR
   ! - ONAME
   ! - INCREMENTS
   ! - OUTPUT
   ! - CUBE
   ! - ATCUBE
   ! - FRAGMENT
   ! - CUTOFFS
   ! - CUTPLOT
   ! - ISORDG
   ! - RHOCUT2
   ! - DGRID
   ! - RANGE
   ! - CG2FG
   ! - INTEGRATE
   !===============================================================================!
   do while (.true.)
      read (uin, '(a)', end=11) line
      line = trim(adjustl(line))
      oline = line
      call upper(line)
      if (line(1:1) == "#") cycle ! skip comments
      if (len(trim(line)) < 1) cycle ! skip blank lines
      idx = index(line, ' ')
      word = line(1:idx - 1)
      line = line(idx:)
      oline = oline(idx:)
      select case (trim(word))

      case ("RTHRES")       ! extra box limits
         read (line, *) rthres
         rthres = rthres/bohrtoa ! change to bohr

      case ("LIGAND")
         ligand = .true.              ! ligand option
         inter = .true.               ! intermolecular option automatically on
         read (line, *) udat0, rthres ! system in (udat0) as center
         rthres = rthres/bohrtoa      ! change to bohr

      case ("INTERMOLECULAR")
         inter = .true.              ! intermolecular option

      case ("RADIUS")
         autor = .false.
         read (line, *) x, rdum      ! center of the box
         xinit = (x - rdum)/bohrtoa  ! box limits
         xmax = (x + rdum)/bohrtoa

      case ("ONAME")         ! output name
         read (oline, '(a)') oname
         oname = trim(adjustl(oname))
         oname = oname(1:index(oname, ' '))

      case ("OUTPUT")
         read (line, *) noutput          ! number of output files

      case ("CUBE")         ! defining cube limits from coordinates. Example:
         autor = .false.    !CUBE x0,y0,z0,x1,y1,z1 format
         read (line, *) xinit, xmax

      case ("ATCUBE")          ! defining cube limits from atoms. Example:
         autor = .false.       !ATCUBE
         xinit = 1d40          !ifile atom1,atom2,...,atomn
         xmax = -1d40          !END
         natommax = 1
         do i = 1, nfiles
            natommax = max(natommax, m(i)%n)
         enddo

         allocate (group(ntotal, natommax))
         allocate (xmaxat(ntotal, 3))
         allocate (xinitat(ntotal, 3))
         allocate (nstepat(ntotal, 3))

         igroup = 0
         do while (.true.)
            igroup = igroup + 1
            read (uin, '(a)') line
            line = trim(adjustl(line))
            call upper(line)
            if (line(1:1) == "#") cycle ! skip comments
            if (len(trim(line)) < 1) cycle ! skip blank lines
            lp = 1
            idx = index(line, ' ')
            word = line(1:idx - 1)
            if (trim(word) /= "END" .and. trim(word) /= "ENDATCUBE") then
               ok = isinteger(ifile, line, lp)
               if (.not. ok) call error('nciplot', 'bad atcube syntax', faterr)
               if (ifile < 1 .or. ifile > nfiles) call error('nciplot', 'atcube: wrong file number', faterr)
               ok = isinteger(n1, line, lp)
               natom = 0
               do while (ok)
                  if (n1 < 1 .or. n1 > m(ifile)%n) call error('nciplot', 'atcube: wrong atom number', faterr)
                  natom = natom + 1
                  group(igroup, natom) = n1
                  xinit = min(xinit, m(ifile)%x(:, n1))
                  xmax = max(xmax, m(ifile)%x(:, n1))
                  ok = isinteger(n1, line, lp)
               end do
            else
               exit
            end if

            xinitat(igroup, :) = xinit - rthres
            xmaxat(igroup, :) = xmax + rthres
            nstepat(igroup, :) = abs(ceiling((xmax - xinit)/xinc))
         end do

      case ("FRAGMENT")          !defining fragments. Example:
         if (autofrag) then      !FRAGMENT
            nfrag = 0            !ifile atom1, atom2,...,atomn
            do ifile = 1, nfiles !END
               do i = 1, m(ifile)%n
                  m(ifile)%ifrag(i) = 0
               end do
            end do
         end if
         autofrag = .false.
         inter = .true.

         nfrag = nfrag + 1
         do while (.true.)
            read (uin, '(a)') line
            line = trim(adjustl(line))
            call upper(line)
            if (line(1:1) == "#") cycle ! skip comments
            if (len(trim(line)) < 1) cycle ! skip blank lines
            lp = 1
            idx = index(line, ' ')
            word = line(1:idx - 1)
            if (trim(word) /= "END" .and. trim(word) /= "ENDFRAGMENT") then
               ok = isinteger(ifile, line, lp)
               if (.not. ok) call error('nciplot', 'bad fragment syntax', faterr)
               if (ifile < 1 .or. ifile > nfiles) call error('nciplot', 'fragment: wrong file number', faterr)
               ok = isinteger(n1, line, lp)
               do while (ok)
                  if (n1 < 1 .or. n1 > m(ifile)%n) call error('nciplot', 'fragment: wrong atom number', faterr)
                  m(ifile)%ifrag(n1) = nfrag
                  ok = isinteger(n1, line, lp)
               end do
            else
               exit
            end if
         end do

      case ("INCREMENTS")   ! grid increments
         read (line, *) xinc
         xinc = max(xinc, 1d-4)
         xinc = xinc/bohrtoa ! transforming to angstrom to bohr

      case ("CUTOFFS")          ! density and RDG cutoffs
         read (line, *) rhocut, dimcut

      case ("CUTPLOT")          ! density cutoff used in the VMD script
         read (line, *) rhoplot, isordg

      case ("ISORDG")           !RDG isosurface used in the RDG script
         read (line, *) isordg

      case ("RHOCUT2")          ! cutoffs for intermolecularity definition
         read (line, *) rhoparam, rhoparam2

      case ("DGRID")            ! using grids for promolecular densities
         do i = 1, nfiles
            if (m(i)%ifile == ifile_xyz) m(i)%ifile = ifile_grd
         end do

      case ("CG2FG") ! coarse grid to fine grid multi-level
         read (line, *) ng ! number of multi-level grids
         if (allocated(fginc)) then
            deallocate (fginc)
         end if
         allocate (fginc(ng)) ! factors of grid increments. Example:
         read (line(:), *) ng, fginc ! CG2FG 4 8 4 2 1

      case ("RANGE")  ! range integration
         read (line, *) nranges ! number of ranges
         if (nranges .le. 0) then
            call error('nciplot', 'No ranges were given', faterr)
         else
            dorange = .true.
            allocate (srhorange(nranges, 2), stat=istat)
            if (istat /= 0) call error('nciplot', 'could not allocate memory for range intervals', faterr)
            do i = 1, nranges
               read (uin, *) srhorange(i, :)
               do j = 1, 2
                  if (abs(srhorange(i, j)) .lt. 1d-15) then
                     srhorange(i, j) = srhorange(i, j) + 1d-15
                  endif
               enddo
            enddo
         endif

      case ("INTEGRATE")  ! integration
         dointeg = .true.              ! intermolecular option

      case default ! something else is read
         call error('nciplot', 'Don''t know what to do with '//trim(word)//' keyword', faterr)
      end select

   enddo
11 continue

   !===============================================================================!
   ! Defining box in detail now.
   !===============================================================================!
   if (autor) then ! automatically build grid, it has not been done
      if (ligand) then ! ligand mode is enabled
         xinit = m(udat0)%x(:, 1)
         xmax = m(udat0)%x(:, 1)
         do j = 1, m(udat0)%n
            xinit = min(xinit, m(udat0)%x(:, j))
            xmax = max(xmax, m(udat0)%x(:, j))
         end do
      end if
      xinit = xinit - rthres
      xmax = xmax + rthres
   end if
   nstep = abs(ceiling((xmax - xinit)/xinc)) !number of grid steps

   !===============================================================================!
   ! Information for user and output files. Default logical units first.
   !===============================================================================!
   lugc = -1 ! RDG logical unit
   ludc = -1 ! Density logical unit
   luvmd = -1 ! VMD logical unit
   write (uout, 131)
   if (inter) then
      write (uout, 132) ! tell user intermolecular mode is on
      if (nfrag .eq. 1) then ! not enough fragments
         call error('nciplot', 'not enough fragments for intermolecular', faterr)
      end if
   end if
   if (ligand) write (uout, 130) trim(m(udat0)%name)
   if (ispromol) write (uout, 133) ! tell user promolecular mode is on
   write (uout, 120)
   write (uout, 110) 'RHO  THRESHOLD   (au):', rhocut
   write (uout, 110) 'RDG  THRESHOLD   (au):', dimcut
   if (inter) write (uout, 110) 'DISCARDING RHO PARAM :', rhoparam
   if (ligand) write (uout, 110) 'RADIAL THRESHOLD  (A):', rthres
   write (uout, *)
!  write(uout,121) xinit, xmax, xinc, nstep ! this is currently not used because it will be adapted
   if (noutput >= 2) then    ! number of outputs --> 2: Only .CUBE files
      lugc = 9
      ludc = 10
      luvmd = 11
      open (lugc, file=trim(oname)//"-grad.cube")    ! RDG cube file
      open (ludc, file=trim(oname)//"-dens.cube")    ! Density cube file
      open (luvmd, file=trim(oname)//".vmd")         ! VMD script
   endif

   if (noutput == 1 .or. noutput == 3) then
      ludat = 16
      open (ludat, file=trim(oname)//".dat")         ! RDG vs sign(lambda2) file
   else
      ludat = -1
   endif

   !===============================================================================!
   ! Write output files.
   !===============================================================================!
   write (uout, 122)
   if (noutput == 1 .or. noutput == 3) then
      write (uout, 123) trim(oname)//".dat"
   end if

   if (noutput >= 2) then
      write (uout, 124) trim(oname)//"-grad.cube", &
         trim(oname)//"-dens.cube", &
         trim(oname)//".vmd"
   end if
   if (lugc > 0) call write_cube_header(lugc, 'grad_cube', '3d plot, reduced density gradient')
   if (ludc > 0) call write_cube_header(ludc, 'dens_cube', '3d plot, density')

   !===============================================================================!
   ! Start run, using multi-level grids.
   !===============================================================================!
   ind_g = 1  ! index of the multi-level grids
   xinc_init = xinc ! initial coarse grid
   allocate (rho_n(1:nfiles))

12 continue    ! set grids from coarse to fine
   xinc = fginc(ind_g)*xinc_init
   nstep = ceiling((xmax - xinit)/xinc)

   write (uout, *)    ! punch info
   write (uout, 121) ind_g, xinit, xmax, xinc, nstep

   !===============================================================================!
   ! Allocate memory for density and gradient.
   !===============================================================================!
   if (ind_g .eq. 1) then
      allocate (crho(0:nstep(1) - 1, 0:nstep(2) - 1, 0:nstep(3) - 1), stat=istat)
      if (istat /= 0) call error('nciplot', 'could not allocate memory for density cube', faterr)
      allocate (crho_n(0:nstep(1) - 1, 0:nstep(2) - 1, 0:nstep(3) - 1, 1:nfiles), stat=istat)
      if (istat /= 0) call error('nciplot', 'could not allocate memory for density cube', faterr)
      allocate (cheig(0:nstep(1) - 1, 0:nstep(2) - 1, 0:nstep(3) - 1), stat=istat)
      if (istat /= 0) call error('nciplot', 'could not allocate memory for density cube', faterr)
      allocate (cgrad(0:nstep(1) - 1, 0:nstep(2) - 1, 0:nstep(3) - 1), stat=istat)
      if (istat /= 0) call error('nciplot', 'could not allocate memory for grad', faterr)
   end if

   !===============================================================================!
   ! Allocate memory for coarse grid.
   !===============================================================================!
   if (ind_g .gt. 1) then
      if (allocated(tmp_crho)) deallocate (tmp_crho)
      allocate (tmp_crho(0:nstep(1) - 1, 0:nstep(2) - 1, 0:nstep(3) - 1), stat=istat)
      call move_alloc(tmp_crho, crho)
      if (allocated(tmp_crho_n)) deallocate (tmp_crho_n)
      allocate (tmp_crho_n(0:nstep(1) - 1, 0:nstep(2) - 1, 0:nstep(3) - 1, 1:nfiles), stat=istat)
      call move_alloc(tmp_crho_n, crho_n)
      if (allocated(tmp_cheigs)) deallocate (tmp_cheigs)
      allocate (tmp_cheigs(0:nstep(1) - 1, 0:nstep(2) - 1, 0:nstep(3) - 1), stat=istat)
      call move_alloc(tmp_cheigs, cheig)
      if (allocated(tmp_cgrad)) deallocate (tmp_cgrad)
      allocate (tmp_cgrad(0:nstep(1) - 1, 0:nstep(2) - 1, 0:nstep(3) - 1), stat=istat)
      call move_alloc(tmp_cgrad, cgrad)
   end if

   if (ispromol) then    ! promolecular densities
      call system_clock(count=c1)
      !$omp parallel do private (x,rho,grad,hess,heigs,hvecs,wk1,wk2,istat,grad2,&
      !$omp dimgrad,intra,rhom,flag,indx,i0,j0,k0) schedule(dynamic)
      do k = 0, nstep(3) - 1
         do j = 0, nstep(2) - 1
            do i = 0, nstep(1) - 1
               x = xinit + (/i, j, k/)*xinc
               if (.not. firstgrid) then
                  flag = .false.
                  do i0 = max(0, i - 1), i
                     do j0 = max(0, j - 1), j
                        do k0 = max(0, k - 1), k
                           ! For each x, look for i, j, k indexes in the previous coarser grid
                           indx = floor(((/i0, j0, k0/)*xinc)/xinc_coarse)
                           indx = (/min(nstep_coarse(1) - 2, indx(1)), min(nstep_coarse(2) - 2, indx(2)), &
                                    min(nstep_coarse(3) - 2, indx(3))/)

                           if ((.not. flag) .and. (.not. (rmbox_coarse(indx(1), indx(2), indx(3))))) then
                              flag = .true.
                              goto 20
                           end if
                        end do
                     end do
                  end do

                  if (.not. flag) then
                     crho(i, j, k) = 100d0
                     cgrad(i, j, k) = 100d0
                     cheig(i, j, k) = 0d0
                     cycle
                  end if
20                continue

               end if

               ! calculate properties at x
               rho_n = 0d0
               call calcprops_pro(x, m, nfiles, rho, rho_n, rhom(1:nfrag), nfrag, autofrag, &
                                  grad, hess, deltag)
               call rs(3, 3, hess, heigs, 0, hvecs, wk1, wk2, istat)
               rho = max(rho, 1d-30)
               grad2 = dot_product(grad, grad)
               dimgrad = sqrt(grad2)/(const*rho**(4.D0/3.D0))
               intra = inter .and. ((any(rhom(1:nfrag) >= sum(rhom(1:nfrag))*rhoparam)) .or. &
                                    (sum(rhom(1:nfrag)) < rhoparam2*rho))
               if (intra) dimgrad = -dimgrad
               !$omp critical (cubewrite)
               crho(i, j, k) = sign(rho, heigs(2))*100.D0
               cgrad(i, j, k) = dimgrad
               !crho_n variable
               do i0 = 1, nfiles
                  crho_n(i, j, k, i0) = rhom(i0)
               enddo
               !$omp end critical (cubewrite)

            end do !i=0,nstep(1)-1
         end do ! j=0,nstep(2)-1
      end do  ! k=0,nstep(3)-1
      !$omp end parallel do
      call system_clock(count=c2)
      write (*, "(A, F6.2, A)") ' Time for computing density & RDG = ', real(dble(c2 - c1)/dble(cr), kind=4), ' secs'

   else  ! wavefunction densities
      call system_clock(count=c1)
      call calcprops_wfn(xinit, xinc, nstep, m, nfiles, crho, cgrad, cheig)
      !$omp parallel do private (x,rho,grad,hess,heigs,hvecs,wk1,wk2,istat,grad2,&
      !$omp dimgrad,intra,rhom) schedule(dynamic)
      do k = 0, nstep(3) - 1
         do j = 0, nstep(2) - 1
            do i = 0, nstep(1) - 1
               x = xinit + (/i, j, k/)*xinc
               if (.not. firstgrid) then
                  ! check if x is used, not removed
                  flag = .false.
                  do i0 = max(0, i - 1), i
                     do j0 = max(0, j - 1), j
                        do k0 = max(0, k - 1), k
                           ! For each x, look for i, j, k indexes in the previous coarser grid
                           indx = floor(((/i0, j0, k0/)*xinc)/xinc_coarse)
                           indx = (/min(nstep_coarse(1) - 2, indx(1)), min(nstep_coarse(2) - 2, indx(2)), &
                                    min(nstep_coarse(3) - 2, indx(3))/)
                           if ((.not. flag) .and. (.not. (rmbox_coarse(indx(1), indx(2), indx(3))))) then
                              flag = .true.
                              goto 21
                           end if
                        end do
                     end do
                  end do

                  if (.not. flag) then
                     crho(i, j, k) = 100d0
                     cgrad(i, j, k) = 100d0
                     cheig(i, j, k) = 0d0
                     cycle
                  end if
21                continue

               end if
            end do
         end do
      end do
      ! Checking if interatomic
      if (inter) then
         !$omp parallel do private (x,rho,grad,hess,intra,rhom) schedule(dynamic)
         do k = 0, nstep(3) - 1
            do j = 0, nstep(2) - 1
               do i = 0, nstep(1) - 1
                  x = xinit + (/i, j, k/)*xinc
                  call calcprops_pro(x, m, nfiles, rho, rho_n, rhom(1:nfrag), nfrag, autofrag, &
                                     grad, hess, deltag)
                  intra = ((any(rhom(1:nfrag) >= sum(rhom(1:nfrag))*rhoparam)) .or. &
                           (sum(rhom(1:nfrag)) < rhoparam2*rho))
                  !$omp critical (cubewrite)
                  if (intra) cgrad(i, j, k) = -abs(cgrad(i, j, k))
                  do i0 = 1, nfiles
                     crho_n(i, j, k, i0) = rhom(i0)
                  enddo
                  !$omp end critical (cubewrite)
               enddo
            enddo
         enddo
         !$omp end parallel do
      endif !interatomic run
      call system_clock(count=c2)
      write (*, "(A, F6.2, A)") ' Time for computing density & RDG = ', real(dble(c2 - c1)/dble(cr), kind=8), ' secs'
   endif !ispromol
   if ((ind_g .le. ng) .or. (ng .eq. 1)) then
      xinc_coarse = xinc ! record increments of the previous coarse grid
      nstep_coarse = nstep
      firstgrid = .false.
      if (allocated(rmbox_coarse)) then
         deallocate (rmbox_coarse)
         allocate (tmp_rmbox(0:nstep(1) - 2, 0:nstep(2) - 2, 0:nstep(3) - 2), stat=istat)
         if (istat /= 0) call error('nciplot', 'could not allocate memory for tmp_rmbox', faterr)
         call build_rmbox_coarse(rhocut, dimcut, ng, ind_g, fginc, tmp_rmbox, crho, cgrad, nstep_coarse)
         call move_alloc(tmp_rmbox, rmbox_coarse)
      else
         allocate (rmbox_coarse(0:nstep(1) - 2, 0:nstep(2) - 2, 0:nstep(3) - 2), stat=istat)
         if (istat /= 0) call error('nciplot', 'could not allocate memory for rmbox_coarse', faterr)
         call build_rmbox_coarse(rhocut, dimcut, ng, ind_g, fginc, rmbox_coarse, crho, cgrad, nstep_coarse)
      end if
      if (allocated(tmp_rmbox)) then
         deallocate (tmp_rmbox)
      end if
   end if

! loop over multi-level grids
   ind_g = ind_g + 1
   if (ind_g .le. ng) then
      goto 12 ! shameful goto to end multilevel grids.
   end if

   !===============================================================================!
   ! Output of .mesh and .sol files.
   !===============================================================================!
   call system_clock(count=c3)
   if (noutput .eq. 4) then
      allocate (vert_use(0:nstep(1) - 1, 0:nstep(2) - 1, 0:nstep(3) - 1))
      if (noutput .eq. 4) then
         lumesh = 20
         lusol = 21
      else
         lumesh = -1
         lusol = -1
      endif
      if (noutput .eq. 3) then
         open (lumesh, file=trim(oname)//".mesh")
         open (lusol, file=trim(oname)//".sol")
         call write_mesh_file(lumesh, lusol, xinit, xinc, nstep, cgrad, xinc_coarse, rmbox_coarse, &
                              nstep_coarse, vert_use)
         close (lumesh)
         close (lusol)
      endif
   end if

   !===============================================================================!
   ! Write .dat file.
   !===============================================================================!
   do k = 0, nstep(3) - 1
      do j = 0, nstep(2) - 1
         do i = 0, nstep(1) - 1
            ! fragments for the wfn case
            intra = (cgrad(i, j, k) < 0)
            cgrad(i, j, k) = abs(cgrad(i, j, k))
            dimgrad = cgrad(i, j, k)
            rho = crho(i, j, k)/100d0
            ! write the dat file
            if (ludat > 0 .and. .not. intra .and. (abs(rho) < rhocut) .and. (dimgrad < dimcut) .and. &
                abs(rho) > 1d-30) then
               write (ludat, '(1p,E18.10,E18.10)') rho, dimgrad
            end if ! rhocut/dimcut

            ! write the cube files
            if ((abs(rho) > rhoplot) .or. (dimgrad > dimcut) .or. intra) then
               cgrad(i, j, k) = 100d0
            endif !rho cutoff
         end do
      end do
   end do

   !===============================================================================!
   ! Write cube files.
   !===============================================================================!
   if (ludc > 0) call write_cube_body(ludc, nstep, crho)          ! density
   if (lugc > 0) call write_cube_body(lugc, nstep, cgrad)         ! RDG

   call system_clock(count=c4)
   write (*, "(A, F6.2, A)") ' Time for writing outputs = ', real(dble(c4 - c3)/dble(cr), kind=8), ' secs'

   !===============================================================================!
   ! Integration for promolecular systems. Box removal.
   !===============================================================================!
   if (dointeg) then    ! dointeg
      if (ispromol) then
         !$omp parallel do private (x,rho,grad,hess,heigs,hvecs,wk1,wk2,istat,grad2,&
         !$omp dimgrad,intra,rhom) schedule(dynamic)
         do k = 0, nstep(3) - 2
            do j = 0, nstep(2) - 2
               do i = 0, nstep(1) - 2
                  x = xinit + (/i, j, k/)*xinc
                  call calcprops_pro(x, m, nfiles, rho, rho_n, rhom(1:nfrag), nfrag, autofrag, &
                                     grad, hess, deltag)
                  call rs(3, 3, hess, heigs, 0, hvecs, wk1, wk2, istat)
                  rho = max(rho, 1d-30)
                  grad2 = dot_product(grad, grad)
                  dimgrad = sqrt(grad2)/(const*rho**(4.D0/3.D0))
                  if (inter) then
                     if ((any(rhom(1:nfrag) >= sum(rhom(1:nfrag))*rhoparam)) .or. &
                         (sum(rhom(1:nfrag)) < rhoparam2*rho)) then
                        rmbox_coarse(i, j, k) = .true. !inactive
                     end if
                  end if
                  if (((dimgrad > dimcut) .and. .not. rmbox_coarse(i, j, k))) then
                     rmbox_coarse(i, j, k) = .true. !inactive if dimgrad > dimcut
                  endif
               end do !i = 0, nstep(3) - 2
            end do !j = 0, nstep(3) - 2
         end do !k = 0, nstep(3) - 2
         !$omp end parallel do
         percent = real(count(rmbox_coarse), kind=8)/(real(size(rmbox_coarse), kind=8))*100d0
         write (*, '(F6.2, A)') percent, '% of small boxes removed for promolecular integration'
      endif  ! ispromol
   endif !dointeg

   !===============================================================================!
   ! Integration for non-promolecular systems. Box removal.
   !===============================================================================!
   if (dointeg) then
      if (.not. ispromol) then
         do k = 0, nstep(3) - 2
            do j = 0, nstep(2) - 2
               do i = 0, nstep(1) - 2
                  x = xinit + (/i, j, k/)*xinc
                  rho = abs(crho(i, j, k))/100d0
                  dimgrad = abs(cgrad(i, j, k))
                  if (inter) then
                     if ((any(rhom(1:nfrag) >= sum(rhom(1:nfrag))*rhoparam)) .or. &
                         (sum(rhom(1:nfrag)) < rhoparam2*rho)) then
                        rmbox_coarse(i, j, k) = .true.
                     end if
                  end if
                  if (((dimgrad > dimcut) .and. .not. rmbox_coarse(i, j, k))) then
                     rmbox_coarse(i, j, k) = .true. !inactive
                  endif ! rhocut/dimcut
               enddo  !k = 0,nstep(3)-1
            enddo !j = 0,nstep(2)-1
         enddo  !i = 0,nstep(1)-1
         percent = real(count(rmbox_coarse), kind=8)/(real(size(rmbox_coarse), kind=8))*100d0
         write (*, '(F6.2, A)') percent, '% of small boxes removed for density integration'
      endif !not ispromol
   endif !dointeg

   !===============================================================================!
   ! Integration and printing, uses the defined rmbox_coarse.
   !===============================================================================!
   ! compute geometry date in the region enclosed by the RDG isosurface:
   ! integral of rho^n (n = 1, 1.5, 2, 2.5, 3, 4/3, 5/3, 0) and rho1*rho2, respectively over the volume and the surface: sum_rhon_vol, sum_rhon_area
   if (dointeg) then
      call dataGeom(sum_rhon_vol, sum_rhon_area, sum_signrhon_vol, xinc, nstep, crho, crho_n, rmbox_coarse, nfiles)
      write (uout, 117) sum_rhon_vol, sum_signrhon_vol, sum_rhon_area
      call system_clock(count=c5)
      write (*, "(A, F6.2, A)") ' Time for integration = ', real(dble(c5 - c4)/dble(cr), kind=8), ' secs'
   end if

   !===============================================================================!
   ! Range integration. Also uses the previously defined rmbox_coarse.
   !===============================================================================!
   if (dorange) then
      allocate (tmp_rmbox(0:nstep(1) - 2, 0:nstep(2) - 2, 0:nstep(3) - 2), stat=istat)
      if (istat /= 0) call error('nciplot', 'could not allocate memory for tmp_rmbox', faterr)
      write (uout, 134)
      call DoRangeWarning() ! warning: change in active box criterion
      do i = 1, nranges
         if (srhorange(i, 1) .lt. srhorange(i, 2)) then
            upperbound = srhorange(i, 2)
            lowerbound = srhorange(i, 1)
         else
            upperbound = srhorange(i, 1)
            lowerbound = srhorange(i, 2)
         endif
         if ((upperbound .gt. rhocut) .or. (abs(lowerbound) .gt. rhocut)) then
            call DoRangeWarning2() ! warning: range outside rhocut
         endif

         tmp_rmbox = rmbox_coarse
         do k = 0, nstep(3) - 2
            do j = 0, nstep(2) - 2
               do l = 0, nstep(1) - 2
                  if (.not. (tmp_rmbox(l, j, k))) then
                     l1 = (/l, l + 1/)
                     j1 = (/j, j + 1/)
                     k1 = (/k, k + 1/)
                     flag_range = (((crho(l1, j1, k1)/100d0) .lt. lowerbound) .or. ((crho(l1, j1, k1)/100d0) .gt. upperbound))

                     if (count(flag_range) .eq. 0) then
                        tmp_rmbox(l, j, k) = .false.

                     else
                        tmp_rmbox(l, j, k) = .true.

                     endif
                  endif
               enddo
            enddo
         enddo
         percent = real(count(tmp_rmbox), kind=8)/(real(size(rmbox_coarse), kind=8))*100d0
         write (*, '(F6.2, A,I2)') percent, '% of small boxes removed for integration of interval', i
         call dataGeom(sum_rhon_vol, sum_rhon_area, sum_signrhon_vol, xinc, nstep, crho, crho_n, tmp_rmbox, nfiles)
         write (uout, 135) lowerbound, upperbound, sum_rhon_vol, sum_signrhon_vol

      enddo
   endif
   call system_clock(count=c6)

   !===============================================================================!
   ! Deallocate grids.
   !===============================================================================!
   if (allocated(tmp_rmbox)) deallocate (tmp_rmbox)
   if (allocated(srhorange)) deallocate (srhorange)
   if (allocated(rmbox_coarse)) deallocate (rmbox_coarse)

   !===============================================================================!
   ! Deallocate grids, close files.
   !===============================================================================!
   if (allocated(crho)) deallocate (crho)
   if (allocated(cheig)) deallocate (cheig)
   if (allocated(cgrad)) deallocate (cgrad)
   if (allocated(crho_n)) deallocate (crho_n)
   if (ludat > 0) close (ludat)

   !===============================================================================!
   ! Write VMD script.
   !===============================================================================!
   if (ligand) then
      nn0 = sum(m(1:udat0 - 1)%n) + 1
      nnf = sum(m(1:udat0 - 1)%n) + m(udat0)%n
   else
      nn0 = 1
      nnf = ntotal
   end if

   if (luvmd > 0) then
      write (luvmd, 114) trim(oname)//"-dens.cube"
      write (luvmd, 115) trim(oname)//"-grad.cube"
      write (luvmd, 116) nn0 - 1, nnf - 1, isordg, 2, 2, 2, -rhoplot*100D0, rhoplot*100D0, 2, 2
      close (luvmd)
   end if

   !===============================================================================!
   ! Deallocate arrays, call clock end, close output and input files.
   !===============================================================================!
   if (allocated(group)) deallocate (group)
   if (allocated(xmaxat)) deallocate (xmaxat)
   if (allocated(xinitat)) deallocate (xinitat)
   if (allocated(nstepat)) deallocate (nstepat)
   if (allocated(m)) deallocate (m)
   if (allocated(tmp_crho)) deallocate (tmp_crho)
   if (allocated(tmp_crho_n)) deallocate (tmp_crho_n)
   if (allocated(tmp_cheigs)) deallocate (tmp_cheigs)
   if (allocated(tmp_cgrad)) deallocate (tmp_cgrad)
   if (allocated(fginc)) deallocate (fginc)

   call tictac('End')
   if (uin /= stdin) close (uin)
   if (uout /= stdout) close (uout)

   !===============================================================================!
   ! Formats used by NCIPLOT.
   !===============================================================================!
110 format(A, F5.2)

   ! VMD script
114 format('#!/usr/local/bin/vmd', /, &
          '# VMD script written by save_state $Revision: 1.10 $', /, &
          '# VMD version: 1.8.6            ', /, &
          'set viewplist            ', /, &
          'set fixedlist            ', /, &
          '# Display settings            ', /, &
          'display projection   Orthographic            ', /, &
          'display nearclip set 0.000000            ', /, &
          '# load new molecule         ', /, &
          'mol new ', a, ' type cube first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all')
115 format('mol addfile ', a, ' type cube first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all')
116 format('#', /, &
          '# representation of the atoms', /, &
          'mol delrep 0 top', /, &
          'mol representation Lines 1.00000', /, &
          'mol color Name', /, &
          'mol selection {all}', /, &
          'mol material Opaque', /, &
          'mol addrep top', /, &
          'mol representation CPK 1.000000 0.300000 118.000000 131.000000', /, &
          'mol color Name', /, &
          'mol selection {index ', i5, ' to ', i5, ' }', /, &
          'mol material Opaque', /, &
          'mol addrep top', /, &
          '#', /, &
          '# add representation of the surface', /, &
          'mol representation Isosurface ', f8.5, ' 1 0 0 1 1', /, &
          'mol color Volume 0', /, &
          'mol selection {all}', /, &
          'mol material Opaque', /, &
          'mol addrep top', /, &
          'mol selupdate ', i1, ' top 0', /, &
          'mol colupdate ', i1, ' top 0', /, &
          'mol scaleminmax top ', i1, ' ', f7.4, f7.4, /, &
          'mol smoothrep top ', i1, ' 0', /, &
          'mol drawframes top ', i1, ' {now}', /, &
          'color scale method BGR', /, &
          'set colorcmds { {color Name {C} gray} }', /, &
          '#some more',/)
117 format('                                                     ', / &
          '----------------------------------------------------------------------', /, &
          '                                                    ', / &
          '                     INTEGRATION DATA                        ', /, &
          '                                                     ', / &
          '----------------------------------------------------------------------', /, &
          ' Integration  over the volumes of rho^n                               '/, &
          '----------------------------------------------------------------------', /, &
          ' n=1.0           :', 3X, F15.8, /, &
          ' n=1.5           :', 3X, F15.8, /, &
          ' n=2.0           :', 3X, F15.8, /, &
          ' n=2.5           :', 3X, F15.8, /, &
          ' n=3.0           :', 3X, F15.8, /, &
          ' n=4/3           :', 3X, F15.8, /, &
          ' n=5/3           :', 3X, F15.8, /, &
          ' Volume          :', 3X, F15.8, /, &
          ' rho-sum_i rho_i :', 3X, F15.8, / &
          '                  ', / &
          '---------------------------------------------------------------------', /, &
          ' Integration  over the volumes of sign(lambda2)(rho)^n             '/, &
          '---------------------------------------------------------------------', /, &
          ' n=1.0           :', 3X, F15.8, /, &
          ' n=1.5           :', 3X, F15.8, /, &
          ' n=2.0           :', 3X, F15.8, /, &
          ' n=2.5           :', 3X, F15.8, /, &
          ' n=3.0           :', 3X, F15.8, /, &
          ' n=4/3           :', 3X, F15.8, /, &
          ' n=5/3           :', 3X, F15.8, /, &
          '----------------------------------------------------------------------', /, &
          ' Integration  over the areas of rho^n              '/, &
          '----------------------------------------------------------------------', /, &
          ' n=1.0           :', 3X, F15.8, /, &
          ' n=1.5           :', 3X, F15.8, /, &
          ' n=2.0           :', 3X, F15.8, /, &
          ' n=2.5           :', 3X, F15.8, /, &
          ' n=3.0           :', 3X, F15.8, /, &
          ' n=4/3           :', 3X, F15.8, /, &
          ' n=5/3           :', 3X, F15.8, /, &
          ' Area            :', 3X, F15.8, /, &
          ' rho-sum_i rho_i :', 3X, F15.8, /, &
          '----------------------------------------------------------------------', /, &
          '                   ',/)

120 format(/'-----------------------------------------------------'/ &
           '      Calculation details:'/ &
           '-----------------------------------------------------')
121 format(/, '-----------------------------------------------------'/ &
           '      Operating grid and increments: Grid-', I1/ &
           '-----------------------------------------------------'/ &
           ' x0,y0,z0  = ', f10.4, ' ', f10.4, ' ', f10.4/ &
           ' x1,y1,z1  = ', f10.4, ' ', f10.4, ' ', f10.4/ &
           ' ix,iy,iz  = ', f5.2, '   ', f5.2, '   ', f5.2/ &
           ' nx,ny,nz  = ', i4, '    ', i4, '    ', i4/)

122 format('-----------------------------------------------------'/ &
          '      Writing output in the following units:'/ &
          '-----------------------------------------------------'/)

   ! noutput=1 .or. noutput =3
123 format(' Sign(lambda2)xDensity x Reduced Density Gradient    = ', a,/)

   ! noutput >=2
124 format(' Reduced Density Gradient,RDG      = ', a, / &
          ' Sign(lambda2)xDensity,LS          = ', a, / &
          ' VMD script                        = ', a,/)

130 format('      Using ', a40, ' as LIGAND')
131 format('-----------------------------------------------------'/ &
          '      INPUT INFORMATION:'/ &
          '-----------------------------------------------------')
132 format(/'      MIND YOU'/ &
           '      ONLY ANALYZING INTERMOLECULAR INTERACTIONS     '/)

133 format(/'      MIND YOU'/ &
           '      RUNNING IN PROMOLECULAR MODE     '/)

134 format('                                                                      ', / &
          '----------------------------------------------------------------------', / &
          '                                                                      ', / &
          '               RANGE INTEGRATION DATA                                 ', /, &
          '----------------------------------------------------------------------', /, &
          '                                                                      ')

135 format('----------------------------------------------------------------------', /, &
          ' Interval        :', 2(3X, F15.8) '                     ', /, &
          '                               ', / &
          '----------------------------------------------------------------------', /, &
          ' Integration  over the volumes of rho^n                               '/, &
          '----------------------------------------------------------------------', /, &
          ' n=1.0           :', 3X, F15.8, /, &
          ' n=1.5           :', 3X, F15.8, /, &
          ' n=2.0           :', 3X, F15.8, /, &
          ' n=2.5           :', 3X, F15.8, /, &
          ' n=3.0           :', 3X, F15.8, /, &
          ' n=4/3           :', 3X, F15.8, /, &
          ' n=5/3           :', 3X, F15.8, /, &
          ' Volume          :', 3X, F15.8, /, &
          ' rho-sum_i rho_i :', 3X, F15.8, / &
          '                  ', / &
          '---------------------------------------------------------------------', /, &
          ' Integration  over the volumes of sign(lambda2)(rho)^n             '/, &
          '---------------------------------------------------------------------', /, &
          ' n=1.0           :', 3X, F15.8, /, &
          ' n=1.5           :', 3X, F15.8, /, &
          ' n=2.0           :', 3X, F15.8, /, &
          ' n=2.5           :', 3X, F15.8, /, &
          ' n=3.0           :', 3X, F15.8, /, &
          ' n=4/3           :', 3X, F15.8, /, &
          ' n=5/3           :', 3X, F15.8, /, &
          '---------------------------------------------------------------------', /, &
          '                         ',/)

contains

   !===============================================================================!
   ! Subroutines for .cube, .mesh and .sol writing.
   !===============================================================================!

   subroutine write_cube_header(lu, l1, l2)

      integer, intent(in) :: lu
      character*(*), intent(in) :: l1, l2

      integer :: i, j

      write (lu, *) trim(l1)
      write (lu, *) trim(l2)
      write (lu, '(I5,3(F12.6))') ntotal, xinit

      write (lu, '(I5,3(F12.6))') nstep(1), xinc(1), 0d0, 0d0
      write (lu, '(I5,3(F12.6))') nstep(2), 0d0, xinc(2), 0d0
      write (lu, '(I5,3(F12.6))') nstep(3), 0d0, 0d0, xinc(3)
      do i = 1, nfiles
         do j = 1, m(i)%n
            write (lu, '(I4,F5.1,F11.6,F11.6,F11.6)') m(i)%z(j), 0d0, m(i)%x(:, j)
         end do
      enddo

   end subroutine write_cube_header

   subroutine write_cube_body(lu, n, c)

      integer, intent(in) :: lu
      integer, intent(in) :: n(3)
      real*8, intent(in) :: c(0:n(1) - 1, 0:n(2) - 1, 0:n(3) - 1)

      integer :: i, j

      do i = 0, n(1) - 1
         do j = 0, n(2) - 1
            write (lu, '(6(1x,e12.5))') (c(i, j, k), k=0, n(3) - 1)
         enddo
      enddo
      close (lu)

   end subroutine write_cube_body

   ! write the .mesh file
   subroutine write_mesh_file(lumesh, lusol, xinit, xinc, nstep, cgrad, xinc_coarse, rmbox_coarse, nstep_coarse, vert_use)
      integer, intent(in) :: lumesh, lusol
      real*8, intent(in) :: xinit(3), xinc(3), xinc_coarse(3)
      integer, intent(in) :: nstep(3), nstep_coarse(3)
      real*8, intent(in) :: cgrad(0:nstep(1) - 1, 0:nstep(2) - 1, 0:nstep(3) - 1)
      logical, intent(in) :: rmbox_coarse(0:nstep_coarse(1) - 2, 0:nstep_coarse(2) - 2, 0:nstep_coarse(3) - 2)
      integer :: i, j, k, i0, j0, k0, m, count_vert, count_cube, c
      integer :: i1(2), j1(2), k1(2)
      integer :: tetra_cube(1:6, 1:4), ind_cube(1:8), indx(3)
      logical, intent(inout) :: vert_use(0:nstep(1) - 1, 0:nstep(2) - 1, 0:nstep(3) - 1)
      integer :: ind_vert(0:nstep(1) - 1, 0:nstep(2) - 1, 0:nstep(3) - 1)
      logical :: flag_use, ohexa, otetra
      ohexa = .false.
      otetra = .true.
      count_vert = 0
      count_cube = 0
      vert_use = .false.
      do i = 0, nstep(1) - 1
         do j = 0, nstep(2) - 1
            do k = 0, nstep(3) - 1
               ! eight neighbor boxes
               flag_use = .false.
               do i0 = max(0, i - 1), i
                  do j0 = max(0, j - 1), j
                     do k0 = max(0, k - 1), k
                        indx = floor(((/i0, j0, k0/)*xinc)/xinc_coarse)
                        indx = (/min(nstep_coarse(1) - 2, indx(1)), min(nstep_coarse(2) - 2, indx(2)), &
                                 min(nstep_coarse(3) - 2, indx(3))/)
                        if (.not. (rmbox_coarse(indx(1), indx(2), indx(3)))) then
                           vert_use(i, j, k) = .true.
                           flag_use = .true.
                           goto 30
                        end if
                     end do
                  end do
               end do
30             continue
               if (flag_use) then
                  count_vert = count_vert + 1
                  ind_vert(i, j, k) = count_vert
               end if
               indx = floor(((/i, j, k/)*xinc)/xinc_coarse)
               indx = (/min(nstep_coarse(1) - 2, indx(1)), min(nstep_coarse(2) - 2, indx(2)), &
                        min(nstep_coarse(3) - 2, indx(3))/)
               if (.not. (rmbox_coarse(indx(1), indx(2), indx(3)))) then
                  count_cube = count_cube + 1
               end if
            end do
         enddo
      enddo

      ! write data
      ! .mesh file
      write (lumesh, "(A)") 'MeshVersionFormatted 2'
      write (lumesh, "(A,/)") 'Dimension 3'
      write (lumesh, "(A)") 'Vertices'
      write (lumesh, '(I20)') count(vert_use)
      ! .sol file
      write (lusol, "(A)") 'MeshVersionFormatted 2'
      write (lusol, "(A,/)") 'Dimension 3'
      write (lusol, "(A)") 'SolAtVertices'
      write (lusol, '(I20)') count(vert_use)
      write (lusol, '(I16, I16)') 1, 1
      do i = 0, nstep(1) - 1
         do j = 0, nstep(2) - 1
            do k = 0, nstep(3) - 1
               if (vert_use(i, j, k)) then
                  write (lumesh, '(1x,3(e16.6), I8)') xinit + (/i, j, k/)*xinc, 1
                  write (lusol, '(1x, 1(f16.6))') cgrad(i, j, k)
               end if
            end do
         enddo
      enddo

      if (ohexa) then
         write (lumesh, "(/,A)") 'Hexaedra '
         write (lumesh, "(I20)") count_cube
         do i = 0, nstep(1) - 2
            do j = 0, nstep(2) - 2
               do k = 0, nstep(3) - 2
                  i1 = (/i, i + 1/)
                  j1 = (/j, j + 1/)
                  k1 = (/k, k + 1/)
                  c = count(vert_use(i1, j1, k1))
                  if (c .eq. 8) then
                     ind_cube = (/ind_vert(i, j, k), ind_vert(i, j, k + 1), ind_vert(i, j + 1, k + 1), ind_vert(i, j + 1, k), &
                        ind_vert(i + 1, j, k), ind_vert(i + 1, j, k + 1), ind_vert(i + 1, j + 1, k + 1), ind_vert(i + 1, j + 1, k)/)
                     write (lumesh, "(1x, 8(I12), I12)") ind_cube, 1
                  end if
               end do
            end do
         end do
      end if

      ! save tetras, by diving each cube into six tetras
      if (otetra) then
         write (lumesh, "(/,A)") 'Tetrahedra'
         write (lumesh, "(I20)") 6*count_cube
         do i = 0, nstep(1) - 2
            do j = 0, nstep(2) - 2
               do k = 0, nstep(3) - 2
                  i1 = (/i, i + 1/)
                  j1 = (/j, j + 1/)
                  k1 = (/k, k + 1/)
                  c = count(vert_use(i1, j1, k1))
                  if (c .eq. 8) then
                     ind_cube = (/ind_vert(i, j, k), ind_vert(i, j, k + 1), ind_vert(i, j + 1, k + 1), ind_vert(i, j + 1, k), &
                        ind_vert(i + 1, j, k), ind_vert(i + 1, j, k + 1), ind_vert(i + 1, j + 1, k + 1), ind_vert(i + 1, j + 1, k)/)
                     call cube2tetra(tetra_cube, ind_cube)
                     do m = 1, 6
                        write (lumesh, "(1x, 4(I16), I16)") tetra_cube(m, :), 1
                     end do
                  end if
               end do
            end do
         end do
      end if

      write (lumesh, '(A)') 'END'
      write (lusol, '(A)') 'END'
   end subroutine write_mesh_file

   !===============================================================================!
   ! Subroutines for multi-level grids.
   !===============================================================================!

   subroutine build_rmbox_coarse(rhocut, dimcut, ng, ind_g, fginc, rmbox_coarse, crho, cgrad, nstep)
      integer, intent(in) :: nstep(3), ng, ind_g
      real*8, intent(in) :: crho(0:nstep(1) - 1, 0:nstep(2) - 1, 0:nstep(3) - 1), &
                            cgrad(0:nstep(1) - 1, 0:nstep(2) - 1, 0:nstep(3) - 1)
      real*8, intent(in) :: rhocut, dimcut
      real*8, intent(in) :: fginc(ng)
      real*8 :: rhocut0, dimcut0
      integer :: i, j, k
      integer :: i1(2), j1(2), k1(2)
      real*8 :: percent
      logical, intent(out) :: rmbox_coarse(0:nstep(1) - 2, 0:nstep(2) - 2, 0:nstep(3) - 2)
      logical :: flag_grad(2, 2, 2)
      rmbox_coarse = .true.
      rhocut0 = rhocut*fginc(ind_g)
      dimcut0 = dimcut*fginc(ind_g)
      do i = 0, nstep(1) - 2
         do j = 0, nstep(2) - 2
            do k = 0, nstep(3) - 2
               i1 = (/i, i + 1/)
               j1 = (/j, j + 1/)
               k1 = (/k, k + 1/)
               if ((ind_g .eq. ng) .and. (count(crho(i1, j1, k1) .gt. 99) .gt. 0)) then
                  ! small box with removed vertices
                  cycle
               end if

               flag_grad = (((abs(crho(i1, j1, k1))/100) .gt. rhocut0) .or. (cgrad(i1, j1, k1) .gt. dimcut0) .or. &
                            (cgrad(i1, j1, k1) .lt. 0))

               if (count(flag_grad) .lt. 8) then
                  rmbox_coarse(i, j, k) = .false.  ! active box
               endif

            end do
         end do
      end do
      percent = real(count(rmbox_coarse), kind=8)/(real(size(rmbox_coarse), kind=8))*100d0
      write (*, '(F6.2, A)') percent, '% of small boxes are removed.'
   end subroutine build_rmbox_coarse

! write six faces of a cube
   subroutine cube2quad(quad_cube, ind_cube)
      integer, intent(in) :: ind_cube(1:8)
      integer, intent(inout) :: quad_cube(1:6, 1:4)
      integer :: i1, i2, i3, i4, i5, i6, i7, i8
      i1 = ind_cube(1)
      i2 = ind_cube(2)
      i3 = ind_cube(3)
      i4 = ind_cube(4)
      i5 = ind_cube(5)
      i6 = ind_cube(6)
      i7 = ind_cube(7)
      i8 = ind_cube(8)
      quad_cube = reshape((/i1, i2, i3, i4, &
                            i5, i6, i7, i8, &
                            i1, i2, i6, i5, &
                            i4, i3, i7, i8, &
                            i1, i5, i8, i4, &
                            i2, i6, i7, i3/), (/6, 4/), order=(/2, 1/))
   end subroutine cube2quad

   ! divide a cube into six triangles
   subroutine cube2tetra(tetra_cube, ind_cube)
      integer, intent(in) :: ind_cube(1:8)
      integer, intent(inout) :: tetra_cube(1:6, 1:4)
      integer :: i1, i2, i3, i4, i5, i6, i7, i8
      i1 = ind_cube(1)
      i2 = ind_cube(2)
      i3 = ind_cube(3)
      i4 = ind_cube(4)
      i5 = ind_cube(5)
      i6 = ind_cube(6)
      i7 = ind_cube(7)
      i8 = ind_cube(8)
      tetra_cube = reshape((/i1, i2, i3, i7, &
                             i1, i2, i6, i7, &
                             i1, i4, i3, i7, &
                             i1, i4, i8, i7, &
                             i1, i5, i6, i7, &
                             i1, i5, i8, i7/), (/6, 4/), order=(/2, 1/))
   end subroutine cube2tetra

! compute NCI geometry of the region enclosed by RDG isosurface
   subroutine dataGeom(sum_rhon_vol, sum_rhon_area, sum_signrhon_vol, xinc, nstep, crho, crho_n, rmbox_coarse, nfiles)
      real*8, intent(inout) :: sum_rhon_vol(9)
      real*8, intent(inout) :: sum_rhon_area(9)
      real*8, intent(inout) :: sum_signrhon_vol(7)
      real*8, intent(in) :: xinc(3)
      integer, intent(in) :: nstep(3)
      real*8, intent(in) :: crho(0:nstep(1) - 1, 0:nstep(2) - 1, 0:nstep(3) - 1), &
                            crho_n(0:nstep(1) - 1, 0:nstep(2) - 1, 0:nstep(3) - 1, 1:nfiles)
      logical, intent(in) :: rmbox_coarse(0:nstep(1) - 2, 0:nstep(2) - 2, 0:nstep(3) - 2)
      integer :: i, j, k, n
      integer, intent(in) :: nfiles
      integer :: i1(2), j1(2), k1(2), negative, positive
      real*8 :: sum_signrho, signlambda_2

      sum_rhon_vol = 0
      sum_signrhon_vol = 0
      negative = 0
      positive = 0
      ! integral of rho^n over the volume of cubes
      do i = 0, nstep(1) - 2
         do j = 0, nstep(2) - 2
            do k = 0, nstep(3) - 2
               if (.not. rmbox_coarse(i, j, k)) then
                  i1 = (/i, i + 1/)
                  j1 = (/j, j + 1/)
                  k1 = (/k, k + 1/)
                  ! n = 1, 1.5, 2, 2.5, 3, 4/3, 5/3: sum of rho^n
                  sum_rhon_vol(1) = sum_rhon_vol(1) + sum(abs(crho(i1, j1, k1)/100))*xinc(1)*xinc(2)*xinc(3)/8
                  sum_rhon_vol(2) = sum_rhon_vol(2) + sum(abs(crho(i1, j1, k1)/100)**1.5) &
                                    *xinc(1)*xinc(2)*xinc(3)/8
                  sum_rhon_vol(3) = sum_rhon_vol(3) + sum(abs(crho(i1, j1, k1)/100)**2) &
                                    *xinc(1)*xinc(2)*xinc(3)/8
                  sum_rhon_vol(4) = sum_rhon_vol(4) + sum(abs(crho(i1, j1, k1)/ &
                                                              100)**2.5)*xinc(1)*xinc(2)*xinc(3)/8
                  sum_rhon_vol(5) = sum_rhon_vol(5) + sum(abs(crho(i1, j1, k1)/100)**3) &
                                    *xinc(1)*xinc(2)*xinc(3)/8
                  sum_rhon_vol(6) = sum_rhon_vol(6) + sum(abs(crho(i1, j1, k1)/100)**(1.333)) &
                                    *xinc(1)*xinc(2)*xinc(3)/8
                  sum_rhon_vol(7) = sum_rhon_vol(7) + sum(abs(crho(i1, j1, k1)/100)**(1.666)) &
                                    *xinc(1)*xinc(2)*xinc(3)/8
                  ! n = 0: volume of cubes
                  sum_rhon_vol(8) = sum_rhon_vol(8) + xinc(1)*xinc(2)*xinc(3)
                  ! sum of rho1*rho2
                  do n = 1, nfiles
                     sum_rhon_vol(9) = sum_rhon_vol(9) + sum(crho_n(i1, j1, k1, n)/100* &
                                                      (abs(crho(i1, j1, k1)) - crho_n(i1, j1, k1, n))/100)*xinc(1)*xinc(2)*xinc(3)/8
                  end do

                  !sign_lambda2 x rho
                  sum_signrho = sum(crho(i1, j1, k1))
                  signlambda_2 = sign(1d0, sum_signrho)

                  sum_signrhon_vol(1) = sum_signrhon_vol(1) + sum(crho(i1, j1, k1)/100)*xinc(1)*xinc(2)*xinc(3)/8
                  sum_signrhon_vol(2) = sum_signrhon_vol(2) + signlambda_2*sum(abs(crho(i1, j1, k1)/100)**1.5) &
                                        *xinc(1)*xinc(2)*xinc(3)/8
                  sum_signrhon_vol(3) = sum_signrhon_vol(3) + signlambda_2*sum(abs(crho(i1, j1, k1)/100)**2) &
                                        *xinc(1)*xinc(2)*xinc(3)/8
                  sum_signrhon_vol(4) = sum_signrhon_vol(4) + signlambda_2*sum(abs(crho(i1, j1, k1)/ &
                                                                                   100)**2.5)*xinc(1)*xinc(2)*xinc(3)/8
                  sum_signrhon_vol(5) = sum_signrhon_vol(5) + signlambda_2*sum(abs(crho(i1, j1, k1)/100)**3) &
                                        *xinc(1)*xinc(2)*xinc(3)/8
                  sum_signrhon_vol(6) = sum_signrhon_vol(6) + signlambda_2*sum(abs(crho(i1, j1, k1)/100)**(1.333)) &
                                        *xinc(1)*xinc(2)*xinc(3)/8
                  sum_signrhon_vol(7) = sum_signrhon_vol(7) + signlambda_2*sum(abs(crho(i1, j1, k1)/100)**(1.666)) &
                                        *xinc(1)*xinc(2)*xinc(3)/8
               end if
            end do
         end do
      end do

      ! integral of rho^n over the surface of cubes
      sum_rhon_area = 0
      do i = 0, nstep(1) - 2
         do j = 0, nstep(2) - 2
            do k = 0, nstep(3) - 2
               if (.not. rmbox_coarse(i, j, k)) then
                  if (i .eq. 0) then
                     i1 = (/i, i/)
                     j1 = (/j, j + 1/)
                     k1 = (/k, k + 1/)
                     call compArea(sum_rhon_area, nstep, crho, crho_n, nfiles, i1, j1, k1, xinc(2), xinc(3))
                  elseif (rmbox_coarse(i - 1, j, k)) then
                     i1 = (/i, i/)
                     j1 = (/j, j + 1/)
                     k1 = (/k, k + 1/)
                     call compArea(sum_rhon_area, nstep, crho, crho_n, nfiles, i1, j1, k1, xinc(2), xinc(3))
                  elseif ((i .eq. nstep(1) - 2) .or. (rmbox_coarse(i, j, k))) then
                     i1 = (/i + 1, i + 1/)
                     j1 = (/j, j + 1/)
                     k1 = (/k, k + 1/)
                     call compArea(sum_rhon_area, nstep, crho, crho_n, nfiles, i1, j1, k1, xinc(2), xinc(3))
                  end if
                  if (j .eq. 0) then
                     i1 = (/i, i + 1/)
                     j1 = (/j, j/)
                     k1 = (/k, k + 1/)
                     call compArea(sum_rhon_area, nstep, crho, crho_n, nfiles, i1, j1, k1, xinc(1), xinc(3))
                  elseif (rmbox_coarse(i, j - 1, k)) then
                     i1 = (/i, i + 1/)
                     j1 = (/j, j/)
                     k1 = (/k, k + 1/)
                     call compArea(sum_rhon_area, nstep, crho, crho_n, nfiles, i1, j1, k1, xinc(1), xinc(3))
                  elseif ((j .eq. nstep(2) - 2) .or. (rmbox_coarse(i, j, k))) then
                     i1 = (/i, i + 1/)
                     j1 = (/j + 1, j + 1/)
                     k1 = (/k, k + 1/)
                     call compArea(sum_rhon_area, nstep, crho, crho_n, nfiles, i1, j1, k1, xinc(1), xinc(3))
                  end if
                  if (k .eq. 0) then
                     i1 = (/i, i + 1/)
                     j1 = (/j, j + 1/)
                     k1 = (/k, k/)
                     call compArea(sum_rhon_area, nstep, crho, crho_n, nfiles, i1, j1, k1, xinc(1), xinc(3))
                  elseif (rmbox_coarse(i, j, k - 1)) then
                     i1 = (/i, i + 1/)
                     j1 = (/j, j + 1/)
                     k1 = (/k, k/)
                     call compArea(sum_rhon_area, nstep, crho, crho_n, nfiles, i1, j1, k1, xinc(1), xinc(2))
                  elseif ((k .eq. nstep(3) - 2) .or. (rmbox_coarse(i, j, k))) then
                     i1 = (/i, i + 1/)
                     j1 = (/j, j + 1/)
                     k1 = (/k + 1, k + 1/)
                     call compArea(sum_rhon_area, nstep, crho, crho_n, nfiles, i1, j1, k1, xinc(1), xinc(2))
                  end if
               end if
            end do
         end do
      end do

   end subroutine dataGeom

   subroutine compArea(sum_rhon_area, nstep, crho, crho_n, nfiles, i1, j1, k1, a, b)
      ! compute integrals over the boundary surface of the union of cubes
      real*8, intent(inout) :: sum_rhon_area(9)
      real*8, intent(in) :: a, b
      integer, intent(in) :: nstep(3)
      real*8, intent(in) :: crho(0:nstep(1) - 1, 0:nstep(2) - 1, 0:nstep(3) - 1), &
                            crho_n(0:nstep(1) - 1, 0:nstep(2) - 1, 0:nstep(3) - 1, 1:nfiles)
      integer :: n
      integer, intent(in) :: nfiles
      integer, intent(in) :: i1(2), j1(2), k1(2)
      ! face sides (a, b)
      ! n = 1, 1.5, 2, 2.5, 3, 4/3, 5/3: sum of rho^n
      sum_rhon_area(1) = sum_rhon_area(1) + sum(abs(crho(i1, j1, k1)/100))*a*b/8
      sum_rhon_area(2) = sum_rhon_area(2) + sum(abs(crho(i1, j1, k1)/100)**1.5)*a*b/8
      sum_rhon_area(3) = sum_rhon_area(3) + sum(abs(crho(i1, j1, k1)/100)**2)*a*b/8
      sum_rhon_area(4) = sum_rhon_area(4) + sum(abs(crho(i1, j1, k1)/100)**2.5)*a*b/8
      sum_rhon_area(5) = sum_rhon_area(5) + sum(abs(crho(i1, j1, k1)/100)**3)*a*b/8
      sum_rhon_area(6) = sum_rhon_area(6) + sum(abs(crho(i1, j1, k1)/100)**(1.3333333)) &
                         *a*b/8
      sum_rhon_area(7) = sum_rhon_area(7) + sum(abs(crho(i1, j1, k1)/100)**(1.6666666)) &
                         *a*b/8
      ! area of iso-surface
      sum_rhon_area(8) = sum_rhon_area(8) + a*b
      ! sum of rho1*rho2
      do n = 1, nfiles
         sum_rhon_area(9) = sum_rhon_area(9) + sum(crho_n(i1, j1, k1, n)/100* &
                                                   (abs(crho(i1, j1, k1)) - crho_n(i1, j1, k1, n))/100)*a*b/8
      end do
   end subroutine compArea

end program
