! Julia Conteras-Garcia <julia.contreras.garcia@gmail.com>, 
! Erin R. Johnson <ejohnson29@ucmerced.edu>,
! A. Otero-de-la-Roza <aoterodelaroza@ucmerced.edu>
! Weitao Yang <weitao.yang@duke.edu>, 
! Roberto A. Boto <robalboto@gmail.com>, 
! Chaoyou Quan  <quanchaoyu@gmail.com>
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
!
! Integration algorithm modified. We made the comparation
! at the final curves
!
! This is NCIPLOT Ver. 4.0
! Simplified version without ELF,Laplacian, kinetic energy and energy densities
! and IGM.
! Computation of NCI volumns and charges with several references
! xyz, wfn available.
! C.Quan Multigrid implementation  
! Integration keyworkd suppresed.
! Line 1476 
! A cube is selected only if rho and redgrad values at its vertexes are within rho and redgrad ranges. 
! C.Quan original code rejected a cube only if rho and redgrad values at the vertexes are  outside 
! rho and redgrad ranges. 


program nciplot
  use param
  use tools_io
  use tools_math
  use reader
  use props

  implicit none

  integer, parameter :: mfiles = 100 !< only applies to rhom().

  character*(mline) :: argv(2), oname
  integer :: argc, nfiles, ifile, idx, istat, ntotal
  integer :: i, j, k, nn0, nnf, lp, n1
  character*(mline) :: filein, line, oline, word, wx, wc
  logical :: ok, ispromol
  real*8 :: rdum
  ! the molecular info
  type(molecule), allocatable :: m(:),mref(:)
  ! logical units
  integer :: lugc, ludc, luvmd, ludat, luelf, luxc, luchk
  logical :: lchk
  ! cubes
  real*8, allocatable, dimension(:,:,:) :: crho, cgrad, celf, cxc
  real*8, allocatable, dimension(:,:,:) :: ctp,cvirial,cedensity,clap
  ! ligand, intermolecular keywords
  logical :: ligand, inter, intra
  real*8 :: rthres
  integer :: udat0
  ! radius and cube keywords
  logical :: autor
  real*8 :: x(3), xinit(3), xmax(3), xinc(3),xrefinc(3)
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
  real*8 :: rho, grad(3), dimgrad, grad2, hess(3,3), elf, exc
  integer, parameter :: mfrag = 100
  real*8 :: rhom(mfrag)
  ! eispack
  real*8 :: wk1(3), wk2(3), heigs(3), hvecs(3,3)
  ! elf
  logical :: doelf
  logical :: dosplines
  ! xc
  integer :: ixc(2)
  ! fragments
  integer :: nfrag
  logical :: autofrag
  ! chk file
  logical :: alcrho, alcgrad, alcelf, alcxc, alctp, alcvirial, alclap, alcedensity
  real*8 :: xinit0(3), xinc0(3),xinc_init(3)
  integer :: nstep0(3)
  ! Modification
  logical ::  dointeg
  integer :: nrefs,splineref,nsteps,nelec,filemin,filemax,iref,ii
  real*8 :: dimgradref,step,volcube,rhocube,tp,rhoref
  real*8 :: gradm2(mfrag),gradm(mfrag)
  real*8, allocatable, dimension(:) :: xspl,yspl,ye2
  real*8, allocatable, dimension(:,:,:,:) :: crefgrad,crhow,cfilegrad,cfilerho

  real*8, allocatable, dimension(:,:,:,:) :: cheigs,refxinit,refnstep

  real*8 :: chargetot,chargepos,chargeneg,voltot,volpos,volneg,yy1,ypn
  real*8 :: dy(nspl)
  real*8 :: lap,virial,edensity
  logical :: dokinetic,dovirial,dointerpol
  real*8 :: egauss 

 ! Independent Gradient Model
  real*8 :: deltag
  real*8, allocatable, dimension(:,:,:) :: cdeltag
  logical :: doigm

 ! ATCUBE 
  integer :: natommax,igroup,natom
  integer, allocatable, dimension(:,:) :: group
  real*8 , allocatable, dimension(:,:) :: xinitat,xmaxat,nstepat


  integer :: nfile,nmol,ncpoints
  real*8  ::  integval,refval
  integer :: lulap, luvir, luctp, lueden
  type (reference), allocatable :: ref(:)
  real*8:: v1(3),v2(3),origen(3),isovalue,w(3),center(3) 
  ! initialize 
  ! multi-level grids, C Quan
  integer :: ind_g, ng
  integer :: indx(3), nstep_coarse(3)
  real*8, allocatable :: fginc(:)
  real*8 :: xinc_coarse(3)
  real*8, allocatable, dimension(:,:,:) :: crho_coarse, cgrad_coarse, tmp_crho, tmp_cgrad
  logical :: flag, firstgrid
  logical, allocatable :: rmbox_coarse(:,:,:), tmp_rmbox(:,:,:)
  integer :: cr, c0, c1, c2, c3, c4, c5, c6
  integer :: i0, j0, k0
  integer :: lumesh, lusol, ludata_geom
  logical, allocatable :: vert_use(:,:,:)
  real*8 :: sum_rhon_vol(9),sum_signrhon_vol(7)
  real*8 :: sum_rhon_area(9)
  real*8, allocatable ::  rho_n(:)
  real*8, allocatable :: crho_n(:,:,:,:), cheig(:,:,:) 
  !RAB: crho_n is a variable created by CQUAN. It is equivalen to rhom (Density per fragment)
  real*8, allocatable :: tmp_crho_n(:,:,:,:), tmp_cheigs(:,:,:) 
  real*8 :: percent 

 !Range integration:RAB 
  logical :: dorange,flag_range(2,2,2) 
  integer :: nrange,nranges,l,l1(2),j1(2),k1(2)
  real*8, allocatable :: srhorange(:,:)
  real*8 :: signlambda_2,upperbound,lowerbound

  firstgrid = .true.

  ! record run time
  call system_clock(count_rate=cr)
  call system_clock(count=c0)


  call param_init()
   
  ! I/O units, process arguments
  call getargs(argc,argv)
  if (argc >= 1) then
     open (uin,file=argv(1),status='old')
     if (argc == 2) then
        open(uout,file=argv(2),status='unknown')
     endif
  endif
  
  ! header and time
  call header()
  call tictac(' # Start')
 
  ! read files, define fragments
  read (uin,*) nfiles   ! number of files to read
  if (nfiles > mfiles) call error('nciplot','too many files, increase mfiles',faterr)
  allocate(m(nfiles),stat=istat)
  if (istat /= 0) call error('nciplot','could not allocate memory for molecules',faterr)

  ! reading files
  do ifile = 1, nfiles
     read (uin,'(a)') filein !wfn or xyz file
     filein = trim(adjustl(filein))
     inquire(file=filein,exist=ok)
     if(.not.ok) &
        call error('nciplot','requested file does not exist: '//trim(filein),faterr)
     m(ifile) = readfile(filein)  !reading wave function(wfn file), wfx file or structure (xyz file)   
     if (ifile == 1) oname = filein(1:index(filein,'.',.true.)-1) 
     ! assigning atoms to fragments
     ! every atom in ifile is assigned to the fragment ifile
     do i = 1, m(ifile)%n   
        m(ifile)%ifrag(i) = ifile 
     end do
  enddo

  nfrag = nfiles ! By default each file defines a fragment
  autofrag = .true.
  if (nfrag > mfrag) then ! Too many fragments?
     call error('nciplot','too many fragments. Increase mfrag',faterr)
  end if
  
  ! Stop if xyz and wfn are mixed
  if (any(m(:)%ifile == ifile_xyz) .and.any(m(:)%ifile == ifile_wfn)) then
     call error('nciplot','mixing xyz and wfn  not allowed',faterr) 
  end if 
  if (any(m(:)%ifile == ifile_xyz) .and. any(m(:)%ifile == ifile_wfx)) then     !ARIAS
     call error('nciplot','mixing xyz and wfx not allowed',faterr)              !ARIAS
  end if                                                                        !ARIAS
  if (any(m(:)%ifile == ifile_wfn) .and. any(m(:)%ifile == ifile_wfx)) then     !ARIAS
     call error('nciplot','mixing wfn and wfx not allowed',faterr)              !ARIAS
  end if                                                                        !ARIAS

  ispromol = .not.(all(m(:)%ifile == ifile_wfn).or.(all(m(:)%ifile == ifile_wfx))) !Is promolecular? 

  ! read density grids (props)  
  call init_rhogrid(m,nfiles) 

  ! by default, use density grids for heavier or charged atoms.
  do i = 1, nfiles 
     if (m(i)%ifile == ifile_xyz .and. (any(m(i)%z > atomic_zmax) .or. any(m(i)%q > 0))) then
        m(i)%ifile = ifile_grd
        do j = 1, m(i)%n
           if (m(i)%z(j) > atomic_zmax .and..not.grd(iztype(m(i)%z(j),m(i)%q(j)))%init) then 
              call error('nciplot','Some atomic density grids for heavy atoms are needed but not initialized',faterr)
           end if
        end do
     end if
  end do

  ! default values
  rhocut = 0.2d0 ! Density cutoff
  xinc = 0.1d0   ! Grid step 
  if (any(m(:)%ifile == ifile_wfn)) then
     dimcut = 2.0d0  ! RDG cutoff
     isordg = 0.5d0  ! RDG isosurface
     rhoplot= 0.05d0 ! Density isosurface
  else
     dimcut = 1.0d0  ! RDG cutoff
     isordg = 0.3d0  ! RDG isosurface 
     rhoplot= 0.07d0 ! Density isosurface
  end if 
  if (any(m(:)%ifile == ifile_wfx)) then        !ARIAS
    dimcut = 1.0d0                              !ARIAS, change 2.0d0 to 1.0d0
    isordg=0.5d0                                !ARIAS
    rhoplot=0.05d0                              !ARIAS
  else                                          !ARIAS
     dimcut = 1.0d0                             !ARIAS
     isordg=0.3d0                               !ARIAS
     rhoplot=0.07d0                             !ARIAS
  end if                                        !ARIAS
  rhoparam = 0.95d0  !
  rhoparam2 = 0.75d0 !
  noutput = 3        ! Number of outputs
  udat0 = 1          
  autor = .true. 
  ligand = .false.   ! Ligand keyword  
  inter = .false.    ! Intermolecular keyword
  rthres = 2.d0      ! Box limits around the molecule
  dosplines=.false.   ! Doing spline analysis?
  dointeg =.false.   ! Integrating properties? 

  ! Defining box around the molecule
  ! read the rest of (optional) keywords 
  ! xinit and xmax for the main system
  xinit = m(1)%x(:,1)
  xmax = m(1)%x(:,1) 
   do i = 1, nfiles 
      do j = 1, m(i)%n 
         xinit = min(xinit,m(i)%x(:,j))
         xmax = max(xmax,m(i)%x(:,j)) 
     enddo
   enddo 

  !Compute the total number of atoms   
  ntotal=0 
  do i=1,nfiles          
      ntotal=ntotal+m(i)%n   
  enddo


 !Reading optional keywords   
  do while (.true.)
     read (uin,'(a)',end=11) line
     line = trim(adjustl(line))
     oline = line
     call upper(line)
     if (line(1:1) == "#") cycle ! skip comments
     if (len(trim(line)) < 1) cycle ! skip blank lines
     idx = index(line,' ')
     word = line(1:idx-1)
     line = line(idx:)
     oline = oline(idx:)
     select case(trim(word)) 
     case("SPLINE")   !Doing splines and computing properties withing NCI regions
         dosplines=.true.  
         dointeg=.true.  
         ! allocating spline arrays
         allocate(xspl(nspl),stat=istat)
         if (istat /= 0) call error('nciplot','could not allocate memory for splines',faterr) 
         allocate(yspl(nspl),stat=istat)
         if (istat /= 0) call error('nciplot','could not allocate memory for splines',faterr) 
         allocate(ye2(nspl),stat=istat)
         if (istat /= 0) call error('nciplot','could not allocate memory for splines',faterr) 
         !Generating spline points of rho
         step=rhocut/nspl !spline steps
         xspl=step*(/(i,i=1,nspl)/) ! spline points
         !Initialiting spline variables
         yspl=5.0
         dy=0.d0
         ypn=0d0
 
     case("NREFS")  !References for computing properties
         read(line,*) nrefs !Number of reference files
         if (nrefs.le.0)  then 
             call error('nciplot','No references were given',faterr) 
         else     
             allocate(mref(nrefs),stat=istat) 
             if (istat /= 0) call error('nciplot','could not allocate memory for references molecules',faterr)
             allocate(ref(nrefs),stat=istat)
             if (istat /= 0) call error('nciplot','could not allocate memory for references data',faterr)
             do iref=1,nrefs                   
                read (uin,'(a)') filein
                filein = adjustl(filein)
                inquire(file=filein,exist=ok)
                if(.not.ok) &
                   call error('nciplot','reference requested file does not exist: '//trim(filein),faterr)
                   mref(iref) = readfile(filein)   
                   do i = 1, mref(iref)%n
                      mref(iref)%ifrag(i) = iref
                   end do !i
             enddo !iref  
             ! Defining box around the reference files
             do i = 1, nrefs
                    ref(i)%xinit = mref(i)%x(:,1) ! initial values
                    ref(i)%xmax  = mref(i)%x(:,1)  ! initial values 
                    do j = 1, mref(i)%n       
                      ref(i)%xinit= min(ref(i)%xinit,mref(i)%x(:,j))
                      ref(i)%xmax= max(ref(i)%xmax,mref(i)%x(:,j))
                     end do !m(i)%n              
                     ref(i)%xinit=ref(i)%xinit-rthres
                     ref(i)%xmax=ref(i)%xmax+rthres 
               ! number of steps
                     ref(i)%nstep=ceiling((ref(i)%xmax-ref(i)%xinit)/xinc) 
             enddo !nrefs
         endif  ! nrefs
   
     ! Extra box limits
     case ("RTHRES")
        read (line,*) rthres
        rthres = rthres / bohrtoa !Change to bohr
     
     case ("LIGAND")
        ligand = .true.              !Ligand option
        inter = .true.               !Intermolecular option
        read (line,*) udat0, rthres  !System in (udat0) as center
        rthres = rthres / bohrtoa    !Change to bohr

! 12 continue
     case ("INTERMOLECULAR")
        inter = .true.              !Intermolecular option 

     case ("RADIUS")                      
        autor=.false.
        read(line,*) x, rdum          !Center of the box 
        !Box limits
        xinit = (x - rdum) / bohrtoa  
        xmax = (x + rdum) / bohrtoa

     !Output name
     case ("ONAME")
        read(oline,'(a)') oname
        oname = trim(adjustl(oname))
        oname = oname(1:index(oname,' '))
     
     !Number of output files
     case ("OUTPUT")
        read(line,*) noutput

     !Defining cube limits from coordinates
     !CUBE x0,y0,z0,x1,y1,z1
     case ("CUBE")
        autor = .false.
        read(line,*) xinit, xmax 

     !Defining cube limits from atoms  
     !CUBE 
     !ifile atom1,atom2,...,atomn
     !END 
     case ("ATCUBE")
        autor = .false.
        xinit = 1d40
        xmax = -1d40
        natommax=1
        do i=1,nfiles
           natommax=max(natommax,m(i)%n) 
        enddo

        allocate(group(ntotal,natommax)) 
        allocate(xmaxat(ntotal,3)) 
        allocate(xinitat(ntotal,3)) 
        allocate(nstepat(ntotal,3)) 

        igroup=0 
        do while (.true.) 
           igroup=igroup+1
           read (uin,'(a)') line
           line = trim(adjustl(line))
           call upper(line)
           if (line(1:1) == "#") cycle ! skip comments
           if (len(trim(line)) < 1) cycle ! skip blank lines
           lp = 1
           idx = index(line,' ')
           word = line(1:idx-1)
           if (trim(word) /= "END" .and. trim(word) /= "ENDATCUBE") then
              ok = isinteger(ifile,line,lp)
              if (.not.ok) call error('nciplot','bad atcube syntax',faterr)
              if (ifile < 1.or. ifile > nfiles) call error('nciplot','atcube: wrong file number',faterr)
              ok = isinteger(n1,line,lp) 
              natom=0 
              do while (ok) 
                 if (n1<1.or.n1>m(ifile)%n) call error('nciplot','atcube: wrong atom number',faterr)          
                 natom=natom+1
                 group(igroup,natom) = n1
                 xinit = min(xinit,m(ifile)%x(:,n1))
                 xmax = max(xmax,m(ifile)%x(:,n1))
                 ok = isinteger(n1,line,lp)
              end do
           else
              exit
           end if 
          
          xinitat(igroup,:) = xinit - rthres
          xmaxat(igroup,:) = xmax + rthres
          nstepat(igroup,:) = abs(ceiling((xmax - xinit) / xinc)) 
        end do

     !Defining fragments
     !FRAGMENT
     !ifile atom1, atom2,...,atomn
     !END 
     case ("FRAGMENT")
        if (autofrag) then
           nfrag = 0
           do ifile = 1, nfiles
              do i = 1, m(ifile)%n
                 m(ifile)%ifrag(i) = 0
              end do
           end do
        end if
        autofrag = .false.
        inter = .true.

        nfrag = nfrag + 1
        do while (.true.)
           read (uin,'(a)') line
           line = trim(adjustl(line))
           call upper(line)
           if (line(1:1) == "#") cycle ! skip comments
           if (len(trim(line)) < 1) cycle ! skip blank lines
           lp = 1
           idx = index(line,' ')
           word = line(1:idx-1)
           if (trim(word) /= "END" .and. trim(word) /= "ENDFRAGMENT") then
              ok = isinteger(ifile,line,lp)
              if (.not.ok) call error('nciplot','bad fragment syntax',faterr)
              if (ifile < 1.or. ifile > nfiles) call error('nciplot','fragment: wrong file number',faterr)
              ok = isinteger(n1,line,lp)
              do while (ok)
                 if (n1<1.or.n1>m(ifile)%n) call error('nciplot','fragment: wrong atom number',faterr)
                 m(ifile)%ifrag(n1) = nfrag
                 ok = isinteger(n1,line,lp)
              end do
           else
              exit
           end if
        end do 

     !Grid increments  
     case ("INCREMENTS")
         read(line,*) xinc
         xinc = xinc / bohrtoa !Transforming to angstrom to bohr

     !Density and RDG cutoffs
     case ("CUTOFFS")
        read(line,*) rhocut, dimcut

     !Density cutoff used in the VMD script
     case ("CUTPLOT")
        read(line,*) rhoplot,isordg

     !RDG isosurface used in the RDG script
     case ("ISORDG")
        read(line,*) isordg
   
     case ("RHOCUT2")
        read(line,*) rhoparam, rhoparam2 

     ! Using grids for promolecular densities
     case ("DGRID")
        do i = 1, nfiles
           if (m(i)%ifile == ifile_xyz) m(i)%ifile = ifile_grd
        end do 
     case ("CG2FG") ! coarse grid to fine grid, C Quan
        read (line,*) ng ! number of multi-level grids
        allocate(fginc(ng))
        read (line(:),*) ng, fginc ! factors of grid increments 
     case("RANGE")  !Range integration
         read(line,*) nranges !Number of ranges
         if (nranges.le.0)  then 
             call error('nciplot','No ranges were given',faterr) 
         else     
             dorange=.true.
             allocate(srhorange(nranges,2),stat=istat)
             if (istat /= 0) call error('nciplot','could not allocate memory for range intervals',faterr)
             do i=1,nranges                   
                read (uin,*) srhorange(i,:)
             enddo  
         endif  ! nrange

     case default
        call error('nciplot','Don''t know what to do with '//trim(word)//' keyword',faterr)
     end select 

  enddo
11 continue 

! set fginc if fginc is not allocated. Default setting CG2FG 1 1 
  if (.not. allocated(fginc)) then 
    ng=1
    allocate(fginc(ng))
    fginc = 1
  end if


  !Defining Box 
  ! set grid limits and npts
  if (autor) then
     if (ligand) then
        xinit = m(udat0)%x(:,1)
        xmax = m(udat0)%x(:,1)
        do j = 1, m(udat0)%n
           xinit = min(xinit,m(udat0)%x(:,j))
           xmax = max(xmax,m(udat0)%x(:,j))
        end do 
     end if  
     xinit = xinit - rthres
     xmax = xmax + rthres 
  end if

  nstep = abs(ceiling((xmax - xinit) / xinc)) !number of grid steps

  
  
  ! punch info
  write(uout,131)
  if (inter) write(uout,132) 
  if (ligand) write(uout,130) trim(m(udat0)%name) 
  if (dosplines) then
     write(uout,133) 
     do iref=1,nrefs   ! Running over the references files
        write(uout,*) trim(adjustl(mref(iref)%name)) 
     enddo  
  end if 

  write(uout,120)
  write(uout,110) 'RHO  THRESHOLD   (au):', rhocut
  write(uout,110) 'RDG  THRESHOLD   (au):', dimcut
  if (inter) write(uout,110) 'DISCARDING RHO PARAM :',rhoparam
  if (ligand) write(uout,110) 'RADIAL THRESHOLD  (A):',rthres
  write (uout,*)
!  write(uout,121) xinit, xmax, xinc, nstep

  ! open output files 

  !Default option
  lugc  = -1 ! RDG logical unit
  ludc  = -1 ! Density logical unit
  luvmd = -1 ! VMD logical unit

  !Number of outputs --> 2: Only .CUBE files 
  if (noutput >= 2) then
     lugc  = 9
     ludc  = 10 
     luvmd = 11 
     open(lugc,file=trim(oname)//"-grad.cube")    ! RDG cube file
     open(ludc,file=trim(oname)//"-dens.cube")    ! Density cube file
     open(luvmd,file=trim(oname)//".vmd")         ! VMD script
  endif

  if (noutput == 1 .or. noutput == 3) then
      ludat = 16
      open(ludat,file=trim(oname)//".dat")         ! RDG vs sign(lambda2) file
  else
        ludat = -1
  endif  
 


  ! Printing ouputs files 
  
  write(uout,122)

  if (noutput==1 .or. noutput==3) then
     write(uout,123) trim(oname)//".dat"
  end if 

  if (noutput >=2) then
     write(uout,124) trim(oname)//"-grad.cube",&
                     trim(oname)//"-dens.cube",&
                     trim(oname)//".vmd"
                  
 end if
 
  ! write cube headers
  if (lugc > 0) call write_cube_header(lugc,'grad_cube','3d plot, reduced density gradient')
  if (ludc > 0) call write_cube_header(ludc,'dens_cube','3d plot, density')


!!!!!!!!!!!!    BEGIN  Multi-level grids, C Quan   !!!!!!!!!!!!!
  ind_g = 1  ! index of the multi-level grids 
  xinc_init= xinc ! initial coarse grid
  allocate( rho_n(1:nfiles))

12 continue
 ! set grids from coarse to fine
  xinc = fginc(ind_g)*xinc_init
  nstep = ceiling((xmax - xinit) / xinc)

 ! punch info
  write (uout,*)
  write(uout,121) ind_g, xinit, xmax, xinc, nstep
 
  ! allocate memory for density and gradient, C Quan 
  
  if (ind_g.eq.1) then
      allocate(crho(0:nstep(1)-1,0:nstep(2)-1,0:nstep(3)-1),stat=istat)
      if (istat /= 0) call error('nciplot','could not allocate memory for density cube',faterr) 
      allocate(crho_n(0:nstep(1)-1,0:nstep(2)-1,0:nstep(3)-1,1:nfiles),stat=istat)
      if (istat /= 0) call error('nciplot','could not allocate memory for density cube',faterr) 
      allocate(cheig(0:nstep(1)-1,0:nstep(2)-1,0:nstep(3)-1),stat=istat)
      if (istat /= 0) call error('nciplot','could not allocate memory for density cube',faterr)
      allocate(cgrad(0:nstep(1)-1,0:nstep(2)-1,0:nstep(3)-1),stat=istat)
      if (istat /= 0) call error('nciplot','could not allocate memory for grad',faterr)  
     
      if (doelf) then  
         allocate(celf(0:nstep(1)-1,0:nstep(2)-1,0:nstep(3)-1),stat=istat)
         if (istat /= 0) call error('nciplot','could not allocate memory for elf',faterr)  
      endif 
  end if


  ! Allocation for a coarser grid

  if (ind_g .gt. 1) then
      if (allocated(tmp_crho)) deallocate(tmp_crho)
      allocate(tmp_crho(0:nstep(1)-1,0:nstep(2)-1,0:nstep(3)-1),stat=istat)
      call move_alloc( tmp_crho, crho )
      if (allocated(tmp_crho_n)) deallocate(tmp_crho_n)
      allocate(tmp_crho_n(0:nstep(1)-1,0:nstep(2)-1,0:nstep(3)-1, 1:nfiles),stat=istat)
      call move_alloc( tmp_crho_n, crho_n )
      if (allocated(tmp_cheigs)) deallocate(tmp_cheigs)
      allocate(tmp_cheigs(0:nstep(1)-1,0:nstep(2)-1,0:nstep(3)-1),stat=istat)
      call move_alloc( tmp_cheigs, cheig)
      if (allocated(tmp_cgrad)) deallocate(tmp_cgrad)
      allocate(tmp_cgrad(0:nstep(1)-1,0:nstep(2)-1,0:nstep(3)-1),stat=istat)
      call move_alloc( tmp_cgrad, cgrad ) 
  end if 
 
    
  !Starting properties calculation 

  !if promolecular
  !   if compute splines
  !      if compute properties
  ! else --> wfn   
  !   if compute splines
  !     if compute properties
  !   if interatomic    
  !end   


  if (ispromol) then    ! Promolecular densities 
     call system_clock(count=c1) 
        !$omp parallel do private (x,rho,grad,hess,heigs,hvecs,wk1,wk2,istat,grad2,&
        !$omp dimgrad,intra,rhom,elf,exc,flag,indx,i0,j0,k0) schedule(dynamic)
        do k = 0, nstep(3)-1
           do j = 0, nstep(2)-1
              do i = 0, nstep(1)-1
                 x = xinit + (/i,j,k/) * xinc 
                 if (.not. firstgrid) then
                    ! check if x is used, not removed
                    flag = .false.
                    do i0 = max(0,i-1),i
                        do j0 = max(0,j-1),j
                           do k0 = max(0,k-1),k
                               ! For each x, look for i, j, k indexes in the previous coarser grid
                               indx = floor( ((/i0,j0,k0/)*xinc)/xinc_coarse ) 
                               indx = (/ min(nstep_coarse(1)-2, indx(1)), min(nstep_coarse(2)-2, indx(2)), &
                                    min(nstep_coarse(3)-2, indx(3))/) 
                              
                               if ( (.not. flag) .and. (.not.(rmbox_coarse(indx(1), indx(2), indx(3)))) ) then
                                    flag = .true.
                                    goto 20
                               end if                          
                            end do
                        end do
                    end do 
                  
                    if (.not. flag) then
                        crho(i,j,k) = 100d0
                        cgrad(i,j,k) = 100d0
                        cheig(i,j,k) = 0d0
                        cycle
                    end if
20                  continue 
                
                 end if
  
                 ! calculate properties at x  
                 rho_n= 0d0  
                 call calcprops_pro(x,m,nfiles,rho,rho_n,rhom(1:nfrag),nfrag,autofrag,&
                    grad,hess,doelf,elf,ixc,exc,deltag)
                 call rs(3,3,hess,heigs,0,hvecs,wk1,wk2,istat) 
                 rho = max(rho,1d-30) 
                 grad2 = dot_product(grad,grad)
                 dimgrad = sqrt(grad2) / (const*rho**(4.D0/3.D0))           
                 intra = inter .and. ((any(rhom(1:nfrag) >= sum(rhom(1:nfrag))*rhoparam)) .or. &
                    (sum(rhom(1:nfrag)) < rhoparam2 * rho))
                 if (intra) dimgrad = -dimgrad 
                 !$omp critical (cubewrite)
                 crho(i,j,k) = sign(rho,heigs(2))*100.D0
                 cgrad(i,j,k) = dimgrad  
                 !crho_n variable 
                 do i0=1,nfiles 
                    crho_n(i,j,k,i0)=rhom(i0) 
                 enddo
                  !$omp end critical (cubewrite) 

               
              end do !i=0,nstep(1)-1 
           end do ! j=0,nstep(2)-1
        end do  ! k=0,nstep(3)-1
       !$omp end parallel do
        call system_clock(count=c2)
   write(*,"(A, F6.2, A)") ' Time for computing density & RDG = ', real(dble(c2-c1)/dble(cr),kind=4), ' secs'
else   
   
       call calcprops_wfn(xinit,xinc,nstep,m,nfiles,crho,cgrad,doelf,celf,ixc,cxc,ctp,cheig) 


       call system_clock(count=c1)
        !$omp parallel do private (x,rho,grad,hess,heigs,hvecs,wk1,wk2,istat,grad2,&
        !$omp dimgrad,intra,rhom,elf,exc) schedule(dynamic) 
        do k = 0, nstep(3)-1
           do j = 0, nstep(2)-1
              do i = 0, nstep(1)-1
                 x = xinit + (/i,j,k/) * xinc  
                 if (.not. firstgrid) then
                    ! check if x is used, not removed
                    flag = .false.
                    do i0 = max(0,i-1),i
                        do j0 = max(0,j-1),j
                            do k0 = max(0,k-1),k
                               ! For each x, look for i, j, k indexes in the previous coarser grid
                               indx = floor( ((/i0,j0,k0/)*xinc)/xinc_coarse )
                               indx = (/ min(nstep_coarse(1)-2, indx(1)), min(nstep_coarse(2)-2, indx(2)), &
                                    min(nstep_coarse(3)-2, indx(3))/)
                               if ( (.not. flag) .and. (.not.(rmbox_coarse(indx(1), indx(2), indx(3)))) ) then
                                    flag = .true.
                                    goto 21
                               end if
                            end do
                        end do
                    end do 
                  
                    if (.not. flag) then
                        crho(i,j,k) = 100d0
                        cgrad(i,j,k) = 100d0
                        cheig(i,j,k) = 0d0
                        cycle
                    end if
21               continue 
                
                 end if
              end do 
            end do
         end do 

        ! Checking if interatomic
        if (inter) then
           !$omp parallel do private (x,rho,grad,hess,intra,rhom,elf,exc) schedule(dynamic)
           do k = 0, nstep(3)-1
              do j = 0, nstep(2)-1
                 do i = 0, nstep(1)-1
                    x = xinit + (/i,j,k/) * xinc
                    call calcprops_pro(x,m,nfiles,rho,rho_n,rhom(1:nfrag),nfrag,autofrag,&
                       grad,hess,doelf,elf,ixc,exc,deltag)
                       intra = ((any(rhom(1:nfrag) >= sum(rhom(1:nfrag))*rhoparam)) .or. &
                       (sum(rhom(1:nfrag)) < rhoparam2 * rho))
                    !$omp critical (cubewrite)
                    if (intra) cgrad(i,j,k) = -abs(cgrad(i,j,k))  
                    do i0=1,nfiles 
                       crho_n(i,j,k,i0)=rhom(i0) 
                    enddo
                    !$omp end critical (cubewrite) 
                 enddo
              enddo
           enddo
         !$omp end parallel do
         endif !inter
        call system_clock(count=c2)
        write(*,"(A, F6.2, A)") ' Time for computing density & RDG = ', real(dble(c2-c1)/dble(cr),kind=8), ' secs'
endif !ispromol
  ! record coarse grid info, C Quan
   !if ( (ind_g .lt. ng) .or. (ng.eq.1) ) then
      xinc_coarse = xinc ! record increments of the previous coarse grid
      nstep_coarse = nstep
      firstgrid = .false.
      if (allocated(rmbox_coarse)) then
          allocate(tmp_rmbox(0:nstep(1)-2, 0:nstep(2)-2, 0:nstep(3)-2),stat=istat)
          if (istat /= 0) call error('nciplot','could not allocate memory for tmp_rmbox',faterr)
          call build_rmbox_coarse( rhocut, dimcut, ng, ind_g, fginc, tmp_rmbox, xinit, crho, cgrad, xinc, &
                nstep)
          call move_alloc(tmp_rmbox, rmbox_coarse)
      else
         allocate(rmbox_coarse(0:nstep(1)-2, 0:nstep(2)-2, 0:nstep(3)-2),stat=istat)
         if (istat /= 0) call error('nciplot','could not allocate memory for rmbox_coarse',faterr)
         call build_rmbox_coarse( rhocut, dimcut, ng, ind_g, fginc, rmbox_coarse, xinit, crho, cgrad, xinc, &
                nstep)
      end if
      if (allocated(tmp_rmbox)) then
        deallocate(tmp_rmbox)
      end if
   !end if

! loop over multi-level grids
  ind_g = ind_g+1
  if (ind_g .le. ng) then
    goto 12
  end if

!!!!!!!!!!!!    END  Multi-level grids, C Quan   !!!!!!!!!!!!! 

!!!!!!!!!!!!   STAR REMOVE BOXES ABOVE  REFERENCE LINE. Roberto A. Boto  !!!!!!!!!!!!!!!!!!!!
                 


  ! Starting integration with the finest grid. 
    !if promolecular
  !   if compute splines
  !      if compute properties
  ! else --> wfn   
  !   if compute splines
  !     if compute properties
  !   if interatomic    
  !end 

 ! Allocate memory for reference                           
  if (dosplines) then
     do i=1,nrefs
        allocate(ref(i)%cgrad(0:ref(i)%nstep(1)-1,0:ref(i)%nstep(2)-1,0:ref(i)%nstep(3)-1),stat=istat)
        if (istat /= 0) call error('nciplot','could not allocate memory for reference gradients',faterr)
        allocate(ref(i)%crho(0:ref(i)%nstep(1)-1,0:ref(i)%nstep(2)-1,0:ref(i)%nstep(3)-1),stat=istat)
        if (istat /= 0) call error('nciplot','could not allocate memory for reference densities',faterr)
        allocate(ref(i)%celf(0:ref(i)%nstep(1)-1,0:ref(i)%nstep(2)-1,0:ref(i)%nstep(3)-1),stat=istat)
        if (istat /= 0) call error('nciplot','could not allocate memory for reference elf',faterr)
        allocate(ref(i)%cxc(0:ref(i)%nstep(1)-1,0:ref(i)%nstep(2)-1,0:ref(i)%nstep(3)-1),stat=istat)
        if (istat /= 0) call error('nciplot','could not allocate memory for reference xc',faterr)
        allocate(ref(i)%ctp(0:ref(i)%nstep(1)-1,0:ref(i)%nstep(2)-1,0:ref(i)%nstep(3)-1),stat=istat)
        if (istat /= 0) call error('nciplot','could not allocate memory for reference kinetic energy ',faterr)
        allocate(ref(i)%cheigs(3,0:ref(i)%nstep(1)-1,0:ref(i)%nstep(2)-1,0:ref(i)%nstep(3)-1),stat=istat)
        if (istat /= 0) call error('nciplot','could not allocate memory for reference densities',faterr)
     enddo
  endif

  if (dointeg) then    ! dointeg

!  1) Compute RHO and RedGrad in reference systems 
!  2) Use that value to draw a spline line 
!-----------------------------------------------------------------------------------------------  
      if (ispromol) then ! Do spline for  computing properties 
         if (dosplines) then  
           !$omp parallel do private (iref,x,rho,grad,hess,heigs,hvecs,wk1,wk2,istat,grad2,&
           !$omp dimgrad) schedule(dynamic) 
           do iref=1,nrefs ! Running over the number of referencesa    
              do k = 0, ref(iref)%nstep(3)-1
                 do j = 0, ref(iref)%nstep(2)-1
                    do i = 0, ref(iref)%nstep(1)-1
                       x = ref(iref)%xinit + (/i,j,k/) * xinc
                      ! calculate properties at the point x for the references files               
                      call calcprops_proref(x,mref,nrefs,iref,rho,dimgrad,nfrag,autofrag,&
                      grad,hess,doelf,elf,ixc,exc)
                      call rs(3,3,hess,heigs,0,hvecs,wk1,wk2,istat)
                      rho = max(rho,1d-30)
                      grad2 = dot_product(grad,grad)
                      ref(iref)%dimgrad = sqrt(grad2) / (const*rho**(4.D0/3.D0))
                      ! compute splines                      
                      call spline(ref(iref)%rho,ref(iref)%dimgrad,xspl,yspl,ye2)
                    enddo !i=0,ref(iref)%nstep(1)-1 
                 enddo ! j=0,ref(iref)%nstep(2)-1
              enddo  ! k=0,ref(iref)%nstep(3)-1 
           enddo !nrefs
          !$omp end parallel do 

          call nrspline(xspl,yspl,nspl,dy(1),ypn,ye2)
        endif !splines  
!-----------------------------------------------------------------------------------------------   
! 3)  At each point, compute RHO and RedGrad for the main system. For every RHO, interpolate 
!     the reference values by cubic splines and get a reference value. RedGrad0. If RedGrad < 
!     RedGrad0, do not remove the box.             
!-----------------------------------------------------------------------------------------------   
        !$omp parallel do private (x,rho,grad,hess,heigs,hvecs,wk1,wk2,istat,grad2,&
        !$omp dimgrad,intra,rhom,elf,exc) schedule(dynamic) 
        do k = 0,nstep(3)-1
           do j = 0,nstep(2)-1
              do i = 0,nstep(1)-1
                 x = xinit + (/i,j,k/) * xinc
                 call calcprops_pro(x,m,nfiles,rho,rho_n,rhom(1:nfrag),nfrag,autofrag,&
                      grad,hess,doelf,elf,ixc,exc,deltag)
                 call rs(3,3,hess,heigs,0,hvecs,wk1,wk2,istat)
                 rho = max(rho,1d-30)
                 grad2 = dot_product(grad,grad)
                 dimgrad = sqrt(grad2) / (const*rho**(4.D0/3.D0))
                 intra = inter .and. ((any(rhom(1:nfrag) >= sum(rhom(1:nfrag))*rhoparam)) .or. &
                    (sum(rhom(1:nfrag)) < rhoparam2 * rho))
                 if (intra) dimgrad = -dimgrad  
                 do i0=1,nfiles
                    crho_n(i,j,k,i0)=rhom(i0) 
                 enddo
                 rhocube=rhocube+abs(rho)
                 if (dosplines) then
                     CALL NRSPLINT(XSPL,YSPL,YE2,NSPL,RHO,YY1)
                     refval=yy1
                 ! else
                 !    refval=integval
                 endif 

                 !  RAB: Intra, rhocut, dimcut filter are not needed.  
                 !  RAB: rmbox_coarse points array already accounts for them.
                 !  set to true. 
                 cgrad(i,j,k)=100d0 
                 crho(i,j,k)=100d0
                 if (.not.intra.and.(abs(rho)< rhocut).and.(dimgrad < dimcut).and. .not. rmbox_coarse(i,j,k)) then
                
                        if (dimgrad.lt.refval-0.01) then
                            rmbox_coarse(i,j,k)=.false.  
                            crho(i,j,k) = sign(rho,heigs(2))*100.D0 
                            cgrad(i,j,k) = dimgrad
                        else 
                            rmbox_coarse(i,k,k)=.true. 
                        !    crho(i,j,k)=100d0
                        !    cgrad(i,j,k)=100d0
                        endif   
       
                    endif
              end do
           end do
        end do

        !$omp end parallel do

!-----------------------------------------------------------------------------------------------   
! Boxes removed after integration  
        percent = real(count(rmbox_coarse),kind=8)/(real(size(rmbox_coarse),kind=8))*100d0
        write(*,'(F6.2, A)') percent, '% small boxes are removed after reference is applied.XYZ' 
     endif  ! ispromol 
  endif !dointeg


!-----------------------------------------------------------------------------------------------   
! Integration wfn files 

  if (dointeg) then
!  1) Compute RHO and RedGrad in reference systems 
!  2) Use that value to draw a spline line 
!-----------------------------------------------------------------------------------------------  
      if (.not. ispromol) then  
        if (dosplines) then
           do iref=1,nrefs   ! Running over the references files 
              nmol=1
              call calcprops_wfn(ref(iref)%xinit,xinc,ref(iref)%nstep,&
              mref(iref),nmol,ref(iref)%crho,ref(iref)%cgrad,doelf,&
              ref(iref)%celf,ixc,ref(iref)%cxc,ref(iref)%ctp,ref(iref)%cheigs)
              do k = 0, ref(iref)%nstep(3)-1
                 do j = 0, ref(iref)%nstep(2)-1
                   do i = 0, ref(iref)%nstep(1)-1
                         x = ref(i)%xinit + (/i,j,k/) * xinc
                         rhoref=abs(ref(iref)%crho(i,j,k))/100d0
                         dimgradref=abs(ref(iref)%cgrad(i,j,k))
                         call spline(rhoref,dimgradref,xspl,yspl,ye2)
                    enddo  !k = 0,nstep(3)-1
                  enddo !j = 0,nstep(2)-1
              enddo  !i = 0,nstep(1)-1  
            enddo !nrefs               
            call nrspline(xspl,yspl,nspl,dy(1),ypn,ye2)
        endif ! do splines
!-----------------------------------------------------------------------------------------------   

! 3)  At each point, compute RHO and RedGrad for the main system. For every RHO, interpolate 
!     the reference values by cubic splines and get a reference value. RedGrad0. If RedGrad < 
!     RedGrad0, do not remove the box.             
!-----------------------------------------------------------------------------------------------   
     ! Compute properties in the interacting system. Not needed. Already computed  
     !   nmol=1 
     !   call calcprops_wfn(xinit,xinc,nstep,m,nmol,crho,cgrad,doelf,celf,ixc, &
     !   cxc,ctp,cheigs)   
    
  !      if (dointeg) then
            do k = 0, nstep(3)-1
               do j = 0, nstep(2)-1
                  do i = 0, nstep(1)-1
                     x = xinit + (/i,j,k/) * xinc
                     rho=abs(crho(i,j,k))/100d0
                     if (dosplines) then
                        CALL NRSPLINT(XSPL,YSPL,YE2,NSPL,RHO,YY1)
                        refval=yy1
                    ! else 
                    !    refval=integval
                    endif
                    dimgrad = abs(cgrad(i,j,k))
                    !rhocube=rhocube+abs(rho) 
                 !  RAB: Intra, rhocut, dimcut filter are not needed.  
                 !  RAB: rmbox_coarse points array already accounts for them.
                 !  set to true. 
                    cgrad(i,j,k)=100d0
                    crho(i,j,k)=100d0 
                    if (((abs(rho) < rhocut) .and. (dimgrad < dimcut).and. .not.rmbox_coarse(i,j,k))) then 
                    !Computing properties
                         if (dimgrad.lt.refval-0.01) then
                              rmbox_coarse(i,j,k)=.false. 
                           else 
                              rmbox_coarse(i,j,k)=.true.
                    !          cgrad(i,j,k)=100d0
                    !          crho(i,j,k)=100d0
                         endif  !integration    
                    endif ! rhocut/dimcut 
                  enddo  !k = 0,nstep(3)-1
               enddo !j = 0,nstep(2)-1
            enddo  !i = 0,nstep(1)-1    

! Boxes removed after integration  
  percent = real(count(rmbox_coarse),kind=8)/(real(size(rmbox_coarse),kind=8))*100d0 
 ! write(*,'(F6.2, A)') percent, '% small boxes are removed after reference is applied. WFN'
 endif !wfn
endif !dointeg 

!!!!!!!!!!!!!!!!!! END REMOVE BOXES ABOVE REFERENCE LINE. Roberto A. Boto !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 allocate( vert_use(0:nstep(1)-1, 0:nstep(2)-1, 0:nstep(3)-1) )
 ! output .mesh and .sol file, C Quan
  call system_clock(count=c5)
  if (noutput .eq. 3) then
    lumesh = 20
    lusol = 21
  else
    lumesh = -1
    lusol = -1
  endif
  if (noutput .eq. 3) then
      open(lumesh,file=trim(oname)//".mesh")
      open(lusol,file=trim(oname)//".sol")
      call write_mesh_file(lumesh, lusol, xinit, xinc, nstep, crho, cgrad, xinc_coarse, rmbox_coarse, &
      nstep_coarse, vert_use)
      close(lumesh)
      close(lusol)
  endif
  call system_clock(count=c6) 
 
 ! Print integration data

 ! compute geometry date in the region enclosed by the RDG isosurface: C Quan
 ! integral of rho^n (n = 1, 1.5, 2, 2.5, 3, 4/3, 5/3, 0) and rho1*rho2, respectively over the volume and the surface: sum_rhon_vol, sum_rhon_area
    call dataGeom( sum_rhon_vol, sum_rhon_area, sum_signrhon_vol, xinc, nstep, crho, crho_n, cgrad, cheig, rmbox_coarse,nfiles)  
    
    write(uout,117) sum_rhon_vol,sum_signrhon_vol,sum_rhon_area 


 ! Temporal rmbox array for range integration 
  allocate(tmp_rmbox(0:nstep(1)-2, 0:nstep(2)-2, 0:nstep(3)-2),stat=istat)
  if (istat /= 0) call error('nciplot','could not allocate memory for tmp_rmbox',faterr) 

  if (dorange) then 
     write(uout,134)     
     call DoRangeWarning() !Warning: Change in active box criterion
     do i=1,nranges 
       if (srhorange(i,1).lt.srhorange(i,2)) then 
            upperbound=srhorange(i,2)
            lowerbound=srhorange(i,1) 
       else 
            upperbound=srhorange(i,1)
            lowerbound=srhorange(i,2) 
       endif 
       
       !Warning call: Range boundaries exceed rhocut. 

       if ((upperbound.gt.rhocut).or.(abs(lowerbound).gt.rhocut)) then 
           call DoRangeWarning2() 
       endif        


        tmp_rmbox=rmbox_coarse 
        do k = 0, nstep(3)-2
             do j = 0, nstep(2)-2
                 do l = 0, nstep(1)-2  
                    if (.not. (tmp_rmbox(l,j,k))) then      
                        l1 = (/l, l+1/)
                        j1 = (/j, j+1/)
                        k1 = (/k, k+1/)  
                        flag_range = ( ((crho(l1,j1,k1)/100d0) .lt. lowerbound) .or. ((crho(l1,j1,k1)/100d0) .gt. upperbound))  
    
                        if (count(flag_range) .eq. 0) then 
                            tmp_rmbox(l,j,k)=.false. 
                        
                        else 
                            tmp_rmbox(l,j,k)=.true. 
                                               
                        endif     
                    endif
                 enddo
             enddo
       enddo  
       percent = real(count(tmp_rmbox),kind=8)/(real(size(rmbox_coarse),kind=8))*100d0 
       write(*,'(F6.2, A,I2)') percent, '% small boxes are removed after reference is applied in interval', i
      call dataGeom(sum_rhon_vol, sum_rhon_area,sum_signrhon_vol ,xinc, nstep, crho, crho_n, cgrad, cheig, tmp_rmbox, nfiles)   

     write(uout,135) lowerbound,upperbound, sum_rhon_vol, sum_signrhon_vol

    enddo
  endif 
  
  if (allocated(tmp_rmbox)) deallocate(tmp_rmbox)
  if (allocated(srhorange)) deallocate(srhorange)

!  ludata_geom = 22
!  open(ludata_geom,position='Append',file="data_geom.txt")  ! or REWIND
!    write(ludata_geom, '(18(f20.10))') sum_rhon_vol, sum_rhon_area
!  close(ludata_geom)

!   write(uout, 127)
!   write(*,'(A,/,7(f16.6))') ' Sum of  rho^n (n = 1, 1.5, 2, 2.5, 3, 4/3, 5/3) =', sum_rhon_vol(1:7)
!   write(*,'(A,f16.6)') ' Volume of cubes (n = 0) = ', sum_rhon_vol(8)
!   write(*,'(A,f16.6)') ' Sum of rho1*rho2 =', sum_rhon_vol(9)



  if (allocated(rmbox_coarse)) deallocate(rmbox_coarse)


  ! apply cutoffs and writting .dat files 
  do k = 0, nstep(3)-1
     do j = 0, nstep(2)-1
        do i = 0, nstep(1)-1
           ! fragments for the wfn case
           intra = (cgrad(i,j,k) < 0)
           cgrad(i,j,k) = abs(cgrad(i,j,k)) 
           dimgrad = cgrad(i,j,k) 
           rho = crho(i,j,k) / 100d0 
           ! write the dat file
           if (ludat>0 .and. .not.intra .and. (abs(rho) < rhocut) .and. (dimgrad < dimcut) .and.& 
           abs(rho)>1d-30) then 
                 write(ludat,'(1p,E18.10,E18.10)') rho, dimgrad  
           end if ! rhocut/dimcut
           
           ! write the cube files
           if ((abs(rho) > rhoplot) .or. (dimgrad > dimcut) .or. intra)  then 
                  cgrad(i,j,k) = 100d0 
           ! Which values to discard ctp,clap,cvirial and cedensity ?
           endif !rho cutoff
        end do
     end do
  end do 


  ! write cubes
   if (ludc > 0) call write_cube_body(ludc,nstep,crho)          ! Density
   if (lugc > 0) call write_cube_body(lugc,nstep,cgrad)         ! RDG
  

  ! deallocate grids and close files
  if (allocated(crho)) deallocate(crho)
  if (allocated(cheig)) deallocate(cheig)
  if (allocated(cgrad)) deallocate(cgrad)

  if (ludat > 0) close(ludat) 


  ! write vmd script
  if (ligand) then
     nn0 = sum(m(1:udat0-1)%n) + 1
     nnf = sum(m(1:udat0-1)%n) + m(udat0)%n
  else
     nn0 = 1
     nnf = ntotal
  end if 
  
   if (luvmd > 0) then
     write (luvmd,114) trim(oname)//"-dens.cube"
     write (luvmd,115) trim(oname)//"-grad.cube"
     write (luvmd,116) nn0-1,nnf-1,isordg,2,2,2,-rhoplot*100D0,rhoplot*100D0,2,2
     close(luvmd) 
  end if 


  ! deallocate property arrays for references
    if (dosplines) then 
     if (allocated(mref)) deallocate(mref)
     do i=1,nrefs
       if (allocated(ref(i)%crho)) deallocate(ref(i)%crho)
       if (allocated(ref(i)%cgrad)) deallocate(ref(i)%cgrad)
       if (allocated(ref(i)%celf)) deallocate(ref(i)%celf)
       if (allocated(ref(i)%ctp)) deallocate(ref(i)%ctp)
       if (allocated(ref(i)%cheigs)) deallocate(ref(i)%cheigs)
       if (allocated(ref(i)%cxc)) deallocate(ref(i)%cxc)
     enddo 
    endif
    !if (allocated(refgrad)) deallocate(refgrad)
    !if (allocated(cheigs)) deallocate(cheigs)

    if (allocated(ye2)) deallocate(ye2) 
    
    ! Deallocate ATCUBE arrays
    if (allocated(group)) deallocate(group)    
    if (allocated(xmaxat)) deallocate(xmaxat)    
    if (allocated(xinitat)) deallocate(xinitat)    
    if (allocated(nstepat)) deallocate(nstepat)    
        
  ! end
  call tictac('End')

  ! close files
  if (uin /= stdin) close(uin)
  if (uout /= stdout) close(uout) 


  ! Formats

110 format (A,F5.2)

  ! VMD script
114 format ('#!/usr/local/bin/vmd',/,&
     '# VMD script written by save_state $Revision: 1.10 $',/,&
     '# VMD version: 1.8.6            ',/,&
     'set viewplist            ',/,&
     'set fixedlist            ',/,&
     '# Display settings            ',/,&
     'display projection   Orthographic            ',/,&
     'display nearclip set 0.000000            ',/,&
     '# load new molecule         ',/,&
     'mol new ',a,' type cube first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all')
115 format ('mol addfile ',a,' type cube first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all')
116 format ('#',/,&
     '# representation of the atoms',/,&
     'mol delrep 0 top',/,&
     'mol representation Lines 1.00000',/,&
     'mol color Name',/,&
     'mol selection {all}',/,&
     'mol material Opaque',/,&
     'mol addrep top',/,&
     'mol representation CPK 1.000000 0.300000 118.000000 131.000000',/,&
     'mol color Name',/,&
     'mol selection {index ', i5, ' to ', i5, ' }',/,&
     'mol material Opaque',/,&
     'mol addrep top',/,&
     '#',/,&
     '# add representation of the surface',/,&
     'mol representation Isosurface ',f8.5,' 1 0 0 1 1',/,&
     'mol color Volume 0',/,&
     'mol selection {all}',/,&
     'mol material Opaque',/,&
     'mol addrep top',/,&
     'mol selupdate ',i1,' top 0',/,&
     'mol colupdate ',i1,' top 0',/,&
     'mol scaleminmax top ',i1,' ', f7.4, f7.4,/,&
     'mol smoothrep top ', i1,' 0',/,&
     'mol drawframes top ',i1,' {now}',/,&
     'color scale method BGR',/,&
     'set colorcmds {{color Name {C} gray}}',/,&
     '#some more',/)
117  format('                                                     ',/&
            '----------------------------------------------------------------------',/,& 
            '                                                    ',/&
            '                     INTEGRATION DATA                        ',/,&
            '                                                     ',/& 
            '----------------------------------------------------------------------',/,& 
            ' Integration  over the volumes of rho^n                               '/,& 
            '----------------------------------------------------------------------',/,&  
            ' n=1.0           :',3X,F15.8,/,&
            ' n=1.5           :',3X,F15.8,/,&
            ' n=2.0           :',3X,F15.8,/,&
            ' n=2.5           :',3X,F15.8,/,&
            ' n=3.0           :',3X,F15.8,/,&
            ' n=4/3           :',3X,F15.8,/,&
            ' n=5/3           :',3X,F15.8,/,&  
            ' Volume          :',3X,F15.8,/,&  
            ' rho-sum_i rho_i :',3X,F15.8,/&   
            '                  ',/&      
            '---------------------------------------------------------------------',/,&  
            ' Integration  over the volumes of sign(lambda2)(rho)^n             '/,& 
            '---------------------------------------------------------------------',/,&  
            ' n=1.0           :',3X,F15.8,/,&
            ' n=1.5           :',3X,F15.8,/,&
            ' n=2.0           :',3X,F15.8,/,&
            ' n=2.5           :',3X,F15.8,/,&
            ' n=3.0           :',3X,F15.8,/,&
            ' n=4/3           :',3X,F15.8,/,&
            ' n=5/3           :',3X,F15.8,/,&  
            '----------------------------------------------------------------------',/,&  
            ' Integration  over the areas of rho^n              '/,& 
            '----------------------------------------------------------------------',/,&  
            ' n=1.0           :',3X,F15.8,/,&
            ' n=1.5           :',3X,F15.8,/,&
            ' n=2.0           :',3X,F15.8,/,&
            ' n=2.5           :',3X,F15.8,/,&
            ' n=3.0           :',3X,F15.8,/,&
            ' n=4/3           :',3X,F15.8,/,&
            ' n=5/3           :',3X,F15.8,/,&  
            ' Area            :',3X,F15.8,/,&  
            ' rho-sum_i rho_i :',3X,F15.8,/,& 
            '----------------------------------------------------------------------',/,&  
            '                   '         ,/) 

118  format(' Kinetic Energy Density :',3X,F15.8,/)
119  format(' Laplacian              :',3X,F15.8,/,&
            ' Local Virial Field     :',3X,F15.8,/,&
            ' Energy Density         :',3X,F15.8,/,&
             '-----------------------------------------------------',/)
120 format(/'-----------------------------------------------------'/&
            '      Calculation details:'/&
            '-----------------------------------------------------')
121 format(/,'-----------------------------------------------------'/&
             '      Operating grid and increments: Grid-', I1/&
             '-----------------------------------------------------'/&
             ' x0,y0,z0  = ',f10.4,' ',f10.4,' ',f10.4/&
             ' x1,y1,z1  = ',f10.4,' ',f10.4,' ',f10.4/&
             ' ix,iy,iz  = ',f5.2,'   ',f5.2,'   ',f5.2/&
             ' nx,ny,nz  = ',i4,'    ', i4,'    ', i4/) 


122 format('-----------------------------------------------------'/&
           '      Writing output in the following units:'/&
           '-----------------------------------------------------'/)

          ! noutput=1 .or. noutput =3 
123 format(' Sign(lambda2)xDensity x Reduced Density Gradient    = ',a,/)

          ! noutput >=2 
124 format(' Reduced Density Gradient,RDG      = ',a,/&
           ' Sign(lambda2)xDensity,LS          = ',a,/&
           ' VMD script                        = ',a,/)

          ! doelf
125 format(' ELF cube file                     = ',a,/)

          ! EXC
126 format(' XC energy density cube file       = ',a,/)

           !dokinetic
127 format( ' Kinetic energy density cube file  = ',a,/)
      
           ! dovirial
128 format(' Local Virial field cube file      = ',a,/&
           ' Laplacian cube file               = ',a,/&
           ' Energy density cube file          = ',a,/) 
           

130 format('      Using ',a40,' as LIGAND')
131 format('-----------------------------------------------------'/&
           '      INPUT INFORMATION:'/&
           '-----------------------------------------------------')
132 format(/'      MIND YOU'/&
            '      ONLY ANALYZING INTERMOLECULAR INTERACTIONS     '/)

133 format(/,'      Reading reference from : ',/)  



134  format('                                                                      ',/&
            '----------------------------------------------------------------------',/&
            '                                                                      ',/&
            '               RANGE INTEGRATION DATA                                 ',/,&
            '----------------------------------------------------------------------',/,& 
            '                                                                      ')

135  format('----------------------------------------------------------------------',/,& 
            ' Interval        :',2(3X,F15.8)'                     ',/,&
            '                               ',/&      
            '----------------------------------------------------------------------',/,& 
            ' Integration  over the volumes of rho^n                               '/,& 
            '----------------------------------------------------------------------',/,&  
            ' n=1.0           :',3X,F15.8,/,&
            ' n=1.5           :',3X,F15.8,/,&
            ' n=2.0           :',3X,F15.8,/,&
            ' n=2.5           :',3X,F15.8,/,&
            ' n=3.0           :',3X,F15.8,/,&
            ' n=4/3           :',3X,F15.8,/,&
            ' n=5/3           :',3X,F15.8,/,&  
            ' Volume          :',3X,F15.8,/,&  
            ' rho-sum_i rho_i :',3X,F15.8,/&   
            '                  '         ,/&       
            '---------------------------------------------------------------------',/,&  
            ' Integration  over the volumes of sign(lambda2)(rho)^n             '/,& 
            '---------------------------------------------------------------------',/,&  
            ' n=1.0           :',3X,F15.8,/,&
            ' n=1.5           :',3X,F15.8,/,&
            ' n=2.0           :',3X,F15.8,/,&
            ' n=2.5           :',3X,F15.8,/,&
            ' n=3.0           :',3X,F15.8,/,&
            ' n=4/3           :',3X,F15.8,/,&
            ' n=5/3           :',3X,F15.8,/,&  
            '---------------------------------------------------------------------',/,&  
            '                         ',/) 

contains 

  subroutine write_cube_header(lu,l1,l2)

    integer, intent(in) :: lu
    character*(*), intent(in) :: l1, l2

    integer :: i, j, k

    write(lu,*) trim(l1)
    write(lu,*) trim(l2) 
    write(lu,'(I5,3(F12.6))') ntotal, xinit 
 
    write(lu,'(I5,3(F12.6))') nstep(1), xinc(1), 0d0, 0d0
    write(lu,'(I5,3(F12.6))') nstep(2), 0d0, xinc(2), 0d0
    write(lu,'(I5,3(F12.6))') nstep(3), 0d0, 0d0, xinc(3) 
    do i=1,nfiles
      do j = 1, m(i)%n
         write(lu,'(I4,F5.1,F11.6,F11.6,F11.6)') m(i)%z(j), 0d0, m(i)%x(:,j)
      end do
    enddo 
   
  end subroutine write_cube_header

  subroutine write_cube_body(lu,n,c)
    
    integer, intent(in) :: lu
    integer, intent(in) :: n(3)
    real*8, intent(in) :: c(0:n(1)-1,0:n(2)-1,0:n(3)-1)

    integer :: i, j
   
    do i = 0, n(1)-1 
       do j = 0, n(2)-1 
          write (lu,'(6(1x,e12.5))') (c(i,j,k),k=0,n(3)-1)
       enddo
    enddo
    close(lu)

  end subroutine write_cube_body

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutines added for multi-grids, .mesh file, geometry data
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! build the rmbox_coarse which collects the removing information, C Quan
  subroutine build_rmbox_coarse( rhocut, dimcut, ng, ind_g, fginc, rmbox_coarse, xinit, crho, cgrad, &
  xinc, nstep)
  real*8, intent(in) ::  xinit(3), xinc(3)
  integer :: nstep(3), ng, ind_g
  real*8, intent(in) :: crho(0:nstep(1)-1,0:nstep(2)-1,0:nstep(3)-1), &
        cgrad(0:nstep(1)-1,0:nstep(2)-1,0:nstep(3)-1)
  real*8, intent(in) :: rhocut, dimcut
  real*8, intent(in) :: fginc(ng) 
  real*8 :: rhocut0, dimcut0
  integer :: i, j, k, cou
  integer :: i1(2), j1(2), k1(2)
  logical :: rm ! remove or not
  real*8 :: percent
  logical, intent(out) :: rmbox_coarse(0:nstep(1)-2, 0:nstep(2)-2, 0:nstep(3)-2)
  logical :: flag_grad(2,2,2) 
  rmbox_coarse = .true.
  rhocut0 = rhocut*fginc(ind_g)
  dimcut0 = dimcut*fginc(ind_g)
  do i = 0, nstep(1)-2
    do j = 0, nstep(2)-2
        do k = 0, nstep(3)-2
            i1 = (/i, i+1/)
            j1 = (/j, j+1/)
            k1 = (/k, k+1/) 
            if ( (ind_g .eq. ng) .and. (count(crho(i1,j1,k1) .gt. 99) .gt. 0) ) then
                ! small box with removed vertices
                cycle
            end if  

            flag_grad = (  ((abs(crho(i1,j1,k1))/100) .gt. rhocut0) .or. (cgrad(i1,j1,k1) .gt. dimcut0) .or. &
                                (cgrad(i1,j1,k1) .lt. 0))   

            if ( count( flag_grad ) .lt. 8 ) then
                 rmbox_coarse(i,j,k) = .false.  !Active box
            endif
       
        end do
    end do
  end do
  percent = real(count(rmbox_coarse),kind=8)/(real(size(rmbox_coarse),kind=8))*100d0
  write(*,'(F6.2, A)') percent, '% small boxes are removed'
  end subroutine build_rmbox_coarse

  ! write the .mesh file
  subroutine write_mesh_file(lumesh, lusol, xinit, xinc, nstep, crho, cgrad, xinc_coarse, rmbox_coarse, nstep_coarse, vert_use)
  integer, intent(in) :: lumesh, lusol
  real*8, intent(in) :: xinit(3), xinc(3), xinc_coarse(3)
  integer, intent(in) :: nstep(3), nstep_coarse(3)
  real*8, intent(in) :: crho(0:nstep(1)-1,0:nstep(2)-1,0:nstep(3)-1), &
        cgrad(0:nstep(1)-1,0:nstep(2)-1,0:nstep(3)-1)
  logical, intent(in) :: rmbox_coarse(0:nstep_coarse(1)-2, 0:nstep_coarse(2)-2, 0:nstep_coarse(3)-2)
  integer :: i, j, k, i0, j0, k0, m, n, count_vert, count_cube, c
  integer :: i1(2), j1(2), k1(2)
  integer :: tetra_cube(1:6,1:4), quad_cube(1:6,1:4), ind_cube(1:8), indx(3), vec(4)
  logical, intent(inout) :: vert_use(0:nstep(1)-1, 0:nstep(2)-1, 0:nstep(3)-1)
  integer :: ind_vert(0:nstep(1)-1, 0:nstep(2)-1, 0:nstep(3)-1)
  logical :: flag_use, ohexa, otetra
  ohexa = .false.
  otetra = .true.
  count_vert = 0
  count_cube = 0
  vert_use = .false.
do i = 0, nstep(1)-1
    do j = 0, nstep(2)-1
      do k = 0, nstep(3)-1
        ! eight neighbor boxes
        flag_use = .false.
        do i0 = max(0,i-1),i
            do j0 = max(0,j-1),j
                do k0 = max(0,k-1),k
                    indx = floor( ((/i0,j0,k0/)*xinc)/xinc_coarse )
                    indx = (/ min(nstep_coarse(1)-2, indx(1)), min(nstep_coarse(2)-2, indx(2)), &
                        min(nstep_coarse(3)-2, indx(3))/)
                    if ( .not.(rmbox_coarse(indx(1), indx(2), indx(3))) )  then
                        vert_use(i, j, k) = .true.
                        flag_use = .true.
                        goto 30
                    end if
                end do
            end do
        end do
30    continue
        if (flag_use) then
            count_vert = count_vert+1
            ind_vert(i,j,k) = count_vert
        end if
        indx = floor( ((/i,j,k/)*xinc)/xinc_coarse )
        indx = (/ min(nstep_coarse(1)-2, indx(1)), min(nstep_coarse(2)-2, indx(2)), &
            min(nstep_coarse(3)-2, indx(3))/)
        if ( .not.(rmbox_coarse(indx(1), indx(2), indx(3))) ) then
            count_cube = count_cube+1
        end if
      end do
    enddo
  enddo

  ! write data
  ! .mesh file
  write(lumesh,"(A)") 'MeshVersionFormatted 2'
  write(lumesh,"(A,/)") 'Dimension 3'
  write(lumesh,"(A)") 'Vertices'
  write(lumesh,'(I20)') count(vert_use)
  ! .sol file
  write(lusol,"(A)") 'MeshVersionFormatted 2'
  write(lusol,"(A,/)") 'Dimension 3'
  write(lusol,"(A)") 'SolAtVertices'
  write(lusol,'(I20)') count(vert_use)
  write(lusol,'(I16, I16)') 1, 1
  do i = 0, nstep(1)-1
    do j = 0, nstep(2)-1
      do k = 0, nstep(3)-1
        if (vert_use(i, j, k)) then
            write(lumesh,'(1x,3(e16.6), I8)') xinit+(/i,j,k/)*xinc, 1
            write(lusol,'(1x, 1(f16.6))') cgrad(i,j,k)
        end if
      end do
    enddo
  enddo

if ( ohexa ) then
      write(lumesh,"(/,A)") 'Hexaedra '
      write(lumesh,"(I20)") count_cube
      do i = 0, nstep(1)-2
        do j = 0, nstep(2)-2
          do k = 0, nstep(3)-2
            i1 = (/i, i+1/)
            j1 = (/j, j+1/)
            k1 = (/k, k+1/)
            c = count(vert_use(i1, j1, k1))
            if (c .eq. 8) then
              ind_cube = (/ind_vert(i,j,k), ind_vert(i,j,k+1), ind_vert(i,j+1,k+1), ind_vert(i,j+1,k), &
                ind_vert(i+1,j,k), ind_vert(i+1,j,k+1), ind_vert(i+1,j+1,k+1), ind_vert(i+1,j+1,k) /)
              write(lumesh,"(1x, 8(I12), I12)") ind_cube, 1
            end if
          end do
        end do
      end do
  end if

  ! save tetras, by diving each cube into six tetras
  if ( otetra ) then
      write(lumesh,"(/,A)") 'Tetrahedra'
      write(lumesh,"(I20)") 6*count_cube
      do i = 0, nstep(1)-2
        do j = 0, nstep(2)-2
          do k = 0, nstep(3)-2
            i1 = (/i, i+1/)
            j1 = (/j, j+1/)
            k1 = (/k, k+1/)
            c = count(vert_use(i1, j1, k1))
            if (c .eq. 8) then
              ind_cube = (/ind_vert(i,j,k), ind_vert(i,j,k+1), ind_vert(i,j+1,k+1), ind_vert(i,j+1,k), &
                ind_vert(i+1,j,k), ind_vert(i+1,j,k+1), ind_vert(i+1,j+1,k+1), ind_vert(i+1,j+1,k) /)
              call cube2tetra(tetra_cube, ind_cube)
              do m = 1,6
                write(lumesh,"(1x, 4(I16), I16)") tetra_cube(m,:), 1
              end do
            end if
          end do
        end do
      end do
  end if

write(lumesh, '(A)') 'END'
  write(lusol, '(A)') 'END'
  end subroutine write_mesh_file

! write six faces of a cube
  subroutine cube2quad(quad_cube, ind_cube)
  integer, intent(in) :: ind_cube(1:8)
  integer, intent(inout) :: quad_cube(1:6, 1:4)
  integer :: i, i1, i2, i3, i4, i5, i6, i7, i8
  i1 = ind_cube(1)
  i2 = ind_cube(2)
  i3 = ind_cube(3)
  i4 = ind_cube(4)
  i5 = ind_cube(5)
  i6 = ind_cube(6)
  i7 = ind_cube(7)
  i8 = ind_cube(8)
  quad_cube = reshape( (/i1, i2, i3, i4, &
                                    i5, i6, i7, i8, &
                                    i1, i2, i6, i5, &
                                    i4, i3, i7, i8, &
                                    i1, i5, i8, i4, &
                                    i2, i6, i7, i3 /), (/6,4/), order = (/2,1/) )
  end subroutine cube2quad

  ! divide a cube into six triangles, C Quan
  subroutine cube2tetra(tetra_cube, ind_cube)
  integer, intent(in) :: ind_cube(1:8)
  integer, intent(inout) :: tetra_cube(1:6, 1:4)
  integer :: i, i1, i2, i3, i4, i5, i6, i7, i8
  i1 = ind_cube(1)
  i2 = ind_cube(2)
  i3 = ind_cube(3)
  i4 = ind_cube(4)
  i5 = ind_cube(5)
  i6 = ind_cube(6)
  i7 = ind_cube(7)
  i8 = ind_cube(8)
  tetra_cube = reshape( (/i1, i2, i3, i7, &
                                    i1, i2, i6, i7, &
                                    i1, i4, i3, i7, &
                                    i1, i4, i8, i7, &
                                    i1, i5, i6, i7, &
                                    i1, i5, i8, i7 /), (/6,4/), order = (/2,1/) )
  end subroutine cube2tetra

! compute NCI geometry of the region enclosed by RDG isosurface, C Quan
  subroutine dataGeom( sum_rhon_vol, sum_rhon_area, sum_signrhon_vol, xinc, nstep, crho, crho_n, cgrad, cheig, rmbox_coarse,nfiles)
  real*8, intent(inout) :: sum_rhon_vol(9)
  real*8, intent(inout) :: sum_rhon_area(9)  
  real*8, intent(inout) :: sum_signrhon_vol(7) 
  real*8, intent(in) :: xinc(3)
  integer, intent(in) :: nstep(3)
  real*8, intent(in) :: crho(0:nstep(1)-1,0:nstep(2)-1,0:nstep(3)-1), &
        crho_n(0:nstep(1)-1,0:nstep(2)-1,0:nstep(3)-1,1:nfiles), &
        cgrad(0:nstep(1)-1,0:nstep(2)-1,0:nstep(3)-1), &
        cheig(0:nstep(1)-1,0:nstep(2)-1,0:nstep(3)-1)
  logical, intent(in) :: rmbox_coarse(0:nstep(1)-2, 0:nstep(2)-2, 0:nstep(3)-2)
  integer :: i, j, k, i0, j0, k0, n
  integer, intent(in) :: nfiles
  integer :: i1(2), j1(2), k1(2),negative,positive
  real*8 :: sum_signrho,signlambda_2
 

  sum_rhon_vol = 0 
  sum_signrhon_vol = 0
  negative=0
  positive=0
  ! integral of rho^n over the volume of cubes
  do i = 0, nstep(1)-2
    do j = 0, nstep(2)-2
      do k = 0, nstep(3)-2
        if ( .not. rmbox_coarse(i,j,k) ) then
            i1 = (/i, i+1/)
            j1 = (/j, j+1/)
            k1 = (/k, k+1/)
            ! n = 1, 1.5, 2, 2.5, 3, 4/3, 5/3: sum of rho^n
            sum_rhon_vol(1) = sum_rhon_vol(1) + sum(abs(crho(i1,j1,k1)/100)) *xinc(1)*xinc(2)*xinc(3) /8
            sum_rhon_vol(2) = sum_rhon_vol(2) + sum(abs(crho(i1,j1,k1)/100)**1.5) &
            *xinc(1)*xinc(2)*xinc(3) /8
            sum_rhon_vol(3) = sum_rhon_vol(3) + sum(abs(crho(i1,j1,k1)/100)**2) &
            *xinc(1)*xinc(2)*xinc(3) /8
            sum_rhon_vol(4) = sum_rhon_vol(4) + sum(abs(crho(i1,j1,k1)/ &
            100)**2.5)*xinc(1)*xinc(2)*xinc(3) /8
            sum_rhon_vol(5) = sum_rhon_vol(5) + sum(abs(crho(i1,j1,k1)/100)**3) &
            *xinc(1)*xinc(2)*xinc(3) /8
            sum_rhon_vol(6) = sum_rhon_vol(6) + sum(abs(crho(i1,j1,k1)/100)**(1.333)) &
            *xinc(1)*xinc(2)*xinc(3) /8
            sum_rhon_vol(7) = sum_rhon_vol(7) + sum(abs(crho(i1,j1,k1)/100)**(1.666)) &
            *xinc(1)*xinc(2)*xinc(3) /8
            ! n = 0: volume of cubes
            sum_rhon_vol(8) = sum_rhon_vol(8) + xinc(1)*xinc(2)*xinc(3) 
            ! sum of rho1*rho2
            do n = 1,nfiles
                sum_rhon_vol(9) = sum_rhon_vol(9) + sum(crho_n(i1,j1,k1,n)/100*&
                    (abs(crho(i1,j1,k1)) - crho_n(i1,j1,k1,n))/100) *xinc(1)*xinc(2)*xinc(3) /8
            end do 
             
            !sign_lambda2 RHO
            sum_signrho=sum(crho(i1,j1,k1))
            signlambda_2=sign(1d0,sum_signrho) 

        
            sum_signrhon_vol(1) = sum_signrhon_vol(1) + sum(crho(i1,j1,k1)/100) *xinc(1)*xinc(2)*xinc(3) /8
            sum_signrhon_vol(2) = sum_signrhon_vol(2) + signlambda_2*sum(abs(crho(i1,j1,k1)/100)**1.5) &
            *xinc(1)*xinc(2)*xinc(3) /8
            sum_signrhon_vol(3) = sum_signrhon_vol(3) + signlambda_2*sum(abs(crho(i1,j1,k1)/100)**2) &
            *xinc(1)*xinc(2)*xinc(3) /8
            sum_signrhon_vol(4) = sum_signrhon_vol(4) + signlambda_2*sum(abs(crho(i1,j1,k1)/ &
            100)**2.5)*xinc(1)*xinc(2)*xinc(3) /8
            sum_signrhon_vol(5) = sum_signrhon_vol(5) + signlambda_2*sum(abs(crho(i1,j1,k1)/100)**3) &
            *xinc(1)*xinc(2)*xinc(3) /8
            sum_signrhon_vol(6) = sum_signrhon_vol(6) + signlambda_2*sum(abs(crho(i1,j1,k1)/100)**(1.333)) &
            *xinc(1)*xinc(2)*xinc(3) /8
            sum_signrhon_vol(7) = sum_signrhon_vol(7) + signlambda_2*sum(abs(crho(i1,j1,k1)/100)**(1.666)) &
            *xinc(1)*xinc(2)*xinc(3) /8
        end if
      end do
    end do
  end do


  ! integral of rho^n over the surface of cubes
  sum_rhon_area = 0
  do i = 0, nstep(1)-2
    do j = 0, nstep(2)-2
      do k = 0, nstep(3)-2
        if ( .not. rmbox_coarse(i,j,k) ) then
            if ( (i .eq. 0) .or. (rmbox_coarse(i-1,j,k)) )then
                 i1 = (/i, i/)
                 j1 = (/j, j+1/)
                 k1 = (/k, k+1/)
                 call compArea(sum_rhon_area, nstep, crho, crho_n, nfiles, i1, j1, k1, xinc(2), xinc(3))
            elseif ( (i .eq. nstep(1)-2) .or. (rmbox_coarse(i+1,j,k)) ) then
                 i1 = (/i+1, i+1/)
                 j1 = (/j, j+1/)
                 k1 = (/k, k+1/)
                 call compArea(sum_rhon_area, nstep, crho, crho_n, nfiles, i1, j1, k1, xinc(2), xinc(3))
            end if
            if ( (j .eq. 0) .or. (rmbox_coarse(i,j-1,k)) )then
                 i1 = (/i, i+1/)
                 j1 = (/j, j/)
                 k1 = (/k, k+1/)
                 call compArea(sum_rhon_area, nstep, crho, crho_n, nfiles, i1, j1, k1, xinc(1), xinc(3))
            elseif ( (j .eq. nstep(2)-2) .or. (rmbox_coarse(i,j+1,k)) ) then
                 i1 = (/i, i+1/)
                 j1 = (/j+1, j+1/)
                 k1 = (/k, k+1/)
                 call compArea(sum_rhon_area, nstep, crho, crho_n, nfiles, i1, j1, k1, xinc(1), xinc(3))
            end if
            if ( (k .eq. 0) .or. (rmbox_coarse(i,j,k-1)) )then
                 i1 = (/i, i+1/)
                 j1 = (/j, j+1/)
                 k1 = (/k, k/)
                 call compArea(sum_rhon_area, nstep, crho, crho_n, nfiles, i1, j1, k1, xinc(1), xinc(2))
            elseif ( (k .eq. nstep(3)-2) .or. (rmbox_coarse(i,j,k+1)) ) then
                 i1 = (/i, i+1/)
                 j1 = (/j, j+1/)
                 k1 = (/k+1, k+1/)
                 call compArea(sum_rhon_area, nstep, crho, crho_n, nfiles, i1, j1, k1, xinc(1), xinc(2))
            end if
        end if
      end do
    end do
  end do

  end subroutine dataGeom

  subroutine compArea(sum_rhon_area, nstep, crho, crho_n, nfiles, i1, j1, k1, a, b)
  ! compute integrals over the boundary surface of the union of cubes, C Quan
  real*8, intent(inout) :: sum_rhon_area(9)
  real*8, intent(in) :: a, b
  integer, intent(in) :: nstep(3)
  real*8, intent(in) :: crho(0:nstep(1)-1,0:nstep(2)-1,0:nstep(3)-1), &
        crho_n(0:nstep(1)-1,0:nstep(2)-1,0:nstep(3)-1,1:nfiles)
  integer :: n
  integer, intent(in) :: nfiles
  integer, intent(in) :: i1(2), j1(2), k1(2)
  ! face sides (a, b)
  ! n = 1, 1.5, 2, 2.5, 3, 4/3, 5/3: sum of rho^n
    sum_rhon_area(1) = sum_rhon_area(1) + sum(abs(crho(i1,j1,k1)/100)) *a*b /8
    sum_rhon_area(2) = sum_rhon_area(2) + sum(abs(crho(i1,j1,k1)/100)**1.5) *a*b /8
    sum_rhon_area(3) = sum_rhon_area(3) + sum(abs(crho(i1,j1,k1)/100)**2) *a*b /8
    sum_rhon_area(4) = sum_rhon_area(4) + sum(abs(crho(i1,j1,k1)/100)**2.5)*a*b /8
    sum_rhon_area(5) = sum_rhon_area(5) + sum(abs(crho(i1,j1,k1)/100)**3) *a*b /8
    sum_rhon_area(6) = sum_rhon_area(6) + sum(abs(crho(i1,j1,k1)/100)**(1.333)) &
        *a*b /8
    sum_rhon_area(7) = sum_rhon_area(7) + sum(abs(crho(i1,j1,k1)/100)**(1.666)) &
        *a*b /8
    ! area of iso-surface
    sum_rhon_area(8) = sum_rhon_area(8) + a*b
    ! sum of rho1*rho2
    do n = 1,nfiles
        sum_rhon_area(9) = sum_rhon_area(9) + sum(crho_n(i1,j1,k1,n)/100*&
            (abs(crho(i1,j1,k1)) - crho_n(i1,j1,k1,n))/100) *a*b /8
    end do
  end subroutine compArea


end program
