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

! props: Computing electron density, its derivative and DFT functional Computing electron density, its derivative and DFT functionalss 

module props
  use reader, only: grid1
  implicit none

  public

  private :: promolecular_grid, index0, pri012, phi012, propwfn, promolecular
  private :: grid1_interp
  private :: gascor, perdc

  integer :: iztype(120,0:20)
  type(grid1), allocatable :: grd(:)

contains

   subroutine calcprops_wfn(xinit,xinc,n,mol,nmol,rho,grad,doelf,elf,ixc,xc,ctp,cheig)
    use reader
    use tools_io
    use tools_math
    use param

    real*8, intent(in) :: xinit(3), xinc(3)
    integer, intent(in) :: n(3),nmol
    type(molecule) :: mol(nmol) 
    real*8, dimension(n(1),n(2),n(3)), intent(out) :: rho, grad 
    logical, intent(in) :: doelf
    real*8, dimension(n(1),n(2),n(3)), intent(out) :: elf, xc 
    integer, intent(in) :: ixc(2) 
!   real*8, dimension(3,n(1),n(2),n(3)), intent(out) :: cheig
    
    real*8, dimension(n(1),n(2),n(3)), intent(out) :: ctp,cheig
    logical :: doxc
    real*8, allocatable, dimension(:,:) :: dx, dy, dz, d2, gg
    real*8, allocatable, dimension(:) :: tp, maxc, rhoaux
    integer :: nmcent, istat
    real*8, allocatable :: chi(:,:), phi(:,:), hess(:,:,:)
    logical, allocatable :: ldopri(:,:)

    integer :: i, j, k, m, iat, ip, jp, kp, nn, l(3)
    integer :: ityp, ipri, ipria, ix, imo
    real*8 :: ex, xl2, xl(3,0:2), x0(3), al
    real*8 :: wk1(3), wk2(3), hvecs(3,3), grad2,heigs(3)
    real*8 :: rho53, ebose, df, eelf, eexc

    real*8, parameter :: cutoff_pri = 1d-10
    real*8, parameter :: fothirds = 4d0/3d0
    real*8, parameter :: m53 = -5d0/3d0
    real*8, parameter :: ctelf = 10d0*3d0**(-5d0/3d0)*pi**(-4d0/3d0)
    

    rho = 0d0 

    grad = 0d0 

   ! ctp= 0d0  
   !cheig = 0d0 
   !RAB: Uncomment cheigh does error. WHY? 
    cheig = 0d0  

 
    if (doelf) elf = 0d0
   ! doxc = any(ixc /= 0)
   ! if (doxc) xc = 0d0

    ! distance matrix, hessian and tp allocates 

    nmcent = 0
    do i = 1,nmol
      
       nmcent = max(nmcent,mol(i)%n) 
    enddo 
    




    allocate(dx(n(1),nmcent),dy(n(1),nmcent),dz(n(1),nmcent),d2(n(1),nmcent),stat=istat)
    if (istat /= 0) call error('calcprops','error allocating distance matrix',faterr)
    allocate(hess(n(1),3,3),stat=istat)
    if (istat /= 0) call error('calcprops','error allocating hessian',faterr)
    allocate(tp(n(1)),rhoaux(n(1)),stat=istat)
    if (istat /= 0) call error('calcprops','error allocating rhoaux',faterr)
    allocate(gg(n(1),3),stat=istat)
    if (istat /= 0) call error('calcprops','error allocating gg',faterr)

    ! run over y and z
    !$omp parallel do private (ip,jp,kp,hess,tp,gg,chi,phi,maxc,ldopri,dx,dy,dz,d2,&
    !$omp nn,ipri,iat,al,ex,x0,l,xl,xl2,grad2,rho53,ebose,df,hvecs,wk1,wk2,&
    !$omp istat,rhoaux,eelf,eexc) schedule(dynamic) 

    do j = 0, n(2)-1
       jp = j + 1
       do k = 0, n(3)-1
          kp = k + 1
          ! zero in-line hessian and tp
          hess = 0d0
          rhoaux = 0d0
          tp = 0d0
          gg = 0d0

          ! run over molecules 
          ! Original 
           do m = 1, nmol 

             ! allocate primitive evaluation array
             allocate(chi(mol(m)%npri,10),phi(mol(m)%nmo,10))

             ! identify the max coefficient 
             allocate(maxc(mol(m)%npri),ldopri(mol(m)%npri,10))
             maxc = 0d0
             do imo = 1, mol(m)%nmo 
                do ipri = 1, mol(m)%npri 
                   maxc(ipri) = max(maxc(ipri),abs(mol(m)%c(imo,ipri)))                   
                enddo
             enddo

             ! calculate distances  
         
             do iat = 1, mol(m)%n
                do i = 0, n(1)-1 
                   ip = i + 1
                   dx(ip,iat) = xinit(1) + i * xinc(1) - mol(m)%x(1,iat)
                   dy(ip,iat) = xinit(2) + j * xinc(2) - mol(m)%x(2,iat) 
                   dz(ip,iat) = xinit(3) + k * xinc(3) - mol(m)%x(3,iat)

                   d2(ip,iat) = dx(ip,iat)*dx(ip,iat)+dy(ip,iat)*dy(ip,iat)+dz(ip,iat)*dz(ip,iat)
                enddo
             enddo
             
             
             ! calculate primitives at the points 
             do i = 0, n(1)-1
                ip = i + 1
                nn = 0
                do ityp = 1, mol(m)%maxntyp
                   do ipria = nn+1,nn+mol(m)%ntyp(ityp) 
                      ipri = mol(m)%intyp(ipria)
                      iat = mol(m)%icenter(ipri)
                      al = mol(m)%e(ipri)
                      ex = exp(-al * d2(ip,iat))
            
                      x0 = (/ dx(ip,iat), dy(ip,iat), dz(ip,iat) /)
            
                      call index0(ityp,l)
                      do ix = 1, 3
                         if (l(ix) == 0) then
                            xl(ix,0) = 1d0
                            xl(ix,1) = 0d0
                            xl(ix,2) = 0d0
                         else if (l(ix) == 1) then
                            xl(ix,0) = x0(ix)
                            xl(ix,1) = 1d0
                            xl(ix,2) = 0d0
                         else if (l(ix) == 2) then
                            xl(ix,0) = x0(ix) * x0(ix)
                            xl(ix,1) = 2d0 * x0(ix)
                            xl(ix,2) = 2d0
                         else if (l(ix) == 3) then
                            xl(ix,0) = x0(ix) * x0(ix) * x0(ix)
                            xl(ix,1) = 3d0 * x0(ix) * x0(ix)
                            xl(ix,2) = 6d0 * x0(ix)
                         else if (l(ix) == 4) then
                            xl2 = x0(ix) * x0(ix)
                            xl(ix,0) = xl2 * xl2
                            xl(ix,1) = 4d0 * xl2 * x0(ix)
                            xl(ix,2) = 12d0 * xl2
                         else
                            call error('pri012','power of L not supported',faterr)
                         end if
                      end do
                      
                      chi(ipri,1) = xl(1,0)*xl(2,0)*xl(3,0)*ex
                      chi(ipri,2) = (xl(1,1)-2*al*x0(1)**(l(1)+1))*xl(2,0)*xl(3,0)*ex
                      chi(ipri,3) = (xl(2,1)-2*al*x0(2)**(l(2)+1))*xl(1,0)*xl(3,0)*ex
                      chi(ipri,4) = (xl(3,1)-2*al*x0(3)**(l(3)+1))*xl(1,0)*xl(2,0)*ex
                      chi(ipri,5) = (xl(1,2)-2*al*(2*l(1)+1)*xl(1,0)&
                         +4*al*al*x0(1)**(l(1)+2))*xl(2,0)*xl(3,0)*ex
                      chi(ipri,6) = (xl(2,2)-2*al*(2*l(2)+1)*xl(2,0)&
                         +4*al*al*x0(2)**(l(2)+2))*xl(3,0)*xl(1,0)*ex
                      chi(ipri,7) = (xl(3,2)-2*al*(2*l(3)+1)*xl(3,0)&
                         +4*al*al*x0(3)**(l(3)+2))*xl(1,0)*xl(2,0)*ex
                      chi(ipri,8) = (xl(1,1)-2*al*x0(1)**(l(1)+1))*&
                         (xl(2,1)-2*al*x0(2)**(l(2)+1))*xl(3,0)*ex
                      chi(ipri,9) = (xl(1,1)-2*al*x0(1)**(l(1)+1))*&
                         (xl(3,1)-2*al*x0(3)**(l(3)+1))*xl(2,0)*ex
                      chi(ipri,10)= (xl(3,1)-2*al*x0(3)**(l(3)+1))*&
                         (xl(2,1)-2*al*x0(2)**(l(2)+1))*xl(1,0)*ex

                      do ix = 1, 10
                         ldopri(ipri,ix) = (abs(chi(ipri,ix))*maxc(ipri) > cutoff_pri)
                      enddo
                   enddo ! ipria = nn+1, nn+ntyp

                   nn = nn + mol(m)%ntyp(ityp)
                enddo ! ityp = 1, maxntyp
                
                ! build the MO avlues at the point 
                phi = 0d0
                do ix = 1, 10
                   do ipri = 1, mol(m)%npri
                      if (.not.ldopri(ipri,ix)) cycle
                      do imo = 1, mol(m)%nmo
                         phi(imo,ix) = phi(imo,ix) + mol(m)%c(imo,ipri)*chi(ipri,ix)
                      enddo
                   enddo
                enddo
                
                ! contribution to the density, etc.
                do imo = 1, mol(m)%nmo
                   rhoaux(ip) = rhoaux(ip) + mol(m)%occ(imo) * phi(imo,1) * phi(imo,1)
                   gg(ip,1) = gg(ip,1) + 2 * mol(m)%occ(imo) * phi(imo,1) * phi(imo,2)
                   gg(ip,2) = gg(ip,2) + 2 * mol(m)%occ(imo) * phi(imo,1) * phi(imo,3)
                   gg(ip,3) = gg(ip,3) + 2 * mol(m)%occ(imo) * phi(imo,1) * phi(imo,4)
                   hess(ip,1,1) = hess(ip,1,1) + 2 * mol(m)%occ(imo) * (phi(imo,1)*phi(imo,5)+phi(imo,2)**2)
                   hess(ip,2,2) = hess(ip,2,2) + 2 * mol(m)%occ(imo) * (phi(imo,1)*phi(imo,6)+phi(imo,3)**2)
                   hess(ip,3,3) = hess(ip,3,3) + 2 * mol(m)%occ(imo) * (phi(imo,1)*phi(imo,7)+phi(imo,4)**2)
                   hess(ip,1,2) = hess(ip,1,2) + 2 * mol(m)%occ(imo) * (phi(imo,1)*phi(imo,8)+phi(imo,2)*phi(imo,3))
                   hess(ip,1,3) = hess(ip,1,3) + 2 * mol(m)%occ(imo) * (phi(imo,1)*phi(imo,9)+phi(imo,2)*phi(imo,4))
                   hess(ip,2,3) = hess(ip,2,3) + 2 * mol(m)%occ(imo) * (phi(imo,1)*phi(imo,10)+phi(imo,3)*phi(imo,4))
                   tp(ip) = tp(ip) + mol(m)%occ(imo) * (phi(imo,2)*phi(imo,2)+phi(imo,3)*phi(imo,3)+phi(imo,4)*phi(imo,4))
                enddo
             enddo ! i = 0, n-1

             deallocate(chi,phi,maxc,ldopri)

          enddo ! m = 1, nmol 
          
        

          ! accumulate intermediate variables
          do ip = 1, n(1)
             ! rho and grad
             rhoaux(ip) = max(rhoaux(ip),1d-30)
             grad2 = gg(ip,1)**2 + gg(ip,2)**2 + gg(ip,3)**2

             ! tau and elf
             tp(ip) = 0.5d0*tp(ip)
             if (doelf) then
                rho53 = ctelf * rhoaux(ip)**m53
                ebose = 0.125d0 * grad2 / rhoaux(ip)
                df = rho53 * (tp(ip)-ebose)
                eelf = 1d0 / (1d0 + df*df)
             endif

             ! xc, ixc
            ! if (doxc) then
            !    call calcexc(ixc,eexc,rhoaux(ip),gg(ip,:),tp(ip))
            ! endif

             ! hessian eigenvalue and sign of rho
             hess(ip,2,1) = hess(ip,1,2)
             hess(ip,3,1) = hess(ip,1,3)
             hess(ip,3,2) = hess(ip,2,3)
             call rs(3,3,hess(ip,:,:),heigs,0,hvecs,wk1,wk2,istat)

             !$omp critical (writeshared)
             rho(ip,jp,kp) = sign(rhoaux(ip),heigs(2)) * 100d0
             grad(ip,jp,kp) = sqrt(grad2) / (const * rhoaux(ip)**fothirds) 

             ! Roberto kinetic energy density
             ! ctp(ip,jp,kp)=tp(ip)
             ! Roberto save eigenvalues
         !!    cheig(1,ip,jp,kp)=heigs(1)
         !!     cheig(2,ip,jp,kp)=heigs(2)
              cheig(ip,jp,kp)=heigs(2)
            ! cheig(3,ip,jp,kp)=heigs(3)
             if (doelf) elf(ip,jp,kp) = eelf
            ! if (doxc) xc(ip,jp,kp) = eexc
             !$omp end critical (writeshared)
          enddo
       enddo ! k = 0, n(3)-1
    enddo ! j = 0, n(2)-1
    !$omp end parallel do

    deallocate(dx,dy,dz,d2,hess,tp,gg)
   

  end subroutine calcprops_wfn

  subroutine init_rhogrid(m,nm)
    use reader
    use tools_io
    use param

    type(molecule), intent(in) :: m(nm)
    integer, intent(in) :: nm

    integer :: i, j, iz, iq, nat, istat
    logical :: isat(120,0:20)
    
    isat = .false. 


    do i = 1, nm
       do j = 1, m(i)%n  
          iz = m(i)%z(j) 
          iq = m(i)%q(j)   
          if (iz > size(isat,1)) call error('init_rhogrid','isat size exceeded',faterr)
          if (iq > size(isat,2)) call error('init_rhogrid','isat size exceeded',faterr)
          isat(iz,iq) = .true.
       end do
    end do
    

    nat = count(isat)
    allocate(grd(nat),stat=istat)
    if (istat /= 0) call error('nciplot','could not allocate memory for atomic density grids',faterr)

    nat = 0
    do i = 1, size(isat,1)
       do j = 0, size(isat,2)-1
          if (isat(i,j)) then
             nat = nat + 1
             iztype(i,j) = nat 
             call readgrid(grd(nat),i-j)
          end if
       end do
    end do

  end subroutine init_rhogrid

  subroutine calcprops_pro(x,m,nm,rho,rho_n,rhom,nfr,autofr,grad,hess,doelf,elf,ixc,exc,deltag)
    use reader
    use tools_io
    use param

    real*8, intent(in) :: x(3)
    integer, intent(in) :: nm, nfr
    logical, intent(in) :: autofr
    type(molecule) :: m(nm)
    real*8, intent(out) :: rho, rho_n(1:nm), rhom(nfr),grad(3), hess(3,3)
    logical, intent(in) :: doelf
    integer, intent(in) :: ixc(2)
    real*8, intent(out) :: elf, exc,deltag

    integer :: i,j
    real*8 :: hh(3,3), gg(3), rr, tpp, rho53, grad2, ebose, df, rhomr(nfr),tp
    real*8 :: igmgg(3), gradigm(3), gradigm2

    real*8, parameter :: fthird = -5d0/3d0
    real*8, parameter :: ctelf = 10d0*3d0**(-5d0/3d0)*pi**(-4d0/3d0)

    rho = 0d0
    rhom = 0d0
    grad = 0d0
    gradigm = 0d0
    hess = 0d0
    tp=0d0
    do i = 1, nm
       select case (m(i)%ifile)
       case (ifile_xyz)
          ! promolecular densities
          call promolecular(x,m(i),rr,gg,hh,rhomr,nfr,autofr,igmgg)
          tpp = 0d0
       case (ifile_grd,ifile_wfn,ifile_wfx)
          call promolecular_grid(x,m(i),rr,gg,hh,rhomr,nfr,autofr,.false.)
          tpp = 0d0
       case default
          call error('calcprops','unknown molecule type',faterr)
       end select
       rho = rho + rr 
       rho_n(i)= rr
       if (autofr) then
          rhom(i) = rhom(i) + rr 
       else
          rhom = rhom + rhomr
       end if
       grad = grad + gg 
       hess = hess + hh  

       ! IGM model
       gradigm  = gradigm +igmgg
    end do 

    ! IGM model
    gradigm2 = gradigm(1)*gradigm(1)+gradigm(2)*gradigm(2)+gradigm(3)*gradigm(3) 
    grad2 = grad(1)*grad(1)+grad(2)*grad(2)+grad(3)*grad(3)
    deltag = sqrt(gradigm2)-sqrt(grad2)  

    tp=tp+tpp
    exc = 0d0
    if (all(ixc /= 0)) call calcexc(ixc,exc,rho,grad,tp)
    
    elf = 0d0
    if (doelf) then
       rho53 = ctelf * rho**fthird
       grad2 = grad(1)*grad(1)+grad(2)*grad(2)+grad(3)*grad(3)
       ebose = 0.125d0 * grad2 / rho
       df = rho53 * (tp-ebose)
       elf = 1d0 / (1d0 + df*df)
    end if

  end subroutine calcprops_pro
!< Calculate promolecular densities for indiviual files
  subroutine calcprops_proref(x,m,nm,ifile,rho,dimgrad,nfr,autofr,grad,hess,doelf,elf,ixc,exc)
    use reader
    use tools_io
    use param

    real*8, intent(in) :: x(3)
    integer, intent(in) :: nm,ifile,nfr
    logical, intent(in) :: autofr
    type(molecule) :: m(nm)
    real*8, intent(inout) :: rho,dimgrad,grad(3), hess(3,3)
    logical, intent(in) :: doelf
    integer, intent(in) :: ixc(2)
    real*8, intent(inout) :: elf, exc

    integer :: i
    real*8 :: hh(3,3), gg(3), rr, tpp, rho53, grad2, ebose, df, rhomr(nfr),tp,igmgg(3)

    real*8, parameter :: fthird = -5d0/3d0
    real*8, parameter :: ctelf = 10d0*3d0**(-5d0/3d0)*pi**(-4d0/3d0)

    rho = 0d0
    grad = 0d0
    hess = 0d0
    tp=0d0
    rr=0d0
    gg=0d0
    hh=0d0
    tp=0.d0
 
!    do i = 1, nm
       select case (m(ifile)%ifile)
       case (ifile_xyz)
          ! promolecular densities
          call promolecular(x,m(ifile),rr,gg,hh,rhomr,nfr,autofr,igmgg)
!          tpp = 0d0
 !      case (ifile_grd,ifile_wfn)
 !         call promolecular_grid(x,m(ifile),rr,gg,hh,rhomr,nfr,autofr,.false.)
 !         tpp = 0d0
       case default
          call error('calcprops','unknown molecule type',faterr)
       end select
       rho = rho + rr
       if (autofr) then
          rho = rho + rr
          grad= grad+ sum(gg(:))
 !      else
 !         rhom = rhom + rhomr
       end if
       grad = grad + gg
       hess = hess + hh
!    end do 
 !   tp=tp+tpp
    exc = 0d0
    if (all(ixc /= 0)) call calcexc(ixc,exc,rho,grad,tp)

    elf = 0d0
    if (doelf) then
       rho53 = ctelf * rho**fthird
       grad2 = grad(1)*grad(1)+grad(2)*grad(2)+grad(3)*grad(3)
       ebose = 0.125d0 * grad2 / rho
       df = rho53 * (tp-ebose)
       elf = 1d0 / (1d0 + df*df)
    end if

  end subroutine calcprops_proref 

  !> Calculate the density and gradient using atomic densities
  subroutine promolecular(x,m,rho,grad,hess,rhom,nfr,autofr,gradigm)
    use reader
    use tools_io
    use param

    real*8, intent(in) :: x(3)
    type(molecule), intent(in) :: m
    integer, intent(in) :: nfr
    real*8, intent(out) :: rho, grad(3), hess(3,3), rhom(nfr)
    logical, intent(in) :: autofr
    real*8, intent(out) :: gradigm(3)

    ! The C and ZETA data blocks are the coefficients and exponents
    ! for the atoms in the first three row of the periodic table, H-Ar.
    !    rho_atom = sum_i c_i * exp(-r / z_i) where i = 1..3
    ! The index runs over atomic numbers
    real*8, parameter :: c1(18) = (/&
       0.2815D0, 2.437D0, 11.84D0, 31.34D0, 67.82D0, 120.2D0, 190.9D0,&
       289.5D0,  406.3D0, 561.3D0, 760.8D0, 1016.D0, 1319.D0, 1658.D0,&
       2042.D0, 2501.D0, 3024.D0, 3625.D0/)
    real*8, parameter :: c2(18) = (/&
         0.D0,    0.D0, 0.06332D0, 0.3694D0, 0.8527D0, 1.172D0, 2.247D0,&
      2.879D0, 3.049D0,   6.984D0,  22.42D0,  37.17D0, 57.95D0, 87.16D0,&
      115.7D0, 158.0D0,   205.5D0,  260.0D0/)
    real*8, parameter :: c3(18) = (/&
       0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0,&
       0.06358D0, 0.3331D0, 0.8878D0, 0.7888D0, 1.465D0, 2.170D0,&
       3.369D0, 5.211D0/)
    real*8, parameter :: zeta1(18) = (/&
       0.5288D0, 0.3379D0, 0.1912D0, 0.1390D0, 0.1059D0, 0.0884D0,&
       0.0767D0, 0.0669D0, 0.0608D0, 0.0549D0, 0.0496D0, 0.0449D0,&
       0.0411D0, 0.0382D0, 0.0358D0, 0.0335D0, 0.0315D0, 0.0296D0/)
    real*8, parameter :: zeta2(18) = (/&
       1.D0, 1.D0, 0.9992D0, 0.6945D0, 0.5300D0, 0.5480D0,&
       0.4532D0, 0.3974D0, 0.3994D0, 0.3447D0, 0.2511D0, 0.2150D0,&
       0.1874D0, 0.1654D0, 0.1509D0, 0.1369D0, 0.1259D0, 0.1168D0/)
    real*8, parameter :: zeta3(18) = (/&
       1.D0, 1.D0, 1.D0, 1.D0, 1.D0, 1.D0, 1.D0, 1.D0, 1.D0, 1.D0,&
       1.0236D0, 0.7753D0, 0.5962D0, 0.6995D0, 0.5851D0, 0.5149D0,&
       0.4974D0, 0.4412D0/)
    real*8, parameter :: czeta1(18) = c1 / zeta1
    real*8, parameter :: czeta2(18) = c2 / zeta2
    real*8, parameter :: czeta3(18) = c3 / zeta3
    real*8, parameter :: czzeta1(18) = czeta1 / zeta1
    real*8, parameter :: czzeta2(18) = czeta2 / zeta2
    real*8, parameter :: czzeta3(18) = czeta3 / zeta3

    integer :: i, iz, j, k
    real*8 :: exp1, exp2, exp3, r, fac0, fac1, fac2, xd(3), xu(3), xu2(3), r1, r2
    real*8 :: fac3

    rho = 0d0
    rhom = 0d0
    grad = 0d0 
    gradigm = 0d0
    hess = 0d0
    do i = 1, m%n
       iz = m%z(i)
       if (iz < 1 .or. iz > atomic_zmax) &
          call error('promolecular','Atomic densities for high Z are not implemented (check atomic_zmax).',faterr)
       xd = x - m%x(:,i)
       r = sqrt(dot_product(xd,xd))
       r = max(r,1d-10)
       r1 = 1d0 / r
       r2 = r1 * r1
       xu = xd * r1
       xu2 = xu * xu

       exp1 = exp(-r/zeta1(iz))
       exp2 = exp(-r/zeta2(iz))
       exp3 = exp(-r/zeta3(iz))
       fac0 = c1(iz) * exp1 + c2(iz) * exp2 + c3(iz) * exp3
       fac1 = czeta1(iz)*exp1 + czeta2(iz)*exp2 + czeta3(iz)*exp3
       fac2 = czzeta1(iz)*exp1 + czzeta2(iz)*exp2 + czzeta3(iz)*exp3
       fac3 = fac2 + fac1 * r1

       rho = rho + fac0
       if (.not.autofr .and. m%ifrag(i)>0) rhom(m%ifrag(i)) = rhom(m%ifrag(i)) + fac0
       grad = grad - fac1 * xu 
       do j=1,3 
           gradigm(j) = gradigm(j) + sqrt(fac1*fac1*xu(j)*xu(j))
       enddo
       ! hessian diagonal elements
       do j = 1, 3
          hess(j,j) = hess(j,j) + fac3 * xu2(j) - fac1 * r1
          do k = j+1, 3
             hess(j,k) = hess(j,k) + r2 * xd(j) * xd(k) * fac3
          end do
       end do
    end do
    do j = 1, 3
       do k = j+1, 3
          hess(k,j) = hess(j,k)
       end do
    end do

  end subroutine promolecular

  subroutine propwfn(x,m,rho,g,hess,tp)
    use reader

    type(molecule), intent(in) :: m
    real*8, intent(in) :: x(3)
    real*8, intent(out) :: rho, g(3), hess(3,3), tp

    real*8 :: psi(m%nmo), f(10), store(m%nmo,10)
    integer :: i, imo

    call phi012(x,m,store)

    rho = 0d0
    tp = 0d0
    g = 0d0
    hess = 0d0
    do imo = 1, m%nmo
       f = store(imo,1:10)
       ! Value of orbital
       psi(imo)=f(1)
       ! Density
       rho=rho+m%occ(imo)*f(1)*f(1)
       ! Gradient
       g = g + 2 * m%occ(imo) * f(1) * f(2:4)
       ! Diagonal Hessian
       do i = 1, 3
          hess(i,i) = hess(i,i) + 2 * m%occ(imo) * (f(1)*f(4+i)+f(1+i)**2)
       end do
       ! Off-diagonal Hessian
       hess(1,2) = hess(1,2) + 2 * m%occ(imo) * (f(1)*f(8)+f(2)*f(3))
       hess(1,3) = hess(1,3) + 2 * m%occ(imo) * (f(1)*f(9)+f(2)*f(4))
       hess(2,3) = hess(2,3) + 2 * m%occ(imo) * (f(1)*f(10)+f(3)*f(4))
       ! Kinetic Engery
       tp = tp + m%occ(imo) * (f(2)*f(2)+f(3)*f(3)+f(4)*f(4))
    end do
    tp = 0.5d0 * tp
    hess(2,1) = hess(1,2)
    hess(3,1) = hess(1,3)
    hess(3,2) = hess(2,3)

  end subroutine propwfn

  subroutine phi012(x,m,phi)
    use reader

    type(molecule), intent(in) :: m
    real*8, intent(in) :: x(3)
    real*8, intent(out) :: phi(m%nmo,10)

    real*8 :: pri(m%npri,10)
    integer :: i, ipri, imo

    ! values of primitives
    call pri012(x,m,pri)

    ! sum
    phi = 0d0
    do i = 1, 10
       do ipri = 1, m%npri
          do imo = 1, m%nmo
             phi(imo,i) = phi(imo,i) + m%c(imo,ipri)*pri(ipri,i)
          end do
       end do
    end do

  end subroutine phi012

  subroutine pri012(x,m,pri)
    use reader
    use tools_io
    use param

    type(molecule), intent(in) :: m
    real*8, intent(in) :: x(3)
    real*8, intent(out) :: pri(m%npri,10)

    integer :: l(3), ipri, ic, it, i
    real*8 :: al, x0(3), x2, ex, xl(3,0:2), xl2

    do ipri = 1, m%npri
       ic = m%icenter(ipri)
       it = m%itype(ipri)

       al = m%e(ipri)
       x0 = x - m%x(:,ic)
       x2 = dot_product(x0,x0)
       ex = exp(-al * x2)

       call index0(it,l)
       do i = 1, 3
          if (l(i) == 0) then
             xl(i,0) = 1d0
             xl(i,1) = 0d0
             xl(i,2) = 0d0
          else if (l(i) == 1) then
             xl(i,0) = x0(i)
             xl(i,1) = 1d0
             xl(i,2) = 0d0
          else if (l(i) == 2) then
             xl(i,0) = x0(i) * x0(i)
             xl(i,1) = 2d0 * x0(i)
             xl(i,2) = 2d0
          else if (l(i) == 3) then
             xl(i,0) = x0(i) * x0(i) * x0(i)
             xl(i,1) = 3d0 * x0(i) * x0(i)
             xl(i,2) = 6d0 * x0(i)
          else if (l(i) == 4) then
             xl2 = x0(i) * x0(i)
             xl(i,0) = xl2 * xl2
             xl(i,1) = 4d0 * xl2 * x0(i)
             xl(i,2) = 12d0 * xl2
          else
             call error('pri012','power of L not supported',faterr)
          end if
       end do

       pri(ipri,1) = xl(1,0)*xl(2,0)*xl(3,0)*ex
       pri(ipri,2) = (xl(1,1)-2*al*x0(1)**(l(1)+1))*xl(2,0)*xl(3,0)*ex
       pri(ipri,3) = (xl(2,1)-2*al*x0(2)**(l(2)+1))*xl(1,0)*xl(3,0)*ex
       pri(ipri,4) = (xl(3,1)-2*al*x0(3)**(l(3)+1))*xl(1,0)*xl(2,0)*ex
       pri(ipri,5) = (xl(1,2)-2*al*(2*l(1)+1)*xl(1,0)&
          +4*al*al*x0(1)**(l(1)+2))*xl(2,0)*xl(3,0)*ex
       pri(ipri,6) = (xl(2,2)-2*al*(2*l(2)+1)*xl(2,0)&
          +4*al*al*x0(2)**(l(2)+2))*xl(3,0)*xl(1,0)*ex
       pri(ipri,7) = (xl(3,2)-2*al*(2*l(3)+1)*xl(3,0)&
          +4*al*al*x0(3)**(l(3)+2))*xl(1,0)*xl(2,0)*ex
       pri(ipri,8) = (xl(1,1)-2*al*x0(1)**(l(1)+1))*&
          (xl(2,1)-2*al*x0(2)**(l(2)+1))*xl(3,0)*ex
       pri(ipri,9) = (xl(1,1)-2*al*x0(1)**(l(1)+1))*&
          (xl(3,1)-2*al*x0(3)**(l(3)+1))*xl(2,0)*ex
       pri(ipri,10)= (xl(3,1)-2*al*x0(3)**(l(3)+1))*&
          (xl(2,1)-2*al*x0(2)**(l(2)+1))*xl(1,0)*ex
    end do

   end subroutine pri012

   subroutine index0(i,l)
     use tools_io
     use param

     integer, intent(in) :: i
     integer, intent(out) :: l(3)

     integer, parameter :: LI(3,35) = reshape((/&
        0,0,0, 1,0,0, 0,1,0, 0,0,1, 2,0,0, 0,2,0, 0,0,2, 1,1,0, 1,0,1, 0,1,1,&
        3,0,0, 0,3,0, 0,0,3, 2,1,0, 2,0,1, 0,2,1, 1,2,0, 1,0,2, 0,1,2, 1,1,1,& ! f primitives
        0,0,4, 0,1,3, 0,2,2, 0,3,1, 0,4,0, 1,0,3, 1,1,2, 1,2,1, 1,3,0, 2,0,2,& ! g primitives
        2,1,1, 2,2,0, 3,0,1, 3,1,0, 4,0,0/),shape(li))

     if (i < 1 .or. i > 35) then  
             write(*,*) 'itype',i 
             call error('index0','type not allowed',faterr) 
     end if        
     l = li(:,i)

   end subroutine index0

  !> Calculate the density and gradient using atomic densities
  subroutine promolecular_grid(x,m,rho,grad,hess,rhom,nfr,autofr,fronly)
    use reader
    use tools_io
    use param

    real*8, intent(in) :: x(3)
    type(molecule), intent(in) :: m
    integer, intent(in) :: nfr
    real*8, intent(out) :: rho, grad(3), hess(3,3), rhom(nfr)
    logical, intent(in) :: autofr, fronly

    integer :: i, j, k, iz, iq, ityp
    real*8 :: xd(3), r, r1, r2, f, fp, fpp, rfac, radd

    rho = 0d0
    rhom = 0d0
    grad = 0d0
    hess = 0d0
    do i = 1, m%n
       iz = m%z(i)
       iq = m%q(i)
       ityp = iztype(iz,iq)
       if (ityp == 0) call error('promolecular_grid','Atom type not initialized',faterr)
       if (.not.grd(ityp)%init) call error('promolecular_grid','Atom type not initialized',faterr)

       xd = x - m%x(:,i)
       r = sqrt(dot_product(xd,xd))
       r = max(r,1d-10)
       r1 = 1d0 / r
       r2 = r1 * r1
       call grid1_interp(grd(ityp),r,f,fp,fpp)

       if (.not.autofr.and.m%ifrag(i)>0) rhom(m%ifrag(i)) = rhom(m%ifrag(i)) + f
       if (.not.fronly) then
          rho = rho + f
          grad = grad + fp * xd * r1
          rfac = (fpp - fp * r1)
          do j = 1, 3
             hess(j,j) = hess(j,j) + fp * r1 + rfac * r2 * xd(j) * xd(j)
             do k = j+1, 3
                radd = rfac * r2 * xd(j) * xd(k)
                hess(j,k) = hess(j,k) + radd
             end do
          end do
       end if
    end do
    do j = 1, 3
       do k = j+1, 3
          hess(k,j) = hess(j,k)
       end do
    end do

  end subroutine promolecular_grid

  !> Interpolate the radial grid g at distance r0, and obtain the value,
  !> first derivative and second derivative.
  subroutine grid1_interp(g,r0,f,fp,fpp)
    use reader, only: grid1

    type(grid1), intent(in) :: g !< The radial grid.
    real*8, intent(in) :: r0 !< Value of the radial coordinate.
    real*8, intent(out) :: f !< Interpolated value
    real*8, intent(out) :: fp !< Interpolated first derivative
    real*8, intent(out) :: fpp !< Interpolated second derivative

    integer :: ir, i, j, ii
    real*8 :: r, prod, rr(4), dr1(4), x1dr12(4,4)

    f = 0d0
    fp = 0d0
    fpp = 0d0

    if (.not.g%init) return
    if (r0 >= g%rmax) return

    ! careful with grid limits.
    if (r0 <= g%r(1)) then
       ir = 1
       r = g%r(1)
    else
       ir = 1 + floor(log(r0/g%a)/g%b)
       r = r0
    end if

    x1dr12 = 0d0
    do i = 1, 4
       ii = min(max(ir,2),g%ngrid-2) - 2 + i
       rr(i) = g%r(ii)
       dr1(i) = r - rr(i)
       do j = 1, i-1
          x1dr12(i,j) = 1d0 / (rr(i) - rr(j))
          x1dr12(j,i) = -x1dr12(i,j)
       end do
    end do

    ! interpolate, lagrange 3rd order, 4 nodes
    do i = 1, 4
       ii = min(max(ir,2),g%ngrid-2) - 2 + i
       prod = 1.d0
       do j = 1 ,4
          if (i == j) cycle
          prod = prod * dr1(j) * x1dr12(i,j)
       end do
       f = f + g%f(ii) * prod
       fp = fp + g%fp(ii) * prod
       fpp = fpp + g%fpp(ii) * prod
    end do

  end subroutine grid1_interp

  subroutine checkden()
    use reader
    use tools_io
    use param

    type(molecule) :: m

    integer, parameter :: mpts = 1000
    real*8, parameter :: rmax = 8d0

    integer :: i, iz
    real*8 :: r, x(3), rho1, rho2, grad1(3), grad2(3), hess1(3,3), hess2(3,3),igmgg(3)
    real*8 :: rhom(100)

    integer, parameter :: lu = 10

    ! define a hydrogen atom
    m%name = "h.dummy"
    m%ifile = ifile_xyz
    m%n = 1
    allocate(m%x(3,1))
    allocate(m%z(1))
    allocate(m%q(1))
    m%x = 0d0
    m%z = 1
    m%q = 0

    open(unit=lu,file="checkden.dat")
    do iz = 1, 18
       m%z = iz
       write (lu,'("# Z = ",I3)') iz
       do i = 1, mpts
          r = rmax * real(i,8) / real(mpts,8)
          x = (/0d0,0d0,r/)
          call promolecular(x,m,rho1,grad1,hess1,rhom,100,.true.,igmgg)
          call promolecular_grid(x,m,rho2,grad2,hess2,rhom,100,.true.,.false.)
          write (lu,*) i, r, rho1, rho2
       end do
       write (lu,*)
       write (lu,*)
    end do
    close(lu)

    stop 1

  end subroutine checkden

  subroutine calcexc(ixc,exc,rho,grad,tau)
    use tools_io
    use param
    implicit none

    integer, intent(in) :: ixc(2)
    real*8, intent(out) :: exc
    real*8, intent(in) :: rho, grad(3), tau

    real*8, parameter :: pol_2_13  = 2d0**(1d0/3d0)
    real*8, parameter :: pol_2_23  = 2d0**(2d0/3d0)
    real*8, parameter :: cslater   = 3d0/2d0 * (3d0/4d0/pi)**(1d0/3d0)
    real*8, parameter :: cslater_u = cslater / 2**(1d0/3d0)
    real*8, parameter :: b_b86     = 0.00336d0
    real*8, parameter :: b_b86_u   = b_b86 / 2**(1d0/3d0)
    real*8, parameter :: g_b86     = 0.00449d0
    real*8, parameter :: b_b88     = 0.0042d0
    real*8, parameter :: b_b88_u   = b_b88 / 2**(1d0/3d0)
    real*8 :: ex, ec, dum1, dum2, ucor0
    real*8 :: rho43, rdg, rdg2, asnh, rho2, grad2, tau2

    exc = 0d0
    ! Exchange
    select case (ixc(1))
    case (1)
       ! LDA exchange, Slater with alpha = 4/3, total rho
       ex = -cslater_u * rho**(4d0/3d0)
    case (2)
       ! PBE exchange (similar to B86), total rho
       rho43 = rho**(4d0/3d0)
       rdg2 = dot_product(grad,grad) / rho43**2 * pol_2_23
       ex = (-cslater_u -b_b86_u * rdg2 / (1d0 + g_b86 * rdg2)) * rho43
    case (3)
       ! Becke88 exchange, total rho
       rho43 = rho**(4d0/3d0)
       rdg = sqrt(dot_product(grad,grad)) / rho43 * pol_2_13
       rdg2 = rdg * rdg
       asnh = log(rdg+sqrt(rdg2+1d0))
       ex = (-cslater_u -b_b88_u * rdg2 / (1d0 + 6d0*b_b88*rdg*asnh)) * rho43
    case (-1)
       ex = 0d0
    case (99)
       ! test
       ! rs = (3d0 / (4d0*pi*rho))**(1d0/3d0)
       !call pbex(rho,dot_product(grad,grad),1,ex,dum1,dum2)
       !ex = -cslater_u * rho**(4d0/3d0) + ex
       ex = 0
    case default
       call error('calcprops','unknown ixc option',faterr)
    end select

    ! Correlation
    select case (ixc(2))
    case (1)
       ! Perdew-Wang 88 parametrization of HEG eps_c
       call gascor(0.5d0*rho,0.5d0*rho,ec,dum1,dum2,ucor0)
    case (2)
       ! PBE correlation (similar to B86), total rho
       call perdc(0.5d0*rho,0.5d0*rho,sqrt(dot_product(grad,grad)),0d0,ec)
    case (3)
       ! Becke88 correlation, total rho
       rho2 = 0.5d0 * rho
       grad2 = 0.25d0 * dot_product(grad,grad)
       tau2 = 0.5d0 * tau
       call b88_corr(rho2,rho2,grad2,grad2,tau2,tau2,ec)
    case (-1)
       ec = 0d0
    case (99)
       ! test
       !call gascor(0.5d0*rho,0.5d0*rho,ec,dum1,dum2,ucor0)
       !call pbec(rho,dot_product(grad,grad),1,ucor0,dum1,dum2)
       !ec = ec + ucor0
       ec = 0
    case default
       call error('calcprops','unknown ixc option',faterr)
    end select
    exc = ex + ec

  end subroutine calcexc

  ! Erin Johnson private communication, Sep 14th 2011
  subroutine gascor(RHO1,RHO2,ECOR,VC1,VC2,UCOR)
    use param, only: pi, eps
    IMPLICIT REAL*8(A-H,O-Z)

    DATA PGAM,PFZZ/0.5198421D0,1.709921D0/
    THRD=1.D0/3.D0
    THRD4=4.D0*THRD
    !**** PW91 CORRELATION LSDA,
    PRS = (0.75D0/PI/(RHO1+RHO2+EPS))**THRD
    PZET = (RHO1-RHO2)/(RHO1+RHO2+EPS)/(1.D0+EPS)
    PF = ((1.D0+PZET)**THRD4+(1.D0-PZET)**THRD4-2.D0)/PGAM
    PA = 0.0310907D0
    PA1 = 0.21370D0
    PB1 = 7.5957D0
    PB2 = 3.5876D0
    PB3 = 1.6382D0
    PB4 = 0.49294D0
    PP = 1.00D0
    PP1 = PP + 1.D0
    PQ0 = -2.D0*PA*(1.D0+PA1*PRS)
    PRS12 = DSQRT(PRS)
    PRS32 = PRS12**3
    PRSP = PRS**PP
    PQ1 = 2.D0*PA*(PB1*PRS12+PB2*PRS+PB3*PRS32+PB4*PRS*PRSP)
    PQ2 = DLOG(1.D0+1.D0/PQ1)
    PEU = PQ0*PQ2
    PQ3 = PA*(PB1/PRS12+2.D0*PB2+3.D0*PB3*PRS12+2.D0*PB4*PP1*PRSP)
    PEURS = -2.D0*PA*PA1*PQ2-PQ0*PQ3/(PQ1**2+PQ1)
    PA = 0.01554535D0
    PA1 = 0.20548D0
    PB1 = 14.1189D0
    PB2 = 6.1977D0
    PB3 = 3.3662D0
    PB4 = 0.62517D0
    PP = 1.00D0
    PP1 = PP + 1.D0
    PQ0 = -2.D0*PA*(1.D0+PA1*PRS)
    PRSP = PRS**PP
    PQ1 = 2.D0*PA*(PB1*PRS12+PB2*PRS+PB3*PRS32+PB4*PRS*PRSP)
    PQ2 = DLOG(1.D0+1.D0/PQ1)
    PEP = PQ0*PQ2
    PQ3 = PA*(PB1/PRS12+2.D0*PB2+3.D0*PB3*PRS12+2.D0*PB4*PP1*PRSP)
    PEPRS = -2.D0*PA*PA1*PQ2-PQ0*PQ3/(PQ1**2+PQ1)
    PA = 0.0168869D0
    PA1 = 0.11125D0
    PB1 = 10.357D0
    PB2 = 3.6231D0
    PB3 = 0.88026D0
    PB4 = 0.49671D0
    PP = 1.00D0
    PP1 = PP + 1.D0
    PQ0 = -2.D0*PA*(1.D0+PA1*PRS)
    PRSP = PRS**PP
    PQ1 = 2.D0*PA*(PB1*PRS12+PB2*PRS+PB3*PRS32+PB4*PRS*PRSP)
    PQ2 = DLOG(1.D0+1.D0/PQ1)
    PALFM = PQ0*PQ2
    PQ3 = PA*(PB1/PRS12+2.D0*PB2+3.D0*PB3*PRS12+2.D0*PB4*PP1*PRSP)
    ALFRSM = -2.D0*PA*PA1*PQ2-PQ0*PQ3/(PQ1**2+PQ1)
    PALFC = -PALFM
    PZ4 = PZET**4
    PEC = PEU*(1.D0-PF*PZ4)+PEP*PF*PZ4-PALFM*PF*(1.D0-PZ4)/PFZZ
    ECOR = (RHO1+RHO2) * PEC
    PECRS = PEURS*(1.D0-PF*PZ4)+PEPRS*PF*PZ4-ALFRSM*PF*(1.D0-PZ4)/PFZZ
    PFZ = THRD4*((1.D0+PZET)**THRD-(1.D0-PZET)**THRD)/PGAM
    PECZET = 4.D0*(PZET**3)*PF*(PEP-PEU+PALFM/PFZZ)+PFZ*(PZ4*PEP-PZ4*PEU-(1.D0-PZ4)*PALFM/PFZZ)
    PCOMM = PEC -PRS*PECRS/3.D0-PZET*PECZET
    VC1 = PCOMM + PECZET
    VC2 = PCOMM - PECZET
    !**** END PW91 CORRELATION LSDA.
    PTC = -4.D0*PEC + 1.5D0*((1.D0+PZET)* VC1+(1.D0-PZET)* VC2)
    PTCOR = (RHO1+RHO2) * PTC
    UCOR =   ECOR-PTCOR
    !**** END PW91 CORRELATION KINETIC AND POTENTIAL ENERGIES.

  end subroutine gascor

  ! Erin Johnson private communication, Sep 14th 2011
  subroutine perdc(RHO1,RHO2,DROT,ZET,CPBE)
    use param, only: pi
    IMPLICIT REAL*8(A-H,O-Z)

    ROT=RHO1+RHO2
    CALL GASCOR(RHO1,RHO2,ECOR,VC1,VC2,UCOR)
    THRD=1.D0/3.D0
    THRD2=2.D0/3.D0
    FERMIK=(3.D0*PI*PI)**(1.D0/3.D0)
    PHI=0.5D0*((1.D0+ZET)**THRD2+(1.D0-ZET)**THRD2)
    PKF=FERMIK*ROT**THRD
    PKS=DSQRT(4.D0*PKF/PI)
    T=DROT/(2.D0*PHI*PKS*ROT)
    A=2.14612D0/(DEXP(-ECOR/ROT/0.031091D0/PHI**3)-1.D0)
    TSTUFF=(1.D0+A*T**2)/(1.D0+A*T**2+A**2*T**4)
    H=0.031091D0*PHI**3*DLOG(1.D0+2.14612D0*T**2*TSTUFF)
    CPBE=ECOR+ROT*H

  end subroutine perdc

  ! Erin Johnson private communication, Sep 14th 2011
  subroutine b88_corr(rho1,rho2,grad1,grad2,tau1,tau2,ec)
    use param, only: pi
    implicit none

    real*8, intent(in) :: rho1, rho2 !< rho, 2 spins
    real*8, intent(in) :: grad1, grad2 !< nablarho * nablarho, 2 spins
    real*8, intent(in) :: tau1, tau2 !< kinetic energy density, 2 spins
    real*8, intent(out) :: ec

    real*8, parameter :: thrd4 = 4d0/3d0
    real*8, parameter :: cex = 1.5D0*(0.75D0/pi)**(1.D0/3.D0)

    real*8 :: drhosq, drho, weizs1, weizs2, rho43, dlss, xl, asnh
    real*8 :: x88, rfg1, rfg2, z11, z22, z12, c8811, c8812, c8822
    real*8 :: dsigs1, dsigs2

    ! SPIN 1:
    weizs1 = 0.25d0 * grad1 / rho1
    drhosq = 4.D0 * rho1 * weizs1
    drho = sqrt(drhosq)
    rho43 = rho1**thrd4
    dlss = drho / rho43
    xl = -cex * rho43
    asnh = log(dlss+sqrt(dlss**2+1.D0))
    x88 = xl - 0.0042D0*drho*dlss/(1.D0+6.D0*0.0042D0*dlss*asnh)
    rfg1 = -0.5D0*rho1/x88
    z11 = 0.88D0*2.D0*rfg1
    dsigs1 = 2d0 * tau1 - weizs1
    c8811 = -0.01D0*rho1*dsigs1*z11**4*(1.D0-2.D0/z11*log(1.D0+0.5D0*z11))

    ! SPIN 2:
    weizs2 = 0.25d0 * grad2 / rho2
    drhosq = 4.D0 * rho2 * weizs2
    drho = sqrt(drhosq)
    rho43 = rho2**thrd4
    dlss = drho/rho43
    xl = -cex*rho43
    asnh = log(dlss+dsqrt(dlss**2+1.d0))
    x88 = xl - 0.0042d0*drho*dlss/(1.d0+6.d0*0.0042d0*dlss*asnh)
    rfg2 = -0.5d0*rho2/x88
    z22 = 0.88d0*2.d0*rfg2
    dsigs2 = 2d0 * tau2 - weizs2
    c8822 = -0.01d0*rho2*dsigs2*z22**4*(1.d0-2.d0/z22*log(1.d0+0.5d0*z22))

    ! opp spins:
    z12=0.63d0*(rfg1+rfg2)
    c8812=-0.8d0*rho1*rho2*z12**2*(1.d0-log(1.d0+z12)/z12)

    ec = c8811 + c8822 + c8812

  end subroutine b88_corr

 ! Computing splines for many files
   
  subroutine spline(rhomref,dimgradmref,xspl,yspl,ye2)
    use reader
    use tools_io
    use tools_math
    use param

    implicit none
    real*8, intent(in) :: rhomref,dimgradmref
    real*8, intent(inout) :: xspl(nspl),yspl(nspl)
    real*8, intent(inout) :: ye2(nspl) 
    real*8 :: low,high
    integer :: i


          


    do i=1,nspl
       low=xspl(i)-0.0005D0
       high=xspl(i)+0.0005D0 
       if ((rhomref.gt.low).and.(rhomref.lt.high)) then

          if (dimgradmref.lt.yspl(i)) yspl(i)=dimgradmref
       endif 
    enddo  

   end subroutine spline     

 !splines


subroutine DoSpline(nrefs,mref,ref,xinc,nfrag,xspl,yspl,ye2,dy,yy1,ypn,autofrag,doelf,ixc) 

  use reader
  use tools_io
  use tools_math
  use param

  implicit none
  integer, intent(in):: nrefs
  type(molecule),intent(inout) :: mref(nrefs)
  type(reference),intent(inout) :: ref(nrefs)
  real*8, intent(in) :: xinc(3)
  integer,intent(in) :: nfrag
  real*8,dimension(nspl),intent(inout) :: xspl,yspl,ye2,dy
  real*8, intent(inout) :: yy1,ypn
  logical, intent(in) :: autofrag, doelf
  integer, intent(in) :: ixc(2)

  real*8 :: rho, grad(3), dimgrad, grad2, hess(3,3), elf, exc
  real*8 :: wk1(3), wk2(3), heigs(3), hvecs(3,3)

  integer :: iref,k,j,i,istat
  real*8  :: x(3) 
  


!$omp parallel do private (iref,x,rho,grad,hess,heigs,hvecs,wk1,wk2,istat,grad2,&
!$omp dimgrad) schedule(dynamic)
do iref=1,nrefs ! Running over the number of references
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
         end do !i=0,ref(iref)%nstep(1)-1 
      end do ! j=0,ref(iref)%nstep(2)-1
   end do  ! k=0,ref(iref)%nstep(3)-1 
enddo !nrefs
!$omp end parallel do 

call nrspline(xspl,yspl,nspl,dy(1),ypn,ye2)     

end subroutine Dospline

end module props
