program main

  use mtmod

  implicit none

  integer, parameter :: MAXFRAME = 2000000000
  character (len=255) :: infile, outfile, ndxfile, pdbfile, line
  character (len=255), allocatable :: xtcfile(:)
  integer, allocatable :: indxrep(:,:), idstep(:), indx(:)
  integer, allocatable :: sd(:), step(:)
  integer :: narg, ndihedrange, i, j, k, l, ii, jj
  integer :: imc, ndim, idrep, nreplica, ncomplex
  integer :: molid(2), hassign(4), nneighbor, numat, natom, magic, ierr, neq, ntot
  integer :: nrow, idr, idn, idi, NMC, nmcwrt, nblock, ntot_block, nn1, nn2
  character (len=80), allocatable :: ARGM(:)
  integer :: nrestype, rngseed, numconf
  character(len=3), dimension(:), allocatable :: restype, restypetmp
  integer, dimension(:), allocatable :: idres
  integer, dimension(:,:), allocatable :: Nw, Nv
  character(len=3) :: aa
  real(8), dimension(:,:), allocatable :: lnw, lnv
  real(8), dimension(:), allocatable :: wave, vave, wstd, vstd
  integer, dimension(:), allocatable :: idconf
  real(4), allocatable :: xt(:,:)
  real(8), allocatable :: x0(:,:)
  real(8) :: boxsize(3), cm(3), dr(3), r
  real(4) :: time, prec, box(9)
  real(8), allocatable :: HelixProp(:)
  real(8) :: Dihedral, r14, lnL, lnL_new, xx
  real(8), allocatable :: dihedrange(:,:)
  real(8) :: r14range(2)
  integer, allocatable :: H_Dihedral(:)
  logical, allocatable :: HelixPart(:)
  real(8), dimension(:), allocatable :: hfrac
  real(8), dimension(:,:), allocatable :: hfracsim
  real(8), parameter :: dv = 0.1d0
  real(8), parameter :: dw = 0.1d0
  real(8), parameter :: lwmax = 10.d0
  real(8), parameter :: lvmax = 10.d0
  real(8), parameter :: lwmin = -10.d0
  real(8), parameter :: lvmin = -10.d0
  real(8) :: kSpr 
  integer, parameter :: Ncol = 10
  real(8), dimension(:,:), allocatable :: lnwcol, lnvcol
  logical :: helical

  narg = iargc()
  if (narg < 1) then
     call print_usage()
     stop
  else
     allocate(ARGM(narg))
     do i = 1, narg
        call getarg(i,ARGM(i))
     end do
     infile = ''
     outfile = ''
     ndxfile = ''
     pdbfile = ''
     idrep = 0
     molid = 0
     hassign = 0
     nneighbor = 0
     ndihedrange = 0
     r14range(1) = 0.0
     r14range(2) = 100.0
     neq = 0
     rngseed = 1234567
     NMC = 100000
     nblock = 1 ! default
     do i = 1, narg-1
        select case (trim(ARGM(i)))
        case ('-x','-xtc')
           infile = trim(ARGM(i+1))
        case ('-o','-out')
           outfile = trim(ARGM(i+1))
        case ('-n','-ndx')
           ndxfile = trim(ARGM(i+1))
        case ('-p','-pdb')
           pdbfile = trim(ARGM(i+1))
        case ('-assign')
           do j = 1, 4
              if (ARGM(i+1)(j:j) == '1') hassign(j) = 1
           end do
        case ('-dihed','-dihedral')
           ndihedrange = 0
           do j = i+1, narg
              read(ARGM(j), *, iostat=ierr) r
              if (ierr /= 0) exit
              ndihedrange = ndihedrange + 1
           end do
           if (ndihedrange == 0) stop 'Invalid entry in -dihed option!'
           if (mod(ndihedrange,2) /= 0) stop 'Invalid entry in -dihed option!'
           ndihedrange = ndihedrange/2

           allocate(dihedrange(2,ndihedrange))
           do j = 1, ndihedrange 
              read(ARGM(i+2*(j-1)+1),*,iostat=ierr) dihedrange(1,j)
              if (ierr /= 0) stop 'Invalid entry in -t option!'
              read(ARGM(i+2*j),*,iostat=ierr) dihedrange(2,j)
              if (ierr /= 0) stop 'Invalid entry in -t option!'
           end do
        case ('-r14')
           read(ARGM(i+1),*,iostat=ierr) r14range(1)
           if (ierr /= 0) stop 'Invalid entry in -r14 option!'
           read(ARGM(i+2),*,iostat=ierr) r14range(2)
           if (ierr /= 0) stop 'Invalid entry in -r14 option!'
        case ('-r','-rep')
           read(ARGM(i+1),*,iostat=ierr) idrep
           if (ierr /= 0) stop 'Invalid entry in -r option!'
        case ('-m','-mol')
           read(ARGM(i+1),*,iostat=ierr) molid(1)
           if (ierr /= 0) stop 'Invalid entry in -m option!'
           read(ARGM(i+2),*,iostat=ierr) molid(2)
        case ('-e','-eq')
           read(ARGM(i+1),*,iostat=ierr) neq
           if (ierr /= 0) stop 'Invalid entry in -e option!'
        case ('-nn','-nneighbor')
           read(ARGM(i+1),*,iostat=ierr) nneighbor
           if (ierr /= 0) stop 'Invalid entry in -nn option!'
        case ('-s','-seed')
           read(ARGM(i+1),*,iostat=ierr) rngseed
           if (ierr /= 0) stop 'Invalid entry in -s option!'
        case ('-nmc')
           read(ARGM(i+1),*,iostat=ierr) NMC
           if (ierr /= 0) stop 'Invalid entry in -nmc option!'
        case ('-block')
           read(ARGM(i+1),*,iostat=ierr) nblock
           if (ierr /= 0) stop 'Invalid entry in -block option!'
        case ('-h','-help')
           call print_usage()
           stop
        end select
     end do
  end if

  if (rngseed < 1) stop 'Wrong RNG seed!'
  if (nblock < 1) nblock = 1

! determine the number of replicas
  if (len_trim(ndxfile) > 0) then
     nreplica = NumberColumn(ndxfile) - 1
     if (nreplica > 1) then
        nrow = NumberRow(ndxfile)
        allocate(indxrep(nreplica,nrow))
        allocate(idstep(nrow))
        idstep = 0
        indxrep = 0
        open(10, file=trim(ndxfile), status='old', iostat=ierr)
        if (ierr /= 0) stop 'Cannot open the ndx file!'
        do i = 1, nrow
           read(10,*) idstep(i), indxrep(:,i)
           if (minval(indxrep(:,i)) < 1) indxrep(:,i) = indxrep(:,i) + 1
        end do
        close(10)
     end if
  else
     nreplica = 1
  end if
  allocate(xtcfile(nreplica))
  allocate(sd(nreplica))

  if (nreplica > 1) then
     open(10, file=trim(infile), status='old', iostat=ierr)
     if (ierr /= 0) stop 'Cannot open the infile!'
     do i = 1, nreplica
        read(10,'(a)',iostat=ierr) xtcfile(i)
        if (ierr /= 0) stop 'Invalid entry in the infile!'
     end do
     close(10)
  else
     i = len_trim(infile)
     if (infile(i-3:i) == '.xtc') then
        xtcfile(1) = trim(infile)
     else
        open(10, file=trim(infile), status='old', iostat=ierr)
        read(10,'(a)',iostat=ierr) xtcfile(1)
        if (ierr /= 0) stop 'Invalid entry in the infile!'
        close(10)
     end if
  end if

! Read the trajectory
  allocate(step(nreplica))
  allocate(indx(nreplica))
  call xdrfopen(sd(1), trim(xtcfile(1)), "r", ierr)
  if (ierr /= 1) then
     stop 'Cannot open the xtc file!'
  end if
  call xtcheader(sd(1), magic, natom, step(1), time, ierr)
  call xdrfclose(sd(1), ierr)
  allocate(xt(3*natom,nreplica))
  allocate(x0(3,natom))
  write(0,*) 'Number of atoms = ', natom

  if (all(molid > 0)) then
     numat = molid(2) - molid(1) + 1
  else
     numat = natom
     molid(1) = 1
     molid(2) = natom
  end if
  if (numat < 1) stop 'Invalid entry in -m option!'

! determine the number of residue types
  allocate(idres(numat))
  open(10, file=trim(pdbfile), status='old', iostat=ierr)
  if (ierr /= 0) stop 'Cannot open the pdb file!'
  nrestype = 0
  allocate(restypetmp(100))
  k = 0
  l = 0
  idres = 0
  restypetmp = ''
  do while (ierr == 0)
     read(10, '(a)', iostat=ierr) line
     if (ierr == 0 .and. line(1:6) == 'ATOM  ' .and. line(14:15) == 'CA') then
        k = k + 1
        if (k >= molid(1) .and. k <= molid(2)) then
           l = l + 1
           if (l > 0 .and. l < (numat + 1)) then
              aa = line(18:20)
              j = 0
              do i = 1, nrestype
                 if (aa == restypetmp(i)) j = i
              end do
              if (j == 0) then
                 nrestype = nrestype + 1
                 restypetmp(nrestype) = aa
                 idres(l) = nrestype
              else
                 idres(l) = j
              end if
           end if
        end if
     end if
  end do
  close(10)

  if (nrestype < 1) stop 'No residue type!'
  allocate(restype(nrestype))
  restype = restypetmp(1:nrestype)
  deallocate(restypetmp)

! determine the number of complexes
  ncomplex = 0
  call xdrfopen(sd(1), trim(xtcfile(1)), "r", ierr)
  do imc = 1, MAXFRAME
     call readxtc(sd(1), natom, step(1), time, box, xt(:,1), prec, ierr)
     if (ierr /= 1) exit
     if (step(1) > neq) ncomplex = ncomplex + 1
  end do
  call xdrfclose(sd(1), ierr)
  write(0,*) 'Number of complexes = ', ncomplex

  ntot_block = ncomplex/nblock
  write(0,*) 'Number of complexes for each block = ', ntot_block

  do i = 1, nreplica
     call xdrfopen(sd(i), trim(xtcfile(i)), "r", ierr)
     if (ierr /= 1) stop 'Cannot open the xtc file!'
  end do

  allocate(H_Dihedral(numat))
  allocate(HelixPart(numat))
  allocate(idconf(numat))
  allocate(hfracsim(numat,nblock+1))

  numconf = nrestype
  idconf = idres

  allocate(Nw(numconf,nblock+1))
  allocate(Nv(numconf,nblock+1))
  Nw = 0
  Nv = 0
  hfracsim = 0.d0

  ntot = 0 
  idi = 1
  outer: do imc = 1, MAXFRAME
     do i = 1, nreplica
        call readxtc(sd(i), natom, step(i), time, box, xt(:,i), prec, ierr)
        if (ierr /= 1) exit outer
     end do
     if (any(step /= step(1))) stop 'Steps mismatch!'
     if (step(1) > neq) then
        if (nreplica > 1) then
           idn = 0
           do i = idi, nrow
              if (idstep(i) == step(1)) then
                 idn = i
                 exit
              end if
           end do
           if (idn == 0) exit outer
           idi = idn 
           indx = indxrep(:,idn)
           if (minval(indx) < 1) indx = indx + 1
           idr = 0
           do j = 1, nreplica
              if (indx(j) == idrep) idr = j
           end do
           if (idr == 0) stop 'No replica!'
        else
           idr = 1
        end if
        boxsize(1) = dble(box(1)*10.0)
        boxsize(2) = dble(box(5)*10.0)
        boxsize(3) = dble(box(9)*10.0)
        do i = 1, natom
           cm = dble(xt(3*i-2:3*i,idr)*10.0)
           if (i == 1) then
              x0(:,i) = cm
           else
              cm = cm - x0(:,i-1)
              do j = 1, 3
                 cm(j) = cm(j) - boxsize(j)*ANINT(cm(j)/boxsize(j))
              end do
              x0(:,i) = x0(:,i-1) + cm
           end if
        end do

        H_Dihedral = 0
        do i = 1, natom-3
           Dihedral = ComputeDihedral(x0(:,i),x0(:,i+1),x0(:,i+2),x0(:,i+3))
           r14 = sqrt(dot_product(x0(:,i)-x0(:,i+3),x0(:,i)-x0(:,i+3)))
           helical = .TRUE.
           do j = 1, ndihedrange
              if (Dihedral < dihedrange(1,j) .or. Dihedral > dihedrange(2,j)) helical = .FALSE.
           end do
           if (r14 < r14range(1) .or. r14 > r14range(2)) helical = .FALSE.
           if (helical) then
              do j = 1, 4
                 if (hassign(j) == 1) H_Dihedral(i+j-1) = 1
              end do
           end if
        end do

        HelixPart = .TRUE.
        do i = 1, numat
           do j = i-nneighbor, i+nneighbor
              if (j >= 1 .and. j <= numat) then
                 if (H_Dihedral(j) == 0) HelixPart(i) = .FALSE.
              end if
           end do
        end do

        ntot = ntot + 1
        j = (ntot-1)/ntot_block + 1

!        do i = 2, numat-1
!           if (HelixPart(i) .and. HelixPart(i-1) .and. HelixPart(i+1)) then
!              hfracsim(i,j) = hfracsim(i,j) + 1.d0
!           end if
!        end do
        do i = 1, numat
           if (HelixPart(i)) hfracsim(i,j) = hfracsim(i,j) + 1.d0
        end do

        if (j <= nblock) then
           do i = 1, numat
              if (HelixPart(i)) then
                 if (i == 1 .or. i == numat) then  !  end residues; weight = v
                    Nv(idconf(i),j) = Nv(idconf(i),j) + 1
                 else
                    if (Helixpart(i-1) .and. HelixPart(i+1)) then ! helix segment; weight = w
                       Nw(idconf(i),j) = Nw(idconf(i),j) + 1
                    else
                       Nv(idconf(i),j) = Nv(idconf(i),j) + 1
                    end if
                 end if
              end if
           end do
        end if
     end if
  end do outer
  Nw(:,nblock+1) = sum(Nw(:,1:nblock),dim=2)
  Nv(:,nblock+1) = sum(Nv(:,1:nblock),dim=2)
  hfracsim(:,1:nblock) = hfracsim(:,1:nblock)/dble(ntot_block)
  hfracsim(:,nblock+1) = sum(hfracsim(:,1:nblock),dim=2)/dble(nblock)

  do i = 1, nreplica
     call xdrfclose(sd(i), ierr)
  end do

  call sgrnd(rngseed)

! Maximum likelihood estimation of w and v
  allocate(lnw(numconf,nblock+1))
  allocate(lnv(numconf,nblock+1))
  allocate(lnwcol(numconf,Ncol))
  allocate(lnvcol(numconf,Ncol))
  allocate(hfrac(numat))
  allocate(wave(numconf))
  allocate(vave(numconf))
  nmcwrt = NMC/10

  do j = 1, nblock+1
     write(0,*) "Block = ", j
     if (j <= nblock) then
        ncomplex = ntot_block
     else
        ncomplex = ntot_block*nblock
     end if
     kSpr = dble(ncomplex)*10.d0

     do ii = 1, Ncol

        do i = 1, numconf
           lnwcol(i,ii) = 2.d0*(grnd() - 0.5d0)
           lnvcol(i,ii) = 2.d0*(grnd() - 0.5d0)
        end do

        lnL = LogLikelihood(ncomplex,numat,numconf,Nv(:,j),Nw(:,j), &
           lnvcol(:,ii),lnwcol(:,ii),idconf)

        do imc = 1, NMC

           i = int(grnd()*numconf) + 1
           i = min(i,numconf)
           xx = lnwcol(i,ii)
           lnwcol(i,ii) = lnwcol(i,ii) + dw*(grnd() - 0.5d0)
           if (lnwcol(i,ii) > lwmax) lnwcol(i,ii) = 2.d0*lwmax - lnwcol(i,ii)
           if (lnwcol(i,ii) < lwmin) lnwcol(i,ii) = 2.d0*lwmin - lnwcol(i,ii)
           lnL_new = LogLikelihood(ncomplex,numat,numconf,Nv(:,j),Nw(:,j), &
               lnvcol(:,ii),lnwcol(:,ii),idconf)
           if (lnL_new > lnL) then
              lnL = lnL_new
           else if (grnd() < exp(lnL_new - lnL)) then
              lnL = lnL_new
           else
              lnwcol(i,ii) = xx
           end if
    
           i = int(grnd()*numconf) + 1
           i = min(i,numconf)
           xx = lnvcol(i,ii)
           lnvcol(i,ii) = lnvcol(i,ii) + dv*(grnd() - 0.5d0)
           if (lnvcol(i,ii) > lvmax) lnvcol(i,ii) = 2.d0*lvmax - lnvcol(i,ii)
           if (lnvcol(i,ii) < lvmin) lnvcol(i,ii) = 2.d0*lvmin - lnvcol(i,ii)
           lnL_new = LogLikelihood(ncomplex,numat,numconf,Nv(:,j),Nw(:,j), &
              lnvcol(:,ii),lnwcol(:,ii),idconf)
           if (lnL_new > lnL) then
              lnL = lnL_new
           else if (grnd() < exp(lnL_new - lnL)) then
              lnL = lnL_new
           else
              lnvcol(i,ii) = xx
           end if

!           if (mod(imc,nmcwrt) == 0)  write(0,*) imc, lnL, lnwcol(:,ii)
        end do

        do imc = 1, NMC

           i = int(grnd()*(numconf)) + 1
           i = min(i,numconf)
           xx = lnwcol(i,ii)
           lnwcol(i,ii) = lnwcol(i,ii) + 0.2*dw*(grnd() - 0.5d0)
           if (lnwcol(i,ii) > lwmax) lnwcol(i,ii) = 2.d0*lwmax - lnwcol(i,ii)
           if (lnwcol(i,ii) < lwmin) lnwcol(i,ii) = 2.d0*lwmin - lnwcol(i,ii)
           lnL_new = LogLikelihood(ncomplex,numat,numconf,Nv(:,j),Nw(:,j), &
              lnvcol(:,ii),lnwcol(:,ii),idconf)
           if (lnL_new > lnL) then
              lnL = lnL_new
           else
              lnwcol(i,ii) = xx
           end if
    
           i = int(grnd()*(numconf)) + 1
           i = min(i,numconf)
           xx = lnvcol(i,ii)
           lnvcol(i,ii) = lnvcol(i,ii) + 0.2*dv*(grnd() - 0.5d0)
           if (lnvcol(i,ii) > lvmax) lnvcol(i,ii) = 2.d0*lvmax - lnvcol(i,ii)
           if (lnvcol(i,ii) < lvmin) lnvcol(i,ii) = 2.d0*lvmin - lnvcol(i,ii)
           lnL_new = LogLikelihood(ncomplex,numat,numconf,Nv(:,j),Nw(:,j), &
              lnvcol(:,ii),lnwcol(:,ii),idconf)
           if (lnL_new > lnL) then
              lnL = lnL_new
           else
              lnvcol(i,ii) = xx
           end if

!           if (mod(imc,nmcwrt) == 0)  write(0,*) imc, lnL, lnwcol(:,ii)
        end do
      end do

      do i = 1, numconf
         lnw(i,j) = sum(lnwcol(i,:))/dble(Ncol)
         lnv(i,j) = sum(lnvcol(i,:))/dble(Ncol)
      end do
  end do

  allocate(wstd(numconf))
  allocate(vstd(numconf))

  wave = exp(lnw(:,nblock+1))
  vave = exp(lnv(:,nblock+1))
  wstd = 0.0
  vstd = 0.0
  do i = 1, nblock
     wstd = wstd + (exp(lnw(:,i)) - wave)**2
     vstd = vstd + (exp(lnv(:,i)) - vave)**2
  end do
  if (nblock > 1) then
     wstd = sqrt(wstd/dble(nblock-1))
     vstd = sqrt(vstd/dble(nblock-1))
  end if

!  compute the helical fraction for each residue
  do i = 1, numat
     hfrac(i) = HelixFraction(numat,numconf,wave,vave,idconf,i)
  end do
 
  open(10, file=trim(outfile), status='unknown', iostat=ierr)
  if (ierr /= 0) stop 'Cannot open the outfile!'
  do i = 1, numat
     write(10,'(i6,1x,a,1x,4ES16.6,2f12.6)') i, restype(idres(i)), wave(idconf(i)), &
        wstd(idconf(i)), vave(idconf(i)), vstd(idconf(i)), hfrac(i), hfracsim(i,nblock+1)
  end do
  close(10)

contains

  subroutine print_usage()

    implicit none

    write (0, '(a)') "Usage: qvalue OPTIONS"
    write (0, '(a)') "       OPTIONS:::"
    write (0, '(a)') "       -x(-xtc): xtc file"
    write (0, '(a)') "       -p(-pdb): pdb file"
    write (0, '(a)') "       -o(-out): output file"
    write (0, '(a)') "       -dihed(-dihedral): dihedral angle range"
    write (0, '(a)') "       -r14: r14 distance range"
    write (0, '(a)') "       -assign: h' assignment rule (0100, 0110, etc.)"
    write (0, '(a)') "       -n(-ndx): ndx file"
    write (0, '(a)') "       -r(-rep): replica"
    write (0, '(a)') "       -m(-mol): molecule id (start end)"
    write (0, '(a)') "       -e(-eq): equilibration step"
    write (0, '(a)') "       -s(-seed): RNG seed"
    write (0, '(a)') "       -nmc: number of MC steps"
    write (0, '(a)') "       -block: number of blocks"
    write (0, '(a)') "       -h(-help): help"

  end subroutine print_usage

  function HelixFraction(n, m, w, v, id, i) result(h)

    implicit none

    integer, intent(in) :: n, m, i
    real(8), dimension(m), intent(in) :: w, v
    integer, dimension(n), intent(in) :: id
    real(8) :: h
    real(8), dimension(3,3) :: Mtot, M1, M2
    real(8), dimension(1,3) :: x1
    real(8), dimension(3,1) :: x2
    real(8), dimension(1,1) :: xx, yy
    integer :: k, n1, n2

    h = 0.d0
    n1 = 1
    n2 = n
    if (i < n1 .or. i > n2) return

    Mtot = 0.d0
    Mtot(1,1) = 1.d0
    Mtot(2,2) = 1.d0
    Mtot(3,3) = 1.d0
    do k = n1, n2
       M1 = 0.d0
       M1(1,1) = w(id(k))
       M1(1,2) = v(id(k))
       M1(2,3) = 1.d0
       M1(3,1) = v(id(k))
       M1(3,2) = v(id(k))
       M1(3,3) = 1.d0
       M2 = matmul(Mtot,M1)
       Mtot = M2
    end do

    x1 = 0.d0
    x1(1,3) = 1.d0
    x2 = 1.d0
    x2(1,1) = 0.d0

    xx = matmul(x1,matmul(Mtot,x2))

    Mtot = 0.d0
    Mtot(1,1) = 1.d0
    Mtot(2,2) = 1.d0
    Mtot(3,3) = 1.d0
    do k = n1, n2
       M1 = 0.d0
       if (k == i) then
          M1(1,1) = w(id(k))
       else
          M1(1,1) = w(id(k))
          M1(1,2) = v(id(k))
          M1(2,3) = 1.d0
          M1(3,1) = v(id(k))
          M1(3,2) = v(id(k))
          M1(3,3) = 1.d0
       end if
       M2 = matmul(Mtot,M1)
       Mtot = M2
    end do
    
    x1 = 0.d0
    x1(1,3) = 1.d0
    x2 = 1.d0
    x2(1,1) = 0.d0

    yy = matmul(x1,matmul(Mtot,x2))

    h = yy(1,1)/xx(1,1)
 
  end function HelixFraction

  function LogLikelihood(nc, n, m, nv, nw, lv, lw, id) result(y)

    implicit none

    integer, intent(in) :: nc, n, m
    integer, dimension(m), intent(in) :: nv, nw
    real(8), dimension(m), intent(in) :: lv, lw
    integer, dimension(n), intent(in) :: id
    real(8) :: y
    integer :: i, n1, n2
    real(8), dimension(1,3) :: x1
    real(8), dimension(3,1) :: x2
    real(8), dimension(3,3) :: M1, M2, Mtot
    real(8) :: x(1,1), lnZ
    real(8), dimension(m) :: v, w

    v = exp(lv)
    w = exp(lw)

    n1 = 1
    n2 = n

    Mtot = 0.d0
    Mtot(1,1) = 1.d0
    Mtot(2,2) = 1.d0
    Mtot(3,3) = 1.d0
    do i = n1, n2
       M1 = 0.d0
       M1(1,1) = w(id(i))
       M1(1,2) = v(id(i))
       M1(2,3) = 1.d0
       M1(3,1) = v(id(i))
       M1(3,2) = v(id(i))
       M1(3,3) = 1.d0
       M2 = matmul(Mtot,M1)
       Mtot = M2
    end do

    x1 = 0.d0
    x1(1,3) = 1.d0
    x2 = 1.d0
    x2(1,1) = 0.d0
    x = matmul(x1,matmul(Mtot,x2))
    lnZ = log(x(1,1))

    y = - dble(nc)*lnZ
    do i = 1, m
       y = y + dble(nv(i))*lv(i) + dble(nw(i))*lw(i)
    end do

  end function LogLikelihood

  function ComputeDihedral(x1, x2, x3, x4) result(d)

    implicit none

    real(8), intent(in) :: x1(3), x2(3), x3(3), x4(3)
    real(8) :: d
    real(8) :: b1(3), b2(3), b3(3), n1(3), n2(3), m(3)
    real(8) :: r, x, y

    b1 = x2 - x1
    b2 = x2 - x3
    b3 = x4 - x3

    call CrossProduct(b1,b2,n1)
    call CrossProduct(b2,b3,n2)
    call CrossProduct(n1,b2,m)

    r = sqrt(dot_product(n1,n1))
    n1 = n1/r
    r = sqrt(dot_product(n2,n2))
    n2 = n2/r
    r = sqrt(dot_product(m,m))
    m = m/r

    x = dot_product(n1,n2)
    y = dot_product(m,n2)

    d = atan2(y,x)

  end function ComputeDihedral

  subroutine CrossProduct(a, b, c)

    implicit none

    real(8), intent(in) :: a(3), b(3)
    real(8), intent(out) :: c(3)

    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)

  end subroutine CrossProduct

  integer function NumberColumn(filename)

    implicit none

    character (len=*), intent(in) :: filename
    character (len=255) :: line, string
    integer :: i

    NumberColumn = 0
    open(10, file=trim(filename), status='old', iostat=ierr)
    if (ierr /= 0) stop 'Cannot open the file!'
    read(10, '(a)') line
    i = 1
    call read_lineblock(i, line, string)
    do while (len_trim(string) > 0)
       NumberColumn = NumberColumn + 1
       call read_lineblock(i, line, string)
    end do
    close(10)

  end function NumberColumn

  integer function NumberRow(filename)

    implicit none

    character (len=*), intent(in) :: filename
    character (len=255) :: line
    integer :: ierr

    open(10, file=trim(filename), status='old', iostat=ierr)
    if (ierr /= 0) stop 'Cannot open the file!'
    NumberRow = 0
    do while (ierr == 0)
       read(10, '(a)', iostat=ierr) line
       if (ierr == 0 .and. len_trim(line) > 0) NumberRow = NumberRow + 1
    end do
    close(10)

  end function NumberRow

  subroutine read_lineblock(i, line, bl)

    implicit none

    integer, intent(inout) :: i
    character (len=*), intent(in) :: line
    character (len=*), intent(out) :: bl
    integer :: i1, i2

    bl = ''

    i1 = i
    do while (line(i1:i1) == ' ' .and. i1 <= len_trim(line))
       i1 = i1 + 1
    end do
    i2 = i1
    do while (line(i2:i2) /= ' ' .and. i2 <= len_trim(line))
       i2 = i2 + 1
    end do
    read (line(i1:i2), '(a)') bl
    i = i2
    
  end subroutine read_lineblock

end program main
