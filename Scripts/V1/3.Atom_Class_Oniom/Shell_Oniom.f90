  implicit none

  character datafile*24, outfile*15, outfile1*31, outfile2*31, outfile3*31, outfile4*31,  outfile5*31, outfile6*31
  character title*50, dummy, molec_name*5, atom_type*1, num*3, So_name, Sv_name
  character Elem_So(20)*2, Type_So(20)*2, q1_So(20)*9, q2_So(20)*9, q3_So(20)*9
  character Elem_Sv(20)*2, Type_Sv(20)*2, q1_Sv(20)*9, q2_Sv(20)*9, q3_Sv(20)*9
  character Head1_So(20)*17, Head2_So(20)*17, Head3_So(20)*17, Head1_Sv(20)*17, Head2_Sv(20)*17, Head3_Sv(20)*17, intra(100)*50
  real Sigma_So(20), Eps_So(20), Sigma_Sv(20), Eps_Sv(20)
  integer bonds_So(20), link_So(20,20), bonds_Sv(20), link_Sv(20,20), is_dummy(20), Ndum
  real Lx, Ly, Lz, x, y, z, Rcut1, rx, ry, rz, rrr, COM(3), sumx, sumy, sumz, largest
  integer i, j, N, atom_nr, molec_nr, So_atoms, Sv_atoms, Sv_mols, ii, jj, So_ctr, Sv_ctr, cnt, N_conf, cf, index, extra
  real, allocatable :: So_xyz(:,:), Sv_xyz(:,:), NSo_xyz(:,:), NSv_xyz(:,:), flag(:)

  open(9, file="oniom.inp", status='old')

  read(9,*) N_conf !Number of configurations (< 1000)
  read(9,*) outfile !Output file motif (max 15 chars)
  read(9,*) Rcut1 !Cutoff radius for shell

  Ndum = 0
  read(9,*)  !Solute block
  read(9,*) So_atoms, So_ctr  !Number of solute atoms, Central atom
  do i=1,So_atoms !notice fixed format for input
     read(9,'(a2,x,a2,2f7.4,3a9,x,i1)') Elem_So(i), Type_So(i), Sigma_So(i), Eps_So(i), q1_So(i), q2_So(i), q3_So(i), is_dummy(i)
     Head1_So(i) = trim(adjustl(Elem_So(i)))//"-"//trim(adjustl(Type_So(i)))//"-"//trim(adjustl(q1_So(i)))//"  0"
     Head2_So(i) = trim(adjustl(Elem_So(i)))//"-"//trim(adjustl(Type_So(i)))//"-"//trim(adjustl(q2_So(i)))//"  0"
     Head3_So(i) = trim(adjustl(Elem_So(i)))//"-"//trim(adjustl(Type_So(i)))//"-"//trim(adjustl(q3_So(i)))//"  0"
     Sigma_So(i) = Sigma_So(i)*0.561231 !convert to Amber distance
     Eps_So(i) = Eps_So(i)/4.184 !convert to kcal/mol
     if (is_dummy(i) == 1) Ndum = Ndum+1
     print *,is_dummy(i)
  end do
  print *,Ndum

  read(9,*)  !Solvent block
  read(9,*) Sv_atoms, Sv_ctr !Number of solvent atoms, central atom
  do i=1,Sv_atoms !notice fixed format for input
     read(9,'(a2,x,a2,2f7.4,3a9)') Elem_Sv(i), Type_Sv(i), Sigma_Sv(i), Eps_Sv(i), q1_Sv(i), q2_Sv(i), q3_Sv(i)
     Head1_Sv(i) = trim(adjustl(Elem_Sv(i)))//"-"//trim(adjustl(Type_Sv(i)))//"-"//trim(adjustl(q1_Sv(i)))//" -1"
     Head2_Sv(i) = trim(adjustl(Elem_Sv(i)))//"-"//trim(adjustl(Type_Sv(i)))//"-"//trim(adjustl(q2_Sv(i)))//" -1"
     Head3_Sv(i) = trim(adjustl(Elem_Sv(i)))//"-"//trim(adjustl(Type_Sv(i)))//"-"//trim(adjustl(q3_Sv(i)))//" -1"
     Sigma_Sv(i) = Sigma_Sv(i)*0.561231 !convert to Amber distance
     Eps_Sv(i) = Eps_Sv(i)/4.184 !convert to kcal/mol
  end do
  read(9,*) Sv_mols !Number of water molecules

  do cf=1,N_conf !Loop over all configurations

     if (cf < 10) then
        write (num,"(I1)") cf
     elseif (cf < 100) then
        write (num,"(I2)") cf
     else
        write (num,"(I3)") cf
     endif

     datafile = "conf_" // trim(num) // ".gro" !Hard-coded file name for configurations
     print *, datafile

     outfile1 = trim(outfile) // "_c" // trim(num) // '_q1.inp'
     outfile2 = trim(outfile) // "_c" // trim(num) // '_q2.inp'
     outfile3 = trim(outfile) // "_c" // trim(num) // '_q3.inp'
     outfile4 = trim(outfile) // "_c" // trim(num) // '_q1_chg.inp'
     outfile5 = trim(outfile) // "_c" // trim(num) // '_q2_chg.inp'
     outfile6 = trim(outfile) // "_c" // trim(num) // '_q3_chg.inp'

     open(10, file=datafile, status='old')
     open(11, file=outfile1, status='replace')
     open(12, file=outfile2, status='replace')
     open(13, file=outfile3, status='replace')
     open(14, file=outfile4, status='replace')
     open(15, file=outfile5, status='replace')
     open(16, file=outfile6, status='replace')

     !!Basic settings

     read (10,*) title
     read (10,*) N

     allocate (So_xyz(So_atoms,3))
     allocate (Sv_xyz(Sv_mols*Sv_atoms,3))
     allocate (NSo_xyz(So_atoms,3))
     allocate (NSv_xyz(Sv_mols*Sv_atoms,3))
     allocate (flag(Sv_mols))

     ! Check for consistency
     if (N .ne. So_atoms+Sv_atoms*Sv_mols ) then
        print*, 'ERROR: number of atoms does not match!'
        stop
     endif

     ! Read solute coordinates
     do i=1,So_atoms
        read (10,'(i5,2a5,i5,3f8.3,a)') molec_nr, molec_name, So_name, atom_nr, x, y, z, dummy
        So_xyz(i,1) = x
        So_xyz(i,2) = y
        So_xyz(i,3) = z
     end do

     ! Read solvent coordinates
     do i=1,Sv_mols
        flag(i)=0 ! reset flag
        do j=1,Sv_atoms
           read (10,'(i5,2a5,i5,3f8.3,a)') molec_nr, molec_name, Sv_name, atom_nr, x, y, z, dummy
           Sv_xyz((i-1)*Sv_atoms+j,1) = x
           Sv_xyz((i-1)*Sv_atoms+j,2) = y
           Sv_xyz((i-1)*Sv_atoms+j,3) = z
        end do
     end do

     ! Read Box dimensions
     read (10,*) Lx, Ly, Lz

     ! Write solute coordinates
     do i=1,So_atoms
        ! Translate coordinates to place central solute atom at origin; then convert units to Angstrom
        rx = So_xyz(i,1) - So_xyz(So_ctr,1)
        ! PBC
        if (rx > 0.5*Lx) then
           NSo_xyz(i,1) = (rx - Lx)*10.
        elseif (-rx > 0.5*Lx) then
           NSo_xyz(i,1) = (rx + Lx)*10.
        else
           NSo_xyz(i,1) = rx*10.
        end if
        ry = So_xyz(i,2) - So_xyz(So_ctr,2)
        if (ry > 0.5*Ly) then
           NSo_xyz(i,2) = (ry - Ly)*10.
        elseif (-ry > 0.5*Ly) then
           NSo_xyz(i,2) = (ry + Ly)*10.
        else
           NSo_xyz(i,2) = ry*10.
        end if
        rz = So_xyz(i,3) - So_xyz(So_ctr,3)
        if (rz > 0.5*Lz) then
           NSo_xyz(i,3) = (rz - Lz)*10.
        elseif (-rz > 0.5*Lz) then
           NSo_xyz(i,3) = (rz + Lz)*10.
        else
           NSo_xyz(i,3) = rz*10.
        end if
        if (is_dummy(i) == 0) then ! do not write coordinates of dummy atoms
           write (11,'(a17,x,3f8.3,a2)') Head1_So(i), NSo_xyz(i,1), NSo_xyz(i,2), NSo_xyz(i,3), " H"
           write (12,'(a17,x,3f8.3,a2)') Head2_So(i), NSo_xyz(i,1), NSo_xyz(i,2), NSo_xyz(i,3), " H"
           write (13,'(a17,x,3f8.3,a2)') Head3_So(i), NSo_xyz(i,1), NSo_xyz(i,2), NSo_xyz(i,3), " H"
        end if
     end do

     cnt=1
     largest=0.0
     ! Find solvent within cutoff radius
     do i=1,Sv_mols

        ! coordinates of central solvent atom
        j=(i-1)*Sv_atoms+Sv_ctr
        ! distance to coordinates of solute COM
        rx = Sv_xyz(j,1) - So_xyz(So_ctr,1)
        ! PBC
        if (rx > 0.5*Lx) then
           rx = rx - Lx
        elseif (-rx > 0.5*Lx) then
           rx = rx + Lx
        end if
        ry = Sv_xyz(j,2) - So_xyz(So_ctr,2)
        if (ry > 0.5*Ly) then
           ry = ry - Ly
        elseif (-ry > 0.5*Ly) then
           ry = ry + Ly
        end if
        rz = Sv_xyz(j,3) - So_xyz(So_ctr,3)
        if (rz > 0.5*Lz) then
           rz = rz - Lz
        elseif (-rz > 0.5*Lz) then
           rz = rz + Lz
        end if
        ! Total distance
        rrr = sqrt(rx*rx + ry*ry + rz*rz)

        ! If distance is within cutoff, write coordinates
        if (rrr <= Rcut1) then
           cnt = cnt+1
           flag(i)=1
           if (rrr>largest) then
              largest=rrr
              index=i
           end if
        end if
     end do

     print*, largest, "  ", index
     if (mod(cnt,2) == 0) then
        flag(index)=0 ! remove furthest molecule if count is even (Gaussian bug)
        cnt = cnt-1
        print*, "REMOVED"
     end if

     ! write molecules
     do i=1,Sv_mols
        ! If distance is within cutoff, write coordinates
        if (flag(i) == 1) then
           do ii = 1,Sv_atoms
              jj = (i-1)*Sv_atoms+ii
              ! distance to coordinates of solute COM
              rx = Sv_xyz(jj,1) - So_xyz(So_ctr,1)
              ! PBC
              if (rx > 0.5*Lx) then
                 NSv_xyz(jj,1) = (rx - Lx)*10.
              elseif (-rx > 0.5*Lx) then
                 NSv_xyz(jj,1) = (rx + Lx)*10.
              else
                 NSv_xyz(jj,1) = rx*10.
              end if
              ry = Sv_xyz(jj,2) - So_xyz(So_ctr,2)
              if (ry > 0.5*Ly) then
                 NSv_xyz(jj,2) = (ry - Ly)*10.
              elseif (-ry > 0.5*Ly) then
                 NSv_xyz(jj,2) = (ry + Ly)*10.
              else
                 NSv_xyz(jj,2) = ry*10.
              end if
              rz = Sv_xyz(jj,3) - So_xyz(So_ctr,3)
              if (rz > 0.5*Lz) then
                 NSv_xyz(jj,3) = (rz - Lz)*10.
              elseif (-rz > 0.5*Lz) then
                 NSv_xyz(jj,3) = (rz + Lz)*10.
              else
                 NSv_xyz(jj,3) = rz*10.
              end if
              write (11,'(a17,x,3f8.3,a2)') Head1_Sv(ii), NSv_xyz(jj,1), NSv_xyz(jj,2), NSv_xyz(jj,3), " L"
              write (12,'(a17,x,3f8.3,a2)') Head2_Sv(ii), NSv_xyz(jj,1), NSv_xyz(jj,2), NSv_xyz(jj,3), " L"
              write (13,'(a17,x,3f8.3,a2)') Head3_Sv(ii), NSv_xyz(jj,1), NSv_xyz(jj,2), NSv_xyz(jj,3), " L"
              write (14,'(3f8.3,x,a10)') NSv_xyz(jj,1), NSv_xyz(jj,2), NSv_xyz(jj,3), q1_Sv(ii)
              write (15,'(3f8.3,x,a10)') NSv_xyz(jj,1), NSv_xyz(jj,2), NSv_xyz(jj,3), q2_Sv(ii)
              write (16,'(3f8.3,x,a10)') NSv_xyz(jj,1), NSv_xyz(jj,2), NSv_xyz(jj,3), q3_Sv(ii)
           end do
        end if
     end do
     write (11,*)
     write (12,*)
     write (13,*)
     write (14,*)
     write (15,*)
     write (16,*)

     ! connectivity for first molecule (with dummies)
     do i=1,So_atoms
        if (is_dummy(i) == 0) then ! do not write connectivity for dummy atoms
           write(11,'(i5)',advance='no') i
           write(12,'(i5)',advance='no') i
           write(13,'(i5)',advance='no') i
           write (11,*)
           write (12,*)
           write (13,*)
        end if
     end do

     ! write connectivity information for solvent
     do i=2,cnt
        do ii=1,Sv_atoms
           write(11,'(i5)',advance='no') (i-2)*Sv_Atoms+So_Atoms+ii-Ndum
           write(12,'(i5)',advance='no') (i-2)*Sv_Atoms+So_Atoms+ii-Ndum
           write(13,'(i5)',advance='no') (i-2)*Sv_Atoms+So_Atoms+ii-Ndum
           write (11,*)
           write (12,*)
           write (13,*)
        end do
     end do
     write (11,*)
     write (12,*)
     write (13,*)

     do i=1,So_atoms
        write (11,*) "VDW ", Type_So(i), " ", Sigma_So(i), " ", Eps_So(i)
        write (12,*) "VDW ", Type_So(i), " ", Sigma_So(i), " ", Eps_So(i)
        write (13,*) "VDW ", Type_So(i), " ", Sigma_So(i), " ", Eps_So(i)        
     end do
     do i=1,Sv_atoms
        write (11,*) "VDW ", Type_Sv(i), " ", Sigma_Sv(i), " ", Eps_Sv(i)
        write (12,*) "VDW ", Type_Sv(i), " ", Sigma_Sv(i), " ", Eps_Sv(i)
        write (13,*) "VDW ", Type_Sv(i), " ", Sigma_Sv(i), " ", Eps_Sv(i)        
     end do
     write (11,*)
     write (12,*)
     write (13,*)

     print*, 'Total solvent molecules = ', cnt
     print*, 'Total atoms = ', cnt*Sv_atoms + So_atoms
     !print*, 'Dont forget to change the number of atoms in output file!'

     ! Write Box dimensions
     !write (11,*) Lx, Ly, Lz

     close (9)
     close (10)
     close (11)
     close (12)
     close (13)
     close (14)
     close (15)
     close (16)

     deallocate (So_xyz)
     deallocate (Sv_xyz)
     deallocate (NSo_xyz)
     deallocate (NSv_xyz)
     deallocate (flag)

  end do

end program
