c-----------------------------------------------------------------------
c
c       Mrg13  -  Merge W1 and W3 template-fit photometry for Option 1b
c                 (some cloning from mrgad.f)
c
c  version 1.0  B71225: initial version
c
c-----------------------------------------------------------------------
c
                     Integer*4  MaxFld
                     Parameter (MaxFld = 1000)
c
      Character*5000 Line
      Character*500  InFNam, OutFNam
      Character*25   Field(MaxFld)      
      Character*11   Vsn, Fmt, IFmt
      Character*8    CDate, CTime
      Character*3    Flag, Flag0
      Real*8         R8tmp1, R8tmp2, v1, v2
      Integer*4      IArgC, LNBlnk, NArgs, NArg, nSrc, Access, nIn,
     +               nOut, IFa(MaxFld), IFb(MaxFld), NF, k, Itmp1,
     +               Itmp2, w1M, w3M
      Logical*4      dbg, GotIn, GotOut, GotHdr, Good1, Good2      
c
      Data Vsn/'1.0  B71225'/, nSrc/0/,dbg,GotIn,GotOut/3*.false./,
     +     nIn,nOut/2*0/, GotHdr/.false./ 
c
      Common / VDT / CDate, CTime, Vsn
c
c-----------------------------------------------------------------------
c
      NArgs = IArgC()
      If (NArgs .lt. 4) then
        print *,'Mrg13 vsn ', Vsn
        print *
        print *,'Usage: mrg13 <flags specifications>'
        print *
        print *,'Where the REQUIRED flags and specifications are:'
        print *,'    -i  name of the Option-1b mdex file'
        print *,'    -o  name of the W1/W3 merged output file'
        print *
        print *,'The OPTIONAL flag is:'
        print *,'    -d  turn on debug prints'
        call exit(32)
      end if
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      call signon('mrg13')
c
      NArg = 0
5     NArg = NArg + 1
      call GetArg(NArg,Flag)
      Flag0 = Flag
      call UpCase(Flag)
c                                      ! input mdex file
      If (Flag .eq. '-I') then
        call NextNarg(NArg,Nargs)
        call GetArg(NArg,InFNam)
        if (Access(InFNam(1:lnblnk(InFNam)),' ') .ne. 0) then
          print *,'File not found: ',InFNam(1:lnblnk(InFNam))
          call exit(64)
        end if
        GotIn = .true.
        if (dbg) print *,'option-1b mdex file: '
     +           //InFNam(1:lnblnk(InFNam))
c                                      ! Turn debug prints on
      else if (Flag .eq. '-D') then
        dbg = .true.
        print *,'Debug prints enabled'
c                                      ! output   file
      else if (Flag .eq. '-O') then
        call NextNarg(NArg,Nargs)
        call GetArg(NArg,OutFNam)
        GotOut = .true.
        if (dbg) print *,'W1/W3 merged output file: '
     +           //OutFNam(1:lnblnk(OutFNam))
c        
      Else
        print *,'ERROR: unrecognized command-line specification: '
     +          //Flag0
      end if
c 
      If (NArg .lt. NArgs) Go to 5
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c                                      ! Verify input
      If (.not.GotIn) then
        print *,'ERROR: Input mdex file not specified'
        call exit(64)
      end if
c
      If (.not.GotOut) then
        print *,'ERROR: Output filename stem not specified'
        call exit(64)
      end if
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c                                      ! Get source count
      open (10, file = InFNam)
      read (10, '(a)', end=3000) Line
      if (Line(1:7) .ne. '\Nsrc =') then
        print *,
     +    'ERROR: "\Nsrc" not found on first header line of input file:'
        print *,'       '//InFNam(1:lnblnk(InFNam))
        print *,'       Line found: ',Line(1:lnblnk(Line))
        print *,'Execution aborted'
        call exit(64)
      end if
      read (Line(8:lnblnk(Line)), *, err=3001) nSrc
      if (dbg) print *,'No. of sources: ', nSrc      
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c                                      ! Create output table header
      open (12, file = OutFNam)
20    write(12,'(a)') Line(1:lnblnk(Line))
25      read (10, '(a)', end=3000) Line
      if (Line(1:1) .eq. '\') go to 20
      if (Line(1:1) .eq. '|') then
        if (GotHdr) go to 20
        GotHdr = .true.
        call GetFlds(Line,Field,IFa,IFb,NF)
        if (dbg) print *,'No. fields returned by GetFlds:', NF
c                                          ! verify fields in mdex file
        call ChkFld(Field(3),'ra',3)
        call ChkFld(Field(4),'dec',4)
        call ChkFld(Field(5),'sigra',5)
        call ChkFld(Field(6),'sigdec',6)
        call ChkFld(Field(7),'sigradec',7)
        go to 20
      end if
c
      nIn = nIn + 1
      Good1 = index(Line(IFA(32):IFB(32)),'null') .eq. 0  ! w1snr
      Good2 = index(Line(IFA(34):IFB(34)),'null') .eq. 0  ! w3snr
      k = 32
      if (Good1) read (Line(IFA(32):IFB(32)), *, err=3006) R8tmp1
      k = 34
      if (Good2) read (Line(IFA(34):IFB(34)), *, err=3006) R8tmp2
      if (Good1 .and. Good2) then
        R8tmp1 = sqrt(R8tmp1**2 + R8tmp2**2)
        write (Line(IFA(32):IFB(32)),'(F7.1)') R8tmp1
      else if (Good2) then
        Line(IFA(32):IFB(32)) = Line(IFA(34):IFB(34))
      end if
c
      Good1 = index(Line(IFA(36):IFB(37)),'null') .eq. 0  ! w1flux,w1sigflux
      Good2 = index(Line(IFA(40):IFB(41)),'null') .eq. 0  ! w3flux,w3sigflux
      if (Good1 .and. Good2) then
        k = 37
        read (Line(IFA(k):IFB(k)), *, err=3006) R8tmp1    ! w1sigflux
        k = 41
        read (Line(IFA(k):IFB(k)), *, err=3006) R8tmp2    ! w3sigflux
        v1 = R8tmp1**2
        v2 = R8tmp2**2
        k = 36
        Read(Line(IFA(k):IFB(k)), *, err = 3006) R8tmp1   ! w1flux
        k = 40
        Read(Line(IFA(k):IFB(k)), *, err = 3006) R8tmp2   ! w3flux
        R8tmp1 = (v2*R8tmp1 + v1*R8tmp2)/(v1+v2)
        write (Line(IFA(36):IFB(36)), '(1pe12.4)') R8tmp1
        R8tmp1 = dsqrt(v1*v2/(v1+v2))
        write (Line(IFA(37):IFB(37)), '(1pe14.4)') R8tmp1
      else if (Good2) then
        Line(IFA(36):IFB(37)) = Line(IFA(40):IFB(41))
      end if
c
      Good1 = index(Line(IFA(44):IFB(46)),'null') .eq. 0  ! w1mpro-w1rchi2
      Good1 = Good1 .and. (index(Line(IFA(184):IFB(185)),'null') .eq. 0) ! w1NM,w1M
      Good2 = index(Line(IFA(50):IFB(52)),'null') .eq. 0  ! w3mpro-w3rchi2
      Good2 = Good2 .and. (index(Line(IFA(206):IFB(207)),'null') .eq. 0) ! w3NM,w3M
      if (Good1 .and. Good2) then
        k = 45
        read (Line(IFA(k):IFB(k)), *, err=3006) R8tmp1    ! w1sigmpro
        k = 51
        read (Line(IFA(k):IFB(k)), *, err=3006) R8tmp2    ! w3sigmpro
        v1 = R8tmp1**2
        v2 = R8tmp2**2
        k = 44
        Read(Line(IFA(k):IFB(k)), *, err = 3006) R8tmp1   ! w1mpro
        k = 50
        Read(Line(IFA(k):IFB(k)), *, err = 3006) R8tmp2   ! w3mpro
        R8tmp1 = (v2*R8tmp1 + v1*R8tmp2)/(v1+v2)
        write (Line(IFA(44):IFB(44)), '(F7.3)') R8tmp1
        R8tmp1 = dsqrt(v1*v2/(v1+v2))
        write (Line(IFA(45):IFB(45)), '(F10.3)') R8tmp1
        k = 46
        Read(Line(IFA(k):IFB(k)), *, err = 3006) R8tmp1   ! w1rchi2
        k = 52
        Read(Line(IFA(k):IFB(k)), *, err = 3006) R8tmp2   ! w3rchi2
        k = 185
        Read(Line(IFA(k):IFB(k)), *, err = 3006) w1M      ! w1M
        k = 207
        Read(Line(IFA(k):IFB(k)), *, err = 3006) w3M      ! w3M
        Itmp1 = w1M + w3M
        k = 46
        if (Itmp1 .gt. 0) then
          R8tmp1 = (dfloat(w1M)*R8tmp1 + dfloat(w3M)*R8tmp2)
     +   /dfloat(Itmp1)
          write (Line(IFA(k):IFB(k)),'(1pe11.3)') R8tmp1
        end if                         ! else leave as "null"
        k = 185
        write (Line(IFA(k):IFB(k)),'(I6)') Itmp1
        k = 206
        Read(Line(IFA(k):IFB(k)), *, err = 3006) Itmp1    ! w3NM
        k = 184
        Read(Line(IFA(k):IFB(k)), *, err = 3006) Itmp2    ! w1NM
        Itmp1 = Itmp1 + Itmp2
        write (Line(IFA(k):IFB(k)),'(I7)') Itmp1
      else if (Good2) then
        Line(IFA(44):IFB(46))   = Line(IFA(50):IFB(52))
        Line(IFA(184):IFB(185)) = Line(IFA(206):IFB(207))
      end if
c
      Good1 = index(Line(IFA(192):IFB(194)),'null') .eq. 0 ! w1mJDmin-w1mJDmean
      Good2 = index(Line(IFA(214):IFB(216)),'null') .eq. 0 ! w3mJDmin-w3mJDmean
      if (Good1 .and. Good2) then
        k = 214
        read(Line(IFA(k):IFB(k)), *, err = 3006) R8tmp2  !  w3mJDmin
        k = 192
        read(Line(IFA(k):IFB(k)), *, err = 3006) R8tmp1  !  w1mJDmin
        if (R8tmp2 .lt. R8tmp1) then
          Line(IFA(192):IFB(192)) = Line(IFA(214):IFB(214))
          v1 = R8tmp2
        else
          v1 = R8tmp1
        end if
        k = 215
        read(Line(IFA(k):IFB(k)), *, err = 3006) R8tmp2  !  w3mJDmax
        k = 193
        read(Line(IFA(k):IFB(k)), *, err = 3006) R8tmp1  !  w1mJDmax
        if (R8tmp2 .gt. R8tmp1) then
          Line(IFA(193):IFB(193)) = Line(IFA(215):IFB(215))
          v2 = R8tmp2
        else
          v2 = R8tmp1
        end if
        R8tmp1 = (v1+v2)/2.0d0        
        write(Line(IFA(194):IFB(194)),'(F18.8)') R8tmp1
      else if (Good2) then
        Line(IFA(192):IFB(194)) = Line(IFA(214):IFB(216))
      end if
c
      Good1 = index(Line(IFA(274):IFB(274)),'null') .eq. 0  ! w1snr_pm
      Good2 = index(Line(IFA(276):IFB(276)),'null') .eq. 0  ! w3snr_pm
      if (Good1 .and. Good2) then
        k = 276
        read (Line(IFA(k):IFB(k)), *, err = 3006) R8tmp2  ! w3snr_pm
        k = 274
        read (Line(IFA(k):IFB(k)), *, err = 3006) R8tmp1  ! w1snr_pm
        R8tmp1 = dsqrt(R8tmp1**2 + R8tmp2**2)
        write (Line(IFA(k):IFB(k)), '(F9.1)') R8tmp1
      else if (Good2) then
        Line(IFA(274):IFB(274)) = Line(IFA(276):IFB(276))
      end if
c      
      Good1 = index(Line(IFA(278):IFB(279)),'null') .eq. 0  ! w1flux_pm-w1sigflux_pm
      Good2 = index(Line(IFA(282):IFB(283)),'null') .eq. 0  ! w3flux_pm-w3sigflux_pm
      if (Good1 .and. Good2) then
        k = 283
        read (Line(IFA(k):IFB(k)), *, err = 3006) R8tmp2  ! w3sigflux_pm
        k = 279
        read (Line(IFA(k):IFB(k)), *, err = 3006) R8tmp1  ! w1sigflux_pm
        v1 = R8tmp1**2
        v2 = R8tmp2**2
        k = 282
        read (Line(IFA(k):IFB(k)), *, err = 3006) R8tmp2  ! w3flux_pm
        k = 278
        read (Line(IFA(k):IFB(k)), *, err = 3006) R8tmp1  ! w1flux_pm
        R8tmp1 = (v2*R8tmp1+v1*R8tmp2)/(v1+v2)
        write (Line(IFA(k):IFB(k)), '(1pE12.4)') R8tmp1
        R8tmp1 = dsqrt(v1*v2/(v1+v2))
        k = 279
        write (Line(IFA(k):IFB(k)), '(1pE14.4)') R8tmp1
      else if (Good2) then
        Line(IFA(278):IFB(279)) = Line(IFA(282):IFB(283))
      end if
c
      Good1 = index(Line(IFA(286):IFB(288)),'null') .eq. 0  ! w1mpro_pm-w1rchi2_pm
      Good2 = index(Line(IFA(292):IFB(294)),'null') .eq. 0  ! w3mpro_pm-w3rchi2_pm
      if (Good1 .and. Good2) then
        k = 293
        read (Line(IFA(k):IFB(k)), *, err = 3006) R8tmp2  ! w3sigmpro_pm
        k = 287
        read (Line(IFA(k):IFB(k)), *, err = 3006) R8tmp1  ! w1sigmpro_pm
        v1 = R8tmp1**2
        v2 = R8tmp2**2
        k = 292
        read (Line(IFA(k):IFB(k)), *, err = 3006) R8tmp2  ! w3mpro_pm
        k = 286
        read (Line(IFA(k):IFB(k)), *, err = 3006) R8tmp1  ! w1mpro_pm
        R8tmp1 = (v2*R8tmp1+v1*R8tmp2)/(v1+v2)
        write (Line(IFA(k):IFB(k)), '(F10.3)') R8tmp1
        R8tmp1 = dsqrt(v1*v2/(v1+v2))
        k = 287
        write (Line(IFA(k):IFB(k)), '(F13.3)') R8tmp1
        k = 294
        read (Line(IFA(k):IFB(k)), *, err = 3006) R8tmp2  ! w3rchi2_pm
        k = 288
        read (Line(IFA(k):IFB(k)), *, err = 3006) R8tmp1  ! w1rchi2_pm
        Itmp1 = w1M + w3M
        if (Itmp1 .gt. 0) then
          R8tmp1 = (dfloat(w1M)*R8tmp1 + dfloat(w3M)*R8tmp2)
     +   /dfloat(Itmp1)
          write (Line(IFA(k):IFB(k)),'(1pe11.3)') R8tmp1
        end if                         ! else leave as "null"
      else if (Good2) then
        Line(IFA(156):IFB(158)) = Line(IFA(322):IFB(324))
      end if
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c                                      ! Output data row
c
      write (12,'(a)') Line(1:lnblnk(Line))
      nOut = nOut + 1
      if (nOut .lt. nSrc) go to 25
c
      print *,'No. rows read and output:', nIn, nOut
c
      call signoff('mrg13')
      stop
c
3000  print *,'ERROR: Unexpected Eof-of-file encountered in'
     +      //' mdex input file'
      print *,'No. data rows read and output:', nIn, nOut
      go to 3333
c
3001  print *,'ERROR: read error encountered in mdex file line:'
      print *,Line(1:lnblnk(Line))
      go to 3333
c
3006  print *,'ERROR: read error encountered in mdex input file '
     +      //'on data line no.',nIn
      print *,'       Data field no.',k,', column name "',
     +                Field(k)(1:lnblnk(Field(k))),'"'
      print *,'        Numeric field: "',Line(IFA(k):IFB(k)),'"'
c
3333  print *,'Execution aborted'      
      call exit(64)
      stop
      end
c      
c=======================================================================
c
      subroutine SignOn(pgmnam)
c
c *** signon- routine which provides sign-on and sign-off messages
c             (orig by John Fowler- mod by Howard McCallon-041214-SIRTF)
c
c     inputs:  pgmnam = program name                                 [call arg]
c
c     outputs: message to stdout
c
      character*(*) pgmnam
      character vsn*11,cdate*8,ctime*8,Fmt*11,FLen*4
      integer*4 onoff,jdate(3),jtime(3),lnblnk
      real*4    dummyt,second(2),etime
c
      common /vdt/ cdate,ctime,vsn
c##
      onoff = 1
c
c         i. obtain date
c
100   cdate = '00-00-00'
      call idate(jdate)    ! Linux call
c
      jdate(3) = mod(jdate(3), 100)
      write(cdate(1:2), '(i2)') jdate(2)
      write(cdate(4:5), '(i2)') jdate(1)
      write(cdate(7:8), '(i2)') jdate(3)
c
      if(cdate(4:4) .eq. ' ') cdate(4:4) = '0'
      if(cdate(7:7) .eq. ' ') cdate(7:7) = '0'
c
c         ii. obtain time
c
      ctime = '00:00:00'
      call itime(jtime)
      write(ctime(1:2), '(i2)') jtime(1)
      write(ctime(4:5), '(i2)') jtime(2)
      write(ctime(7:8), '(i2)') jtime(3)
c
      if(ctime(4:4) .eq. ' ') ctime(4:4) = '0'
      if(ctime(7:7) .eq. ' ') ctime(7:7) = '0'
c
c         iii. set up format for pgmnam
c
      write(Flen,'(I4)') lnblnk(pgmnam)
      Fmt = '(A'//Flen//'$)'
c
c         iv. write out results
c
      write(*,Fmt) pgmnam
      if(onoff .eq. 1) then                      ! sign on
        write(*,301) vsn,cdate,ctime
      else                                       ! sign off
        dummyt = etime(second)
        write(*,302) vsn,cdate,ctime,second
      endif
  301 format(' version: ',a11,' - execution begun on ',a8,' at ',a8)
  302 format(' version: ',a11,' - execution ended on ',a8,' at ',a8
     *    /1x,f9.2,' cpu seconds used;',f8.2,' system seconds used.')
c
      return
c
      entry SignOff(pgmnam)
      OnOff = 2
      go to 100
c
      end
c      
c=======================================================================
c
      subroutine GetFlds(ColNam,Field,IFa,IFb,NF)
c-----------------------------------------------------------------------
c
c  Get fields in a table-file header line
c
c-----------------------------------------------------------------------
                     Integer*4  MaxFld
                     Parameter (MaxFld = 1000)
c
      character*5000 ColNam
      Character*300  Line
      character*25   Field(MaxFld)
      integer*4      IFa(MaxFld), IFb(MaxFld), NF, N, M, L, K, LNBlnk,
     +               LastErr
c
c-----------------------------------------------------------------------
c
      N = 0
      K = 0
      LastErr = 0
      do 100 M = 1, LNBlnk(ColNam)
        if (ColNam(M:M) .eq. '|') then
          N = N + 1
          NF = N - 1
          if (N .gt. 1) IFb(N-1) = M-1
          if (N .gt. MaxFld) return
          IFa(N) = M
          do 10 L = 1, 25
            Field(N)(L:L) = ' '
10        continue
          K = 0
        else
          if (ColNam(M:M) .ne. ' ') then
            K = K + 1
            if (K .le. 25) then
              Field(N)(K:K) = ColNam(M:M)
            else
              if (LastErr .ne. N) then
                write(Line,*) N
                Line = 'GetFlds - Table column name no. '
     +               //Line(1:lnblnk(Line))//' longer than 25 '
     +               //'characters: '//Field(N)//'....; excess ignored'
                print *,Line(1:lnblnk(line))
                LastErr = N
              end if
            end if
          end if
        end if
100   continue
c
      return
      end
c
c=======================================================================
c
      Subroutine NextNarg(NArg,NArgs)
c
      integer NArg, NArgs
c
c-----------------------------------------------------------------------
c
      if (NArg .lt. NArgs) then
        NArg = NArg + 1
        return
      else
        print *,'ERROR: expected another argument but none found'
        call exit(64)
      end if
      return
      end
c
c=======================================================================
c
      subroutine upcase(string)
      character*(*) string
      integer*4 j, lnblnk
c
      do 10 j = 1,lnblnk(string)
         if(string(j:j) .ge. "a" .and. string(j:j) .le. "z") then
            string(j:j) = achar(iachar(string(j:j)) - 32)
         end if
10    continue
      return
      end
c
c=======================================================================
c
      subroutine ChkFld(Fld1, Fld2, k)
c      
      character*(*) Fld1, Fld2
      integer*4     k, lnblnk      
c
      if (Fld1 .ne. Fld2) then
        print *,'ERROR: input field no.',k,' expected to be ',
     +           Fld2(1:lnblnk(Fld2)),'; got ',Fld1(1:lnblnk(Fld2))
        call exit(64)
      end if
c      
      return
c      
      end
