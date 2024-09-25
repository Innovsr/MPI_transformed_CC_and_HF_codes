c   dimension
#define   MNS   11
#define MAXVAL  10
#define MNOCC   30
#define MNBAS  700
#define   MNB  700
#define    MN  750
c   files
#define  STDIN   5
#define STDOUT   4
#define IPUNCH   7
#define   IWFN   9
c   variables
#define  TINY  1.0D-10
#define  MASTER  0

c---------------------------------------------------------------------------
c         Relavistic SCF code (based on gaussian basis) 
c  This code takes three different types of basis 
c (i)    Universal Basis set (UB)
c (ii)   Well Tempered (WT)
c (iii)  Even Tempered (ET)
c  geetha 11/7/2000
c---------------------------------------------------------------------------
c             sample 'inputfile'
c  ba.dat
c  ba.out
c-------------------------------------------------------------------------
c 1st card: title
c 2nd card: name of the atom in captial
c 3rd card: total number symetry,  total s,  total p_1/2,  total p_3/2, ....
c 4th card:                         occ  s,    occ p_1/2,   occ  P_3/2,.....
c 5th card:    Kappa values
c 6th card:    J values
c 7th card:    Which method to use 
c   UB(Universal Basis)  WT (Well Tempered)  ET (Even Tempered)
c 8th card:
c  if(UB)   read alpha1,beta1
c  if(WT)   read alpha2,beta2,gamma,delta
c  if(ET)   read alpha3(isym) for each symmetry 
c           read beta3(isym)  for each symmetry
c 9th card:    charge
c 10th card:    number iteration, accuracy (give the integervalue)
c 11th card:    printing
c 12th card:   point nuclei (0) ferminuclei (1)
c 13th card:   Atomic weight and atomic number 
c-------------------------------------------------------------------------
c   Also one has to give a filename 'inputfile'
c   this contains two lines only
c 1st card:   name of the inputfile
c 2nd card:   name of the outputfile name
c-------------------------------------------------------------------------
      program main
      implicit real*8(a-h,o-z)
      include 'mpif.h'
      character*2 method 
   

      common/quan1/kap(MNS),nak(MNOCC),np(MNOCC),nsy(MNS),nkk(MNS)
      common /grid/steplin,r(MN),rp(MN),rpor(MN),rnt,h,hp,n
     :      /adat/occj(MNOCC),occl(MNOCC),unocl(MNS),unocj(MNS)
     :  /UB/alpha1,beta1
     :  /WT/alpha2,beta2,gamma,delta
     :  /ET/alpha3(MNS),beta3(MNS)

      data  (nkk(i),i=1,11)/1,2,2,3,3,4,4,5,5,6,6/
      data  (kap(i),i=1,11)/-1,1,-2,2,-3,3,-4,4,-5,5,-6/
      data  (unocj(i),i=1,11)/0.5d00,0.5d00,1.5d00,1.5d00,2.5d00,
     :                       2.5d00,3.5d00,3.5d00,4.5d00,4.5d00,
     :                                       5.5d00/

      call MPI_INIT(ierr)
      START_TIME=MPI_WTIME()
      call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)

      call fzero
      call readinp(mocc,method)
      call setgrd
      call nucparm
      call setqic
      call nucpot
      call gbasis(nbasis,method)
      call rinti(nbasis)
      call scfiter(nbasis)
      call output(nbasis)
      END_TIME=MPI_WTIME()
       if (iproc .eq. 0) then
        print *, 'Time spent (in seconds) is: ',END_TIME-START_TIME
        endif
       if(iproc.eq.0)then
       write(STDOUT,*)'Time spent (in seconds) is: ',
     :  END_TIME-START_TIME
      endif
       call MPI_FINALIZE(ierr)

      stop
      end

      subroutine simpson(tmp,result)
c    this is the simpson integration 

      implicit real*8(a-h,o-z)
      include 'mpif.h'
      real*8 tmp,result
      dimension tmp(MN)
      common /grid/steplin,r(MN),rp(MN),rpor(MN),rnt,h,hp,n
     :  /cons/zero,half,tenth,one,two,three,ten
c
      result=zero
      do l=2,n,2
      result=result+tmp(l)+tmp(l)+tmp(l+1)
      result1=result+tmp(l)+tmp(l)+tmp(l+1)
      enddo
      result=h*(result+result+tmp(1)-tmp(n))/three
      result=result+half*tmp(1)
      return
      end      

      subroutine rinti(nbasis)
c-------------------------------------------------------------------------
c     this generates overlap (sss and sll) and one-electron fock matrix,
c     where sss and sll  stands for small and large component part of the
c     overlap matrix and stored in the array emt.
c     the off-diagonal part of the fock matrix contains the kinetic energy
c     part while the diagonal part contains the nuclear contribution and
c     is stored in the array df_single. If nucopt.ne.0 then finite nucleus
c     will be considered else point nuclues.
c--------------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      include 'mpif.h'
      common/basis/gl(MN,MNBAS),gs(MN,MNBAS),alpha(MNBAS,MNS),mg(MNBAS)
     :  /wave/pf(MN,MNBAS),qf(MN,MNBAS),mf(MNBAS)
     :  /grid/steplin,r(MN),rp(MN),rpor(MN),rnt,h,hp,n
     :  /ncc/c1,c2,c3,c4,cnc5c(5,5),c5num(5,5),pi,c
     :  /overlap/sss(MNB,MNB),sll(MNB,MNB)
     :  /dumm/tmp1(MN),tmp2(MN),tmp3(MN),tmp4(MN)
     :  /quan1/kap(MNS),nak(MNOCC),np(MNOCC),nsy(MNS),nkk(MNS)
     :  /cons/zero,half,tenth,one,two,three,ten
     :  /datt/nsym,nocc,nbas(MNS),ind
     :  /single/df_single(MNB,MNB)
     :  /def1/amass,z
     :  /npot/zz(MN)
     :  /skip/nskip(MNS)
     :  /ovlp/emt(MNB,MNB)
     :  /print/jprint,iout
     :  /nucop/nucopt,isms
     :  /info2/qqf(MN,MNBAS),qql(MN,MNBAS),qqs(MN,MNBAS)
     :  /info3/df_t(MNB,MNB)

      dimension atemp(MNB,MNB)
      do i=1,MNB
      do j=1,MNB
      atemp(MNB,MNB)=0.d0
      enddo
      enddo

      spi = dsqrt(pi)
      four=two*two
      sq2=dsqrt(two)
      
      if(nucopt.eq.0)then
      do i=1,n
      zz(i)=z
      enddo
      endif

      do isym=1,nsym
      do ia=1,nbas(isym)
      ja=ia+nskip(isym)
      do ib=1,nbas(isym)
      jb=ib+nskip(isym)

      ka=ia+2*nskip(isym)
      index1=ia+nbas(isym)+2*nskip(isym)

      kb=ib+2*nskip(isym)
      index2=ib+nbas(isym)+2*nskip(isym)

      do la = 1,n
      tmp1(la)=gl(la,ja)*gl(la,jb)*rp(la)
      tmp2(la)=gs(la,ja)*gs(la,jb)*rp(la)
      tmp3(la)=gs(la,jb)*qqf(la,ja)*rp(la)
      enddo
      call simpson(tmp1,result1)
      call simpson(tmp2,result2)
      call simpson(tmp3,result3)
      sll(ja,jb)=result1
      sss(ja,jb)=result2
      atemp(ja,jb)=atemp(ja,jb)+c*result3
      emt(ka,kb)=result1
      emt(index1,index2) = result2
      enddo
      enddo

      do ia=1,nbas(isym)
      ja=ia+2*nskip(isym)
      ka=ia+nskip(isym)
      index1=ia+nbas(isym)+2*nskip(isym)
      do ib=1,nbas(isym)
      jb=ib+2*nskip(isym)
      kb=ib+nskip(isym)
      index2=ib+nbas(isym)+2*nskip(isym)
      df_single(index1,index2)=df_single(index1,index2)-
     :                         two*c*c*sss(ka,kb)
      enddo
      enddo

      do ia=1,nbas(isym)
      ja=ia+2*nskip(isym)
      ka=ia+nskip(isym)
      do ib=1,nbas(isym)
      kb=ib+nskip(isym)
      index2=ib+nbas(isym)+2*nskip(isym)
      df_single(ja,index2)=df_single(ja,index2)+atemp(ka,kb)
      df_single(index2,ja)=df_single(ja,index2)
      df_t(ja,index2)=df_single(ja,index2)
      df_t(index2,ja)=df_single(ja,index2)
      enddo
      enddo

c     nuclear repulsion energy

      do ia=1,nbas(isym)
      ka=ia+nskip(isym)
      ja=ia+2*nskip(isym)
      index1=ia+nbas(isym)+2*nskip(isym)
      do ib= 1,nbas(isym)
      kb=ib+nskip(isym)
      jb=ib+2*nskip(isym)
      index2=ib+nbas(isym)+2*nskip(isym)
      do la = 1,n
      tmp1(la)=gl(la,ka)*gl(la,kb)*rpor(la)*zz(la)
      tmp2(la)=gs(la,ka)*gs(la,kb)*rpor(la)*zz(la)
      enddo
      call simpson(tmp1,result1)
      call simpson(tmp2,result2)
c     matrix df_single
      df_single(ja,jb)=df_single(ja,jb)-result1
      df_single(index1,index2)=df_single(index1,index2)
     :                         -result2
      enddo
      enddo
      enddo
      mbasis=2*nbasis
      if(jprint.ne.0)then
      call matout(df_single,MNB,mbasis,mbasis,mbasis,4hH   )
      endif
      return
      end
c----------------------------------------------------------------------
      subroutine quad (ta,mtp,result)
c------------------------------------------------------------------------
c   the integral in the interval zero to r(2) is computed by use of an 
c   analytical fit                                                    
c                                sigma                                
c                      f(r) = a r                                     
c   a five-point  closed  newton-cotes  formula (cf. f b hildebrand,  
c   introduction to numerical analysis, second edition, mcgraw-hill,  
c   new york, 1974, p 93)  is  used  to  compute the integral in the  
c   interval  r(2:mtp).                                              
c------------------------------------------------------------------------
      implicit real*8(a-h, o-z)
      include 'mpif.h'
      dimension ta(MN)
      common/cons/zero,half,tenth,one,two,three,ten
     : /grid/steplin,r(MN),rp(MN),rpor(MN),rnt,h,hp,n
     : /ncc/c1,c2,c3,c4,cnc5c(5,5),c5num(5,5),pi,c
c------------------------------------------------------------------------
c   find first values that will permit computation of exponent
c------------------------------------------------------------------------
      result = zero
      mtpm1 = mtp - 1
      do 3 i = 2,mtpm1
        tai = ta(i)
        if (dabs (tai) .gt. zero) then
          ip1   = i+1
          taip1 = ta(ip1)
          quott = taip1/tai
          if (quott .gt. zero) then
c------------------------------------------------------------------------
c                          exponent from fit
c-----------------------------------------------------------------------
               frip1 = taip1/rp(ip1)
               fri   = tai  /rp(i  )
               ratio = frip1/fri
               rip1  = r (ip1)
               ri    = r (i  )
               sigma = dlog (ratio)/dlog (rip1/ri)
c------------------------------------------------------------------------------
c       analytical integration and error estimate for interval r(1:i)
c------------------------------------------------------------------------------
               fri    = ri*fri
               result = fri/(sigma+one)
c------------------------------------------------------------------------------
c                       set the tail to zero
c------------------------------------------------------------------------------
               do 1 loc = 1,3
                  ta(mtp+loc) = zero
    1          continue
c-------------------------------------------------------------------------
c              newton-cotes quadature for the remainder
c------------------------------------------------------------------------
               result = result+c1*tai
               do 2 loc = ip1,mtp,4
                  result = result+c2*(ta(loc  )+ta(loc+2))
     :                           +c3* ta(loc+1)
     :                           +c4* ta(loc+3)
    2          continue
               if (mod (mtp-i,4) .eq. 0) result = result-c1*ta(mtp)
               goto 4
            endif
         endif
    3 continue
c------------------------------------------------------------------------
c          no value which will permit computation of exponent
c-------------------------------------------------------------------------
      result = zero
    4 return
      end
      
      subroutine factt(gam)
c---------------------------------------------------------------------
c calculates the logs  of factorials required by the racah coefficient 
c routine  dracah. written by n.s. scott.
c---------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      include 'mpif.h'
      dimension gam(500)
      common /cons/zero,half,tenth,one,two,three,ten
      data thirty/3.0d 01/
*
      gam(1) = one
      gam(2) = one
      x = two
      do i = 3,30
         gam(i) = gam(i-1)*x
         x = x+one
      enddo
      do  i = 1,30
         gam(i) = dlog(gam(i))
      enddo
      x = thirty
      do i = 31,500
         gam(i) = gam(i-1)+dlog(x)
         x = x+one
      enddo
      return
*
      end
      subroutine init(ii) 
c----------------------------------------------------------------------
c   this subroutine will initialiaze the arrays and the constants in  *
c   dummy blocks dumm and datt                                        *
c----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      include 'mpif.h'
      common /dumm/tmp1(MN),tmp2(MN),tmp3(MN),tmp4(MN)
     :  /cons/zero,half,tenth,one,two,three,ten
     :  /datt/nsym,nocc,nbas(MNS),ind
      do 10 ia=1,MN
      tmp1(ia)=zero
      tmp2(ia)=zero
      tmp3(ia)=zero
      tmp4(ia)=zero
 10   continue
      cn2=zero
      cn3=zero
      cn4=zero
      return
      end
      subroutine scf(nbasis)
c------------------------------------------------------------------
c     This performs the scf iteration.
c------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      include 'mpif.h'
      logical set
      character*3 nh_i
      real*8 df_mat,swap
      dimension can(MNB,MNB)
      dimension  drs(MNB,MNB),sa(MNB),df_mat(MNB,MNB)
      common/basis/gl(MN,MNBAS),gs(MN,MNBAS),alpha(MNBAS,MNS),mg(MNBAS)
      common /grid/steplin,r(MN),rp(MN),rpor(MN),rnt,h,hp,n
     :  /ncc/c1,c2,c3,c4,cnc5c(5,5),c5num(5,5),pi,c
     :  /cons/zero,half,tenth,one,two,three,ten
     :  /overlap/sss(MNB,MNB),sll(MNB,MNB)
     :  /adat/occj(MNOCC),occl(MNOCC),unocl(MNS),unocj(MNS)
     :  /quan1/kap(MNS),nak(MNOCC),np(MNOCC),nsy(MNS),nkk(MNS)
     :  /diag/eigv(MNB,MNB),eng(MNB)
     :  /datt/nsym,nocc,nbas(MNS),ind
     :  /single/df_single(MNB,MNB)
     :  /double/df_two(MNB,MNB)
     :  /wave/pf(MN,MNBAS),qf(MN,MNBAS),mf(MNBAS)
     :  /conv/crit,maxit,npower
     :  /print/jprint,iout
     :  /jv/occ(MNOCC),orbe(MNB),iq(MNB)
     :  /skip/nskip(MNS)
     :  /dens/p(MNB,MNB)
     :  /ovlp/emt(MNB,MNB)
     :  /info1/occk(MNOCC,MNS)
     :  /info/nocorb(MNOCC),nucorb(MNB)
     :  /vector/vec(MNB,MNB),oe(MNB,MNS)
     :  /info3/df_t(MNB,MNB)

      dimension a(MNB,MNB),st(MNB,MNB),b(MNB,MNB)
      dimension pold(MNB,MNB),stemp(MNB,MNB)
      dimension nh_i(MNS),scr(MNB,MNB)

      data nh_i/' S ',' P-',' P ',' D-',' D ',' F-',' F ',
     :                  ' G-',' G ','H-','H'/

c--------------------------------------------------------------------
c       initialize arrays and constants required for the computation.
c--------------------------------------------------------------------
      tcsqr=-two*c*c
      set = .true.
      iter=0
      en=zero
      eold=zero
      delta=zero

      do ia=1,MNB
      sa(ia)=zero
      do ib=1,MNB
      can(ib,ia)=zero
      drs(ib,ia)=zero
      p(ia,ib)=zero
      enddo
      enddo

      mbasis=2*nbasis

      do isym=1,nsym
      ndim=nbas(isym)
      ndim2=2*ndim
 
      do ia=1,ndim2
      ja=ia+2*nskip(isym)
      do ib=1,ndim2
      jb=ib+2*nskip(isym)
      stemp(ia,ib)=emt(ja,jb)
      enddo
      enddo
     
      if(jprint.ne.0)call matout(stemp,MNB,ndim2,ndim2,ndim2,4hS   )

      do i=1,ndim2
      do j=1,ndim2
      drs(i,j)=stemp(i,j)
      enddo
      enddo

      call eigen (drs,scr,ndim2,sa,ndim2,bnorm,MNB)
c-------------------------------------------------------------------
c       canonical transformation
c------------------------------------------------------------------
      do ia=1,ndim2
      do ib=1,ndim2
      if(dabs(sa(ib)).gt.TINY) then
      can(ia,ib)=drs(ia,ib)/dsqrt(dabs(sa(ib)))
      st(ib,ia)=can(ia,ib)
      else
      endif
      enddo
      enddo

      if(jprint.ne.0)call matout(st,MNB,ndim2,ndim2,ndim2,4hXT  )

      if(jprint.ne.0)
     :          call matout(df_single,MNB,mbasis,mbasis,mbasis,4hH   )

      call twoe(isym,nbasis)

      do ic=1,ndim2
      jc=ic+2*nskip(isym)
      do id=1,ndim2
      jd=id+2*nskip(isym)
      a(ic,id)=df_single(jc,jd)+df_two(jc,jd)
      enddo
      enddo

      if(jprint.ne.0)call matout(a,MNB,ndim2,ndim2,ndim2,4hFOCK )

      if(jprint.ne.0)
     :call matout(df_two,MNB,mbasis,mbasis,mbasis,4hVHF )

      call mult(a,can,b,MNB,ndim2)
      call mult(st,b,df_mat,MNB,ndim2)
      if(jprint.ne.0)
     :      call matout(df_mat,MNB,ndim2,ndim2,ndim2,6hFRIME )
      call eigen (df_mat,scr,ndim2,eng,ndim2,bnorm,MNB)
      do i=1,ndim2
      do j=1,ndim2
      drs(i,j)=df_mat(i,j)
      enddo
      enddo

      if(jprint.ne.0)
     : call matout(drs,MNB,ndim2,ndim2,ndim2,4hC   )
      call mult(can,drs,eigv,MNB,ndim2)
  
      loc = 1
      last=ndim-1
      dowhile(loc.le.last)
      samin = eng(loc)
      locus = loc
      do  ia = loc,ndim
      if (eng(ia).lt.samin) then
      samin = eng(ia)
      locus = ia
      endif
      enddo
      if (locus .gt. loc) then
      swap = eng(loc)
      eng(loc) = eng(locus)
      eng(locus) = swap
      do  ib = 1,ndim2
      swap = eigv(ib,loc)
      eigv(ib,loc) = eigv(ib,locus)
      eigv(ib,locus) = swap
      enddo
      endif
      loc = loc + 1
      enddo
c-----------------------------------------------------------------
      do ia=1,ndim2
      do ib=1,ndim2
      b(ia,ib)=zero
      enddo
      b(ia,ia)=eng(ia)
      oe(ia,isym)=eng(ia)
      if(isym.eq.2.and.ia.eq.1)print*,'energy =',oe(ia,isym)
      enddo

      if(jprint.ne.0)call matout(b,MNB,ndim2,ndim2,ndim2,4hE   )

      if(jprint.ne.0)call matout(eigv,MNB,ndim2,ndim2,ndim2,4hCP  )

      do ia=1,ndim2
      ja=ia+2*nskip(isym)
      do ib=1,ndim2
      jb=ib+2*nskip(isym)
      pold(ja,jb)=p(ja,jb)
      p(ja,jb)=0.0d0
      nocc=nocorb(isym)
      if(nocc.ne.0)then
      do ibas=1,nocc
      twojp1=two*occ(isym)+one
      p(ja,jb)=p(ja,jb)+twojp1*eigv(ia,ibas)*eigv(ib,ibas)
      enddo
      endif
      vec(ja,jb)=eigv(ia,ib)
      enddo
      enddo

      enddo

      if(jprint.ne.0)call matout(p,MNB,mbasis,mbasis,mbasis,4hP   )
      
      call integral(nbasis)
      eold=en
      en=0.d0
      do i=1,mbasis
      do j=1,mbasis
      en=en+0.5d0*p(i,j)*(two*df_single(i,j)+df_two(i,j))
      enddo
      enddo
      
      write(STDOUT,204)en    
      do jsym=1,nsym
      do ibas=1,nbas(jsym)
      write(STDOUT,206)ibas,nh_i(jsym),oe(ibas,jsym),
     :                 oe(ibas+nbas(jsym),jsym)
      enddo
      enddo

      return
 199  format(/,1x,'CONVERGENCE DATA',/,
     1        ' MAXIMUM NO. OF ITERATION=',i6,/,
     2        ' CONVERGENCE CRITERION   = 1.0D-',i2,//)
 200  format(/,1x,'CYCLE      DENSITY CONV',/)
 201  format(i3,1x,(2x,d20.7))
 202  format(///,4x,33hSCF FAILS TO CONVERGES AT CYCLE  ,i4)
 203  format(///,4x,22hSCF CONVERGES AT CYCLE, i4)
 204  format(///,4x,20hELECTRONIC ENERGY=  ,d20.11)
 205  format(///,4x,32h ORBITAL ENERGIES (+ve and -ve) ,i4,///)
 206  format(/,4x,I2,A,2D20.11)
      end            
c---------------------------------------------------------------------
      subroutine output(nbasis)
c----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      include 'mpif.h'
      common/wave/pf(MN,MNBAS),qf(MN,MNBAS),mf(MNBAS)
     : /basis/gl(MN,MNBAS),gs(MN,MNBAS),alpha(MNBAS,MNS),mg(MNBAS)
     : /dumm/tmp1(MN),tmp2(MN),tmp3(MN),tmp4(MN)
     : /quan1/kap(MNS),nak(MNOCC),np(MNOCC),nsy(MNS),nkk(MNS)
     : /grid/steplin,r(MN),rp(MN),rpor(MN),rnt,h,hp,n
     : /datt/nsym,nocc,nbas(MNS),ind
     : /vector/eigv(MNB,MNB),eng(MNB,MNS)
     :  /jv/occ(MNOCC),orbe(MNB),iq(MNB)
     :  /skip/nskip(MNS)
     :  /print/jprint,iout
     :  /kapp/kappa(MNB),kappc(MNB),iqc(MNOCC)
     : /freez/enmin(MNS),enmax(MNS)
 
      write(STDOUT,208)
c-----------------------------------------------------------------------
c            now construct the gaussian orbitals on a grid.
c-----------------------------------------------------------------------
      write(IWFN)h,n
      write(IWFN)(r(i),i = 1,n),(rp(i),i = 1,n),(rpor(i),i=1,n)
      do isym=1,nsym

c     write only the relevant orbitals i.e skip writting the
c     negavtive energy orbitals
******************************************
*   if you like to have expansion coeff
*****************************************
c         write(IWFN)nbas(isym),nskip(isym)
*****************************************
      do  ia=1,nbas(isym)
c   initialise 
               do ic=1,n
                       pf(ic,ia)=0.0d0
                       qf(ic,ia)=0.0d0
               enddo
       ka=ia+nskip(isym)
       ja=ia+2*nskip(isym)
 	do ib =1,nbas(isym)
        jb=ib+2*nskip(isym)
         kb=ib+nskip(isym)
          do  ic=1,n
            pf(ic,ia)=pf(ic,ia)+eigv(jb,ja)*gl(ic,kb)
            qf(ic,ia)=qf(ic,ia)+eigv(jb+nbas(isym),ja)*gs(ic,kb)
	  enddo
******************************************
*   if you like to have expansion coeff
c         write(IWFN)eigv(jb,ja),eigv(jb+nbas(isym),ja)
*******************************************          
        enddo
      if (isym.eq.1.and.ia.eq.1) print*,'2p- radial',
c       if (iout.eq.0) print*,'2p- radial',
     :    eng(ia,isym),pf(100,ia) 

      if (iout.eq.0) then
      write(15,*)ia,eng(ia,isym)
      write(IWFN)eng(ia,isym)
      write(IWFN)(pf(ii,ia),ii=1,n),(qf(ii,ia),ii=1,n)
      do ib =1,ia
      kb=ib+nskip(isym)
      result = zero
      do ic = 1,n
      tmp1(ic) = (pf(ic,ib)*pf(ic,ia)+qf(ic,ib)*
     :                     qf(ic,ia))*rp(ic)
      enddo
c  integrate tmp1 for orthogonality
      call quad(tmp1,n,result)
      r0=rint(ia,ib,0)
      r1=0
      r2=0
      r3=0
      if(ia.eq.ib)r1=rint(ia,ib,1)
      if(ia.eq.ib)r2=rint(ia,ib,2)
      if(ia.eq.ib)r3=rint(ia,ib,-1)
      write(STDOUT,'(2i4,2x,4(2x,D15.9))')ka,kb,r0,r1,r2,r3
      enddo
      else

      if(eng(ia,isym).le.0.0d0.and.eng(ia,isym).
     c                         lt.enmin(isym))go to 8888
      if(eng(ia,isym).gt.0.0d0.and.eng(ia,isym).
     c                         gt.enmax(isym))go to 8888
      write(15,*)ia,eng(ia,isym),enmax(isym)
      write(IWFN)eng(ia,isym)
      write(IWFN)(pf(ii,ia),ii=1,n),(qf(ii,ia),ii=1,n)
      do ib =1,ia
      kb=ib+nskip(isym)
      result = zero
      do ic = 1,n
      tmp1(ic) = (pf(ic,ib)*pf(ic,ia)+qf(ic,ib)*
     :                     qf(ic,ia))*rp(ic)
      enddo
c  integrate tmp1 for orthogonality
      call quad(tmp1,n,result)
      r0=rint(ia,ib,0)
      r1=0
      r2=0
      r3=0
      if(ia.eq.ib)r1=rint(ia,ib,1)
      if(ia.eq.ib)r2=rint(ia,ib,2)
      if(ia.eq.ib)r3=rint(ia,ib,-1)
      write(STDOUT,'(2i4,2x,4(2x,D15.9))')ka,kb,r0,r1,r2,r3
      enddo
c     endif
 8888  continue
      endif
      enddo
      enddo

c     rewind(IWFN)
c     read(IWFN)h,n
c     read(IWFN)(r(i),i = 1,n),(rp(i),i = 1,n),(rpor(i),i=1,n)
c     do i=1,nbasis
c     read(IWFN)orben
c     read(IWFN)(pf(ii,i),ii=1,n),(qf(ii,i),ii=1,n)
c     enddo

      return
 208  format(2x,'THESE ARE THE ORTHOGONALITY OF THE  ',
     :  1x,'ORBITALS',//,2x,'Sl. No.','    <A|B>      ',
     :  '    <A|R|B>     ','    <A|R*R|B>   ',
     :  '     <A|(1/R)|B>',//)
 210  format(///,2x,'THESE ARE THE DIPOLE MOMENTS OF THE  ',
     :  1x,'ORBITALS',//,2x,'Sl. No.',5x,'DIPOLE MOMENTS'/)
      end
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      real*8 function clrx (kappaa,k,kappab)
c----------------------------------------------------------------------
c   the value of clrx is the 3-j symbol:
c
c                    ( ja        k        jb  )
c                    ( 1/2       0       -1/2 )
c 
c the  k's are kappa angular quantum numbers. the formula is taken from d m
c brink and g r satchler, <angular momentum>, second edition (oxford: 
c clarendon press, 1968), p 138.   the logarithms of the first  mfact  
c factorials must be available in  common/facts/ for this program to function 
c correctly. note that  n!  is stored in fact(n+1).
c----------------------------------------------------------------------
      implicit real*8(a-h, o-z)
      include 'mpif.h'
      dimension gam(500)
      common/cons/zero,half,tenth,one,two,three,ten

       call factt(gam) 
c----------------------------------------------------------------------
c               determine the absolute values of the kappas
c----------------------------------------------------------------------
      ka =iabs (kappaa)
      kb =iabs (kappab)
c----------------------------------------------------------------------
c                   perform the triangularity check
c----------------------------------------------------------------------
      if ((iabs(ka-kb) .le. k) .and. (ka+kb-1 .ge. k)) then
c----------------------------------------------------------------------
c   triangularity satisfied; compute the 3j coefficient
c   begin with the logarithm of the square of the leading term
c----------------------------------------------------------------------
         exptrm = -dlog (dble (ka*kb))
c----------------------------------------------------------------------
c    compute the logarithm of the square root of the leading term and the 
c    factorial part that doesn't depend on the parity of ka+kb+k (the delta 
c    factor)
c----------------------------------------------------------------------
         kapkb = ka+kb
         kabkp = kapkb+k
         kamkb = ka-kb
         kbmka = kb-ka
         exptrm = half*(exptrm+gam(kapkb-k  )+gam(kamkb+k+1)
     :                        +gam(kbmka+k+1)-gam(kabkp  +1) )
c----------------------------------------------------------------------
c   the remainder depends on the parity of ka+kb+k
c----------------------------------------------------------------------
         if (mod (kabkp,2) .eq. 0) then
c----------------------------------------------------------------------
c   computation for even parity case
c   include the phase factor: a minus sign if necessary
c----------------------------------------------------------------------
            if (mod (3*kabkp/2,2) .eq. 0) then
               clrx =  one
            else
               clrx = -one
            endif
c----------------------------------------------------------------------
c   include the contribution from the factorials
c----------------------------------------------------------------------
            exptrm = exptrm+gam((kabkp  +2)/2)-gam((kapkb-k  )/2)
     :                     -gam((kamkb+k+2)/2)-gam((kbmka+k+2)/2)
*
         else
c----------------------------------------------------------------------
c                   computation for odd parity case
c         include the phase factor: a minus sign if necessary
c----------------------------------------------------------------------
            if (mod ((3*kabkp-1)/2,2) .eq. 0) then
               clrx =  one
            else
               clrx = -one
            endif
c----------------------------------------------------------------------
c   include the contribution from the factorials
c----------------------------------------------------------------------
            exptrm = exptrm+gam((kabkp  +1)/2)-gam((kapkb-k+1)/2)
     :                     -gam((kamkb+k+1)/2)-gam((kbmka+k+1)/2)
*
         endif
c----------------------------------------------------------------------
c   final assembly
c----------------------------------------------------------------------
         clrx = clrx*dexp (exptrm)
*
      else
c----------------------------------------------------------------------
c   triangularity violated; set the coefficient to zero
c----------------------------------------------------------------------
         clrx = zero
*
      endif
*
      return
      end
      subroutine setqic
c--------------------------------------------------------------------------
c this  subroutine sets up the coefficients for subroutines 
c quad, rinti, yzk 
c--------------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      include 'mpif.h'
      dimension b13(13,13)
      common/cons/zero,half,tenth,one,two,three,ten
     :  /grid/steplin,r(MN),rp(MN),rpor(MN),rnt,h,hp,n
     :  /lic13/a(13,13)
     :  /ncc/c1,c2,c3,c4,cnc5c(5,5),c5num(5,5),pi,c
c------------------------------------------------------------------------
c   five-point newton-cotes coefficients for closed integration. 
c   expressed as  rational numbers
c-------------------------------------------------------------------------
      data (c5num(i,2),i = 1,5)/ 251.0d 00, 646.0d 00,
     :       -264.0d 00, 106.0d 00, -19.0d 00/
      data (c5num(i,3),i = 1,5)/ 232.0d 00, 992.0d 00,
     :        192.0d 00,  32.0d 00,  -8.0d 00/
      data (c5num(i,4),i = 1,5)/ 243.0d 00, 918.0d 00,
     :        648.0d 00, 378.0d 00, -27.0d 00/
      data (c5num(i,5),i = 1,5)/ 224.0d 00,1024.0d 00,
     :        384.0d 00,1024.0d 00, 224.0d 00/
      data c5den/ 720.0d 00/
c------------------------------------------------------------------------
*----------------------------------------------------------------------*
*                                                                      *
*   Six-point Newton-Cotes coefficients for closed integration.  Ex-   *
*   pressed as rational numbers                                        *
*
c     data (c6num(i,2),i = 1,6)/ 475.0d 00,1427.0d 00,
c    :                          -798.0d 00, 482.0d 00,
c    :                          -173.0d 00,  27.0d 00/
c     data (c6num(i,3),i = 1,6)/ 448.0d 00,2064.0d 00,
c    :                           224.0d 00, 224.0d 00,
c    :                           -96.0d 00,  16.0d 00/
c     data (c6num(i,4),i = 1,6)/ 459.0d 00,1971.0d 00,
c    :                          1026.0d 00,1026.0d 00,
c    :                          -189.0d 00,  27.0d 00/
c     data (c6num(i,5),i = 1,6)/ 448.0d 00,2048.0d 00,
c    :                           768.0d 00,2048.0d 00,
c    :                           448.0d 00,   0.0d 00/
c     data (c6num(i,6),i = 1,6)/ 475.0d 00,1875.0d 00,
c    :                          1250.0d 00,1250.0d 00,
c    :                          1875.0d 00, 475.0d 00/
*
c     data c6den/1440.0d 00/
*
c--------------------------------------------------------------------------
      data (b13( 1,i),i = 1,13) /            -1486442880.0d 00,
     :            5748019200.0d 00, -15807052800.0d 00,
     :           35126784000.0d 00, -59276448000.0d 00,
     :           75873853440.0d 00, -73766246400.0d 00,
     :           54195609600.0d 00, -29638224000.0d 00,
     :           11708928000.0d 00,  -3161410560.0d 00,
     :             522547200.0d 00,    -39916800.0d 00/
      data (b13( 2,i),i = 1,13) /              -39916800.0d 00,
     :            -967524480.0d 00,   2634508800.0d 00,
     :           -4390848000.0d 00,   6586272000.0d 00,
     :           -7903526400.0d 00,   7376624640.0d 00,
     :           -5269017600.0d 00,   2822688000.0d 00,
     :           -1097712000.0d 00,    292723200.0d 00,
     :             -47900160.0d 00,      3628800.0d 00/
      data (b13( 3,i),i = 1,13) /                3628800.0d 00,
     :             -87091200.0d 00,   -684478080.0d 00,
     :            1596672000.0d 00,  -1796256000.0d 00,
     :            1916006400.0d 00,  -1676505600.0d 00,
     :            1149603840.0d 00,   -598752000.0d 00,
     :             228096000.0d 00,    -59875200.0d 00,
     :               9676800.0d 00,      -725760.0d 00/
      data (b13( 4,i),i = 1,13) /                -725760.0d 00,
     :              13063680.0d 00,   -143700480.0d 00,
     :            -476910720.0d 00,   1077753600.0d 00,
     :            -862202880.0d 00,    670602240.0d 00,
     :            -431101440.0d 00,    215550720.0d 00,
     :             -79833600.0d 00,     20528640.0d 00,
     :              -3265920.0d 00,       241920.0d 00/
      data (b13( 5,i),i = 1,13) /                 241920.0d 00,
     :              -3870720.0d 00,     31933440.0d 00,
     :            -212889600.0d 00,   -303937920.0d 00,
     :             766402560.0d 00,   -447068160.0d 00,
     :             255467520.0d 00,   -119750400.0d 00,
     :              42577920.0d 00,    -10644480.0d 00,
     :               1658880.0d 00,      -120960.0d 00/
      data (b13( 6,i),i = 1,13) /                -120960.0d 00,
     :               1814400.0d 00,    -13305600.0d 00,
     :              66528000.0d 00,   -299376000.0d 00,
     :            -148262400.0d 00,    558835200.0d 00,
     :            -239500800.0d 00,     99792000.0d 00,
     :             -33264000.0d 00,      7983360.0d 00,
     :              -1209600.0d 00,        86400.0d 00/
      data (b13( 7,i),i = 1,13) /                  86400.0d 00,
     :              -1244160.0d 00,      8553600.0d 00,
     :             -38016000.0d 00,    128304000.0d 00,
     :            -410572800.0d 00,            0.0d 00,
     :             410572800.0d 00,   -128304000.0d 00,
     :              38016000.0d 00,     -8553600.0d 00,
     :               1244160.0d 00,       -86400.0d 00/
      data (b13( 8,i),i = 1,13) /                 -86400.0d 00,
     :               1209600.0d 00,     -7983360.0d 00,
     :              33264000.0d 00,    -99792000.0d 00,
     :             239500800.0d 00,   -558835200.0d 00,
     :             148262400.0d 00,    299376000.0d 00,
     :             -66528000.0d 00,     13305600.0d 00,
     :              -1814400.0d 00,       120960.0d 00/
      data (b13( 9,i),i = 1,13) /                 120960.0d 00,
     :              -1658880.0d 00,     10644480.0d 00,
     :             -42577920.0d 00,    119750400.0d 00,
     :            -255467520.0d 00,    447068160.0d 00,
     :            -766402560.0d 00,    303937920.0d 00,
     :             212889600.0d 00,    -31933440.0d 00,
     :               3870720.0d 00,      -241920.0d 00/
      data (b13(10,i),i = 1,13) /                -241920.0d 00,
     :               3265920.0d 00,    -20528640.0d 00,
     :              79833600.0d 00,   -215550720.0d 00,
     :             431101440.0d 00,   -670602240.0d 00,
     :             862202880.0d 00,  -1077753600.0d 00,
     :             476910720.0d 00,    143700480.0d 00,
     :             -13063680.0d 00,       725760.0d 00/
      data (b13(11,i),i = 1,13) /                 725760.0d 00,
     :              -9676800.0d 00,     59875200.0d 00,
     :            -228096000.0d 00,    598752000.0d 00,
     :           -1149603840.0d 00,   1676505600.0d 00,
     :           -1916006400.0d 00,   1796256000.0d 00,
     :           -1596672000.0d 00,    684478080.0d 00,
     :              87091200.0d 00,     -3628800.0d 00/
      data (b13(12,i),i = 1,13) /               -3628800.0d 00,
     :              47900160.0d 00,   -292723200.0d 00,
     :            1097712000.0d 00,  -2822688000.0d 00,
     :            5269017600.0d 00,  -7376624640.0d 00,
     :            7903526400.0d 00,  -6586272000.0d 00,
     :            4390848000.0d 00,  -2634508800.0d 00,
     :             967524480.0d 00,     39916800.0d 00/
      data (b13(13,i),i = 1,13) /               39916800.0d 00,
     :            -522547200.0d 00,   3161410560.0d 00,
     :          -11708928000.0d 00,  29638224000.0d 00,
     :          -54195609600.0d 00,  73766246400.0d 00,
     :          -75873853440.0d 00,  59276448000.0d 00,
     :          -35126784000.0d 00,  15807052800.0d 00,
     :           -5748019200.0d 00,   1486442880.0d 00/
*
      data b13den/479001600.0d 00/
*
c------------------------------------------------------------------------
*   THIRTEEN-POINT COEFFICIENTS FOR DPBDT
*
         factor = 1/b13den
         do j = 1,13
            do i = 1,13
               a(i,j) = b13(i,j)*factor
c            print*,'i,j,a(i,j)',i,j,a(i,j)
           enddo
          enddo

c------------------------------------------------------------------------
c                  newton-cotes coefficients for yzk and quad
c-----------------------------------------------------------------------
      factor = h/c5den
      do j = 2,4
         do i = 2,5
            cnc5c(i,j) = factor*c5num(i,j)
         enddo
       enddo
*
      c1 = factor*c5num(1,5)
      c2 = factor*c5num(2,5)
      c3 = factor*c5num(3,5)
      c4 = c1+c1
*
      do j = 2,5
         do i = 2,5
            c5num(i,j) = c5num(i,j)/c5den
         enddo
      enddo
      pi = 3.141592653589790d 000
      c = 1.3703598950000d 002
c     c = 1.00000d006
      zero=0.0d 00
      half=0.5d 00 
      tenth=0.1d 00
      one=1.0d 00
      two=2.0d 00
      three=3.0d 00
      ten=10.0d 00
      return
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer function idbf(n,i)
c     calculates double factorial
c-----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      include 'mpif.h'
      if(n.lt.-1) then
      write(*,*)'n.lt.-1',n
      stop
      else
        if((n.eq.0).or.(n.eq.-1)) then
          idbf=1
        else
          idbf=1
          do j=n,1,-i
            idbf=j*idbf
          enddo
        endif
      endif
      return
      end
c----------------------------------------------------------------------
      subroutine nucpot
c----------------------------------------------------------------------
c evaluate the nuclear potential for point and fermi models.
c subroutines called: es.
c-----------------------------------------------------------------------
      implicit doubleprecision (a-h, o-z)
      include 'mpif.h'
      logical set
*
      common /def1/amass,z
     :  /npar/parm(2),nparm
     :  /npot/zz(MN)
     :  /cons/zero,half,tenth,one,two,three,ten
     :  /grid/steplin,r(MN),rp(MN),rpor(MN),rnt,h,hp,n

      data pi/3.141592653589793d0/

c-----------------------------------------------------------------------
c                                 point nucleus
c-----------------------------------------------------------------------
      if (nparm.eq.0) then
         do 1 i = 1,MN
            zz(i) = z
    1    continue
c-----------------------------------------------------------------------
c                               fermi distribution
c-----------------------------------------------------------------------
      elseif (nparm.eq.2) then
         c = parm(1)
         a = parm(2)
         abc = a/c
         tabc = two*abc
         abc2 = abc*abc
         thabc2 = three*abc2
         abc3 = abc2*abc
         cba = c/a
         pi2 = pi*pi
         hpiac2 = half*pi2*abc2
         six = two*three
         h3 = half*three
         h3php = h3+hpiac2
         call es (-cba,s2mcba,s3mcba)
         sabc3 = six*abc3
         dmsas = -sabc3*s3mcba
         en = one + abc2*pi2 + dmsas
         zbn = z/en
*
         do 2 i = 1,MN
            ri = r(i)
            rmc = ri-c
            rmcba = rmc/a
            rbc = ri/c
            if (rbc.le.one) then
               call es (rmcba,s2rcba,s3rcba)
               zz(i) = zbn*( dmsas + sabc3*s3rcba
     :                      +rbc*( h3php-thabc2*s2rcba
     :                            -half*rbc*rbc) )
            else
               if (.not. set) then
                  nnuc = i
                  set = .true.
             endif
               call es (-rmcba,s2rcba,s3rcba)
               zz(i) = z * ( one
     :                        +thabc2 * ( rbc *s2rcba
     :                                   +tabc*s3rcba ) / en )
            endif
    2    continue
      endif
      return
      end
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine es (f,s2f,s3f)
c------------------------------------------------------------------------
c                        evaluate the sum of the series
c
c          k(f)=sum_(0)^(infinity)(-1)^n  exp(n*f)/ (n^k)
c
c       for k = 2, 3 to machine precision. this is a utility 
c------------------------------------------------------------------------
      implicit doubleprecision (a-h, o-z)
      include 'mpif.h'
      common /cons/zero,half,tenth,one,two,three,ten
*
      n = 0
      s2f = zero
      s3f = zero
      fase = one
    1 n = n+1
         en = dble (n)
         obn = one/en
         fase = -fase
         enf = dexp (en*f)
         term2 = fase*enf*obn*obn
         term3 = term2*obn
         s2last = s2f
         s2f = s2f+term2
         s3f = s3f+term3
         if (dabs (s2f) .ne. dabs (s2last)) goto 1
      return
*
      end
c------------------------------------------------------------
      subroutine mult(a,b,c,im,m)
c    matrice multiplication
c------------------------------------------------------------
      implicit double precision (a-h, o-z)
      include 'mpif.h'
      dimension a(im,im),b(im,im),c(im,im)
      do i=1,m
          do j=1,m
          c(i,j)=0.0d0
            do k=1,m
            c(i,j)=c(i,j)+a(i,k)*b(k,j)
            end do
          end do
       end do
       return
       end
c------------------------------------------------------------
      subroutine fzero
c    initialise
c------------------------------------------------------------
      implicit real*8(a-h,o-z)
      include 'mpif.h'
      common /grid/steplin,r(MN),rp(MN),rpor(MN),rnt,h,hp,n
     :  /wave/pf(MN,MNBAS),qf(MN,MNBAS),mf(MNBAS)
     :  /basis/gl(MN,MNBAS),gs(MN,MNBAS),alpha(MNBAS,MNS),mg(MNBAS)
     :  /ncc/c1,c2,c3,c4,cnc5c(5,5),c5num(5,5),pi,c
     :  /quan1/kap(MNS),nak(MNOCC),np(MNOCC),nsy(MNS),nkk(MNS)
     :  /cons/zero,half,tenth,one,two,three,ten
     :  /dumm/tmp1(MN),tmp2(MN),tmp3(MN),tmp4(MN)
     :  /adat/occj(MNOCC),occl(MNOCC),unocl(MNS),unocj(MNS)
     :  /overlap/sss(MNB,MNB),sll(MNB,MNB)
     :  /datt/nsym,nocc,nbas(MNS),ind
     :  /single/df_single(MNB,MNB)
     :  /npar/parm(2),nparm
      common /def1/amass,z
     :  /npot/zz(MN)
      do i=1,MN
                 r(i)=zero
                 rp(i)=zero
                 rpor(i)=zero
                 tmp1(i)=zero
                 tmp2(i)=zero
                 tmp3(i)=zero
                 tmp4(i)=zero
                 zz(i)=zero
           do j=1,MNBAS
                 pf(i,j)=zero
                 qf(i,j)=zero
                 gl(i,j)=zero
                 gs(i,j)=zero
           enddo
      enddo
      do i=1,MNS
                 nsy(i)=0
                 unocl(i)=zero
                 unocj(i)=zero
         do j=1,MNBAS
                 alpha(j,i)=zero
         enddo
      enddo
      do i=1,MNB
         do j=1,MNB
                 sss(i,j)=zero
                 sll(i,j)=zero
                 df_single(i,j)=zero
         enddo
      enddo
      return
      end
c***********************************************************************
      subroutine slater(ia,ib,ic,id,k,rkll,rkss,rkls,rksl)

*   the value of this  function is the  slater integral  as normally   *
*   defined in terms of the four sets of quantum numbers  a,b,c,d.     *
*                                                                      *
*   subroutines called: quad, yzk.                                     *
*                                                                      *
c***********************************************************************
*
      implicit real*8(a-h, o-z)
      include 'mpif.h'
c     dimension rhop(MN),rttk(MN),wk(MN),temp(MN),zk(MN)
      common /cons/zero,half,tenth,one,two,three,ten
     :  /grid/steplin,r(MN),rp(MN),rpor(MN),rnt,h,hp,n
     :  /ncc/c1,c2,c3,c4,cnc5c(5,5),c5num(5,5),pi,c
     :  /dumm/tmp1(MN),tmp2(MN),tmp3(MN),tmp4(MN)
     :  /wave/pf(MN,MNBAS),qf(MN,MNBAS),mf(MNBAS)
     :  /basis/gl(MN,MNBAS),gs(MN,MNBAS),alpha(MNBAS,MNS),mg(MNBAS)

      mtp=n
      do i=1,mtp
         tmp1(i)=0
         tmp2(i)=0
         tmp3(i)=0
         tmp4(i)=0
      enddo

c     call yzk(k,ib,id)
      do i = 2,mtp
         tmp3(i)=gl(i,ia)*gl(i,ic)*rpor(i)*tmp1(i)
      enddo
      call quad (tmp3,mtp,rkll)

      do i = 2,mtp
         tmp3(i)=gs(i,ia)*gs(i,ic)*rpor(i)*tmp2(i)
      enddo

      call quad (tmp3,mtp,rkss)

      do i = 2,mtp
         tmp3(i)=gs(i,ia)*gs(i,ic)*rpor(i)*tmp1(i)
      enddo

      call quad (tmp3,mtp,rksl)

      do i = 2,mtp
         tmp4(i)=gl(i,ia)*gl(i,ic)*rpor(i)*tmp2(i)
      enddo

      call quad (tmp4,mtp,rkls)

      return
      end

      subroutine matout(a,im,in,m,n,label)
      implicit double precision (a-h,o-z)
      include 'mpif.h'
      character*6 label
      dimension a(im,im)
      ihigh=0
  10  low=ihigh+1
      ihigh=ihigh+5
      ihigh=min0(ihigh,n)
              write(STDOUT,20)label,(i,i=low,ihigh)
  20  format(///,3x,5h THE ,a6,6h ARRAY,/15x,5(10x,i3,6x)//)
      do i=1,m
              write(STDOUT,40) i,(a(i,j),j=low,ihigh)
       print*,a(i,j),i
      enddo
  40  format(i8,5x,8(1x,d12.6))
      if(n-ihigh)50,50,10
  50  return
      end

      subroutine scfout
      write(STDOUT,9990)
      write(STDOUT,9999)
      write(STDOUT,9998)
      write(STDOUT,9998)
      write(STDOUT,9998)
      write(STDOUT,9997)
      write(STDOUT,9998)
      write(STDOUT,9996)
      write(STDOUT,9998)
      write(STDOUT,9995)
      write(STDOUT,9998)
      write(STDOUT,9998)
      write(STDOUT,9998)
      write(STDOUT,9999)

9999  format(20x,'**************************************')
9998  format(20x,'*                                    *')
9997  format(20x,'*        Relativistic SCF Run        *')
9996  format(20x,'*            written by              *')
9995  format(20x,'*   P. K.  Panda and R. K. Chaudhuri *')
9990  format(//)
      return
      end
c*****************************************************************
c     generate the basis functions on a grid point. It generates
c     (at present) geometrical basis but can be easily transformed
c     to a general one and also to contracted gaussian type.
c*****************************************************************
      subroutine gbasis(nbasis,method)
      implicit real*8(a-h,o-z)
      include 'mpif.h'
      real*8 norm_large,norm_small
      character*2 method
      common /grid/steplin,r(MN),rp(MN),rpor(MN),rnt,h,hp,n
     :  /wave/pf(MN,MNBAS),qf(MN,MNBAS),mf(MNBAS)
     :  /basis/gl(MN,MNBAS),gs(MN,MNBAS),alpha(MNBAS,MNS),mg(MNBAS)
     :  /ncc/c1,c2,c3,c4,cnc5c(5,5),c5num(5,5),pi,c
     :  /quan1/kap(MNS),nak(MNOCC),np(MNOCC),nsy(MNS),nkk(MNS)
     :  /cons/zero,half,tenth,one,two,three,ten
     :  /dumm/tmp1(MN),tmp2(MN),tmp3(MN),tmp4(MN)
     :  /adat/occj(MNOCC),occl(MNOCC),unocl(MNS),unocj(MNS)
     :  /overlap/sss(MNB,MNB),sll(MNB,MNB)
     :  /datt/nsym,nocc,nbas(MNS),ind
     :  /single/df_single(MNB,MNB)
     :  /npar/parm(2),nparm
      common /def1/amass,z
     :  /npot/zz(MN)
     :  /jvalue/orbj(MNOCC,MNS),iocc(MNBAS)
     :  /skip/nskip(MNS)
     :  /kapp/kappa(MNB),kappc(MNB),iqc(MNOCC)
     :  /info1/occk(MNOCC,MNS)
     :  /jv/occ(MNOCC),orbe(MNB),iq(MNB)
     :  /info/nocorb(MNOCC),nucorb(MNB)
     :  /info2/qqf(MN,MNBAS),qql(MN,MNBAS),qqs(MN,MNBAS)
     :  /UB/alpha1,beta1
     :  /WT/alpha2,beta2,gamma,delta
     :  /ET/alpha3(MNS),beta3(MNS)
     :  /tab/tgl(MN),tgs(MN)

      dimension rnk(MN)
     :    ,dtgl(MN,MNBAS),dtgs(MN,MNBAS)
      spi=dsqrt(pi)
      four = two*two
      onehalf = one + half
      nbasis=0
      nskip(1)=0
      do isym=1,nsym
      nbasis=nbasis+nbas(isym)
      nskip(isym+1)=nskip(isym)+nbas(isym)
      enddo
  
      ii=0
      do isym=1,nsym
      do ibas=1,nbas(isym)
      ii=ii+1
      kappa(ii)=nak(isym)
      orbe(ii)=occ(isym)
      if(kappa(ii).eq.-1)iq(ii)=1
      if(kappa(ii).eq. 1)iq(ii)=-1
      if(kappa(ii).eq.-2)iq(ii)= 1
      if(kappa(ii).eq. 2)iq(ii)=-1
      if(kappa(ii).eq.-3)iq(ii)= 1
      if(kappa(ii).eq. 3)iq(ii)=-1
      if(kappa(ii).eq.-4)iq(ii)= 1
      if(kappa(ii).eq. 4)iq(ii)=-1
      if(kappa(ii).eq.-5)iq(ii)= 1
      if(kappa(ii).eq. 5)iq(ii)=-1
      if(kappa(ii).eq.-6)iq(ii)= 1
      enddo
      enddo
      
      ii=0
      do isym=1,nsym
      nn=nocorb(isym)
      if(nn.ne.0)then
      do ibas=1,nocorb(isym)
      ii=ii+1
      occl(ii)=occ(isym)
      occk(ibas,isym)=occ(isym)
      kappc(ii)=nak(isym)
      if(kappc(ii).eq.-1)iqc(ii)=1
      if(kappc(ii).eq. 1)iqc(ii)=-1
      if(kappc(ii).eq.-2)iqc(ii)= 1
      if(kappc(ii).eq. 2)iqc(ii)=-1
      if(kappc(ii).eq.-3)iqc(ii)= 1
      if(kappc(ii).eq. 3)iqc(ii)=-1
      if(kappc(ii).eq.-4)iqc(ii)= 1
      if(kappc(ii).eq. 4)iqc(ii)=-1
      if(kappc(ii).eq.-5)iqc(ii)= 1
      if(kappc(ii).eq. 5)iqc(ii)=-1
      if(kappc(ii).eq.-6)iqc(ii)= 1
      enddo
      endif
      enddo

      nocc=ii
             write(STDOUT,200)nbasis
             write(STDOUT,201)
      do isym=1,nsym
             write(STDOUT,202)isym,nbas(isym),nskip(isym)
      enddo

      do i=1,nbasis
            mf(i)=n
             mg(i)=n
      enddo

      if(method.eq.'UB')then
      iopt=1
      endif
      if(method.eq.'WT')then
      iopt=2
      endif
      if(method.eq.'ET')then
      iopt=3
      endif
 
      sq2=dsqrt(two)
      call init(1)
c-----------------------------------------------------------------
c              generate the gaussians on a grid
c         one thing to be noted: these are normalized  guassians. 
c-----------------------------------------------------------------
c  here we give different alpha3 for different symmetries. So alpha3
c  is also given thro  the common block instead thro arguments.
      do isym=1,nsym

c----------------------------------------------------------------
c  here the value of the exponential factor for the 
c                 gaussians are computed.
c----------------------------------------------------------------

      do ib=1,n
      rnk(ib) = r(ib)**nkk(isym)
      enddo
      do ic=1,nbas(isym)
      jc=ic+nskip(isym)

C here the alpha is calculated using alpha(i) and beta(i) for the 
C three different types of basis given by UB,WT and ET.

      goto (1,2,3),iopt

   1  alpha(jc,isym) = alpha1*(beta1**(ic-1))
      goto 4

   2  alpha(jc,isym) = alpha2*(beta2**(ic-1))*(1+ gamma *
     :               (ic/nbas(isym))**delta)
      goto 4

   3  alpha(jc,isym) = alpha3(isym)*(beta3(isym)**(ic-1))
   4  continue
      enddo

      dkap=dble(kap(isym))
      dnkk=dble(nkk(isym))
      ifact=idbf(2*nkk(isym)-1,2)
      power = (two*dnkk+one)/two
      do id = 1,nbas(isym)
      jd=id+nskip(isym)
      cl=two**(two*dnkk+onehalf)*(alpha(jd,isym)**power)/(ifact*spi)
      cs =dsqrt((two*dnkk-one)/(alpha(jd,isym)*
     :           (four*(dkap*dkap+dkap+dnkk)-one)))
      norm_large=dsqrt(cl/two)
      norm_small=norm_large*cs
      gl(1,jd)=0
      gs(1,jd)=0
      do ie = 2,n
      expon= dexp(-alpha(jd,isym)*r(ie)*r(ie))*rnk(ie)
      gl(ie,jd) = norm_large*expon
      gs(ie,jd) = norm_small*((dnkk+dkap)/r(ie)-
     :                   two*alpha(jd,isym)*r(ie))*expon
      qqf(ie,jd) = norm_large*((dnkk+dkap)/r(ie)-
     :                   two*alpha(jd,isym)*r(ie))*expon
      qql(ie,jd) = norm_large*((dnkk)/r(ie)-
     :                   two*alpha(jd,isym)*r(ie))*expon
c      if (jd.eq.1.and.ie.eq.100)then
c       write(18,*)norm_large, expon, r(ie),dnkk,alpha(jd,isym)
c      endif
      qqs(ie,jd) = norm_small*((dnkk+dkap)/r(ie)*(two*alpha(jd,isym)*
     :   r(ie)+(dnkk-1)/r(ie))+two*alpha(jd,isym)*(two*alpha(jd,isym)*
     :              r(ie)*r(ie)-(dnkk+1)))*expon
      
c       enddo
c        call diff2 (jd)
c       print*,'----------jd----------',jd
c       if(ie.eq.100)then
c       print*,'----------gbasis----------',jd
c       print*,'gl,gs',gl(ie,jd),gs(ie,jd)
c       endif
c       print*,'qql,qqs',qql(100,jd),qqs(100,jd)
c       print*,'tgl,tgs',tgl(100),tgs(100)
c       
c       print*,'**************************************'
      enddo
      enddo
      enddo

 200  format(/,2x,31hTOTAL NUMBER OF BASIS FUNCTION=i4,//)
 201  format(2x,'SYMMETRY   NO. OF BASIS FUNCTION   NSKIP',/)
 202  format(2x,i4,12x,i4,12x,i4)
      return
      end
c------------------------------------------------------------------
c     read input data
      subroutine readinp(mocc,method)
c------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      include 'mpif.h'
      character*80 title
      character*14 SCFIN,OUTPUT
      character*2 atom,method
      common /datt/nsym,nocc,nbas(MNS),ind
     :  /quan1/kap(MNS),nak(MNOCC),np(MNOCC),nsy(MNS),nkk(MNS)
     :  /cons/zero,half,tenth,one,two,three,ten
     :  /adat/occj(MNOCC),occl(MNOCC),unocl(MNS),unocj(MNS)
      common /def1/amass,z
     :  /jvalue/orbj(MNOCC,MNS),iocc(MNBAS)
     :  /conv/crit,maxit,npower
     :  /jv/occ(MNOCC),orbe(MNB),iq(MNB)
     :  /print/jprint,iout
     :  /info/nocorb(MNOCC),nucorb(MNB)
     :  /nucop/nucopt,isms
     :  /UB/alpha1,beta1
     :  /WT/alpha2,beta2,gamma,delta
     :  /ET/alpha3(MNS),beta3(MNS)
     :  /freez/enmin(MNS),enmax(MNS)

      open(unit=IPUNCH,status='old',file='inputfile')
      read(IPUNCH,'(a)')SCFIN
      read(IPUNCH,'(a)')OUTPUT
      open (unit=STDIN,file=SCFIN,status='old')
      open (unit=STDOUT,file=OUTPUT,status='unknown')
      read(STDIN,'(a)')title
      read(STDIN,'(a)')atom
      read(STDIN,*)nsym,(nbas(i),i=1,nsym)
      read(STDIN,*)(nocorb(i),i=1,nsym)
      read(STDIN,*)(nak(i),i=1,nsym)
      read(STDIN,*)(occ(i),i=1,nsym)
      read(STDIN,'(a)')method

      if(method.eq.'UB')then
      iopt=1
      endif
      if(method.eq.'WT')then
      iopt=2
      endif
      if(method.eq.'ET')then
      iopt=3
      endif
      print*,method
      goto  (1, 2, 3),iopt
  1   read(STDIN,*)alpha1,beta1
      goto 4
  2   read(STDIN,*)alpha2,beta2,gamma,delta 
      goto 4
  3   read(STDIN,*)(alpha3(i),i=1,nsym)
      read(STDIN,*)(beta3(i),i=1,nsym)
  4   continue
      read(STDIN,*)icharg
      read(STDIN,*)maxit,npower
*  iout=0 is output without energy limit, otherwise with limit
      read(STDIN,*)jprint,iout
      read(STDIN,*)nucopt
      read(STDIN,*)amass,z
c     read(STDIN,*)isms
      isms=0

      if (iout.ne.0) then
      read(STDIN,*)(enmin(i),i=1,nsym)
      read(STDIN,*)(enmax(i),i=1,nsym)
      endif

      crit=(0.1d0)**npower
      iele=0
      mocc=0
      do isym=1,nsym
      nn=nocorb(isym)
      if(nn.ne.0)then
      do ibas=1,nocorb(isym)
      iele=iele+2*occ(isym)+1
      mocc=mocc+1
      enddo
      endif
      enddo
      call scfout
      write(STDOUT,9994)title
      write(STDOUT,301)
      call nucleus(atom,nnuc)
      write(STDOUT,302)nsym
      write(STDOUT,*)'BASIS   --> ',(nbas(i),i=1,nsym)
      write(STDOUT,*)'OCC ORB --> ',(nocorb(i),i=1,nsym)
      write(STDOUT,*)'KAPPA   --> ',(nak(i),i=1,nsym)
      write(STDOUT,*)'J VALUE --> ',(occ(i),i=1,nsym)
      write(STDOUT,9995)method
      goto (5,6,7),iopt
  5   write(STDOUT,310)alpha1,beta1
      go to 8
  6   write(STDOUT,311)alpha2,beta2,gamma,delta
      go to 8
  7   write(STDOUT,*)'ALPHA3  --> ',(alpha3(i),i=1,nsym)
      write(STDOUT,*)'BETA3   --> ',(beta3(i),i=1,nsym)
  8   continue
      write(STDOUT,303)icharg
      write(STDOUT,304)maxit, crit
      write(STDOUT,305)jprint
      write(STDOUT,306)nucopt
      write(STDOUT,309)
  301 format('>>>>>>>>>>>>>>> INPUT DATA >>>>>>>>>>>>>>>>>>>>>>>>')
  302 format('NUMBER OF SYMMETRY = ',i2)
  304 format('MAX ITER = ',i2,8x,'CONV. =',d10.4)
  303 format('CHARGE =',i2)
  305 format('JPRINT =',i2)
  306 format('NUCLEAR OPTIONAL = ',i2)
  309 format('<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<',/)
  310 format('ALPHA1 = ',d10.4,8x,'BETA1 =',d10.4)
  311 format('ALPHA2 = ',d10.4,3x,'BETA2 =',d10.4,3x,'GAMMA = ',
     :   d10.4,3x,'DELTA =',d10.4)
9994  format(//,27x,a,//)
9995  format(//,25x,a,//)
      if(iele+icharg.ne.nnuc)then
      write(STDOUT,*)'>>>>>>>>> PROGRAMME TERMINATES DUE TO IO ERROR'
      write(STDOUT,*)' CHECK ATOMIC NUMBER, CHARGE AND OCCUPANCY'
      stop
      else
      z=dble(nnuc)
      endif
      open (unit=IWFN,file='wfn.dat',form='unformatted',
     :          status='unknown')
c     open (unit=IWFN,file='wfn.dat',status='old')

      return
      end
************************************************************************
*                                                                      *
      subroutine diff1(J)
*                                                                      *
*   ----------------   SECTION 10   SUBPROGRAM 08   ----------------   *
*                                                                      *
*   THIS SUBROUTINE COMPUTES H TIMES THE DERIVATIVE, WITH RESPECT TO   *
*   THE INTERNAL GRID, OF THE LARGE AND SMALL COMPONENTS OF THE WAVE   *
*   FUNCTION WITH INDEX  J .  THESE  ARE TABULATED, RESPECTIVELY, IN   *
*   ARRAYS  TA  AND  TB  IN  COMMON  BLOCK  /TATB/ .                   *
*                                                                      *
*   A  THIRTEEN-POINT  LAGRANGE FORMAULA IS USED FOR THE CALCULATION   *
*   OF DERIVATIVES.                                                    *
*                                                                      *
*   NO SUBROUTINES CALLED.                                             *
*                                                                      *
*   WRITTEN BY FARID A PARPIA, AT OXFORD   LAST UPDATED: 25 MAR 1988   *
*                                                                      *
************************************************************************
*
      implicit real*8 (a-h, o-z)
      include 'mpif.h'
*
      COMMON/cons/zero,half,tenth,one,two,three,ten
     :      /grid/steplin,r(MN),rp(MN),rpor(MN),rnt,h,hp,n
     :      /lic13/a(13,13)
     :      /tatb/tpf(MN),tqf(MN)
     :      /wave/pf(MN,MNBAS),qf(MN,MNBAS),mf(MNBAS)
*
      equivalence (a1,a(7,1)),(a2,a(7,2)),(a3,a(7,3)),
     :            (a4,a(7,4)),(a5,a(7,5)),(a6,a(7,6))
*
*   COMPUTE DERIVATIVE IN THREE SEPARATE REGIONS
*
*   FIRST, POINTS 1 TO 6
*     
      do i = 1,6
         hdpbdt = 0.0d 00
         hdqbdt = 0.0d 00
         do k = 1,13
            aik = a(i,k)
            hdpbdt = hdpbdt+aik*pf(k,J)
            hdqbdt = hdqbdt+aik*qf(k,J)
         enddo
         tpf(i) = hdpbdt
         tqf(i) = hdqbdt
      enddo
*
*   NEXT, POINTS 7 TO N-6
*
*   SPECIAL TREATMENT FOR THIS REGION BECAUSE OF THE SYMMETRY OF
*   THE DIFFERENTIATION FORMULA
*
      do i = 7,n-6
         tpf(i) =  a1*(pf(i-6,J)-pf(i+6,J))+a2*(pf(i-5,J)-pf(i+5,J))
     :           +a3*(pf(i-4,J)-pf(i+4,J))+a4*(pf(i-3,J)-pf(i+3,J))
     :           +a5*(pf(i-2,J)-pf(i+2,J))+a6*(pf(i-1,J)-pf(i+1,J))
         tqf(i) =  a1*(qf(i-6,J)-qf(i+6,J))+a2*(qf(i-5,J)-qf(i+5,J))
     :           +a3*(qf(i-4,J)-qf(i+4,J))+a4*(qf(i-3,J)-qf(i+3,J))
     :           +a5*(qf(i-2,J)-qf(i+2,J))+a6*(qf(i-1,J)-qf(i+1,J))
       enddo
*
*   LAST, POINTS N-5 TO N
*
      do i = n-5,n
         irow = i-n+13
         hdpbdt = 0.0d 00
         hdqbdt = 0.0d 00
         do k = 1,13
            aik = a(irow,k)
            loc = n-13+k
            hdpbdt = hdpbdt+aik*pf(loc,J)
            hdqbdt = hdqbdt+aik*qf(loc,J)
         enddo
         tpf(i) = hdpbdt
         tqf(i) = hdqbdt
       enddo
       return
      end
c------------------------------------------------------------------
      subroutine diff2(J)

c------------------------------------------------------------------
      implicit real*8 (a-h, o-z)
      include 'mpif.h'
*
      COMMON/cons/zero,half,tenth,one,two,three,ten
     :      /grid/steplin,r(MN),rp(MN),rpor(MN),rnt,h,hp,n
     :      /lic13/a(13,13)
     :      /tab/tgl(MN),tgs(MN)
     :      /basis/gl(MN,MNBAS),gs(MN,MNBAS),alpha(MNBAS,MNS),mg(MNBAS)
     :  /info2/qqf(MN,MNBAS),qql(MN,MNBAS),qqs(MN,MNBAS)
*
      equivalence (a1,a(7,1)),(a2,a(7,2)),(a3,a(7,3)),
     :            (a4,a(7,4)),(a5,a(7,5)),(a6,a(7,6))
*
*   COMPUTE DERIVATIVE IN THREE SEPARATE REGIONS
*
*   FIRST, POINTS 1 TO 6
*
      do i = 1,6
         hdpbdt = 0.0d 00
         hdqbdt = 0.0d 00
         do k = 1,13
            aik = a(i,k)
            hdpbdt = hdpbdt+aik*gl(k,J)
            hdqbdt = hdqbdt+aik*gs(k,J)
         enddo
         tgl(i) = hdpbdt
         tgs(i) = hdqbdt
      enddo
*
*   NEXT, POINTS 7 TO N-6
*
*   SPECIAL TREATMENT FOR THIS REGION BECAUSE OF THE SYMMETRY OF
*   THE DIFFERENTIATION FORMULA
*
      do i = 7,n-6
         tgl(i) =  a1*(gl(i-6,J)-gl(i+6,J))+a2*(gl(i-5,J)-gl(i+5,J))
     :           +a3*(gl(i-4,J)-gl(i+4,J))+a4*(gl(i-3,J)-gl(i+3,J))
     :           +a5*(gl(i-2,J)-gl(i+2,J))+a6*(gl(i-1,J)-gl(i+1,J))
         tgs(i) =  a1*(gs(i-6,J)-gs(i+6,J))+a2*(gs(i-5,J)-gs(i+5,J))
     :           +a3*(gs(i-4,J)-gs(i+4,J))+a4*(gs(i-3,J)-gs(i+3,J))
     :           +a5*(gs(i-2,J)-gs(i+2,J))+a6*(gs(i-1,J)-gs(i+1,J))
      enddo
*
*   LAST, POINTS N-5 TO N
*
      do i = n-5,n
         irow = i-n+13
         hdpbdt = 0.0d 00
         hdqbdt = 0.0d 00
         do k = 1,13
            aik = a(irow,k)
            loc = n-13+k
            hdpbdt = hdpbdt+aik*gl(loc,J)
            hdqbdt = hdqbdt+aik*gs(loc,J)
         enddo
         tgl(i) = hdpbdt
         tgs(i) = hdqbdt
       enddo
      return
      end
c---------------------------------------------------------------
      subroutine nucparm
c---------------------------------------------------------------

          implicit real*8 (a-h,o-z)
      include 'mpif.h'
      common /grid/steplin,r(MN),rp(MN),rpor(MN),rnt,h,hp,n
     :  /cons/zero,half,tenth,one,two,three,ten
     :  /npar/parm(2),nparm
     :  /def1/amass,z
     :  /wave/pf(MN,MNBAS),qf(MN,MNBAS),mf(MNBAS)

      nparm=2
      pi = 3.141592653589790d 000
      tfm=2.3d0
      ainfcm = 0.529177249d-08
      fmtoau = 1.0d-13/ainfcm
      t=tfm*fmtoau
      const= 4.0d0* log (3.0d0)
      parm(2)=t/const
      rrmsfm=0.836d0*amass**(1.d0/3.d0)+0.57d0
      rrms = rrmsfm*fmtoau
      factor = rrms**2-(7.0d 00/5.0d 00)
     :                    *(pi**2)*(parm(2)**2)
      parm(1) =   sqrt (5.0d 00/3.0d 00)* dsqrt (factor)
      return
      end
c-----------------------------------------------------------
      subroutine findk(orba,orbb,orbp,orbq,ia,ib,ip,iq,
     $jeven,kkk)
c     this subroutine finds the common k value
c-----------------------------------------------------------
      implicit real*8 (a-h,o-z)
      include 'mpif.h'
      real*8 jj
      dimension ja(MAXVAL),jb(MAXVAL),jp(MAXVAL),jq(MAXVAL),
     :          jeven(MAXVAL)

      do i=1,MAXVAL
          jeven(i)=0
      enddo

      iab=ia*ib
      ipq=ip*iq

      jminab=idint(dmax1(orba,orbb)-dmin1(orba,orbb))
      jminpq=idint(dmax1(orbp,orbq)-dmin1(orbp,orbq))

      jmaxab=idint(orba+orbb)
      jmaxpq=idint(orbp+orbq)

      ll=0
      kk=0
      do i=jminab,jmaxab
          jtot=jmaxab+i
          jj=float(jtot)-2.0d0*(float(jtot/2))
          if(jj.ne.0.and.iab.gt.0)then
          ll=ll+1
          ja(ll)=i
      endif
      if(jj.eq.0.and.iab.lt.0)then
      kk=kk+1
      jb(kk)=i
      endif
      enddo
      mm=0
      nn=0
      do i=jminpq,jmaxpq
      jtot=jmaxpq+i
      jj=float(jtot)-2.0*(float(jtot/2))
      if(jj.ne.0.and.ipq.gt.0)then
      mm=mm+1
      jp(mm)=i
      endif
      if(jj.eq.0.and.ipq.lt.0)then
      nn=nn+1
      jq(nn)=i
      endif
      enddo
      kkk=0
      do i=1,ll
      do j=1,mm
      if(ja(i).eq.jp(j))then
      kkk=kkk+1
      jeven(kkk)=ja(i)
c      print*,'jeven',jeven(kkk)
      endif
      enddo
      enddo

      do i=1,ll
      do j=1,nn
      if(ja(i).eq.jq(j))then
      kkk=kkk+1
      jeven(kkk)=ja(i)
c      print*,'jeven',jeven(kkk)
      endif
      enddo
      enddo
      do i=1,kk
      do j=1,mm
      if(jb(i).eq.jp(j))then
      kkk=kkk+1
      jeven(kkk)=jb(i)
c      print*,'jeven',jeven(kkk)
      endif
      enddo
      enddo
      do i=1,kk
      do j=1,nn
      if(jb(i).eq.jq(j))then
      kkk=kkk+1
      jeven(kkk)=jb(i)
c      print*,'jeven',jeven(kkk)
      endif
      enddo
      enddo
c      print*,'kkk,jeven',kkk,jeven(kkk)
      return
      end

      subroutine integral(nbasis)
c-------------------------------------------------------------------
c This evaluates the direct and exchange matrix elements numericaly.
c It computes submatrices skll,skss and c skls. the submatrix skls 
c is not computed as the fock matrix is symmetric. P(i,j) is the
c density matrix.
c-------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      include 'mpif.h'
      dimension momj(MAXVAL),jq(MNS),mmomj(MAXVAL),moomj(MAXVAL)
c     :   ,kapdac(MNB),kappac(MNB)
c     :   kapc(MNB),kapa(MNB),kapcc(MNB),kapaa(MNB)
      common/cons/zero,half,tenth,one,two,three,ten
     : /grid/steplin,r(MN),rp(MN),rpor(MN),rnt,h,hp,n
     : /dumm/tmp1(MN),tmp2(MN),tmp3(MN),tmp4(MN)
     : /wave/pf(MN,MNBAS),qf(MN,MNBAS),mf(MNBAS)
     : /basis/gl(MN,MNBAS),gs(MN,MNBAS),alpha(MNBAS,MNS),mg(MNBAS)
     :  /diag/eigv(MNB,MNB),eng(MNB)
     : /adat/occj(MNOCC),occl(MNOCC),unocl(MNS),unocj(MNS)
     : /double/df_two(MNB,MNB)
     : /quan1/kap(MNS),nak(MNOCC),np(MNOCC),nsy(MNS),nkk(MNS)
     : /datt/nsym,nocc,nbas(MNS),ind
     :  /jvalue/orbj(MNOCC,MNS),iocc(MNBAS)
     :  /info/nocorb(MNOCC),nucorb(MNB)
     :  /skip/nskip(MNS)
     :  /dens/p(MNB,MNB)
c     :  /kapp/kappa(MNB),kappc(MNB),iqc(MNOCC)
     :  /jv/occ(MNOCC),orbe(MNB),iq(MNB)
     :  /info1/occk(MNOCC,MNS)
     :  /kapp/kappa(MNB),kappc(MNB),iqc(MNOCC)
     :  /nucop/nucopt,isms
     :  /ncc/c1,c2,c3,c4,cnc5c(5,5),c5num(5,5),pi,c
     :  /def1/amass,z

       data jq/1,-1,1,-1,1,-1,1,-1,1,-1,1/
      do ia = 1,MNB
      do ib = 1,MNB
      df_two(ib,ia) = zero
      enddo
      enddo

c      option = 1
      aiso33=zero
      aiso333=zero
      aiso33_3=zero
      aiso3_33=zero
      mtp=n
      do isym=1,nsym
      do ia=1,nbas(isym)
      ka=ia+nskip(isym)
      ja=ia+2*nskip(isym)
      index1=ia+nbas(isym)+2*nskip(isym)
      do ib=1,ia
      index2=ib+nbas(isym)+2*nskip(isym)
      jb=ib+2*nskip(isym)
      kb=ib+nskip(isym)

      df_two(ja,jb)=zero
      df_two(index1,index2)=zero

      
      do jsym=1,nsym
      nn=nocorb(jsym)
      if(nn.ne.0)then
      do ic=1,nocorb(jsym)
      orbc=occ(jsym)
      orba=orbe(ka)
      orbb=orbe(kb)
      iaa=iq(ka)
      ibb=iq(kb)
c     icc=jq(jsym)
      kc=ic+nskip(jsym)
      icc=iq(kc)
c---------------------------------------------------------------
c     direct part
c----------------------------------------------------------------
c     find k
      call findk(orba,orbb,orbc,orbc,iaa,ibb,icc,icc,momj,lmax)
      if(lmax.ne.0) then
      do  lk=1,lmax
      ll=momj(lk)
      if(ll.eq.0)then
      call slater1(ka,kc,kb,kc,ll,rkll,rkss)

      df_two(ja,jb)=df_two(ja,jb)+rkll*(2.d0*orbc+1.d0)

      df_two(jb,ja)=df_two(ja,jb)

      df_two(index1,index2)=df_two(index1,index2)+rkss*(2.d0*orbc+1.d0)
      
      df_two(index2,index1)=df_two(index1,index2)
      endif
      enddo
      endif
c------------------------------------------------------------------
c     exchange part
c------------------------------------------------------------------
      call findk(orba,orbc,orbc,orbb,iaa,icc,icc,ibb,momj,kmax)
      if(kmax.ne.0)then
      do ik=1,kmax
      kk=momj(ik)
      fact=(2.d0*orbc+1.d0)*clrx (kappa(ka),kk,kappa(kc))**2

      call slater2(ka,kc,kc,kb,kk,rkll,rkss,rkls,rksl)
      df_two(ja,jb)=df_two(ja,jb)-rkll*fact
      df_two(jb,ja)=df_two(ja,jb)
      df_two(index1,index2)=df_two(index1,index2)-rkss*fact
      df_two(index2,index1)=df_two(index1,index2)
      df_two(ja,index2)=df_two(ja,index2)-rkls*fact
      df_two(index2,ja)=df_two(ja,index2)
      df_two(index1,jb)=df_two(index1,jb)-rksl*fact
      df_two(jb,index1)=df_two(index1,jb)
      enddo
      endif
      
       if(isms.eq.1)then
c-----------------------------------------------------------------
c     exchange part of isotope shift
c------------------------------------------------------------------
c       print*,'ka,kc',ka,kc,kb
       orbac=orba+orbc+1
       onein=(-1)**orbac
       iaacc=iaa*icc
***********************************************
*   Parity selection rule
*   1/2[1-\eta_a*\eta_c*(-1)^{j_a+j_c+k}]
*   same value for both <a||c> and <c||b>
*******************************************
       iac=half*(1-iaacc*onein)
        if(iac.eq.1) then

c        if(orba.eq.orbc.or.orbc.eq.orba+1.or.orbc.eq.orba-1) then
c       jminac=idint(dmax1(orba,orbc)-dmin1(orba,orbc))
c       jmaxac=idint(orba+orbc)
c       do kk=jminac,jmaxac
c        if (kk.eq.1) then


        call findk(orba,orbc,orbc,orbb,iaa,icc,icc,ibb,momj,kmax)
       if(kmax.ne.0)then
       do ik=1,kmax
       kk=momj(ik)
       if (kk.eq.1)then

       kacl=((kappa(ka)*(kappa(ka)+1))-(kappa(kc)*(kappa(kc)+1)))*half
       kacs=((-kappa(ka)*(-kappa(ka)+1))+(kappa(kc))*(-kappa(kc)+1))
     : *half
       kcal=((kappa(kc)*(kappa(kc)+1))-(kappa(ka)*(kappa(ka)+1)))*half
       kcas=((-kappa(kc)*(-kappa(kc)+1))+(kappa(ka)*(-kappa(ka)+1)))
     : *half

c       print*,'kkka,kkkc,kkkb',ka,kc,kb
       call gen_iso(ka,kc,kc,kb,kacl,kacs,kcal,kcas,rtll,
     :     rtss,rtls,rtsl)
       ama=z*1836.15267247+(amass-z)*1838.6836605
c       ama=(amass*1822.888)
c       ama=amass
c      fact1=((-1)**(orba+orbc+1))
      am=one/(ama)
      fact=(2.d0*orbc+1.d0)*am*
     : clrx (kappa(ka),kk,kappa(kc))**2

      df_two(ja,jb)=df_two(ja,jb)-rtll*fact
      df_two(jb,ja)=df_two(ja,jb)

c       if(ja.eq.jb.and.ja.eq.3)then
c        aiso33=aiso33-rtll*fact
c       endif


      df_two(index1,index2)=df_two(index1,index2)-rtss*fact
      df_two(index2,index1)=df_two(index1,index2)

c       if(index1.eq.33.and.index2.eq.33)then
c        aiso3333=aiso3333-rtss*fact
c       endif


      df_two(ja,index2)=df_two(ja,index2)-rtls*fact
      df_two(index2,ja)=df_two(ja,index2)

c       if(ja.eq.3.and.index2.eq.33)then
c        aiso3_33=aiso3_33-rtls*fact
c       endif
      df_two(index1,jb)=df_two(index1,jb)-rtsl*fact
      df_two(jb,index1)=df_two(index1,jb)

c       if(index1.eq.33.and.jb.eq.3)then
c        aiso33_3=aiso33_3-rtsl*fact
c       endif

      endif

c      if((ia.eq.ib).and.(ia.eq.3)) then
c      print*,'Iso_integ=',isym,jsym,ia,ib,ic,aiso33,aiso3333,fact
c     :   aiso33_3,aiso3333,aiso3_33,fact
c      endif

      enddo
      endif
      endif
      endif


*  end of SMS
      enddo
      endif
      enddo

      enddo
      enddo
      enddo
c      print*,'Iso_integ_f=',aiso33,aiso3_33,aiso33_3,aiso3333
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine nucleus(tex,nmas)
      parameter (kerne = 105)
      character*2 text,tex
      dimension text(kerne)

      data text/'  ',' H','HE','LI','BE',' B',' C',' N',' O',' F','NE',
     &'NA','MG','AL','SI',' P',' S','CL','AR',' K','CA','SC','TI',' V',
     &'CR','MN','FE','CO','NI','CU','ZN','GA','GE','AS','SE','BR','KR',
     &'RB','SR',' Y','ZR','NB','MO','TC','RU','RH','PD','AG','CD','IN',
     &'SN','SB','TE',' I','XE','CS','BA','LA','CE','PR','ND','PM','SM',
     &'EU','GD','TB','DY','HO','ER','TM','YB','LU','HF','TA',' W','RE',
     &'OS','IR','PT','AU','HG','TL','PB','BI','PO','AT','RN','FR','RA',
     &'AC','TH','PA',' U','NP','PU','AM','CM','BK','CF','ES','FM','MD',
     &'NO','LR','04'/

      do 1 nl=1,kerne
      if  (tex.eq.text(nl)) goto 2
    1 continue
      write(STDOUT,201) tex,nmas
  201 format(//,' ATOM ',a2,i4,'  UNKNOWN')
      stop
    2 nmas=nl-1
      write(STDOUT,202) tex,nmas
  202 format(//,'    ATOM =   ',a2,8x,' ATOMIC NUMBER = ',i3)
      return
      end
c------------------------------------------------------------------
      subroutine twoe(isym,nbasis)
c-------------------------------------------------------------------
c This evaluates the direct and exchange matrix elements numericaly.
c It computes submatrices skll,skss and c skls. the submatrix skls 
c is not computed as the fock matrix is symmetric. P(i,j) is the
c density matrix.
c------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      include 'mpif.h'
      real*8 df_two_local
      dimension momj(MAXVAL),jq(MNS),mmomj(MAXVAL),moomj(MAXVAL)
     :,df_two_local(MNB,MNB)
c     :    ,kapdac(MNB),kappac(MNB)
c     :   kapc(MNB),kapa(MNB),kapcc(MNB),kapaa(MNB)
      common/cons/zero,half,tenth,one,two,three,ten
     : /grid/steplin,r(MN),rp(MN),rpor(MN),rnt,h,hp,n
     : /dumm/tmp1(MN),tmp2(MN),tmp3(MN),tmp4(MN)
     : /wave/pf(MN,MNBAS),qf(MN,MNBAS),mf(MNBAS)
     : /basis/gl(MN,MNBAS),gs(MN,MNBAS),alpha(MNBAS,MNS),mg(MNBAS)
     :  /diag/eigv(MNB,MNB),eng(MNB)
     : /adat/occj(MNOCC),occl(MNOCC),unocl(MNS),unocj(MNS)
     : /double/df_two(MNB,MNB)
     : /quan1/kap(MNS),nak(MNOCC),np(MNOCC),nsy(MNS),nkk(MNS)
     : /datt/nsym,nocc,nbas(MNS),ind
     :  /jvalue/orbj(MNOCC,MNS),iocc(MNBAS)
     :  /info/nocorb(MNOCC),nucorb(MNB)
     :  /skip/nskip(MNS)
     :  /dens/p(MNB,MNB)
     :  /kapp/kappa(MNB),kappc(MNB),iqc(MNOCC)
     :  /jv/occ(MNOCC),orbe(MNB),iq(MNB)
     :  /info1/occk(MNOCC,MNS)
     :  /nucop/nucopt,isms
     :  /rest/rll(MN),rss(MN),rls(MN),rsl(MN)
     :  /ncc/c1,c2,c3,c4,cnc5c(5,5),c5num(5,5),pi,c
     :  /def1/amass,z

       data jq/1,-1,1,-1,1,-1,1,-1,1,-1,1/
      call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)

      do ib = 1,MNB
      do ia = 1,MNB
      df_two(ib,ia) = zero
      df_two_local(ib,ia) = zero
      enddo
      enddo
c      option = 1
      mtp=n
      do ia=iproc+1,nbas(isym),nproc
      ka=ia+nskip(isym)
      ja=ia+2*nskip(isym)
      index1=ia+nbas(isym)+2*nskip(isym)
      do ib=1,ia
      index2=ib+nbas(isym)+2*nskip(isym)
      jb=ib+2*nskip(isym)
      kb=ib+nskip(isym)

      df_two(ja,jb)=0
      df_two(index1,index2)=0

      do jsym=1,nsym
      nn=nocorb(jsym)
      if(nn.ne.0)then
      do ic=1,nocorb(jsym)
      orbc=occ(jsym)
      orba=orbe(ka)
      orbb=orbe(kb)
      iaa=iq(ka)
      ibb=iq(kb)
*    Modified by Sonjoy
c     icc=jq(jsym)
      kc=ic+nskip(jsym)
      icc=iq(kc)
c
c     direct part
c
c     find k
      call findk(orba,orbb,orbc,orbc,iaa,ibb,icc,icc,momj,lmax)
      if(lmax.ne.0) then
      do  lk=1,lmax
      ll=momj(lk)

      if(ll.eq.0)then
      call slater1(ka,kc,kb,kc,ll,rkll,rkss)
      df_two_local(ja,jb)=df_two_local(ja,jb)+rkll*(2.d0*orbc+1.d0)

      df_two_local(jb,ja)=df_two_local(ja,jb)

      df_two_local(index1,index2)=df_two_local(index1,index2)+
     :  rkss*(2.d0*orbc+1.d0)
      
      df_two_local(index2,index1)=df_two_local(index1,index2)

      endif
      enddo
      endif
c
c     exchange part
c
c     find k
      call findk(orba,orbc,orbc,orbb,iaa,icc,icc,ibb,momj,kmax)
      if(kmax.ne.0)then
      do ik=1,kmax
      kk=momj(ik)

      fact=(2.d0*orbc+1.d0)*clrx (kappa(ka),kk,kappa(kc))**2

      call slater2(ka,kc,kc,kb,kk,rkll,rkss,rkls,rksl)

      df_two_local(ja,jb)=df_two_local(ja,jb)-rkll*fact
      df_two_local(jb,ja)=df_two_local(ja,jb)
      df_two_local(index1,index2)=df_two_local(index1,index2)
     :   -rkss*fact
      df_two_local(index2,index1)=df_two_local(index1,index2)
      df_two_local(ja,index2)=df_two_local(ja,index2) -rkls*fact
      df_two_local(index2,ja)=df_two_local(ja,index2)
      df_two_local(index1,jb)=df_two_local(index1,jb) -rksl*fact
      df_two_local(jb,index1)=df_two_local(index1,jb)

      enddo
      endif

       if (isms.eq.1) then
c--------------------------------------------------------------

        orbac=orba+orbc+1
       onein=(-1)**orbac
       iaacc=iaa*icc
***********************************************
*   Parity selection rule
*   1/2[1-\eta_a*\eta_c*(-1)^{j_a+j_c+k}]
*   same value for both <a||c> and <c||b>
*******************************************
       iac=half*(1-iaacc*onein)
        if(iac.eq.1) then


c       if(orba.eq.orbc.or.orbc.eq.orba+1.or.orbc.eq.orba-1) then
c       jminac=idint(dmax1(orba,orbc)-dmin1(orba,orbc))
c       jmaxac=idint(orba+orbc)
c       do kk=jminac,jmaxac
c        if (kk.eq.1) then


       call findk(orba,orbc,orbc,orbb,iaa,icc,icc,ibb,momj,kmax)
       if(kmax.ne.0)then
       do ik=1,kmax
       kk=momj(ik)
       if (kk.eq.1)then

       kacl=((kappa(ka)*(kappa(ka)+1))-(kappa(kc)*(kappa(kc)+1)))*half
       kacs=((-kappa(ka)*(-kappa(ka)+1))+(kappa(kc)*(-kappa(kc)+1)))
     : *half
       kcal=((kappa(kc)*(kappa(kc)+1))-(kappa(ka)*(kappa(ka)+1)))*half
       kcas=((-kappa(kc)*(-kappa(kc)+1))+(kappa(ka)*(-kappa(ka)+1)))
     : *half
       
c      print*,'kacl,kacs,kcal,kcas,ka,kc',kacl,kacs,kcal,kcas,kappa(ka),
c     : kappa(kc)
      call gen_iso(ka,kc,kc,kb,kacl,kacs,kcal,kcas,rtll,
     :     rtss,rtls,rtsl)
       ama=z*1836.15267247+(amass-z)*1838.6836605
c      ama=(amass*1822.888)
c       fact1=((-1)**(orba+orbc+1))
c      ama=amass
      am=one/(ama)
      fact=(2.d0*orbc+1.d0)*am*
     :  clrx (kappa(ka),kk,kappa(kc))**2
c       print*,'fact,ka,kc,kb',fact,ka,kc,kb      
      df_two_local(ja,jb)=df_two_local(ja,jb)-rtll*fact
      df_two_local(jb,ja)=df_two_local(ja,jb)
      df_two_local(index1,index2)=df_two_local(index1,index2)-
     : rtss*fact
      df_two_local(index2,index1)=df_two_local(index1,index2)
      df_two_local(ja,index2)=df_two_local(ja,index2)-rtls*fact
      df_two_local(index2,ja)=df_two_local(ja,index2)
      df_two_local(index1,jb)=df_two_local(index1,jb)-rtsl*fact
      df_two_local(jb,index1)=df_two_local(index1,jb)
      endif
*   end of selection iac.eq.1.and.kk.eq.1 
      enddo
*   end of kk loop
      endif
*   end of kb selection
      endif
      endif
*   end of SMS calculation

      enddo
      endif
      enddo
      enddo
      enddo
             call MPI_ALLREDUCE(df_two_local,df_two,MNB*MNB,MPI_REAL8
     : ,MPI_SUM,MPI_COMM_WORLD,ierr)

      do i= 1,MNB
       do j=1,MNB
       df_two_local(i,j)=0.0
       enddo
      enddo


      return
      end
c***************************************************************
      subroutine eigen (a,b,nsub,valu,msub,anorm,nnmax)
      parameter(nmax=1600)
c     eigenvalues and eigenvectors of a real symmetric matrix
c     a ----- secular determinant (input)
c             eigenvetors  column-wise (output)
c             a is destroyed
c     b ----- nsub x nsub array for intermidiate storage
c     nsub ----- order of a (input)
c                must be equal to or greater than 1
c     valu ----- eigenvalues in decreasing order (output)
c     msub -----  number of eigenvalues (input)
c                 if equal to -1, eigenvectors are not computed
c     anorm ----- matrix norm which is a sqrt of sum of a(i,j)**2
c                                        over all i & j (output)
c     nmax ----- row dimension of a & b
c
      implicit real*8 (a-h,o-z)
      include 'mpif.h'
      integer index
      dimension a(nnmax,nnmax), b(nnmax,nnmax), valu(nnmax)
      dimension diag(nmax), superd(nmax), wvec(nmax),pvec(nmax),
     1          qvec(nmax), vall(nmax), q(nmax), u(nmax),
     2          index(nmax), factor(nmax), v(nmax), t(nmax,3)
      common /f202/ diag, superd, wvec, pvec, t
      equivalence (wvec(1),vall(1),factor(1),u(1)),
     1            (pvec(1),qvec(1),q(1),v(1))
      data e1/2.3d-15/
c     householder similarity transformation to co-diagonal form
c     ----------- ---------- -------------- -- ----------- ----
      lpout = 0
      n=nsub
      m=msub
      do 10 i = 1,n
      valu(i) = 0.0d+00
   10 continue
      if (n.lt.0) go to 220
      if (m.le.0.or.n.le.2) go to 50
      do 40 i=2,n
      do 40 j=2,n
      b(j,i)=0.0d0
      if (i.eq.j) b(j,i)=1.0d0
   40 continue
   50 do 200 i=1,n
      i1=i+1
      i2=i1+1
      if (i2.gt.n) go to 160
      sum=0.0d0
      do 70 j=i2,n
   70 sum=sum+a(j,i)**2
      if (sum.eq.0.0d0) go to 160
      j=i1
      temp=a(j,i)
      sum=dsqrt(sum+ temp  **2)
      a(j,i)=-dsign(sum, temp  )
      wvec(j )=dsqrt(     1.0d0+dabs( temp  )/sum)
      div=dsign(    wvec(j )*sum, temp  )
      do 85 j=i2,n
   85 wvec(j)=a(j,i)/div
      scalar=0.0d0
      do 95 j=i1,n
      pvec(j)=0.0d0
      do 90 k=i1,n
   90 pvec(j)=pvec(j)+a(k,j)*wvec(k)
      scalar=scalar+pvec(j)*wvec(j)
   95 continue
      scalar=scalar/2.0d0
      do 120 j=i1,n
      qvec(j)=pvec(j)-scalar*wvec(j)
      do 120 k=i1,j
      a(k,j)=a(k,j)-(wvec(k)*qvec(j)+wvec(j)*qvec(k))
      a(j,k)=a(k,j)
  120 continue
      if (m.le.0) go to 160
      do 150 k=2,n
      temp=0.0d0
      do 140 j=i1,n
  140 temp=temp+wvec(j)*b(j,k)
      do 150 j=i1,n
      b(j,k)=b(j,k)-wvec(j)*temp
  150 continue
  160 j=i
      diag(i)=a(j,i)
      if (i.ne.n) superd(i)=a(j+1,i)
  200 continue
c     givens eigenvalue iteration from sturm chain of co-diagonal minors
  220 n=iabs (n)
      m=iabs (m)
c     calculate norm of matrix and initialize eigenvalue bounds
      anorm2=diag(1)**2
      if (n.eq.1) go to 235
      do 230 l=2,n
      q(l-1)=superd(l-1)**2
      anorm2=diag(l)**2+q(l-1)+q(l-1)+anorm2
  230 continue
  235 anorm=dsqrt(anorm2)
      if (m.eq.0) return
      do 240 l=1,m
      valu(l)=anorm
      vall(l)=-anorm
  240 continue
      eps1=anorm*e1
      if (eps1.eq.0.0d0) return
  250 do 570 l=1,m
  260 tau=(valu(l)+vall(l))/2.0d0
      if (    (tau-vall(l)).le.e1*(dabs(tau)+e1)) go to 570
      match=0
      t2=0.0d0
      t1=1.0d0
      do 450 l1=1,n
      p=diag(l1)-tau
      if (t2.ne.0.0d0) go to 330
  300 t1=dsign(1.0d0,t1)
  330 if (t1.ne.0.0d0) go to 400
      t0=-dsign(1.0d0,t2)
      t2=0.0d0
      if (q(l1-1).ne.0.0d0) go to 410
      go to 300
  400 t0=p-q(l1-1)*t2/t1
      t2=1.0d0
  410 if (t0) 440, 420, 430
  420 t2=t1
      if (t2) 440, 430, 430
  430 match=match+1
  440 t1=t0
  450 continue
      do 530 l1=l,m
      if (l1.le.match) go to 500
      if (valu(l1).le.tau) go to 260
      valu(l1)=tau
      go to 530
  500 vall(l1)=tau
  530 continue
      go to 260
  570 continue
c     eigenvectors of co-diagonal symmetric matrix--inverse iteration
      m=msub
      if (m.lt.0) return
      do 970 i=1,m
      if (i.eq.1) go to 725
      if (valu(i-1)-valu(i).lt.1.0e+4*eps1) go to 730
  725 i1=-1
  730 i1=i1+1
      do 760 l=1,n
      v(l)=eps1*(i1*l+1)
      t(l,2)=diag(l)-valu(i)
      if (l.lt.n) go to 740
      t(l,3)=0.0d0
      go to 760
  740 t(l,3)=superd(l)
      if (t(l,3).eq.0.0d0) t(l,3)=eps1
      t(l+1,1)=t(l,3)
  760 continue
      do 820 j=1,n
      t(j,1)=t(j,2)
      t(j,2)=t(j,3)
      t(j,3)=0.0d0
      go to 780
  770 t(j,1)=eps1
c**** replace function name  in next      statement  in single precision
  780 vtemp=dabs(t(j,1))
      if (vtemp.lt.eps1) go to 770
      if (j.eq.n) go to 820
      index(j)=0
      if (dabs(t(j+1,1)).le.vtemp) go to 810
      index(j)=1
      do 800 k=1,3
      vtemp=t(j,k)
      t(j,k)=t(j+1,k)
      t(j+1,k)=vtemp
  800 continue
  810 vtemp   =t(j+1,1)/t(j,1)
      factor(j)=vtemp
      t(j+1,2)=t(j+1,2)-  vtemp  *t(j,2)
      t(j+1,3)=t(j+1,3)-  vtemp  *t(j,3)
  820 continue
      iter=1
      if (i1.gt.0) go to 920
  860 do 870 l1=1,n
      l=n+1-l1
      v( l )=(v( l )-t(l,2)*v( l+1 )-t(l,3)*v( l+2 ))/t(l,1)
  870 continue
      go to (875,920), iter
  875 iter=2
      if (n.eq.1) go to 860
  880 do 910 l=2,n
      if (index(l-1).eq.0) go to 900
      vtemp=v(l-1)
      v(l-1)=v(l)
      v(l)=vtemp
  900 v(l)=v(l)-factor(l-1)*v(l-1)
  910 continue
      go to 860
  920 if (i1.eq.0) go to 945
      do 940 l1=1,i1
      k=i-l1
      vtemp=0.0d0
      do 930 j=1,n
  930 vtemp=vtemp+a(j,k)*v(j)
      do 940 j=1,n
  940 v(j)=v(j)-a(j,k)*vtemp
      go to (880,945), iter
  945 vnorm2=0.0d0
      do 950 l=1,n
  950 vnorm2=vnorm2+v( l )**2
      vnorm=dsqrt(vnorm2)
      do 960 j=1,n
  960 a(j,i)=v( j )/vnorm
  970 continue
      n=nsub
      if (n.le.2) return
      do 990 i=1,m
      do 980 k=2,n
      u(k)=0.0d0
      do 980 j=2,n
  980 u(k)=u(k)+b(j,k)*a(j,i)
      do 990 j=2,n
  990 a(j,i)=u(j)
      return
      end
c***********************************************************************
      subroutine slater1(ia,ib,ic,id,k,rkll,rkss)

*   the value of this  function is the  slater integral  as normally   *
*   defined in terms of the four sets of quantum numbers  a,b,c,d.     *
*                                                                      *
*   subroutines called: quad, yzk.                                     *
*                                                                      *
c***********************************************************************
*
      implicit real*8(a-h, o-z)
      include 'mpif.h'
c     dimension rhop(MN),rttk(MN),wk(MN),temp(MN),zk(MN)
      common /cons/zero,half,tenth,one,two,three,ten
     :  /grid/steplin,r(MN),rp(MN),rpor(MN),rnt,h,hp,n
     :  /ncc/c1,c2,c3,c4,cnc5c(5,5),c5num(5,5),pi,c
     :  /dumm/tmp1(MN),tmp2(MN),tmp3(MN),tmp4(MN)
     :  /wave/pf(MN,MNBAS),qf(MN,MNBAS),mf(MNBAS)
     :  /basis/gl(MN,MNBAS),gs(MN,MNBAS),alpha(MNBAS,MNS),mg(MNBAS)

      mtp=n
      do i=1,mtp
         tmp1(i)=0
         tmp3(i)=0
      enddo

      call yzk1(k,ib,id)
      do i = 2,mtp
         tmp3(i)=gl(i,ia)*gl(i,ic)*rpor(i)*tmp1(i)
      enddo
      call quad (tmp3,mtp,rkll)

      do i = 2,mtp
         tmp3(i)=gs(i,ia)*gs(i,ic)*rpor(i)*tmp1(i)
      enddo

      call quad (tmp3,mtp,rkss)

      return
      end
c---------------------------------------------------------------------
      subroutine slater2(ia,ib,ic,id,k,rkll,rkss,rkls,rksl)

*   the value of this  function is the  slater integral  as normally   *
*   defined in terms of the four sets of quantum numbers  a,b,c,d.     *
*                                                                      *
*   subroutines called: quad, yzk.                                     *
*                                                                      *
c***********************************************************************
*
      implicit real*8(a-h, o-z)
      include 'mpif.h'
c     dimension rhop(MN),rttk(MN),wk(MN),temp(MN),zk(MN)
      common /cons/zero,half,tenth,one,two,three,ten
     :  /grid/steplin,r(MN),rp(MN),rpor(MN),rnt,h,hp,n
     :  /ncc/c1,c2,c3,c4,cnc5c(5,5),c5num(5,5),pi,c
     :  /dumm/tmp1(MN),tmp2(MN),tmp3(MN),tmp4(MN)
     :  /wave/pf(MN,MNBAS),qf(MN,MNBAS),mf(MNBAS)
     :  /basis/gl(MN,MNBAS),gs(MN,MNBAS),alpha(MNBAS,MNS),mg(MNBAS)

      mtp=n
      do i=1,mtp
         tmp1(i)=0
         tmp2(i)=0
         tmp3(i)=0
         tmp4(i)=0
      enddo

      call yzk2(k,ib,id)
      do i = 2,mtp
         tmp3(i)=gl(i,ia)*pf(i,ic)*rpor(i)*tmp1(i)
      enddo
      call quad (tmp3,mtp,rkll)

      do i = 2,mtp
         tmp3(i)=gs(i,ia)*qf(i,ic)*rpor(i)*tmp2(i)
      enddo

      call quad (tmp3,mtp,rkss)
 
      do i = 2,mtp
         tmp3(i)=gs(i,ia)*qf(i,ic)*rpor(i)*tmp1(i)
      enddo

      call quad (tmp3,mtp,rksl)

      do i = 2,mtp
         tmp4(i)=gl(i,ia)*pf(i,ic)*rpor(i)*tmp2(i)
      enddo

      call quad (tmp4,mtp,rkls)

      return
      end
c***********************************************************************
      subroutine slater_iso(ia,ib,ic,id,k,rkll,rkss,rkls,rksl)

c***********************************************************************
*
      implicit real*8(a-h, o-z)
      include 'mpif.h'
c     dimension rhop(MN),rttk(MN),wk(MN),temp(MN),zk(MN)
      common /cons/zero,half,tenth,one,two,three,ten
     :  /grid/steplin,r(MN),rp(MN),rpor(MN),rnt,h,hp,n
     :  /ncc/c1,c2,c3,c4,cnc5c(5,5),c5num(5,5),pi,c
     :  /dumm/tmp1(MN),tmp2(MN),tmp3(MN),tmp4(MN)
     :  /wave/pf(MN,MNBAS),qf(MN,MNBAS),mf(MNBAS)
     :  /basis/gl(MN,MNBAS),gs(MN,MNBAS),alpha(MNBAS,MNS),mg(MNBAS)

      mtp=n
      do i=1,mtp
         tmp1(i)=0
         tmp2(i)=0
         tmp3(i)=0
         tmp4(i)=0
      enddo

      call yzk2(k,ib,id)
      do i = 2,mtp
         tmp3(i)=gl(i,ia)*pf(i,ic)*rpor(i)*tmp1(i)
      enddo
      call quad (tmp3,mtp,rkll)
      do i = 2,mtp
         tmp3(i)=gs(i,ia)*qf(i,ic)*rpor(i)*tmp2(i)
      enddo

      call quad (tmp3,mtp,rkss)

      do i = 2,mtp
         tmp3(i)=gs(i,ia)*qf(i,ic)*rpor(i)*tmp1(i)
      enddo

      call quad (tmp3,mtp,rksl)

      do i = 2,mtp
         tmp4(i)=gl(i,ia)*pf(i,ic)*rpor(i)*tmp2(i)
      enddo

      call quad (tmp4,mtp,rkls)

      return
      end
c***********************************************************************
      subroutine gen_iso(ia,ib,ic,id,kkc,kkkc,kkd,kkkd,rll,rss,
     :      rls,rsl)
c***********************************************************************
       implicit real*8(a-h, o-z)
      include 'mpif.h'
       dimension gemat1(MN),gemat2(MN),gemat3(MN),gemat4(MN)
     :  ,dgl(MN,MNBAS),dgs(MN,MNBAS),dpf(MN,MNBAS),dqf(MN,MNBAS)
     :  ,gem1(MN),gem2(MN),gem3(MN),gem4(MN)
       common /cons/zero,half,tenth,one,two,three,ten
     :  /grid/steplin,r(MN),rp(MN),rpor(MN),rnt,h,hp,n
     :  /ncc/c1,c2,c3,c4,cnc5c(5,5),c5num(5,5),pi,c
     :  /dumm/tmp1(MN),tmp2(MN),tmp3(MN),tmp4(MN)
     :  /wave/pf(MN,MNBAS),qf(MN,MNBAS),mf(MNBAS)
     :  /basis/gl(MN,MNBAS),gs(MN,MNBAS),alpha(MNBAS,MNS),mg(MNBAS)
     :  /info2/qqf(MN,MNBAS),qql(MN,MNBAS),qqs(MN,MNBAS)
     :  /tatb/tpf(MN),tqf(MN)
     :  /tab/tgl(MN),tgs(MN)
   
          

       mtp=n
       call diff1 (ic)
       call diff2 (id)

       do i=2,mtp
       dpf(i,ic)=tpf(i)
       dqf(i,ic)=tqf(i)
       enddo

       do i=2,mtp
       dgl(i,id)=tgl(i)
       dgs(i,id)=tgs(i)
       enddo
c       print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       do i=1,mtp
       gemat1(i)=zero
       gem1(i)=zero
       gemat2(i)=zero
       gem2(i)=zero
       gemat3(i)=zero
       gem3(i)=zero
       gemat4(i)=zero
       gem4(i)=zero
       enddo
c       print*,'kkc,kkkc,kkd,kkkd,ka,kc',kkc,kkkc,kkd,kkkd,ia,ic
       do i=2,mtp
***********************************
* gemat* are the integran with 1st derivative of PF and QF
**********************************
       gemat1(i)=gl(i,ia)*dpf(i,ic)
       gem1(i)=(kkc/r(i))*gl(i,ia)*pf(i,ic)*rp(i)
       gemat2(i)=gs(i,ia)*dqf(i,ic)
       gem2(i)=(kkkc/r(i))*gs(i,ia)*qf(i,ic)*rp(i)
       gemat3(i)=pf(i,ib)*dgl(i,id)
       gem3(i)=(kkd/r(i))*pf(i,ib)*gl(i,id)*rp(i)
       gemat4(i)=qf(i,ib)*dgs(i,id)
       gem4(i)=(kkkd/r(i))*qf(i,ib)*gs(i,id)*rp(i)
       enddo

       call quad (gemat1,mtp,rll1)
       rlll1=rll1/h
       call quad (gem1,mtp,rrl1)
       call quad (gemat2,mtp,rss2)
       rsss2=rss2/h
       call quad (gem2,mtp,rrs2)
       call quad (gemat3,mtp,rll3)
       rlll3=rll3/h
       call quad (gem3,mtp,rrl3)
       call quad (gemat4,mtp,rss4)
        rsss4=rss4/h
       call quad (gem4,mtp,rrs4)

c      print*,'ia,ib,ic,id',ia,ib,ic,id
      

      rllll1=rlll1-rrl1
      rssss2=rsss2-rrs2
      rllll3=rlll3-rrl3
      rssss4=rsss4-rrs4


       rll=rllll1*rllll3
       rss=rssss2*rssss4
       rls=rllll1*rssss4
       rsl=rssss2*rllll3
c       print*,'rll,rss,rls,rsl,ia,ib,ic,id',rll,rss,rls,rsl,ia,ib,ic,id
       return
       end
c********************************************************************
       subroutine gen2_iso(ia,ib,ic,id,rll,rss,rls,rsl)
c***********************************************************************
       implicit real*8(a-h, o-z)
      include 'mpif.h'
       dimension gemat1(MN),gemat2(MN),gemat3(MN),gemat4(MN)
     :  ,dgl(MN,MNBAS),dgs(MN,MNBAS),dpf(MN,MNBAS),dqf(MN,MNBAS)
       common /cons/zero,half,tenth,one,two,three,ten
     :  /grid/steplin,r(MN),rp(MN),rpor(MN),rnt,h,hp,n
     :  /ncc/c1,c2,c3,c4,cnc5c(5,5),c5num(5,5),pi,c
     :  /dumm/tmp1(MN),tmp2(MN),tmp3(MN),tmp4(MN)
     :  /wave/pf(MN,MNBAS),qf(MN,MNBAS),mf(MNBAS)
     :  /basis/gl(MN,MNBAS),gs(MN,MNBAS),alpha(MNBAS,MNS),mg(MNBAS)
     :  /info2/qqf(MN,MNBAS),qql(MN,MNBAS),qqs(MN,MNBAS)
     :  /tatb/tpf(MN),tqf(MN)
     :  /tab/tgl(MN),tgs(MN)

       mtp=n
       call diff1 (ic)
       call diff2 (id)

       do i=2,mtp
       dpf(i,ic)=tpf(i)
       dqf(i,ic)=tqf(i)
       enddo

       do i=2,mtp
       dgl(i,id)=tgl(i)
       dgs(i,id)=tgs(i)
       enddo

       do i=1,mtp
       gemat1(i)=zero
       gemat2(i)=zero
       gemat3(i)=zero
       gemat4(i)=zero
       enddo

       do i=2,mtp
***********************************
* gemat* are the integran with 1st derivative of PF and QF
**********************************
       gemat1(i)=gl(i,ia)*dgl(i,id)*rp(i)
       gemat2(i)=gs(i,ia)*dgs(i,id)*rp(i)
       gemat3(i)=pf(i,ib)*dpf(i,ic)*rp(i)
       gemat4(i)=qf(i,ib)*dqf(i,ic)*rp(i)
       enddo
       call quad (gemat1,mtp,rll1)
       rll1=rll1/h
       call quad (gemat2,mtp,rss2)
       rss2=rss2/h
       call quad (gemat3,mtp,rll3)
       rll3=rll3/h
       call quad (gemat4,mtp,rss4)
        rss4=rss4/h


      rlll1=rll1
      rsss2=rss2
      rlll3=rll3
      rsss4=rss4


       rll=rlll1*rlll3
       rss=rsss2*rsss4
       rls=rlll1*rsss4
       rsl=rsss2*rlll3

       return
       end


c***********************************************************************
      subroutine gen3_iso(ia,ib,ic,id,rll,rss,rls,rsl)
c***********************************************************************
       implicit real*8(a-h, o-z)
      include 'mpif.h'
       dimension gemat1(MN),gemat2(MN),gemat3(MN),gemat4(MN)
     :  ,dgl(MN,MNBAS),dgs(MN,MNBAS),dpf(MN,MNBAS),dqf(MN,MNBAS)
     :  ,gem1(MN),gem2(MN),gem3(MN),gem4(MN)
       common /cons/zero,half,tenth,one,two,three,ten
     :  /grid/steplin,r(MN),rp(MN),rpor(MN),rnt,h,hp,n
     :  /ncc/c1,c2,c3,c4,cnc5c(5,5),c5num(5,5),pi,c
     :  /dumm/tmp1(MN),tmp2(MN),tmp3(MN),tmp4(MN)
     :  /wave/pf(MN,MNBAS),qf(MN,MNBAS),mf(MNBAS)
     :  /basis/gl(MN,MNBAS),gs(MN,MNBAS),alpha(MNBAS,MNS),mg(MNBAS)
     :  /info2/qqf(MN,MNBAS),qql(MN,MNBAS),qqs(MN,MNBAS)
     :  /tatb/tpf(MN),tqf(MN)
     :  /tab/tgl(MN),tgs(MN)

       mtp=n
       call diff1 (ic)
       call diff2 (id)

       do i=2,mtp
       dpf(i,ic)=tpf(i)
       dqf(i,ic)=tqf(i)
       enddo

       do i=2,mtp
       dgl(i,id)=tgl(i)
       dgs(i,id)=tgs(i)
       enddo
       do i=1,mtp
       gemat1(i)=zero
c       gem1(i)=zero
       gemat2(i)=zero
c       gem2(i)=zero
       gemat3(i)=zero
c       gem3(i)=zero
       gemat4(i)=zero
c       gem4(i)=zero
       enddo
c       print*,'kkc,kkkc,kkd,kkkd',kkc,kkkc,kkd,kkkd
       do i=2,mtp
***********************************
* gemat* are the integran with 1st derivative of PF and QF
**********************************
c       print*,'i,ia,ib',i,ia,ib
c        print*,'gl,gsia,pf,qfib',gl(i,ia),gs(i,ia),pf(i,ib),qf(i,ib)
       gemat1(i)=gl(i,ia)*dgl(i,id)*rp(i)
c       gem1(i)=(kkc/(2*r(i)))*gl(i,ia)*gl(i,id)*rp(i)
       gemat2(i)=gs(i,ia)*dgs(i,id)*rp(i)
c       gem2(i)=(kkkc/(2*r(i)))*gs(i,ia)*gs(i,id)*rp(i)
       gemat3(i)=pf(i,ib)*dpf(i,ic)*rp(i)
c       gem3(i)=(kkd/(2*r(i)))*pf(i,ib)*pf(i,ic)*rp(i)
       gemat4(i)=qf(i,ib)*dqf(i,ic)*rp(i)
c       gem4(i)=(kkkd/(2*r(i)))*qf(i,ib)*qf(i,ic)*rp(i)
c       print*,'gen',gem1(i),gem2(i),gem3(i),gem4(i)
       enddo

      
       call quad (gemat1,mtp,rll1)
       rlll1=rll1/h
c       call quad (gem1,mtp,rrl1)
       call quad (gemat2,mtp,rss2)
       rsss2=rss2/h
c       call quad (gem2,mtp,rrs2)
       call quad (gemat3,mtp,rll3)
       rlll3=rll3/h
c       call quad (gem3,mtp,rrl3)
       call quad (gemat4,mtp,rss4)
        rsss4=rss4/h
c       call quad (gem4,mtp,rrs4)

c      print*,'rll1,rrl1',rll1,rrl1
c      print*,'rss2,rrs2',rss2,rrs2
c      print*,'rll3,rrl3',rll3,rrl3
c      print*,'rss4,rrs4',rss4,rrs4

c      rlll1=rll1
c      rsss2=rss2
c      rlll3=rll3
c      rsss4=rss4

c        rlll1=rll1
c        rsss2=rss2
c        rlll3=rll3
c        rsss4=rss4

       rll=rlll1*rlll3
       rss=rsss2*rsss4
       rls=rlll1*rsss4
       rsl=rsss2*rlll3

       return
       end



c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       subroutine iso_mat(ia,ib,ic,id,rll,rss,rls,rsl)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
        implicit real*8(a-h, o-z)
      include 'mpif.h'
c       dimension clg(MN,MNBAS),csm(MN,MNBAS),dlg(MN,MNBAS),dsm(MN,MNBAS)
      dimension gemat1(MN),gemat2(MN),gemat3(MN),gemat4(MN)
c     :  ,dcsm(MN,MNBAS),ddlg(MN,MNBAS),ddsm(MN,MNBAS),dclg(MN,MNBAS)
       common /cons/zero,half,tenth,one,two,three,ten
     :  /grid/steplin,r(MN),rp(MN),rpor(MN),rnt,h,hp,n
     :  /ncc/c1,c2,c3,c4,cnc5c(5,5),c5num(5,5),pi,c
     :  /dumm/tmp1(MN),tmp2(MN),tmp3(MN),tmp4(MN)
     :  /wave/pf(MN,MNBAS),qf(MN,MNBAS),mf(MNBAS)
     :  /basis/gl(MN,MNBAS),gs(MN,MNBAS),alpha(MNBAS,MNS),mg(MNBAS)
c     :  /rest/rll(MN),rss(MN),rls(MN),rsl(MN)
c     :  /info2/qqf(MN,MNBAS),qql(MN,MNBAS),qqs(MN,MNBAS)
c     :  /tatb/tpf(MN),tqf(MN)
c     :  /tab/tgl(MN),tgs(MN)


       mtp=n
        do i=1,mtp
       gemat1(i)=0
       gemat2(i)=0
       gemat3(i)=0
       gemat4(i)=0
       enddo
       do i=2,mtp
       gemat1(i)=(gl(i,ia)*qf(i,ic))*rp(i)
       enddo
       call quad (gemat1,mtp,rll1)
       do i=2,mtp
       gemat2(i)=(gs(i,ia)*pf(i,ic))*rp(i)
       enddo
       call quad (gemat2,mtp,rss2)
       do i=2,mtp
       gemat3(i)=(pf(i,ib)*gs(i,id))*rp(i)
       enddo
       call quad (gemat3,mtp,rss3)
       do i=2,mtp
       gemat4(i)=(qf(i,ib)*gl(i,id))*rp(i)
       enddo
       call quad (gemat4,mtp,rll4)
       rll=rll1*rll4
       rss=rss2*rss3
       rls=rll1*rss3
       rsl=rss2*rll4
       return
       end
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
      subroutine yzk1(k,i,j)
c-----------------------------------------------------------------------
c   this subroutine evaluates hartree y- and z-functions:              *
c                                                                      *
c               (k)            (k)           (k)                       *
c              y   (i,j;r) =  z   (i,j;r) + w   (i,j;r)                *
c                                                                      *
c   where                                                              *
c                                                                      *
c    (k)                                                               *
c   z   (i,j;r) =  i ( (s/r)   (p (s)*p (s) + q (s)*q (s)) ; 0 - r )   *
c                                i     j       i     j                 *
c                                                                      *
c   where                                                              *
c                                                                      *
c    (k)                    k+1                                        *
c   w   (i,j;r) =  i ( (r/s)   (p (s)*p (s) + q (s)*q (s)) ; r -       *
c                                i     j       i     j    infinity )   *
c                                                                      *
c   where  i ( g(r,s) ; range )  denotes the integral of g(r,s) over   *
c   range  in  s .  the y-function is tabulated in  common/tatb/  in   *
c   array  tb , the z-function in array ta .                           *
c-----------------------------------------------------------------------
      implicit real*8(a-h, o-z)
      include 'mpif.h'
      dimension rhop(MN),rttk(MN),wk(MN),temp(MN),zk(MN)
      common /cons/zero,half,tenth,one,two,three,ten
     :  /grid/steplin,r(MN),rp(MN),rpor(MN),rnt,h,hp,n
     :  /ncc/c1,c2,c3,c4,cnc5c(5,5),c5num(5,5),pi,c
     :  /dumm/tmp1(MN),tmp2(MN),tmp3(MN),tmp4(MN)
     :  /wave/pf(MN,MNBAS),qf(MN,MNBAS),mf(MNBAS)
     :  /basis/gl(MN,MNBAS),gs(MN,MNBAS),alpha(MNBAS,MNS),mg(MNBAS)
      accy=1.0d-12
c----------------------------------------------------------------------
c   for  k > 0  compute  r^k   and store in  rttk
c----------------------------------------------------------------------
      if (k .gt. 0) then
         do 1 ii = 2,n
            rttk(ii) = r(ii)**k
    1    continue
      endif
c---------------------------------------------------------------------
c     determine maximum tabulation point as location beyond which rhop (see 
c comment statements below) would be zero; determine other important locations
c--------------------------------------------------------------------
      mtp = min(mf(i),mf(j))
      mtpp1 = mtp+1
      mtpp3 = mtp+3
      mtpp4 = mtp+4
c--------------------------------------------------------------------
c   compute rp(s)*(p_i (s)*p_j (s)+q_i (s)*q_j (s)) and store in rhop
c--------------------------------------------------------------------
      do 2 ii = 2,mtp
         rhop(ii) = rp(ii)*(pf(ii,i)*pf(ii,j)+qf(ii,i)*qf(ii,j))
    2 continue
c--------------------------------------------------------------------
c                   fill array temp with r**k * rhop
c--------------------------------------------------------------------
      temp(1) = zero
      if (k .eq. 0) then
         do 3 ii = 2,mtp
            temp(ii) = rhop(ii)
    3    continue
      else
         do 4 ii = 2,mtp
            temp(ii) = rttk(ii)*rhop(ii)
    4    continue
      endif
c---------------------------------------------------------------------
c                 set an additional four points to zero
c---------------------------------------------------------------------
      do 5 ii = mtpp1,mtpp4
         temp(ii) = zero
    5 continue
c-----------------------------------------------------------------------
c           compute the first few values of  r^k  * zk  using semi-open
c                            newton-cotes formulae
c-----------------------------------------------------------------------
      zk(1) = zero
      do 7 ii = 2,4
         sum = zero
         do 6 kk = 2,5
            sum = sum+cnc5c(kk,ii)*temp(kk)
    6    continue
         zk(ii) = sum
    7 continue
c-------------------------------------------------------------------
c   compute remainder of r^k  * zk: march out to mtp+3; use closed
c                       newton-cotes formula
c-------------------------------------------------------------------
      do 8 ii = 5,mtpp3
         zk(ii) = zk(ii-4)+c1*(temp(ii-4)+temp(ii  ))
     :                    +c2*(temp(ii-3)+temp(ii-1))
     :                    +c3* temp(ii-2)
    8 continue
c------------------------------------------------------------------------------
c determine the asymptotic value of  r^k * z^(k) correction to z^(0): in the 
c      manner of  c froese fischer,the hartree-fock method for atoms, 
c                  john wiley & sons,new york, 1977, p 235.
c------------------------------------------------------------------------------
      if (k .eq. 0) then
         if (i .eq. j) then
            zklim = one
         else
            zklim = zero
         endif
         do 10 kk = mtpp3,mtp,-1
            dif = zk(kk)-zklim
            if (abs (dif) .gt. accy) then
               do 9 ii = kk,2,-4
                  zk(ii) = zk(ii)-dif
    9          continue
            endif
   10    continue
      else
         zklim = zk(mtpp3)
      endif
c------------------------------------------------------------------------------
c                tabulate  zk  for entire internal lgrid
c------------------------------------------------------------------------------
      if (k .eq. 0) then
         do 11 ii = mtpp4,n
            zk(ii) = zklim
   11    continue
      else
         do 12 ii = 2,mtpp3
            zk(ii) = zk(ii)/rttk(ii)
   12    continue
         do 13 ii = mtpp4,n
            zk(ii) = zklim/rttk(ii)
   13    continue
      endif
c------------------------------------------------------------------------------
c                       start array wk / r**(k+1)
c------------------------------------------------------------------------------
      np4 = n+4
      do 14 ii = np4,mtpp1,-1
         wk(ii) = zero
   14 continue
c------------------------------------------------------------------------------
c                   compute  r^(k+1) and store in rttk
c------------------------------------------------------------------------------
      if (k .gt. 0) then
         do 15 ii = 2,n
            rttk(ii) = rttk(ii)*r(ii)
   15    continue
      endif
c------------------------------------------------------------------------------
c        fill array temp with rhop / r**(k+1) ; set temp(1) = zero
c                         to avoid 0/0 case
c------------------------------------------------------------------------------
      temp(1) = zero
      if (k .eq. 0) then
         do 16 ii = 2,mtp
            temp(ii) = rhop(ii)/r(ii)
   16    continue
      else
         do 17 ii = 2,mtp
            temp(ii) = rhop(ii)/rttk(ii)
   17    continue
      endif
c------------------------------------------------------------------------------
c         compute remainder of wk / r**(k+1): march in to the origin
c------------------------------------------------------------------------------
      do 18 ii = mtp,2,-1
         wk(ii) = wk(ii+4)+c1*(temp(ii  )+temp(ii+4))
     :                    +c2*(temp(ii+1)+temp(ii+3))
     :                    +c3*(temp(ii+2))
   18 continue
      wk(1) = zero
c------------------------------------------------------------------------------
c                               compute wk
c------------------------------------------------------------------------------
      if (k .eq. 0) then
         do 19 ii = 2,mtp
            wk(ii) = wk(ii)*r(ii)
   19    continue
      else
         do 20 ii = 2,mtp
            wk(ii) = wk(ii)*rttk(ii)
   20    continue
      endif
c------------------------------------------------------------------------------
c                           assemble solution
c------------------------------------------------------------------------------
      tmp1(1) = zero
      do 21 ii = 2,n
         tmp1(ii) = zk(ii)+wk(ii)
   21 continue
*
      end
c---------------------------------------------------------------------------
      subroutine yzk_iso1(k,i,j)
c--------------------------------------------------------------------------
      implicit real*8(a-h, o-z)
      include 'mpif.h'
      dimension rhop(MN),rttk(MN),wk(MN),temp(MN),zk(MN)
      common /cons/zero,half,tenth,one,two,three,ten
     :  /grid/steplin,r(MN),rp(MN),rpor(MN),rnt,h,hp,n
     :  /ncc/c1,c2,c3,c4,cnc5c(5,5),c5num(5,5),pi,c
     :  /dumm/tmp1(MN),tmp2(MN),tmp3(MN),tmp4(MN)
     :  /wave/pf(MN,MNBAS),qf(MN,MNBAS),mf(MNBAS)
     :  /basis/gl(MN,MNBAS),gs(MN,MNBAS),alpha(MNBAS,MNS),mg(MNBAS)
      accy=1.0d-12
c----------------------------------------------------------------------
c   for  k > 0  compute  r^k   and store in  rttk
c----------------------------------------------------------------------
      if (k .gt. 0) then
         do 1 ii = 2,n
            rttk(ii) = r(ii)**k
    1    continue
      endif
c---------------------------------------------------------------------
c     determine maximum tabulation point as location beyond which rhop (see 
c comment statements below) would be zero; determine other important locations
c--------------------------------------------------------------------
      mtp = min(mf(i),mf(j))
      mtpp1 = mtp+1
      mtpp3 = mtp+3
      mtpp4 = mtp+4
c--------------------------------------------------------------------
c   compute rp(s)*(p_i (s)*p_j (s)+q_i (s)*q_j (s)) and store in rhop
c--------------------------------------------------------------------
      do 2 ii = 2,mtp
         rhop(ii) = rp(ii)*(pf(ii,i)*pf(ii,j)+qf(ii,i)*qf(ii,j))
    2 continue
c--------------------------------------------------------------------
c                   fill array temp with r**k * rhop
c--------------------------------------------------------------------
      temp(1) = zero
      if (k .eq. 0) then
         do 3 ii = 2,mtp
            temp(ii) = rhop(ii)
    3    continue
      else
         do 4 ii = 2,mtp
            temp(ii) = rttk(ii)*rhop(ii)
    4    continue
      endif
c---------------------------------------------------------------------
c                 set an additional four points to zero
c---------------------------------------------------------------------
      do 5 ii = mtpp1,mtpp4
         temp(ii) = zero
    5 continue
c-----------------------------------------------------------------------
c           compute the first few values of  r^k  * zk  using semi-open
c                            newton-cotes formulae
c-----------------------------------------------------------------------
      zk(1) = zero
      do 7 ii = 2,4
         sum = zero
         do 6 kk = 2,5
            sum = sum+cnc5c(kk,ii)*temp(kk)
    6    continue
         zk(ii) = sum
    7 continue
c-------------------------------------------------------------------
c   compute remainder of r^k  * zk: march out to mtp+3; use closed
c                       newton-cotes formula
c-------------------------------------------------------------------
      do 8 ii = 5,mtpp3
         zk(ii) = zk(ii-4)+c1*(temp(ii-4)+temp(ii  ))
     :                    +c2*(temp(ii-3)+temp(ii-1))
     :                    +c3* temp(ii-2)
    8 continue
c------------------------------------------------------------------------------
c determine the asymptotic value of  r^k * z^(k) correction to z^(0): in the 
c      manner of  c froese fischer,the hartree-fock method for atoms, 
c                  john wiley & sons,new york, 1977, p 235.
c------------------------------------------------------------------------------
      if (k .eq. 0) then
         if (i .eq. j) then
            zklim = one
         else
            zklim = zero
         endif
         do 10 kk = mtpp3,mtp,-1
            dif = zk(kk)-zklim
            if (abs (dif) .gt. accy) then
               do 9 ii = kk,2,-4
                  zk(ii) = zk(ii)-dif
    9          continue
            endif
   10    continue
      else
         zklim = zk(mtpp3)
      endif
c------------------------------------------------------------------------------
c                tabulate  zk  for entire internal lgrid
c------------------------------------------------------------------------------
      if (k .eq. 0) then
         do 11 ii = mtpp4,n
            zk(ii) = zklim
   11    continue
      else
         do 12 ii = 2,mtpp3
            zk(ii) = zk(ii)/rttk(ii)
   12    continue
         do 13 ii = mtpp4,n
            zk(ii) = zklim/rttk(ii)
   13    continue
      endif
c------------------------------------------------------------------------------
c                       start array wk / r**(k+1)
c------------------------------------------------------------------------------
      np4 = n+4
      do 14 ii = np4,mtpp1,-1
         wk(ii) = zero
   14 continue
c------------------------------------------------------------------------------
c                   compute  r^(k+1) and store in rttk
c------------------------------------------------------------------------------
      if (k .gt. 0) then
         do 15 ii = 2,n
            rttk(ii) = rttk(ii)*r(ii)
   15    continue
      endif
c------------------------------------------------------------------------------
c        fill array temp with rhop / r**(k+1) ; set temp(1) = zero
c                         to avoid 0/0 case
c------------------------------------------------------------------------------
      temp(1) = zero
      if (k .eq. 0) then
         do 16 ii = 2,mtp
            temp(ii) = rhop(ii)/r(ii)
   16    continue
      else
         do 17 ii = 2,mtp
            temp(ii) = rhop(ii)/rttk(ii)
   17    continue
      endif
c------------------------------------------------------------------------------
c         compute remainder of wk / r**(k+1): march in to the origin
c------------------------------------------------------------------------------
      do 18 ii = mtp,2,-1
         wk(ii) = wk(ii+4)+c1*(temp(ii  )+temp(ii+4))
     :                    +c2*(temp(ii+1)+temp(ii+3))
     :                    +c3*(temp(ii+2))
   18 continue
      wk(1) = zero
c------------------------------------------------------------------------------
c                               compute wk
c------------------------------------------------------------------------------
      if (k .eq. 0) then
         do 19 ii = 2,mtp
            wk(ii) = wk(ii)*r(ii)
   19    continue
      else
         do 20 ii = 2,mtp
            wk(ii) = wk(ii)*rttk(ii)
   20    continue
      endif
c------------------------------------------------------------------------------
c                           assemble solution
c------------------------------------------------------------------------------
      tmp1(1) = zero
      do 21 ii = 2,n
         tmp1(ii) = zk(ii)+wk(ii)
   21 continue
*
      end

C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine yzk2(k,i,j)
c------------------------------------------------------------------------------
c   this subroutine evaluates hartree y- and z-functions:              *
c               (k)            (k)           (k)                       *
c              y#   (i,j;r) =  z#   (i,j;r) + w#   (i,j;r)             *
c                                                                      * 
c   where # can be p or q and                                          *
c    (k)                                                               *
c   z#   (i,j;r) =  i ( (s/r)   # (s)*# (s) ; 0 - r )                  *
c                                i     j                               *
c    (k)                    k+1                                        *
c   w#   (i,j;r) =  i ( (r/s)    # (s)*# (s) ; r - infinity)           *
c                                 i     j                              *
c                                                                      *
c   where  i ( g(r,s) ; range )  denotes the integral of g(r,s) over   *
c   range  in  s .  the yp-function is tabulated in  common/tatb/  in  *
c   array  ta , the yq-function in array tb . for the direct part the  *
c   two araays yp and yq are summed and tabulated in array ta.         *
c------------------------------------------------------------------------------
      implicit real*8(a-h, o-z)
      include 'mpif.h'
      dimension rhopp(MN),wkp(MN+5),tempp(MN),zkp(MN),
     :          rhopq(MN),wkq(MN+5),tempq(MN),zkq(MN),rttk(MN)
      common /cons/zero,half,tenth,one,two,three,ten
     :  /grid/steplin,r(MN),rp(MN),rpor(MN),rnt,h,hp,n
     :  /ncc/c1,c2,c3,c4,cnc5c(5,5),c5num(5,5),pi,c
     :  /dumm/tmp1(MN),tmp2(MN),tmp3(MN),tmp4(MN)
     :  /wave/pf(MN,MNBAS),qf(MN,MNBAS),mf(MNBAS)
     :  /basis/gl(MN,MNBAS),gs(MN,MNBAS),alpha(MNBAS,MNS),mg(MNBAS)
c------------------------------------------------------------------------------
c           for  k > 0  compute  r^k   and store in  rttk
c------------------------------------------------------------------------------
      if (k .gt. 0) then
        do 1 ii = 2,n
          rttk(ii) = r(ii)**k
    1   continue
      endif
c------------------------------------------------------------------------------
c   determine maximum tabulation point as location beyond which rhop  (see 
c     comment statements below) would be zero; determine other important 
c                               locations
c------------------------------------------------------------------------------
      mtp = mf(i)
      mtpp1 = mtp+1
      mtpp3 = mtp+3
      mtpp4 = mtp+4
c------------------------------------------------------------------------------
c   compute rp(s)* p_i (s)*gl_j (s) and store in rhopp
c                  q_i     gs_j
c------------------------------------------------------------------------------
      do 2 ii = 2,mtp
         rhopp(ii) = rp(ii)*pf(ii,i)*gl(ii,j)
         rhopq(ii) = rp(ii)*qf(ii,i)*gs(ii,j)
    2 continue
c------------------------------------------------------------------------------
c                   fill array tempp with r**k * rhopp
c                              tempq             rhopq
c------------------------------------------------------------------------------
      tempp(1) = zero
      tempq(1) = zero
      if (k .eq. 0) then
        do 3 ii = 2,mtp
          tempp(ii) = rhopp(ii)
          tempq(ii) = rhopq(ii)
    3   continue
      else
        do 4 ii = 2,mtp
          tempp(ii) = rttk(ii)*rhopp(ii)
          tempq(ii) = rttk(ii)*rhopq(ii)
    4   continue
      endif
c------------------------------------------------------------------------------
c                   set an additional four points to zero
c------------------------------------------------------------------------------
      do 5 ii = mtpp1,mtpp4
        tempp(ii) = zero
        tempq(ii) = zero
    5 continue
c------------------------------------------------------------------------------
c compute first few values of  r^k*zk using semi-open newton-cotes formulae
c------------------------------------------------------------------------------
      zkp(1) = zero
      zkq(1) = zero
      do 7 ii = 2,4
        sump = zero
        sumq = zero
        do 6 kk = 2,5
          sump = sump+cnc5c(kk,ii)*tempp(kk)
          sumq = sumq+cnc5c(kk,ii)*tempq(kk)
    6   continue
        zkp(ii) = sump
        zkq(ii) = sumq
    7 continue
c------------------------------------------------------------------------------
c        compute remainder of r^k*zk: march out to mtp+3; use closed
c                         newton-cotes formula
c------------------------------------------------------------------------------
      do 8 ii = 5,mtpp3
         zkp(ii) = zkp(ii-4)+c1*(tempp(ii-4)+tempp(ii  ))
     :                      +c2*(tempp(ii-3)+tempp(ii-1))
     :                      +c3* tempp(ii-2)
         zkq(ii) = zkq(ii-4)+c1*(tempq(ii-4)+tempq(ii  ))
     :                      +c2*(tempq(ii-3)+tempq(ii-1))
     :                      +c3* tempq(ii-2)
    8 continue
c------------------------------------------------------------------------------
c   determine the asymptotic value of  r^k*z^(k) correction to  z^(k)   
c         : in the manner of  c froese fischer,the hartree-fock method for 
c            atoms, john wiley & sons, new york, 1977, p 235.
c------------------------------------------------------------------------------
      if (k .eq. 0) then
        zklimp = zero
        zklimq = zero
      else
        zklimp = zkp(mtpp3)
        zklimq = zkq(mtpp3)
      endif
c------------------------------------------------------------------------------
c                     tabulate  zk  for entire internal lgrid
c------------------------------------------------------------------------------
      if (k .eq. 0) then
        do 12 ii = mtpp4,n
          zkp(ii) = zklimp
          zkq(ii) = zklimq
   12   continue
      else
        do 13 ii = 2,mtpp3
          zkp(ii) = zkp(ii)/rttk(ii)
          zkq(ii) = zkq(ii)/rttk(ii)
   13   continue
        do 14 ii = mtpp4,n
          zkp(ii) = zklimp/rttk(ii)
          zkq(ii) = zklimq/rttk(ii)
   14   continue
      endif
c------------------------------------------------------------------------------
c                          start array wk / r**(k+1)
c------------------------------------------------------------------------------
         np4 = n+4
         do 15 ii = np4,mtpp1,-1
            wkp(ii) = zero
            wkq(ii) = zero
   15    continue
c------------------------------------------------------------------------------
c               compute  r^(k+1)  and store in rttk
c------------------------------------------------------------------------------
         if (k .gt. 0) then
         do 16 ii = 2,n
            rttk(ii) = rttk(ii)*r(ii)
   16    continue
         endif
c------------------------------------------------------------------------------
c        fill array temp with rhop / r**(k+1) ; set temp(1) = zero
c                          to avoid 0/0 case
c------------------------------------------------------------------------------
         tempp(1) = zero
         tempq(1) = zero
         if (k .eq. 0) then
         do 17 ii = 2,mtp
            tempp(ii) = rhopp(ii)/r(ii)             
            tempq(ii) = rhopq(ii)/r(ii)
   17    continue
         else
         do 18 ii = 2,mtp
            tempp(ii) = rhopp(ii)/rttk(ii)
            tempq(ii) = rhopq(ii)/rttk(ii)
   18    continue
         endif
c------------------------------------------------------------------------------
c         compute remainder of wk / r**(k+1): march in to the origin
c------------------------------------------------------------------------------
      do 19 ii = mtp,2,-1
         wkp(ii) = wkp(ii+4)+c1*(tempp(ii  )+tempp(ii+4))
     :                    +c2*(tempp(ii+1)+tempp(ii+3))
     :                    +c3*(tempp(ii+2))
         wkq(ii) = wkq(ii+4)+c1*(tempq(ii  )+tempq(ii+4))
     :                    +c2*(tempq(ii+1)+tempq(ii+3))
     :                    +c3*(tempq(ii+2))
   19 continue
      wkp(1) = zero
      wkq(1) = zero
c------------------------------------------------------------------------------
c                                 compute wk
c------------------------------------------------------------------------------
      if (k .eq. 0) then
         do 20 ii = 2,mtp
            wkp(ii) = wkp(ii)*r(ii)
            wkq(ii) = wkq(ii)*r(ii)
   20    continue
      else
         do 21 ii = 2,mtp
            wkp(ii) = wkp(ii)*rttk(ii)
            wkq(ii) = wkq(ii)*rttk(ii)
   21    continue
      endif
c-----------------------------------------------------------
c                             assemble solution
c----------------------------------------------------------
      tmp1(1) = zero
      tmp2(1) = zero
      do 22 ii = 2,n
        tmp1(ii) = zkp(ii)+wkp(ii)
        tmp2(ii) = zkq(ii)+wkq(ii) 
   22 continue
      return
*
      end
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine scfiter(nbasis)
c------------------------------------------------------------------
c     This performs the scf iteration.
c------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      include 'mpif.h'
      logical set
      character*3 nh_i
      real*8 df_mat,swap
      dimension can(MNB,MNB)
      dimension  drs(MNB,MNB),sa(MNB),df_mat(MNB,MNB)
      common/basis/gl(MN,MNBAS),gs(MN,MNBAS),alpha(MNBAS,MNS),mg(MNBAS)
      common /grid/steplin,r(MN),rp(MN),rpor(MN),rnt,h,hp,n
     :  /ncc/c1,c2,c3,c4,cnc5c(5,5),c5num(5,5),pi,c
     :  /cons/zero,half,tenth,one,two,three,ten
     :  /overlap/sss(MNB,MNB),sll(MNB,MNB)
     :  /adat/occj(MNOCC),occl(MNOCC),unocl(MNS),unocj(MNS)
     :  /quan1/kap(MNS),nak(MNOCC),np(MNOCC),nsy(MNS),nkk(MNS)
     :  /diag/eigv(MNB,MNB),eng(MNB)
     :  /datt/nsym,nocc,nbas(MNS),ind
     :  /single/df_single(MNB,MNB)
     :  /double/df_two(MNB,MNB)
     :  /wave/pf(MN,MNBAS),qf(MN,MNBAS),mf(MNBAS)
     :  /conv/crit,maxit,npower
     :  /print/jprint,iout
     :  /jv/occ(MNOCC),orbe(MNB),iq(MNB)
     :  /skip/nskip(MNS)
     :  /dens/p(MNB,MNB)
     :  /ovlp/emt(MNB,MNB)
     :  /info1/occk(MNOCC,MNS)
     :  /info/nocorb(MNOCC),nucorb(MNB)
     :  /vector/vec(MNB,MNB),oe(MNB,MNS)
     :  /info3/df_t(MNB,MNB)

      dimension a(MNB,MNB),st(MNB,MNB),b(MNB,MNB)
      dimension pold(MNB,MNB),stemp(MNB,MNB)
      dimension nh_i(MNS),scr(MNB,MNB)

      data nh_i/' S ',' P-',' P ',' D-',' D ',' F-',' F ',
     :                     ' G-',' G ','H-','H'/

c--------------------------------------------------------------------
c       initialize arrays and constants required for the computation.
c--------------------------------------------------------------------
      tcsqr=-two*c*c
      set = .true.
      iter=0
      en=zero
      eold=zero
      delta=zero

      do ia=1,MNB
      sa(ia)=zero
      do ib=1,MNB
      can(ib,ia)=zero
      drs(ib,ia)=zero
      p(ia,ib)=zero
      enddo
      enddo

      mbasis=2*nbasis
      write(STDOUT,199)maxit,npower
      write(STDOUT,200)    
      write(STDOUT,201)iter,delta
      write(*,201)iter,delta
 100  iter=iter+1
      do isym=1,nsym
      ndim=nbas(isym)
      ndim2=2*ndim
 
      do ia=1,ndim2
      ja=ia+2*nskip(isym)
      do ib=1,ndim2
      jb=ib+2*nskip(isym)
      stemp(ia,ib)=emt(ja,jb)
      enddo
      enddo
c***********************************************************
c      matout is wrirting the energy in the output         *
c***********************************************************   
     
      if(jprint.ne.0)call matout(stemp,MNB,ndim2,ndim2,ndim2,4hS   )

      do i=1,ndim2
      do j=1,ndim2
      drs(i,j)=stemp(i,j)
      enddo
      enddo

c**************************************************************
c    'stemp' or 'drs' is basically the overlap matrix         *
c  so at first we compute the eigen value of overlap matrix 
c  'sa' is the eigen value and then use calay-hamilton theorem to
calculate the inverse of overlap matrix                       *
c                  our equation is
c      < |H| > = <overlap mat> E or, E= < |H| > (<overlap mat>^-1)  
c**************************************************************

      call eigen (drs,scr,ndim2,sa,ndim2,bnorm,MNB)
c-------------------------------------------------------------------
c       canonical transformation
c    calculating the inverse of overlap mat
c------------------------------------------------------------------
      do ia=1,ndim2
      do ib=1,ndim2
      if(dabs(sa(ib)).gt.TINY) then
      can(ia,ib)=drs(ia,ib)/dsqrt(dabs(sa(ib)))
      st(ib,ia)=can(ia,ib)
      else
      write(*,*)sa(ib)
      endif
      enddo
      enddo

      if(jprint.ne.0)call matout(st,MNB,ndim2,ndim2,ndim2,4hXT  )

      if(jprint.ne.0)
     :          call matout(df_single,MNB,mbasis,mbasis,mbasis,4hH   )

      if(iter.ne.1)call twoe(isym,nbasis)

      do ic=1,ndim2
      jc=ic+2*nskip(isym)
      do id=1,ndim2
      jd=id+2*nskip(isym)
      if(iter.ne.1)then
      a(ic,id)=df_single(jc,jd)+df_two(jc,jd)
      else
      a(ic,id)=df_single(jc,jd)
      endif

      enddo
      enddo

      if(jprint.ne.0)call matout(a,MNB,ndim2,ndim2,ndim2,4hFOCK )

      if(jprint.ne.0)
     :call matout(df_two,MNB,mbasis,mbasis,mbasis,4hVHF )

      call mult(a,can,b,MNB,ndim2)
      call mult(st,b,df_mat,MNB,ndim2)
      if(jprint.ne.0)
     :      call matout(df_mat,MNB,ndim2,ndim2,ndim2,6hFRIME )
      call eigen (df_mat,scr,ndim2,eng,ndim2,bnorm,MNB)
      do i=1,ndim2
      do j=1,ndim2
      drs(i,j)=df_mat(i,j)
      enddo
      enddo

      if(jprint.ne.0)
     : call matout(drs,MNB,ndim2,ndim2,ndim2,4hC   )
      call mult(can,drs,eigv,MNB,ndim2)
  
      loc = 1
      last=ndim-1
      dowhile(loc.le.last)
      samin = eng(loc)
      locus = loc
      do  ia = loc,ndim
      if (eng(ia).lt.samin) then
      samin = eng(ia)
      locus = ia
      endif
      enddo
      if (locus .gt. loc) then
      swap = eng(loc)
      eng(loc) = eng(locus)
      eng(locus) = swap
      do  ib = 1,ndim2
      swap = eigv(ib,loc)
      eigv(ib,loc) = eigv(ib,locus)
      eigv(ib,locus) = swap
      enddo
      endif
      loc = loc + 1
      enddo
c-----------------------------------------------------------------
      do ia=1,ndim2
      do ib=1,ndim2
      b(ia,ib)=zero
      enddo
      b(ia,ia)=eng(ia)
      oe(ia,isym)=eng(ia)
      enddo

      if(jprint.ne.0)call matout(b,MNB,ndim2,ndim2,ndim2,4hE   )

      if(jprint.ne.0)call matout(eigv,MNB,ndim2,ndim2,ndim2,4hCP  )

      do ia=1,ndim2
      ja=ia+2*nskip(isym)
      do ib=1,ndim2
      jb=ib+2*nskip(isym)
      pold(ja,jb)=p(ja,jb)
      p(ja,jb)=0.0d0
      nocc=nocorb(isym)
      if(nocc.ne.0)then
      do ibas=1,nocc
      twojp1=two*occ(isym)+one
      p(ja,jb)=p(ja,jb)+twojp1*eigv(ia,ibas)*eigv(ib,ibas)
      enddo
      endif
      vec(ja,jb)=eigv(ia,ib)
      enddo
      enddo

c   initialise 
       do ia =1,nbas(isym)
       ka=ia+nskip(isym)
       ja=ia+2*nskip(isym)
       do ic=1,n
       pf(ic,ka)=0.0d0
       qf(ic,ka)=0.0d0
       enddo
       ja=ia+2*nskip(isym)
       do ib =1,nbas(isym)
       jb=ib+2*nskip(isym)
       kb=ib+nskip(isym)
       do  ic=1,n
       pf(ic,ka)=pf(ic,ka)+vec(jb,ja)*gl(ic,kb)
       qf(ic,ka)=qf(ic,ka)+vec(jb+nbas(isym),ja)*gs(ic,kb)
       enddo
       enddo
       enddo

       enddo

      if(jprint.ne.0)call matout(p,MNB,mbasis,mbasis,mbasis,4hP   )

c     calculate delta

      delta=0.0d0
      do i=1,mbasis
      do j=1,mbasis
      delta=delta+(p(i,j)-pold(i,j))**2
      enddo
      enddo
      delta=dsqrt(delta/4.0d0)
      write(STDOUT,201)iter,delta
      write(*,201)iter,delta
      if(delta.lt.crit)go to 160
      if(iter.lt.maxit)go to 100
      write(STDOUT,202)
      stop
  160 continue

      call integral(nbasis)
      eold=en
      en=zero
c     print*,'density=',p(3,3),p(3,33),p(33,3),p(33,33)
      do i=1,mbasis
      do j=1,mbasis
      en=en+0.5d0*p(i,j)*(two*df_single(i,j)+df_two(i,j))
      enddo
      enddo
      
      write(STDOUT,203)iter
      write(STDOUT,204)en
      write(STDOUT,205)
      do jsym=1,nsym
      do ibas=1,nbas(jsym)
      write(STDOUT,206)ibas,nh_i(jsym),oe(ibas,jsym),
     :                 oe(ibas+nbas(jsym),jsym)
      enddo
      enddo

      return
 199  format(/,1x,'CONVERGENCE DATA',/,
     1        ' MAXIMUM NO. OF ITERATION=',i6,/,
     2        ' CONVERGENCE CRITERION   = 1.0D-',i2,//)
 200  format(/,1x,'CYCLE      DENSITY CONV',/)
 201  format(i3,1x,(2x,d20.7))
 202  format(///,4x,33hSCF FAILS TO CONVERGES AT CYCLE , i4)
 203  format(///,4x,22hSCF CONVERGES AT CYCLE, i4)
 204  format(///,4x,20hELECTRONIC ENERGY=  ,d20.11)
 205  format(///,4x,32h ORBITAL ENERGIES (+ve and -ve) ,i4,///)
 206  format(/,4x,I2,A,2D20.11)
      end            
************************************************************************
*                                                                      *
      double precision function rint (i,j,k)
*                                                                      *
*                                                                      *
*   the value of this rint is an approximation to:                     *
*                                                                      *
*              k                                                       *
*         i ( r  *  ( p (r)*p (r) + q (r)*q (r) ; 0 to infinity)       *
*                      i     j       i     j                           *
*                                                                      *
*   where   i ( g(r) ; range )  denotes  the  integral  of g(r) over   *
*   range.                                                             *
*                                                                      *
*   subroutines called: quad.                                          *
************************************************************************
*
      implicit real*8 (a-h,o-z)
      include 'mpif.h'
      common/cons/zero,half,tenth,one,two,three,ten
     :  /grid/steplin,r(MN),rp(MN),rpor(MN),rnt,h,hp,n
     :  /wave/pf(MN,MNBAS),qf(MN,MNBAS),mf(MNBAS)
      dimension ta(MN)
*
      mtp = n
*
*   tabulate integrand as required for subroutine quad
*   value at first tabulation point is arbitrary
*
      ta(1) = zero
      do 1 l = 2,mtp
      ta(l)=(r(l)**k)*(pf(l,i)*pf(l,j)+qf(l,i)*qf(l,j))*rp(l)
    1 continue
*
*   perform integration
*
      call quad (ta,mtp,result)
      rint = result
*
      return
      end
c
      subroutine setgrd
      implicit real*8 (a-h,o-z)
      include 'mpif.h'
      common /grid/steplin,r(MN),rp(MN),rpor(MN),rnt,h,hp,n
      n = MN-10
      mtp=n
      h = 5.0d-02
      rnt = 2.0d-06
      r(1) = 0.0d0
      rp(1) = rnt
      eph = exp (h)
      ett = 1.d0
*   set up the arrays r, rp, rpor
      do i = 2,MN
      ett = eph*ett
      ettm1 = ett-1.0d0
      r(i) = rnt*ettm1
      rp(i) = rnt*ett
      rpor(i) = ett/ettm1
      enddo
      return
      end

