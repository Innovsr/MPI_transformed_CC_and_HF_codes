CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC                                                                  CC
CC               IBRA = 0   for only couloumb interaction           CC
CC                    = 1   for all parity channel with coulomb int.CC
CC                    = 2   for both coulomb and breit interaction  CC
CC                    = 3   for only breit interaction              CC
CC                                                                  CC
CC                                                                  CC
CC               IOPT = 0   for CCD                                 CC
CC                    = 1   for CCSD                                CC
CC                                                                  CC
CC             NONLIN = 0   for linear CC                           CC
CC                    = 1   for non-linear CC                       CC
CC                                                                  CC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC                                                                  CC
CC      INPUT FILE NAME: input                                      CC
CC    READ:                                                         CC
CC   1. a) Total no. of symm.  b) Active orbitals for each symm.    CC
CC   2. a) Active occupied orbitals                                 CC
CC   3. a) iopt  b) nonlin  c) isolv d) nprint  e) iall             CC
CC                                                                  CC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


C CPP VARIABLES
C
C DIMENSION VARIABLES
#define     NHO      23
#define     NPO      82
#define     MXV       9
#define    MDIM   2000000

c
#define      MN     750
#define   MNSYM      11
#define     MNS       6
#define    MXVR    (MXV+1)/2
#define   MNBAS    NHO+NPO
#define   MNOCC    NHO+NPO
#define   MNEXC    NHO+NPO
#define   NTEMP   NPO*(NPO+1)/2
#define   MULLI   NTEMP*(NTEMP+1)/2
#define   IMSHPPH  (NPO**2)*(NHO**2)*(MXV+1)
#define   IMSHHHH  (NHO**4)*(MXV+1)


#define    NEV    57               !Total no. of even parity orbitals
#define    NEE    NEV*(NEV+1)/2
C#define          NEEEE=NEE*(NEE+1)/2
#define    NEEEE  1367031

#define    NOD    48                   !Total no. of  odd parity orbitals
#define    NOO    NOD*(NOD+1)/2
C#define          NOOOO=NOO*(NOO+1)/2
#define    NOOOO  692076

#define    NOE    NEV*NOD
C#define   NOEOE  NOE*(NOE+1)/2
#define    NOEOE  3744216

c#define   N2INT   2*(NOOOO+NEEEE+NOEOE)
#define   N2INT   (NOOOO+NEEEE+NOEOE+2000000)

C CONSTANTS
#define  STDIN    5
#define STDOUT    6
#define  WFNIN    7
#define  NTWOE    8
#define NTFILE   16
#define NITMAX   30
#define  NPMAX   128
#define  SMALL   1.0d-6
#define MASTER   0
     
      implicit real*8(a-h,o-z)
      include 'mpif.h'
      integer drow,rowstart,rowend,extrarow
      integer error,np
      common/scratch/b(MDIM),ad(MDIM),baux(MDIM)
      common/scratch2/h2d(MDIM),raux(MDIM),ap(MDIM)                 
      common/iallparity/iall
      common/skip/nnskip(0:MXV,2),n2skip,nsing
c     local dimension
      integer dspl(NPMAX),delrow(NPMAX)
      dimension nbas(MNSYM),iocc(MNSYM)
      dimension told(MDIM),t(MDIM)

C     THIS IS THE MAIN ROUTINE



c     iproc=0
c     nproc = 1
c     error = 0

      call MPI_INIT(ierr)
      START_TIME=MPI_WTIME()
      call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)

      call daopen
      call symmcal
      call setgrd
      print*,'setgrd is over'

      call readinp(nsym,nbas,iocc,iopt,nonlin,isolv,nprint,itwoe)
      print*,'readinp is over'
      call setup(nsym,nbasis,nbas,iocc,iproc)
      print*,'setup is over'
      call symm(iopt,ntmax,isolv,nprint,iproc)

      print*,'symm is over'

      if(iproc.eq.0)then

      write(STDOUT,*)'<<<<<<<<<<<<PLEASE CHECK>>>>>>>>>>>'
      write(STDOUT,*)'<<**********WHETHER**************>>'
      write(STDOUT,*)'<<*********YOU WANT THIS*********>>'
      write(STDOUT,*)'<</                              />>'
      if(iall.eq.0)write(STDOUT,*)'<< Only Coulomb interaction      >>'
      if(iall.eq.1)write(STDOUT,*)'<< All Parity Channel with Coul. >>'
      if(iall.eq.2)write(STDOUT,*)'<< Both Coul. and Breit int.     >>'
      if(iall.eq.3)write(STDOUT,*)'<< Only Breit intearction        >>'
      if(iopt.eq.0)write(STDOUT,*)'<< CCD  is considered            >>'
      if(iopt.eq.1)write(STDOUT,*)'<< CCSD is considered            >>'
      if(nonlin.eq.0)write(STDOUT,*)'<<Linearlised CC                >>'
      if(nonlin.eq.1)write(STDOUT,*)'<<Non-linearlised CC            >>'
      write(STDOUT,*)'<</                              />>'
      write(STDOUT,*)'<<<<<<<<<<<<OOOOO>>>>>>>>>>>>>>>>>>'

      write(STDOUT,*)'____________________________________________'
      write(STDOUT,*)'____________________________________________'
      write(STDOUT,*)'isolv value =  ',isolv
      write(STDOUT,*)'nprint value = ',nprint 
      write(STDOUT,*)'Total no. single excitations = ',nsing
      write(STDOUT,*)'Total no. double excitations = ',ntmax-nsing
      write(STDOUT,*)'____________________________________________'
      write(STDOUT,*)'Total no. excitations =        ',ntmax
      write(STDOUT,*)'____________________________________________'

      endif

      call grasprd(nbasis,iproc)
      print*,'finish grasprd'
      call setqic
      print*,'finish setqic'
      call mInit() 
      print*,'finish mInit'

      call twoint(nbasis,itwoe)
      print*,'finish twoint'

c     call a matrix

      call amat(iopt,ntmax,isolv,iproc)

c     initial guess
      do i =1,ntmax
      t(i)=b(i)/ad(i)
      told(i)=t(i)
      enddo

c     compute second order correlation energy
      call cizec2(t,iopt,iproc)

      drow = ntmax/nproc
      do i = 1, nproc
      delrow(i) = drow
      end do
      extrarow = mod(ntmax,nproc)
      do i = 1, extrarow
      delrow(i) = delrow(i) + 1
      end do
      dspl(1) = 0
      do i = 2, nproc
      dspl(i) = dspl(i - 1) + delrow(i - 1)
      end do
      rowstart = dspl(iproc + 1) + 1
      drow = delrow(iproc + 1)
      rowend = rowstart + drow - 1

      if (isolv.eq.0) then
      call jacobi1(told,t,iopt,isolv,ntmax,nonlin,
     :            rowstart,rowend,iproc,nproc,dspl,delrow)

      if(nonlin.ne.0)then
      call jacobi2(told,t,iopt,isolv,ntmax,nonlin,
     :            rowstart,rowend,iproc,nproc,dspl,delrow)
      endif

      else

      call jacobi1(told,t,iopt,isolv,ntmax,nonlin,
     :            rowstart,rowend,iproc,nproc,dspl,delrow)

      if(nonlin.ne.0)then
      call jacobi2(told,t,iopt,isolv,ntmax,nonlin,
     :            rowstart,rowend,iproc,nproc,dspl,delrow)
      endif

      endif

      call daclose
      END_TIME=MPI_WTIME()
       if (iproc .eq. 0) then
        print *, 'Time spent (in seconds) is: ',END_TIME-START_TIME
        endif
       if(iproc.eq.0)then
       write(STDOUT,*)'Time spent (in seconds) is: ',
     :  END_TIME-START_TIME
      endif
       call MPI_FINALIZE(ierr)


c      call MPI_FINALIZE(error)

      stop
      end


      integer function Calc4Binary(i, j, k, l)
c----------------------------------------------------------------------
      implicit none
      include 'mpif.h'
      integer i, j, k, l
                       
      Calc4Binary = i*8 + j*4 + k*2 + l
      return
      end


      subroutine findk2(orba,orbb,orbp,orbq,ia,ib,ip,iq,
     $jeven,kkk)

C     THIS SUBROUTINE FINDS THE ALLOWED MULTIPOLE MOMENT
C     VARIBALES "ORBA, ORBB ETC. ARE J-VALUE OF THE A, B
C     ORBITAL, IA, IB ETC. ARE 'IA" VALUE. 
C     
C     ON EXIT KKK--> NO. OF ALLOWED MOMENT AND JEVEN
C     STORES THE MOMENT VALUE    

      implicit real*8 (a-h,o-z)
      include 'mpif.h'
      real*8 jj
      dimension ja(MXV),jb(MXV),jp(MXV),jq(MXV)
      dimension jeven(MXV+1)

      do i=1,MXV
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
      jj=dfloat(jtot)-2.0*(dfloat(jtot/2))
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
      jj=dfloat(jtot)-2.0*(dfloat(jtot/2))
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
      endif
      enddo
      enddo

      do i=1,ll
      do j=1,nn
      if(ja(i).eq.jq(j))then
      kkk=kkk+1
      jeven(kkk)=ja(i)
      endif
      enddo
      enddo

      do i=1,kk
      do j=1,mm
      if(jb(i).eq.jp(j))then
      kkk=kkk+1
      jeven(kkk)=jb(i)
      endif
      enddo
      enddo

      do i=1,kk
      do j=1,nn
      if(jb(i).eq.jq(j))then
      kkk=kkk+1
      jeven(kkk)=jb(i)
      endif
      enddo
      enddo
      return
      end
 
      subroutine init(a)

C     Initialize the ARRYA A
c     include 'mpif.h'

      real*8 a(MDIM)
      do i=1,MDIM
      a(i)=0
      enddo
      return
      end

      subroutine iloc(i,j,k,l,ij,kl)
c     pack two index to one index
      implicit integer (a-z)
      include 'mpif.h'
      ij=j*100+i
      kl=l*100+k
      return
      end

      subroutine sixj(rj1,rj2,rj3,rl1,rl2,rl3,s6j)
C
C     THIS SUBROUTINE COMPUTES SIXJ VALUE
C
      implicit real*8 (a-h,o-z)
      include 'mpif.h'
        common/factor/ fct(0:50),mfd
!comment function s6j(j1,j2,j3,l1,l2,l3) calculates the 6j symbol
!        of its arguments.  it returns an error message if the
!        sum of angular momenta involved in a triangular
!        relationship is not an integer;
!for all j1,j2,j3,l1,l2,l3 such that
!        fixp(2*j1) and fixp(2*j2) and fixp(2*j3) and
!        fixp(2*l1) and fixp(2*l2) and fixp(2*l3)
! let s6j(j1,j2,j3,l1,l2,l3) =
!    begin scalar ws6j,delprod,zmin,zmax,tria,trib,tric,trid,s;

        if(mfd.eq.33)goto 704
        mfd=33
        fct(0)=1.
        fct(1)=1.
        do 20 i=2,33
        fct(i)=fct(i-1)*i
20      continue

704    tria=tri6j(rj1,rj2,rj3)
       trib=tri6j(rl1,rj2,rl3)
       tric=tri6j(rj1,rl2,rl3)
       trid=tri6j(rl1,rl2,rj3)
c       type *,tria,trib,tric,trid
       if(tria.eq.-1..or.trib.eq.-1..or.tric.eq.-1..or.trid.eq.-1.)then
c      type *,'inconsistent arguments to 6j symbol'
c      type *,tria,trib,tric,trid
c      type *,rj1,rj2,rj3
c      type *,rl1,rl2,rl3
        return
        endif
       if(tria.eq.0..or.trib.eq.0..or.tric.eq.0..or.trid.eq.0.)then
c         type *,'triangle condition'
          s6j=0
        return
        endif
!    comment apply formula in edmonds [(6.3.7) p.99]

       delprod=tridel(rj1,rj2,rj3)*tridel(rj1,rl2,rl3)*
     *         tridel(rl1,rj2,rl3)*tridel(rl1,rl2,rj3)
c       type *,'delprod=',delprod
c       type *,'args to min',(rj1+rj2+rl1+rl2),(rj2+rj3+rl2+rl3),
c     *                       (rj3+rj1+rl3+rl1)
       izmin=dmax1((rj1+rj2+rj3),(rj1+rl2+rl3),(rl1+rj2+rl3)
     *           ,(rl1+rl2+rj3))
       izmax=dmin1((rj1+rj2+rl1+rl2),(rj2+rj3+rl2+rl3)
     *           ,(rj3+rj1+rl3+rl1))
        ws6j=0
c       type *,'izmin,izmax=',izmin,izmax
        do 10 iz=izmin,izmax
        z=iz
         ws6j=ws6j+((-1)**iz)*fct(int(z+1.))/(
     *           fct(int(z-rj1-rj2-rj3))*
     *           fct(int(z-rj1-rl2-rl3))*
     *           fct(int(z-rl1-rj2-rl3))*
     *           fct(int(z-rl1-rl2-rj3))*
     *           fct(int(rj1+rj2+rl1+rl2-z))*
     *           fct(int(rj2+rj3+rl2+rl3-z))*
     *           fct(int(rj3+rj1+rl3+rl1-z)))
10      continue
       s6j=delprod*ws6j
       return
       end

        function tridel(rm1,rm2,rm3)
        implicit real*8 (a-h,o-z)
        include 'mpif.h'
        common/factor/ fct(0:50),mfd

c       type *,'+++',rm1,rm2,rm3
!comment function tridel evaluates the delta symbol defined on
!        p.99 of edmonds;
!for all m1,m2,m3 such that fixp(2*m1) and fixp(2*m2) and fixp(2*m3)
c       type *,m1+m2-m3,m1-m2+m3,-m1+m2+m3,m1+m2+m3+1

        i1=rm1+rm2-rm3
        i2=rm1-rm2+rm3
        i3=-rm1+rm2+rm3
        i4=rm1+rm2+rm3+1

        tridel=sqrt(fct(i1)*fct(i2)*fct(i3)/fct(i4))

c       type *,'***',rm1,rm2,rm3

        return
        end
        function tri6j(rj1,rj2,rj3)
        implicit real*8 (a-h,o-z)
        include 'mpif.h'
        sumj=rj1+rj2+rj3
        isumj=int(sumj)
        if(abs(sumj-isumj).gt.1e-3)then !if j1+j2+j3 is not an integer
          tri6j=-1                      !then tri6j=-1
        goto 999
        endif
        if((2.*abs(rj1-rj2).le.2.*rj3).and.
     *        (2.*(rj1+rj2).ge.2.*rj3))then
          tri6j=1                               !triangle condition is satisfied
        goto 999
        endif

        tri6j=0                         !j1+j2+j3 is an integer
c       type *,rj1,rj2,rj3
999     return
        end
      
      function dk(cj1,cj11,cj2,cj22,cm1,ak)

      implicit real*8(a-h,o-z)
        include 'mpif.h'
      real*8 cj1,cj2,x,y,w3j1,w3j2,cj11,cj22
      integer j1,j2


         j1=idint(2.0*cj1)
         j2=idint(2.0*cj2)
         x=dsqrt(dfloat(j1)+1)
         y=dsqrt(dfloat(j2)+1)

* loop over the m1(j1+1) and m2(j2+1) values. j1=2*cj1,j2=2*cj2

* w3j value is calculated here

              w3j1=dr(cj1,cj2,ak,0.5d0,-0.5d0,0.0d0)
              w3j2=dr(cj1,cj2,ak,-cj11,cj22,cm1)

* phase part is calculated here

              cpower1=j2-j1
              ipower1=idint(cpower1)
              cpower2=0.5+cj22
              ipower2=idint(cpower2)

* the dk coefficient part is calculted here

              dk=x*y*(-1)**(ipower1)*(-1)**(ipower2)*w3j1*w3j2
       return
      end
*-----------------------------------------------------------------
*    programe to calculate the 3j symbols
*-----------------------------------------------------------------
      function dr(cj1,cj2,j3,cm1,cm2,m3)
      implicit real*8(a-h,o-z)
        include 'mpif.h'
      real*8 m3,j3
      ddr=0.d0
      if(abs(cm1+cm2+m3).gt.0.1d0)go to 80
      if(abs(cm1).gt.cj1+0.1d0)go to 80
      if(abs(cm2).gt.cj2+0.1d0)go to 80
      if(abs(cm1+cm2).gt.j3+0.1d0)go to 80
      if((cj1+cj2).lt.j3-0.1d0)go to 80
      if(abs(cj1-cj2).gt.j3+0.1d0)go to 80
      xnum=fact(j3+cj1-cj2,j3+cj2+cj1+1)*fact(j3-cj1+cj2,cj1-cm1)*
     1fact(cj1+cj2-j3,cj1+cm1)*fact(j3+m3,cj2-cm2)*
     2fact(j3-m3,cj2+cm2)
      xnum=sqrt(xnum)
      a=-cj1+cm1
      b=cj2-cj1-m3
      c=cj2+j3+cm1
      d=-cj1+cj2+j3
      e=j3-m3
      numin=0
      if(b.gt.0.d0)numin=b+0.1d0
      numax=c+0.1d0
      if(c.gt.d+0.1d0)numax=d+0.1d0
      f=numax
      if(f.gt.e+0.1d0)numax=e+0.1d0
      cj=cj1+cm2-m3
      kmi=numin+1
      kma=numax+1
      do 1 k=kmi,kma
      ck=k-1
      term=fact(1.d0,d-ck)*fact(c-ck,e-ck)*fact(ck-a,ck-b)*fact(1.d0,ck)
      ci=dmod(cj+ck,2.d0)
      if(dabs(ci).gt.0.5d0)term=-term
      ddr=ddr+term
   1  continue
      ddr=ddr*xnum
  80  dr=ddr
      return
      end
*************************************************************************
      function fact(cj1,cj2)
      implicit real*8(a-h,o-z)
      include 'mpif.h'
      n=abs(cj1-cj2)
      r=1.d0
      if(n.eq.0)go to 4
      k=1
      if(cj1.lt.cj2)k=-1
      cjx=dmax1(cj1,cj2)
      cjy=dmin1(cj1,cj2)
      hr=cjy
      do 1 i=1,n
      hr=hr+1.d0
      r=r*hr
  1   continue
      if(k.lt.0)r=1.d0/r
  4   fact=r
      return
      end

************************************************************************
      real*8 function slater (ia,ib,ic,id,k)

*   the value of this  function is the  slater integral  as normally   *
*   defined in terms of the four sets of quantum numbers  a,b,c,d.     *
*                                                                      *
*   subroutines called: quad, yzk.                                     *
*                                                                      *
************************************************************************
*
      implicit real*8 (a-h, o-z)
      include 'mpif.h'
*
      common/cons/zero,half,tenth,one,two,three,ten
     :/grid/r(MN),rp(MN),rpor(MN),rnt,h,hp,n
     :/tatb/ta(MN),tb(MN),mtp
     :/wave/pz(MNBAS),pf(MN,MNBAS),qf(MN,MNBAS),
     :mf(MNBAS)
*
      call yzk (k,ib,id)
*
*   multiply by second term, and obtain result by integration
*
      ta(1) = zero
      do 1 i = 2,mtp
         ta(i) =  (pf(i,ia)*pf(i,ic)+qf(i,ia)*qf(i,ic))
     :           *rpor(i)*tb(i)
    1 continue
*
      call quad (result)
      slater = result
*
*   debug printout
*
      return
*
* 300 format (/'  (',1i1,')',
*    :        /' r   (',1i2,1a2,',',1i2,1a2,';',1i2,1a2,',',1i2,1a2,') '
*    :        ,'= ',1pd21.14)
*
      end
************************************************************************
*                                                                      *
      subroutine quad (result)
*                                                                      *
*   the argument result is an approximation  to the integral of f(r)   *
*   from  zero  to  infinity,  where the values of rp(i)*f(r(i)) are   *
*   tabulated in the array  ta(i). the integral in the interval zero   *
*   to r(2) is computed by use of an analytical fit                    *
*                                                                      *
*                                sigma                                 *
*                      f(r) = a r                                      *
*                                                                      *
*   a five-point  closed  newton-cotes  formula (cf. f b hildebrand,   *
*   introduction to numerical analysis, second edition, mcgraw-hill,   *
*   new york, 1974, p 93)  is  used  to  compute the integral in the   *
*   interval  r(2:mtp).  the  contribution  from  the  tail  of  the   *
*   function beyond the last  tabular  point  (mtp) is assumed to be   *
*   negligible. the method uses  mtp+3  tabulation points. array  ta   *
*   should therefore be dimensioned to at least  n+4 .                 *
*                                                                      *
*   no subroutines called.                                             *
*                                                                      *
*   written by farid a parpia, at oxford   last updated: 10 nov 1989   *
*                                                                      *
************************************************************************
*
      implicit real*8 (a-h, o-z)
        include 'mpif.h'
*
      common/cons/zero,half,tenth,one,two,three,ten
     :      /def4/accy,nscf,nsic,nsolv
     :      /grid/r(MN),rp(MN),rpor(MN),rnt,h,hp,n
     :      /tatb/ta(MN),tb(MN),mtp
     :      /ncc/c1,c2,c3,c4
*
*   find first values that will permit computation of exponent
*
      mtpm1 = mtp-1
      do 3 i = 2,mtpm1
*
         tai = ta(i)
         if (abs (tai) .gt. zero) then
*
            ip1   = i+1
            taip1 = ta(ip1)
            quott = taip1/tai
*
            if (quott .gt. zero) then
*
*   exponent from fit
*
               frip1 = taip1/rp(ip1)
               fri   = tai  /rp(i  )
               ratio = frip1/fri
               rip1  = r (ip1)
               ri    = r (i  )
               sigma = log (ratio)/log (rip1/ri)
*
*   analytical integration and error estimate for interval r(1:i)
*
               fri    = ri*fri
               result = fri/(sigma+one)
*
*   set the tail to zero
*
               do 1 loc = 1,3
                  ta(mtp+loc) = zero
    1          continue
*
*   newton-cotes quadature for the remainder
*
               result = result+c1*tai
               do 2 loc = ip1,mtp,4
                  result = result+c2*(ta(loc  )+ta(loc+2))
     :                           +c3* ta(loc+1)
     :                           +c4* ta(loc+3)
    2          continue
               if (mod (mtp-i,4) .eq. 0) result = result-c1*ta(mtp)
*
*   test of result's accuracy
*
*              if (itc(18) .eq. 1) then
*                 estder = ten*ri*fri
*                 ratio = abs (estder/result)
*                 if (ratio .gt. accy) write (*,300) ratio
*              endif
*
               goto 4
*
            endif
*
         endif
*
    3 continue
*
*   no value which will permit computation of exponent
*
      result = zero
*
    4 return
*
* 300 format (/' ***** warning in subroutine quad *****'
*    :       //' estimated accuracy is ',1pd10.3,
*    :        /' decrease rnt or improve input data conditioning to'
*    :        ,' ameliorate.'/)
 
      end
************************************************************************
************************************************************************
*                                                                      *
      subroutine yzk (k,i,j)
*                                                                      *
*   ----------------   section 09   subprogram 32   ----------------   *
*                                                                      *
*   this subroutine evaluates hartree y- and z-functions:              *
*                                                                      *
*               (k)            (k)           (k)                       *
*              y   (i,j;r) =  z   (i,j;r) + w   (i,j;r)                *
*                                                                      *
*   where                                                              *
*                                                                      *
*    (k)                                                               *
*   z   (i,j;r) =  i ( (s/r)   (p (s)*p (s) + q (s)*q (s)) ; 0 - r )   *
*                                i     j       i     j                 *
*                                                                      *
*   where                                                              *
*                                                                      *
*    (k)                    k+1                                        *
*   w   (i,j;r) =  i ( (r/s)   (p (s)*p (s) + q (s)*q (s)) ; r -       *
*                                i     j       i     j    infinity )   *
*                                                                      *
*   where  i ( g(r,s) ; range )  denotes the integral of g(r,s) over   *
*   range  in  s .  the y-function is tabulated in  common/tatb/  in   *
*   array  tb , the z-function in array ta .                           *
*                                                                      *
*   subroutines called: draw                                           *
*                                                                      *
*   written by farid a parpia, at oxford   last updated: 13 nov 1989   *
*                                                                      *
************************************************************************
*
      implicit real*8 (a-h, o-z)
        include 'mpif.h'

      dimension rhop(MN),rttk(MN),wk(MN),temp(MN),
     :          yk(MN),zk(MN)
*
      common/cnc5/cnc5c(2:5,2:4)
     :      /cons/zero,half,tenth,one,two,three,ten
     :      /def4/accy,nscf,nsic,nsolv
     :      /grid/r(MN),rp(MN),rpor(MN),rnt,h,hp,n
     :      /ncc/c1,c2,c3,c4
     :      /tatb/ta(MN),tb(MN),mtp
     :      /wave/pz(MNBAS),pf(MN,MNBAS),qf(MN,MNBAS),
     :                                      mf(MNBAS)
*
      equivalence (ta(1),zk(1)),(tb(1),yk(1))
*
*                         k
*   for  k > 0  compute  r   and store in  rttk
*
      if (k .gt. 0) then
         do 1 ii = 2,n+4
            rttk(ii) = r(ii)**k
    1    continue
      endif
*
*   determine maximum tabulation point as location beyond which
*   rhop  (see comment statements below) would be zero; determine
*   other important locations
*
c     mtp = min (mf(i),mf(j))
      mtpp1 = mtp+1
      mtpp3 = mtp+3
      mtpp4 = mtp+4
*
*   compute rp(s)*(p (s)*p (s)+q (s)*q (s)) and store in rhop
*                   i     j     i     j
*
      do 2 ii = 2,mtp
         rhop(ii) = rp(ii)*(pf(ii,i)*pf(ii,j)+qf(ii,i)*qf(ii,j))
    2 continue
*
*   fill array temp with r**k * rhop
*
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
*
*   set an additional four points to zero
*
      do 5 ii = mtpp1,mtpp4
         temp(ii) = zero
    5 continue
*
*                                     k
*   compute the first few values of  r  * zk  using semi-open
*   newton-cotes formulae
*
      zk(1) = zero
      do 7 ii = 2,4
         sum = zero
         do 6 kk = 2,5
            sum = sum+cnc5c(kk,ii)*temp(kk)
    6    continue
         zk(ii) = sum
    7 continue
*                         k
*   compute remainder of r  * zk: march out to mtp+3; use closed
*   newton-cotes formula
*
      do 8 ii = 5,mtpp3
         zk(ii) = zk(ii-4)+c1*(temp(ii-4)+temp(ii  ))
     :                    +c2*(temp(ii-3)+temp(ii-1))
     :                    +c3* temp(ii-2)
    8 continue
*
*                                       k   (k)
*   determine the asymptotic value of  r * z
*
*                   (0)
*   correction to  z   : in the manner of  c froese fischer,
*   the hartree-fock method for atoms, john wiley & sons,
*   new york, 1977, p 235.
*
      if (k .eq. 0) then
*
         if (i .eq. j) then
            zklim = one
         else
            zklim = zero
         endif
*
         do 10 kk = mtpp3,mtp,-1
            dif = zk(kk)-zklim
            if (abs (dif) .gt. accy) then
               do 9 ii = kk,2,-4
                  zk(ii) = zk(ii)-dif
    9          continue
            endif
   10    continue
*
      else
*
         zklim = zk(mtpp3)
*
      endif
*
*   tabulate  zk  for entire internal grid
*
      if (k .eq. 0) then
*
         do 11 ii = mtpp4,n
            zk(ii) = zklim
   11    continue
*
      else
*
         do 12 ii = 2,mtpp3
            zk(ii) = zk(ii)/rttk(ii)
   12    continue
*
         do 13 ii = mtpp4,n
            zk(ii) = zklim/rttk(ii)
   13    continue
*
      endif
*
*   start array wk / r**(k+1)
*
      np4 = n+4
      do 14 ii = np4,mtpp1,-1
         wk(ii) = zero
   14 continue
*
*             k+1
*   compute  r       and store in rttk
*
      if (k .gt. 0) then
         do 15 ii = 2,n
            rttk(ii) = rttk(ii)*r(ii)
   15    continue
      endif
*
*   fill array temp with rhop / r**(k+1) ; set temp(1) = zero
*   to avoid 0/0 case
*
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
*
*   compute remainder of wk / r**(k+1): march in to the origin
*
      do 18 ii = mtp,2,-1
         wk(ii) = wk(ii+4)+c1*(temp(ii  )+temp(ii+4))
     :                    +c2*(temp(ii+1)+temp(ii+3))
     :                    +c3*(temp(ii+2))
   18 continue
      wk(1) = zero
*
*   compute wk
*
      if (k .eq. 0) then
         do 19 ii = 2,mtp
            wk(ii) = wk(ii)*r(ii)
   19    continue
      else
         do 20 ii = 2,mtp
            wk(ii) = wk(ii)*rttk(ii)
   20    continue
      endif
*
*   assemble solution
*
      yk(1) = zero
      do 21 ii = 2,n
         yk(ii) = zk(ii)+wk(ii)
   21 continue
*  the yk factors are multiplied by alpha when appropriate
*     if(iyk.ne.0)then
*     do 997 iii=1,n
*        yk(iii)=yk(iii)*alpha(k+1)
* 997 continue
*     endif
*
*   debug printout
*
*     if (itc(13) .ne. 0) then
*        write (ipd,300) k,np(i),nh(i),np(j),nh(j)
*        nb3 = n/3
*        if (3*nb3 .eq. n) then
*           nrows = nb3
*        else
*           nrows = nb3+1
*        endif
*        do 22 ii = 1,nrows
*           ii1 = ii
*           ii2 = ii1+nrows
*           ii3 = ii2+nrows
*           if (ii3 .le. n) then
*              write (ipd,301) r(ii1),yk(ii1),r(ii2),yk(ii2),
*    :                         r(ii3),yk(ii3)
*           elseif (ii2 .le. n) then
*              write (ipd,301) r(ii1),yk(ii1),r(ii2),yk(ii2)
*           else
*              write (ipd,301) r(ii1),yk(ii1)
*           endif
*  22    continue
*     endif
*
      return
*
* 300 format (///'  (',1i1,')'
*    :          /' y   (',1i2,1a2,',',1i2,1a2,';r) :'
*    :         //3(' --------- r --------- ------- y (r) -------'))
* 301 format (1p,6(1x,1d21.14))
*
      end

      subroutine ninej(rj11,rj12,rj13,rj21,rj22,rj23,rj31,rj32,
     : rj33,r9j)
      implicit real*8(a-h,o-z)
        include 'mpif.h'
      common/factor/ fct(0:50),mfd

      if(mfd.eq.33)goto 101
      mfd=33                          !calculate factorials if this
      fct(0)=1
      fct(1)=1                        !hasn't already been done
      do 20 i=2,mfd
      fct(i)=fct(i-1)*i
20    continue

101   tria=tri1(rj11,rj12,rj13)
      trib=tri1(rj21,rj22,rj23)
      tric=tri1(rj31,rj32,rj33)
      trid=tri1(rj11,rj21,rj31)
      trie=tri1(rj12,rj22,rj32)
      trif=tri1(rj13,rj23,rj33)

      if(tria.eq.-1.or.trib.eq.-1.or.tric.eq.-1.or.
     *  trid.eq.-1.or.trie.eq.-1.or.trif.eq.-1)then
c     type *,'inconsistent arguments to 9j symbol'
      endif

      if(tria.eq.0.or.trib.eq.0.or.tric.eq.0.or.
     *  trid.eq.0.or.trie.eq.0.or.trif.eq.0)then
      s9j=0
      return
      endif

!    comment apply formula in edmonds [(6.4.3) p.101]
      rkmin=rmax3(abs(2*(rj11-rj33)),abs(2*(rj32-rj21)),
     *          abs(2*(rj12-rj23)))/2
      rkmax=rmin3(2*(rj11+rj33),2*(rj32+rj21),2*(rj12+rj23))/2
c     if(rkmax.lt.rkmin)type *,'something wrong with k'
      sum=0
      rk=rkmin
99    if(rk.gt.rkmax)goto 22

      call sixj(rj11,rj21,rj31,rj32,rj33,rk,sja6)
      call sixj(rj12,rj22,rj32,rj21,rk,rj23,sjb6)
      call sixj(rj13,rj23,rj33,rk,rj11,rj12,sjc6)
      sum=sum+((-1)**(int(2*rk)))*(2*rk+1)*
     x sja6*sjb6*sjc6
      rk=rk+1
      goto 99
22    r9j=sum
      return
      end

      function tri1(rj1,rj2,rj3)
      implicit real*8 (a-h,o-z)
        include 'mpif.h'
      sumj=rj1+rj2+rj3
      isumj=int(sumj)
      if(sumj.ne.isumj)then           !if j1+j2+j3 is not an integer
      tri1=-1                       !then tri1=-1
      goto 999
      endif
      if((2*abs(rj1-rj2).le.2*rj3).and.
     *        (2*(rj1+rj2).ge.2*rj3))then
      tri1=1                        !triangle condition is satisfied
      goto 999                      !so tri=1
      endif

      tri1=0                          !j1+j2+j3 is an integer

999   return
      end
      function rmax3(x,y,z)
      implicit real*8 (a-h,o-z)
        include 'mpif.h'
      rmax3=x
      if(x.gt.y.and.x.gt.z)rmax3=x
      if(y.gt.x.and.y.gt.z)rmax3=y
      if(z.gt.x.and.z.gt.y)rmax3=z
      return
      end

      function rmin3(x,y,z)
      implicit real*8 (a-h,o-z)
        include 'mpif.h'
      rmin3=x
      if(x.lt.y.and.x.lt.z)rmin3=x
      if(y.lt.x.and.y.lt.z)rmin3=y
      if(z.lt.x.and.z.lt.y)rmin3=z
      return
      end

      subroutine symmcal
      implicit real*8 (a-h,o-z)
      include 'mpif.h'
      common/symmetry/mtbl(MNS,MNS)

C     THIS SUBROUTINE STORES THE PRODUCT PARITY
C     TWO ORBITALS OF SAME OR OPPOSITE PARITY
c     HERE 1 STANDS FOR EVEN AND 2 STANDS FOR ODD PARITY

c     even parity
      mtbl(1,1)=1
      mtbl(2,2)=1
      mtbl(3,3)=1
      mtbl(4,4)=1
      mtbl(5,5)=1
      mtbl(6,6)=1

      mtbl(1,3)=1
      mtbl(3,1)=1

      mtbl(3,5)=1
      mtbl(5,3)=1

      mtbl(1,5)=1
      mtbl(5,1)=1

      mtbl(2,6)=1
      mtbl(6,2)=1

      mtbl(6,4)=1
      mtbl(4,6)=1

      mtbl(2,4)=1
      mtbl(4,2)=1



c     odd parity

      mtbl(1,2)=2
      mtbl(2,1)=2

      mtbl(1,4)=2
      mtbl(4,1)=2

      mtbl(2,3)=2
      mtbl(3,2)=2

      mtbl(3,4)=2
      mtbl(4,3)=2

      mtbl(2,5)=2
      mtbl(5,2)=2

      mtbl(5,4)=2
      mtbl(4,5)=2

      mtbl(1,6)=2
      mtbl(6,1)=2

      mtbl(3,6)=2
      mtbl(6,3)=2

      mtbl(6,5)=2
      mtbl(5,6)=2

      return
      end

      subroutine setgrd
C
C     THIS SUBROUTINE SETS UP THE GRID USED 
C     IN SLATER, YZK AND QUAD
C
      implicit real*8 (a-h,o-z)
      include 'mpif.h'
      common/grid/r(MN),rp(MN),rpor(MN),rnt,h,hp,n
     :      /tatb/ta(MN),tb(MN),mtp
     :      /cons/zero,half,tenth,one,two,three,ten
      n = MN-10
      mtp=n
      h = 6.25d-02
      rnt = 2.0d-06
      r(1) = 0.
      rp(1) = rnt
      eph = exp (h)
      ett = 1.
*   set up the arrays r, rp, rpor
      do i = 2,MN
      ett = eph*ett
      ettm1 = ett-1.
      r(i) = rnt*ettm1
      rp(i) = rnt*ett
      rpor(i) = ett/ettm1
      enddo
c     costant
      zero=0.0d0
      half=0.50d0
      one=1.00d0
      two=2.00d0
      three=3.00d0
      ten=10.00d0
      tenth=one/ten
      return
      end
      
      subroutine daopen
C
C     OPEN THE FILES
C
      implicit real*8 (a-h,o-z)
      include 'mpif.h'

       call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)
      if(iproc .eq. MASTER)then

      open(STDIN,file='input',form='formatted',status='old')
      open(STDOUT,file='ccsd_all.out',form='formatted')
      open(WFNIN,file='wfn.dat',form='unformatted',status='old')
      endif
      return
      end

c****************************************************************
c   This subroutine  read the input file                        *
c****************************************************************

      subroutine readinp(nsym,nbas,iocc,iopt,nonlin,isolv,
     :nprint,itwoe)
C
C     I/O ROUTINE
      implicit real*8(a-h,o-z)
      include 'mpif.h'

      common/iallparity/iall
      dimension nbas(MNSYM),iocc(MNSYM)

      call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)

      if(iproc .eq. MASTER)then

c     read number of symmetry   
      read(STDIN,*)nsym,(nbas(i),i=1,nsym)
c     occupancy of each orbital
      read(STDIN,*)(iocc(i),i=1,nsym)
c     option for ccsd and ccd (if iopt=1 then ccsd else ccd)
      read(STDIN,*)iopt,nonlin,iall,itwoe,isolv,nprint
      endif
      call MPI_BCAST(nsym,1,MPI_INTEGER,MASTER,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(nbas,MNSYM,MPI_INTEGER,MASTER,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(iocc,MNSYM,MPI_INTEGER,MASTER,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(iopt,1,MPI_INTEGER,MASTER,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(iall,1,MPI_INTEGER,MASTER,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(itwoe,1,MPI_INTEGER,MASTER,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(isolv,1,MPI_INTEGER,MASTER,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(nprint,1,MPI_INTEGER,MASTER,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(nonlin,1,MPI_INTEGER,MASTER,MPI_COMM_WORLD,ierr)

      return
      end

      subroutine grasprd(nbasis,iproc)
C
C     THIS SUBROUTINE READS THE VARIABLE
C     PARAMETERS FROM TAPE "wfn.dat" CREATED IN
C     THE "GAUSSIAN BASIS CODE"
C
C     THE VAIRABLES R,RP, RPOR, MTP ETC. ARE
C     GIVEN ELSEWHERE (GRASP CODE)
C
C     THE VARIABLE "eorb" IS THE ORBITAL ENERGY
C
      implicit real*8 (a-h,o-z)
      include 'mpif.h'
      common/wave/pz(MNBAS),pf(MN,MNBAS),qf(MN,MNBAS),
     :       mf(MNBAS)
     :      /grid/r(MN),rp(MN),rpor(MN),rnt,h,hp,n
     :      /tatb/ta(MN),tb(MN),mtp
      common/orbital_energy/eorb(MNBAS)
c     local dimension 
      
      call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)

      if(iproc .eq. MASTER)then

      rewind (WFNIN)
      read(WFNIN)h,n
      mtp=n
      read(WFNIN)(r(i),i = 1,n),(rp(i),i = 1,n),(rpor(i),i=1,n)
      do i=1,nbasis
      read(WFNIN)eorb(i)
      mf(i)=n
      mfj=mf(i)
      read(WFNIN)(pf(j,i),j=1,mfj),(qf(j,i),j=1,mfj)
      enddo
      rewind WFNIN
       endif
      call MPI_BCAST(h,1,MPI_REAL8,MASTER,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(n,1,MPI_INTEGER,MASTER,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(r,MN,MPI_REAL8,MASTER,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(rp,MN,MPI_REAL8,MASTER,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(rpor,MN,MPI_REAL8,MASTER,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(pf,MNBAS*MN,MPI_REAL8,MASTER,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(qf,MNBAS*MN,MPI_REAL8,MASTER,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(eorb,MNBAS,MPI_REAL8,MASTER,MPI_COMM_WORLD,ierr)

      
      if(iproc.eq.0)then
         write(STDOUT,9999)
         write(STDOUT,9998)(eorb(i),i=1,nbasis)
      endif
      mtp=n
9999  format(//,10x,17h ORBITAL ENERGIES,/)
9998  format(4d20.11)
      return
      end
     

      subroutine setup(nsym,nbasis,nbas,iocc,iproc)
C
C     THIS SUBROUTINE SETS THE J, IA, PARITY ETC.
C     IT ALSO DISECTS THE ORBITALS INTO 'OCCUPIED
C     AND UNOCCUPIED' CATEGORY.
C
C     VARIABLES:
C                    NOCC==== NO. OF OCCUPIED ORBITALS
C                    NEXCIT== NO. OF EXCITED ORBTALS
C
C                    FOR OCCUPIED ORBITALS;
C                                           KC= ACTUAL INDEX
C                                           ORBC= J-VALES
C                                           ISYMC=PARITY
C                                           IQC= IA VALUE
C                    FOR UNOCCUPIED ORBITALS
C                                           KE= ACTUAL INDEX
C                                           ORBE=J-VALUE
C                                           ISYME=PARITY
C                                           IQE=IA VALUE
C
      implicit real*8 (a-h,o-z)
      include 'mpif.h'
      character*1 orbtyp,orbsym
      character*5 qj1,qj
      common/jvalue/orbc(MNOCC),orbe(MNEXC)
      common/index/ke(MNEXC),kc(MNOCC)
      common/info/nocc,nexcit,iType(MNBAS),iTCount(MNBAS)
      common/parity/isymc(MNOCC),isyme(MNEXC)
      common/symmetry/mtbl(MNS,MNS)
      common/kpavlaue/iqc(MNOCC),iqe(MNEXC)
      common/syminfo/orbj(MNBAS),iiq(MNBAS),iparity(MNBAS)
      common/maxj/jmax

c     local dimension
      dimension kap(MNSYM),ia(MNSYM),iorbsm(MNSYM)
      dimension orbtyp(MNS),nka(MNBAS)
      dimension orb(MNSYM),np(MNBAS)
      dimension qj(MNBAS),qj1(MNS)
      dimension nprimc(MNOCC),nkprimc(MNOCC)
      dimension nprime(MNEXC),nkprime(MNEXC)
      dimension nmc(MNBAS),nme(MNBAS)
      dimension nbas(MNSYM),iocc(MNSYM)


      data kap/-1, 1,-2, 2,-3, 3,-4,4,-5,5,-6/
      data orb/0.5,0.5,1.5,1.5,2.5,2.5,3.5,3.5,
     :         4.5,4.5,5.5/
      data ia /1,-1,1,-1,1,-1,1,-1,1,-1,1/
      data iorbsm/1,2,2,3,3,4,4,5,5,6,6/
      data qj1/' 1/2',' 3/2',' 5/2',' 7/2',' 9/2',
     :         ' 11/2'/
      data orbtyp/'s','p','d','f','g','h'/

      nocc=0
      nexcit=0
      nbasis=0
      do i=1,MNBAS
         isyme(i)=0
         isymc(i)=0
         iType(i)=0
      enddo

      jmax=2*orb(nsym)-1
      ii=0
      do isym=1,nsym              ! sum over symmetries s, p, d, ...
         ibasis=nbas(isym)
         if(ibasis.ne.0) then
            nbasis = nbasis + ibasis
            do jbas=1,ibasis      ! sum over all orbitals of given symm.
               ii = ii + 1
c-----------------------------------------------------------------------
c              ii is a running counter for the orbitals irrespective of
c              their symmetry.
c-----------------------------------------------------------------------
               nka(ii)=kap(isym)
               orbj(ii)=orb(isym)
               np(ii)=jbas
               iparity(ii)=iorbsm(isym)
               if(nka(ii).eq.kap(isym)) then
                  iiq(ii)=ia(isym)
               endif
               if(nka(ii).eq.kap(isym).and.
     &            jbas.le.iocc(isym))then
                  iType(ii) = 0
                  nocc=nocc+1
                  iTCount(ii) = nocc
                  kc(nocc)=ii
                  nkprimc(nocc)=kap(isym)
                  nprimc(nocc)=np(ii)
                  orbc(nocc)=orbj(ii)
                  nmc(nocc)=idint(2*(orbc(nocc))) +1
                  if(nkprimc(nocc).eq.kap(isym)) then
                     isymc(nocc)=iorbsm(isym)
                  endif
                  if(nkprimc(nocc).eq.kap(isym)) then
                     iqc(nocc)=ia(isym)
                  endif
               else
                  iType(ii) = 1
                  nexcit=nexcit+1
                  iTCount(ii) = nexcit
                  ke(nexcit)=ii
                  nprime(nexcit)=np(ii)
                  nkprime(nexcit)=nka(ii)
                  orbe(nexcit)=orb(isym)
                  nme(nexcit)=idint(2*(orbj(ii))) +1
                  if(nkprime(nexcit).eq.kap(isym)) then
                     isyme(nexcit)=iorbsm(isym)
                  endif
                  if(nkprime(nexcit).eq.kap(isym)) then
                     iqe(nexcit)=ia(isym)
                  endif
               endif
            enddo
         endif
      enddo

c     write occupied and unoccupied orbitals
      if(iproc.eq.0)then
         write(STDOUT,101)
         write(STDOUT,100)
      endif
      do i=1,nocc
         if(orbc(i).eq.0.5)qj(i)=qj1(1)
         if(orbc(i).eq.1.5)qj(i)=qj1(2)
         if(orbc(i).eq.2.5)qj(i)=qj1(3)
         if(orbc(i).eq.3.5)qj(i)=qj1(4)
         if(orbc(i).eq.4.5)qj(i)=qj1(5)
         if(orbc(i).eq.5.5)qj(i)=qj1(6)
         if(isymc(i).eq.1)orbsym=orbtyp(1)
         if(isymc(i).eq.2)orbsym=orbtyp(2)
         if(isymc(i).eq.3)orbsym=orbtyp(3)
         if(isymc(i).eq.4)orbsym=orbtyp(4)
         if(isymc(i).eq.5)orbsym=orbtyp(5)
         if(isymc(i).eq.6)orbsym=orbtyp(6)
         if(isymc(i).eq.3)isymc(i)=1
         if(isymc(i).eq.4)isymc(i)=2
         if(isymc(i).eq.5)isymc(i)=1
         if(isymc(i).eq.6)isymc(i)=2
         if(iproc.eq.0) then
            write(STDOUT,200)orbsym,nprimc(i),nmc(i),kc(i), qj(i)
         endif
      enddo

      if(iproc.eq.0)write(STDOUT,102)
      if(iproc.eq.0)write(STDOUT,100)
      do i=1,nexcit
         if(orbe(i).eq.0.5)qj(i)=qj1(1)
         if(orbe(i).eq.1.5)qj(i)=qj1(2)
         if(orbe(i).eq.2.5)qj(i)=qj1(3)
         if(orbe(i).eq.3.5)qj(i)=qj1(4)
         if(orbe(i).eq.4.5)qj(i)=qj1(5)
         if(orbe(i).eq.5.5)qj(i)=qj1(6)
         if(isyme(i).eq.1)orbsym=orbtyp(1)
         if(isyme(i).eq.2)orbsym=orbtyp(2)
         if(isyme(i).eq.3)orbsym=orbtyp(3)
         if(isyme(i).eq.4)orbsym=orbtyp(4)
         if(isyme(i).eq.5)orbsym=orbtyp(5)
         if(isyme(i).eq.6)orbsym=orbtyp(6)
         if(isyme(i).eq.3)isyme(i)=1
         if(isyme(i).eq.4)isyme(i)=2
         if(isyme(i).eq.5)isyme(i)=1
         if(isyme(i).eq.6)isyme(i)=2
         if(iproc.eq.0) then
            write(STDOUT,200)orbsym,nprime(i),nme(i),ke(i),qj(i)
         endif
      enddo
 100  format(/,1x,'TYPE',2x,' N ',2x,'2J+1',2x,'ORB. NO.',
     $3x,' J ',/)
 101  format(/,10x,' OCCUPIED ORBITALS',/)
 102  format(/,10x,' EXCITED ORBITALS',/)
 200  format(2x,a,2x,3(i4,2x),6x,a4)
      return
      end



************************************************************************
*                                                                      *
      subroutine setqic
*                                                                      *
*   this  subroutine sets up the coefficients for subroutines dpbdt,   *
*   quad, rinti, start, yzk, zkf.                                      *

*   no subroutines called.                                             *
*                                                                      *
*   written by farid a parpia, at oxford    last update: 14 nov 1989   *
*                                                                      *
************************************************************************
*
      implicit real*8 (a-h, o-z)
      include 'mpif.h'
      logical first
*
      dimension b13(13,13),cg(6),c5num(1:5,2:5),c6num(1:6,2:6)
*
      common/cnc5/cnc5c(2:5,2:4)
     :      /cons/zero,half,tenth,one,two,three,ten
     :      /grid/r(MN),rp(MN),rpor(MN),rnt,h,hp,n
     :      /cnc6/cnc6c(1:6,2:6)
     :      /lic13/a13(13,13)
     :      /ncc/c1,c2,c3,c4
     :      /sbc/c(6)
*
*----------------------------------------------------------------------*
*                                                                      *
*                                                                      *
*   THIRTEEN-POINT  LAGRANGE  INTERPOLATION  COEFFICIENTS FOR FIRST    *
*   DERIVATIVE                                                         *
*
      DATA (B13( 1,I),I = 1,13) /            -1486442880.0D 00,
     :            5748019200.0D 00, -15807052800.0D 00,
     :           35126784000.0D 00, -59276448000.0D 00,
     :           75873853440.0D 00, -73766246400.0D 00,
     :           54195609600.0D 00, -29638224000.0D 00,
     :           11708928000.0D 00,  -3161410560.0D 00,
     :             522547200.0D 00,    -39916800.0D 00/
      DATA (B13( 2,I),I = 1,13) /              -39916800.0D 00,
     :            -967524480.0D 00,   2634508800.0D 00,
     :           -4390848000.0D 00,   6586272000.0D 00,
     :           -7903526400.0D 00,   7376624640.0D 00,
     :           -5269017600.0D 00,   2822688000.0D 00,
     :           -1097712000.0D 00,    292723200.0D 00,
     :             -47900160.0D 00,      3628800.0D 00/
      DATA (B13( 3,I),I = 1,13) /                3628800.0D 00,
     :             -87091200.0D 00,   -684478080.0D 00,
     :            1596672000.0D 00,  -1796256000.0D 00,
     :            1916006400.0D 00,  -1676505600.0D 00,
     :            1149603840.0D 00,   -598752000.0D 00,
     :             228096000.0D 00,    -59875200.0D 00,
     :               9676800.0D 00,      -725760.0D 00/
      DATA (B13( 4,I),I = 1,13) /                -725760.0D 00,
     :              13063680.0D 00,   -143700480.0D 00,
     :            -476910720.0D 00,   1077753600.0D 00,
     :            -862202880.0D 00,    670602240.0D 00,
     :            -431101440.0D 00,    215550720.0D 00,
     :             -79833600.0D 00,     20528640.0D 00,
     :              -3265920.0D 00,       241920.0D 00/
      DATA (B13( 5,I),I = 1,13) /                 241920.0D 00,
     :              -3870720.0D 00,     31933440.0D 00,
     :            -212889600.0D 00,   -303937920.0D 00,
     :             766402560.0D 00,   -447068160.0D 00,
     :             255467520.0D 00,   -119750400.0D 00,
     :              42577920.0D 00,    -10644480.0D 00,
     :               1658880.0D 00,      -120960.0D 00/
      DATA (B13( 6,I),I = 1,13) /                -120960.0D 00,
     :               1814400.0D 00,    -13305600.0D 00,
     :              66528000.0D 00,   -299376000.0D 00,
     :            -148262400.0D 00,    558835200.0D 00,
     :            -239500800.0D 00,     99792000.0D 00,
     :             -33264000.0D 00,      7983360.0D 00,
     :              -1209600.0D 00,        86400.0D 00/
      DATA (B13( 7,I),I = 1,13) /                  86400.0D 00,
     :              -1244160.0D 00,      8553600.0D 00,
     :             -38016000.0D 00,    128304000.0D 00,
     :            -410572800.0D 00,            0.0D 00,
     :             410572800.0D 00,   -128304000.0D 00,
     :              38016000.0D 00,     -8553600.0D 00,
     :               1244160.0D 00,       -86400.0D 00/
      DATA (B13( 8,I),I = 1,13) /                 -86400.0D 00,
     :               1209600.0D 00,     -7983360.0D 00,
     :              33264000.0D 00,    -99792000.0D 00,
     :             239500800.0D 00,   -558835200.0D 00,
     :             148262400.0D 00,    299376000.0D 00,
     :             -66528000.0D 00,     13305600.0D 00,
     :              -1814400.0D 00,       120960.0D 00/
      DATA (B13( 9,I),I = 1,13) /                 120960.0D 00,
     :              -1658880.0D 00,     10644480.0D 00,
     :             -42577920.0D 00,    119750400.0D 00,
     :            -255467520.0D 00,    447068160.0D 00,
     :            -766402560.0D 00,    303937920.0D 00,
     :             212889600.0D 00,    -31933440.0D 00,
     :               3870720.0D 00,      -241920.0D 00/
      DATA (B13(10,I),I = 1,13) /                -241920.0D 00,
     :               3265920.0D 00,    -20528640.0D 00,
     :              79833600.0D 00,   -215550720.0D 00,
     :             431101440.0D 00,   -670602240.0D 00,
     :             862202880.0D 00,  -1077753600.0D 00,
     :             476910720.0D 00,    143700480.0D 00,
     :             -13063680.0D 00,       725760.0D 00/
      DATA (B13(11,I),I = 1,13) /                 725760.0D 00,
     :              -9676800.0D 00,     59875200.0D 00,
     :            -228096000.0D 00,    598752000.0D 00,
     :           -1149603840.0D 00,   1676505600.0D 00,
     :           -1916006400.0D 00,   1796256000.0D 00,
     :           -1596672000.0D 00,    684478080.0D 00,
     :              87091200.0D 00,     -3628800.0D 00/
      DATA (B13(12,I),I = 1,13) /               -3628800.0D 00,
     :              47900160.0D 00,   -292723200.0D 00,
     :            1097712000.0D 00,  -2822688000.0D 00,
     :            5269017600.0D 00,  -7376624640.0D 00,
     :            7903526400.0D 00,  -6586272000.0D 00,
     :            4390848000.0D 00,  -2634508800.0D 00,
     :             967524480.0D 00,     39916800.0D 00/
      DATA (B13(13,I),I = 1,13) /               39916800.0D 00,
     :            -522547200.0D 00,   3161410560.0D 00,
     :          -11708928000.0D 00,  29638224000.0D 00,
     :          -54195609600.0D 00,  73766246400.0D 00,
     :          -75873853440.0D 00,  59276448000.0D 00,
     :          -35126784000.0D 00,  15807052800.0D 00,
     :           -5748019200.0D 00,   1486442880.0D 00/
*
      DATA B13DEN/479001600.0D 00/
*
*----------------------------------------------------------------------*
*   COEFFICIENTS FOR SIENKIEWICZ-BAYLIS ALGORITHM                      *
*
      DATA CG / 1771.0D 00,
     :          9235.0D 00,
     :          5890.0D 00,
     :          4610.0D 00,
     :            35.0D 00,
     :            59.0D 00/
*
      DATA DENOM /5760.0D 00/
*
*----------------------------------------------------------------------*
*                                                                      *
*   FIVE-POINT NEWTON-COTES COEFFICIENTS FOR CLOSED INTEGRATION. EX-   *
*   PRESSED AS RATIONAL NUMBERS                                        *
*
      DATA (C5NUM(I,2),I = 1,5)/ 251.0D 00, 646.0D 00,
     :       -264.0D 00, 106.0D 00, -19.0D 00/
      DATA (C5NUM(I,3),I = 1,5)/ 232.0D 00, 992.0D 00,
     :        192.0D 00,  32.0D 00,  -8.0D 00/
      DATA (C5NUM(I,4),I = 1,5)/ 243.0D 00, 918.0D 00,
     :        648.0D 00, 378.0D 00, -27.0D 00/
      DATA (C5NUM(I,5),I = 1,5)/ 224.0D 00,1024.0D 00,
     :        384.0D 00,1024.0D 00, 224.0D 00/
*
      DATA C5DEN/ 720.0D 00/
*
*----------------------------------------------------------------------*
*                                                                      *
*   SIX-POINT NEWTON-COTES COEFFICIENTS FOR CLOSED INTEGRATION.  EX-   *
*   PRESSED AS RATIONAL NUMBERS                                        *
*
      DATA (C6NUM(I,2),I = 1,6)/ 475.0D 00,1427.0D 00,
     :                          -798.0D 00, 482.0D 00,
     :                          -173.0D 00,  27.0D 00/
      DATA (C6NUM(I,3),I = 1,6)/ 448.0D 00,2064.0D 00,
     :                           224.0D 00, 224.0D 00,
     :                           -96.0D 00,  16.0D 00/
      DATA (C6NUM(I,4),I = 1,6)/ 459.0D 00,1971.0D 00,
     :                          1026.0D 00,1026.0D 00,
     :                          -189.0D 00,  27.0D 00/
      DATA (C6NUM(I,5),I = 1,6)/ 448.0D 00,2048.0D 00,
     :                           768.0D 00,2048.0D 00,
     :                           448.0D 00,   0.0D 00/
      DATA (C6NUM(I,6),I = 1,6)/ 475.0D 00,1875.0D 00,
     :                          1250.0D 00,1250.0D 00,
     :                          1875.0D 00, 475.0D 00/
*
      DATA C6DEN/1440.0D 00/
*
*----------------------------------------------------------------------*
*
      DATA FIRST /.TRUE./
*
      DATA H/6.250D-2/

*   LAGRANGE INTERPOLATION COEFFICIENTS
*
*   DO THIS INITIALIZATION ONCE PER RUN ONLY
*
      IF (FIRST) THEN
*
*   THIRTEEN-POINT COEFFICIENTS FOR DPBDT
*
         FACTOR = ONE/B13DEN
         DO 2 J = 1,13
            DO 1 I = 1,13
               A13(I,J) = B13(I,J)*FACTOR
    1       CONTINUE
    2    CONTINUE
*
         FIRST = .FALSE.
*
      ENDIF
*
*   SIENKIEWICZ-BAYLIS COEFFICIENTS FOR SBSTEP
*
      C(1) = CG(1)/DENOM
      FACTOR = H/DENOM
      DO 3 I = 2,6
         C(I) = CG(I)*FACTOR
    3 CONTINUE
*
*   NEWTON-COTES COEFFICIENTS FOR YZK AND QUAD
*
      FACTOR = H/C5DEN
      DO 5 J = 2,4
         DO 4 I = 2,5
            CNC5C(I,J) = FACTOR*C5NUM(I,J)
    4    CONTINUE
    5 CONTINUE
*
      C1 = FACTOR*C5NUM(1,5)
      C2 = FACTOR*C5NUM(2,5)
      C3 = FACTOR*C5NUM(3,5)
      C4 = C1+C1
*
*   NEWTON-COTES COEFFICIENTS FOR START
*
      FACTOR = H/C6DEN
      DO 7 J = 2,6
         DO 6 I = 1,6
            CNC6C(I,J) = FACTOR*C6NUM(I,J)
    6    CONTINUE
    7 CONTINUE
*
      RETURN
      END

      subroutine daclose
      implicit real*8(a-h,o-z)
      include 'mpif.h'
       call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)
      if(iproc .eq. MASTER)then

      close (STDIN)
      close (STDOUT)
      close (WFNIN)
      endif
      return
      end


      subroutine amat(iopt,ntmax,isolv,iproc)

c     THIS SUBROUTINE CONSTRUCTS THE A MATRIX OF
C     CC EQUATION OF THE FORM B(T)*T + A =0
C     FOR DETAILS CHECK THE DOC. FILE
C
      implicit real*8 (a-h,o-z)
      include 'mpif.h'
c     common block
      common/jvalue/orbc(MNOCC),orbe(MNEXC)
      common/scratch/b(MDIM),ad(MDIM),baux(MDIM)
      common/index/ke(MNEXC),kc(MNOCC)
      common/info/nocc,nexcit,iType(MNBAS),iTCount(MNBAS)
      common/parity/isymc(MNOCC),isyme(MNEXC)
      common/symmetry/mtbl(MNS,MNS)
      common/kpavlaue/iqc(MNOCC),iqe(MNEXC)
      common/orbital_energy/eorb(MNBAS)
      common/idra1/idra(MNBAS,MNBAS,0:MXV)
      common/idpa1/idpa(MNBAS,MNBAS)

c     local dimension
      dimension kgot(MXV+1)
      
c     data
      data two,half,onep5,one/2.0,0.5,1.5,1.0d0/
      data zero/0.0d0/

      nbasis=nocc+nexcit

      if(iopt.ne.0)then

      do 1 ia=1,nocc
      iasym=isymc(ia)
      do 2 ip=1,nexcit
      ipsym=isyme(ip)
      if(iasym.ne.ipsym)go to 2
      if(orbc(ia).ne.orbe(ip))go to 2
      iia=kc(ia)
      iip=ke(ip)
      neqn=idpa(iip,iia)
      b(neqn)=zero
      ad(neqn)=eorb(iip)-eorb(iia)
  2   continue
  1   continue

      endif

c     construction of a-matrix

      do 50 ia=1,nocc
      iasym=isymc(ia)
      do 60 ip=1,nexcit
      ipsym=isyme(ip)
      iapsym=mtbl(ipsym,iasym)
      do 70 ib=1,nocc
      ibsym=isymc(ib)
      do 80 iq=1,nexcit
      iqsym=isyme(iq)
      ibqsym=mtbl(iqsym,ibsym)
      if(ibqsym.ne.iapsym)go to 80
      orba=orbc(ia)
      orbb=orbc(ib)
      orbp=orbe(ip)
      orbq=orbe(iq)
      iaa=iqc(ia)
      iab=iqc(ib)
      iap=iqe(ip)
      iaq=iqe(iq)
      iia=kc(ia)
      iib=kc(ib)
      iip=ke(ip)
      iiq=ke(iq)
      call findk(orba,orbp,orbb,orbq,iaa,iap,iab,iaq,kgot,
     $kmax)
      if(kmax.eq.0)go to 80
      do 90 nloop=1,kmax
      kk=kgot(nloop)
      ak=dfloat(kk)
      nap=idra(iip,iia,kk)
      nbq=idra(iiq,iib,kk)
      call iloc(iip,iia,iiq,iib,ipa,iqb)
      if(ipa.lt.iqb)go to 90
      neqn=ntloc(nap,nbq,iapsym,kk)
      ipqdis=idis(iip,iiq)
      iabloc=indxr(iia)+iib
      if(iiq.gt.iip)iabloc=indxr(iib)+iia
c      rk=slater(iip,iiq,iia,iib,kk)
      f1=dr(orbp,ak,orba,half,zero,-half)
      f2=dr(orbq,ak,orbb,half,zero,-half)
      f3=dsqrt((two*orba+1)*(two*orbb+1)*(two*orbp+1)*(two*orbq+1))
      f4=(-one)**(kk+orbp+orbq+1)
      rk = vint(iip,iiq,iia,iib,kk)
c     b(neqn)=-rk*f1*f2*f3*f4
      b(neqn)=-rk*f4
      ad(neqn)=eorb(iip)-eorb(iia)+eorb(iiq)-eorb(iib)
c     if ((ip.ge.7.and.ip.le.9).and.(iq.ge.7.and.iq.le.9)) then
c     print*,'Radial integral = ',iip,iiq,iia,iib,rk,f1,f2,f3,f4
c     print*,'energy denominator = ',ip,iq,ia,ib,neqn,ad(neqn)
c     print*,'energies = ',eorb(iip),eorb(iia),eorb(iiq),eorb(iib)
c     endif
      baux(neqn)=rk/(f1*f2*f3)
  90  continue
  80  continue
  70  continue
  60  continue
  50  continue

      if(isolv.ne.0) then
         write(STDOUT,*)'please use Jacobi method, isolv = 0'
         isolv = 0
      endif

9999  format(/,2x,25h COUPLED CLUSTER A MATRIX,/)
9998  format(/,2x,26h COUPLED CLUSTER BD MATRIX,/)
9997  format(4f20.11)
      return
      end




c******************************************************************
c                                                                 *
c                 SOUBROUTINE SINGLES                             *
c                                                                 *
c     THIS SUBROUTINE SETS UP THE B- MATRIX FOR LINEARIZED CC     *
c     EQUATIONS T1-T1 AND T1-T2 BLOCK                             *
c                                                                 *
c     FOR DETAILS CHECK THE DOC. FILE                             *
c******************************************************************

      subroutine singles(raux,t,ntmax,isolv,rowstart,rowend)


      implicit real*8 (a-h,o-z)
      include 'mpif.h'
      integer row,rowstart,rowend
      real*8 raux_local,baux_local
c     common block
      common/jvalue/orbc(MNOCC),orbe(MNEXC)
      common/scratch/b(MDIM),ad(MDIM),baux(MDIM)
      common/index/ke(MNEXC),kc(MNOCC)
      common/info/nocc,nexcit,iType(MNBAS),iTCount(MNBAS)
      common/parity/isymc(MNOCC),isyme(MNEXC)
      common/symmetry/mtbl(MNS,MNS)
      common/orbital_energy/eorb(MNBAS)
      common/kpavlaue/iqc(MNOCC),iqe(MNEXC)
      common/idra1/idra(MNBAS,MNBAS,0:MXV)
      common/idpa1/idpa(MNBAS,MNBAS)
      common/skip/nnskip(0:MXV,2),n2skip,nsing
c     local dimension
      dimension t(MDIM),raux(MDIM),baux_local(MDIM),
     :raux_local(MDIM)

      dimension kgot(MXV+1),jgot(MXV+1),igot(MXV+1)

      data zero,half,one,two/0.0d0,0.50d0,1.0d0,2.0d0/
      call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)

      do i= 1,MDIM
       raux_local(i)=0.0
      enddo
      do i= 1,MDIM
       baux_local(i)=0.0
      enddo


c-------------------------------------------------------------------------
c     start with core orbital 'A'  and  virtual orbital 'P'
c-------------------------------------------------------------------------

      do 10 ia=1,nocc
      iasym=isymc(ia)
      orba=orbc(ia)
      iia=kc(ia)
      iaa=iqc(ia)

      do 20 ip=iproc+1,nexcit,nproc
      ipsym=isyme(ip)
      orbp=orbe(ip)
      iip=ke(ip)
      iap=iqe(ip)

      if(iasym.ne.ipsym)go to 20
      if(orba.ne.orbp)go to 20
                                                                                
      call init(baux_local)
      neqn=idpa(iip,iia)
                                                                                
      if(isolv.ne.0)then
      baux_local(neqn)=baux_local(neqn)+eorb(iip)-eorb(iia)
      endif
                                                                                
      row=idpa(iip,iia)
      irow=row+1-rowstart
      if(irow.le.0)go to 20
      if(irow.gt.rowend-rowstart+1)go to 20

c---------------------------------------------------------------------------
c      start with core orbital 'B' and virtual orbital 'Q'
c---------------------------------------------------------------------------
                                                                                
      do 30 ib=1,nocc
      ibsym=isymc(ib)
      orbb=orbc(ib)
      iib=kc(ib)
      iab=iqc(ib)

      do 40 iq=1,nexcit
      iqsym=isyme(iq)
      iiq=ke(iq)
      iaq=iqe(iq)
      orbq=orbe(iq)
                                                                                
      if(iqsym.ne.ibsym)go to 50
      if(orbb.ne.orbq)go to 50
                                                                                
      in=idpa(iiq,iib)

c-------------------------------------------------------------------------
c       diagram  <pb|v|aq><q|t1|b>
c-------------------------------------------------------------------------  
                                                                                                                                                               
      rk=vint(iip,iib,iia,iiq,0)
      f4=(-1)**(orbb+orbp+1+0)
      twoe=rk*f4
      fact=dsqrt((two*orbb+1)/(two*orba+1))
      baux_local(in)=baux_local(in)+twoe*fact

c-------------------------------------------------------------------------
c      diagram <pb|v|qa><q|t1|b> 
c-------------------------------------------------------------------------                                                                               

      call findk(orba,orbb,orbp,orbq,iaa,iab,iap,iaq,kgot,
     $kmax)
      if(kmax.eq.0) go to 50
      do  60 nloop=1,kmax
      ak=dfloat(kgot(nloop))
      kk=kgot(nloop)
      rk=vint(iip,iib,iiq,iia,kk)
      f4=(-1)**(kk+orbb+orbp+1)
      twoe=rk*f4
      fact=one/(two*orba+one)
      isign=(-one)**(orbb+orba+kk)
      baux_local(in)=baux_local(in)+twoe*fact*isign

  60  continue
  50  continue

c---------------------------------------------------------------------------
c     start with virtual orbital 'R'
c---------------------------------------------------------------------------
                                                                                

      do 70 ir=1,nexcit
      iir=ke(ir)
      iar=iqe(ir)
      irsym=isyme(ir)
      orbr=orbe(ir)
                                                                                
      ibrsym=mtbl(ibsym,irsym)
      ibqsym=mtbl(ibsym,iqsym)
      iaqsym=mtbl(iasym,iqsym)
      ipqsym=mtbl(ipsym,iqsym)

      if(ipqsym.ne.ibrsym)go to 70

c-------------------------------------------------------------------------
c               diagram <pb|v|qr><qr|t2|ab>
c--------------------------------------------------------------------------

      call findk(orbp,orbq,orbb,orbr,iap,iaq,iab,iar,igot,
     $imax)
      if(imax.eq.0)go to 70
      do 80 mloop=1,imax
      kk=igot(mloop)
      ak=dfloat(kk)
      ipbdis=idis(iip,iib)
      iqrloc=indxr(iiq)+iir
      if(iib.gt.iip)iqrloc=indxr(iir)+iiq
      rk=vint(iip,iib,iiq,iir,kk)
      f4=(-1)**(orbp+orbb+kk+1)
      twoe=rk*f4
      fact=one/((two*ak+one)*(two*orba+one))
      isign=(-one)**(orbp+orbq+orbb+orbr)
      naq=idra(iiq,iia,kk)
      nbr=idra(iir,iib,kk)
      in=ntloc(naq,nbr,iaqsym,kk)
      baux_local(in)=baux_local(in)+twoe*fact*isign

c----------------------------------------------------------------------------
c     diagram  <pb|v|qr><qr|t2|ba>
c----------------------------------------------------------------------------
                                                                                
      call findk(orbb,orbq,orbr,orba,iab,iaq,iar,iaa,jgot,
     $jmax)
      if(jmax.eq.0)go to 80
      do 90 kloop=1,jmax
      ll=jgot(kloop)
      al=dfloat(ll)
      nbq=idra(iiq,iib,ll)
      nar=idra(iir,iia,ll)
      in=ntloc(nbq,nar,ibqsym,ll)
      call sixj(ak,orbr,orbb,al,orbq,orba,s6j)
      baux_local(in)=baux_local(in)+twoe*s6j*isign/(two*orba+1)

   90 continue
   80 continue
   70 continue


c----------------------------------------------------------------------
c      start with core orbital 'G'
c-----------------------------------------------------------------------
                                                                                

      do 100 ig=1,nocc
      iig=kc(ig)
      iag=iqc(ig)
      igsym=isymc(ig)
      orbg=orbc(ig)
                                                                                
      iabsym=mtbl(iasym,ibsym)
      iqgsym=mtbl(iqsym,igsym)
      ipgsym=mtbl(ipsym,igsym)

      if(iabsym.ne.iqgsym)go to 100

c--------------------------------------------------------------------------
c          diagram <bg|v|aq><pq|t2|bg> 
c--------------------------------------------------------------------------

      call findk(orba,orbb,orbq,orbg,iaa,iab,iaq,iag,igot,
     $imax)
      if(imax.eq.0)go to 110
      do 120 mloop=1,imax
      kk=igot(mloop)
      ak=dfloat(kk)
      nbp=idra(iip,iib,kk)
      ngq=idra(iiq,iig,kk)
      in=ntloc(nbp,ngq,iqgsym,kk)
      rk=vint(iib,iig,iia,iiq,kk)
      f4=(-one)**(orbg+orbb+one+ak)
      twoe=rk*f4
      isign=(-one)**(orbg+orbb+orbq+orba)
      fact=one/((two*ak+one)*(two*orba+one))
      baux_local(in)=baux_local(in)-twoe*fact*isign

c----------------------------------------------------------------------------
c                diagram <bg|v|aq><pq|t2|gb>
c----------------------------------------------------------------------------


      call findk(orbp,orbg,orbq,orbb,iap,iag,iaq,iab,jgot,
     $jmax)
      if(jmax.eq.0)go to 120
      do 130 kloop=1,jmax
      ll=jgot(kloop)
      al=dfloat(ll)
      ngp=idra(iip,iig,ll)
      nbq=idra(iiq,iib,ll)
      in=ntloc(ngp,nbq,ipgsym,ll)
      call sixj(ak,orbp,orbb,al,orbq,orbg,s6j)
      baux_local(in)=baux_local(in)-twoe*s6j*isign/(two*orba+one)

  130 continue
  120 continue
  110 continue
  100 continue
   40 continue
   30 continue
                                                                                
      xsum=0
      do in=1,ntmax
      xsum=xsum+baux_local(in)*t(in)
      enddo
      raux_local(neqn)=xsum
     
  20  continue
  10  continue
             call MPI_ALLREDUCE(raux_local,raux,MDIM,MPI_REAL8,MPI_SUM,
     :                   MPI_COMM_WORLD,ierr)

      do i= 1,MDIM
       raux_local(i)=0.0
      enddo

             call MPI_ALLREDUCE(baux_local,baux,MDIM,MPI_REAL8,MPI_SUM,
     :                   MPI_COMM_WORLD,ierr)

      do i= 1,MDIM
       baux_local(i)=0.0
      enddo


c     print*,'single = raux,baux',raux


      return
      end

c****************************************************************************
c                                                                           *
c       subroutine doubles(raux,t,ntmax,iopt,isolv,rowstart,rowend)         *
c                                                                           *   
c     THIS SUBROUTINE SETS UP THE B- MATRIX FOR LINEARIZED CC               *
c     EQUATIONS T2-T2 AND T2-T1 BLOCK                                       *
c                                                                           *
c     FOR DETAILS CHECK THE DOC. FILE                                       *
c                                                                           *
c****************************************************************************

      subroutine doubles(raux,t,ntmax,iopt,isolv,rowstart,rowend)
      implicit real*8 (a-h,o-z)

      include 'mpif.h'
      integer row,rowstart,rowend
c     real*8 raux_local,baux_local
c     common block
      common/jvalue/orbc(MNOCC),orbe(MNEXC)
      common/index/ke(MNEXC),kc(MNOCC)
      common/scratch/b(MDIM),ad(MDIM),baux(MDIM)
      common/info/nocc,nexcit,iType(MNBAS),iTCount(MNBAS)
      common/orbital_energy/eorb(MNBAS)
      common/parity/isymc(MNOCC),isyme(MNEXC)
      common/symmetry/mtbl(MNS,MNS)
      common/kpavlaue/iqc(MNOCC),iqe(MNEXC)
      common/idra1/idra(MNBAS,MNBAS,0:MXV)
      common/idpa1/idpa(MNBAS,MNBAS)
      common/skip/nnskip(0:MXV,2),n2skip,nsing
c     local dimension
      dimension t(MDIM),raux(MDIM)
     : ,raux_local(MDIM),baux_local(MDIM)
      dimension kgot(MXV+1),jgot(MXV+1),igot(MXV+1)

      data zero,half,one,onep5,two/0.0,0.5,1.0,1.5,2.0/
      data four/4.0d0/

c     call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
c     call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)

      do i=0,MXV
      kgot(i)=0
      jgot(i)=0
      igot(i)=0
      enddo

c     do i= 1,MDIM
c      raux_local(i)=0.0
c     enddo
c     do i= 1,MDIM
c      baux_local(i)=0.0
c     enddo



      do 10 ia=1,nocc
      iasym=isymc(ia)

c     do 20 ip=iproc+1,nexcit,nproc
      do 20 ip=1,nexcit
      ipsym=isyme(ip)

      iapsym=mtbl(ipsym,iasym)

      do 30 ib=1,nocc
      ibsym=isymc(ib)

      do 40 iq=1,nexcit
      iqsym=isyme(iq)

      ibqsym=mtbl(ibsym,iqsym)

      if(ibqsym.ne.iapsym)go to 40

      orba=orbc(ia)
      orbb=orbc(ib)
      orbp=orbe(ip)
      orbq=orbe(iq)

      iaa=iqc(ia)
      iab=iqc(ib)
      iap=iqe(ip)
      iaq=iqe(iq)

      iia=kc(ia)
      iib=kc(ib)
      iip=ke(ip)
      iiq=ke(iq)

      call findk(orba,orbp,orbb,orbq,iaa,iap,iab,iaq,kgot,
     $kmax)
      if(kmax.eq.0)go to 40
      do 45 jloop=1,kmax
      l1=kgot(jloop)
      al1=dfloat(l1)
      call iloc(iip,iia,iiq,iib,ipa,iqb)
      if(ipa.lt.iqb)go to 45
      nap=idra(iip,iia,l1)
      nbq=idra(iiq,iib,l1)
      neqn=ntloc(nap,nbq,iapsym,l1)
      row=neqn
      irow=neqn+1-rowstart
      if(irow.le.0)go to 45
      if(irow.gt.rowend-rowstart+1)go to 45
c     call init(baux_local)
      call init(baux)
      mult=2
      mloop=1
      if(iip.eq.iiq.and.iia.eq.iib)go to 50
      mult=1
      mloop=2

  50  do 60 nloop=1,mloop
      go to (70,80),nloop

  70  jp=iip
      ja=iia
      jq=iiq
      jb=iib

      jasym=iasym
      jbsym=ibsym
      jpsym=ipsym
      jqsym=iqsym

      jaa=iaa
      jab=iab
      jap=iap
      jaq=iaq

      orbap=orba
      orbbp=orbb
      orbpp=orbp
      orbqp=orbq

      go to 90

   80 jp=iiq
      ja=iib
      jq=iip
      jb=iia

      jasym=ibsym
      jbsym=iasym
      jpsym=iqsym
      jqsym=ipsym

      jaa=iab
      jab=iaa
      jap=iaq
      jaq=iap

      orbap=orbb
      orbbp=orba
      orbpp=orbq
      orbqp=orbp

   90 continue
                                                                                
      if(isolv.ne.0)then
c     baux_local(neqn)=baux_local(neqn)+(eorb(jp)-eorb(ja))*mult
      baux(neqn)=baux(neqn)+(eorb(jp)-eorb(ja))*mult
      endif
                                                                                

      do 800 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      iir=ke(ir)

      do 900 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      iig=kc(ig)

      irgsym=mtbl(irsym,igsym)
      iqgsym=mtbl(jqsym,igsym)
      jqbsym=mtbl(jqsym,jbsym)

      if(irgsym.ne.jqbsym)go to 910

c---------------------------------------------------------------------------
c     diagram  <pg|v|ar><qr|t(2)|bg>
c--------------------------------------------------------------------------

      f4=(-1)**(l1+orbg+orbpp+1)
      rk=vint(jp,iig,ja,iir,l1)
      twoe=rk*f4
      ngr=idra(iir,iig,l1)
      nbq=idra(jq,jb,l1)
      in=ntloc(ngr,nbq,jqbsym,l1)
      fact=(-1)**(orbr+orbg+al1)
c     baux_local(in)=baux_local(in)-twoe*fact*mult/(two*l1+1)
      baux(in)=baux(in)-twoe*fact*mult/(two*l1+1)

c---------------------------------------------------------------------------
c     diagram  <pg|v|ar><qr|t(2)|gb>
c---------------------------------------------------------------------------

      call findk(orbg,orbqp,orbbp,orbr,iag,jaq,jab,iar,igot,imax)
      if(imax.eq.0)go to 920
      do 930 mloop=1,imax
      l3=igot(mloop)
      al3=dfloat(l3)
      nbr=idra(iir,jb,l3)
      ngq=idra(jq,iig,l3)
      in=ntloc(nbr,ngq,iqgsym,l3)
      call sixj(al1,orbr,orbg,al3,orbqp,orbbp,s6j)
c     baux_local(in)=baux_local(in)-twoe*fact*s6j*mult
      baux(in)=baux(in)-twoe*fact*s6j*mult

  930 continue
  920 continue

c---------------------------------------------------------------------------
c       diagram  <pg|v|ra><qr|t2|bg>
c---------------------------------------------------------------------------

      nrg=idra(iir,iig,l1)
      nqb=idra(jq,jb,l1)
      in=ntloc(nrg,nqb,irgsym,l1)
      call findk(orbg,orbap,orbpp,orbr,iag,jaa,jap,iar,jgot,jmax)
      if (jmax.eq.0) go to 940
      do 950 kloop=1,jmax
      l2=jgot(kloop)
      al2=dfloat(l2)
      rk=vint(jp,iig,iir,ja,l2)
      f4=(-1)**(l2+orbg+orbpp+1)
      twoe=rk*f4
      fact=(-one)**(orbr+orbg+l1)
      call sixj(al1,orbg,orbr,al2,orbpp,orbap,s6j)
c     baux_local(in)=baux_local(in)-twoe*fact*mult*s6j
      baux(in)=baux(in)-twoe*fact*mult*s6j

  950 continue
  940 continue
  910 continue

      irasym=mtbl(irsym,jasym)
      iqgsym=mtbl(jqsym,igsym)

      if(irasym.ne.iqgsym)go to 960

c----------------------------------------------------------------------------
c           diagram  <pg|v|rb><rq|t2|qg>
c----------------------------------------------------------------------------

      call findk(orbg,orbbp,orbpp,orbr,iag,jab,jap,iar,jgot,jmax)
      if (jmax.eq.0) go to 970
      do 980 kloop=1,jmax
      l2=jgot(kloop)
      al2=dfloat(l2)
      f4=(-1)**(l2+orbg+orbpp+1)
      rk=vint(jp,iig,iir,jb,l2)
      twoe=rk*f4

      call findk(orbg,orbqp,orbap,orbr,iag,jaq,jaa,iar,igot,imax)
      if(imax.eq.0)go to 980
      do 990 mloop=1,imax
      l3=igot(mloop)
      al3=dfloat(l3)
      nra=idra(iir,ja,l3)
      nqg=idra(jq,iig,l3)
      in=ntloc(nra,nqg,irasym,l3)
      fact=(two*al1+1)
      call sixj(al1,al2,al3,orbg,orbqp,orbbp,s6j1)
      call sixj(al1,al3,al2,orbr,orbpp,orbap,s6j2)
      f5=(-one)**(orba+orbp+orbb+orbq+l1+l2+l3)
c     baux_local(in)=baux_local(in)-twoe*mult*f5*s6j1*s6j2*fact
      baux(in)=baux(in)-twoe*mult*f5*s6j1*s6j2*fact

  990 continue
  980 continue
  970 continue
  960 continue


  900 continue

      do 1000 is=1,nexcit
      issym=isyme(is)
      orbs=orbe(is)
      ias=iqe(is)
      iis=ke(is)

      isbsym=mtbl(issym,jbsym)
      if(isbsym.ne.irasym)go to 1000

c------------------------------------------------------------------------
c          diagram  <pq|v|rs><rs|t2|ab>
c------------------------------------------------------------------------

      call findk(orbr,orbpp,orbs,orbqp,iar,jap,ias,jaq,jgot,jmax)
      if (jmax.eq.0) go to 1000
      do 1100 kloop=1,jmax
      l2=jgot(kloop)
      al2=dfloat(l2)
      rk=vint(jp,jq,iir,iis,l2)
      f4=(-1)**(l2+orbpp+orbqp+1)
      twoe=rk*f4

      call findk(orbr,orbap,orbs,orbbp,iar,jaa,ias,jab,igot,imax)
      if (imax.eq.0)go to 1100
      do 1200 mloop=1,imax
      l3=igot(mloop)
      al3=dfloat(l3)
      f5=(-1)**(orba+orbb+orbp+orbq)
      nra=idra(iir,ja,l3)
      nsb=idra(iis,jb,l3)
      in=ntloc(nra,nsb,irasym,l3)
      fact=two*al1+1
      call sixj(al1,al3,al2,orbr,orbpp,orbap,s6j1)
      call sixj(al1,al3,al2,orbs,orbqp,orbbp,s6j2)
c     baux_local(in)=baux_local(in)+twoe*half*mult*fact*s6j1*s6j2*f5
      baux(in)=baux(in)+twoe*half*mult*fact*s6j1*s6j2*f5

 1200 continue
 1100 continue


 1000 continue

      if(iopt.eq.0)go to 800

c----------------------------------------------------------------------------
c            diagram  <pq|v|ar><r|t1|b>
c----------------------------------------------------------------------------

      if(irsym.ne.jbsym)go to 1300
      iir=ke(ir)
      iar=iqe(ir)
      orbr=orbe(ir)
      if(orbbp.ne.orbr)go to 1300
      in=idpa(iir,jb)
      rk=vint(jp,jq,ja,iir,l1)
      f4=(-1)**(l1+orbpp+orbqp+1)
c     baux_local(in)=baux_local(in)+rk*mult*f4
      baux(in)=baux(in)+rk*mult*f4

 1300 continue


  800 continue

                                                                                
      do 1400 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      iig=kc(ig)

      ipgsym=mtbl(jpsym,igsym)

      do 1500 id=1,nocc
      idsym=isymc(id)
      orbd=orbc(id)
      iad=iqc(id)
      iid=kc(id)

      iqdsym=mtbl(jqsym,idsym)

      if(iqdsym.ne.ipgsym)go to 1500

c------------------------------------------------------------------------------
c            diagram  <gd|v|ab><pq|t2|gd>
c------------------------------------------------------------------------------

      call findk(orbd,orbbp,orbg,orbap,iad,jab,iag,jaa,jgot,jmax)
      if (jmax.eq.0)go to 1500
      do 1600 kloop=1,jmax
      l2=jgot(kloop)
      al2=dfloat(l2)
      rk=vint(iig,iid,ja,jb,l2)
      f4=(-1)**(l2+orbg+orbd+1)
      twoe=rk*f4

      call findk(orbd,orbqp,orbg,orbpp,iad,jaq,iag,jap,igot,imax)
      if(imax.eq.0) go to 1600
      do 1700 mloop=1,imax
      l3=igot(mloop)
      al3=dfloat(l3)
      npg=idra(jp,iig,l3)
      nqd=idra(jq,iid,l3)
      in=ntloc(npg,nqd,ipgsym,l3)
      fact=(two*al1+1)
      call sixj(al1,al2,al3,orbg,orbpp,orbap,s6j1)
      call sixj(al1,al2,al3,orbd,orbqp,orbbp,s6j2)
      f5=(-1)**(orba+orbp+orbb+orbq)
c     baux_local(in)=baux_local(in)+twoe*mult*half*fact*s6j1*s6j2*f5
      baux(in)=baux(in)+twoe*mult*half*fact*s6j1*s6j2*f5

 1700 continue
 1600 continue


 1500 continue

c-----------------------------------------------------------------------------
c          diagram 7 <pg|v|ab><q|t1|g>
c-----------------------------------------------------------------------------

      if(iopt.eq.0)go to 1400
      if(igsym.ne.jqsym)go to 1800
      iig=kc(ig)
      iag=iqc(ig)
      orbg=orbc(ig)
      if(orbg.ne.orbqp)go to 1800
      in=idpa(jq,iig)
      rk=vint(jp,iig,ja,jb,l1)
      f4=(-1)**(l1+orbpp+orbg+1)
c     baux_local(in)=baux_local(in)-rk*mult*f4
      baux(in)=baux(in)-rk*mult*f4

 1800 continue
 1400 continue
                                                                                
                                                                                
  60  continue
      xsum=0
      do in=1,ntmax
c     xsum=xsum+baux_local(in)*t(in)
      xsum=xsum+baux(in)*t(in)

      enddo
c     raux_local(neqn)=xsum
      raux(neqn)=xsum

  45  continue

  40  continue
  30  continue
  20  continue
  10  continue
c       call MPI_ALLREDUCE(raux_local,raux,MDIM,MPI_REAL8,MPI_SUM,
c    :                   MPI_COMM_WORLD,ierr)

c     do i= 1,MDIM
c      raux_local(i)=0.0
c     enddo

c       call MPI_ALLREDUCE(baux_local,baux,MDIM,MPI_REAL8,MPI_SUM,
c    :                   MPI_COMM_WORLD,ierr)

c     do i= 1,MDIM
c      baux_local(i)=0.0
c     enddo


c     print*,'double = raux,baux',raux,baux
      return
      end

      subroutine cizec(t,iopt,iproc)
      implicit real*8 (a-h,o-z)
      real*8 xnergy_local,ynergy_local,xnergy1_local,ynergy1_local
C
C     THIS SUBROUTINE COMPUTES THE CORRELATION ENERGY
C
      include 'mpif.h'
c     common block
      common/jvalue/orbc(MNOCC),orbe(MNEXC)
      common/index/ke(MNEXC),kc(MNOCC)
      common/info/nocc,nexcit,iType(MNBAS),iTCount(MNBAS)
      common/orbital_energy/eorb(MNBAS)
      common/parity/isymc(MNOCC),isyme(MNEXC)
      common/symmetry/mtbl(MNS,MNS)
      common/kpavlaue/iqc(MNOCC),iqe(MNEXC)
      common/idra1/idra(MNBAS,MNBAS,0:MXV)
      common/idpa1/idpa(MNBAS,MNBAS)
      common/skip/nnskip(0:MXV,2),n2skip,nsing
c     local dimension
      dimension kgot(MXV+1),igot(MXV+1),t(MDIM)
c     data
      data zero,half,one,two/0.0d0,0.50d0,1.0d0,2.0d0/
      call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)



      xnergy=0
      ynergy=0
      xnergy1=0
      ynergy1=0
      xnergy_local=0
      ynergy_local=0
      xnergy1_local=0
      ynergy1_local=0

      do 10 ia=1,nocc
      iasym=isymc(ia)

      do 20 ip=iproc+1,nexcit,nproc
      ipsym=isyme(ip)

      iapsym=mtbl(ipsym,iasym)

      do 30 ib=1,nocc
      ibsym=isymc(ib)

      do 40 iq=1,nexcit
      iqsym=isyme(iq)

      ibqsym=mtbl(ibsym,iqsym)
      iaqsym=mtbl(iasym,iqsym)

      orba=orbc(ia)
      orbb=orbc(ib)
      orbp=orbe(ip)
      orbq=orbe(iq)

      iaa=iqc(ia)
      iab=iqc(ib)
      iap=iqe(ip)
      iaq=iqe(iq)

      iia=kc(ia)
      iib=kc(ib)
      iip=ke(ip)
      iiq=ke(iq)

      if(ibqsym.ne.iapsym)go to 40

c------------------------------------------------------------------------
c      energy contributed fron <pq|v|ab><pq|t2|ab>
c------------------------------------------------------------------------

      call findk(orba,orbp,orbb,orbq,iaa,iap,iab,iaq,kgot,
     $kmax)
      if(kmax.eq.0)go to 50
      do 60 kloop=1,kmax
      l1=kgot(kloop)
      al1=dfloat(l1)
      rk=vint(iia,iib,iip,iiq,l1)
      f4=(-1)**(l1+orba+orbb+1)
      nap=idra(iip,iia,l1)
      nbq=idra(iiq,iib,l1)
      in=ntloc(nap,nbq,iapsym,l1)
      f5=(-1)**(orba+orbb+orbp+orbq)
c     xnergy=xnergy+rk*f5*t(in)*half*f4/(two*l1+1)
      xnergy_local=xnergy_local+rk*f5*t(in)*half*f4/(two*l1+1)


c---------------------------------------------------------------------------
c        energy contributed from <pq|v|ab><pq|t2|ba>
c---------------------------------------------------------------------------

      call findk(orba,orbq,orbb,orbp,iaa,iaq,iab,iap,igot,
     $imax)
      if(imax.eq.0)go to 60
      do 70 nloop=1,imax
      l2=igot(nloop)
      al2=dfloat(l2)
      naq=idra(iiq,iia,l2)
      nbp=idra(iip,iib,l2)
      in=ntloc(naq,nbp,iaqsym,l2)
      call sixj(al1,orbp,orba,al2,orbq,orbb,s6j)
c     ynergy=ynergy+rk*t(in)*s6j*half*f5*f4
      ynergy_local=ynergy_local+rk*t(in)*s6j*half*f5*f4

  70  continue
  60  continue
  50  continue

      if(iopt.ne.0)then

      if(orba.ne.orbp)go to 80 
      if(orbb.ne.orbq)go to 80
      if(iasym.ne.ipsym)go to 80
      if(ibsym.ne.iqsym)go to 80

      in1=idpa(iip,iia)
      in2=idpa(iiq,iib)

c----------------------------------------------------------------------------
c
c----------------------------------------------------------------------------

      rk=vint(iia,iib,iip,iiq,0)
      f4=(-1)**(orba+orbb+1)
      f5=dsqrt((two*orba+1)*(two*orbb+1))
c     xnergy1=xnergy1+rk*f4*f5*half*t(in1)*t(in2)
      xnergy1_local=xnergy1_local+rk*f4*f5*half*t(in1)*t(in2)

c---------------------------------------------------------------------------
c
c---------------------------------------------------------------------------

      call findk(orba,orbq,orbb,orbp,iaa,iaq,iab,iap,kgot,
     $kmax)
      if(kmax.eq.0)go to 80
      do 90 kloop=1,kmax
      kk=kgot(kloop)
      ak=dfloat(kk)
      rk=vint(iia,iib,iiq,iip,kk)
      f4=(-1)**(kk+orba+orbb+1)
      f5=(-1)**(orba+orbb+kk)
c     ynergy1=ynergy1+rk*f4*f5*half*t(in1)*t(in2)
      ynergy1_local=ynergy1_local+rk*f4*f5*half*t(in1)*t(in2)

 90   continue
 80   continue

      endif

  40  continue
  30  continue
  20  continue
  10  continue
      call MPI_ALLREDUCE(xnergy_local,xnergy,1,MPI_REAL8,MPI_SUM,
     :                   MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(ynergy_local,ynergy,1,MPI_REAL8,MPI_SUM,
     :                   MPI_COMM_WORLD,ierr)

       xnergy_local=0.0
       ynergy_local=0.0
      call MPI_ALLREDUCE(xnergy1_local,xnergy1,1,MPI_REAL8,MPI_SUM,
     :                   MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(ynergy1_local,ynergy1,1,MPI_REAL8,MPI_SUM,
     :                   MPI_COMM_WORLD,ierr)

       xnergy1_local=0.0
       ynergy1_local=0.0

 
      energy=xnergy+ynergy+xnergy1+ynergy1

      if (iproc .eq. 0) then
      write(STDOUT,9999)
      write(STDOUT,9996)xnergy,ynergy
      write(STDOUT,9995)xnergy+ynergy
      open(42,file='energy.out',form='formatted')
      write(42,9999)
      write(42,9996)xnergy,ynergy
      write(42,9995)xnergy+ynergy
      close(42)
      if(iopt.ne.0)then
      write(STDOUT,9998)
      write(STDOUT,9996)xnergy1,ynergy1
      write(STDOUT,9995)xnergy1+ynergy1
      write(STDOUT,9997)
      write(STDOUT,9995)energy

      endif
      endif

 9999 format(//,2x,'correlation energy (from t2) bt=-a',/)
 9998 format(//,2x,'correlation energy (from t1) bt=-a',/)
 9997 format(//,2x,'correlation energy (t1+t2) bt=-a',/)
 9996 format(/,2x,'direct   part                  =',f20.11,/
     :         2x,'exchange part                  =',f20.11)
 9995 format(/,2x,'correlation energy(au)         =',f20.11)

      return
      end


c***************************************************************************
c                                                                          *
c      SUBROUTINE   SYMM                                                   *
c                                                                          *
c     THIS SUBROUTINE SETS UP THE T1 AND T2 EQUATION INDEX                 *
c                                                                          *
c***************************************************************************

     
      subroutine symm(iopt,ntmax,isolv,nprint,iproc)

      implicit real*8 (a-h,o-z)
      include 'mpif.h'
      common/iallparity/iall
      common/jvalue/orbc(MNOCC),orbe(MNEXC)
      common/index/ke(MNEXC),kc(MNOCC)
      common/info/nocc,nexcit,iType(MNBAS),iTCount(MNBAS)
      common/parity/isymc(MNOCC),isyme(MNEXC)
      common/symmetry/mtbl(MNS,MNS)
      common/kpavlaue/iqc(MNOCC),iqe(MNEXC)
      common/maxj/jmax
      common/idra1/idra(MNBAS,MNBAS,0:MXV)
      common/idpa1/idpa(MNBAS,MNBAS)
      common/skip/nnskip(0:MXV,2),n2skip,nsing

c     local dimension
      dimension kgot(MXV+1),ieven(0:MXV),iodd(0:MXV)
      dimension indx(MDIM),indx1(MDIM)
      dimension ncheck(MDIM)

      do i=1, MDIM
         ncheck(i)=0
      enddo

c---------------------------------------------------------------------------
c     set the equation index for T1
c---------------------------------------------------------------------------

      if(iopt.ne.0)then
      if (iproc .eq. 0) then
         write(STDOUT,9999)
         write(STDOUT,9997)
      endif

      do i=1,MNBAS
      do j=1,MNBAS
      do k=0,MXV
      idra(i,j,k)=0

      enddo
      enddo
      enddo

      nsing=0

      do 1 ia=1,nocc
      iasym=isymc(ia)

      do 1 ip=1,nexcit
      ipsym=isyme(ip)
      iapsym=mtbl(ipsym,iasym)

      orba=orbc(ia)
      orbp=orbe(ip)

      iaa=iqc(ia)
      iap=iqe(ip)

      iia=kc(ia)
      iip=ke(ip)

      if(iapsym.ne.1)go to 1
      if(orba.ne.orbp)go to 1

      nsing=nsing+1

      idpa(iip,iia)=nsing

      indx(nsing)=nsing
      indx1(nsing)=nsing


   1  continue

      endif

      if (iproc .eq. 0) write(STDOUT,9998)

c----------------------------------------------------------------------------
c    set the skip information
c----------------------------------------------------------------------------

      do 2 i=0,MXV
      ieven(i)=0
      iodd(i)=0
      nnskip(i,1)=0
      nnskip(i,2)=0
    2 continue



      do 3 ia=1,nocc
      iasym=isymc(ia)

      do 3 ip=1,nexcit
      ipsym=isyme(ip)

      iapsym=mtbl(ipsym,iasym)

      orba=orbc(ia)
      orbp=orbe(ip)

      iaa=iqc(ia)
      iap=iqe(ip)

      iia=kc(ia)
      iip=ke(ip)

      call findk(orba,orbp,orba,orbp,iaa,iap,iaa,iap,kgot,
     $kmax)
      if(kmax.eq.0)go to 3
      do 4 nloop=1,kmax
      kk=kgot(nloop)
      if(iapsym.eq.1)then
        ieven(kk)=ieven(kk)+1
        idra(iip,iia,kk)=ieven(kk)
      else
        iodd(kk)=iodd(kk)+1
        idra(iip,iia,kk)=iodd(kk)
      endif

  4   continue
  3   continue


      n2skip=0

      if(iall.eq.0)then

      do 5 i=2,MXV,2
      nnskip(i,2)=nnskip(i-2,2)+ieven(i-2)*(ieven(i-2)+1)/2
  5   continue

      do  i=0,JMAX,2
      n2skip=n2skip+ieven(i)*(ieven(i)+1)/2
      enddo

      do 6 i=3,MXV-2,2
      nnskip(i,1)=nnskip(i-2,1)+iodd(i-2)*(iodd(i-2)+1)/2
  6   continue

        else

      do 9 i=1,MXV
      nnskip(i,2)=nnskip(i-1,2)+ieven(i-1)*(ieven(i-1)+1)/2
  9   continue
      
      do  i=0,MXV 
      n2skip=n2skip+ieven(i)*(ieven(i)+1)/2
      enddo

      do 10 i=1,MXV
      nnskip(i,1)=nnskip(i-1,1)+iodd(i-1)*(iodd(i-1)+1)/2
 10   continue

       endif

      nn=nsing

      if (iproc .eq. 0) write(STDOUT,9996)

c--------------------------------------------------------------------------
c    set the equation index for T2
c--------------------------------------------------------------------------

      do 15 ia=1,nocc
      iasym=isymc(ia)

      do 20 ip=1,nexcit
      ipsym=isyme(ip)

      iapsym=mtbl(ipsym,iasym)

      do 30 ib=1,nocc
      ibsym=isymc(ib)

      do 40 iq=1,nexcit
      iqsym=isyme(iq)

      ibqsym=mtbl(ibsym,iqsym)

      if(ibqsym.ne.iapsym)go to 40

      orba=orbc(ia)
      orbb=orbc(ib)
      orbp=orbe(ip)
      orbq=orbe(iq)

      iaa=iqc(ia)
      iab=iqc(ib)
      iap=iqe(ip)
      iaq=iqe(iq)

      iia=kc(ia)
      iib=kc(ib)
      iip=ke(ip)
      iiq=ke(iq)

      call findk(orba,orbp,orbb,orbq,iaa,iap,iab,iaq,kgot,
     $kmax)
      if(kmax.eq.0)go to 40
      do 45 nloop=1,kmax
      kk=kgot(nloop)
      nap=idra(iip,iia,kk)
      nbq=idra(iiq,iib,kk)
      call iloc(iip,iia,iiq,iib,ipa,iqb)
      if(ipa.lt.iqb)go to 45
      nn=nn+1
      neqn=ntloc(nap,nbq,iapsym,kk)
      if (ncheck(neqn).eq.1) then
       if(iproc.eq.0)write(STDOUT,*)'double count: neqn=', neqn
      stop
      endif

      ncheck(neqn)=1

      indx(nn)=nn
      indx1(nn)=neqn

  45  continue
  40  continue
  30  continue
  20  continue
  15  continue

      ntmax=nn

      if(isolv.ne.0) then
         isolv = 0
      endif

9999  format(//,'SINGLE EXCITATIONS',/)
9998  format(//,'DOUBLE EXCITATIONS',/)
9997  format(3x,'IN','    IA ','  IP ','J(IA)',' J(IP)',/)
9996  format(3x,'IN','    IA ','  IB ','  IP ',' IQ' ,3x,'JA',
     $         4x,'JB',5x,'JP',5x,'JQ',5x,'ASM',2x,'BSM',2x,
     $         'PSM',2x,'QSM',2x,'KMAX',/)
9995  format(3i5,2(f5.1,2x))
9994  format(5i5,4(f5.1,2x),8(i4,1x))

      return
      end


      function ntloc(iap,ibq,iapsym,kk)
      implicit real*8 (a-h,o-z)
      include 'mpif.h'
c     common block
      common/skip/nnskip(0:MXV,2),n2skip,nsing
      mm=max0(iap,ibq)
      nn=min0(iap,ibq)
      if(iapsym.eq.1)then
      ntloc=mm*(mm-1)/2+nn+nnskip(kk,2)+nsing
      else
      ntloc=mm*(mm-1)/2+nn+nnskip(kk,1)+n2skip+nsing
      endif
      return
      end

      subroutine mInit()
c----------------------------------------------------------------------
c     The integrals vanish for each 2'nd multipole moment, so that
c     the tables just store the nonvanishing contributions. Here,
c     the index is determined, into which each of the multipole moments
c     is mapped.
c----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      include 'mpif.h'
      common/mind/mIndex(0:20)
      integer icount, m
      icount = 1
      do m = 0, 19, 2
         mIndex(m) = icount
         mIndex(m+1) = icount
         icount = icount + 1
      enddo
      return
      end
c***************************************************************************
c                                                                          *
c      SUBROUTINE TWOINT                                                   *
c                                                                          *
c     TWO ELECTRON INTEGRALS ARE GENERATED AND STORED                      *              
c     IN MULLIKEN ORDER (I.GE.J , K.GE.L AND IJ.GE.KL                      *
c***************************************************************************


      subroutine twoint(nbasis,itwoe)

      IMPLICIT REAL*8 (A-H,O-Z)
C     include 'mpif.h'
C     COMMON BLOCK
      COMMON/INFO/NOCC,NEXCIT,Itype(MNBAS),Itcount(MNBAS)
      COMMON/SYMMETRY/MTBL(MNS,MNS)
      COMMON/MIND/MINDEX(0:20)
      COMMON/SYMINFO/ORB(MNBAS),IIQ(MNBAS),IPARITY(MNBAS)

      COMMON/INTGRL/V2E(N2INT,0:MXV)

c     common block
       common/iallparity/iall

c     local dimension
      dimension igot(MXV+1),jgot(MXV+1)
c     data
      data zero,half,one,two/0.0d0,0.50d0,1.0d0,2.0d0/
c     initialize the vint array
 
      print*,'I am here 1'
      DO I=1,N2INT
      DO J=0,MXV
      IGOT(J)=0
      JGOT(J)=0
      V2E(I,J)= ZERO
      ENDDO
      ENDDO

      print*,'I am here 2: calling Info_2e'
      CALL INFO_2E
      print*,'I am here 3: came out from Info_2e'

      NBASIS=NOCC+NEXCIT

      if(itwoe.eq.1)goto 11
      if(itwoe.eq.2)goto 11
          GOTO 12

 11   CONTINUE

      open(NTWOE,file='twoea.dat',form='unformatted',
     &     status='unknown')

 12   CONTINUE

      print*,'I am here 4'

      DO 100 I=1,NBASIS
      DO 100 J=1,I
      Do 100 K=1,NBASIS
      Do 100 L=1,NBASIS

      ISYM=IPARITY(I)
      JSYM=IPARITY(J)
      KSYM=IPARITY(K)
      LSYM=IPARITY(L)

      IJPAR=MTBL(ISYM,JSYM)
      KLPAR=MTBL(KSYM,LSYM)

      IF(IJPAR.NE.KLPAR)GO TO 100

      ORBI=ORB(I)
      ORBJ=ORB(J)
      ORBK=ORB(K)
      ORBL=ORB(L)

      II=IIQ(I)
      JJ=IIQ(J)
      KK=IIQ(K)
      LL=IIQ(L)

      CALL FINDK(ORBI,ORBK,ORBJ,ORBL,II,KK,JJ,LL,IGOT,IMAX)
      IF(IMAX.EQ.0)GO TO 200
      DO 210 NLOOP=1,IMAX
      L1=IGOT(NLOOP)

      INDX=IGET2(I,K,J,L)

      IF(INDX.LE.0)GOTO 130
      IF(INDX.GT.N2INT) GOTO 130
            GOTO 140

 130  CONTINUE

      WRITE(*,*)'INVALID INDEX: CHECK IGET2',INDX
      WRITE(*,*)I,J,K,L
      STOP

 140  CONTINUE

      if(itwoe.eq.0)goto 141
      if(itwoe.eq.1)goto 141
           goto 142
 141  continue

      RK=0.0
      VAL1=0.0  
      VAL2=0.0
      VAL3=0.0

c-----------------------------------------------------------------------------
c     Coulomb interaction part
c-----------------------------------------------------------------------------

      RK=SLATER_NEW(I,J,K,L,ORBI,ORBJ,ORBK,ORBL,
     :II,JJ,KK,LL,L1)

      IF(IALL.EQ.2) GOTO 220
      IF(IALL.EQ.3) GOTO 220
          GOTO 230

 220  CONTINUE

C       CALL BREIT(I,J,K,L,ORBI,ORBJ,ORBK,ORBL,
C     : II,JJ,KK,LL,L1,VAL,IALL)



c-----------------------------------------------------------------------------
c      Breit interaction part with n=k
c-----------------------------------------------------------------------------

       VAL1=SLATER_NEWBR1(I,J,K,L,ORBI,ORBJ,ORBK,ORBL,
     : II,JJ,KK,LL,L1)

c----------------------------------------------------------------------------
c      Breit interaction part with n=k-1
c----------------------------------------------------------------------------

        VAL2=SLATER_NEWBR2(I,J,K,L,ORBI,ORBJ,ORBK,ORBL,
     : II,JJ,KK,LL,L1)

c----------------------------------------------------------------------------
c      Breit interaction part with n=k+1
c----------------------------------------------------------------------------

        VAL3=SLATER_NEWBR3(I,J,K,L,ORBI,ORBJ,ORBK,ORBL,
     : II,JJ,KK,LL,L1)

 
  230 CONTINUE

      V2E(INDX,L1)=RK+VAL1+VAL2+VAL3


      if(itwoe.eq.1)then
      write(NTWOE)V2E(INDX,L1)
      endif

 142  continue

      if(itwoe.eq.2)then
      read(NTWOE)V2E(INDX,L1)
      endif

 210  CONTINUE
 200  CONTINUE

 100  CONTINUE

      RETURN
      END

      real*8 function vint(i,j,k,l,m)
      implicit real*8 (a-h,o-z)
      include 'mpif.h'

      COMMON/INFO/NOCC,NEXCIT,ITYPE(MNBAS),ITCOUNT(MNBAS)
      COMMON/MIND/MINDEX(0:20)
      COMMON/SYMMETRY/MTBL(MNS,MNS)
      COMMON/SYMINFO/ORB(MNBAS),IIq(MNBAS),IPARITY(MNBAS)
      COMMON/INTGRL/V2E(N2INT,0:MXV)

      DIMENSION KGOT(MXV+1)

c-------------------------------------------------------------------
c     check validity of the given multipole moment
c-------------------------------------------------------------------

      ISYM=IPARITY(I)
      JSYM=IPARITY(J)
      KSYM=IPARITY(K)
      LSYM=IPARITY(L)

      IJPAR=MTBL(ISYM,JSYM)
      KLPAR=MTBL(KSYM,LSYM)

      IF(IJPAR.NE.KLPAR)THEN
      VINT=0.0D0

      ELSE

      INDX=IGET2(I,K,J,L)

      IF(INDX.LE.0) GOTO 100
      IF(INDX.GT.N2INT) GOTO 100
          GOTO 110

 100  CONTINUE

      WRITE(*,*)'INVALID INDEX: CHECK IGET2',INDX
      WRITE(*,*)I,J,K,L
      STOP

 110  CONTINUE

      VINT = V2E(INDX,M)

      ENDIF

      RETURN
      END
       
      SUBROUTINE INFO_2E
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/INFO/NOCC,NEXCIT,ITYPE(MNBAS),ITCOUNT(MNBAS)
      COMMON/PARITY/ISYMC(MNBAS),ISYME(MNBAS)
      COMMON/SYMINFO/ORBJ(MNBAS),IIQ(MNBAS),IPARITY(MNBAS)
      COMMON/INFO1/NORBSM(MNBAS),LORBCM(MNBAS,2),NUMORB(2)
      COMMON/INFO2/INT2AL(0:204)
      DIMENSION MULT(2,2)
      DATA MULT/1,2,2,1/

      NBASIS=NOCC+NEXCIT

      DO 100 I=1,NBASIS
      NORBSM(I)=IPARITY(I)
      IF(NORBSM(I).EQ.3)NORBSM(I)=1
      IF(NORBSM(I).EQ.4)NORBSM(I)=2
      IF(NORBSM(I).EQ.5)NORBSM(I)=1
      IF(NORBSM(I).EQ.6)NORBSM(I)=2
  100 CONTINUE

      DO 110 I=1,2
      NUMORB(I)=0
  110 CONTINUE

      DO 200 J=1,2
      DO 210 I=1,NBASIS
      II=NORBSM(I)
      IF(II.NE.J)GO TO 210
      NUMORB(II)=NUMORB(II)+1
  210 CONTINUE
  200 CONTINUE

      JJ=0
      DO 300 I=1,2
      II=0
      DO 310 J=1,NBASIS
      IF(NORBSM(J).NE.I)GOTO 310
      JJ=JJ+1
      II=II+1
      LORBCM(J,1)=JJ
      LORBCM(J,2)=II
  310 CONTINUE
  300 CONTINUE

      DO 280 L1=1,2
      IND=(L1-1)*L1*(2*L1-1)/6
      DO 280 L2=1,L1
      DO 280 IS1=1,L1
      INDX=IND+L1*(IS1-1)+L2
      IS2=MULT(L1,L2)
      IS2=MULT(IS1,IS2)
      IF((L1.NE.L2.AND.L2.GE.IS2).OR.(L1.EQ.L2.AND.IS2.LE.IS1))THEN
      NH1=NUMORB(L1)
      IF(L1.EQ.IS1)THEN
      J1=(NH1+1)*NH1/2
      IND1=(L1+1)*L1/2
      ELSE
      NH2=NUMORB(IS1)
      J1=NH1*NH2
      IND1=(L1-1)*L1/2+IS1
      END IF
      NH3=NUMORB(L2)
      IF(L2.EQ.IS2)THEN
      J2=(NH3+1)*NH3/2
      IND2=(L2+1)*L2/2
      ELSE
      NH4=NUMORB(IS2)
      IND2=(L2-1)*L2/2+IS2
      J2=NH3*NH4
      END IF
      IF(IND1.EQ.IND2)THEN
      INT2AL(INDX)=(J1+1)*J1/2
      ELSE
      INT2AL(INDX)=J1*J2
      END IF
      ELSE
      INT2AL(INDX)=0
      END IF

      IF(INDX.GT.204)THEN
      WRITE(*,*)
     :'INCREASE THE DIMENSION OF INT2AL FROM 204 TO ',INDX
      STOP
      ENDIF

  280 CONTINUE

      INT2AL(0)=0
      DO 290 I=1,204
      INT2AL(I)=INT2AL(I)+INT2AL(I-1)
  290 CONTINUE
      RETURN
      END

      INTEGER FUNCTION IGET2 (I,J,K,L)
      IMPLICIT INTEGER (A-Z)
      COMMON/INFO2/INT2AL(0:204)
      COMMON/INFO1/NORBSM(MNBAS),LORBCM(MNBAS,2),NUMORB(2)
      DIMENSION MULT(2,2)
      DATA MULT/1,2,2,1/

      IGET2=0

      ISM1=NORBSM(I)
      ISM2=NORBSM(J)
      ISM3=NORBSM(K)
      ISM4=NORBSM(L)
      M1=MULT(ISM1,ISM2)
      M2=MULT(ISM3,ISM4)
      IF(M1.NE.M2)GO TO 99
      I1=LORBCM(I,1)
      II1=LORBCM(I,2)
      I2=LORBCM(J,1)
      II2=LORBCM(J,2)
      IF(I1.LT.I2)THEN
      ITP=I1
      I1=I2
      I2=ITP
      ITP=II1
      II1=II2
      II2=ITP
      ITP=ISM1
      ISM1=ISM2
      ISM2=ITP
      END IF
      I3=LORBCM(K,1)
      II3=LORBCM(K,2)
      I4=LORBCM(L,1)
      II4=LORBCM(L,2)
      IF(I3.LT.I4)THEN
      ITP=I3
      I3=I4
      I4=ITP
      ITP=II3
      II3=II4
      II4=ITP
      ITP=ISM3
      ISM3=ISM4
      ISM4=ITP
      END IF
      IF(I1.GT.I3)GOTO 50
      IF(I1.EQ.I3.AND.I2.GT.I4)GOTO 50
      ITP=II1
      II1=II3
      II3=ITP
      ITP=ISM1
      ISM1=ISM3
      ISM3=ITP
      ITP=II2
      II2=II4
      II4=ITP
      ITP=ISM2
      ISM2=ISM4
      ISM4=ITP
 50   CONTINUE
      IOFSET=INT2AL((ISM1-1)*ISM1
     &      *(2*ISM1-1)/6+ISM1*(ISM2-1)+ISM3-1)
      IF(ISM1.EQ.ISM2)THEN
      J1=(II1-1)*II1/2+II2
      J2=(II3-1)*II3/2+II4
      IF(ISM1.EQ.ISM3)THEN
      IGET2=(J1-1)*J1/2+J2+IOFSET
      ELSE
      IGET2=(J1-1)*(NUMORB(ISM3)+1)*NUMORB(ISM3)/2+J2+IOFSET
      END IF
      ELSE
      J1=(II1-1)*NUMORB(ISM2)+II2
      J2=(II3-1)*NUMORB(ISM4)+II4
      IF(ISM1.EQ.ISM3)THEN
      IGET2=(J1-1)*J1/2+J2+IOFSET
      ELSE
      IGET2=(J1-1)*NUMORB(ISM3)*NUMORB(ISM4)+J2+IOFSET
      END IF
      END IF
 99   CONTINUE
c     print*, I,J,K,L,IGET2
      RETURN
      END

      SUBROUTINE PACK(I,J,K,L,IJ,KL)
      implicit real*8 (a-h,o-z)
      include 'mpif.h'
      common/info/nocc,nexcit,iType(MNBAS),iTCount(MNBAS)
      IJ=I*(I-1)/2+J
      KL=(NOCC+NEXCIT)*(K-1)+L
      RETURN
      END

 

      subroutine jacobi1(told,t,iopt,isolv,ntmax,nonlin,
     :                  irows,irowe,iproc,nproc,dspl,delrow)
      implicit real*8 (a-h,o-z)
      include 'mpif.h'
      common/scratch/b(MDIM),ad(MDIM),baux(MDIM)
      common/scratch2/h2d(MDIM),raux(MDIM),ap(MDIM)                 
C
C     THIS SUBROITINE SOLVES THE BT=A FOR LINEARIZED
C     COUPLED CLUSTER EQUATION ON FLY USING JACOBI
C     SCHEME
C
c     local dimension
      integer dspl(NPMAX),delrow(NPMAX),drow,error
      dimension told(MDIM),t(MDIM),caux(MDIM)
      start_time = 0
      end_time = 0

c      start_time = MPI_WTIME()
      if (iproc .eq. 0) write(STDOUT,1000)
c--------------------------------------------------------------------
c     Start linear iteration
c--------------------------------------------------------------------
      nit=0
   30 nit=nit+1

      if(iopt.ne.0)then
      call singles(caux,told,ntmax,isolv,irows,irowe)
      endif

      call doubles(caux,told,ntmax,iopt,isolv,irows,irowe)

c     check convergence
      do i=irows,irowe
      t(i)=(b(i)-caux(i))/ad(i)
      enddo

c      call MPI_ALLGATHERV(t(irows),drow,MPI_DOUBLE_PRECISION,t,
c     >delrow,dspl,MPI_DOUBLE_PRECISION, MPI_COMM_WORLD,error)

      eps=0
      do i=1,ntmax
      eps=eps+dabs(told(i)-t(i))
      told(i)=t(i)
      enddo

      conv=eps/ntmax
      if(iproc.eq.0)then
      write(STDOUT,9997)nit,conv
      endif

      if(conv.le.SMALL)go to 40
      if(nit.ge.NITMAX)go to 50
      go to 30

c-----------------------------------------------------------------------
c     End linear iteration
c-----------------------------------------------------------------------

 50   if (iproc .eq. 0) write(STDOUT,9999)
      go to 60
 40   if (iproc .eq. 0) write(STDOUT,9998)nit
 60   continue

c      end_time= MPI_WTIME()

      if(iproc.eq.0)then
      open(NTFILE,file='to.dat',form='unformatted',
     &     status='unknown')
      rewind(NTFILE)
      write(NTFILE)ntmax,(t(i),i=1,ntmax)
      endif

      if (iproc .eq. 0) write(STDOUT,1002)end_time-start_time
   42 continue

      call cizec(t,iopt,iproc)

 998  format(/,1x,'diagonal part of b-matrix',//)
 999  format(4f20.11)
1000  format(/,1x,'Jacobi iteration begins',//)
1001  format(/,1x,'solution vectors',//)
1002  format(/,1x,'Time elapsed=',2f20.11//)
9990  format(//,1x,'End of linear CC',//)
9991  format(//,1x,'Non-linear Iteration begins',//)
9996  format(//,1x,'Iteration begins',//)
9997  format(1x,'Iteration no=',i4,5x,'eps=',f20.11)
9998  format(1x,' Iteration converges at',I4,' cycle')
9999  format(1x,' Convergence failure; Programme terminates')

      return
      end



      subroutine cizec2(t,iopt,iproc)
      implicit real*8 (a-h,o-z)
      include 'mpif.h'
      real*8 xnergy_local,ynergy_local,xnergy1_local,ynergy1_local
c     common block

C     THIS SUBROUTINE COMPUTES THE
C     SECOND ORDER CORRELATION ENERGY
C     I.E. MP2
C
      common/jvalue/orbc(MNOCC),orbe(MNEXC)
      common/index/ke(MNEXC),kc(MNOCC)
      common/info/nocc,nexcit,iType(MNBAS),iTCount(MNBAS)
      common/orbital_energy/eorb(MNBAS)
      common/parity/isymc(MNOCC),isyme(MNEXC)
      common/symmetry/mtbl(MNS,MNS)
      common/kpavlaue/iqc(MNOCC),iqe(MNEXC)
      common/idra1/idra(MNBAS,MNBAS,0:MXV)
      common/idpa1/idpa(MNBAS,MNBAS)
      common/skip/nnskip(0:MXV,2),n2skip,nsing
c     local dimension
      dimension kgot(MXV+1),igot(MXV+1),t(MDIM)
c     data
      data zero,half,one,two/0.0d0,0.50d0,1.0d0,2.0d0/
      call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)


      xnergy=0
      ynergy=0
      xnergy1=0
      ynergy1=0
      xnergy_local=0
      ynergy_local=0
      xnergy1_local=0
      ynergy1_local=0

      do 10 ia=1,nocc
      iasym=isymc(ia)
      do 20 ip=iproc+1,nexcit,nproc
      ipsym=isyme(ip)
      iapsym=mtbl(ipsym,iasym)
      do 30 ib=1,nocc
      ibsym=isymc(ib)
      do 40 iq=1,nexcit
      iqsym=isyme(iq)
      ibqsym=mtbl(ibsym,iqsym)
      iaqsym=mtbl(iasym,iqsym)
      if(ibqsym.ne.iapsym)go to 40
      orba=orbc(ia)
      orbb=orbc(ib)
      orbp=orbe(ip)
      orbq=orbe(iq)
      iaa=iqc(ia)
      iab=iqc(ib)
      iap=iqe(ip)
      iaq=iqe(iq)
      iia=kc(ia)
      iib=kc(ib)
      iip=ke(ip)
      iiq=ke(iq)
      call findk(orba,orbp,orbb,orbq,iaa,iap,iab,iaq,kgot,
     $kmax)
      if(kmax.eq.0)go to 40
      do 50 kloop=1,kmax
      l1=kgot(kloop)
      al1=dfloat(l1)
      rk=vint(iia,iib,iip,iiq,l1)
      f4=(-1)**(l1+orba+orbb+1)
      nap=idra(iip,iia,l1)
      nbq=idra(iiq,iib,l1)
      in=ntloc(nap,nbq,iapsym,l1)
      f5=(-1)**(orba+orbb+orbp+orbq)
c     xnergy=xnergy+rk*f4*f5*t(in)*half/(two*l1+1)
      xnergy_local=xnergy_local+rk*f4*f5*t(in)*half/(two*l1+1)
      call findk(orba,orbq,orbb,orbp,iaa,iaq,iab,iap,igot,
     $imax)
      if(imax.eq.0)go to 60
      do 70 nloop=1,imax
      l2=igot(nloop)
      al2=dfloat(l2)
      naq=idra(iiq,iia,l2)
      nbp=idra(iip,iib,l2)
      in=ntloc(naq,nbp,iaqsym,l2)
      call sixj(al1,orbp,orba,al2,orbq,orbb,s6j)
c      ynergy=ynergy+rk*f4*t(in)*s6j*half*f5
C  The sign is changed by Bijaya Sahoo
c     ynergy=ynergy-rk*f4*t(in)*s6j*half*f5
      ynergy_local=ynergy_local-rk*f4*t(in)*s6j*half*f5
  70  continue
  60  continue
  50  continue
  40  continue
  30  continue
  20  continue
  10  continue
      call MPI_ALLREDUCE(xnergy_local,xnergy,1,MPI_REAL8,MPI_SUM,
     :                   MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(ynergy_local,ynergy,1,MPI_REAL8,MPI_SUM,
     :                   MPI_COMM_WORLD,ierr)

       xnergy_local=0.0
       ynergy_local=0.0

      if (iproc .eq. 0) then
         write(STDOUT,9999)
         write(STDOUT,9998)xnergy,ynergy
         write(STDOUT,9997)xnergy+ynergy
      endif
 9999 format(//,2x,'second order correlation energy ',/)
 9998 format(/,2x,'direct   part                  =',f20.11,/
     :         2x,'exchange part                  =',f20.11)
 9997 format(/,2x,'correlation energy(au)         =',f20.11)
      return
      end

      integer function idis(i,j)
c     implicit real*8 (a-h,o-z)
      include 'mpif.h'
c     common block
      idis=max0(i,j)*(max0(i,j)-1)/2+min0(i,j)
      return
      end

      integer function indxr(i)
      implicit real*8 (a-h,o-z)
      include 'mpif.h'
c     common block
      common/info/nocc,nexcit,iType(MNBAS),iTCount(MNBAS)
      nbasis=nocc+nexcit
      indxr=nbasis*(i-1)
      return
      end

      integer function imslocHPPH(ia,ip,iq,ib,l)
      implicit none
      include 'mpif.h'
      integer ia,ip,iq,ib,l,iab,ipq
c     pack the holes:
      iab=(ia-1)*NHO + ib  

c     pack the particles:
      ipq=(ip-1)*NPO + iq

c     pack both together with the angular part:
      imslocHPPH=((iab-1)*NPO*NPO+(ipq-1))*(MXV+1)+l+1
      return
      end


      integer function imslocHHHH(ia,ib,ig,id,l)
      implicit none
      include 'mpif.h'
      integer ia,ib,ig,id,l,iab,igd
c     pack the first pair
      iab=(ia-1)*NHO + ib
      
c     pack the second pair
      igd=(ig-1)*NHO + id
      
c     pack both together with the angular part:
      imslocHHHH=((iab-1)*NHO*NHO+(igd-1))*(MXV+1)+l+1
      return
      end

 
      subroutine vbar(t,iopt)
      implicit real*8 (a-h,o-z)
      include 'mpif.h'
      real*8 aims40_local,aims41_local,aims42_local,aims43_local
     :,aims31_local,aims32_local,aims32,aims39_local,aims39

C
C     THIS SUBROUTINE COMPUTES THE EFFECTIVE 
C     TWO-ELECTRON INTEGRALS USED IN CMAT2
C     [THIS IS THE INTERMEDIATE STORAGE FOR
C      FOR THE TWO-BODY EFFECTIVE INTEGRALS
C     USED TO SOLVE CC-EQUATION EFFICIENTLY.]
C
      common/jvalue/orbc(MNOCC),orbe(MNEXC)
      common/index/ke(MNEXC),kc(MNOCC)
      common/info/nocc,nexcit,iType(MNBAS),iTCount(MNBAS)
      common/scratch2/h2d(MDIM),raux(MDIM),ap(MDIM)                 
      common/orbital_energy/eorb(MNBAS)
      common/parity/isymc(MNOCC),isyme(MNEXC)
      common/symmetry/mtbl(MNS,MNS)
      common/kpavlaue/iqc(MNOCC),iqe(MNEXC)
      common/idra1/idra(MNBAS,MNBAS,0:MXV)
      common/idpa1/idpa(MNBAS,MNBAS)
      common/skip/nnskip(0:MXV,2),n2skip,nsing
c
      common/vb/aims31(IMSHPPH),
     :aims32(IMSHPPH),
     :aims39(IMSHPPH),
     :aims40(IMSHHHH),
     :aims41(IMSHPPH),
     :aims42(IMSHPPH),
     :aims43(IMSHPPH)
c
c     local dimension
      dimension jgot(MXV+1),lgot(MXV+1),mgot(MXV+1),kgot(MXV+1)
      dimension igot(MXV+1),t(MDIM)
      dimension aims31_local(IMSHPPH),aims32_local(IMSHPPH),
     :aims39_local(IMSHPPH),aims41_local(IMSHPPH),
     :aims42_local(IMSHPPH),aims43_local(IMSHPPH),
     :aims40_local(IMSHPPH)


      data zero,half,one,onep5,two/0.0,0.5,1.0,1.5,2.0/
      data four/4.0d0/

      call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)

      do 1010 i=0,MXV
      kgot(i)=0
      jgot(i)=0
      mgot(i)=0
      lgot(i)=0
 1010 igot(i)=0
        do i=1,IMSHPPH
         aims31_local(i)=0.0
        enddo

        do i=1,IMSHPPH
         aims32_local(i)=0.0
        enddo
        do i=1,IMSHPPH
         aims39_local(i)=0.0
        enddo
        do i=1,IMSHPPH
         aims41_local(i)=0.0
        enddo
        do i=1,IMSHPPH
         aims42_local(i)=0.0
        enddo
        do i=1,IMSHPPH
         aims43_local(i)=0.0
        enddo
        do i=1,IMSHPPH
         aims40_local(i)=0.0
        enddo


C     <QG|Vbar|BR> AND <QG|Vbar|RB>
C     <PH||HP> AND <PH||PH> Block
                                                                                
      do 1020 ib=1,nocc
      ibsym=isymc(ib)
      orbb=orbc(ib)
      iab=iqc(ib)
      iib=kc(ib)
                                                                                
      do 1030 iq=iproc+1,nexcit,nproc
      iqsym=isyme(iq)
      ibqsym=mtbl(ibsym,iqsym)
      orbq=orbe(iq)
      iaq=iqe(iq)
      iiq=ke(iq)
c------------------------------------------------------------------------
c         generate the ims arrays for each L value 0....MXV
c------------------------------------------------------------------------
      do 1040 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iig=kc(ig)
      iag=iqc(ig)
                                                                                
      do 1050 ir=1,nexcit
      irsym=isyme(ir)
      iir=ke(ir)
      iar=iqe(ir)
      orbr=orbe(ir)
      irgsym=mtbl(irsym,igsym)
      if(irgsym.ne.ibqsym)go to 1050
      call findk(orbq,orbb,orbr,orbg,iaq,iab,iar,iag,kgot,kmax)
      if(kmax.eq.0)go to 1050
      do 1060 kloop=1,kmax
      kk=kgot(kloop)
      ak=float(kk)
      rk=vint(iiq,iig,iib,iir,kk)
      f4=(-1)**(orbq+orbg+kk+1)
      vqgbr=rk*f4
      sum=zero
                                                                                
      imsloc=imslocHPPH(ib,ir,iq,ig,kk)
      aims31_local(imsloc) = vqgbr
      aims32_local(imsloc) = 0.0
                                                                                
c------------------------------------------------------------------
c                 Begin of inner loops
c------------------------------------------------------------------
      do is=1,nexcit
      issym=isyme(is)
      iis=ke(is)
      orbs=orbe(is)
      ias=iqe(is)
      if(iopt.ne.0)then
      if(orbs.eq.orbb)then
      if(issym.eq.ibsym)then
      in2=idpa(iis,iib)
      rk=vint(iiq,iig,iis,iir,kk)
      d4=(-1)**(orbg+orbq+kk+1)
      twoe=rk*d4
      aims31_local(imsloc)=aims31_local(imsloc)+twoe*t(in2)
      endif
      endif
      endif
      enddo
                                                                                
      do 1070 id=1,nocc
      idsym=isymc(id)
      orbd=orbc(id)
      iad=iqc(id)
      iid=kc(id)
c     Diagrams with T1
      if(iopt.ne.0)then
      if(orbq.eq.orbd)then
      if(iqsym.eq.idsym)then
      in1=idpa(iiq,iid)
      rk=vint(iig,iid,iir,iib,kk)
      d4=(-1)**(orbg+orbd+kk+1)
      twoe=rk*d4
      aims31_local(imsloc) = aims31_local(imsloc)-twoe*t(in1)
      endif
      endif
      endif
                                                                                
      do 1080 is=1,nexcit
      issym=isyme(is)
      iis=ke(is)
      orbs=orbe(is)
      ias=iqe(is)
                                                                                
c     Diagrams with T1 and T1**2
      if(iopt.ne.0)then
      if(orbs.eq.orbb)then
      if(issym.eq.ibsym)then
      in2=idpa(iis,iib)
      if(orbq.eq.orbd)then
      if(iqsym.eq.idsym)then
      in1=idpa(iiq,iid)
      rk=vint(iig,iid,iir,iis,kk)
      d4=(-1)**(orbg+orbd+kk+1)
      twoe=rk*d4
      aims31_local(imsloc) = aims31_local(imsloc)-twoe*t(in1)*t(in2)
      endif
      endif
      endif
      endif
                                                                                
      endif
                                                                                
      isdsym=mtbl(issym,idsym)
      if(isdsym.ne.irgsym)go to 1090
      if(isdsym.ne.ibqsym)go to 1090
c     part of the diagram 31
      nqb=idra(iiq,iib,kk)
      nsd=idra(iis,iid,kk)
      in=ntloc(nqb,nsd,ibqsym,kk)
      rk=vint(iir,iis,iig,iid,kk)
      d4=(-1)**(orbg+orbd+kk+1)
      vn=rk*d4
      f1=one/(two*ak+1.0)
      f2=(-1)**(orbs+orbd+kk)
      aims31_local(imsloc) = aims31_local(imsloc)-vn*t(in)*f1*f2*half
c     part of the diagram 33
      call findk(orbb,orbs,orbq,orbd,iab,ias,iaq,iad,jgot,jmax)
      if(jmax.eq.0)go to 2000
      do 2010 jloop=1,jmax
      ll=jgot(jloop)
      al=float(ll)
      nqd=idra(iiq,iid,ll)
      nsb=idra(iis,iib,ll)
      isbsym=mtbl(issym,ibsym)
      in=ntloc(nqd,nsb,isbsym,ll)
      f1=(-1)**(orbs+orbd+kk)
      call sixj(ak,orbs,orbd,al,orbq,orbb,s6j)
      aims31_local(imsloc) =aims31_local(imsloc)-vn*t(in)*f1*s6j*half
 2010 continue
 2000 continue
 1090 continue
c     part of the diagram 32
      irdsym=mtbl(irsym,idsym)
      isgsym=mtbl(issym,igsym)
      if(irdsym.ne.isgsym)go to 2020
      if(isdsym.ne.ibqsym)go to 2020
      call findk(orbd,orbr,orbs,orbg,iad,iar,ias,iag,jgot,jmax)
      if(jmax.eq.0)go to 2030
      do 2040 jloop=1,jmax
      ll=jgot(jloop)
      al=float(ll)
      nqb=idra(iiq,iib,kk)
      nsd=idra(iis,iid,kk)
      isdsym=mtbl(issym,idsym)
      in=ntloc(nqb,nsd,ibqsym,kk)
      rk=vint(iir,iis,iid,iig,ll)
      d4=(-1)**(orbd+orbg+ll+1)
      vn=rk*d4
      f1=(-1)**(orbs+orbd+kk)
      call sixj(ak,orbs,orbd,al,orbr,orbg,s6j)
      aims32_local(imsloc) =aims32_local(imsloc)-vn*t(in)*f1*s6j*half
 2040 continue
 2030 continue
 2020 continue
 1080 continue
 1070 continue
 1060 continue
 1050 continue
                                                                                
c     part of the diagram 41
      do 2050 ir=1,nexcit
      irsym=isyme(ir)
      iir=ke(ir)
      iar=iqe(ir)
      orbr=orbe(ir)
      irgsym=mtbl(irsym,igsym)
      if(irgsym.ne.ibqsym)go to 2050
      call findk(orbq,orbr,orbb,orbg,iaq,iar,iab,iag,kgot,kmax)
      if(kmax.eq.0)go to 2050
      do 2060 kloop=1,kmax
      kk=kgot(kloop)
      ak=float(kk)
      rk=vint(iiq,iig,iir,iib,kk)
      d4=(-1)**(orbq+orbg+kk+1)
      vqgrb=rk*d4
                                                                                
      imsloc=imslocHPPH(ib,iq,ir,ig,kk)
      aims41_local(imsloc)=vqgrb
      aims42_local(imsloc)=zero
      aims43_local(imsloc)=zero
                                                                                
c     This part aims42 will go into two places
c     (1) sum(r,g)<qr|T2|ga><pg|V|rb>
c     (2) sum(r,g)<qr|T2|bg><pg|V|rb>
                                                                                
      do is=1,nexcit
      issym=isyme(is)
      iis=ke(is)
      orbs=orbe(is)
      ias=iqe(is)
                                                                                
      if(iopt.ne.0)then
      if(orbs.eq.orbb)then
      if(issym.eq.ibsym)then
      in2=idpa(iis,iib)
      rk=vint(iiq,iig,iir,iis,kk)
      d4=(-1)**(orbg+orbq+kk+1)
      twoe=rk*d4
      aims41_local(imsloc)=aims41_local(imsloc)+twoe*t(in2)

      aims42_local(imsloc)=aims42_local(imsloc)+twoe*t(in2)
      aims43_local(imsloc)=aims43_local(imsloc)+twoe*t(in2)
      endif
      endif
      endif
      enddo
                                                                                
      do 2070 id=1,nocc
      idsym=isymc(id)
      orbd=orbc(id)
      iad=iqc(id)
      iid=kc(id)
                                                                                
c     Diagram with T1
      if(iopt.ne.0)then
      if(orbd.eq.orbq)then
      if(idsym.eq.iqsym)then
      in1=idpa(iiq,iid)
      rk=vint(iid,iig,iir,iib,kk)
      d4=(-1)**(orbg+orbd+kk+1)
      twoe=rk*d4
      aims41_local(imsloc)=aims41_local(imsloc)-twoe*t(in1)
      aims42_local(imsloc)=aims42_local(imsloc)-twoe*t(in1)
      aims43_local(imsloc)=aims43_local(imsloc)-twoe*t(in1)
      endif
      endif
      endif
                                                                                
      do 2080 is=1,nexcit
      issym=isyme(is)
      iis=ke(is)
      orbs=orbe(is)
      ias=iqe(is)
                                                                                
c     Diagram with T1 and T1**2
      if(iopt.ne.0)then
      if(orbs.eq.orbb)then
      if(issym.eq.ibsym)then
      in2=idpa(iis,iib)
      if(orbd.eq.orbq)then
      if(idsym.eq.iqsym)then
      in1=idpa(iiq,iid)
      rk=vint(iid,iig,iir,iis,kk)
      d4=(-1)**(orbg+orbd+kk+1)
      twoe=rk*d4
      aims41_local(imsloc)=aims41_local(imsloc)-twoe*t(in1)*t(in2)
      aims42_local(imsloc)=aims42_local(imsloc)-twoe*t(in1)*t(in2)
     : *half
      aims43_local(imsloc)=aims43_local(imsloc)-twoe*t(in1)*t(in2)
      endif
      endif
                                                                                
                                                                                
      endif
      endif
                                                                                
      endif
                                                                                
      isdsym=mtbl(issym,idsym)
      if(isdsym.ne.irgsym)go to 2090
      if(isdsym.ne.ibqsym)go to 2090
      call findk(orbd,orbr,orbs,orbg,iad,iar,ias,iag,jgot,jmax)
      if(jmax.eq.0)go to 2090
      do 3000 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      rk=vint(iir,iis,iid,iig,l2)
      d4=(-1)**(orbg+orbd+l2+1)
      vn=rk*d4
      call findk(orbd,orbq,orbs,orbb,iad,iaq,ias,iab,igot,imax)
      if(imax.eq.0)go to 3000
      do 3010 iloop=1,imax
      l4=igot(iloop)
      al4=float(l4)
      nsb=idra(iis,iib,l4)
      nqd=idra(iiq,iid,l4)
      isbsym=mtbl(issym,ibsym)
      in=ntloc(nsb,nqd,isbsym,l4)
      call sixj(al2,orbg,orbs,orbb,al4,ak,s6j1)
      call sixj(al2,orbr,orbd,orbq,al4,ak,s6j2)
      f1=(-1)**(kk+l2+l4+orbb+orbg+orbr+orbq)
      f2=s6j1*s6j2*(2*ak+1)
      aims41_local(imsloc)=aims41_local(imsloc)-vn*t(in)*f1*f2*half
 3010 continue
 3000 continue
 2090 continue
 2080 continue
 2070 continue
 2060 continue
 2050 continue
                                                                                
c     part of the diagram 39
      do 3050 ir=1,nexcit
      irsym=isyme(ir)
      iir=ke(ir)
      iar=iqe(ir)
      orbr=orbe(ir)
      irgsym=mtbl(irsym,igsym)
      ibqsym=mtbl(ibsym,iqsym)
      if(irgsym.ne.ibqsym)go to 3050
      call findk(orbq,orbb,orbr,orbg,iaq,iab,iar,iag,kgot,kmax)
      if(kmax.eq.0)go to 3050
      do 3060 kloop=1,kmax
      kk=kgot(kloop)
      ak=float(kk)
                                                                                
      imsloc=imslocHPPH(ib,ir,iq,ig,kk)
      aims39_local(imsloc)=0.0
                                                                                
      do 3070 id=1,nocc
      idsym=isymc(id)
      orbd=orbc(id)
      iad=iqc(id)
      iid=kc(id)
      do 3080 is=1,nexcit
      issym=isyme(is)
      iis=ke(is)
      orbs=orbe(is)
      ias=iqe(is)
      isdsym=mtbl(issym,idsym)
      if(isdsym.ne.irgsym)go to 3090
      if(isdsym.ne.ibqsym)go to 3090
      call findk(orbd,orbr,orbs,orbg,iad,iar,ias,iag,jgot,jmax)
      if(jmax.eq.0)go to 3090
      do 4000 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      rk=vint(iir,iis,iid,iig,l2)
      d4=(-1)**(orbg+orbd+l2+1)
      vn=rk*d4
      call findk(orbd,orbq,orbs,orbb,iad,iaq,ias,iab,igot,imax)
      if(imax.eq.0)go to 4000
      do 4010 iloop=1,imax
      l4=igot(iloop)
      al4=float(l4)
      nsb=idra(iis,iib,l4)
      nqd=idra(iiq,iid,l4)
      isbsym=mtbl(issym,ibsym)
      in=ntloc(nsb,nqd,isbsym,l4)
      call sixj(ak,orbq,orbb,al4,orbs,orbd,s6j1)
      call sixj(ak,orbg,orbr,al2,orbd,orbs,s6j2)
      f1=(-1)**(orbs+orbd+orbr+orbg)
      f2=s6j1*s6j2
      aims39_local(imsloc)=aims39_local(imsloc)+vn*t(in)*f1*f2
 4010 continue
 4000 continue
 3090 continue
 3080 continue
 3070 continue
 3060 continue
 3050 continue
 1040 continue
                                                                                
                                                                                
 1030 continue
 1020 continue
                                                                                
        call MPI_ALLREDUCE(aims39_local,aims39,IMSHPPH,MPI_REAL8,
     : MPI_SUM,MPI_COMM_WORLD,ierr)

        call MPI_ALLREDUCE(aims32_local,aims32,IMSHPPH,MPI_REAL8,
     : MPI_SUM,MPI_COMM_WORLD,ierr)

        call MPI_ALLREDUCE(aims31_local,aims31,IMSHPPH,MPI_REAL8,
     : MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(aims42_local,aims42,IMSHPPH,MPI_REAL8,
     : MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(aims43_local,aims43,IMSHPPH,MPI_REAL8,
     : MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(aims41_local,aims41,IMSHPPH,MPI_REAL8,
     : MPI_SUM,MPI_COMM_WORLD,ierr)

        do i=1,IMSHPPH
         aims31_local(i)=0.0
        enddo

        do i=1,IMSHPPH
         aims32_local(i)=0.0
        enddo
        do i=1,IMSHPPH
         aims39_local(i)=0.0
        enddo
        do i=1,IMSHPPH
         aims41_local(i)=0.0
        enddo
        do i=1,IMSHPPH
         aims42_local(i)=0.0
        enddo
        do i=1,IMSHPPH
         aims43_local(i)=0.0
        enddo

c     <AB|Vbar|GD> ==== <HH|HH> Block
                                                                                
      do 4020 ia=1,nocc
      iasym=isymc(ia)
      orba=orbc(ia)
      iaa=iqc(ia)
      iia=kc(ia)
                                                                                
      do 4030 ib=1,nocc
      ibsym=isymc(ib)
      orbb=orbc(ib)
      iab=iqc(ib)
      iib=kc(ib)
                                                                                
      do 4040 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      iig=kc(ig)
                                                                                
      do 4040 id=1,nocc
      idsym=isymc(id)
      orbd=orbc(id)
      iad=iqc(id)
      iid=kc(id)
                                                                                
      iabsym=mtbl(iasym,ibsym)
      igdsym=mtbl(igsym,idsym)
      if(iabsym.ne.igdsym)go to 4040
      call findk(orba,orbg,orbb,orbd,iaa,iag,iab,iad,kgot,kmax)
      if(kmax.eq.0)go to 4040
      do 4050 kloop=1,kmax
      kk=kgot(kloop)
      ak=float(kk)
      rk=vint(iia,iib,iig,iid,kk)
      d4=(-1)**(orbg+orbd+kk+1)
      vabgd=rk*d4
      imsloc=imslocHHHH(ia,ig,ib,id,kk)
      aims40(imsloc)=vabgd
                                                                                
c     part of the diagram 40
      do ip=1,nexcit
      ipsym=isyme(ip)
      iip=ke(ip)
      iap=iqe(ip)
      orbp=orbe(ip)
                                                                                
c     diagram with T1
      if(iopt.ne.0)then
      if(orbp.eq.orbb)then
      if(ipsym.eq.ibsym)then
      in=idpa(iip,iib)
      rk=vint(iig,iid,iia,iip,kk)
      d4=(-1)**(orbg+orbd+kk+1)
      twoe=rk*d4
      aims40(imsloc)=aims40(imsloc)+twoe*t(in)*two
      endif
      endif
      endif
      enddo
                                                                                
      do 4060 ip=1,nexcit
      ipsym=isyme(ip)
      iip=ke(ip)
      iap=iqe(ip)
      orbp=orbe(ip)
      do 4070 iq=1,nexcit
      iqsym=isyme(iq)
      iiq=ke(iq)
      iaq=iqe(iq)
      orbq=orbe(iq)
                                                                                
c     diagram with T1
      if(iopt.ne.0)then
      if(orbp.eq.orbb)then
      if(ipsym.eq.ibsym)then
      if(orbq.eq.orba)then
      if(iqsym.eq.iasym)then
      in1=idpa(iip,iib)
      in2=idpa(iiq,iia)
      rk=vint(iig,iid,iiq,iip,kk)
      d4=(-1)**(orbg+orbd+kk+1)
      twoe=rk*d4
      aims40(imsloc)=aims40(imsloc)+twoe*t(in1)*t(in2)
      endif
      endif
      endif
      endif
                                                                                
      endif
                                                                                
c     diagram with T2
      ipqsym=mtbl(ipsym,iqsym)
      ipasym=mtbl(ipsym,iasym)
      iqbsym=mtbl(iqsym,ibsym)
      if(iabsym.ne.ipqsym)go to 4070
      if(ipasym.ne.iqbsym)go to 4070
      call findk(orbp,orbg,orbq,orbd,iap,iag,iaq,iad,jgot,jmax)
      if(jmax.eq.0)go to 4070
      do 4080 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      rk=vint(iig,iid,iip,iiq,l2)
      d4=(-1)**(orbg+orbd+l2+1)
      vn=rk*d4
      call findk(orbp,orba,orbq,orbb,iap,iaa,iaq,iab,igot,imax)
      if(imax.eq.0)go to 4080
      do 4090 iloop=1,imax
      l4=igot(iloop)
      al4=float(l4)
      npa=idra(iip,iia,l4)
      nqb=idra(iiq,iib,l4)
      in=ntloc(npa,nqb,ipasym,l4)
      call sixj(al2,orbg,orbp,orba,al4,ak,s6j1)
      call sixj(al2,orbd,orbq,orbb,al4,ak,s6j2)
      f1=(-1)**(orba+orbb+orbg+orbd)
      f2=s6j1*s6j2*(2*ak+1)
      aims40(imsloc)=aims40(imsloc)+vn*t(in)*f1*f2
                                                                                
 4090 continue
 4080 continue
 4070 continue
 4060 continue
 4050 continue
 4040 continue
 4030 continue
 4020 continue
      
      return
      end      

      function slater_new(ii,jj,kk,ll,orbi,orbj,orbk,orbl,
     :iai,iaj,iak,ial,l1)
      implicit real*8 (a-h,o-z)
       include 'mpif.h'
      data zero,half,one,two/0.0d0,0.50d0,1.0d0,2.0d0/


      slater_new=0.0
      fact=0.0

      fact=select1(iai,iak,orbi,orbk,l1)*
     : select1(iaj,ial,orbj,orbl,l1)

      if(fact.ne.zero)then

      ak=dfloat(l1)
      f1=dr(orbi,ak,orbk,half,zero,-half)
      f2=dr(orbj,ak,orbl,half,zero,-half)
      f3=dsqrt((two*orbi+1)*(two*orbj+1)*(two*orbk+1)*(two*orbl+1))
c     f4=(-one)**(l1+orbi+orbj+1)
c       if(iall.eq.0)then
      slater_new=f1*f2*f3*slater(ii,jj,kk,ll,l1)
  

c     else

c      slater_new1=f1*f2*f3*slater(ii,jj,kk,ll,l1)
c      slater_new2=f1*f2*f3*slater_breitkm1(ii,jj,kk,
c     : ll,orbi,orbj,orbk,orbl,iai,iaj,iak,ial,l1)
c      slater_new3=f1*f2*f3*slater_breitkp1(ii,jj,kk,
c     : ll,orbi,orbj,orbk,orbl,iai,iaj,iak,ial,l1)

c       slater_new=slater_new1+slater_new2+slater_new3

c      endif

      endif

      return
      end

      integer function idelta(ja,jp,k)
      implicit none
      include 'mpif.h'
      real*8 ja,jp,eps
      integer k,kmax,kmin
      eps=1.0d-10
      kmax = int(ja + jp + eps)
      kmin = int(abs(ja - jp)+eps)
      if (k.le.kmax.and.k.ge.kmin) then
         idelta = 1   ! triangle condition satisfied
      else
         idelta = 0   ! not satisfied
      endif
      return
      end


      subroutine findk(orba,orbb,orbp,orbq,ia,ib,ip,iq,jeven,kkk)
      implicit real*8 (a-h,o-z)
      include 'mpif.h'
       common/iallparity/iall

      dimension jeven(MXV+1)

        if(iall.eq.0)then

       call findk2(orba,orbb,orbp,orbq,ia,ib,ip,iq,jeven,kkk)

       else

       call findk_all(orba,orbb,orbp,orbq,ia,ib,ip,iq,jeven,kkk)

       endif

      return
      end

      subroutine findk_all(orba,orbb,orbp,orbq,ia,ib,ip,iq,jeven,
     :kkk)
      implicit real*8 (a-h,o-z)
      include 'mpif.h'

      dimension jeven(MXV+1)

      eps=1.0d-10
      
      kkk=0
      do k=0, MXV
         if(idelta(orba,orbb,k).eq.1.and.
     &      idelta(orbp,orbq,k).eq.1) then
            kkk=kkk+1
            jeven(kkk)=k
         endif
      enddo

      return
      end

      subroutine fbar(t,iopt)
c----------------------------------------------------------------
c     THIS SUBROUTINE IS A PART OF INTERMEDIATE 
C     STORAGE. IT COMPUTES THE EFFECTIVE ONE-BODY
C     INTEGRALS USED IN CMAT2 AND CMAT2.
c----------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      include 'mpif.h'
       real*8 sum_local

      common/jvalue/orbc(MNOCC),orbe(MNEXC)
      common/index/ke(MNEXC),kc(MNOCC)
      common/info/nocc,nexcit,iType(MNBAS),iTCount(MNBAS)
      common/scratch2/h2d(MDIM),raux(MDIM),ap(MDIM)
      common/orbital_energy/eorb(MNBAS)
      common/parity/isymc(MNOCC),isyme(MNEXC)
      common/symmetry/mtbl(MNS,MNS)
      common/kpavlaue/iqc(MNOCC),iqe(MNEXC)
      common/idra1/idra(MNBAS,MNBAS,0:MXV)
      common/idpa1/idpa(MNBAS,MNBAS)
      common/skip/nnskip(0:MXV,2),n2skip,nsing
c
      common/ims1/f(MNBAS,MNBAS)
c     local dimension
      dimension jgot(MXV+1),lgot(MXV+1),mgot(MXV+1),kgot(MXV+1)
      dimension igot(MXV+1),t(MDIM)
     

      data zero,half,one,onep5,two/0.0,0.5,1.0,1.5,2.0/
      data four/4.0d0/
      call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)

      
      do 1 i=0,MXV
      kgot(i)=0
      jgot(i)=0
      mgot(i)=0
      lgot(i)=0
   1  igot(i)=0
      
      do i=1,MNBAS
      do j=1,MNBAS
      f(i,j)=0
      enddo
      enddo

      do 2210 ip=1,nexcit
      ipsym=isyme(ip)
      orbp=orbe(ip)
      iap=iqe(ip)
      iip=ke(ip)
                                                                                
      do 2220 ir=1,nexcit
      irsym=isyme(ir)
      iir=ke(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
                                                                                
      if(orbr.ne.orbp)go to 2220
      if(irsym.ne.ipsym)go to 2220
                                                                                
      sum_local=zero
                                                                                
      do 2230 ig=1,nocc
      igsym=isymc(ig)
      iig=kc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      irgsym=mtbl(irsym,igsym)
                                                                                
      do 2240 is=iproc+1,nexcit,nproc
      issym=isyme(is)
      iis=ke(is)
      ias=iqe(is)
      orbs=orbe(is)
                                                                                
c     Fbar with T1
      if(iopt.ne.0)then
                                                                                
      if(orbs.eq.orbg)then
      if(issym.eq.igsym)then
      in=idpa(iis,iig)
      rk=vint(iip,iig,iir,iis,0)
      d4=(-1)**(orbp+orbg+1)
      twoe=rk*d4
      fact=dsqrt((2*orbg+1)/(2*orbp+1))
      sum_local=sum_local+twoe*t(in)*fact
      call findk(orbp,orbs,orbg,orbr,iap,ias,iag,iar,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      kk=igot(iloop)
      rk=vint(iip,iig,iis,iir,kk)
      d4=(-1)**(orbp+orbg+kk+1)
      twoe=rk*d4
      fact=(-1)**(orbg+orbp+kk)/(2*orbp+1)
      sum_local=sum_local+twoe*t(in)*fact
      enddo
      endif
                                                                                
      endif
      endif
                                                                                
      endif
                                                                                
      do 2250 id=1,nocc
c--------------------------------------------------------------
c                   diagram 35 -<gd|v|rs><qr|t2|ba><ps|t2|gd>
c--------------------------------------------------------------
      idsym=isymc(id)
      orbd=orbc(id)
      iad=iqc(id)
      iid=kc(id)
                                                                                
c     Fbar with T1*T1
      if(iopt.ne.0)then
                                                                                
      if(orbs.eq.orbg)then
      if(issym.eq.igsym)then
      if(orbp.eq.orbd)then
      if(ipsym.eq.idsym)then
      in1=idpa(iis,iig)
      in2=idpa(iip,iid)
      rk=vint(iid,iig,iir,iis,0)
      d4=(-1)**(orbd+orbg+1)
      twoe=rk*d4
      fact=dsqrt((2*orbg+1)/(2*orbp+1))
      sum_local=sum_local-twoe*t(in1)*t(in2)*fact*half
      call findk(orbd,orbs,orbg,orbr,iad,ias,iag,iar,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      kk=igot(iloop)
      rk=vint(iid,iig,iis,iir,kk)
      d4=(-1)**(orbd+orbg+kk+1)
      twoe=rk*d4
      fact=(-1)**(orbg+orbp+kk)/(2*orbp+1)
      sum_local=sum_local-twoe*t(in1)*t(in2)*fact*half
      enddo
      endif
                                                                                
      endif
      endif
      endif
      endif
                                                                                
      endif
                                                                                
c     Fbar with T2
      isdsym=mtbl(issym,idsym)
      ipgsym=mtbl(ipsym,igsym)
      if(irgsym.ne.isdsym)go to 2250
      call findk(orbg,orbr,orbd,orbs,iag,iar,iad,ias,igot,imax)
      if(imax.eq.0)go to 2260
      do 2270 kloop=1,imax
      l2=igot(kloop)
      al2=dfloat(l2)
      npg=idra(iip,iig,l2)
      nsd=idra(iis,iid,l2)
      in=ntloc(npg,nsd,ipgsym,l2)
      rk=vint(iir,iis,iig,iid,l2)
      f1=one/(two*orbp+1)*(one/(two*l2+1))
      f2=(-one)**(orbg+orbd+orbr+orbs)
      d4=(-1)**(orbg+orbd+l2+1)
      vn=rk*d4*f1*f2
      sum_local=sum_local - vn*t(in)
c-------------------------------------------------------------------
c                   diagram 37 <gd|v|sr><ps|t|gd><qr|t2|ba>
c-------------------------------------------------------------------
      call findk(orbg,orbs,orbd,orbp,iag,ias,iad,iap,jgot,jmax)
      if(jmax.eq.0)go to 2270
      do 2280 lloop=1,jmax
      l4=jgot(lloop)
      al4=dfloat(l4)
      npd=idra(iip,iid,l4)
      nsg=idra(iis,iig,l4)
      isgsym=mtbl(issym,igsym)
      in=ntloc(npd,nsg,isgsym,l4)
      call sixj(al2,orbs,orbd,al4,orbp,orbg,s6j)
      f1=one/(two*orbp+1)*s6j
      f2=(-one)**(orbg+orbd+orbr+orbs)
      vn=rk*f1*f2*d4
      sum_local=sum_local- vn*t(in)
 2280 continue
 2270 continue
 2260 continue
 2250 continue
 2240 continue
 2230 continue
       call MPI_ALLREDUCE(sum_local,sum,1,MPI_REAL8,MPI_SUM,
     :                   MPI_COMM_WORLD,ierr)
       sum_local=0

      f(iip,iir)=sum
 2220 continue
 2210 continue
c     print*,'fbar1 = f',f
                                                                                
c-------------------------------------------------------------
c           diagram 36 -<gd|v|rs><pq|t2|db><rs|t2|ag>
c-------------------------------------------------------------
      do 2290 ia=1,nocc
      iasym=isymc(ia)
      iia=kc(ia)
      orba=orbc(ia)
      iaa=iqc(ia)
                                                                                
      do 2300 id=1,nocc
      idsym=isymc(id)
      iid=kc(id)
      orbd=orbc(id)
      iad=iqc(id)
                                                                                
      if(orba.ne.orbd)go to 2300
      if(iasym.ne.idsym)go to 2300
                                                                                
      sum_local=zero
                                                                                
      do 2310 ig=1,nocc
      igsym=isymc(ig)
      iig=kc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
                                                                                
      do 2320 is=iproc+1,nexcit,nproc
      issym=isyme(is)
      iis=ke(is)
      ias=iqe(is)
      orbs=orbe(is)
                                                                                
c     Fbar with T1
      if(iopt.ne.0)then
                                                                                
      if(orbs.eq.orbg)then
      if(issym.eq.igsym)then
      in=idpa(iis,iig)
      rk=vint(iig,iid,iis,iia,0)
      d4=(-1)**(orbg+orbd+1)
      twoe=rk*d4
      fact=dsqrt((2*orbg+1)/(2*orba+1))
      sum_local=sum_local+twoe*t(in)*fact
      call findk(orba,orbg,orbs,orbd,iaa,iag,ias,iad,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      kk=igot(iloop)
      rk=vint(iig,iid,iia,iis,kk)
      d4=(-1)**(orbg+orbd+kk+1)
      twoe=rk*d4
      fact=(-1)**(orbs+orba+kk)/(2*orba+1)
      sum_local=sum_local+twoe*t(in)*fact
      enddo
      endif
                                                                                
      endif
      endif
                                                                                
      endif
                                                                                
      do 2330 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      iir=ke(ir)
                                                                                
c     Fbar with T1
      if(iopt.ne.0)then
                                                                                
      if(orbs.eq.orbg)then
      if(issym.eq.igsym)then
      if(orbr.eq.orba)then
      if(irsym.eq.iasym)then
      in1=idpa(iis,iig)
      in2=idpa(iir,iia)
      rk=vint(iig,iid,iis,iir,0)
      d4=(-1)**(orbg+orbd+1)
      twoe=rk*d4
      fact=dsqrt((2*orbg+1)/(2*orba+1))
      sum_local=sum_local+twoe*t(in1)*t(in2)*fact*half
      call findk(orbr,orbg,orbs,orbd,iar,iag,ias,iad,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      kk=igot(iloop)
      rk=vint(iig,iid,iir,iis,kk)
      d4=(-1)**(orbg+orbd+kk+1)
      twoe=rk*d4
      fact=(-1)**(orbs+orba+kk)/(2*orba+1)
      sum_local=sum_local+twoe*t(in1)*t(in2)*fact*half
      enddo
      endif
                                                                                
      endif
      endif
      endif
      endif
                                                                                
      endif
                                                                                
c     Fbar with T2
      isgsym=mtbl(issym,igsym)
      irdsym=mtbl(irsym,idsym)
      irgsym=mtbl(irsym,igsym)
      if(irdsym.ne.isgsym)go to 2330
      call findk(orbg,orbs,orbd,orbr,iag,ias,iad,iar,igot,imax)
      if(imax.eq.0)go to 2330
      do 2340 kloop=1,imax
      l2=igot(kloop)
      al2=dfloat(l2)
      nra=idra(iir,iia,l2)
      nsg=idra(iis,iig,l2)
      isgsym=mtbl(issym,igsym)
      in=ntloc(nra,nsg,isgsym,l2)
      d4=(-1)**(orbg+orbd+l2+1)
      rk=vint(iis,iir,iig,iid,l2)
      f1=one/(two*orba+1)*(one/(two*l2+1))
      f2=(-one)**(orbg+orbd+orbr+orbs)
      vn=rk*d4*f1*f2
      sum_local=sum_local+ vn*t(in)
c---------------------------------------------------------------------
c                   diagram 38 <gd|v|rs><rs|t2|ag><pq|t2|db>
c---------------------------------------------------------------------
      call findk(orbg,orbr,orba,orbs,iag,iar,iaa,ias,jgot,jmax)
      if(jmax.eq.0)go to 2340
      do 2350 lloop=1,jmax
      l4=jgot(lloop)
      al4=dfloat(l4)
      nrg=idra(iir,iig,l4)
      nsa=idra(iis,iia,l4)
      isasym=mtbl(issym,iasym)
      in=ntloc(nrg,nsa,isasym,l4)
      call sixj(al2,orbs,orbg,al4,orbr,orba,s6j)
      f1=s6j/(two*orba+1)
      f2=(-1)**(orbg+orbr+orbd+orbs)
      vn=rk*d4*f1*f2
      sum_local=sum_local+vn*t(in)
 2350 continue
 2340 continue
 2330 continue
 2320 continue
 2310 continue
       call MPI_ALLREDUCE(sum_local,sum,1,MPI_REAL8,MPI_SUM,
     :                   MPI_COMM_WORLD,ierr)
       sum_local=0

      f(iid,iia)=sum
 2300 continue
 2290 continue
c     print*,'fbar2 = f',f
                                                                                
c     Fbar (a,p)
      if(iopt.ne.0)then
                                                                                
      do 2490 ia=1,nocc
      iasym=isymc(ia)
      iia=kc(ia)
      orba=orbc(ia)
      iaa=iqc(ia)
      do 2500 ip=1,nexcit
      ipsym=isyme(ip)
      orbp=orbe(ip)
      iap=iqe(ip)
      iip=ke(ip)
                                                                                
      if(orbp.ne.orba)go to 2500
      if(ipsym.ne.iasym)go to 2500
      sum_local=zero
                                                                                
      do 2510 ig=1,nocc
      igsym=isymc(ig)
      iig=kc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
                                                                                
      do 2520 is=iproc+1,nexcit,nproc
      issym=isyme(is)
      iis=ke(is)
      ias=iqe(is)
      orbs=orbe(is)
                                                                                
      if(orbs.ne.orbg)go to 2520
      if(issym.ne.igsym)go to 2520
      in=idpa(iis,iig)
      rk=vint(iig,iia,iis,iip,0)
      d4=(-1)**(orbg+orba+1)
      twoe=rk*d4
      fact=dsqrt((2*orbg+1)/(2*orba+1))
      sum_local=sum_local+twoe*t(in)*fact
      call findk(orba,orbs,orbg,orbp,iaa,ias,iag,iap,igot,imax)
      if(imax.eq.0)go to 2520
      do 2530 iloop=1,imax
      kk=igot(iloop)
      rk=vint(iig,iia,iip,iis,kk)
      d4=(-1)**(orbg+orba+kk+1)
      twoe=rk*d4
      fact=(-1)**(orbg+orba+kk)/(2*orba+1)
      sum_local=sum_local+twoe*t(in)*fact
 2530 continue
 2520 continue
 2510 continue
       call MPI_ALLREDUCE(sum_local,sum,1,MPI_REAL8,MPI_SUM,
     :                   MPI_COMM_WORLD,ierr)
       sum_local=0

      f(iia,iip)=sum
 2500 continue
 2490 continue
      endif

c     print*,'fbar3 = f',f
      return
      end

      subroutine jacobi2(told,t,iopt,isolv,ntmax,nonlin,
     :          irows,irowe,iproc,nproc,dspl,delrow)
      implicit real*8 (a-h,o-z)
      include 'mpif.h'
      common/scratch/b(MDIM),ad(MDIM),baux(MDIM)
      common/scratch2/h2d(MDIM),raux(MDIM),ap(MDIM)                 
C
C     THIS SUBROITINE SOLVES THE BT=A FOR THE
C     COUPLED CLUSTER EQUATION ON FLY USING JACOBI
C     SCHEME
C
c     local dimension
      integer dspl(NPMAX),delrow(NPMAX),drow,error
      dimension told(MDIM),t(MDIM),caux(MDIM)
      start_time = 0
      end_time = 0

c     rewind(NTFILE)
c     read(NTFILE)ntmax,(told(i),i=1,ntmax)

c      start_time = MPI_WTIME()

      if (iproc .eq. 0) write(STDOUT,1000)

c--------------------------------------------------------------------
c     Start linear iteration
c--------------------------------------------------------------------

      nit=0
   30 nit=nit+1

      call fbar(told,iopt)
      call vbar(told,iopt)

      if(iopt.ne.0)then
      call cmat1(caux,told,ntmax,isolv,irows,irowe)
      endif

      call cmat2(caux,told,ntmax,iopt,isolv,irows,irowe)

c     check convergence

      do i=irows,irowe
      t(i)=(b(i)-caux(i))/ad(i)
      enddo

      eps=0

c      call MPI_ALLGATHERV(t(irows),drow,MPI_DOUBLE_PRECISION,t,
c     >delrow,dspl,MPI_DOUBLE_PRECISION, MPI_COMM_WORLD,error)

      do i=1,ntmax
      eps=eps+dabs(told(i)-t(i))
      told(i)=t(i)
      enddo

      conv=eps/ntmax

      call cizec(t,iopt,iproc)

      if(iproc.eq.0)then
      write(STDOUT,9997)nit,conv
      endif

      if(conv.le.SMALL)go to 40
      if(nit.ge.NITMAX)go to 50
      go to 30

c-----------------------------------------------------------------------
c     End linear iteration
c-----------------------------------------------------------------------

 50   if (iproc .eq. 0) write(STDOUT,9999)
      go to 60

 40   if (iproc .eq. 0) write(STDOUT,9998)nit
 60   continue

c      end_time= MPI_WTIME()

      if(iproc.eq.0)then
      open(NTFILE,file='to.dat',form='unformatted',
     &     status='unknown')

      rewind(NTFILE)
      write(NTFILE)ntmax,(t(i),i=1,ntmax)
      endif

      if (iproc .eq. 0) write(STDOUT,1002)end_time-start_time
   42 continue

      call cizec(t,iopt,iproc)

 998  format(/,1x,'diagonal part of b-matrix',//)
 999  format(4f20.11)
1000  format(/,1x,'Jacobi iteration begins',//)
1001  format(/,1x,'solution vectors',//)
1002  format(/,1x,'Time elapsed=',2f20.11//)
9990  format(//,1x,'End of linear CC',//)
9991  format(//,1x,'Non-linear Iteration begins',//)
9996  format(//,1x,'Iteration begins',//)
9997  format(1x,'Iteration no=',i4,5x,'eps=',f20.11)
9998  format(1x,' Iteration converges at',I4,' cycle')
9999  format(1x,' Convergence failure; Programme terminates')

      return
      end


      subroutine cmat2(raux,t,ntmax,iopt,isolv,rowstart,rowend)
      implicit real*8 (a-h,o-z)
c
c     THIS SUBROUTINE IS SOMEWHAT IDENTICAL
C     TO SUBROUTINE DOUBLES USED FOR LINEARIZED CC EQUATIONS.
C     HERE IT CONTAINS THE NON-LINEAR PART TROUGH EFFECTIVE
C     ONE (FB) AND TWO-BODY (VIA AIMS31,AIMS32,AIMS39,AIMS40,
C     AIMS41,AIMS42 IN ADDITION TO SOME BARE TWO ELECTRON
C     INTEGRALS)
C
C
C     SEE DOC. FOR DETAILS. SOME DETAILS OF CMAT1 IS AVAILABLE
C     ELSE WHERE

      include 'mpif.h'
      integer row,rowstart,rowend
c     common block
      common/jvalue/orbc(MNOCC),orbe(MNEXC)
      common/index/ke(MNEXC),kc(MNOCC)
      common/scratch/b(MDIM),ad(MDIM),baux(MDIM)
      common/info/nocc,nexcit,iType(MNBAS),iTCount(MNBAS)
      common/orbital_energy/eorb(MNBAS)
      common/parity/isymc(MNOCC),isyme(MNEXC)
      common/symmetry/mtbl(MNS,MNS)
      common/kpavlaue/iqc(MNOCC),iqe(MNEXC)
      common/idra1/idra(MNBAS,MNBAS,0:MXV)
      common/idpa1/idpa(MNBAS,MNBAS)
      common/skip/nnskip(0:MXV,2),n2skip,nsing
c
      common/ims1/f(MNBAS,MNBAS)
c
      common/vb/aims31(IMSHPPH),
     :aims32(IMSHPPH),
     :aims39(IMSHPPH),
     :aims40(IMSHHHH),
     :aims41(IMSHPPH),
     :aims42(IMSHPPH),
     :aims43(IMSHPPH)
c
c     local dimension
      dimension t(MDIM),raux(MDIM)
      dimension kgot(MXV+1),jgot(MXV+1),igot(MXV+1)

      data zero,half,one,onep5,two/0.0,0.5,1.0,1.5,2.0/
      data four/4.0d0/

      do i=0,MXV
      kgot(i)=0
      jgot(i)=0
      igot(i)=0
      enddo

      do 10 ia=1,nocc
      iasym=isymc(ia)
      do 20 ip=1,nexcit
      ipsym=isyme(ip)
      iapsym=mtbl(ipsym,iasym)
      do 30 ib=1,nocc
      ibsym=isymc(ib)
      do 40 iq=1,nexcit
      iqsym=isyme(iq)
      ibqsym=mtbl(ibsym,iqsym)
      if(ibqsym.ne.iapsym)go to 40
      orba=orbc(ia)
      orbb=orbc(ib)
      orbp=orbe(ip)
      orbq=orbe(iq)
      iaa=iqc(ia)
      iab=iqc(ib)
      iap=iqe(ip)
      iaq=iqe(iq)
      iia=kc(ia)
      iib=kc(ib)
      iip=ke(ip)
      iiq=ke(iq)
      call findk(orba,orbp,orbb,orbq,iaa,iap,iab,iaq,kgot,
     $kmax)
      if(kmax.eq.0)go to 40
      do 45 jloop=1,kmax
      l1=kgot(jloop)
      al1=dfloat(l1)
      call iloc(iip,iia,iiq,iib,ipa,iqb)
      if(ipa.lt.iqb)go to 45
      nap=idra(iip,iia,l1)
      nbq=idra(iiq,iib,l1)
      neqn=ntloc(nap,nbq,iapsym,l1)
      row=neqn
      irow=neqn+1-rowstart
      if(irow.le.0)go to 45
      if(irow.gt.rowend-rowstart+1)go to 45
      call init(baux)
      mult=2
      mloop=1
      if(iip.eq.iiq.and.iia.eq.iib)go to 50
      mult=1
      mloop=2
  50  do 60 nloop=1,mloop
      go to (70,80),nloop
  70  jp=iip
      ja=iia
      jq=iiq
      jb=iib
      kp=ip
      kq=iq
      ka=ia
      kb=ib
      jasym=iasym
      jbsym=ibsym
      jpsym=ipsym
      jqsym=iqsym
      jaa=iaa
      jab=iab
      jap=iap
      jaq=iaq
      orbap=orba
      orbbp=orbb
      orbpp=orbp
      orbqp=orbq
      go to 90
   80 jp=iiq
      ja=iib
      jq=iip
      jb=iia
      kp=iq
      kq=ip
      ka=ib
      kb=ia
      jasym=ibsym
      jbsym=iasym
      jpsym=iqsym
      jqsym=ipsym
      jaa=iab
      jab=iaa
      jap=iaq
      jaq=iap
      orbap=orbb
      orbbp=orba
      orbpp=orbq
      orbqp=orbp
   90 continue
                                                                                
c     ph-hp diagrams
      do 800 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      iir=ke(ir)
                                                                                
      do 900 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      iig=kc(ig)
      irgsym=mtbl(irsym,igsym)
      iqgsym=mtbl(jqsym,igsym)
      jqbsym=mtbl(jqsym,jbsym)
                                                                                
c
      if(iopt.ne.0)then
                                                                                
      if(orbg.eq.orbr)then
      if(igsym.eq.irsym)then
                                                                                
      if(orbqp.eq.orbg)then
      if(jqsym.eq.igsym)then
      in1=idpa(jq,iig)
      npa=idra(jp,ja,l1)
      nrb=idra(iir,jb,l1)
      irbsym=mtbl(irsym,jbsym)
      in=ntloc(npa,nrb,irbsym,l1)
      baux(in)=baux(in)-f(iig,iir)*t(in1)*mult*half
      endif
      endif
                                                                                
      if(orbbp.eq.orbr)then
      if(jbsym.eq.irsym)then
      in1=idpa(iir,jb)
      npa=idra(jp,ja,l1)
      nqg=idra(jq,iig,l1)
      igqsym=mtbl(igsym,jqsym)
      in=ntloc(npa,nqg,igqsym,l1)
      baux(in)=baux(in)-f(iig,iir)*t(in1)*mult*half
      endif
      endif
                                                                                
      endif
      endif
c     T1*T1
                                                                                
      if(orbg.eq.orbqp)then
      if(igsym.eq.jqsym)then
      if(orbr.eq.orbap)then
      if(irsym.eq.jasym)then
      in1=idpa(jq,iig)
      in2=idpa(iir,ja)
      imsloc=imslocHPPH(kb,kp,ir,ig,l1)
      vn=aims42(imsloc)
      rk=vint(jp,iig,iir,jb,l1)
      f4=(-1)**(l1+orbpp+orbg+1)
      twoe=rk*f4
      baux(in2)=baux(in2)-(twoe+vn)*t(in1)*mult
      endif
      endif
      endif
      endif
                                                                                
      if(orbg.eq.orbqp)then
      if(igsym.eq.jqsym)then
      if(orbr.eq.orbbp)then
      if(irsym.eq.jbsym)then
      in1=idpa(jq,iig)
      in2=idpa(iir,jb)
      rk=vint(jp,iig,ja,iir,l1)
      f4=(-1)**(l1+orbpp+orbg+1)
      twoe=rk*f4
      baux(in2)=baux(in2)-twoe*t(in1)*mult
      endif
      endif
      endif
      endif
                                                                                
      endif
                                                                                
                                                                                
      if(irgsym.ne.jqbsym)go to 910
c
      imsloc=imslocHPPH(ka,ir,kp,ig,l1)
      vn1=aims31(imsloc)
      vn2=aims32(imsloc)
      vn3=aims39(imsloc)
c
c     diagram 1 -sum(r,g) <pg|v|ar><qr|t(2)|bg>
                                                                                
      ngr=idra(iir,iig,l1)
      nbq=idra(jq,jb,l1)
      in=ntloc(ngr,nbq,jqbsym,l1)
      fact=(-1)**(orbr+orbg+al1)
                                                                                
      baux(in)=baux(in)-vn1*fact*mult/(two*l1+1)
      baux(in)=baux(in)-vn2*fact*mult/(two*l1+1)
      baux(in)=baux(in)+vn3*mult
c     diagram 2 -sum(r,g) <pg|v|ar><qr|t(2)|gb>
      call findk(orbg,orbqp,orbbp,orbr,iag,jaq,jab,iar,igot,imax)
      if(imax.eq.0)go to 920
      do 930 mloop=1,imax
      l3=igot(mloop)
      al3=dfloat(l3)
      nbr=idra(iir,jb,l3)
      ngq=idra(jq,iig,l3)
      in=ntloc(nbr,ngq,iqgsym,l3)
      call sixj(al1,orbr,orbg,al3,orbqp,orbbp,s6j)
      baux(in)=baux(in)-vn1*fact*s6j*mult
  930 continue
  920 continue
c     ph-ph  diagram 3= (-)sum(r,g) <pg|v|ra><qr|t2|bg>
      nrg=idra(iir,iig,l1)
      nqb=idra(jq,jb,l1)
      in=ntloc(nrg,nqb,irgsym,l1)
c     find allowed k for <pg|v|ra>
      call findk(orbg,orbap,orbpp,orbr,iag,jaa,jap,iar,jgot,jmax)
      if (jmax.eq.0) go to 940
      do 950 kloop=1,jmax
      l2=jgot(kloop)
      al2=dfloat(l2)
                                                                                
      imsloc=imslocHPPH(ka,kp,ir,ig,l2)
      vn=aims43(imsloc)
                                                                                
      rk=vint(jp,iig,iir,ja,l2)
      f4=(-1)**(l2+orbg+orbpp+1)
      twoe=rk*f4
      fact=(-one)**(orbr+orbg+l1)
      call sixj(al1,orbg,orbr,al2,orbpp,orbap,s6j)
      baux(in)=baux(in)-twoe*fact*mult*s6j
      baux(in)=baux(in)-vn*fact*mult*s6j
  950 continue
  940 continue
  910 continue
      irasym=mtbl(irsym,jasym)
      iqgsym=mtbl(jqsym,igsym)
      if(irasym.ne.iqgsym)go to 960
c     find allowed k for t(2) amplitudes
      call findk(orbg,orbbp,orbpp,orbr,iag,jab,jap,iar,jgot,jmax)
      if (jmax.eq.0) go to 970
      do 980 kloop=1,jmax
      l2=jgot(kloop)
      al2=dfloat(l2)
      f4=(-1)**(l2+orbg+orbpp+1)
      rk=vint(jp,iig,iir,jb,l2)
      twoe=rk*f4
      imsloc=imslocHPPH(kb,kp,ir,ig,l2)
      vn=aims41(imsloc)
      call findk(orbg,orbqp,orbap,orbr,iag,jaq,jaa,iar,igot,imax)
      if(imax.eq.0)go to 980
      do 990 mloop=1,imax
      l3=igot(mloop)
      al3=dfloat(l3)
      nra=idra(iir,ja,l3)
      nqg=idra(jq,iig,l3)
      in=ntloc(nra,nqg,irasym,l3)
      fact=(two*al1+1)
      call sixj(al1,al2,al3,orbg,orbqp,orbbp,s6j1)
      call sixj(al1,al3,al2,orbr,orbpp,orbap,s6j2)
      f5=(-one)**(orba+orbp+orbb+orbq+l1+l2+l3)
      baux(in)=baux(in)-vn*mult*f5*s6j1*s6j2*fact
  990 continue
  980 continue
  970 continue
  960 continue
                                                                                
  900 continue
                                                                                
                                                                                
c     diagram 4 =sum(r,s) <pq|v|rs><rs|t2|ab>
      do 1000 is=1,nexcit
      issym=isyme(is)
      orbs=orbe(is)
      ias=iqe(is)
      iis=ke(is)
                                                                                
c     T1*T1
      if(iopt.ne.0)then
                                                                                
      if(orbbp.eq.orbs)then
      if(jbsym.eq.issym)then
      if(orbr.eq.orbap)then
      if(irsym.eq.jasym)then
      in1=idpa(iis,jb)
      in2=idpa(iir,ja)
      rk=vint(jp,jq,iir,iis,l1)
      f4=(-1)**(l1+orbpp+orbqp+1)
      twoe=rk*f4
      baux(in2)=baux(in2)+twoe*t(in1)*mult*half
      endif
      endif
      endif
      endif
                                                                                
      endif
                                                                                
      isbsym=mtbl(issym,jbsym)
      if(isbsym.ne.irasym)go to 1000
c     find allowed k for <pq|v|rs> for two electron integrals
      call findk(orbr,orbpp,orbs,orbqp,iar,jap,ias,jaq,jgot,jmax)
      if (jmax.eq.0) go to 1000
      do 1100 kloop=1,jmax
      l2=jgot(kloop)
      al2=dfloat(l2)
      rk=vint(jp,jq,iir,iis,l2)
      f4=(-1)**(l2+orbpp+orbqp+1)
      twoe=rk*f4
                                                                                
c     include T1 and T1*T1
      if(iopt.ne.0)then
                                                                                
      do ig=1,nocc
      iig=kc(ig)
      iag=iqc(ig)
      orbg=orbc(ig)
      igsym=isymc(ig)
      if(orbg.eq.orbqp)then
      if(igsym.eq.jqsym)then
      in1=idpa(jq,iig)
      rk=vint(jp,iig,iir,iis,l2)
      f4=(-1)**(orbpp+orbg+l2+1)
      twoe=twoe-rk*f4*t(in1)*two
      do id=1,nocc
      iid=kc(id)
      iad=iqc(id)
      orbd=orbc(id)
      idsym=isymc(id)
      if(orbd.eq.orbpp)then
      if(idsym.eq.jpsym)then
      in2=idpa(jp,iid)
      rk=vint(iig,iid,iis,iir,l2)
      f4=(-1)**(orbg+orbd+l2+1)
      twoe=twoe+rk*f4*t(in1)*t(in2)
      endif
      endif
      enddo
      endif
      endif
      enddo
                                                                                
      endif
                                                                                
      call findk(orbr,orbap,orbs,orbbp,iar,jaa,ias,jab,igot,imax)
      if (imax.eq.0)go to 1100
      do 1200 mloop=1,imax
      l3=igot(mloop)
      al3=dfloat(l3)
      f5=(-1)**(orba+orbb+orbp+orbq)
      nra=idra(iir,ja,l3)
      nsb=idra(iis,jb,l3)
      in=ntloc(nra,nsb,irasym,l3)
      fact=two*al1+1
      call sixj(al1,al3,al2,orbr,orbpp,orbap,s6j1)
      call sixj(al1,al3,al2,orbs,orbqp,orbbp,s6j2)
      baux(in)=baux(in)+twoe*half*mult*fact*s6j1*s6j2*f5
 1200 continue
 1100 continue
 1000 continue
                                                                                
      if(iopt.eq.0)go to 1300
c     diagram 8 <pq|v|ar><r|t1|b>
      iir=ke(ir)
      iar=iqe(ir)
      orbr=orbe(ir)
      if(irsym.ne.jbsym)go to 1300
      if(orbbp.ne.orbr)go to 1300
      in=idpa(iir,jb)
      rk=vint(jp,jq,ja,iir,l1)
      f4=(-1)**(l1+orbpp+orbqp+1)
      twoe=rk*f4
      baux(in)=baux(in)+twoe*mult
 1300 continue
      if(irsym.eq.jqsym)then
      if(orbr.eq.orbqp)then
      npa=idra(jp,ja,l1)
      nrb=idra(iir,jb,l1)
      ipasym=mtbl(jpsym,jasym)
      in=ntloc(npa,nrb,ipasym,l1)
      baux(in)=baux(in)+f(jq,iir)*mult
      endif
      endif
  800 continue
                                                                                
c     diagram 5 = sum(g,d) <gd|v|ab><pq|t2|gd>
      do 1400 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      iig=kc(ig)
      ipgsym=mtbl(jpsym,igsym)
                                                                                
      do 1500 id=1,nocc
      idsym=isymc(id)
      orbd=orbc(id)
      iad=iqc(id)
      iid=kc(id)
                                                                                
      if(iopt.ne.0)then
                                                                                
      if(orbpp.eq.orbg)then
      if(jpsym.eq.igsym)then
      if(orbqp.eq.orbd)then
      if(jqsym.eq.idsym)then
      in1=idpa(jq,iid)
      in2=idpa(jp,iig)
      rk=vint(iig,iid,ja,jb,l1)
      d4=(-1)**(orbg+orbd+l1+1)
      twoe=rk*d4
      baux(in2)=baux(in2)+twoe*t(in1)*mult*half
      endif
      endif
      endif
      endif
                                                                                
      endif
                                                                                
      iqdsym=mtbl(jqsym,idsym)
      if(iqdsym.ne.ipgsym)go to 1500
c     find allowd k for <gd|v|ab>
      call findk(orbd,orbbp,orbg,orbap,iad,jab,iag,jaa,jgot,jmax)
      if (jmax.eq.0)go to 1500
      do 1600 kloop=1,jmax
      l2=jgot(kloop)
      al2=dfloat(l2)
      rk=vint(iig,iid,ja,jb,l2)
      f4=(-1)**(l2+orbg+orbd+1)
      twoe=rk*f4
c
      imsloc=imslocHHHH(ka,ig,kb,id,l2)
      vn=aims40(imsloc)
c
c     find allowed k for <pq|t(2)|gd>
      call findk(orbd,orbqp,orbg,orbpp,iad,jaq,iag,jap,igot,imax)
      if(imax.eq.0) go to 1600
      do 1700 mloop=1,imax
      l3=igot(mloop)
      al3=dfloat(l3)
      npg=idra(jp,iig,l3)
      nqd=idra(jq,iid,l3)
      in=ntloc(npg,nqd,ipgsym,l3)
      fact=(two*al1+1)
      call sixj(al1,al2,al3,orbg,orbpp,orbap,s6j1)
      call sixj(al1,al2,al3,orbd,orbqp,orbbp,s6j2)
      f5=(-1)**(orba+orbp+orbb+orbq)
      baux(in)=baux(in)+vn*mult*fact*s6j1*s6j2*f5*half
 1700 continue
 1600 continue
 1500 continue
c     diagram 7 -sum(g)<pg|v|ab><q|t1|g>
      if(iopt.eq.0)go to 1800
      iig=kc(ig)
      iag=iqc(ig)
      orbg=orbc(ig)
      if(igsym.ne.jqsym)go to 1800
      if(orbg.ne.orbqp)go to 1800
      in=idpa(jq,iig)
      rk=vint(jp,iig,ja,jb,l1)
      f4=(-1)**(l1+orbpp+orbg+1)
      twoe=rk*f4
      baux(in)=baux(in)-twoe*mult
 1800 continue
      if(igsym.eq.jbsym)then
      if(orbg.eq.orbbp)then
      npa=idra(jp,ja,l1)
      nqg=idra(jq,iig,l1)
      ipasym=mtbl(jpsym,jasym)
      in=ntloc(npa,nqg,ipasym,l1)
      baux(in)=baux(in)-f(iig,jb)*mult
      endif
      endif
 1400 continue
                                                                                
                                                                                
  60  continue
      xsum=0
      do in=1,ntmax
      xsum=xsum+baux(in)*t(in)
      enddo
      raux(neqn)=xsum
  45  continue
  40  continue
  30  continue
  20  continue
  10  continue

      return
      end

      subroutine cmat1(raux,t,ntmax,isolv,rowstart,rowend)
C
C     THIS SUBROUTINE SETS UP THE B- MATRIX FOR LINEARIZED CC
C     EQUATIONS T1-T1 AND T1-T2 BLOCK   
C
C     FOR DETAILS CHECK THE DOC. FILE
C
C
      implicit real*8 (a-h,o-z)
      include 'mpif.h'
      integer row,rowstart,rowend
c     common block
      common/jvalue/orbc(MNOCC),orbe(MNEXC)
      common/scratch/b(MDIM),ad(MDIM),baux(MDIM)
      common/index/ke(MNEXC),kc(MNOCC)
      common/info/nocc,nexcit,iType(MNBAS),iTCount(MNBAS)
      common/parity/isymc(MNOCC),isyme(MNEXC)
      common/symmetry/mtbl(MNS,MNS)
      common/orbital_energy/eorb(MNBAS)
      common/kpavlaue/iqc(MNOCC),iqe(MNEXC)
      common/idra1/idra(MNBAS,MNBAS,0:MXV)
      common/idpa1/idpa(MNBAS,MNBAS)
      common/skip/nnskip(0:MXV,2),n2skip,nsing
c
      common/ims1/f(MNBAS,MNBAS)
c
c     local dimension
      dimension t(MDIM),raux(MDIM)
      dimension kgot(MXV+1),jgot(MXV+1),igot(MXV+1)

      data zero,half,one,two/0.0d0,0.50d0,1.0d0,2.0d0/

      do 10 ia=1,nocc
      iasym=isymc(ia)
      orba=orbc(ia)
      iia=kc(ia)
      iaa=iqc(ia)
      do 20 ip=1,nexcit
      ipsym=isyme(ip)
      orbp=orbe(ip)
      iip=ke(ip)
      iap=iqe(ip)
      if(iasym.ne.ipsym)go to 20
      if(orba.ne.orbp)go to 20
                                                                                
      call init(baux)
      neqn=idpa(iip,iia)
                                                                                
      if(isolv.ne.0)then
      baux(neqn)=baux(neqn)+eorb(iip)-eorb(iia)
      endif
                                                                                
      row=idpa(iip,iia)
      irow=row+1-rowstart
      if(irow.le.0)go to 20
      if(irow.gt.rowend-rowstart+1)go to 20
                                                                                
      do ib=1,nocc
      ibsym=isymc(ib)
      orbb=orbc(ib)
      iib=kc(ib)
      iab=iqc(ib)
                                                                                
      if(orbb.eq.orba)then
      if(iasym.eq.ibsym)then
      in=idpa(iip,iib)
      baux(in)=baux(in)-f(iib,iia)
      endif
      endif
      enddo
                                                                                
      do iq=1,nexcit
      iqsym=isyme(iq)
      iiq=ke(iq)
      iaq=iqe(iq)
      orbq=orbe(iq)
                                                                                
      if(orbp.eq.orbq)then
      if(ipsym.eq.iqsym)then
      in=idpa(iiq,iia)
      baux(in)=baux(in)+f(iip,iiq)
      endif
      endif
      enddo
                                                                                
      do 30 ib=1,nocc
      ibsym=isymc(ib)
      orbb=orbc(ib)
      iib=kc(ib)
      iab=iqc(ib)
                                                                                
      do 40 iq=1,nexcit
      iqsym=isyme(iq)
      iiq=ke(iq)
      iaq=iqe(iq)
      orbq=orbe(iq)
                                                                                
      if(iqsym.ne.ibsym)go to 45
      if(orbb.ne.orbq)go to 45
c     one body f(iib,iiq)
      npa=idra(iip,iia,0)
      nqb=idra(iiq,iib,0)
      ipasym=mtbl(ipsym,iasym)
      in=ntloc(npa,nqb,ipasym,0)
      fact=dsqrt((two*orbb+1)/(two*orba+1))
      baux(in)=baux(in)+f(iib,iiq)*fact
      call findk(orba,orbq,orbp,orbb,iaa,iaq,iap,iab,kgot,kmax)
      if(kmax.ne.0)then
      do kloop=1,kmax
      kk=kgot(kloop)
      npb=idra(iip,iib,kk)
      nqa=idra(iiq,iia,kk)
      ipbsym=mtbl(ipsym,ibsym)
      in=ntloc(npb,nqa,ipbsym,kk)
      fact=(-1)**(orbp+orbb+kk)/(2*orba+1)
      baux(in)=baux(in)+f(iib,iiq)*fact
      enddo
      endif
 45   continue
                                                                                
      if(iqsym.ne.ibsym)go to 50
      if(orbb.ne.orbq)go to 50
                                                                                
      in=idpa(iiq,iib)
                                                                                
c     <pb|aq><q|t1|b>
                                                                                
      rk=vint(iip,iib,iia,iiq,0)
      f4=(-1)**(orbb+orbp+1+0)
      twoe=rk*f4
      fact=dsqrt((two*orbb+1)/(two*orba+1))
      baux(in)=baux(in)+twoe*fact
                                                                                
c     direct part <pb|v|qa><q|t1|b>
      call findk(orba,orbb,orbp,orbq,iaa,iab,iap,iaq,kgot,
     $kmax)
      if(kmax.eq.0) go to 50
      do  60 nloop=1,kmax
      ak=dfloat(kgot(nloop))
      kk=kgot(nloop)
      rk=vint(iip,iib,iiq,iia,kk)
      f4=(-1)**(kk+orbb+orbp+1)
      twoe=rk*f4
      fact=one/(two*orba+one)
      isign=(-one)**(orbb+orba+kk)
      baux(in)=baux(in)+twoe*fact*isign
  60  continue
  50  continue
                                                                                
c     insert <pb|v|qr><qr|t2|ab>
      do 70 ir=1,nexcit
      iir=ke(ir)
      iar=iqe(ir)
      irsym=isyme(ir)
      orbr=orbe(ir)
                                                                                
      ibrsym=mtbl(ibsym,irsym)
      ibqsym=mtbl(ibsym,iqsym)
      iaqsym=mtbl(iasym,iqsym)
      ipqsym=mtbl(ipsym,iqsym)
      if(ipqsym.ne.ibrsym)go to 70
c     find allowed k for <pb|v|qr>
c     diagram <pb|v|qr><qr|t2|ab>
      call findk(orbp,orbq,orbb,orbr,iap,iaq,iab,iar,igot,
     $imax)
      if(imax.eq.0)go to 70
      do 80 mloop=1,imax
      kk=igot(mloop)
      ak=dfloat(kk)
      ipbdis=idis(iip,iib)
      iqrloc=indxr(iiq)+iir
      if(iib.gt.iip)iqrloc=indxr(iir)+iiq
      rk=vint(iip,iib,iiq,iir,kk)
      f4=(-1)**(orbp+orbb+kk+1)
      twoe=rk*f4
      fact=one/((two*ak+one)*(two*orba+one))
      isign=(-one)**(orbp+orbq+orbb+orbr)
      naq=idra(iiq,iia,kk)
      nbr=idra(iir,iib,kk)
      in=ntloc(naq,nbr,iaqsym,kk)
      baux(in)=baux(in)+twoe*fact*isign
                                                                                
c     exchange diagram <pb|v|qr><qr|t2|ba>
      call findk(orbb,orbq,orbr,orba,iab,iaq,iar,iaa,jgot,
     $jmax)
      if(jmax.eq.0)go to 80
      do 90 kloop=1,jmax
      ll=jgot(kloop)
      al=dfloat(ll)
      nbq=idra(iiq,iib,ll)
      nar=idra(iir,iia,ll)
      in=ntloc(nbq,nar,ibqsym,ll)
      call sixj(ak,orbr,orbb,al,orbq,orba,s6j)
      baux(in)=baux(in)+twoe*s6j*isign/(two*orba+1)
   90 continue
   80 continue
   70 continue
                                                                                
c     insert <bg|v|aq><pq|t2|bg> (direct part)
      do 100 ig=1,nocc
      iig=kc(ig)
      iag=iqc(ig)
      igsym=isymc(ig)
      orbg=orbc(ig)
                                                                                
      iabsym=mtbl(iasym,ibsym)
      iqgsym=mtbl(iqsym,igsym)
      ipgsym=mtbl(ipsym,igsym)
      if(iabsym.ne.iqgsym)go to 100
c     allowed k for <bg|v|aq)
c     insert <bg|v|aq><pq|t2|bg> (direct part)
      call findk(orba,orbb,orbq,orbg,iaa,iab,iaq,iag,igot,
     $imax)
      if(imax.eq.0)go to 110
      do 120 mloop=1,imax
      kk=igot(mloop)
      ak=dfloat(kk)
      nbp=idra(iip,iib,kk)
      ngq=idra(iiq,iig,kk)
      in=ntloc(nbp,ngq,iqgsym,kk)
      rk=vint(iib,iig,iia,iiq,kk)
      f4=(-one)**(orbg+orbb+one+ak)
      twoe=rk*f4
      isign=(-one)**(orbg+orbb+orbq+orba)
      fact=one/((two*ak+one)*(two*orba+one))
      baux(in)=baux(in)-twoe*fact*isign
c     <bg|v|aq><pq|t2|gb>
      call findk(orbp,orbg,orbq,orbb,iap,iag,iaq,iab,jgot,
     $jmax)
      if(jmax.eq.0)go to 120
      do 130 kloop=1,jmax
      ll=jgot(kloop)
      al=dfloat(ll)
      ngp=idra(iip,iig,ll)
      nbq=idra(iiq,iib,ll)
      in=ntloc(ngp,nbq,ipgsym,ll)
      call sixj(ak,orbp,orbb,al,orbq,orbg,s6j)
      baux(in)=baux(in)-twoe*s6j*isign/(two*orba+one)
  130 continue
  120 continue
  110 continue
  100 continue
   40 continue
   30 continue

      xsum=0
      do in=1,ntmax
      xsum=xsum+baux(in)*t(in)
      enddo
      raux(neqn)=xsum
  20  continue
  10  continue

      return
      end


      subroutine breit(ia,ib,ic,id,orba,orbb,orbc,orbd,
     :iaa,iab,iac,iad,l1,brval,ibra)
      implicit real*8(a-h,o-z)
**************************************************************************
** Coded by Bijaya K. Sahoo on 20.01.2005                             ***|
**************************************************************************
       include 'mpif.h'
      common /cons/zero,half,tenth,one,two,three,ten

      kapa=-int((orba+half)*iaa)
      kapb=-int((orbb+half)*iab)
      kapc=-int((orbc+half)*iac)
      kapd=-int((orbd+half)*iad)

      dl1=dfloat(l1)
      brval=zero

      d1=dr(orba,dl1,orbc,half,zero,-half)
      d2=dr(orbb,dl1,orbd,half,zero,-half)

      factor=dsqrt((2*orba+one)*(2*orbb+one)*
     :(2*orbc+one)*(2*orbd+one))

      do 50 k=int(abs(l1-1)),l1+1

*  Radial part of Breit interaction

      fact=0.0
      fact=select(iaa,iac,orba,orbc,k)*
     : select(iab,iad,orbb,orbd,k)

      if(fact.ne.zero)then

      call slater_br1(ia,ib,ic,id,k,val1,val2,val3,val4)

      brval=brval+factor*d1*d2*fact*
     : (rk(kapa,kapb,kapc,kapd,l1,k,1)*val1
     :+ rk(kapa,kapb,kapc,kapd,l1,k,2)*val2
     :+ rk(kapa,kapb,kapc,kapd,l1,k,3)*val3
     :+ rk(kapa,kapb,kapc,kapd,l1,k,4)*val4)

*   Derivative terms

      if(ibra.ne.3)goto 60

      if(k.eq.l1)goto 60

      brval=brval-factor*d1*d2*fact*
     : (sk(kapa,kapb,kapc,kapd,l1,k,1)*val1
     :+ sk(kapa,kapb,kapc,kapd,l1,k,2)*val2
     :+ sk(kapa,kapb,kapc,kapd,l1,k,3)*val3
     :+ sk(kapa,kapb,kapc,kapd,l1,k,4)*val4)

      if(k.eq.l1-1)then

      call slater_br2(ia,ib,ic,id,k,val5,val6,val7,val8)
      call slater_br2(ib,ia,id,ic,k,val9,val10,val11,val12)

      brval=brval+factor*d1*d2*fact*
     : (gk(kapa,kapb,kapc,kapd,l1,5)*val5
     :+ gk(kapa,kapb,kapc,kapd,l1,6)*val6
     :+ gk(kapa,kapb,kapc,kapd,l1,7)*val7
     :+ gk(kapa,kapb,kapc,kapd,l1,8)*val8
     :+ gk(kapb,kapa,kapd,kapc,l1,5)*val9
     :+ gk(kapb,kapa,kapd,kapc,l1,6)*val10
     :+ gk(kapb,kapa,kapd,kapc,l1,7)*val11
     :+ gk(kapb,kapa,kapd,kapc,l1,8)*val12)

      endif

  60  continue

      endif

  50  continue

      return
      end

      function rk(kapa,kapb,kapc,kapd,l1,k,i)
      implicit real*8(a-h,o-z)
       include 'mpif.h'
      common /cons/zero,half,tenth,one,two,three,ten

      rk=0.0

      kp1=(kapc-kapa)
      kp2=(kapd-kapb)

      if(k.eq.l1-1)then
      A=one/dfloat(l1*(2*l1-1))

      if(i.eq.1)then
      rk=A*(k1+l1)*(k2+l1)
      endif

      if(i.eq.2)then
      rk=A*(k1-l1)*(k2-l1)
      endif

      if(i.eq.3)then
      rk=A*(k1+l1)*(k2-l1)
      endif

      if(i.eq.4)then
      rk=A*(k1-l1)*(k2+l1)
      endif

      endif

      if(k.eq.l1)then
      A=-dfloat((kapa+kapc)*(kapb+kapd))/dfloat(l1*(l1+1))
      rk=A
      endif

      if(k.eq.l1+1)then
      A=one/dfloat((l1+1)*(2*l1+3))

      if(i.eq.1)then
      rk=A*(k1-l1-1)*(k2-l1-1)
      endif

      if(i.eq.2)then
      rk=A*(k1+l1+1)*(k2+l1+1)
      endif

      if(i.eq.3)then
      rk=A*(k1-l1-1)*(k2+l1+1)
      endif

      if(i.eq.4)then
      rk=A*(k1+l1+1)*(k2-l1-1)
      endif

      endif

      return
      end

      function sk(kapa,kapb,kapc,kapd,l1,k,i)
      implicit real*8(a-h,o-z)
       include 'mpif.h'
      common /cons/zero,half,tenth,one,two,three,ten

      sk=0.0

      kp1=(kapc-kapa)
      kp2=(kapd-kapb)

      if(k.eq.l1-1)then

      if(i.eq.1)then
      sk=dfloat((kp1-k-1)*(kp2-k-1))/dfloat((2*k+3)*(2*k+1))
      endif

      if(i.eq.2)then
      sk=dfloat((kp1-k-1)*(kp2+k+1))/dfloat((2*k+3)*(2*k+1))
      endif

      if(i.eq.3)then
      sk=dfloat((kp1+k+1)*(kp2-k-1))/dfloat((2*k+3)*(2*k+1))
      endif

      if(i.eq.4)then
      sk=dfloat((kp1+k+1)*(kp2+k+1))/dfloat((2*k+3)*(2*k+1))
      endif

      endif

      if(k.eq.l1+1)then

      if(i.eq.1)then
      sk=dfloat((kp1+k)*(kp2+k))/dfloat((2*k-1)*(2*k+1))
      endif

      if(i.eq.2)then
      sk=dfloat((kp1+k)*(kp2-k))/dfloat((2*k-1)*(2*k+1))
      endif

      if(i.eq.3)then
      sk=dfloat((kp1-k)*(kp2+k))/dfloat((2*k-1)*(2*k+1))
      endif

      if(i.eq.4)then
      sk=dfloat((kp1-k)*(kp2-k))/dfloat((2*k-1)*(2*k+1))
      endif

      endif

      return
      end

      function gk(kapa,kapb,kapc,kapd,k,i)
      implicit real*8(a-h,o-z)
       include 'mpif.h'
      common /cons/zero,half,tenth,one,two,three,ten

      gk=0.0

      kp1=(kapc-kapa)
      kp2=(kapd-kapb)

      if(i.eq.1)then
      gk=dfloat((kp1+k+1)*(kp2-k))/dfloat(2*(2*k+1))
      endif

      if(i.eq.2)then
      gk=dfloat((kp1+k+1)*(kp2+k))/dfloat(2*(2*k+1))
      endif

      if(i.eq.3)then
      gk=dfloat((kp1-k-1)*(kp2-k))/dfloat(2*(2*k+1))
      endif

      if(i.eq.4)then
      gk=dfloat((kp1-k-1)*(kp2+k))/dfloat(2*(2*k+1))
      endif

      return
      end


c***********************************************************************
      subroutine slater_br1(ia,ib,ic,id,k,val1,val2,val3,val4)
c***********************************************************************
      implicit real*8(a-h, o-z)
       include 'mpif.h'
      common/grid/r(MN),rp(MN),rpor(MN),rnt,h,hp,n
     :      /ncc/c1,c2,c3,c4
      common /cons/zero,half,tenth,one,two,three,ten
     :/wave/pz(MNBAS),pf(MN,MNBAS),qf(MN,MNBAS),mf(MNBAS)
     :/tatb/ta(MN),tb(MN),mtp
     :/tatb1/ta1(MN),tb1(MN)
     :/tatb2/ta2(MN),tb2(MN)

      mtp=n

      do i=1,mtp
      ta(i)=zero
      tb1(i)=zero
      tb2(i)=zero
      enddo

      val1=0.0
      val2=0.0
      val3=0.0
      val4=0.0

      call yzk_br(k,ib,id)

      do i = 2,mtp
         ta(i)=qf(i,ia)*pf(i,ic)*rpor(i)*tb2(i)
      enddo

      call quad (val1)

      do i=1,mtp
      ta(i)=zero
      enddo

      do i = 2,mtp
         ta(i)=qf(i,ia)*pf(i,ic)*rpor(i)*tb1(i)
      enddo
      call quad (val2)

      do i=1,mtp
      ta(i)=zero
      enddo

      do i = 2,mtp
         ta(i)=pf(i,ia)*qf(i,ic)*rpor(i)*tb2(i)
      enddo
      call quad (val3)

      do i=1,mtp
      ta(i)=zero
      enddo

      do i = 2,mtp
         ta(i)=pf(i,ia)*qf(i,ic)*rpor(i)*tb1(i)
      enddo
      call quad (val4)

      return
      end

c***********************************************************************
      subroutine slater_br2(ia,ib,ic,id,k,val1,val2,val3,val4)
c***********************************************************************
*
      implicit real*8(a-h, o-z)
       include 'mpif.h'
      common/grid/r(MN),rp(MN),rpor(MN),rnt,h,hp,n
     :      /ncc/c1,c2,c3,c4
      common /cons/zero,half,tenth,one,two,three,ten
     :/wave/pz(MNBAS),pf(MN,MNBAS),qf(MN,MNBAS),mf(MNBAS)
     :/tatb/ta(MN),tb(MN),mtp
     :/tatb1/ta1(MN),tb1(MN)
     :/tatb2/ta2(MN),tb2(MN)

      mtp=n

      do i=1,mtp
      ta(i)=zero
      ta1(i)=zero
      ta2(i)=zero
      enddo

       val1=0.0
       val11=0.0
       val12=0.0
       val2=0.0
       val21=0.0
       val22=0.0
       val3=0.0
       val31=0.0
       val32=0.0
       val4=0.0
       val41=0.0
       val42=0.0

      call yzk_br(k+2,ib,id)

      do i = 2,mtp
         ta(i)=qf(i,ia)*pf(i,ic)*rpor(i)*ta2(i)
      enddo
      call quad (val11)

      do i = 2,mtp
         ta(i)=qf(i,ia)*pf(i,ic)*rpor(i)*ta1(i)
      enddo
      call quad (val21)

      do i = 2,mtp
         ta(i)=pf(i,ia)*qf(i,ic)*rpor(i)*ta2(i)
      enddo
      call quad (val31)

      do i = 2,mtp
         ta(i)=pf(i,ia)*qf(i,ic)*rpor(i)*ta1(i)
      enddo
      call quad (val41)

      do i=1,mtp
      ta(i)=zero
      ta1(i)=zero
      ta2(i)=zero
      enddo

      call yzk_br(k,ib,id)

      do i = 2,mtp
         ta(i)=qf(i,ia)*pf(i,ic)*rpor(i)*ta2(i)
      enddo
      call quad (val12)

      do i = 2,mtp
         ta(i)=qf(i,ia)*pf(i,ic)*rpor(i)*ta1(i)
      enddo
      call quad (val22)

      do i = 2,mtp
         ta(i)=pf(i,ia)*qf(i,ic)*rpor(i)*ta2(i)
      enddo
      call quad (val32)

      do i = 2,mtp
         ta(i)=pf(i,ia)*qf(i,ic)*rpor(i)*ta1(i)
      enddo
      call quad (val42)

      val1=(val11-val12)
      val2=(val21-val22)
      val3=(val31-val32)
      val4=(val41-val42)

      return
      end

************************************************************************
      subroutine yzk_br(k,i,j)
*   ----------------   section 09   subprogram 32   ----------------   *
*                                                                      *
*   this subroutine evaluates hartree y- and z-functions:              *
*                                                                      *
*               (k)            (k)           (k)                       *
*              y   (i,j;r) =  z   (i,j;r) + w   (i,j;r)                *
*                                                                      *
*   where                                                              *
*                                                                      *
*    (k)                                                               *
*   z   (i,j;r) =  i ( (s/r)  p (s)*q (s)  ; 0 - r )                   *
*                             i     j                                  *
*                                                                      *
*   where                                                              *
*                                                                      *
*    (k)                    k+1                                        *
*   w   (i,j;r) =  i ( (r/s)    p (s)*q (s)  ; r -                     *
*                                i     j      infinity )               *
*                                                                      *
*   where  i ( g(r,s) ; range )  denotes the integral of g(r,s) over   *
*   range  in  s .  the y-function is tabulated in  common/tatb/  in   *
*   array  tb , the z-function in array ta .                           *
*                                                                      *
************************************************************************
*
      implicit doubleprecision (a-h, o-z)
       include 'mpif.h'

      dimension rhop(MN),rhoq(MN),rttk(MN),wk1(MN),temp1(MN),
     :yk1(MN),zk1(MN),wk2(MN),temp2(MN),yk2(MN),zk2(MN)
*
      common/cnc5/cnc5c(2:5,2:4)
     :      /cons/zero,half,tenth,one,two,three,ten
     :      /def4/accy,nscf,nsic,nsolv
     :      /grid/r(MN),rp(MN),rpor(MN),rnt,h,hp,n
     :      /ncc/c1,c2,c3,c4
     :      /tatb/ta(MN),tb(MN),mtp
     :      /tatb1/ta1(MN),tb1(MN)
     :      /tatb2/ta2(MN),tb2(MN)
     :      /wave/pz(MNBAS),pf(MN,MNBAS),qf(MN,MNBAS),
     :                                      mf(MNBAS)
*
      equivalence (ta1(1),zk1(1)),(tb1(1),yk1(1))
      equivalence (ta2(1),zk2(1)),(tb2(1),yk2(1))
*                         k
*   for  k > 0  compute  r   and store in  rttk
      if (k .gt. 0) then
         do 1 ii = 2,n+4
            rttk(ii) = r(ii)**k
    1    continue
      endif

*   determine maximum tabulation point as location beyond which
*   rhop  (see comment statements below) would be zero; determine
*   other important locations
      mtpp1 = mtp+1
      mtpp3 = mtp+3
      mtpp4 = mtp+4
*   compute rp(s)* P(s)*Q (s) and store in rhop
*                  i    j     i     j

      do 2 ii = 2,mtp
      rhop(ii) = rp(ii)*pf(ii,i)*qf(ii,j)
      rhoq(ii) = rp(ii)*qf(ii,i)*pf(ii,j)
    2 continue

*   fill array temp with r**k * rhop
      temp1(1) = zero
      temp2(1) = zero

      if (k .eq. 0) then

         do 3 ii = 2,mtp
            temp1(ii) = rhop(ii)
            temp2(ii) = rhoq(ii)
    3    continue

      else

         do 4 ii = 2,mtp
            temp1(ii) = rttk(ii)*rhop(ii)
            temp2(ii) = rttk(ii)*rhoq(ii)
    4    continue

      endif

*   set an additional four points to zero

      do 5 ii = mtpp1,mtpp4
         temp1(ii) = zero
         temp2(ii) = zero
    5 continue

*   compute the first few values of  r  * zk  using semi-open
*   newton-cotes formulae

      zk1(1) = zero
      zk2(1) = zero

      do 7 ii = 2,4
         sum1 = zero
         sum2 = zero

         do 6 kk = 2,5
            sum1 = sum1+cnc5c(kk,ii)*temp1(kk)
            sum2 = sum2+cnc5c(kk,ii)*temp2(kk)
    6    continue

         zk1(ii) = sum1
         zk2(ii) = sum2

    7 continue

*                         k
*   compute remainder of r  * zk: march out to mtp+3; use closed
*   newton-cotes formula

      do 8 ii = 5,mtpp3

         zk1(ii) = zk1(ii-4)+c1*(temp1(ii-4)+temp1(ii))
     :                    +c2*(temp1(ii-3)+temp1(ii-1))
     :                    +c3*temp1(ii-2)

         zk2(ii) = zk2(ii-4)+c1*(temp2(ii-4)+temp2(ii))
     :                    +c2*(temp2(ii-3)+temp2(ii-1))
     :                    +c3*temp2(ii-2)

    8 continue

*                                       k   (k)
*   determine the asymptotic value of  r * z
*                   (0)
*   correction to  z   : in the manner of  c froese fischer,
*   the hartree-fock method for atoms, john wiley & sons,
*   new york, 1977, p 235.

      if (k .eq. 0) then

         if (i .eq. j) then
            zklim1 = one
            zklim2 = one
         else
            zklim1 = zero
            zklim2 = zero
         endif

         do 10 kk = mtpp3,mtp,-1
            dif1 = zk1(kk)-zklim1
            dif2 = zk2(kk)-zklim2

            if (abs (dif1) .gt. accy) then
               do 9 ii = kk,2,-4
                  zk1(ii) = zk1(ii)-dif1
    9          continue
            endif

            if (abs (dif2) .gt. accy) then
               do ii = kk,2,-4
                  zk2(ii) = zk2(ii)-dif2
               enddo
            endif

   10    continue

      else

         zklim1 = zk1(mtpp3)
         zklim2 = zk2(mtpp3)

      endif

*   tabulate  zk  for entire internal grid

      if (k .eq. 0) then

         do 11 ii = mtpp4,n
            zk1(ii) = zklim1
            zk2(ii) = zklim2
   11    continue

      else

         do 12 ii = 2,mtpp3
            zk1(ii) = zk1(ii)/rttk(ii)
            zk2(ii) = zk2(ii)/rttk(ii)
   12    continue

         do 13 ii = mtpp4,n
            zk1(ii) = zklim1/rttk(ii)
            zk2(ii) = zklim2/rttk(ii)
   13    continue

      endif

*   start array wk / r**(k+1)

      np4 = n+4

      do 14 ii = np4,mtpp1,-1
         wk1(ii) = zero
         wk2(ii) = zero
   14 continue

*             k+1
*   compute  r       and store in rttk

      if (k .gt. 0) then

         do 15 ii = 2,n
            rttk(ii) = rttk(ii)*r(ii)
   15    continue

      endif

*   fill array temp with rhop / r**(k+1) ; set temp(1) = zero
*   to avoid 0/0 case

      temp1(1) = zero
      temp2(1) = zero

      if (k .eq. 0) then

         do 16 ii = 2,mtp
            temp1(ii) = rhop(ii)/r(ii)
            temp2(ii) = rhoq(ii)/r(ii)
   16    continue

      else

         do 17 ii = 2,mtp
            temp1(ii) = rhop(ii)/rttk(ii)
            temp2(ii) = rhoq(ii)/rttk(ii)
   17    continue

      endif

*   compute remainder of wk / r**(k+1): march in to the origin

      do 18 ii = mtp,2,-1

         wk1(ii) = wk1(ii+4)+c1*(temp1(ii)+temp1(ii+4))
     :                    +c2*(temp1(ii+1)+temp1(ii+3))
     :                    +c3*(temp1(ii+2))

         wk2(ii) = wk2(ii+4)+c1*(temp2(ii)+temp2(ii+4))
     :                    +c2*(temp2(ii+1)+temp2(ii+3))
     :                    +c3*(temp2(ii+2))

   18 continue

      wk1(1) = zero
      wk2(1) = zero

*   compute wk

      if (k .eq. 0) then

         do 19 ii = 2,mtp
            wk1(ii) = wk1(ii)*r(ii)
            wk2(ii) = wk2(ii)*r(ii)
   19    continue

      else

         do 20 ii = 2,mtp
            wk1(ii) = wk1(ii)*rttk(ii)
            wk2(ii) = wk2(ii)*rttk(ii)
   20    continue

      endif

*   assemble solution

      yk1(1) = zero
      yk2(1) = zero

      do 21 ii = 2,n
         yk1(ii) = zk1(ii)+wk1(ii)
         yk2(ii) = zk2(ii)+wk2(ii)
   21 continue

      return
      end


      function select1(ia,ib,orba,orbb,k)
      implicit real*8 (a-h,o-z)
       include 'mpif.h'

      select1=0.0

      dk1=dfloat(k)
      da=dfloat(ia)
      db=dfloat(ib)
      select1=(1.0d0-(da*db*(-1)**(orba+orbb+dk1)))/2.0d0

      return
      end

      function select(ia,ib,orba,orbb,k)
      implicit real*8 (a-h,o-z)
       include 'mpif.h'

      select=0.0

      dk1=dfloat(k)
      da=dfloat(ia)
      db=dfloat(ib)
      select=(1.0d0+(da*db*(-1)**(orba+orbb+dk1)))/2.0d0

      return
      end

c*******************************************************************************
c                                                                              *
c                   SUBROUTINE  YZK_BREIT(I,J,K)                               *
c                                                                              *
c      This subroutine evaluates hartree y and z functions                     *
c        (k)             (k)             (k)                                   *
c       y   #$(i,j;r) = z   #$(i,j;r) + w   #$(i,j;r)                          *
c                                                                              *
c       where #$ can be either p*q or q*p                                      *
c                                                                              *
c       (k)                  k                                                 *
c      z   #$(i,j;r)=i ((s/r)  # (s) * $ (s)  ; 0-r)                           *
c                               i       j                                      *
c       (k)                  k+1                                               *
c      w   #$(i,j;r)=i ((r/s)  # (s) * $ (s)  ; r-infinity)                    *
c                               i       j                                      *
c                                                                              *
c     where i(g(r,s);range) denotes the integral of g(r,s)                     *
c     over range in s                                                          *
c*******************************************************************************

         subroutine yzk_breit(k,i,j)

      implicit real*8(a-h, o-z)
      include 'mpif.h'
      dimension rhoppq(MN),wkpq(MN+5),temppq(MN),zkpq(MN),
     :          rhopqp(MN),wkqp(MN+5),tempqp(MN),zkqp(MN),rttk(MN)
      common /cons/zero,half,tenth,one,two,three,ten
     :  /cnc5/cnc5c(2:5,2:4)
     :  /grid/r(MN),rp(MN),rpor(MN),rnt,h,hp,n
     :  /ncc/c1,c2,c3,c4
     :  /dumm/tmp1(MN),tmp2(MN),tmp3(MN),tmp4(MN)
     :  /wave/pz(MNBAS),pf(MN,MNBAS),qf(MN,MNBAS),mf(MNBAS)
     :  /basis/gl(MN,MNBAS),gs(MN,MNBAS),alpha(MNBAS,MNS),mg(MNBAS)

c------------------------------------------------------------------------------
c           for  k > 0  compute  r^k   and store in  rttk
c------------------------------------------------------------------------------
        if (k.gt.0) then
          do 121 ii = 2,n
            rttk(ii) = r(ii)**k
 121   continue
        endif
c------------------------------------------------------------------------------
c
c    determine maximum tabulation point as location beyond which rhop  (see
c     comment statements below) would be zero; determine other important
c                               locations
c------------------------------------------------------------------------------

        mtp = mf(i)
        mtpp1 = mtp+1
        mtpp3 = mtp+3
        mtpp4 = mtp+4
c------------------------------------------------------------------------------
c     compute rp(s)* p_i (s)*q_j (s) and store in rhoppq
c                    q_i (s) p_j (s)              rhopqp
c------------------------------------------------------------------------------
         do 221 ii = 2,mtp
         rhoppq(ii) = rp(ii)*pf(ii,i)*qf(ii,j)
         rhopqp(ii) = rp(ii)*qf(ii,i)*pf(ii,j)
  221    continue
c------------------------------------------------------------------------------
c                   fill array temppq with r**k * rhoppq
c                              tempqp             rhopqp
c------------------------------------------------------------------------------
         temppq(1) = zero
         tempqp(1) = zero
         if (k.eq.0) then
           do 231 ii = 2,mtp
             temppq(ii) = rhoppq(ii)
             tempqp(ii) = rhopqp(ii)
  231     continue
         else
        do 241 ii = 2,mtp
         temppq(ii) = rttk(ii)*rhoppq(ii)
         tempqp(ii) = rttk(ii)*rhopqp(ii)
 241    continue
         endif
c------------------------------------------------------------------------------
c                   set an additional four points to zero
c------------------------------------------------------------------------------
        do 251 ii = mtpp1,mtpp4
           temppq(ii) = zero
           tempqp(ii) = zero
 251     continue
c------------------------------------------------------------------------------
c  compute first few values of  r^k*zk using semi-open newton-cotes formulae
c------------------------------------------------------------------------------
         zkpq(1) = zero
         zkqp(1) = zero
       do 271 ii = 2,4
         sumpq = zero
         sumqp = zero
       do 261 kk = 2,5
          sumpq = sumpq+cnc5c(kk,ii)*temppq(kk)
          sumqp = sumqp+cnc5c(kk,ii)*tempqp(kk)
 261    continue
         zkpq(ii) = sumpq
         zkqp(ii) = sumqp
 271    continue
c------------------------------------------------------------------------------
c        compute remainder of r^k*zk: march out to mtp+3; use closed
c                         newton-cotes formula
c------------------------------------------------------------------------------
       do 281 ii = 5,mtpp3
          zkpq(ii) = zkpq(ii-4)+c1*(temppq(ii-4)+temppq(ii  ))
     :                       +c2*(temppq(ii-3)+temppq(ii-1))
     :                       +c3* temppq(ii-2)
          zkqp(ii) = zkqp(ii-4)+c1*(tempqp(ii-4)+tempqp(ii  ))
     :                       +c2*(tempqp(ii-3)+tempqp(ii-1))
     :                       +c3* tempqp(ii-2)
 281     continue
c------------------------------------------------------------------------------
c      determine the asymptotic value of  r^k*z^(k) correction to  z^(k)
c         : in the manner of  c froese fischer,the hartree-fock method for
c            atoms, john wiley & sons, new york, 1977, p 235.
c------------------------------------------------------------------------------
        if (k .eq. 0) then
         zklimpq = zero
         zklimqp = zero
      else
         zklimpq = zkpq(mtpp3)
         zklimqp = zkqp(mtpp3)
      endif
c------------------------------------------------------------------------------
c                     tabulate  zk  for entire internal lgrid
c------------------------------------------------------------------------------
       if (k .eq. 0) then
        do 212 ii = mtpp4,n
          zkpq(ii) = zklimpq
          zkqp(ii) = zklimqp
 212    continue
      else
         do 213 ii = 2,mtpp3
          zkpq(ii) = zkpq(ii)/rttk(ii)
          zkqp(ii) = zkqp(ii)/rttk(ii)
 213   continue
        do 214 ii = mtpp4,n
          zkpq(ii) = zklimpq/rttk(ii)
          zkqp(ii) = zklimqp/rttk(ii)
 214   continue
      endif
c------------------------------------------------------------------------------
c                          start array wk / r**(k+1)
c------------------------------------------------------------------------------

          np4 = n+4
         do 215 ii = np4,mtpp1,-1
            wkpq(ii) = zero
            wkqp(ii) = zero
 215    continue
c------------------------------------------------------------------------------
c               compute  r^(k+1)  and store in rttk
c------------------------------------------------------------------------------
         if (k .gt. 0) then
         do 216 ii = 2,n
            rttk(ii) = rttk(ii)*r(ii)
 216    continue
         endif

c------------------------------------------------------------------------------
c        fill array temp with rhop / r**(k+1) ; set temp(1) = zero
c                          to avoid 0/0 case
c------------------------------------------------------------------------------
         temppq(1) = zero
         tempqp(1) = zero
         if (k .eq. 0) then
         do 217 ii = 2,mtp
            temppq(ii) = rhoppq(ii)/r(ii)
            tempqp(ii) = rhopqp(ii)/r(ii)
 217    continue
         else
         do 218 ii = 2,mtp
            temppq(ii) = rhoppq(ii)/rttk(ii)
            tempqp(ii) = rhopqp(ii)/rttk(ii)
 218    continue
         endif
c------------------------------------------------------------------------------
c         compute remainder of wk / r**(k+1): march in to the origin
c------------------------------------------------------------------------------
        do 219 ii = mtp,2,-1
         wkpq(ii) = wkpq(ii+4)+c1*(temppq(ii  )+temppq(ii+4))
     :                    +c2*(temppq(ii+1)+temppq(ii+3))
     :                    +c3*(temppq(ii+2))
         wkqp(ii) = wkqp(ii+4)+c1*(tempqp(ii  )+tempqp(ii+4))
     :                    +c2*(tempqp(ii+1)+tempqp(ii+3))
     :                    +c3*(tempqp(ii+2))
 219   continue
        wkpq(1) = zero
        wkqp(1) = zero
c------------------------------------------------------------------------------
c                   compute wk
c-------------------------------------------------------------------------------

         if (k .eq. 0) then
         do 220 ii = 2,mtp
            wkpq(ii) = wkpq(ii)*r(ii)
            wkqp(ii) = wkqp(ii)*r(ii)
 220    continue
        else
         do 421 ii = 2,mtp
            wkpq(ii) = wkpq(ii)*rttk(ii)
            wkqp(ii) = wkqp(ii)*rttk(ii)
 421    continue
        endif
c-----------------------------------------------------------
c                             assemble solution
c----------------------------------------------------------
           tmp1(1) = zero
           tmp2(1) = zero
        do 422 ii = 2,n
           tmp1(ii) = zkpq(ii)+wkpq(ii)
           tmp2(ii) = zkqp(ii)+wkqp(ii)
 422  continue
       return

       end



c ****************************************************************************
c                                                                            *
c        SUBROUTINE slater_breit                                             *
c                                                                            *
c                                                                            *
c   The value of the function is the slater integral as normally defined in  *
c   terms  of the four sets of quantum numbers a,b,c,d of breit interaction  *
c                                                                            *
c                                                                            *
c              subroutines called : quad,yzk_breit                           *
c                                                                            *
c ****************************************************************************

          subroutine slater_breit(ia,ib,ic,id,k)

      implicit real*8(a-h, o-z)
      include  'mpif.h'

      common /cons/zero,half,tenth,one,two,three,ten
     :  /grid/r(MN),rp(MN),rpor(MN),rnt,h,hp,n
     :  /ncc/c1,c2,c3,c4
     :  /dumm/tmp1(MN),tmp2(MN),tmp3(MN),tmp4(MN)
     :  /wave/pz(MNBAS),pf(MN,MNBAS),qf(MN,MNBAS),mf(MNBAS)
     :  /basis/gl(MN,MNBAS),gs(MN,MNBAS),alpha(MNBAS,MNS),mg(MNBAS)
     : /brresult/result1,result2,result3,result4


       mtp=n
       do i=1,mtp
       tmp1(i)=0
       tmp2(i)=0
       tmp3(i)=0
       tmp4(i)=0
         enddo

             call yzk_breit(k,ib,id)

         

c-----------------------------------------------------------------
c          compute PQPQ component of the slater integral
c-----------------------------------------------------------------

           do   i=2,mtp
         tmp3(i)=pf(i,ia)*qf(i,ic)*rpor(i)*tmp1(i)
           enddo
         call quad1 (tmp3,mtp,result1)

c------------------------------------------------------------------
c          compute QPQP component of the slater integral
c------------------------------------------------------------------

           do   i=2,mtp
         tmp3(i)=qf(i,ia)*pf(i,ic)*rpor(i)*tmp2(i)
           enddo
        call quad1 (tmp3,mtp,result2)

c-------------------------------------------------------------------
c          compute PQQP component of the slater integral
c-------------------------------------------------------------------

           do   i=2,mtp
         tmp3(i)=pf(i,ia)*qf(i,ic)*rpor(i)*tmp2(i)
           enddo
         call quad1 (tmp3,mtp,result3)

c-------------------------------------------------------------------
c         compute QPPQ component of the slater integral
c-------------------------------------------------------------------

       do   i=2,mtp
         tmp3(i)=qf(i,ia)*pf(i,ic)*rpor(i)*tmp1(i)
           enddo
         call quad1 (tmp3,mtp,result4)

           return
           end

c*******************************************************************************
c                                                                              *
c                                                                              *
c  Slater interegals with angular coefficients computed for breit              *
c  interaction for n=k                                                         *
c                                                                              *
c                                                                              *
c*******************************************************************************


      function slater_breitk(ii,jj,kk,ll,orbi,orbj,orbk,orbl,
     :iai,iaj,iak,ial,l1)
      implicit real*8 (a-h,o-z)
      include 'mpif.h'
      common /cons/zero,half,tenth,one,two,three,ten
     : /brresult/result1,result2,result3,result4

        slater_breitk=0.0     

c---------------------------------------------------------------------------
c      kappa values are determined
c---------------------------------------------------------------------------
      
      kapi=-dble((orbi+0.5)*iai)
      kapj=-dble((orbj+0.5)*iaj)
      kapk=-dble((orbk+0.5)*iak)
      kapl=-dble((orbl+0.5)*ial)




c--------------------------------------------------------------------------
c    angular coefficients are calculated
c--------------------------------------------------------------------------

      x=dble(l1*(l1+1))
      p=1/x
      k1=kapi+kapk
      k2=kapj+kapl
      r=-p*dble(k1*k2)

c--------------------------------------------------------------------------
c   angular coefficients are multiplied with the slater integrals
c--------------------------------------------------------------------------


      call slater_breit(ii,jj,kk,ll,l1)
      slater_breitk=(result1+result2+result3+result4)*r

      

      return
      end
c*******************************************************************************
c                                                                              *
c                                                                              *
c  Slater integrals with angular coefficients have been evaluated for breit    *
c  interaction for n=k-1                                                       *
c                                                                              *
c                                                                              *
c*******************************************************************************
      function slater_breitkm1(ii,jj,kk,ll,orbi,orbj,orbk,orbl,
     :iai,iaj,iak,ial,l1)
      implicit real*8 (a-h,o-z)
      include 'mpif.h'
       common /cons/zero,half,tenth,one,two,three,ten
     : /brresult/result1,result2,result3,result4


       slater_breitkm1=0.0

c-----------------------------------------------------------------------------
c      kappa values are determined
c-----------------------------------------------------------------------------
     
      
      kapi=-dble((orbi+0.5)*iai)
      kapj=-dble((orbj+0.5)*iaj)
      kapk=-dble((orbk+0.5)*iak)
      kapl=-dble((orbl+0.5)*ial)

c-----------------------------------------------------------------------------
c     angular coefficients are calculated
c-----------------------------------------------------------------------------

      x=dble(l1*(2*l1-1))
      p=1/x
      k1=kapk-kapi
      k2=kapl-kapj
      r1=p*dble((l1+k1)*(l1+k2))
      r2=p*dble((l1-k1)*(l1-k2))
      r3=p*dble((k1+l1)*(k2-l1))
      r4=p*dble((k1-l1)*(k2+l1))

c-----------------------------------------------------------------------------
c     angular coefficients are multiplied with slater integrals
c-----------------------------------------------------------------------------
      
      l11=l1-1
      call slater_breit(ii,jj,kk,ll,l11)
      slater_breitkm1=result1*r1+result2*r2+result3*r3+result4*r4

      
      return
      end
c***********************************************************************
c                                                                      *
c                                                                      *
c  Slater integrals with angular coefficients  have been evaluated     *
c  for breit  intearction for n=k+1                                    *
c                                                                      *
c***********************************************************************
      function slater_breitkp1(ii,jj,kk,ll,orbi,orbj,orbk,orbl,
     :iai,iaj,iak,ial,l1)
      implicit real*8 (a-h,o-z)
      include 'mpif.h'
       common /cons/zero,half,tenth,one,two,three,ten
     : /brresult/result1,result2,result3,result4


       slater_breitkp1=0.0

c-----------------------------------------------------------------------------
c      kappa values are determined
c-----------------------------------------------------------------------------

      kapi=-dble((orbi+0.5)*iai)
      kapj=-dble((orbj+0.5)*iaj)
      kapk=-dble((orbk+0.5)*iak)
      kapl=-dble((orbl+0.5)*ial)

c-----------------------------------------------------------------------------
c     angular coefficients are calculated
c-----------------------------------------------------------------------------

    
      x=dble((l1+1)*(2*l1+3))
      p=1/x
      k1=kapk-kapi
      k2=kapl-kapj
      r1=p*dble((k1-l1-1)*(k2-l1-1))
      r2=p*dble((k1+l1+1)*(k2+l1+1))
      r3=p*dble((k1-l1-1)*(k2+l1+1))
      r4=p*dble((k1+l1+1)*(k2-l1-1))

c-----------------------------------------------------------------------------
c     angular coefficients are multiplied with slater integrals
c-----------------------------------------------------------------------------

      
      l11=l1+1

      call slater_breit(ii,jj,kk,ll,l11)

      slater_breitkp1=result1*r1+result2*r2+result3*r3+result4*r4


      return
      end

c*************************************************************************
c    Operator strength   have been  calculated  without phase factor     *
c     for n=k part of the breit  interaction                             *
c                                                                        *
c                                                                        *
c*************************************************************************

      function slater_newbr1(ii,jj,kk,ll,orbi,orbj,orbk,orbl,
     :iai,iaj,iak,ial,l1)
      implicit real*8 (a-h,o-z)
       include 'mpif.h'
      data zero,half,one,two/0.0d0,0.50d0,1.0d0,2.0d0/

      slater_newbr1=0.0
      fact=0.0

c--------------------------------------------------------------------------
c      multipoles (k's) are selected for the odd operators
c--------------------------------------------------------------------------

      fact=select(iai,iak,orbi,orbk,l1)*
     : select(iaj,ial,orbj,orbl,l1)



      if(fact.ne.zero)then
      if(l1.ne.0)then


      ak=dfloat(l1)
      f1=dr(orbi,ak,orbk,half,zero,-half)
      f2=dr(orbj,ak,orbl,half,zero,-half)
      f3=dsqrt((two*orbi+1)*(two*orbj+1)*(two*orbk+1)*(two*orbl+1))
c     f4=(-one)**(l1+orbi+orbj+1)
      slater_newbr1=f1*f2*f3*slater_breitk(ii,jj,kk,ll,orbi,
     :orbj,orbk,orbl,iai,iaj,iak,ial,l1)

       
      endif
      endif

      return
      end

c*************************************************************************
c    Operator strength   have been  calculated  without phase factor     *
c     for n=k-1 part of the breit  interaction                           *
c                                                                        *
c                                                                        *
c*************************************************************************

      function slater_newbr2(ii,jj,kk,ll,orbi,orbj,orbk,orbl,
     :iai,iaj,iak,ial,l1)
      implicit real*8 (a-h,o-z)
      include 'mpif.h'
      data zero,half,one,two/0.0d0,0.50d0,1.0d0,2.0d0/

      slater_newbr2=0.0
      fact=0.0

c--------------------------------------------------------------------------
c      multipoles (k's) are selected for the odd operators
c--------------------------------------------------------------------------
      fact=select1(iai,iak,orbi,orbk,l1)*
     : select1(iaj,ial,orbj,orbl,l1)

      if(fact.ne.zero)then
      if(l1.ne.0)then


      ak=dfloat(l1)
      f1=dr(orbi,ak,orbk,half,zero,-half)
      f2=dr(orbj,ak,orbl,half,zero,-half)
      f3=dsqrt((two*orbi+1)*(two*orbj+1)*(two*orbk+1)*(two*orbl+1))
c     f4=(-one)**(l1+orbi+orbj+1)
      slater_newbr2=f1*f2*f3*slater_breitkm1(ii,jj,kk,ll,orbi,
     :orbj,orbk,orbl,iai,iaj,iak,ial,l1)


      endif
      endif

      return
      end


c*************************************************************************
c    Operator strength   have been  calculated  without phase factor     *
c     for n=k+1 part of the breit  interaction                           *
c                                                                        *
c                                                                        *
c*************************************************************************

      function slater_newbr3(ii,jj,kk,ll,orbi,orbj,orbk,orbl,
     :iai,iaj,iak,ial,l1)
      implicit real*8 (a-h,o-z)
      include 'mpif.h'
      data zero,half,one,two/0.0d0,0.50d0,1.0d0,2.0d0/

      slater_newbr3=0.0
      fact=0.0
c--------------------------------------------------------------------------
c      multipoles (k's) are selected for the odd operators
c--------------------------------------------------------------------------

      fact=select1(iai,iak,orbi,orbk,l1)*
     : select1(iaj,ial,orbj,orbl,l1)
       if(fact.ne.zero)then


      ak=dfloat(l1)
      f1=dr(orbi,ak,orbk,half,zero,-half)
      f2=dr(orbj,ak,orbl,half,zero,-half)
      f3=dsqrt((two*orbi+1)*(two*orbj+1)*(two*orbk+1)*(two*orbl+1))
c     f4=(-one)**(l1+orbi+orbj+1)
      slater_newbr3=f1*f2*f3*slater_breitkp1(ii,jj,kk,ll,orbi,
     :orbj,orbk,orbl,iai,iaj,iak,ial,l1)

      
      endif

      return
      end


c------------------------------------------------------------------------------
c   the integral in the interval zero to r(2) is computed by use of an
c   analytical fit
c                                sigma
c                      f(r) = a r
c   a five-point  closed  newton-cotes  formula (cf. f b hildebrand,
c   introduction to numerical analysis, second edition, mcgraw-hill,
c   new york, 1974, p 93)  is  used  to  compute the integral in the
c   interval  r(2:mtp).
c------------------------------------------------------------------------------

       subroutine quad1 (ta,mtp,result)

      implicit real*8(a-h, o-z)
      include 'mpif.h'
      dimension ta(MN)
      common/cons/zero,half,tenth,one,two,three,ten
     : /grid/r(MN),rp(MN),rpor(MN),rnt,h,hp,n
     : /ncc/c1,c2,c3,c4
c------------------------------------------------------------------------------
c   find first values that will permit computation of exponent
c------------------------------------------------------------------------------
      result = zero
      mtpm1 = mtp - 1
      do 3 i = 2,mtpm1
        tai = ta(i)
        if (dabs (tai) .gt. zero) then
          ip1   = i+1
          taip1 = ta(ip1)
          quott = taip1/tai
          if (quott .gt. zero) then
c------------------------------------------------------------------------------

c                          exponent from fit
c------------------------------------------------------------------------------
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
c------------------------------------------------------------------------------
c              newton-cotes quadature for the remainder
c------------------------------------------------------------------------------
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
c------------------------------------------------------------------------------
c          no value which will permit computation of exponent
c------------------------------------------------------------------------------
      result = zero
    4 return
      end


  

