CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC                                                                    CC
CC               IBRA = 0   for only couloumb interaction             CC
CC                    = 1   for all parity channel with coulomb int.  CC
CC                    = 2   for both coulomb and breit interaction    CC
CC                    = 3   for only breit interaction                CC
CC                                                                    CC
CC                                                                    CC
CC               IOPT = 0   for CCD                                   CC
CC                    = 1   for CCSD                                  CC
CC                                                                    CC
CC             NONLIN = 0   for linear CC                             CC
CC                    = 1   for non-linear CC                         CC
CC                                                                    CC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC                                                                    CC
CC     INPUT FILE NAME:   input.ea                                    CC
CC    READ:                                                           CC
CC   1.  a) No. of symm  b) No. of active orbitals for each symm      CC
CC   2.  a) No. of occupied orbitals for each symm                    CC
CC   3.  a) iopt b) ibra c) nprint d) simopt e) t3opt                 CC
C   4.  a) No. of maximum open-shell will be obtained b) orbital no.sCC
CC   5.  a) open-shell name for each valence orbital                  CC
CC                                                                    CC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


c dimension variables

#define     NHO     23
#define     NPO     82
#define     MXV      9
#define    MDIM   1000000

#define      MN     750
#define   MNSYM      11
#define     MNS       6
#define    MXVR    (MXV+1)/2
#define   MNBAS    NHO+NPO
#define   MNOCC    NHO+NPO
#define   MNEXC    NHO+NPO
#define   NTEMP    NPO*(NPO+1)/2
#define   MULLI    NTEMP*(NTEMP+1)/2
#define   N2ORB    (NHO+NPO)*((NHO+NPO)+1)/2
#define   N2BAS    (NHO+NPO)*(NHO+NPO)


#define    NEV    57                   !Total no. of even parity orbitals^M
#define    NEE    NEV*(NEV+1)/2
C#define          NEEEE=NEE*(NEE+1)/2
#define    NEEEE  1367031 
 
#define    NOD    48                   !Total no. of  odd parity orbitals^M
#define    NOO    NOD*(NOD+1)/2
C#define          NOOOO=NOO*(NOO+1)/2
#define    NOOOO  692076

#define    NOE    NEV*NOD
C#define   NOEOE  NOE*(NOE+1)/2
#define    NOEOE  3744216

c#define   N2INT   2*(NOOOO+NEEEE+NOEOE)
#define   N2INT   (NOOOO+NEEEE+NOEOE+2000000)


c constants
#define NFVBAR    3
#define MFVBAR    4
#define  STDIN    5
#define STDOUT    6
#define  WFNIN    7
#define  NTWOE    8
#define NTFILE   16
#define NSFILE   17
#define LINEAR   15
#define   NBFL   12
#define SIMPOT    1
#define NPRINT    1
#define   ITOL    6
#define NITMAX   20
#define  ITMAX   30
#define  NPMAX   32
#define  SMALL  1.0d-5
#define  NVAL    50
#define  MASTER  0
      
      implicit real*8(a-h,o-z)
      include 'mpif.h'

      integer drow,rowstart,rowend,extrarow,row
      integer error,dspl,delrow
      integer simopt,t3opt
      integer valorb(NVAL)            ! List of valence orbitals
      character*14 sfile(NVAL)        ! List of so- file names
      character*14 soname             ! Present name for so- file
c     local dimension 
      dimension nbas(MNSYM),iocc(MNSYM),t(MDIM)
      dimension dspl(NPMAX),delrow(NPMAX)

      common/scratch1/b(MDIM),ad(MDIM)
      common/scratch2/raux(MDIM),caux(MDIM)
      common/iallparity/ibra
      common/skip3/ncis
      dimension s(MDIM),sold(MDIM)


c     nproc=1
c     iproc=0
c     error=0


      call MPI_INIT(ierr)
      START_TIME=MPI_WTIME()
      call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)

      call daopen
      call symmcal
      call setgrd

      call readinp(nsym,nbas,iocc,maxvalence,valorb,
     &  iopt,itwoe,ivbar,nprint,simopt,t3opt,sfile)

      call setup(nsym,nbasis,nbas,iocc,iproc)
      call grasprd(nbasis,iproc)
      call setqic
      call symm(iopt,nprint,ntmax,iproc)

      if(iproc.eq.0)then

      write(STDOUT,*)'<<<<<<<<<<<<PLEASE CHECK>>>>>>>>>>>>>>'
      write(STDOUT,*)'<<*************WHETHER**************>>'
      write(STDOUT,*)'<<************YOU WANT THIS*********>>'
      write(STDOUT,*)'<<************FOR OPEN-SHELL********>>'
      write(STDOUT,*)'<</                                 />>'
      if(ibra.eq.0)write(STDOUT,*)'<< Only Coulomb interaction       >>'
      if(ibra.eq.1)write(STDOUT,*)'<< All parity channel with Coulomb>>'
      if(ibra.eq.2)write(STDOUT,*)'<< Both Coulomb and Breit         >>'
      if(ibra.eq.3)write(STDOUT,*)'<< Only Breit interaction         >>'
      if(iopt.eq.0)write(STDOUT,*)'<< CCD  is considered             >>'
      if(iopt.eq.1)write(STDOUT,*)'<< CCSD is considered             >>'
      if(t3opt.eq.0)write(STDOUT,*)'<< Partial triples neglected     >>'
      if(t3opt.eq.1)write(STDOUT,*)'<< Partial triples considered    >>'
      write(STDOUT,*)'<</                                 />>'
      write(STDOUT,*)'<<<<<<<<<<<<OOOOO>>>>>>>>>>>>>>>>>>>>>'
      write(STDOUT,*)'____________________________________________'
      write(STDOUT,*)'____________________________________________'
      write(STDOUT,*)'simopt value = ',simopt,' nprint value = ',
     : nprint
      write(STDOUT,*)'Total no. excitations from Closed-shell = ',ntmax
      write(STDOUT,*)'____________________________________________'
      write(STDOUT,*)'____________________________________________'

      endif

      call mInit()
c     begin one and two-electron part
      call twoint(nbasis,itwoe)

      call fbar(t,ivbar)
      call vbar(t,ivbar)

c     end of one and two-electron generation

c---------------------------------------------------------------------
c     Start the valence-orbital loop
c---------------------------------------------------------------------
      do 10 ivalence = 1, maxvalence 
         do i=1, MDIM
            s(i)=0.0
            b(i)=0.0
            ad(i)=0.0
            raux(i)=0.0
            caux(i)=0.0
         enddo

         iorb=valorb(ivalence)

         soname=sfile(ivalence)

         call ealoc(iorb,nsmax,nprint,simopt,iopt,iproc)

       if(iproc.eq.0)then
      write(STDOUT,*)'____________________________________________'
      write(STDOUT,*)'____________________________________________'
      write(STDOUT,*)'Total no. single excitations in open-shell = ',
     : ncis
      write(STDOUT,*)'Total no. double excitations in open-shell = ',
     : nsmax-ncis
      write(STDOUT,*)'____________________________________________'
      write(STDOUT,*)'Total no. excitations in open-shell =        ',
     : nsmax
      write(STDOUT,*)'____________________________________________'
       endif
    
         call eamata(iorb,nsmax,simopt,iopt,iproc)

c        initial guess
         do i =1,nsmax
            s(i)=b(i)/ad(i)

            if(dabs(ad(i)) .lt. 1.0d-30)then
            if(iproc.eq.0)write(STDOUT,*)'Boooiiiinnnng' 
            stop
            endif

            sold(i)=s(i)
         enddo

         drow = nsmax/nproc
         do i = 1, nproc
         delrow(i) = drow
         end do
         extrarow = mod(nsmax,nproc)
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

         simopt = 0              ! Only Jacobi is enabled


         call jacobi(itol,nprint,rowstart,rowend,simopt,
     :   nsmax,iopt,iproc,iorb,t3opt,t,dspl,delrow,drow,sold,s,
     &   soname)

 10   continue

c---------------------------------------------------------------------
c     End of valence orbital loop
c---------------------------------------------------------------------

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

c     call MPI_FINALIZE(error)

9993  format(/,1x,'LINEAR RUN BEGINS',/)
9994  format(/,1x,'SOLUTION VECTORS',/)
9995  format(4f20.11)
9996  format(/,1x,'NON-LINEAR RUN BEGINS',/)
9997  format(/,1x,'ITERATION NO=',i4,5x,'EPS=',f20.11,/)
9998  format(/,1x,'CONVERGENCE FAILURE; PROGRAMM TERMINATES') 
9999  format(/,1x,'ITERATION CONVERGES AT ',I4,' CYCLES')
      stop
      end

c***********************************************************************
c        SUBROUTINE  PUTVBAR                                           *
c                                                                      *
c        This subroutine put the Vbar value in a particular location   *
c***********************************************************************

      subroutine PutVbar(i,j,m,val)
      implicit real*8 (a-h,o-z)
      include 'mpif.h'
      common/iallparity/iall
      common/v_bar/
     :      vb(N2ORB,N2BAS,0:MXV)
      common/mind/mIndex(0:20)      

      if(iall.eq.0)then
      m1 = mIndex(m)
      vb(i,j,m1) = val
      else
      vb(i,j,m) = val
      endif

      return
      end

c*********************************************************************
c       FUNCTION GETVBAR                                             *
c                                                                    *
c   This function is used to extract the vbar values in 'B' matrixes *
c   from the subroutine PutVbar                                      *
c*********************************************************************

      real*8 function GetVbar(i,j,m)
      implicit real*8 (a-h,o-z)
      include 'mpif.h'
      common/iallparity/iall
      common/v_bar/
     :       vb(N2ORB,N2BAS,0:MXV)
      common/mind/mIndex(0:20)
      
      if(iall.eq.0)then
      m1 = mIndex(m)
      GetVbar = vb(i,j,m1)
      else
      GetVbar = vb(i,j,m)
      endif

      return
      end

c****************************************************************************
c                       SUBROUTINE  EAMATB                                  *
c                                                                           *
c          This subroutine calculate the 'B' matrixes of Open shell CC      *
c          equations                                                        *
c                                                                           *
c****************************************************************************


      subroutine eamatb(iopt,iorb,nsmax,isolv,rows,rowe,s,t)

      implicit real*8 (a-h,o-z)
      include 'mpif.h'
      integer row,rows,rowe
c     real*8 raux_local,bpar_local
      common/jvalue/orbc(MNOCC),orbe(MNEXC)
      common/index/ke(MNEXC),kc(MNOCC)
      common/info/nocc,nexcit,iType(MNBAS),iTCount(MNBAS)
      common/parity/isymc(MNOCC),isyme(MNEXC)
      common/symmetry/mtbl(MNS,MNS)
      common/kpavlaue/iqc(MNOCC),iqe(MNEXC)
      common/skip1/npbps(0:MXV,2)
      common/skip2/mmskip(0:MXV,2),m2skip
      common/skip3/ncis
      common/idra1/idra(MNBAS,MNBAS,0:MXV)
      common/idpa1/idpa(MNBAS,MNBAS)
      common/jdra1/jdra(MNBAS,MNBAS,0:MXV)
      common/jdpa1/jdpa(MNBAS,MNBAS)
      common/f_bar/f(MNBAS,MNBAS)
      common/scratch2/raux(MDIM),caux(MDIM)
      dimension bpar(MDIM),kgot(MXV+1),s(MDIM),t(MDIM)
      dimension igot(MXV+1),jgot(MXV+1)
c    :,bpar_local(MDIM),raux_local(MDIM)


      data zero,half,one,two/0.0,0.5,1.0,2.0/
      real*8 GetVbar      
c     call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
c     call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)


      do i=0,MXV
      kgot(i)=0
      jgot(i)=0
      igot(i)=0
      enddo


      do i = 1, MDIM
      bpar(i) = 0.0
c     bpar_local(i) = 0.0
c     raux_local(i) = 0.0
      enddo

C---------------S1-S1 and S1-S2 Block-----------------------------------

      if(iopt.ne.0)then

c     do 1 ip=iproc+1,nexcit,nproc
      do 1 ip=1,nexcit
      ipsym=isyme(ip)
      iip=ke(ip)
      orbp=orbe(ip)
      iap=iqe(ip)

      do 1 iq=1,nexcit
      iqsym=isyme(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      iiq=ke(iq)

      if(iip.ne.iorb)go to 1
      if(iiq.eq.iip)go to 1
      if(orbp.ne.orbq)go to 1
      if(ipsym.ne.iqsym)go to 1
      neqn=jdpa(iip,iiq)
                                                                                
      row=neqn
      irow=row+1-rows
      if(irow.le.0)go to 1
      if(irow.gt.rowe-rows+1)go to 1
                                                                                
c     call init(bpar_local,nsmax)
      call init(bpar,nsmax)
                                                                                
      if(isolv.ne.0)then
c     bpar_local(neqn)=bpar_local(neqn)-f(iip,iip)
      bpar(neqn)=bpar(neqn)-f(iip,iip)
      endif
                                                                                
c------------diagram 1     <q|f|r><r|s|p>-----------------------------------

      do 7 ir=1,nexcit
      irsym=isyme(ir)
      iir=ke(ir)
      iar=iqe(ir)
      orbr=orbe(ir)

      if(irsym.ne.iqsym)go to 2
      if(orbr.ne.orbq)go to 2
      if(iir.eq.iip)go to 2
      in=jdpa(iip,iir)

      if(isolv.ne.0)then
c     bpar_local(in)=bpar_local(in)+f(iiq,iir)
      bpar(in)=bpar(in)+f(iiq,iir)
      else
c     if(iiq.ne.iir)bpar_local(in)=bpar_local(in)+f(iiq,iir)
      if(iiq.ne.iir)bpar(in)=bpar(in)+f(iiq,iir)
      endif

  2   continue
                                                                                
                                                                                
      do 7 ig=1,nocc
      igsym=isymc(ig)
      iig=kc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)

      if(orbr.ne.orbg)go to 4
      if(irsym.ne.igsym)go to 4

c--------------------- diagram 2 <g|f|r><qr|s|pg>-----------------------------
                                                                                
      call findk(orbp,orbq,orbr,orbg,iap,iaq,iar,iag,igot,
     $imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      ipq=jdra(iip,iiq,l1)
      irg=idra(iir,iig,l1)
      ipqsym=mtbl(ipsym,iqsym)
      in=nsloc(ipq,irg,ipqsym,l1)
      amult=dsqrt((two*orbg+one)/(two*orbp+one))
c     bpar_local(in)=bpar_local(in)+f(iig,iir)*amult
      bpar(in)=bpar(in)+f(iig,iir)*amult

      endif
      enddo
      endif
                                                                                
c---------------- diagram  3  <g|f|r><pq|s2|rg>-------------------------

      call findk(orbp,orbr,orbg,orbq,iap,iar,iag,iaq,igot,
     $imax)
      if(imax.eq.0)go to 4
      do 5 mloop=1,imax
      ll=igot(mloop)
      al=float(ll)
      ipr=jdra(iip,iir,ll)
      iqg=idra(iiq,iig,ll)
      iprsym=mtbl(ipsym,irsym)
      in=nsloc(ipr,iqg,iprsym,ll)
      amult=one/(two*orbp+one)
      isign=(-1)**(orbq+orbr+ll)
c     bpar_local(in)=bpar_local(in)+f(iig,iir)*amult*isign
      bpar(in)=bpar(in)+f(iig,iir)*amult*isign

   5  continue
   4  continue
   3  continue
                                                                                
      iqgsym=mtbl(iqsym,igsym)
      iqgrsm=mtbl(iqgsym,irsym)

      do 8 is=1,nexcit
      issym=isyme(is)
      iis=ke(is)
      ias=iqe(is)
      orbs=orbe(is)

      if(issym.ne.iqgrsm)go to 8

c----------------------diagram  4   <qs|v|rg><ps|s2|rg>------------------------

      call findk(orbq,orbr,orbg,orbs,iaq,iar,iag,ias,igot,
     $imax)
      if(imax.eq.0)go to 8
      do 9 iloop=1,imax
      kk=igot(iloop)
      ak=float(kk)
      iqg=idis(iiq,iig)
      irs=indxr(iir)+iis
      if(iig.gt.iiq)irs=indxr(iis)+iir
      rk=GetVbar(iqg,irs,kk)
      f1=(-1)**(orbr+orbg+orbq+orbs)
      f2=one/((two*orbp+one)*(two*ak+one))
      npr=jdra(iip,iir,kk)
      nsg=idra(iis,iig,kk)
      iprsym=mtbl(ipsym,irsym)
      in=nsloc(npr,nsg,iprsym,kk)
c     bpar_local(in)=bpar_local(in)+rk*f1*f2
      bpar(in)=bpar(in)+rk*f1*f2
                                                                                
c-------------------diagram 5     <rs|v|qg><rp|s2|gs>---------------------------

      call findk(orbp,orbs,orbg,orbr,iap,ias,iag,iar,jgot,
     $jmax)
      if(jmax.eq.0)go to 9
      do 10 jloop=1,jmax
      ll=jgot(jloop)
      al=float(ll)
      nps=jdra(iip,iis,ll)
      nrg=idra(iir,iig,ll)
      ipssym=mtbl(ipsym,issym)
      in=nsloc(nps,nrg,ipssym,ll)
      call sixj(ak,orbs,orbg,al,orbr,orbq,s6j)
      f2=s6j*one/(two*orbp+1)
c     bpar_local(in)=bpar_local(in)+rk*f1*f2
      bpar(in)=bpar(in)+rk*f1*f2

 10   continue
  9   continue
  8   continue
  7   continue

c----------------------- Sum of S1-S1  and S1-S2  block-------------------------
                                                                                
      xsum=0.0d0
      do in=1,nsmax
c     xsum=xsum+bpar_local(in)*s(in)
      xsum=xsum+bpar(in)*s(in)
      enddo
                                                                                
c     raux_local(neqn)=xsum
      raux(neqn)=xsum
                                                                                
  1   continue
c            call MPI_ALLREDUCE(raux_local,raux,MDIM,MPI_REAL8,MPI_SUM,
c    :                   MPI_COMM_WORLD,ierr)

c     do i= 1,MDIM
c      raux_local(i)=0.0
c     enddo

c            call MPI_ALLREDUCE(bpar_local,bpar,MDIM,MPI_REAL8,MPI_SUM,
c    :                   MPI_COMM_WORLD,ierr)

c     do i= 1,MDIM
c      bpar_local(i)=0.0
c     enddo


      endif
     
c-------------------- S2-S2 and S2-S1 Block-------------------------------------

c     do 15 ip=iproc+1,nexcit,nproc
      do 15 ip=1,nexcit
      ipsym=isyme(ip)

      do 20 iq=1,nexcit
      iqsym=isyme(iq)

      ipqsym=mtbl(ipsym,iqsym)

      do 25 ir=1,nexcit
      irsym=isyme(ir)

      ipqrsm=mtbl(ipqsym,irsym)

      do 30 ig=1,nocc
      igsym=isymc(ig)

      irgsym=mtbl(irsym,igsym)

      if(ipqrsm.ne.igsym)go to 30

      orbp=orbe(ip)
      orbq=orbe(iq)
      orbr=orbe(ir)
      orbg=orbc(ig)

      iap=iqe(ip)
      iaq=iqe(iq)
      iar=iqe(ir)
      iag=iqc(ig)

      iip=ke(ip)
      iiq=ke(iq)
      iir=ke(ir)
      iig=kc(ig)

      if(iip.ne.iorb)go to 30

      call findk(orbp,orbq,orbr,orbg,iap,iaq,iar,iag,kgot,
     $kmax)
      if(kmax.eq.0)go to 30
      do 35 nloop=1,kmax
      l1=kgot(nloop)
      al1=float(l1)
      npq=jdra(iip,iiq,l1)
      nrg=idra(iir,iig,l1)
      ipqsym=mtbl(ipsym,iqsym)
      neqn=nsloc(npq,nrg,ipqsym,l1)
                                                                                
      row=neqn
      irow=row+1-rows
      if(irow.le.0)go to 35
      if(irow.gt.rowe-rows+1)go to 35
                                                                                
c     call init(bpar_local,nsmax)
      call init(bpar,nsmax)
      if(isolv.ne.0)then
c     bpar_local(neqn)=bpar_local(neqn)-f(iip,iip)
      bpar(neqn)=bpar(neqn)-f(iip,iip)
      endif
                                                                                
c------------------------- diagram 8      <dp|v|gq><rp|s2|dq>--------------------

      do 40 id=1,nocc
      idsym=isymc(id)
      iid=kc(id)
      orbd=orbc(id)
      iad=iqc(id)

      if(orbd.ne.orbg)go to 40
      if(idsym.ne.igsym)go to 40

      npq=jdra(iip,iiq,l1)
      nrd=idra(iir,iid,l1)
      ipqsym=mtbl(ipsym,iqsym)
      in=nsloc(npq,nrd,ipqsym,l1)
                                                                                
      if(isolv.ne.0)then
c     bpar_local(in)=bpar_local(in)-f(iid,iig)
      bpar(in)=bpar(in)-f(iid,iig)
      else
c     if(iid.ne.iig)bpar_local(in)=bpar_local(in)-f(iid,iig)
      if(iid.ne.iig)bpar(in)=bpar(in)-f(iid,iig)
      endif
                                                                                
   40 continue
                                                                                
      do 110 is=1,nexcit
      issym=isyme(is)
      iis=ke(is)
      orbs=orbe(is)
      ias=iqe(is)
                                                                                
      if(orbr.eq.orbs) then
      if(irsym.eq.issym)then

      npq=jdra(iip,iiq,l1)
      nsg=idra(iis,iig,l1)
      ipqsym=mtbl(ipsym,iqsym)
      in=nsloc(npq,nsg,ipqsym,l1)

c------------------------ diagram 6   <ps|v|qr><ps|s2|qg>-----------------------
                                                                                
      if(isolv.ne.0)then
c     bpar_local(in)=bpar_local(in)+f(iir,iis)
      bpar(in)=bpar(in)+f(iir,iis)
      else
c     if(iir.ne.iis)bpar_local(in)=bpar_local(in)+f(iir,iis)
      if(iir.ne.iis)bpar(in)=bpar(in)+f(iir,iis)
      endif
                                                                                
      endif
      endif
                                                                                
      if(iqsym.eq.issym)then
      if(orbs.eq.orbq) then

      nps=jdra(iip,iis,l1)
      nrg=idra(iir,iig,l1)
      ipssym=mtbl(ipsym,issym)
      in=nsloc(nps,nrg,ipssym,l1)

c-------------------------diagram   7   <qr|v|sg><pr|s2|sg>----------------------
                                                                                
      if(isolv.ne.0)then
c     bpar_local(in)=bpar_local(in)+f(iiq,iis)
      bpar(in)=bpar(in)+f(iiq,iis)
      else
c     if(iiq.ne.iis)bpar_local(in)=bpar_local(in)+f(iiq,iis)
      if(iiq.ne.iis)bpar(in)=bpar(in)+f(iiq,iis)
      endif
                                                                                
      endif
      endif
                                                                                
   45 continue
                                                                                                                                                               
      if(iis.eq.iorb)go to 60
      if(issym.ne.ipsym) go to 60
      if(orbs.ne.orbp) go to 60

c------------------------diagram   9    <qr|v|sg><p|s1|s>------------------------

      ias=iqe(is)
      in=jdpa(iip,iis)
      iqr=idis(iiq,iir)
      isg=indxr(iis)+iig
      if(iir.gt.iiq)isg=indxr(iig)+iis
      rk=GetVbar(iqr,isg,l1)
c     bpar_local(in)=bpar_local(in)+rk
      bpar(in)=bpar(in)+rk

   60 continue
                                                                                
                                                                                
      do 120 id=1,nocc
      idsym=isymc(id)
      iid=kc(id)
      orbd=orbc(id)
      iad=iqc(id)
                                                                                
      igdsym=mtbl(igsym,idsym)
      irgdsm=mtbl(irsym,igdsym)
                                                                                
      if(issym.eq.irgdsm)then

c--------------------- diagram  10 <rs|v|gd><ps|s2|qd>---------------------------
                                                                                
      ird=idis(iir,iid)
      igs=indxr(iig)+iis
      if(iid.gt.iir)igs=indxr(iis)+iig
      rk=GetVbar(ird,igs,l1)
      npq=jdra(iip,iiq,l1)
      nsd=idra(iis,iid,l1)
      isdsym=mtbl(issym,idsym)                                                                                
      in=nsloc(npq,nsd,isdsym,l1)
      f1=(-1)**(orbs+orbd+l1)
      f2=1.0/(2*l1+1)
c     bpar_local(in)=bpar_local(in)-rk*f1*f2
      bpar(in)=bpar(in)-rk*f1*f2
                                                                                
c-------------------- diagram   11  <rs|v|gd><pq|s2|sd>--------------------------

      call findk(orbd,orbq,orbs,orbp,iad,iaq,ias,iap,igot,imax)
      if(imax.ne.0)then
      do  iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      nps=jdra(iip,iis,l2)
      nqd=idra(iiq,iid,l2)
      iqdsym=mtbl(iqsym,idsym)
      in=nsloc(nps,nqd,iqdsym,l2)
      call sixj(al1,orbq,orbp,al2,orbs,orbd,s6j)
      f2=s6j
c     bpar_local(in)=bpar_local(in)-rk*f1*f2
      bpar(in)=bpar(in)-rk*f1*f2

      enddo
      endif
                                                                                
c------------------- diagram   13    <rd|v|sg><ps|s2|qd>-------------------------

      call findk(orbd,orbg,orbr,orbs,iad,iag,iar,ias,igot,
     $imax)
      if(imax.ne.0)then
      do mloop=1,imax
      l2=igot(mloop)
      al2=float(l2)
      ird=idis(iir,iid)
      isg=indxr(iis)+iig
      if(iid.gt.iir)isg=indxr(iig)+iis
      rk=GetVbar(ird,isg,l2)
      npq=jdra(iip,iiq,l1)
      nsd=idra(iis,iid,l1)
      isdsym=mtbl(issym,idsym)
      in=nsloc(npq,nsd,isdsym,l1)
      call sixj(al1,orbd,orbs,al2,orbr,orbg,s6j)
      f1=s6j*(-1)**(orbd+orbs+l1)
c     bpar_local(in)=bpar_local(in)-rk*f1
      bpar(in)=bpar(in)-rk*f1

      enddo
      endif
                                                                                
      endif
                                                                                
      igdsym=mtbl(igsym,idsym)
      iqgdsm=mtbl(iqsym,igdsym)
                                                                                
      if(issym.eq.iqgdsm)then
                                                                                
c------------------ diagram 14   <qd|v|sg><sr|s2|pd>-----------------------------

      call findk(orbd,orbg,orbq,orbs,iad,iag,iaq,ias,igot,
     $imax)
      if(imax.ne.0)then
      do mloop=1,imax
      l2=igot(mloop)
      al2=float(l2)
      iqd=idis(iiq,iid)
      isg=indxr(iis)+iig
      if(iid.gt.iiq)isg=indxr(iig)+iis
      rk=GetVbar(iqd,isg,l2)

      call findk(orbp,orbs,orbd,orbr,iap,ias,iad,iar,jgot,
     $jmax)
      if(jmax.ne.0)then
      do iloop=1,jmax
      l3=jgot(iloop)
      al3=float(l3)
      nps=jdra(iip,iis,l3)
      nrd=idra(iir,iid,l3)
      ipssym=mtbl(ipsym,issym)
      in=nsloc(nps,nrd,ipssym,l3)

      call sixj(al1,al2,al3,orbd,orbr,orbg,s6j1)
      call sixj(al1,al3,al2,orbs,orbq,orbp,s6j2)
      f1=(two*al1+1)*s6j1*s6j2
      f2=(-1)**(orbp+orbq+orbg+orbr+l1+l2+l3)
c     bpar_local(in)=bpar_local(in)-rk*f1*f2
      bpar(in)=bpar(in)-rk*f1*f2

      enddo
      endif
      enddo
      endif
                                                                                
      endif
                                                                                
 120  continue
                                                                                
c---------------------- diagram 12   <qt|v|sr><pt|s2|sg>------------------------

      do 130 it=1,nexcit
      itsym=isyme(it)
      iit=ke(it)
      iat=iqe(it)
      orbt=orbe(it)
                                                                                
      isqsym=mtbl(issym,iqsym)
      isqrsm=mtbl(isqsym,irsym)
      if(itsym.eq.isqrsm)then
                                                                                
      call findk(orbt,orbr,orbs,orbq,iat,iar,ias,iaq,igot,
     $imax)
      if(imax.ne.0)then
      do mloop=1,imax
      l2=igot(mloop)
      al2=float(l2)
      iqr=idis(iiq,iir)
      ist=indxr(iis)+iit
      if(iir.gt.iiq)ist=indxr(iit)+iis
      rk=GetVbar(iqr,ist,l2)

      call findk(orbp,orbs,orbt,orbg,iap,ias,iat,iag,jgot,
     $jmax)
      if(jmax.ne.0)then
      do iloop=1,jmax
      l3=jgot(iloop)
      al3=float(l3)
      nps=jdra(iip,iis,l3)
      ntg=idra(iit,iig,l3)
      ipssym=mtbl(ipsym,issym)
      in=nsloc(nps,ntg,ipssym,l3)
      call sixj(al1,al3,al2,orbs,orbq,orbp,s6j1)
      call sixj(al1,al3,al2,orbt,orbr,orbg,s6j2)
      f1=s6j1*s6j2*(two*al1+1)
      f2=(-1)**(orbp+orbq+orbg+orbr)
c     bpar_local(in)=bpar_local(in)+rk*f1*f2
      bpar(in)=bpar(in)+rk*f1*f2

      enddo
      endif
      enddo
      endif
                                                                                
      endif
c----------------three body vt2 contribution to S equation----------------------

      do 170 ia=1,nocc
      iia=kc(ia)
      iaa=iqc(ia)
      orba=orbc(ia)
      iasym=isymc(ia)

      do 170 ib=1,nocc
      iib=kc(ib)
      iab=iqc(ib)
      orbb=orbc(ib)
      ibsym=isymc(ib)

      iabdis=idis(iia,iib)
      iabsym=mtbl(iasym,ibsym)
      istsym=mtbl(issym,itsym)
      iptsym=mtbl(ipsym,itsym)
      if(istsym.ne.iabsym)go to 170
      itsloc=indxr(iit)+iis
      if(iib.gt.iia)itsloc=indxr(iis)+iit
      if(orba.ne.orbp)go to 170
      if(iasym.ne.ipsym)go to 170
      nrg=idra(iir,iig,l1)
      nqa=idra(iiq,iia,l1)
      irgsym=mtbl(irsym,igsym)
      in1=ntloc(nrg,nqa,irgsym,l1)

      call findk(orbs,orbb,orbt,orbp,ias,iab,iat,iap,igot,imax)
      if(imax.eq.0)go to 170
      do 180 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      npt=jdra(iip,iit,l2)
      nsb=idra(iis,iib,l2)
      in=nsloc(npt,nsb,iptsym,l2)
      rk=vint(iia,iib,iit,iis,l2)
      d4=(-1)**(orba+orbb+l2+1)
      twoe=rk*d4

      fact=one/((two*orbp+1)*(two*al2+1))
      isign=(-1)**(orbb+orba+orbs+orbt)
c     bpar_local(in)=bpar_local(in)-twoe*fact*isign*t(in1)
      bpar(in)=bpar(in)-twoe*fact*isign*t(in1)

      call findk(orbs,orbp,orbb,orbt,ias,iap,iab,iat,jgot,jmax)
      if(jmax.eq.0)go to 180
      do 190 jloop=1,jmax
      l3=jgot(jloop)
      al3=float(l3)
      nps=jdra(iip,iis,l3)
      ntb=idra(iit,iib,l3)
      ipssym=mtbl(ipsym,issym)
      in=nsloc(nps,ntb,ipssym,l3)
      call sixj(al2,orbp,orbt,al3,orbb,orbs,s6j)
      isign=(-1)**(orbb+orba+orbs+orbt+2*orbt)
      fact=s6j/(two*orbp+1)
c     bpar_local(in)=bpar_local(in)+twoe*fact*isign*t(in1)
      bpar(in)=bpar(in)+twoe*fact*isign*t(in1)

 190  continue
 180  continue
 170  continue
                                                                                
 999  continue
                                                                                
 130  continue
 110  continue

c------------  Calculate sum of S2-S1 and S2-S2 block--------------------------
                                                                                
      xsum=0.0d0
      do in=1,nsmax
c     xsum=xsum+bpar_local(in)*s(in)
      xsum=xsum+bpar(in)*s(in)
      enddo
c     raux_local(neqn)=xsum
      raux(neqn)=xsum
                                                                                
  35  continue
  30  continue
  25  continue
  20  continue
  15  continue
c            call MPI_ALLREDUCE(raux_local,raux,MDIM,MPI_REAL8,MPI_SUM,
c    :                   MPI_COMM_WORLD,ierr)

c     do i= 1,MDIM
c      raux_local(i)=0.0
c     enddo

c            call MPI_ALLREDUCE(bpar_local,bpar,MDIM,MPI_REAL8,MPI_SUM,
c    :                   MPI_COMM_WORLD,ierr)

c     do i= 1,MDIM
c      bpar_local(i)=0.0
c     enddo


      return
      end

c********************************************************************
c                                                                   *
c       INTEGER FUNCTION IDELTA                                     *
c                                                                   *
c       Select multipoles according to traingular condition         *
c                                                                   *
c********************************************************************

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

c************************************************************************
c                                                                       *
c     SUBROUTINE FINDK                                                  *
c                                                                       *
c    Select the multipoles either for all parity channel or for even    *
c    parity channel                                                     *
c                                                                       *
c************************************************************************


      subroutine findk(orba,orbb,orbp,orbq,ia,ib,ip,iq,jeven,kkk)

      implicit real*8 (a-h,o-z)
      include 'mpif.h'

      common/iallparity/iall

      dimension jeven(MXV+1)

      if(iall.eq.0)then

c-----------------------select even parity channel--------------------------


      call findk2(orba,orbb,orbp,orbq,ia,ib,ip,iq,
     $jeven,kkk)

      else

c------------------------ select all parity channel-----------------------


      call findk_all(orba,orbb,orbp,orbq,ia,ib,ip,iq,jeven,kkk)

      endif
  
      return
      end


c**********************************************************************
c                                                                     *
c      SUBROUTINE   FINDK_ALL                                         *
c                                                                     *
c      This subroutine select the multipoles according to traingular  *
c      condition                                                      *
c                                                                     *
c**********************************************************************


      subroutine findk_all(orba,orbb,orbp,orbq,ia,ib,ip,iq,jeven,kkk)

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


c***************************************************************************
c                                                                          *
c                      SUBROUTINE    FINDK2                                *
c                                                                          *
c             This subroutine select the multipole for even operator       *
c                                                                          *
c***************************************************************************


      subroutine findk2(orba,orbb,orbp,orbq,ia,ib,ip,iq,
     $jeven,kkk)
c     this subroutine finds the common k value
      implicit real*8 (a-h,o-z)
      include 'mpif.h'
      real*8 jj
      dimension ja(MXV),jb(MXV),
     :          jp(MXV),jq(MXV),jeven(MXV)

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
      jj=float(jtot)-2.0*(float(jtot/2))
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
 
      subroutine init(a,n)
c     initialization
      real*8 a(MDIM)
      do i=1,n
      a(i)=0
      enddo
      return
      end


c*************************************************************************
c                           SUBROUTINE  ILOC                             *  
c                                                                        *
c      Pack two index to one index                                       *
c*************************************************************************

      subroutine iloc(i,j,k,l,ij,kl)
c     pack two index to one index
      implicit integer (a-z)
      include 'mpif.h'
      ij=j*1000+i
      kl=l*1000+k
      return
      end

        subroutine sixj(rj1,rj2,rj3,rl1,rl2,rl3,s6j)
        implicit real*8 (a-h,o-z)
      include 'mpif.h'
        common/factor/fct(0:50),mfd
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
       print *,'inconsistent arguments to 6j symbol'
       print *,tria,trib,tric,trid
       print *,rj1,rj2,rj3
       print *,rl1,rl2,rl3
       stop
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
         ws6j=ws6j+((-1)**iz)*fct(int(z+1.))/(fct(z-rj1-rj2-rj3)*
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
        common/factor/fct(0:50),mfd

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
      
* ------------------------------------------------------------------
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
      double precision function fact(cj1,cj2)
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
      function slater(ia,ib,ic,id,k)

*   the value of this  function is the  slater integral  as normally   *
*   defined in terms of the four sets of quantum numbers  a,b,c,d.     *
*                                                                      *
*   subroutines called: quad, yzk.                                     *
*                                                                      *
************************************************************************
*
      implicit doubleprecision (a-h, o-z)
      include 'mpif.h'

      common/cons/zero,half,tenth,one,two,three,ten
     :/grid/r(MN),rp(MN),rpor(MN),rnt,h,hp,n
     :/tatb/ta(MN),tb(MN),mtp
     :/wave/pz(MNBAS),pf(MN,MNBAS),qf(MN,MNBAS),
     :mf(MNBAS)

      call yzk (k,ib,id)
      mtp=n
      ta(1) = zero
      do 1 i = 2,mtp
      ta(i) =(pf(i,ia)*pf(i,ic)+qf(i,ia)*qf(i,ic))
     :       *rpor(i)*tb(i)
    1 continue
      call quad (result)
      slater = result
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
      slater_new=f1*f2*f3*slater(ii,jj,kk,ll,l1)


      endif

      return
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
*   introduction to numerical analysis,  edition, mcgraw-hill,   *
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
      implicit doubleprecision (a-h, o-z)
      include 'mpif.h'

      common/cons/zero,half,tenth,one,two,three,ten
     :      /def4/accy,nscf,nsic,nsolv
     :      /grid/r(MN),rp(MN),rpor(MN),rnt,h,hp,n
     :      /tatb/ta(MN),tb(MN),mtp
     :      /ncc/c1,c2,c3,c4

*   find first values that will permit computation of exponent
      mtpm1 = mtp-1
      do 3 i = 2,mtpm1
         tai = ta(i)
         if (abs (tai) .gt. zero) then
            ip1   = i+1
            taip1 = ta(ip1)
            quott = taip1/tai
            if (quott .gt. zero) then
*   exponent from fit
               frip1 = taip1/rp(ip1)
               fri   = tai  /rp(i  )
               ratio = frip1/fri
               rip1  = r (ip1)
               ri    = r (i  )
               sigma = log (ratio)/log (rip1/ri)
*   analytical integration and error estimate for interval r(1:i)
               fri    = ri*fri
               result = fri/(sigma+one)
*   set the tail to zero
               do 1 loc = 1,3
                  ta(mtp+loc) = zero
    1          continue
*   newton-cotes quadature for the remainder
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
      result = zero
    4 return
      end

************************************************************************
*                                                                      *
      subroutine yzk(k,i,j)
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
      implicit doubleprecision (a-h, o-z)
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
*   compute rp(s)*(p (s)*p (s)+q (s)*q (s)) and store in rhop
*                   i     j     i     j
      do 2 ii = 2,mtp
      rhop(ii) = rp(ii)*(pf(ii,i)*pf(ii,j)+qf(ii,i)*qf(ii,j))
    2 continue
*   fill array temp with r**k * rhop
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
*   set an additional four points to zero
      do 5 ii = mtpp1,mtpp4
         temp(ii) = zero
    5 continue
*                                     k
*   compute the first few values of  r  * zk  using semi-open
*   newton-cotes formulae
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
      do 8 ii = 5,mtpp3
         zk(ii) = zk(ii-4)+c1*(temp(ii-4)+temp(ii  ))
     :                    +c2*(temp(ii-3)+temp(ii-1))
     :                    +c3* temp(ii-2)
    8 continue
*                                       k   (k)
*   determine the asymptotic value of  r * z
*                   (0)
*   correction to  z   : in the manner of  c froese fischer,
*   the hartree-fock method for atoms, john wiley & sons,
*   new york, 1977, p 235.
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
*   tabulate  zk  for entire internal grid
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
*   start array wk / r**(k+1)
      np4 = n+4
      do 14 ii = np4,mtpp1,-1
         wk(ii) = zero
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
*   compute remainder of wk / r**(k+1): march in to the origin
      do 18 ii = mtp,2,-1
         wk(ii) = wk(ii+4)+c1*(temp(ii  )+temp(ii+4))
     :                    +c2*(temp(ii+1)+temp(ii+3))
     :                    +c3*(temp(ii+2))
   18 continue
      wk(1) = zero
*   compute wk
      if (k .eq. 0) then
         do 19 ii = 2,mtp
            wk(ii) = wk(ii)*r(ii)
   19    continue
      else
         do 20 ii = 2,mtp
            wk(ii) = wk(ii)*rttk(ii)
   20    continue
      endif
*   assemble solution
      yk(1) = zero
      do 21 ii = 2,n
         yk(ii) = zk(ii)+wk(ii)
   21 continue
      return
      end

      subroutine ninej(rj11,rj12,rj13,rj21,rj22,rj23,rj31,rj32,
     : rj33,r9j)
      implicit real*8(a-h,o-z)
      include 'mpif.h'
      common/factor/fct(0:50),mfd

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
      print *,'inconsistent arguments to 9j symbol'
      stop
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
      if(rkmax.lt.rkmin) then
         print *,'something wrong with k'
         stop
      endif
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
c     1 stands for even and 2 stands for odd parity
c     even parity
      mtbl(1,1)=1
      mtbl(2,2)=1
      mtbl(3,3)=1
      mtbl(4,4)=1
      mtbl(5,5)=1

      mtbl(1,3)=1
      mtbl(3,1)=1

      mtbl(2,4)=1
      mtbl(4,2)=1

      mtbl(1,5)=1
      mtbl(5,1)=1

      mtbl(3,5)=1
      mtbl(5,3)=1

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

      return
      end

      subroutine setgrd
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
      implicit real*8 (a-h,o-z)
      include 'mpif.h'


      call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)
      if(iproc .eq. MASTER)then

      open(STDIN,file='input.ea',form='formatted',status='old')
      open(STDOUT,file='ccea_all.out',form='formatted',status='unknown')
      open(WFNIN,file='wfn.dat',form='unformatted',status='old')
      open(NTFILE,file='to.dat',form='unformatted',status='old')
      endif
      return
      end

      subroutine readinp(nsym,nbas,iocc,maxvalence,valorb,
     &  iopt,itwoe,ivbar,nprint,simopt,t3opt,sfile)
      implicit real*8(a-h,o-z)
      include 'mpif.h'
      integer simopt,t3opt
      integer valorb(NVAL)
      character*14 sfile(NVAL)
      common/iallparity/iall
      dimension nbas(MNSYM),iocc(MNSYM)

      call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)

      if(iproc .eq. MASTER)then


c     read number of symmetry
      read(STDIN,*)nsym,(nbas(i),i=1,nsym)
c     occupancy of each orbital
      read(STDIN,*)(iocc(i),i=1,nsym)
c     truncation scheme and option flag
c iopt =!0 excited to excited bloch is taken
c t3opt= S_1V_2 and S_2V_1 operator considered only for energy
c iall = 0, only even rank operator considered
c itwoe = 1, twoea.dat file is used
c ivar =! two body part of open shell considered
c nprint =!0 intermediate print
c simopt =0 for jacobi diagonalization. It is fored to be zero
      read(STDIN,*)iopt,t3opt,iall,itwoe,ivbar,nprint,simopt
c     information of active orbital
      read(STDIN,*) maxvalence, (valorb(i),i=1,maxvalence)
      read(STDIN,'(a)') (sfile(i),i=1,maxvalence)
      endif

       call MPI_BCAST(sfile,20,MPI_CHARACTER,MASTER,
     :  MPI_COMM_WORLD,ierr)
      call MPI_BCAST(nsym,1,MPI_INTEGER,MASTER,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(nbas,nsym,MPI_INTEGER,MASTER,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(maxvalence,1,MPI_INTEGER,MASTER
     :  ,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(valorb,maxvalence,MPI_INTEGER,MASTER
     : ,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(iocc,nsym,MPI_INTEGER,MASTER,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(iopt,1,MPI_INTEGER,MASTER,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(iall,1,MPI_INTEGER,MASTER,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(itwoe,1,MPI_INTEGER,MASTER,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(ivbar,1,MPI_INTEGER,MASTER,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(simopt,1,MPI_INTEGER,MASTER,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(nprint,1,MPI_INTEGER,MASTER,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(t3opt,1,MPI_INTEGER,MASTER,MPI_COMM_WORLD,ierr)

      return
      end

      subroutine grasprd(nbasis,iproc)
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
      call MPI_BCAST(mtp,1,MPI_INTEGER,MASTER,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(r,MN,MPI_REAL8,MASTER,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(rp,MN,MPI_REAL8,MASTER,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(rpor,MN,MPI_REAL8,MASTER,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(pf,MNBAS*MN,MPI_REAL8,MASTER,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(qf,MNBAS*MN,MPI_REAL8,MASTER,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(eorb,MNBAS,MPI_REAL8,MASTER,MPI_COMM_WORLD,ierr)


      if(iproc.eq.0)write(STDOUT,9999)
      if(iproc.eq.0)write(STDOUT,9998)(eorb(i),i=1,nbasis)

      mtp=n

9999  format(//,10x,17h ORBITAL ENERGIES,/)
9998  format(4d20.11)
      return
      end
     

      subroutine setup(nsym,nbasis,nbas,iocc,iproc)
      implicit real*8 (a-h,o-z)
      include 'mpif.h'
      character*1 orbtyp,orbsym
      character*4 qj1,qj
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
      dimension orbtyp(MNS),nkj(MNBAS),nka(MNBAS)
      dimension orb(MNSYM),np(MNBAS)
      dimension qj(MNBAS),qj1(MNS)
      dimension nprimc(MNOCC),nkprimc(MNOCC)
      dimension nprime(MNEXC),nkprime(MNEXC)
      dimension nmc(MNBAS),nme(MNBAS),nm(MNBAS)
      dimension nbas(MNSYM),iocc(MNSYM)
      data kap/-1, 1,-2, 2,-3, 3,-4,4,-5,5,-6/
      data orb/0.5,0.5,1.5,1.5,2.5,2.5,3.5,3.5,4.5,4.5,5.5/
      data ia /1,-1,1,-1,1,-1,1,-1,1,-1,1/
      data iorbsm/1,2,2,3,3,4,4,5,5,6,6/
      data qj1/'1/2','3/2','5/2','7/2','9/2','11/2'/
      data orbtyp/'s','p','d','f','g','h'/

      nocc=0
      nexcit=0
      nbasis=0
      jmax=2*orb(nsym)-1
      ii=0
      do isym=1,nsym
      ibasis=nbas(isym)
      if(ibasis.ne.0) then
      nbasis = nbasis + ibasis
      do jbas=1,ibasis
      ii = ii + 1
      nka(ii)=kap(isym)
      orbj(ii)=orb(isym)
      np(ii)=jbas
      iparity(ii)=iorbsm(isym)
      if(nka(ii).eq.kap(isym))iiq(ii)=ia(isym)
      if(nka(ii).eq.kap(isym).and.jbas.le.iocc(isym))then
      nocc=nocc+1
      iType(ii) = 0
      iTCount(ii) = nocc
      kc(nocc)=ii
      nkprimc(nocc)=kap(isym)
      nprimc(nocc)=np(ii)
      orbc(nocc)=orbj(ii)
      nmc(nocc)=idint(2*(orbc(nocc))) +1
      if(nkprimc(nocc).eq.kap(isym))isymc(nocc)=iorbsm(isym)
      if(nkprimc(nocc).eq.kap(isym))iqc(nocc)=ia(isym)
      else
      nexcit=nexcit+1
      iType(ii) = 1
      iTCount(ii) = nexcit
      ke(nexcit)=ii
      nprime(nexcit)=np(ii)
      nkprime(nexcit)=nka(ii)
      orbe(nexcit)=orb(isym)
      nme(nexcit)=idint(2*(orbj(ii))) +1
      if(nkprime(nexcit).eq.kap(isym))isyme(nexcit)=iorbsm(isym)
      if(nkprime(nexcit).eq.kap(isym))iqe(nexcit)=ia(isym)
      nup=nup+1
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
      if(isymc(i).eq.1)orbsym=orbtyp(1)
      if(isymc(i).eq.2)orbsym=orbtyp(2)
      if(isymc(i).eq.3)orbsym=orbtyp(3)
      if(isymc(i).eq.4)orbsym=orbtyp(4)
      if(isymc(i).eq.5)orbsym=orbtyp(5)
      if(isymc(i).eq.3)isymc(i)=1
      if(isymc(i).eq.4)isymc(i)=2
      if(isymc(i).eq.5)isymc(i)=1
      if(iproc.eq.0)write(STDOUT,200)orbsym,nprimc(i),
     &nmc(i),kc(i),qj(i)
      enddo

      if(iproc.eq.0)write(STDOUT,102)
      if(iproc.eq.0)write(STDOUT,100)
      do i=1,nexcit
      if(orbe(i).eq.0.5)qj(i)=qj1(1)
      if(orbe(i).eq.1.5)qj(i)=qj1(2)
      if(orbe(i).eq.2.5)qj(i)=qj1(3)
      if(orbe(i).eq.3.5)qj(i)=qj1(4)
      if(orbe(i).eq.4.5)qj(i)=qj1(5)
      if(isyme(i).eq.1)orbsym=orbtyp(1)
      if(isyme(i).eq.2)orbsym=orbtyp(2)
      if(isyme(i).eq.3)orbsym=orbtyp(3)
      if(isyme(i).eq.4)orbsym=orbtyp(4)
      if(isyme(i).eq.5)orbsym=orbtyp(5)
      if(isyme(i).eq.3)isyme(i)=1
      if(isyme(i).eq.4)isyme(i)=2
      if(isyme(i).eq.5)isyme(i)=1
      if(iproc.eq.0)write(STDOUT,200)orbsym,nprime(i),
     &nme(i),ke(i),qj(i)
      enddo
 100  format(/,1x,'TYPE',2x,' N ',2x,'2J+1',2x,'ORB. NO.',
     $3x,' J ',/)
 101  format(/,10x,' OCCUPIED ORBITALS',/)
 102  format(/,10x,' EXCITED ORBITALS',/)
 200  format(2x,a,2x,3(i4,2x),6x,a3)
      return
      end



      subroutine setqic
      implicit real*8 (a-h,o-z)
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
 
C     NUMERICAL CONSTANTS

      ZERO=0.0D0
      ONE=1.0D0
      TWO=2.0D0
      THREE=3.0D0
      TEN=10.0D0
      TENTH=ONE/TEN
      HALF=ONE/TWO

      RETURN
      END

      subroutine daclose
      implicit real*8 (a-h,o-z)
      include 'mpif.h'
      
      call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)

      if(iproc .eq. MASTER)then

      close (STDIN)
      close (STDOUT)
      close (WFNIN)
      endif
      return
      end


c************************************************************************
c                                                                       *
c     SUBROUTINE EAMATA                                                 *
c                                                                       *
c   This subroutine calculate the 'A' matrixes                          *
c                                                                       *
c************************************************************************



      subroutine eamata(iorb,nsmax,isolv,iopt,iproc)

      implicit real*8 (a-h,o-z)
      include 'mpif.h'
      real*8 ad_local,baux_local
c     common block
      common/jvalue/orbc(MNOCC),orbe(MNEXC)
      common/index/ke(MNEXC),kc(MNOCC)
      common/info/nocc,nexcit,iType(MNBAS),iTCount(MNBAS)
      common/parity/isymc(MNOCC),isyme(MNEXC)
      common/symmetry/mtbl(MNS,MNS)
      common/kpavlaue/iqc(MNOCC),iqe(MNEXC)
      common/maxj/jmax
      common/skip1/npbps(0:MXV,2)
      common/skip2/mmskip(0:MXV,2),m2skip
      common/skip3/ncis
      common/idra1/idra(MNBAS,MNBAS,0:MXV)
      common/idpa1/idpa(MNBAS,MNBAS)
      common/jdra1/jdra(MNBAS,MNBAS,0:MXV)
      common/jdpa1/jdpa(MNBAS,MNBAS)
      common/f_bar/f(MNBAS,MNBAS)
      common/scratch1/baux(MDIM),ad(MDIM)
c      common/v_bar/vb(N2ORB,N2BAS,0:MXV)
c      common/mind/mIndex(0:20)


      dimension kgot(MXV+1),ad_local(MDIM),baux_local(MDIM)

      data zero,half,one,two/0.0d0,0.50d0,1.0d0,2.0d0/
      call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)

      do i= 1,MDIM
       ad_local(i)=0.0
      enddo

      do i= 1,MDIM
       baux_local(i)=0.0
      enddo



      nbasis=nocc+nexcit
     
      if(iopt.ne.0)then

c----------------------Single excitation  <P|S2|Q>------------------------


      do 1 ip=iproc+1,nexcit,nproc
      ipsym=isyme(ip)
      iip=ke(ip)
      orbp=orbe(ip)

      do 2 iq=1,nexcit
      iqsym=isyme(iq)
      iiq=ke(iq)
      orbq=orbe(iq)

      if(iip.ne.iorb)go to 2
      if(iip.eq.iiq)go to 2
      if(orbp.ne.orbq)go to 2
      if(ipsym.ne.iqsym)go to 2
      in=jdpa(iip,iiq)
      baux_local(in)=-f(iiq,iip)
      if(isolv.eq.0)then
      ad_local(in)=f(iiq,iiq)-f(iip,iip)
      endif

  2   continue
  1   continue
c            call MPI_ALLREDUCE(ad_local,ad,MDIM,MPI_REAL8,MPI_SUM,
c    :                   MPI_COMM_WORLD,ierr)

c     do i= 1,MDIM
c      ad_local(i)=0.0
c     enddo

c            call MPI_ALLREDUCE(baux_local,baux,MDIM,MPI_REAL8,MPI_SUM,
c    :                   MPI_COMM_WORLD,ierr)

c     do i= 1,MDIM
c      baux_local(i)=0.0
c     enddo


      endif

c---------------Double excitation     <PQ|S2|RA>------------------------------
      
      do 10 ip=iproc+1,nexcit,nproc
      ipsym=isyme(ip)

      do 20 iq=1,nexcit
      iqsym=isyme(iq)

      ipqsym=mtbl(ipsym,iqsym)

      do 30 ir=1,nexcit
      irsym=isyme(ir)

      ipqrsm=mtbl(ipqsym,irsym)

      do 40 ig=1,nocc
      igsym=isymc(ig)

      if(ipqrsm.ne.igsym)go to 40

      orbg=orbc(ig)
      orbp=orbe(ip)
      orbq=orbe(iq)
      orbr=orbe(ir)

      iap=iqe(ip)
      iag=iqc(ig)
      iar=iqe(ir)
      iaq=iqe(iq)

      iip=ke(ip)
      iiq=ke(iq)
      iir=ke(ir)
      iig=kc(ig)

      if(iip.ne.iorb)go to 40

      call findk(orbp,orbq,orbr,orbg,iap,iaq,iar,iag,kgot,
     $kmax)
      if(kmax.eq.0)go to 40
      do 50 nloop=1,kmax
      kk=kgot(nloop)
      ak=float(kk)
      npq=jdra(iip,iiq,kk)
      nrg=idra(iir,iig,kk)
      in=nsloc(npq,nrg,ipqsym,kk)
c     rk=vb(iiq,iir,iip,iig,kk+1)
      iqr=idis(iiq,iir)
      ipg=indxr(iip)+iig
      if(iir.gt.iiq)ipg=indxr(iig)+iip
c      rk=vb(iqr,ipg,kk)
      rk=GetVbar(iqr,ipg,kk)
      baux_local(in)=-rk
      if(isolv.eq.0)then
      ad_local(in)=f(iiq,iiq)+f(iir,iir)-f(iig,iig)-f(iip,iip)
      endif

  50  continue
  40  continue
  30  continue
  20  continue
  10  continue
             call MPI_ALLREDUCE(ad_local,ad,MDIM,MPI_REAL8,MPI_SUM,
     :                   MPI_COMM_WORLD,ierr)

      do i= 1,MDIM
       ad_local(i)=0.0
      enddo

             call MPI_ALLREDUCE(baux_local,baux,MDIM,MPI_REAL8,MPI_SUM,
     :                   MPI_COMM_WORLD,ierr)

      do i= 1,MDIM
       baux_local(i)=0.0
      enddo


 9999 format(//,' A Matrix of A+BS=S*Heff',//)
 9998 format(4f20.11)
 9997 format(//,' Diagonal part of B (MP)',//)

      return
      end

c********************************************************************
c                                                                   *
c    SUBROUTINE  HEFF                                               *
c                                                                   *
c    This subroutine calculate the electron affinity                *
c                                                                   *
c********************************************************************

      subroutine heff(iorb,t,b,iproc)

      implicit real*8 (a-h,o-z)
      include 'mpif.h'
      real*8 sum_local
c     common block
      common/jvalue/orbc(MNOCC),orbe(MNEXC)
      common/index/ke(MNEXC),kc(MNOCC)
      common/info/nocc,nexcit,iType(MNBAS),iTCount(MNBAS)
      common/parity/isymc(MNOCC),isyme(MNEXC)
      common/symmetry/mtbl(MNS,MNS)
      common/kpavlaue/iqc(MNOCC),iqe(MNEXC)
      common/jdra1/jdra(MNBAS,MNBAS,0:MXV)
      common/idra1/idra(MNBAS,MNBAS,0:MXV)
      common/idpa1/idpa(MNBAS,MNBAS)
      common/jdpa1/jdpa(MNBAS,MNBAS)
      common/f_bar/f(MNBAS,MNBAS)
      common/hefen/hef

      dimension t(MDIM)
      dimension kgot(MXV+1),igot(MXV+1),jgot(MXV+1)
      dimension a(MNEXC,MNEXC),b(MNEXC,MNEXC)


      data zero,half,one,two/0.0,0.5,1.0,2.0/
      real*8 GetVbar

      call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)

      do i=0,MXV
      kgot(i)=0
      jgot(i)=0
      igot(i)=0
      enddo

      if(iproc.eq.0)write(STDOUT,9999)
      if(iproc.eq.0)write(STDOUT,9998)

      do i=1,nexcit
      do j=1,nexcit
      a(i,j)=0.0d0
      b(i,j)=0.0d0
      enddo
      enddo

      do 1 ip=1,nexcit
      ipsym=isyme(ip)
      iip=ke(ip)
      orbp=orbe(ip)
      iap=iqe(ip)

      if(iip.ne.iorb)go to 1

      sum1_local=0
      sum2_local=0

      do 4 id=1,nocc
      idsym=isymc(id)
      iid=kc(id)
      orbd=orbc(id)
      iad=iqc(id)

      ipdsym=mtbl(ipsym,idsym)

      do 5 ir=iproc+1,nexcit,nproc
      irsym=isyme(ir)
      iir=ke(ir)
      iar=iqe(ir)
      orbr=orbe(ir)

      ipdrsm=mtbl(ipdsym,irsym)

      do 6 is=1,nexcit
      issym=isyme(is)
      iis=ke(is)
      ias=iqe(is)
      orbs=orbe(is)

      if(issym.ne.ipdrsm)go to 6

      call findk(orbp,orbr,orbs,orbd,iap,iar,ias,iad,jgot,
     $jmax)
      if(jmax.eq.0)go to 6
      do 7 iloop=1,jmax
      kk=jgot(iloop)
      ak=float(kk)
      ipd=idis(iip,iid)
      irs=indxr(iir)+iis
      if(iid.gt.iip)irs=indxr(iis)+iir
      rk=GetVbar(ipd,irs,kk)

c--------------diagram   1  <ps|v|rd><ps|t2|rd>-------------------------------

      npr=jdra(iip,iir,kk)
      nsd=idra(iis,iid,kk)
      isdsym=mtbl(issym,idsym)
      in=nsloc(npr,nsd,isdsym,kk)
      f1=(-1)**(orbp+orbs+orbd+orbr)
      f2=1/((two*kk+1)*(two*orbp+1))
      sum1_local=sum1_local+rk*f1*f2*t(in)

c---------------diagram    2  <ps|v|rd><pr|t2|sd>-------------------------------

      call findk(orbp,orbs,orbr,orbd,iap,ias,iar,iad,igot,
     $imax)
      if(imax.eq.0)go to 7
      do 8 mloop=1,imax
      ll=igot(mloop)
      al=float(ll)
      nps=jdra(iip,iis,ll)
      nrd=idra(iir,iid,ll)
      irdsym=mtbl(irsym,idsym)
      in=nsloc(nps,nrd,irdsym,ll)
      f2=1/(two*orbp+1)
      call sixj(ak,orbs,orbd,al,orbr,orbp,s6j)
      sum2_local=sum2_local+rk*t(in)*f1*f2*s6j

  8   continue
  7   continue
  6   continue
  5   continue
  4   continue
       call MPI_ALLREDUCE(sum1_local,sum1,1,MPI_REAL8,MPI_SUM,
     :                   MPI_COMM_WORLD,ierr)
       sum1_local=0
       call MPI_ALLREDUCE(sum2_local,sum2,1,MPI_REAL8,MPI_SUM,
     :                   MPI_COMM_WORLD,ierr)
       sum2_local=0



      sum_local3=0

      do 15 iq=iproc+1,nexcit,nproc
      iqsym=isyme(iq)
      orbq=orbe(iq)
      iiq=ke(iq)

      if(iiq.eq.iip)go to 15
      if(ipsym.ne.iqsym)go to 15
      if(orbq.ne.orbp)go to 15

c-----------------diagram  3   <p|f|q><p|t1|q>--------------------------------

      in=jdpa(iip,iiq)
      sum3_local=sum3_local+f(iip,iiq)*t(in)

  15  continue
       call MPI_ALLREDUCE(sum3_local,sum3,1,MPI_REAL8,MPI_SUM,
     :                   MPI_COMM_WORLD,ierr)
       sum3_local=0


      sum4_local=0
      sum5_local=0

      do 16 ig=1,nocc
      igsym=isymc(ig)
      iig=kc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)

      do 16 iq=iproc+1,nexcit,nproc
      iqsym=isyme(iq)
      iiq=ke(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)

      if(orbq.ne.orbg)go to 16
      if(iqsym.ne.igsym)go to 16

c-------------------- diagram   4   <pq|v|pg><q|t1|g>-----------------------------
                                                                                
      call findk(orbp,orbp,orbq,orbg,iap,iap,iaq,iag,igot,
     $imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      ipp=jdra(iip,iip,l1)
      iqg=idra(iiq,iig,l1)
      iqgsym=mtbl(iqsym,igsym)
      in=nsloc(ipp,iqg,iqgsym,l1)
      amult=dsqrt((two*orbg+one)/(two*orbp+one))
      sum4_local=sum4_local+f(iig,iiq)*amult*t(in)

      endif
      enddo
      endif

c---------------------diagram   5  <pp|v|qg><q|t1|g>-----------------------------

      call findk(orbp,orbq,orbp,orbg,iap,iaq,iap,iag,igot,
     $imax)
      if(imax.eq.0)go to 17
      do 18 mloop=1,imax
      kk=igot(mloop)
      ipq=jdra(iip,iiq,kk)
      ipg=idra(iip,iig,kk)
      ipgsym=mtbl(ipsym,igsym)
      in=nsloc(ipq,ipg,ipgsym,kk)
      amult=one/(two*orbp+one)
      sum5_local=sum5_local+f(iig,iiq)*amult*t(in)*(-1)**(orbg+orbp+kk)

  18  continue
  17  continue
  16  continue
       call MPI_ALLREDUCE(sum4_local,sum4,1,MPI_REAL8,MPI_SUM,
     :                   MPI_COMM_WORLD,ierr)
       sum4_local=0
      call MPI_ALLREDUCE(sum5_local,sum5,1,MPI_REAL8,MPI_SUM,
     :                   MPI_COMM_WORLD,ierr)
       sum5_local=0


c------------------- Calculate total electron affinity------------------------

      a(ip,ip)=a(ip,ip)+sum1+sum2+sum3+sum4+sum5
      b(iip,iip)=a(ip,ip)
      if(iproc.eq.0)then
      write(STDOUT,9997)sum1,sum2,a(ip,ip),a(ip,ip)+f(iip,iip)
      endif

  1   continue
9996  format(6f12.8)
9997  format(4f20.11)
9998  format(8x,'DIRECT PART',7x,'EXCHANGE PART',10x,
     :       'SUM',17x,'TOTAL'/)
9999  format(//,1x,'OSCC ELECTRON AFFINITY (AU)',//)

      return
      end

c     this subroutine constructs the nonlinear a matrix for oscc-ip
      subroutine nlheff(iorb,a,s,t,nsmax,t3opt,isolv,rows,rowe,
     $iopt,iproc)
      implicit real*8 (a-h,o-z)
      include 'mpif.h'
c     real*8 sum_local
c     common block
      integer t3opt,rows,rowe,row
      common/jvalue/orbc(MNOCC),orbe(MNEXC)
      common/index/ke(MNEXC),kc(MNOCC)
      common/info/nocc,nexcit,iType(MNBAS),iTCount(MNBAS)
      common/parity/isymc(MNOCC),isyme(MNEXC)
      common/jdra1/jdra(MNBAS,MNBAS,0:MXV)
      common/idra1/idra(MNBAS,MNBAS,0:MXV)
      common/idpa1/idpa(MNBAS,MNBAS)
      common/jdpa1/jdpa(MNBAS,MNBAS)
      common/symmetry/mtbl(MNS,MNS)
      common/skip2/mmskip(0:MXV,2),m2skip
      common/kpavlaue/iqc(MNOCC),iqe(MNEXC)
      common/f_bar/f(MNBAS,MNBAS)
      common/scratch2/raux(MDIM),caux(MDIM)

      dimension kgot(MXV+1),a(MNEXC,MNEXC)
      dimension s(MDIM),baux(MDIM),t(MDIM)

      data zero,half,one,two/0.0d0,0.50d0,1.0d0,2.0d0/
c     call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
c     call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)

      nbasis=nocc+nexcit
      t3cont=0

      if(t3opt.ne.0)call ccsdtea(iorb,s,t,t3cont,iproc)
   
      if(iopt.ne.0)then

      do 1 ip=1,nexcit
      ipsym=isyme(ip)
      iip=ke(ip)
      orbp=orbe(ip)
      do 2 iq=1,nexcit
      iqsym=isyme(iq)
      iiq=ke(iq)
      orbq=orbe(iq)
      if(iip.ne.iorb)go to 2
      if(iip.eq.iiq)go to 2
      if(orbp.ne.orbq)go to 2
      if(ipsym.ne.iqsym)go to 2
      in=jdpa(iip,iiq)

      row=in
      irow=row+1-rows
      if(irow.le.0)go to 2
      if(irow.gt.rowe-rows+1)go to 2

      caux(in)=a(iip,iip)*s(in)+t3cont*s(in)
  2   continue
  1   continue

      endif
      
      do 10 ip=1,nexcit
      ipsym=isyme(ip)
      do 20 iq=1,nexcit
      iqsym=isyme(iq)
      ipqsym=mtbl(ipsym,iqsym)
      do 30 ir=1,nexcit
      irsym=isyme(ir)
      ipqrsm=mtbl(ipqsym,irsym)
      do 40 ig=1,nocc
      igsym=isymc(ig)
      if(ipqrsm.ne.igsym)go to 40
      orbg=orbc(ig)
      orbp=orbe(ip)
      orbq=orbe(iq)
      orbr=orbe(ir)
      iap=iqe(ip)
      iag=iqc(ig)
      iar=iqe(ir)
      iaq=iqe(iq)
      iip=ke(ip)
      iiq=ke(iq)
      iir=ke(ir)
      iig=kc(ig)
      if(iip.ne.iorb)go to 40
      call findk(orbp,orbq,orbr,orbg,iap,iaq,iar,iag,kgot,
     $kmax)
      if(kmax.eq.0)go to 40
      do 50 nloop=1,kmax
      kk=kgot(nloop)
      ak=float(kk)
      ipq=jdra(iip,iiq,kk)
      irg=idra(iir,iig,kk)
      ipqsym=mtbl(ipsym,iqsym)
      in=nsloc(ipq,irg,ipqsym,kk)

      row=in
      irow=row+1-rows
      if(irow.le.0)go to 50
      if(irow.gt.rowe-rows+1)go to 50

      caux(in)=a(iip,iip)*s(in)+t3cont*s(in)
  50  continue
  40  continue
  30  continue
  20  continue
  10  continue

      return
      end

C
C     Subroutine for triples
C
      subroutine ccsdtea(iorb,s,t,sum,iproc)
      implicit real*8 (a-h,o-z)
      include 'mpif.h'
      real*8 sum_local
c     common block
      common/jvalue/orbc(MNOCC),orbe(MNEXC)
      common/index/ke(MNEXC),kc(MNOCC)
      common/info/nocc,nexcit,iType(MNBAS),iTCount(MNBAS)
      common/parity/isymc(MNOCC),isyme(MNEXC)
      common/symmetry/mtbl(MNS,MNS)
      common/kpavlaue/iqc(MNOCC),iqe(MNEXC)
      common/jdra1/jdra(MNBAS,MNBAS,0:MXV)
      common/idra1/idra(MNBAS,MNBAS,0:MXV)
      common/idpa1/idpa(MNBAS,MNBAS)
      common/jdpa1/jdpa(MNBAS,MNBAS)
      common/orbital_energy/eorb(MNBAS)

      dimension s(MDIM),t(MDIM)
      dimension igot(MXV+1),jgot(MXV+1),kgot(MXV+1)
      data zero,half,one,two/0.0,0.5,1.0,2.0/
      call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)

c     VS3=VVS2 and VS3=VVT2 approximate contribution of triples
c     to effective Hamiltonian.
c
      do 1 ip=1,nexcit
      ipsym=isyme(ip)
      iip=ke(ip)
      orbp=orbe(ip)
      iap=iqe(ip)
      if(iip.ne.iorb)go to 1
                                                                                
      sum_local=zero
c     Diagrams A1-A5, A10-A14, A16 and A23-A26. This is from VS2=S3 part
                                                                                
c     external loop
      do 2 ia=1,nocc
      iia=kc(ia)
      iaa=iqc(ia)
      orba=orbc(ia)
      iasym=isymc(ia)
      do 2 ib=1,nocc
      iib=kc(ib)
      iab=iqc(ib)
      orbb=orbc(ib)
      ibsym=isymc(ib)
      do 2 ir=iproc+1,nexcit,nproc
      irsym=isyme(ir)
      iir=ke(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      do 2 is=1,nexcit
      issym=isyme(is)
      iis=ke(is)
      orbs=orbe(is)
      ias=iqe(is)
      irasym=mtbl(iasym,irsym)
      isbsym=mtbl(issym,ibsym)
      if(irasym.ne.isbsym)go to 2
      call findk(orba,orbr,orbb,orbs,iaa,iar,iab,ias,kgot,kmax)
      if(kmax.eq.0)go to 2
      do 3 kloop=1,kmax
      k1=kgot(kloop)
      ak1=float(k1)
                                                                                
      rk=vint(iia,iib,iir,iis,k1)
      d4=(-1)**(k1+orba+orbb+1)
      vabrs=rk*d4
      denom= (eorb(iia)+eorb(iib)-eorb(iir)-eorb(iis))
                                                                                
c     diagram evaluation
                                                                                
      do 4 iq=1,nexcit
      iiq=ke(iq)
      iqsym=isyme(iq)
      iaq=iqe(iq)
      orbq=orbe(iq)
                                                                                
      ipqsym=mtbl(ipsym,iqsym)
      iqbsym=mtbl(iqsym,ibsym)
      iqasym=mtbl(iqsym,iasym)
      iarsym=mtbl(iasym,irsym)
                                                                                
      if(ipqsym.ne.iarsym)go to 4
      if(ipqsym.ne.isbsym)go to 4
                                                                                
                                                                                
c     diagram A1
      rk=vint(iip,iir,iiq,iia,k1)
      d8=(-1)**(orbp+orbr+k1+1)
      vprqa=rk*d8
                                                                                
      call findk(orbp,orbq,orbb,orbs,iap,iaq,iab,ias,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      kk=jgot(jloop)
      if(kk.eq.k1)then
      npq=jdra(iip,iiq,k1)
      nsb=idra(iis,iib,k1)
      in=nsloc(npq,nsb,isbsym,k1)
      fact=((2*ak1+1)**2)*(2*orbp+1)
      fact=one/fact
      isign=(-1)**(orba+orbb+orbs+orbr+orbp+orbq+k1)
      sum_local=sum_local-vabrs*vprqa*s(in)*fact*isign/denom
      endif
      enddo
      endif
                                                                                
c     diagram A2
      call findk(orba,orbs,orbp,orbq,iaa,ias,iap,iaq,jgot,jmax)
      if(jmax.eq.0)go to 5
      do 6 jloop=1,jmax
      k2=jgot(jloop)
      ak2=float(k2)
                                                                                
      rk=vint(iip,iis,iiq,iia,k2)
      d8=(-1)**(orbp+orbs+k2+1)
      vpsqa=rk*d8
      npq=jdra(iip,iiq,k2)
      nrb=idra(iir,iib,k2)
      in=nsloc(npq,nrb,ipqsym,k2)
      call sixj(ak1,orbs,orbb,ak2,orbr,orba,s6j)
      fact=one/((2*orbp+1)*(2*ak2+1))*s6j
      isign=(-1)**(orbp+orbq+orba+orbs+orbr+orbb+k2)
      sum_local=sum_local-vabrs*vpsqa*s(in)*fact*isign/denom
 6    continue
 5    continue
                                                                               
c     diagram A3
      call findk(orba,orbp,orbq,orbr,iaa,iap,iaq,iar,jgot,jmax)
      if(jmax.eq.0)go to 7
      do 8 jloop=1,jmax
      k2=jgot(jloop)
      ak2=float(k2)
                                                                                
      rk=vint(iip,iir,iia,iiq,k2)
      d8=(-1)**(orbp+orbr+k2+1)
      vpraq=rk*d8
      call findk(orbp,orbq,orbb,orbs,iap,iaq,iab,ias,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      kk=igot(iloop)
      if(kk.eq.k1)then
      npq=jdra(iip,iiq,k1)
      nsb=idra(iis,iib,k1)
      in=nsloc(npq,nsb,ipqsym,k1)
      call sixj(ak2,orbr,orbq,ak1,orbp,orba,s6j)
      fact=one/((2*orbp+1)*(2*ak1+1))*s6j
      isign=(-1)**(orbp+orbq+orba+orbs+orbr+orbb+k1)
      sum_local=sum_local-vabrs*vpraq*s(in)*fact*isign/denom
      endif
      enddo
      endif
 8    continue
 7    continue
                                                                                
c     diagram A4
      call findk(orba,orbp,orbq,orbr,iaa,iap,iaq,iar,jgot,jmax)
      if(jmax.eq.0)go to 9
      do 10 jloop=1,jmax
      k2=jgot(jloop)
      ak2=float(k2)
                                                                                
      rk=vint(iip,iir,iia,iiq,k2)
      d8=(-1)**(orbp+orbr+k2+1)
      vpraq=rk*d8
      call findk(orbp,orbs,orbq,orbb,iap,ias,iaq,iab,igot,imax)
      if(imax.eq.0)go to 10
      do 11 iloop=1,imax
      k3=igot(iloop)
      ak3=float(k3)
      call sixj(orbs,orbp,ak3,orbq,orbb,ak1,s6j1)
      call sixj(orba,ak2,orbp,orbq,ak1,orbr,s6j2)
      isign=(-1)**(orba+orbb+orbp+orbq+orbr+orbs+k1)
      fact=one/(two*orbp+1)*s6j1*s6j2
      nps=jdra(iip,iis,k3)
      nqb=idra(iiq,iib,k3)
      in=nsloc(nps,nqb,iqbsym,k3)
      sum_local=sum_local-vabrs*vpraq*s(in)*fact*isign/denom
 11   continue
 10   continue
  9   continue
                                                                                
c     diagram A5
      call findk(orbp,orbs,orbq,orbb,iap,ias,iaq,iab,jgot,jmax)
      if(jmax.eq.0)go to 12
      do 13 jloop=1,jmax
      k2=jgot(jloop)
      ak2=float(k2)
      call sixj(ak1,orbb,orbs,ak2,orbp,orbq,s6j)
      nps=jdra(iip,iis,k2)
      nqb=idra(iiq,iib,k2)
      in=nsloc(nps,nqb,iqbsym,k2)
      fact=(two*ak1+1)*(two*orbp+1)
      fact=s6j/fact
      isign=(-1)**(orba+orbb+orbr+orbp+orbq+orbs+k1)
      jsign=(-1)**(2*orbs+2*orbr)
      sum_local=sum_local-vabrs*vprqa*s(in)*fact*isign/denom
 13   continue
 12   continue
                                                                                
c     diagram A10
      if(orbq.eq.orbb)then
      if(iqsym.eq.ibsym)then
      rk=vint(iir,iis,iia,iiq,k1)
      d8=(-1)**(orbs+orbr+k1+1)
      vrsaq=rk*d8
                                                                                
      call findk(orbp,orbq,orbp,orbb,iap,iaq,iap,iab,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k2=jgot(jloop)
      if(k2.eq.0)then
      npp=jdra(iip,iip,0)
      nqb=idra(iiq,iib,0)
      in=nsloc(npp,nqb,iqbsym,0)
      fact=dsqrt((2*orbq+1)/(2*orbp+1))/((two*ak1+1)*(two*orbq+1))
      isign=(-1)**(orba+orbb+orbs+orbr)
      sum_local=sum_local+vabrs*vrsaq*s(in)*fact*isign/denom
      endif
      enddo
      endif
                                                                                
c     diagram A11
      call findk(orbp,orbq,orbp,orbb,iap,iaq,iap,iab,jgot,jmax)
      if(jmax.eq.0)go to 14
      do 15 jloop=1,jmax
      k2=jgot(jloop)
      ak2=float(k2)
      npq=jdra(iip,iiq,k2)
      npb=idra(iip,iib,k2)
      in=nsloc(npq,npb,ipqsym,k2)
      fact=((2*orbp+1)*(2*orbq+1)*(2*ak1+1))
      fact=one/fact
      isign=(-1)**(orbp+orbq+orba+orbb+orbr+orbs+k2+2*orbp)
      sum_local=sum_local-vabrs*vrsaq*s(in)*fact*isign/denom
 15   continue
 14   continue
                                                                                
      endif
      endif
                                                                                
c     diagram A12
      if(orbq.ne.orba)go to 16
      if(iqsym.ne.iasym)go to 16
      call findk(orbs,orbq,orbr,orbb,ias,iaq,iar,iab,igot,imax)
      if(imax.eq.0) go to 16
      do 17 iloop=1,imax
      k2=jgot(iloop)
      ak2=float(k2)
                                                                                
      rk=vint(iir,iis,iib,iiq,k2)
      d8=(-1)**(orbs+orbr+k2+1)
      vrsbq=rk*d8
                                                                                
      call findk(orbp,orbp,orbq,orbb,iap,iap,iaq,iaa,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k3=jgot(jloop)
      if(k3.eq.0)then
      npp=jdra(iip,iip,0)
      nqa=idra(iiq,iia,0)
      in=nsloc(npp,nqa,iqasym,0)
      call sixj(ak1,orbb,orbs,ak2,orba,orbr,s6j)
      fact=dsqrt((2*orbq+1)/(2*orbp+1))/(two*orbq+1)
      fact=s6j*fact
      isign=(-1)**(orba+orbb+orbs+orbr+2*orbs)
      sum_local=sum_local-vabrs*vrsbq*s(in)*fact*isign/denom
      endif
      enddo
      endif
                                                                                
c     diagram A13
      call findk(orbp,orbq,orbp,orba,iap,iaq,iap,iaa,jgot,jmax)
      if(jmax.eq.0)go to 17
      do 18 jloop=1,jmax
      k3=jgot(jloop)
      ak3=float(k3)
      npq=jdra(iip,iiq,k3)
      npa=idra(iip,iia,k3)
      in=nsloc(npq,npa,ipqsym,k3)
      fact=((2*orbp+1)*(2*orbq+1))
      call sixj(ak1,orbb,orbs,ak3,orba,orbr,s6j)
      fact=s6j/fact
      isign=(-1)**(orbp+orbq+orba+orbb+orbr+orbs+k3)
      sum_local=sum_local+vabrs*vrsbq*s(in)*fact*isign/denom
 18   continue
 17   continue
 16   continue
                                                                                
c     diagram A16
      call findk(orbr,orbq,orbp,orbb,iar,iaq,iap,iab,igot,imax)
      if(imax.eq.0) go to 19
      do 20 iloop=1,imax
      k2=igot(iloop)
      ak2=float(k2)
                                                                                
      rk=vint(iip,iir,iib,iiq,k2)
      d8=(-1)**(orbp+orbr+k2+1)
      vprbq=rk*d8
      call findk(orbp,orbq,orbs,orba,iap,iaq,ias,iaa,jgot,jmax)
      if(jmax.eq.0)go to 20
      do 21 jloop=1,jmax
      k3=jgot(jloop)
      ak3=float(k3)
      npq=jdra(iip,iiq,k3)
      nsa=idra(iis,iia,k3)
      in=nsloc(npq,nsa,ipqsym,k3)
      fact=one/(two*orbp+1)
      call sixj(orbs,ak1,orbb,orbr,ak3,orba,s6j1)
      call sixj(orbp,orbb,ak2,orbr,orbq,ak3,s6j2)
      fact=fact*s6j1*s6j2
      isign=(-1)**(orbp+orbq+orba+orbs+orbb+orbr+k3+2*orbr)
      sum_local=sum_local+vabrs*vprbq*s(in)*fact*isign/denom
 21   continue
 20   continue
 19   continue
                                                                                
c     diagram A23
      if(orbr.eq.orbq)then
      if(irsym.eq.iqsym)then
                                                                                
      call findk(orbs,orbq,orbr,orbb,ias,iaq,iar,iab,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k3=jgot(jloop)
      if(k3.eq.0)then
      rk=vint(iip,iir,iip,iiq,0)
      d8=(-1)**(orbp+orbr+1)
      vprpq=rk*d8
      nqa=idra(iiq,iia,k1)
      nsb=idra(iis,iib,k1)
      in=ntloc(nqa,nsb,iqasym,k1)
      fact=dsqrt((2*orbq+1)/(2*orbp+1))/((2*ak1+1)*(2*orbq+1))
      isign=(-1)**(orba+orbb+orbr+orbs)
      sum_local=sum_local+vabrs*vprpq*t(in)*fact*isign/denom
                                                                                
c     diagram A24
      call findk(orbq,orbb,orbs,orba,iaq,iab,ias,iaa,igot,imax)
      if(imax.eq.0)go to 22
      do 23 iloop=1,imax
      k2=igot(iloop)
      ak2=float(k2)
      nqb=idra(iiq,iib,k2)
      nsa=idra(iis,iia,k2)
      in=ntloc(nqb,nsa,iqbsym,k2)
      fact=dsqrt((2*orbq+1)/(2*orbp+1))
      call sixj(ak1,orbb,orbs,ak2,orba,orbr,s6j)
      fact=fact*s6j/(two*orbq+1)
      isign=(-1)**(orba+orbb+orbr+orbs+2*orbs)
      sum_local=sum_local-vabrs*vprpq*t(in)*fact*isign/denom
 23   continue
 22   continue
                                                                                
      endif
      enddo
      endif
                                                                                
      endif
      endif
                                                                                
c     diagram A25
      call findk(orbp,orba,orbq,orbr,iap,iaa,iaq,iar,igot,imax)
      if(imax.eq.0)go to 24
      do 25 iloop=1,imax
      k2=igot(iloop)
      ak2=float(k2)
                                                                                
      rk=vint(iip,iir,iia,iiq,k2)
      d8=(-1)**(orbp+orbr+k2+1)
      vpraq=rk*d8
      call findk(orbp,orbq,orbs,orbb,iap,iaq,ias,iab,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      kk=jgot(jloop)
      if(kk.eq.k1)then
      npq=jdra(iip,iiq,k1)
      nsb=idra(iis,iib,k1)
      in=nsloc(npq,nsb,ipqsym,k1)
      call sixj(ak1,orbq,orbr,ak2,orba,orbp,s6j)
      fact=(two*ak1+1)*(two*orbp+1)
      fact=s6j/fact
      isign=(-1)**(orbs+orbb+orbp+orbq+orba+orbr+k1)
      sum_local=sum_local-vabrs*vpraq*s(in)*fact*isign/denom
      endif
      enddo
      endif
                                                                                
c     diagram A26
      call findk(orbp,orbs,orbq,orbb,iap,ias,iaq,iab,jgot,jmax)
      if(jmax.eq.0)go to 25
      do 26 jloop=1,jmax
      k3=jgot(jloop)
      ak3=float(k3)
      nps=jdra(iip,iis,k3)
      nqb=idra(iiq,iib,k3)
      in=nsloc(nps,nqb,iqbsym,k3)
      call sixj(orba,ak3,orbp,orbq,ak1,orbr,s6j1)
      call sixj(orbs,orbp,ak2,orbq,orbb,ak1,s6j2)
      isign=(-1)**(orbr+orba+orbp+orbq+orbs+orbb+k1+2*orbq)
      fact=s6j1*s6j2/(two*orbp+1)
      sum_local=sum_local+vabrs*vpraq*s(in)*fact*isign/denom
 26   continue
 25   continue
 24   continue
  4   continue
                                                                                
c     diagrams A6-A9, A14-A15, A17-A22, A27-A29
c     diagram A6
      do 30 ig=1,nocc
      iig=kc(ig)
      iag=iqc(ig)
      orbg=orbc(ig)
      igsym=isymc(ig)
                                                                                
      ipgsym=mtbl(ipsym,igsym)
      isgsym=mtbl(issym,igsym)
      iabsym=mtbl(iasym,ibsym)
      irgsym=mtbl(irsym,igsym)
                                                                                
      if(iabsym.ne.irgsym)go to 30
      if(orbg.ne.orbs)go to 31
                                                                                
      call findk(orbp,orbp,orbs,orbg,iap,iap,ias,iag,igot,imax)
      if (imax .ne. 0) then
      do iloop=1,imax
      k2=igot(iloop)
      if(k2.eq.0)then
                                                                                
      rk=vint(iir,iig,iia,iib,k1)
      d8=(-1)**(orbr+orbg+k1+1)
      vrgab=rk*d8
      fact=dsqrt((2*orbs+1)/(2*orbp+1))/(2*orbg+1)
      fact=fact/(two*ak1+1)
      isign=(-1)**(orba+orbr+orbs+orbb)
      npp=jdra(iip,iip,0)
      nsg=idra(iis,iig,0)
      in=nsloc(npp,nsg,isgsym,0)
      sum_local=sum_local-vabrs*vrgab*s(in)*fact*isign/denom
                                                                                
      endif
      enddo
      endif
                                                                                
c     diagram A7
      call findk(orbp,orbs,orbp,orbg,iap,ias,iap,iag,jgot,jmax)
      if (jmax.gt.0) then
      do 32 jloop=1,jmax
      k2=jgot(jloop)
      ak2=float(k2)
      nps=jdra(iip,iis,k2)
      npg=idra(iip,iig,k2)
      in=nsloc(nps,npg,ipgsym,k2)
      fact=one/((two*orbs+1)*(two*orbp+1)*(two*ak1+1))
      isign=(-1)**(orba+orbb+orbs+orbr+orbp+orbg+k2+2*orbp)
      sum_local=sum_local+vabrs*vrgab*s(in)*fact*isign/denom
 32   continue
      endif
                                                                                
c     diagram A8
      call findk(orbr,orbb,orba,orbg,iar,iab,iaa,iag,jgot,jmax)
      if(jmax.eq.0) goto 31
      do 33 jloop=1,jmax
      k2=jgot(jloop)
      ak2=float(k2)
                                                                                
      rk=vint(iir,iig,iib,iia,k2)
      d8=(-1)**(orbr+orbg+k2+1)
      vrgba=rk*d8
                                                                                
      call findk(orbp,orbp,orbs,orbg,iap,iap,ias,iag,igot,imax)
      if (imax.ne.0) then
      do iloop=1,imax
      k3=igot(iloop)
      if(k3.eq.0)then
                                                                                
      npp=jdra(iip,iip,0)
      nsg=idra(iis,iig,0)
      in=nsloc(npp,nsg,isgsym,0)
      fact=dsqrt((2*orbs+1)/(2*orbp+1))/(two*orbs+1)
      call sixj(ak1,orbr,orba,ak2,orbs,orbb,s6j)
      fact=s6j*fact
      isign=(-1)**(orba+orbb+orbs+orbr)
      sum_local=sum_local-vabrs*vrgba*s(in)*fact*isign/denom
      endif
      enddo
      endif
                                                                                
c     diagram A9
      call findk(orbp,orbs,orbp,orbg,iap,ias,iap,iag,igot,imax)
      if(imax.eq.0)go to 33
      do 34 iloop=1,imax
      k3=igot(iloop)
      ak3=float(k3)
      nps=jdra(iip,iis,k3)
      npg=idra(iip,iig,k3)
      in=nsloc(nps,npg,ipgsym,k3)
      call sixj(ak1,orbr,orba,ak2,orbs,orbb,s6j)
      fact=one/((two*orbp+1)*(two*orbs+1))
      fact=fact*s6j
      isign=(-1)**(orbp+orbg+orbb+orbr+orba+orbs+k3)
      sum_local=sum_local-vabrs*vrgba*s(in)*fact*isign/denom
  34  continue
  33  continue
  31  continue
                                                                                
c     diagram A15
      call findk(orbp,orba,orbb,orbg,iap,iaa,iab,iag,igot,imax)
      if(imax.eq.0)go to 35
      do 36 iloop=1,imax
      k2=igot(iloop)
      ak2=float(k2)
                                                                                
      rk=vint(iip,iig,iia,iib,k2)
      d8=(-1)**(orbp+orbg+k2+1)
      vpgab=rk*d8
      call findk(orbp,orbr,orbs,orbg,iap,iar,ias,iag,jgot,jmax)
      if(jmax.eq.0)go to 36
      do 37 jloop=1,jmax
      k3=jgot(jloop)
      ak3=float(k3)
      npr=jdra(iip,iir,k3)
      nsg=idra(iis,iig,k3)
      in=nsloc(npr,nsg,isgsym,k3)
      call sixj(orba,orbp,ak2,ak3,ak1,orbr,s6j1)
      call sixj(orbb,ak2,orbg,ak3,orbs,ak1,s6j2)
      fact=s6j1*s6j2/(two*orbp+1)
      sum_local=sum_local+vabrs*vpgab*s(in)*fact/denom
  37  continue
  36  continue
  35  continue
                                                                                
c     diagram A17
      call findk(orbp,orbb,orba,orbg,iap,iab,iaa,iag,igot,imax)
      if(imax.eq.0)go to 38
      do 39 iloop=1,imax
      k2=igot(iloop)
      ak2=float(k2)
                                                                                
      rk=vint(iip,iig,iib,iia,k2)
      d8=(-1)**(orbp+orbg+k2+1)
      vpgba=rk*d8
      call findk(orbp,orbr,orbs,orbg,iap,iar,ias,iag,jgot,jmax)
      if(jmax.eq.0)go to 39
      do 40 jloop=1,jmax
      k3=jgot(jloop)
      ak3=float(k3)
      npr=jdra(iip,iir,k3)
      nsg=idra(iis,iig,k3)
      in=nsloc(npr,nsg,isgsym,k3)
      call ninej(orba,ak2,orbg,orbr,orbp,ak3,ak1,orbb,orbs,s9j)
      fact=s9j/(two*orbp+1)
      sum_local=sum_local+vabrs*vpgba*s(in)*fact/denom
  40  continue
  39  continue
  38  continue
                                                                                
c     diagram A18
      call findk(orbs,orbb,orbp,orbg,ias,iab,iap,iag,igot,imax)
      if(imax .ne. 0) then
      do iloop=1,imax
      k2=igot(iloop)
      if(k2.eq.k1)then
                                                                                
      rk=vint(iir,iig,iia,iip,k1)
      d8=(-1)**(orbr+orbg+k1+1)
      vrgap=rk*d8
      nsb=idra(iis,iib,k1)
      npg=idra(iip,iig,k1)
      in=ntloc(npg,nsb,ipgsym,k1)
      fact=(two*ak1+1)**2*(two*orbp+1)
      fact=one/fact
      isign=(-1)**(orba+orbr+orbb+orbs+orbp+orbg+k1)
      sum_local=sum_local+vabrs*vrgap*t(in)*fact*isign/denom
      endif
      enddo
      endif
                                                                                
c     diagram A19
      call findk(orbp,orbg,orbr,orbb,iap,iag,iar,iab,jgot,jmax)
      if(jmax.eq.0)go to 41
      do 42 jloop=1,jmax
      k2=jgot(jloop)
      ak2=float(k2)
                                                                                
      rk=vint(iis,iig,iia,iip,k2)
      d8=(-1)**(orbs+orbg+k2+1)
      vsgap=rk*d8
      call sixj(ak1,orbs,orbb,ak2,orbr,orba,s6j)
      fact=(two*ak2+1)*(two*orbp+1)
      fact=s6j/fact
      nrb=idra(iir,iib,k2)
      npg=idra(iip,iig,k2)
      in=ntloc(npg,nrb,ipgsym,k2)
      isign=(-1)**(orbp+orbg+orba+orbb+orbr+orbs+k2)
      sum_local=sum_local+vabrs*vsgap*t(in)*fact*isign/denom
 42   continue
 41   continue
                                                                                
c     diagram A20
      call findk(orbp,orbb,orbs,orbg,iap,iab,ias,iag,jgot,jmax)
      if(jmax.eq.0)go to 43
      do 44 jloop=1,jmax
      k2=jgot(jloop)
      ak2=float(k2)
      call sixj(ak1,orbb,orbs,ak2,orbg,orbp,s6j)
      fact=(two*ak1+1)*(two*orbp+1)
      fact=s6j/fact
      npb=idra(iip,iib,k2)
      nsg=idra(iis,iig,k2)
      in=ntloc(npb,nsg,isgsym,k2)
      isign=(-1)**(orbp+orbg+orba+orbb+orbr+orbs+k1)
      sum_local=sum_local+vabrs*vrgap*t(in)*fact*isign/denom
 44   continue
 43   continue
                                                                                
c     diagram A21
      call findk(orbp,orbr,orba,orbg,iap,iar,iaa,iag,jgot,jmax)
      if(jmax.eq.0)go to 45
      do 46 jloop=1,jmax
      k2=jgot(jloop)
      ak2=float(k2)
      rk=vint(iir,iig,iip,iia,k2)
      d8=(-1)**(orbr+orbg+k2+1)
      vrgpa=rk*d8
      call sixj(ak1,orbg,orbp,ak2,orbr,orba,s6j)
      fact=(two*ak2+1)*(two*orbp+1)
      fact=s6j/fact
      npg=idra(iip,iig,k1)
      nsb=idra(iis,iib,k1)
      in=ntloc(npg,nsb,ipgsym,k1)
      isign=(-1)**(orbp+orbg+orba+orbb+orbr+orbs+k1)
      sum_local=sum_local+vabrs*vrgpa*t(in)*fact*isign/denom
 46   continue
 45   continue
                                                                                
c     diagram A22
      call findk(orbp,orbr,orba,orbg,iap,iar,iaa,iag,jgot,jmax)
      if(jmax.eq.0)go to 47
      do 48 jloop=1,jmax
      k2=jgot(jloop)
      ak2=float(k2)
      irgdis=idis(iir,iig)
      ipaloc=indxr(iip)+iia
      if(iig.gt.iir)ipaloc=indxr(iia)+iip
      rk=vint(iir,iig,iip,iia,k2)
      d8=(-1)**(orbr+orbg+k2+1)
      vrgpa=rk*d8
      call findk(orbp,orbb,orbs,orbg,iap,iab,ias,iag,igot,imax)
      if(imax.eq.0)go to 48
      do 49 iloop=1,imax
      k3=igot(iloop)
      ak3=float(k3)
      npb=idra(iip,iib,k3)
      nsg=idra(iis,iig,k3)
      in=ntloc(npb,nsg,isgsym,k3)
      call sixj(orbb,ak3,orbp,orbg,ak1,orbs,s6j1)
      call sixj(orbr,orbp,ak2,orbg,orba,ak1,s6j2)
      fact=s6j1*s6j2/(two*orbp+1)
      isign=(-1)**(orbp+orbg+orba+orbb+orbr+orbs+k1)
      sum_local=sum_local-vabrs*vrgpa*t(in)*fact*isign/denom
 49   continue
 48   continue
 47   continue
                                                                                
c     diagram A27
      if(orbg.eq.orba)then
      if(igsym.eq.iasym)then
                                                                                
      nrg=idra(iir,iig,k1)
      nsb=idra(iis,iib,k1)
      in=ntloc(nrg,nsb,isbsym,k1)
      call findk(orbr,orbb,orbs,orbg,iar,iab,ias,iag,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k3=jgot(jloop)
      if(k3.eq.0)then
      rk=vint(iip,iig,iip,iia,0)
      d8=(-1)**(orbp+orbg+1)
      vpgpa=rk*d8
      fact=dsqrt((2*orba+1)/(2*orbp+1))/((two*orbg+1)*(two*ak1+1))
      isign=(-1)**(orba+orbb+orbs+orbr)
      sum_local=sum_local-vabrs*vpgpa*t(in)*fact*isign/denom
                                                                                
c     diagram A28
      call findk(orbr,orbb,orbs,orbg,iar,iab,ias,iag,igot,imax)
      if(imax.eq.0)go to 50
      do 51 iloop=1,imax
      k2=igot(iloop)
      ak2=float(k2)
      nrb=idra(iir,iib,k2)
      nsg=idra(iis,iig,k2)
      in=ntloc(nrb,nsg,isgsym,k2)
      call sixj(ak1,orbs,orbb,ak2,orbr,orba,s6j)
      fact=s6j*dsqrt((2*orba+1)/(2*orbp+1))/(two*orbg+1)
      isign=(-1)**(orba+orbb+orbs+orbr+2*orbs)
      sum_local=sum_local+vabrs*vpgpa*t(in)*fact*isign/denom
 51   continue
 50   continue
                                                                                
      endif
      enddo
      endif
                                                                                
      endif
      endif
                                                                                
c     diagram A29
                                                                                
      call findk(orbr,orbp,orba,orbg,iar,iap,iaa,iag,igot,imax)
      if(imax.eq.0)go to 30
      do 52 iloop=1,imax
      k2=igot(iloop)
      ak2=float(k2)
      nsb=idra(iis,iib,k1)
      npg=idra(iip,iig,k1)
      ipgsym=mtbl(ipsym,igsym)
      in=ntloc(npg,nsb,ipgsym,k1)
      rk=vint(iir,iig,iip,iia,k2)
      d4=(-1)**(orbg+orbr+k2+1)
      vrgpa=rk*d4
      call sixj(ak2,orbg,orba,ak1,orbr,orbp,s6j)
      fact=one/((two*orbp+1)*(two*ak1+1))
      isign=(-1)**(orbs+orbb+k1+orbr+orba+orbg+orbp)
      fact=s6j*fact
      sum_local=sum_local+vabrs*vrgpa*t(in)*fact*isign/denom
                                                                                
c     diagram A30
      call findk(orbb,orbp,orbs,orbg,iab,iap,ias,iag,jgot,jmax)
      if(jmax.eq.0)go to 52
      do 53 jloop=1,jmax
      k3=jgot(jloop)
      ak3=float(k3)
      npb=idra(iip,iib,k3)
      nsg=idra(iis,iig,k3)
      isgsym=mtbl(issym,igsym)
      in=ntloc(npb,nsg,isgsym,k3)
      fact=one/(two*orbp+1)
      isign=(-1)**(orba+orbb+orbr+orbs+orbp+orbg+2*orbp+k1)
      call sixj(orbg,orba,ak2,orbr,orbp,ak1,s6j1)
      call sixj(orbg,ak3,orbs,orbb,ak1,orbp,s6j2)
      fact=fact*s6j1*s6j2
      sum_local=sum_local-vabrs*vrgpa*t(in)*fact*isign/denom
 53   continue
 52   continue
 30   continue
  3   continue
  2   continue
       call MPI_ALLREDUCE(sum_local,sum,1,MPI_REAL8,MPI_SUM,
     :                   MPI_COMM_WORLD,ierr)
       sum_local=0

                                                                                
      write(STDOUT,9999)sum
  1   continue
 9999 format(/,1x,'Contribution of approx. Triples to Heff=',F20.11,/)

      return
      end


c***************************************************************************
c                                                                          *
c      SUBROUTINE   SYMM                                                   *
c                                                                          *
c     THIS SUBROUTINE SETS UP THE T1 AND T2 EQUATION INDEX                 *
c                                                                          *
c***************************************************************************

      subroutine symm(iopt,nprint,ntmax,iproc)

      implicit real*8 (a-h,o-z)
      include 'mpif.h'

c     common block
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
      common/skip1/npbps(0:MXV,2)

c     local dimension
      dimension kgot(MXV+1),ieven(0:MXV),iodd(0:MXV)
      dimension indx(MDIM),indx1(MDIM)
      dimension ncheck(MDIM)

      do i=1, MDIM
         ncheck(i)=0
      enddo

c------------------ set the equation index for T1-----------------------------

      nsing=0

      if(iopt.ne.0)then

      if(iproc.eq.0)write(STDOUT,9999)
      if(iproc.eq.0)write(STDOUT,9997)

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

c------------------------Set the skip information--------------------------------

      do 2 i=0,MXV
      ieven(i)=0
      iodd(i)=0
      npbps(i,1)=0
      npbps(i,2)=0
      nnskip(i,1)=0
   2  nnskip(i,2)=0


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

c-----------------------iall=0 for only even rank----------------------------

      do 5 i=2,MXV,2
      nnskip(i,2)=nnskip(i-2,2)+ieven(i-2)*(ieven(i-2)+1)/2
  5   continue

      do  i=0,JMAX,2
      n2skip=n2skip+ieven(i)*(ieven(i)+1)/2
      enddo

      do 6 i=3,MXV-2,2
      nnskip(i,1)=nnskip(i-2,1)+iodd(i-2)*(iodd(i-2)+1)/2
  6   continue

      do 7 i=0,MXV,2
      npbps(i,2)=ieven(i)
  7   continue

      do 8 i=1,MXV,2
      npbps(i,1)=iodd(i)
  8   continue

      else

      do 13 i=1,MXV
      nnskip(i,2)=nnskip(i-1,2)+ieven(i-1)*(ieven(i-1)+1)/2
 13   continue

      do  i=0,MXV
      n2skip=n2skip+ieven(i)*(ieven(i)+1)/2
      enddo

      do 14 i=1,MXV
      nnskip(i,1)=nnskip(i-1,1)+iodd(i-1)*(iodd(i-1)+1)/2
 14   continue

      do 15 i=0,MXV
      npbps(i,2)=ieven(i)
 15   continue

      do 16 i=0,MXV
      npbps(i,1)=iodd(i)
 16   continue

      endif

      nn=nsing

      if(iproc.eq.0)write(STDOUT,9998)
      if(iproc.eq.0)write(STDOUT,9996)

c------------------------ set the equation index for T2--------------------------

      do 17 ia=1,nocc
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
      indx(nn)=nn
      indx1(nn)=neqn


      if(iproc.eq.0)then
      if (ncheck(neqn).eq.1) then
         write(STDOUT,*)'double count in symm =', neqn
      stop
      endif
      endif

      ncheck(neqn)=1

  45  continue
  40  continue
  30  continue
  20  continue
  17  continue

      ntmax=nn

9999  format(//,'SINGLE EXCITATIONS',/)
9998  format(//,'DOUBLE EXCITATIONS',/)
9997  format(3x,'IN','    IA ','  IP ','J(IA)',' J(IP)',/)
9996  format(3x,'IN','    IA ','  IB ','  IP ',' IQ' ,3x,'JA',
     $         4x,'JB',5x,'JP',5x,'JQ',5x,'ASM',2x,'BSM',2x,
     $         'PSM',2x,'QSM',2x,'KMAX',/)
9995  format(3i5,2(f5.1,2x))
9994  format(i7,4i5,4(f5.1,2x),8(i4,1x))

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

      DO I=1,N2INT
      DO J=0,MXV
      IGOT(J)=0
      JGOT(J)=0
      V2E(I,J)= ZERO

      ENDDO
      ENDDO

      CALL INFO_2E

      NBASIS=NOCC+NEXCIT

      if(itwoe.eq.1)goto 11
      if(itwoe.eq.2)goto 11
          GOTO 12

 11   CONTINUE

      open(NTWOE,file='twoea.dat',form='unformatted',
     &     status='unknown')

 12   CONTINUE

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
      VAL=0.0
      VAL1=0.0
      VAL2=0.0
      VAL3=0.0

c-------------------------- Coulomb interaction part----------------------------

      RK=SLATER_NEW(I,L,K,J,ORBI,ORBL,ORBK,ORBJ,
     :II,LL,KK,JJ,L1)

      IF(IALL.EQ.2) GOTO 220
      IF(IALL.EQ.3) GOTO 220
          GOTO 230

  220 CONTINUE

C      CALL BREIT(I,J,K,L,ORBI,ORBJ,ORBK,ORBL,
C     :II,JJ,KK,LL,L1,VAL,IALL)


c-----------------------Breit interaction part with n=k--------------------------


        VAL1=SLATER_NEWBR1(I,J,K,L,ORBI,ORBJ,ORBK,ORBL,
     : II,JJ,KK,LL,L1)


c-----------------------Breit interaction part with n=k-1------------------------


         VAL2=SLATER_NEWBR2(I,J,K,L,ORBI,ORBJ,ORBK,ORBL,
     : II,JJ,KK,LL,L1)

c-----------------------Breit interaction part with n=k+1------------------------


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

c*************************************************************************
c                                                                        *
c               FUNCTION VINT                                            *
c                                                                        *
c     This function stores the radial integral values for Coulomb and    *
c     Breit interaction                                                  *
c*************************************************************************      

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
      RETURN
      END



      integer function Calc4Binary(i, j, k, l)
c----------------------------------------------------------------------
      implicit none
      include 'mpif.h'
      integer i, j, k, l

      Calc4Binary = i*8 + j*4 + k*2 + l
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


      SUBROUTINE PACK(I,J,K,L,IJ,KL)
      IMPLICIT INTEGER (A-Z)
      common/info/nocc,nexcit,iType(MNBAS),iTCount(MNBAS)
      IJ=I*(I-1)/2+J
      KL=(NOCC+NEXCIT)*(K-1)+L
      RETURN
      END


      function nsloc(ipq,irg,ipar,kk)
      implicit real*8 (a-h,o-z)
      include 'mpif.h'
c     common block
      common/skip1/npbps(0:MXV,2)
      common/skip2/mmskip(0:MXV,2),m2skip
      common/skip3/nsing
      if(ipar.eq.1)then
      nsloc=(ipq-1)*npbps(kk,2)+irg+mmskip(kk,2)+nsing
      else
      nsloc=(ipq-1)*npbps(kk,1)+irg+mmskip(kk,1)+m2skip
      endif
      return
      end

      integer function idis(iip,iiq)
c     include 'mpif.h'
      ip=max0(iip,iiq)
      iq=min0(iip,iiq)
      idis=ip*(ip-1)/2+iq
      return
      end

      integer function indxr(ir)
c     include 'mpif.h'
      common/info/nocc,nexcit,iType(MNBAS),iTCount(MNBAS)
      nbasis=nocc+nexcit
      indxr=nbasis*(ir-1)
      return
      end


c------------------------------------------------------------------------------


      subroutine vbar(t,ivbar)


c----------------This subroutine calculate diagram associated with 'T'------
c-----------------amplitudes of the PPPP,PHPP,PHHP and PPHP type------------

      implicit real*8 (a-h,o-z)
      include 'mpif.h'
      real*8 sum_local
c     common block
      common/jvalue/orbc(MNOCC),orbe(MNEXC)
      common/index/ke(MNEXC),kc(MNOCC)
      common/info/nocc,nexcit,iType(MNBAS),iTCount(MNBAS)
      common/parity/isymc(MNOCC),isyme(MNEXC)
      common/symmetry/mtbl(MNS,MNS)
      common/kpavlaue/iqc(MNOCC),iqe(MNEXC)
      common/skip/nnskip(0:MXV,2),n2skip,nsing
      common/idra1/idra(MNBAS,MNBAS,0:MXV)
      common/idpa1/idpa(MNBAS,MNBAS)

      dimension t(MDIM)
      dimension kgot(MXV+1),igot(MXV+1),jgot(MXV+1)

      data zero,half,one,two/0.0,0.5,1.0,2.0/
      call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)


      do 1 i=0,MXV
      igot(i)=0
      jgot(i)=0
      kgot(i)=0
    1 continue

      if(ivbar.eq.1)goto 871
      if(ivbar.eq.2)goto 871
          goto 872

 871  continue

c----------------------Open MFVBAR file------------------------------------

      open(MFVBAR,file='vfbar.dat',form='unformatted',
     :status='unknown')
 872  continue


      if(ivbar.eq.0)goto 881
      if(ivbar.eq.1)goto 881

          goto 882

 881   continue

c----------------PP-PP BLOCK <pq|v|rs>--------------------------------------

      do 10 ip=1,nexcit
      ipsym=isyme(ip)
      iip=ke(ip)
      iap=iqe(ip)
      orbp=orbe(ip)

      do 80 iq=1,ip
      iqsym=isyme(iq)
      iiq=ke(iq)
      iaq=iqe(iq)
      orbq=orbe(iq)

      ipq=idis(iip,iiq)

      do 80 ir=1,nexcit
      irsym=isyme(ir)
      iir=ke(ir)
      iar=iqe(ir)
      orbr=orbe(ir)

      do 80 is=1,nexcit
      issym=isyme(is)
      iis=ke(is)
      ias=iqe(is)
      orbs=orbe(is)

      irs=indxr(iir)+iis
      ipqsym=mtbl(ipsym,iqsym)
      irssym=mtbl(irsym,issym)

      if(ipqsym.ne.irssym)go to 80

c------------------ diagram 1   bare  <pq|v|rs>-----------------------------

      call findk(orbp,orbr,orbq,orbs,iap,iar,iaq,ias,kgot,kmax)
      if(kmax.eq.0)go to 80
      do 90 kloop=1,kmax
      kk=kgot(kloop)
      ak=float(kk)
      rk=vint(iip,iiq,iir,iis,kk)
      d4=(-1)**(kk+orbp+orbq+1)
      vpqrs=rk*d4

      sum_local=zero

      do 100 ig=iproc+1,nocc,nproc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      iig=kc(ig)

      if(orbp.ne.orbg)go to 105
      if(ipsym.ne.igsym)go to 105
      in1=idpa(iip,iig)

c----------------------diagram   2    <qr|v|sg><p|t1|g>----------------------

      rk=vint(iiq,iig,iis,iir,kk)
      iqg=idis(iiq,iig)
      isr=indxr(iis)+iir
      if(iig.gt.iiq)isr=indxr(iir)+iis
      d4=(-1)**(kk+orbg+orbq+1)
      sum_local=sum_local-rk*d4*t(in1)

105   continue

c--------------------------Nonlinear diagram--------------------------------

      do 110 id=1,nocc
      idsym=isymc(id)
      orbd=orbc(id)
      iad=iqc(id)
      iid=kc(id)

      if(orbp.ne.orbg)go to 115
      if(ipsym.ne.igsym)go to 115
      in1=idpa(iip,iig)

      if(orbq.ne.orbd)go to 115
      if(iqsym.ne.idsym)go to 115
      in2=idpa(iiq,iid)

      rk=vint(iig,iid,iir,iis,kk)
      d4=(-1)**(kk+orbg+orbd+1)
      sum_local=sum_local+rk*d4*t(in1)*t(in2)

 115  continue

      igdsym=mtbl(igsym,idsym)
      if(irssym.ne.igdsym)go to 120

c----------------------- diagram 3    <rs|v|gd><pq|t2|gd>----------------------

      call findk(orbr,orbg,orbs,orbd,iar,iag,ias,iad,igot,imax)
      if(imax.eq.0)go to 120
      do 130 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      rk=vint(iig,iid,iir,iis,l1)
      d4=(-1)**(l1+orbg+orbd+1)

      call findk(orbp,orbg,orbq,orbd,iap,iag,iaq,iad,jgot,jmax)
      if(jmax.eq.0)go to 130
      do 140 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      npg=idra(iip,iig,l2)
      nqd=idra(iiq,iid,l2)
      ipgsym=mtbl(ipsym,igsym)
      in=ntloc(npg,nqd,ipgsym,l2)

      call sixj(ak,al1,al2,orbg,orbp,orbr,s6j1)
      call sixj(ak,al1,al2,orbd,orbq,orbs,s6j2)
      d5=(two*ak+1)*s6j1*s6j2
      d6=(-1)**(orbp+orbq+orbs+orbr)
      sum_local=sum_local+rk*d4*d5*d6*t(in)

 140  continue
 130  continue
 120  continue

 110  continue
 100  continue
       call MPI_ALLREDUCE(sum_local,sum,1,MPI_REAL8,MPI_SUM,
     :                   MPI_COMM_WORLD,ierr)
       sum_local=0

c---------------------- Calculate   <pq|V|rs> and write in MFBAR file------------

      val = sum+vpqrs
c     vb(ipq,irs,kk) = sum+vpqrs
      call PutVbar(ipq,irs,kk,val)
      if(ivbar.eq.1)then
      write(MFVBAR)val
      endif

  90  continue
  80  continue

c---------------------PH-PP BLOCK <pa|v|qr>---------------------------------

      do 210 ib=1,nocc
      ibsym=isymc(ib)
      iib=kc(ib)
      iab=iqc(ib)
      orbb=orbc(ib)

      ipb=idis(iip,iib)
      ipbsym=mtbl(ipsym,ibsym)

      do 210 iq=1,nexcit
      iqsym=isyme(iq)
      iiq=ke(iq)
      iaq=iqe(iq)
      orbq=orbe(iq)

      do 150 ir=1,nexcit
      irsym=isyme(ir)
      iir=ke(ir)
      iar=iqe(ir)
      orbr=orbe(ir)

      iqrsym=mtbl(iqsym,irsym)
      iqr=indxr(iiq)+iir
      if(iib.gt.iip)iqr=indxr(iir)+iiq
      if(ipbsym.ne.iqrsym)go to 150

c----------------------------- diagram 4      bare <pb|v|qr>----------------------

      call findk(orbp,orbq,orbr,orbb,iap,iaq,iar,iab,kgot,kmax)
      if(kmax.eq.0)go to 150
      do 160 kloop=1,kmax
      kk=kgot(kloop)
      ak=float(kk)
      rk=vint(iip,iib,iiq,iir,kk)
      d4=(-1)**(kk+orbp+orbb+1)
      vpbqr=rk*d4

      sum_local=zero

      do 170 ig=iproc+1,nocc,nproc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      iig=kc(ig)

      if(orbp.ne.orbg)go to 170
      if(ipsym.ne.igsym)go to 170
      in=idpa(iip,iig)

c-----------------------------diagram 5   <qb|v|gr><p|t1|g>----------------------

      rk=vint(iig,iib,iiq,iir,kk)
      d4=(-1)**(kk+orbg+orbb+1)
      sum_local=sum_local-rk*d4*t(in)

 170  continue
 161  continue
       call MPI_ALLREDUCE(sum_local,sum,1,MPI_REAL8,MPI_SUM,
     :                   MPI_COMM_WORLD,ierr)
       sum_local=0

c---------------------Calculate <pb|qr> and write in MFVBAR file----------------

      val = sum+vpbqr
c     vb(ipb,iqr,kk) =  sum+vpbqr
      call PutVbar(ipb,iqr,kk,val)
      if(ivbar.eq.1)then
      write(MFVBAR)val
      endif

 160  continue
 150  continue

c---------------------PH-PH BLOCK <pb|v|qa>------------------------------------

      do 210 ia=1,nocc
      iasym=isymc(ia)
      iia=kc(ia)
      iaa=iqc(ia)
      orba=orbc(ia)

      iqasym=mtbl(iqsym,iasym)
      iqa=indxr(iiq)+iia
      if(iib.gt.iip)iqa=indxr(iia)+iiq
      if(ipbsym.ne.iqasym)go to 215

c---------------------- diagram   6     bare <pa|v|qb>--------------------------

      call findk(orbp,orbq,orbb,orba,iap,iaq,iab,iaa,kgot,kmax)
      if(kmax.eq.0)go to 215
      do 220 kloop=1,kmax
      kk=kgot(kloop)
      ak=float(kk)
      rk=vint(iip,iib,iiq,iia,kk)
      ipb=idis(iip,iib)
      iqa=indxr(iiq)+iia
      if(iib.gt.iip)iqa=indxr(iia)+iiq
      d4=(-1)**(kk+orbp+orbb+1)
      vpbqa=rk*d4
      sum_local=zero


      do 230 ig=iproc+1,nocc,nproc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      iig=kc(ig)

      if(orbp.ne.orbg)go to 230
      if(ipsym.ne.igsym)go to 230
      in1=idpa(iip,iig)

c------------------------diagram  7    <qa|v|gb><p|t1|g>-----------------------

      rk=vint(iib,iig,iia,iiq,kk)
      d4=(-1)**(kk+orbg+orbb+1)
      sum_local=sum_local-rk*d4*t(in1)

230   continue

      do 250 ir=iproc+1,nexcit,nproc
      irsym=isyme(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      iir=ke(ir)

      if(orbr.eq.orba)then
      if(irsym.eq.iasym)then
      in2=idpa(iir,iia)

c-----------------------diagram   8   <pr|v|qb><r|t1|a>------------------------

      rk=vint(iip,iib,iiq,iir,kk)
      d4=(-1)**(kk+orbp+orbb+1)
      sum_local=sum_local+rk*d4*t(in2)
      endif
      endif

      do 250 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      iig=kc(ig)

c----------------------------- Non linear diagram----------------------------

      if(orbr.ne.orba)go to 240
      if(irsym.ne.iasym)go to 240
      in2=idpa(iir,iia)

      if(orbp.ne.orbg)go to 240
      if(ipsym.ne.igsym)go to 240
      in1=idpa(iip,iig)

      rk=vint(iib,iig,iir,iiq,kk)
      d4=(-1)**(kk+orbg+orbb+1)
      sum_local=sum_local-rk*d4*t(in1)*t(in2)

 240  continue

c---------------------------diagram  9   <qr|v|gb><pr|t2|ga>-------------------

      igrsym=mtbl(igsym,irsym)
      ipasym=mtbl(ipsym,iasym)
      if(ipasym.ne.igrsym)go to 250

      call findk(orbq,orbg,orbr,orbb,iaq,iag,iar,iab,igot,imax)
      if(imax.eq.0)go to 250
      do 260 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      rk=vint(iib,iig,iir,iiq,l1)
      d4=(-1)**(l1+orbb+orbg+1)

      call findk(orbp,orbg,orbr,orba,iap,iag,iar,iaa,jgot,jmax)
      if(jmax.eq.0)go to 260
      do 270 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      npg=idra(iip,iig,l2)
      nra=idra(iir,iia,l2)
      ipgsym=mtbl(ipsym,igsym)
      in=ntloc(npg,nra,ipgsym,l2)

      call sixj(ak,al1,al2,orbg,orbp,orbq,s6j1)
      call sixj(ak,al2,al1,orbr,orbb,orba,s6j2)
      d5=(two*ak+1)*s6j1*s6j2
      d6=(-1)**(orba+orbp+orbq+orbb+kk+l1+l2)
      sum_local=sum_local-rk*d4*d5*d6*t(in)

 270  continue
 260  continue
 250  continue
       call MPI_ALLREDUCE(sum_local,sum,1,MPI_REAL8,MPI_SUM,
     :                   MPI_COMM_WORLD,ierr)
       sum_local=0


c-------------Calculate total <pa|V|qb> and write in MFBAR-------------------

      val = sum+vpbqa
c     vb(ipb,iqa,kk) = sum+vpbqa
      call PutVbar(ipb,iqa,kk,val)
      if(ivbar.eq.1)then
      write(MFVBAR)val
      endif

 220  continue
 215  continue

c-------------------PH-HP BLOCK <pb|v|aq>-----------------------------------

      ipbsym=mtbl(ipsym,ibsym)
      iaqsym=mtbl(iasym,iqsym)
      iaqloc=indxr(iia)+iiq
      if(iib.gt.iip)iaqloc=indxr(iiq)+iia
      if(ipbsym.ne.iaqsym)go to 280

c------------------------- diagram   10   bare  <pq|v|ab>--------------------

      call findk(orbp,orba,orbq,orbb,iap,iaa,iaq,iab,kgot,kmax)
      if(kmax.eq.0)go to 280
      do 290 kloop=1,kmax
      kk=kgot(kloop)
      ak=float(kk)
      rk=vint(iip,iib,iia,iiq,kk)
      d4=(-1)**(kk+orbp+orbb+1)
      vpbaq=rk*d4

      sum_local=zero

      do 300 ir=iproc+1,nexcit,nproc
      irsym=isyme(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      iir=ke(ir)

   
c----------------------- diagram   11    <pq|v|rb><r|t1|a>--------------------

      if(orbr.ne.orba)go to 300
      if(irsym.ne.iasym)go to 300
      in1=idpa(iir,iia)
      rk=vint(iip,iib,iir,iiq,kk)
      d4=(-1)**(kk+orbp+orbb+1)
      sum_local=sum_local+rk*d4*t(in1)

300   continue


      do 330 ig=iproc+1,nocc,nproc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      iig=kc(ig)

c----------------------diagram    12    <qa|v|bg><p|t1|g>---------------------

      if(orbg.eq.orbp)then
      if(ipsym.eq.igsym)then
      in2=idpa(iip,iig)
      rk=vint(iib,iig,iiq,iia,kk)
      d4=(-1)**(kk+orbb+orbg+1)
      sum_local=sum_local-rk*d4*t(in2)

      endif
      endif

      do 330 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      iir=ke(ir)

c-----------------Non linear diagram---------------------------------------

      if(orbg.ne.orbp)go to 310
      if(ipsym.ne.igsym)go to 310
      in2=idpa(iip,iig)

      if(orbr.ne.orba)go to 310
      if(irsym.ne.iasym)go to 310
      in1=idpa(iir,iia)

      rk=vint(iib,iig,iiq,iir,kk)
      d4=(-1)**(kk+orbg+orbb+1)
      sum_local=sum_local-rk*d4*t(in1)*t(in2)

 310  continue



      irgsym=mtbl(irsym,igsym)
      ipasym=mtbl(ipsym,iasym)
      if(ipasym.ne.irgsym)go to 330

c----------------------diagram   13    <rb|v|gq><ra|t2|gp>------------------

      rk=vint(iib,iig,iiq,iir,kk)
      d4=(-1)**(kk+orbb+orbg+1)
      d5=(-1)**(orbg+orbr+kk)
      npa=idra(iip,iia,kk)
      nrg=idra(iir,iig,kk)
      in=ntloc(npa,nrg,ipasym,kk)
      d6=one/(two*kk+1)
      sum_local=sum_local-rk*d4*d5*d6*t(in)

c---------------------- diagram    14   <rb|v|gq><pa|t2|gr>-----------------

      call findk(orbp,orbg,orbr,orba,iap,iag,iar,iaa,jgot,jmax)
      if(jmax.eq.0)go to 340
      do 350 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      npg=idra(iip,iig,l2)
      nra=idra(iir,iia,l2)
      ipgsym=mtbl(ipsym,igsym)
      in=ntloc(npg,nra,ipgsym,l2)
      call sixj(ak,orbp,orba,al2,orbr,orbg,s6j)
      d6=s6j
      sum_local=sum_local-rk*d4*d5*d6*t(in)

 350  continue
 340  continue

c----------------------diagram  15     <rq|v|bg><pa|t2|gr>---------------------

      call findk(orbr,orbb,orbq,orbg,iar,iab,iaq,iag,jgot,jmax)
      if(jmax.eq.0)go to 360
      do 370 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      rk=vint(iib,iig,iir,iiq,l2)
      d4=(-1)**(l2+orbb+orbg+1)
      npa=idra(iip,iia,kk)
      nrg=idra(iir,iig,kk)
      in=ntloc(npa,nrg,ipasym,kk)
      call sixj(ak,orbg,orbr,al2,orbb,orbq,s6j)
      d6=s6j
      sum_local=sum_local-rk*d4*d5*d6*t(in)

 370  continue
 360  continue
 330  continue

       call MPI_ALLREDUCE(sum_local,sum,1,MPI_REAL8,MPI_SUM,
     :                   MPI_COMM_WORLD,ierr)
       sum_local=0

c--------------Calculate total <pq|V|ab>  and write in MFVBAR-------------------

      val = sum+vpbaq
c     vb(ipb,iaqloc,kk) = sum+vpbaq
      call PutVbar(ipb,iaqloc,kk,val)
      if(ivbar.eq.1)then
      write(MFVBAR)val
      endif

 290  continue
 280  continue

 210  continue


c--------------------------PP-HP BLOCK <pr|v|qa>-------------------------------

      do 380 ir=1,nexcit
      irsym=isyme(ir)
      iir=ke(ir)
      iar=iqe(ir)
      orbr=orbe(ir)

      ipr=idis(iip,iir)

      do 380 iq=1,nexcit
      iqsym=isyme(iq)
      iiq=ke(iq)
      iaq=iqe(iq)
      orbq=orbe(iq)

      do 380 ia=1,nocc
      iasym=isymc(ia)
      iia=kc(ia)
      iaa=iqc(ia)
      orba=orbc(ia)

      iprsym=mtbl(ipsym,irsym)
      iqasym=mtbl(iqsym,iasym)

      iqaloc=indxr(iiq)+iia
      if(iir.gt.iip)iqaloc=indxr(iia)+iiq
      if(iprsym.ne.iqasym)go to 380

c---------------diagram   16      bare <pr|v|qa>------------------------------

      call findk(orbp,orbq,orbr,orba,iap,iaq,iar,iaa,kgot,kmax)
      if(kmax.eq.0)go to 380
      do 390 kloop=1,kmax
      kk=kgot(kloop)
      ak=float(kk)
      rk=vint(iip,iir,iiq,iia,kk)
      d4=(-1)**(kk+orbp+orbr+1)
      vprqa=rk*d4

c     sum=zero
       sum_local=0

      do 430 is=iproc+1,nexcit,nproc
      issym=isyme(is)
      orbs=orbe(is)
      ias=iqe(is)
      iis=ke(is)

      do 430 ib=1,nocc
      ibsym=isymc(ib)
      orbb=orbc(ib)
      iab=iqc(ib)
      iib=kc(ib)

      irasym=mtbl(irsym,iasym)
      isbsym=mtbl(issym,ibsym)
      if(irasym.ne.isbsym)go to 400

c------------------diagram   17   <ps|v|qb><rb|t2|as>------------------------

      rk=vint(iip,iib,iiq,iis,kk)
      d4=(-1)**(kk+orbp+orbb+1)
      d5=(-1)**(orbb+orbs+kk)
      nra=idra(iir,iia,kk)
      nsb=idra(iis,iib,kk)
      irasym=mtbl(irsym,iasym)
      in=ntloc(nra,nsb,irasym,kk)
      d6=one/(two*kk+1)
      sum_local=sum_local-rk*d4*d5*d6*t(in)


c--------------------- diagram   18  <ps|v|qb><ra|t2|bs>---------------------- 

      call findk(orbr,orbb,orbs,orba,iar,iab,ias,iaa,igot,imax)
      if(imax.eq.0)go to 400
      do 410 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      nrb=idra(iir,iib,l2)
      nsa=idra(iis,iia,l2)
      irbsym=mtbl(irsym,ibsym)
      in=ntloc(nrb,nsa,irbsym,l2)
      call sixj(ak,orbs,orbb,al2,orbr,orba,s6j)
      d6=s6j
      sum_local=sum_local-rk*d4*d5*d6*t(in)

 410  continue

c--------------------- diagram   19   <pq|v|sb><rs|t2|ab>--------------------

      call findk(orbq,orbb,orbp,orbs,iaq,iab,iap,ias,igot,imax)
      if(imax.eq.0)go to 400
      do 420 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      nra=idra(iir,iia,kk)
      nsb=idra(iis,iib,kk)
      irasym=mtbl(irsym,iasym)
      in=ntloc(nra,nsb,irasym,kk)
      rk=vint(iip,iib,iis,iiq,l2)
      d4=(-1)**(l2+orbp+orbb+1)
      d5=(-1)**(orbs+orbb+kk)
      call sixj(ak,orbb,orbs,al2,orbp,orbq,s6j)
      d6=s6j
      sum_local=sum_local-rk*d4*d5*d6*t(in)

 420  continue
 400  continue

      ipasym=mtbl(ipsym,iasym)
      isbsym=mtbl(issym,ibsym)
      if(ipasym.ne.isbsym)go to 425

c--------------------diagram 20    <rq|v|sb><ps|t2|ba>------------------------

      call findk(orbr,orbs,orbb,orbq,iar,ias,iab,iaq,igot,imax)
      if(imax.eq.0)go to 425
      do 440 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      rk=vint(iir,iib,iis,iiq,l2)
      irb=idis(iir,iib)
      d4=(-1)**(l2+orbr+orbb+1)

      call findk(orba,orbs,orbb,orbp,iaa,ias,iab,iap,jgot,jmax)
      if(jmax.eq.0)go to 440
      do 450 jloop=1,jmax
      l3=jgot(jloop)
      al3=float(l3)
      npb=idra(iip,iib,l3)
      nsa=idra(iis,iia,l3)
      isasym=mtbl(issym,iasym)
      in=ntloc(npb,nsa,isasym,l3)

      call sixj(ak,al3,al2,orbs,orbr,orba,s6j1)
      call sixj(ak,al2,al3,orbb,orbp,orbq,s6j2)
      d5=(-1)**(orba+orbp+orbq+orbr+kk+l2+l3)
      d6=s6j1*s6j2*(2*kk+1)
      sum_local=sum_local-rk*d4*d5*d6*t(in)

 450  continue
 440  continue
 425  continue

c---------------------------- Non linear diagram----------------------------

      do ig=1,nocc
      igsym=isymc(ig)
      iig=kc(ig)
      iag=iqc(ig)
      orbg=orbc(ig)

      if(ipsym.eq.igsym)then
      if(orbp.eq.orbg)then
      in1=idpa(iip,iig)

      irasym=mtbl(irsym,iasym)
      isbsym=mtbl(issym,ibsym)
      if(irasym.eq.isbsym)then
      nra=idra(iir,iia,kk)
      nsb=idra(iis,iib,kk)
      in=ntloc(nra,nsb,isbsym,kk)
      rk=vint(iib,iig,iis,iiq,kk)
      d4=(-1)**(orbb+orbg+kk+1)
      fact=(-1)**(orbb+orbs+kk)/(2*ak+1)
      sum_local=sum_local+rk*d4*fact*t(in)*t(in1)

      irbsym=mtbl(irsym,ibsym)

      call findk(orbb,orbr,orbs,orba,iab,iar,ias,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nsa=idra(iis,iia,l1)
      nrb=idra(iir,iib,l1)
      in=ntloc(nsa,nrb,irbsym,l1)
      call sixj(ak,orbs,orbb,al1,orbr,orba,s6j)
      fact=(-1)**(orbb+orbs+kk)*s6j
      sum_local=sum_local+rk*d4*fact*t(in)*t(in1)

      enddo
      endif
c
      nra=idra(iir,iia,kk)
      nsb=idra(iis,iib,kk)
      in=ntloc(nra,nsb,isbsym,kk)
      call findk(orbb,orbq,orbs,orbg,iab,iaq,ias,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      rk=vint(iib,iig,iiq,iis,l1)
      d4=(-1)**(orbb+orbg+l1+1)
      call sixj(ak,orbb,orbs,al1,orbp,orbq,s6j)
      fact=s6j*(-1)**(orbs+orbb+kk)
      sum_local=sum_local+rk*d4*fact*t(in)*t(in1)

      enddo
      endif
      endif
      endif
      endif

c
      if(irsym.eq.igsym)then
      if(orbr.eq.orbg)then
      in1=idpa(iir,iig)

      ipbsym=mtbl(ipsym,ibsym)
      isasym=mtbl(issym,iasym)
      if(ipbsym.eq.isasym)then

      call findk(orbb,orbp,orbs,orba,iab,iap,ias,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npb=idra(iip,iib,l1)
      nsa=idra(iis,iia,l1)
      in=ntloc(npb,nsa,ipbsym,l1)

      call findk(orbb,orbq,orbs,orbg,iab,iaq,ias,iag,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      rk=vint(iib,iig,iiq,iis,l2)
      d4=(-1)**(orbb+orbg+l2+1)

      call sixj(al2,al1,ak,orbp,orbq,orbb,s6j1)
      call sixj(al2,al1,ak,orba,orbr,orbs,s6j2)
      d5=(-1)**(orba+orbp+orbq+orbr+kk+l2+l1)
      d6=s6j1*s6j2*(2*kk+1)
      sum_local=sum_local+rk*d4*d5*d6*t(in1)*t(in)

      enddo
      endif
      enddo
      endif

      endif
      endif
      endif
c
      if(issym.eq.iasym)then
      if(orbs.eq.orba)then
      in1=idpa(iis,iia)

      ipbsym=mtbl(ipsym,ibsym)
      irgsym=mtbl(irsym,igsym)
      if(ipbsym.eq.irgsym)then

      call findk(orbb,orbp,orbr,orbg,iab,iap,iar,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npb=idra(iip,iib,l1)
      nrg=idra(iir,iig,l1)
      in=ntloc(npb,nrg,ipbsym,l1)

      call findk(orbb,orbq,orbs,orbg,iab,iaq,ias,iag,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      rk=vint(iib,iig,iiq,iis,l2)
      d4=(-1)**(orbb+orbg+l2+1)
      call sixj(al2,al1,ak,orbp,orbq,orbb,s6j1)
      call sixj(al2,al1,ak,orbr,orba,orbg,s6j2)
      d5=s6j1*s6j2*(2*kk+1)
      d6=(-1)**(orbp+orbq+orbr+orba)
      sum_local=sum_local+rk*d4*d5*d6*t(in)*t(in1)

      enddo
      endif
      enddo
      endif

      endif
      endif
      endif

      enddo

 430  continue

      do 460 ib=iproc+1,nocc,nproc
      ibsym=isymc(ib)
      orbb=orbc(ib)
      iab=iqc(ib)
      iib=kc(ib)


      if(orbb.eq.orbr)then
      if(ibsym.eq.irsym)then
      in1=idpa(iir,iib)

c--------------------------- diagram   21    <pb|v|qa><r|f|b>-----------------

      d4=(-1)**(kk+orbp+orbb+1)
      rk=vint(iip,iib,iiq,iia,kk)
      sum_local=sum_local-rk*d4*t(in1)
      endif
      endif

c---------------------------diagram    22   <qr|v|ba><p|t1|b>------------------

      if(ipsym.eq.ibsym)then
      if(orbp.eq.orbb)then
      in1=idpa(iip,iib)
      rk=vint(iir,iib,iia,iiq,kk)
      d4=(-1)**(kk+orbr+orbb+1)
      sum_local=sum_local-rk*d4*t(in1)
      endif
      endif

460   continue


c     diagram 7
      do 470 is=iproc+1,nexcit,nproc
      issym=isyme(is)
      orbs=orbe(is)
      ias=iqe(is)
      iis=ke(is)

c---------------------- diagram    23    <ps|v|qr><s|t1|a>-------------------

      if(orbs.eq.orba)then
      if(issym.eq.iasym)then
      in2=idpa(iis,iia)
      rk=vint(iip,iir,iiq,iis,kk)
      d4=(-1)**(kk+orbp+orbr+1)
      sum_local=sum_local+rk*d4*t(in2)
      endif
      endif

      do 470 ib=1,nocc
      ibsym=isymc(ib)
      orbb=orbc(ib)
      iab=iqc(ib)
      iib=kc(ib)

      if(orbs.eq.orba)then
      if(issym.eq.iasym)then
      in2=idpa(iis,iia)

c------------------------- Non linear diagram---------------------------

      if(orbb.eq.orbr)then
      if(ibsym.eq.irsym)then
      in1=idpa(iir,iib)

      rk=vint(iip,iib,iiq,iis,kk)
      d4=(-1)**(kk+orbp+orbb+1)
      sum_local=sum_local-rk*d4*t(in1)*t(in2)

      endif
      endif
      endif
      endif

      if(orbs.eq.orba)then
      if(issym.eq.iasym)then
      in2=idpa(iis,iia)

      if(orbb.eq.orbp)then
      if(ibsym.eq.ipsym)then
      in1=idpa(iip,iib)

      rk=vint(iir,iib,iis,iiq,kk)
      d4=(-1)**(kk+orbr+orbb+1)
      sum_local=sum_local-rk*d4*t(in1)*t(in2)

      endif
      endif
      endif
      endif

 470  continue


      do 480 ib=iproc+1,nocc,nproc
      ibsym=isymc(ib)
      orbb=orbc(ib)
      iab=iqc(ib)
      iib=kc(ib)

      do 480 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      iig=kc(ig)

      if(ipsym.eq.ibsym)then
      if(irsym.eq.igsym)then
      if(orbr.eq.orbg)then
      if(orbp.eq.orbb)then
      in1=idpa(iip,iib)
      in2=idpa(iir,iig)

      rk=vint(iib,iig,iiq,iia,kk)
      d4=(-1)**(kk+orbb+orbg+1)
      sum_local=sum_local+rk*d4*t(in1)*t(in2)

      endif
      endif
      endif
      endif

      ipbsym=mtbl(ipsym,ibsym)
      irgsym=mtbl(irsym,igsym)
      if(ipbsym.ne.irgsym)go to 485

c---------------------diagram  24    <qa|v|bg><pr|t2|bg>---------------------

      call findk(orba,orbg,orbb,orbq,iaa,iag,iab,iaq,igot,imax)
      if(imax.eq.0) go to 485
      do 490 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      rk=vint(iib,iig,iiq,iia,l2)
      d4=(-1)**(l2+orbb+orbg+1)

      call findk(orbr,orbg,orbb,orbp,iar,iag,iab,iap,jgot,jmax)
      if(jmax.eq.0) go to 490
      do 500 jloop=1,jmax
      l3=jgot(jloop)
      al3=float(l3)
      npb=idra(iip,iib,l3)
      nrg=idra(iir,iig,l3)
      in=ntloc(npb,nrg,irgsym,l3)

      call sixj(ak,al2,al3,orbg,orbr,orba,s6j1)
      call sixj(ak,al2,al3,orbb,orbp,orbq,s6j2)
      d5=s6j1*s6j2*(2*kk+1)
      d6=(-1)**(orbp+orbq+orbr+orba)
      sum_local=sum_local+rk*d4*d5*d6*t(in)

 500  continue
 490  continue
 485  continue

      go to 999

      do 510 is=1,nexcit
      issym=isyme(is)
      orbs=orbe(is)
      ias=iqe(is)
      iis=ke(is)

c     new diagram (13)

      if(ipsym.eq.ibsym)then
      if(irsym.eq.igsym)then
      if(issym.eq.iasym)then
      if(orbr.eq.orbg)then
      if(orbp.eq.orbb)then
      if(orbs.eq.orba)then
      in1=idpa(iip,iib)
      in2=idpa(iir,iig)
      in3=idpa(iis,iia)
      rk=vint(iib,iig,iiq,iis,kk)
      d4=(-1)**(kk+orbb+orbg+1)
      sum_local=sum_local+rk*d4*t(in1)*t(in2)*t(in3)

      endif
      endif
      endif
      endif
      endif
      endif

 510  continue

 999  continue

 480  continue
       call MPI_ALLREDUCE(sum_local,sum,1,MPI_REAL8,MPI_SUM,
     :                   MPI_COMM_WORLD,ierr)
       sum_local=0


c---------------Calculate total <pr|V|qa> and write in MFVBAR------------------

      val = sum+vprqa
      call PutVbar(ipr,iqaloc,kk,val)
      if(ivbar.eq.1)then
      write(MFVBAR)val
      endif

 390  continue
 380  continue
  10  continue

 882  continue

      if(ivbar.eq.2)then

c-------------------  Read  MFBAR  file----------------------------

      do 1000 ip=1,nexcit
      ipsym=isyme(ip)
      iip=ke(ip)
      iap=iqe(ip)
      orbp=orbe(ip)

      do 1080 iq=1,ip
      iqsym=isyme(iq)
      iiq=ke(iq)
      iaq=iqe(iq)
      orbq=orbe(iq)

      ipq=idis(iip,iiq)

      do 1080 ir=1,nexcit
      irsym=isyme(ir)
      iir=ke(ir)
      iar=iqe(ir)
      orbr=orbe(ir)

      do 1080 is=1,nexcit
      issym=isyme(is)
      iis=ke(is)
      ias=iqe(is)
      orbs=orbe(is)

      irs=indxr(iir)+iis
      ipqsym=mtbl(ipsym,iqsym)
      irssym=mtbl(irsym,issym)

      if(ipqsym.ne.irssym)go to 1080

      call findk(orbp,orbr,orbq,orbs,iap,iar,iaq,ias,kgot,kmax)
      if(kmax.eq.0)go to 1080
      do 1090 kloop=1,kmax
      kk=kgot(kloop)

      read(MFVBAR)val
      call PutVbar(ipq,irs,kk,val)

 1090 continue
 1080 continue

      do 1210 ib=1,nocc
      ibsym=isymc(ib)
      iib=kc(ib)
      iab=iqc(ib)
      orbb=orbc(ib)

      ipb=idis(iip,iib)

      ipbsym=mtbl(ipsym,ibsym)

      do 1210 iq=1,nexcit
      iqsym=isyme(iq)
      iiq=ke(iq)
      iaq=iqe(iq)
      orbq=orbe(iq)

      do 1150 ir=1,nexcit
      irsym=isyme(ir)
      iir=ke(ir)
      iar=iqe(ir)
      orbr=orbe(ir)

      iqrsym=mtbl(iqsym,irsym)

      iqr=indxr(iiq)+iir
      if(iib.gt.iip)iqr=indxr(iir)+iiq
      if(ipbsym.ne.iqrsym)go to 1150

      call findk(orbp,orbq,orbr,orbb,iap,iaq,iar,iab,kgot,kmax)
      if(kmax.eq.0)go to 1150
      do 1160 kloop=1,kmax
      kk=kgot(kloop)

      read(MFVBAR)val
      call PutVbar(ipb,iqr,kk,val)

 1160 continue
 1150 continue

      do 1210 ia=1,nocc
      iasym=isymc(ia)
      iia=kc(ia)
      iaa=iqc(ia)
      orba=orbc(ia)

      iqasym=mtbl(iqsym,iasym)
      iqa=indxr(iiq)+iia
      if(iib.gt.iip)iqa=indxr(iia)+iiq
      if(ipbsym.ne.iqasym)go to 1215

      call findk(orbp,orbq,orbb,orba,iap,iaq,iab,iaa,kgot,kmax)
      if(kmax.eq.0)go to 1215
      do 1220 kloop=1,kmax
      kk=kgot(kloop)

      read(MFVBAR)val
      call PutVbar(ipb,iqa,kk,val)

 1220 continue
 1215 continue

      ipbsym=mtbl(ipsym,ibsym)
      iaqsym=mtbl(iasym,iqsym)
      iaqloc=indxr(iia)+iiq
      if(iib.gt.iip)iaqloc=indxr(iiq)+iia
      if(ipbsym.ne.iaqsym)go to 1280
      call findk(orbp,orba,orbq,orbb,iap,iaa,iaq,iab,kgot,kmax)
      if(kmax.eq.0)go to 1280
      do 1290 kloop=1,kmax
      kk=kgot(kloop)

      read(MFVBAR)val
      call PutVbar(ipb,iaqloc,kk,val)

 1290 continue
 1280 continue

 1210 continue


c     PP-HP BLOCK <pr|v|qa>
      do 1380 ir=1,nexcit
      irsym=isyme(ir)
      iir=ke(ir)
      iar=iqe(ir)
      orbr=orbe(ir)
      ipr=idis(iip,iir)
      do 1380 iq=1,nexcit
      iqsym=isyme(iq)
      iiq=ke(iq)
      iaq=iqe(iq)
      orbq=orbe(iq)
      do 1380 ia=1,nocc
      iasym=isymc(ia)
      iia=kc(ia)
      iaa=iqc(ia)
      orba=orbc(ia)
      iprsym=mtbl(ipsym,irsym)
      iqasym=mtbl(iqsym,iasym)

      iqaloc=indxr(iiq)+iia
      if(iir.gt.iip)iqaloc=indxr(iia)+iiq

      if(iprsym.ne.iqasym)go to 1380
      call findk(orbp,orbq,orbr,orba,iap,iaq,iar,iaa,kgot,kmax)
      if(kmax.eq.0)go to 1380
      do 1390 kloop=1,kmax
      kk=kgot(kloop)

      read(MFVBAR)val
      call PutVbar(ipr,iqaloc,kk,val)

 1390 continue
 1380 continue
 1000 continue

      endif

      return
      end


      subroutine fbar(t,ivbar)

c***************This subroutine calculate diagrmas associated with **************
c***************the contraction between H and cluster operator T ****************
c***************to produce PP,HH and PH type diagrams **************************

      implicit real*8 (a-h,o-z)
      include 'mpif.h'
      real*8 sum_local
      common/jvalue/orbc(MNOCC),orbe(MNEXC)
      common/index/ke(MNEXC),kc(MNOCC)
      common/info/nocc,nexcit,iType(MNBAS),iTCount(MNBAS)
      common/parity/isymc(MNOCC),isyme(MNEXC)
      common/symmetry/mtbl(MNS,MNS)
      common/kpavlaue/iqc(MNOCC),iqe(MNEXC)
      common/f_bar/f(MNBAS,MNBAS)
      common/orbital_energy/eorb(MNBAS)
      common/idra1/idra(MNBAS,MNBAS,0:MXV)
      common/idpa1/idpa(MNBAS,MNBAS)
      common/skip/nnskip(0:MXV,2),n2skip,nsing

      dimension kgot(MXV+1),igot(MXV+1),jgot(MXV+1)
      dimension t(MDIM)

      data zero,half,one,two/0.0,0.5,1.0,2.0/
      call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)

      do i=0,MXV
      kgot(i)=0
      igot(i)=0
      enddo

      nbasis=nocc+nexcit

      do i=1,nbasis
      do j=1,nbasis
      f(i,j)=zero
      enddo
      enddo

      rewind NTFILE
      read(NTFILE)n,(t(i),i=1,n)                                        

      if(ivbar.eq.1)goto 871
      if(ivbar.eq.2)goto 871
          goto 872

 871  continue

c-------------------Open the NFVBAR file-----------------------------


      open(NFVBAR,file='fvbar.dat',form='unformatted',
     :status='unknown')

 872  continue

      if(ivbar.eq.0)goto 881
      if(ivbar.eq.1)goto 881

          goto 882

 881   continue

c-------------------H-H type <b|f|a>--------------------------------

      do 1 ia=1,nocc
      iasym=isymc(ia)
      iia=kc(ia)
      orba=orbc(ia)
      iaa=iqc(ia)

      do 1 ib=1,nocc
      ibsym=isymc(ib)
      iib=kc(ib)
      iab=iqc(ib)
      orbb=orbc(ib)

      if(ibsym.ne.iasym)go to 1
      if(orba.ne.orbb)go to 1
      if(iia.eq.iib)then

c------------------------diagram 1   bare <b|f|a>---------------------

      f(iib,iia)=f(iib,iia)+eorb(iia)
      else
      f(iib,iia)=0
      endif

      sum_local=zero


      do 7 ip=iproc+1,nexcit,nproc
      ipsym=isyme(ip)
      iip=ke(ip)
      orbp=orbe(ip)
      iap=iqe(ip)

      do 7 ig=1,nocc
      igsym=isymc(ig)
      iig=kc(ig)
      iag=iqc(ig)
      orbg=orbc(ig)

      if(orbp.ne.orbg)go to 2
      if(ipsym.ne.igsym)go to 2
      in=idpa(iip,iig)

c----------------diagram  2   <pa|v|gb><p|t1|g>----------------------


      rk=vint(iib,iig,iia,iip,0)
      ibgdis=idis(iib,iig)
      iaploc=indxr(iia)+iip
      if(iig.gt.iib)iaploc=indxr(iip)+iia
      d4=(-1)**(orbb+orbg+1)
      d5=dsqrt((two*orbg+one)/(two*orba+one))
      sum_local=sum_local+rk*t(in)*d4*d5

c ------------------diagram 3    <pa|v|bg><p|t1|g>-------------------

      call findk (orba,orbg,orbp,orbb,iaa,iag,iap,iab,igot,imax)
      if(imax.eq.0)go to 2
      do 3 iloop=1,imax
      kk=igot(iloop)
      ak=float(kk)
      rk=vint(iib,iig,iip,iia,kk)
      ibgdis=idis(iib,iig)
      ipaloc=indxr(iip)+iia
      if(iig.gt.iib)ipaloc=indxr(iia)+iip
      d4=(-1)**(kk+orbb+orbg+1)
      d5=(-1)**(orba+orbp+kk)
      d6=one/(two*orba+one)
      sum_local=sum_local+rk*d4*d5*d6*t(in)

   3  continue
   2  continue


      do 7 iq=1,nexcit
      iqsym=isyme(iq)
      iiq=ke(iq)
      iaq=iqe(iq)
      orbq=orbe(iq)

      ipgsym=mtbl(ipsym,igsym)
      ipgbsm=mtbl(ipgsym,ibsym)
      if(ipgbsm.ne.iqsym)go to 4

      call findk(orbp,orbb,orbq,orbg,iap,iab,iaq,iag,jgot,
     $jmax)
      if(jmax.eq.0)go to 4
      do 5 jloop=1,jmax
      kk=jgot(jloop)
      ak=float(kk)
      rk=vint(iib,iig,iip,iiq,kk)
      ibgdis=idis(iib,iig)
      ipqloc=indxr(iip)+iiq
      if(iig.gt.iib)ipqloc=indxr(iiq)+iip
      d4=(-1)**(orbb+orbg+kk+1)
      twoe=rk*d4

c-------------------------diagram 4  <pq|v|bg><pq|t2|ag>----------------------

      iapsym=mtbl(iasym,ipsym)
      npa=idra(iip,iia,kk)
      nqg=idra(iiq,iig,kk)
      in=ntloc(npa,nqg,iapsym,kk)
      d5=(-1)**(orbb+orbp+orbq+orbg)
      d6=one/((two*orba+one)*(two*ak+one))
      sum_local=sum_local+twoe*d5*d6*t(in)

c-------------------------diagram 5    <pq|v|bg><pq|t2|ga>--------------------

      call findk(orbp,orbg,orbq,orba,iap,iag,iaq,iaa,igot,
     $imax)
      if(imax.eq.0)go to 5
      do 6 iloop=1,imax
      ll=igot(iloop)
      al=float(ll)
      iaqsym=mtbl(iasym,iqsym)
      npg=idra(iip,iig,ll)
      nqa=idra(iiq,iia,ll)
      in=ntloc(npg,nqa,iaqsym,ll)
      call sixj(ak,orbq,orbg,al,orbp,orba,s6j)
      d6=one/(two*orba+1)*s6j
      sum_local=sum_local+twoe*d5*d6*t(in)

  6   continue
  5   continue
  4   continue

c --------------------------t1*t1 contribution---------------------------

      if(orbg.ne.orbp)go to 7
      if(orba.ne.orbq)go to 7
      if(iasym.ne.iqsym)go to 7
      if(igsym.ne.ipsym)go to 7
      in1=idpa(iiq,iia)

c     new diagram
      in2=idpa(iip,iig)
      rk=vint(iib,iig,iiq,iip,0)
      d4=(-1)**(orbg+orbb+1)
      twoe=rk*d4
      fact=dsqrt((2*orbg+1)/(2*orba+1))
      sum_local=sum_local+twoe*fact*t(in1)*t(in2)

      call findk(orbb,orbp,orbg,orbq,iab,iap,iag,iaq,igot,imax)
      if(imax.eq.0)go to 7
      do 8 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      rk=vint(iib,iig,iip,iiq,l1)
      d4=(-1)**(orbg+orbb+l1+1)
      twoe=rk*d4
      in1=idpa(iiq,iia)
      in2=idpa(iip,iig)
      fact=one/(two*orba+1)
      isign=(-1)**(orba+orbp+l1)
      sum_local=sum_local+twoe*fact*isign*t(in1)*t(in2)

  8   continue
  7   continue
       call MPI_ALLREDUCE(sum_local,sum,1,MPI_REAL8,MPI_SUM,
     :                   MPI_COMM_WORLD,ierr)
       sum_local=0


c------ Calculate total <b|f|a> and  write in NFVBAR file---------------------


      f(iib,iia)=f(iib,iia)+sum
      if(ivbar.eq.1)then
      write(NFVBAR)f(iib,iia)
      endif

  1   continue

c-----------------P-P type <p|f|q>----------------------------------

      do 10 ip=1,nexcit
      ipsym=isyme(ip)
      iip=ke(ip)
      orbp=orbe(ip)
      iap=iqe(ip)

      do 10 iq=1,nexcit
      iqsym=isyme(iq)
      iiq=ke(iq)
      iaq=iqe(iq)
      orbq=orbe(iq)

      if(iqsym.ne.ipsym)go to 10
      if(orbp.ne.orbq)go to 10
      if(iip.eq.iiq)then

c------------------diagram 6     bare  <p|f|q>---------------------

      f(iip,iiq)=f(iip,iiq)+eorb(iip)
      else
      f(iip,iiq)=0
      endif

      sum_local=zero

      do 70 ir=iproc+1,nexcit,nproc
      irsym=isyme(ir)
      iir=ke(ir)
      orbr=orbe(ir)
      iar=iqe(ir)

      do 70 ia=1,nocc
      iasym=isymc(ia)
      iia=kc(ia)
      iaa=iqc(ia)
      orba=orbc(ia)

      if(orbr.ne.orba)go to 20
      if(irsym.ne.iasym)go to 20
      in=idpa(iir,iia)

c-------------------diargam 7  <pr|v|qa><r|t1|a>----------------------

      rk=vint(iip,iia,iiq,iir,0)
      d4=(-1)**(orbp+orba+1)
      d5=dsqrt((two*orba+one)/(two*orbp+one))
      sum_local=sum_local+rk*t(in)*d4*d5

c------------------- diagram 8  <pq|v|ra><r|t1|a>---------------------

      call findk(orbp,orbr,orbq,orba,iap,iar,iaq,iaa,igot,imax)
      if(imax.eq.0)go to 20
      do 30 iloop=1,imax
      kk=igot(iloop)
      ak=float(kk)
      rk=vint(iip,iia,iir,iiq,kk)
      d4=(-1)**(kk+orbp+orba+1)
      d5=one/(two*orbp+one)
      d6=(-1)**(orbp+orba+kk)
      sum_local=sum_local+rk*t(in)*d4*d5*d6

 30   continue
 20   continue

      do 70 ib=1,nocc
      ibsym=isymc(ib)
      iib=kc(ib)
      iab=iqc(ib)
      orbb=orbc(ib)

      irasym=mtbl(irsym,iasym)
      iraqsm=mtbl(irasym,iqsym)
      if(iraqsm.ne.ibsym)go to 40

      call findk(orbq,orba,orbr,orbb,iaq,iaa,iar,iab,igot,
     $imax)
      if(imax.eq.0)go to 40
      do 50 iloop=1,imax
      kk=igot(iloop)
      ak=float(kk)
      rk=vint(iia,iib,iiq,iir,kk)
      iabdis=idis(iia,iib)
      iqrloc=indxr(iiq)+iir
      if(iib.gt.iia)iqrloc=indxr(iir)+iiq
      d4=(-1)**(orba+orbb+kk+1)
      twoe=rk*d4

c------------------ diagram 9 <qr|v|ab><pr|t2|ab>-----------------------

      d5=(-1)**(orbq+orba+orbb+orbr)
      d6=one/((two*orbp+one)*(two*ak+one))
      npa=idra(iip,iia,kk)
      nrb=idra(iir,iib,kk)
      iapsym=mtbl(iasym,ipsym)
      in=ntloc(npa,nrb,iapsym,kk)
      sum_local=sum_local-twoe*d5*d6*t(in)

c------------------diagram 10 <qr|v|ab><pr|t2|ba>-----------------------

      call findk(orbp,orbb,orbr,orba,iap,iab,iar,iaa,jgot,
     $jmax)
      if(jmax.eq.0)go to 50
      do 60 jloop=1,jmax
      ll=jgot(jloop)
      al=float(ll)
      npb=idra(iip,iib,ll)
      nra=idra(iir,iia,ll)
      ibpsym=mtbl(ibsym,ipsym)
      in=ntloc(npb,nra,ibpsym,ll)
      call sixj(ak,orbr,orbb,al,orbp,orba,s6j)
      d6=(one)/(two*orbp+1)*s6j
      sum_local=sum_local-twoe*d5*d6*t(in)

  60  continue
  50  continue
  40  continue

c--------------------------t1*t1 contribution--------------------------

      if(orbb.ne.orbp)go to 70
      if(orba.ne.orbr)go to 70
      if(ibsym.ne.ipsym)go to 70
      if(iasym.ne.irsym)go to 70

c     new diagram
      in1=idpa(iir,iia)
      in2=idpa(iip,iib)
      rk=vint(iia,iib,iir,iiq,0)
      d4=(-1)**(orba+orbb+1)
      fact=dsqrt((2*orba+1)/(2*orbp+1))
      twoe=rk*d4
      sum_local=sum_local-twoe*fact*t(in1)*t(in2)

      call findk(orbb,orbr,orba,orbq,iab,iar,iaa,iaq,igot,imax)
      if(imax.eq.0)go to 70
      do 75 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      rk=vint(iia,iib,iiq,iir,l1)
      d4=(-1)**(orba+orbb+l1+1)
      twoe=rk*d4
      in1=idpa(iip,iib)
      in2=idpa(iir,iia)
      fact=one/(two*orbp+1)
      isign=(-1)**(orba+orbp+l1)
      sum_local=sum_local-twoe*fact*isign*t(in1)*t(in2)

  75  continue
  70  continue
       call MPI_ALLREDUCE(sum_local,sum,1,MPI_REAL8,MPI_SUM,
     :                   MPI_COMM_WORLD,ierr)
       sum_local=0


c--------- Calculate total <p|f|q> and write in NFVBAR file-------------------

      f(iip,iiq)=f(iip,iiq)+sum
      if(ivbar.eq.1)then
      write(NFVBAR)f(iip,iiq)
      endif

  10  continue

c--------------------P-H type <a|f|p>---------------------------

      do 80 ia=1,nocc
      iasym=isymc(ia)
      iia=kc(ia)
      orba=orbc(ia)
      iaa=iqc(ia)

      do 80 ip=1,nexcit
      ipsym=isyme(ip)
      iip=ke(ip)
      orbp=orbe(ip)
      iap=iqe(ip)

      if(ipsym.ne.iasym)go to 80
      if(orbp.ne.orba)go to 80

      sum_local=0

      do 90 iq=iproc+1,nexcit,nproc
      iqsym=isyme(iq)
      iiq=ke(iq)
      iaq=iqe(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)

      do 90 ib=1,nocc
      ibsym=isymc(ib)
      iib=kc(ib)
      iab=iqc(ib)
      orbb=orbc(ib)
      iab=iqc(ib)

      if(iqsym.ne.ibsym)go to 90
      if(orbq.ne.orbb)go to 90
      in=idpa(iiq,iib)

c--------------------------diagram  11  <pq|v|ab><q|t1|b>-------------------

      rk=vint(iia,iib,iip,iiq,0)
      d4=(-1)**(orba+orbb+1)
      d5=dsqrt((two*orbb+1)/(two*orba+1))
      sum_local=sum_local+rk*d4*d5*t(in)

c-------------------------diagram   12   <pq|v|ba><q|t1|b>------------------

      call findk(orbq,orba,orbp,orbb,iaq,iaa,iap,iab,igot,imax)
      if(imax.eq.0)go to 90
      do 100 iloop=1,imax
      kk=igot(iloop)
      ak=float(kk)
      rk=vint(iia,iib,iiq,iip,kk)
      iabdis=idis(iia,iib)
      iqploc=indxr(iiq)+iip
      if(iib.gt.iia)iqploc=indxr(iip)+iiq
      d4=(-1)**(orba+orbb+1+kk)
      d5=one/(two*orba+1)
      d6=(-1)**(orba+orbq+kk)
      sum_local=sum_local+rk*d4*d5*d6*t(in)

  100 continue
   90 continue
       call MPI_ALLREDUCE(sum_local,sum,1,MPI_REAL8,MPI_SUM,
     :                   MPI_COMM_WORLD,ierr)
       sum_local=0


c----------------Calculate total <a|f|p> and write in NFVBAR file-------------

      f(iia,iip)=sum
      if(ivbar.eq.1)then
      write(NFVBAR)f(iia,iip)
      endif

   80 continue

 882  continue

c---------------------------- Read NFVBAR file------------------------------

      if(ivbar.eq.2)then
      do 101 ia=1,nocc
      iasym=isymc(ia)
      iia=kc(ia)
      orba=orbc(ia)
      iaa=iqc(ia)

      do 101 ib=1,nocc
      ibsym=isymc(ib)
      iib=kc(ib)
      iab=iqc(ib)
      orbb=orbc(ib)

      if(ibsym.ne.iasym)go to 101
      if(orba.ne.orbb)go to 101
      read(NFVBAR)f(iib,iia)

 101  continue

      do 110 ip=1,nexcit
      ipsym=isyme(ip)
      iip=ke(ip)
      orbp=orbe(ip)
      iap=iqe(ip)

      do 110 iq=1,nexcit
      iqsym=isyme(iq)
      iiq=ke(iq)
      iaq=iqe(iq)
      orbq=orbe(iq)

      if(iqsym.ne.ipsym)go to 110
      if(orbp.ne.orbq)go to 110
      read(NFVBAR)f(iip,iiq)

  110 continue

      do 180 ia=1,nocc
      iasym=isymc(ia)
      iia=kc(ia)
      orba=orbc(ia)
      iaa=iqc(ia)

      do 180 ip=1,nexcit
      ipsym=isyme(ip)
      iip=ke(ip)
      orbp=orbe(ip)
      iap=iqe(ip)

      if(ipsym.ne.iasym)go to 180
      if(orbp.ne.orba)go to 180
      read(NFVBAR)f(iia,iip)

  180 continue
      endif

c     write(STDOUT,9999)
c     write(STDOUT,9998)((f(i,j),i=1,nbasis),j=1,nbasis)
 9999 format(//,1x,' FBAR ',//)
 9998 format(4F20.11)

      return
      end


      subroutine jacobi(itol,nprint,irows,irowe,isolv,nsmax,
     :iopt,iproc,iorb,t3opt,t,dspl,delrow,drow,sold,s,soname)
      implicit real*8 (a-h,o-z)
      include 'mpif.h'
      integer dspl,delrow,error,drow
      integer t3opt
      character*14 soname
      common/scratch1/b(MDIM),ad(MDIM)
      common/scratch2/raux(MDIM),caux(MDIM)
      
      dimension t(MDIM),a(MNEXC,MNEXC)
      dimension dspl(NPMAX),delrow(NPMAX)
      dimension sold(MDIM),s(MDIM)


      open(42,file='ccea.eps',form='formatted',status='unknown')
      open(NSFILE,file=soname,form='unformatted',status='unknown')

      if(iproc.eq.0) then
      write(STDOUT,*)'------------------------------------------------'
      write(STDOUT,*)'        Valence orbital No.', iorb
      write(STDOUT,*)'------------------------------------------------'
      endif
      nit=0
c------------------------------------------------------------------------
c     Begin of iteration
c------------------------------------------------------------------------
   10 nit=nit+1
c     start_time = MPI_WTIME()

      call eamatb(iopt,iorb,nsmax,isolv,irows,irowe,sold,t)
      call heff(iorb,sold,a,iproc)
      call nlheff(iorb,a,sold,t,nsmax,t3opt,isolv,irows,irowe,
     $iopt,iproc)

c     check convergence

      do 20 i=irows,irowe
   20 s(i)=(b(i)-raux(i)+caux(i))/ad(i)

c     call MPI_ALLGATHERV(s(irows),drow,MPI_DOUBLE_PRECISION,s,
c    >delrow,dspl,MPI_DOUBLE_PRECISION, MPI_COMM_WORLD,error)

      eps=0
      do 30 i=1,nsmax
      eps=eps+dabs(sold(i)-s(i))
   30 sold(i)=s(i)      
      conv=eps/nsmax

      rewind(NSFILE)
      if(iproc.eq.0)write(NSFILE)nsmax,(s(i),i=1,nsmax)

C     close(NSFILE)

      if(iproc.eq.0)write(42,*)nit,conv
      if(iproc.eq.0)write(STDOUT,9997)nit,conv

c     end_time= MPI_WTIME()
c     if(iproc.eq.0)write(STDOUT,1002)end_time-start_time

      if(conv.le.SMALL)go to 40
      if(nit.ge.NITMAX)go to 50
      go to 10

c------------------------------------------------------------------------
c     End of iteration
c------------------------------------------------------------------------

   50 if(iproc.eq.0)write(STDOUT,9999)
      go to 60
   40 if(iproc.eq.0)write(STDOUT,9998)nit
   60 continue
      
      call heff(iorb,s,a,iproc)

      if(t3opt.ne.0)then
c     start_time = MPI_WTIME()
      call ccsdtea(iorb,s,t,sum,iproc)

c     end_time= MPI_WTIME() 
      if(iproc.eq.0)write(STDOUT,1002)end_time-start_time
      endif   


 998  format(/,1x,'diagonal part of b-matrix',//)
 999  format(4f20.11)
1000  format(/,1x,'Jacobi iteration begins',//)
1001  format(/,1x,'solution vectors',//)
1002  format(/,1x,'Time elapsed=',2f20.11//)
9990  format(//,1x,'End of linear CC',//)
9994  format(/,1x,'Total energy(including triples)=',F20.11/)
9995  format(4f20.11)
9991  format(//,1x,'Non-linear Iteration begins',//)
9996  format(//,1x,'Iteration begins',//)
9997  format(1x,'Iteration no=',i4,5x,'eps=',f20.11)
9998  format(1x,' Iteration converges at',I4,' cycle')
9999  format(1x,' Convergence failure; Programme terminates')

      return
      end


c***********************************************************************
c              SUBROUTINE   EALOC                                      *
c                                                                      *
c     This subroutine sets the equation  index for the open shell      *
c      cluster amplitudes                                               *
c                                                                      *
c***********************************************************************


      subroutine ealoc(iorb,nsmax,nprint,isolv,iopt,iproc)

      implicit real*8 (a-h,o-z)
      include 'mpif.h'
      common/iallparity/iall
      common/jvalue/orbc(MNOCC),orbe(MNEXC)
      common/index/ke(MNEXC),kc(MNOCC)
      common/info/nocc,nexcit,iType(MNBAS),iTCount(MNBAS)
      common/parity/isymc(MNOCC),isyme(MNEXC)
      common/symmetry/mtbl(MNS,MNS)
      common/kpavlaue/iqc(MNOCC),iqe(MNEXC)
      common/idra1/idra(MNBAS,MNBAS,0:MXV)
      common/idpa1/idpa(MNBAS,MNBAS)
      common/jdra1/jdra(MNBAS,MNBAS,0:MXV)
      common/jdpa1/jdpa(MNBAS,MNBAS)
      common/maxj/jmax
      common/skip1/npbps(0:MXV,2)
      common/skip2/mmskip(0:MXV,2),m2skip
      common/skip3/ncis
c     local dimension
      dimension kgot(MXV+1),ieven(0:MXV),iodd(0:MXV)
      dimension npqbps(0:MXV,2)
      dimension indx(MDIM),indx1(MDIM)
      dimension ncheck(MDIM)
      data zero/0.0d0/

      nbasis=nocc+nexcit

c      if(iproc.eq.0)write(STDOUT,9991)
c      if(iproc.eq.0)write(STDOUT,9990)
     
      do i=1, MDIM
         ncheck(i)=0
      enddo

      do i=0,MXV
        kgot(i)=0
        ieven(i)=0
        iodd(i)=0
        npqbps(i,1)=0
        mmskip(i,1)=0
        npqbps(i,2)=0
        mmskip(i,2)=0
      enddo

c--------------------------Set S1 equation index--------------------------------

      ncis=0

      if(iopt.ne.0)then

      do 1 ip=1,nexcit
      ipsym=isyme(ip)
      iip=ke(ip)
      orbp=orbe(ip)

      do 2 iq=1,nexcit
      iqsym=isyme(iq)
      iiq=ke(iq)
      orbq=orbe(iq)

      if(iip.ne.iorb)go to 2
      if(iip.eq.iiq)go to 2
      if(orbp.ne.orbq)go to 2
      if(ipsym.ne.iqsym)go to 2

      ncis=ncis+1
      jdpa(iip,iiq)=ncis
      indx(ncis)=ncis
      indx1(ncis)=ncis

  2   continue
  1   continue

      endif

c----------------------Set the skip information--------------------------------

      do 3 ip=1,nexcit
      ipsym=isyme(ip)
      iip=ke(ip)
      orbp=orbe(ip)
      iap=iqe(ip)

      do 4 iq=1,nexcit
      iqsym=isyme(iq)
      iiq=ke(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)

      if(iip.ne.iorb)go to 4

      ipqsym=mtbl(ipsym,iqsym)

      call findk(orbp,orbq,orbp,orbq,iap,iaq,iap,iaq,kgot,
     $kmax)
      if(kmax.eq.0)go to 4
      do 5 kloop=1,kmax
      kk=kgot(kloop)
      if(ipqsym.eq.1)then
      ieven(kk)=ieven(kk)+1
      jdra(iip,iiq,kk)=ieven(kk)
      else
      iodd(kk)=iodd(kk)+1
      jdra(iip,iiq,kk)=iodd(kk)
      endif

  5   continue
  4   continue
  3   continue

      m2skip=0

      if(iall.eq.0)then

      do 6 i=2,MXV,2
      npqbps(i,2)=npqbps(i-2,2)+ieven(i-2)*npbps(i-2,2)
      mmskip(i,2)=npqbps(i,2)
  6   continue

      do  i=0,JMAX,2
      m2skip=m2skip+ieven(i)*npbps(i,2)
      enddo

      do 7 i=3,MXV-2,2
      npqbps(i,1)=npqbps(i-2,1)+iodd(i-2)*npbps(i-2,1)
      mmskip(i,1)=npqbps(i,1)
  7   continue

      else

      do 10 i=1,MXV
      npqbps(i,2)=npqbps(i-1,2)+ieven(i-1)*npbps(i-1,2)
      mmskip(i,2)=npqbps(i,2)
 10   continue

      do  i=0,MXV
      m2skip=m2skip+ieven(i)*npbps(i,2)
      enddo

      do 11 i=1,MXV
      npqbps(i,1)=npqbps(i-1,1)+iodd(i-1)*npbps(i-1,1)
      mmskip(i,1)=npqbps(i,1)
 11   continue

      endif

      m2skip=m2skip+ncis
      
c      if(iproc.eq.0)write(STDOUT,9999)
c      if(iproc.eq.0)write(STDOUT,9998)

c-----------------------Set S2 equation index----------------------------------

      nn=ncis

      do 15 ip=1,nexcit
      ipsym=isyme(ip)

      do 20 iq=1,nexcit
      iqsym=isyme(iq)

      ipqsym=mtbl(ipsym,iqsym)

      do 30 ir=1,nexcit
      irsym=isyme(ir)

      ipqrsm=mtbl(ipqsym,irsym)

      do 40 ig=1,nocc
      igsym=isymc(ig)

      if(ipqrsm.ne.igsym)go to 40

      orbg=orbc(ig)
      orbp=orbe(ip)
      orbq=orbe(iq)
      orbr=orbe(ir)

      iap=iqe(ip)
      iag=iqc(ig)
      iaq=iqe(iq)
      iar=iqe(ir)

      iip=ke(ip)
      iiq=ke(iq)
      iir=ke(ir)
      iig=kc(ig)

      if(iip.ne.iorb)go to 40

      call findk(orbp,orbq,orbr,orbg,iap,iaq,iar,iag,kgot,
     $kmax)
      if(kmax.eq.0)go to 40
      do 50 nloop=1,kmax
      kk=kgot(nloop)
      npq=jdra(iip,iiq,kk)
      nrg=idra(iir,iig,kk)
      nn=nn+1
      neqn=nsloc(npq,nrg,ipqsym,kk)
      indx(nn)=nn
      indx1(nn)=neqn

      

      if (ncheck(neqn).eq.1) then
        if(iproc.eq.0)write(STDOUT,*)'double count in ealoc =', neqn
      stop
      endif

      ncheck(neqn)=1

  50  continue
  40  continue
  30  continue
  20  continue
  15  continue

      nsmax=nn

9999  format(//,'DOUBLE EXCITATIONS',/)
9998  format(3x,'IN','   IIP ',' IIQ ',' IIS ',' IIA',
     $1x,'JP',3X,'JQ',3x,'JR',3x,'JA',3x,'PSM',3x,'QSM',3x,
     $'RSM',3x,'ASM',3x,'KMAX',/)
9997  format(5i5,4f5.1,6i6)
9996  format(/,2x,'COUPLED CLUSTER A MATRIX (IP)',/)
9995  format(4f20.11)
9991  format(//,'SINGLE EXCITATIONS',/)
9990  format(3x,'IN','   IIP ',' IIQ',' J(IIP)',' J(IIQ)'/)
9989  format(3i5,2(f5.1,2x))

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

      if(i.eq.1)then
      rk=dfloat((kp1-k-1)*(kp2-k-1))/dfloat((k+1)*(2*k+1))
      endif

      if(i.eq.2)then
      rk=dfloat((kp1-k-1)*(kp2+k+1))/dfloat((k+1)*(2*k+1))
      endif

      if(i.eq.3)then
      rk=dfloat((kp1+k+1)*(kp2-k-1))/dfloat((k+1)*(2*k+1))
      endif

      if(i.eq.4)then
      rk=dfloat((kp1+k+1)*(kp2+k+1))/dfloat((k+1)*(2*k+1))
      endif

      endif

      if(k.eq.l1)then

      rk=-dfloat((kapa+kapc)*(kapb+kapd))/dfloat(k*(k+1))

      endif

      if(k.eq.l1+1)then

      if(i.eq.1)then
      rk=dfloat((kp1+k)*(kp2+k))/dfloat(k*(2*k+1))
      endif

      if(i.eq.2)then
      rk=dfloat((kp1+k)*(kp2-k))/dfloat(k*(2*k+1))
      endif

      if(i.eq.3)then
      rk=dfloat((kp1-k)*(kp2+k))/dfloat(k*(2*k+1))
      endif

      if(i.eq.4)then
      rk=dfloat((kp1-k)*(kp2-k))/dfloat(k*(2*k+1))
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

c********************************************************************
c                                                                   *
c              FUNCTION  SELECT1                                    *
c                                                                   *
c   This function select the multipoles of even operator            *
c                                                                   *
c********************************************************************


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

c********************************************************************
c                                                                   *
c                FUNCTION SELECT                                    *
c                                                                   *
c     This subroutine select the multipoles of odd operators        *
c                                                                   *
c********************************************************************


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


