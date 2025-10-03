C
C  
C
C  Program for computing breast cancer penetrance using family data.
C  TP53 data for  Clare Turnbull.
C  This version assumes  fixed background incidences independent of birth cohort.
C  Breast Cancer incidences for England and Wales


      PARAMETER(LENC=100,LENI=100000,LENL=100,LENR=1000000,MXLOCI=10
     :,NEXTRA=100)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION EXTRA(NEXTRA),RARRAY(LENR)
      INTEGER IARRAY(LENI),COND
      CHARACTER*8 CARRAY(LENC),LNAME(MXLOCI),BASE,BATFIL*40,LOCFIL*40
     :,MUTLOC,OUTFIL*40,PEDFIL*40,TITLE*40,TRAVEL,XXSIGN,XYSIGN
      LOGICAL LARRAY(LENL),ASYCV,ECHO,LUMP,ORDERD,PMODE,STAND
C
C     DEFAULT VALUES FOR SOME VARIABLES.  TO RUN THE PROGRAM IN
C     BATCH MODE FILL IN THE NAME OF THE BATCH FILE BATFIL.  SET
C     PMODE TO TRUE TO DO CALCULATIONS IN PRODUCT MODE AND LUMP
C     TO TRUE TO AMALGAMATE ALLELES AT EACH LOCUS.
C
      DATA MXTWIN,ABSENT/10,-1.0D20/
      DATA CONV,NCONV,MXSTEP,DP/1.0D-7,4,3,1.0D-7/
      DATA BATFIL,LUMP,ORDERD,PMODE/' ',.FALSE.,.FALSE.,.FALSE./
C
C     DEFAULT VALUES FOR THE PROBLEM MENU.
C
      DATA TITLE,LOCFIL/'TP53_data','locusbr.dat'/
      DATA PEDFIL,OUTFIL/'pedigree.txt','single_out.dat'/
      DATA ECHO,XXSIGN,XYSIGN,NVAR/.FALSE.,'F','M',7/
      DATA EXTRA/100*0.0D0/
      DATA LNAME/'MAJOR',9*' '/
      DATA MUTLOC,XXRATE,XYRATE,COND/' ',0.0001199d0,0.0001199d0,2/
      DATA BASE,STAND,TRAVEL,NPOINT/'E',.FALSE.,'SEARCH',1/
      DATA NPAR,NCNSTR,ASYCV,MXITER/3,0,.TRUE.,599/
C
C     THE ARRAYS AND VARIABLES BEGIN A LONG DESCENT INTO THE
C     PROGRAM.  SAY GOODBYE TO THEM AND WISH THEM LUCK.
C
      CALL MENDEL(EXTRA,RARRAY,IARRAY,CARRAY,LNAME,LARRAY,ABSENT
     :,CONV,DP,XXRATE,XYRATE,COND,LENC,LENI,LENL,LENR,MXITER,MXLOCI
     :,MXSTEP,MXTWIN,NCNSTR,NCONV,NEXTRA,NPAR,NPOINT,NVAR,BASE,BATFIL
     :,LOCFIL,MUTLOC,OUTFIL,PEDFIL,TITLE,TRAVEL,XXSIGN,XYSIGN,ASYCV
     :,ECHO,LUMP,ORDERD,PMODE,STAND)
      END


      SUBROUTINE INITAL(ALLFRQ,CNSTR,CVALUE,EXTRA,GRID,PAR,PARMAX
     1,PARMIN,PNAME,XLINK,XXRATE,XYRATE,MAXALL,MUTATE,NCNSTR,NEXTRA
     2,NLOCI,NPAR,NPOINT,NVAR,PROBLM,UNIT3,TRAVEL)
C
C     IN THIS SUBROUTINE THE USER SHOULD DEFINE THE INITIAL
C     PARAMETER VALUES, THE PARAMETER BOUNDS, AND THE LINEAR
C     EQUALITY CONSTRAINTS FOR A LIKELIHOOD SEARCH.  WHEN A
C     GRID OF LIKELIHOOD VALUES IS DESIRED, THEN ONLY DEFINE
C     THE ARRAY GRID.  PARAMETER NAMES CAN BE OPTIONALLY INPUT.
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION ALLFRQ(NLOCI,MAXALL),CNSTR(NCNSTR,NPAR)
     1,CVALUE(NCNSTR),EXTRA(NEXTRA),GRID(NPOINT,NPAR),PAR(NPAR)
     2,PARMAX(NPAR),PARMIN(NPAR)
      INTEGER PROBLM,UNIT3
      CHARACTER*8 PNAME(NPAR),TRAVEL
      LOGICAL XLINK(NLOCI)

      double precision   popbr(0:89)

      common /popincid/popbr



C-----------------------------------------------------------------------
C The log-likelihood is parameterised in terms of the log relative risks
C Assume all individuals are censored at 70
C-----------------------------------------------------------------------



      pname(1)='RR'

      CALL RANDOM_SEED()

      DO I=1,NPAR
        par(I) = 0.5d0
        parmin(I) =  0.0d0
        parmax(I) =  5.0d0
        IF ( par(I) .EQ. 0.0d0 ) THEN
            CALL RANDOM_NUMBER(HARVEST)
            par(I) = parmin(I) + HARVEST * (parmax(I) - parmin(I))
        ENDIF
      END DO


C----------------------------------------------------------------
C  Read in the population incidence rates for England and Wales 1993-97
c  per 100000 population.
C----------------------------------------------------------------
      
      data  popbr/0.0010, 0.0050, 0.0090, 0.0220, 0.0430,
     :0.0650, 0.0860, 0.1160, 0.1990, 0.2910,
     :0.3830, 0.4750, 0.5730, 0.7090, 0.8510,
     :0.9930, 1.1350, 1.2730, 1.3820, 1.4870,
     :1.5920, 1.6960, 1.8240, 2.0860, 2.3720,
     :2.6570, 2.9430, 3.2620, 3.7800, 4.3320,
     :4.8840, 5.4360, 6.0360, 6.9270, 7.8660,
     :8.8060, 9.7450, 10.7850, 12.4340, 14.1830,
     :15.9330, 17.6820, 19.5540, 22.1620, 24.8920,
     :27.6220, 30.3510, 33.5010, 39.1800, 45.2780,
     :51.3740, 57.4690, 62.9320, 64.6190, 65.6750,
     :66.7300, 67.7830, 69.2050, 72.8420, 76.8460,
     :80.8460, 84.8440, 89.0500, 94.5240, 100.2040,
     :105.8790, 111.5490, 117.1350, 122.2510, 127.2830,
     :132.3070, 137.3240, 142.7000, 150.2750, 158.2030,
     :166.1160, 174.0110, 181.8510, 189.4490, 196.9850,
     :204.4950, 211.9770, 219.2250, 225.2180, 230.9740,
     :236.6940, 242.3780, 246.8160, 243.9830, 239.9320 /
      



      END

      SUBROUTINE OUTPUT(EXTRA,PAR,SCORE,PNAME,LOGLIK,FINAL,ITER,MAXPAR
     1,NEXTRA,NPAR,NSTEP,UNIT3,BASE,TRAVEL,STAND,UMOVE)
C
C     THIS SUBROUTINE OUTPUTS THE LOGLIKELIHOOD AND PARAMETERS
C     AT EACH ITERATION.
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION EXTRA(NEXTRA),PAR(MAXPAR),SCORE(MAXPAR),LOGLIK
      INTEGER FINAL,UNIT3
      CHARACTER*8 PNAME(MAXPAR),BASE,TRAVEL
      LOGICAL STAND,UMOVE
C

      SAVE BEST,IBEST,START
C
      IF (ITER.EQ.1) THEN
      BEST=LOGLIK
      START=LOGLIK
      write(*,10) (pname(i),i=1,npar)
      WRITE(UNIT3,10) (PNAME(I),I=1,NPAR)
 10   FORMAT(/,' ITER  NSTEP  LOGLIKELIHOOD',(T28,4(4X,A8),:))
      END IF
      IF (LOGLIK.GE.BEST) THEN
      BEST=LOGLIK
      IBEST=ITER
      END IF
      IF (STAND) LOGLIK=LOGLIK-START
      IF (BASE.EQ.'10') LOGLIK=LOG10(EXP(1.0D0))*LOGLIK
      write(*,20) iter,nstep,loglik,(par(i),i=1,npar)
      WRITE(UNIT3,20) ITER,NSTEP,LOGLIK,(PAR(I),I=1,NPAR)
 20   FORMAT(/,I4,3X,I3,3X,D14.7,(T28,4(1X,D11.4),:))


      IF (ITER.EQ.FINAL) WRITE(UNIT3,30) IBEST
 30   FORMAT(/,' THE MAXIMUM LOGLIKELIHOOD OCCURS AT ITERATION',I4,'.')
      END



      SUBROUTINE NEWLIK(EXTRA,PAR,LOGLIK,COND,ITER,NEXTRA,NPAR,NPED,PED)
C
C     THIS SUBROUTINE ALLOWS THE USER TO MODIFY THE LOGLIKELIHOOD
C     OF A PARTICULAR PEDIGREE.  FOR INSTANCE, IN GENETIC COUNSELING
C     PROBLEMS IT IS MAY BE NECESSARY TO FORM CONDITIONAL PROBABILITIES.
C     COND TELLS WHICH PEDIGREE, IF ANY, TO CONDITION ON.
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION EXTRA(NEXTRA),PAR(NPAR),LOGLIK
      INTEGER COND,PED
C
C     THE LOGLIKELIHOOD FOR EACH PEDIGREE GIVEN FREE RECOMBINATION
C     IS STORED IN THE ARRAY EXTRA.  THE LIKELIHOOD IS THEN ADJUSTED
C     TO BE A CONVEX COMBINATION OF THIS LIKELIHOOD AND THE CURRENT
C     LIKELIHOOD WITH THE RECOMBINATION PARAMETERS TAKING VALUES LESS
C     THAN .5.



c Ascertainment part should follow as the second half of the predigrees:

C     if(ped.gt.nped/2) loglik=-loglik


      END






      SUBROUTINE APEN(EXTRA,PAR,PEN,VAR,GENES,XLINK,ABSENT,XXRATE
     1,XYRATE,FIRST,LAST,MUTATE,NEXTRA,NGTYPE,NLOCI,NPAR,NVAR,PED
     2,PER,MALE,like)
C
C     PEN SUPPLIES THE PENETRANCE PROBABILITIES FOR EACH PERSON
C     IN A PEDIGREE.  THE CURRENT VERSION IS VALID ONLY FOR SIMPLE
C     EITHER/OR TRAITS.
C

      IMPLICIT NONE

      double precision absent
      double precision xxrate, xyrate
      integer          first,last
      integer          nextra
      integer          ngtype
      integer          nloci
      integer          npar
      integer          nvar
      integer          ped, per
      logical          male
      integer          mutate, like

      DOUBLE PRECISION EXTRA(NEXTRA),PAR(NPAR),PEN(NGTYPE),VAR(NVAR)
      INTEGER GENES(FIRST:LAST,2,NGTYPE)
      LOGICAL XLINK(NLOCI)

      double precision  rrbr(0:89)
      double precision  cumbrrisk
      double precision  popbr(0:89)
      double precision  ffbr(0:90)
      double precision  ffncbr(0:90)
      double precision  cumncbr
      double precision  p1, p2
      double precision  lambda(0:89,0:1)


      integer i
      integer idis
      integer agelfu, agebc, age, agebc2, ageoc, agedeath
      integer ageother, ageother2
      integer isex
      integer imut
      integer is
      integer iage





      common /popincid/popbr
      common /risk/ffbr



c---------------------------------------------------------------------
c The risks will correspond to the risk of a first BC
c Censoring occurs at the first of bc dx, oc dx, mastectomy or agelfu
c Note that this is the minimum of the age at last follow-up
c All censored at age 80.
c Ignore oophorectomy [CAN BE CHANGED LATER]

c Define the disease status, idis(=0 if unaffected, =1 if BC)

c Define the censoring variable in the data. If icens=1 then the individual
c is censored at birth.
c---------------------------------------------------------------------

       agebc  = var(1)
       agebc2 = var(2)
       ageoc  = var(3)
       ageother=var(4)
       ageother2=var(5)
       agelfu = var(6)
       agedeath= var(7)



c By default everyone censored at 0, unless age information is available:

       age = 0
       idis=0

	if (ageother.eq. 999) then
	 	if(agelfu.gt.0 .and. agedeath.gt.0) then
			ageother=min(agelfu,agedeath)
		elseif(agelfu.gt.0 .and. agedeath.eq.0) then
			ageother=agelfu
		elseif(agelfu.eq.0 .and. agedeath.gt.0) then
			ageother=agedeath
		else
			ageother=0
		endif
	endif


	if (agebc.eq. 999) then
	 	if(agelfu.gt.0 .and. agedeath.gt.0) then
			agebc=min(agelfu,agedeath)
		elseif(agelfu.gt.0 .and. agedeath.eq.0) then
			agebc=agelfu
		elseif(agelfu.eq.0 .and. agedeath.gt.0) then
			agebc=agedeath
		else
			agebc=0
		endif
	endif

	if (ageoc.eq. 999) then
	 	if(agelfu.gt.0 .and. agedeath.gt.0) then
			ageoc=min(agelfu,agedeath)
		elseif(agelfu.gt.0 .and. agedeath.eq.0) then
			ageoc=agelfu
		elseif(agelfu.eq.0 .and. agedeath.gt.0) then
			ageoc=agedeath
		else
			ageoc=0
		endif
	endif

C Censor at the first cancer.

	if(agebc.gt.0 .and. ageoc.gt.0 .and. ageother.gt.0) then
	 age=min(agebc, ageoc, ageother)
	elseif(agebc.gt.0 .and. ageoc.gt.0 .and. ageother.eq.0) then
	 age=min(agebc, ageoc)
	elseif(agebc.gt.0 .and. ageoc.eq.0 .and. ageother.gt.0) then
	 age=min(agebc, ageother)
	elseif(agebc.eq.0 .and. ageoc.gt.0 .and. ageother.gt.0) then
	 age=min(ageoc, ageother)
	elseif(agebc.gt.0 .and. ageoc.eq.0 .and. ageother.eq.0) then
	 age=agebc
	elseif(agebc.eq.0 .and. ageoc.gt.0 .and. ageother.eq.0) then
	 age=ageoc
	elseif(agebc.eq.0 .and. ageoc.eq.0 .and. ageother.gt.0) then
	 age=ageother
	elseif(agelfu.gt.0 .and. agedeath.eq.0) then
		age=agelfu
	elseif(agelfu.eq.0 .and. agedeath.gt.0) then
	        age=agedeath
	elseif(agelfu.gt.0 .and. agedeath.gt.0) then
		age=min(agelfu, agedeath)
	endif

	if (age.eq.agebc .and. agebc.gt.0) idis=1

c	write(*,*) agebc,ageoc,ageother,agelfu,agedeath,age,idis


c---------------------------------------------------------------------
c Censor at age 90
c---------------------------------------------------------------------

       if(age.ge.90) then
          idis=0
          age=90
       endif


c---------------------------------------------------------------------
c For carriers we estimate the incidence rates using the
c relative risks. In the fixed incidence version, no Rel Risk applies
c to non-carriers.
c---------------------------------------------------------------------

      
       do i=0,89
          if(i.lt.30) then
            rrbr(i) = exp(par(1))
          elseif(i.lt.60) then
            rrbr(i) = exp(par(2))
          else
            rrbr(i) = exp(par(3))
          endif
       end do
      


c---------------------------------------------------------------------
c This code constraints overall incidence to the population rates to obtain
c the female non carriers have incidence rates.
c---------------------------------------------------------------------

c non carrier frequency
       p1=(1-0.0004)**2

c carrier frequency:
       p2=1-p1

c initialise survival probabilities:

c non-carriers:

       ffncbr(0)=1.0d0

c carriers:

       ffbr(0)=1.0d0

c initialise sums:

       cumncbr=0.0d0
       cumbrrisk=0.0d0


       do iage=0,89

          lambda(iage,0)=(popbr(iage)/100000.0d0)*
     :             (p1*ffncbr(iage)+p2*ffbr(iage))/
     :              (p1*ffncbr(iage)+p2*rrbr(iage)*ffbr(iage))

          lambda(iage,1)=lambda(iage,0)*rrbr(iage)



c---------------------------------------------------------------------
c Cpmpute the survivor function for non carriers
c incidence rates
c---------------------------------------------------------------------

        cumncbr = cumncbr + lambda(iage,0)

        ffncbr(iage+1) = exp(-cumncbr)


c---------------------------------------------------------------------
c Compute the survivor function for carriers,
c---------------------------------------------------------------------


          cumbrrisk = cumbrrisk+lambda(iage,1)

          ffbr(iage+1) = exp(-cumbrrisk)


       end do



c---------------------------------------------------------------------
c Code the sex of each member of the pedigree
c---------------------------------------------------------------------

      if(male) then
       isex=1
      else
       isex=2
      end if




c---------------------------------------------------------------------
c Assume a  dominant model. So carriers denoted
c by is=2 and non carriers by  is=1
c---------------------------------------------------------------------

      DO 10 I=1,NGTYPE

       if(genes(1,1,I).eq.1.and.genes(1,2,I).eq.1) then
        is=1
       else
        is=2
       endif




c---------------------------------------------------------------------
c Assume that males do not develop either cancer and specify the penetrance
c for each male as 1*fact
c---------------------------------------------------------------------



c---------------------------------------------------------------------
c Non-Carriers develop the cancers according to population specific rates.
c---------------------------------------------------------------------

c=====================================================================
         if(is.eq.1) then
c=====================================================================


c---------------------------------------------------------------------
c Compute the penetrance for females
c---------------------------------------------------------------------


c---------------------------------------------------------------------
c For disease free females:
c---------------------------------------------------------------------

            if(idis.eq.0) then
              pen(i) = ffncbr(age)


c---------------------------------------------------------------------
c If she develops  breast cancer
c---------------------------------------------------------------------
            elseif(idis.eq.1) then
              pen(i) = ffncbr(age)*lambda(age,0)


            endif


c---------------------------------------------------------------------
c Carriers develop the cancers according to the fixed background
c incidence times the relative risks.
c---------------------------------------------------------------------

c=====================================================================
         else
c=====================================================================


c---------------------------------------------------------------------
c For disease free carriers:
c---------------------------------------------------------------------

            if(idis.eq.0) then
               pen(i) = ffbr(age)

c---------------------------------------------------------------------
c If carriers  develops  breast cancer
c---------------------------------------------------------------------
            elseif(idis.eq.1) then
               pen(i) = ffbr(age)*lambda(age,1)

            endif

c=====================================================================
         endif
c=====================================================================





c	write(*,*) ped, per, age, pen(i)

10    continue
      return
      END




      SUBROUTINE APRIOR(ALLFRQ,EXTRA,PAR,PRIOR,VAR,GENES,XLINK
     1,ABSENT,XXRATE,XYRATE,FIRST,LAST,MAXALL,MUTATE,NEXTRA,NGTYPE
     2,NLOCI,NPAR,NVAR,PED,PER,MALE)
C
C     PRIOR SUPPLIES THE PRIOR PROBABILITIES FOR EACH FOUNDER IN
C     A PEDIGREE.  THE CURRENT VERSION IS VALID IF HARDY-WEINBERG
C     AND LINKAGE EQUILIBRIUM HOLD AT ALL LOCI.  LOCI MAY BE AUTOSOMAL
C     OR X-LINKED.
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION ALLFRQ(NLOCI,MAXALL),EXTRA(NEXTRA),PAR(NPAR)
     1,PRIOR(NGTYPE),VAR(NVAR)
      INTEGER FIRST,GENES(FIRST:LAST,2,NGTYPE),PED,PER
      LOGICAL XLINK(NLOCI),MALE
C

      allfrq(1,2)=0.0004d0
      allfrq(1,1)=1.0d0-allfrq(1,2)



      DO 10 I=1,NGTYPE
      P=1.0D0
      DO 20 LOCUS=FIRST,LAST
      IG1=GENES(LOCUS,1,I)
      IF (XLINK(LOCUS).AND.MALE) THEN
      P=P*ALLFRQ(LOCUS,IG1)
      ELSE
      IG2=GENES(LOCUS,2,I)
      P=P*ALLFRQ(LOCUS,IG1)*ALLFRQ(LOCUS,IG2)
      END IF
 20   CONTINUE
 10   PRIOR(I)=P
      END





      SUBROUTINE ATRANS(EXTRA,PAR,TRANS,VARI,VARJ,GAMETE,GENES,XLINK
     1,ABSENT,XXRATE,XYRATE,FIRST,LAST,MUTATE,NEXTRA,NGTYPE,NLOCI
     2,NPAR,NVAR,PED,PERI,PERJ,MALEI,MALEJ)
C
C     TRANS SUPPLIES THE TRANSMISSION PROBABILITIES FOR EACH
C     PARENT-OFFSPRING PAIR IN A PEDIGREE.  THE SUFFIX I INDICATES
C     THE PARENT AND THE SUFFIX J THE CHILD.  SEVERAL MULTIPLE
C     LOCUS GENOTYPES OF THE PARENT ARE PASSED VIA THE ARRAY GENES.
C     THE ARRAY GAMETE REPRESENTS ONE OF THE TWO GAMETES MAKING UP
C     A MULTIPLE LOCUS GENOTYPE OF THE CHILD.


      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION EXTRA(NEXTRA),PAR(NPAR),TRANS(NGTYPE)
     1,VARI(NVAR),VARJ(NVAR)
      INTEGER FIRST,GAMETE(FIRST:LAST),GENES(FIRST:LAST,2,NGTYPE)
     1,PED,PERI,PERJ
      LOGICAL XLINK(NLOCI),MALEI,MALEJ

      DO 10 I=1,NGTYPE

         TRANS(I)=1.0

         T=0.0D0
         IF (GENES(1,1,I).EQ.GAMETE(1)) T=T+0.5D0
         IF (GENES(1,2,I).EQ.GAMETE(1)) T=T+0.5D0
         TRANS(I)=TRANS(I)*T



10    CONTINUE

      END
