C	Program StrucFunc


	subroutine StrucFunc(dfile_id, dW2, dQ2, dT, dYoniIndex, dIPOL, dIA1, dSFChoice, dAsymChoice, dAnswer, F1, F2, G1, G2 )
C	Program to generate f1(x,Q2), f2(x,Q2) etc. using proper low-Q2 behavior
C	
C	4/2/97 SEK
C       
C	Use new NMC parametrization for F2 and
C	and SLAC 1998 fit for R
C	
C	Rewritten to call subroutines in "models.f" for all
C	calculations
c       
C       Final H version 10/26/1999 SEK
C       
C	Added necessary stuff for deuteron 12/2/1999 SEK
C       
C 	Changed to allow calling of new Hall C parametrization for F1, Fl 23-Feb-2004 SEK
C       
C	Updated to refer only to the newest models 15-Jun-2009 SEK
C       
C	*********************************************************
C 	Modified by C. Ayerbe to be used as a subroutine wrapped
C	trough a C++ code.
C 	The subroutine needs to receive W2, Q2, and the options
C	for the different model to use.
C	it should return pol/unpol SF and more depending of the arguments
C	*********************************************************

	
	implicit none

	include 'instruct.inc'
	
C       sk	REAL*8 VSTRT(37),STP(37),MY_LO(37),MY_HIGH(37),XVAL(37)
C       sk	INTEGER NPRM(37)
C       sk	CHARACTER*10 PNAM(37)

	real*8 a, b

	real*8 dQ2, dW2
	integer  dYoniIndex, dIPOL, dIA1, dSFChoice, dAsymChoice, dfile_id
	character*1 dT, dAnswer


	
	real*8 Q2, Nu, X, R, F2, F1, qsq, W2, DR,E
	real*8 TH, G1, G2, A1, A2
	real*8 fn,fnerr,xval1(41),xvall(41),w2inc,temp(4)
        real*8 mp,mp2,pi,alpha

	character*1 Answer, T, nT /'N'/, pT /'P'/ !Target Type
	logical dummy, suppress	!, newhallc
	integer model/12/, i, imax, npts, sf, YoniIndex
	integer pSFChoice, nSFChoice, pAsymChoice, nAsymChoice
	integer pIPOL, nIPOL, pIA1, nIA1
		
	common /YoniDeut/ pSFChoice, nSFChoice, pAsymChoice, nAsymChoice, 
     *		pIPOL, nIPOL, pIA1, nIA1

C	print *, "options: Target: ",  dT," Yoni: ", dYoniIndex, " IPOL: ",dIPOL, " IA1: ",
C     *	dIA1," SFCHOICE: ", dSFChoice, " Asym: ", dAsymChoice, "Answer :" ,dAnswer

	SFChoice = dSFChoice
	if ((T.eq.'n').or.(T.eq.'N')) then
	   nSFChoice = dSFChoice
	elseif ((T.eq.'p').or.(T.eq.'P')) then
	   pSFChoice = dSFChoice
	endif

	Q2 = dQ2
	W2 = dW2
	YoniIndex = dYoniIndex
	IPOL =  dIPOL
	IA1 = dIA1
	AsymChoice = dAsymChoice
	Answer =  dAnswer	
	T = dT
	
C	print*, W2, Q2, T, IPOL	
	   
	NMC = .TRUE.
CC	suppress = .FALSE.
	BODEK1 = .FALSE.
	NORES = .FALSE.
	ERROR = .FALSE.
CC	IPOL = 1
c       sk	IPOL = 1: standard version of A2 in the DIS region
c       sk	IPOL = 2: extra term g2tw3 in A2 in DIS
c       sk	IPOL = 4: g2 = 0 in DIS
c       sk	IPOL = 5: A2 = 0 in DIS
c       sk	IPOL = 7: A2 = Soffer bound in DIS
CC	IPOLRES = 1
CC	IA1 = 4			! Use to determine specific A1/A2 DIS model used
c       sk     IA1 = 4: sek_me (2006/7 combined Doaa/sek fit to world data including EG1)
c       sk     IA1 = 5: sek_me; add correlated fit error to A1
c       sk     IA1 = 6: sek_me; subtract correlated fit error to A1
c       sk     IA1 = 7: sek_2021 (New version by Pushpa)
c       sk     IA1 = 8: sek_2021 + dA1
c       sk     IA1 = 9: sek_2021 - dA1

CC	AsymChoice = 11		! use to determine specific A1/A2 models used in the RR
c       sk    11: New Standard Resonance Model 2008-9 (c) Guler/Kuhn
c       sk    12: Preliminary version v1 of A1 model - 2008-9 (c) Guler/Kuhn 
c       sk    13: Old Standard for A1 2008-9 (c) Guler/Kuhn
c       sk    14: Alternative A2 resonance model: A2_MAID
c       sk    15: Old version of A2 resonance model 

CC	SFChoice = 23
c       sk    10: 2007 version of F1n/F1p/F1A by Peter Bosted/Eric Christie (c) 2007, HERMES
c       sk    11: 2009 version of F1n/F1p/F1A by Peter Bosted/Eric Christie (c) 2009, HERMES
c       sk    12: Same version as 11, but with errors added to F2 (and proportionally F1)
c       sk    13: Same version as 11, but with errors subtracted from R (F2 unchanged)
c       sk    14: Same version as 11, but without tabulated R values in RR substituted
c       sk    15-16: Analog to 12-13, but without tabulated R values in RR substituted
c       sk    17: Newest kludge for Rd, F2d from Eric Christy December 2011. Should be only for D
c       sk    18-19: Analog to 15-16 for new kludge
c       sk	20: New version by Christie/Kalantarians 2014
c       sk	21-22 corresponding error estimates
c       sk 23-25: New version F1F221
	
C	open (unit=22, file='f1.in',status='old')
C	open (unit=11, file='f1.out',status='new')



c$$$ 1009	write(6,*) ' Enter Target type (P,N,3 or D - please capitalize):'
c$$$	read(5,111) T
c$$$ 111	format(a1)
c$$$	if ((T .ne. 'P').and.(T .ne. 'D').and.(T. ne. 'N')
c$$$     *		.and.(T .ne. '3')) goto 1009
c$$$C       ********************************************************
c$$$
c$$$	
c$$$C       THIS LINE OPTION--> ASK XIAOCHAO
c$$$1523	write(6,*) ' Enter smearing preference: ',
c$$$     *		'0 = none, 1 = quasielastic only, 2 = all, 3 = inelastic only.'
c$$$      	read(5,112) YoniIndex
c$$$
c$$$ 112	format(I1)
c$$$	if (YoniIndex.gt.3 .or. YoniIndex.lt.0) goto 1523
c$$$
c$$$
c$$$
c$$$ 	
c$$$	if ((T.eq.pT).or.(T.eq.nT).or.(YoniIndex.lt.1)) then
c$$$ 1245	   write(6,*) ' What value of IPOL do you want use? 1-7'
c$$$	   read(5,*) IPOL
c$$$	   if (IPOL.lt.1 .or. IPOL.gt.7) goto 1245
c$$$	   
c$$$	   
c$$$	   
c$$$ 1246	   write(6,*) ' What value of IA1 do you want use? 4-9'
c$$$	   read(5,*) IA1
c$$$	   if (IA1.lt.4 .or. IA1.gt.9) goto 1246
c$$$	   
c$$$	   
c$$$	   
c$$$ 1247	   continue
c$$$	   if (T.eq. 'D' .and. (YoniIndex.lt.1)) then
c$$$	      write(6,*) ' What value for SFChoice? 10-19'
c$$$	      read(5,*) SFChoice
c$$$	      if (SFChoice.lt.10 .or. SFChoice .gt. 19) goto 1247
c$$$	   else
c$$$	      write(6,*) ' What value for SFChoice? 11-25'
c$$$	      read(5,*) SFChoice
c$$$	      if (SFChoice.lt.10 .or. SFChoice .gt. 25) goto 1247
c$$$	   endif
c$$$	   
c$$$	   
c$$$ 1248	   write(6,*) ' What value for AsymChoice? 11-13'
c$$$	   read(5,*) AsymChoice
c$$$	   if (AsymChoice.lt.11 .or. AsymChoice .gt. 13) goto 1248
c$$$	   
c$$$	else
c$$$	   
c$$$	   write(6,*) ' Enter model parameters for proton'
c$$$ 2245	   write(6,*) ' What value of IPOL do you want use? 1-7'
c$$$	   read(5,*) pIPOL
c$$$	   if (pIPOL.lt.1 .or. pIPOL.gt.7) goto 2245
c$$$	   
c$$$ 2246	   write(6,*) ' What value of IA1 do you want use? 4-9'
c$$$	   read(5,*) pIA1
c$$$	   if (pIA1.lt.4 .or. pIA1.gt.9) goto 2246
c$$$	   
c$$$ 2247	   write(6,*) ' What value for SFChoice? 10-25'
c$$$	   read(5,*) pSFChoice
c$$$	   if (pSFChoice.lt.10 .or. pSFChoice .gt. 25) goto 2247
c$$$	   
c$$$ 2248	   write(6,*) ' What value for AsymChoice? 11-13'
c$$$	   read(5,*) pAsymChoice
c$$$	   if (pAsymChoice.lt.11 .or. pAsymChoice .gt. 13) goto 2248
c$$$	   
c$$$	   write(6,*) ' Enter model parameters for neutron'
c$$$ 3245	   write(6,*) ' What value of IPOL do you want use? 1-7'
c$$$	   read(5,*) nIPOL
c$$$	   if (nIPOL.lt.1 .or. nIPOL.gt.7) goto 3245
c$$$	   
c$$$ 3246	   write(6,*) ' What value of IA1 do you want use? 4-6'
c$$$	   read(5,*) nIA1
c$$$	   if (nIA1.lt.4 .or. nIA1.gt.6) goto 3246
c$$$	   
c$$$ 3247	   write(6,*) ' What value for SFChoice? 10-25'
c$$$	   read(5,*) nSFChoice
c$$$	   if (nSFChoice.lt.10 .or. nSFChoice .gt. 25) goto 3247
c$$$	   
c$$$ 3248	   write(6,*) ' What value for AsymChoice? 11-13'
c$$$	   read(5,*) nAsymChoice
c$$$	   if (nAsymChoice.lt.11 .or. nAsymChoice .gt. 13) goto 3248
c$$$	   
c$$$	endif
c$$$	
c$$$	write(6,*) ' Include Resonances? [Y]/n'
c$$$	read(5,111) Answer
c$$$	if ((Answer.eq.'N').or.(Answer.eq.'n')) then
c$$$	   NORES = .TRUE.
c$$$C	this makes h2mod use only the nonres. = background terms for sigma
c$$$	   suppress = .TRUE.
c$$$c       this skips h2mod entirely, using a DIS extrapolation instead.
c$$$	endif
c$$$
c$$$
c$$$
c$$$
c$$$	
c	do i = 1, 200000
c	   read(22,*,end=99) W2, Q2
	   Nu = (W2-0.939**2+Q2)/2.0/0.939
	   qsq = Q2 + Nu*Nu
	   X = Q2/2.0/0.939/Nu
	   if (YoniIndex .eq. 0)then
C	      print *, 'going to F1F2new'
	      call F1F2new (x, Q2, T, F1, F2, R, suppress)
C	      print*,"FROM F1F2new: ",F1, Q2
C	R is here only used as placeholder for sign/sigp
	      TH = 25.0 D00	! Arbitrary, not really needed
	      call G1G2new(x, Q2, TH, T, G1, G2, A1, A2)
	      if (F1 .gt. 0. .and. Nu .gt. 0. .and. Q2 .gt. 0.) then
		 R = F2/F1/2.0/x*(1+Q2/Nu/Nu)-1.0
	      else
		 R = 0.0
	      endif
	   else if (T .ne. 'P' .and. T .ne. 'N') then
	      print *, 'going to DSF ', T
	      call DSFs(X, Q2, YoniIndex, F1, F2, R, G1, G2, A1, A2, suppress)
	   else
C	      print *, 'going to Yoni ', T
	      call Yoni(X, Q2, T, YoniIndex, F1, F2, R, G1, G2, A1, A2, suppress)
	   endif
C	   write(11,1234) Q2,W2,X,F1,F2,R,A1,A2,G1,G2

 1234	   format(' ',10f12.6)
C	   print*,  Q2,W2,F1,F2,G1,G2
c	enddo
c 99	continue
c	close (22)
c	close (11)
c	Stop
c$$$

C       THIS IS JUST TO CONTROL FLUX OF FILES
C	print*, "file_id: ", dfile_id
	return
	end

