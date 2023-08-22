CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC                                                                      CCC
CCC   F1F221 version 0.997  -  April 10, 2021                            CCC
CCC   Collection of subroutines to calculate inclusive cross sections    CCC
CCC   for range of nuclei.  For A > 2 the parameterization is based on   CCC
CCC   by M. E. Christy, T. Gautam, and A Bodek to 12C, 27Al, 56Fe and    CCC
CCC   64Cu.  However, the fit scales relatively well with 'A' and        CCC
CCC   should be good for all nuclei with 10 < A < 80.                    CCC
CCC   Also included is the proton cross section fit and a preliminary    CCC
CCC   deuteron/neutron fit by M. E. Christy, N. Kalantarians, J. Either  CCC
CCC   and W. Melnitchouk (to be published) based on both inclusive       CCC
CCC   deuteron and tagged deuteron data on n/d from BONuS.               CCC
CCC   Range of validity is W^2 < 32,  Q^2 < 32.0                         CCC
CCC   New data included in deuteron fit, including photoproduction at    CCC
CCC   Q^2 = 0.                                                           CCC
CCC                                                                      CCC
CCC                                                                      CCC
CCC   A > 2 nuclei fit are 4He, 12C, 27Al, and 56Fe.  More nuclei        CCC
CCC   will be fit in the next version.                                   CCC
CCC                                                                      CCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CSK
CSK REDUCED SET OF SUBROUTINES NEEDED FOR "newSFs.f"
CSK ONLY valid for p, n
CSK Modified by S. Kuhn June 2021
CSK
      
      SUBROUTINE RESCSP21(w2,q2,sigtp,siglp)
CCCC  Returns proton transverse and longitudinal cross sections   CCCC
CCCC  February 4, 2021                                            CCCC
      
      IMPLICIT none

      real*8 w2,q2,sigtp,siglp
      real*8 xvalp(100),xval1(50),xvalL(50)
      Integer i
      real*8 sigtdis,sigLdis,w2max,w2min

      data xvalp / 
     & 0.12287E+01,0.15194E+01,0.15044E+01,0.17010E+01,0.16743E+01,
     & 0.14468E+01,0.12480E+00,0.23000E+00,0.91261E-01,0.87852E-01,
     & 0.77402E-01,0.37857E+00,0.76997E+01,0.45404E+01,0.42604E+01,
     & 0.18299E+01,0.68745E+01,0.19713E-07,0.36429E+04,0.44334E+01,
     & 0.59193E+00,0.18984E+02,0.62651E-01,0.24619E+01,0.76417E+00,
     & 0.67506E+01,0.18840E+00,0.18956E+01,0.23118E+01,0.49609E+01,
     & 0.35804E+01,0.26219E+01,0.11076E-02,0.10000E+05,0.80725E-01,
     & 0.23272E+01,0.25357E+00,0.73390E+00,0.20111E+01,0.92868E-01,
     & 0.12817E-01,0.40981E+00,0.19435E+01,0.44089E+01,0.14325E-01,
     & 0.17597E+00,0.19923E+01,0.55000E+00,0.42998E+01,0.42338E+00,
     & 0.99122E+00,0.99085E+00,0.99798E+00,0.10028E+01,0.98145E+00,
     & 0.10163E+01,0.10226E+01,0.10165E+01,0.00000E+00,0.00000E+00,
     & 0.00000E+00,0.00000E+00,0.11353E-05,0.10605E+02,0.19705E+01,
     & 0.00000E+00,0.19429E-09,0.75543E+02,0.60997E+01,0.00000E+00,
     & 0.73839E+01,0.59820E-05,0.19247E+01,0.00000E+00,0.23239E+01,
     & 0.16384E+01,0.14220E+01,0.00000E+00,0.10015E-03,0.57993E+00,
     & 0.64963E+00,0.00000E+00,0.10093E-03,0.41844E+00,0.42011E+00,
     & 0.00000E+00,0.39197E+00,0.85540E+01,0.10000E+00,0.10000E+01,
     & -.84009E+01,0.50000E+02,0.53511E+02,0.29984E+01,0.11295E-02,
     & 0.48118E+00,0.00000E+00,0.00000E+00,0.00000E+00,0.18796E+02 / 

      do i=1,50
        xval1(i) = xvalp(i)
        xvalL(i) = xvalp(50+i) 
        if(i.LE.12) xvalL(i) = xval1(i)
        if(i.EQ.47.OR.i.EQ.48) xvalL(i) = xval1(i)
      enddo


      w2max = 34.0  !!! max for resonance fit
      w2min = 36.0  !!! min for dis fit

      call resmodp21(1,w2,q2,xval1,sigTp)
      call resmodp21(2,w2,q2,xvalL,sigLp)

c      if(w2.GT.w2min) then 
c        call disp(w2,q2,sigTdis,sigLdis)
c        if(w2.LE.w2max) then
c         sigTp = (w2max-w2)*sigTp+(w2-w2min)*sigTdis
c         sigLp = (w2max-w2)*sigLp+(w2-w2min)*sigLdis
c         sigTp = sigTp/2.0
c         sigLp = sigLp/2.0
c         else
c          sigTp = sigTdis
c          sigLp = sigLdis
c        endif
c      endif

      return
      end


      SUBROUTINE RESCSN21(w2,q2,sigtn,sigln)
CCCC  Returns neutron transverse and longitudinal cross sections   CCCC
CCCC  January 15, 2021 version                                     CCCC      
      IMPLICIT none

      real*8 w2,q2,sigtn,sigln
      real*8 xvaln(100),xval1(50),xvalL(50)
      integer i
      data xvaln / 
     & 0.12287E+01,0.15195E+01,0.15042E+01,0.17062E+01,0.16787E+01,
     & 0.14474E+01,0.12502E+00,0.23000E+00,0.90468E-01,0.86480E-01,
     & 0.75000E-01,0.38012E+00,0.72000E+01,0.34406E+01,0.14715E+01,
     & 0.21795E+01,0.28269E+01,0.14249E+01,0.18085E+03,0.11604E+01,
     & 0.29603E+01,0.92847E+01,0.21512E+02,0.22213E+01,0.51569E+03,
     & 0.82835E+03,0.85402E+02,0.64661E+05,0.26191E+00,0.22097E+02,
     & 0.28088E+00,0.23009E+01,0.27367E+01,0.69913E+04,0.85719E+04,
     & 0.32685E+02,0.33579E+00,0.16534E+01,0.26841E+02,0.11457E+00,
     & 0.81967E-01,0.52000E+00,0.22614E+01,0.49169E+01,-.78524E-01,
     & 0.10827E-01,0.19782E+01,0.54078E+00,0.26027E+01,0.19898E-01,
     & 0.10186E+01,0.97493E+00,0.99258E+00,0.98113E+00,0.10425E+01,
     & 0.10154E+01,0.98069E+00,0.98176E+00,0.10402E+01,0.10018E+01,
     & 0.10146E+01,0.10114E+01,0.13775E+02,0.67322E+02,0.58419E+01,
     & 0.25048E+01,0.64847E+00,0.21008E+02,0.30625E+01,0.29919E+00,
     & 0.48625E-03,0.72055E+01,0.22848E+01,0.77741E+02,0.14788E+01,
     & 0.12655E+02,0.29691E+01,0.98035E+00,0.11913E+02,0.23614E+02,
     & 0.76568E+01,0.25731E+01,0.29579E-03,0.24162E+01,0.83270E+00,
     & 0.99999E+00,0.86995E-01,0.29264E+01,0.20862E+01,0.10000E+01,
     & 0.36439E+02,0.51543E+01,0.91411E+01,0.77592E+01,0.68505E-01,
     & 0.84495E+00,0.19410E+01,0.44998E+00,0.10011E+01,0.37921E+01  /
      
          
      do i=1,50
        xval1(i) = xvaln(i)
        xvalL(i) = xvaln(50+i) 
        if(i.LE.12) xvalL(i) = xval1(i)  !!! use same masses for L as T !!!
        if(i.EQ.47.OR.i.EQ.48) xvalL(i) = xval1(i)
      enddo

      call resmodn21(1,w2,q2,xval1,sigtn)
      call resmodn21(2,w2,q2,xvalL,sigLn)

      return
      end
     


CCC-----------------

      SUBROUTINE RESMODP21(sf,w2,q2,xval,sig)
CCC  Returns proton transverse (sf=1) and longitudinal (sf=2 Cross sections CCC      
CCC  Version from February 21, 2021  -  Author:  M.E. Christy               CCC
CCC  This routine returns proton photo-absorbtion cross sections            CCC
CCC  for either transverse or longitudinal photons in units of ub/Sr/Gev.   CCC
CCC                                                                         CCC
CCC  Fit form is empirical.  Interpret physics from it at your own risk.    CCC
CCC  replaced 2-pi threshold with eta                                       CCC
      
      IMPLICIT NONE
      REAL*8 W,w2,q2,mp,mp2,mpi2,xb,sig,xval(50),mass(7),width(7)
      REAL*8 height(7),rescoef(6,4)
      REAL*8 nr_coef(3,4),sigr(7),wdif(2),sig_nr,pi2,alpha
      REAL*8 mpi,meta,intwidth(7),k,kcm,kr(7),kcmr(7),ppicm,ppi2cm
      REAL*8 petacm,ppicmr(7),ppi2cmr(7),petacmr(7),epicmr(7),epi2cmr(7)
      REAL*8 eetacmr(7),epicm,epi2cm,eetacm,br(7,3),ang(7)
      REAL*8 pgam(7),pwid(7,3),x0(7),dip,mon,q20,A0
      REAL*8 sig_res,xpr(2),t1,t2
      INTEGER i,j,num,sf
      
      mp = 0.9382727
      mpi = 0.134977
      mpi2 = mpi*mpi
      meta = 0.547862
      mp2 = mp*mp
      alpha = 1./137.036
      W = sqrt(w2)
      wdif(1) = w - (mp + mpi)
      wdif(2) = w - (mp + meta)
      
      q20 = 0.05
      q20= xval(50)

       
CCCC   single pion branching ratios  CCCC

      br(1,1) = 1.00       !!!  P33(1232)       
      br(2,1) = 0.45      !!!  S11(1535)   
      br(3,1) = 0.60      !!!  D13(1520)
      br(4,1) = 0.65      !!!  F15(1680)
      br(5,1) = 0.60       !!!  S11(1650)
      br(6,1) = 0.65     !!!  P11(1440) roper 
      br(7,1) = 0.60      !!!  F37(1950)

CCCC  eta branching ratios   CCCC

      br(1,3) = 0.0       !!!  P33(1232)
      br(2,3) = 0.40      !!!  S11(1535) 
      br(3,3) = 0.08      !!!  D13(1520)
      br(4,3) = 0.0       !!!  F15(1680)
      br(5,3) = 0.20      !!!  S11(1650)
      br(6,3) = 0.0       !!!  P11(1440) roper   
      br(7,3) = 0.0       !!!  F37(1950)

CCCC  2-pion branching ratios  CCCC

      do i=1,7
        br(i,2) = 1.-br(i,1)-br(i,3)
      enddo

CCCC   Meson angular momentum   CCCC

      ang(1) = 1.       !!!  P33(1232)
      ang(2) = 0.       !!!  S11(1535)
      ang(3) = 2.       !!!  D13(1520)
      ang(4) = 3.       !!!  F15(1680)
      ang(5) = 0.       !!!  S15(1650)
      ang(6) = 1.       !!!  P11(1440) roper   
      ang(7) = 3.       !!!  F37(1950)

      do i=1,7     !!!  resonance damping parameter  !!!
         x0(i) = 0.178   !!! 
      enddo

c      x0(1) = 0.125
      x0(1) = 0.142
      
      do i=1,7
        br(i,2) = 1.-br(i,1)-br(i,3)
      enddo
               
      dip = 1./(1.+q2/1.05)**2.             !!!  Dipole parameterization  !!!
      mon = 1./(1.+q2/1.05)**1.
           
      xb = q2/(q2+w2-mp2)
      xpr(1) = 1.+(w2-(mp+mpi)**2)/(q2+q20)
      xpr(1) = 1./xpr(1)
      xpr(2) = 1.+(w2-(mp+meta)**2)/(q2+q20)
      xpr(2) = 1./xpr(2)
      if(w.LE.(mp+mpi)) xpr(1) = 1.0      
      if(w.LE.(mp+meta)) xpr(2) = 1.0
     

CCC    Calculate kinematics needed for threshold Relativistic B-W  CCC

      k = (w2 - mp2)/2./mp
      kcm = (w2-mp2)/2./w

      epicm = (W2 + mpi**2 -mp2 )/2./w
      ppicm = SQRT(MAX(0.0,(epicm**2 - mpi**2)))
      epi2cm = (W2 + (2.*mpi)**2 -mp2 )/2./w
      ppi2cm = SQRT(MAX(0.0,(epi2cm**2 - (2.*mpi)**2)))
      eetacm = (W2 + meta*meta -mp2 )/2./w
      petacm =  SQRT(MAX(0.0,(eetacm**2 - meta**2)))

      num = 0

      do i=1,6              !!!  Read in resonance masses     !!!
        num = num + 1
        mass(i) = xval(i)
      enddo
      do i=1,6              !!!  Read in resonance widths     !!!
        num = num + 1
        intwidth(i) = xval(num)
        width(i) = intwidth(i)
      enddo

      mass(7) = xval(47)
      intwidth(7) = xval(48)
      width(7) = intwidth(7) 

      do i=1,7
        kr(i) = (mass(i)**2-mp2)/2./mp
        kcmr(i) = (mass(i)**2.-mp2)/2./mass(i)
        epicmr(i) = (mass(i)**2 + mpi**2 -mp2 )/2./mass(i)
        ppicmr(i) = SQRT(MAX(0.0,(epicmr(i)**2 - mpi**2)))
        epi2cmr(i) = (mass(i)**2 + (2.*mpi)**2 -mp2 )/2./mass(i)
        ppi2cmr(i) = SQRT(MAX(0.0,(epi2cmr(i)**2 - (2.*mpi)**2)))
        eetacmr(i) = (mass(i)**2 + meta*meta -mp2 )/2./mass(i)
        petacmr(i) =  SQRT(MAX(0.0,(eetacmr(i)**2 - meta**2)))

CCC   Calculate partial widths   CCC

        pwid(i,1) = intwidth(i)*(ppicm/ppicmr(i))**(2.*ang(i)+1.)
     &           *((ppicmr(i)**2+x0(i)**2)/(ppicm**2+x0(i)**2))**ang(i)
c         !!!  1-pion decay mode


        pwid(i,2) = intwidth(i)*(ppi2cm/ppi2cmr(i))**(2.*ang(i)+4.)
     &         *((ppi2cmr(i)**2+x0(i)**2)/(ppi2cm**2+x0(i)**2))
     &         **(ang(i)+2)   !!!  2-pion decay mode

        pwid(i,2) = W/mass(i)*pwid(i,2)


        pwid(i,3) = 0.          !!!  eta decay mode


        if(i.EQ.2.OR.i.EQ.5) then
          pwid(i,3) =  intwidth(i)*(petacm/petacmr(i))**(2.*ang(i)+1.)
     &          *((petacmr(i)**2+x0(i)**2)/(petacm**2+x0(i)**2))**ang(i)
c         !!!  eta decay only for S11's 
        endif 


        pgam(i) = (kcm/kcmr(i))**2*
     &                   (kcmr(i)**2+x0(i)**2)/(kcm**2+x0(i)**2)

        pgam(i) = intwidth(i)*pgam(i)

        width(i) = br(i,1)*pwid(i,1)+br(i,2)*pwid(i,2)+br(i,3)*pwid(i,3)

      enddo

CCC    End resonance kinematics and Widths calculations   CCC


CCC    Begin resonance Q^2 dependence calculations   CCC
           

      do i=1,6
        do j=1,4
          num = num + 1
          rescoef(i,j)=xval(num)
        enddo

        if(sf.EQ.1) then

          height(i) = rescoef(i,1)*
     &          (1.+rescoef(i,2)*q2/(1.+rescoef(i,3)*q2))*
     &          mon**rescoef(i,4)


        else

           height(i) = (rescoef(i,1)+rescoef(i,2)*q2)
     &           *exp(-1.*rescoef(i,3)*q2)

          
        endif
 
        height(i) = height(i)*height(i)

      enddo


      if(sf.EQ.2) then      !!!  4th resonance region  !!!

        height(7) = (xval(44)+xval(45)*q2)*exp(-1.0*xval(46)*q2)
        
      else
        height(7) = xval(49)*mon    
      endif
      
      height(7) = height(7)*height(7)


CCC    End resonance Q^2 dependence calculations   CCC

     
      do i=1,3               !!!  Non-Res coefficients  !!!
        do j=1,4
          num = num + 1
          nr_coef(i,j)=xval(num)
        enddo
      enddo


CCC   Calculate Breit-Wigners for all resonances   CCC

      sig_res = 0.0

      do i=1,7

       sigr(i) = width(i)*pgam(i)/((W2 - mass(i)**2.)**2. 
     &              + (mass(i)*width(i))**2.)
        sigr(i) = height(i)*kr(i)/k*kcmr(i)/kcm*sigr(i)/intwidth(i)
        sig_res = sig_res + sigr(i)   
      enddo

      sig_res = sig_res*w
      if(sf.EQ.2) sig_res = sig_res*q2


CCC    Finish resonances / start non-res background calculation   CCC

 
      sig_nr = 0.

      if(sf.EQ.1.and.xpr(1).LT.1.0) then
        A0 = (1.0+xval(44)*q2)/(1.+q2/xval(42))**xval(43)  !!! overall amplitude
        t1 = xval(38)*log(2.0+q2/xval(39))                 !!! exponent of (1-xpr)
        t2 = xval(40)*log(2.0+q2/xval(41))+xval(45)        !!! exponent of xpr


        if(xpr(1).LE.1.0) then                             !!! 1-pi threshold
          sig_nr = xval(37)*389.4*A0*(1.-xpr(1))**t1*xpr(1)**t2
        endif
          
        if(xpr(2).LE.1.0) then                             !!! eta threshold 
          sig_nr = sig_nr+xval(46)*389.4*A0*(1.-xpr(2))**t1*xpr(2)**t2
         endif
       
      elseif(sf.EQ.2.and.xpr(1).LT.1.0) then
         
         t1 = xval(38)*log(1.001+q2/q20)                   !!! exponent of (1-xpr)
         t2 = xval(41)/(1.001+q2/xval(42))**xval(43)       !!! exponent of xpr

        
        if(xpr(1).LE.1.0) then
          sig_nr = sig_nr + xval(37)*389.4*
     &       xb*(1.-xpr(1))**t1*xpr(1)**t2
        endif
        sig_nr =sig_nr/log(1.001+q2/xval(39))
        
      endif
    
      sig = sig_res + sig_nr

      if((w-mp).LT.wdif(1)) sig = 0.0    


 1000  format(8f12.5)
 1001  format(7f12.3)

      RETURN 
      END 



      SUBROUTINE RESMODN21(sf,w2,q2,xval,sig)
CCC  Returns proton transverse (sf=1) and longitudinal (sf=2 Cross sections CCC      
CCC  Version from February 21, 2021  -  Author:  M.E. Christy               CCC
CCC  This routine returns proton photo-absorbtion cross sections            CCC
CCC  for either transverse or longitudinal photons in units of ub/Sr/Gev.   CCC
CCC                                                                         CCC
CCC  Fit form is empirical.  Interpret physics from it at your own risk.    CCC
CCC  replaced 2-pi threshold with eta                                       CCC
      
      IMPLICIT NONE
      REAL*8 W,w2,q2,mp,mp2,mpi2,xb,sig,xval(50),mass(7),width(7)
      REAL*8 height(7),rescoef(6,4)
      REAL*8 nr_coef(3,4),sigr(7),wdif(2),sig_nr,pi2,alpha
      REAL*8 mpi,meta,intwidth(7),k,kcm,kr(7),kcmr(7),ppicm,ppi2cm
      REAL*8 petacm,ppicmr(7),ppi2cmr(7),petacmr(7),epicmr(7),epi2cmr(7)
      REAL*8 eetacmr(7),epicm,epi2cm,eetacm,br(7,3),ang(7)
      REAL*8 pgam(7),pwid(7,3),x0(7),dip,mon,q20,A0
      REAL*8 sig_res,xpr(2),t1,t2
      INTEGER i,j,num,sf
      
      mp = 0.939565
      mpi = 0.134977
      mpi2 = mpi*mpi
      meta = 0.547862
      mp2 = mp*mp
      alpha = 1./137.036
      W = sqrt(w2)
      wdif(1) = w - (mp + mpi)
      wdif(2) = w - (mp + meta)
      
      q20 = 0.05
      q20= xval(50)

      
CCCC   single pion branching ratios  CCCC

      br(1,1) = 1.00       !!!  P33(1232)       
      br(2,1) = 0.45      !!!  S11(1535)   
      br(3,1) = 0.60      !!!  D13(1520)
      br(4,1) = 0.65      !!!  F15(1680)
      br(5,1) = 0.60       !!!  S11(1650)
      br(6,1) = 0.65     !!!  P11(1440) roper 
      br(7,1) = 0.60      !!!  F37(1950)

CCCC  eta branching ratios   CCCC

      br(1,3) = 0.0       !!!  P33(1232)
      br(2,3) = 0.40      !!!  S11(1535) 
      br(3,3) = 0.08      !!!  D13(1520)
      br(4,3) = 0.0       !!!  F15(1680)
      br(5,3) = 0.20      !!!  S11(1650)
      br(6,3) = 0.0       !!!  P11(1440) roper   
      br(7,3) = 0.0       !!!  F37(1950)

CCCC  2-pion branching ratios  CCCC

      do i=1,7
        br(i,2) = 1.-br(i,1)-br(i,3)
      enddo

CCCC   Meson angular momentum   CCCC

      ang(1) = 1.       !!!  P33(1232)
      ang(2) = 0.       !!!  S11(1535)
      ang(3) = 2.       !!!  D13(1520)
      ang(4) = 3.       !!!  F15(1680)
      ang(5) = 0.       !!!  S15(1650)
      ang(6) = 1.       !!!  P11(1440) roper   
      ang(7) = 3.       !!!  F37(1950)

      do i=1,7     !!!  resonance damping parameter  !!!
         x0(i) = 0.178   !!! 
      enddo

c      x0(1) = 0.125
      x0(1) = 0.142
      
      do i=1,7
        br(i,2) = 1.-br(i,1)-br(i,3)
      enddo
               
      dip = 1./(1.+q2/1.05)**2.             !!!  Dipole parameterization  !!!
      mon = 1./(1.+q2/1.1)**1.

      xb = q2/(q2+w2-mp2)
      xpr(1) = 1.+(w2-(mp+mpi)**2)/(q2+q20)
      xpr(1) = 1./xpr(1)
      xpr(2) = 1.+(w2-(mp+meta)**2)/(q2+q20)
      xpr(2) = 1./xpr(2)
      if(w.LE.(mp+mpi)) xpr(1) = 1.0      
      if(w.LE.(mp+meta)) xpr(2) = 1.0

CCC    Calculate kinematics needed for threshold Relativistic B-W  CCC

      k = (w2 - mp2)/2./mp
      kcm = (w2-mp2)/2./w

      epicm = (W2 + mpi**2 -mp2 )/2./w
      ppicm = SQRT(MAX(0.0,(epicm**2 - mpi**2)))
      epi2cm = (W2 + (2.*mpi)**2 -mp2 )/2./w
      ppi2cm = SQRT(MAX(0.0,(epi2cm**2 - (2.*mpi)**2)))
      eetacm = (W2 + meta*meta -mp2 )/2./w
      petacm =  SQRT(MAX(0.0,(eetacm**2 - meta**2)))

      num = 0

      do i=1,6              !!!  Read in resonance masses     !!!
        num = num + 1
        mass(i) = xval(i)
      enddo
      do i=1,6              !!!  Read in resonance widths     !!!
        num = num + 1
        intwidth(i) = xval(num)
        width(i) = intwidth(i)
      enddo

      mass(7) = xval(47)
      intwidth(7) = xval(48)
      width(7) = intwidth(7) 

      do i=1,7
        kr(i) = (mass(i)**2-mp2)/2./mp
        kcmr(i) = (mass(i)**2.-mp2)/2./mass(i)
        epicmr(i) = (mass(i)**2 + mpi**2 -mp2 )/2./mass(i)
        ppicmr(i) = SQRT(MAX(0.0,(epicmr(i)**2 - mpi**2)))
        epi2cmr(i) = (mass(i)**2 + (2.*mpi)**2 -mp2 )/2./mass(i)
        ppi2cmr(i) = SQRT(MAX(0.0,(epi2cmr(i)**2 - (2.*mpi)**2)))
        eetacmr(i) = (mass(i)**2 + meta*meta -mp2 )/2./mass(i)
        petacmr(i) =  SQRT(MAX(0.0,(eetacmr(i)**2 - meta**2)))

CCC   Calculate partial widths   CCC

        pwid(i,1) = intwidth(i)*(ppicm/ppicmr(i))**(2.*ang(i)+1.)
     &           *((ppicmr(i)**2+x0(i)**2)/(ppicm**2+x0(i)**2))**ang(i)
c         !!!  1-pion decay mode


        pwid(i,2) = intwidth(i)*(ppi2cm/ppi2cmr(i))**(2.*ang(i)+4.)
     &         *((ppi2cmr(i)**2+x0(i)**2)/(ppi2cm**2+x0(i)**2))
     &         **(ang(i)+2)   !!!  2-pion decay mode

        pwid(i,2) = W/mass(i)*pwid(i,2)


        pwid(i,3) = 0.          !!!  eta decay mode


        if(i.EQ.2.OR.i.EQ.5) then
          pwid(i,3) =  intwidth(i)*(petacm/petacmr(i))**(2.*ang(i)+1.)
     &          *((petacmr(i)**2+x0(i)**2)/(petacm**2+x0(i)**2))**ang(i)
c         !!!  eta decay only for S11's 
        endif 


        pgam(i) = (kcm/kcmr(i))**2*
     &                   (kcmr(i)**2+x0(i)**2)/(kcm**2+x0(i)**2)

        pgam(i) = intwidth(i)*pgam(i)

        width(i) = br(i,1)*pwid(i,1)+br(i,2)*pwid(i,2)+br(i,3)*pwid(i,3)

      enddo

CCC    End resonance kinematics and Widths calculations   CCC


CCC    Begin resonance Q^2 dependence calculations   CCC
           

      do i=1,6
        do j=1,4
          num = num + 1
          rescoef(i,j)=xval(num)
        enddo

        if(sf.EQ.1) then

          height(i) = rescoef(i,1)*
     &          (1.+rescoef(i,2)*q2/(1.+rescoef(i,3)*q2))*
     &          mon**rescoef(i,4)


        else

           height(i) = (rescoef(i,1)+rescoef(i,2)*q2)
     &           *exp(-1.*rescoef(i,3)*q2)


        endif
                  
        height(i) = height(i)*height(i)

      enddo


      if(sf.EQ.2) then      !!!  4th resonance region  !!!

        height(7) = (xval(44)+xval(45)*q2)*exp(-1.0*xval(46)*q2)

      else
        height(7) = xval(49)*mon
      endif

      height(7) = height(7)*height(7)


CCC    End resonance Q^2 dependence calculations   CCC

     
      do i=1,3               !!!  Non-Res coefficients  !!!
        do j=1,4
          num = num + 1
          nr_coef(i,j)=xval(num)
        enddo
      enddo


CCC   Calculate Breit-Wigners for all resonances   CCC

      sig_res = 0.0

      do i=1,7

       sigr(i) = width(i)*pgam(i)/((W2 - mass(i)**2.)**2. 
     &              + (mass(i)*width(i))**2.)
        sigr(i) = height(i)*kr(i)/k*kcmr(i)/kcm*sigr(i)/intwidth(i)
        sig_res = sig_res + sigr(i)   
      enddo

      sig_res = sig_res*w
      if(sf.EQ.2) sig_res = sig_res*q2


CCC    Finish resonances / start non-res background calculation   CCC

 
      sig_nr = 0.

      if(sf.EQ.1.and.xpr(1).LT.1.0) then
        A0 = (1.0+xval(44)*q2)/(1.+q2/xval(42))**xval(43)  !!! overall amplitude
        
        t1 = xval(38)*log(2.51+q2/xval(39))                !!! exponent of (1-xpr)
        t2 = xval(40)*log(2.51+q2/xval(41))+xval(45)       !!! exponent of xpr
        
        
        if(xpr(1).LE.1.0) then                             !!! 1-pi threshold
          sig_nr = xval(37)*389.4*A0*(1.-xpr(1))**t1*xpr(1)**t2
        endif
          
        if(xpr(2).LE.1.0) then                             !!! eta threshold 
          sig_nr = sig_nr+xval(46)*389.4*A0*(1.-xpr(2))**t1*xpr(2)**t2
         endif
       
      elseif(sf.EQ.2.and.xpr(1).LT.1.0) then
         
         t1 = xval(38)*log(1.001+q2/q20)                   !!! exponent of (1-xpr)
         t2 = xval(41)/(1.001+q2/xval(42))**xval(43)       !!! exponent of xpr

        
        if(xpr(1).LE.1.0) then
          sig_nr = sig_nr + xval(37)*389.4*
     &       xb*(1.-xpr(1))**t1*xpr(1)**t2
        endif
        sig_nr =sig_nr/log(1.001+q2/xval(39))


      endif
    
      sig = sig_res + sig_nr

      if((w-mp).LT.wdif(1)) sig = 0.0    


 1000  format(8f12.5)
 1001  format(7f12.3)

      RETURN 
      END 


      
      SUBROUTINE SF21(w2,q2,F1p,FLp,F2p,F1n,FLn,F2n)
CCCC   Converts reduced cross sections to structure functions for protons and neutrons  CCCCC 

      IMPLICIT none

      real*8 w2,q2,x,sigtp,siglp,sigtn,sigln,f1p,f2p,fLp
      real*8 f1n,f2n,fLn,pi,pi2,alpha,mp,mp2

      mp = 0.938272
      mp2 = mp*mp
      pi = 3.14159
      pi2 = pi*pi
      alpha = 1/137.03599 
      x = q2/(q2+w2-mp2)

   
      call rescsp21(w2,q2,sigTp,sigLp)
      call rescsn21(w2,q2,sigTn,sigLn)

      f1p = sigTp/0.3894e3/pi2/alpha/8.0*abs(w2-mp2)
      f1n = sigTn/0.3894e3/pi2/alpha/8.0*abs(w2-mp2)
      fLp = sigLp*2.0*x/0.3894e3/pi2/alpha/8.0*abs(w2-mp2)
      fLn = sigLn*2.0*x/0.3894e3/pi2/alpha/8.0*abs(w2-mp2)
      f2p = (2.*x*f1p+fLp)/(1.+4.*mp2*x*x/q2)
      f2n = (2.*x*f1n+fLn)/(1.+4.*mp2*x*x/q2)

      return
      end

      
CCC   -----------------
