      Program SimpleExample
C Simple example of using of UpperLim.f
C Example of the compilation instruction with GNU Fortran f77=g77:

C  f77 -o SimpleExample SimpleExample.f UpperLim.f y_vs_CLf.f CMaxinf.f ConfLev.f Cinf.f CERN_Stuff.f
C Or, for gfortran,
C  gfortran -frecord-marker=4 -o SimpleExample SimpleExample.f UpperLim.f y_vs_CLf.f CMaxinf.f ConfLev.f Cinf.f CERN_Stuff.f
C
C The data will consist of just two events described later in this example.
      Implicit None
      Integer If,N,Iflag
      Real sigma,mu0,sigma0,muB,CL
      Real UpperLim,FC(0:10) ! FC holds the data.  It can store up to 9 events.
C FC(0) is reserved for UpperLim use, as is FC(10) with 9 events.
C
C Assume the distribution of events is expected to be uniform from 5-100 keV,
C and suppose there are N=2 events, one at 10 keV and one at 50 keV.  Then
C given an event, it has probability (10-5)/(100-5) = 0.05263 of being
C below 10 keV and probability (50-5)/(100-5) = 0.47368 of being below
C 50 keV.  Put these values into FC in order of lowest to highest energy.
      N=2 ! Two events
      FC(1)=0.05263 ! Probability of event being below 10 keV
      FC(2)=0.47368 ! Probability of event being below 50 keV
C Here we could put FC(0)=0.0 and FC(3)=1.0, but UpperLim will do it for us.
C 
      CL=0.9 ! 90% Confidence level
      If=1 ! fmin = 0.
      muB=0. ! Don't subtract background
C Assume cross section sigma0=1.E-42 cm^2 would produce expected number of
C events mu0= 1.85.  Then the factor converting the upper limit expected number
C of events into the upper limit cross section is sigma0/mu0.
      sigma0=1.E-42
      mu0=1.85
      sigma=(sigma0/mu0)*UpperLim(CL,If,N,FC,muB,FC,Iflag)
C sigma is now the upper limit cross section.  Iflag is a status; 0 is good.
      Write(6,10) sigma,Iflag
 10   Format('Upper limit cross section =',E9.4,'cm^2 with status flag',
     1 I5)
      Stop
      End
