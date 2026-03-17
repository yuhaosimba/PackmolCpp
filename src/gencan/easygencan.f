      subroutine easygcf(n,x,l,u,m,lambda,rho,epsgpsn,maxit,maxfc,
     +trtype,iprint,ncomp,f,g,gpsupn,iter,fcnt,gcnt,cgcnt,inform,wi,wd,
     +delmin)

      implicit none

C     SCALAR ARGUMENTS
      integer cgcnt,fcnt,gcnt,m,maxfc,maxit,n,ncomp,inform,iprint,iter
      double precision epsgpsn,f,gpsupn

C     ARRAY ARGUMENTS
      integer wi(n)
      double precision g(n),l(n),lambda(m),rho(m),u(n),wd(8*n),x(n)

C     This subroutine aims to simplify the use of GENCAN. For this 
C     purpose it gives values to most of the GENCAN arguments and 
C     leaves to the user those arguments which he/she may would like to 
C     set by him/herself.
C
C     The arguments of EASYGENCAN are the input and output arguments of 
C     GENCAN that are supposed to be useful for a common user. The input 
C     arguments are mostly related to basic problem information, like 
C     dimension and bounds, and the initial point. There are also input 
C     arguments related to simple stopping criteria (like norm of the 
C     projected gradient, and maximum number of iterations and 
C     functional evaluations). There are also two input arguments 
C     related to control the amount of information written into the 
C     screen. The output arguments are related to information of the 
C     solution and some few performance measurements. Basically, on 
C     return, EASYGENCAN gives to the user the solution, the objective 
C     functional value and its gradient at the solution, Euclidian and 
C     sup-norm of the projected gradient at the solution, the number of 
C     iterations, functional and gradient evaluations, and Conjugate 
C     Gradient iterations used to reach the solution, and, finally, a 
C     flag that indicates the stopping criterion that was satisfied.
C
C     All the other arguments of GENCAN are setted with its default 
C     values by EASYGENCAN. EASYGENCAN divides the arguments of GENCAN 
C     in two sets. Those that are related to the behaviour of GENCAN are 
C     declared as Fortran parameters (constants). The other arguments of 
C     GENCAN, most of them related to alternative stopping criteria, and 
C     that may depend of, for example, maxit, are declared as local 
C     variables of EASYGENCAN.
C
C     GENCAN arguments that are defined as Fortran parameters in this 
C     subroutine are GENCAN arguments that should not be modified by a 
C     common user. They are arguments that modify the behaviour of 
C     GENCAN and whos values were selected because they are classical 
C     values in some cases or because some numerical experiments seemed 
C     to indicate that they are the best choices.
C
C     GENCAN arguments that are declared as local variables in this 
C     subroutine are GENCAN arguments that may be modified if, with 
C     their suggested values, GENCAN does not give the desired result. 
C     Most of them are related to Conjugate Gradients or to disabled 
C     stopping criteria that may be useful in bad-scaled problems or 
C     problems with not trustable derivatives.
C
C     Finally, this subroutine declares as local variables some 
C     arguments of GENCAN which in fact are output arguments. Most of 
C     them are related to quantities that can be used for statistics 
C     related to the GENCAN performance, like number Spectral Projected 
C     Gradient iterations, Truncated Newton iterations, Conjugate 
C     Gradient iterations, etc. As we assume that this values are not 
C     useful for the common user, this subroutine throw all of them 
C     away.
C
C     We describe below the meaning of the arguments of the EASYGENCAN
C     subroutine. More detailed descriptions as well as the descriptions 
C     of all the other GENCAN arguments that are not arguments of 
C     EASYGENCAN are also described at the begining of the GENCAN 
C     subroutine.
C     
C     On entry:
C
C     n        integer 
C              number of variables
C
C     x        double precision x(n)
C              initial estimation of the solution
C
C     l        double precision l(n)
C              lower bounds on the variables
C
C     u        double precision u(n)
C              upper bounds on the variables
C
C     m        integer
C     lambda   double precision lambda(m)
C     rho      double precision rho(m)
C              These three parameters are not used nor modified by 
C              GENCAN and they are passed as arguments to the user-
C              defined subroutines evalal and evalnal to compute the 
C              objective function and its gradient, respectively. 
C              Clearly, in an Augmented Lagrangian context, if GENCAN is 
C              being used to solve the bound-constrainted subproblems, m 
C              would be the number of constraints, lambda the Lagrange 
C              multipliers approximation and rho the penalty parameters
C
C     epsgpsn  double precision
C              GENCAN stops declaring convergence if it finds a point 
C              whos projected gradient sup-norm is smaller than or equal 
C              to epsgpsn
C
C     maxit    integer
C              GENCAN stops declaring ''maximum number of iteration 
C              achieved'' if the number of iterations exceeds maxit
C
C     maxfc    integer
C              the same as before but with the number of functional 
C              evaluations
C
C     iprint   integer
C              indicates the degree of details of the output generated 
C              by GENCAN. Setting iprint to a value smaller than 2 will 
C              make GENCAN to generate no output at all. An iprint value 
C              greater than or equal to 2 will generate information of 
C              every GENCAN iteration. An iprint value greater than or 
C              equal to 3 will also show information of the Conjugate 
C              Gradient iterations (used to compute the Truncated Newton 
C              direction) and also information related to the line 
C              search procedures in the Spectral Projected Gradient 
C              direction and the Truncated Newton direction.
C
C     ncomp    integer
C              Sometimes, vectors like the current point x, the gradient 
C              of the objective function g, or the search directions 
C              (Spectral Projected Gradient direction or Truncated 
C              Newton direction), among other vector, are showed in the 
C              screen. In such cases, if the problem dimension is large, 
C              to show just a few elements of these vectors may be 
C              preferable. Argument ncomp can be used to indicate how 
C              many array elements must be displayed.
C
C     wi       integer wi(n)
C              integer working space
C
C     wd       double precision wd(8*n)
C              double precision working space
C
C     On return:
C
C     x        double precision x(n)
C              estimation of the solution
C
C     f        double precision
C              objective function value at the solution
C
C     g        double precision g(n)
C              gradient of the objective function at the solution
C
C     gpsupn   double precision
C              sup-norm of the continuous projected gradient
C
C     iter     integer
C              number of iterations used to reach the solution
C
C     fcnt     integer
C              number of functional evaluations
C
C     gcnt     integer
C              number of gradient evaluations
C
C     cgcnt    integer
C              number of Conjugate Gradient iterations
C
C     inform   integer
C              termination criteria. inform equal to 1 means that 
C              GENCAN converged with the sup-norm of the continuous 
C              projected gradient stopping criterion (inform equal to 0
C              means the same but with the Euclidian norm). Other 
C              positive values means that GENCAN stopped by a may be not  
C              successful stopping criteria. A negative value means that 
C              there was an error in the user-defined subroutines that 
C              computes the objective function (subroutine evalal), the 
C              gradient (subroutine evalnal), or the Hessian-vector
C              product (subroutine evalhd). See the GENCAN description 
C              for more details.

C     HERE STARTS THE DESCRIPTION OF SOME GENCAN ARGUMENTS THAT ARE 
C     BEING SETTED INSIDE EASYGENCAN. THE FIRST SET OF ARGUMENTS ARE 
C     THOSE ARGUMENTS THAT WE WILL CALL ''CONSTANTS'' AND THAT, AS THEIR 
C     VALUES ALTER THE BEHAVIOUR OF GENCAN, SHOULD NOT BE MODIFIED BY A 
C     COMMON USER.

C     CONSTANTS FOR GENERAL USES

C     Steps: h = max( steabs, sterel * abs( x ) ) should be a number 
C     such that h is small ( relatively to x ) and x + h is different 
C     from x. So, h is something that can be used a a step for a finite 
C     differences approximation of a partial derivative relative to x.

C     Epsilons: something smaller than max( epsabs, epsrel * abs( x ) )
C     should be considered as ``zero'' when compared with x. It is used, 
C     for example, to detect that a step taken during a line search is 
C     too small.

C     Infinitys: infrel is a big number that may appear in the 
C     calculations. infabs is a number that should never be reached in 
C     the calculations and is used the represent ``infinite''. Detailed 
C     explanations of how are they used are rather cumbersome.

      double precision steabs,sterel,epsabs,epsrel,infabs,infrel
      parameter ( steabs    = 1.0d-10 )
      parameter ( sterel    = 1.0d-07 )
      parameter ( epsabs    = 1.0d-20 )
      parameter ( epsrel    = 1.0d-10 )
      parameter ( infabs    = 1.0d+99 )
      parameter ( infrel    = 1.0d+20 )

C     CONSTANTS FOR CLASSICAL LINE-SEARCH CONDITIONS

C     beta is the constant for the ''beta condition''. We use this 
C     condition to test whether is promising to extrapolate or not.

C     gamma is the constant for the sufficient decrease ''Armijo 
C     condition''.

C     theta is the constant for the ''angle condition''.

C     sigma1 and sigma2 are the constants for the safeguarding quadratic 
C     interpolations. We use them in a rather unusual way. Instead of 
C     discarding a new step anew if it does not belong to the interval 
C     [ sigma1 * aprev, sigma2 * aprev ], we discard it if it does not 
C     belong to the interval [ sigma1, sigma2 * aprev ]. In such a case 
C     we take something similar to ''anew = aprev / 2''.

      double precision beta,gamma,theta,sigma1,sigma2
      parameter ( beta   =   0.5d0 )
      parameter ( gamma  = 1.0d-04 )
      parameter ( theta  = 1.0d-06 )
      parameter ( sigma1 =   0.1d0 )
      parameter ( sigma2 =   0.9d0 )

C     CONSTANTS FOR SPECIFIC PROCEDURES (NOT SO CLASSICAL)

C     In line searches, when interpolating, the step may become so 
C     small that we should declare a line search failure indicating that 
C     direction may not be a descent direction. This decision is never 
C     take before doing at least mininterp interpolations.

C     In line searches, the beta condition (see above) may recommend to
C     extrapolate. We never do more than maxextrap extrapolations.

C     In the line searches, when we need to interpolate and the result 
C     of the quadratic interpolation is rejected, the new step is 
C     computed as anew = aprev / nint. When the beta condition 
C     recommends to extrapolate, we compute anew = aprev * next.

C     When computing the Newton direction by Conjugate Gradients we 
C     never go further an artificial ''trust region''. This ''trust 
C     radius'' is never smaller than delmin.

C     In active set strategies, constants eta is used to decide whether 
C     the current face should be abandoned or not. In particular, the 
C     current face is abandoned when the norm of the internal to face 
C     component of the continuous projected gradient is smaller than 
C     ( 1 - eta ) times the norm of the continuous projected gradient. 
C     In this way, values of eta near 1 makes the method to work hard 
C     inside the faces and values of eta near 0 makes the method to 
C     abandon the faces very quickly.

C     We always use as a first step in a line search procedure along a
C     first order direction the spectral steplength. This steplength 
C     must belong to the interval [lspgmi,lspgma].

      integer maxextrap,mininterp
      parameter ( maxextrap = 100 )
      parameter ( mininterp =   4 )

      double precision nint,next,delmin,eta,lspgma,lspgmi
      parameter ( nint    =   2.0d0 )
      parameter ( next    =   2.0d0 )
c      parameter ( delmin  =   1.d4 )
      parameter ( eta     =   0.9d0 )
      parameter ( lspgma  = 1.0d+10 )
      parameter ( lspgmi  = 1.0d-10 )

C     DIMENSIONS FOR SOME WORKING SPACES

C     In non-monotone line searches, given p, the last p objective 
C     functional values must be stored. For this reason we declare a 
C     vector with pmax double precision elements. So p must be less than
C     or equal to pmax.

C     Sometimes, is the problem is bad scaled, to request a small 
C     gradient norm at the solution may be inadequate. For this reason, 
C     a test to verify if this norm is not decreasing during maxitngp 
C     (MAXimum of ITerations with No Gradient Progress) consecutive 
C     iterations then we stop the method with a warning. As it is not 
C     expected a monotone decreasing of the gradient norm, again, the 
C     norm of the last maxitngp iterations must be saved. For this 
C     purpose, we declare a vector of tmax elements. So maxitngp must 
C     be less than or equal to tmax.
     
      integer tmax
      parameter ( tmax = 1000 )

C     HERE STARTS THE DESCRIPTION OF THE OTHER ARGUMENTS OF GENCAN BEING 
C     SETTED BY EASYGENCAN. THESE ARGUMENTS MAY BE MODIFIED BY A COMMON 
C     USER IF, WITH THEIR SUGGESTED VALUES, GENCAN DOES NOT GIVE THE 
C     EXPECTED RESULT.

C     GENCAN INPUT ARGUMENTS THAT WILL BE SETTED BELOW

      logical nearlyq

      integer cgmaxit,cgscre,gtype,htvtype,maxitnfp,maxitngp,maxitnqmp,
     +        trtype

      double precision cgepsf,cgepsi,cggpnf,delta0,epsgpen,epsnfp,
     +        epsnqmp,fmin

C     GENCAN OUTPUT ARGUMENTS THAT WILL BE DISCARDED

      integer spgfcnt,spgiter,tnexbcnt,tnexgcnt,tnexbfe,tnexgfe,tnfcnt,
     +        tnintcnt,tnintfe,tniter,tnstpcnt

      double precision gpeucn2

C     GENCAN WORKING VECTORS (WHICH DIMENSION IS NOT RELATED TO THE 
C     PROBLEM DIMENSION)

      double precision lastgpns(tmax)

C     ARGUMENTS RELATED TO DERIVATIVES CALCULATIONS

C     gtype indicates in which way the gradient of the objective 
C     function will be computed. If the user have been implemented the 
C     user-supplied evalnal subroutine to compute the gradient of the 
C     objective function then gtype argument must be set to 0 (ZERO) and 
C     the user-supplied evalnal subroutine will be called by GENCAN any 
C     time the gradient would be required.
C
C     The prototype of the evalnal subroutine must be:
C
C         subroutine evalnal(n,x,m,lambda,rho,nal,flag)
C
C         SCALAR ARGUMENTS
C         integer n,m,flag
C
C         ARRAY ARGUMENTS
C         double precision x(n),lambda(m),rho(m),nal(n)
C
C         ''Here must be written the subroutine body that calculates the 
C         n-dimensional gradient vector of the objective function 
C         evaluated at x and saves it in nal. It also must set flag to 0 
C         (ZERO) if the gradient was successfully computed and to any 
C         other value if the gradient vector is not well defined at the 
C         required point x. If GENCAN is been used stand-alone to solve 
C         a unique bound-constrained problem then m, lambda and rho are 
C         dummy arguments. On the other hand, if GENCAN is been used in 
C         an Augmented Lagrangian framework then these arguments should 
C         be used for the number of constraints, the Lagrange 
C         multipliers approximation and the penalty parameters, 
C         respectively.''
C
C         end
C
C     If, on the other hand, the user is not able to provide evalnal 
C     subroutine, gtype argument must be set to 1 (ONE). In this case, 
C     every time GENCAN needs to compute the gradient of the objective 
C     function, an internal subroutine that approximates it by finite-
C     differences will be used (be aware that it maybe very time 
C     consuming). Moreover, note that the evalnal subroutine must still 
C     be present (with an empty body).

      gtype     =        0

C     htvtype indicates in which way the product of the Hessian of the
C     objective function times an arbitrary vector will be computed. If 
C     the user has not been implemented the user-supplied evalhd 
C     subroutine to do this task then htvtype argument must be set to 1 
C     (ONE). In this case an internal subroutine that approximates this 
C     product by incremental quotients will be used. Note that, even in 
C     this case, evalhd subroutine must be present (with an empty body). 
C     This is the default option and the empty-body subroutine follows:
C
C         subroutine evalhd(nind,ind,n,x,m,lambda,rho,d,hd,flag)
C
C         SCALAR ARGUMENTS
C         integer nind,n,m,flag
C
C         ARRAY ARGUMENTS
C         integer ind(nind)
C         double precision d(n),hd(n),lambda(m),rho(m),x(n)
C
C         flag = - 1
C
C         end
C
C     If, on the other hand, the user prefers to implement his/her own 
C     evalhd subroutine then htvtype argument must be set to 0 (ZERO). 
C     In this case, the product of the Hessian times vector d (input 
C     argument of evalhd subroutine) must be saved in vector hd (output 
C     argument of evalhd subroutine). The other arguments description as
C     well as some hints on how to implement your own evalhd subroutine 
C     can be found in the GENCAN arguments description.

C     When ALGENCAN uses GENCAN to solve the subproblems in the classical
C     Augmented Lagrangian framework, ALGENCAN uses its own evalhd
C     subroutine to overcome the lack of continuity of the second 
C     derivatives. So, when GENCAN is being used toghether with ALGENCAN,
C     htvtype must be equal to 0 (ZERO). On the other hand, if GENCAN is
C     being used stand-alone, just set htvtype equal to 1 (ONE) and add 
C     the empty-body subroutine described above.

      htvtype   =        1

C     ARGUMENTS RELATED TO STOPPING CRITERIA

C     Besides the stopping criterion related to the sup-norm of the 
C     continuous projected gradient, there is another stopping criterion 
C     related to its Euclidian norm. So, GENCAN stops the process if it 
C     finds a point at which the Euclidian norm of the continuous 
C     projected gradient is smaller than epsgpen.

      epsgpen   =    0.0d0

C     For an explanation of maxitngp see above the explanation of tmax 
C     in ''DIMENSIONS FOR SOME WORKING SPACES''. Just note that the 
C     value of maxitngp must be less than or equal to tmax.

      maxitngp  =     tmax

C     maxitnfp means MAXimum of allowed number of iterations with No 
C     Progress in the objective functional value. ''Progress'' from one 
C     iteration to the next one refers to ( fnew - fprev ). Since the 
C     begining of the algorithm we save the ''best progress'' and 
C     consider that there was no progress in an iteration if the 
C     progress of this iterations was smaller than epsnfp times the best 
C     progress. Finally, the algorithm stops if there was no progress 
C     during maxitnfp consecutive iterations.

      maxitnfp  =    maxit
      epsnfp    =    0.0d0

C     There is a stopping criterion that stops the method if a point 
C     with a functional value smaller than fmin is found. The idea 
C     behind this stopping criterion is to stop the method if the 
C     objective function is not bounded from below.

      fmin      = 1.0d-05

C     ARGUMENTS RELATED TO CONJUGATE GRADIENTS

C     When computing the Truncated Newton direction by Conjugate 
C     Gradients there is something similar to a ''trust-region radius''. 
C     This trust radius is updated from iteration to iteration depending 
C     on the agreement of the objective function and its quadratic 
C     model. But an initial value for the trust radius is required. If 
C     the user has a good guess for this initial value then it should be 
C     passed to GENCAN using the delta0 arguments. On the other hand, if 
C     delta0 is set to -1, a default value depending on the norm of the 
C     current point will be used.

      delta0    =  - 1.0d0
      delmin  =  1.d-2
c      delta0 = delmin

C     The ''trust-region'' can be like a ball (using Euclidian norm) or 
C     like a box (using sup-norm). This choice can be made using trtype 
C     (TRust region TYPE) argument. trtype equal to 0 means Euclidian 
C     norm and trtype equal to 1 means sup-norm.

      trtype    =        1

C     When the method is far from the solution, it may be not useful to 
C     do a very large effort in computing the Truncated Newton direction 
C     precisely. To avoid it, a fixed maximum number of iterations for 
C     Conjugate Gradients can be given to GENCAN. If the user would like 
C     to choose this maximum number of iterations for Conjugate 
C     Gradient then it should use the cgmaxit arguments. On the other 
C     hand he/she prefers to leave this task to GENCAN then he/she 
C     should set cgmaxit to -1.
 
      cgmaxit   =   -1

C     If the task of deciding the accuracy for computing the Truncated 
C     Newton direction is leaved to GENCAN then a default strategy based 
C     on increasing accuracies will be used. The proximity to the 
C     solution is estimated observing the norm of the projected gradient 
C     at the current point and locating it between that norm at the 
C     initial point and the expected value of that norm at the solution. 
C     Then the accuracy for the Truncated Newton direction of the 
C     current iteration will be computed taking a precision located in 
C     the same relative position with respect to two given values for 
C     the accuracies for the first and the last Truncated Newton 
C     direction calculations. These two accuracies (cgepsi and cgepsf, 
C     respectively) must be given by the user. Moreover, the expected 
C     value of the projected gradient norm at the solution (cggpnf) must
C     also be given by the user who must indicate setting argument 
C     cgscre to 1 or 2 if that norm is the Euclidian or the sup-norm.
      
      cggpnf    =  max( 1.0d-04, max( epsgpen, epsgpsn ) ) 
      cgscre    =        2
      cgepsi    =  1.0d-01
      cgepsf    =  1.0d-05

C     The next two arguments are used for an alternative stopping 
C     criterion for Conjugate Gradients. Conjugate Gradients method is 
C     stopped if the quadratic model makes no progress during maxitnqmp 
C     (MAXimum of ITerations with No Quadratic Model Progress) 
C     consecutive iterations. In this context, ''no progress'' means 
C     that the progress is smaller than epsnqmp (EPSilon to measure the 
C     No Quadratic Model Progress) times the best progress obtained 
C     during the previous iterations.

      epsnqmp   =  1.0d-04
      maxitnqmp =        5

C     Depending on how much the objective function seems to be a 
C     quadratic, function, Conjugate Gradients may take different 
C     decision. So, if the objective function is a quadratic function or 
C     is very similar to a quadratic function then the nearlyq argument 
C     should be set to TRUE, else, it should be set to FALSE. However, 
C     the option with nearlyq equal TRUE never showed good results. 
C     Regarding this unexpected no good performance, rather recently it 
C     was found a bug that affected the behaviour of GENCAN just in this 
C     case (See the April 1st, 2003 modifications report at the end of 
C     this file). So, new experiments setting nearlyq equal TRUE should 
C     be made. 

      nearlyq   =   .false.

C     FINALLY, CALL GENCAN

      call gencan(n,x,l,u,m,lambda,rho,epsgpen,epsgpsn,maxitnfp,epsnfp,
     +maxitngp,fmin,maxit,maxfc,delta0,cgmaxit,cgscre,cggpnf,cgepsi,
     +cgepsf,epsnqmp,maxitnqmp,nearlyq,nint,next,mininterp,maxextrap,
     +gtype,htvtype,trtype,iprint,ncomp,f,g,gpeucn2,gpsupn,iter,fcnt,
     +gcnt,cgcnt,spgiter,spgfcnt,tniter,tnfcnt,tnstpcnt,tnintcnt,
     +tnexgcnt,tnexbcnt,tnintfe,tnexgfe,tnexbfe,inform,wd(1),wd(n+1),
     +wd(2*n+1),wi,lastgpns,wd(3*n+1),eta,delmin,lspgma,lspgmi,theta,
     +gamma,beta,sigma1,sigma2,sterel,steabs,epsrel,epsabs,infrel,
     +infabs)

      end

C     ******************************************************************
C     ******************************************************************

C     Last update of GENCAN or any of its dependencies: 
C
C     February 18th, 2005.
C
C     See report of modifications at the end of this file.
