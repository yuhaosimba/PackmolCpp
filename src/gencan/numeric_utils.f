      double precision function norm2s(n,x)

      implicit none

C     SCALAR ARGUMENTS
      integer n

C     ARRAY ARGUMENTS
      double precision x(n)

C     This subroutine computes the squared Euclidian norm of an 
C     n-dimensional vector.
C
C     On entry:
C
C     n        integer
C              dimension
C
C     x        double precision x(n)
C              vector
C
C     On return:
C
C     The function return the squared Euclidian norm of the 
C     n-dimensional vector x.

      external hsldnrm2
      double precision hsldnrm2

      norm2s = hsldnrm2(n,x,1) ** 2

      return 

      end

C     ******************************************************************
C     ******************************************************************

      DOUBLE PRECISION FUNCTION HSLDNRM2(N,DX,INCX)
      DOUBLE PRECISION ZERO,ONE
      PARAMETER (ZERO=0.0D0,ONE=1.0D0)
      DOUBLE PRECISION CUTLO,CUTHI
      PARAMETER (CUTLO=8.232D-11,CUTHI=1.304D19)
      INTEGER INCX,N
      DOUBLE PRECISION DX(*)
      DOUBLE PRECISION HITEST,SUM,XMAX
      INTEGER I,J,NN
      INTRINSIC DABS,DSQRT,FLOAT
      IF (N.GT.0) GO TO 10
      HSLDNRM2 = ZERO
      GO TO 300
   10 CONTINUE
      SUM = ZERO
      NN = N*INCX
      I = 1
   20 CONTINUE
      IF (DABS(DX(I)).GT.CUTLO) GO TO 85
      XMAX = ZERO
      IF (DX(I).EQ.ZERO) GO TO 200
      IF (DABS(DX(I)).GT.CUTLO) GO TO 85
      GO TO 105
  100 I = J
      SUM = (SUM/DX(I))/DX(I)
  105 XMAX = DABS(DX(I))
      GO TO 115
      IF (DABS(DX(I)).GT.CUTLO) GO TO 75
      IF (DABS(DX(I)).LE.XMAX) GO TO 115
      SUM = ONE + SUM* (XMAX/DX(I))**2
      XMAX = DABS(DX(I))
      GO TO 200
  115 SUM = SUM + (DX(I)/XMAX)**2
      GO TO 200
   75 SUM = (SUM*XMAX)*XMAX
   85 HITEST = CUTHI/DFLOAT(N)
      DO 95 J = I,NN,INCX
        IF (DABS(DX(J)).GE.HITEST) GO TO 100
        SUM = SUM + DX(J)**2
   95 END DO
      HSLDNRM2 = DSQRT(SUM)
      GO TO 300
  200 CONTINUE
      I = I + INCX
      IF (I.LE.NN) GO TO 20
      HSLDNRM2 = XMAX*DSQRT(SUM)
  300 CONTINUE
      RETURN
      END

C     ******************************************************************
C     ******************************************************************
C
C     Report of modifications.
C
C     February 18th, 2005.
C
C     1) An unsed format statement, previously used to automaticaly
C     generates some tables, was deleted.
C
C     2) An unmateched parenthesis was corrected in the format
C     statement used to stop GENCAN due to a small step in a line search.
C
C     February 16th, 2005.
C
C     1) The evalhd subroutine used by default in GENCAN is now the one
C     implemented in calchddiff, which approximates the Hessian-vector
C     product by incremental quotients. The implementation used to 
C     overcome the non twice continuously differentiability of the 
C     classical (PHR) Augmented Lagrangian function is now part of 
C     ALGENCAN (and not GENCAN). So, to use GENCAN inside ALGENCAN, 
C     htvtype argument must be set equal to 0 (ZERO).
C
C     2) The commented version of the empty function evalhd that must
C     be added when GENCAN is beinf used stand-alone was wrong. The
C     arguments declarations had been copied from evalnal. It was 
C     corrected.
C
C     November 10th, 2004.
C
C     1) After several test, all references to nonmontone line search
C     schemes were deleted.
C
C     September 28th, 2004.
C
C     1) Subroutines were checked an some absent arguments explanations
C     were added
C
C     2) Some calling sequences were modified to group related arguments
C
C     3) Arguments and local variables declarations were reordered in
C     alphabetical order.
C
C     3) Shrink and expand subroutines were modified to deal with just
C     one vector at a time. In this way, they are now being called from
C     calc* subroutines.
C
C     September 27th, 2004.
C
C     1) All comments were arranged to fit into the 72-columns format
C
C     2) Unused variable goth, which was prepared to indicate whether 
C     the Hessian matrix have been evaluated at the current point, was 
C     deleted from CG subroutine.
C
C     3) A spell check was used to correct the comments
C
C     September 21th, 2004.
C
C     1) In the stopping criterion where the progress in the objective 
C     function is verified, ''itnfp .ge. maxitnfp'' was changed for 
C     ''itnfp .gt. maxitnfp'', to make the choice maxitnfp equal to 1 
C     sounds reasonable.
C
C     2) Moreover, the previous chance came from the addition in the 
C     comments of GENCAN of the ''constraints'' information which makes 
C     clear to the user the values each argument may assume.
C
C     3) In the calculations of the first ''trust-radius'' for Conjugate 
C     Gradients, ''if( udelta0 .lt. 0.d0 ) then'' was changed by ''if 
C     ( udelta0 .le. 0.0d0 ) then'' to also make the default GENCAN 
C     choice of this initial trust-radius in the case of the user have 
C     been setted udelta = 0 by mistake.
C
C     4) The same for ucgmaxit.
C
C     5) In the line search subroutines spgls and tnls, ''if ( interp 
C     .gt. mininterp .and. samep ) then'' was changes by ''.ge.''.
C
C     6) Some comments of GENCAN arguments were re-written.
C
C     September 16th, 2004.
C
C     1) With the reconfiguration of the calc* subroutines (see (1) 
C     below) there were a number of redundant parameters in calchd and 
C     evalhd subroutines. These parameters were eliminated.
C
C     September 13th, 2004.
C
C     1) Subroutines named calc* that work in the reduced space always
C     call the corresponding eval* subroutine. As it was, calcg (that
C     computes the gradient in the reduced space) called evalg or 
C     evalgdiff depending on gtype parameter. The same was for calchd. 
C     Now, calcg calls evalg, calchd calls evalhd, and calchddiff (new) 
C     approximates the Hessian times vector product by incremental 
C     quotients calling calcg or calcgdiff depending on gtype parameter.
C     An improvement of this modification is that calcg does not call 
C     evalg or evalgdiff (both work in the full space) any more but it 
C     approximates the gradient vector in the reduced space (by central 
C     finite differences) calling 2 * nind times evalf subroutine.
C
C     2) Some comments were added inside evalg and evalhd user supplied
C     subroutines alerting about the relation of these subroutines and
C     the parameters gtype and htvtype, respectively.
C
C     3) Description of tnls subroutine was slightly modified.
C
C     4) The description of htvtype parameter in gencan was again 
C     slightly modified.
C
C     5) With the introduction of the parameter lambda (that in the
C     context of Augmented Lagrangians is used to store the 
C     approximation of the Lagrange multipliers) the name of the 
C     variable used for spectral steplength was changed from lambda to 
C     lamspg. In addition, lammax was changed to lspgma and lammin to 
C     lspgmi.
C
C     6) Modifications introduced in June 15th, 2004 and May 5th, 2004
C     were, in fact, made in this version on September 13th, 2004.
C
C     June 15th, 2004.
C
C     1) The fmin stopping criterion and the maximum number of
C     functional evaluation stopping criterion were erroneously being 
C     tested before the main loop. It was just redundant and, for this 
C     reason, deleted.
C
C     May 5th, 2004.
C
C     1) Incorporated into an Augmented Lagrangian framework.
C
C     a) evalf and evalg were renamed as evalal and evalnal, 
C        respectively.
C
C     b) m,lambda,rho were added as parameters of the subroutines evalal 
C        and evalnal, and, as a consequence, as parameters of almost all 
C        the other subroutines.
C
C     2) The comment of htvtype parameter of gencan was in portuguese
C     and it was translated into english.
C
C     3) A nonmonotone version of gencan is starting to be studied.
C     Parameters p and lastfv(0:p-1) were added to gencan, spgls, and
C     tnls to allow a nonmonotone line search. Array lastfv is now 
C     been updated for saving the last p functional values and the 
C     nonmonotone line searches are been done in a SPG or a 
C     Truncated Newton direction. p = 1 means monotone line search 
C     and is recommended until this study finish.
C
C     April 13th, 2004.
C
C     1) The modifications introduced in the occasion of the IRLOC 
C     development and re-development (October 21th, 2003 and February 
C     19th, 2003, respectively) were in fact made in this version on 
C     April 13th, 2004. The motivation to do this was to unify two 
C     parallel and different version of GENCAN (created, obviously, by 
C     mistake).
C
C     2) The complete reference of the GENCAN paper was finally added.
C
C     May 14th, 2003.
c
C     1) The way amax2 and amax2n were being computing may caused a 
C     segmentation fault. Its initialization was changed from infty and
C     -infty to 1.0d+99 and -1.0d+99, respectively. Using infty, when
C     combined with a big trust region radius, the final value of amax2
C     or amax2n may cause the impression that a bound is being attained, 
C     when it is not. "Redundant" ifs inside the amax2 and anax2n 
C     calculation were deleted. It should considered the possibility of 
C     using two constants, namely, bignum = 1.0d+20 and infty = 1.0d+99, 
C     instead of just infty. 
C
C     Modification introduced in October 21, 2003 in occasion of the
C     IRLOC re-development:
C
C     1) The stooping criteria related to functional value smaller than
C     fmin and exhaustion of maximum allowed number of functional 
C     evaluations have been done after the line search. And the 
C     questions were done as "if line search flag is equal to 4" or "if 
C     line search flag is equal to 8". But it was wrong in the case, for 
C     example, inside the line search, a functional value such that f <= 
C     fmin and the Armijo criterion was satisfied. In such case, the 
C     line search flag was being setted to 0 and not to 4. And gencan 
C     did not stop by the fmin criterion. Now, both stooping criteria 
C     are tested at the begining of the main gencan loop and just the 
C     stooping criteria by small line search step is tested after the 
C     line search.
C
C     Modification introduced in February 19, 2003 in occasion of the
C     IRLOC development:
C
C     1) The description of epsnfp parameter of GENCAN was modified. It
C     was written that to inhibit the related stopping criterion (lack
C     of function progress) it was necessary just set epsnfp = 0 when
C     it is also necessary to set maxitnfp = maxit. it was added in the
C     explanation.
C
C     2) In the explanation at the beginning of GENCAN it was written 
C     that cgscre parameter should be double precision. This comment was 
C     wrong. The correct type for cgscre parameter is integer.
C
C     Modifications introduced near April 1st 2003 in occasion of the 
C     PHR and inequality-constraints Augmented Lagrangian methods 
C     development:
C
C     1) The use of iprint was redefined and iprint2 was deleted.
C
C     2) The way to detect no progress in the log of the projected 
C     gradient norm was changed. As it was, ''no progress'' means no
C     reduction in the projected gradient norm over M iterations.
C     But this criterion implicitly assumed that the projected
C     gradient norm must decrease monotonously. Is it is clearly not
C     true, the criterion was changed by a non-monotone decrease
C     criterion. Now, progress means that the projected gradient
C     norm is, at each iteration, smaller than the maximum over the
C     last M iterations. And "no progress" means the it does not 
C     occurs during  not smaller than the 
C
C     3 ) The computation of qamaxn inside cg subroutine was in the 
C     wrong place (it was being used before computed) and it may was 
C     the reason for which the option nearlyq = .true. never worked 
C     properly. With this correction this option should be tested again.
C
C     On September 29th, 2004, we did a new test using the 41 bound
C     constrained problems with quadratic objective function from the
C     CUTE collection. The behaviour of GENCAN setting nearly equal
C     to true or false was indistinguishable. The test did not
C     include the different choices for the maximum number of CG
C     iterations being restricted to evaluate the different
C     alternatives for the case of finding a direction d such that
C     d^t H d <= 0. As a conclusion of this experiment we continue
C     recommending as a default choice to set nearlyq equal to false.
C
C     Modifications introduced from March 1st to March 21th of 2002
C     in occasion of the ISPG development:
C
C     1) Comments of some new parameters introduced in the previous
C     modification
C
C     2) As it was, in the first iteration of GENCAN (when kappa takes
C     value equal 1) and for one-dimensional faces, cgmaxit(the maximum 
C     number of Conjugate Gradient iterations to compute the internal to
C     the face truncated-Newton direction) was being 0. As it is 
C     obviously wrong, we add a max between what was being computed and 
C     one to allow at least one CG iteration.
C
C     3) Parameter inform in subroutines evalf, evalg and evalhd 
C     supplied by the user was added
C
C     Modifications introduced from May 31th to November 2nd of 2001
C     in occasion of the ALGENCAN development:
C
C     Fixed bugs:
C
C     1) The first spectral steplength was not been projected in the
C     [lspgmi,lspgma] interval.
C
C     2) The conjugate gradients accuracy (cgeps) which is linearly
C     dependent of the Euclidian norm of the projected gradient, was
C     also not been projected in the interval [cgepsi,cgepsf].
C
C     3) Conjugate gradients said that it was being used an Euclidian
C     norm trust region when it has really being used an infinite norm
C     trust region and viceversa.
C
C     4) Sometimes, the analytic gradient has been used although the
C     user choose the finite differences option.
C
C     Modifications:
C
C     1) To avoid roundoff errors, an explicit detection of at least one
C     variable reaching its bound when a maximum step is being made was
C     added.
C
C     2) The way in which two points were considered very similar in, 
C     for example, the interpolations and the extrapolations (which was 
C     dependent of the infinity norm of the points) showed to be very
C     scale dependent. A new version which test the difference 
C     coordinate to coordinate was done. In this was the calculus of the 
C     current point x and the descent direction sup-norm is not done any
C     more.
C
C     3) The same constants epsrel and epsabs were used as small 
C     relative and absolute values for, for example, detecting similar
C     points and for finite differences. Now, epsrel and epsabs are used 
C     for detecting similar points (and the recommended values are 
C     10^{-10} and 10^{-20}, respectively) and new constants sterel and 
C     steabs were introduced for finite differences (and the recommended 
C     values are 10^{-7} and 10^{-10}, respectively).
C
C     4) Two new stopping criteria for CG were added: (i) we stop if
C     two consecutive iterates are too  close; and (ii) we also
C     stop if there is no enough quadratic model progress during
C     maxitnqmp iterations.
C
C     5) The linear relation between the conjugate gradient accuracy
C     and the norm of the projected gradient can be computed using
C     the Euclidian- and the sup-norm of the projected gradient (only
C     Euclidian norm version was present in the previous version. The
C     linear relation is such that the CG accuracy is cgepsi when the
C     projected gradient norm value is equal to the value corresponding
C     to the initial guess and the CG accuracy is cgepsf when the
C     projected gradient norm value is cgrelf).
C
C     6) Inside Conjugate Gradients, the Euclidian-norm is been computed 
C     using an algorithm developed by C.L.LAWSON, 1978 JAN 08. Numerical 
C     experiments showed that the performance of GENCAN depends 
C     basically on the conjugate gradients performance and stopping
C     criteria and that the conjugate gradients depends on the way the
C     Euclidian-norm is been computed. These things deserve further 
C     research.
C
C     7) In the Augmented Lagrangian algorithm ALGENCAN, which uses
C     GENCAN to solve the bounded constrained subproblems, the maximum
C     number of Conjugate Gradients iterations (cgmaxit), which in this
C     version is linearly dependent of the projected gradient norm, was 
C     set to 2 * (# of free variables). As CG is not using restarts we 
C     do not know very well what this means. On the other hand, the 
C     accuracy (given by cgeps) continues being more strict when we are 
C     near to the solution and less strict when we ar far from the 
C     solution. 
c
