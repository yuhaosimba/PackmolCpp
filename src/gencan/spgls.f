      subroutine spglsf(n,x,m,lambda,rho,f,g,l,u,lamspg,nint,mininterp,
     +fmin,maxfc,iprint,fcnt,inform,xtrial,d,gamma,sigma1,sigma2,sterel,
     +steabs,epsrel,epsabs,infrel,infabs) 

      implicit none

C     SCALAR ARGUMENTS
      integer fcnt,m,maxfc,mininterp,n,inform,iprint
      double precision epsabs,epsrel,f,fmin,gamma,infrel,infabs,lamspg,
     +        nint,sigma1,sigma2,steabs,sterel

C     ARRAY ARGUMENTS
      double precision d(n),g(n),l(n),lambda(m),rho(m),u(n),x(n),
     +        xtrial(n)
 
C     Safeguarded quadratic interpolation, used in the Spectral 
C     Projected Gradient directions.
C
C     On Entry
C
C     n        integer
C              the order of the x
C
C     x        double precision x(n)
C              current point
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
C     f        double precision
C              function value at the current point
C
C     g        double precision g(n)
C              gradient vector at the current point
C
C     l        double precision l(n)
C              lower bounds
C
C     u        double precision u(n)
C              upper bounds
C
C     lamspg   double precision
C              spectral steplength
C
C     nint     double precision
C              constant for the interpolation. See the description of
C              sigma1 and sigma2 above. Sometimes we take as a new 
C              trial step the previous one divided by nint
C
C              RECOMMENDED: nint = 2.0
C
C     mininterp integer
C              constant for testing if, after having made at least 
C              mininterp interpolations, the steplength is so small. In 
C              that case failure of the line search is declared (may be 
C              the direction is not a descent direction due to an error 
C              in the gradient calculations) 
C
C              RECOMMENDED: mininterp = 4
C
C     fmin     double precision
C              functional value for the stopping criterion f <= fmin
C
C     maxfc    integer
C              maximum number of functional evaluations
C
C     iprint   integer
C              Commands printing. Nothing is printed if iprint is 
C              smaller than 2. If iprint is greater than or equal to 
C              2, GENCAN iterations information is printed. If iprint 
C              is greater than or equal to 3, line searches and 
C              Conjugate Gradients information is printed.
C
C              RECOMMENDED: iprint = 2
C
C              CONSTRAINTS: allowed values are just 2 or 3.
C
C     xtrial   double precision xtrial(n)
C     d        double precision d(n)
C              working vectors
C
C     gamma    double precision
C              constant for the Armijo criterion
C              f(x + alpha d) <= f(x) + gamma * alpha * <\nabla f(x),d>
C
C              RECOMMENDED: gamma = 10^{-4}
C
C     sigma1   double precision
C     sigma2   double precision
C              constant for the safeguarded interpolation
C              if alpha_new \notin [sigma1, sigma*alpha] then we take
C              alpha_new = alpha / nint
C
C              RECOMMENDED: sigma1 = 0.1 and sigma2 = 0.9
C
C     sterel   double precision
C     steabs   double precision
C              this constants mean a ``relative small number'' and ``an 
C              absolute small number'' for the increments in finite
C              difference approximations of derivatives
C
C              RECOMMENDED: epsrel = 10^{-7}, epsabs = 10^{-10} 
C
C     epsrel   double precision
C     epsabs   double precision
C     infrel   double precision
C     infabs   double precision
C              this constants mean a ``relative small number'', ``an 
C              absolute small number'', and ``infinite or a very big
C              number''. Basically, a quantity A is considered 
C              negligible with respect to another quantity B if |A| < 
C              max ( epsrel * |B|, epsabs ) 
C
C              RECOMMENDED: epsrel = 10^{-10}, epsabs = 10^{-20}, 
C                           infrel = 10^{+20}, infabs = 10^{+99}
C
C     On Return
C
C     x        double precision
C              final estimation of the solution
C
C     f        double precision
C              functional value at the final estimation 
C
C     fcnt     integer
C              number of functional evaluations used in the line search   
C
C     inform   integer
C              This output parameter tells what happened in this 
C              subroutine, according to the following conventions:
C
C              0 = convergence with an Armijo-like criterion
C                  (f(xnew) <= f(x) + gamma * alpha * <g,d>);
C
C              4 = the algorithm stopped because the functional value
C                  is smaller than fmin;
C
C              6 = too small step in the line search. After having made 
C                  at least mininterp interpolations, the steplength 
C                  becames small. ''small steplength'' means that we are 
C                  at point x with direction d and step alpha, and, for 
C                  all i, 
C
C                  | alpha * d(i) | <= max ( epsrel * |x(i)|, epsabs ). 
C 
C                  In that case failure of the line search is declared 
C                  (maybe the direction is not a descent direction due
C                  to an error in the gradient calculations). Use
C                  mininterp > maxfc to inhibit this criterion;
C
C              8 = it was achieved the maximum allowed number of
C                  function evaluations (maxfc);
C
C            < 0 = error in evalf subroutine.

C     LOCAL SCALARS
      logical samep
      integer i,interp
      double precision alpha,atmp,ftrial,gtd

C     Print presentation information

      if ( iprint .ge. 4 ) then
          write(*, 980) lamspg
          write(10,980) lamspg
      end if

C     Initialization

      interp = 0

C     Compute first trial point, spectral projected gradient direction, 
C     and directional derivative <g,d>.

      alpha = 1.0d0

      gtd = 0.0d0
      do i = 1,n
          xtrial(i) = min( u(i), max( l(i), x(i) - lamspg * g(i) ) )
          d(i)      = xtrial(i) - x(i)
          gtd       = gtd + g(i) * d(i)
      end do

      call evalal(n,xtrial,m,lambda,rho,ftrial,inform)
      fcnt = fcnt + 1

      if ( inform .lt. 0 ) then

          if ( iprint .ge. 4 ) then
              write(*, 1000) inform
              write(10,1000) inform
          end if

          return

      end if

C     Print information of the first trial

      if ( iprint .ge. 4 ) then
          write(*, 999) alpha,ftrial,fcnt
          write(10,999) alpha,ftrial,fcnt
      end if

C     Main loop

 100  continue

C     Test Armijo stopping criterion

      if ( ftrial .le. f + gamma * alpha * gtd ) then

          f = ftrial

          do i = 1,n
              x(i) = xtrial(i)
          end do

          inform = 0

          if ( iprint .ge. 4 ) then
              write(*, 990) inform
              write(10,990) inform
          end if

          go to 500

      end if

C     Test whether f is very small

      if ( ftrial .le. fmin ) then

          f = ftrial

          do i = 1,n
              x(i) = xtrial(i)
          end do

          inform = 4

          if ( iprint .ge. 4 ) then
              write(*, 994) inform
              write(10,994) inform
          end if

          go to 500

      end if

C     Test whether the number of functional evaluations is exhausted

      if ( fcnt .ge. maxfc ) then

          if ( ftrial .lt. f ) then

              f = ftrial

              do i = 1,n
                  x(i) = xtrial(i)
              end do

          end if

          inform = 8

          if ( iprint .ge. 4 ) then
              write(*, 998) inform
              write(10,998) inform
          end if

          go to 500

      end if

C     Compute new step (safeguarded quadratic interpolation)

      interp = interp + 1

      if ( alpha .lt. sigma1 ) then
          alpha = alpha / nint      

      else
          atmp = ( - gtd * alpha ** 2 ) / 
     +           ( 2.0d0 * ( ftrial - f - alpha * gtd ) )

          if ( atmp .lt. sigma1 .or. atmp .gt. sigma2 * alpha ) then
              alpha = alpha / nint

          else
              alpha = atmp
          end if
      end if

C     Compute new trial point

      do i = 1,n
          xtrial(i) = x(i) + alpha * d(i)
      end do

      call evalal(n,xtrial,m,lambda,rho,ftrial,inform)
      fcnt = fcnt + 1

      if ( inform .lt. 0 ) then

          if ( iprint .ge. 4 ) then
              write(*, 1000) inform
              write(10,1000) inform
          end if

          return

      end if

C     Print information of the current trial

      if ( iprint .ge. 4 ) then
          write(*, 999) alpha,ftrial,fcnt
          write(10,999) alpha,ftrial,fcnt
      end if

C     Test whether at least mininterp interpolations were made and two 
C     consecutive iterates are close enough

      samep = .true.
      do i = 1,n
         if ( abs( alpha * d(i) ) .gt. 
     +        max( epsrel * abs( x(i) ), epsabs ) ) then
             samep = .false.
         end if
      end do

      if ( interp .ge. mininterp .and. samep ) then

          if ( ftrial .lt. f ) then

              f = ftrial

              do i = 1,n
                  x(i) = xtrial(i)
              end do

          end if

          inform = 6

          if ( iprint .ge. 4 ) then
              write(*, 996) inform
              write(10,996) inform
          end if
  
          go to 500

      end if

C     Iterate

      go to 100

C     Return

 500  continue

      return

C     Non-executable statements

 980  format(/,6x,'SPG (spectral steplength ',1PD11.4,')',/,/,
     *         6x,'SPG Line search')
 999  format(6x,'Alpha= ',1PD11.4,' F= ',1PD11.4,' FE= ',I5)
 990  format(6x,'Flag of SPG Line search = ',I3,
     *          ' (Convergence with an Armijo-like criterion)')
 994  format(6x,'Flag of SPG Line search = ',I3,
     *          ' (Small functional value, smaller than ',/,
     *       6X,'parameter fmin)')
 996  format(6x,'Flag of SPG Line search = ',I3,
     *          ' (Too small step in the interpolation)')
 998  format(6x,'Flag of SPG Line search = ',I3,
     *          ' (Too many functional evaluations)')
 1000 format(6x,'Flag of SPG Line search = ',I3,' Fatal Error')

      end

C     ******************************************************************
C     ******************************************************************
