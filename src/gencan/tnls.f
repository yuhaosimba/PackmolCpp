      subroutine tnlsf(nind,ind,n,x,m,lambda,rho,l,u,f,g,d,amax,rbdtype,
     +rbdind,nint,next,mininterp,maxextrap,fmin,maxfc,gtype,iprint,fcnt,
     +gcnt,intcnt,exgcnt,exbcnt,inform,xplus,xtmp,xbext,gamma,beta,
     +sigma1,sigma2,sterel,steabs,epsrel,epsabs,infrel,infabs)

      implicit none

C     SCALAR ARGUMENTS
      integer exbcnt,exgcnt,fcnt,gcnt,gtype,inform,intcnt,iprint,m,
     +        maxextrap,maxfc,mininterp,n,nind,rbdind,rbdtype
      double precision amax,beta,epsabs,epsrel,f,fmin,gamma,infabs,
     +        infrel,next,nint,sigma1,sigma2,steabs,sterel

C     ARRAY ARGUMENTS
      integer ind(nind)
      double precision d(n),g(n),l(n),lambda(m),rho(m),u(n),x(n),
     +        xbext(n),xplus(n),xtmp(n)

C     This subroutine implements the line search used in the Truncated
C     Newton direction.
C
C     On Entry
C
C     nind     integer
C              number of free variables (this is thee dimension in 
C              which this subroutine will work)
C
C     ind      integer ind(n)
C              array which contains, in the first nind positions, the
C              identifiers of the free variables
C
C     n        integer
C              dimension of the full space
C
C     x        double precision x(n)
C              current point
C
C              The first nind positions of x contains the free variables 
C              x_ind(1), x_ind(2), ..., x_ind(nind).
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
C     l        double precision l(nind)
C              lower bounds on x. It components are ordered in the
C              same way as x and g.
C
C     u        double precision u(nind)
C              upper bounds on x. It components are ordered in the
C              same way as x, g and l.
C
C     f        double precision
C              functional value at x
C
C     g        double precision g(n)
C              gradient vector at x
C
C              It also contains in the first nind positions the 
C              components g_ind(1), g_ind(2), ..., g_ind(nind).
C    
C              IMPORTANT: the linear algebra of this subroutine lies in 
C              a space of dimension nind. The value of the full 
C              dimension n, the non-free variables (which are at the end 
C              of array x) and its gradient components (which are at the 
C              end of array g) are also used and updated any time the 
C              gradient is being computed.
C
C     d        double precision d(nind)
C              descent direction 
C    
C     amax     double precision
C
C     rbdtype  integer
C
C     rbdind   integer
C
C     nint     double precision
C              constant for the interpolation. See the description of
C              sigma1 and sigma2 above. Sometimes we take as a new 
C              trial step the previous one divided by nint
C
C              RECOMMENDED: nint = 2.0
C
C     next     double precision
C              constant for the extrapolation
C              when extrapolating we try alpha_new = alpha * next
C
C              RECOMMENDED: next = 2.0
C
C     mininterp integer
C              constant for testing if, after having made at least 
C              mininterp interpolations, the steplength is so small. 
C              In that case failure of the line search is declared (may 
C              be the direction is not a descent direction due to an 
C              error in the gradient calculations) 
C
C              RECOMMENDED: mininterp = 4
C
C     maxextrap integer
C              constant to limit the number of extrapolations
C
C              RECOMMENDED: maxextrap = 1000 (a big number)
C
C     fmin     double precision
C              functional value for the stopping criteria f <= fmin
C
C     maxfc    integer
C              maximum number of functional evaluations
C
C     gtype    integer
C              type of gradient calculation
C              gtype = 0 means user suplied evalg subroutine,
C              gtype = 1 means central difference approximation.
C
C              RECOMMENDED: gtype = 0
C
C              (provided you have the evalg subroutine)
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
C     xplus    double precision xplus(nind)
C     xtmp     double precision xtmp(nind)
C     xbext    double precision xbext(nind) 
C              working vectors
C
C     gamma    double precision
C              constant for the Armijo criterion
C              f(x + alpha d) <= f(x) + gamma * alpha * <\nabla f(x),d>
C
C              RECOMMENDED: gamma = 10^{-4}
C
C     beta     double precision
C              constant for the beta condition <dk, g(xk + dk)>  <  beta 
C              * <dk,gk>. If (xk + dk) satisfies the Armijo condition 
C              but does not satisfy the beta condition then the point is 
C              accepted, but if it satisfied the Armijo condition and 
C              also satisfies the beta condition then we know that there 
C              is the possibility for a successful extrapolation
C
C              RECOMMENDED: beta = 0.5
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
C              negligible with respect to another quantity B if 
C              |A| < max ( epsrel * |B|, epsabs ) 
C
C              RECOMMENDED: epsrel = 10^{-10}, epsabs = 10^{-20}, 
C                           infrel = 10^{+20}, infabs = 10^{+99}
C
C     On Return
C
C     x        double precision x(n)
C              new current point
C
C     f        double precision
C              functional value at x
C
C     g        double precision g(n)
C              gradient vector at x
C
C     fcnt     integer
C              number of functional evaluations used in this line search
C
C     gcnt     integer
C              number of gradient evaluations used in this line search
C
C     intcnt   integer
C              number of interpolations
C
C     exgcnt   integer
C              number of good extrapolations
C
C     exbcnt   integer
C              number of bad extrapolations
C
C     inform   integer
C              This output parameter tells what happened in this 
C              subroutine, according to the following conventions:
C
C              0 = convergence with an Armijo-like criterion
C                  (f(xnew) <= f(x) + 1.0d-4 * alpha * <g,d>);
C
C              4 = the algorithm stopped because the functional value
C                  is very small (f <= fmin);
C
C              6 = so small step in the line search. After having made 
C                  at least mininterp interpolations, the steplength 
C                  becames small. ``small steplength'' means that we are 
C                  at point x with direction d and step alpha, and, for 
C                  all i, 
C
C                  |alpha * d(i)| .le. max ( epsrel * |x(i)|, epsabs ). 
C 
C                  In that case failure of the line search is declared 
C                  (may be the direction is not a descent direction 
C                  due to an error in the gradient calculations). Use
C                  mininterp > maxfc for inhibit this criterion;
C
C              8 = it was achieved the maximum allowed number of
C                  function evaluations (maxfc);
C
C            < 0 = error in evalf or evalg subroutines.

C     LOCAL SCALARS
      logical samep
      integer extrap,i,interp
      double precision alpha,atmp,fbext,fplus,ftmp,gptd,gtd

C     ==================================================================
C     Initialization 
C     ==================================================================

C     ==================================================================
C     Compute directional derivative
C     ==================================================================

      gtd = 0.0d0
      do i = 1,nind
          gtd = gtd + g(i) * d(i)
      end do

C     ==================================================================
C     Compute first trial
C     ==================================================================

      alpha = min( 1.0d0, amax )

      do i = 1,nind
          xplus(i) = x(i) + alpha * d(i)
      end do

      if ( alpha .eq. amax ) then
          if ( rbdtype .eq. 1 ) then
              xplus(rbdind) = l(rbdind)
          else ! if (rbdtype.eq.2) then
              xplus(rbdind) = u(rbdind)
          end if
      end if

      call calcf(nind,ind,xplus,n,x,m,lambda,rho,fplus,inform)
      fcnt = fcnt + 1

      if ( inform .lt. 0 ) then
 
          if ( iprint .ge. 4 ) then
              write(*, 1000) inform
              write(10,1000) inform
          end if

          return

      end if

C     Print initial information

      if ( iprint .ge. 4 ) then
          write(*, 980) amax
          write(*, 999) alpha,fplus,fcnt

          write(10,980) amax
          write(10,999) alpha,fplus,fcnt
      end if

C     ==================================================================
C     Test Armijo and beta-condition and decide for accepting the trial 
C     point, interpolate or extrapolate.
C     ==================================================================

      if ( amax .gt. 1.0d0 ) then

C         x + d belongs to the interior of the feasible set
          if ( iprint .ge. 4 ) then
              write(*, *) '     x+d belongs to int of the feasible set'
              write(10,*) '     x+d belongs to int of the feasible set'
          end if

C         Verify Armijo

          if ( fplus .le. f + gamma * alpha * gtd ) then

C             Armijo condition holds  
              if ( iprint .ge. 4 ) then
                  write(*, *) '     Armijo condition holds' 
                  write(10,*) '     Armijo condition holds' 
              end if

              if ( gtype .eq. 0 ) then
                  call calcg(nind,ind,xplus,n,x,m,lambda,rho,g,inform)
              else if ( gtype .eq. 1 ) then
                  call calcgdiff(nind,ind,xplus,n,x,m,lambda,rho,g,
     +            sterel,steabs,inform)
              end if
              gcnt = gcnt + 1

              if ( inform .lt. 0 ) then

                  if ( iprint .ge. 4 ) then
                      write(*, 1000) inform
                      write(10,1000) inform
                  end if

                  return

              end if

              gptd = 0.0d0
              do i = 1,nind
                  gptd = gptd + g(i) * d(i)
              end do

C             Verify directional derivative (beta condition)

              if ( gptd .lt. beta * gtd ) then

C                 Extrapolate
                  if ( iprint .ge. 4 ) then
                      write(*, *)'     The beta-condition does not hold'
                      write(*, *)'     We will extrapolate'
                      write(10,*)'     The beta-condition does not hold'
                      write(10,*)'     We will extrapolate'
                  end if

C                 f and x before extrapolation
                  fbext = fplus

                  do i = 1,nind
                      xbext(i) = xplus(i)
                  end do

                  go to 100

              else

C                 Step = 1 was ok, finish the line search
                  if ( iprint .ge. 4 ) then
                      write(*, *) '     The beta condition is also true'
                      write(*, *) '     Line search is over'
                      write(10,*) '     The beta condition is also true'
                      write(10,*) '     Line search is over'
                  end if

                  f = fplus

                  do i = 1,nind
                      x(i) = xplus(i)
                  end do

                  inform = 0

                  if ( iprint .ge. 4 ) then
                      write(*, 990) inform
                      write(10,990) inform
                  end if

                  go to 500

              end if

          else 

C             Interpolate
              if ( iprint .ge. 4 ) then
                  write(*, *) '     Armijo does not hold'
                  write(*, *) '     We will interpolate'
                  write(10,*) '     Armijo does not hold'
                  write(10,*) '     We will interpolate'
              end if

              go to 200

          end if

      else

C         x + d does not belong to the feasible set (amax <= 1)
          if ( iprint .ge. 4 ) then
              write(*, *) '     x+d does not belong to box-interior'
              write(10,*) '     x+d does not belong to box-interior'
          end if

          if ( fplus .lt. f ) then

C             Extrapolate
              if ( iprint .ge. 4 ) then
                  write(*, *) '     f(x+d) < f(x)'
                  write(*, *) '     We will extrapolate'
                  write(10,*) '     f(x+d) < f(x)'
                  write(10,*) '     We will extrapolate'
              end if

C             f and x before extrapolation
              fbext = fplus

              do i = 1,nind
                  xbext(i) = xplus(i)
              end do

              go to 100

          else 

C             Interpolate
              if ( iprint .ge. 4 ) then
                  write(*, *) '     f(x+d) >= f(x)'
                  write(*, *) '     We will interpolate'
                  write(10,*) '     f(x+d) >= f(x)'
                  write(10,*) '     We will interpolate'
              end if

              go to 200

          end if

      end if


C     ==================================================================
C     Extrapolation
C     ==================================================================

 100  continue

      extrap = 0

C     Test f going to -inf

 120  if ( fplus .le. fmin ) then

C         Finish the extrapolation with the current point

          f = fplus

          do i = 1,nind
              x(i) = xplus(i)
          end do

          if ( extrap .ne. 0 .or. amax .le. 1.0d0 ) then

              if ( gtype .eq. 0 ) then
                  call calcg(nind,ind,x,n,x,m,lambda,rho,g,inform)
              else if ( gtype .eq. 1 ) then
                  call calcgdiff(nind,ind,x,n,x,m,lambda,rho,g,sterel,
     +            steabs,inform)
              end if
              gcnt = gcnt + 1

              if ( inform .lt. 0 ) then
 
                  if ( iprint .ge. 4 ) then
                      write(*, 1000) inform
                      write(10,1000) inform
                  end if

                  return

              end if

              if ( f .lt. fbext ) then
                  exgcnt = exgcnt + 1
              else
                  exbcnt = exbcnt + 1
              end if

          end if

          inform = 4

          if ( iprint .ge.3 ) then
              write(*, 994) inform
              write(10,994) inform
          end if

          go to 500

      end if

C     Test maximum number of functional evaluations

      if ( fcnt .ge. maxfc ) then

C         Finish the extrapolation with the current point

          f = fplus

          do i = 1,nind
              x(i) = xplus(i)
          end do

C         If extrap=0 and amax>1 the gradient was computed for testing 
C         the beta condition and it is not necessary to compute it again
          if ( extrap .ne. 0 .or. amax .le. 1.0d0 ) then

              if ( gtype .eq. 0 ) then
                  call calcg(nind,ind,x,n,x,m,lambda,rho,g,inform)
              else if ( gtype .eq. 1 ) then
                  call calcgdiff(nind,ind,x,n,x,m,lambda,rho,g,sterel,
     +            steabs,inform)
              end if
              gcnt = gcnt + 1

              if ( inform .lt. 0 ) then

                  if ( iprint .ge. 4 ) then
                      write(*, 1000) inform
                      write(10,1000) inform
                  end if

                  return

              end if

              if ( f .lt. fbext ) then
                  exgcnt = exgcnt + 1
              else
                  exbcnt = exbcnt + 1
              end if

          end if

          inform = 8

          if ( iprint .ge. 4 ) then
              write(*, 998) inform
              write(10,998) inform
          end if

          go to 500

      end if

C     Test if the maximum number of extrapolations was exceeded

      if ( extrap .ge. maxextrap ) then

C         Finish the extrapolation with the current point

          f = fplus

          do i = 1,nind
              x(i) = xplus(i)
          end do

C         If extrap=0 and amax>1 the gradient was computed for testing 
C         the beta condition and it is not necessary to compute it again
          if ( extrap .ne. 0 .or. amax .le. 1.0d0 ) then

              if ( gtype .eq. 0 ) then
                  call calcg(nind,ind,x,n,x,m,lambda,rho,g,inform)
              else if ( gtype .eq. 1 ) then
                  call calcgdiff(nind,ind,x,n,x,m,lambda,rho,g,sterel,
     +            steabs,inform)
              end if
              gcnt = gcnt + 1

              if ( inform .lt. 0 ) then

                  if ( iprint .ge. 4 ) then
                      write(*, 1000) inform
                      write(10,1000) inform
                  end if

                  return

              end if

              if ( f .lt. fbext ) then
                  exgcnt = exgcnt + 1
              else
                  exbcnt = exbcnt + 1
              end if

          end if

          inform = 7

          if ( iprint .ge. 4 ) then
              write(*, 997) inform
              write(10,997) inform
          end if

          go to 500

      end if

C     Chose new step 

      if ( alpha .lt. amax .and. next * alpha .gt. amax ) then
          atmp = amax
      else
          atmp = next * alpha
      end if

C     Compute new trial point

      do i = 1,nind
          xtmp(i) = x(i) + atmp * d(i)
      end do

      if ( atmp .eq. amax ) then
          if ( rbdtype .eq. 1 ) then
              xtmp(rbdind) = l(rbdind)
          else ! if ( rbdtype .eq. 2 ) then
              xtmp(rbdind) = u(rbdind)
          end if
      end if

C     Project

      if ( atmp .gt. amax ) then
          do i = 1,nind
              xtmp(i) = max( l(i), min( xtmp(i), u(i) ) )
          end do
      end if

C     Test if this is not the same point as the previous one.
C     This test is performed only when alpha > amax.

      if( alpha .gt. amax ) then

          samep = .true.
          do i = 1,nind
              if ( abs( xtmp(i) - xplus(i) ) .gt.
     +             max( epsrel * abs( xplus(i) ), epsabs ) ) then
                  samep = .false.
              end if
          end do

          if ( samep ) then

C             Finish the extrapolation with the current point

              f = fplus

              do i = 1,nind
                  x(i) = xplus(i)
              end do

C             If extrap=0 and amax>1 the gradient was computed for 
C             testing the beta condition and it is not necessary to 
C             compute it again
              if ( extrap .ne. 0 .or. amax .le. 1.0d0 ) then

                  if ( gtype .eq. 0 ) then
                      call calcg(nind,ind,x,n,x,m,lambda,rho,g,inform)
                  else if ( gtype .eq. 1 ) then
                      call calcgdiff(nind,ind,x,n,x,m,lambda,rho,g,
     +                sterel,steabs,inform)
                  end if
                  gcnt = gcnt + 1

                  if ( inform .lt. 0 ) then

                      if ( iprint .ge. 4 ) then
                          write(*, 1000) inform
                          write(10,1000) inform
                      end if

                      return

                  end if

                  if ( f .lt. fbext ) then
                      exgcnt = exgcnt + 1
                  else
                      exbcnt = exbcnt + 1
                  end if

              end if

              inform = 0

              if ( iprint .ge. 4 ) then
                  write(*, 990) inform
                  write(10,990) inform
              end if

              go to 500

          end if

      end if

C     Evaluate function

      call calcf(nind,ind,xtmp,n,x,m,lambda,rho,ftmp,inform)
      fcnt = fcnt + 1

      if ( inform .lt. 0 ) then

C         if ( iprint .ge. 4 ) then
C             write(*, 1000) inform
C             write(10,1000) inform
C         end if

C         return

C         If the objective function is not well defined in an 
C         extrapolated point, we discard all the extrapolated points
C         and return to a safe region (where the point before
C         starting the extrapolations is)

          f = fbext

          do i = 1,nind
              x(i) = xbext(i)
          end do

C         If extrap=0 and amax>1 the gradient was computed for testing 
C         the beta condition and it is not necessary to compute it again
          if ( extrap .ne. 0 .or. amax .le. 1.0d0 ) then

              if ( gtype .eq. 0 ) then
                  call calcg(nind,ind,x,n,x,m,lambda,rho,g,inform)
              else if ( gtype .eq. 1 ) then
                  call calcgdiff(nind,ind,x,n,x,m,lambda,rho,g,sterel,
     +            steabs,inform)
              end if
              gcnt = gcnt + 1

              if ( inform .lt. 0 ) then

                  if ( iprint .ge. 4 ) then
                      write(*, 1000) inform
                      write(10,1000) inform
                  end if

                  return

              end if

              exbcnt = exbcnt + 1

          end if

          inform = 0

          if ( iprint .ge. 4 ) then
              write(*, 1010) inform
              write(10,1010) inform
          end if

          go to 500

      end if

C     Print information of this iteration

      if ( iprint .ge. 4 ) then
          write(*, 999) atmp,ftmp,fcnt
          write(10,999) atmp,ftmp,fcnt
      end if

C     If the functional value decreases then set the current point and 
C     continue the extrapolation

      if ( ftmp .lt. fplus ) then

          alpha = atmp

          fplus = ftmp

          do i = 1,nind
              xplus(i) = xtmp(i)
          end do

          extrap = extrap + 1

          go to 120

C     If the functional value does not decrease then discard the last 
C     trial and finish the extrapolation with the previous point

      else

          f = fplus

          do i = 1,nind
              x(i) = xplus(i)
          end do

C         If extrap=0 and amax>1 the gradient was computed for testing 
C         the beta condition and it is not necessary to compute it again
          if ( extrap .ne. 0 .or. amax .le. 1.0d0 ) then

              if ( gtype .eq. 0 ) then
                  call calcg(nind,ind,x,n,x,m,lambda,rho,g,inform)
              else if ( gtype .eq. 1 ) then
                  call calcgdiff(nind,ind,x,n,x,m,lambda,rho,g,sterel,
     +            steabs,inform)
              end if
              gcnt = gcnt + 1

              if ( inform .lt. 0 ) then

                  if ( iprint .ge. 4 ) then
                      write(*, 1000) inform
                      write(10,1000) inform
                  end if

                  return

              end if

              if ( f .lt. fbext ) then
                  exgcnt = exgcnt + 1
              else
                  exbcnt = exbcnt + 1
              end if

          end if

          inform = 0

          if ( iprint .ge.3 ) then
              write(*, 990) inform
              write(10,990) inform
          end if

          go to 500

      end if
C     ==================================================================
C     End of extrapolation
C     ==================================================================

C     ==================================================================
C     Interpolation
C     ==================================================================

 200  continue

      intcnt = intcnt + 1

      interp = 0

 210  continue

C     Test f going to -inf

      if ( fplus .le. fmin ) then

C         Finish the interpolation with the current point

          f = fplus

          do i = 1,nind
              x(i) = xplus(i)
          end do

          if ( gtype .eq. 0 ) then
              call calcg(nind,ind,x,n,x,m,lambda,rho,g,inform)
          else if ( gtype .eq. 1 ) then
              call calcgdiff(nind,ind,x,n,x,m,lambda,rho,g,sterel,
     +        steabs,inform)
          end if
          gcnt = gcnt + 1

          if ( inform .lt. 0 ) then

              if ( iprint .ge. 4 ) then
                  write(*, 1000) inform
                  write(10,1000) inform
              end if

              return

          end if

          inform = 4

          if ( iprint .ge. 4 ) then
              write(*, 994) inform
              write(10,994) inform
          end if

          go to 500

      end if

C     Test maximum number of functional evaluations

      if ( fcnt .ge. maxfc ) then

C         As this is an abrupt termination then the current point of the 
C         interpolation may be worst than the initial one

C         If the current point is better than the initial one then
C         finish the interpolation with the current point else discard
C         all we did inside this line search and finish with the initial
C         point

          if ( fplus .lt. f ) then

              f = fplus

              do i = 1,nind
                  x(i) = xplus(i)
              end do

              if ( gtype .eq. 0 ) then
                  call calcg(nind,ind,x,n,x,m,lambda,rho,g,inform)
              else if ( gtype .eq. 1 ) then
                  call calcgdiff(nind,ind,x,n,x,m,lambda,rho,g,sterel,
     +            steabs,inform)
              end if
              gcnt = gcnt + 1

              if ( inform .lt. 0 ) then

                  if ( iprint .ge. 4 ) then
                      write(*, 1000) inform
                      write(10,1000) inform
                  end if

                  return

              end if

          end if

          inform = 8

          if ( iprint .ge. 4 ) then
              write(*, 998) inform
              write(10,998) inform
          end if

          go to 500

      end if

C     Test Armijo condition

      if ( fplus .le. f + gamma * alpha * gtd ) then

C         Finish the line search
      
          f = fplus

          do i = 1,nind
              x(i) = xplus(i)
          end do

          if ( gtype .eq. 0 ) then
              call calcg(nind,ind,x,n,x,m,lambda,rho,g,inform)
          else if ( gtype .eq. 1 ) then
              call calcgdiff(nind,ind,x,n,x,m,lambda,rho,g,sterel,
     +        steabs,inform)
          end if
          gcnt = gcnt + 1

          if ( inform .lt. 0 ) then

              if ( iprint .ge. 4 ) then
                  write(*, 1000) inform
                  write(10,1000) inform
              end if

              return

          end if

          inform = 0

          if ( iprint .ge. 4 ) then
              write(*, 990) inform
              write(10,990) inform
          end if

          go to 500

      end if

C     Compute new step

      interp = interp + 1

      if ( alpha .lt. sigma1 ) then
          alpha = alpha / nint      

      else
          atmp = ( - gtd * alpha **2 ) / 
     +           (2.0d0 * ( fplus - f - alpha * gtd ) )

          if ( atmp .lt. sigma1 .or. atmp .gt. sigma2 * alpha ) then
              alpha = alpha / nint

          else
              alpha = atmp
          end if
      end if

C     Compute new trial point

      do i = 1,nind
          xplus(i) = x(i) + alpha * d(i)
      end do

      call calcf(nind,ind,xplus,n,x,m,lambda,rho,fplus,inform)
      fcnt = fcnt + 1

      if ( inform .lt. 0 ) then

          if ( iprint .ge. 4 ) then
              write(*, 1000) inform
              write(10,1000) inform
          end if

          return

      end if

C     Print information of this iteration

      if ( iprint .ge. 4 ) then
          write(*, 999) alpha,fplus,fcnt
          write(10,999) alpha,fplus,fcnt
      end if

C     Test whether at least mininterp interpolations were made and two 
C     consecutive iterates are much close

      samep = .true.
      do i = 1,nind
         if ( abs( alpha * d(i) ) .gt. 
     +        max( epsrel * abs( x(i) ), epsabs ) ) then
             samep = .false.
         end if
      end do

      if ( interp .ge. mininterp .and. samep ) then

C         As this is an abrupt termination then the current point of the 
C         interpolation may be worst than the initial one

C         If the current point is better than the initial one then
C         finish the interpolation with the current point else discard 
C         all we did inside this line search and finish with the initial 
C         point

C         if ( fplus .lt. f ) then

C             f = fplus

C             do i = 1,nind
C                 x(i) = xplus(i)
C             end do

C             if ( gtype .eq. 0 ) then
C                 call calcg(nind,ind,x,n,x,m,lambda,rho,g,inform)
C             else if ( gtype .eq. 1 ) then
C                 call calcgdiff(nind,ind,x,n,x,m,lambda,rho,g,
c    +            sterel,steabs,inform)
C             end if
C             gcnt = gcnt + 1

C             if ( inform .lt. 0 ) then 

C                 if ( iprint .ge. 4 ) then
C                     write(*, 1000) inform
C                     write(10,1000) inform
C                 end if

C                 return

C             end if

C         end if

C         The previous lines were commented because, as it is been used, 
C         this subroutine must return with the initial point in case of 
C         finding a very small interpolation step. From that initial 
C         point, something different will be tried.

          inform = 6

          if ( iprint .ge. 4 ) then
              write(*, 996) inform
              write(10,996) inform
          end if
  
          go to 500

      end if

C     Else, iterate

      go to 210
C     ==================================================================
C     End of interpolation
C     ==================================================================

 500  continue

C     ==================================================================
C     Return
C     ==================================================================

      return

C     Non-executable statements

 980  format(/,6X,'TN Line search (alphamax= ',1PD11.4,')')
 999  format(6X,'Alpha= ',1PD11.4,' F= ',1PD11.4,' FE= ',I5)
 990  format(6X,'Flag of TN Line search= ',I3,
     +          ' (Convergence with an Armijo-like criterion)')
 994  format(6X,'Flag of TN Line search= ',I3,
     +          ' (Small functional value, smaller than ',/,
     +       6X,'parameter fmin)')
 996  format(6X,'Flag of TN Line search= ',I3,
     +          ' (Too small step in the interpolation)')
 997  format(6X,'Flag of TN Line search= ',I3,
     +          ' (Too many extrapolations)')
 998  format(6X,'Flag of TN Line search= ',I3,
     +          ' (Too many functional evaluations)')
 1000 format(6X,'Flag of TN Line search = ',I3,' Fatal Error')
 1010 format(6X,'Flag of TN Line search= ',I3,
     +          ' (Fatal Error in an extrapolated point)')

      end

C     ******************************************************************
C     ******************************************************************
