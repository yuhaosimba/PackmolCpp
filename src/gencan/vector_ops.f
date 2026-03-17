      subroutine shrink(nind,ind,n,v)

      implicit none

C     SCALAR ARGUMENTS
      integer n,nind

C     ARRAY ARGUMENTS
      integer ind(nind)
      double precision v(n)

C     This subroutine shrinks vector v from the full dimension space 
C     (dimension n) to the reduced space (dimension nind).
C
C     On entry:
C
C     nind     integer 
C              dimension of the reduced space
C
C     ind      integer ind(nind)
C              components ind(1)-th, ..., ind(nind)-th are the
C              components that belong to the reduced space
C
C     n        integer
C              dimension of the full space
C
C     v        double precision v(n)
C              vector to be shrinked
C
C     On Return
C
C     v        double precision v(n)
C              shrinked vector

C     LOCAL SCALARS
      integer i,indi
      double precision tmp

      do i = 1,nind
           indi = ind(i)
           if ( i .ne. indi ) then
               tmp     = v(indi)
               v(indi) = v(i)
               v(i)    = tmp
          end if
      end do

      return

      end

C     ******************************************************************
C     ******************************************************************

      subroutine expand(nind,ind,n,v)

      implicit none

C     SCALAR ARGUMENTS
      integer n, nind

C     ARRAY ARGUMENTS
      integer ind(nind)
      double precision v(n)

C     This subroutine expands vector v from the reduced space 
C     (dimension nind) to the full space (dimension n).
C
C     On entry:
C
C     nind     integer 
C              dimension of the reduced space
C
C     ind      integer ind(nind)
C              components ind(1)-th, ..., ind(nind)-th are the
C              components that belong to the reduced space
C
C     n        integer
C              dimension of the full space
C
C     v        double precision v(n)
C              vector to be expanded
C
C     On Return
C
C     v        double precision v(n)
C              expanded vector

C     LOCAL SCALARS
      integer i,indi
      double precision tmp

      do i = nind,1,- 1
          indi = ind(i)
          if ( i .ne. indi ) then
              tmp     = v(indi)
              v(indi) = v(i)
              v(i)    = tmp
          end if
      end do
     
      return

      end
 
C     ******************************************************************
C     ******************************************************************

      subroutine evalnaldiff(n,x,m,lambda,rho,g,sterel,steabs,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer n,m,inform
      double precision sterel,steabs

C     ARRAY ARGUMENTS
      double precision x(n),lambda(m),rho(m),g(n)

C     Approximates the gradient vector g(x) of the objective function by 
C     central finite differences. This subroutine, which works in the 
C     full space, is prepared to replace the subroutine evalnal (to 
C     evaluate the gradient vector) in the case of the lastest have not
C     being provided by the user.
C
C     On entry:
C
C     n        integer
C              number of variables
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
C     sterel   double precision
C              See below
C
C     steabs   double precision
C              This constants mean a ''relative small number'' and ''an 
C              absolute small number'' for the increments in finite 
C              difference approximations of derivatives
C
C              RECOMMENDED: epsrel = 1.0d-07 and epsabs = 1.0d-10 
C
C              CONSTRAINTS: sterel >= steabs > 0
C
C     On Return
C
C     g        double precision g(n)
C              approximation of the gradient vector at x
C
C     inform   integer
C              0 = no errors,
C            < 0 = there was an error in the gradient calculation.

C     LOCAL SCALARS
      integer j
      double precision tmp,step,fplus,fminus

      inform = 0

      do j = 1,n
          tmp  = x(j)
          step = max( steabs, sterel * abs( tmp ) )

          x(j) = tmp + step
          call evalal(n,x,m,lambda,rho,fplus,inform)
          if ( inform .lt. 0 ) then
              return
          end if

          x(j) = tmp - step
          call evalal(n,x,m,lambda,rho,fminus,inform)
          if ( inform .lt. 0 ) then
              return
          end if

          g(j) = ( fplus - fminus ) / ( 2.0d0 * step )
          x(j) = tmp
      end do

      return

      end

C     *****************************************************************
C     *****************************************************************

