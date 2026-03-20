C     *****************************************************************
C     *****************************************************************

      subroutine evalal(n,x,m,lambda,rho,f,flag)

C     This subroutine computes the objective function when GENCAN is
C     being used stand-alone to solve a unique bound-constrained problem. 
C     When GENCAN is being used in an Augmented Lagrangian framework, 
C     this subroutine must compute the Augmented Lagrangian function.
C
C     On Entry:
C
C     n     integer,
C           number of variables,
C
C     x     double precision x(n),
C           current point,
C
C     m     integer,
C           number of constraints (equalities plus inequalities),
C
C     lambda double precision lambdae(m),
C           current estimation of the Lagrange multipliers,
C
C     rho   double precision rho(m)
C           penalty parameters,
C
C     NOTE: arguments m, lambda and rho are useful when GENCAN is being used
C     for solving the box-constrained subproblems of an Augmented Lagrangian
C     framework. When GENCAN is being used stand-alone for solving a bound-
C     constrained problem, these arguments are dummy arguments.
C
C     On Return
C
C     f     double precision,
C           objective function value at x,
C
C     flag  integer
C           0 means "no errors",
C           1 means "some error occurs in the objective funtion evaluation".

      implicit none

C     SCALAR ARGUMENTS
      integer flag,m,n
      double precision f

C     ARRAY ARGUMENTS
      double precision lambda(m),rho(m),x(n)

C     LOCAL SCALARS

      flag = 0

      call computef(n,x,f)

      end

C     *****************************************************************
C     *****************************************************************

      subroutine evalnal(n,x,m,lambda,rho,g,flag)

C     This subroutine computes the gradient of the objective function 
C     when GENCAN is being used stand-alone to solve a unique bound-
C     constrained problem. When GENCAN is being used in an Augmented
C     Lagrangian framework, this subroutine must compute the gradient of
C     Augmented Lagrangian.
C
C     On Entry:
C
C     n     integer,
C           number of variables,
C
C     x     double precision x(n),
C           current point,
C
C     m     integer,
C           number of constraints (equalities plus inequalities),
C
C     lambda double precision lambdae(m),
C           current estimation of the Lagrange multipliers,
C
C     rho   double precision rho(m)
C           penalty parameters,
C
C     NOTE: arguments m, lambda and rho are useful when GENCAN is being used
C     for solving the box-constrained subproblems of an Augmented Lagrangian
C     framework. When GENCAN is being used stand-alone for solving a bound-
C     constrained problem, these arguments are dummy arguments.
C
C     On Return
C
C     g     double precision g(n),
C           gradient of the objective function at x,
C
C     flag  integer
C           0 means "no errors",
C           1 means "some error occurs in the gradient evaluation".

      implicit none

C     SCALAR ARGUMENTS
      integer flag,m,n

C     ARRAY ARGUMENTS
      double precision g(n),lambda(m),rho(m),x(n)

C     LOCAL SCALARS

      flag = 0

      call computeg(n,x,g)

      end

C     *****************************************************************
C     *****************************************************************

c Modified by L. Martinez (there was an error on the number of
c parameters when calling this subroutine). This subroutine does
c nothing.
c      subroutine evalhd(nind,ind,n,x,m,lambda,rho,d,hd,flag)

      subroutine evalhd(n)

C     This subroutine computes the product of the Hessian matrix times
C     the input vector argument d. If GENCAN is being used stand-alone 
C     to solve a bound-constrained problem, the ''Hessian matrix'' must
C     be the Hessian matrix of the objective function. On the other hand,
C     if GENCAN is being used to solve the bound-constrained subproblems
C     in an Augmented Lagrangian framework, the Hessian matrix must be
C     the Hessian of the Augmented Lagrangian function.
C
C     IMPORTANT: This subroutine does not need to be coded if the user
C     prefers to approximate the Hessian-vector product by incremental 
C     quotients. In this case, it is enough to set the GENCAN input
C     argument htvtype equal to 1 and an internal GENCAN subroutine will
C     be used to compute the approximation. In fact, this is the default
C     GENCAN option. See the GENCAN and EASYGENCAN arguments descriptions
C     for details.
C
C     On Entry:
C
C     nind  integer
C           number of component of the Hessian-vector product that
C           must be computed,
C
C     ind   integer ind(nind)
C           the component that must be computed are ind(1)-th ... ind(nind)-th,
C
C     n     integer,
C           number of variables,
C
C     x     double precision x(n),
C           current point,
C
C     m     integer,
C           number of constraints (equalities plus inequalities),
C
C     lambda double precision lambdae(m),
C           current estimation of the Lagrange multipliers,
C
C     rho   double precision rho(m)
C           penalty parameters,
C
C     d     double precision d(n)
C           vector of the Hessian-vector product.
C
C     NOTE: arguments m, lambda and rho are useful when GENCAN is being used
C     for solving the box-constrained subproblems of an Augmented Lagrangian
C     framework. When GENCAN is being used stand-alone for solving a bound-
C     constrained problem, these arguments are dummy arguments.
C
C     On Return
C
C     hd    double precision g(n),
C           Hessian-vector product,
C
C     flag  integer
C           0 means "no errors",
C           1 means "some error occurs in the gradient evaluation".

      implicit none

C     SCALAR ARGUMENTS
c      integer flag,m,n,nind
      integer n

C     ARRAY ARGUMENTS
c      integer ind(nind)
c      double precision d(n),hd(n),lambda(m),rho(m),x(n)

c      flag = - 1

      end
 
C**************************************************************************

C     Last update of EASYGENCAN: February 18th, 2005.

