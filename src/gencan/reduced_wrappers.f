      subroutine calcf(nind,ind,x,n,xc,m,lambda,rho,f,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer nind,n,m,inform
      double precision f

C     ARRAY ARGUMENTS
      integer ind(nind)
      double precision x(n),xc(n),lambda(m),rho(m)

C     This subroutines computes the objective function. 
C
C     It is called from the reduced space (dimension nind), expands the
C     point x where the function will be evaluated and call the 
C     subroutine evalf to compute the objective function Finally, 
C     shrinks vector x to the reduced space. 
C
C     About subroutines named calc[something]. The subroutines whos 
C     names start with ``calc'' work in (are called from) the reduced 
C     space. Their tasks are (i) expand the arguments to the full space, 
C     (ii) call the corresponding ``eval'' subroutine (which works in 
C     the full space), and (iii) shrink the parameters again and also 
C     shrink a possible output of the ``eval'' subroutine. Subroutines
C     of this type are: calcf, calcg, calchd, calcgdiff and calchddiff. 
C     The corresponding subroutines in the full space are the user 
C     defined subroutines evalf, evalg and evalhd.

C     LOCAL SCALARS
      integer i

C     Complete x

      do i = nind + 1,n
          x(i) = xc(i)
      end do

C     Expand x to the full space

      call expand(nind,ind,n,x)

C     Compute f calling the user supplied subroutine evalf

      call evalal(n,x,m,lambda,rho,f,inform)

C     Shrink x to the reduced space

      call shrink(nind,ind,n,x)

      return

      end

C     ******************************************************************
C     ******************************************************************

      subroutine calcg(nind,ind,x,n,xc,m,lambda,rho,g,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer nind,n,m,inform

C     ARRAY ARGUMENTS
      integer ind(nind)
      double precision x(n),xc(n),lambda(m),rho(m),g(n)

C     This subroutine computes the gradient vector g of the objective 
C     function. 
C
C     It is called from the reduced space (dimension nind), expands the
C     point x where the gradient will be evaluated and calls the user 
C     supplied subroutine evalg to compute the gradient vector. Finally, 
C     shrinks vectors x and g to the reduced space. 
C
C     About subroutines named calc[something]. The subroutines whos 
C     names start with ``calc'' work in (are called from) the reduced 
C     space. Their tasks are (i) expand the arguments to the full space, 
C     (ii) call the corresponding ``eval'' subroutine (which works in 
C     the full space), and (iii) shrink the parameters again and also 
C     shrink a possible output of the ``eval'' subroutine. Subroutines
C     of this type are: calcf, calcg, calchd, calcgdiff and calchddiff. 
C     The corresponding subroutines in the full space are the user 
C     defined subroutines evalf, evalg and evalhd.

C     LOCAL SCALARS
      integer i

C     Complete x

      do i = nind + 1,n
          x(i) = xc(i)
      end do

C     Expand x to the full space

      call expand(nind,ind,n,x)

C     Compute the gradient vector calling the user supplied subroutine 
C     evalg

      call evalnal(n,x,m,lambda,rho,g,inform)

C     Shrink x and g to the reduced space

      call shrink(nind,ind,n,x)
      call shrink(nind,ind,n,g)

      return

      end

C     ******************************************************************
C     ******************************************************************

      subroutine calcgdiff(nind,ind,x,n,xc,m,lambda,rho,g,sterel,steabs,
     +inform)

      implicit none

C     SCALAR ARGUMENTS
      integer nind,n,m,inform
      double precision sterel,steabs

C     ARRAY ARGUMENTS
      integer ind(nind)
      double precision x(n),xc(n),lambda(m),rho(m),g(n)

C     This subroutine approximates the gradient vector g of the 
C     objective function in the reduced space using central finite 
C     differences.
C
C     It is called from the reduced space (dimension nind), expands the
C     point x where the gradient will be estimated and calls evalf
C     subroutine (to evaluate the objective function) 2 * nind times.
C     Finally, shrinks vectors x and g to the reduced space. 
C
C     About subroutines named calc[something]. The subroutines whos 
C     names start with ``calc'' work in (are called from) the reduced 
C     space. Their tasks are (i) expand the arguments to the full space, 
C     (ii) call the corresponding ``eval'' subroutine (which works in 
C     the full space), and (iii) shrink the parameters again and also 
C     shrink a possible output of the ``eval'' subroutine. Subroutines
C     of this type are: calcf, calcg, calchd, calcgdiff and calchddiff. 
C     The corresponding subroutines in the full space are the user 
C     defined subroutines evalf, evalg and evalhd.

C     LOCAL SCALARS
      integer i,indi
      double precision fminus,fplus,step,tmp

C     Complete x

      do i = nind + 1,n
          x(i) = xc(i)
      end do

C     Expand x to the full space

      call expand(nind,ind,n,x)
      
C     Approximate the gradient vector by central finite differences

      do i = 1,nind
          indi = ind(i)

          tmp  = x(indi)
          step = max( steabs, sterel * abs( tmp ) )

          x(indi) = tmp + step
          call evalal(n,x,m,lambda,rho,fplus,inform)
          if ( inform .lt. 0 ) then
              return
          end if

          x(indi) = tmp - step
          call evalal(n,x,m,lambda,rho,fminus,inform)
          if ( inform .lt. 0 ) then
              return
          end if

          g(indi) = ( fplus - fminus ) / ( 2.0d0 * step )
          x(indi) = tmp
      end do

C     Shrink x and g to the reduced space

      call shrink(nind,ind,n,x)
      call shrink(nind,ind,n,g)

      return

      end


C     ******************************************************************
C     ******************************************************************

      subroutine calchd(nind,ind,x,d,g,n,xc,m,lambda,rho,hd,xtmp,sterel,
     +steabs,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,m,n,nind
      double precision steabs,sterel

C     ARRAY ARGUMENTS
      integer ind(nind)
      double precision d(n),g(n),hd(n),lambda(m),rho(m),x(n),xc(n),
     +        xtmp(n)        

C     This subroutine computes the product Hessian times vector d. As it
C     is called from the reduced space, it expands vectors x and d,  
C     calls the user supplied subroutine evalhd to compute the Hessian 
C     times vector d product, and shrinks vectors x, d and hd. 
C
C     About subroutines named calc[something]. The subroutines whos 
C     names start with ``calc'' work in (are called from) the reduced 
C     space. Their tasks are (i) expand the arguments to the full space, 
C     (ii) call the corresponding ``eval'' subroutine (which works in 
C     the full space), and (iii) shrink the parameters again and also 
C     shrink a possible output of the ``eval'' subroutine. Subroutines
C     of this type are: calcf, calcg, calchd, calcgdiff and calchddiff. 
C     The corresponding subroutines in the full space are the user 
C     defined subroutines evalf, evalg and evalhd.

C     LOCAL SCALARS
      integer i

C     Complete d with zeroes

      do i = nind + 1,n
          d(i) = 0.0d0
      end do

C     Complete x

      do i = nind + 1,n
          x(i) = xc(i)
      end do

C     Expand x and d to the full space

      call expand(nind,ind,n,x)
      call expand(nind,ind,n,d)
      call expand(nind,ind,n,g)

C     Compute the Hessian times vector d product calling the user 
C     supplied subroutine evalhd

      call evalhd(n)

C     Shrink x, d and hd to the reduced space

      call shrink(nind,ind,n,x)
      call shrink(nind,ind,n,d)
      call shrink(nind,ind,n,g)
      call shrink(nind,ind,n,hd)
      
      end

C     ******************************************************************
C     ******************************************************************

      subroutine calchddiff(nind,ind,x,d,g,n,xc,m,lambda,rho,gtype,hd,
     +xtmp,sterel,steabs,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer gtype,inform,m,n,nind
      double precision steabs,sterel

C     ARRAY ARGUMENTS
      integer ind(nind)
      double precision d(n),g(n),hd(n),lambda(m),rho(m),x(n),xc(n),
     +        xtmp(n)

C     This subroutine computes the Hessian times vector d product by
C     means of a ``directional finite difference''. The idea is that, at 
C     the current point x, the product H d is the limit of 
C 
C     [ Gradient(x + t d) - Gradient(x) ] / t
C
C     In this implementation we use
C
C     t = max(steabs, sterel ||x||_\infty) / ||d||_\infty
C
C     provided that d is not equal 0, of course. 
C
C     So, we evaluate the Gradient at the auxiliary point x + t d and
C     use the quotient above to approximate H d. To compute the gradient 
C     vector at the auxiliary point it is used evalg or evalgdiff 
C     depending on gtype parameter.
C
C     About subroutines named calc[something]. The subroutines whos 
C     names start with ``calc'' work in (are called from) the reduced 
C     space. Their tasks are (i) expand the arguments to the full space, 
C     (ii) call the corresponding ``eval'' subroutine (which works in 
C     the full space), and (iii) shrink the parameters again and also 
C     shrink a possible output of the ``eval'' subroutine. Subroutines
C     of this type are: calcf, calcg, calchd, calcgdiff and calchddiff. 
C     The corresponding subroutines in the full space are the user 
C     defined subroutines evalf, evalg and evalhd.

C     On Entry
C
C     n        integer
C              order of the x
C
C     x        double precision x(n)
C              point for which Hessian(x) times d will be approximated
C
C     d        double precision d(n)
C              vector for which the Hessian times vetor product will
C              be approximated
C
C     g        double precision g(n)
C              gradient at x
C
C     xtmp     double precision xtmp(n)
C              working vector
C
C     sterel   double precision
C     steabs   double precision
C              these constants mean a ``relative small number'' and 
C              ``an absolute small number''
C
C     On Return
C
C     hd       double precision hd(n)
C              approximation of H d

C     LOCAL SCALARS
      integer flag,i,indi
      double precision dsupn,step,tmp,xsupn

      inform = 0

C     Compute incremental quotients step

      xsupn = 0.0d0
      dsupn = 0.0d0
      do i = 1,nind
          xsupn = max( xsupn, abs( x(i) ) )
          dsupn = max( dsupn, abs( d(i) ) )
      end do

c Safeguard added by LM
      if(dsupn.lt.1.d-20) dsupn = 1.d-20

      step = max( sterel * xsupn, steabs ) / dsupn 

C     Set the point at which the gradient will be evaluated

      do i = 1,nind
          xtmp(i) = x(i) + step * d(i)
      end do

C     Evaluate the gradient at xtmp = x + step * d

      if ( gtype .eq. 0 ) then

C         Complete xtmp

          do i = nind + 1,n
              xtmp(i) = xc(i)
          end do

C         Expand xtmp to the full space

          do i = nind,1,-1
              indi = ind(i)
              if ( i .ne. indi ) then
                  tmp        = xtmp(indi)
                  xtmp(indi) = xtmp(i)
                  xtmp(i)    = tmp
              end if
          end do

c         Compute the gradient at xtmp = x + step * d

          call evalnal(n,xtmp,m,lambda,rho,hd,flag)

C         Shrink hd to the reduced space

          do i= 1, nind
              indi= ind(i)
              if (i.ne.indi) then
                  tmp      = hd(indi)
                  hd(indi) = hd(i)
                  hd(i)    = tmp
              end if    
          end do

      else if ( gtype .eq. 1 ) then

          call calcgdiff(nind,ind,xtmp,n,xc,m,lambda,rho,hd,sterel,
     +    steabs,inform)

      end if

C     Compute incremental quotients

      do i = 1,nind
          hd(i) = ( hd(i) - g(i) ) / step
      end do 

      return

      end


C     ******************************************************************
C     ******************************************************************

