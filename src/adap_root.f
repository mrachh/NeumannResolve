ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
c
        subroutine adap_root_getroots(k,nn,coefs,iposits,w,
     1      ab,roots,nroots,epsrt)
        implicit real *8 (a-h,o-z)
        dimension iposits(1),w(1),roots(1),ab(2,1),
     1      recoef(1000),coefim(1000),rtemp(1000)
        complex *16 coefs(1),ctemp(1000)

        nroots = 0

        do 3000 ind=1,nn

        do 1500 i=1,k
        recoef(i) = real(coefs(iposits(ind)+i-1))
        coefim(i) = imag(coefs(iposits(ind)+i-1))
 1500 continue

ccc        call prin2('recoef = *',recoef,k)

        nord = k-1
        call adap_root_find(nord,
     1      rtemp,errs,nr,recoef,w,epsrt)

        nt = 0
            
        do 1800 i=1,nr

ccc            call prin2('root fnd = *',roots(nroots+1),1)
            r = rtemp(i)
            call chebexev(r,val1,recoef,k-1)
ccc            call prin2('val = *',val,1)
            r = rtemp(i)
            call chebexev(r,val2,coefim,k-1)
ccc            call prin2('val = *',val,1)
ccc            stop

            err = sqrt(val1*val1+val2*val2)
ccc            call prin2('err here= *',err,1)
ccc            call prin2('val1 = *',val1,1)
ccc            call prin2('val2 = *',val2,1)
ccc    
        if (sqrt(val1*val1+val2*val2) .lt. epsrt) then

            nt = nt + 1
            roots(nroots + nt) = r
cc            call prinf('found ! *',0,0)
ccc            call prin2('rtemp = *',rtemp,nr)
ccc            call prinf('nr = *',nr,1)

            err = sqrt(val1*val1+val2*val2)
ccc            call prin2('err = *',err,1)


        end if     

        
ccc        call prinf('nroots found = *',nn,1)
ccc        call prin2('roots = *',roots,nn)

 1800 continue

        nr = nt

        do 2000 i=1,nr
        roots(nroots+i) = (roots(nroots+i)+1)/2*
     1      (ab(2,ind)-ab(1,ind))+ab(1,ind)
 2000 continue

ccc        call prin2('ab = *',ab,2)
ccc        call prin2_long('roots = *',roots,nroots)
ccc        call prinf('here *',0,0)
        nroots = nroots + nr

 3000 continue

c
c           .   .   .   remove duplicates
c

        if (nroots .eq. 0) return

        nr = 2

        do 4000 i=2,nroots
    
        if (abs(roots(i)-roots(nr-1)) .gt. epsrt) then

            roots(nr) = roots(i)
            nr = nr + 1

        end if

 4000 continue

        nroots = nr - 1

        call prin2('roots = *',roots,nroots)

        return
        end
c
c
c
c
c
         subroutine adap_root_disc(ier,a,b,fun,ipar1,par1,par2,
     1    k,ab,lab,eps,t,u,v,whts,x,f,coefs,nn,lused,ifinit7,
     1    lencfs,iposits)
        implicit real *8 (a-h,o-z)
        dimension t(1),u(k,k),v(k,k),whts(1),x(1),
     1      par1(1),par2(1),ab(2,1),iposits(1)
        complex *16 coefs(1),f(1)
        equivalence (rea,intint)
        external fun
c
c       this subroutine finds a subdivision of the interval [a,b],
c       such that on every subinterval, the user-provided function
c       fun is integrated to precision eps by a k-point Chebyshev 
c       quadrature. 
c
c                          input parameters:
c
c  a,b - the ends of the interval to be subdivided
c  fun - the user-supplied function to be integrated
c      The calling sequence of funget must be:
c
c      funget(x,par1,par2,f),
c  
c      where x is the input argument where the function is to be 
c      evaluated, f  (output) is the function, and par1,par2 are 
c      both input parameters, carrying whatever information funget 
c      might need.
c
c  par1,par2 - whatever the user-supplied subroutine funget needs
c  k - the number of Chebyshev nodes on each subinterval
c  lab - the length of the user-provided array ab (see below)
c  eps - the accuracy to which the chebyshev expansion must represent 
c     the user-supplied function for the subroutine to decide that 
c     the representation is accurate
c  t - the Chebyshev nodes on the interval [-1,1]. This is only an 
c     input parameter is ifinit7 (see below) has been set to 0;
c     otherwise, it is an output parameter.
c  u - the 2*k * 2*k  matrix converting the  values at of a polynomial 
c     of order 2*k-1 at 2*k chebyshev nodes into the coefficients of its 
c     chebyshev expansion. This is only an input parameter is ifinit7
c     (see below) has been set to 0; otherwise, it is an output parameter.
c  v - the 2*k * 2*k matrix converting the coefficients
c     of an 2*k-term chebyshev expansion into its values at 2*k chebyshev 
c     nodes (note that v is the inverse of u). This is only an input
c     parameter is ifinit7 (see below) has been set to 0; otherwise, it
c     is an output parameter.
c  whts - the Chebyshev quadrature weights (2*k of them).  This is only 
c     an input parameter is ifinit7 (see below) has been set to 0; 
c     otherwise, its an output parameter.
c  ifinit7 - the index telling the subroutine whether it should construct
c     the parameters t,u,v,whts (by calling the subroutine legeexps),
c     or view them as input parameters. In the latter case, these 
c     must have been created by a prior call to this subroutine (or to 
c     subroutine legeexps, for an adventurous and creative user).
c     ifinit7=1 means that the parameters t,u,v,whts will be created. 
c     ifinit7=0 means that the parameters t,u,v,whts will be viewed 
c     as input parameters. 
c
c                     output parameters:
c
c  ier - error return code. 
c         ier=0 means successful execution
c         ier=4 means that the space allocated in the array ab is
c         insufficient
c  ab - the array of ends of subintervals created. For each i=1,2,...,nn,
c     ab(1,i) is the left end of the i-th subinterval; ab(2,i) is
c     the right end of the i-th subinterval
c  t,u,v,whts - are only output is the parameter ifinit (see above)
c     has been set by the user to 1.
c  x - the Chebyshev nodes on the interval [a,b]
c  f - the values of the user-specified function fun at
c     the nodes x
c  coefs - the coefficients of the Chebyshev expansion of the 
c     function fun on the interval [a,b]
c  nn - the number of subintervals created on the interval [a,b];
c     also the number of entries in the array ab (see above)
c  lused - the total amount of space used by the subroutine in the
c     array ab; note that normally, it  is somewhat greater than 
c     the expected nn*2; The difference is not large, but should
c     be taken into account when allocating memory.
c
c       . . . start the recursive subdivision process
c

        ncused = 1
        nints  = 0

        ab(1,1)=a
        ab(2,1)=b
        nn=0
        ifinit=ifinit7
        istore=1
        ier=0
c
c       recursively subdivide the intervals till none need 
c       to be subdivided
c
        lused=0
        do 2000 i=1,1 000 000
c
        nn=nn+1
c
c       check if any more intervals need to be subdivided
c
        if(nn .gt. istore) goto 2200        
c
c       check if this interval needs to be subdivided
c
        d=ab(2,nn)-ab(1,nn)

c
        call adap_root_ifsplit(ier,ab(1,nn),ab(2,nn),k,fun,
     1      ipar1,par1,par2,eps/d,t,u,v,whts,x(ncused),ifinit,
     1      f(ncused),coefs(ncused))

ccc        call prinf('ier = *',ier,1)
c
        ifinit=0
c
        if(ier .eq. 0) then
        
            nints = nints + 1
            iposits(nints) = ncused
            ncused = ncused + k

            if (ncused .gt. lencfs) then
                call prinf('insuff. work array *',0,0)
                return
            end if        

            goto 2000

        end if
c
c       this interval needs to be subdivided. act accordingly
c
        if(lused .lt. istore*2+4) lused=istore*2+4
        if(istore*2 +4 .lt. lab) goto 1800
        ier=4
        return
 1800 continue
c
        istore=istore+1
        ab(1,istore)=ab(1,nn)
        ab(2,istore)=(ab(1,nn)+ab(2,nn))/2
c
        istore=istore+1
        ab(1,istore)=ab(2,istore-1)
        ab(2,istore)=ab(2,nn)
c
        intint=9654236
        ab(1,nn)=rea
        ab(2,nn)=1.0d35
 2000 continue
 2200 continue
c
        nn=nn-1
c
c       scan the intervals we have created, discarding those that have kids
c
        ii=0
        do 2400 i=1,nn
        rea=ab(1,i)
c       
        if( (ab(2,i) .gt. 1.0d34 ) .and. 
     1      (intint .eq. 9654236) ) goto 2400
        ii=ii+1
        ab(1,ii)=ab(1,i)
        ab(2,ii)=ab(2,i)
 2400 continue
        nn=ii
c
c        bubble-sort the array ab
c
        call adap_root_bubble(ab,iposits,nn)

ccc        call prinf('iposits = *',iposits,nints)

        return
        end  

c
c
c
c
c
        subroutine adap_root_ifsplit(ier,a,b,k,fun,ipar1,par1,
     1    par2,eps,t,u,v,whts,x,ifinit,f,coefs)
        implicit real *8 (a-h,o-z)
        dimension t(1),u(k,k),v(k,k),whts(1),x(1),
     1      par1(1),par2(1)
        complex *16 f(1),coefs(1)
c
c       this subroutine determines whether a k-point chebyshev
c       expansion of the user-supplied function f on the interval
c       [a,b] represents it with accuracy eps.
c
c                          input parameters:
c
c  a,b - the ends of the interval on which the approximation is checked
c  k - the order of the approximation for which the accuracy is checked
c  funget - the user-supplied function for which the accuracy is checked
c      The calling sequence of funget must be:
c
c      fun(x,i,par1,par2,f),
c  
c      where x is the input argument where the function is to be 
c      evaluated, f  (output) is the function, and par1,par2 are 
c      both input parameters, carrying whatever information funget 
c      might need.
c
c  par1,par2 - whatever the user-supplied subroutine funget needs
c  eps - the accuracy to which the chebyshev expansion must represent 
c     the user-supplied function for the subroutine to decide that 
c     the representation is accurate
c  t - the Chebyshev nodes on the interval [-1,1]. This is only an 
c     input parameter is ifinit (see below) has been set to 0;
c     otherwise, it is an output parameter.
c  u - the 2*k * 2*k  matrix converting the  values at of a polynomial 
c     of order 2*k-1 at 2*k chebyshev nodes into the coefficients of its 
c     chebyshev expansion. This is only an input parameter is ifinit 
c     (see below) has been set to 0; otherwise, it is an output parameter.
c  v - the 2*k * 2*k matrix converting the coefficients
c     of an 2*k-term chebyshev expansion into its values at 2*k chebyshev 
c     nodes (note that v is the inverse of u). This is only an input
c     parameter is ifinit (see below) has been set to 0; otherwise, it
c     is an output parameter.
c  whts - the Chebyshev quadrature weights (2*k of them).  This is only 
c     an input parameter is ifinit (see below) has been set to 0; 
c     otherwise, its an output parameter.
c  ifinit - the index telling the subroutine whether it sholud construct
c     the parameters t,u,v,whts (by calling the subroutine legeexps),
c     or view them as input parameters. In the latter case, these 
c     must have been created by a prior call to this subroutine (or to 
c     subroutine legeexps, for an adventurous and creative user).
c     ifinit=1 means that the parameters t,u,v,whts will be created. 
c     ifinit=0 means that the parameters t,u,v,whts will be viewed 
c     as input parameters. 
c
c                     output parameters:
c
c  t,u,v,whts - are only output is the parameter ifinit (see above)
c     has been set by the user to 1.
c  x - the Chebyshev nodes on the interval [a,b]
c  f - the values of the user-specified function funget at
c     the nodes x
c  coefs - the coefficients of the Chebyshev expansion of the 
c     function f on the interval [a,b]
c
c       . . . if the need be, construct the gaussian nodes on the 
c             interval [-1,1], and the matrix conevrting the values
c             of the function at these nodes into coefficients of 
c             its chebyshev expansion
c
        if(ifinit .eq. 0) goto 1200
c      
        itype=2
        call chebexps(itype,k,t,u,v,whts)
c
 1200 continue
c
c        discretize the interval and evaluate the
c        user-supplied function at its nodes
c     
        alpha=(b-a)/2
        beta=(b+a)/2
c
        do 1400 i=1,k
        x(i)=alpha*t(i)+beta       
cc        call prin2('x = *',x(i),1)
cccc        call fun(x(i),ipar1,par1,par2,par3,f(i) )
cc        call prin2('par1=*',par1,24)
cc        call prin2('par2=*',par2,13)
cc        call prinf('ipar1=*',ipar1,1)
        call fun(x(i),ipar1,par1,par2,f(i) )

ccc        call prin2('f val = *',f(i),2)
 1400 continue
cccc        call prin2('in ifsplit, x=*',x,k*2)
cccc        call prin2('in ifsplit, f=*',f,k*2)
c
c       construct the chebyshev expansion of the function on
c       the interval
c

        call adap_root_qrma2(u,f,coefs,k,1)

        call prin2('coefs=  *',coefs,2*k)
        call prin2('f = *',f,2*k)

ccc        call prinf('k = *',k,1)

c


cccc        call prin2('whts=*',whts,k*2)

        ddth=0
        do 1600 i=1,k
c
        ddth=ddth+abs(f(i))**2*whts(i)
 1600 continue
c
        ddth=sqrt(ddth)  
c
c       scan the coefficients and see if the last k 
c       of them are small
c
        ier=0
        do 1800 i=k-2,k
        d=abs(coefs(i))
c
        if(ddth .eq. 0) goto 1800

        if(d .gt. ddth*eps) then
            ier=4
            goto 2200
        endif
c
 1800 continue
c
 2200 continue
c
        if(ier .eq. 4) call prin2('ddth=*',ddth,1)
        if(ier .eq. 4) call prin2('and d=*',d,1)
c
        return
        end
c
c
c
c
c
        subroutine adap_root_qrma2(a,x,y,n,n0)
        implicit real *8 (a-h,o-z)
        dimension a(n,n)
        complex *16 d,x(1),y(1)
c
        do 2000 i=n0,n
        d=0
        do 1400 j=1,n
        d=d+a(i,j)*x(j)
 1400 continue
        y(i)=d
 2000 continue
        return
        end
c
c
c
c
c
        subroutine adap_root_bubble(ab,ipts,nn)
        implicit real *8 (a-h,o-z)
        dimension ab(2,1),ipts(1)
c
c       sort the array ab with respect to the first coordinate
c
        do 2000 i=1,nn
        do 1200 j=1,i-1
        if(ab(1,i) .ge. ab(1,j) ) goto 1200

        d=ab(1,i)
        ab(1,i)=ab(1,j)
        ab(1,j)=d
c
        d=ab(2,i)
        ab(2,i)=ab(2,j)
        ab(2,j)=d

        ival = ipts(i)
        jval = ipts(j)
        ipts(i) = jval
        ipts(j) = ival

 1200 continue
 2000 continue
        return
        end
c
c
c
c
c
c
        subroutine adap_root_pltcmplx(xs,cvs,nvs,iw)
        implicit real *8 (a-h,o-z)
        dimension xs(1),yre(100 000),yim(100 000)
        complex *16 cvs(1)


        do 2000 i=1,nvs
        yre(i) = real(cvs(i))
        yim(i) = imag(cvs(i))
 2000 continue        

        itype1 = 2
        itype2 = 2
        call quagraph2(iw,xs,yre,nvs,itype1,
     1      xs,yim,nvs,itype2,'*')

        return
        end
c
c
c
c
c
        subroutine adap_root_find(n,
     1      roots,errs,nn,coefs,w2,epsuser)
        implicit real *8 (a-h,o-z)
        save
        external fun_eval
        dimension roots(1),errs(1),coefs(1),par1(1),par2(1)
        complex *16 ima,w2(1)
        data ima/(0.0d0,1.0d0)/
c
c       This subroutine finds the roots of the user-specified 
c       function fun_eval on the interval [-1,1]. It is assumed 
c       that this function is well approximated by the Chebychev
c       expansion of order n, where n is also specified by the 
c       user. 
c
c                    Input parameters:
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc   
c
c  n - the order of the Chebychev polynomial with which the subroutine 
c       will approximate on the interval [-1,1] the user-supplied 
c       function fun_eval
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                   Output Parameters:
c
c  roots - the zeroes of the function fun_eval as found by the subroutine
c  errs - the values of fun_eval evaluated at the points roots
c  nn - the number of zeroes found by the subroutine
c  coefs - the coefficients of the Chebychev expansion of the function
c       fun_eval on the interval [-1,1] (n+1 of them things). Please note 
c       that this parameter is returned entirely for the user's edification    
c
c                   Work arrays:
c
c  w2 - must be ??? real *8 locations long.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       . . . allocate memory
c
        idiag=1
        ldiag=n+4
c
c
        ivals=idiag + ldiag
        lvals=n/2+10
c
        iw=ivals+lvals
        lw=n*6+100
 
cc        lw=lw*2
c
        ip=iw+lw
        lp=n+10+n
c
        iq=ip+lp
        lq=n+10
c
        isupdiag=iq+lq
        lsupdiag=n+10
c
        ltot=isupdiag+lsupdiag

c
        call adap_root_find0(n,w2,
     1      coefs,w2(ivals),w2(iw),w2(ip),
     2      w2(iq),w2(isupdiag),zero_mach)
c
c       find the roots that are on the interval [-1,1]
c
        nn=0
        eps=zero_mach*10000
        eps=epsuser
ccc        call prin2('eps=*',eps,1)
        do 2000 i=1,n
c
        d=ima*w2(i)
        if(abs(d) .gt. eps) goto 2000
        d=w2(i)
        if(d .gt. 1+eps) goto 2000
        if(d .lt. -1-eps) goto 2000
c
        nn=nn+1
        roots(nn)=d
 2000 continue
c
        call adap_root_bubble2(roots,nn)
c  
        return
        end

c
c
c
c 
c
        subroutine adap_root_find0(n,
     1      diag,vals,
     2      coefs,w,p,q,supdiag,zero_mach)
        implicit real *8 (a-h,o-z)
        dimension w(1),vals(1),coefs(1),
     1      par1(1),par2(1)
        complex *16 q(1),p(1),diag(1),supdiag(1)
        data ifcalled/0/
        save
c
        done=1
        pi=atan(done)*4        
c
        do 1400 i=1,n+1
c
        coefs(i) = vals(i)
 1400 continue
c
        if(ifcalled .eq. 0) call adap_root_mach_zero(zero_mach)
        ifcalled=1
        d=0
        do 1600 i=1,n+1
c
        d=d+coefs(i)**2
 1600 continue
c
        dd=sqrt(d)
        d=abs(coefs(n+1))/dd
c
        if(d .lt. zero_mach) coefs(n+1)=dd*zero_mach
ccc        call prin2('and coefs=*',coefs,n)
c
c       construct the matrix
c
        done=1
c
        do 2200 i=1,n-1
c
        supdiag(i)=done/2
 2200 continue
c
        supdiag(1)=1
c
        do 2400 i=1,n
c
        w(i)=-coefs(i)/coefs(n+1)/2
        diag(i)=0
 2400 continue
c
c       scale the first row/column
c
        d=2
        d=sqrt(d)
c
        supdiag(1)=supdiag(1)/d
        w(1)=w(1)*d
c
c       scale the last column
c
        d=0
        do 2800 i=1,n-1
c
        d=d+w(i)**2
 2800 continue
c
        d=1/sqrt(sqrt(d))
c
        do 3200 i=1,n-1
c
        w(i)=w(i)*d
 3200 continue
c
        supdiag(n-1)=supdiag(n-1)/d
c
c       evaluate the roots via the sparse subroutine
c
        diag(n)=w(n)
c
        do 4600 i=1,n
c
        p(i)=w(i)
        q(i)=0
 4600 continue
c
        q(n)=1
        p(n)=0
        p(n-1)=p(n-1)-supdiag(n-1)
c

        call trid_p_rank1(ier,diag,supdiag,q,p,n,
     1      w,niter_max,aver_niter)

ccc        call prin2('diag=*',diag,2*n)
cc        call adap_root_trid_p_rank1(ier,diag,supdiag,q,p,n,
cc     1      w,niter_max,aver_niter)
c
        return
        end
c
c
c
c 
c
        subroutine adap_root_bubble2(cs,n)
        implicit real *8 (a-h,o-z)
        real *8 cs(n)
c
        do 1400 i=1,n
c
        do 1200 j=1,n-1
c
        d=cs(j)
        dd=cs(j+1)
        if(d .lt. dd) goto 1200

        cc=cs(j)
        cs(j)=cs(j+1)
        cs(j+1)=cc
 1200 continue
 1400 continue
c
        return
        end
c
c 
c 
c 
c 
        subroutine adap_root_mach_zero(zero_mach)
        implicit real *8 (a-h,o-z)
c
        zero_mach=100       
c
        d1=1.1d0
        d=1.11d0
        do 1200 i=1,10000
c
        d=d/2
        d2=d1+d
        d4=sin(sin(d1)-sin(d2))
c
        if(d4 .eq. 0) goto 1400
c
 1200 continue
 1400 continue
c
        zero_mach=d
        return
        end


