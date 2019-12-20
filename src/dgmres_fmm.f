ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       This is the end of the debugging routines and the beginning 
c       of the complex GMRES routines
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
c
        subroutine dgmres_fmm(ier,n,multa,ifdn,ifinout,
     1     xys,dxys,qwts,nedges,ncint,lns,rns,icl,icr,icsgnl,icsgnr,
     2     xsgnl,xsgnr,ncorner,xmatc,xmatcsub,y,eps,numit,
     1     x,niter,errs,ngmrec,w)
        implicit real *8 (a-h,o-z)
        dimension errs(1)
        real *8 x(1),y(1),w(1)

        real *8 xys(2,n),dxys(2,n),qwts(n)
        integer ncint,nedges
        integer lns(nedges),rns(nedges),icl(ncint),icr(ncint)
        integer icsgnl(ncint),icsgnr(ncint)
        integer ncorner
        real *8 xsgnl(ncint),xsgnr(ncint)
        real *8 xmatc(ncorner,ncorner,ncint)
        real *8 xmatcsub(ncorner,ncorner,ncint)
        
        
c
        external multa
c
c     This subroutine solves a complex linear system Ax=y by means
c     of GMRES algorithm. This is a memory management routine for 
c     cgmres1 which performs the actual work.
c
c     Input parameters:
c
c     n - the dimensionality of the linear system
c     y - the right hand side
c     a - the matrix of the system (or whatever other parameter)
c     eps - the required accuracy
c     multa - the user-defined matrix-vector multiplication subroutine
c
c     the calling sequence for multa must be
c
c     multa(x,y,n,ifdn,ifinout,xys,dxys,qwts,nedges,ncint,
c         lns,rns,icl,icr,icsgnl,icsgnr,xsgnl,xsgnr,ncorner,
c         xmatc)
c
c       in (1), a is a matrix the system or whatever other parameter,
c       par1, par2 are whatever other parameters, x is an input vector,
c       y is a product Ax, and n is the dimensionality of a, x, and y.
c
c     numit - the maximum number of iteration permitted
c     ngmrec - the maximum number of iteration 
c            which GMRES algorithm needs to be restarted
c
c     w - must be at least (ngmrec*2+4)*n complex *16 elements long
c
c     Output parameters:
c
c     ier - error return code
c        ier=0 normal execution of the subroutine
c        ier=4 means that the maximum number iterations numit
c           has been reached without achieving the required accuracy eps
c        ier=8 means that the errors failed to decrease before the maximum
c           number of iterations has been reached or the required accuracy
c           eps has benn reached. 
c
c     x - the solution of the system
c     niter - the number of iterations performed 
c     errs - the array of errors produced by the algorithm. 
c        errs(i)=||y-Ax_i||, where x_i is a solution obtained on the i-th
c        GMRES iteration.
c
        ie=1
        lie=n*ngmrec
c
        iae=ie+lie
        liae=n*ngmrec
c
        iz=iae+liae
        lz=n
c
        iw=iz+lz
        lw=n
c
        ixr=iw+lw
        lxr=n
c
        izr=ixr+lxr
        lzr=n
c
        call dgmres_fmm1(ier,n,multa,ifdn,ifinout,
     1     xys,dxys,qwts,nedges,ncint,lns,rns,icl,icr,
     2     icsgnl,icsgnr,xsgnl,xsgnr,ncorner,xmatc,xmatcsub,
     2     y,eps,numit,x,niter,errs,
     1     ngmrec,w(ie),w(iae),w(iz),w(iw),w(ixr),w(izr))

        return
        end
c
c
c
c
c
        subroutine dgmres_fmm1(ier,n,multa,ifdn,ifinout,xys,
     1     dxys,qwts,nedges,ncint,lns,rns,icl,icr,icsgnl,icsgnr,
     2     xsgnl,xsgnr,ncorner,xmatc,xmatcsub,y,eps,numit,
     1     x,niter,errs,ngmrec,e,ae,z,w,xr,zr)
        implicit real *8 (a-h,o-z)
c
        dimension errs(1)
        real *8 x(1),y(1),e(n,1),ae(n,1),
     1       z(1),w(1),xr(1),zr(1),d
        real *8 xys(2,n),dxys(2,n),qwts(n)
        integer ncint,nedges
        integer lns(nedges),rns(nedges),icl(ncint),icr(ncint)
        integer icsgnl(ncint),icsgnr(ncint)
        integer ncorner
        real *8 xsgnl(ncint),xsgnr(ncint)
        real *8 xmatc(ncorner,ncorner,ncint)
        real *8 xmatcsub(ncorner,ncorner,ncint)

        
c
        external multa
c
        ier=0
        niter=0
c
        do 1000 i=1,n
ccc     x(i)=0
c
c     set x=y for cgmrespp  
c
        x(i)=y(i)
 1000 continue
c
c
        dold=1d+100
c
        m=0
c
        do 6000 iter=1,numit
c     
c       ... restart GMRES
c
        if( m .eq. ngmrec ) then
           m=0
        endif
c
        m=m+1
        niter=niter+1
c
        if(mod(niter,10).eq.0) 
     1     print *, "cgmres iter=",niter,", err=",dold 
cccc        call prinf('in cgmres, m=*',m,1)
c
        nsteps=5
c
        ifrec=0
        if( m .eq. 1 ) ifrec=1
        if( mod(niter,nsteps) .eq. 1 ) ifrec=1
        if( ifrec .eq. 0 ) goto 1600
c
cccc        call prinf('in cgmres, recalculate the residual at iter=*'
cccc     1     ,iter,1)
c
c       ... this is the first iteration, calculate the residual z=y-Ax
c       also, recalculate the residual every 5 steps
c
        call multa(x,z,n,ifdn,ifinout,xys,dxys,qwts,nedges,ncint,
     1    lns,rns,icl,icr,icsgnl,icsgnr,xsgnl,xsgnr,ncorner,xmatc,
     2    xmatcsub)
c
        do 1200 i=1,n
        z(i)=y(i)-z(i)
 1200 continue
c
 1600 continue
c
c
c       ... compute the error
c
        d=0
        do 2000 i=1,n
        d=d+z(i)*z(i)
 2000 continue
c
        d=sqrt(d)
        errs(niter)=dble(d)
c
cccc        call prin2('in cgmres, d=*',d,1)
c
        if( abs(d) .lt. eps ) return
c
        if( abs(d) .gt. dold ) then
c
c       ... the errors stopped decreasing, abort
c
        niter=niter-1
        do 2200 i=1,n
        x(i)=xr(i)
 2200 continue 
c
        ier=8
        return
        endif
c
 2400 continue
c
        dold=d
c
c       ... compute the new direction w=Az
c
        call multa(z,w,n,ifdn,ifinout,xys,dxys,qwts,nedges,ncint,
     1    lns,rns,icl,icr,icsgnl,icsgnr,xsgnl,xsgnr,ncorner,xmatc,
     2    xmatcsub)
c
c       ... orthogonalize w to all preceeding ae
c
        do 3000 i=1,n
        zr(i)=z(i)
 3000 continue
c
        do 3600 j=1,m-1
        d=0
        do 3200 i=1,n
        d=d+ae(i,j)*w(i)
 3200 continue
        do 3400 i=1,n
        w(i)=w(i)-d*ae(i,j)
        zr(i)=zr(i)-d*e(i,j)
 3400 continue
 3600 continue
c       
c       ... normalize the current direction w
c
        d=0
        do 3800 i=1,n
        d=d+w(i)*w(i)
 3800 continue
        d=1/sqrt(d)
c
c     ... store e and ae 
c
        do 4000 i=1,n
        ae(i,m)=w(i)*d
        e(i,m)=zr(i)*d
 4000 continue
c
c
c       ... double orthogonalize the current ae to all preceeding ae
c
        do 4030 j=1,m-1
        d=0
        do 4010 i=1,n
        d=d+ae(i,j)*ae(i,m)
 4010 continue
        do 4020 i=1,n
        ae(i,m)=ae(i,m)-d*ae(i,j)
        e(i,m)=e(i,m)-d*e(i,j)
 4020 continue
 4030 continue
c
        d=0
        do 4040 i=1,n
        d=d+ae(i,m)*ae(i,m)
 4040 continue
        d=1/sqrt(d)
c
cccc        call prin2('d=*',d-1,1)
c
c     ... store e and ae 
c
        do 4050 i=1,n
        ae(i,m)=ae(i,m)*d
        e(i,m)=e(i,m)*d
 4050 continue
c
c
c       ... update the solution x and the residual z
c
        d=0
        do 4200 i=1,n
        d=d+ae(i,m)*z(i)
 4200 continue
c
        do 4400 i=1,n
        z(i)=z(i)-d*ae(i,m)
 4400 continue
c
        do 4600 i=1,n
        xr(i)=x(i)
        x(i)=x(i)+d*e(i,m)
 4600 continue
c
 6000 continue
c
c       ... the maximum number of iterations has been reached, abort
c
        ier=4
c
        return
        end
c
c
c
c
