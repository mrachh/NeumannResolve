c
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       This is the end of the debugging code, and the beginning of 
c       the code for the finding of all roots of a smooth function
c       on the interval [-1,1]
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c
        subroutine roots_m1p1_find(n,
     1      roots,errs,nn,coefs,w2)
        implicit real *8 (a-h,o-z)
        save
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
c       will approximate on the interval [-1,1] 
c  coefs - the coefficients of the Chebychev expansion of the function
c       fun_eval on the interval [-1,1] (n+1 of them things). Please note 
c       that this parameter is returned entirely for the user's edification    
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                   Output Parameters:
c
c  roots - the zeroes of the function fun_eval as found by the subroutine
c  errs - the values of fun_eval evaluated at the points roots
c  nn - the number of zeroes found by the subroutine
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
        ixs=idiag+ldiag
        lxs=n/2+10
c
        iws=ixs+lxs
        lws=n/2+10
c
        iu=iws+lws
        lu=(n+1)**2/2+10
c
        iv=iu+lu
        lv=(n+1)**2/2+10
c
        ivals=iv+lv
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
ccc        call prin2('in roots *',coefs,2*(n+1))

        call roots_m1p1_find0(n,w2,w2(ixs),w2(iws),
     1      w2(iu),w2(iv),w2(ivals),coefs,w2(iw),w2(ip),
     2      w2(iq),w2(isupdiag),zero_mach)
c
c       find the roots that are on the interval [-1,1]
c
        call prin2('w2=*',w2,2*n) 
        nn=0
        eps=zero_mach*100000
ccc        call prin2('eps=*',eps,1)
        eps = 1.0d-4

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
        call roots_m1p1_bubble2(roots,nn)
c  
        return
        end

c
c
c
c 
c
        subroutine roots_m1p1_find0(n,
     1      diag,xs,ws,u,v,vals,
     2      coefs,w,p,q,supdiag,zero_mach)
        implicit real *8 (a-h,o-z)
        save
        dimension xs(1),ws(1),u(1),v(1),vals(1),
     1      par1(1),par2(1)
        complex *16 q(1),p(1),diag(1),supdiag(1),
     1      coefs(1),d,w(1)
        data ifcalled/0/
c
c       construct the test polynomial
c
        itype=2
        call chebexps(itype,n+1,xs,u,v,ws)
c
        done=1
        pi=atan(done)*4        

c       ... decompose
c
ccc        call roots_m1p1_matvec(u,n+1,f,coefs)

ccc        call prin2('coefs = *',coefs,2*(n+1))
c
c       . . . make sure that the last coefficient is not 
c             equal to zero
c
        if(ifcalled .eq. 0) call roots_m1p1_mach_zero(zero_mach)
        ifcalled=1
        d=0
        do 1600 i=1,n+1
c
        d=d+abs(coefs(i))**2
 1600 continue
c
        dd=sqrt(abs(d))
        d=abs(coefs(n+1))/dd
c
        if(abs(d) .lt. zero_mach) coefs(n+1)=dd*zero_mach
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

cc        call prin2('diags=  *',diag,2*n)

        return
        end
c
c
c
c 
c
        subroutine roots_m1p1_bubble2(cs,n)
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
        subroutine roots_m1p1_matvec(a,n,x,y)
        implicit real *8 (a-h,o-z)
        dimension a(n,n)
        complex *16 cd,x(n),y(n)
c 
c        apply the matrix a to the vector x obtaining y
c 
        do 1400 i=1,n
        cd=0
        do 1200 j=1,n
        cd=cd+a(i,j)*x(j)
 1200 continue
        y(i)=cd
 1400 continue
c
        return
        end  
c 
c 
c 
c 
c  
        subroutine roots_m1p1_mach_zero(zero_mach)
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


