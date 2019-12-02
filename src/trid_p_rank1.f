c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       This is the end of the debugging code, and the beginning of 
c       the diagonalizing code proper. This file contains one (1) 
c       user-callable subbroutine: trid_p_rank1. Its description 
c       can be found in the body of the subroutine.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
        subroutine trid_p_rank1(ier,diag,supdiag,p,q,n,
     1      w,niter_max,aver_niter)
        implicit real *8 (a-h,o-z)
        save
        complex *16 diag(1),supdiag(1),p(1),q(1),w(1)
        data ifcalled/0/
c
c       This subroutine finds the spectrum of the upper Hessenberg 
c       matrix B of the form 
c
c                    B = A + P \circ Q^*,                               (1)
c
c       where A is a tridiagonal Hermitian matrix, and P, Q are a 
c       pair of vectors. The claim to fame of this subroutine is 
c       the fact that it costs O(n^2) operations, and is numerically
c       stable.
c
c       PLEASE NOTE that this is a very specialized subroutine, in 
c       in that the matrix B is BOTH of the form (1) and upper 
c       Hessenberg, so that 
c
c                    B(i,j)=0                                           (2)
c
c       for all 
c
c                    j > i+1.                                           (3)
c
c
c                    Input parameters:
c
c  diag - the diagonal of the matrix A in (1);
c       destroyed by this subroutine utterly
c  supdiag - the superdiagonal of the matrix A in (1);
c       destroyed by this subroutine utterly
c  p - the vector P in (1); destroyed by this subroutine utterly
c  q - the vector Q in (1); destroyed by this subroutine utterly
c  n - the dimensionality of all matrices and vectors involved
c
c
c                    Output parameters:
c
c  ier - error return code:
c    ier=0 means successful execution
c    ier=128 means that at some point the QR process failed to comverge;
c      this should not happen, and most likely indicates a user bug.
c  diag - the eigenvaalues of the matrix B in (2), in no particular
c      order
c  niter_max - the maximum number of iterations taken by the QR process
c      on any iteration; provided solely for a curious user's edification
c  aver_niter - the average number of iterations taken by the QR process;
c      provided solely for a curious user's edification
c       
c                    Work array:
c
c  w - must be at least 8*n+100 real *8 locations long
c
c       . . . find the machine zero
c
        if(ifcalled .eq. 0) call trid_p_rank1_mach_zero(zero_mach)
        ifcalled=1
c
        dm=p(n)*conjg(p(n))+q(n)*conjg(q(n))+diag(n)*conjg(diag(n))
        do 1200 i=1,n-1
c
        dm=dm+p(i)*conjg(p(i))+q(i)*conjg(q(i))+
     1      diag(i)*conjg(diag(i))+supdiag(i)*conjg(supdiag(i))
 1200 continue
c
        eps=sqrt(dm)*zero_mach
c
c       . . . find the first n-2 eigenvalues
c
        ier=0
        niter_max=0
c
        aver_niter=0
c
c       augment the user-supplied tridiagonal matrix
c
        do 1300 i=1,n
c
        diag(i)=diag(i)+p(i)*conjg(q(i))
 1300 continue
c
        do 1400 i=1,n-1
c
        w(i)=conjg(supdiag(i))+p(i+1)*conjg(q(i))
        supdiag(i)=supdiag(i)+p(i)*conjg(q(i+1))
 1400 continue
c
        iw=n+3
c
        do 2000 i=1,n-1
c
        nn=n-i+1
        call trid_p_rank1_one_dimen(jer,diag(i),w(i),supdiag(i),
     1      p(i),q(i),nn,eps,w(iw),niter)
c
        aver_niter=aver_niter+niter
c
        if(jer .ne. 0) then
            ier=128
            return
        endif
c
        if(niter_max .lt. niter) niter_max=niter
ccc        call prinf('niter=*',niter,1)
ccc        call prin2('supdiag(i)=*',supdiag(i),n*2)
 2000 continue
c
        aver_niter=aver_niter/n
c
        return
        end
c
c
c
c 
c
        subroutine trid_p_rank1_one_dimen(ier,diag,subdiag,supdiag,
     1      p,q,n,eps,w,niter)
        implicit real *8 (a-h,o-z)
        complex *16 diag(1),subdiag(1),supdiag(1),p(1),q(1),
     1      w(1),aa,bb,cc,clam,clam1,clam2,clamsum
c
c       conduct QR iterations
c
        ier=16
        ifout=0
        clamsum=0
        do 2000 ijk=1,1000
c
c       . . . find the shift
c
        aa=1
        bb=-(diag(1)+diag(2))
        cc=diag(1)*diag(2) - supdiag(1)*subdiag(1)
c
        clam1=(-bb+sqrt(bb**2-4*aa*cc)) / (2*aa)
        clam2=(-bb-sqrt(bb**2-4*aa*cc)) / (2*aa)
c
        clam=clam1
        if( abs(clam2-diag(1)) .lt. abs(clam1-diag(1)) ) clam=clam2
c
        clamsum=clamsum+clam
c
        do 1400 i=1,n
c
        diag(i)=diag(i)-clam
 1400 continue
c
c       conduct a QR ireration
c
        iw=1
        lw=n+4
c
        iuus=iw+lw
        luus=n*4+16
c
        ltot=iuus+luus
c
        call trid_p_rank1_iter(n,p,q,diag,subdiag,supdiag,w(iw),w(iuus))
c
c       check if the time has come to terminate the iterations
c
        d=abs(supdiag(1))
        if(d .lt. eps) ifout=ifout+1
        if(ifout .eq. 2) goto 2400
c
 2000 continue
c
        return
c
 2400 continue
c
        ier=0
        niter=ijk
c
        do 2600 i=1,n
        diag(i)=diag(i)+clamsum
 2600 continue
c
        return
        end
c
c
c
c 
c
        subroutine trid_p_rank1_iter(n,p,q,diag,subdiag,supdiag,w,uus)
        implicit real *8 (a-h,o-z)
        complex *16 xy(2),uu(2,2),uus(2,2,1),p(1),q(1),
     1      diag(1),subdiag(1),supdiag(1),w(1)
c
        ii=0
c
        zero=0
        do 1200 i=1,n
c
        w(i)=zero
 1200 continue
c       
        do 2000 k=n,2,-1
c
        ii=ii+1
        xy(1)=supdiag(k-1)
        xy(2)=diag(k)
c
        call trid_p_rank1_rotfnd4(xy,uus(1,1,ii))
c
        if(k .ne. 2)
     1      w(k-2)=p(k)*conjg(q(k-2))-q(k)*conjg(p(k-2))
c
        call trid_p_rank1_two_elems_rotate01(uus(1,1,ii),
     1      supdiag(k-1),diag(k))
c
        call trid_p_rank1_two_elems_rotate(uus(1,1,ii),
     1      diag(k-1),subdiag(k-1))
        call trid_p_rank1_two_elems_rotate(uus(1,1,ii),
     1      subdiag(k-2),w(k-2))
c
        call trid_p_rank1_two_elems_rotate(uus(1,1,ii),q(k-1),q(k))
        call trid_p_rank1_two_elems_rotate(uus(1,1,ii),p(k-1),p(k))
c
 2000 continue
c
c       restore the matrix to its tridiagonal form
c
        iii=0
        do 3000 k=n,2,-1
c
        iii=iii+1
c
        uu(2,1)=conjg(uus(2,1,iii))
        uu(2,2)=conjg(uus(2,2,iii))
        uu(1,1)=conjg(uus(1,1,iii))
        uu(1,2)=conjg(uus(1,2,iii))
c
        call trid_p_rank1_two_elems_rotate(uu,subdiag(k-1),diag(k))
        call trid_p_rank1_two_elems_rotate10(uu,diag(k-1),supdiag(k-1))
c
        if(k .le. n-1) 
     1      call trid_p_rank1_two_elems_rotate01(uu,w(k-1),subdiag(k))
 3000 continue
c
        return
        end
c
c
c
c 
c
        subroutine trid_p_rank1_iter_readable
     1      (n,p,q,diag,subdiag,supdiag,w,uus)
        implicit real *8 (a-h,o-z)
        complex *16 xy(2),uu(2,2),uus(2,2,1),p(1),q(1),
     1      diag(1),subdiag(1),supdiag(1),w(1)
c
        ii=0
c
        do 1200 i=1,n
c
        w(i)=0
 1200 continue
c       
        do 2000 k=n,2,-1
c
        xy(1)=supdiag(k-1)
        xy(2)=diag(k)
c
        call trid_p_rank1_rotfnd4(xy,uu)
c
        if(k .ne. 2)
     1      w(k-2)=p(k)*conjg(q(k-2))-q(k)*conjg(p(k-2))
c
        call trid_p_rank1_two_elems_rotate(uu,supdiag(k-1),diag(k))
        call trid_p_rank1_two_elems_rotate(uu,diag(k-1),subdiag(k-1))
        call trid_p_rank1_two_elems_rotate(uu,subdiag(k-2),w(k-2))
c
        call trid_p_rank1_two_elems_rotate(uu,q(k-1),q(k))
        call trid_p_rank1_two_elems_rotate(uu,p(k-1),p(k))
c
        ii=ii+1
        call trid_p_rank1_rarrcopy(uu,uus(1,1,ii),8)
 2000 continue
c
c       restore the matrix to its tridiagonal form
c
        iii=0
        do 3000 k=n,2,-1
c
        iii=iii+1
        call trid_p_rank1_rarrcopy(uus(1,1,iii),uu,8)
c
        uu(2,1)=conjg(uu(2,1))
        uu(2,2)=conjg(uu(2,2))
        uu(1,1)=conjg(uu(1,1))
        uu(1,2)=conjg(uu(1,2))
c
        call trid_p_rank1_two_elems_rotate(uu,subdiag(k-1),diag(k))
        call trid_p_rank1_two_elems_rotate(uu,diag(k-1),supdiag(k-1))
c
        if(k .le. n-1) 
     1      call trid_p_rank1_two_elems_rotate(uu,w(k-1),subdiag(k))
 3000 continue
c
        return
        end
c 
c 
c 
c 
c 
        subroutine trid_p_rank1_two_elems_rotate(u,x,y)
        implicit complex *16 (a-h,o-z)
        dimension u(2,2)
c
        d1=u(1,1)*x+u(1,2)*y
        y=u(2,1)*x+u(2,2)*y
c 
        x=d1
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine trid_p_rank1_two_elems_rotate01(u,x,y)
        implicit complex *16 (a-h,o-z)
        dimension u(2,2)
        data zero/0.0d0/
c
        y=u(2,1)*x+u(2,2)*y
        x=zero
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine trid_p_rank1_two_elems_rotate10(u,x,y)
        implicit complex *16 (a-h,o-z)
        dimension u(2)
c
        y=u(2)*x
        x=u(1)*x
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine trid_p_rank1_rotfnd4(a,u)
        implicit complex *16 (a-h,o-z)
        dimension a(2),u(2,2)
        real *8 d,d2,done
        data done/1.0d0/
c 
        alpha=a(1)
        beta=a(2)
c 
        d=(a(1)*conjg(a(1))+a(2)*conjg(a(2)))
        d=sqrt(d)
c 
        if(d .eq. 0) then
c 
            u(2,2)=1
            u(1,2)=0
            u(1,1)=1
            u(2,1)=0
            return
        endif
c
        d2=done/d
        u(1,1)=beta*d2
        u(2,2)=conjg(u(1,1))
        u(1,2)=-alpha*d2
        u(2,1)=-conjg(u(1,2))
c
        return
        end
c 
c 
c 
c 
c 
        subroutine trid_p_rank1_rarrcopy(x,y,n)
        implicit real *8 (a-h,o-z)
        real *8 x(1),y(1)
c 
        do 2400 i=1,n
        y(i)=x(i)
 2400 continue
        return
        end  
c 
c 
c 
c 
c  
        subroutine trid_p_rank1_mach_zero(zero_mach)
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
