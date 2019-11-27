      implicit real *8 (a-h,o-z)
      real *8, allocatable :: verts(:,:)
      real *8 xyin(2),xyout(2),dpars(2)
      complex *16 zpars,zk

      character *100 fname

      external helm_slp

      call prini(6,13)

      nvmax = 100000
      allocate(verts(2,nvmax))

      igeom = 1

      call loadverts_demos(igeom,nverts,verts,xyin,xyout)

c
c
c     This subroutine solves the interior/exterior
c        Dirichlet/Neumann Helmholtz boundary value
c        problem using a double layer as the representation
c        for the dirichlet problem, and single layer as the
c        representation for the neumann problem.
c         (Beware of spurious resonances)
c
c     The right hand side is a user supplied subroutine 
c       fker whose calling sequence is of the form
c
c      fker(y,dpars,zpars,ipars,f)
c       
c       y contains the source info 
c        y(1:2) - location of source
c        y(3:4) - dydt for the source point (note not normalized)
c
c
c     The output is the potential at the user specified targets
c     and also stores the following information in a file
c     if requested
c      line 1: nverts - number of vertices
c      line 2: n - number of discretization points
c      line 3: nch - number of panels
c      line 4: ntarg - targets for plotting in interior/exterior
c          (auto generated and potential computed to 13 digits)
c      line 4-3+nverts- xy locations of vertices
c      line 4+nverts-3+n+nverts: 
c        y1,y2, dy1/dt,dy2/dt,real(rhs),imag(rhs),real(soln),imag(soln)
c      line 4+n+nverts-3+n+nverts+nch: ixyzs,ks,iscorn
c      line 4+n+nverts+nch -3+n+nverts+nch+ntarg:
c         targ1,targ2,real(pot),imag(pot)
c
c       ifdn = dirichlet/Neumann flag
c          ifdn = 0 => Dirichlet problem
c          ifdn = 1 => Neumann problem     
c 
c       ifinout = interior/exterior problem flag
c          ifinout = 0 => interior problem
c          ifinout = 1 => exterior problem
c
c       iffast = flag for solving linear system
c          iffast = 0 => direct inversion using blas
c          iffast = 1 => iterative solver using blas for matvec
c          iffast = 2 => iterative solver using fmm for matvec
c          
 
      ifdn = 1 
      ifinout = 0
      iffast = 0

      zk = 2.1d0

      ifwrite = 1
      fname = 'equi_interior_test.dat'

      call helm_solver(zk,nverts,verts,ifdn,ifinout,iffast,fker,
     1   dpars,zpars,ipars,ntarg,xyin,pot,ifwrite,fname)

      
      stop
      end
c
c
c
c
c 
      subroutine loadverts_demos(igeom,nverts,verts,xyin,xyout)
      implicit real *8 (a-h,o-z)
      real *8 verts(2,*),xyin(2),xyout(2)

      if(igeom.eq.1) then

        nverts = 3
        verts(1,1) = 0
        verts(2,1) = 0

        verts(1,2) = 0.5d0
        verts(2,2) = -sqrt(3.0d0)/2.0d0-0.1d0
        
        verts(1,3) = 0.57d0
        verts(2,3) = sqrt(3.0d0)/2.0d0
        
        xyout(1) = 1.0d0
        xyout(2) = 0.1d0

        xyin(1) = 0.3d0
        xyin(2) = 0.2d0
      endif
      


      return
      end
c
c
c
c       main subroutine
c


      subroutine helm_solver(zk,nverts,verts,ifdn,ifinout,iffast,fker,
     1   dpars,zpars,ipars,ntarg,xytarg,pot,ifwrite,fname)
      implicit real *8 (a-h,o-z)

c
c       calling sequence variables
c 

      
      real *8 verts(2,nverts)
      real *8 dpars(*)
      complex *16 zpars(*),zk
      integer ipars(*)

      real *8 xytarg(2,ntarg)
      complex *16 pot(ntarg)

      character (len=*) fname
c
c
c      temporary variables
c
      real *8, allocatable :: xys(:,:),dxys(:,:),pl(:),pr(:),qwts(:)
      real *8, allocatable :: xs(:),ys(:)
      integer, allocatable :: el(:),er(:),imid(:)
      integer, allocatable :: iscorn(:),ixys(:),ks(:),nepts(:)
      integer, allocatable :: lns(:),rns(:)

      real *8, allocatable :: vtmp(:,:),rpan(:),rmid(:),angs(:)
      
      external fker

      done = 1
      pi = atan(done)*4

c
c      find lenght of each of the edges
c
      allocate(vtmp(2,nverts+2),rpan(nverts+1))

      do i=1,nverts+2
        ii = i
        if(ii.gt.nverts) ii = ii-nverts
        vtmp(1,i) = verts(1,ii)
        vtmp(2,i) = verts(2,ii)
      enddo

      do i=1,nverts+1
        rr = (vtmp(1,i+1) - vtmp(1,i))**2 + (vtmp(2,i+1)-vtmp(2,i))**2
        rpan(i) = sqrt(rr)
      enddo

      allocate(pl(nverts+1),pr(nverts+1))


c
c       compute lengths of corner panels
c
      do i=1,nverts
        rmin = 0.15d0*min(rpan(i+1),rpan(i))
        rr = pi/real(zk)
        if(rmin.gt.rr) rmin = rr
        pr(i) = rmin
        pl(i+1) = rmin
      enddo
      pl(1) = pl(nverts+1)
      pr(nverts+1) = pr(1)

      allocate(rmid(nverts+1),angs(nverts+1),imid(nverts+1))

      do i=2,nverts+1
        dx1 = vtmp(1,i-1)-vtmp(1,i)
        dy1 = vtmp(2,i-1)-vtmp(2,i)

        dx2 = vtmp(1,i+1)-vtmp(1,i)
        dy2 = vtmp(2,i+1)-vtmp(2,i)
        
        angs(i) = atan2(dy2,dx2)-atan2(dy1,dx1)
      enddo
      angs(1) = angs(nverts+1)

      nch = 0
      npts = 0
      kmid = 16

      call prin2('angs=*',angs,nverts+1)


      do i=1,nverts
        rmid(i) = rpan(i) -pl(i)-pr(i)
        hmid = 2.0d0*min(pl(i),pr(i))
        hmida = 2.0d0*min(pl(i)*abs(tan(angs(i))),
     1      pr(i)*abs(tan(angs(i+1))))
        print *, i,hmid,hmida
        if(hmid.gt.hmida) hmid = hmida
        imid(i) = rmid(i)/hmid + 1

        nch = nch + imid(i)
        npts = npts  + imid(i)*kmid
      enddo

      allocate(el(nverts),er(nverts))
      do i=1,nverts
        el(i) = i
        er(i) = i+1
      enddo
      er(nverts) = 1

      nedges = nverts

      nch = nch + 2*nedges
      npts = npts + 72*nedges

      allocate(xys(2,npts),dxys(2,npts),qwts(npts))
      allocate(nepts(nedges),lns(nedges),rns(nedges))
      allocate(ixys(nch),ks(nch),iscorn(nch))

      call getedgedis(nverts,verts,nedges,el,er,pl,pr,imid,kmid,
     1       npts,xys,dxys,qwts,lns,rns,nepts,nch,ixys,ks,iscorn)
      
      allocate(xs(npts),ys(npts))
      do i=1,npts
        xs(i) = xys(1,i)
        ys(i) = xys(2,i)
      enddo


      call pyplot2(12,xverts,yverts,nverts,3,xs,ys,npts,
     1    1,'a*')

      call prin2('verts=*',verts,2*nverts)
      call prinf('el=*',el,nverts)
      call prinf('er=*',er,nverts)
      call prin2('pl=*',pl,nverts)
      call prin2('pr=*',pr,nverts)
      call prinf('imid=*',imid,nverts)
      call prinf('kmid=*',kmid,1)
      call prinf('npts=*',npts,1)
      call prinf('ixys=*',ixys,nch)
      call prinf('ks=*',ks,nch)
      call prinf('iscorn=*',iscorn,nch)
      return
      end

c
c
c

       subroutine getcornstruct(nr0,ts,wts,u,v)
       implicit real *8 (a-h,o-z)
       integer nr0
       real *8 ts(nr0),wts(nr0),u(nr0,nr0),v(nr0,nr0)

       itype = 1
       call lapdisc(ts,wts,u,v,nr0,itype)


       return
       end

       subroutine getedgedis(nverts,verts,nedges,el,er,pl,pr,
     1    imid,kmid,n,xys,dxyds,qwts,lns,rns,nepts,npan,
     2    ichstart,ks,iscorn)
       
       implicit real *8 (a-h,o-z)

c        this subroutine computes the discretization for the given
c        geometry
c  
c       input:
c       nverts - number of vertices
c       verts(2,nverts) - locations of vertices
c       nedges - number of edges
c       el(i)  - denotes the left end vertex on the edge 
c       er(i)  - denotes the right end vertex on the edge 
c       pl(i) - length of panel corresponding to corner discretization
c               on left end vertex of the edge
c       pr(i) - length of panel corresponding to corner discretization
c               on right end vertex of the edge
c      imid(i) - denotes the number of gauss-legendre panels in
c                the middle of the edge
c      kmid - is the number of nodes on each of the gauss-legendre
c             panels
c      n  - total number of discretization nodes
c            = \sum_{i=1}^{nedges} imid(i)*kmid + 72*nedges
c      npan - total number of panels
c            = \sum_{i=1}^{nedges} imid(i) + 2*nedges
c      
c      output:
c       xys - x and y coordinates of discretization nodes
c       dxyds - x and y components of tangents at the
c                 discretization nodes
c       qwts - quadrature weights at the discretization nodes
c       lns(i) - denotes the starting location in the discretization 
c               nodes array where the corner discretization nodes
c               corresponding to the left end vertex start
c             
c       rns(i) - denotes the starting location in the discretization 
c               nodes array where the corner discretization nodes
c               corresponding to the right end vertex start
c
c       nepts(i) - total number of discretization nodes on edge i
c       ichstart(i) - index in xys array where nodes on panel i start
c      
c       ks(i) - number of points on the panel 
c       iscorn(i) - whether panel is a corner panel
c 
c
       real *8 verts(2,nverts),pl(nedges),pr(nedges)
       real *8, allocatable :: ts(:),wts(:),umatt(:,:),vmatt(:,:)
       integer el(nedges),er(nedges),imid(nedges)
       real *8 xys(2,n),dxyds(2,n),qwts(n)
       real *8, allocatable :: rlen(:)
       
       integer lns(nedges),rns(nedges),nepts(nedges),ichstart(npan)
       integer ks(npan),iscorn(npan)
       real *8, allocatable :: ts0(:),wts0(:),uk(:,:),vk(:,:)

       allocate(ts0(kmid),wts0(kmid),uk(kmid,kmid),vk(kmid,kmid))

       allocate(rlen(nedges))
     
       n = 0
       itype = 1
       call legeexps(itype,kmid,ts0,uk,vk,wts0)

       ncorner = 36

       allocate(ts(ncorner),wts(ncorner),umatt(ncorner,ncorner))
       allocate(vmatt(ncorner,ncorner))


       itype = 2
       call lapdisc(ts,wts,umatt,vmatt,ncorner,itype)

       ipan = 0

       do ie=1,nedges 

         lns(ie) = n+1
         ipan = ipan + 1
         ichstart(ipan) = n+1
         ks(ipan) = ncorner
         iscorn(ipan) = 1

         dxt = verts(1,er(ie))-verts(1,el(ie))
         dyt = verts(2,er(ie))-verts(2,el(ie))


         dst = sqrt(dxt**2+dyt**2)

         rlen(ie) = dst
         dxt = dxt/dst
         dyt = dyt/dst

         do i=1,ncorner
           xys(1,n+i) = verts(1,el(ie)) + dxt*pl(ie)*ts(i)
           xys(2,n+i) = verts(2,el(ie)) + dyt*pl(ie)*ts(i)
           dxyds(1,n+i) = dxt*pl(ie)
           dxyds(2,n+i) = dyt*pl(ie)

           qwts(n+i) = pl(ie)*wts(i)
         enddo

         n = n+ ncorner
         rmid = rlen(ie) - pl(ie)-pr(ie)
       
         rpanel =  rmid/imid(ie)
         rstart = pl(ie)


         ii = 1

         do i=1,imid(ie)
           ipan = ipan + 1
           ichstart(ipan) = n+ii
           ks(ipan) = kmid
           iscorn(ipan) = 0
           do j=1,kmid
             xys(1,n+ii) = verts(1,el(ie)) + 
     1           dxt*(rstart+(ts0(j)+1)/2*rpanel)
             xys(2,n+ii) = verts(2,el(ie)) + 
     2           dyt*(rstart+(ts0(j)+1)/2*rpanel)
             dxyds(1,n+ii) = dxt*rpanel/2
             dxyds(2,n+ii) = dyt*rpanel/2
             qwts(n+ii) = wts0(j)*rpanel/2
             ii = ii+1
           enddo
           rstart = rstart + rpanel 
         enddo
         n = n + imid(ie)*kmid

         ipan = ipan+1
         rns(ie) = n+1
         ichstart(ipan) = n+1
         ks(ipan) = ncorner
         iscorn(ipan) = 1

         do i=1,ncorner
           xys(1,n+i) = verts(1,er(ie)) - dxt*ts(i)*pr(ie)
           xys(2,n+i) = verts(2,er(ie)) - dyt*ts(i)*pr(ie)
           dxyds(1,n+i) = dxt*pr(ie)
           dxyds(2,n+i) = dyt*pr(ie)

           qwts(n+i) = pr(ie)*wts(i)
         enddo
         nepts(ie) = 2*ncorner + imid(ie)*kmid
         n = n+ ncorner
       enddo


       return
       end

c
c
c
c
c
c
        subroutine helmcornmat(alpha,zk,rpan,x0s,w0s,u,v,nr0,coefs,
     1     lcoefs,akern)
        implicit real *8 (a-h,o-z)
        complex *16 zk, akern(nr0,nr0)
        complex *16, allocatable :: work(:)
        real *8 coefs(lcoefs),t0,pi
        real *8, allocatable :: rquads(:),wquads(:)
        real *8 x0s(nr0),w0s(nr0),u(nr0,nr0),v(nr0,nr0)
        data pi/3.141592653589793d0/

        allocate(work(20000),rquads(1000),wquads(1000))

        t0 = alpha/pi

        if(t0.lt.0) t0 = t0+2
    
        iw = 1
        iw2 = iw+ nr0*nr0*3+20

        do ipt=1,nr0
          call lap_load_quads(t0,rquads,wquads,nquads,ipt)
          call helmkern(x0s(ipt),w0s(ipt),nr0,t0,rquads,
     1      wquads,nquads,v,
     2      work(iw),work(iw2),zk,rpan,akern,ipt,coefs,lcoefs)
        enddo

        return
        end
c
c
c
c
c
        subroutine helmkern(x0,w0,nr0,ang,rq,wq,nq,umat,
     1      xmat,work,zk,rpan,amat,ipt,coefs,lcoefs)
        implicit real *8 (a-h,o-z)
        real *8 pi
        dimension coefs(lcoefs),rq(nq),wq(nq),umat(nr0,nr0),
     1      xmat(nq,nr0)
        complex *16 zk,amat(nr0,nr0),work(*),z,h0,h1,imainv4,val
        integer ifexpon
        data imainv4/(0.0d0,0.25d0)/
        data pi/3.141592653589793d0/

c                           Input parameters
c
c       x0: the target point
c       w0: the weight associated with the target point x0
c       nr0: the number of discretization nodes
c       rq: the quadrature nodes
c       wq: the quadrature weights
c       umat: the val to coef interp matrix for the discretization nodes
c
c                           Output parameters
c
c       amat: the row in the kernel matrix corresponding to the point x0


        ifexpon = 1
        done = 1
        pi = atan(done)*4

        do 1000 i=1,nq
    
        top = x0*sin(pi*ang)*rpan
        bottom = sqrt((cos(pi*ang)*x0-rq(i))**2+
     1      (x0*sin(pi*ang))**2)*rpan
         
          
         z = zk*bottom
         call hank103(z,h0,h1,ifexpon)
         work(i) = -imainv4*h1*zk*top/bottom*wq(i)*rpan
 1000 continue

        do 2000 i=1,nq

        do 1500 j=1,nr0

        call lapnestev2(ier,coefs,rq(i),j,xmat(i,j))
    
 1500 continue
 2000 continue



        do 3000 i=1,nr0

        val = 0

        do 2500 j=1,nq

        val = val+work(j)*xmat(j,i)

 2500 continue


        amat(ipt,i) = val

 3000 continue

        do 4000 i=1,nr0

        val = 0

        do 3500 j=1,nr0 

        val=val+amat(ipt,j)*umat(j,i)

 3500 continue

        work(i)=val

 4000 continue

        do 4100 i=1,nr0

        amat(ipt,i) = sqrt(w0)*work(i)

 4100 continue

        return
        end

c
c
c
c
c
c
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c           The following are various methods for constructing
c           and returning quadratures, weights, and kernels for
c           the Laplace equation with corners.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
c
        subroutine lapinterps(x,rhs,val)
        implicit real *8 (a-h,o-z)
        dimension rquads(100),wquads(200),x0s(36),w0s(36),rhs(1),
     1      u(36,36),v(36,36),work(9000),cfs(36)
c
c       This subroutine interpolates the value of a function near a
c       corner given the value of the function at a set of
c       discretization nodes.
c
c                           Input parameters
c
c       x: the value at which the function is to be interpolated
c       rhs: the values of the function sampled at the discretization
c               nodes (not scaled by the square roots of the weights).
c
c                           Output parameters
c
c       val: the interpolated value of the function at the point x.

        save
        
        data ifirst/1/

        lww = 4000

        if (ifirst .eq. 1) then

        itype = 1
        call lapdisc(x0s,w0s,u,v,nr0,itype)

        ifirst = 0

        endif


        do 5000 i=1,nr0

        val = 0
    
        do 4500 j=1,nr0

        val = val+v(i,j)*sqrt(w0s(j))*rhs(j)

 4500 continue

        cfs(i) = val
        
 5000 continue

        do 5500 j=1,nr0

        call lapeval(x,j,vt)
        if (ier .ne. 0) then
            call prinf('Bombing in check interp nestev*',0,0)
            call prinf('ier = *',ier,1)
            stop
        endif

        val = val+vt*cfs(j)

 5500 continue

        return
        end
c
c
c
c
c
        subroutine lapeval(x,ifun,val)
        implicit real *8 (a-h,o-z)
        dimension work(9000),coefs(200 000)
c
c       This subroutine evaluates the singular functions for the powers
c       x^mu for 0.5<mu<40.
c
c                           Input parameters
c
c       x: the point at which to evaluate the function
c       ifun: the index of the singular function
c
c                           Output parameters
c
c       val: the value of the ifun singular function at the point x.
c
        save

        data isvd/32/
        data ifirst/1/

        lww = 4000

        if (ifirst .eq. 1) then

        rewind(isvd)
        
 1500 format(i8)
        
        read(isvd,1500) keep
        call laparrread(isvd,coefs,keep)

        ifirst = 0

        endif

        call lapnestev(ier,coefs,x,ifun,val)
        if (ier .ne. 0) then
            call prinf('Bombing in check interp nestev*',0,0)
            call prinf('ier = *',ier,1)
            stop
        endif

        return
        end
c
c
c
c
c
        subroutine laparrread(ir,arr,n)
        implicit real *8 (a-h,o-z)
        dimension arr(n)
c
 1200 format(6x,d38.32)
        read(ir,1200) (arr(i),i=1,n)
c
        return
        end
c
c
c
c
c
        subroutine lapnestev(ier,coefs,x,i,val)
        implicit real *8 (a-h,o-z)
        dimension coefs(1),pjcoefs1(300),pjcoefs2(300)

c
c       This subroutine is a copy of nestev2 for evaluating functions
c       generated via allsvcmb. It has no stand-alone purpose.
c

        save
        data ifcalled/0/
c
c       decode the beginning of the array coefs
c
        k=coefs(6)
        nn=coefs(5)
        iab=coefs(1)
        n=coefs(7)
        icoefs=coefs(8)
c
c       . . . find the subinterval in which the point x lives
c
        ier=0
        call findinte(ier,x,coefs(iab),nn,intnum)
        if(ier .ne. 0) return
c
c       the point x lives in the interval number intnum.
c       evaluate the expansion at this point
c
        iii=iab+intnum*2-1
        u=2/(coefs(iii)-coefs(iii-1))
        v=1-coefs(iii)*u
c
        t=u*x+v
c
        jj=(intnum-1)*k+1
c
        ninit=290
        if(ifcalled .eq. 1) ninit=0
        ifcalled=1
c
        iijj=(i-1)*nn*k+jj +icoefs -1
        call legeexe2(t,VAL,coefs(iijj),k-1,
     1      pjcoefs1,pjcoefs2,ninit)
c
        return
        end
c
c
c
c
c
c
        subroutine lapnestev2(ier,coefs,x,i,val)
        implicit real *8 (a-h,o-z)
        dimension coefs(*)

c
c       This subroutine is a copy of nestev2 for evaluating functions
c       generated via allsvcmb. It has no stand-alone purpose.
c

c
c       decode the beginning of the array coefs
c
        k=coefs(6)
        nn=coefs(5)
        iab=coefs(1)
        n=coefs(7)
        icoefs=coefs(8)
c
c       . . . find the subinterval in which the point x lives
c
        ier=0
        call findinte(ier,x,coefs(iab),nn,intnum)
        if(ier .ne. 0) return
c
c       the point x lives in the interval number intnum.
c       evaluate the expansion at this point
c
        iii=iab+intnum*2-1
        u=2/(coefs(iii)-coefs(iii-1))
        v=1-coefs(iii)*u
c
        t=u*x+v
c
        jj=(intnum-1)*k+1
c
c
        iijj=(i-1)*nn*k+jj +icoefs -1
        call legeexev(t,val,coefs(iijj),k-1)
c
        return
        end
c
c
c
        subroutine findinte(ier,x,ab,nn,intnum)
        implicit real *8 (a-h,o-z)
        dimension ab(2,nn)
c
        data intold/-10/,ithresh/10/
c
c       check if the point is on the subinterval as the preceding one
c
        ier=0
        if(intold .lt. 0) goto 2000
c
        intnum=intold
        if( (x .ge. ab(1,intnum) ) .and. (x .le. ab(2,intnum) ) ) return
c
 2000 continue
c
c      the point is not on the same subinterval as the preceding one.
c      if nn is less than ithresh, use direct scan to find the proper 
c      interval
c
        
       if(nn .gt. ithresh) goto 3000
c
       if(x .lt. ab(1,1)) then
           ier=4
           return
       endif
c
        do 2200 j=1,nn
c
           intnum=j
c
        if(ab(2,j) .ge. x) goto 2400
 2200 continue
c
        ier=4
        return
 2400 continue
c
        intold=intnum
        return
c
 3000 continue
c
c      The point is not on the same subinterval as the preceding one,
c      and nn is greater than ithresh; use bisection to find the proper 
c      interval
c
       i1=1
       i2=nn
       i3=(i1+i2)/2
c
cccc       nsteps=0
       do 3400 i=1,100
c
       if(x .ge. ab(1,i3)) i1=i3
       if(x .le. ab(2,i3)) i2=i3
c
       if(i2 .eq. i1) goto 3600
c
       i3=(i1+i2)/2
 3400 continue
c
 3600 continue

       if(x .lt. ab(1,i3)) i3=i3-1
       if(x .gt. ab(2,i3)) i3=i3+1

       intnum=i3
       intold=intnum
c       
        return
        end
c
c
c
c
c
c
        subroutine lapdisc(roots,weights,u,v,npts,itype)
        implicit real *8 (a-h,o-z)
        dimension roots(1), weights(1), u(1), v(1),rts(36),wts(36),
     1      work(16000)

c
c       This subroutine returns a set of discretization nodes and
c       weights for Laplace near corners. Depending on the value of
c       itype it will also output matrices mapping the coefficients of
c       the singular function to their values at the discretization
c       nodes (and back).
c
c                           Input parameters
c
c       itype: if itype is greater than or equal to one then the
c               interpolation matrices are returned. Otherwise not.
c
c                           Output parameters
c
c       roots: the locations of the discretization nodes.
c       weights: the corresponding weights for the discretization nodes.
c       npts: the number of discretization nodes.
c       u: the coefficient-to-value 'interpolation' matrix.
c       v: the value-to-coefficient 'interpolation' matrix.
c

        nr0  = 36
        npts = 36

        data rts/
     1  0.7153241747307978D-09,0.3377836743170287D-07,
     2  0.5130805470357669D-06,0.4302091983344161D-05,
     3  0.2438570399508877D-04,0.1038302657440851D-03,
     4  0.3542151321928258D-03,0.1011404847825894D-02,
     5  0.2494729075160877D-02,0.5444867944769177D-02,
     6  0.1071566026812782D-01,0.1930823611072878D-01,
     7  0.3225634017695757D-01,0.5048916587188894D-01,
     8  0.7470463868724012D-01,0.1052804775045622D+00,
     9  0.1422368424680723D+00,0.1852494478758348D+00,
     *  0.2337009173467016D+00,0.2867532037701939D+00,
     1  0.3434246385591752D+00,0.4026595369943845D+00,
     2  0.4633838423882242D+00,0.5245451728362305D+00,
     3  0.5851389017555692D+00,0.6442234580910580D+00,
     4  0.7009282778474498D+00,0.7544573296534130D+00,
     5  0.8040903461004613D+00,0.8491831219651994D+00,
     6  0.8891676311778808D+00,0.9235522974976429D+00,
     7  0.9519225051333568D+00,0.9739413233554091D+00,
     8  0.9893505102238157D+00,0.9979739287518251D+00/

c        
c        Weights:
c        
        data wts/
     1  0.3488387406713975D-08,0.1060713246829122D-06,
     2  0.1219072309361280D-05,0.8208292861991787D-05,
     3  0.3849614123472857D-04,0.1380261910723107D-03,
     4  0.4010043380017863D-03,0.9825669945713395D-03,
     5  0.2090967113128153D-02,0.3952737459803072D-02,
     6  0.6757680951805507D-02,0.1060213136740334D-01,
     7  0.1545112262562957D-01,0.2113236530646019D-01,
     8  0.2736215158818752D-01,0.3379263575171562D-01,
     9  0.4006514738343413D-01,0.4585591975580118D-01,
     *  0.5090603137005028D-01,0.5503375475249415D-01,
     1  0.5813223680645907D-01,0.6015788774245568D-01,
     2  0.6111489952499047D-01,0.6104007704184891D-01,
     3  0.5999032684117640D-01,0.5803363810654483D-01,
     4  0.5524331217551819D-01,0.5169474359179677D-01,
     5  0.4746392056021400D-01,0.4262695578774503D-01,
     6  0.3726012465109946D-01,0.3144009320827304D-01,
     7  0.2524415159573378D-01,0.1875041493753455D-01,
     8  0.1203832814441368D-01,0.5196613268519835D-02/

        do 1000 i=1,npts

        roots(i) = rts(i)
        weights(i) = wts(i)

 1000 continue

        if (itype .gt. 0) then

        lww = 16000
        call lapcreainterp(u,v,npts,npts,wts,
     1      rts,work,lww)

        endif

        return
        end
c
c
c
c
c
        subroutine lapcreainterp(xmat,vmat,nroots,nfuns,weights,
     1      roots,work,lw)
        implicit real *8 (a-h,o-z)
        dimension xmat(nroots,nfuns),vmat(nroots,nfuns),weights(1),
     1      work(1),roots(1)

c
c       This subroutine creates the interpolation matrices taking the
c       coefficients of functions in the basis of singular functions to
c       their values at the discretization nodes (and the inverse of the
c       corresponding matrix). Both matrices are scaled by the square
c       root of the weights of the discretization nodes.
c
c                               Input parameters
c       
c       roots: the positions of the discretization nodes.
c       weights: the corresponding weights.
c       nfuns and nroots: notice that these will always be the same.
c       work: a work array sufficiently large for orthom.
c
c                               Output parameters
c       
c       xmat: coefficient-to-value interpolation matrix.
c       vmat: value-to-coefficient interpolation matrix.


        eps = 1.0d-16

        if (nfuns .ne. nroots) then
        call prinf('bombing in creainterp, incorrect dims*',0,0)
        stop
        end if

        par = 0 

        do 2000 i=1,nroots
        do 1000 j=1,nfuns

        call lapeval(roots(i),j,xmat(i,j))

        xmat(i,j) = xmat(i,j)*sqrt(weights(i))
        vmat(i,j) = xmat(i,j)

 1000 continue
 2000 continue
c
        call orthom(vmat,nfuns,work,cond)
c
        return
        end
c
c
c
c
c
        subroutine lap_load_quads(ang,roots,weights,nqs,ipt)
        implicit real *8 (a-h,o-z)
        dimension roots(1), weights(1)

        alpha = ang
        if (ang .gt. 1) then
        alpha = 2-ang
ccc        call prin2('ang = *',alpha,1)
        endif

        if (ipt .eq. 1) then
        call lap_load_quads1(alpha,roots,weights,nqs)
        endif

        if (ipt .eq. 2) then
        call lap_load_quads2(alpha,roots,weights,nqs)
        endif

        if (ipt .eq. 3) then
        call lap_load_quads3(alpha,roots,weights,nqs)
        endif

        if (ipt .eq. 4) then
        call lap_load_quads4(alpha,roots,weights,nqs)
        endif

        if (ipt .eq. 5) then
        call lap_load_quads5(alpha,roots,weights,nqs)
        endif

        if (ipt .eq. 6) then
        call lap_load_quads6(alpha,roots,weights,nqs)
        endif

        if (ipt .eq. 7) then
        call lap_load_quads7(alpha,roots,weights,nqs)
        endif

        if (ipt .eq. 8) then
        call lap_load_quads8(alpha,roots,weights,nqs)
        endif

        if (ipt .eq. 9) then
        call lap_load_quads9(alpha,roots,weights,nqs)
        endif

        if (ipt .eq. 10) then
        call lap_load_quads10(alpha,roots,weights,nqs)
        endif

        if (ipt .eq. 11) then
        call lap_load_quads11(alpha,roots,weights,nqs)
        endif

        if (ipt .eq. 12) then
        call lap_load_quads12(alpha,roots,weights,nqs)
        endif

        if (ipt .eq. 13) then
        call lap_load_quads13(alpha,roots,weights,nqs)
        endif

        if (ipt .eq. 14) then
        call lap_load_quads14(alpha,roots,weights,nqs)
        endif

        if (ipt .eq. 15) then
        call lap_load_quads15(alpha,roots,weights,nqs)
        endif

        if (ipt .eq. 16) then
        call lap_load_quads16(alpha,roots,weights,nqs)
        endif

        if (ipt .eq. 17) then
        call lap_load_quads17(alpha,roots,weights,nqs)
        endif

        if (ipt .eq. 18) then
        call lap_load_quads18(alpha,roots,weights,nqs)
        endif

        if (ipt .eq. 19) then
        call lap_load_quads19(alpha,roots,weights,nqs)
        endif

        if (ipt .eq. 20) then
        call lap_load_quads20(alpha,roots,weights,nqs)
        endif

        if (ipt .eq. 21) then
        call lap_load_quads21(alpha,roots,weights,nqs)
        endif

        if (ipt .eq. 22) then
        call lap_load_quads22(alpha,roots,weights,nqs)
        endif

        if (ipt .eq. 23) then
        call lap_load_quads23(alpha,roots,weights,nqs)
        endif

        if (ipt .eq. 24) then
        call lap_load_quads24(alpha,roots,weights,nqs)
        endif

        if (ipt .eq. 25) then
        call lap_load_quads25(alpha,roots,weights,nqs)
        endif

        if (ipt .eq. 26) then
        call lap_load_quads26(alpha,roots,weights,nqs)
        endif

        if (ipt .eq. 27) then
        call lap_load_quads27(alpha,roots,weights,nqs)
        endif

        if (ipt .eq. 28) then
        call lap_load_quads28(alpha,roots,weights,nqs)
        endif

        if (ipt .eq. 29) then
        call lap_load_quads29(alpha,roots,weights,nqs)
        endif

        if (ipt .eq. 30) then
        call lap_load_quads30(alpha,roots,weights,nqs)
        endif

        if (ipt .eq. 31) then
        call lap_load_quads31(alpha,roots,weights,nqs)
        endif

        if (ipt .eq. 32) then
        call lap_load_quads32(alpha,roots,weights,nqs)
        endif

        if (ipt .eq. 33) then
        call lap_load_quads33(alpha,roots,weights,nqs)
        endif

        if (ipt .eq. 34) then
        call lap_load_quads34(alpha,roots,weights,nqs)
        endif

        if (ipt .eq. 35) then
        call lap_load_quads35(alpha,roots,weights,nqs)
        endif

        if (ipt .eq. 36) then
        call lap_load_quads36(alpha,roots,weights,nqs)
        endif

        return
        end
c
c
c
c
c
        subroutine lap_load_quads1(alpha,roots,weights,nquads)
        implicit real *8 (a-h,o-z)
        dimension r1(35),r2(35),r3(35),r4(35),r5(35),r6(31),
     1          w1(35),w2(35),w3(35),w4(35),w5(35),w6(31),nqs(6),
     2          roots(1),weights(1)

        data nqs/35,35,35,35,35,31/

        data r1/
     1  0.5366556771387501D-11,0.1647617655595256D-10,
     2  0.6111992295278005D-10,0.4111737415339524D-09,
     3  0.5981481753613897D-09,0.6530850552669246D-09,
     4  0.6800311286538681D-09,0.6977085376956023D-09,
     5  0.7119322624298200D-09,0.7257621260784202D-09,
     6  0.7420599755680398D-09,0.7656160423659988D-09,
     7  0.8101087718434004D-09,0.9515471494487160D-09,
     8  0.2532393009736539D-08,0.3273696204083351D-07,
     9  0.6625069631590972D-06,0.1151591155196457D-04,
     *  0.1248127356206399D-03,0.8378735650232596D-03,
     1  0.3793167448778285D-02,0.1257372887111414D-01,
     2  0.3259314922435407D-01,0.6954798487483777D-01,
     3  0.1271880771100305D+00,0.2058047030989370D+00,
     4  0.3021781650911163D+00,0.4106465365107982D+00,
     5  0.5244702648605223D+00,0.6368896586962888D+00,
     6  0.7417221730571566D+00,0.8336252929439057D+00,
     7  0.9081972194612087D+00,0.9620244968515138D+00,
     *  0.9927267918020293D+00/

        data w1/
     1  0.2292199769624283D-10,-.2004980523866348D-10,
     2  0.1718323502850508D-09,0.3577576415685633D-09,
     3  0.8625274197473594D-10,0.3549012474210544D-10,
     4  0.2087797955286745D-10,0.1533583368273358D-10,
     5  0.1358863536312142D-10,0.1452006538418858D-10,
     6  0.1879563100801639D-10,0.3017732414593796D-10,
     7  0.6683673146640176D-10,0.2971410573713252D-09,
     8  0.5011115196447901D-08,0.9491441960674334D-07,
     9  0.1999419067301139D-05,0.3038610982543869D-04,
     *  0.2660209915294601D-03,0.1417759349491064D-02,
     1  0.5090079031328400D-02,0.1340125553962105D-01,
     2  0.2762817456363217D-01,0.4693359606195215D-01,
     3  0.6839487775012053D-01,0.8828165246169287D-01,
     4  0.1034973626370962D+00,0.1122982308570652D+00,
     5  0.1142130607484336D+00,0.1095913709918301D+00,
     6  0.9918181415587997D-01,0.8389009008551759D-01,
     7  0.6468421249574790D-01,0.4257215354601376D-01,
     *  0.1862580214714210D-01/
        
        data r2/
     1  0.5209659790520027D-12,0.1715906146393313D-10,
     2  0.1579567933515668D-09,0.4280001296316204D-09,
     3  0.5611003865462505D-09,0.6249999971344377D-09,
     4  0.6647239566559766D-09,0.6948522546911944D-09,
     5  0.7221746917887031D-09,0.7519477545053957D-09,
     6  0.7912757715557787D-09,0.8561366000436423D-09,
     7  0.1005997951397069D-08,0.1682143071930841D-08,
     8  0.8668780155991793D-08,0.1062310565820858D-06,
     9  0.1641819084540212D-05,0.2113434310470096D-04,
     *  0.1846248453088460D-03,0.1081774054932928D-02,
     1  0.4494394644383011D-02,0.1409245356338990D-01,
     2  0.3520880040652434D-01,0.7328988879705800D-01,
     3  0.1317938721829506D+00,0.2108225399517790D+00,
     4  0.3071253626327841D+00,0.4151318377928457D+00,
     5  0.5282446266281935D+00,0.6398417679476488D+00,
     6  0.7438500344643960D+00,0.8350057617332579D+00,
     7  0.9089617354864849D+00,0.9623413174866205D+00,
     *  0.9927875171122753D+00/

        data w2/
     1  0.2486359881070530D-11,0.4619330079089234D-10,
     2  0.2644313430812053D-09,0.2008634984228429D-09,
     3  0.8606923722372180D-10,0.4811992373787808D-10,
     4  0.3347018680049364D-10,0.2783693283560273D-10,
     5  0.2763809379160480D-10,0.3300467349251262D-10,
     6  0.4796185216572196D-10,0.8929466911914309D-10,
     7  0.2528030001587562D-09,0.1603582799582197D-08,
     8  0.1922326376435237D-07,0.2868866090179167D-06,
     9  0.4436931557026539D-05,0.5018562777604112D-04,
     *  0.3620135607794045D-03,0.1715025242622284D-02,
     1  0.5724283199068594D-02,0.1438865865474966D-01,
     2  0.2879046334229536D-01,0.4797251666544075D-01,
     3  0.6905122740504005D-01,0.8844502582151465D-01,
     4  0.1032092720155450D+00,0.1116871809363262D+00,
     5  0.1134255533788581D+00,0.1087525912015391D+00,
     6  0.9838511843798887D-01,0.8320098525538202D-01,
     7  0.6414747674467172D-01,0.4221735940237304D-01,
     *  0.1847031730284352D-01/

        data r3/
     1  0.1312569502711661D-13,0.9665957112285794D-11,
     2  0.4710757860974359D-10,0.1132434549575401D-09,
     3  0.3184834745330229D-09,0.4750004279498002D-09,
     4  0.5723354392463968D-09,0.6417614999502461D-09,
     5  0.7007547229570332D-09,0.7617076909786532D-09,
     6  0.8399076551438375D-09,0.9677884818646642D-09,
     7  0.1258709146452518D-08,0.2401801578369711D-08,
     8  0.1117245968198632D-07,0.1185720593600601D-06,
     9  0.1744005334121912D-05,0.2210507638742527D-04,
     *  0.1911883174590722D-03,0.1110983712515486D-02,
     1  0.4584893953886909D-02,0.1430100427978574D-01,
     2  0.3558684505625448D-01,0.7385394834146841D-01,
     3  0.1325124181175333D+00,0.2116274161050390D+00,
     4  0.3079366984772393D+00,0.4158803582108978D+00,
     5  0.5288830555907543D+00,0.6403462655633777D+00,
     6  0.7442164746185503D+00,0.8352448294179939D+00,
     7  0.9090946557163836D+00,0.9623965415749393D+00,
     *  0.9927981164105578D+00/

        data w3/
     1  0.4656137137289190D-12,0.3011055915353166D-10,
     2  0.1320862574322122D-10,0.1692614327491604D-09,
     3  0.1967466609479557D-09,0.1204443086006678D-09,
     4  0.7957031207041982D-10,0.6200679989417723D-10,
     5  0.5793202511499134D-10,0.6634213063412267D-10,
     6  0.9478382621171203D-10,0.1757272561188046D-09,
     7  0.4807948643566282D-09,0.2449895977073532D-08,
     8  0.2264327104337715D-07,0.3098014999667746D-06,
     9  0.4671524822888308D-05,0.5224330343042511D-04,
     *  0.3731481083127526D-03,0.1753008150307155D-02,
     1  0.5812297710201225D-02,0.1453598987853833D-01,
     2  0.2897597953790274D-01,0.4815055827114680D-01,
     3  0.6917581889103487D-01,0.8849121013358440D-01,
     4  0.1031779795602242D+00,0.1115966836076141D+00,
     5  0.1132998561865894D+00,0.1086137250541226D+00,
     6  0.9825032413182784D-01,0.8308277420527945D-01,
     7  0.6405458531355087D-01,0.4215562322501696D-01,
     *  0.1844318676443172D-01/

        data r4/
     1  0.4134068042510158D-12,0.5351330333705171D-11,
     2  0.3025367485609098D-10,0.1232956215176688D-09,
     3  0.2701084317658861D-09,0.4095049613956526D-09,
     4  0.5276683789273348D-09,0.6342267926906696D-09,
     5  0.7441338737749390D-09,0.8798536929936753D-09,
     6  0.1086797297859527D-08,0.1493538243068581D-08,
     7  0.2628937181142153D-08,0.7746473740818010D-08,
     8  0.4636906545026061D-07,0.4886089467267708D-06,
     9  0.6124981436502247D-05,0.6199975983691773D-04,
     *  0.4340862574943042D-03,0.2118563598480827D-02,
     1  0.7605485965361278D-02,0.2123734537057477D-01,
     2  0.4840416738773842D-01,0.9369616921868880D-01,
     3  0.1590973689617414D+00,0.2432597425475658D+00,
     4  0.3420219515193725D+00,0.4495884696492851D+00,
     5  0.5596832413217629D+00,0.6663279952870817D+00,
     6  0.7642406501288724D+00,0.8489941335120379D+00,
     7  0.9170606813756916D+00,0.9658059223445785D+00,
     *  0.9934636790791986D+00/

        data w4/
     1  0.1649983859061274D-11,0.9333063620580560D-11,
     2  0.5142605023443081D-10,0.1317420043089895D-09,
     3  0.1494765820932035D-09,0.1279629012919057D-09,
     4  0.1101944781576128D-09,0.1053618155540904D-09,
     5  0.1179882819681084D-09,0.1602941232338844D-09,
     6  0.2716869079429530D-09,0.6090194911101528D-09,
     7  0.2026363869175409D-08,0.1107609756076640D-07,
     8  0.9839354200080973D-07,0.1230580369037732D-05,
     9  0.1512108546119964D-04,0.1325213079515256D-03,
     *  0.7631702652101177D-03,0.3012470178826648D-02,
     1  0.8698777227959519D-02,0.1949718086443111D-01,
     2  0.3563447758098515D-01,0.5529258046496931D-01,
     3  0.7527479041747457D-01,0.9233492909691760D-01,
     4  0.1042035265177506D+00,0.1098706181660086D+00,
     5  0.1093163079474507D+00,0.1030905166746063D+00,
     6  0.9199634892715432D-01,0.7692280332243890D-01,
     7  0.5877176922480016D-01,0.3842431712880607D-01,
     *  0.1674642967828948D-01/

        data r5/
     1  0.6478071761871498D-09,0.2096449553122768D-01,
     2  0.4673721487340153D-01,0.6595335611246786D-09,
     3  0.5175646295502460D-11,0.2822001328188289D-12,
     4  0.7788645766904316D-02,0.8946907183294964D-01,
     5  0.1514264737326203D+00,0.2293628804693522D-02,
     6  0.2319877066207031D+00,0.5088162356259250D-03,
     7  0.3278361285787289D+00,0.4338123053411188D+00,
     8  0.1806570932127369D-10,0.5439108091024668D+00,
     9  0.8091877161728405D-04,0.6520682856774026D+00,
     *  0.6701161757349363D-10,0.7526505126694929D+00,
     1  0.9125054978542732D-09,0.9108369428808025D-05,
     2  0.8407100062267510D+00,0.8230810242513981D-06,
     3  0.1320666106930265D-08,0.4708128489721399D-09,
     4  0.9121189662855137D+00,0.8554827158267787D-07,
     5  0.1698194052819045D-09,0.3093070832465000D-09,
     6  0.1476863804108549D-07,0.2165698057760137D-08,
     7  0.4559018119712391D-08,0.9636492286342817D+00,
     *  0.9930381593056095D+00/

        data w5/
     1  -.1477726984529933D-10,0.1864349429109287D-01,
     2  0.3365683834521095D-01,0.2266073895630594D-09,
     3  0.8313328507165564D-11,0.1300791005684248D-11,
     4  0.8554741378648688D-02,0.5219835592172685D-01,
     5  0.7160864760869978D-01,0.3108526104736665D-02,
     6  0.8894546229302795D-01,0.8472072950043104D-03,
     7  0.1018739881203503D+00,0.1090641601609281D+00,
     8  0.2371157642481558D-10,0.1101144001058027D+00,
     9  0.1631336288204146D-03,0.1052574968618884D+00,
     *  0.7702504067065786D-10,0.9508032970539065D-01,
     1  0.3064362155348750D-09,0.2126203106506305D-04,
     2  0.8034943014183074D-01,0.1982892795758362D-05,
     3  0.5496801415170908D-09,0.1725967577292516D-09,
     4  0.6192947104031722D-01,0.1754338119614575D-06,
     5  0.1249119592363475D-09,0.1513896585817447D-09,
     6  0.2145328478072134D-07,0.1285092941029919D-08,
     7  0.4225807665404617D-08,0.4075233775576936D-01,
     *  0.1782853029170032D-01/

        data r6/
     1  0.4864039469047956D-12,0.5476812527968156D-11,
     2  0.2779496339486087D-10,0.9284497525630514D-10,
     3  0.2243130155345543D-09,0.4489215047608004D-09,
     4  0.8338375965187342D-09,0.1587162953753886D-08,
     5  0.3435942991689221D-08,0.9613575928185074D-08,
     6  0.4037584460283148D-07,0.2776527612037404D-06,
     7  0.2702176952979200D-05,0.2641929117736848D-04,
     8  0.2013851737181476D-03,0.1109886372102987D-02,
     9  0.4495172757551667D-02,0.1396606536797713D-01,
     *  0.3482268471894703D-01,0.7255072168953097D-01,
     1  0.1307051111743677D+00,0.2094851665082472D+00,
     2  0.3056921467791303D+00,0.4137533505371346D+00,
     3  0.5270346053194151D+00,0.6388662713924229D+00,
     4  0.7431315209389860D+00,0.8345324319809254D+00,
     5  0.9086968324851899D+00,0.9622307996521897D+00,
     *  0.9927662586820376D+00/

        data w6/
     1  0.1654138045975191D-11,0.1017808708446708D-10,
     2  0.3916913771463986D-10,0.9477415902707518D-10,
     3  0.1720266403578981D-09,0.2870387761761137D-09,
     4  0.5129456058824634D-09,0.1093868966147377D-08,
     5  0.3012204873340283D-08,0.1160580468350229D-07,
     6  0.6787647402705465D-07,0.5971871869115283D-06,
     7  0.6312516179886936D-05,0.5761645448232835D-04,
     8  0.3767512156714162D-03,0.1716483233799729D-02,
     9  0.5656959230333329D-02,0.1419531222218418D-01,
     *  0.2847119176234911D-01,0.4760291992262654D-01,
     1  0.6873873471704765D-01,0.8826902290030011D-01,
     2  0.1031928360511044D+00,0.1118065070749223D+00,
     3  0.1136349660113104D+00,0.1090048863121153D+00,
     4  0.9864108968179042D-01,0.8343125270588773D-01,
     5  0.6433120752339988D-01,0.4234056444895384D-01,
     *  0.1852470412221541D-01/
c        
c         %%%%%%%%%%%%%%%%%%%%%

ccc        call prin2('alpha = *',alpha,1)
        dtwo = 2
        if (alpha .ge. 1/dtwo) then
        do 1000 i=1,nqs(6)
        roots(i) =r6(i)
        weights(i) = w6(i)
 1000 continue

        nquads = nqs(6)
        endif

        if (alpha .ge. dtwo**(-2)) then
        if (alpha .lt. dtwo**(-1)) then

        do 2000 i=1,nqs(5)
        roots(i) =r5(i)
        weights(i) = w5(i)
 2000 continue

        nquads = nqs(5)

        endif
        endif


        if (alpha .ge. dtwo**(-3)) then
        if (alpha .lt. dtwo**(-2)) then

        do 3000 i=1,nqs(4)
        roots(i) =r4(i)
        weights(i) = w4(i)
 3000 continue

        nquads = nqs(4)

        endif
        endif


        if (alpha .ge. dtwo**(-4)) then
        if (alpha .lt. dtwo**(-3)) then

        do 4000 i=1,nqs(3)
        roots(i) =r3(i)
        weights(i) = w3(i)
 4000 continue

        nquads = nqs(3)

        endif
        endif


        if (alpha .ge. dtwo**(-5)) then
        if (alpha .lt. dtwo**(-4)) then

        do 5000 i=1,nqs(2)
        roots(i) =r2(i)
        weights(i) = w2(i)
 5000 continue

        nquads = nqs(2)

        endif
        endif


        if (alpha .lt. dtwo**(-5)) then

        do 6000 i=1,nqs(1)
        roots(i) =r1(i)
        weights(i) = w1(i)
 6000 continue

        nquads = nqs(1)

        endif

        return
        end
c
c
c
c
c
        subroutine lap_load_quads2(alpha,roots,weights,nquads)
        implicit real *8 (a-h,o-z)
        dimension r1(35),r2(35),r3(35),r4(36),r5(35),r6(32),
     1          w1(35),w2(35),w3(35),w4(36),w5(35),w6(32),nqs(6),
     2          roots(1),weights(1)

        data nqs/35,35,35,36,35,32/


        data r1/
     1  0.8057070639869677D-10,0.1531987497915018D-08,
     2  0.1122986121441717D-07,0.2549185810824127D-07,
     3  0.2980623577969206D-07,0.3153879226335014D-07,
     4  0.3254774060894748D-07,0.3328671659712658D-07,
     5  0.3394127710634526D-07,0.3464143719406773D-07,
     6  0.3555069441345191D-07,0.3702467844349090D-07,
     7  0.4039545179703278D-07,0.5701031710530871D-07,
     8  0.2525112434627803D-06,0.2022444672349308D-05,
     9  0.1341823151095407D-04,0.6819521967621083D-04,
     *  0.3419151085135046D-03,0.1539283323548171D-02,
     1  0.5582485472244185D-02,0.1619437580749949D-01,
     2  0.3857577637164082D-01,0.7788420735495923D-01,
     3  0.1372749875366233D+00,0.2166720059113269D+00,
     4  0.3128149894078884D+00,0.4202452810973353D+00,
     5  0.5325235059695234D+00,0.6431765356659786D+00,
     6  0.7462482201880514D+00,0.8365593218532520D+00,
     7  0.9098213073489613D+00,0.9626973283617665D+00,
     *  0.9928557343220865D+00/

        data w1/
     1  0.3222467552105505D-09,0.3526434826815409D-08,
     2  0.1752114237410321D-07,0.7527327397437417D-08,
     3  0.2445553923758112D-08,0.1249744142076787D-08,
     4  0.8324612597229760D-09,0.6738439507749214D-09,
     5  0.6557552248271072D-09,0.7698422204395129D-09,
     6  0.1100123385171327D-08,0.2014745857768382D-08,
     7  0.5710276446243910D-08,0.4278674736860752D-07,
     8  0.5051827180733641D-06,0.4153031873703729D-05,
     9  0.2298735664211952D-04,0.1094617587072956D-03,
     *  0.5431099201013144D-03,0.2159856336570692D-02,
     1  0.6551406557824064D-02,0.1556668906032753D-01,
     2  0.3009164407154550D-01,0.4907298748177956D-01,
     3  0.6969401959158923D-01,0.8853826376652542D-01,
     4  0.1028163801854956D+00,0.1109549752143442D+00,
     5  0.1125137784338477D+00,0.1077956343038568D+00,
     6  0.9748284956496136D-01,0.8242366476421880D-01,
     7  0.6354337871379040D-01,0.4181852934304492D-01,
     *  0.1829563822399065D-01/

        data r2/
     1  0.6275677350274061D-10,0.1229454231178169D-08,
     2  0.8294008466313078D-08,0.2041139726395772D-07,
     3  0.2652938400860658D-07,0.2952001662365732D-07,
     4  0.3139011603265726D-07,0.3281169865872725D-07,
     5  0.3410186863784092D-07,0.3550740434686462D-07,
     6  0.3736124242294133D-07,0.4040605014138623D-07,
     7  0.4734841924605678D-07,0.7676790288576285D-07,
     8  0.3215933213765933D-06,0.2518474366111572D-05,
     9  0.1979034270785539D-04,0.1139994204632077D-03,
     *  0.4646899491720300D-03,0.1692275024712113D-02,
     1  0.5647059615998554D-02,0.1598442532490391D-01,
     2  0.3789202980143180D-01,0.7661825804700863D-01,
     3  0.1354736659507074D+00,0.2145178244498937D+00,
     4  0.3105505603774698D+00,0.4180964695275734D+00,
     5  0.5306546576246804D+00,0.6416793545199993D+00,
     6  0.7451501539662155D+00,0.8358380418860303D+00,
     7  0.9094184081187592D+00,0.9625294386708561D+00,
     *  0.9928234602795552D+00/
        
        data w2/
     1  0.2530957523401563D-09,0.2820658849340026D-08,
     2  0.1216447359494860D-07,0.9100506911426233D-08,
     3  0.4009083994729115D-08,0.2261275360064146D-08,
     4  0.1578081561373073D-08,0.1314212425450684D-08,
     5  0.1305126439408426D-08,0.1557388515311440D-08,
     6  0.2257974261894510D-08,0.4178040090852691D-08,
     7  0.1158155742092173D-07,0.6664606472447601D-07,
     8  0.6077267351834886D-06,0.5361552471756615D-05,
     9  0.3842923683590507D-04,0.1766738368972282D-03,
     *  0.6114295009244828D-03,0.2141971552514835D-02,
     1  0.6379662698952728D-02,0.1518619835436756D-01,
     2  0.2954049719966102D-01,0.4848720653343714D-01,
     3  0.6923266490144084D-01,0.8830386559595661D-01,
     4  0.1028268706266933D+00,0.1111648696969748D+00,
     5  0.1128513854200541D+00,0.1081905692660233D+00,
     6  0.9787787365711752D-01,0.8277624812289895D-01,
     7  0.6382343309505520D-01,0.4200584163915762D-01,
     *  0.1837821875829001D-01/

        data r3/
     1  0.2716734502692850D-10,0.7095316226801927D-09,
     2  0.4447376357415601D-08,0.1201292351129975D-07,
     3  0.1930360420205505D-07,0.2438253087003122D-07,
     4  0.2793811748464048D-07,0.3070614190090232D-07,
     5  0.3316589066724979D-07,0.3571150155470600D-07,
     6  0.3882648549852718D-07,0.4338809729460722D-07,
     7  0.5164769047087436D-07,0.7164362562876159D-07,
     8  0.1429216310780054D-06,0.5246922071375863D-06,
     9  0.3342057619693829D-05,0.2648449122064023D-04,
     *  0.1892213072737604D-03,0.1037743392489309D-02,
     1  0.4252821901838976D-02,0.1339537002654093D-01,
     2  0.3378898134782234D-01,0.7102358759782665D-01,
     3  0.1287852254604662D+00,0.2073637750154753D+00,
     4  0.3035806751765986D+00,0.4118267451628158D+00,
     5  0.5254063788419043D+00,0.6375890499240876D+00,
     6  0.7422090963832467D+00,0.8339332077067524D+00,
     7  0.9083646861532496D+00,0.9620930810344998D+00,
     *  0.9927398545659463D+00/

        data w3/
     1  0.1209315702935356D-09,0.1697508705439435D-08,
     2  0.6044025039372632D-08,0.8162306377102517D-08,
     3  0.6148932908096626D-08,0.4167244963800613D-08,
     4  0.3065938843360345D-08,0.2546759540908490D-08,
     5  0.2435832995523554D-08,0.2731802955566038D-08,
     6  0.3632014252670645D-08,0.5816215464387641D-08,
     7  0.1177833379005776D-07,0.3307115060162525D-07,
     8  0.1385984207652364D-06,0.8495708118524635D-06,
     9  0.6761911943602060D-05,0.5456166983213934D-04,
     *  0.3494491636277766D-03,0.1612663134720743D-02,
     1  0.5411231771788879D-02,0.1378806114411231D-01,
     2  0.2797085403727152D-01,0.4713951235087272D-01,
     3  0.6843223549792381D-01,0.8817556790862587D-01,
     4  0.1032994593539636D+00,0.1120589370032717D+00,
     5  0.1139689680540657D+00,0.1093647483621362D+00,
     6  0.9898497493141289D-01,0.8372973035872934D-01,
     7  0.6456416330396812D-01,0.4249473397419885D-01,
     *  0.1859230647930362D-01/

        data r4/
     1  0.7105235957609258D-08,0.5057962167597516D-03,
     2  0.1076186898979171D-03,0.1995397371584887D-02,
     3  0.1929907600516001D-04,0.9122010162329947D-08,
     4  0.6559968525577829D-02,0.3190324546335021D-05,
     5  0.1791581656466906D-01,0.2498617793627381D-10,
     6  0.1900194910720483D-06,0.9350356252962257D-07,
     7  0.5616278013332875D-09,0.4115472543795059D-01,
     8  0.1532595546239986D-07,0.6253917300637021D-07,
     9  0.8125222299052793D-01,0.3249723361804282D-08,
     *  0.3167873542248459D+00,0.2207860873001322D+00,
     1  0.4238049552134485D+00,0.2094774680084402D-07,
     2  0.5355006253401945D+00,0.4864372187654672D-07,
     3  0.6454981503596820D+00,0.1411876083517896D+00,
     4  0.7479194650363703D+00,0.2581683393119519D-07,
     5  0.8376430670404712D+00,0.6249435916467596D-06,
     6  0.4062589846360984D-07,0.9104214482506999D+00,
     7  0.3032272559287880D-07,0.3499969784359711D-07,
     8  0.9629460435173859D+00,0.9929034086034010D+00/
        
        data w4/
     1  0.8330585410015175D-09,0.7385431434527547D-03,
     2  0.1760634433364155D-03,0.2560148443061556D-02,
     3  0.3451125062741129D-04,0.5776070597623021D-08,
     4  0.7193475273438003D-02,0.5637184633459150D-05,
     5  0.1639382386298983D-01,0.1061058127803228D-09,
     6  0.1768374502927206D-06,0.4884885703062911D-07,
     7  0.1283148925042053D-08,0.3094771446406178D-01,
     8  0.6051334439324155D-08,0.1927387830112129D-07,
     9  0.4976440945745604D-01,0.4233758113915214D-08,
     *  0.1025233298757542D+00,0.8856070754141295D-01,
     1  0.1104400635855852D+00,0.5197414659500582D-08,
     2  0.1118799588261987D+00,0.1006268752286475D-07,
     3  0.1071314629344889D+00,0.7007502077259993D-01,
     4  0.9685583196690482D-01,0.4609504399279335D-08,
     5  0.8188245144525946D-01,0.9020406987603631D-06,
     6  0.6486594106077442D-08,0.6312205274815898D-01,
     7  0.4489461042901864D-08,0.4988262919092849D-08,
     8  0.4154002064729297D-01,0.1817357201500115D-01/

        data r5/
     1  0.6521719142155402D-07,0.1202549750969990D-07,
     2  0.7566313400414526D-07,0.1465495449601864D-07,
     3  0.8599929943723271D-12,0.5021323600578882D-10,
     4  0.1994209296313935D-07,0.1164183440133353D-06,
     5  0.6984312471450402D-08,0.6386880890219265D-09,
     6  0.4855839306362314D-07,0.3623146561750795D-07,
     7  0.2816942193302889D-08,0.2724848623593784D-07,
     8  0.2362523834427984D-06,0.3189054633609749D-05,
     9  0.2059552780710950D-04,0.6934595514690132D-06,
     *  0.1407132248509217D-03,0.8007976844767223D-03,
     1  0.3473836347641625D-02,0.1153135012591734D-01,
     2  0.3033524597048090D-01,0.6581519484436763D-01,
     3  0.1221299671914139D+00,0.1999218325295823D+00,
     4  0.2961121738752803D+00,0.4049749684603796D+00,
     5  0.5195957687151465D+00,0.6330213730596652D+00,
     6  0.7389060039762100D+00,0.8317857785505407D+00,
     7  0.9071738197563792D+00,0.9615991751734732D+00,
     *  0.9926451476867222D+00/

        data w5/
     1  0.1359868525898122D-07,0.4219099830368778D-08,
     2  0.1867624798774938D-07,0.2941398660958637D-08,
     3  0.2976889075589546D-11,0.1687953957380264D-09,
     4  0.6668414156656341D-08,0.6591361710067259D-07,
     5  0.4951806808809796D-08,0.1212548655054726D-08,
     6  0.1490574857646397D-07,0.1025683881287604D-07,
     7  0.3221663456305613D-08,0.7981674966022933D-08,
     8  0.2056314207631640D-06,0.5529107026348647D-05,
     9  0.4005329833730193D-04,0.9012637647123848D-06,
     *  0.2613043801998679D-03,0.1287270674416889D-02,
     1  0.4616768345439582D-02,0.1241728475948083D-01,
     2  0.2622213273642515D-01,0.4546252178632747D-01,
     3  0.6727651271051805D-01,0.8777370940544613D-01,
     4  0.1036265588833903D+00,0.1129271041805682D+00,
     5  0.1151452952965909D+00,0.1106442148089739D+00,
     6  0.1002131017012666D+00,0.8479808795010499D-01,
     7  0.6539896995388908D-01,0.4304754160780158D-01,
     *  0.1883477679909480D-01/

        data r6/
     1  0.6304804513609647D-11,0.1426421750302470D-09,
     2  0.7396106438319135D-09,0.1542864720427925D-08,
     3  0.3853823183581938D-08,0.8704572081149007D-08,
     4  0.1673219360776321D-07,0.2947844715472169D-07,
     5  0.5106964530201361D-07,0.9316648419785648D-07,
     6  0.1931807298032595D-06,0.4985804691837017D-06,
     7  0.1750739997791942D-05,0.8492304092302991D-05,
     8  0.4965922991610142D-04,0.2815265329835350D-03,
     9  0.1329715169111711D-02,0.4971331001040423D-02,
     *  0.1479902622577472D-01,0.3603663382760395D-01,
     1  0.7407111266405468D-01,0.1323876946570983D+00,
     2  0.2111696484525078D+00,0.3072462267349789D+00,
     3  0.4150918221548906D+00,0.5281177525001545D+00,
     4  0.6396889697070616D+00,0.7437117987177007D+00,
     5  0.8349030254224299D+00,0.9088998305188319D+00,
     6  0.9623143280232493D+00,0.9927822082244512D+00/

        data w6/
     1  0.2697515949405009D-10,0.3230982811701180D-09,
     2  0.7576582748077717D-09,0.1198290379858687D-08,
     3  0.3512840726840224D-08,0.6288416520486111D-08,
     4  0.9998867017441425D-08,0.1611500086826658D-07,
     5  0.2878433732981644D-07,0.6068368416077373D-07,
     6  0.1588779834454017D-06,0.5434103305093806D-06,
     7  0.2497570466494985D-05,0.1450137517158607D-04,
     8  0.8863882715595172D-04,0.4672225226880401D-03,
     9  0.1911845995113538D-02,0.5971782656166438D-02,
     *  0.1458064427730375D-01,0.2883040708846649D-01,
     1  0.4784454931551144D-01,0.6881903447685163D-01,
     2  0.8819765281487984D-01,0.1030116112957522D+00,
     3  0.1115642419955522D+00,0.1133723423511203D+00,
     4  0.1087504129645590D+00,0.9841311960374674D-01,
     5  0.8324137944604093D-01,0.6418687950260878D-01,
     6  0.4224657571272793D-01,0.1848383023063403D-01/
c        
c         %%%%%%%%%%%%%%%%%%%%%

ccc        call prin2('alpha = *',alpha,1)
        dtwo = 2
        if (alpha .ge. 1/dtwo) then
        do 1000 i=1,nqs(6)
        roots(i) =r6(i)
        weights(i) = w6(i)
 1000 continue

        nquads = nqs(6)
        endif

        if (alpha .ge. dtwo**(-2)) then
        if (alpha .lt. dtwo**(-1)) then

        do 2000 i=1,nqs(5)
        roots(i) =r5(i)
        weights(i) = w5(i)
 2000 continue

        nquads = nqs(5)

        endif
        endif


        if (alpha .ge. dtwo**(-3)) then
        if (alpha .lt. dtwo**(-2)) then

        do 3000 i=1,nqs(4)
        roots(i) =r4(i)
        weights(i) = w4(i)
 3000 continue

        nquads = nqs(4)

        endif
        endif


        if (alpha .ge. dtwo**(-4)) then
        if (alpha .lt. dtwo**(-3)) then

        do 4000 i=1,nqs(3)
        roots(i) =r3(i)
        weights(i) = w3(i)
 4000 continue

        nquads = nqs(3)

        endif
        endif


        if (alpha .ge. dtwo**(-5)) then
        if (alpha .lt. dtwo**(-4)) then

        do 5000 i=1,nqs(2)
        roots(i) =r2(i)
        weights(i) = w2(i)
 5000 continue

        nquads = nqs(2)

        endif
        endif


        if (alpha .lt. dtwo**(-5)) then

        do 6000 i=1,nqs(1)
        roots(i) =r1(i)
        weights(i) = w1(i)
 6000 continue

        nquads = nqs(1)

        endif

        return
        end
c
c
c
c
c
        subroutine lap_load_quads3(alpha,roots,weights,nquads)
        implicit real *8 (a-h,o-z)
        dimension r1(35),r2(35),r3(36),r4(36),r5(34),r6(32),
     1          w1(35),w2(35),w3(36),w4(36),w5(34),w6(32),nqs(6),
     2          roots(1),weights(1)

        data nqs/35,35,36,36,34,32/


        data r1/
     1  0.1805487218555334D-09,0.1212275228359914D-07,
     2  0.1284227126071698D-06,0.1830453848159829D-06,
     3  0.3772979943105942D-06,0.4462086024980011D-06,
     4  0.4740560178738865D-06,0.4900260541125241D-06,
     5  0.5013459758897048D-06,0.5108049427410262D-06,
     6  0.5200609694434561D-06,0.5306983352996112D-06,
     7  0.5451522125795341D-06,0.5692201249626951D-06,
     8  0.6251592829426325D-06,0.8987511938212741D-06,
     9  0.3886849134643068D-05,0.2875492219854735D-04,
     *  0.1939544896166486D-03,0.1029882581459594D-02,
     1  0.4169724158466954D-02,0.1311116402572781D-01,
     2  0.3315680088161359D-01,0.6995874616695880D-01,
     3  0.1273216449067481D+00,0.2056424183842035D+00,
     4  0.3017899169697274D+00,0.4101406109772463D+00,
     5  0.5239491979247005D+00,0.6364277346350874D+00,
     6  0.7413609177898303D+00,0.8333778666681363D+00,
     7  0.9080552100484312D+00,0.9619643238841511D+00,
     *  0.9927151241925502D+00/

        data w1/
     1  0.1057802991305544D-08,0.3807992661897960D-07,
     2  0.1414804981003144D-06,0.1313221981941030D-06,
     3  0.1185222444011825D-06,0.3940115046304848D-07,
     4  0.1996881646681344D-07,0.1299679085531226D-07,
     5  0.1006773670808307D-07,0.9112185716012771D-08,
     6  0.9651232529857960D-08,0.1199828181896397D-07,
     7  0.1774486487685703D-07,0.3318051155839965D-07,
     8  0.9511015206828043D-07,0.6935217714501879D-06,
     9  0.7503592418367999D-05,0.5735545402032385D-04,
     *  0.3493475112135047D-03,0.1579427082833151D-02,
     1  0.5282306245594998D-02,0.1351042928701259D-01,
     2  0.2756343683964651D-01,0.4670269042848771D-01,
     3  0.6809028209363153D-01,0.8800984719845587D-01,
     4  0.1033234700141851D+00,0.1122354008280486D+00,
     5  0.1142405900328938D+00,0.1096765286482062D+00,
     6  0.9929323611327934D-01,0.8400274217161918D-01,
     7  0.6477988428536580D-01,0.4263854160553587D-01,
     *  0.1865559735138674D-01/

        data r2/
     1  0.2498399244559789D-09,0.9574469091217448D-08,
     2  0.8410201161999000D-07,0.2464787381837469D-06,
     3  0.3632970651921083D-06,0.4225419319046052D-06,
     4  0.4569632778101988D-06,0.4807535490847752D-06,
     5  0.4998371453926280D-06,0.5175162267615506D-06,
     6  0.5365386529241978D-06,0.5604461089917226D-06,
     7  0.5960466238664378D-06,0.6621789129867101D-06,
     8  0.8364172976551161D-06,0.1596673765445400D-05,
     9  0.6321752802582027D-05,0.3765174448792888D-04,
     *  0.2249889243090204D-03,0.1119649028058666D-02,
     1  0.4377105162892054D-02,0.1349802636345287D-01,
     2  0.3375513432092290D-01,0.7074856109815952D-01,
     3  0.1282358388463982D+00,0.2065922431824465D+00,
     4  0.3026926772135806D+00,0.4109363753194458D+00,
     5  0.5246046470010458D+00,0.6369322031151307D+00,
     6  0.7417202178008668D+00,0.8336089514701711D+00,
     7  0.9081824117422336D+00,0.9620168300177505D+00,
     *  0.9927251670119255D+00/

        data w2/
     1  0.1204717319299347D-08,0.2664649243706006D-07,
     2  0.1339902911257218D-06,0.1561710619685148D-06,
     3  0.8091702776495329D-07,0.4325071621166525D-07,
     4  0.2772249778319588D-07,0.2076670193512961D-07,
     5  0.1792190322656422D-07,0.1787504779338550D-07,
     6  0.2072680614269314D-07,0.2811683476816226D-07,
     7  0.4574241010145619D-07,0.9630335393096237D-07,
     8  0.3066549782686194D-06,0.1605736878231447D-05,
     9  0.1054903458655302D-04,0.6898026732574348D-04,
     *  0.3856187905172795D-03,0.1664736853446596D-02,
     1  0.5433003104488688D-02,0.1371351003297060D-01,
     2  0.2777411869003385D-01,0.4686674547687584D-01,
     3  0.6817135862860695D-01,0.8800125603758564D-01,
     4  0.1032419573618542D+00,0.1121074865717644D+00,
     5  0.1140916789586455D+00,0.1095262049173790D+00,
     6  0.9915505123322215D-01,0.8388570560574137D-01,
     7  0.6468995230580873D-01,0.4257958575385162D-01,
     *  0.1862987062757630D-01/

        data r3/
     1  0.1086716818114864D-06,0.1200054648972930D-03,
     2  0.5252613985635641D-03,0.5618119315158980D-07,
     3  0.2537485563265096D-04,0.2013159822883746D-02,
     4  0.6521951154530152D-02,0.5624720594412130D-05,
     5  0.1768163846399866D-01,0.8511775731581190D-08,
     6  0.2111853668893954D-06,0.4052117425967484D-01,
     7  0.2867548817579140D-09,0.1786350674163324D-05,
     8  0.8005430210817712D-01,0.1394007287619701D+00,
     9  0.9887075050053395D-06,0.2185495321179273D+00,
     *  0.3126927123317131D-06,0.3143446137058399D+00,
     1  0.7516055555054246D-06,0.4214152527778254D+00,
     2  0.3812738627393588D-06,0.5333727801546942D+00,
     3  0.6461153181402344D-06,0.6437628334017034D+00,
     4  0.4301636856112372D-06,0.7466297681641056D+00,
     5  0.5844082366391510D-06,0.8367877099308710D+00,
     6  0.4691434464891472D-06,0.9099404209007163D+00,
     7  0.5406866183632631D-06,0.9627447194585455D+00,
     8  0.5043033594624691D-06,0.9928646171514977D+00/

        data w3/
     1  0.6209018672270703D-07,0.1831718162587153D-03,
     2  0.7446276734537522D-03,0.6511940069669588D-07,
     3  0.3967554880259226D-04,0.2543369161841548D-02,
     4  0.7083059161168087D-02,0.7842184610701444D-05,
     5  0.1609987496547276D-01,0.2174785027421426D-07,
     6  0.1179751967292685D-06,0.3044926797097564D-01,
     7  0.1292977954992012D-08,0.1524624337689102D-05,
     8  0.4915906200264103D-01,0.6953130644145892D-01,
     9  0.3808612454497425D-06,0.8822289248976031D-01,
     *  0.8313863351400470D-07,0.1024509214147290D+00,
     1  0.1458721186075257D-06,0.1106094421281519D+00,
     2  0.5663729595262236D-07,0.1122207904626886D+00,
     3  0.7692195218469183D-07,0.1075625648055084D+00,
     4  0.4271752468773728D-07,0.9730521689902013D-01,
     5  0.5026007264017589D-07,0.8229365757040197D-01,
     6  0.3621349771636688D-07,0.6345377343518570D-01,
     7  0.3882207473511349D-07,0.4176396928328413D-01,
     8  0.3491561915192619D-07,0.1827277537460141D-01/

        data r4/
     1  0.2290010088858140D-06,0.2376809971718999D-06,
     2  0.5966472821662121D-10,0.1833423098155487D-08,
     3  0.1430163037683537D-06,0.2330093014939201D-01,
     4  0.8761105298967350D-02,0.5185093237126674D-01,
     5  0.2698069142288247D-02,0.9915496760113392D-01,
     6  0.1672515098908588D+00,0.6047628148244620D-07,
     7  0.2546671755136946D+00,0.3198149584665396D-06,
     8  0.6793908992306910D-03,0.3569408151826680D+00,
     9  0.1513700312790496D-07,0.4678664845970772D+00,
     *  0.5807133218570147D+00,0.6890441202509378D+00,
     1  0.7871334706437450D+00,0.1463616950595870D-03,
     2  0.8701456648028855D+00,0.9342224898742861D+00,
     3  0.9766189765085179D+00,0.7522306601221303D-05,
     4  0.2647562020214068D-05,0.9965784235307057D+00,
     5  0.1382248427457567D-05,0.5315521027775555D-06,
     6  0.9417170692275908D-06,0.3053829832256961D-04,
     7  0.3928259133722243D-06,0.4608200513596277D-06,
     8  0.7367039721487189D-06,0.6164916613269864D-06/

        data w4/
     1  0.1657143029753335D-07,0.7340156567664806D-07,
     2  0.2727661212764934D-09,0.4777680392037202D-08,
     3  0.9292674424370482D-07,0.2061490507894380D-01,
     4  0.9421626785200346D-02,0.3731084155756543D-01,
     5  0.3456019763058139D-02,0.5765085736563156D-01,
     6  0.7828611729676773D-01,0.6645938577468285D-07,
     7  0.9577899432650368D-01,0.7756435292763212D-07,
     8  0.9989501820691382D-03,0.1077127493372904D+00,
     9  0.2572162796885087D-07,0.1130032493409744D+00,
     *  0.1116110417942550D+00,0.1040937180635673D+00,
     1  0.9127780631650564D-01,0.2311396165896312D-03,
     2  0.7410197098902930D-01,0.5358690512903867D-01,
     3  0.3100979191995219D-01,0.9379257216074231D-05,
     4  0.2209694292824870D-05,0.9793510975172180D-02,
     5  0.6783447365016343D-06,0.7542876619075980D-07,
     6  0.2811744105972579D-06,0.4643752612395183D-04,
     7  0.6939325638055487D-07,0.6786358211909448D-07,
     8  0.1500720903041710D-06,0.9771185714081539D-07/

        data r5/
     1  0.1291004426496013D-09,0.2654127081078098D-08,
     2  0.1361063656420815D-07,0.3508007747291853D-07,
     3  0.7868473082564474D-07,0.1483889388784878D-06,
     4  0.2337043982959294D-06,0.3272331721538856D-06,
     5  0.4283629492829166D-06,0.5449807013489431D-06,
     6  0.7000528118749445D-06,0.9432115599483349D-06,
     7  0.1392935681837077D-05,0.2405166130854322D-05,
     8  0.5300505867709876D-05,0.1598378323777745D-04,
     9  0.6390503271039933D-04,0.2902255381619946D-03,
     *  0.1249899246512469D-02,0.4569179659285092D-02,
     1  0.1368813916856916D-01,0.3383498638184307D-01,
     2  0.7062590556044459D-01,0.1278857878373189D+00,
     3  0.2060617427676864D+00,0.3020683929169170D+00,
     4  0.4103071603359495D+00,0.5240385006037359D+00,
     5  0.6364695864852530D+00,0.7413768952891221D+00,
     6  0.8333818140205832D+00,0.9080549819835483D+00,
     7  0.9619635944151491D+00,0.9927149203599683D+00/

        data w5/
     1  0.5367928426889987D-09,0.5806856063565804D-08,
     2  0.1592237263980103D-07,0.2968873207853766D-07,
     3  0.5804285700762720D-07,0.7929000200528030D-07,
     4  0.9008815572257575D-07,0.9685858922955332D-07,
     5  0.1065859615807331D-06,0.1303405874483019D-06,
     6  0.1875405587714759D-06,0.3163670193085985D-06,
     7  0.6337627943918141D-06,0.1569745271027755D-05,
     8  0.4987927994341353D-05,0.2016937971071181D-04,
     9  0.9470455942861733D-04,0.4389994219560997D-03,
     *  0.1736572234910641D-02,0.5473766989599823D-02,
     1  0.1366050344401106D-01,0.2760949705197200D-01,
     2  0.4663850260022482D-01,0.6795744097995343D-01,
     3  0.8786094466114486D-01,0.1031947746616217D+00,
     4  0.1121410401503538D+00,0.1141795177909663D+00,
     5  0.1096412689582469D+00,0.9927544104676593D-01,
     6  0.8399547954535481D-01,0.6477810459886066D-01,
     7  0.4263888189618991D-01,0.1865608152418344D-01/

        data r6/
     1  0.3953214119924202D-03,0.1003981549953298D-03,
     2  0.4640704735242578D-08,0.1511510995107804D-02,
     3  0.2782481227627916D-04,0.1718072286607394D-07,
     4  0.5109428804130664D-02,0.9192207082426758D-05,
     5  0.1460831537115056D-01,0.3514767128471333D-01,
     6  0.4647936578562926D-07,0.3726089353983933D-05,
     7  0.6073683374398020D-09,0.7223334012281243D-01,
     8  0.1296208814325718D+00,0.1036849792697204D-06,
     9  0.2077538099360583D+00,0.1805527421766401D-05,
     *  0.3035898769732535D+00,0.5834967593749685D-06,
     1  0.4115871249263999D+00,0.1555194020013263D-10,
     2  0.1992415036366161D-06,0.5250535146727502D+00,
     3  0.3479143307714730D-06,0.9919305366025207D-06,
     4  0.6372277945409344D+00,0.7419047008144393D+00,
     5  0.8337155530359967D+00,0.9082364842682121D+00,
     6  0.9620379245298399D+00,0.9927290773157608D+00/

        data w6/
     1  0.5435882689032238D-03,0.1350937926275206D-03,
     2  0.7305889850068755D-08,0.1951632905994262D-02,
     3  0.3354599729166260D-04,0.1912644509730883D-07,
     4  0.5812813687175159D-02,0.9236221529037789D-05,
     5  0.1406561809211463D-01,0.2796980572397013D-01,
     6  0.4144404301402482D-07,0.3003701230609629D-05,
     7  0.1591810269529538D-08,0.4685572347868969D-01,
     8  0.6799534650330006D-01,0.7460963256177968D-07,
     9  0.8774479234564085D-01,0.1177385189705925D-05,
     *  0.1029795814321038D+00,0.3004393145810986D-06,
     1  0.1118812055449611D+00,0.8080470437601588D-10,
     2  0.1187352136950522D-06,0.1139146521851929D+00,
     3  0.1837679132610744D-06,0.5509464789088935D-06,
     4  0.1093954748580578D+00,0.9906202009062036D-01,
     5  0.8382164834754727D-01,0.6464798675139636D-01,
     6  0.4255497463069616D-01,0.1861978000822164D-01/
c        
c         %%%%%%%%%%%%%%%%%%%%%

ccc        call prin2('alpha = *',alpha,1)
        dtwo = 2
        if (alpha .ge. 1/dtwo) then
        do 1000 i=1,nqs(6)
        roots(i) =r6(i)
        weights(i) = w6(i)
 1000 continue

        nquads = nqs(6)
        endif

        if (alpha .ge. dtwo**(-2)) then
        if (alpha .lt. dtwo**(-1)) then

        do 2000 i=1,nqs(5)
        roots(i) =r5(i)
        weights(i) = w5(i)
 2000 continue

        nquads = nqs(5)

        endif
        endif


        if (alpha .ge. dtwo**(-3)) then
        if (alpha .lt. dtwo**(-2)) then

        do 3000 i=1,nqs(4)
        roots(i) =r4(i)
        weights(i) = w4(i)
 3000 continue

        nquads = nqs(4)

        endif
        endif


        if (alpha .ge. dtwo**(-4)) then
        if (alpha .lt. dtwo**(-3)) then

        do 4000 i=1,nqs(3)
        roots(i) =r3(i)
        weights(i) = w3(i)
 4000 continue

        nquads = nqs(3)

        endif
        endif


        if (alpha .ge. dtwo**(-5)) then
        if (alpha .lt. dtwo**(-4)) then

        do 5000 i=1,nqs(2)
        roots(i) =r2(i)
        weights(i) = w2(i)
 5000 continue

        nquads = nqs(2)

        endif
        endif


        if (alpha .lt. dtwo**(-5)) then

        do 6000 i=1,nqs(1)
        roots(i) =r1(i)
        weights(i) = w1(i)
 6000 continue

        nquads = nqs(1)

        endif

        return
        end
c
c
c
c
c
        subroutine lap_load_quads4(alpha,roots,weights,nquads)
        implicit real *8 (a-h,o-z)
        dimension r1(35),r2(35),r3(35),r4(35),r5(35),r6(31),
     1          w1(35),w2(35),w3(35),w4(35),w5(35),w6(31),nqs(6),
     2          roots(1),weights(1)

        data nqs/35,35,35,35,35,31/

        data r1/
     1  0.1059993813945641D-08,0.3916818682885566D-07,
     2  0.3971613082630646D-06,0.1756312324159999D-05,
     3  0.3226314820110193D-05,0.3755034685929197D-05,
     4  0.3980347275630697D-05,0.4111451179505097D-05,
     5  0.4204920932585694D-05,0.4283238813508983D-05,
     6  0.4359905988861207D-05,0.4447832190807734D-05,
     7  0.4566733122859630D-05,0.4762656865798449D-05,
     8  0.5203235766352740D-05,0.7013331546928513D-05,
     9  0.1942313624963768D-04,0.8766602824008714D-04,
     *  0.4105822617506400D-03,0.1694771599514038D-02,
     1  0.5806939617604720D-02,0.1636284314937582D-01,
     2  0.3848409218597422D-01,0.7735600154229800D-01,
     3  0.1362606819108401D+00,0.2152675687276333D+00,
     4  0.3112067205000019D+00,0.4186339427231879D+00,
     5  0.5310706918137091D+00,0.6419837555163668D+00,
     6  0.7453584883965373D+00,0.8359680349187027D+00,
     7  0.9094884117630064D+00,0.9625579175348178D+00,
     *  0.9928288647821826D+00/

        data w1/
     1  0.4985557074462531D-08,0.1104658460793608D-06,
     2  0.7459275318607946D-06,0.1867098228032338D-05,
     3  0.8686519971172985D-06,0.3149639569758498D-06,
     4  0.1632313045522593D-06,0.1070953187284464D-06,
     5  0.8327944092206894D-07,0.7549174851093128D-07,
     6  0.7989211556249258D-07,0.9900129464393497D-07,
     7  0.1454820418505127D-06,0.2679277540520277D-06,
     8  0.7280839220986226D-06,0.4023351863867031D-05,
     9  0.2676102066683701D-04,0.1373337509275373D-03,
     *  0.6145796465121726D-03,0.2252853566186016D-02,
     1  0.6576923161504502D-02,0.1541441100484985D-01,
     2  0.2972844812868183D-01,0.4858584588668965D-01,
     3  0.6923457169266309D-01,0.8823280310729238D-01,
     4  0.1027161249971412D+00,0.1110419680056961D+00,
     5  0.1127334430462473D+00,0.1080861495626856D+00,
     6  0.9779048608483493D-01,0.8270704158299527D-01,
     7  0.6377267391013965D-01,0.4197354815774579D-01,
     *  0.1836434875661838D-01/

        data r2/
     1  0.7446273523878573D-09,0.3369049352122097D-07,
     2  0.3657404275255330D-06,0.8973354553309214D-06,
     3  0.1914140440584365D-05,0.2997410013193958D-05,
     4  0.3525687266246351D-05,0.3823726220971193D-05,
     5  0.4027234941708568D-05,0.4189611308169932D-05,
     6  0.4339765834001274D-05,0.4501531146613825D-05,
     7  0.4705797502150559D-05,0.5013184303784516D-05,
     8  0.5599223558252155D-05,0.7268119995999740D-05,
     9  0.1584076634837391D-04,0.6727528965393195D-04,
     *  0.3284666313773889D-03,0.1426249609362276D-02,
     1  0.5106054456741301D-02,0.1489685324952129D-01,
     2  0.3597467023981040D-01,0.7374395861871383D-01,
     3  0.1317681103779638D+00,0.2103189419680338D+00,
     4  0.3062789360900211D+00,0.4141288338434832D+00,
     5  0.5272543657006729D+00,0.6389834887715631D+00,
     6  0.7431875946120303D+00,0.8345557034983718D+00,
     7  0.9087047187706956D+00,0.9622327390215185D+00,
     *  0.9927664956194669D+00/

        data w2/
     1  0.3695715678409841D-08,0.1002750102415574D-06,
     2  0.6406046006333814D-06,0.4213582875582112D-06,
     3  0.1412515663858830D-05,0.7369436213276440D-06,
     4  0.3779603384669745D-06,0.2381141664977551D-06,
     5  0.1770282536667355D-06,0.1522757716680802D-06,
     6  0.1518341229360948D-06,0.1765436263454500D-06,
     7  0.2410931274051740D-06,0.3980378213492158D-06,
     8  0.8705153998795342D-06,0.3100995947234053D-05,
     9  0.1882174600634935D-04,0.1067561898920751D-03,
     *  0.5084736596605367D-03,0.1964576869561904D-02,
     1  0.5983875858533123D-02,0.1448579207431350D-01,
     2  0.2860870807777073D-01,0.4754988455020094D-01,
     3  0.6854443349977631D-01,0.8801904201465953D-01,
     4  0.1029578034285470D+00,0.1116220176042343D+00,
     5  0.1135073798319795D+00,0.1089252905562921D+00,
     6  0.9859617075414107D-01,0.8340873689039988D-01,
     7  0.6432162751116040D-01,0.4233739125521081D-01,
     *  0.1852401783618523D-01/

        data r3/
     1  0.7169935857529262D-10,0.1620063862577103D-07,
     2  0.2273360966836847D-06,0.7941143495813348D-06,
     3  0.1299287883991965D-05,0.2244174016467081D-05,
     4  0.2945773313091079D-05,0.3422827825009347D-05,
     5  0.3781081265280639D-05,0.4084187100506913D-05,
     6  0.4376376462879391D-05,0.4701594013461139D-05,
     7  0.5123964571982582D-05,0.5775761024940094D-05,
     8  0.7036479081324257D-05,0.1046989607014254D-04,
     9  0.2487794844063623D-04,0.9955320282891636D-04,
     *  0.4645682675514940D-03,0.1955862316047836D-02,
     1  0.6776664487091988D-02,0.1906662264459815D-01,
     2  0.4436111968994979D-01,0.8774909413617102D-01,
     3  0.1517834067552197D+00,0.2354539247157524D+00,
     4  0.3346070674786194D+00,0.4432175741147756D+00,
     5  0.5546944689886045D+00,0.6627714901718117D+00,
     6  0.7619564853849661D+00,0.8477021984549095D+00,
     7  0.9164455120411763D+00,0.9655852471940818D+00,
     *  0.9934254856416695D+00/

        data w3/
     1  0.7983579719877429D-09,0.5561543585064123D-07,
     2  0.4468457616457022D-06,0.4217028831065741D-06,
     3  0.8638331947767411D-06,0.8547484082455534D-06,
     4  0.5677210861773482D-06,0.4041999959303056D-06,
     5  0.3224403194717018D-06,0.2908515484368089D-06,
     6  0.3005578265882897D-06,0.3599240180354502D-06,
     7  0.5047012410151333D-06,0.8522554028487949D-06,
     8  0.1868379848449827D-05,0.6061445648708709D-05,
     9  0.2913305953324194D-04,0.1518590374671077D-03,
     *  0.7035350548439830D-03,0.2633976027379589D-02,
     1  0.7705816051043109D-02,0.1783041934909292D-01,
     2  0.3364527117178830D-01,0.5357415424272354D-01,
     3  0.7431731828881288D-01,0.9230805164975573D-01,
     4  0.1049711921624970D+00,0.1111365454496536D+00,
     5  0.1107658506901630D+00,0.1044703549515382D+00,
     6  0.9314125532269521D-01,0.7775671880967020D-01,
     7  0.5929744138616287D-01,0.3870013300138509D-01,
     *  0.1684679827281716D-01/

        data r4/
     1  0.5960204081759053D-09,0.1600079946222516D-07,
     2  0.1072583256806995D-06,0.3363643078199049D-06,
     3  0.8035802562452548D-06,0.1490236229591095D-05,
     4  0.2191079789374611D-05,0.2814607038589473D-05,
     5  0.3369860713563083D-05,0.3895875628798055D-05,
     6  0.4445215061743258D-05,0.5093739165846994D-05,
     7  0.5973256307165460D-05,0.7365740514694280D-05,
     8  0.1002199698375403D-04,0.1643052002749496D-04,
     9  0.3640777334023690D-04,0.1116024989548059D-03,
     *  0.4127948871046316D-03,0.1543183459943612D-02,
     1  0.5167140523929174D-02,0.1470283472212973D-01,
     2  0.3528224359666048D-01,0.7240119776054615D-01,
     3  0.1298079335125832D+00,0.2079436281176629D+00,
     4  0.3037679550329249D+00,0.4117431422761551D+00,
     5  0.5251817568604317D+00,0.6373264957287603D+00,
     6  0.7419750593933174D+00,0.8337608529977570D+00,
     7  0.9082614441391893D+00,0.9620482349195940D+00,
     *  0.9927310502014603D+00/

        data w4/
     1  0.2632584237733635D-08,0.3934128358909442D-07,
     2  0.1529478907534533D-06,0.3247034834902775D-06,
     3  0.6088620459840252D-06,0.7225729107323680D-06,
     4  0.6656046422538303D-06,0.5841555270660910D-06,
     5  0.5330349410653626D-06,0.5275876075388303D-06,
     6  0.5832616735672008D-06,0.7344609971711408D-06,
     7  0.1067430217565926D-05,0.1827041159237842D-05,
     8  0.3832353038331015D-05,0.1031538768279799D-04,
     9  0.3528828035077502D-04,0.1389743765269133D-03,
     *  0.5529962759582973D-03,0.1971369891813471D-02,
     1  0.5845008248491401D-02,0.1410570463018419D-01,
     2  0.2800805861380143D-01,0.4688278360414000D-01,
     3  0.6800630288545189D-01,0.8773972523340606D-01,
     4  0.1029618979896856D+00,0.1118555530571275D+00,
     5  0.1138854169987014D+00,0.1093661195260665D+00,
     6  0.9903503365164804D-01,0.8379875344475088D-01,
     7  0.6463036165498072D-01,0.4254340369929523D-01,
     *  0.1861472655993423D-01/

        data r5/
     1  0.2262430269781617D-02,0.9386862444230902D-11,
     2  0.8588271364128832D-03,0.5867085753247661D-02,
     3  0.1483276534428340D-08,0.2910708260177847D-03,
     4  0.1965307032708896D-07,0.1501542181355715D-01,
     5  0.1254635301906868D-05,0.9909871178739790D-07,
     6  0.9845867964596219D-04,0.3482518710805700D-01,
     7  0.7092393967608478D-01,0.1941584342148552D-04,
     8  0.3032300071731556D-06,0.1167723606007495D-04,
     9  0.1938713077724199D-05,0.1273390800257244D+00,
     *  0.6901031419641947D-06,0.2047776289151140D+00,
     1  0.3911914065696030D-04,0.3003259838098989D+00,
     2  0.8008446765847068D-05,0.4580469825547848D-05,
     3  0.2702693036761467D-05,0.4084236361094997D+00,
     4  0.5939967498718950D-05,0.5222739789619451D+00,
     5  0.6349910159873866D+00,0.3558353638119591D-05,
     6  0.7402615310378139D+00,0.8326358964481178D+00,
     7  0.9076335348856710D+00,0.9617867558396507D+00,
     *  0.9926808054990950D+00/

        data w5/
     1  0.2127210894048255D-02,0.8313371176425357D-10,
     2  0.8795898002926717D-03,0.5647153337238339D-02,
     3  0.4952675234125894D-08,0.3245172704639595D-03,
     4  0.3921275573204346D-07,0.1352208115382500D-01,
     5  0.6343994758745780D-06,0.1303095154746884D-06,
     6  0.1007891352335563D-03,0.2707467241225418D-01,
     7  0.4582553825691528D-01,0.1155041608892026D-04,
     8  0.2891570445688625D-06,0.5053321850438725D-05,
     9  0.7267282450802848D-06,0.6712694395775582D-01,
     *  0.4834425532811359D-06,0.8724593595906486D-01,
     1  0.3178348501393442D-04,0.1028987116270952D+00,
     2  0.2645931589079479D-05,0.1151015719594823D-05,
     3  0.8030979410076135D-06,0.1121442829599567D+00,
     4  0.1625522974730212D-05,0.1143980911657633D+00,
     5  0.1099794267257950D+00,0.9202226650279670D-06,
     6  0.9965186491599589D-01,0.8434948052566986D-01,
     7  0.6506729060776358D-01,0.4283527625705690D-01,
     *  0.1874330173857002D-01/

        data r6/
     1  0.2978624732214353D-09,0.6720053212458587D-08,
     2  0.4561734003435004D-07,0.1696624602648334D-06,
     3  0.4457328611405440D-06,0.9447646188655679D-06,
     4  0.1744178663321435D-05,0.2963765755694934D-05,
     5  0.4864216435204008D-05,0.8072969180505095D-05,
     6  0.1418484092977766D-04,0.2763965996312870D-04,
     7  0.6232984521314316D-04,0.1664765753393147D-03,
     8  0.5133642874713985D-03,0.1677591744165506D-02,
     9  0.5228732647177122D-02,0.1444972342950120D-01,
     *  0.3441182252322664D-01,0.7071805944162654D-01,
     1  0.1273477535197926D+00,0.2049590453357277D+00,
     2  0.3006121958171207D+00,0.4087465333241848D+00,
     3  0.5225810695454728D+00,0.6352497738761548D+00,
     4  0.7404570950515785D+00,0.8327667418687377D+00,
     5  0.9077074579715845D+00,0.9618177677966631D+00,
     *  0.9926867872914786D+00/

        data w6/
     1  0.1246600832056045D-08,0.1585547411844744D-07,
     2  0.7100859521251142D-07,0.1883885202257766D-06,
     3  0.3756015689871911D-06,0.6348278386100083D-06,
     4  0.9824252040760656D-06,0.1496891637562176D-05,
     5  0.2400194920347445D-05,0.4257849851529539D-05,
     6  0.8619684897519731D-05,0.2028325370550476D-04,
     7  0.5590621129629236D-04,0.1771276762849948D-03,
     8  0.6022642815665083D-03,0.1973690213537773D-02,
     9  0.5672259315143491D-02,0.1363699803697145D-01,
     *  0.2726090038639708D-01,0.4604571299308095D-01,
     1  0.6732729860171451D-01,0.8738668341114303D-01,
     2  0.1029677315666777D+00,0.1121514407508686D+00,
     3  0.1143627561037299D+00,0.1099210678627585D+00,
     4  0.9958604934993043D-01,0.8428739571004177D-01,
     5  0.6501655361868886D-01,0.4280082894989821D-01,
     *  0.1872800773145546D-01/
c        
c        
c         %%%%%%%%%%%%%%%%%%%%%

ccc        call prin2('alpha = *',alpha,1)
        dtwo = 2
        if (alpha .ge. 1/dtwo) then
        do 1000 i=1,nqs(6)
        roots(i) =r6(i)
        weights(i) = w6(i)
 1000 continue

        nquads = nqs(6)
        endif

        if (alpha .ge. dtwo**(-2)) then
        if (alpha .lt. dtwo**(-1)) then

        do 2000 i=1,nqs(5)
        roots(i) =r5(i)
        weights(i) = w5(i)
 2000 continue

        nquads = nqs(5)

        endif
        endif


        if (alpha .ge. dtwo**(-3)) then
        if (alpha .lt. dtwo**(-2)) then

        do 3000 i=1,nqs(4)
        roots(i) =r4(i)
        weights(i) = w4(i)
 3000 continue

        nquads = nqs(4)

        endif
        endif


        if (alpha .ge. dtwo**(-4)) then
        if (alpha .lt. dtwo**(-3)) then

        do 4000 i=1,nqs(3)
        roots(i) =r3(i)
        weights(i) = w3(i)
 4000 continue

        nquads = nqs(3)

        endif
        endif


        if (alpha .ge. dtwo**(-5)) then
        if (alpha .lt. dtwo**(-4)) then

        do 5000 i=1,nqs(2)
        roots(i) =r2(i)
        weights(i) = w2(i)
 5000 continue

        nquads = nqs(2)

        endif
        endif


        if (alpha .lt. dtwo**(-5)) then

        do 6000 i=1,nqs(1)
        roots(i) =r1(i)
        weights(i) = w1(i)
 6000 continue

        nquads = nqs(1)

        endif

        return
        end
c
c
c
c
c
        subroutine lap_load_quads5(alpha,roots,weights,nquads)
        implicit real *8 (a-h,o-z)
        dimension r1(35),r2(35),r3(36),r4(35),r5(34),r6(31),
     1          w1(35),w2(35),w3(36),w4(35),w5(34),w6(31),nqs(6),
     2          roots(1),weights(1)

        data nqs/35,35,36,35,34,31/


        data r1/
     1  0.4250430297776588D-08,0.1493769854414862D-06,
     2  0.1581445398005698D-05,0.7826712781719706D-05,
     3  0.1638508918911356D-04,0.2041263736481108D-04,
     4  0.2209251395025545D-04,0.2299598554518875D-04,
     5  0.2359624076641904D-04,0.2406559541522347D-04,
     6  0.2449198906539893D-04,0.2494301284139143D-04,
     7  0.2550008781290888D-04,0.2631259806747105D-04,
     8  0.2778066326188823D-04,0.3150166281571988D-04,
     9  0.4703238410677645D-04,0.1300179326108574D-03,
     *  0.4908085530911444D-03,0.1802524826426396D-02,
     1  0.5818289093172258D-02,0.1599135431495479D-01,
     2  0.3736573240297302D-01,0.7524500536181848D-01,
     3  0.1331818989652865D+00,0.2115124755490829D+00,
     4  0.3072037647487513D+00,0.4147994481587772D+00,
     5  0.5277153927624592D+00,0.6392853565909898D+00,
     6  0.7433746519378933D+00,0.8346629800179082D+00,
     7  0.9087587582264081D+00,0.9622537067437443D+00,
     *  0.9927703698411021D+00/

        data w1/
     1  0.1961270579822455D-07,0.4201344265560530D-06,
     2  0.3117935713550381D-05,0.9243479690306325D-05,
     3  0.6257667757054421D-05,0.2424841047479262D-05,
     4  0.1164596371112347D-05,0.7101536786431700D-06,
     5  0.5163443791608880D-06,0.4361560829610499D-06,
     6  0.4274136673038034D-06,0.4874710035161098D-06,
     7  0.6493454309199104D-06,0.1032486795283552D-05,
     8  0.2106870739924559D-05,0.6435160773515331D-05,
     9  0.3219120133189398D-04,0.1642067170129862D-03,
     *  0.6586808039747206D-03,0.2244522140335590D-02,
     1  0.6362742478413792D-02,0.1484771826439830D-01,
     2  0.2882198385695904D-01,0.4755517661593269D-01,
     3  0.6837638936852477D-01,0.8776140237606168D-01,
     4  0.1026886715287850D+00,0.1113875604443492D+00,
     5  0.1133234043807823D+00,0.1087897906903006D+00,
     6  0.9850043379511652D-01,0.8334348062248825D-01,
     7  0.6427934706561353D-01,0.4231283539796075D-01,
     *  0.1851401258139531D-01/

        data r2/
     1  0.4683822828198871D-08,0.1047340153614535D-07,
     2  0.2627239995592630D-06,0.2393336686308890D-05,
     3  0.8787669560003065D-05,0.1532498785341567D-04,
     4  0.1898712983322485D-04,0.2101624189225421D-04,
     5  0.2233134176884613D-04,0.2332093046082991D-04,
     6  0.2417682027172423D-04,0.2503122652664824D-04,
     7  0.2602181744210555D-04,0.2736394753043622D-04,
     8  0.2954034663591405D-04,0.3407123084861979D-04,
     9  0.4773475180114918D-04,0.1070259354726575D-03,
     *  0.3779213106579535D-03,0.1450337305054721D-02,
     1  0.4968975529534347D-02,0.1434280415000069D-01,
     2  0.3472139109086060D-01,0.7164432701503754D-01,
     3  0.1289087489374163D+00,0.2069860494031291D+00,
     4  0.3028383806480954D+00,0.4109096455735304D+00,
     5  0.5244860836486618D+00,0.6367856936251278D+00,
     6  0.7415870283701960D+00,0.8335099567786214D+00,
     7  0.9081228244126260D+00,0.9619908787770233D+00,
     *  0.9927200657965406D+00/

        data w2/
     1  0.9098768249254281D-08,0.2878689666640698D-07,
     2  0.7102919496084977D-06,0.4189045612798429D-05,
     3  0.7618017735753816D-05,0.4990023956967169D-05,
     4  0.2623706673055990D-05,0.1578771406997162D-05,
     5  0.1110794654014574D-05,0.8984084974555421D-06,
     6  0.8344738294965576D-06,0.8962685461382373D-06,
     7  0.1118329535995043D-05,0.1637359548487174D-05,
     8  0.2924450298151778D-05,0.6990936262361831D-05,
     9  0.2498413022487342D-04,0.1169913840780591D-03,
     *  0.5120620884846629D-03,0.1892786583543690D-02,
     1  0.5711094148047351D-02,0.1391881731575884D-01,
     2  0.2780100010921746D-01,0.4670640238880410D-01,
     3  0.6790324238270406D-01,0.8772626484302514D-01,
     4  0.1030280759622944D+00,0.1119769773464374D+00,
     5  0.1140354641297681D+00,0.1095226289075478D+00,
     6  0.9918183520739228D-01,0.8392470594649367D-01,
     7  0.6472794944187470D-01,0.4260770181445468D-01,
     *  0.1864285710567657D-01/

        data r3/
     1  0.2576504748745078D+00,0.3451770714731933D+00,
     2  0.1790369369448084D+00,0.4418387055595634D+00,
     3  0.1126935617589391D+00,0.5457843021692586D+00,
     4  0.6289075865308832D-01,0.6510353730657077D+00,
     5  0.3050817874096912D-01,0.7507566313272438D+00,
     6  0.1264270328889074D-01,0.4445347755992486D-02,
     7  0.1355130926594161D-02,0.3894808986947493D-03,
     8  0.8390023168560877D+00,0.2875095599931166D-05,
     9  0.1284100263687461D-03,0.7863814685361336D-06,
     *  0.9110036392274737D+00,0.6161390302143396D-04,
     1  0.7098752657502829D-05,0.1207751135705030D-04,
     2  0.2990856297760982D-08,0.9631442221890465D+00,
     3  0.4200765905692085D-04,0.9700025196336043D-07,
     4  0.1591064890587050D-04,0.3420752622071206D-04,
     5  0.9929371536846508D+00,0.1864269920617513D-04,
     6  0.3012994897767879D-04,0.2070951676802512D-04,
     7  0.2753783360648634D-04,0.2242020547002323D-04,
     8  0.2399206361857497D-04,0.2561940482256129D-04/

        data w3/
     1  0.8325959317730900D-01,0.9198762987214637D-01,
     2  0.7332049674811129D-01,0.1010522371971943D+00,
     3  0.5857436942723385D-01,0.1057661718293547D+00,
     4  0.4088481540233655D-01,0.1035727210611099D+00,
     5  0.2439872966528156D-01,0.9487663898665144D-01,
     6  0.1218123211713668D-01,0.4996081991457436D-02,
     7  0.1679043660385561D-02,0.4759895920714401D-03,
     8  0.8082079616801113D-01,0.3081580551914307D-05,
     9  0.1218170911508259D-03,0.1264770490998429D-05,
     *  0.6257960844135950D-01,0.3257633302780905D-04,
     1  0.5089329431600612D-05,0.4509389880564734D-05,
     2  0.1375215183210535D-07,0.4128625640682546D-01,
     3  0.1135134918444411D-04,0.2572657290079106D-06,
     4  0.3204496158840159D-05,0.5340054968172330D-05,
     5  0.1808469057567189D-01,0.2336401757974379D-05,
     6  0.3136496948350525D-05,0.1847445930543967D-05,
     7  0.2170481450829847D-05,0.1609009983867337D-05,
     8  0.1566121296257141D-05,0.1726310258179095D-05/

        data r4/
     1  0.1539330935972923D-08,0.5239700420079496D-07,
     2  0.4820921255906951D-06,0.1975623797523771D-05,
     3  0.4113248764790821D-05,0.6868432230882420D-05,
     4  0.1067712607400962D-04,0.1436184845752877D-04,
     5  0.1763535320997093D-04,0.2063498695473463D-04,
     6  0.2361019675488914D-04,0.2690294683686638D-04,
     7  0.3105096013102666D-04,0.3708106853130529D-04,
     8  0.4740218517628337D-04,0.6895538896710321D-04,
     9  0.1254684758847792D-03,0.3034691341485716D-03,
     *  0.8998394636892727D-03,0.2792940919962116D-02,
     1  0.8064380945636416D-02,0.2047430289338212D-01,
     2  0.4505951435500846D-01,0.8661612644460944D-01,
     3  0.1478466559590279D+00,0.2282978216006795D+00,
     4  0.3245385685752964D+00,0.4311691136793360D+00,
     5  0.5419685733065455D+00,0.6507449882174135D+00,
     6  0.7518146424478912D+00,0.8402283737405358D+00,
     7  0.9118772613794628D+00,0.9635560791979949D+00,
     *  0.9930210403013454D+00/

        data w4/
     1  0.7137182718756107D-08,0.1428640958056207D-06,
     2  0.8588962026803042D-06,0.2068474810127461D-05,
     3  0.2137657936541478D-05,0.3508554081885053D-05,
     4  0.3864910166306574D-05,0.3472834829429308D-05,
     5  0.3101062898067725D-05,0.2939773548616153D-05,
     6  0.3066060582248771D-05,0.3604903779920366D-05,
     7  0.4850201297806502D-05,0.7570859876794294D-05,
     8  0.1407197601363658D-04,0.3236994958038464D-04,
     9  0.9294672423260586D-04,0.3066069118336629D-03,
     *  0.1021234943283489D-02,0.3101807517838105D-02,
     1  0.8073580688054505D-02,0.1761324112605740D-01,
     2  0.3239062000295955D-01,0.5119184639989259D-01,
     3  0.7118941874545461D-01,0.8912236171925548D-01,
     4  0.1024389797491660D+00,0.1097710478227247D+00,
     5  0.1107882465254575D+00,0.1058138635045045D+00,
     6  0.9549900223337813D-01,0.8064264582152518D-01,
     7  0.6212006016803702D-01,0.4086205129614485D-01,
     *  0.1787280198331713D-01/

        data r5/
     1  0.2470734352249801D-09,0.1912349940143634D-07,
     2  0.2144047417286105D-06,0.9860393439032848D-06,
     3  0.1829116497047657D-05,0.3088864194815432D-05,
     4  0.5964695544468922D-05,0.9560085173965030D-05,
     5  0.1356144120523476D-04,0.1792191888833115D-04,
     6  0.2286337084486753D-04,0.2895543793605528D-04,
     7  0.3737002602446689D-04,0.5059803807724020D-04,
     8  0.7458387147373786D-04,0.1256089763519875D-03,
     9  0.2539658860980783D-03,0.6248853816720818D-03,
     *  0.1762658131079181D-02,0.5104704970129487D-02,
     1  0.1375841420872468D-01,0.3273430573080307D-01,
     2  0.6781594135423240D-01,0.1233302738643669D+00,
     3  0.2002426244200857D+00,0.2957362887414446D+00,
     4  0.4041924227739895D+00,0.5186776661959515D+00,
     5  0.6321621123035077D+00,0.7382161185547859D+00,
     6  0.8313067318750690D+00,0.9068968123829182D+00,
     7  0.9614813280894640D+00,0.9926222536350561D+00/

        data w5/
     1  0.1633818662777671D-08,0.5852634341769963D-07,
     2  0.4076935108969995D-06,0.1149878197328796D-05,
     3  0.4265347888785221D-06,0.2322040061904949D-05,
     4  0.3311513460661938D-05,0.3827797915725162D-05,
     5  0.4167200271987156D-05,0.4589592919623515D-05,
     6  0.5384544772994073D-05,0.6984755581905904D-05,
     7  0.1023002096999087D-04,0.1712397915080187D-04,
     8  0.3324551189475233D-04,0.7607854233906545D-04,
     9  0.2045663249953435D-03,0.6154170588893888D-03,
     *  0.1885172262617024D-02,0.5310234619982778D-02,
     1  0.1285167278962882D-01,0.2610537032388498D-01,
     2  0.4481404042744382D-01,0.6638203583842924D-01,
     3  0.8695430592658117D-01,0.1030697492688891D+00,
     4  0.1126661626996049D+00,0.1151214825319993D+00,
     5  0.1107714020369892D+00,0.1004131150534191D+00,
     6  0.8501143141193097D-01,0.6558408601525370D-01,
     7  0.4317722963464720D-01,0.1889321600881537D-01/

        data r6/
     1  0.4950573639015966D-10,0.1067742502707840D-07,
     2  0.1349999122992512D-06,0.6772866411281184D-06,
     3  0.2031017301137529D-05,0.4446868201631449D-05,
     4  0.7999707612636393D-05,0.1283840900147100D-04,
     5  0.1963037589619494D-04,0.2993445722752328D-04,
     6  0.4703133630954248D-04,0.7875645376069558D-04,
     7  0.1460195666770755D-03,0.3101259291270945D-03,
     8  0.7610499022448045D-03,0.2071004336433264D-02,
     9  0.5744361607003832D-02,0.1493227291984373D-01,
     *  0.3460152999414036D-01,0.7038841120350478D-01,
     1  0.1264368416970531D+00,0.2035883325985839D+00,
     2  0.2990062874844741D+00,0.4071347393158657D+00,
     3  0.5211368394577259D+00,0.6340741720543648D+00,
     4  0.7395875272291478D+00,0.8321929812574004D+00,
     5  0.9073862084790718D+00,0.9616837434441179D+00,
     *  0.9926610097840707D+00/

        data w6/
     1  0.5978092458213636D-09,0.3496701807539413D-07,
     2  0.2666439965775811D-06,0.8884467319274009D-06,
     3  0.1861316482494428D-05,0.2978769232105274D-05,
     4  0.4142989959999330D-05,0.5640714024274475D-05,
     5  0.8192938071876456D-05,0.1292086838982508D-04,
     6  0.2246139414967061D-04,0.4410214973287509D-04,
     7  0.9941248696032155D-04,0.2566637487612200D-03,
     8  0.7312994140828527D-03,0.2125537101651863D-02,
     9  0.5741074898005514D-02,0.1348233539382162D-01,
     *  0.2683558305522989D-01,0.4546200522553063D-01,
     1  0.6678085971195463D-01,0.8703091725915216D-01,
     2  0.1028531195170530D+00,0.1122441810642368D+00,
     3  0.1145926693506494D+00,0.1102177057744669D+00,
     4  0.9989358117875699D-01,0.8456603028800368D-01,
     5  0.6523932347615969D-01,0.4295025132535783D-01,
     *  0.1879395793456611D-01/
c        
c         %%%%%%%%%%%%%%%%%%%%%

ccc        call prin2('alpha = *',alpha,1)
        dtwo = 2
        if (alpha .ge. 1/dtwo) then
        do 1000 i=1,nqs(6)
        roots(i) =r6(i)
        weights(i) = w6(i)
 1000 continue

        nquads = nqs(6)
        endif

        if (alpha .ge. dtwo**(-2)) then
        if (alpha .lt. dtwo**(-1)) then

        do 2000 i=1,nqs(5)
        roots(i) =r5(i)
        weights(i) = w5(i)
 2000 continue

        nquads = nqs(5)

        endif
        endif


        if (alpha .ge. dtwo**(-3)) then
        if (alpha .lt. dtwo**(-2)) then

        do 3000 i=1,nqs(4)
        roots(i) =r4(i)
        weights(i) = w4(i)
 3000 continue

        nquads = nqs(4)

        endif
        endif


        if (alpha .ge. dtwo**(-4)) then
        if (alpha .lt. dtwo**(-3)) then

        do 4000 i=1,nqs(3)
        roots(i) =r3(i)
        weights(i) = w3(i)
 4000 continue

        nquads = nqs(3)

        endif
        endif


        if (alpha .ge. dtwo**(-5)) then
        if (alpha .lt. dtwo**(-4)) then

        do 5000 i=1,nqs(2)
        roots(i) =r2(i)
        weights(i) = w2(i)
 5000 continue

        nquads = nqs(2)

        endif
        endif


        if (alpha .lt. dtwo**(-5)) then

        do 6000 i=1,nqs(1)
        roots(i) =r1(i)
        weights(i) = w1(i)
 6000 continue

        nquads = nqs(1)

        endif

        return
        end
c
c
c
c
c
        subroutine lap_load_quads6(alpha,roots,weights,nquads)
        implicit real *8 (a-h,o-z)
        dimension r1(34),r2(35),r3(35),r4(34),r5(33),r6(30),
     1          w1(34),w2(35),w3(35),w4(34),w5(33),w6(30),nqs(6),
     2          roots(1),weights(1)

        data nqs/34,35,35,34,33,30/

        data r1/
     1  0.1819963154508990D-07,0.6470293158994573D-06,
     2  0.6924377077149701D-05,0.3421730496343683D-04,
     3  0.7039295388478593D-04,0.8706925175016186D-04,
     4  0.9411327140835799D-04,0.9793152693370735D-04,
     5  0.1004768437277980D-03,0.1024701014388531D-03,
     6  0.1042818985800451D-03,0.1061979094411251D-03,
     7  0.1085618533722361D-03,0.1120011172070331D-03,
     8  0.1181785677719659D-03,0.1335956105135498D-03,
     9  0.1963261002842855D-03,0.5206292746637409D-03,
     *  0.1795930035390451D-02,0.5756539148846519D-02,
     1  0.1585564075519776D-01,0.3716730392164989D-01,
     2  0.7503350645330144D-01,0.1330199598043523D+00,
     3  0.2114409691304642D+00,0.3072266202425685D+00,
     4  0.4148920445000214D+00,0.5278423360383994D+00,
     5  0.6394146341905714D+00,0.7434842533268175D+00,
     6  0.8347419374129067D+00,0.9088055574439528D+00,
     7  0.9622739301525131D+00,0.9927743310327933D+00/

        data w1/
     1  0.8412938495749568D-07,0.1828221392194316D-05,
     2  0.1369878689112013D-04,0.3999199258800120D-04,
     3  0.2588257611760556D-04,0.1011275428790987D-04,
     4  0.4909644225582280D-05,0.3007912966576928D-05,
     5  0.2191652780200199D-05,0.1853028034453724D-05,
     6  0.1816155544847717D-05,0.2070088725208464D-05,
     7  0.2753127669749609D-05,0.4362345203721740D-05,
     8  0.8826609302697347D-05,0.2640065721941203D-04,
     9  0.1291284505012731D-03,0.6248723620164419D-03,
     *  0.2201272175930178D-02,0.6295492891514712D-02,
     1  0.1477246142162952D-01,0.2877889076702911D-01,
     2  0.4757445403870925D-01,0.6845220023746759D-01,
     3  0.8785985914328331D-01,0.1027740336654906D+00,
     4  0.1114399714851734D+00,0.1133403902958316D+00,
     5  0.1087792343063884D+00,0.9847351863162610D-01,
     6  0.8331069346227020D-01,0.6424899941900426D-01,
     7  0.4229081711586171D-01,0.1850392044793807D-01/

        data r2/
     1  0.1032741921655046D-07,0.3707228851360952D-06,
     2  0.3600842303543439D-05,0.1334791063268988D-04,
     3  0.3299389060001187D-04,0.6112548492255608D-04,
     4  0.7831164444782268D-04,0.8760914221839177D-04,
     5  0.9348098880046874D-04,0.9777181507496800D-04,
     6  0.1013293986960909D-03,0.1046684012626427D-03,
     7  0.1082334238548843D-03,0.1125823010976320D-03,
     8  0.1187022731808560D-03,0.1289875658098374D-03,
     9  0.1515600387091602D-03,0.2261620843057826D-03,
     *  0.5511774914499785D-03,0.1775820708873975D-02,
     1  0.5559999608352206D-02,0.1526415185567494D-01,
     2  0.3593892763588215D-01,0.7303387668345754D-01,
     3  0.1303146947720349D+00,0.2082798117154714D+00,
     4  0.3039444561702827D+00,0.4118014769386927D+00,
     5  0.5251692802412838D+00,0.6372822716695400D+00,
     6  0.7419255079841665D+00,0.8337206455165880D+00,
     7  0.9082361281611982D+00,0.9620369390703489D+00,
     *  0.9927288038146624D+00/

        data w2/
     1  0.4808694018693344D-07,0.1038474477146325D-05,
     2  0.6397165333545741D-05,0.1276089081721498D-04,
     3  0.2796605800203315D-04,0.2342979169439348D-04,
     4  0.1217161670500576D-04,0.7132997677426405D-05,
     5  0.4890846947075770D-05,0.3822786097824933D-05,
     6  0.3374386730885468D-05,0.3375057145200323D-05,
     7  0.3843562363144264D-05,0.5006210956159275D-05,
     8  0.7583975422063614D-05,0.1407833949975028D-04,
     9  0.3594131808935419D-04,0.1412072299346236D-03,
     *  0.6079980731362196D-03,0.2104015623976964D-02,
     1  0.6022739921588829D-02,0.1425014683625131D-01,
     2  0.2804610187527237D-01,0.4680015580956264D-01,
     3  0.6784742664266031D-01,0.8756709591298852D-01,
     4  0.1028202235170870D+00,0.1117616614126436D+00,
     5  0.1138361344416252D+00,0.1093497165486034D+00,
     6  0.9903879989991911D-01,0.8381215583825024D-01,
     7  0.6464567989226725D-01,0.4255544856764772D-01,
     *  0.1862043039168630D-01/

        data r3/
     1  0.2430909783458630D-08,0.1467941430461446D-06,
     2  0.1903107588548290D-05,0.1045680502331800D-04,
     3  0.2942587486166218D-04,0.5067332747623360D-04,
     4  0.6715975075133592D-04,0.7902036734538396D-04,
     5  0.8798065328388669D-04,0.9536662030323206D-04,
     6  0.1021337602723724D-03,0.1091384301912572D-03,
     7  0.1174198292765242D-03,0.1286699683829335D-03,
     8  0.1464567891083981D-03,0.1801987318640321D-03,
     9  0.2593024862113864D-03,0.4769476372897521D-03,
     *  0.1101240024198094D-02,0.2907596152271531D-02,
     1  0.7840954729242885D-02,0.1958155200192869D-01,
     2  0.4326798995757451D-01,0.8398089756715249D-01,
     3  0.1447157743617252D+00,0.2251427464656276D+00,
     4  0.3217535732355256D+00,0.4289656428165374D+00,
     5  0.5403811560114305D+00,0.6496950031397275D+00,
     6  0.7511773888233231D+00,0.8398784058923029D+00,
     7  0.9117103197931966D+00,0.9634945756265461D+00,
     *  0.9930100663203175D+00/

        data w3/
     1  0.1337028259784817D-07,0.4600821745519124D-06,
     2  0.3963732994971434D-05,0.1410664460982781D-04,
     3  0.2204424516517081D-04,0.1923735649469143D-04,
     4  0.1388420188038352D-04,0.1014571036756945D-04,
     5  0.7992199322494419D-05,0.6934305667957228D-05,
     6  0.6738407967459203D-05,0.7437349424998864D-05,
     7  0.9390578625469774D-05,0.1365175850587723D-04,
     8  0.2330369971863427D-04,0.4842048610941767D-04,
     9  0.1237319614831141D-03,0.3515025914768492D-03,
     *  0.1013912297330747D-02,0.2911643951255836D-02,
     1  0.7571546379212921D-02,0.1679341706960490D-01,
     2  0.3146377304244875D-01,0.5048350982249875D-01,
     3  0.7092742964794998D-01,0.8932119494838864D-01,
     4  0.1029478006126447D+00,0.1103948458660603D+00,
     5  0.1113777282161569D+00,0.1062919860814793D+00,
     6  0.9584662541692972D-01,0.8087339668408589D-01,
     7  0.6225995873575503D-01,0.4093688338156986D-01,
     *  0.1790138916435571D-01/

        data r4/
     1  0.4885764316407881D-08,0.1520225964528051D-06,
     2  0.1335165322233526D-05,0.5748791664838276D-05,
     3  0.1551048185795070D-04,0.3026046632383547D-04,
     4  0.4665345068045615D-04,0.6191743872604438D-04,
     5  0.7552855189069966D-04,0.8807928939167394D-04,
     6  0.1005685020944581D-03,0.1143774020316154D-03,
     7  0.1316646402048454D-03,0.1564557171047365D-03,
     8  0.1977916165889639D-03,0.2799030617448179D-03,
     9  0.4768462534203541D-03,0.1021325217699369D-02,
     *  0.2584643318539502D-02,0.6818081054672245D-02,
     1  0.1697865787651098D-01,0.3796219926666934D-01,
     2  0.7511924206758543D-01,0.1322288321936045D+00,
     3  0.2098809750352147D+00,0.3051920611981826D+00,
     4  0.4127228315752884D+00,0.5258207375605747D+00,
     5  0.6377236897157739D+00,0.7422091799245446D+00,
     6  0.8338890835648701D+00,0.9083235469215212D+00,
     7  0.9620716243408423D+00,0.9927352965895457D+00/

        data w4/
     1  0.2208686596255040D-07,0.4041051050610863D-06,
     2  0.2355588292959093D-05,0.6867423665781630D-05,
     3  0.1262444143769827D-04,0.1622517193264960D-04,
     4  0.1607689060619193D-04,0.1439209109547516D-04,
     5  0.1293950702008432D-04,0.1232809574994895D-04,
     6  0.1287449302325233D-04,0.1508745167092962D-04,
     7  0.2011442589793034D-04,0.3083794436372751D-04,
     8  0.5540339394783858D-04,0.1195893957397663D-03,
     9  0.3087276204901906D-03,0.8828660812977013D-03,
     *  0.2506695323473882D-02,0.6500438219280337D-02,
     1  0.1465736137263092D-01,0.2823887128969420D-01,
     2  0.4673455578606027D-01,0.6758714201736460D-01,
     3  0.8721902032105655D-01,0.1024725328820486D+00,
     4  0.1114614251420155D+00,0.1135969054189073D+00,
     5  0.1091674621283119D+00,0.9890387672017388D-01,
     6  0.8371521770043009D-01,0.6457959617338203D-01,
     7  0.4251544898826669D-01,0.1860371430870008D-01/

        data r5/
     1  0.1346842015479817D-08,0.4820099868342430D-07,
     2  0.5044387084270718D-06,0.2529846304217983D-05,
     3  0.7497181505327395D-05,0.1531098950220943D-04,
     4  0.2587793160225312D-04,0.3987045453705892D-04,
     5  0.5659889646382703D-04,0.7560345012238795D-04,
     6  0.9773030886290021D-04,0.1255818997341992D-03,
     7  0.1647676828149750D-03,0.2275228071477827D-03,
     8  0.3431677323324708D-03,0.5894180938372727D-03,
     9  0.1184891365515578D-02,0.2734207602435245D-02,
     *  0.6737784414160010D-02,0.1626225258806900D-01,
     1  0.3613894258744480D-01,0.7191802844917882D-01,
     2  0.1277663371441660D+00,0.2046207495258095D+00,
     3  0.2997396591162084D+00,0.4076211850725444D+00,
     4  0.5214425743137922D+00,0.6342574218054045D+00,
     5  0.7396918459199535D+00,0.8322483145245032D+00,
     6  0.9074122735435033D+00,0.9616933528948093D+00,
     *  0.9926627322293233D+00/

        data w5/
     1  0.6195889487713051D-08,0.1364915746643803D-06,
     2  0.9754315996575286D-06,0.3348752586112896D-05,
     3  0.6539863937283406D-05,0.9032636310520015D-05,
     4  0.1229223940729079D-04,0.1552270890559462D-04,
     5  0.1784595616857044D-04,0.2029932512455493D-04,
     6  0.2437012534148811D-04,0.3221618993026859D-04,
     7  0.4804741100693540D-04,0.8191105936222769D-04,
     8  0.1610047266044324D-03,0.3642018086747294D-03,
     9  0.9189753230586053D-03,0.2416930690871239D-02,
     *  0.6095556323580597D-02,0.1377635000896792D-01,
     1  0.2694093120899405D-01,0.4534570177309825D-01,
     2  0.6651422065870629D-01,0.8671990842245450D-01,
     3  0.1025749609388368D+00,0.1120304904354038D+00,
     4  0.1144433417760990D+00,0.1101194967752351D+00,
     5  0.9983162925513587D-01,0.8452831154959129D-01,
     6  0.6521736743210880D-01,0.4293859904768289D-01,
     *  0.1878947745775130D-01/

        data r6/
     1  0.3403193665910710D-08,0.9340825980753985D-07,
     2  0.7440594338986110D-06,0.3067315039690537D-05,
     3  0.8453165722882082D-05,0.1818814386061695D-04,
     4  0.3371170029722435D-04,0.5714627760972298D-04,
     5  0.9241255063209328D-04,0.1480003436151646D-03,
     6  0.2432732532746448D-03,0.4237920209376221D-03,
     7  0.8004168643026135D-03,0.1648117570713680D-02,
     8  0.3639510656755988D-02,0.8304918345581734D-02,
     9  0.1862872597190215D-01,0.3923114928068218D-01,
     *  0.7543044981625047D-01,0.1312887493789146D+00,
     1  0.2078065415699217D+00,0.3023957565589766D+00,
     2  0.4096998796145435D+00,0.5229865930639233D+00,
     3  0.6353477288895687D+00,0.7404166453234609D+00,
     4  0.8326903932151478D+00,0.9076463869032013D+00,
     5  0.9617875315575062D+00,0.9926804960211762D+00/

        data w6/
     1  0.1493171475493908D-07,0.2375057446540038D-06,
     2  0.1257243539023226D-05,0.3629847431257922D-05,
     3  0.7349563629607994D-05,0.1235264008750804D-04,
     4  0.1903235031497996D-04,0.2844728191519853D-04,
     5  0.4340630832811839D-04,0.7077190079519109D-04,
     6  0.1267925800760420D-03,0.2512495540925551D-03,
     7  0.5442259719059915D-03,0.1255576339213560D-02,
     8  0.2970925042973352D-02,0.6852051263547999D-02,
     9  0.1458129350965842D-01,0.2754586871632906D-01,
     *  0.4556305412992587D-01,0.6632907374833464D-01,
     1  0.8625922071157628D-01,0.1020012358926921D+00,
     2  0.1114643463081669D+00,0.1139461026621979D+00,
     3  0.1097102050750125D+00,0.9950884325756516D-01,
     4  0.8428426019168496D-01,0.6504435646595075D-01,
     5  0.4283099478645813D-01,0.1874382421913766D-01/
c        
c         %%%%%%%%%%%%%%%%%%%%%

ccc        call prin2('alpha = *',alpha,1)
        dtwo = 2
        if (alpha .ge. 1/dtwo) then
        do 1000 i=1,nqs(6)
        roots(i) =r6(i)
        weights(i) = w6(i)
 1000 continue

        nquads = nqs(6)
        endif

        if (alpha .ge. dtwo**(-2)) then
        if (alpha .lt. dtwo**(-1)) then

        do 2000 i=1,nqs(5)
        roots(i) =r5(i)
        weights(i) = w5(i)
 2000 continue

        nquads = nqs(5)

        endif
        endif


        if (alpha .ge. dtwo**(-3)) then
        if (alpha .lt. dtwo**(-2)) then

        do 3000 i=1,nqs(4)
        roots(i) =r4(i)
        weights(i) = w4(i)
 3000 continue

        nquads = nqs(4)

        endif
        endif


        if (alpha .ge. dtwo**(-4)) then
        if (alpha .lt. dtwo**(-3)) then

        do 4000 i=1,nqs(3)
        roots(i) =r3(i)
        weights(i) = w3(i)
 4000 continue

        nquads = nqs(3)

        endif
        endif


        if (alpha .ge. dtwo**(-5)) then
        if (alpha .lt. dtwo**(-4)) then

        do 5000 i=1,nqs(2)
        roots(i) =r2(i)
        weights(i) = w2(i)
 5000 continue

        nquads = nqs(2)

        endif
        endif


        if (alpha .lt. dtwo**(-5)) then

        do 6000 i=1,nqs(1)
        roots(i) =r1(i)
        weights(i) = w1(i)
 6000 continue

        nquads = nqs(1)

        endif

        return
        end
c
c
c
c
c
        subroutine lap_load_quads7(alpha,roots,weights,nquads)
        implicit real *8 (a-h,o-z)
        dimension r1(35),r2(35),r3(35),r4(34),r5(33),r6(30),
     1          w1(35),w2(35),w3(35),w4(34),w5(33),w6(30),nqs(6),
     2          roots(1),weights(1)

        data nqs/35,35,35,34,33,30/


        data r1/
     1  0.2487312624346449D-07,0.9133202026070931D-06,
     2  0.9976656537866171D-05,0.4792153281024741D-04,
     3  0.1037285516284718D-03,0.2174899557203055D-03,
     4  0.2890633309714642D-03,0.3167425663350694D-03,
     5  0.3309125842208205D-03,0.3400174785627511D-03,
     6  0.3468748260541815D-03,0.3527835486139950D-03,
     7  0.3585935822418601D-03,0.3651281165608433D-03,
     8  0.3735685388827754D-03,0.3863389237468704D-03,
     9  0.4105154907712407D-03,0.4784476356781711D-03,
     *  0.8130137724554206D-03,0.2270015014816349D-02,
     1  0.6635063715611449D-02,0.1729881962113916D-01,
     2  0.3919005876899809D-01,0.7747971393460707D-01,
     3  0.1356302462891007D+00,0.2139569196019337D+00,
     4  0.3094623251999870D+00,0.4167515048358774D+00,
     5  0.5293015959037671D+00,0.6404948793477572D+00,
     6  0.7442306335421216D+00,0.8352111088111143D+00,
     7  0.9090595974281260D+00,0.9623776586324545D+00,
     *  0.9927940552580398D+00/

        data w1/
     1  0.1158172799029205D-06,0.2606737975326503D-05,
     2  0.1983764164719252D-04,0.5126565718414318D-04,
     3  0.7926319523444888D-04,0.1110601112783682D-03,
     4  0.4122978257708696D-04,0.1858546668982478D-04,
     5  0.1093532348895645D-04,0.7694565943788833D-05,
     6  0.6222454731707402D-05,0.5730661548204018D-05,
     7  0.6021120182770183D-05,0.7233075113303711D-05,
     8  0.1001169207638380D-04,0.1651262162774582D-04,
     9  0.3571794575177724D-04,0.1257034911392052D-03,
     *  0.6923494267896844D-03,0.2498674028437909D-02,
     1  0.6797544672278146D-02,0.1537430532491179D-01,
     2  0.2930586174210915D-01,0.4787656268152568D-01,
     3  0.6847910116835024D-01,0.8765737756558669D-01,
     4  0.1024316057230349D+00,0.1110420322556268D+00,
     5  0.1129453266602239D+00,0.1084201057401849D+00,
     6  0.9816672875802843D-01,0.8306391544506636D-01,
     7  0.6406595014972754D-01,0.4217348713103086D-01,
     *  0.1845332416561690D-01/

        data r2/
     1  0.2152879838597372D-07,0.8068323434640118D-06,
     2  0.8919004828605090D-05,0.4597113084002858D-04,
     3  0.9767457152664714D-04,0.1695453985730217D-03,
     4  0.2444209424925005D-03,0.2857966572272453D-03,
     5  0.3101645568181449D-03,0.3268546894748186D-03,
     6  0.3398900632644080D-03,0.3513954198562589D-03,
     7  0.3629062300275764D-03,0.3760233370591469D-03,
     8  0.3931363823250430D-03,0.4191626644232500D-03,
     9  0.4678671811765823D-03,0.5931866201173689D-03,
     *  0.1054192119043335D-02,0.2725150388501871D-02,
     1  0.7536787030098970D-02,0.1900758148259203D-01,
     2  0.4206514427359989D-01,0.8168555584493655D-01,
     3  0.1410051651681112D+00,0.2200459266407508D+00,
     4  0.3156787122386415D+00,0.4225525909023477D+00,
     5  0.5342958708198886D+00,0.6444700247596822D+00,
     6  0.7471335587210025D+00,0.8371123917967729D+00,
     7  0.9101195680968827D+00,0.9628188109846972D+00,
     *  0.9928788050143136D+00/

        data w2/
     1  0.1008155600111533D-06,0.2314501694345112D-05,
     2  0.1787417103734666D-04,0.5657592249574973D-04,
     3  0.4606164776992708D-04,0.9065492394383071D-04,
     4  0.5550323379629892D-04,0.3058380856675903D-04,
     5  0.1958559815877118D-04,0.1441456426543736D-04,
     6  0.1199167089080563D-04,0.1126542104016127D-04,
     7  0.1201263071672320D-04,0.1459641958858509D-04,
     8  0.2037180517669125D-04,0.3363375882619733D-04,
     9  0.7079100044797524D-04,0.2143735664208604D-03,
     *  0.8491314349563087D-03,0.2801242242222658D-02,
     1  0.7410536424523438D-02,0.1637673652469148D-01,
     2  0.3060245790083032D-01,0.4918397865847987D-01,
     3  0.6945859931046545D-01,0.8808398335099525D-01,
     4  0.1022684455699250D+00,0.1104012038737383D+00,
     5  0.1120031879231427D+00,0.1073507393456483D+00,
     6  0.9711220804268050D-01,0.8212999310371625D-01,
     7  0.6332751661049370D-01,0.4168089385523902D-01,
     *  0.1823644036785493D-01/

        data r3/
     1  0.1617505256036589D-07,0.5750555495419473D-06,
     2  0.5920989287288817D-05,0.2842095551337739D-04,
     3  0.6943772665437422D-04,0.1234076473707860D-03,
     4  0.1877069832228859D-03,0.2385903795474658D-03,
     5  0.2753847556431821D-03,0.3035334513093420D-03,
     6  0.3270511418529014D-03,0.3487816071662975D-03,
     7  0.3712535812516342D-03,0.3974973596154514D-03,
     8  0.4323500643033748D-03,0.4857791639510525D-03,
     9  0.5838396549839779D-03,0.8111971590028980D-03,
     *  0.1468275875638804D-02,0.3447539204166504D-02,
     1  0.8776369645318737D-02,0.2110944603376731D-01,
     2  0.4546799853142869D-01,0.8677794025107749D-01,
     3  0.1479115578458122D+00,0.2284833427435924D+00,
     4  0.3249924010011095D+00,0.4319025851137210D+00,
     5  0.5428784349748027D+00,0.6516783920098862D+00,
     6  0.7526306894326636D+00,0.8408356093675163D+00,
     7  0.9122470341633437D+00,0.9637189809234820D+00,
     *  0.9930532953967668D+00/

        data w3/
     1  0.7498518284057079D-07,0.1615080320641557D-05,
     2  0.1137725443357638D-04,0.3457967977774898D-04,
     3  0.4484446118774821D-04,0.6415242440948269D-04,
     4  0.5930368424076492D-04,0.4290275841341366D-04,
     5  0.3165904305368075D-04,0.2529269509234362D-04,
     6  0.2220230608121293D-04,0.2166766265539036D-04,
     7  0.2376239569478060D-04,0.2948550203652984D-04,
     8  0.4173854120766254D-04,0.6892799806627672D-04,
     9  0.1391101661509388D-03,0.3590611101004453D-03,
     *  0.1094835040869629D-02,0.3191079661162222D-02,
     1  0.8078540639354489D-02,0.1745084349395908D-01,
     2  0.3212349243198247D-01,0.5099671619482059D-01,
     3  0.7120257453661928D-01,0.8933734806777863D-01,
     4  0.1027362587076900D+00,0.1100132469438454D+00,
     5  0.1108906947930340D+00,0.1057609424155467D+00,
     6  0.9532586681561557D-01,0.8040890452250391D-01,
     7  0.6188891063211116D-01,0.4068715301486289D-01,
     *  0.1779083434013859D-01/

        data r4/
     1  0.5195286302645521D-08,0.2070170002416046D-06,
     2  0.2305910978342840D-05,0.1249749633320935D-04,
     3  0.4033650426856944D-04,0.8693774623421572D-04,
     4  0.1410760990874595D-03,0.1928375506381782D-03,
     5  0.2393725368477825D-03,0.2817469267509323D-03,
     6  0.3225191762414025D-03,0.3652364486849407D-03,
     7  0.4150203307584360D-03,0.4803272466362484D-03,
     8  0.5774787764890410D-03,0.7432871046291173D-03,
     9  0.1072677506998347D-02,0.1828822640072925D-02,
     *  0.3716003130604670D-02,0.8378097264772706D-02,
     1  0.1896643952924760D-01,0.4020660812292555D-01,
     2  0.7733657333093855D-01,0.1341583968629448D+00,
     3  0.2113835311065410D+00,0.3062582395082249D+00,
     4  0.4134235894269877D+00,0.5262523773953903D+00,
     5  0.6379740725756161D+00,0.7423454221115036D+00,
     6  0.8339574746854390D+00,0.9083539082756880D+00,
     7  0.9620822319380696D+00,0.9927371308604571D+00/

        data w4/
     1  0.2483215206999479D-07,0.5998856128928521D-06,
     2  0.4614574546293406D-05,0.1761412191378932D-04,
     3  0.3837581888919684D-04,0.5260408762907020D-04,
     4  0.5395190680471483D-04,0.4917462704400757D-04,
     5  0.4411173848886132D-04,0.4106934478638149D-04,
     6  0.4105674566783525D-04,0.4520519738794150D-04,
     7  0.5570662359955919D-04,0.7744965249626026D-04,
     8  0.1224629882368178D-03,0.2233767298428400D-03,
     9  0.4749491784421907D-03,0.1146160598579189D-02,
     *  0.2891267158824772D-02,0.6954991030433000D-02,
     1  0.1502750133594252D-01,0.2836172580897009D-01,
     2  0.4656221716732254D-01,0.6720593011655616D-01,
     3  0.8676885909778736D-01,0.1020633328398263D+00,
     4  0.1111435870934096D+00,0.1133747238007943D+00,
     5  0.1090234760142884D+00,0.9881604689340663D-01,
     6  0.8366459742537783D-01,0.6455224633267807D-01,
     7  0.4250208870828602D-01,0.1859890052397644D-01/

        data r5/
     1  0.7422573724251981D-08,0.2289628607252402D-06,
     2  0.2038885459143855D-05,0.9207717242237511D-05,
     3  0.2661821256284761D-04,0.5652570596790915D-04,
     4  0.9736726856514205D-04,0.1462301036866216D-03,
     5  0.2012182275213653D-03,0.2625370331717857D-03,
     6  0.3332367509734553D-03,0.4204446123436248D-03,
     7  0.5380952646910298D-03,0.7133440214118608D-03,
     8  0.1002106448777624D-02,0.1528000892721307D-02,
     9  0.2581266649114568D-02,0.4853632504222448D-02,
     *  0.9890841748689841D-02,0.2066264486186020D-01,
     1  0.4172747842715965D-01,0.7831336497631083D-01,
     2  0.1343996217597281D+00,0.2109398601824313D+00,
     3  0.3053468691954372D+00,0.4123055293372156D+00,
     4  0.5251453273390502D+00,0.6370202428489838D+00,
     5  0.7416149840298420D+00,0.8334648142471463D+00,
     6  0.9080741912875423D+00,0.9619645421231603D+00,
     *  0.9927143966206169D+00/

        data w5/
     1  0.3343776973936726D-07,0.6086453191546617D-06,
     2  0.3658767432643485D-05,0.1157318651896279D-04,
     3  0.2364658322118712D-04,0.3585234648476272D-04,
     4  0.4530484570321406D-04,0.5208989218918825D-04,
     5  0.5790983272804093D-04,0.6521769237899955D-04,
     6  0.7732862834253766D-04,0.9931375910088523D-04,
     7  0.1402834130097920D-03,0.2189615789078704D-03,
     8  0.3777011454076792D-03,0.7188145002858004D-03,
     9  0.1495371075059363D-02,0.3296179640837584D-02,
     *  0.7268337707104577D-02,0.1505031348903336D-01,
     1  0.2798596688290803D-01,0.4588183841932359D-01,
     2  0.6645818400801289D-01,0.8617528090545985D-01,
     3  0.1017280356361513D+00,0.1110572790093211D+00,
     4  0.1134698795896860D+00,0.1092225717638008D+00,
     5  0.9905445308250718D-01,0.8389527522997607D-01,
     6  0.6474316374379081D-01,0.4263256446348303D-01,
     *  0.1865700709874396D-01/

        data r6/
     1  0.6790755136555082D-09,0.7316371048755115D-07,
     2  0.9645513843288480D-06,0.5312596566392021D-05,
     3  0.1757629163117028D-04,0.4249502559432317D-04,
     4  0.8474845628926390D-04,0.1500753036657554D-03,
     5  0.2477943700504340D-03,0.3971110530724728D-03,
     6  0.6411824045062535D-03,0.1079685071666443D-02,
     7  0.1948876908305696D-02,0.3811419256812634D-02,
     8  0.7934593026067321D-02,0.1678533858926773D-01,
     9  0.3406994749251143D-01,0.6313907392805667D-01,
     *  0.1032280212703197D+00,0.1526619778078553D+00,
     1  0.2205115829163587D+00,0.3097633896263971D+00,
     2  0.4140417394932700D+00,0.5255505822130547D+00,
     3  0.6368430551608213D+00,0.7412633748179417D+00,
     4  0.8331430839693296D+00,0.9078626772545887D+00,
     5  0.9618684371478063D+00,0.9926951424986748D+00/

        data w6/
     1  0.5096430962956880D-08,0.2372080048602673D-06,
     2  0.1984433612980797D-05,0.7491472744126793D-05,
     3  0.1782978197543711D-04,0.3276722784201036D-04,
     4  0.5266162609201969D-04,0.7944688723442530D-04,
     5  0.1189790676311627D-03,0.1863660498960043D-03,
     6  0.3172007452646768D-03,0.5962579795273653D-03,
     7  0.1229603446044321D-02,0.2697179684175778D-02,
     8  0.5956587287426441D-02,0.1238587639277225D-01,
     9  0.2281471611935271D-01,0.3522697055902178D-01,
     *  0.4418737321870515D-01,0.5683220714602500D-01,
     1  0.7924899045144862D-01,0.9805813803395595D-01,
     2  0.1091839705385862D+00,0.1125856756076581D+00,
     3  0.1088869459213706D+00,0.9900815809163680D-01,
     4  0.8398115958608792D-01,0.6486546667559914D-01,
     5  0.4273393473867691D-01,0.1870581892520034D-01/
c        
c         %%%%%%%%%%%%%%%%%%%%%

ccc        call prin2('alpha = *',alpha,1)
        dtwo = 2
        if (alpha .ge. 1/dtwo) then
        do 1000 i=1,nqs(6)
        roots(i) =r6(i)
        weights(i) = w6(i)
 1000 continue

        nquads = nqs(6)
        endif

        if (alpha .ge. dtwo**(-2)) then
        if (alpha .lt. dtwo**(-1)) then

        do 2000 i=1,nqs(5)
        roots(i) =r5(i)
        weights(i) = w5(i)
 2000 continue

        nquads = nqs(5)

        endif
        endif


        if (alpha .ge. dtwo**(-3)) then
        if (alpha .lt. dtwo**(-2)) then

        do 3000 i=1,nqs(4)
        roots(i) =r4(i)
        weights(i) = w4(i)
 3000 continue

        nquads = nqs(4)

        endif
        endif


        if (alpha .ge. dtwo**(-4)) then
        if (alpha .lt. dtwo**(-3)) then

        do 4000 i=1,nqs(3)
        roots(i) =r3(i)
        weights(i) = w3(i)
 4000 continue

        nquads = nqs(3)

        endif
        endif


        if (alpha .ge. dtwo**(-5)) then
        if (alpha .lt. dtwo**(-4)) then

        do 5000 i=1,nqs(2)
        roots(i) =r2(i)
        weights(i) = w2(i)
 5000 continue

        nquads = nqs(2)

        endif
        endif


        if (alpha .lt. dtwo**(-5)) then

        do 6000 i=1,nqs(1)
        roots(i) =r1(i)
        weights(i) = w1(i)
 6000 continue

        nquads = nqs(1)

        endif

        return
        end
c
c
c
c
c
        subroutine lap_load_quads8(alpha,roots,weights,nquads)
        implicit real *8 (a-h,o-z)
        dimension r1(35),r2(35),r3(34),r4(33),r5(32),r6(29),
     1          w1(35),w2(35),w3(34),w4(33),w5(32),w6(29),nqs(6),
     2          roots(1),weights(1)

        data nqs/35,35,34,33,32,29/

        data r1/
     1  0.2909868432906089D-07,0.1094285193236835D-05,
     2  0.1224955517930444D-04,0.6676633502651332D-04,
     3  0.2035297664383023D-03,0.4759796097208765D-03,
     4  0.7556968332856968D-03,0.8740913531484868D-03,
     5  0.9277608479386993D-03,0.9591973175431588D-03,
     6  0.9812617335233496D-03,0.9990865679552013D-03,
     7  0.1015499450130764D-02,0.1032750682769910D-02,
     8  0.1053498486348582D-02,0.1082297564335431D-02,
     9  0.1130106451359104D-02,0.1234996168794656D-02,
     *  0.1602179963837567D-02,0.3234487243718297D-02,
     1  0.8198562060114765D-02,0.1980975556032667D-01,
     2  0.4279656791920486D-01,0.8205387082702390D-01,
     3  0.1408064985594572D+00,0.2192657975106145D+00,
     4  0.3144714009844764D+00,0.4211493358184530D+00,
     5  0.5329164721240650D+00,0.6432742364033822D+00,
     6  0.7462092327058331D+00,0.8364833693859415D+00,
     7  0.9097598496283932D+00,0.9626666950416958D+00,
     *  0.9928493381082172D+00/

        data w1/
     1  0.1361611503174099D-06,0.3148622340476582D-05,
     2  0.2477544148907435D-04,0.9195739344187543D-04,
     3  0.1894323506446835D-03,0.3408944938048777D-03,
     4  0.1838693056717150D-03,0.7408026290375124D-04,
     5  0.3919727622020000D-04,0.2554883618403004D-04,
     6  0.1936592356195086D-04,0.1673111537400612D-04,
     7  0.1645551573879601D-04,0.1846928331834191D-04,
     8  0.2372810050095614D-04,0.3546936669586141D-04,
     9  0.6518908027387611D-04,0.1685419209346170D-03,
     *  0.7222899248273074D-03,0.2882219341167907D-02,
     1  0.7594864767680459D-02,0.1643834534535328D-01,
     2  0.3038465623388236D-01,0.4869170239770245D-01,
     3  0.6885080072076555D-01,0.8755821089589722D-01,
     4  0.1019525016357679D+00,0.1103225874324594D+00,
     5  0.1121188694952635D+00,0.1075899804671502D+00,
     6  0.9740513762889157D-01,0.8241937514058386D-01,
     7  0.6357082457515102D-01,0.4184893816882760D-01,
     *  0.1831170537837907D-01/

        data r2/
     1  0.3657675448517922D-08,0.3171001136520413D-06,
     2  0.5257528149721599D-05,0.3891257786672552D-04,
     3  0.1633032045725587D-03,0.4034807506255823D-03,
     4  0.6338645977725506D-03,0.7754386097904803D-03,
     5  0.8583701662971716D-03,0.9128857388253575D-03,
     6  0.9534816003325337D-03,0.9874496735711221D-03,
     7  0.1019434951008638D-02,0.1053521946026189D-02,
     8  0.1094789475865220D-02,0.1151923936435152D-02,
     9  0.1244609443553884D-02,0.1430939906551449D-02,
     *  0.1924131169791879D-02,0.3474298901591026D-02,
     1  0.7919937439971997D-02,0.1857558210990539D-01,
     2  0.4023304675236405D-01,0.7804970984276181D-01,
     3  0.1356155731169730D+00,0.2134177075137444D+00,
     4  0.3085753790376568D+00,0.4157230213142543D+00,
     5  0.5283040656535130D+00,0.6396422571691991D+00,
     6  0.7435795147463060D+00,0.8347722143532910D+00,
     7  0.9088103558985246D+00,0.9622727503331817D+00,
     *  0.9927737842325511D+00/

        data w2/
     1  0.2244156017550048D-07,0.1074753729228778D-05,
     2  0.1229664450163157D-04,0.6619006835407558D-04,
     3  0.1909725400142208D-03,0.2623226724330584D-03,
     4  0.1850875869536747D-03,0.1054946850030823D-03,
     5  0.6532729481946928D-04,0.4595934476703163D-04,
     6  0.3637987293046394D-04,0.3229976317528779D-04,
     7  0.3232868503042175D-04,0.3665472839436794D-04,
     8  0.4722942161869009D-04,0.6997916761804927D-04,
     9  0.1236320966460232D-03,0.2784476015415972D-03,
     *  0.8216515410811828D-03,0.2579897952348259D-02,
     1  0.6860296012543487D-02,0.1526926432527046D-01,
     2  0.2894093300284114D-01,0.4731870375546810D-01,
     3  0.6789809419977785D-01,0.8720923298501672D-01,
     4  0.1021886562361210D+00,0.1109955885138530D+00,
     5  0.1130434863424297D+00,0.1086020501546644D+00,
     6  0.9838014220333187D-01,0.8326944238412400D-01,
     7  0.6423597567395679D-01,0.4228975410592116D-01,
     *  0.1850513124216035D-01/

        data r3/
     1  0.3745513478956037D-07,0.1307752015720087D-05,
     2  0.1329848874965524D-04,0.6601506041239569D-04,
     3  0.1907600428316124D-03,0.3705672483335856D-03,
     4  0.5486865140556914D-03,0.6868682231412372D-03,
     5  0.7887168905083890D-03,0.8677345434067856D-03,
     6  0.9342816086972613D-03,0.9960097631748612D-03,
     7  0.1059886092291089D-02,0.1134293603601513D-02,
     8  0.1232443793291350D-02,0.1380829139582241D-02,
     9  0.1645475663932381D-02,0.2223847859947565D-02,
     *  0.3737631918902896D-02,0.7806814764001373D-02,
     1  0.1759149128267783D-01,0.3787329216070206D-01,
     2  0.7405719603366946D-01,0.1301649411199933D+00,
     3  0.2070622783572330D+00,0.3020212000845751D+00,
     4  0.4096021726515179D+00,0.5230528355170264D+00,
     5  0.6354831012876386D+00,0.7405570869719825D+00,
     6  0.8328010119188879D+00,0.9077150043365137D+00,
     7  0.9618179034702559D+00,0.9926865126902283D+00/

        data w3/
     1  0.1729063932577517D-06,0.3650748629142365D-05,
     2  0.2551740793705111D-04,0.8648505215300845D-04,
     3  0.1597171467275302D-03,0.1893550324912147D-03,
     4  0.1600570486447845D-03,0.1177055374721222D-03,
     5  0.8838252074492842D-04,0.7135045393782826D-04,
     6  0.6298440223283361D-04,0.6160092200338710D-04,
     7  0.6749905255662964D-04,0.8339962357721379D-04,
     8  0.1169647612061744D-03,0.1895963371034772D-03,
     9  0.3682577494260882D-03,0.8813471541090358D-03,
     *  0.2406707262598883D-02,0.6254689448419759D-02,
     1  0.1412819966259835D-01,0.2737337246470065D-01,
     2  0.4569735427315035D-01,0.6666946886280703D-01,
     3  0.8665349714906096D-01,0.1023325210111355D+00,
     4  0.1116839146681506D+00,0.1140578944588558D+00,
     5  0.1097420283844447D+00,0.9949208775361490D-01,
     6  0.8424515336275139D-01,0.6500209725164424D-01,
     7  0.4279843626957719D-01,0.1872853385914390D-01/

        data r4/
     1  0.2599842544971318D-07,0.8494306701316040D-06,
     2  0.8209492695050897D-05,0.4032230657565318D-04,
     3  0.1220125135656392D-03,0.2538534046653272D-03,
     4  0.4059455923361654D-03,0.5521275655774009D-03,
     5  0.6842405952797647D-03,0.8048887106593335D-03,
     6  0.9211296189133454D-03,0.1042976612150137D-02,
     7  0.1184968341282635D-02,0.1371073179337435D-02,
     8  0.1647208823622824D-02,0.2115128539357322D-02,
     9  0.3027315184472808D-02,0.5032942612691536D-02,
     *  0.9666242465886667D-02,0.1993500382085214D-01,
     1  0.4046911920031969D-01,0.7658557721258343D-01,
     2  0.1323465213183668D+00,0.2087533896988658D+00,
     3  0.3032189127837326D+00,0.4103894889676689D+00,
     4  0.5235387816505388D+00,0.6357661379714870D+00,
     5  0.7407120628089200D+00,0.8328794535177804D+00,
     6  0.9077501620338456D+00,0.9618303004414538D+00,
     *  0.9926886702124254D+00/

        data w4/
     1  0.1182337724345132D-06,0.2318814075827526D-05,
     2  0.1545895913890279D-04,0.5358706525324345D-04,
     3  0.1100905847724412D-03,0.1477636861818303D-03,
     4  0.1518819170037042D-03,0.1393051393907771D-03,
     5  0.1254596790250328D-03,0.1170317781599372D-03,
     6  0.1170931364080437D-03,0.1289507374752554D-03,
     7  0.1588457032347783D-03,0.2205329264169845D-03,
     8  0.3473311791046361D-03,0.6268213595848586D-03,
     9  0.1297374584530850D-02,0.2958386193510396D-02,
     *  0.6805106721062895D-02,0.1451954364197426D-01,
     1  0.2747054209550707D-01,0.4547275852995020D-01,
     2  0.6622483820671394D-01,0.8614168693393461D-01,
     3  0.1018718540257425D+00,0.1113276640341831D+00,
     4  0.1138092308150591D+00,0.1095807235768361D+00,
     5  0.9939331629752433D-01,0.8418780874034130D-01,
     6  0.6497076925031228D-01,0.4278292461907400D-01,
     *  0.1872288083474442D-01/

        data r5/
     1  0.3021449378145195D-08,0.1906510697956316D-06,
     2  0.2431116909823128D-05,0.1386821460223004D-04,
     3  0.4714318248611254D-04,0.1120341160729673D-03,
     4  0.2095257312930590D-03,0.3356893619633423D-03,
     5  0.4848955855891783D-03,0.6535498146021062D-03,
     6  0.8452300263143084D-03,0.1074800827797684D-02,
     7  0.1374761969036442D-02,0.1810611524202475D-02,
     8  0.2519350799253785D-02,0.3804777670657434D-02,
     9  0.6357684110308956D-02,0.1167882406205949D-01,
     *  0.2263081370561138D-01,0.4365064442332184D-01,
     1  0.7993411494067510D-01,0.1355586556888996D+00,
     2  0.2116305338009640D+00,0.3056707372163092D+00,
     3  0.4123953357818417D+00,0.5251134608135862D+00,
     4  0.6369426172095347D+00,0.7415350555116479D+00,
     5  0.8334030927966845D+00,0.9080365736793052D+00,
     6  0.9619481001129724D+00,0.9927111622933196D+00/

        data w5/
     1  0.1715466372968170D-07,0.5952030656809595D-06,
     2  0.5047713197338570D-05,0.2010881970746587D-04,
     3  0.4820566062966824D-04,0.8161643404412075D-04,
     4  0.1126841715897414D-03,0.1386387348872926D-03,
     5  0.1590793596652936D-03,0.1787674733413261D-03,
     6  0.2069788354327297D-03,0.2573522545710151D-03,
     7  0.3529167211937245D-03,0.5401904515242916D-03,
     8  0.9244962709291812D-03,0.1754403299630836D-02,
     9  0.3591206027973985D-02,0.7522544945992783D-02,
     *  0.1513412149828110D-01,0.2780591738280740D-01,
     1  0.4547499432081612D-01,0.6596894313832530D-01,
     2  0.8574631590964058D-01,0.1014280208150534D+00,
     3  0.1108851311153289D+00,0.1133922778723603D+00,
     4  0.1092031766058171D+00,0.9906545005818977D-01,
     5  0.8391830294436373D-01,0.6476690647701042D-01,
     6  0.4265035206158086D-01,0.1866524026838476D-01/

        data r6/
     1  0.2088919255414004D-07,0.5977996813208510D-06,
     2  0.5000126309467296D-05,0.2160525725091777D-04,
     3  0.6070047194097136D-04,0.1246151589442633D-03,
     4  0.2134778981527818D-03,0.3523482165962350D-03,
     5  0.5658732621364966D-03,0.8818668102740808D-03,
     6  0.1362408475679668D-02,0.2144870123623909D-02,
     7  0.3530441944826803D-02,0.6180337854995348D-02,
     8  0.1149701853332058D-01,0.2216380165760432D-01,
     9  0.4245368571144502D-01,0.7759001820787705D-01,
     *  0.1319373153650328D+00,0.2069958736694590D+00,
     1  0.3005587846263130D+00,0.4073902375976383D+00,
     2  0.5206767732373346D+00,0.6333490332507481D+00,
     3  0.7388834458912847D+00,0.8316557572243466D+00,
     4  0.9070589389669738D+00,0.9615403961286748D+00,
     *  0.9926327619346728D+00/

        data w6/
     1  0.9243481347882964D-07,0.1547429183580932D-05,
     2  0.8701164062020763D-05,0.2641824491022184D-04,
     3  0.5225075560304162D-04,0.7454361936115005D-04,
     4  0.1085055916062382D-03,0.1729519580229318D-03,
     5  0.2582493468691934D-03,0.3834866430796127D-03,
     6  0.5991538400948774D-03,0.1012920584078221D-02,
     7  0.1862140350386370D-02,0.3661534710980662D-02,
     8  0.7411790313515511D-02,0.1464361649307942D-01,
     9  0.2683960571486250D-01,0.4419788836499407D-01,
     *  0.6476165415024620D-01,0.8497404265837461D-01,
     1  0.1012525214831221D+00,0.1112519177901471D+00,
     2  0.1141298958918244D+00,0.1101220533263850D+00,
     3  0.1000078012040994D+00,0.8476850751882794D-01,
     4  0.6544536406549184D-01,0.4310501276500615D-01,
     *  0.1886583158697211D-01/
c        
c         %%%%%%%%%%%%%%%%%%%%%

ccc        call prin2('alpha = *',alpha,1)
        dtwo = 2
        if (alpha .ge. 1/dtwo) then
        do 1000 i=1,nqs(6)
        roots(i) =r6(i)
        weights(i) = w6(i)
 1000 continue

        nquads = nqs(6)
        endif

        if (alpha .ge. dtwo**(-2)) then
        if (alpha .lt. dtwo**(-1)) then

        do 2000 i=1,nqs(5)
        roots(i) =r5(i)
        weights(i) = w5(i)
 2000 continue

        nquads = nqs(5)

        endif
        endif


        if (alpha .ge. dtwo**(-3)) then
        if (alpha .lt. dtwo**(-2)) then

        do 3000 i=1,nqs(4)
        roots(i) =r4(i)
        weights(i) = w4(i)
 3000 continue

        nquads = nqs(4)

        endif
        endif


        if (alpha .ge. dtwo**(-4)) then
        if (alpha .lt. dtwo**(-3)) then

        do 4000 i=1,nqs(3)
        roots(i) =r3(i)
        weights(i) = w3(i)
 4000 continue

        nquads = nqs(3)

        endif
        endif


        if (alpha .ge. dtwo**(-5)) then
        if (alpha .lt. dtwo**(-4)) then

        do 5000 i=1,nqs(2)
        roots(i) =r2(i)
        weights(i) = w2(i)
 5000 continue

        nquads = nqs(2)

        endif
        endif


        if (alpha .lt. dtwo**(-5)) then

        do 6000 i=1,nqs(1)
        roots(i) =r1(i)
        weights(i) = w1(i)
 6000 continue

        nquads = nqs(1)

        endif

        return
        end
c
c
c
c
c
        subroutine lap_load_quads9(alpha,roots,weights,nquads)
        implicit real *8 (a-h,o-z)
        dimension r1(35),r2(34),r3(35),r4(33),r5(32),r6(28),
     1          w1(35),w2(34),w3(35),w4(33),w5(32),w6(28),nqs(6),
     2          roots(1),weights(1)

        data nqs/35,34,35,33,32,28/


        data r1/
     1  0.9089343354011086D-03,0.4152422948145573D-03,
     2  0.2158062766486727D-01,0.1011580465499302D-01,
     3  0.1835949972156984D-05,0.4355649686015757D-01,
     4  0.5593966115352519D-07,0.1992522893549072D-04,
     5  0.8119383991266169D-01,0.4943840438516259D-02,
     6  0.1187487718830618D-03,0.1615921492579861D-02,
     7  0.1382904510922805D+00,0.2155258144411689D+00,
     8  0.3101453308355006D+00,0.4168431965803397D+00,
     9  0.3307943371084802D-02,0.5290753575820512D+00,
     *  0.2050330165510863D-02,0.6401541697263276D+00,
     1  0.7439032196793086D+00,0.2881810814421040D-02,
     2  0.8349620543461184D+00,0.2234630424320695D-02,
     3  0.9089079858286648D+00,0.2484696450521689D-02,
     4  0.2718273573731626D-02,0.9623112487138695D+00,
     5  0.2525413976424215D-02,0.2332044978486982D-02,
     6  0.2443301755599264D-02,0.2395367866240711D-02,
     7  0.2629987809562208D-02,0.2571143404978491D-02,
     *  0.9927809664645985D+00/

        data w1/
     1  0.6262563420154605D-03,0.4006823724502621D-03,
     2  0.1587356897762351D-01,0.7773409511438750D-02,
     3  0.5120045887981730D-05,0.2898109548923656D-01,
     4  0.2520487606157455D-06,0.4084891255742309D-04,
     5  0.4695784245929820D-01,0.2979118560699491D-02,
     6  0.1815241329759648D-03,0.6398924228311387D-03,
     7  0.6735307703128387D-01,0.8664169992695115D-01,
     8  0.1016894809890427D+00,0.1105968802478447D+00,
     9  0.7421759482416010D-03,0.1127422442438819D+00,
     *  0.2674977064349998D-03,0.1083813457538382D+00,
     1  0.9822159054883047D-01,0.2374349212642193D-03,
     2  0.8315815438851339D-01,0.1266272963806633D-03,
     3  0.6416142904735472D-01,0.4016574268814105D-04,
     4  0.1134353519541150D-03,0.4224517410326901D-01,
     5  0.4218045669539678D-04,0.7576649951647022D-04,
     6  0.4356049273507543D-04,0.5368982208723902D-04,
     7  0.6960186523709478D-04,0.5055075733700257D-04,
     *  0.1848662558284245D-01/

        data r2/
     1  0.7904333628766886D-07,0.2715451965120571D-05,
     2  0.2794399244188884D-04,0.1460307181184369D-03,
     3  0.4585732367125379D-03,0.9881188598586699D-03,
     4  0.1548455577952813D-02,0.1906649025636927D-02,
     5  0.2114829448508237D-02,0.2250646020137861D-02,
     6  0.2351369247545632D-02,0.2435473773380889D-02,
     7  0.2514610836435058D-02,0.2598984459894436D-02,
     8  0.2701310704130999D-02,0.2843510819474922D-02,
     9  0.3075979302992468D-02,0.3551089351312430D-02,
     *  0.4840113184575644D-02,0.8786673145804008D-02,
     1  0.1887163705917894D-01,0.3982522956399818D-01,
     2  0.7683429064410439D-01,0.1336569979160823D+00,
     3  0.2109552261025660D+00,0.3059294031073719D+00,
     4  0.4131900066936738D+00,0.5260960827321595D+00,
     5  0.6378747789044287D+00,0.7422857533661384D+00,
     6  0.8339242284313348D+00,0.9083375713658046D+00,
     7  0.9620760108650619D+00,0.9927359940012299D+00/

        data w2/
     1  0.3628334429630030D-06,0.7568523458170105D-05,
     2  0.5447288229787550D-04,0.2016378084686487D-03,
     3  0.4273122561534262D-03,0.6031020796580102D-03,
     4  0.4672130901951899D-03,0.2660410756112497D-03,
     5  0.1632131793295353D-03,0.1142044997902057D-03,
     6  0.9014369438737166D-04,0.7992851145824830D-04,
     7  0.7998821611827340D-04,0.9078385820475037D-04,
     8  0.1172669596287420D-03,0.1746295585825657D-03,
     9  0.3117590251205485D-03,0.7174032560741590D-03,
     *  0.2158053397798778D-02,0.6325159468361092D-02,
     1  0.1464274205446552D-01,0.2816504573054045D-01,
     2  0.4651004104165237D-01,0.6725137761402728D-01,
     3  0.8686174184128765D-01,0.1021642753920660D+00,
     4  0.1112309675440959D+00,0.1134416424823007D+00,
     5  0.1090711096366479D+00,0.9884838513185011D-01,
     6  0.8368572529953047D-01,0.6456538940360359D-01,
     7  0.4250946847093240D-01,0.1860184418285991D-01/

        data r3/
     1  0.5855618913687272D-03,0.9919168323451511D-03,
     2  0.2711091736209459D-01,0.1425089190799911D-01,
     3  0.7876108050935163D-04,0.1437889047494125D-04,
     4  0.5080326400236026D-01,0.7979600676381761D-02,
     5  0.1392534652176430D-02,0.2627980926020268D-03,
     6  0.1335932070355482D-05,0.8994714935860374D-01,
     7  0.3687700653052151D-07,0.1478638709124072D+00,
     8  0.5201772176972463D-02,0.2250833265258933D+00,
     9  0.1712722957797390D-02,0.3189702860281871D+00,
     *  0.3982790436864876D-02,0.4244603282107592D+00,
     1  0.5352544702106863D+00,0.1954573472919523D-02,
     2  0.2457062696610847D-02,0.3380182423102961D-02,
     3  0.6448578295373512D+00,0.7472268486161212D+00,
     4  0.3029569946095874D-02,0.2144963934822719D-02,
     5  0.8370876041597583D+00,0.2612765201045873D-02,
     6  0.2793425385250298D-02,0.9100734377091236D+00,
     7  0.2306592229671140D-02,0.9627911127624146D+00,
     *  0.9928726272197995D+00/

        data w3/
     1  0.3757765366810281D-03,0.4217211635688813D-03,
     2  0.1749179928863646D-01,0.8947801901669280D-02,
     3  0.1131210505330065D-03,0.2864355594645729D-04,
     4  0.3069405293217780D-01,0.4108182767908465D-02,
     5  0.3654241851232421D-03,0.2590366687601277D-03,
     6  0.3788454453422766D-05,0.4817471720959798D-01,
     7  0.1713931588833163D-06,0.6775305565674611D-01,
     8  0.1768824539890940D-02,0.8623250518749364D-01,
     9  0.2767809308871604D-03,0.1006758332923216D+00,
     *  0.8175884400807035D-03,0.1092366104532044D+00,
     1  0.1112586559309166D+00,0.2117614150564407D-03,
     2  0.1502489899305451D-03,0.4424233097277890D-03,
     3  0.1069358992105749D+00,0.9692111907819305D-01,
     4  0.2795983780353636D-03,0.1727687474501744D-03,
     5  0.8207162393083180D-01,0.1643333975149428D-03,
     6  0.2018226120826634D-03,0.6333375809230131D-01,
     7  0.1533531844191815D-03,0.4170545263954191D-01,
     *  0.1825174547458373D-01/

        data r4/
     1  0.2839183492260374D-07,0.1080203522264097D-05,
     2  0.1185571211240386D-04,0.6534001054124397D-04,
     3  0.2216612312189372D-03,0.5110765602603504D-03,
     4  0.8772557197778932D-03,0.1241466044280207D-02,
     5  0.1570510894404269D-02,0.1865951606643447D-02,
     6  0.2142133481358509D-02,0.2418737954245417D-02,
     7  0.2721803533303353D-02,0.3089963608318338D-02,
     8  0.3589289612716637D-02,0.4349947336929746D-02,
     9  0.5661545471531676D-02,0.8212969787052357D-02,
     *  0.1357517566829095D-01,0.2476469256672281D-01,
     1  0.4630081393955310D-01,0.8325950314231439D-01,
     2  0.1394671010752846D+00,0.2158165869613994D+00,
     3  0.3097586968652392D+00,0.4160700191170886D+00,
     4  0.5281790284616896D+00,0.6393205590077180D+00,
     5  0.7432365539994302D+00,0.8345003390838807D+00,
     6  0.9086414673218391D+00,0.9621979899840039D+00,
     *  0.9927589767267397D+00/

        data w4/
     1  0.1338985003590562D-06,0.3097796496411414D-05,
     2  0.2373110071623082D-04,0.9429233783660232D-04,
     3  0.2249446445390967D-03,0.3428485748419561D-03,
     4  0.3750645222253333D-03,0.3482843497604610D-03,
     5  0.3104552572073315D-03,0.2829165064156525D-03,
     6  0.2727033890725885D-03,0.2847485998522326D-03,
     7  0.3276006658451218D-03,0.4191959020928405D-03,
     8  0.5997595303058065D-03,0.9661741437408617D-03,
     9  0.1763592437037624D-02,0.3589931581233157D-02,
     *  0.7634959149957889D-02,0.1551316763462380D-01,
     1  0.2843791486038346D-01,0.4614831578469851D-01,
     2  0.6642430812516924D-01,0.8583485449600141D-01,
     3  0.1011547735895342D+00,0.1103527770920165D+00,
     4  0.1127263431512951D+00,0.1085087830832556D+00,
     5  0.9841698238044911D-01,0.8336476331482324D-01,
     6  0.6433991040477756D-01,0.4236990760489514D-01,
     *  0.1854276409039960D-01/

        data r5/
     1  0.2303104052495996D-08,0.1767027743516306D-06,
     2  0.2704195305142447D-05,0.1816791867608605D-04,
     3  0.7179121748726823D-04,0.1948850474050289D-03,
     4  0.4030924772116903D-03,0.6881951686516666D-03,
     5  0.1031323073503100D-02,0.1418727851803798D-02,
     6  0.1851120087655276D-02,0.2349799931863291D-02,
     7  0.2965401211224440D-02,0.3797202146456464D-02,
     8  0.5037925972531993D-02,0.7075540073014392D-02,
     9  0.1070909783446126D-01,0.1753252635247839D-01,
     *  0.3039434030149948D-01,0.5350712627990933D-01,
     1  0.9164053559164112D-01,0.1484215703680135D+00,
     2  0.2246955064430061D+00,0.3179954698299600D+00,
     3  0.4232553045955496D+00,0.5340807264258503D+00,
     4  0.6438664786675856D+00,0.7464813034984659D+00,
     5  0.8365921223944067D+00,0.9097952234661776D+00,
     6  0.9626749322885072D+00,0.9928502766194785D+00/

        data w5/
     1  0.1347828271832848D-07,0.5858423044885678D-06,
     2  0.6079452536771712D-05,0.2918759888055604D-04,
     3  0.8373396651966299D-04,0.1652200649963707D-03,
     4  0.2495334841964934D-03,0.3172483916875787D-03,
     5  0.3665721941005470D-03,0.4082671547748069D-03,
     6  0.4599387271355644D-03,0.5455697354795503D-03,
     7  0.7015020371920778D-03,0.9926634411586619D-03,
     8  0.1550030645683875D-02,0.2651872176935781D-02,
     9  0.4872178690002095D-02,0.9245879818582398D-02,
     *  0.1719311765529562D-01,0.2985843356486028D-01,
     1  0.4704323494618001D-01,0.6666877736616152D-01,
     2  0.8545402220557842D-01,0.1002784747927422D+00,
     3  0.1091562721961012D+00,0.1113825704618447D+00,
     4  0.1071621088805028D+00,0.9717679446997276D-01,
     5  0.8231003034095338D-01,0.6352623212973185D-01,
     6  0.4183492770141805D-01,0.1830892638820728D-01/

        data r6/
     1  0.3676351620791584D-08,0.2936554849888773D-06,
     2  0.3944433458382113D-05,0.2336660850926393D-04,
     3  0.8417259493977365D-04,0.2187786550802334D-03,
     4  0.4577795148384161D-03,0.8310462104001935D-03,
     5  0.1380143183588664D-02,0.2183520256210183D-02,
     6  0.3404455265754079D-02,0.5384677666276548D-02,
     7  0.8829537314053595D-02,0.1514416558251101D-01,
     8  0.2689906792681724D-01,0.4811335076070653D-01,
     9  0.8372925068173950D-01,0.1379848237725986D+00,
     *  0.2124701896770285D+00,0.3051857963880947D+00,
     1  0.4110914387487570D+00,0.5234978679837470D+00,
     2  0.6353939909806947D+00,0.7402758283833846D+00,
     3  0.8325223435209523D+00,0.9075251159185730D+00,
     4  0.9617299834007639D+00,0.9926687385072134D+00/

        data w6/
     1  0.2325717713922233D-07,0.9411211766765315D-06,
     2  0.8334062651925948D-05,0.3491324376361668D-04,
     3  0.9228287751961663D-04,0.1820051335719297D-03,
     4  0.3007690339253542D-03,0.4522953259609724D-03,
     5  0.6582232712407467D-03,0.9742916331823673D-03,
     6  0.1521618850028754D-02,0.2550830637479922D-02,
     7  0.4564646711820240D-02,0.8486526878091470D-02,
     8  0.1570072340897343D-01,0.2758079897008195D-01,
     9  0.4439363191033817D-01,0.6439950698154941D-01,
     *  0.8422678128534074D-01,0.1003387167490729D+00,
     1  0.1103345865056547D+00,0.1132962086802046D+00,
     2  0.1094063534046012D+00,0.9941870146368020D-01,
     3  0.8430571857284005D-01,0.6510706287034207D-01,
     4  0.4289000523395039D-01,0.1877350192577945D-01/
c        
c         %%%%%%%%%%%%%%%%%%%%%

ccc        call prin2('alpha = *',alpha,1)
        dtwo = 2
        if (alpha .ge. 1/dtwo) then
        do 1000 i=1,nqs(6)
        roots(i) =r6(i)
        weights(i) = w6(i)
 1000 continue

        nquads = nqs(6)
        endif

        if (alpha .ge. dtwo**(-2)) then
        if (alpha .lt. dtwo**(-1)) then

        do 2000 i=1,nqs(5)
        roots(i) =r5(i)
        weights(i) = w5(i)
 2000 continue

        nquads = nqs(5)

        endif
        endif


        if (alpha .ge. dtwo**(-3)) then
        if (alpha .lt. dtwo**(-2)) then

        do 3000 i=1,nqs(4)
        roots(i) =r4(i)
        weights(i) = w4(i)
 3000 continue

        nquads = nqs(4)

        endif
        endif


        if (alpha .ge. dtwo**(-4)) then
        if (alpha .lt. dtwo**(-3)) then

        do 4000 i=1,nqs(3)
        roots(i) =r3(i)
        weights(i) = w3(i)
 4000 continue

        nquads = nqs(3)

        endif
        endif


        if (alpha .ge. dtwo**(-5)) then
        if (alpha .lt. dtwo**(-4)) then

        do 5000 i=1,nqs(2)
        roots(i) =r2(i)
        weights(i) = w2(i)
 5000 continue

        nquads = nqs(2)

        endif
        endif


        if (alpha .lt. dtwo**(-5)) then

        do 6000 i=1,nqs(1)
        roots(i) =r1(i)
        weights(i) = w1(i)
 6000 continue

        nquads = nqs(1)

        endif

        return
        end
c
c
c
c
c
        subroutine lap_load_quads10(alpha,roots,weights,nquads)
        implicit real *8 (a-h,o-z)
        dimension r1(34),r2(34),r3(34),r4(33),r5(31),r6(28),
     1          w1(34),w2(34),w3(34),w4(33),w5(31),w6(28),nqs(6),
     2          roots(1),weights(1)

        data nqs/34,34,34,33,31,28/


        data r1/
     1  0.2148693787899062D-06,0.7682249827984919D-05,
     2  0.8161607302122335D-04,0.4263760955451694D-03,
     3  0.1265411104497456D-02,0.2718124835186353D-02,
     4  0.4102495373007676D-02,0.4713032720052750D-02,
     5  0.4997024107935637D-02,0.5164847263684975D-02,
     6  0.5283065309137730D-02,0.5378727179909268D-02,
     7  0.5466864830383844D-02,0.5559481043650584D-02,
     8  0.5670734301531807D-02,0.5824735818101985D-02,
     9  0.6078746707774891D-02,0.6625520618483421D-02,
     *  0.8407575325260914D-02,0.1514207752704323D-01,
     1  0.3148754622201863D-01,0.5750940900204241D-01,
     2  0.8695223174664487D-01,0.1349773818496158D+00,
     3  0.2088371967440173D+00,0.3025783977437335D+00,
     4  0.4096364880862224D+00,0.5228843668265689D+00,
     5  0.6352674025656860D+00,0.7403666390593902D+00,
     6  0.8326638085174944D+00,0.9076344137025830D+00,
     7  0.9617833825671719D+00,0.9926797878689644D+00/

        data w1/
     1  0.9950698103522735D-06,0.2170361569570108D-04,
     2  0.1609135601903773D-03,0.5751471519875607D-03,
     3  0.1115114317611964D-02,0.1697149303816893D-02,
     4  0.9299318717098406D-03,0.3892702292015368D-03,
     5  0.2086858991538265D-03,0.1367161593818033D-03,
     6  0.1038712104107695D-03,0.8983288749062405D-04,
     7  0.8836806781885281D-04,0.9911528817455395D-04,
     8  0.1271121205725515D-03,0.1892835842487576D-03,
     9  0.3446244431044859D-03,0.8650154529103255D-03,
     *  0.3332160800883945D-02,0.1098305692876000D-01,
     1  0.2192712569424202D-01,0.2805982070451387D-01,
     2  0.3480820305129330D-01,0.6182476241838920D-01,
     3  0.8483837967534808D-01,0.1015554767409448D+00,
     4  0.1113590525404784D+00,0.1139514273107574D+00,
     5  0.1097408082102864D+00,0.9953668114189385D-01,
     6  0.8430308709559745D-01,0.6505514124271325D-01,
     7  0.4283627643608558D-01,0.1874568977452161D-01/

        data r2/
     1  0.9111653475559301D-07,0.3223830817867696D-05,
     2  0.3504173355603469D-04,0.1999067946877757D-03,
     3  0.7070644821093017D-03,0.1671490310783921D-02,
     4  0.2894738363366409D-02,0.3855552689168288D-02,
     5  0.4426891984743373D-02,0.4781685878879136D-02,
     6  0.5030701916628412D-02,0.5227483251521095D-02,
     7  0.5402088074587509D-02,0.5576901559794492D-02,
     8  0.5775363773684952D-02,0.6031805340771168D-02,
     9  0.6413878591179010D-02,0.7095413281052894D-02,
     *  0.8640577943214550D-02,0.1294072648918396D-01,
     1  0.2405571516193902D-01,0.4695398851328137D-01,
     2  0.8634685804024754D-01,0.1452804497085980D+00,
     3  0.2238117545824914D+00,0.3188736777677175D+00,
     4  0.4251642759760683D+00,0.5363445846841475D+00,
     5  0.6459984271850948D+00,0.7482007196777633D+00,
     6  0.8377900951452599D+00,0.9104896430850553D+00,
     7  0.9629708388433671D+00,0.9929078125901788D+00/

        data w2/
     1  0.4199382751840686D-06,0.9118745599232907D-05,
     2  0.7073481940828930D-04,0.2976221406110756D-03,
     3  0.7400302231692215D-03,0.1160964314348237D-02,
     4  0.1172826614090091D-02,0.7390839771286446D-03,
     5  0.4374339517499661D-03,0.2897328559976241D-03,
     6  0.2166739617272488D-03,0.1816648355036434D-03,
     7  0.1711467027608458D-03,0.1822442126255680D-03,
     8  0.2200920844461776D-03,0.3031211540863468D-03,
     9  0.4864313883505097D-03,0.9560388304115503D-03,
     *  0.2421243651086832D-02,0.6896722192854173D-02,
     1  0.1616900963804079D-01,0.3044069931798234D-01,
     2  0.4887764999441248D-01,0.6900048254919976D-01,
     3  0.8753410693932175D-01,0.1016847432854526D+00,
     4  0.1098236820333440D+00,0.1114584200451166D+00,
     5  0.1068575847246547D+00,0.9668498273648563D-01,
     6  0.8177997179540710D-01,0.6306357352446188D-01,
     7  0.4150965635231401D-01,0.1816209046957590D-01/

        data r3/
     1  0.5926934781452889D-07,0.2226778155197293D-05,
     2  0.2494371519799498D-04,0.1439676919385080D-03,
     3  0.5199569753992906D-03,0.1267131338773681D-02,
     4  0.2220239947389421D-02,0.3091720741734446D-02,
     5  0.3768997191003446D-02,0.4282694115247870D-02,
     6  0.4690386891372539D-02,0.5038502447524686D-02,
     7  0.5363546329251955D-02,0.5699919695407390D-02,
     8  0.6089011193053984D-02,0.6594046022955120D-02,
     9  0.7334225375907498D-02,0.8577340270296229D-02,
     *  0.1099956325716577D-01,0.1627877188534515D-01,
     1  0.2773960362240308D-01,0.5007744590903861D-01,
     2  0.8814306892509011D-01,0.1453319339306011D+00,
     3  0.2221889240494678D+00,0.3160472149289582D+00,
     4  0.4217699518207880D+00,0.5329665855035551D+00,
     5  0.6430547680164328D+00,0.7459201439212281D+00,
     6  0.8362366253975177D+00,0.9096009325703931D+00,
     7  0.9625949812109684D+00,0.9928350020759568D+00/

        data w3/
     1  0.2775929990833511D-06,0.6397922603779030D-05,
     2  0.5078476205334128D-04,0.2159310933262274D-03,
     3  0.5611617876481317D-03,0.9019155338467249D-03,
     4  0.9497483570757206D-03,0.7762369887013951D-03,
     5  0.5856379994113510D-03,0.4517988376527851D-03,
     6  0.3712376567460337D-03,0.3309422036569737D-03,
     7  0.3247381092947552D-03,0.3546860783181753D-03,
     8  0.4335415100670223D-03,0.5950814445269635D-03,
     9  0.9261262760445468D-03,0.1662793909439353D-02,
     *  0.3446775198817568D-02,0.7671093965165316D-02,
     1  0.1605464863596047D-01,0.2945761771325742D-01,
     2  0.4726541625010024D-01,0.6720917955234609D-01,
     3  0.8604611996341280D-01,0.1007923956808376D+00,
     4  0.1095697576950107D+00,0.1117157093208743D+00,
     5  0.1074366758594285D+00,0.9740442418797336D-01,
     6  0.8249373342696135D-01,0.6366449983270777D-01,
     7  0.4192485052489966D-01,0.1834806412883448D-01/

        data r4/
     1  0.5815691967958484D-07,0.1914945366447122D-05,
     2  0.1925931787900631D-04,0.1029461456097793D-03,
     3  0.3542999612762989D-03,0.8546586446162509D-03,
     4  0.1555052720589521D-02,0.2322167405990840D-02,
     5  0.3059761096378575D-02,0.3736702049734519D-02,
     6  0.4363604818787073D-02,0.4972893773069054D-02,
     7  0.5613197154026117D-02,0.6355931383980912D-02,
     8  0.7316328099920923D-02,0.8703905872752086D-02,
     9  0.1094149533173170D-01,0.1492237191301570D-01,
     *  0.2240734291327250D-01,0.3619241566406929D-01,
     1  0.5963154452998261D-01,0.9632675691347394D-01,
     2  0.1500678095076770D+00,0.2230553465167924D+00,
     3  0.3139333531412210D+00,0.4180059108198276D+00,
     4  0.5287384445447233D+00,0.6391820336079454D+00,
     5  0.7428505908093476D+00,0.8341227792273006D+00,
     6  0.9083852192492348D+00,0.9620795284483868D+00,
     *  0.9927350423839006D+00/

        data w4/
     1  0.2643355304159976D-06,0.5267679833192903D-05,
     2  0.3743978256991238D-04,0.1479715037490125D-03,
     3  0.3709816611027657D-03,0.6204147317556537D-03,
     4  0.7555627639959635D-03,0.7624308401573094D-03,
     5  0.7080768689309257D-03,0.6479185555204263D-03,
     6  0.6113730654736528D-03,0.6152030092099071D-03,
     7  0.6769874847933639D-03,0.8269306396412064D-03,
     8  0.1126982203161261D-02,0.1714522046047034D-02,
     9  0.2902107584847263D-02,0.5348296072850551D-02,
     *  0.1009826846943429D-01,0.1803865810206287D-01,
     1  0.2943086780619675D-01,0.4462733263911465D-01,
     2  0.6326226060131018D-01,0.8247869940397516D-01,
     3  0.9846619990638613D-01,0.1085651798542630D+00,
     4  0.1117274497655375D+00,0.1080726066368655D+00,
     5  0.9832562203980012D-01,0.8344923767285791D-01,
     6  0.6448216412832418D-01,0.4249339757170195D-01,
     *  0.1860332457299993D-01/

        data r5/
     1  0.3746556192396846D-07,0.1354486688448407D-05,
     2  0.1415149155329930D-04,0.7478006401075975D-04,
     3  0.2502593488731440D-03,0.6005188877666074D-03,
     4  0.1127836801687583D-02,0.1782119908539078D-02,
     5  0.2514505535130330D-02,0.3315775376866185D-02,
     6  0.4215164471329522D-02,0.5272228870826639D-02,
     7  0.6596540872559336D-02,0.8395667899719037D-02,
     8  0.1106577136837922D-01,0.1537356860507526D-01,
     9  0.2278802813788635D-01,0.3592277879710312D-01,
     *  0.5874433768837532D-01,0.9595121404668545D-01,
     1  0.1513742333754278D+00,0.2262298151017043D+00,
     2  0.3183741369063833D+00,0.4228864048351117D+00,
     3  0.5333543877057628D+00,0.6430728102232621D+00,
     4  0.7457967243174150D+00,0.8360978727866666D+00,
     5  0.9095028951798423D+00,0.9625489441286763D+00,
     *  0.9928256430534253D+00/

        data w5/
     1  0.1747002419542037D-06,0.3816737136886250D-05,
     2  0.2764986357860025D-04,0.1054752640231542D-03,
     3  0.2565599155959704D-03,0.4442097155674713D-03,
     4  0.6009804428180009D-03,0.6987413617177930D-03,
     5  0.7646161895674371D-03,0.8430435236988138D-03,
     6  0.9652933592035503D-03,0.1166330264116613D-02,
     7  0.1515500517275796D-02,0.2145920079590992D-02,
     8  0.3316560402891161D-02,0.5535706218723464D-02,
     9  0.9723580766391452D-02,0.1721666120212210D-01,
     *  0.2923852296057251D-01,0.4584615438770383D-01,
     1  0.6521586452200823D-01,0.8412634879136280D-01,
     2  0.9931859654541572D-01,0.1086164618346495D+00,
     3  0.1111896531196458D+00,0.1072009218248874D+00,
     4  0.9734021903522636D-01,0.8251572650243468D-01,
     5  0.6371660071451010D-01,0.4197239491650070D-01,
     *  0.1837171432082113D-01/

        data r6/
     1  0.1176285836190869D-07,0.4118799507960122D-06,
     2  0.4828261562063566D-05,0.2964435003118847D-04,
     3  0.1155549677211353D-03,0.3250824995661310D-03,
     4  0.7249900766085124D-03,0.1375017533329217D-02,
     5  0.2338924489991333D-02,0.3719326408913739D-02,
     6  0.5718820086317951D-02,0.8743302132955400D-02,
     7  0.1358145614835102D-01,0.2169781276042130D-01,
     8  0.3560337848627482D-01,0.5904203310230991D-01,
     9  0.9649284175759624D-01,0.1517446468104675D+00,
     *  0.2261761318986889D+00,0.3178760269750357D+00,
     1  0.4220867795234957D+00,0.5324470720334299D+00,
     2  0.6422229814567059D+00,0.7451093265639065D+00,
     3  0.8356169386899254D+00,0.9092231249703529D+00,
     4  0.9624294317385746D+00,0.9928023733156208D+00/

        data w6/
     1  0.5301060767770548D-07,0.1188381482755490D-05,
     2  0.1012684141813936D-04,0.4632148365239601D-04,
     3  0.1362829030926243D-03,0.2941279008005750D-03,
     4  0.5154540132455330D-03,0.7946486106923884D-03,
     5  0.1149284620875184D-02,0.1643646622541006D-02,
     6  0.2420113890007409D-02,0.3755400388731166D-02,
     7  0.6158244484515427D-02,0.1048846147528033D-01,
     8  0.1795398633892139D-01,0.2969525695186404D-01,
     9  0.4586983878858728D-01,0.6487952209203574D-01,
     *  0.8365584833787688D-01,0.9892806758013971D-01,
     1  0.1084117263739179D+00,0.1111733607462129D+00,
     2  0.1073216406815535D+00,0.9753403073538776D-01,
     3  0.8272660035365416D-01,0.6390214753736105D-01,
     4  0.4210353750675339D-01,0.1843108134879161D-01/
c        
c         %%%%%%%%%%%%%%%%%%%%%

ccc        call prin2('alpha = *',alpha,1)
        dtwo = 2
        if (alpha .ge. 1/dtwo) then
        do 1000 i=1,nqs(6)
        roots(i) =r6(i)
        weights(i) = w6(i)
 1000 continue

        nquads = nqs(6)
        endif

        if (alpha .ge. dtwo**(-2)) then
        if (alpha .lt. dtwo**(-1)) then

        do 2000 i=1,nqs(5)
        roots(i) =r5(i)
        weights(i) = w5(i)
 2000 continue

        nquads = nqs(5)

        endif
        endif


        if (alpha .ge. dtwo**(-3)) then
        if (alpha .lt. dtwo**(-2)) then

        do 3000 i=1,nqs(4)
        roots(i) =r4(i)
        weights(i) = w4(i)
 3000 continue

        nquads = nqs(4)

        endif
        endif


        if (alpha .ge. dtwo**(-4)) then
        if (alpha .lt. dtwo**(-3)) then

        do 4000 i=1,nqs(3)
        roots(i) =r3(i)
        weights(i) = w3(i)
 4000 continue

        nquads = nqs(3)

        endif
        endif


        if (alpha .ge. dtwo**(-5)) then
        if (alpha .lt. dtwo**(-4)) then

        do 5000 i=1,nqs(2)
        roots(i) =r2(i)
        weights(i) = w2(i)
 5000 continue

        nquads = nqs(2)

        endif
        endif


        if (alpha .lt. dtwo**(-5)) then

        do 6000 i=1,nqs(1)
        roots(i) =r1(i)
        weights(i) = w1(i)
 6000 continue

        nquads = nqs(1)

        endif

        return
        end
c
c
c
c
c
        subroutine lap_load_quads11(alpha,roots,weights,nquads)
        implicit real *8 (a-h,o-z)
        dimension r1(34),r2(34),r3(34),r4(32),r5(31),r6(27),
     1          w1(34),w2(34),w3(34),w4(32),w5(31),w6(27),nqs(6),
     2          roots(1),weights(1)

        data nqs/34,34,34,32,31,27/


        data r1/
     1  0.2197637641680728D-06,0.8534289676715256D-05,
     2  0.9953168780288998D-04,0.6050717074884018D-03,
     3  0.2309179460338404D-02,0.5509631050178328D-02,
     4  0.8040994540707118D-02,0.9185653153168038D-02,
     5  0.9747935936903278D-02,0.1008578593921477D-01,
     6  0.1032255584659281D-01,0.1050979478750200D-01,
     7  0.1067512991524097D-01,0.1083827724518636D-01,
     8  0.1101854446060441D-01,0.1124181659977486D-01,
     9  0.1155422995020201D-01,0.1206181172018344D-01,
     *  0.1307068076198911D-01,0.1553338056166040D-01,
     1  0.2193302955435724D-01,0.3864539094589668D-01,
     2  0.7256322334817104D-01,0.1273359521381852D+00,
     3  0.2035853580661939D+00,0.2984388609685339D+00,
     4  0.4063001774409668D+00,0.5202599215708591D+00,
     5  0.6332979382903348D+00,0.7389844315541118D+00,
     6  0.8317827487501965D+00,0.9071520889944158D+00,
     7  0.9615849376018067D+00,0.9926418938882594D+00/

        data w1/
     1  0.1036048216050354D-05,0.2482716676361977D-04,
     2  0.2071790237126108D-03,0.9424827585330428D-03,
     3  0.2586143050230190D-02,0.3319398642314783D-02,
     4  0.1686183535392637D-02,0.7582248323130428D-03,
     5  0.4187667895009561D-03,0.2752889759143582D-03,
     6  0.2062009604937321D-03,0.1726557727721467D-03,
     7  0.1611904456795810D-03,0.1682193695922851D-03,
     8  0.1964353097358767D-03,0.2572658783532014D-03,
     9  0.3834198559012108D-03,0.6758837978553501D-03,
     *  0.1486584650300983D-02,0.3805886088439435D-02,
     1  0.1012959649592465D-01,0.2448353462137012D-01,
     2  0.4395893299505272D-01,0.6570309108192267D-01,
     3  0.8630135691740967D-01,0.1024472400652271D+00,
     4  0.1121031993620629D+00,0.1146349621431549D+00,
     5  0.1103649239253776D+00,0.1000838932170193D+00,
     6  0.8475543519481419D-01,0.6539819084490322D-01,
     7  0.4305957099162823D-01,0.1884279919211760D-01/

        data r2/
     1  0.1904923200798393D-08,0.6892819408330922D-06,
     2  0.1408967473531557D-04,0.1184161106659714D-03,
     3  0.5741956224714236D-03,0.1759425725187921D-02,
     4  0.3590019707261859D-02,0.5799862944246102D-02,
     5  0.7617479339017589D-02,0.8722684511568080D-02,
     6  0.9415004527644154D-02,0.9902742353482114D-02,
     7  0.1028886298925564D-01,0.1063175333588918D-01,
     8  0.1097510131058061D-01,0.1136468630826073D-01,
     9  0.1186736651656825D-01,0.1261401411639415D-01,
     *  0.1393636775594003D-01,0.1687829137705945D-01,
     1  0.2470041417780511D-01,0.4348447192006267D-01,
     2  0.7892262536539506D-01,0.1346109308152137D+00,
     3  0.2111363215466894D+00,0.3056565449006001D+00,
     4  0.4127133519704216D+00,0.5255837171339989D+00,
     5  0.6374227205399836D+00,0.7419376690999008D+00,
     6  0.8336896006145452D+00,0.9082046597952544D+00,
     7  0.9620202102950804D+00,0.9927252296466668D+00/

        data w2/
     1  0.2565930028724337D-07,0.2567810525773352D-05,
     2  0.3497172409403526D-04,0.2160825315604874D-03,
     3  0.7729955053084842D-03,0.1573320668314729D-02,
     4  0.2066093256700222D-02,0.2177719883833347D-02,
     5  0.1420016452915029D-02,0.8511979347366890D-03,
     6  0.5667478614381421D-03,0.4248736758956262D-03,
     7  0.3566543032915989D-03,0.3361589910542849D-03,
     8  0.3578907923835029D-03,0.4318266855691123D-03,
     9  0.5935531232781946D-03,0.9484518427539082D-03,
     *  0.1845047645284165D-02,0.4547279198801142D-02,
     1  0.1220557950038129D-01,0.2631768066828123D-01,
     2  0.4517907304271930D-01,0.6630379343150220D-01,
     3  0.8625738513982133D-01,0.1018488891499000D+00,
     4  0.1111246455048456D+00,0.1134643848829406D+00,
     5  0.1091603031902111D+00,0.9896153861972353D-01,
     6  0.8379602592405473D-01,0.6465642553874473D-01,
     7  0.4257143387509082D-01,0.1862936598474477D-01/

        data r3/
     1  0.5305916079263595D-01,0.3938604113593643D-01,
     2  0.1943225459596273D-05,0.8152146615638532D-01,
     3  0.2641381096019005D-01,0.1377771494254568D-04,
     4  0.8950482087429543D-07,0.1333582447636612D+00,
     5  0.7949459320671332D-04,0.1915100170737853D-01,
     6  0.1232975142289999D-02,0.2923986305943803D-02,
     7  0.2072750028410124D+00,0.1022677186014090D-01,
     8  0.1087701710344863D-01,0.4952705705652955D-02,
     9  0.1560678022930967D-01,0.1158810660808724D-01,
     *  0.9561727736076951D-02,0.3665455949341139D-03,
     1  0.3003355081839819D+00,0.6639723801787216D-02,
     2  0.8808955938821126D-02,0.1246271688939833D-01,
     3  0.4069946114911978D+00,0.1367796071381776D-01,
     4  0.5202714437712971D+00,0.7877851090607878D-02,
     5  0.6330060278786125D+00,0.7386274257196144D+00,
     6  0.8314869475512623D+00,0.9069647318869434D+00,
     7  0.9615012165612803D+00,0.9926252411243852D+00/

        data w3/
     1  0.1700621390803869D-01,0.1417083640201001D-01,
     2  0.4418227327353716D-05,0.4041470026763916D-01,
     3  0.1015611569743239D-01,0.2494614206343680D-04,
     4  0.3777769338608874D-06,0.6312890143171152D-01,
     5  0.1324311417068836D-03,0.4919092984538395D-02,
     6  0.1285279841225031D-02,0.2003245358413280D-02,
     7  0.8422117864124611D-01,0.6459256577013131D-03,
     8  0.6668061786750221D-03,0.1922950390343893D-02,
     9  0.2512730816043104D-02,0.7717948827762406D-03,
     *  0.6959943408047094D-03,0.5051289528160775D-03,
     1  0.1009517351806987D+00,0.1444190006666163D-02,
     2  0.8242569095358061D-03,0.1005213346988420D-02,
     3  0.1111810008763092D+00,0.1482564461647588D-02,
     4  0.1141668431265206D+00,0.1059357571528794D-02,
     5  0.1102021900324522D+00,0.1000977425546407D+00,
     6  0.8485093779643440D-01,0.6551101926005743D-01,
     7  0.4314879755386118D-01,0.1888508228321226D-01/

        data r4/
     1  0.1001515691012306D-06,0.3776404340441685D-05,
     2  0.4095631172968316D-04,0.2234096199569968D-03,
     3  0.7616833722867422D-03,0.1811109672396339D-02,
     4  0.3246830433906712D-02,0.4751984273917607D-02,
     5  0.6102229788771814D-02,0.7248610633064243D-02,
     6  0.8297236741860510D-02,0.9372397816182190D-02,
     7  0.1052868332117368D-01,0.1184547421869972D-01,
     8  0.1348308911778985D-01,0.1573950410423329D-01,
     9  0.1919660268867193D-01,0.2505448170069047D-01,
     *  0.3572295367886055D-01,0.5532001949515284D-01,
     1  0.8911876852281828D-01,0.1416717645143639D+00,
     2  0.2147579926971681D+00,0.3064487119762930D+00,
     3  0.4116831991501432D+00,0.5236993515790918D+00,
     4  0.6354004061508288D+00,0.7402083417326560D+00,
     5  0.8324474575486057D+00,0.9074723583189696D+00,
     6  0.9617052271023921D+00,0.9926637063582434D+00/

        data w4/
     1  0.4715607806161734D-06,0.1079012885385676D-04,
     2  0.8146259021087992D-04,0.3211775772088482D-03,
     3  0.7864816044811622D-03,0.1289740713039209D-02,
     4  0.1522500721143738D-02,0.1449956563595450D-02,
     5  0.1242354087184501D-02,0.1070811717492547D-02,
     6  0.1049268790317350D-02,0.1107549617773501D-02,
     7  0.1217675091891964D-02,0.1441844650726055D-02,
     8  0.1880839926516769D-02,0.2723718035926675D-02,
     9  0.4380096400710744D-02,0.7725976375095417D-02,
     *  0.1430533480091479D-01,0.2580252085907522D-01,
     1  0.4259138568122725D-01,0.6282773126028525D-01,
     2  0.8301167756858863D-01,0.9949993351141345D-01,
     3  0.1098182133143666D+00,0.1130173429300277D+00,
     4  0.1092826344998050D+00,0.9938555109697738D-01,
     5  0.8431798761011186D-01,0.6513517301768796D-01,
     6  0.4291557941858527D-01,0.1878621827798438D-01/

        data r5/
     1  0.2859128971673654D-08,0.4265245261795798D-06,
     2  0.7389830913323887D-05,0.5327046474161631D-04,
     3  0.2239143118105624D-03,0.6452685334302751D-03,
     4  0.1407915514094463D-02,0.2504649013820712D-02,
     5  0.3857885107222307D-02,0.5392501351039978D-02,
     6  0.7082746266607393D-02,0.8970467566999995D-02,
     7  0.1118094042889835D-01,0.1395731345039929D-01,
     8  0.1773557244386633D-01,0.2329451907830208D-01,
     9  0.3202285786634384D-01,0.4628907792157192D-01,
     *  0.6967969032585697D-01,0.1065868619218745D+00,
     1  0.1608370682542428D+00,0.2339578978194848D+00,
     2  0.3242308066037604D+00,0.4270606331634511D+00,
     3  0.5361829290698529D+00,0.6449025048376093D+00,
     4  0.7469193540894631D+00,0.8367367595765190D+00,
     5  0.9098229228832786D+00,0.9626726545076336D+00,
     *  0.9928484551518309D+00/

        data w5/
     1  0.2309646389482349D-07,0.1499742522383645D-05,
     2  0.1722569462932720D-04,0.8896537009535376D-04,
     3  0.2747973195222304D-03,0.5843197874277773D-03,
     4  0.9391362307855523D-03,0.1239968532482206D-02,
     5  0.1453230708745578D-02,0.1611480776242037D-02,
     6  0.1775764167005433D-02,0.2020030305289406D-02,
     7  0.2440272185880984D-02,0.3182757807269505D-02,
     8  0.4499147884533083D-02,0.6843552650327657D-02,
     9  0.1100246243583467D-01,0.1813343750033725D-01,
     *  0.2940425535106507D-01,0.4508009315132463D-01,
     1  0.6369214019813420D-01,0.8225540918161707D-01,
     2  0.9750218364958526D-01,0.1070917026399137D+00,
     3  0.1100233846309524D+00,0.1063594546883756D+00,
     4  0.9675578959940647D-01,0.8212285477139396D-01,
     5  0.6346472365391612D-01,0.4182716194627089D-01,
     *  0.1831277434265037D-01/

        data r6/
     1  0.5231053288019080D-07,0.1868520199303450D-05,
     2  0.1938671322746177D-04,0.1026765595705205D-03,
     3  0.3516799302153685D-03,0.8932840429822604D-03,
     4  0.1843537212726329D-02,0.3299653142748365D-02,
     5  0.5377771702129658D-02,0.8288352369647700D-02,
     6  0.1244816384163602D-01,0.1864936740110818D-01,
     7  0.2831287515348635D-01,0.4379733194849981D-01,
     8  0.6855962328013581D-01,0.1067497669845536D+00,
     9  0.1619680188993182D+00,0.2356482184066546D+00,
     *  0.3261101119367515D+00,0.4288599901162376D+00,
     1  0.5377394663793286D+00,0.6461386193903748D+00,
     2  0.7478183453136460D+00,0.8373231897724054D+00,
     3  0.9101487942998400D+00,0.9628079700304761D+00,
     *  0.9928744176236167D+00/

        data w6/
     1  0.2431333276622805D-06,0.5247406073351016D-05,
     2  0.3780634404799202D-04,0.1458828943154042D-03,
     3  0.3740368265466827D-03,0.7286769741042332D-03,
     4  0.1187367226622101D-02,0.1742956512446200D-02,
     5  0.2446544112991877D-02,0.3441347375885713D-02,
     6  0.5006065750703442D-02,0.7627572152115169D-02,
     7  0.1208874581289703D-01,0.1946207788694541D-01,
     8  0.3077602113954266D-01,0.4623431208095259D-01,
     9  0.6445859801066927D-01,0.8261544167245059D-01,
     *  0.9753804011092632D-01,0.1069140091802911D+00,
     1  0.1097296757488409D+00,0.1060223737617580D+00,
     2  0.9642539907047075D-01,0.8183262339457339D-01,
     3  0.6323679456075195D-01,0.4167582707234212D-01,
     *  0.1824631378740810D-01/
c        
c         %%%%%%%%%%%%%%%%%%%%%

ccc        call prin2('alpha = *',alpha,1)
        dtwo = 2
        if (alpha .ge. 1/dtwo) then
        do 1000 i=1,nqs(6)
        roots(i) =r6(i)
        weights(i) = w6(i)
 1000 continue

        nquads = nqs(6)
        endif

        if (alpha .ge. dtwo**(-2)) then
        if (alpha .lt. dtwo**(-1)) then

        do 2000 i=1,nqs(5)
        roots(i) =r5(i)
        weights(i) = w5(i)
 2000 continue

        nquads = nqs(5)

        endif
        endif


        if (alpha .ge. dtwo**(-3)) then
        if (alpha .lt. dtwo**(-2)) then

        do 3000 i=1,nqs(4)
        roots(i) =r4(i)
        weights(i) = w4(i)
 3000 continue

        nquads = nqs(4)

        endif
        endif


        if (alpha .ge. dtwo**(-4)) then
        if (alpha .lt. dtwo**(-3)) then

        do 4000 i=1,nqs(3)
        roots(i) =r3(i)
        weights(i) = w3(i)
 4000 continue

        nquads = nqs(3)

        endif
        endif


        if (alpha .ge. dtwo**(-5)) then
        if (alpha .lt. dtwo**(-4)) then

        do 5000 i=1,nqs(2)
        roots(i) =r2(i)
        weights(i) = w2(i)
 5000 continue

        nquads = nqs(2)

        endif
        endif


        if (alpha .lt. dtwo**(-5)) then

        do 6000 i=1,nqs(1)
        roots(i) =r1(i)
        weights(i) = w1(i)
 6000 continue

        nquads = nqs(1)


        endif

        return
        end
c
c
c
c
c
        subroutine lap_load_quads12(alpha,roots,weights,nquads)
        implicit real *8 (a-h,o-z)
        dimension r1(33),r2(34),r3(33),r4(32),r5(30),r6(27),
     1          w1(33),w2(34),w3(33),w4(32),w5(30),w6(27),nqs(6),
     2          roots(1),weights(1)

        data nqs/33,34,33,32,30,27/


        data r1/
     1  0.2600147730252100D-06,0.9893767194866311D-05,
     2  0.1119553582561351D-03,0.6455364840170744D-03,
     3  0.2255039484451097D-02,0.5357266001753866D-02,
     4  0.1033004943690986D-01,0.1472341559030423D-01,
     5  0.1675387218678468D-01,0.1773381704418521D-01,
     6  0.1832100099346486D-01,0.1873707012774743D-01,
     7  0.1907467805606827D-01,0.1938604632571373D-01,
     8  0.1971310403841361D-01,0.2010518950321668D-01,
     9  0.2064545315058882D-01,0.2152711661431559D-01,
     *  0.2336818993064484D-01,0.2878923999668141D-01,
     1  0.4581319295656373D-01,0.8164025148207771D-01,
     2  0.1385232607571747D+00,0.2161562437004015D+00,
     3  0.3112694032178272D+00,0.4182938251465517D+00,
     4  0.5306041026482978D+00,0.6415430430124946D+00,
     5  0.7450117467025081D+00,0.8357331578390923D+00,
     6  0.9093554988972847D+00,0.9625022558476003D+00,
     *  0.9928181475515668D+00/

        data w1/
     1  0.1221189179445056D-05,0.2854611802387426D-04,
     2  0.2287139728595638D-03,0.9598962518536035D-03,
     3  0.2316148112591462D-02,0.4013069270803560D-02,
     4  0.5460746045830157D-02,0.3013778382318449D-02,
     5  0.1328981246324718D-02,0.7269719555222687D-03,
     6  0.4801866763913472D-03,0.3662216793991930D-03,
     7  0.3172754215577624D-03,0.3121929713589302D-03,
     8  0.3497748861690827D-03,0.4472619292563685D-03,
     9  0.6617848153174456D-03,0.1186400532852751D-02,
     *  0.2842663845252987D-02,0.9465391493446197D-02,
     1  0.2581960966844349D-01,0.4617387962798307D-01,
     2  0.6755723468736809D-01,0.8714940800178057D-01,
     3  0.1021261310975166D+00,0.1108010154557819D+00,
     4  0.1127005467806523D+00,0.1081568960671963D+00,
     5  0.9789952618398797D-01,0.8281725685325042D-01,
     6  0.6386401154452636D-01,0.4203549391108555D-01,
     *  0.1839176332411823D-01/

        data r2/
     1  0.2087228111740779D-06,0.8050225295765713D-05,
     2  0.9347345974845816D-04,0.5656729546895491D-03,
     3  0.2174854097135160D-02,0.5634404998129640D-02,
     4  0.9895321474655801D-02,0.1306276791022664D-01,
     5  0.1505545147230131D-01,0.1636833593468672D-01,
     6  0.1730380165228322D-01,0.1802570327849999D-01,
     7  0.1863080930787953D-01,0.1918513598504277D-01,
     8  0.1974425412509078D-01,0.2036830902391891D-01,
     9  0.2114116144005768D-01,0.2221012226461575D-01,
     *  0.2389354123285556D-01,0.2702120615476691D-01,
     1  0.3401246066193556D-01,0.5067994120363675D-01,
     2  0.8425333839174181D-01,0.1387171286096997D+00,
     3  0.2143543379472566D+00,0.3081522037771674D+00,
     4  0.4146016909409594D+00,0.5269672699621181D+00,
     5  0.6383959376600106D+00,0.7425849585840681D+00,
     6  0.8340852321266404D+00,0.9084146597251640D+00,
     7  0.9621048455590743D+00,0.9927412109969642D+00/

        data w2/
     1  0.9820099184000724D-06,0.2337974446565902D-04,
     2  0.1941437989294208D-03,0.8793352057577152D-03,
     3  0.2496635382415112D-02,0.4243190791501842D-02,
     4  0.3881096150874529D-02,0.2489235358936984D-02,
     5  0.1585497193199445D-02,0.1088547838500682D-02,
     6  0.8085786951788778D-03,0.6507889931681211D-03,
     7  0.5701104035954823D-03,0.5475742415496243D-03,
     8  0.5803995387569291D-03,0.6811120910644245D-03,
     9  0.8876492259786301D-03,0.1298746151659856D-02,
     *  0.2190848141859388D-02,0.4425798991204425D-02,
     1  0.1054779995062139D-01,0.2411095457333997D-01,
     2  0.4368206520645655D-01,0.6529119388386184D-01,
     3  0.8546701755800158D-01,0.1011874861801312D+00,
     4  0.1105698183121319D+00,0.1130083522905451D+00,
     5  0.1087939427830094D+00,0.9867431908184857D-01,
     6  0.8357854780098238D-01,0.6450168664862705D-01,
     7  0.4247486132282800D-01,0.1858830445909955D-01/

        data r3/
     1  0.9511417291965105D-07,0.3574967264297781D-05,
     2  0.4077748271072033D-04,0.2468905094406644D-03,
     3  0.9812900677106206D-03,0.2794086012769782D-02,
     4  0.5847991693845617D-02,0.9282106986061780D-02,
     5  0.1214190470910288D-01,0.1428261961373018D-01,
     6  0.1591741302649146D-01,0.1725143308050916D-01,
     7  0.1843589248485918D-01,0.1959589594088838D-01,
     8  0.2086235414912917D-01,0.2241148595186042D-01,
     9  0.2454022088944045D-01,0.2784918941035817D-01,
     *  0.3370521857464627D-01,0.4514816228893631D-01,
     1  0.6730858041280735D-01,0.1054301133263018D+00,
     2  0.1624687846444874D+00,0.2384331178759682D+00,
     3  0.3305557594935263D+00,0.4339611007412693D+00,
     4  0.5426219318054269D+00,0.6502495265324402D+00,
     5  0.7509124380000900D+00,0.8393828397950699D+00,
     6  0.9113070183987098D+00,0.9632921587597427D+00,
     *  0.9929676211721240D+00/

        data w3/
     1  0.4449445420529576D-06,0.1029893728618006D-04,
     2  0.8426551870229756D-04,0.3871673624357199D-03,
     3  0.1182884776049496D-02,0.2488173989667373D-02,
     4  0.3453951846995516D-02,0.3235632639352243D-02,
     5  0.2473820380872087D-02,0.1848516297382799D-02,
     6  0.1455391510768357D-02,0.1237291150861850D-02,
     7  0.1151967121454223D-02,0.1189259646307233D-02,
     8  0.1371914866058480D-02,0.1773172264802195D-02,
     9  0.2577443484626179D-02,0.4252020814595366D-02,
     *  0.7951009527084059D-02,0.1583063121613306D-01,
     1  0.2941106119812950D-01,0.4733061837301063D-01,
     2  0.6674559380303624D-01,0.8471011080319116D-01,
     3  0.9870563351662731D-01,0.1070851850713386D+00,
     4  0.1091780953133666D+00,0.1050830221302216D+00,
     5  0.9536888493136689D-01,0.8084497280353168D-01,
     6  0.6243702383330150D-01,0.4113664077468837D-01,
     *  0.1800789915221305D-01/

        data r4/
     1  0.7271916744031931D-07,0.3053167822150274D-05,
     2  0.3740216606486756D-04,0.2309137956007307D-03,
     3  0.8872171608648876D-03,0.2357321779014917D-02,
     4  0.4650907162414731D-02,0.7335965878582434D-02,
     5  0.9971581655072344D-02,0.1238052231984431D-01,
     6  0.1457521095573290D-01,0.1664879162541121D-01,
     7  0.1873552026369185D-01,0.2101636002781949D-01,
     8  0.2375489356351136D-01,0.2737779110271642D-01,
     9  0.3264756272032036D-01,0.4101314022083385D-01,
     *  0.5517073607452680D-01,0.7941858424958117D-01,
     1  0.1187706959378821D+00,0.1767720303225598D+00,
     2  0.2537555353222943D+00,0.3466754714606934D+00,
     3  0.4501477642971291D+00,0.5578146038262166D+00,
     4  0.6633918006757339D+00,0.7612346847643766D+00,
     5  0.8465584140545045D+00,0.9154937171472509D+00,
     6  0.9650901927021133D+00,0.9933190987498249D+00/

        data w4/
     1  0.3492636872189320D-06,0.9081788171134375D-05,
     2  0.7904210444105524D-04,0.3602198893540321D-03,
     3  0.1019099904524499D-02,0.1925725027875751D-02,
     4  0.2579244915134041D-02,0.2712766200237090D-02,
     5  0.2531250769303414D-02,0.2290629079404795D-02,
     6  0.2114829980606236D-02,0.2054579156935743D-02,
     7  0.2148472891739114D-02,0.2455903283513527D-02,
     8  0.3089614505554034D-02,0.4277654132993977D-02,
     9  0.6493329737955148D-02,0.1067803247170329D-01,
     *  0.1836724984178166D-01,0.3100652490411100D-01,
     1  0.4832038646029765D-01,0.6773912689561420D-01,
     2  0.8570052895240799D-01,0.9921763195935239D-01,
     3  0.1066495772330502D+00,0.1076277461187962D+00,
     4  0.1025826682029397D+00,0.9230495514732812D-01,
     5  0.7769766553043043D-01,0.5968126859377867D-01,
     6  0.3917274706717812D-01,0.1711209798979958D-01/

        data r5/
     1  0.8963043014461470D-07,0.3182277326885785D-05,
     2  0.3331444745601668D-04,0.1792777914984625D-03,
     3  0.6198545466322204D-03,0.1554985732906477D-02,
     4  0.3071494818244335D-02,0.5097951948736083D-02,
     5  0.7499972147260623D-02,0.1018778486627711D-01,
     6  0.1315719958139520D-01,0.1650379594165773D-01,
     7  0.2045512346085537D-01,0.2543732970410233D-01,
     8  0.3219797472402577D-01,0.4202066676842272D-01,
     9  0.5703910146333484D-01,0.8050478244809074D-01,
     *  0.1165747116861693D+00,0.1691712646333939D+00,
     1  0.2402348635831980D+00,0.3285227207047031D+00,
     2  0.4297471988005547D+00,0.5377312924125763D+00,
     3  0.6457230793551786D+00,0.7473145509229811D+00,
     4  0.8369052674793012D+00,0.9098835182097546D+00,
     5  0.9626892294066489D+00,0.9928507763344369D+00/

        data w5/
     1  0.4152369835142118D-06,0.8944363150976369D-05,
     2  0.6544212053234574D-04,0.2574638202590100D-03,
     3  0.6597569281912928D-03,0.1225907416117052D-02,
     4  0.1793152216442791D-02,0.2235367665666087D-02,
     5  0.2553198296801805D-02,0.2821603123757454D-02,
     6  0.3132637202782087D-02,0.3597419186146688D-02,
     7  0.4374311347685636D-02,0.5711183690278102D-02,
     8  0.8018548727257621D-02,0.1197689093732358D-01,
     9  0.1860392897782029D-01,0.2903991404268404D-01,
     *  0.4378622302895177D-01,0.6175257663815200D-01,
     1  0.8016405664678737D-01,0.9567534111222770D-01,
     2  0.1057208438284697D+00,0.1091067148458592D+00,
     3  0.1058016434215368D+00,0.9644528271885954D-01,
     4  0.8196667505137636D-01,0.6339601360855776D-01,
     5  0.4180221611606243D-01,0.1830632768327901D-01/

        data r6/
     1  0.6663360204690531D-07,0.2384306435707261D-05,
     2  0.2501985302490578D-04,0.1351654909007072D-03,
     3  0.4751869824537306D-03,0.1242317226069181D-02,
     4  0.2635970667204419D-02,0.4828392512387179D-02,
     5  0.7994002965560885D-02,0.1239691351568418D-01,
     6  0.1852717344435190D-01,0.2728791079412730D-01,
     7  0.4023520180332583D-01,0.5981143595307946D-01,
     8  0.8935672768162304D-01,0.1325502614657307D+00,
     9  0.1921706794832875D+00,0.2687340621927567D+00,
     *  0.3599106499433988D+00,0.4610452075195950D+00,
     1  0.5662816398775671D+00,0.6696074766153348D+00,
     2  0.7655013540895603D+00,0.8492272539427214D+00,
     3  0.9169336534422492D+00,0.9656767375595477D+00,
     *  0.9934304861537098D+00/

        data w6/
     1  0.3094904064398250D-06,0.6711732621843138D-05,
     2  0.4918792670514694D-04,0.1952126234099168D-03,
     3  0.5187483853763517D-03,0.1049324943650311D-02,
     4  0.1766183847015293D-02,0.2646239478847537D-02,
     5  0.3726737857424793D-02,0.5157666680754552D-02,
     6  0.7249282965877012D-02,0.1052715317241132D-01,
     7  0.1577331696434958D-01,0.2394461750816299D-01,
     8  0.3578128566284082D-01,0.5109424540978933D-01,
     9  0.6824321052327057D-01,0.8448503678170352D-01,
     *  0.9707164085207873D-01,0.1042033385476180D+00,
     1  0.1052568850636042D+00,0.1004700930810465D+00,
     2  0.9052695749771746D-01,0.7628265878769139D-01,
     3  0.5863945939844215D-01,0.3850821307220352D-01,
     *  0.1682628174498075D-01/
c        
c         %%%%%%%%%%%%%%%%%%%%%

ccc        call prin2('alpha = *',alpha,1)
        dtwo = 2
        if (alpha .ge. 1/dtwo) then
        do 1000 i=1,nqs(6)
        roots(i) =r6(i)
        weights(i) = w6(i)
 1000 continue

        nquads = nqs(6)
        endif

        if (alpha .ge. dtwo**(-2)) then
        if (alpha .lt. dtwo**(-1)) then

        do 2000 i=1,nqs(5)
        roots(i) =r5(i)
        weights(i) = w5(i)
 2000 continue

        nquads = nqs(5)

        endif
        endif


        if (alpha .ge. dtwo**(-3)) then
        if (alpha .lt. dtwo**(-2)) then

        do 3000 i=1,nqs(4)
        roots(i) =r4(i)
        weights(i) = w4(i)
 3000 continue

        nquads = nqs(4)

        endif
        endif


        if (alpha .ge. dtwo**(-4)) then
        if (alpha .lt. dtwo**(-3)) then

        do 4000 i=1,nqs(3)
        roots(i) =r3(i)
        weights(i) = w3(i)
 4000 continue

        nquads = nqs(3)

        endif
        endif


        if (alpha .ge. dtwo**(-5)) then
        if (alpha .lt. dtwo**(-4)) then

        do 5000 i=1,nqs(2)
        roots(i) =r2(i)
        weights(i) = w2(i)
 5000 continue

        nquads = nqs(2)

        endif
        endif


        if (alpha .lt. dtwo**(-5)) then

        do 6000 i=1,nqs(1)
        roots(i) =r1(i)
        weights(i) = w1(i)
 6000 continue

        nquads = nqs(1)

        endif

        return
        end
c
c
c
c
c
        subroutine lap_load_quads13(alpha,roots,weights,nquads)
        implicit real *8 (a-h,o-z)
        dimension r1(33),r2(33),r3(33),r4(32),r5(29),r6(26),
     1          w1(33),w2(33),w3(33),w4(32),w5(29),w6(26),nqs(6),
     2          roots(1),weights(1)

        data nqs/33,33,33,32,29,26/


        data r1/
     1  0.7209215566090158D-07,0.4353131444611706D-05,
     2  0.6723673513027334D-04,0.5030751124840565D-03,
     3  0.2338816081392145D-02,0.7541882214101832D-02,
     4  0.1680504226982146D-01,0.2445417385845515D-01,
     5  0.2789337643010806D-01,0.2954349767745464D-01,
     6  0.3052327779104748D-01,0.3120250673442777D-01,
     7  0.3172433038175935D-01,0.3215497864859573D-01,
     8  0.3256746026021401D-01,0.3305358407698073D-01,
     9  0.3368436939706362D-01,0.3458260901181740D-01,
     *  0.3606682383563164D-01,0.3916961307518419D-01,
     1  0.4811426656599046D-01,0.7441400298169132D-01,
     2  0.1252680108780317D+00,0.1995804783833475D+00,
     3  0.2936545242574038D+00,0.4015367775134146D+00,
     4  0.5160356966210188D+00,0.6298888136055381D+00,
     5  0.7364790555948378D+00,0.8301373821766920D+00,
     6  0.9062338429707757D+00,0.9612026502131601D+00,
     *  0.9925684484490720D+00/

        data w1/
     1  0.3821466459676226D-06,0.1414024387412836D-04,
     2  0.1557457331996558D-03,0.8790721632351829D-03,
     3  0.3143281058793544D-02,0.7527138713944656D-02,
     4  0.9809146622721153D-02,0.5127758977911274D-02,
     5  0.2243204231153287D-02,0.1219632188156769D-02,
     6  0.7945967049229701D-03,0.5858520605565557D-03,
     7  0.4674571267331089D-03,0.4054147318602075D-03,
     8  0.4365307314793811D-03,0.5457221440751149D-03,
     9  0.7345457726489115D-03,0.1108945906635437D-02,
     *  0.2001917631070112D-02,0.4775223860338012D-02,
     1  0.1531696771658181D-01,0.3844884611118870D-01,
     2  0.6299794448210720D-01,0.8501714224947475D-01,
     3  0.1021177043936670D+00,0.1124275683720295D+00,
     4  0.1153482109091852D+00,0.1112506947693602D+00,
     5  0.1009837780316556D+00,0.8556092586642269D-01,
     6  0.6603737151630824D-01,0.4348638370223986D-01,
     *  0.1903075312982380D-01/

        data r2/
     1  0.1619589088298826D-06,0.6432236875264909D-05,
     2  0.7891963451277533D-04,0.5124363382455956D-03,
     3  0.2157860118743714D-02,0.6408176091949890D-02,
     4  0.1350037892452432D-01,0.2030937441774141D-01,
     5  0.2463287022945692D-01,0.2721777179739140D-01,
     6  0.2891992698988255D-01,0.3016503669747777D-01,
     7  0.3116068724551871D-01,0.3203136088618804D-01,
     8  0.3289592988219220D-01,0.3388573256664930D-01,
     9  0.3514847145382749D-01,0.3693502321018314D-01,
     *  0.3981890861727442D-01,0.4536042144446740D-01,
     1  0.5800992628858655D-01,0.8639231514050280D-01,
     2  0.1370295014908244D+00,0.2103998082719458D+00,
     3  0.3031594037916810D+00,0.4094963065616740D+00,
     4  0.5223702729173533D+00,0.6346506256081247D+00,
     5  0.7398154539463711D+00,0.8322589287139534D+00,
     6  0.9073926884324237D+00,0.9616786638564563D+00,
     *  0.9926592616159116D+00/

        data w2/
     1  0.7636856920435687D-06,0.1896971512084289D-04,
     2  0.1689225076291094D-03,0.8366203470114920D-03,
     3  0.2710987744459343D-02,0.5900508442116161D-02,
     4  0.7650822133118264D-02,0.5571050309212320D-02,
     5  0.3268064249252215D-02,0.2044118719405546D-02,
     6  0.1426900013539215D-02,0.1095910616429863D-02,
     7  0.9147718320165224D-03,0.8461425890206485D-03,
     8  0.9052102095880034D-03,0.1097470006379359D-02,
     9  0.1466778434121492D-02,0.2192983707693563D-02,
     *  0.3804751062761580D-02,0.7963870971086917D-02,
     1  0.1892134315939074D-01,0.3899489098406599D-01,
     2  0.6230207566364127D-01,0.8386554303026393D-01,
     3  0.1006628617320348D+00,0.1108177322430118D+00,
     4  0.1137303619362769D+00,0.1097387523991156D+00,
     5  0.9965523840249151D-01,0.8446656341126713D-01,
     6  0.6521088902104859D-01,0.4295015153080731D-01,
     *  0.1879797919093005D-01/

        data r3/
     1  0.1540132230814852D-06,0.5792028008056074D-05,
     2  0.6609858168980683D-04,0.3984826845499190D-03,
     3  0.1566114938858163D-02,0.4393043691589948D-02,
     4  0.9094355557567792D-02,0.1436683298300094D-01,
     5  0.1867300195006689D-01,0.2187703880489213D-01,
     6  0.2457870997976720D-01,0.2694452918459378D-01,
     7  0.2901346662328384D-01,0.3091062555062023D-01,
     8  0.3279685161285977D-01,0.3486665657538536D-01,
     9  0.3739220211153410D-01,0.4082898668342312D-01,
     *  0.4606978178658083D-01,0.5504486769957789D-01,
     1  0.7177149559526477D-01,0.1025297832367582D+00,
     2  0.1528961213822550D+00,0.2245412382763434D+00,
     3  0.3150204813319166D+00,0.4189820145860322D+00,
     4  0.5296255906334179D+00,0.6399336093780695D+00,
     5  0.7434277132819741D+00,0.8345148030271974D+00,
     6  0.9086092995115114D+00,0.9621742759024508D+00,
     *  0.9927533958101433D+00/

        data w3/
     1  0.7203387575789920D-06,0.1669370631979776D-04,
     2  0.1365110130060792D-03,0.6217837246961715D-03,
     3  0.1865704455969042D-02,0.3849135718030897D-02,
     4  0.5307668807991315D-02,0.4945652674915413D-02,
     5  0.3658961784777398D-02,0.2879976832389254D-02,
     6  0.2535273832299989D-02,0.2201776834823430D-02,
     7  0.1958189907298883D-02,0.1862869927549208D-02,
     8  0.1941303991959239D-02,0.2242086950790566D-02,
     9  0.2881158509652194D-02,0.4131936959915790D-02,
     *  0.6651799960218282D-02,0.1195451426105961D-01,
     1  0.2260514983107370D-01,0.3990538995719889D-01,
     2  0.6110094385103958D-01,0.8176010848674747D-01,
     3  0.9828088418460168D-01,0.1084876495058951D+00,
     4  0.1116175574918933D+00,0.1079135018112249D+00,
     5  0.9814082909504379D-01,0.8326854348122611D-01,
     6  0.6433055195039976D-01,0.4238872447461066D-01,
     *  0.1855644568662485D-01/

        data r4/
     1  0.1070439835863702D-06,0.4230867355969822D-05,
     2  0.4859793796379647D-04,0.2742289718134742D-03,
     3  0.9270499844113371D-03,0.2261766549790440D-02,
     4  0.4719044587595633D-02,0.8415114672300812D-02,
     5  0.1271901564448349D-01,0.1692727521071883D-01,
     6  0.2072901380716620D-01,0.2413097306335995D-01,
     7  0.2730693136966368D-01,0.3050856557389381D-01,
     8  0.3401821125846825D-01,0.3819772879217492D-01,
     9  0.4362266535917122D-01,0.5130486493260436D-01,
     *  0.6309135979568504D-01,0.8223500741519210D-01,
     1  0.1135828091314351D+00,0.1622447571699443D+00,
     2  0.2309907420835131D+00,0.3186424121681459D+00,
     3  0.4204979461222542D+00,0.5298636666703857D+00,
     4  0.6395599052603878D+00,0.7428816080071085D+00,
     5  0.8340384535128865D+00,0.9083005867329045D+00,
     6  0.9620347329548082D+00,0.9927254920048813D+00/

        data w4/
     1  0.5071975042491110D-06,0.1233438698106014D-04,
     2  0.9915003312548693D-04,0.3974679782129834D-03,
     3  0.9424367421818419D-03,0.1812921530440300D-02,
     4  0.3131337648865030D-02,0.4139009577398953D-02,
     5  0.4343643460852527D-02,0.4023977957614707D-02,
     6  0.3583234039034838D-02,0.3250948263275993D-02,
     7  0.3144004928339122D-02,0.3305372861118239D-02,
     8  0.3771975005722917D-02,0.4679546815085700D-02,
     9  0.6332493919103160D-02,0.9329177157953734D-02,
     *  0.1477910302646051D-01,0.2433787781982161D-01,
     1  0.3926062771878276D-01,0.5855413755942957D-01,
     2  0.7873116508042127D-01,0.9576007145843334D-01,
     3  0.1068091121377389D+00,0.1107099345246108D+00,
     4  0.1075583364879107D+00,0.9811803718411203D-01,
     5  0.8340693540271262D-01,0.6451156868785483D-01,
     6  0.4253634059808949D-01,0.1862721281081065D-01/

        data r5/
     1  0.6975732704135508D-07,0.2724672060556362D-05,
     2  0.3146498296446522D-04,0.1847372455500818D-03,
     3  0.6903853901524879D-03,0.1870031992258389D-02,
     4  0.3991744332959530D-02,0.7110521308394891D-02,
     5  0.1105737418410519D-01,0.1560248524061346D-01,
     6  0.2062701005027590D-01,0.2620994627140913D-01,
     7  0.3267453262675420D-01,0.4066127956969033D-01,
     8  0.5127595267161441D-01,0.6633930210156577D-01,
     9  0.8868364519946161D-01,0.1222075256375945D+00,
     *  0.1711444834903357D+00,0.2383636307861715D+00,
     1  0.3236465111924864D+00,0.4232636903959266D+00,
     2  0.5310134808648408D+00,0.6397710409351509D+00,
     3  0.7426761164355024D+00,0.8337423517651078D+00,
     4  0.9080755381096961D+00,0.9619254973119699D+00,
     *  0.9927029518391035D+00/

        data w5/
     1  0.3291172114199244D-06,0.7932037946631080D-05,
     2  0.6475041453636608D-04,0.2802999931329115D-03,
     3  0.7873045071078412D-03,0.1619805163415009D-02,
     4  0.2633097687727276D-02,0.3572985814093939D-02,
     5  0.4279329776828502D-02,0.4789812877627412D-02,
     6  0.5271796152427858D-02,0.5946841866101841D-02,
     7  0.7086698969911974D-02,0.9066866123760720D-02,
     8  0.1246037342149641D-01,0.1813592357948130D-01,
     9  0.2721375300643484D-01,0.4057052498201779D-01,
     *  0.5781435484912517D-01,0.7659486194283990D-01,
     1  0.9332909260050172D-01,0.1048435408583703D+00,
     2  0.1094484260288172D+00,0.1069137764902473D+00,
     3  0.9789558779385157D-01,0.8342059850909402D-01,
     4  0.6462063878710159D-01,0.4264671360840534D-01,
     *  0.1868398304038589D-01/

        data r6/
     1  0.1106486229044502D-06,0.3738986219447368D-05,
     2  0.3689048746579527D-04,0.1872504863824436D-03,
     3  0.6254313507824776D-03,0.1590828042288246D-02,
     4  0.3368942060902585D-02,0.6262422375230378D-02,
     5  0.1058019107188314D-01,0.1669780385512941D-01,
     6  0.2521407539580855D-01,0.3719298422467869D-01,
     7  0.5445171503020546D-01,0.7978475133719613D-01,
     8  0.1168490004817827D+00,0.1693418744867384D+00,
     9  0.2394913038386957D+00,0.3266681314424393D+00,
     *  0.4271085581220086D+00,0.5348435131016724D+00,
     1  0.6430632277340750D+00,0.7451852649135409D+00,
     2  0.8354263679089416D+00,0.9090276974322612D+00,
     3  0.9623249150333097D+00,0.9927799772651563D+00/

        data w6/
     1  0.5077636452936383D-06,0.1030314253213696D-04,
     2  0.7018922365698261D-04,0.2593407699008930D-03,
     3  0.6572440064801905D-03,0.1321681256161580D-02,
     4  0.2285328409421557D-02,0.3552171339136035D-02,
     5  0.5143313516003592D-02,0.7187603119054456D-02,
     6  0.1001726038508472D-01,0.1423903556381431D-01,
     7  0.2074499241690571D-01,0.3054403014443733D-01,
     8  0.4422808854466526D-01,0.6115506630373077D-01,
     9  0.7904374146339887D-01,0.9466478575487431D-01,
     *  0.1051965383157199D+00,0.1091202665984461D+00,
     1  0.1062130573858689D+00,0.9706388149985637D-01,
     2  0.8262426840052313D-01,0.6396730890247242D-01,
     3  0.4220296109745426D-01,0.1848703467675490D-01/
c        
c         %%%%%%%%%%%%%%%%%%%%%

ccc        call prin2('alpha = *',alpha,1)
        dtwo = 2
        if (alpha .ge. 1/dtwo) then
        do 1000 i=1,nqs(6)
        roots(i) =r6(i)
        weights(i) = w6(i)
 1000 continue

        nquads = nqs(6)
        endif

        if (alpha .ge. dtwo**(-2)) then
        if (alpha .lt. dtwo**(-1)) then

        do 2000 i=1,nqs(5)
        roots(i) =r5(i)
        weights(i) = w5(i)
 2000 continue

        nquads = nqs(5)

        endif
        endif


        if (alpha .ge. dtwo**(-3)) then
        if (alpha .lt. dtwo**(-2)) then

        do 3000 i=1,nqs(4)
        roots(i) =r4(i)
        weights(i) = w4(i)
 3000 continue

        nquads = nqs(4)

        endif
        endif


        if (alpha .ge. dtwo**(-4)) then
        if (alpha .lt. dtwo**(-3)) then

        do 4000 i=1,nqs(3)
        roots(i) =r3(i)
        weights(i) = w3(i)
 4000 continue

        nquads = nqs(3)

        endif
        endif


        if (alpha .ge. dtwo**(-5)) then
        if (alpha .lt. dtwo**(-4)) then

        do 5000 i=1,nqs(2)
        roots(i) =r2(i)
        weights(i) = w2(i)
 5000 continue

        nquads = nqs(2)

        endif
        endif


        if (alpha .lt. dtwo**(-5)) then

        do 6000 i=1,nqs(1)
        roots(i) =r1(i)
        weights(i) = w1(i)
 6000 continue

        nquads = nqs(1)

        endif

        return
        end
c
c
c
c
c
        subroutine lap_load_quads14(alpha,roots,weights,nquads)
        implicit real *8 (a-h,o-z)
        dimension r1(33),r2(33),r3(33),r4(31),r5(29),r6(25),
     1          w1(33),w2(33),w3(33),w4(31),w5(29),w6(25),nqs(6),
     2          roots(1),weights(1)

        data nqs/33,33,33,31,29,25/


        data r1/
     1  0.1751896845897842D-06,0.7164529688440898D-05,
     2  0.9076023646706345D-04,0.6114868919618596D-03,
     3  0.2695153514821585D-02,0.8572035694907685D-02,
     4  0.2020284965305542D-01,0.3360853159188718D-01,
     5  0.4134009018203094D-01,0.4492516795615610D-01,
     6  0.4688939205607976D-01,0.4817154099306324D-01,
     7  0.4912854933860693D-01,0.4992905088358163D-01,
     8  0.5067648206253510D-01,0.5145686737644407D-01,
     9  0.5236940255423064D-01,0.5356962599398584D-01,
     *  0.5537778544107236D-01,0.5864911591421632D-01,
     1  0.6631904939424636D-01,0.8801561832366835D-01,
     2  0.1354972297147635D+00,0.2083246735919929D+00,
     3  0.3013264296230475D+00,0.4080732280958521D+00,
     4  0.5213367684793159D+00,0.6339374045797700D+00,
     5  0.7393504370248479D+00,0.8319795206768759D+00,
     6  0.9072463142238751D+00,0.9616201985650346D+00,
     *  0.9926482760526306D+00/

        data w1/
     1  0.8312674090398416D-06,0.2134423085981178D-04,
     2  0.1973732618463792D-03,0.1022841379957444D-02,
     3  0.3530017615623868D-02,0.8651684623079726D-02,
     4  0.1403670440100214D-01,0.1101731266421529D-01,
     5  0.5047313262662572D-02,0.2530596613889099D-02,
     6  0.1536370588204057D-02,0.1082128338102322D-02,
     7  0.8582962363922274D-03,0.7592560818435733D-03,
     8  0.7494907026575613D-03,0.8271877250656611D-03,
     9  0.1022321684479741D-02,0.1426790552618237D-02,
     *  0.2313277872652454D-02,0.4635416511008237D-02,
     1  0.1225763165832567D-01,0.3364572828750559D-01,
     2  0.6090584345795546D-01,0.8388305169682639D-01,
     3  0.1010403777332875D+00,0.1112325015938103D+00,
     4  0.1140878976021863D+00,0.1100218156959796D+00,
     5  0.9987033115399404D-01,0.8462439804560479D-01,
     6  0.6532036236301121D-01,0.4301726793281301D-01,
     *  0.1882623716513073D-01/

        data r2/
     1  0.2012685565608248D-06,0.7973799373470384D-05,
     2  0.9591750074553230D-04,0.6081697027865548D-03,
     3  0.2517384810950892D-02,0.7509911392643628D-02,
     4  0.1666365775882002D-01,0.2745889143562194D-01,
     5  0.3558602070929281D-01,0.4060300936920565D-01,
     6  0.4381497594309237D-01,0.4608748675787245D-01,
     7  0.4786462527999293D-01,0.4939593573918951D-01,
     8  0.5085398450832950D-01,0.5239497997001462D-01,
     9  0.5420714055763670D-01,0.5658513928668334D-01,
     *  0.6010569659997348D-01,0.6612677104884506D-01,
     1  0.7819149438750096D-01,0.1040406010858479D+00,
     2  0.1519320774935046D+00,0.2232761663015636D+00,
     3  0.3142493053668189D+00,0.4187573485683151D+00,
     4  0.5297611873869924D+00,0.6402282975488685D+00,
     5  0.7437373573559755D+00,0.8347597200360120D+00,
     6  0.9087614193002460D+00,0.9622416248671835D+00,
     *  0.9927667371201790D+00/

        data w2/
     1  0.9513857507362981D-06,0.2340663419662569D-04,
     2  0.2030386443073801D-03,0.9793636988843363D-03,
     3  0.3138143450353678D-02,0.7092409875658674D-02,
     4  0.1078680256231811D-01,0.9913990043164484D-02,
     5  0.6359597825760075D-02,0.3920600711903495D-02,
     6  0.2642708444523604D-02,0.1972464603505976D-02,
     7  0.1621618281051614D-02,0.1468849639809810D-02,
     8  0.1472553650725221D-02,0.1639710392611000D-02,
     9  0.2031219813900183D-02,0.2813924438945728D-02,
     *  0.4432437844170916D-02,0.8150901661945523D-02,
     1  0.1731121403357953D-01,0.3596614612860611D-01,
     2  0.5996735673644669D-01,0.8203981342035245D-01,
     3  0.9887133100477642D-01,0.1089551643663305D+00,
     4  0.1118704797828237D+00,0.1079892390476206D+00,
     5  0.9810602119965823D-01,0.8318251638343829D-01,
     6  0.6423695544643047D-01,0.4231652861477236D-01,
     *  0.1852254023167753D-01/

        data r3/
     1  0.3889747693145864D-05,0.1188683130411989D-06,
     2  0.4148785666886104D-04,0.5890776716079254D+00,
     3  0.4857862516782219D+00,0.6881930082941032D+00,
     4  0.3838600187378202D+00,0.7787852612762492D+00,
     5  0.2897508621268742D+00,0.2098935156610979D+00,
     6  0.8574278544115238D+00,0.2607374361140754D-03,
     7  0.1491627884596923D+00,0.3697109822698679D-02,
     8  0.9211564450457666D+00,0.9037652263706187D-02,
     9  0.1088185678229825D+00,0.1143795280404856D-02,
     *  0.1697348231264432D-01,0.9673045203778732D+00,
     1  0.8513782412822083D-01,0.5122591288277622D-01,
     2  0.2536590913633964D-01,0.4824023163846316D-01,
     3  0.7176219313651778D-01,0.4519674611441294D-01,
     4  0.5447609570214674D-01,0.9937274797464544D+00,
     5  0.3233975068499385D-01,0.4178944594938275D-01,
     6  0.6374324188205413D-01,0.5841775498571031D-01,
     *  0.3765899112271794D-01/

        data w3/
     1  0.1075466633703940D-04,0.5385606108727445D-06,
     2  0.8524438145551357D-04,0.1020261419592392D+00,
     3  0.1036285162438563D+00,0.9548553576009850D-01,
     4  0.9913307375908381D-01,0.8513696326344549D-01,
     5  0.8798698475459238D-01,0.7087975604252041D-01,
     6  0.7165707744616949D-01,0.4289497732641831D-03,
     7  0.5033737453970117D-01,0.3810311412327248D-02,
     8  0.5535114191728173D-01,0.6860798970528026D-02,
     9  0.3101512310881865D-01,0.1509834302834414D-02,
     *  0.8601239872774373D-02,0.3658082058301710D-01,
     1  0.1749327717899076D-01,0.3059095660068437D-02,
     2  0.7845181398114679D-02,0.2964629804381282D-02,
     3  0.1008779268244540D-01,0.3172214215812754D-02,
     4  0.3510014016645708D-02,0.1605641336780092D-01,
     5  0.6088811183492342D-02,0.3701172003783161D-02,
     6  0.6373916735092371D-02,0.4483366153600271D-02,
     *  0.4637934281816073D-02/

        data r4/
     1  0.1901160669533623D-06,0.7455790077281967D-05,
     2  0.8657402979689508D-04,0.5150195400193861D-03,
     3  0.1944095874010244D-02,0.5192079725578998D-02,
     4  0.1050388269230827D-01,0.1710878435750039D-01,
     5  0.2388852780462018D-01,0.3018889513560471D-01,
     6  0.3588084152492350D-01,0.4111044925658445D-01,
     7  0.4615100357897653D-01,0.5137476017123091D-01,
     8  0.5728468862659293D-01,0.6460222411873129D-01,
     9  0.7446148373760042D-01,0.8876771181113349D-01,
     *  0.1106841438804415D+00,0.1447812091671879D+00,
     1  0.1958117054040025D+00,0.2661279378813851D+00,
     2  0.3540269959734657D+00,0.4543838705265672D+00,
     3  0.5604207102215968D+00,0.6652416977854589D+00,
     4  0.7626909022373923D+00,0.8476772974098648D+00,
     5  0.9162277449906932D+00,0.9654367690064462D+00,
     *  0.9933907591781611D+00/

        data w4/
     1  0.8988418127742727D-06,0.2172176969410800D-04,
     2  0.1790029762373009D-03,0.7889684505371225D-03,
     3  0.2220184826051275D-02,0.4327440625239330D-02,
     4  0.6152772253757641D-02,0.6854368172531394D-02,
     5  0.6597172990767781D-02,0.5987157388209630D-02,
     6  0.5423421345847831D-02,0.5081635001961612D-02,
     7  0.5060980392411338D-02,0.5469043880003465D-02,
     8  0.6467969563872672D-02,0.8350194428181676D-02,
     9  0.1167688322354061D-01,0.1745715356084203D-01,
     *  0.2715659494527195D-01,0.4188040626008977D-01,
     1  0.6060422391915423D-01,0.7972953767850266D-01,
     2  0.9518135514912901D-01,0.1043763546882697D+00,
     3  0.1065342711196362D+00,0.1020769940178774D+00,
     4  0.9197241884934391D-01,0.7733916331165885D-01,
     5  0.5927889416309411D-01,0.3882089416511334D-01,
     *  0.1693192204135923D-01/

        data r5/
     1  0.1019818374512996D-06,0.3891700528892268D-05,
     2  0.4381034734902877D-04,0.2533861176010431D-03,
     3  0.9435583989721332D-03,0.2563552014898579D-02,
     4  0.5504618288264707D-02,0.9886610171727298D-02,
     5  0.1553030619485740D-01,0.2213869396868865D-01,
     6  0.2951104429183719D-01,0.3767168057207757D-01,
     7  0.4694046420755060D-01,0.5800859764003932D-01,
     8  0.7206065881945817D-01,0.9095092842716907D-01,
     9  0.1173620133422466D+00,0.1547038286388376D+00,
     *  0.2063880314555325D+00,0.2744484933138301D+00,
     1  0.3582617070461799D+00,0.4543070891583438D+00,
     2  0.5570313953436559D+00,0.6600885589352746D+00,
     3  0.7573081457132092D+00,0.8432286836049554D+00,
     4  0.9133262360365807D+00,0.9641005106779065D+00,
     *  0.9931197362409661D+00/

        data w5/
     1  0.4797051926956732D-06,0.1122162521881417D-04,
     2  0.8922491184432919D-04,0.3822486342351504D-03,
     3  0.1076623595328734D-02,0.2232741051427161D-02,
     4  0.3669909938278416D-02,0.5060596049873246D-02,
     5  0.6173420853446359D-02,0.7007888307374630D-02,
     6  0.7739859872545451D-02,0.8633553622011839D-02,
     7  0.1001721639939638D-01,0.1231341630367872D-01,
     8  0.1609654830351140D-01,0.2213170839262941D-01,
     9  0.3127049266369518D-01,0.4400681738973789D-01,
     *  0.5972475307936111D-01,0.7629630132608754D-01,
     1  0.9072674126241679D-01,0.1004161315398263D+00,
     2  0.1039583996289363D+00,0.1011190613153749D+00,
     3  0.9240647404314187D-01,0.7867939502710056D-01,
     4  0.6093174223835057D-01,0.4021027320632300D-01,
     *  0.1761675971365528D-01/

        data r6/
     1  0.9200473787558814D-07,0.3585923999018997D-05,
     2  0.4095530562577882D-04,0.2386036625902381D-03,
     3  0.8938888110251579D-03,0.2463349670272235D-02,
     4  0.5458898201907253D-02,0.1035326728029607D-01,
     5  0.1758950875014900D-01,0.2772564930922864D-01,
     6  0.4168574626267199D-01,0.6105500276629350D-01,
     7  0.8830781410195981D-01,0.1267354772796793D+00,
     8  0.1797506812808203D+00,0.2495256263170486D+00,
     9  0.3356243994428460D+00,0.4345930671322732D+00,
     *  0.5407411799498535D+00,0.6474440134353480D+00,
     1  0.7482208359273624D+00,0.8373390897567349D+00,
     2  0.9100653289034241D+00,0.9627491654954613D+00,
     *  0.9928607108097665D+00/

        data w6/
     1  0.4348314715788823D-06,0.1040624954189502D-04,
     2  0.8383345827161270D-04,0.3612033352672331D-03,
     3  0.1027496606292262D-02,0.2198560131729574D-02,
     4  0.3871240776429678D-02,0.5988842558095937D-02,
     5  0.8572314659904468D-02,0.1184883410814081D-01,
     6  0.1633121707424670D-01,0.2281997597973513D-01,
     7  0.3224541684039903D-01,0.4520311729836208D-01,
     8  0.6121499935884417D-01,0.7826970218770486D-01,
     9  0.9333709147222533D-01,0.1036281106338241D+00,
     *  0.1075455361785595D+00,0.1047703728681365D+00,
     1  0.9582415288174648D-01,0.8162231169313076D-01,
     2  0.6322132512436508D-01,0.4172352840075963D-01,
     *  0.1827997529281559D-01/
c        
c        
c         %%%%%%%%%%%%%%%%%%%%%

ccc        call prin2('alpha = *',alpha,1)
        dtwo = 2
        if (alpha .ge. 1/dtwo) then
        do 1000 i=1,nqs(6)
        roots(i) =r6(i)
        weights(i) = w6(i)
 1000 continue

        nquads = nqs(6)
        endif

        if (alpha .ge. dtwo**(-2)) then
        if (alpha .lt. dtwo**(-1)) then

        do 2000 i=1,nqs(5)
        roots(i) =r5(i)
        weights(i) = w5(i)
 2000 continue

        nquads = nqs(5)

        endif
        endif


        if (alpha .ge. dtwo**(-3)) then
        if (alpha .lt. dtwo**(-2)) then

        do 3000 i=1,nqs(4)
        roots(i) =r4(i)
        weights(i) = w4(i)
 3000 continue

        nquads = nqs(4)

        endif
        endif


        if (alpha .ge. dtwo**(-4)) then
        if (alpha .lt. dtwo**(-3)) then

        do 4000 i=1,nqs(3)
        roots(i) =r3(i)
        weights(i) = w3(i)
 4000 continue

        nquads = nqs(3)

        endif
        endif


        if (alpha .ge. dtwo**(-5)) then
        if (alpha .lt. dtwo**(-4)) then

        do 5000 i=1,nqs(2)
        roots(i) =r2(i)
        weights(i) = w2(i)
 5000 continue

        nquads = nqs(2)

        endif
        endif


        if (alpha .lt. dtwo**(-5)) then

        do 6000 i=1,nqs(1)
        roots(i) =r1(i)
        weights(i) = w1(i)
 6000 continue

        nquads = nqs(1)

        endif

        return
        end
c
c
c
c
c
        subroutine lap_load_quads15(alpha,roots,weights,nquads)
        implicit real *8 (a-h,o-z)
        dimension r1(33),r2(33),r3(33),r4(31),r5(28),r6(24),
     1          w1(33),w2(33),w3(33),w4(31),w5(28),w6(24),nqs(6),
     2          roots(1),weights(1)

        data nqs/33,33,33,31,28,24/


        data r1/
     1  0.5295079827384090D-06,0.2041252519662560D-04,
     2  0.2354718295695405D-03,0.1414169522455729D-02,
     3  0.5498534657048116D-02,0.1539064192183174D-01,
     4  0.3222085113633526D-01,0.4988547619995615D-01,
     5  0.6060387876130613D-01,0.6601700420312188D-01,
     6  0.6908455607842244D-01,0.7109746034994100D-01,
     7  0.7259336634314954D-01,0.7383737537820433D-01,
     8  0.7499549613613046D-01,0.7620719640101376D-01,
     9  0.7763337358468075D-01,0.7952420431387879D-01,
     *  0.8238351058619235D-01,0.8748182568763812D-01,
     1  0.9861408313188701D-01,0.1254167529874101D+00,
     2  0.1761274599479964D+00,0.2473067868705679D+00,
     3  0.3342123527300048D+00,0.4333568952657758D+00,
     4  0.5396912207170284D+00,0.6466782808195544D+00,
     5  0.7477348017693314D+00,0.8370675493266672D+00,
     6  0.9099351227862818D+00,0.9627011436544318D+00,
     *  0.9928521413505088D+00/

        data w1/
     1  0.2492604639955403D-05,0.5920788747262752D-04,
     2  0.4870170009276671D-03,0.2193404011757926D-02,
     3  0.6499067939284576D-02,0.1359101963233669D-01,
     4  0.1902454130433178D-01,0.1461096454892585D-01,
     5  0.7404311316369633D-02,0.3918989612767096D-02,
     6  0.2412326609496229D-02,0.1696186383435396D-02,
     7  0.1337337336663670D-02,0.1177286720097819D-02,
     8  0.1161486901698482D-02,0.1287741830403728D-02,
     9  0.1604012640031192D-02,0.2254781034584875D-02,
     *  0.3651690242880268D-02,0.7099901098447016D-02,
     1  0.1683284818775322D-01,0.3852762831081929D-01,
     2  0.6196162651078401D-01,0.7961443962600232D-01,
     3  0.9368188917338439D-01,0.1037547193797714D+00,
     4  0.1077958031942140D+00,0.1050681869270971D+00,
     5  0.9607589989580711D-01,0.8179881598493717D-01,
     6  0.6333013651051674D-01,0.4178194107344389D-01,
     *  0.1830229856891676D-01/

        data r2/
     1  0.9330069816776123D-07,0.5327080734596219D-05,
     2  0.7893288294102843D-04,0.5692333815860413D-03,
     3  0.2554106886508157D-02,0.8057590626585915D-02,
     4  0.1894608426390735D-01,0.3383251699410158D-01,
     5  0.4731283681910917D-01,0.5646479054627678D-01,
     6  0.6234710272372470D-01,0.6638372513348769D-01,
     7  0.6941278283173804D-01,0.7190337452381135D-01,
     8  0.7415776855583885D-01,0.7642077509360442D-01,
     9  0.7895109283381005D-01,0.8210327652337982D-01,
     *  0.8648799326620257D-01,0.9337138839016661D-01,
     1  0.1057091486383661D+00,0.1299986788041655D+00,
     2  0.1750784030217495D+00,0.2441623561476874D+00,
     3  0.3332037810245296D+00,0.4354695087535625D+00,
     4  0.5437666710470152D+00,0.6512343855452133D+00,
     5  0.7517263142133977D+00,0.8399773037163661D+00,
     6  0.9116664622988446D+00,0.9634499987835750D+00,
     *  0.9929988270130290D+00/

        data w2/
     1  0.4882734133843424D-06,0.1705467067529887D-04,
     2  0.1798438094661238D-03,0.9739619730495464D-03,
     3  0.3348244273080515D-02,0.8028239794437749D-02,
     4  0.1356426290909399D-01,0.1513716474932675D-01,
     5  0.1132919551407218D-01,0.7233155809463682D-02,
     6  0.4774905990583253D-02,0.3432689818235191D-02,
     7  0.2700329977348582D-02,0.2329552355327758D-02,
     8  0.2218951956375982D-02,0.2348947937734485D-02,
     9  0.2768428635994767D-02,0.3631946313140370D-02,
     *  0.5333001071547002D-02,0.8889805101750872D-02,
     1  0.1685816986683378D-01,0.3339149774029169D-01,
     2  0.5731225536672224D-01,0.8008102851132917D-01,
     3  0.9685078920859041D-01,0.1064701004799069D+00,
     4  0.1089746513839503D+00,0.1049338999537466D+00,
     5  0.9517155443230539D-01,0.8060978887968372D-01,
     6  0.6221069505997419D-01,0.4096682730925852D-01,
     *  0.1792857087328952D-01/

        data r3/
     1  0.3840435951943103D+00,0.4628878857324323D+00,
     2  0.3055899252035182D+00,0.5539710788516706D+00,
     3  0.2321351553786520D+00,0.6538492276560936D+00,
     4  0.7116047051272663D-02,0.2774119226499005D-02,
     5  0.1746553479723241D+00,0.7515177476035291D+00,
     6  0.2235971052636181D-01,0.1362282271302786D+00,
     7  0.8390437144322692D+00,0.7807301722322991D-01,
     8  0.3352716437194165D-01,0.1347819778034241D-01,
     9  0.7358672460107079D-01,0.9108793690221961D+00,
     *  0.1129409130653321D+00,0.4410501673447665D-01,
     1  0.8320656520731211D-01,0.6924167142908240D-01,
     2  0.7451519965917604D-03,0.9630572318494646D+00,
     3  0.1271682892805219D-03,0.1124799887508182D-04,
     4  0.6460998637077341D-01,0.9885343456785322D-01,
     5  0.8971874768616508D-01,0.2983187950129386D-06,
     6  0.5261337303362807D-01,0.5924627984728199D-01,
     *  0.9929170705363775D+00/

        data w3/
     1  0.7730444724490451D-01,0.8358238198570267D-01,
     2  0.7821084656241578D-01,0.9745492670287210D-01,
     3  0.6670917545537008D-01,0.1003640566115185D+00,
     4  0.5475751144772356D-02,0.3110242112401522D-02,
     5  0.4776132049285158D-01,0.9368010756048750D-01,
     6  0.1042618787330193D-01,0.2984542895048583D-01,
     7  0.8045911771224931D-01,0.4713184004718013D-02,
     8  0.1132162159414625D-01,0.7357224246291020D-02,
     9  0.4341482134254216D-02,0.6255485792378522D-01,
     *  0.1778645736610742D-01,0.9599853814331200D-02,
     1  0.5671405485289632D-02,0.4418423922245886D-02,
     2  0.1132871919026340D-02,0.4135787268745829D-01,
     3  0.2603688320078640D-03,0.3239563819948927D-04,
     4  0.4917839611970904D-02,0.1109059459604067D-01,
     5  0.7551049930438245D-02,0.1397027742693295D-05,
     6  0.7474259288060625D-02,0.5898714501981332D-02,
     *  0.1813413506657100D-01/

        data r4/
     1  0.1929236652276537D-06,0.7682320621847399D-05,
     2  0.9047262653993485D-04,0.5476451763327760D-03,
     3  0.2118753826806249D-02,0.5862375335642816D-02,
     4  0.1241460611751653D-01,0.2122053764129703D-01,
     5  0.3086600098709591D-01,0.4016525258250673D-01,
     6  0.4864584950925994D-01,0.5636023525664242D-01,
     7  0.6363184486029289D-01,0.7092606055534512D-01,
     8  0.7881546586998496D-01,0.8804792336599148D-01,
     9  0.9972215777462291D-01,0.1155675699821715D+00,
     *  0.1383072985168249D+00,0.1718362386871814D+00,
     1  0.2204699509605455D+00,0.2868720488286188D+00,
     2  0.3701560305174174D+00,0.4659680477076250D+00,
     3  0.5680400455357793D+00,0.6697402186426224D+00,
     4  0.7649874845879243D+00,0.8486222436795989D+00,
     5  0.9164966900543817D+00,0.9654687483943414D+00,
     *  0.9933878931500266D+00/

        data w4/
     1  0.9155381610132737D-06,0.2248810200718057D-04,
     2  0.1884293682819405D-03,0.8500718548320260D-03,
     3  0.2482283250058392D-02,0.5127801945838033D-02,
     4  0.7877677479120556D-02,0.9475421768549927D-02,
     5  0.9610432270838642D-02,0.8914412608859277D-02,
     6  0.8059974602556783D-02,0.7424665542172251D-02,
     7  0.7196785148358898D-02,0.7485938129818069D-02,
     8  0.8414565173985600D-02,0.1022885418877437D-01,
     9  0.1340133383495579D-01,0.1873838626642240D-01,
     *  0.2740066876594350D-01,0.4041533283924385D-01,
     1  0.5733381895011317D-01,0.7531881071163717D-01,
     2  0.9049330241247616D-01,0.1000589438679529D+00,
     3  0.1029649803402803D+00,0.9941638467383509D-01,
     4  0.9021714270440253D-01,0.7636135112526522D-01,
     5  0.5886422880561256D-01,0.3872086995348533D-01,
     *  0.1693372777616109D-01/

        data r5/
     1  0.1552931333595893D-06,0.5804872302953473D-05,
     2  0.6481146176897053D-04,0.3760677008416288D-03,
     3  0.1413680336312489D-02,0.3873670719058315D-02,
     4  0.8343370174487121D-02,0.1495087191362457D-01,
     5  0.2336788962620846D-01,0.3313207678429528D-01,
     6  0.4395813139700819D-01,0.5589370275294879D-01,
     7  0.6939453072172244D-01,0.8540788463825124D-01,
     8  0.1054977948031395D+00,0.1319832718006160D+00,
     9  0.1679430457935045D+00,0.2167763848345668D+00,
     *  0.2810700407314034D+00,0.3611437681657923D+00,
     1  0.4542919902092163D+00,0.5552990573763295D+00,
     2  0.6577176390745332D+00,0.7550607231968016D+00,
     3  0.8415206262846788D+00,0.9122849148495830D+00,
     4  0.9636437145949262D+00,0.9930296263918334D+00/

        data w5/
     1  0.7263923191098918D-06,0.1665317656313051D-04,
     2  0.1318691214004976D-03,0.5701484116697442D-03,
     3  0.1627303579513778D-02,0.3397625894463626D-02,
     4  0.5562912098508680D-02,0.7590922738691022D-02,
     5  0.9160565105108916D-02,0.1031794977464762D-01,
     6  0.1134107321120592D-01,0.1260470890589361D-01,
     7  0.1455269431278949D-01,0.1773207259573070D-01,
     8  0.2283336801043853D-01,0.3065713260396262D-01,
     9  0.4184861746310374D-01,0.5627646315108501D-01,
     *  0.7237852473449820D-01,0.8730772312875685D-01,
     1  0.9809573417044978D-01,0.1028231266617558D+00,
     2  0.1009212189812172D+00,0.9279226638938967D-01,
     3  0.7932399864279166D-01,0.6158581242873239D-01,
     4  0.4070281333906094D-01,0.1784597497625180D-01/

        data r6/
     1  0.1837266581404723D-06,0.6621703257044589D-05,
     2  0.7137730547714612D-04,0.4014351449846677D-03,
     3  0.1476169839252338D-02,0.4024413254791613D-02,
     4  0.8836960869824831D-02,0.1659770040052357D-01,
     5  0.2791250765167691D-01,0.4354660007235751D-01,
     6  0.6476361250405618D-01,0.9361628178514843D-01,
     7  0.1329729007497613D+00,0.1859901792327477D+00,
     8  0.2549135716708809D+00,0.3396841085690593D+00,
     9  0.4372930209167346D+00,0.5423399075992096D+00,
     *  0.6482864777761108D+00,0.7486093199181302D+00,
     1  0.8374891310408979D+00,0.9101094928531647D+00,
     2  0.9627571636332322D+00,0.9928612026915437D+00/

        data w6/
     1  0.8526276533891775D-06,0.1875769055834648D-04,
     2  0.1428143178414253D-03,0.5974681616214515D-03,
     3  0.1677161487849099D-02,0.3554444432165495D-02,
     4  0.6184131255260617D-02,0.9434078891707467D-02,
     5  0.1331708079006538D-01,0.1815637648369950D-01,
     6  0.2461855311266178D-01,0.3357309406017619D-01,
     7  0.4568739679716833D-01,0.6074738423764691D-01,
     8  0.7710168275060508D-01,0.9192701184792176D-01,
     9  0.1023676030766112D+00,0.1066164349632099D+00,
     *  0.1041783197270617D+00,0.9549302401683349D-01,
     1  0.8146249940153088D-01,0.6315890329891557D-01,
     2  0.4170682471824018D-01,0.1827810185299490D-01/
c        
c        
c         %%%%%%%%%%%%%%%%%%%%%

ccc        call prin2('alpha = *',alpha,1)
        dtwo = 2
        if (alpha .ge. 1/dtwo) then
        do 1000 i=1,nqs(6)
        roots(i) =r6(i)
        weights(i) = w6(i)
 1000 continue

        nquads = nqs(6)
        endif

        if (alpha .ge. dtwo**(-2)) then
        if (alpha .lt. dtwo**(-1)) then

        do 2000 i=1,nqs(5)
        roots(i) =r5(i)
        weights(i) = w5(i)
 2000 continue

        nquads = nqs(5)

        endif
        endif


        if (alpha .ge. dtwo**(-3)) then
        if (alpha .lt. dtwo**(-2)) then

        do 3000 i=1,nqs(4)
        roots(i) =r4(i)
        weights(i) = w4(i)
 3000 continue

        nquads = nqs(4)

        endif
        endif


        if (alpha .ge. dtwo**(-4)) then
        if (alpha .lt. dtwo**(-3)) then

        do 4000 i=1,nqs(3)
        roots(i) =r3(i)
        weights(i) = w3(i)
 4000 continue

        nquads = nqs(3)

        endif
        endif


        if (alpha .ge. dtwo**(-5)) then
        if (alpha .lt. dtwo**(-4)) then

        do 5000 i=1,nqs(2)
        roots(i) =r2(i)
        weights(i) = w2(i)
 5000 continue

        nquads = nqs(2)

        endif
        endif


        if (alpha .lt. dtwo**(-5)) then

        do 6000 i=1,nqs(1)
        roots(i) =r1(i)
        weights(i) = w1(i)
 6000 continue

        nquads = nqs(1)

        endif

        return
        end
c
c
c
c
c
        subroutine lap_load_quads16(alpha,roots,weights,nquads)
        implicit real *8 (a-h,o-z)
        dimension r1(33),r2(33),r3(32),r4(30),r5(27),r6(24),
     1          w1(33),w2(33),w3(32),w4(30),w5(27),w6(24),nqs(6),
     2          roots(1),weights(1)

        data nqs/33,33,32,30,27,24/


        data r1/
     1  0.4996742557661396D-06,0.1919495916763985D-04,
     2  0.2192608743320138D-03,0.1293130565665469D-02,
     3  0.4899292373715333D-02,0.1343100656215535D-01,
     4  0.2895183819938661D-01,0.5137334042233771D-01,
     5  0.7364106021160841D-01,0.8711853610596712D-01,
     6  0.9396809844387508D-01,0.9789130215882889D-01,
     7  0.1005023396840622D+00,0.1024695869538229D+00,
     8  0.1041228506064478D+00,0.1056692645296571D+00,
     9  0.1072827390738619D+00,0.1091629856680879D+00,
     *  0.1116177365295181D+00,0.1152594922218987D+00,
     1  0.1216186753872256D+00,0.1352392341953818D+00,
     2  0.1676917012248405D+00,0.2294057469793434D+00,
     3  0.3160489713660002D+00,0.4186873721055752D+00,
     4  0.5289707377191591D+00,0.6392865641382565D+00,
     5  0.7429178616592162D+00,0.8341709226209254D+00,
     6  0.9084154031418082D+00,0.9620932483042569D+00,
     *  0.9927378100288193D+00/

        data w1/
     1  0.2351662834772341D-05,0.5555504272575146D-04,
     2  0.4504822642175805D-03,0.1976405801512732D-02,
     3  0.5656313350788543D-02,0.1176915054080802D-01,
     4  0.1935650622995528D-01,0.2430997147815224D-01,
     5  0.1832988411366188D-01,0.9339465287475442D-02,
     6  0.4983460926019304D-02,0.3107952427551423D-02,
     7  0.2216816019921339D-02,0.1769587452885937D-02,
     8  0.1570154733429341D-02,0.1550778536318736D-02,
     9  0.1708254978041186D-02,0.2100778563858108D-02,
     *  0.2902638671721333D-02,0.4607378219998243D-02,
     1  0.8776297247617912D-02,0.2045606048896148D-01,
     2  0.4666269919247720D-01,0.7566999209299619D-01,
     3  0.9609251904426096D-01,0.1077944856352829D+00,
     4  0.1114998327565505D+00,0.1080102062984031D+00,
     5  0.9830294548633320D-01,0.8343138383845408D-01,
     6  0.6446434831965149D-01,0.4247902077045648D-01,
     *  0.1859632252667728D-01/

        data r2/
     1  0.2318136659976222D-06,0.9572027317699132D-05,
     2  0.1181168848932009D-03,0.7584064058397523D-03,
     3  0.3157892862076491D-02,0.9552303249316230D-02,
     4  0.2229490796395620D-01,0.4115773784206448D-01,
     5  0.6085316881963827D-01,0.7563554811477437D-01,
     6  0.8521742792158959D-01,0.9159125259863445D-01,
     7  0.9620157462747097D-01,0.9985181994393753D-01,
     8  0.1030179799090868D+00,0.1060401160951144D+00,
     9  0.1092300668975651D+00,0.1129603394570100D+00,
     *  0.1177974373573766D+00,0.1247938317166186D+00,
     1  0.1362101704559952D+00,0.1570530777319798D+00,
     2  0.1958351913109983D+00,0.2586636472868279D+00,
     3  0.3431064139784270D+00,0.4421459579227672D+00,
     4  0.5481654506580353D+00,0.6540355792486208D+00,
     5  0.7534262399072994D+00,0.8409357142137357D+00,
     6  0.9121425854263322D+00,0.9636328378960288D+00,
     *  0.9930324060369954D+00/

        data w2/
     1  0.1107453637350373D-05,0.2841401477726120D-04,
     2  0.2522064764330609D-03,0.1227136746086083D-02,
     3  0.3956979573995464D-02,0.9270752105867763D-02,
     4  0.1624192552255286D-01,0.2051982093277397D-01,
     5  0.1775609025656160D-01,0.1186957558490150D-01,
     6  0.7667174689708967D-02,0.5316703328323621D-02,
     7  0.4032173814235738D-02,0.3344655554644877D-02,
     8  0.3042694607424962D-02,0.3052200696824629D-02,
     9  0.3387794172020474D-02,0.4162881431125914D-02,
     *  0.5676002994504264D-02,0.8667337782618272D-02,
     1  0.1496566866211289D-01,0.2824526271138226D-01,
     2  0.5049738698566199D-01,0.7463071765145096D-01,
     3  0.9302358682183489D-01,0.1037735717575597D+00,
     4  0.1070726827807662D+00,0.1036105892551037D+00,
     5  0.9426970843862427D-01,0.8001194006600326D-01,
     6  0.6183235102743543D-01,0.4075117871519638D-01,
     *  0.1784172738784951D-01/

        data r3/
     1  0.2709014631443115D-06,0.1039108417334220D-04,
     2  0.1183283774964105D-03,0.6817961353083182D-03,
     3  0.2376795533741046D-02,0.5776352903968398D-02,
     4  0.1245149110351927D-01,0.2426405564109015D-01,
     5  0.3996183839593585D-01,0.5578827564418930D-01,
     6  0.6903105033611473D-01,0.7938850824746093D-01,
     7  0.8760725565745659D-01,0.9448526964437704D-01,
     8  0.1006784572481389D+00,0.1067674837190724D+00,
     9  0.1133691835885528D+00,0.1212832024433256D+00,
     *  0.1317416382254118D+00,0.1468835279238169D+00,
     1  0.1705133653360781D+00,0.2084403196756051D+00,
     2  0.2662331310151964D+00,0.3448773047224028D+00,
     3  0.4398858137111621D+00,0.5440088093729741D+00,
     4  0.6495859152027034D+00,0.7496336104376778D+00,
     5  0.8382190626640051D+00,0.9105442988621739D+00,
     6  0.9629463222620056D+00,0.9928984148668204D+00/

        data w3/
     1  0.1274044991928472D-05,0.3006853622705753D-04,
     2  0.2421098427070259D-03,0.1011141341327533D-02,
     3  0.2451817477689943D-02,0.4608331017186369D-02,
     4  0.9134441179389750D-02,0.1426959258054180D-01,
     5  0.1642104238107433D-01,0.1475813837644664D-01,
     6  0.1171073856490856D-01,0.9144296112197830D-02,
     7  0.7427866386457663D-02,0.6435982745028052D-02,
     8  0.6045333613277483D-02,0.6233241536578590D-02,
     9  0.7099574010013665D-02,0.8926263123602217D-02,
     *  0.1233353139060838D-01,0.1857217772235581D-01,
     1  0.2968993942084210D-01,0.4718711320843563D-01,
     2  0.6857486206650355D-01,0.8790718309686234D-01,
     3  0.1008621416827698D+00,0.1060907153547177D+00,
     4  0.1038943126485434D+00,0.9521132171999308D-01,
     5  0.8116008207665797D-01,0.6287998860971047D-01,
     6  0.4150198316432317D-01,0.1818339496803007D-01/

        data r4/
     1  0.3492887557284512D-06,0.1302525001757679D-04,
     2  0.1450815026226367D-03,0.8376203351991228D-03,
     3  0.3113893750247026D-02,0.8338955027509692D-02,
     4  0.1722754920651789D-01,0.2896320163906717D-01,
     5  0.4172720216874723D-01,0.5410069460876225D-01,
     6  0.6564561678699850D-01,0.7649794955318201D-01,
     7  0.8695019678582991D-01,0.9741493092839936D-01,
     8  0.1085137759925432D+00,0.1211694379880311D+00,
     9  0.1367528087002015D+00,0.1573313147586720D+00,
     *  0.1859654151768022D+00,0.2267042474306684D+00,
     1  0.2835457705201054D+00,0.3582180975689769D+00,
     2  0.4485076895436095D+00,0.5487744543580960D+00,
     3  0.6518091230838458D+00,0.7504415245091465D+00,
     4  0.8383784112259047D+00,0.9104950192663051D+00,
     5  0.9628898201900457D+00,0.9928839590470912D+00/

        data w4/
     1  0.1632636135025546D-05,0.3734212089973439D-04,
     2  0.2947195794584544D-03,0.1263759106415782D-02,
     3  0.3535012311728417D-02,0.7053512204093771D-02,
     4  0.1057718001081140D-01,0.1256180263871093D-01,
     5  0.1272018195776411D-01,0.1196383857118651D-01,
     6  0.1115647082128197D-01,0.1059644326998585D-01,
     7  0.1037514733711789D-01,0.1065658765085540D-01,
     8  0.1169255961179946D-01,0.1384306130894192D-01,
     9  0.1766260682438838D-01,0.2400246947378829D-01,
     *  0.3395655497032173D-01,0.4822954485827261D-01,
     1  0.6577496631113546D-01,0.8318875983156776D-01,
     2  0.9641485776284359D-01,0.1028894367648192D+00,
     3  0.1019717310866145D+00,0.9424013467670065D-01,
     4  0.8077899430567253D-01,0.6280419194976924D-01,
     5  0.4153822586049163D-01,0.1821827418642780D-01/

        data r5/
     1  0.2394629806012706D-06,0.9089991978211130D-05,
     2  0.1018623224860555D-03,0.5855835525470437D-03,
     3  0.2158104714610989D-02,0.5768313269638481D-02,
     4  0.1213438516652844D-01,0.2136024600150313D-01,
     5  0.3302672893943320D-01,0.4658260046944923D-01,
     6  0.6168498893086425D-01,0.7839471875517791D-01,
     7  0.9729875806555243D-01,0.1196211027520224D+00,
     8  0.1473271361620061D+00,0.1831433564562904D+00,
     9  0.2302778430850043D+00,0.2915599925044619D+00,
     *  0.3680591692863713D+00,0.4579642371680875D+00,
     1  0.5566255863262489D+00,0.6576933383162853D+00,
     2  0.7544779231045556D+00,0.8408811587740456D+00,
     3  0.9118318269537742D+00,0.9634299461745898D+00,
     *  0.9929860380886632D+00/

        data w5/
     1  0.1124922646567085D-05,0.2616857497743707D-04,
     2  0.2069584415043777D-03,0.8788516555989403D-03,
     3  0.2433721196143996D-02,0.4913981066582643D-02,
     4  0.7833282204101046D-02,0.1054100527570815D-01,
     5  0.1269428986984079D-01,0.1435730390275934D-01,
     6  0.1585571876815901D-01,0.1765973971003825D-01,
     7  0.2035044454906395D-01,0.2462096718722985D-01,
     8  0.3124932964573770D-01,0.4093200316748068D-01,
     9  0.5383368180259389D-01,0.6893031985374590D-01,
     *  0.8376519194761915D-01,0.9524837771890751D-01,
     1  0.1009916675324606D+00,0.1000134280545793D+00,
     2  0.9253272704395447D-01,0.7942658013778682D-01,
     3  0.6182438302859597D-01,0.4092271105852650D-01,
     *  0.1795604168365711D-01/

        data r6/
     1  0.1738034326416008D-06,0.6835108616748489D-05,
     2  0.7857081868672412D-04,0.4623867880689919D-03,
     3  0.1753936947666065D-02,0.4885677444088860D-02,
     4  0.1089439094402132D-01,0.2068621187186298D-01,
     5  0.3502327149494967D-01,0.5475347866755926D-01,
     6  0.8114906615951964D-01,0.1161560044348819D+00,
     7  0.1623076692165838D+00,0.2220664595915370D+00,
     8  0.2966419456101960D+00,0.3848860840019453D+00,
     9  0.4830529364395184D+00,0.5855744225293911D+00,
     *  0.6862526678050489D+00,0.7792359177449232D+00,
     1  0.8595885708478710D+00,0.9235746199852068D+00,
     2  0.9687833740630569D+00,0.9940731965266603D+00/

        data w6/
     1  0.8240127749016774D-06,0.1987547794948667D-04,
     2  0.1614594860359048D-03,0.7054565607502984D-03,
     3  0.2037606271856667D-02,0.4403218548673904D-02,
     4  0.7765459781415567D-02,0.1194100830516354D-01,
     5  0.1686832052727821D-01,0.2279938483333690D-01,
     6  0.3031687605004689D-01,0.4013014252748967D-01,
     7  0.5260485707550948D-01,0.6713506973013861D-01,
     8  0.8182898371464977D-01,0.9402020384311163D-01,
     9  0.1013638985815366D+00,0.1026294851321915D+00,
     *  0.9774263381055503D-01,0.8739691191495233D-01,
     1  0.7268356538488463D-01,0.5488725482647525D-01,
     2  0.3534590190345732D-01,0.1521160169976591D-01/
c        
c         %%%%%%%%%%%%%%%%%%%%%

ccc        call prin2('alpha = *',alpha,1)
        dtwo = 2
        if (alpha .ge. 1/dtwo) then
        do 1000 i=1,nqs(6)
        roots(i) =r6(i)
        weights(i) = w6(i)
 1000 continue

        nquads = nqs(6)
        endif

        if (alpha .ge. dtwo**(-2)) then
        if (alpha .lt. dtwo**(-1)) then

        do 2000 i=1,nqs(5)
        roots(i) =r5(i)
        weights(i) = w5(i)
 2000 continue

        nquads = nqs(5)

        endif
        endif


        if (alpha .ge. dtwo**(-3)) then
        if (alpha .lt. dtwo**(-2)) then

        do 3000 i=1,nqs(4)
        roots(i) =r4(i)
        weights(i) = w4(i)
 3000 continue

        nquads = nqs(4)

        endif
        endif


        if (alpha .ge. dtwo**(-4)) then
        if (alpha .lt. dtwo**(-3)) then

        do 4000 i=1,nqs(3)
        roots(i) =r3(i)
        weights(i) = w3(i)
 4000 continue

        nquads = nqs(3)

        endif
        endif


        if (alpha .ge. dtwo**(-5)) then
        if (alpha .lt. dtwo**(-4)) then

        do 5000 i=1,nqs(2)
        roots(i) =r2(i)
        weights(i) = w2(i)
 5000 continue

        nquads = nqs(2)

        endif
        endif


        if (alpha .lt. dtwo**(-5)) then

        do 6000 i=1,nqs(1)
        roots(i) =r1(i)
        weights(i) = w1(i)
 6000 continue

        nquads = nqs(1)

        endif

        return
        end
c
c
c
c
c
        subroutine lap_load_quads17(alpha,roots,weights,nquads)
        implicit real *8 (a-h,o-z)
        dimension r1(33),r2(32),r3(32),r4(30),r5(27),r6(23),
     1          w1(33),w2(32),w3(32),w4(30),w5(27),w6(23),nqs(6),
     2          roots(1),weights(1)

        data nqs/33,32,32,30,27,23/


        data r1/
     1  0.5341644326459973D-06,0.2080259689413553D-04,
     2  0.2413304080734940D-03,0.1450334953705445D-02,
     3  0.5623712145430014D-02,0.1577234986070453D-01,
     4  0.3425242275918581D-01,0.6054428049221170D-01,
     5  0.8955053408265130D-01,0.1113842532636569D+00,
     6  0.1233980697050673D+00,0.1300271885489452D+00,
     7  0.1342059586544216D+00,0.1371995170470130D+00,
     8  0.1395945903714170D+00,0.1417223587701506D+00,
     9  0.1438247238994059D+00,0.1461396844927474D+00,
     *  0.1489830408441599D+00,0.1529017457328475D+00,
     1  0.1590848844839435D+00,0.1706399260795219D+00,
     2  0.1958700476258023D+00,0.2476515679325035D+00,
     3  0.3276140327326253D+00,0.4263576860769740D+00,
     4  0.5341759631235443D+00,0.6427946732174141D+00,
     5  0.7451975542294879D+00,0.8355473288230311D+00,
     6  0.9091414507632254D+00,0.9623849857754185D+00,
     *  0.9927928280455098D+00/

        data w1/
     1  0.2521057416413655D-05,0.6050641143548509D-04,
     2  0.4998952048937986D-03,0.2247175260195686D-02,
     3  0.6631611259724470D-02,0.1410349186455722D-01,
     4  0.2278132006678373D-01,0.2902969565696440D-01,
     5  0.2700044527154793D-01,0.1630855837526407D-01,
     6  0.8624200254099565D-02,0.5098948951295508D-02,
     7  0.3453020452927311D-02,0.2625779877885682D-02,
     8  0.2217145196567627D-02,0.2077437562131360D-02,
     9  0.2165743565342822D-02,0.2514599288463203D-02,
     *  0.3257571424554241D-02,0.4760482535516396D-02,
     1  0.8067988125731637D-02,0.1635565631508475D-01,
     2  0.3668048745863456D-01,0.6705575846465115D-01,
     3  0.9109562765095698D-01,0.1047722992502094D+00,
     4  0.1094887841921184D+00,0.1065802502156748D+00,
     5  0.9725338294268249D-01,0.8266288261773719D-01,
     6  0.6392624543309603D-01,0.4214556266155552D-01,
     *  0.1845492513430021D-01/

        data r2/
     1  0.4863182234993468D-06,0.1914091817912507D-04,
     2  0.2243748745212372D-03,0.1362421684016608D-02,
     3  0.5337205819608887D-02,0.1510761219463131D-01,
     4  0.3287834609887068D-01,0.5692653070376204D-01,
     5  0.8107476785503066D-01,0.1001747011744118D+00,
     6  0.1135112809429491D+00,0.1226809304489003D+00,
     7  0.1293321462704049D+00,0.1345589893337655D+00,
     8  0.1390578101356622D+00,0.1433397484883206D+00,
     9  0.1478799948006080D+00,0.1532530665556582D+00,
     *  0.1603422479795917D+00,0.1707871892011665D+00,
     1  0.1879667947469951D+00,0.2183999803341702D+00,
     2  0.2703238163570705D+00,0.3463322212732643D+00,
     3  0.4407723342661218D+00,0.5449622268955478D+00,
     4  0.6505960113046318D+00,0.7505406900824923D+00,
     5  0.8389008405964837D+00,0.9109598153650800D+00,
     6  0.9631289474752926D+00,0.9929344960269712D+00/

        data w2/
     1  0.2300485593730921D-05,0.5587355492424739D-04,
     2  0.4670953412879619D-03,0.2124989259693544D-02,
     3  0.6347613370373519D-02,0.1362267045767615D-01,
     4  0.2165971559511681D-01,0.2528199237054148D-01,
     5  0.2209309571919086D-01,0.1604184263264308D-01,
     6  0.1093794561134998D-01,0.7682572258789628D-02,
     7  0.5798058664288938D-02,0.4768062060870383D-02,
     8  0.4312427626691387D-02,0.4328583891927721D-02,
     9  0.4844519698801981D-02,0.6042131627725002D-02,
     *  0.8393939916896757D-02,0.1302979179508059D-01,
     1  0.2241846762200929D-01,0.3998532024692086D-01,
     2  0.6434775689270058D-01,0.8660125696518120D-01,
     3  0.1007815651095346D+00,0.1062052237719147D+00,
     4  0.1038734144780081D+00,0.9503457454137300D-01,
     5  0.8090050183946185D-01,0.6261962100237628D-01,
     6  0.4130542509749863D-01,0.1809165049355711D-01/

        data r3/
     1  0.3187917067350701D-06,0.1267212303208947D-04,
     2  0.1507541703522016D-03,0.9319540228182983D-03,
     3  0.3721867321105569D-02,0.1074660664726466D-01,
     4  0.2392816286815917D-01,0.4279115245870409D-01,
     5  0.6372544648670495D-01,0.8266100238638606D-01,
     6  0.9795417177455321D-01,0.1099907210047717D+00,
     7  0.1197146966672260D+00,0.1279873440266831D+00,
     8  0.1355124849305651D+00,0.1429104547681320D+00,
     9  0.1508193025794229D+00,0.1600259664142258D+00,
     *  0.1716687267549349D+00,0.1875845050105675D+00,
     1  0.2108526388180730D+00,0.2462062000200374D+00,
     2  0.2989274948550104D+00,0.3713264489797008D+00,
     3  0.4603614822830430D+00,0.5593249339805803D+00,
     4  0.6605904946649091D+00,0.7571010258432840D+00,
     5  0.8428625282142978D+00,0.9130424951347154D+00,
     6  0.9639626591372223D+00,0.9930913298556694D+00/

        data w3/
     1  0.1510292301544218D-05,0.3714968889449322D-04,
     2  0.3163482008190021D-03,0.1471159690828208D-02,
     3  0.4499064417723045D-02,0.9911528725680553D-02,
     4  0.1639331882090504D-01,0.2067619247904559D-01,
     5  0.2045260476448835D-01,0.1716213498889875D-01,
     6  0.1351717687604208D-01,0.1072222464318198D-01,
     7  0.8869523225099273D-02,0.7791648871616725D-02,
     8  0.7360369589007508D-02,0.7539818089186817D-02,
     9  0.8405840631386279D-02,0.1019164055612330D-01,
     *  0.1339058995461519D-01,0.1894394737353573D-01,
     1  0.2839419383212472D-01,0.4326309507742979D-01,
     2  0.6261573606836887D-01,0.8163404225275824D-01,
     3  0.9526715178750709D-01,0.1013657553235637D+00,
     4  0.9998243858360329D-01,0.9203859655195294D-01,
     5  0.7867930598268827D-01,0.6106783299567204D-01,
     6  0.4034972431874653D-01,0.1768833534620506D-01/

        data r4/
     1  0.2242476623189804D-06,0.8670634169580458D-05,
     2  0.9974953588772308D-04,0.5941405866900183D-03,
     3  0.2290409078095081D-02,0.6459892322123761D-02,
     4  0.1439492995936701D-01,0.2662883003649958D-01,
     5  0.4233099767678421D-01,0.5965288063402789D-01,
     6  0.7687260440012268D-01,0.9312701402571970D-01,
     7  0.1083511784040336D+00,0.1229781808403587D+00,
     8  0.1377664504004976D+00,0.1537822343168577D+00,
     9  0.1724907503018231D+00,0.1959263579639861D+00,
     *  0.2268689364001555D+00,0.2687632523566228D+00,
     1  0.3248843976622307D+00,0.3966354339156060D+00,
     2  0.4821398005173795D+00,0.5764818483732066D+00,
     3  0.6732006524830701D+00,0.7657368886941216D+00,
     4  0.8482517497220246D+00,0.9159459880832410D+00,
     5  0.9651445860029629D+00,0.9933157328800527D+00/

        data w4/
     1  0.1056878790325632D-05,0.2515208535732013D-04,
     2  0.2057335299997433D-03,0.9155951433879135D-03,
     3  0.2696713007556252D-02,0.5871665901193633D-02,
     4  0.1009993704011102D-01,0.1422361353623742D-01,
     5  0.1684928567441659D-01,0.1749457742191923D-01,
     6  0.1680399766280648D-01,0.1569981450779353D-01,
     7  0.1482350041113184D-01,0.1455706687509273D-01,
     8  0.1519561603790800D-01,0.1707720894483274D-01,
     9  0.2067538611908134D-01,0.2665909003881715D-01,
     *  0.3581319390055886D-01,0.4855288885913735D-01,
     1  0.6395125750000532D-01,0.7923240436563154D-01,
     2  0.9092820014062278D-01,0.9665376004395525D-01,
     3  0.9567756401675097D-01,0.8841716147302297D-01,
     4  0.7581164764091203D-01,0.5896422380479641D-01,
     5  0.3901018501187354D-01,0.1711250242629972D-01/

        data r5/
     1  0.1958252377975966D-06,0.7672553509504654D-05,
     2  0.8899433065060761D-04,0.5320572121539478D-03,
     3  0.2050187945886540D-02,0.5755295073621992D-02,
     4  0.1273287990523775D-01,0.2348925075971767D-01,
     5  0.3775893923334644D-01,0.5482187546377149D-01,
     6  0.7401676667760247D-01,0.9511285256328147D-01,
     7  0.1184837509676397D+00,0.1451836151276628D+00,
     8  0.1769810524336688D+00,0.2162920040784177D+00,
     9  0.2658438912013901D+00,0.3278936271904361D+00,
     *  0.4031057353974747D+00,0.4896969337175238D+00,
     1  0.5834978854110283D+00,0.6788741151399449D+00,
     2  0.7698478258549607D+00,0.8509039857305114D+00,
     3  0.9174028250066891D+00,0.9657440476237282D+00,
     *  0.9934301377301074D+00/

        data w5/
     1  0.9260482864598482D-06,0.2233835382769281D-04,
     2  0.1840734725174147D-03,0.8207963380258188D-03,
     3  0.2408617241911023D-02,0.5193188438557198D-02,
     4  0.8849998213386552D-02,0.1261221056335762D-01,
     5  0.1579764932774403D-01,0.1821177293568069D-01,
     6  0.2013698364026590D-01,0.2211633697498098D-01,
     7  0.2480044537222101D-01,0.2889474654815084D-01,
     8  0.3510783051119464D-01,0.4398177050960346D-01,
     9  0.5551849839376699D-01,0.6870797698430011D-01,
     *  0.8141777994737676D-01,0.9105239514766380D-01,
     1  0.9559039969018637D-01,0.9415227737718189D-01,
     2  0.8686953992654837D-01,0.7446710176303408D-01,
     3  0.5792972122832986D-01,0.3833535639860405D-01,
     *  0.1681926865329639D-01/

        data r6/
     1  0.2445392023226351D-06,0.9353256014397869D-05,
     2  0.1053320550860552D-03,0.6103757577034545D-03,
     3  0.2290235229838391D-02,0.6337258182075233D-02,
     4  0.1408472107680578D-01,0.2670755148724773D-01,
     5  0.4516908078930166D-01,0.7044725909150846D-01,
     6  0.1038784192560930D+00,0.1473572776515585D+00,
     7  0.2031020834782488D+00,0.2728209907054232D+00,
     8  0.3565404002988403D+00,0.4518248982370175D+00,
     9  0.5539693033840179D+00,0.6569736534979403D+00,
     *  0.7546310141667457D+00,0.8412760473434945D+00,
     1  0.9121587416397627D+00,0.9635932584112680D+00,
     *  0.9930201044356364D+00/

        data w6/
     1  0.1151315573093624D-05,0.2697739432372816D-04,
     2  0.2145752273187816D-03,0.9232277275342564D-03,
     3  0.2640767320650855D-02,0.5681198312478554D-02,
     4  0.1001047956678486D-01,0.1539153608958869D-01,
     5  0.2168814623949932D-01,0.2908685658400582D-01,
     6  0.3809462423122850D-01,0.4924723944493187D-01,
     7  0.6254863531361428D-01,0.7690831973543967D-01,
     8  0.9012076590347216D-01,0.9964581546582705D-01,
     9  0.1036237413274517D+00,0.1013380449650607D+00,
     *  0.9302666576037917D-01,0.7946957692775363D-01,
     1  0.6168070173676490D-01,0.4076051710802233D-01,
     *  0.1787043630229600D-01/
c        
c         %%%%%%%%%%%%%%%%%%%%%

ccc        call prin2('alpha = *',alpha,1)
        dtwo = 2
        if (alpha .ge. 1/dtwo) then
        do 1000 i=1,nqs(6)
        roots(i) =r6(i)
        weights(i) = w6(i)
 1000 continue

        nquads = nqs(6)
        endif

        if (alpha .ge. dtwo**(-2)) then
        if (alpha .lt. dtwo**(-1)) then

        do 2000 i=1,nqs(5)
        roots(i) =r5(i)
        weights(i) = w5(i)
 2000 continue

        nquads = nqs(5)

        endif
        endif


        if (alpha .ge. dtwo**(-3)) then
        if (alpha .lt. dtwo**(-2)) then

        do 3000 i=1,nqs(4)
        roots(i) =r4(i)
        weights(i) = w4(i)
 3000 continue

        nquads = nqs(4)

        endif
        endif


        if (alpha .ge. dtwo**(-4)) then
        if (alpha .lt. dtwo**(-3)) then

        do 4000 i=1,nqs(3)
        roots(i) =r3(i)
        weights(i) = w3(i)
 4000 continue

        nquads = nqs(3)

        endif
        endif


        if (alpha .ge. dtwo**(-5)) then
        if (alpha .lt. dtwo**(-4)) then

        do 5000 i=1,nqs(2)
        roots(i) =r2(i)
        weights(i) = w2(i)
 5000 continue

        nquads = nqs(2)

        endif
        endif


        if (alpha .lt. dtwo**(-5)) then

        do 6000 i=1,nqs(1)
        roots(i) =r1(i)
        weights(i) = w1(i)
 6000 continue

        nquads = nqs(1)

        endif

        return
        end
c
c
c
c
c
        subroutine lap_load_quads18(alpha,roots,weights,nquads)
        implicit real *8 (a-h,o-z)
        dimension r1(33),r2(32),r3(31),r4(29),r5(26),r6(23),
     1          w1(33),w2(32),w3(31),w4(29),w5(26),w6(23),nqs(6),
     2          roots(1),weights(1)

        data nqs/33,32,31,29,26,23/


        data r1/
     1  0.3274971584867418D-06,0.1340901629237135D-04,
     2  0.1636866709826301D-03,0.1038781252834765D-02,
     3  0.4285702913418621D-02,0.1296254441212478D-01,
     4  0.3089053008348849D-01,0.6047782993366967D-01,
     5  0.9863351890204554D-01,0.1331852877430527D+00,
     6  0.1543091297033439D+00,0.1656313449817537D+00,
     7  0.1723307307863204D+00,0.1768622686722215D+00,
     8  0.1803055041251541D+00,0.1832118566596838D+00,
     9  0.1859348581080027D+00,0.1887740072089964D+00,
     *  0.1920715767595709D+00,0.1963455575186276D+00,
     1  0.2025935598351438D+00,0.2131565452561935D+00,
     2  0.2341336022104770D+00,0.2778292647137391D+00,
     3  0.3512223030537358D+00,0.4456352463934320D+00,
     4  0.5497333325004169D+00,0.6547494385009601D+00,
     5  0.7537426915335622D+00,0.8410665292783035D+00,
     6  0.9121902591428730D+00,0.9636465716634864D+00,
     *  0.9930344581981684D+00/

        data w1/
     1  0.1562234444004189D-05,0.3967302382956944D-04,
     2  0.3476066866643174D-03,0.1669130720621680D-02,
     3  0.5343991971532411D-02,0.1268215328047532D-01,
     4  0.2364215173683311D-01,0.3511905466974327D-01,
     5  0.3898337924065655D-01,0.2812240365546313D-01,
     6  0.1511042238818440D-01,0.8415107240218429D-02,
     7  0.5362615283198268D-02,0.3867740281789248D-02,
     8  0.3105859158008574D-02,0.2763645467362975D-02,
     9  0.2730825168810612D-02,0.3002561189139280D-02,
     *  0.3674580904804789D-02,0.5027617648192427D-02,
     1  0.7821890371628958D-02,0.1423904231735218D-01,
     2  0.2994722029493699D-01,0.5895329127433267D-01,
     3  0.8598618552137062D-01,0.1009163867355099D+00,
     4  0.1058388910648530D+00,0.1030442433494980D+00,
     5  0.9400366087775670D-01,0.7988963527574405D-01,
     6  0.6177992302956087D-01,0.4073141178157255D-01,
     *  0.1783613615591115D-01/

        data r2/
     1  0.4862653097223501D-06,0.1932517429341174D-04,
     2  0.2285644767109912D-03,0.1400148278030076D-02,
     3  0.5541522104595522D-02,0.1592139038637376D-01,
     4  0.3552164643128919D-01,0.6397079271087176D-01,
     5  0.9543533556939502D-01,0.1223129724707963D+00,
     6  0.1417742780997902D+00,0.1552472526452884D+00,
     7  0.1648848050536027D+00,0.1722553283199803D+00,
     8  0.1783670462998191D+00,0.1839190288567851D+00,
     9  0.1894961563930466D+00,0.1957195420477305D+00,
     *  0.2034249161799799D+00,0.2139847426181284D+00,
     1  0.2299741211987644D+00,0.2563028818296132D+00,
     2  0.3004494459383134D+00,0.3677792893429337D+00,
     3  0.4554498830714147D+00,0.5549582180716351D+00,
     4  0.6573267552603772D+00,0.7549139194884009D+00,
     5  0.8415435489547762D+00,0.9123556090404164D+00,
     6  0.9636904512886575D+00,0.9930404662784719D+00/

        data w2/
     1  0.2305214288451888D-05,0.5658971159608191D-04,
     2  0.4778078627557295D-03,0.2196661561085889D-02,
     3  0.6656357374966831D-02,0.1465645114430132D-01,
     4  0.2453286848581820D-01,0.3133147117381405D-01,
     5  0.3017452880520013D-01,0.2314530785849541D-01,
     6  0.1608562957875951D-01,0.1123572688122629D-01,
     7  0.8297883505567172D-02,0.6606605492834276D-02,
     8  0.5730190456483617D-02,0.5468918410412988D-02,
     9  0.5786134106473597D-02,0.6795032051456470D-02,
     *  0.8834045567810023D-02,0.1269432896369083D-01,
     1  0.2008534171429315D-01,0.3389712841269638D-01,
     2  0.5544190848551530D-01,0.7866584929742284D-01,
     3  0.9515920067648926D-01,0.1023392863208763D+00,
     4  0.1011247832646183D+00,0.9302311171994566D-01,
     5  0.7942807896292854D-01,0.6158693872839713D-01,
     6  0.4066417477239125D-01,0.1781935343738901D-01/

        data r3/
     1  0.3360016576838359D-06,0.1371208224751348D-04,
     2  0.1664192936153023D-03,0.1043922669258076D-02,
     3  0.4215123416422653D-02,0.1229864623987639D-01,
     4  0.2776749206212802D-01,0.5068963872921474D-01,
     5  0.7746000182236091D-01,0.1029582132198156D+00,
     6  0.1242523793131195D+00,0.1412392560046872D+00,
     7  0.1549941682784575D+00,0.1666743953511108D+00,
     8  0.1772790494779247D+00,0.1877236381539531D+00,
     9  0.1989842579945119D+00,0.2122867982533222D+00,
     *  0.2293955182511493D+00,0.2530519315237990D+00,
     1  0.2873773138759648D+00,0.3372621765491667D+00,
     2  0.4055227849029588D+00,0.4899499589425621D+00,
     3  0.5841350246819210D+00,0.6804436890171627D+00,
     4  0.7719110996571431D+00,0.8528214162995221D+00,
     5  0.9187334735244774D+00,0.9663768824285873D+00,
     *  0.9935603772854570D+00/

        data w3/
     1  0.1601823925056982D-05,0.4051401048812366D-04,
     2  0.3520648539622732D-03,0.1660850748229501D-02,
     3  0.5137929281303820D-02,0.1148412066850859D-01,
     4  0.1949446787362840D-01,0.2570904285027514D-01,
     5  0.2689876084878507D-01,0.2360406637248656D-01,
     6  0.1900336004473366D-01,0.1516705785902096D-01,
     7  0.1253808208066740D-01,0.1098717242596848D-01,
     8  0.1037230981440515D-01,0.1067654484583163D-01,
     9  0.1204438634252181D-01,0.1484724193108817D-01,
     *  0.1981827166571916D-01,0.2818663703907780D-01,
     1  0.4133112146920674D-01,0.5897280184709356D-01,
     2  0.7715842696540429D-01,0.9056858134292358D-01,
     3  0.9649872333355143D-01,0.9495030467195391D-01,
     4  0.8703235205367219D-01,0.7405007061103622D-01,
     5  0.5723063106844580D-01,0.3769127496295509D-01,
     *  0.1649122829313044D-01/

        data r4/
     1  0.2869971836289289D-06,0.1137064862159775D-04,
     2  0.1341929702769913D-03,0.8203737081126170D-03,
     3  0.3237746710276030D-02,0.9269181321390676D-02,
     4  0.2066457311904016D-01,0.3770158947165001D-01,
     5  0.5868391217952856D-01,0.8103704809800370D-01,
     6  0.1028335612751435D+00,0.1233169858690987D+00,
     7  0.1425900105753523D+00,0.1612475039895278D+00,
     8  0.1802340178256186D+00,0.2008627999434703D+00,
     9  0.2249197343537080D+00,0.2548018830978881D+00,
     *  0.2935575322911517D+00,0.3444837128769622D+00,
     1  0.4098896481625024D+00,0.4894051612103581D+00,
     2  0.5792717225450137D+00,0.6733345743704928D+00,
     3  0.7647214630794006D+00,0.8470697275781447D+00,
     4  0.9150902500933930D+00,0.9647364649454304D+00,
     *  0.9932320853632617D+00/

        data w4/
     1  0.1359485361254747D-05,0.3326780581481740D-04,
     2  0.2802665204329573D-03,0.1285181649959833D-02,
     3  0.3878178817325343D-02,0.8502557164858198D-02,
     4  0.1434562989358731D-01,0.1942363804371834D-01,
     5  0.2208195776845361D-01,0.2229557105168899D-01,
     6  0.2117893630127890D-01,0.1981160356292317D-01,
     7  0.1883567140360211D-01,0.1863624496658503D-01,
     8  0.1955385225357743D-01,0.2199880114785241D-01,
     9  0.2651308377979654D-01,0.3376464734464267D-01,
     *  0.4431727866999839D-01,0.5795260995932898D-01,
     1  0.7281176864806734D-01,0.8557276760171703D-01,
     2  0.9310345577436753D-01,0.9385626052597844D-01,
     3  0.8784851810663270D-01,0.7596428242578038D-01,
     4  0.5940221354064567D-01,0.3942691855259902D-01,
     *  0.1732347723342559D-01/

        data r5/
     1  0.2351806022319534D-06,0.9371601777833066D-05,
     2  0.1107082643478053D-03,0.6739376379504583D-03,
     3  0.2638272003890354D-02,0.7493792665864039D-02,
     4  0.1669269151817645D-01,0.3086571584100454D-01,
     5  0.4957931792755601D-01,0.7182202370163699D-01,
     6  0.9670746318326264D-01,0.1239292354325439D+00,
     7  0.1539341188863238D+00,0.1879585179894866D+00,
     8  0.2279827897596799D+00,0.2765070517192719D+00,
     9  0.3359654710899962D+00,0.4077134559202752D+00,
     *  0.4909503025536599D+00,0.5823044621295881D+00,
     1  0.6764393917382556D+00,0.7672247916500066D+00,
     2  0.8487745480647058D+00,0.9160538269103440D+00,
     3  0.9651389368767911D+00,0.9933094290125829D+00/

        data w5/
     1  0.1115893603802233D-05,0.2745191953993641D-04,
     2  0.2310340575439405D-03,0.1051314436226808D-02,
     3  0.3136526082030550D-02,0.6831872944864686D-02,
     4  0.1167562642786935D-01,0.1658838250943685D-01,
     5  0.2065726418543863D-01,0.2367293021142841D-01,
     6  0.2604445372243788D-01,0.2847308790498027D-01,
     7  0.3174325549431288D-01,0.3663883042479569D-01,
     8  0.4383584078685604D-01,0.5363609068576796D-01,
     9  0.6552672898282823D-01,0.7784382516131720D-01,
     *  0.8805641328317263D-01,0.9374097638099370D-01,
     1  0.9348922449275047D-01,0.8709078455015987D-01,
     2  0.7516572508113051D-01,0.5873770036661820D-01,
     3  0.3897800722006565D-01,0.1712553679382987D-01/

        data r6/
     1  0.2443400632125834D-06,0.9361435985030183D-05,
     2  0.1068439344896282D-03,0.6333547818598942D-03,
     3  0.2444238223002779D-02,0.6963392465380609D-02,
     4  0.1589038666324393D-01,0.3077774780309904D-01,
     5  0.5283390722889294D-01,0.8308464568179150D-01,
     6  0.1227177826327271D+00,0.1732693817638671D+00,
     7  0.2363272868221974D+00,0.3126761998360653D+00,
     8  0.4012856978117024D+00,0.4988527567895812D+00,
     9  0.6002519801650724D+00,0.6995362709306352D+00,
     *  0.7908888522293347D+00,0.8692441973537732D+00,
     1  0.9306588629578503D+00,0.9726883847505298D+00,
     *  0.9950033761158701D+00/

        data w6/
     1  0.1149293312796231D-05,0.2707365365272704D-04,
     2  0.2196415032785897D-03,0.9745451681444365D-03,
     3  0.2890068260986777D-02,0.6439575267218910D-02,
     4  0.1167679010916158D-01,0.1829615353039985D-01,
     5  0.2598050244243536D-01,0.3471498256401027D-01,
     6  0.4481193311687481D-01,0.5657114774826538D-01,
     7  0.6969435844917330D-01,0.8284039726907010D-01,
     8  0.9381946450237123D-01,0.1004420396713743D+00,
     9  0.1013487073152842D+00,0.9624036014336202D-01,
     *  0.8561500092371670D-01,0.7043297893307652D-01,
     1  0.5197797410006782D-01,0.3201641355583844D-01,
     *  0.1296874247892385D-01/
c        
c         %%%%%%%%%%%%%%%%%%%%%

ccc        call prin2('alpha = *',alpha,1)
        dtwo = 2
        if (alpha .ge. 1/dtwo) then
        do 1000 i=1,nqs(6)
        roots(i) =r6(i)
        weights(i) = w6(i)
 1000 continue

        nquads = nqs(6)
        endif

        if (alpha .ge. dtwo**(-2)) then
        if (alpha .lt. dtwo**(-1)) then

        do 2000 i=1,nqs(5)
        roots(i) =r5(i)
        weights(i) = w5(i)
 2000 continue

        nquads = nqs(5)

        endif
        endif


        if (alpha .ge. dtwo**(-3)) then
        if (alpha .lt. dtwo**(-2)) then

        do 3000 i=1,nqs(4)
        roots(i) =r4(i)
        weights(i) = w4(i)
 3000 continue

        nquads = nqs(4)

        endif
        endif


        if (alpha .ge. dtwo**(-4)) then
        if (alpha .lt. dtwo**(-3)) then

        do 4000 i=1,nqs(3)
        roots(i) =r3(i)
        weights(i) = w3(i)
 4000 continue

        nquads = nqs(3)

        endif
        endif


        if (alpha .ge. dtwo**(-5)) then
        if (alpha .lt. dtwo**(-4)) then

        do 5000 i=1,nqs(2)
        roots(i) =r2(i)
        weights(i) = w2(i)
 5000 continue

        nquads = nqs(2)

        endif
        endif


        if (alpha .lt. dtwo**(-5)) then

        do 6000 i=1,nqs(1)
        roots(i) =r1(i)
        weights(i) = w1(i)
 6000 continue

        nquads = nqs(1)

        endif

        return
        end
c
c
c
c
c
        subroutine lap_load_quads19(alpha,roots,weights,nquads)
        implicit real *8 (a-h,o-z)
        dimension r1(33),r2(32),r3(31),r4(29),r5(25),r6(22),
     1          w1(33),w2(32),w3(31),w4(29),w5(25),w6(22),nqs(6),
     2          roots(1),weights(1)

        data nqs/33,32,31,29,25,22/


        data r1/
     1  0.5777498096872038D-06,0.2288474056400583D-04,
     2  0.2702885712592717D-03,0.1657948198035502D-02,
     3  0.6599237565165783D-02,0.1920829806350446D-01,
     4  0.4393828489101669D-01,0.8246913261355575D-01,
     5  0.1293543807970466D+00,0.1702269091164211D+00,
     6  0.1953704500689091D+00,0.2091822858697152D+00,
     7  0.2174879414754571D+00,0.2231533309129341D+00,
     8  0.2274769397068781D+00,0.2311343712133850D+00,
     9  0.2345636736334689D+00,0.2381376340833007D+00,
     *  0.2422808900837322D+00,0.2476290565664344D+00,
     1  0.2553820883573565D+00,0.2682532714704561D+00,
     2  0.2928028229366426D+00,0.3405066041357018D+00,
     3  0.4147791116978645D+00,0.5029100346655221D+00,
     4  0.5921630158444857D+00,0.6796973995527945D+00,
     5  0.7662541714542700D+00,0.8468350534197480D+00,
     6  0.9146092091576833D+00,0.9644665822530418D+00,
     *  0.9931745205987632D+00/

        data w1/
     1  0.2736489878031350D-05,0.6695745641808172D-04,
     2  0.5649719469254307D-03,0.2606329086433771D-02,
     3  0.7984312815384039D-02,0.1802936440283729D-01,
     4  0.3178535399811143D-01,0.4444361444035084D-01,
     5  0.4668283011465065D-01,0.3315270665927609D-01,
     6  0.1824416503120287D-01,0.1037295834600213D-01,
     7  0.6684065327769558D-02,0.4848767465018319D-02,
     8  0.3905423193150608D-02,0.3479893429552158D-02,
     9  0.3439085261417356D-02,0.3777253979094453D-02,
     *  0.4610101784397522D-02,0.6272648926371216D-02,
     1  0.9647764977757842D-02,0.1712487611086181D-01,
     2  0.3412429211694978D-01,0.6212087125858995D-01,
     3  0.8388356940791066D-01,0.9009407918636028D-01,
     4  0.8806048998543070D-01,0.8736879471659959D-01,
     5  0.8477002386343386D-01,0.7520293922146814D-01,
     6  0.5949859745914455D-01,0.3968248018261188D-01,
     *  0.1746768135863911D-01/

        data r2/
     1  0.3873767093867011D-06,0.1579160449391647D-04,
     2  0.1920578234791119D-03,0.1212441649928748D-02,
     3  0.4955463665604445D-02,0.1474329147805521D-01,
     4  0.3419942503427241D-01,0.6423806851895797D-01,
     5  0.9850198235841527D-01,0.1255301063528927D+00,
     6  0.1522497549922609D+00,0.1769836208656419D+00,
     7  0.1948297133155925D+00,0.2074753706339834D+00,
     8  0.2170236463386961D+00,0.2248707155121976D+00,
     9  0.2319627498718751D+00,0.2390741569047036D+00,
     *  0.2470192493681261D+00,0.2568951298189119D+00,
     1  0.2705130187383804D+00,0.2912546331221134D+00,
     2  0.3252502098367053D+00,0.3802536960010034D+00,
     3  0.4585469452481370D+00,0.5531685163324032D+00,
     4  0.6538269800240759D+00,0.7514472018345842D+00,
     5  0.8389032866186268D+00,0.9107549302885457D+00,
     6  0.9629921236553322D+00,0.9929031711457575D+00/

        data w2/
     1  0.1845792692967511D-05,0.4666397799955374D-04,
     2  0.4070714094615389D-03,0.1939860415423601D-02,
     3  0.6117668799537254D-02,0.1411299510670174D-01,
     4  0.2503021932980798D-01,0.3406958971945714D-01,
     5  0.3172715223689434D-01,0.2479320706052417D-01,
     6  0.2749140618443985D-01,0.2124287837186846D-01,
     7  0.1483848751428991D-01,0.1080971817057836D-01,
     8  0.8513435805477226D-02,0.7333818267296676D-02,
     9  0.6976035568551458D-02,0.7378795901477582D-02,
     *  0.8687660362800196D-02,0.1135224613601805D-01,
     1  0.1642304487457717D-01,0.2608171125653275D-01,
     2  0.4334769959061750D-01,0.6713540766328462D-01,
     3  0.8813991530181894D-01,0.9931640367768606D-01,
     4  0.1004839513002412D+00,0.9357476316057453D-01,
     5  0.8043499372320411D-01,0.6260213821577112D-01,
     6  0.4142033754542796D-01,0.1816887755896598D-01/

        data r3/
     1  0.3287247650618073D-06,0.1189834898827308D-04,
     2  0.1257825909549172D-03,0.7040917014561371D-03,
     3  0.2745975111022903D-02,0.8368362132199745D-02,
     4  0.2058249374082620D-01,0.4160691633532044D-01,
     5  0.7055273125661464D-01,0.1028685431873764D+00,
     6  0.1331553452341919D+00,0.1585956603063203D+00,
     7  0.1791898921446799D+00,0.1961046350353544D+00,
     8  0.2106206673782640D+00,0.2238897469806542D+00,
     9  0.2370051712897761D+00,0.2511519685850517D+00,
     *  0.2678072950629314D+00,0.2890301877463561D+00,
     1  0.3178381665645801D+00,0.3583618799687185D+00,
     2  0.4147926278699505D+00,0.4884976859891866D+00,
     3  0.5758390037231845D+00,0.6695502943806675D+00,
     4  0.7616103051370999D+00,0.8449385997849748D+00,
     5  0.9138836353962480D+00,0.9642322439495634D+00,
     *  0.9931351787267855D+00/

        data w3/
     1  0.1531698227654661D-05,0.3358227951560693D-04,
     2  0.2489101054125620D-03,0.1063784455531823D-02,
     3  0.3373909282534618D-02,0.8402888270071051D-02,
     4  0.1644014248946964D-01,0.2547599344322218D-01,
     5  0.3160242272805179D-01,0.3205821488523040D-01,
     6  0.2804878932729951D-01,0.2286665275821813D-01,
     7  0.1853299967523339D-01,0.1551257741540916D-01,
     8  0.1370980161830690D-01,0.1300805114864566D-01,
     9  0.1341708476004827D-01,0.1511795532227976D-01,
     *  0.1852909475113621D-01,0.2441455823986569D-01,
     1  0.3390393712005471D-01,0.4789704387799045D-01,
     2  0.6523996011729819D-01,0.8153200131106366D-01,
     3  0.9188416836854708D-01,0.9417462028395866D-01,
     4  0.8875922273051629D-01,0.7695610976884027D-01,
     5  0.6023329632877842D-01,0.3998921482686120D-01,
     *  0.1757148061238108D-01/

        data r4/
     1  0.3346929676059012D-06,0.1334813938449397D-04,
     2  0.1584571573312020D-03,0.9728811689676323D-03,
     3  0.3848941761790749D-02,0.1103432423724806D-01,
     4  0.2465381935918670D-01,0.4518669369574161D-01,
     5  0.7082177893928679D-01,0.9851813899063944D-01,
     6  0.1257252794524200D+00,0.1512215241440908D+00,
     7  0.1749417770932063D+00,0.1975310418594082D+00,
     8  0.2200523141400822D+00,0.2438946670536177D+00,
     9  0.2708179277118596D+00,0.3030425336286944D+00,
     *  0.3432348548202127D+00,0.3941249921359251D+00,
     1  0.4574969786489961D+00,0.5328363180691191D+00,
     2  0.6167171478941819D+00,0.7035999520687778D+00,
     3  0.7873178740744279D+00,0.8622245843895498D+00,
     4  0.9237228697464194D+00,0.9683893774116790D+00,
     *  0.9939407672393777D+00/

        data w4/
     1  0.1587676995801523D-05,0.3913844860747862D-04,
     2  0.3317680109464139D-03,0.1527329902923168D-02,
     3  0.4616655613551478D-02,0.1013778311402440D-01,
     4  0.1719461588823031D-01,0.2354833532285968D-01,
     5  0.2718152360193348D-01,0.2777236112144208D-01,
     6  0.2644353467552388D-01,0.2454647340280209D-01,
     7  0.2300637340228412D-01,0.2235068831401490D-01,
     8  0.2292241916404158D-01,0.2505376969455530D-01,
     9  0.2916243432837561D-01,0.3573667250311370D-01,
     *  0.4512021278658856D-01,0.5698316372075962D-01,
     1  0.6968753877434371D-01,0.8040439583938462D-01,
     2  0.8640628904098361D-01,0.8631195653364296D-01,
     3  0.8017474827808601D-01,0.6887099176238526D-01,
     4  0.5355595119373512D-01,0.3539713203158041D-01,
     *  0.1551415585228466D-01/

        data r5/
     1  0.3880492306975629D-06,0.1473348723516151D-04,
     2  0.1657912874731596D-03,0.9620802781151550D-03,
     3  0.3602319775559434D-02,0.9853784124630416D-02,
     4  0.2134891244847540D-01,0.3881869200421601D-01,
     5  0.6188661987432842D-01,0.8950362528043401D-01,
     6  0.1206446942080724D+00,0.1548691914353943D+00,
     7  0.1926007424215349D+00,0.2351762513882884D+00,
     8  0.2846865071282982D+00,0.3434996928482982D+00,
     9  0.4133521690576568D+00,0.4941933376671271D+00,
     *  0.5834313690948669D+00,0.6761828395513481D+00,
     1  0.7663626228851638D+00,0.8478789293742985D+00,
     2  0.9154248720183345D+00,0.9648424610820970D+00,
     *  0.9932489519644638D+00/

        data w5/
     1  0.1822426957289209D-05,0.4244426964623663D-04,
     2  0.3379332756651379D-03,0.1455851090362873D-02,
     3  0.4132461734376675D-02,0.8657142705075239D-02,
     4  0.1445847583748723D-01,0.2040997012093030D-01,
     5  0.2553835991540182D-01,0.2951517429345091D-01,
     6  0.3268909102876364D-01,0.3582979046536404D-01,
     7  0.3985847665294078D-01,0.4564853040161904D-01,
     8  0.5378215199669317D-01,0.6415463001055076D-01,
     9  0.7555632557390468D-01,0.8568937033711310D-01,
     *  0.9195242591866360D-01,0.9251512218996632D-01,
     1  0.8681606328789655D-01,0.7532653985814848D-01,
     2  0.5907069861033862D-01,0.3928292062036103D-01,
     *  0.1727822737832249D-01/

        data r6/
     1  0.2744158607502753D-06,0.1066625871275257D-04,
     2  0.1238540957222321D-03,0.7467057561932747D-03,
     3  0.2922822340284041D-02,0.8413264148232129D-02,
     4  0.1932326440893504D-01,0.3754262443367131D-01,
     5  0.6446084090407024D-01,0.1011161883561661D+00,
     6  0.1485391395329319D+00,0.2078775205197720D+00,
     7  0.2799947651469308D+00,0.3646078649564360D+00,
     8  0.4595087342633007D+00,0.5604856350861258D+00,
     9  0.6620091134139910D+00,0.7581907279588116D+00,
     *  0.8435401086693549D+00,0.9133917763991992D+00,
     1  0.9640980903222380D+00,0.9931162062803247D+00/

        data w6/
     1  0.1294019801796107D-05,0.3102200237153032D-04,
     2  0.2567854350612253D-03,0.1160784568226616D-02,
     3  0.3491988303491068D-02,0.7850448090980897D-02,
     4  0.1428966510676444D-01,0.2237746582084955D-01,
     5  0.3162520105250176D-01,0.4185430417728283D-01,
     6  0.5318967375842669D-01,0.6564683154705083D-01,
     7  0.7855528110966096D-01,0.9030436206309475D-01,
     8  0.9878249602631297D-01,0.1022324840881029D+00,
     9  0.9982023208915139D-01,0.9162114236573427D-01,
     *  0.7829713810219845D-01,0.6079705832228091D-01,
     1  0.4019036779790175D-01,0.1762397415275239D-01/
c        
c         %%%%%%%%%%%%%%%%%%%%%

ccc        call prin2('alpha = *',alpha,1)
        dtwo = 2
        if (alpha .ge. 1/dtwo) then
        do 1000 i=1,nqs(6)
        roots(i) =r6(i)
        weights(i) = w6(i)
 1000 continue

        nquads = nqs(6)
        endif

        if (alpha .ge. dtwo**(-2)) then
        if (alpha .lt. dtwo**(-1)) then

        do 2000 i=1,nqs(5)
        roots(i) =r5(i)
        weights(i) = w5(i)
 2000 continue

        nquads = nqs(5)

        endif
        endif


        if (alpha .ge. dtwo**(-3)) then
        if (alpha .lt. dtwo**(-2)) then

        do 3000 i=1,nqs(4)
        roots(i) =r4(i)
        weights(i) = w4(i)
 3000 continue

        nquads = nqs(4)

        endif
        endif


        if (alpha .ge. dtwo**(-4)) then
        if (alpha .lt. dtwo**(-3)) then

        do 4000 i=1,nqs(3)
        roots(i) =r3(i)
        weights(i) = w3(i)
 4000 continue

        nquads = nqs(3)

        endif
        endif


        if (alpha .ge. dtwo**(-5)) then
        if (alpha .lt. dtwo**(-4)) then

        do 5000 i=1,nqs(2)
        roots(i) =r2(i)
        weights(i) = w2(i)
 5000 continue

        nquads = nqs(2)

        endif
        endif


        if (alpha .lt. dtwo**(-5)) then

        do 6000 i=1,nqs(1)
        roots(i) =r1(i)
        weights(i) = w1(i)
 6000 continue

        nquads = nqs(1)

        endif

        return
        end
c
c
c
c
c
        subroutine lap_load_quads20(alpha,roots,weights,nquads)
        implicit real *8 (a-h,o-z)
        dimension r1(33),r2(32),r3(30),r4(28),r5(25),r6(22),
     1          w1(33),w2(32),w3(30),w4(28),w5(25),w6(22),nqs(6),
     2          roots(1),weights(1)

        data nqs/33,32,30,28,25,22/


        data r1/
     1  0.5043469240042000D-06,0.2022529802564272D-04,
     2  0.2413860152807270D-03,0.1492856971125932D-02,
     3  0.5976865471415200D-02,0.1746807167426997D-01,
     4  0.4013688689281063D-01,0.7614687870643071D-01,
     5  0.1232570396327131D+00,0.1736006367127779D+00,
     6  0.2153222182808499D+00,0.2420114060488575D+00,
     7  0.2575732862404726D+00,0.2672790993174669D+00,
     8  0.2740246460502503D+00,0.2792234713484424D+00,
     9  0.2836439931372219D+00,0.2877967641805648D+00,
     *  0.2921202516619875D+00,0.2971114571086958D+00,
     1  0.3035004582148211D+00,0.3126178094814233D+00,
     2  0.3272958819598049D+00,0.3537540955185076D+00,
     3  0.4019171939369476D+00,0.4758234285138982D+00,
     4  0.5671837458118445D+00,0.6647057629168830D+00,
     5  0.7592808324444416D+00,0.8439865803377077D+00,
     6  0.9135709050107688D+00,0.9641592748616685D+00,
     *  0.9931268955921536D+00/

        data w1/
     1  0.2395511148486089D-05,0.5940752262724848D-04,
     2  0.5068235579227714D-03,0.2357240433515388D-02,
     3  0.7259803188206895D-02,0.1646123185828198D-01,
     4  0.2928024533439229D-01,0.4238692623955571D-01,
     5  0.5051500260878420D-01,0.4801785810261307D-01,
     6  0.3420753349081922D-01,0.2007095581198170D-01,
     7  0.1196585663189468D-01,0.7904600001372589D-02,
     8  0.5808677077773558D-02,0.4711243326465429D-02,
     9  0.4212103808894590D-02,0.4164541677379092D-02,
     *  0.4562658345561503D-02,0.5536035844069973D-02,
     1  0.7450343746880757D-02,0.1122375331819504D-01,
     2  0.1913085857501104D-01,0.3562292106049405D-01,
     3  0.6156547226847629D-01,0.8460392694033493D-01,
     4  0.9617177077752788D-01,0.9736379281996965D-01,
     5  0.9064161007139449D-01,0.7789914678260801D-01,
     6  0.6062579430222006D-01,0.4011342126567930D-01,
     *  0.1759604769794813D-01/

        data r2/
     1  0.4974563341414117D-06,0.2001886832081410D-04,
     2  0.2390958881844835D-03,0.1475430020556342D-02,
     3  0.5873602806243705D-02,0.1699211036710178D-01,
     4  0.3841458895971831D-01,0.7119006198739146D-01,
     5  0.1120261621333440D+00,0.1544153028483606D+00,
     6  0.1917426508019797D+00,0.2202846972948088D+00,
     7  0.2406439060745370D+00,0.2553757664120871D+00,
     8  0.2666934038671084D+00,0.2760950715425225D+00,
     9  0.2846364905759139D+00,0.2932018706869713D+00,
     *  0.3027200710896535D+00,0.3144104220252504D+00,
     1  0.3301867655628224D+00,0.3533540162348442D+00,
     2  0.3892920701826527D+00,0.4440089759474484D+00,
     3  0.5186788660901682D+00,0.6069491027756509D+00,
     4  0.6993681046402296D+00,0.7874425576503459D+00,
     5  0.8646682050879784D+00,0.9264705553623897D+00,
     6  0.9700724233196510D+00,0.9943371736274485D+00/

        data w2/
     1  0.2365234254464787D-05,0.5884340263639214D-04,
     2  0.5018512583431694D-03,0.2324085398877767D-02,
     3  0.7089018654724624D-02,0.1579288565200053D-01,
     4  0.2726409995764868D-01,0.3771711430582221D-01,
     5  0.4280908202383678D-01,0.4082109797734989D-01,
     6  0.3317181910928629D-01,0.2407223503388419D-01,
     7  0.1711213952467921D-01,0.1272207390982313D-01,
     8  0.1015722966734332D-01,0.8817575393063460D-02,
     9  0.8409244717128278D-02,0.8872527820024934D-02,
     *  0.1036077583151983D-01,0.1332720938445591D-01,
     1  0.1876038579482407D-01,0.2848367078987036D-01,
     2  0.4450604012809742D-01,0.6518903419887180D-01,
     3  0.8298443620084922D-01,0.9191161637869169D-01,
     4  0.9150024724842718D-01,0.8356755861294907D-01,
     5  0.7012348451741378D-01,0.5301994579897093D-01,
     6  0.3400223891856094D-01,0.1454806715577048D-01/

        data r3/
     1  0.5639559034922784D-06,0.2224783107486844D-04,
     2  0.2608986897755215D-03,0.1581894138439489D-02,
     3  0.6187949959195583D-02,0.1758167398072902D-01,
     4  0.3900553974620439D-01,0.7088111433460775D-01,
     5  0.1093172826553499D+00,0.1478245977450751D+00,
     6  0.1813826353842322D+00,0.2085415376887265D+00,
     7  0.2302069501594223D+00,0.2479975316819424D+00,
     8  0.2635923919837708D+00,0.2784803746921423D+00,
     9  0.2938933664478119D+00,0.3111037089981701D+00,
     *  0.3318062222983273D+00,0.3584444851126535D+00,
     1  0.3944450198376256D+00,0.4438303953478600D+00,
     2  0.5092510098606737D+00,0.5891561924331386D+00,
     3  0.6774261278816685D+00,0.7659706188798603D+00,
     4  0.8471817513400710D+00,0.9149117882802004D+00,
     5  0.9646034741494559D+00,0.9932009060177432D+00/

        data w3/
     1  0.2669485714694070D-05,0.6497447278213920D-04,
     2  0.5429624834054850D-03,0.2464158957065905D-02,
     3  0.7358862479120722D-02,0.1601717914922827D-01,
     4  0.2693211015754213D-01,0.3614115803121110D-01,
     5  0.3956727142215446D-01,0.3657162318649692D-01,
     6  0.3033400513948747D-01,0.2416493139460210D-01,
     7  0.1944383267939462D-01,0.1642094527201387D-01,
     8  0.1501922567557905D-01,0.1495572585947043D-01,
     9  0.1607579434748954D-01,0.1862210444959608D-01,
     *  0.2318515910292760D-01,0.3066573701674547D-01,
     1  0.4203876694928999D-01,0.5723253076505861D-01,
     2  0.7335815595069952D-01,0.8537761189659130D-01,
     3  0.8977070014310679D-01,0.8603643854882538D-01,
     4  0.7536051376017275D-01,0.5934202795568925D-01,
     5  0.3953281146395924D-01,0.1740001180457914D-01/

        data r4/
     1  0.4533596253276711D-06,0.1766165535407989D-04,
     2  0.2051672902985155D-03,0.1233914296109404D-02,
     3  0.4786904912133911D-02,0.1348390624690898D-01,
     4  0.2971391264446005D-01,0.5402060426646134D-01,
     5  0.8451054990309612D-01,0.1179080122026200D+00,
     6  0.1512557767904907D+00,0.1829102628847130D+00,
     7  0.2125563364586749D+00,0.2408080926261436D+00,
     8  0.2688715491356652D+00,0.2983895485861691D+00,
     9  0.3314119172375188D+00,0.3703645077337488D+00,
     *  0.4178222810288905D+00,0.4758557456036951D+00,
     1  0.5449435764354309D+00,0.6230432690590438D+00,
     2  0.7056610704682070D+00,0.7869447901652918D+00,
     3  0.8609709826485927D+00,0.9225900849839667D+00,
     4  0.9677919163880345D+00,0.9938124435404889D+00/

        data w4/
     1  0.2139331635268970D-05,0.5139119525034050D-04,
     2  0.4252362593515086D-03,0.1912197738042885D-02,
     3  0.5651527755013801D-02,0.1217193729571003D-01,
     4  0.2039047793746439D-01,0.2788666102329133D-01,
     5  0.3251494721479038D-01,0.3377009874942354D-01,
     6  0.3265310674183235D-01,0.3061117814784411D-01,
     7  0.2878911169340514D-01,0.2791781740580022D-01,
     8  0.2848310713222798D-01,0.3089458564590772D-01,
     9  0.3555759956304905D-01,0.4278738239278004D-01,
     *  0.5249287207653327D-01,0.6367197723107812D-01,
     1  0.7416217366550323D-01,0.8126743502313083D-01,
     2  0.8296970570767616D-01,0.7859653761634523D-01,
     3  0.6859212787778883D-01,0.5397643010929802D-01,
     4  0.3596585354128050D-01,0.1583438392854573D-01/

        data r5/
     1  0.3054140140626802D-06,0.1220683347503413D-04,
     2  0.1446161578156958D-03,0.8835270342669381D-03,
     3  0.3477911478320047D-02,0.9964581107329006D-02,
     4  0.2247086660793121D-01,0.4218044110918295D-01,
     5  0.6882902339780016D-01,0.1011151746593468D+00,
     6  0.1375809300191166D+00,0.1773343530030530D+00,
     7  0.2203839728179283D+00,0.2676567929809623D+00,
     8  0.3207720062361733D+00,0.3815438969646411D+00,
     9  0.4511915019713791D+00,0.5294416271924578D+00,
     *  0.6139850714189429D+00,0.7006715681284557D+00,
     1  0.7843218831735529D+00,0.8596602658968753D+00,
     2  0.9219929499549210D+00,0.9675743574954894D+00,
     *  0.9937734880121647D+00/

        data w5/
     1  0.1450136312610632D-05,0.3579153714271754D-04,
     2  0.3022383034241838D-03,0.1382101972678410D-02,
     3  0.4158798165572319D-02,0.9188914102203148D-02,
     4  0.1603043851937725D-01,0.2333637598023061D-01,
     5  0.2972742719077615D-01,0.3458939265943572D-01,
     6  0.3818820815905080D-01,0.4131996180473875D-01,
     7  0.4493744215028966D-01,0.4988619775355380D-01,
     8  0.5666016582781558D-01,0.6510299506924722D-01,
     9  0.7415691335165276D-01,0.8195485432144139D-01,
     *  0.8642533313857791D-01,0.8606795471742072D-01,
     1  0.8034203018824652D-01,0.6954666399610485D-01,
     2  0.5449202193285534D-01,0.3623064176942230D-01,
     *  0.1593568725242926D-01/

        data r6/
     1  0.2497083389218221D-06,0.1025285307845646D-04,
     2  0.1240432432216531D-03,0.7713906329357771D-03,
     3  0.3091311592661669D-02,0.9058257592041935D-02,
     4  0.2108176782243518D-01,0.4133818374859219D-01,
     5  0.7135705453527452D-01,0.1120804440873640D+00,
     6  0.1641684306378662D+00,0.2281328808863286D+00,
     7  0.3040200016669326D+00,0.3907561321515361D+00,
     8  0.4856520487743542D+00,0.5845046077320805D+00,
     9  0.6822380323256484D+00,0.7736663029591145D+00,
     *  0.8540529736729893D+00,0.9194139725008649D+00,
     1  0.9666509471750545D+00,0.9936115935239254D+00/

        data w6/
     1  0.1193352606000712D-05,0.3030342176077959D-04,
     2  0.2615838179347206D-03,0.1219909710762728D-02,
     3  0.3756461973819812D-02,0.8591700523564156D-02,
     4  0.1582549409348640D-01,0.2494060312661562D-01,
     5  0.3524792804918054D-01,0.4630585429916175D-01,
     6  0.5796097251521550D-01,0.6998792135501996D-01,
     7  0.8161692575580054D-01,0.9140521402062715D-01,
     8  0.9767416572802719D-01,0.9917145217981025D-01,
     9  0.9542168072284907D-01,0.8664067435397882D-01,
     *  0.7346666162852140D-01,0.5673813582714233D-01,
     1  0.3737595929914940D-01,0.1635920424496587D-01/
c        
c         %%%%%%%%%%%%%%%%%%%%%

ccc        call prin2('alpha = *',alpha,1)
        dtwo = 2
        if (alpha .ge. 1/dtwo) then
        do 1000 i=1,nqs(6)
        roots(i) =r6(i)
        weights(i) = w6(i)
 1000 continue

        nquads = nqs(6)
        endif

        if (alpha .ge. dtwo**(-2)) then
        if (alpha .lt. dtwo**(-1)) then

        do 2000 i=1,nqs(5)
        roots(i) =r5(i)
        weights(i) = w5(i)
 2000 continue

        nquads = nqs(5)

        endif
        endif


        if (alpha .ge. dtwo**(-3)) then
        if (alpha .lt. dtwo**(-2)) then

        do 3000 i=1,nqs(4)
        roots(i) =r4(i)
        weights(i) = w4(i)
 3000 continue

        nquads = nqs(4)

        endif
        endif


        if (alpha .ge. dtwo**(-4)) then
        if (alpha .lt. dtwo**(-3)) then

        do 4000 i=1,nqs(3)
        roots(i) =r3(i)
        weights(i) = w3(i)
 4000 continue

        nquads = nqs(3)

        endif
        endif


        if (alpha .ge. dtwo**(-5)) then
        if (alpha .lt. dtwo**(-4)) then

        do 5000 i=1,nqs(2)
        roots(i) =r2(i)
        weights(i) = w2(i)
 5000 continue

        nquads = nqs(2)

        endif
        endif


        if (alpha .lt. dtwo**(-5)) then

        do 6000 i=1,nqs(1)
        roots(i) =r1(i)
        weights(i) = w1(i)
 6000 continue

        nquads = nqs(1)

        endif

        return
        end
c
c
c
c
c
        subroutine lap_load_quads21(alpha,roots,weights,nquads)
        implicit real *8 (a-h,o-z)
        dimension r1(33),r2(32),r3(30),r4(28),r5(24),r6(21),
     1          w1(33),w2(32),w3(30),w4(28),w5(24),w6(21),
     1          nqs(6),
     2          roots(1),weights(1)

        data nqs/33,32,30,28,24,21/


        data r1/
     1  0.5930034788191595D-06,0.2371194347328105D-04,
     2  0.2824508493581126D-03,0.1744738498692855D-02,
     3  0.6982933665799360D-02,0.2043088799679578D-01,
     4  0.4710955612889593D-01,0.8996533182691136D-01,
     5  0.1468572336014081D+00,0.2080884178352926D+00,
     6  0.2583358911257991D+00,0.2900891448532474D+00,
     7  0.3085716597827789D+00,0.3201294954005386D+00,
     8  0.3281821569933929D+00,0.3343980350360575D+00,
     9  0.3396868383105411D+00,0.3446542980158413D+00,
     *  0.3498197186011733D+00,0.3557666150688066D+00,
     1  0.3633361858052754D+00,0.3740042338479538D+00,
     2  0.3906167144624769D+00,0.4171311876505639D+00,
     3  0.4481645742317984D+00,0.4931153676527576D+00,
     4  0.5704871371531783D+00,0.6633689558119333D+00,
     5  0.7570227698000755D+00,0.8421121598715006D+00,
     6  0.9124163230660221D+00,0.9636557875972951D+00,
     *  0.9930281734875361D+00/

        data w1/
     1  0.2814533166373770D-05,0.6959283965359171D-04,
     2  0.5926066761375095D-03,0.2753530577989591D-02,
     3  0.8483502726425517D-02,0.1929682266117413D-01,
     4  0.3460340116941303D-01,0.5079612147137572D-01,
     5  0.6137675729025317D-01,0.5824908704890267D-01,
     6  0.4085917788107600D-01,0.2382407756634959D-01,
     7  0.1422988465779077D-01,0.9426666695276542D-02,
     8  0.6940949242199596D-02,0.5635617927191408D-02,
     9  0.5039675761529020D-02,0.4979511124665697D-02,
     *  0.5445660883553423D-02,0.6582614884660566D-02,
     1  0.8790057672593394D-02,0.1300412138483243D-01,
     2  0.2102436282812116D-01,0.3129049267765825D-01,
     3  0.3137946924517803D-01,0.6295200825727951D-01,
     4  0.8816994854965607D-01,0.9518047882363423D-01,
     5  0.9061143812055031D-01,0.7854839396523686D-01,
     6  0.6135359738641388D-01,0.4065998427751820D-01,
     *  0.1784757319254375D-01/

        data r2/
     1  0.4947985409297258D-06,0.2006463868720604D-04,
     2  0.2414080213764556D-03,0.1500461452252477D-02,
     3  0.6018212276149268D-02,0.1756407326351609D-01,
     4  0.4017014509015407D-01,0.7562660187297233D-01,
     5  0.1213350873653182D+00,0.1705765994970462D+00,
     6  0.2158060807298128D+00,0.2522202684243293D+00,
     7  0.2791503443592958D+00,0.2987458137012279D+00,
     8  0.3135558816133551D+00,0.3254998477577873D+00,
     9  0.3359270046758860D+00,0.3458913431714587D+00,
     *  0.3563880795185574D+00,0.3685772213097362D+00,
     1  0.3840867554609869D+00,0.4054999280210817D+00,
     2  0.4369448589895058D+00,0.4837307382247493D+00,
     3  0.5489023736829899D+00,0.6289713511084253D+00,
     4  0.7154054983484628D+00,0.7992132858254817D+00,
     5  0.8731668505856812D+00,0.9321149891468435D+00,
     6  0.9729840661371166D+00,0.9950049194917669D+00/

        data w2/
     1  0.2356486547978569D-05,0.5912775999896325D-04,
     2  0.5084502319417312D-03,0.2374370929052781D-02,
     3  0.7311255674763628D-02,0.1649840577624084D-01,
     4  0.2905331490809556D-01,0.4142530442614598D-01,
     5  0.4880785506932029D-01,0.4835178938745144D-01,
     6  0.4128228857696887D-01,0.3146040229976343D-01,
     7  0.2280708491544491D-01,0.1682298162364157D-01,
     8  0.1311380921981325D-01,0.1099405332689506D-01,
     9  0.1003210711137134D-01,0.1005893926251118D-01,
     *  0.1112310342523858D-01,0.1351848516501510D-01,
     1  0.1791941597875397D-01,0.2559360540892120D-01,
     2  0.3824031139496889D-01,0.5593851791541994D-01,
     3  0.7378456479344345D-01,0.8483807173218866D-01,
     4  0.8649963239224692D-01,0.7991010392675012D-01,
     5  0.6713967442982169D-01,0.5023931157062430D-01,
     6  0.3136921900632003D-01,0.1292208587431833D-01/

        data r3/
     1  0.5125952637933203D-06,0.2054067312239327D-04,
     2  0.2443854051891628D-03,0.1502146444953944D-02,
     3  0.5955555794424073D-02,0.1716541956804334D-01,
     4  0.3872574473900515D-01,0.7185608390637876D-01,
     5  0.1136758043223347D+00,0.1579975133290766D+00,
     6  0.1986891522163996D+00,0.2327214613773546D+00,
     7  0.2602962390142118D+00,0.2832005571652813D+00,
     8  0.3032948615143786D+00,0.3219254519091276D+00,
     9  0.3402096715066642D+00,0.3594120516442717D+00,
     *  0.3811448209163303D+00,0.4075397447859065D+00,
     1  0.4414057314167587D+00,0.4860648074505761D+00,
     2  0.5441971105662744D+00,0.6156114586708203D+00,
     3  0.6959390078250844D+00,0.7780223118257994D+00,
     4  0.8543739116626790D+00,0.9186533350159455D+00,
     5  0.9660901229204681D+00,0.9934794197771529D+00/

        data w3/
     1  0.2434811983770311D-05,0.6029624916836424D-04,
     2  0.5120149454658550D-03,0.2360091518438421D-02,
     3  0.7163833873194724D-02,0.1589969758773040D-01,
     4  0.2745809243427677D-01,0.3831029566900497D-01,
     5  0.4423357470992021D-01,0.4332816034964867D-01,
     6  0.3755349653175738D-01,0.3057720429419980D-01,
     7  0.2489680068309196D-01,0.2123203606418312D-01,
     8  0.1917325603279994D-01,0.1826760399023716D-01,
     9  0.1850881173925228D-01,0.2016185418119653D-01,
     *  0.2365571341198405D-01,0.2960315133732231D-01,
     1  0.3870040120162712D-01,0.5110732332915900D-01,
     2  0.6515397379721813D-01,0.7691710938529293D-01,
     3  0.8249819152490459D-01,0.8039673525610229D-01,
     4  0.7124445992194106D-01,0.5652024758607715D-01,
     5  0.3782003801591924D-01,0.1668309956690181D-01/

        data r4/
     1  0.3463210910990835D-06,0.1408153712246907D-04,
     2  0.1698500545283119D-03,0.1056854373299038D-02,
     3  0.4234197925496467D-02,0.1231942490568680D-01,
     4  0.2808316084414944D-01,0.5288466293863767D-01,
     5  0.8564628012941022D-01,0.1232876067283316D+00,
     6  0.1622900434856405D+00,0.2001412273197074D+00,
     7  0.2358185349293957D+00,0.2695204853728257D+00,
     8  0.3022349813330292D+00,0.3354571026401564D+00,
     9  0.3710719888202864D+00,0.4112786771794225D+00,
     *  0.4583751122253226D+00,0.5142200633670160D+00,
     1  0.5793446001044130D+00,0.6521210324191227D+00,
     2  0.7286898651776831D+00,0.8038267780652760D+00,
     3  0.8721301167546901D+00,0.9288813192255965D+00,
     4  0.9704367375633492D+00,0.9943238298473650D+00/

        data w4/
     1  0.1650180850808305D-05,0.4153609563344654D-04,
     2  0.3580759044191125D-03,0.1672447220921128D-02,
     3  0.5134389067987302D-02,0.1152733771618174D-01,
     4  0.2024644456234287D-01,0.2916879082904653D-01,
     5  0.3581128810416590D-01,0.3885852687081762D-01,
     6  0.3872049659459017D-01,0.3681760896227987D-01,
     7  0.3457797960839852D-01,0.3299523916626758D-01,
     8  0.3268455993750469D-01,0.3407470502764656D-01,
     9  0.3752466632999963D-01,0.4328044627598060D-01,
     *  0.5123799537814653D-01,0.6055639612307569D-01,
     1  0.6942323540897700D-01,0.7546997472720305D-01,
     2  0.7677540052323380D-01,0.7258362782047255D-01,
     3  0.6323062469583324D-01,0.4966343364988536D-01,
     4  0.3303535004013739D-01,0.1452777317800130D-01/

        data r5/
     1  0.3686388509001472D-06,0.1480988697998827D-04,
     2  0.1761042199867546D-03,0.1077969921361915D-02,
     3  0.4243407765376845D-02,0.1213681349406332D-01,
     4  0.2728710031932366D-01,0.5103875604458078D-01,
     5  0.8300294393229916D-01,0.1216036700537700D+00,
     6  0.1651156544930132D+00,0.2124670028467239D+00,
     7  0.2635793254876822D+00,0.3193286024717429D+00,
     8  0.3811867176338477D+00,0.4505283165741299D+00,
     9  0.5276868134791959D+00,0.6111299716493767D+00,
     *  0.6972742156370314D+00,0.7810991414245904D+00,
     1  0.8571552688314326D+00,0.9204333971515023D+00,
     2  0.9668800984150858D+00,0.9936354371307942D+00/

        data w5/
     1  0.1752457312563815D-05,0.4349032773142620D-04,
     2  0.3685477698678081D-03,0.1687315555346077D-02,
     3  0.5070324999591180D-02,0.1116201779090115D-01,
     4  0.1937288214378286D-01,0.2805456615710920D-01,
     5  0.3559183137122648D-01,0.4130831227060480D-01,
     6  0.4553169691318218D-01,0.4915873884540545D-01,
     7  0.5321776012840967D-01,0.5853775531240217D-01,
     8  0.6542423406243687D-01,0.7332977977286884D-01,
     9  0.8073545103508291D-01,0.8554213524652287D-01,
     *  0.8589134031657690D-01,0.8083531381664523D-01,
     1  0.7043376177823526D-01,0.5544341183190217D-01,
     2  0.3697138295918237D-01,0.1628619713767357D-01/

        data r6/
     1  0.3853581433876694D-06,0.1504378465271227D-04,
     2  0.1745109375117942D-03,0.1048072741691043D-02,
     3  0.4083403067587969D-02,0.1170329459883503D-01,
     4  0.2677762308111198D-01,0.5182298811736253D-01,
     5  0.8852369922555939D-01,0.1377789359056462D+00,
     6  0.1999787265722530D+00,0.2750409526963107D+00,
     7  0.3619964084362822D+00,0.4584085542138315D+00,
     8  0.5601894398015496D+00,0.6620568895454270D+00,
     9  0.7583362516752517D+00,0.8436735683371256D+00,
     *  0.9134772292265513D+00,0.9641359690983017D+00,
     *  0.9931236806355162D+00/

        data w6/
     1  0.1820201864508535D-05,0.4377254197194644D-04,
     2  0.3613197739325609D-03,0.1624529291510589D-02,
     3  0.4859805734163352D-02,0.1087336281935753D-01,
     4  0.1970322201869741D-01,0.3066657051733492D-01,
     5  0.4287881927659870D-01,0.5569540966511551D-01,
     6  0.6869816793417891D-01,0.8127607066178786D-01,
     7  0.9223252348845576D-01,0.9990187808720369D-01,
     8  0.1027649519860817D+00,0.1000118923339130D+00,
     9  0.9164734578454736D-01,0.7825866797014187D-01,
     *  0.6074528238226302D-01,0.4014963204963891D-01,
     *  0.1760495548124084D-01/
c        
c         %%%%%%%%%%%%%%%%%%%%%

ccc        call prin2('alpha = *',alpha,1)
        dtwo = 2
        if (alpha .ge. 1/dtwo) then
        do 1000 i=1,nqs(6)
        roots(i) =r6(i)
        weights(i) = w6(i)
 1000 continue

        nquads = nqs(6)
        endif

        if (alpha .ge. dtwo**(-2)) then
        if (alpha .lt. dtwo**(-1)) then

        do 2000 i=1,nqs(5)
        roots(i) =r5(i)
        weights(i) = w5(i)
 2000 continue

        nquads = nqs(5)

        endif
        endif


        if (alpha .ge. dtwo**(-3)) then
        if (alpha .lt. dtwo**(-2)) then

        do 3000 i=1,nqs(4)
        roots(i) =r4(i)
        weights(i) = w4(i)
 3000 continue

        nquads = nqs(4)

        endif
        endif


        if (alpha .ge. dtwo**(-4)) then
        if (alpha .lt. dtwo**(-3)) then

        do 4000 i=1,nqs(3)
        roots(i) =r3(i)
        weights(i) = w3(i)
 4000 continue

        nquads = nqs(3)

        endif
        endif


        if (alpha .ge. dtwo**(-5)) then
        if (alpha .lt. dtwo**(-4)) then

        do 5000 i=1,nqs(2)
        roots(i) =r2(i)
        weights(i) = w2(i)
 5000 continue

        nquads = nqs(2)

        endif
        endif


        if (alpha .lt. dtwo**(-5)) then

        do 6000 i=1,nqs(1)
        roots(i) =r1(i)
        weights(i) = w1(i)
 6000 continue

        nquads = nqs(1)

        endif

        return
        end
c
c
c
c
c
        subroutine lap_load_quads22(alpha,roots,weights,nquads)
        implicit real *8 (a-h,o-z)
        dimension r1(32),r2(32),r3(30),r4(27),r5(24),r6(21),
     1          w1(32),w2(32),w3(30),w4(27),w5(24),w6(21),nqs(6),
     2          roots(1),weights(1)

        data nqs/32,32,30,27,24,21/


        data r1/
     1  0.5725413276283860D-06,0.2315611156711776D-04,
     2  0.2785071656241312D-03,0.1733496169725412D-02,
     3  0.6975616794198303D-02,0.2048416163108566D-01,
     4  0.4736991400597671D-01,0.9082112763913564D-01,
     5  0.1495066454095597D+00,0.2159110461802674D+00,
     6  0.2769927314280076D+00,0.3217157475342793D+00,
     7  0.3499591347774660D+00,0.3675975244767268D+00,
     8  0.3793927668568186D+00,0.3880312232780117D+00,
     9  0.3949867161896577D+00,0.4011716791636832D+00,
     *  0.4072839705334115D+00,0.4140157782553249D+00,
     1  0.4222591583596929D+00,0.4334517406061478D+00,
     2  0.4502935769058299D+00,0.4780116740558569D+00,
     3  0.5245364038359090D+00,0.5939525636254902D+00,
     4  0.6788476683395263D+00,0.7668137948469024D+00,
     5  0.8479384776382350D+00,0.9154731502150898D+00,
     6  0.9648828072704216D+00,0.9932596403360068D+00/

        data w1/
     1  0.2724414071488398D-05,0.6820629550269167D-04,
     2  0.5867503867773594D-03,0.2747031919966667D-02,
     3  0.8505411307548072D-02,0.1940978394104225D-01,
     4  0.3493876862846649D-01,0.5177036358265890D-01,
     5  0.6435123476991759D-01,0.6616142248528024D-01,
     6  0.5393055827533841D-01,0.3564790342211937D-01,
     7  0.2194768209592335D-01,0.1412220736076246D-01,
     8  0.9898509211947662D-02,0.7613495269474030D-02,
     9  0.6444219400556727D-02,0.6038828832617456D-02,
     *  0.6297998615672885D-02,0.7309064546166203D-02,
     1  0.9406191329147015D-02,0.1340636921433831D-01,
     2  0.2112761271883303D-01,0.3573693279520621D-01,
     3  0.5820719250225109D-01,0.7920273189475569D-01,
     4  0.8838738429580594D-01,0.8590032995232690D-01,
     5  0.7524255988223095D-01,0.5908027186728974D-01,
     6  0.3925748029909170D-01,0.1725277848691309D-01/

        data r2/
     1  0.4899953403011341D-06,0.1999722714440589D-04,
     2  0.2420697184165262D-03,0.1513358032262974D-02,
     3  0.6105161395919574D-02,0.1793337766862783D-01,
     4  0.4135871077560469D-01,0.7878409041826494D-01,
     5  0.1284576554822554D+00,0.1841097554027559D+00,
     6  0.2374019942179101D+00,0.2820738497316416D+00,
     7  0.3162547418972567D+00,0.3415015493442630D+00,
     8  0.3604703489878106D+00,0.3754512005444621D+00,
     9  0.3881129398184152D+00,0.3997130250735647D+00,
     *  0.4113437619061007D+00,0.4241487017422634D+00,
     1  0.4395634306666574D+00,0.4596571717215855D+00,
     2  0.4875690177760381D+00,0.5275575904657592D+00,
     3  0.5832239867252578D+00,0.6537123610544116D+00,
     4  0.7325117052215977D+00,0.8108560731614858D+00,
     5  0.8809651897972756D+00,0.9370766562267682D+00,
     6  0.9756013282923993D+00,0.9956516420598025D+00/

        data w2/
     1  0.2336819883288121D-05,0.5905404477826066D-04,
     2  0.5112845966166059D-03,0.2403447029909976D-02,
     3  0.7453288551727699D-02,0.1697387185812751D-01,
     4  0.3031804159630804D-01,0.4424956281950997D-01,
     5  0.5401030293505649D-01,0.5582332452239453D-01,
     6  0.4967577951775980D-01,0.3937801600412928D-01,
     7  0.2928878337694418D-01,0.2167194185074028D-01,
     8  0.1664916260660264D-01,0.1358629048702134D-01,
     9  0.1194260737964030D-01,0.1143616173847935D-01,
     *  0.1201209455161829D-01,0.1383146972793744D-01,
     1  0.1733481955201516D-01,0.2337115050043841D-01,
     2  0.3318804541876867D-01,0.4746795495365847D-01,
     3  0.6375097472666069D-01,0.7606455180344085D-01,
     4  0.8001599431541093D-01,0.7536361371650956D-01,
     5  0.6389445258722484D-01,0.4772289428318182D-01,
     6  0.2915581601163874D-01,0.1139291011586660D-01/

        data r3/
     1  0.5501514789227983D-06,0.2196372157640899D-04,
     2  0.2601486355464068D-03,0.1590661738045614D-02,
     3  0.6270490111902275D-02,0.1797607105314628D-01,
     4  0.4041128964049613D-01,0.7502424749278209D-01,
     5  0.1195333558513453D+00,0.1686066833248838D+00,
     6  0.2163677612502772D+00,0.2588703807864688D+00,
     7  0.2948782942576631D+00,0.3250586597887873D+00,
     8  0.3508629756453048D+00,0.3738963870739845D+00,
     9  0.3957552022696601D+00,0.4180833380637550D+00,
     *  0.4426996377409523D+00,0.4717317133566876D+00,
     1  0.5076583161286358D+00,0.5530012305325799D+00,
     2  0.6093081809680398D+00,0.6756332746551345D+00,
     3  0.7478458833997212D+00,0.8197109077452231D+00,
     4  0.8847631687130651D+00,0.9376724231213696D+00,
     5  0.9749490457891534D+00,0.9953291682229628D+00/

        data w3/
     1  0.2611332065019828D-05,0.6438395340011701D-04,
     2  0.5437704745240752D-03,0.2490295421974482D-02,
     3  0.7506647399807306D-02,0.1656469222970677D-01,
     4  0.2857546717729153D-01,0.4025805847669194D-01,
     5  0.4782449916061857D-01,0.4929666560814074D-01,
     6  0.4555542313189689D-01,0.3926116740405005D-01,
     7  0.3288645406370514D-01,0.2772630568328544D-01,
     8  0.2415185418287473D-01,0.2217975554840835D-01,
     9  0.2180963356646054D-01,0.2314737314876106D-01,
     *  0.2643852794267367D-01,0.3204276008818766D-01,
     1  0.4024768491012314D-01,0.5072534661065237D-01,
     2  0.6174579523536927D-01,0.7019445693853178D-01,
     3  0.7315652053693238D-01,0.6947337326226997D-01,
     4  0.5973233956759084D-01,0.4549749211158631D-01,
     5  0.2884669073150596D-01,0.1205395410091390D-01/

        data r4/
     1  0.3297466873412479D-06,0.1381361370064898D-04,
     2  0.1711843929220451D-03,0.1091593420521096D-02,
     3  0.4470013087924967D-02,0.1325406222498140D-01,
     4  0.3069300651696182D-01,0.5852602890940691D-01,
     5  0.9569969951305334D-01,0.1388024702839367D+00,
     6  0.1838592707671512D+00,0.2279772835086871D+00,
     7  0.2698802909812153D+00,0.3096112882345654D+00,
     8  0.3480820079719525D+00,0.3867912657034145D+00,
     9  0.4276795783653907D+00,0.4729667475797029D+00,
     *  0.5248060945240692D+00,0.5846231285945156D+00,
     1  0.6522134748193047D+00,0.7250877419191858D+00,
     2  0.7986702757480167D+00,0.8673324825200674D+00,
     3  0.9255849954449382D+00,0.9688838276189572D+00,
     *  0.9940060961155116D+00/

        data w4/
     1  0.1581622659103543D-05,0.4113699853537371D-04,
     2  0.3653023247678219D-03,0.1752772176421009D-02,
     3  0.5511713808890548D-02,0.1263209256312257D-01,
     4  0.2256056798540039D-01,0.3292163435228208D-01,
     5  0.4082268707206814D-01,0.4469029794006911D-01,
     6  0.4493346087879749D-01,0.4309676768433979D-01,
     7  0.4072462170140529D-01,0.3889382816064716D-01,
     8  0.3829943850540361D-01,0.3944480729933166D-01,
     9  0.4270716505563387D-01,0.4823712419116135D-01,
     *  0.5569544258857516D-01,0.6391736455092232D-01,
     1  0.7084510054325236D-01,0.7412152416916533D-01,
     2  0.7208413966069603D-01,0.6431199868221032D-01,
     3  0.5143006311530243D-01,0.3462815560769749D-01,
     *  0.1532921076124220D-01/

        data r5/
     1  0.3245735122398547D-06,0.1340075350906127D-04,
     2  0.1634045061753075D-03,0.1024043128559201D-02,
     3  0.4122072113422181D-02,0.1204586981978045D-01,
     4  0.2765206690702527D-01,0.5275108736111824D-01,
     5  0.8732155718052196D-01,0.1298384681438392D+00,
     6  0.1783022986527552D+00,0.2311926076164301D+00,
     7  0.2879421592141206D+00,0.3489302665185113D+00,
     8  0.4150899587513199D+00,0.4871973600314703D+00,
     9  0.5650091821381723D+00,0.6466333990446034D+00,
     *  0.7285374524249186D+00,0.8062089958396365D+00,
     1  0.8750602447359000D+00,0.9311539536431448D+00,
     2  0.9715884562650833D+00,0.9945692561510447D+00/

        data w5/
     1  0.1552327210844114D-05,0.3970187969791374D-04,
     2  0.3459354211804380D-03,0.1626303718526055D-02,
     3  0.5015381342825278D-02,0.1133333390081148D-01,
     4  0.2019555032269095D-01,0.2999210363104312D-01,
     5  0.3887530099503887D-01,0.4580845431571964D-01,
     6  0.5085806063574614D-01,0.5482840132888167D-01,
     7  0.5874056358542536D-01,0.6340331364457744D-01,
     8  0.6906079835716869D-01,0.7512513540606213D-01,
     9  0.8018309347622900D-01,0.8246436292751699D-01,
     *  0.8057548724657570D-01,0.7399054706776181D-01,
     1  0.6304944169627073D-01,0.4865343681835463D-01,
     2  0.3191893522298818D-01,0.1391480473169695D-01/

        data r6/
     1  0.3235955210219312D-06,0.1333040356588220D-04,
     2  0.1626555824574076D-03,0.1022234921114556D-02,
     3  0.4139719184293033D-02,0.1224095735244652D-01,
     4  0.2868248257778606D-01,0.5645538566936932D-01,
     5  0.9747637180895232D-01,0.1525141034611335D+00,
     6  0.2214560711648077D+00,0.3033917590531211D+00,
     7  0.3963127582767359D+00,0.4967650174738328D+00,
     8  0.5999414737720175D+00,0.7002860002232959D+00,
     9  0.7922532319342742D+00,0.8709018977717835D+00,
     *  0.9322855267793159D+00,0.9738577189608327D+00,
     *  0.9953495995318043D+00/

        data w6/
     1  0.1546234507645722D-05,0.3948657240411641D-04,
     2  0.3446304483736304D-03,0.1627632575380070D-02,
     3  0.5071674948500139D-02,0.1170901336503211D-01,
     4  0.2168529002309401D-01,0.3418027940472608D-01,
     5  0.4798385518704894D-01,0.6207088440845878D-01,
     6  0.7566727462980076D-01,0.8787663589667691D-01,
     7  0.9738830209488858D-01,0.1027038680141876D+00,
     8  0.1027094942007432D+00,0.9704541086865896D-01,
     9  0.8605949534005668D-01,0.7057447979826621D-01,
     *  0.5175501069017673D-01,0.3131715859289026D-01,
     *  0.1218857670612857D-01/
c        
c         %%%%%%%%%%%%%%%%%%%%%

ccc        call prin2('alpha = *',alpha,1)
        dtwo = 2
        if (alpha .ge. 1/dtwo) then
        do 1000 i=1,nqs(6)
        roots(i) =r6(i)
        weights(i) = w6(i)
 1000 continue

        nquads = nqs(6)
        endif

        if (alpha .ge. dtwo**(-2)) then
        if (alpha .lt. dtwo**(-1)) then

        do 2000 i=1,nqs(5)
        roots(i) =r5(i)
        weights(i) = w5(i)
 2000 continue

        nquads = nqs(5)

        endif
        endif


        if (alpha .ge. dtwo**(-3)) then
        if (alpha .lt. dtwo**(-2)) then

        do 3000 i=1,nqs(4)
        roots(i) =r4(i)
        weights(i) = w4(i)
 3000 continue

        nquads = nqs(4)

        endif
        endif


        if (alpha .ge. dtwo**(-4)) then
        if (alpha .lt. dtwo**(-3)) then

        do 4000 i=1,nqs(3)
        roots(i) =r3(i)
        weights(i) = w3(i)
 4000 continue

        nquads = nqs(3)

        endif
        endif


        if (alpha .ge. dtwo**(-5)) then
        if (alpha .lt. dtwo**(-4)) then

        do 5000 i=1,nqs(2)
        roots(i) =r2(i)
        weights(i) = w2(i)
 5000 continue

        nquads = nqs(2)

        endif
        endif


        if (alpha .lt. dtwo**(-5)) then

        do 6000 i=1,nqs(1)
        roots(i) =r1(i)
        weights(i) = w1(i)
 6000 continue

        nquads = nqs(1)

        endif

        return
        end
c
c
c
c
c
        subroutine lap_load_quads23(alpha,roots,weights,nquads)
        implicit real *8 (a-h,o-z)
        dimension r1(32),r2(31),r3(29),r4(26),r5(23),r6(20),
     1          w1(32),w2(31),w3(29),w4(26),w5(23),w6(20),nqs(6),
     2          roots(1),weights(1)

        data nqs/32,31,29,26,23,20/


        data r1/
     1  0.5436536720474311D-06,0.2206515003875207D-04,
     2  0.2663376772622917D-03,0.1663926293724576D-02,
     3  0.6724138090532208D-02,0.1985433160241590D-01,
     4  0.4626619559206305D-01,0.8967838006834557D-01,
     5  0.1499633983032681D+00,0.2215100956677301D+00,
     6  0.2921788978017841D+00,0.3447349972641436D+00,
     7  0.3751492132847745D+00,0.4018775717484933D+00,
     8  0.4222324833302102D+00,0.4361664238751991D+00,
     9  0.4463250251496274D+00,0.4544560125156733D+00,
     *  0.4616612266323788D+00,0.4687795059926307D+00,
     1  0.4766435813808924D+00,0.4863396923632861D+00,
     2  0.4996526817200052D+00,0.5199807179827140D+00,
     3  0.5536723032237626D+00,0.6083917868489078D+00,
     4  0.6834457361457694D+00,0.7670981079774494D+00,
     5  0.8469780476200676D+00,0.9145698967624621D+00,
     6  0.9644201530656850D+00,0.9931625642082524D+00/

        data w1/
     1  0.2588836086842587D-05,0.6507126808727398D-04,
     2  0.5620750912158650D-03,0.2643251598446899D-02,
     3  0.8230543283481543D-02,0.1894017464974429D-01,
     4  0.3454475735395470D-01,0.5228446587831890D-01,
     5  0.6736879657053555D-01,0.7370108869095690D-01,
     6  0.6457383828908991D-01,0.3885419304490499D-01,
     7  0.2726657870758232D-01,0.2432756063365420D-01,
     8  0.1665801071703344D-01,0.1167659904062894D-01,
     9  0.8923033545077985D-02,0.7516288912566517D-02,
     *  0.7029996145103139D-02,0.7341418744580687D-02,
     1  0.8561050279825679D-02,0.1111422890552695D-01,
     2  0.1605082333652171D-01,0.2566278187860524D-01,
     3  0.4320482423793907D-01,0.6620376041153489D-01,
     4  0.8168509169373355D-01,0.8350770552305157D-01,
     5  0.7485717725333291D-01,0.5943394253419526D-01,
     6  0.3971179757291615D-01,0.1749648537176614D-01/

        data r2/
     1  0.5347424059600907D-06,0.2168641562064174D-04,
     2  0.2612074141322226D-03,0.1626484807540688D-02,
     3  0.6542635046620713D-02,0.1919561929500003D-01,
     4  0.4433859869753007D-01,0.8493598274660809D-01,
     5  0.1399853127142598D+00,0.2037379252495672D+00,
     6  0.2669882773716296D+00,0.3209859413737631D+00,
     7  0.3621054201742005D+00,0.3920628167981924D+00,
     8  0.4143147307662959D+00,0.4317597476974344D+00,
     9  0.4464364141912605D+00,0.4598406276048893D+00,
     *  0.4732474570497724D+00,0.4879706416370231D+00,
     1  0.5056285198455563D+00,0.5284841786800461D+00,
     2  0.5597760554754921D+00,0.6033593715348079D+00,
     3  0.6612979088987668D+00,0.7305003727966398D+00,
     4  0.8033088469276928D+00,0.8714758902349943D+00,
     5  0.9285901600699341D+00,0.9703796269207432D+00,
     *  0.9943226926454300D+00/

        data w2/
     1  0.2546401467106677D-05,0.6392161341893804D-04,
     2  0.5505347593101747D-03,0.2577421445668993D-02,
     3  0.7973137045018563D-02,0.1816814299111572D-01,
     4  0.3265133413446338D-01,0.4839039514645942D-01,
     5  0.6073072227111483D-01,0.6517325421305024D-01,
     6  0.5977168088244264D-01,0.4761853785660820D-01,
     7  0.3498484973029858D-01,0.2554284495408675D-01,
     8  0.1944496141608004D-01,0.1577840594389499D-01,
     9  0.1381887130558261D-01,0.1319793775232485D-01,
     *  0.1383044340827390D-01,0.1587995020433016D-01,
     1  0.1980618263533232D-01,0.2644646614093084D-01,
     2  0.3681974643054297D-01,0.5076579248223170D-01,
     3  0.6458380614667539D-01,0.7246061836734938D-01,
     4  0.7174312987950405D-01,0.6352158272571282D-01,
     5  0.5000270447497672D-01,0.3316281659123926D-01,
     *  0.1453726065049445D-01/

        data r3/
     1  0.4637109410633896D-06,0.1882483721451764D-04,
     2  0.2268525630644963D-03,0.1412661829347920D-02,
     3  0.5680566147925088D-02,0.1665117258853870D-01,
     4  0.3839168436072837D-01,0.7332539004939046D-01,
     5  0.1204094765780418D+00,0.1749025181345342D+00,
     6  0.2301854297400787D+00,0.2807360078116089D+00,
     7  0.3241006649514624D+00,0.3605511871297653D+00,
     8  0.3916010354804756D+00,0.4189987231884792D+00,
     9  0.4444642660979168D+00,0.4697578164090438D+00,
     *  0.4968092073421705D+00,0.5278245601568920D+00,
     1  0.5652793856869928D+00,0.6115515328591286D+00,
     2  0.6678646660657611D+00,0.7328021200904627D+00,
     3  0.8017645965389041D+00,0.8682340072164650D+00,
     4  0.9256903879897440D+00,0.9688304948939511D+00,
     *  0.9939865911831262D+00/

        data w3/
     1  0.2208767884451704D-05,0.5550081999055627D-04,
     2  0.4781879666880731D-03,0.2238404111010505D-02,
     3  0.6919261692265669D-02,0.1573825985328297D-01,
     4  0.2818052222914758D-01,0.4151821310022978D-01,
     5  0.5182742301676792D-01,0.5599976278417271D-01,
     6  0.5360850072531038D-01,0.4709367585476283D-01,
     7  0.3971037595081217D-01,0.3345974172915651D-01,
     8  0.2893443811365984D-01,0.2614586389533206D-01,
     9  0.2507631028780357D-01,0.2583049945300438D-01,
     *  0.2863951210348064D-01,0.3380883840139065D-01,
     1  0.4151485233418050D-01,0.5126224009503741D-01,
     2  0.6114015316067378D-01,0.6794306029119451D-01,
     3  0.6885971797544397D-01,0.6297244514009334D-01,
     4  0.5105082938619050D-01,0.3461740895673497D-01,
     *  0.1537379180429776D-01/

        data r4/
     1  0.4688127703345338D-06,0.1889931874145700D-04,
     2  0.2259453134160366D-03,0.1393342073726484D-02,
     3  0.5535513666132662D-02,0.1599378696513231D-01,
     4  0.3629382232330974D-01,0.6825180030456466D-01,
     5  0.1107453454705166D+00,0.1601467918741902D+00,
     6  0.2120662983795193D+00,0.2631032751268642D+00,
     7  0.3116124028573282D+00,0.3575370482325245D+00,
     8  0.4019144866278211D+00,0.4464629660735291D+00,
     9  0.4932968647589906D+00,0.5446053192572201D+00,
     *  0.6020974112323766D+00,0.6661773020065684D+00,
     1  0.7351785364603813D+00,0.8052567308533034D+00,
     2  0.8711784459466562D+00,0.9275194437286391D+00,
     3  0.9696278167928694D+00,0.9941425199241753D+00/

        data w4/
     1  0.2229858143475170D-05,0.5558216085610178D-04,
     2  0.4743009499775415D-03,0.2192633623885719D-02,
     3  0.6668873099414160D-02,0.1487091740386841D-01,
     4  0.2605514720243678D-01,0.3766548731942377D-01,
     5  0.4668422220077067D-01,0.5135509258101792D-01,
     6  0.5190483425759196D-01,0.4989948351401093D-01,
     7  0.4711857509464291D-01,0.4491050698446113D-01,
     8  0.4413503135201657D-01,0.4531914155952555D-01,
     9  0.4872256625964456D-01,0.5419371308085722D-01,
     *  0.6087234910124282D-01,0.6701228268259598D-01,
     1  0.7033480373931804D-01,0.6892754284471974D-01,
     2  0.6199538448078548D-01,0.4990318805720926D-01,
     3  0.3374992092300540D-01,0.1497618966857789D-01/

        data r5/
     1  0.4797907353965881D-06,0.1914465716737341D-04,
     2  0.2266303962064141D-03,0.1383779000833555D-02,
     3  0.5444139841967654D-02,0.1559340194915226D-01,
     4  0.3517433463659223D-01,0.6609608730388265D-01,
     5  0.1080321881533203D+00,0.1589764083271159D+00,
     6  0.2165129942861067D+00,0.2788686533964804D+00,
     7  0.3453390705858309D+00,0.4161222148782570D+00,
     8  0.4916859820254972D+00,0.5718533022928709D+00,
     9  0.6549859057063883D+00,0.7377483382434832D+00,
     *  0.8156143704808112D+00,0.8837805318056638D+00,
     1  0.9380451747618910D+00,0.9755445755756333D+00,
     *  0.9955275276832503D+00/

        data w5/
     1  0.2276925572610624D-05,0.5611340792180768D-04,
     2  0.4735017642130799D-03,0.2163851110397363D-02,
     3  0.6507380798255731D-02,0.1437866549485875D-01,
     4  0.2511634495929355D-01,0.3666141226126276D-01,
     5  0.4685263980097866D-01,0.5461478917471064D-01,
     6  0.6015598309740834D-01,0.6444103097315074D-01,
     7  0.6854367782123589D-01,0.7311677145228798D-01,
     8  0.7799622733339019D-01,0.8206979881051337D-01,
     9  0.8363004193185020D-01,0.8112396227277093D-01,
     *  0.7379261608598062D-01,0.6182791676546809D-01,
     1  0.4621087716071030D-01,0.2865170386316993D-01,
     *  0.1161241673459846D-01/

        data r6/
     1  0.4756598539520443D-06,0.1897901230600169D-04,
     2  0.2242456657361756D-03,0.1366930190103947D-02,
     3  0.5384749324567682D-02,0.1553986261015340D-01,
     4  0.3564982473585991D-01,0.6888491002357068D-01,
     5  0.1169977153995884D+00,0.1802926044678277D+00,
     6  0.2579326236226577D+00,0.3480482811747926D+00,
     7  0.4475552023412750D+00,0.5520647648075505D+00,
     8  0.6562254027287370D+00,0.7543829958620055D+00,
     9  0.8412163329185847D+00,0.9121552273759357D+00,
     *  0.9635980495461574D+00,0.9930215603771403D+00/

        data w6/
     1  0.2258003817729880D-05,0.5560396758573240D-04,
     2  0.4680313032887041D-03,0.2136531931940944D-02,
     3  0.6457050998509273D-02,0.1451204006057348D-01,
     4  0.2624986807228886D-01,0.4051165996320371D-01,
     5  0.5575801973090138D-01,0.7069020716342719D-01,
     6  0.8428247583626178D-01,0.9543728751138108D-01,
     7  0.1028390200124513D+00,0.1052761518824848D+00,
     8  0.1020917227212276D+00,0.9332974251198392D-01,
     9  0.7957185160761572D-01,0.6170407846925617D-01,
     *  0.4075939173344732D-01,0.1786700651835336D-01/
c        
c         %%%%%%%%%%%%%%%%%%%%%

ccc        call prin2('alpha = *',alpha,1)
        dtwo = 2
        if (alpha .ge. 1/dtwo) then
        do 1000 i=1,nqs(6)
        roots(i) =r6(i)
        weights(i) = w6(i)
 1000 continue

        nquads = nqs(6)
        endif

        if (alpha .ge. dtwo**(-2)) then
        if (alpha .lt. dtwo**(-1)) then

        do 2000 i=1,nqs(5)
        roots(i) =r5(i)
        weights(i) = w5(i)
 2000 continue

        nquads = nqs(5)

        endif
        endif


        if (alpha .ge. dtwo**(-3)) then
        if (alpha .lt. dtwo**(-2)) then

        do 3000 i=1,nqs(4)
        roots(i) =r4(i)
        weights(i) = w4(i)
 3000 continue

        nquads = nqs(4)

        endif
        endif


        if (alpha .ge. dtwo**(-4)) then
        if (alpha .lt. dtwo**(-3)) then

        do 4000 i=1,nqs(3)
        roots(i) =r3(i)
        weights(i) = w3(i)
 4000 continue

        nquads = nqs(3)

        endif
        endif


        if (alpha .ge. dtwo**(-5)) then
        if (alpha .lt. dtwo**(-4)) then

        do 5000 i=1,nqs(2)
        roots(i) =r2(i)
        weights(i) = w2(i)
 5000 continue

        nquads = nqs(2)

        endif
        endif


        if (alpha .lt. dtwo**(-5)) then

        do 6000 i=1,nqs(1)
        roots(i) =r1(i)
        weights(i) = w1(i)
 6000 continue

        nquads = nqs(1)

        endif

        return
        end
c
c
c
c
c
        subroutine lap_load_quads24(alpha,roots,weights,nquads)
        implicit real *8 (a-h,o-z)
        dimension r1(33),r2(31),r3(29),r4(26),r5(22),r6(20),
     1          w1(33),w2(31),w3(29),w4(26),w5(22),w6(20),nqs(6),
     2          roots(1),weights(1)

        data nqs/33,31,29,26,22,20/


        data r1/
     1  0.3741782341806244D-01,0.6778281630653795D-01,
     2  0.1111284726840989D+00,0.1705509832351944D-01,
     3  0.5984039002807298D-02,0.1707724966075256D+00,
     4  0.1507818167525744D-02,0.2431456244280386D+00,
     5  0.3186490950630096D+00,0.2431889635568469D-03,
     6  0.3846871597603363D+00,0.5226775122188226D+00,
     7  0.2016374438081528D-04,0.4324248549453259D+00,
     8  0.4948610942747535D-06,0.4630631081948744D+00,
     9  0.6517785616391574D+00,0.7139398736676397D+00,
     *  0.6076024055454228D+00,0.4827774203110844D+00,
     1  0.7867198054963687D+00,0.5793180905665986D+00,
     2  0.4963961090369411D+00,0.5609651649367198D+00,
     3  0.8586666933875959D+00,0.5301831984517208D+00,
     4  0.5481425703464485D+00,0.5150888886881447D+00,
     5  0.9206914982297994D+00,0.5383585989935004D+00,
     6  0.5066546897001365D+00,0.9668734093648903D+00,
     *  0.9936249176161080D+00/

        data w1/
     1  0.2526487416606227D-01,0.3598614745029968D-01,
     2  0.5144018108155354D-01,0.1546633204525195D-01,
     3  0.7173396325306201D-02,0.6719816798600634D-01,
     4  0.2377137856767293D-02,0.7591338669638119D-01,
     5  0.7288074523100949D-01,0.5123096127986694D-03,
     6  0.5755885734822561D-01,0.7426814321778733D-02,
     7  0.5950009454177037D-04,0.3825645425618768D-01,
     8  0.2359412065703884D-05,0.2416557557621382D-01,
     9  0.5358851195197580D-01,0.6936181926153457D-01,
     *  0.3526585075485977D-01,0.1606469502168779D-01,
     1  0.7417134872890772D-01,0.2238911556647729D-01,
     2  0.1161087899895049D-01,0.1504277092510422D-01,
     3  0.6820978717368881D-01,0.7706989261992448D-02,
     4  0.1100159643603776D-01,0.7874375542613850D-02,
     5  0.5486990304420604D-01,0.8795591646975631D-02,
     6  0.9151905174171488D-02,0.3690449483647350D-01,
     *  0.1630812621389287D-01/

        data r2/
     1  0.5199837149936321D-06,0.2126580153693187D-04,
     2  0.2580590051847552D-03,0.1617329401151716D-02,
     3  0.6542188732142251D-02,0.1929076366430339D-01,
     4  0.4478374553067801D-01,0.8629567801673609D-01,
     5  0.1433416624599801D+00,0.2109373773826889D+00,
     6  0.2806029920791991D+00,0.3432548268998980D+00,
     7  0.3933534109542883D+00,0.4308066624370946D+00,
     8  0.4586409278416898D+00,0.4801047779030571D+00,
     9  0.4976655374530936D+00,0.5131156162764486D+00,
     *  0.5278882121833343D+00,0.5433287819343739D+00,
     1  0.5609333674874611D+00,0.5826070270522576D+00,
     2  0.6109173150879452D+00,0.6489890921542722D+00,
     3  0.6991102755310684D+00,0.7599079681092074D+00,
     4  0.8252822040306169D+00,0.8871625415057942D+00,
     5  0.9386486086643700D+00,0.9752963583940497D+00,
     *  0.9953945895608946D+00/

        data w2/
     1  0.2480753849530676D-05,0.6285045378417664D-04,
     2  0.5456993941789401D-03,0.2572611248052853D-02,
     3  0.8006706710931715D-02,0.1835293463620551D-01,
     4  0.3321713460869176D-01,0.4973871025568184D-01,
     5  0.6350603626052351D-01,0.7022428913613920D-01,
     6  0.6751359620423606D-01,0.5682624230044361D-01,
     7  0.4340906386963016D-01,0.3205581455070409D-01,
     8  0.2416714364785225D-01,0.1916787060811772D-01,
     9  0.1624672760464975D-01,0.1488786772442123D-01,
     *  0.1487775212523712D-01,0.1624844823939555D-01,
     1  0.1927337109277290D-01,0.2450112533944784D-01,
     2  0.3266108717100986D-01,0.4392042521708604D-01,
     3  0.5612006748723254D-01,0.6439863944089048D-01,
     4  0.6494724453553439D-01,0.5764603878524107D-01,
     5  0.4458169599722808D-01,0.2843292330806362D-01,
     *  0.1188740129276658D-01/

        data r3/
     1  0.4414124387235530D-06,0.1814704243635486D-04,
     2  0.2211771505298278D-03,0.1391311818134536D-02,
     3  0.5644803482927635D-02,0.1667852996288096D-01,
     4  0.3874521887300771D-01,0.7458803382500427D-01,
     5  0.1236242340177057D+00,0.1816270982672504D+00,
     6  0.2422174814768390D+00,0.2994195917652340D+00,
     7  0.3497483244747935D+00,0.3925446045806236D+00,
     8  0.4289013577464109D+00,0.4605224372436341D+00,
     9  0.4891855685314719D+00,0.5166580910794063D+00,
     *  0.5447705632499916D+00,0.5754925000658630D+00,
     1  0.6109206365715812D+00,0.6530330142669244D+00,
     2  0.7029995341940723D+00,0.7600975794045960D+00,
     3  0.8209991277261854D+00,0.8803319762546541D+00,
     4  0.9321900476923804D+00,0.9714610829094674D+00,
     *  0.9944838663991281D+00/

        data w3/
     1  0.2108455147971465D-05,0.5371920139145466D-04,
     2  0.4685967187451754D-03,0.2217781005653144D-02,
     3  0.6922963447987391D-02,0.1589135120406388D-01,
     4  0.2873287085368581D-01,0.4286164724060973D-01,
     5  0.5447960906315799D-01,0.6041632378861561D-01,
     6  0.5972974163673378D-01,0.5408923056761796D-01,
     7  0.4648046383471761D-01,0.3931477061179925D-01,
     8  0.3369093148550942D-01,0.2984813409464830D-01,
     9  0.2777114101369765D-01,0.2747757645404914D-01,
     *  0.2907548749208474D-01,0.3272132635956134D-01,
     1  0.3847592183142467D-01,0.4596388933736050D-01,
     2  0.5386365472663558D-01,0.5976735446554416D-01,
     3  0.6111163381175544D-01,0.5654655051710686D-01,
     4  0.4630749346726599D-01,0.3162139117644127D-01,
     *  0.1409633613698818D-01/

        data r4/
     1  0.3098482870968050D-06,0.1321369340011192D-04,
     2  0.1666344237836556D-03,0.1081434422973171D-02,
     3  0.4510101903614672D-02,0.1363960171981574D-01,
     4  0.3228833319121732D-01,0.6309879935530548D-01,
     5  0.1059398245941637D+00,0.1577757377541505D+00,
     6  0.2141014236891052D+00,0.2708373395745529D+00,
     7  0.3254823181502046D+00,0.3772368242001322D+00,
     8  0.4265836133228924D+00,0.4748260651150451D+00,
     9  0.5237381298634736D+00,0.5752306761551474D+00,
     *  0.6308780822806552D+00,0.6912495861792423D+00,
     1  0.7552403369425741D+00,0.8198370381931386D+00,
     2  0.8806106619445518D+00,0.9326999371961476D+00,
     3  0.9717559970315892D+00,0.9945480587293430D+00/

        data w4/
     1  0.1491846696654181D-05,0.3958234184166517D-04,
     2  0.3585002489997921D-03,0.1755732197312285D-02,
     3  0.5646246027655745D-02,0.1328270775748923D-01,
     4  0.2448538781528063D-01,0.3711611444775088D-01,
     5  0.4803288929714857D-01,0.5485312606691172D-01,
     6  0.5709854556765504D-01,0.5595593786510422D-01,
     7  0.5320906940649810D-01,0.5038838841830944D-01,
     8  0.4852794834155573D-01,0.4825774447342569D-01,
     9  0.4989255245838107D-01,0.5336651821392504D-01,
     *  0.5803553024447963D-01,0.6252817146402894D-01,
     1  0.6494220451614825D-01,0.6349509068742511D-01,
     2  0.5722546016747492D-01,0.4621843345491009D-01,
     3  0.3135030617020116D-01,0.1393632050339037D-01/

        data r5/
     1  0.4826379323120386D-06,0.1952366681478281D-04,
     2  0.2341152347666003D-03,0.1446670140832628D-02,
     3  0.5753902649888168D-02,0.1664368364733422D-01,
     4  0.3787764528583240D-01,0.7173540367774196D-01,
     5  0.1180247054756594D+00,0.1745684741209474D+00,
     6  0.2385649561733131D+00,0.3077846989492891D+00,
     7  0.3810741597140075D+00,0.4581648958502997D+00,
     8  0.5389708325664158D+00,0.6226612356257828D+00,
     9  0.7069692433653891D+00,0.7881596222265524D+00,
     *  0.8616487013229401D+00,0.9228515781600421D+00,
     1  0.9678565586798770D+00,0.9938195672855619D+00/

        data w5/
     1  0.2297246374831517D-05,0.5748389748131726D-04,
     2  0.4920550506586070D-03,0.2278734308501111D-02,
     3  0.6937228900589915D-02,0.1550269269310859D-01,
     4  0.2736365989488286D-01,0.4030974702807134D-01,
     5  0.5188041494027717D-01,0.6071901748083797D-01,
     6  0.6689416853364757D-01,0.7135668764381571D-01,
     7  0.7518981938136059D-01,0.7899393732859122D-01,
     8  0.8249223906362054D-01,0.8451127586259330D-01,
     9  0.8347060979959233D-01,0.7813101720661887D-01,
     *  0.6807155535373215D-01,0.5367601587357755D-01,
     1  0.3585647985694771D-01,0.1581286265511875D-01/

        data r6/
     1  0.4255179756283601D-06,0.1722799292630540D-04,
     2  0.2066616330403500D-03,0.1279465578353538D-02,
     3  0.5118779105309694D-02,0.1499282491935620D-01,
     4  0.3486348525789442D-01,0.6815306524191311D-01,
     5  0.1168212549140445D+00,0.1811717992808818D+00,
     6  0.2601147682002324D+00,0.3513894092300778D+00,
     7  0.4515798160587443D+00,0.5561688424288540D+00,
     8  0.6598839077766567D+00,0.7572642857825773D+00,
     9  0.8431999151530676D+00,0.9132972585778837D+00,
     *  0.9640827662799887D+00,0.9931156179011471D+00/

        data w6/
     1  0.2026115919225865D-05,0.5072974173061039D-04,
     2  0.4345128015655345D-03,0.2019774560029457D-02,
     3  0.6216166742109500D-02,0.1421530967698233D-01,
     4  0.2611286947820424D-01,0.4079522700038727D-01,
     5  0.5658596805092533D-01,0.7192930957724335D-01,
     6  0.8558037514816951D-01,0.9639947095284176D-01,
     7  0.1032224533902882D+00,0.1050666369572347D+00,
     8  0.1014440458280793D+00,0.9245901594798987D-01,
     9  0.7867678048984912D-01,0.6093731253179547D-01,
     *  0.4022516152761659D-01,0.1762685348103860D-01/
c        
c         %%%%%%%%%%%%%%%%%%%%%

ccc        call prin2('alpha = *',alpha,1)
        dtwo = 2
        if (alpha .ge. 1/dtwo) then
        do 1000 i=1,nqs(6)
        roots(i) =r6(i)
        weights(i) = w6(i)
 1000 continue

        nquads = nqs(6)
        endif

        if (alpha .ge. dtwo**(-2)) then
        if (alpha .lt. dtwo**(-1)) then

        do 2000 i=1,nqs(5)
        roots(i) =r5(i)
        weights(i) = w5(i)
 2000 continue

        nquads = nqs(5)

        endif
        endif


        if (alpha .ge. dtwo**(-3)) then
        if (alpha .lt. dtwo**(-2)) then

        do 3000 i=1,nqs(4)
        roots(i) =r4(i)
        weights(i) = w4(i)
 3000 continue

        nquads = nqs(4)

        endif
        endif


        if (alpha .ge. dtwo**(-4)) then
        if (alpha .lt. dtwo**(-3)) then

        do 4000 i=1,nqs(3)
        roots(i) =r3(i)
        weights(i) = w3(i)
 4000 continue

        nquads = nqs(3)

        endif
        endif


        if (alpha .ge. dtwo**(-5)) then
        if (alpha .lt. dtwo**(-4)) then

        do 5000 i=1,nqs(2)
        roots(i) =r2(i)
        weights(i) = w2(i)
 5000 continue

        nquads = nqs(2)

        endif
        endif


        if (alpha .lt. dtwo**(-5)) then

        do 6000 i=1,nqs(1)
        roots(i) =r1(i)
        weights(i) = w1(i)
 6000 continue

        nquads = nqs(1)

        endif

        return
        end
c
c
c
c
c
        subroutine lap_load_quads25(alpha,roots,weights,nquads)
        implicit real *8 (a-h,o-z)
        dimension r1(32),r2(31),r3(28),r4(25),r5(22),r6(20),
     1          w1(32),w2(31),w3(28),w4(25),w5(22),w6(20),nqs(6),
     2          roots(1),weights(1)

        data nqs/32,31,28,25,22,20/


        data r1/
     1  0.2244616469603041D-06,0.1110729266201527D-04,
     2  0.1558091721476007D-03,0.1096770092522937D-02,
     3  0.4878404790993447D-02,0.1555778807657180D-01,
     4  0.3855931189906911D-01,0.7855269256777912D-01,
     5  0.1370070486752347D+00,0.2107825315490306D+00,
     6  0.2921301711161437D+00,0.3700014076260994D+00,
     7  0.4341162947071081D+00,0.4811584339806742D+00,
     8  0.5139352293240593D+00,0.5364562579496742D+00,
     9  0.5523891127413427D+00,0.5644071204338016D+00,
     *  0.5742346388714995D+00,0.5830386250963613D+00,
     1  0.5917513424951173D+00,0.6013046871944095D+00,
     2  0.6128696772500947D+00,0.6282213654105074D+00,
     3  0.6503285837716748D+00,0.6837800770457632D+00,
     4  0.7327905245506637D+00,0.7953039528429985D+00,
     5  0.8617290265324502D+00,0.9214859289507708D+00,
     6  0.9669820748133494D+00,0.9936247083637778D+00/

        data w1/
     1  0.1121315293965638D-05,0.3464375266211377D-04,
     2  0.3495548185002944D-03,0.1859339529958003D-02,
     3  0.6390106183771906D-02,0.1591145947167849D-01,
     4  0.3091450453441900D-01,0.4934128446098596D-01,
     5  0.6704274994707400D-01,0.7919824614431574D-01,
     6  0.8159517796574392D-01,0.7229220758428959D-01,
     7  0.5541213041137104D-01,0.3924952560355368D-01,
     8  0.2698136034389853D-01,0.1868699896692735D-01,
     9  0.1362129487721882D-01,0.1069468908909701D-01,
     *  0.9149554606039087D-02,0.8610127596241366D-02,
     1  0.8966914837541899D-02,0.1032920271047967D-01,
     2  0.1308652675437818D-01,0.1809965348289778D-01,
     3  0.2690144730612449D-01,0.4083009212117863D-01,
     4  0.5688783785240325D-01,0.6638513609012761D-01,
     5  0.6464141135694395D-01,0.5361619703527775D-01,
     6  0.3662298547918990D-01,0.1629651777041702D-01/

        data r2/
     1  0.5481818076904533D-06,0.2224998455653179D-04,
     2  0.2678062878856049D-03,0.1662623439474163D-02,
     3  0.6650128895924223D-02,0.1934989858150264D-01,
     4  0.4424530912293649D-01,0.8391413103946084D-01,
     5  0.1375128516633060D+00,0.2012699388318729D+00,
     6  0.2700259252131252D+00,0.3377370386601366D+00,
     7  0.3981122561804112D+00,0.4472080793129268D+00,
     8  0.4850209660385028D+00,0.5139807490939737D+00,
     9  0.5368474091035494D+00,0.5558627865942899D+00,
     *  0.5727593187785283D+00,0.5889872085222634D+00,
     1  0.6059349524998674D+00,0.6251267644940559D+00,
     2  0.6484129966018003D+00,0.6780789687605436D+00,
     3  0.7165194543187614D+00,0.7648590492067397D+00,
     4  0.8208645861390986D+00,0.8785806587346881D+00,
     5  0.9306395618604057D+00,0.9706821044289356D+00,
     *  0.9943220670198187D+00/

        data w2/
     1  0.2611142435537957D-05,0.6558703485945460D-04,
     2  0.5639790602886665D-03,0.2627552443318343D-02,
     3  0.8057232664553417D-02,0.1813220282905962D-01,
     4  0.3211734425270952D-01,0.4706923804775839D-01,
     5  0.5946584470957783D-01,0.6718399082716577D-01,
     6  0.6931888247961849D-01,0.6500028971890645D-01,
     7  0.5505898361599688D-01,0.4316728214077768D-01,
     8  0.3290296084081504D-01,0.2548857593029246D-01,
     9  0.2061663076668787D-01,0.1769867393823965D-01,
     *  0.1633243697774146D-01,0.1635154370331524D-01,
     1  0.1779470681581213D-01,0.2089474891490395D-01,
     2  0.2606211693315092D-01,0.3369159337046954D-01,
     3  0.4341938442769386D-01,0.5289783311635652D-01,
     4  0.5807249022746976D-01,0.5608770334756123D-01,
     5  0.4695677411046247D-01,0.3239726941208302D-01,
     *  0.1450353619991883D-01/

        data r3/
     1  0.6423750325061252D-06,0.2564232988706757D-04,
     2  0.3038697951641652D-03,0.1859419261021993D-02,
     3  0.7338815011341645D-02,0.2109208653141621D-01,
     4  0.4766944178067650D-01,0.8936164930763658D-01,
     5  0.1445441664900204D+00,0.2080649026218318D+00,
     6  0.2733776543875165D+00,0.3349704400279887D+00,
     7  0.3897244422325202D+00,0.4369176251410677D+00,
     8  0.4774370788121656D+00,0.5129030596172675D+00,
     9  0.5451456748135167D+00,0.5760565207454392D+00,
     *  0.6076065770837583D+00,0.6418445541998460D+00,
     1  0.6807410664604861D+00,0.7257268739227911D+00,
     2  0.7768554007337560D+00,0.8319579249242075D+00,
     3  0.8866151375884861D+00,0.9352534979219449D+00,
     4  0.9726016423965452D+00,0.9946883439558595D+00/

        data w3/
     1  0.3048679923440834D-05,0.7517421814093833D-04,
     2  0.6353572140930340D-03,0.2912733733667881D-02,
     3  0.8797794212559291D-02,0.1951315855249191D-01,
     4  0.3405434504549830D-01,0.4905192809321907D-01,
     5  0.6041485714534800D-01,0.6549133210507593D-01,
     6  0.6419237985588111D-01,0.5847502511105579D-01,
     7  0.5092755470106571D-01,0.4362241232014367D-01,
     8  0.3769519920635830D-01,0.3354362816184050D-01,
     9  0.3125613421566675D-01,0.3089380829290080D-01,
     *  0.3254996891283748D-01,0.3626177069808337D-01,
     1  0.4178443193055449D-01,0.4821265282437616D-01,
     2  0.5368456155438433D-01,0.5576317493555221D-01,
     3  0.5259573136984077D-01,0.4378587412548662D-01,
     4  0.3024173693475208D-01,0.1356422584920208D-01/

        data r4/
     1  0.5048420110925455D-06,0.2051827441293548D-04,
     2  0.2471636555397400D-03,0.1534896934407123D-02,
     3  0.6138094898418407D-02,0.1785094256880432D-01,
     4  0.4079490316079493D-01,0.7734675662852448D-01,
     5  0.1267249747798196D+00,0.1852839691741289D+00,
     6  0.2481765515376879D+00,0.3111764757530947D+00,
     7  0.3717061269503468D+00,0.4289408148811338D+00,
     8  0.4834135934377701D+00,0.5365230894477246D+00,
     9  0.5900663098946368D+00,0.6457273451357610D+00,
     *  0.7044413878815518D+00,0.7657197262862685D+00,
     1  0.8272706188906483D+00,0.8852621930484699D+00,
     2  0.9351679826183068D+00,0.9727410917009203D+00,
     *  0.9947324649206918D+00/

        data w4/
     1  0.2405504125915217D-05,0.6050457652261245D-04,
     2  0.5206411867165781D-03,0.2425741886374744D-02,
     3  0.7434514315632020D-02,0.1671701776986850D-01,
     4  0.2959352957773233D-01,0.4337391275727835D-01,
     5  0.5474379591153246D-01,0.6153343155485568D-01,
     6  0.6353008543004557D-01,0.6203971757615284D-01,
     7  0.5888977806511009D-01,0.5567686941344429D-01,
     8  0.5350934204536612D-01,0.5301740835651843D-01,
     9  0.5436045812369703D-01,0.5712962223628971D-01,
     *  0.6022438926327858D-01,0.6194580886826667D-01,
     1  0.6050189681072902D-01,0.5470957237021485D-01,
     2  0.4438220350136158D-01,0.3021592374439434D-01,
     *  0.1346142915449171D-01/

        data r5/
     1  0.4868815805334776D-06,0.1976669429646944D-04,
     2  0.2377018762352264D-03,0.1473126990796648D-02,
     3  0.5880090860546593D-02,0.1708553846673657D-01,
     4  0.3909503242076269D-01,0.7448780474625430D-01,
     5  0.1232884204091743D+00,0.1833076209255841D+00,
     6  0.2514592412899979D+00,0.3250245019634974D+00,
     7  0.4022522449287079D+00,0.4822445043754465D+00,
     8  0.5643392915080508D+00,0.6473255677555328D+00,
     9  0.7289175402567902D+00,0.8057933340068547D+00,
     *  0.8741226825310395D+00,0.9302235151201735D+00,
     1  0.9710500580487535D+00,0.9944468126625850D+00/

        data w5/
     1  0.2319580734366290D-05,0.5825977667793132D-04,
     2  0.5002487941660100D-03,0.2325084421036994D-02,
     3  0.7112978536619656D-02,0.1599966978872222D-01,
     4  0.2847021811611986D-01,0.4231213169141981D-01,
     5  0.5489935334738548D-01,0.6460215852533612D-01,
     6  0.7123800323247772D-01,0.7560311196913217D-01,
     7  0.7871377586295376D-01,0.8117460912890151D-01,
     8  0.8282276358345787D-01,0.8276794949929078D-01,
     9  0.7985251120238742D-01,0.7324898313262025D-01,
     *  0.6278778931205357D-01,0.4890231334687986D-01,
     1  0.3238975209016683D-01,0.1421601506145979D-01/

        data r6/
     1  0.4296176008459255D-06,0.1768038225703905D-04,
     2  0.2153364610088213D-03,0.1351910061746199D-02,
     3  0.5475738510017789D-02,0.1620393963580829D-01,
     4  0.3797396874176939D-01,0.7460106759681787D-01,
     5  0.1281162327221613D+00,0.1984528416928000D+00,
     6  0.2837496441155955D+00,0.3807377652659484D+00,
     7  0.4850012887492322D+00,0.5912821421907574D+00,
     8  0.6939670022268243D+00,0.7876427706982044D+00,
     9  0.8675545778533243D+00,0.9299544080252617D+00,
     *  0.9725082733891164D+00,0.9949906684858660D+00/

        data w6/
     1  0.2052954335420894D-05,0.5234012807109210D-04,
     2  0.4559023348145186D-03,0.2152005303050176D-02,
     3  0.6711955843251462D-02,0.1551094459197863D-01,
     4  0.2868586883071594D-01,0.4491122566639327D-01,
     5  0.6210137776515058D-01,0.7825466325723517D-01,
     6  0.9179042820649375D-01,0.1014451808552837D+00,
     7  0.1061970054798782D+00,0.1054193996784601D+00,
     8  0.9904048712870223D-01,0.8751316399151823D-01,
     9  0.7167979182028826D-01,0.5271880745410152D-01,
     *  0.3233865158482651D-01,0.1301874712545123D-01/
c        
c         %%%%%%%%%%%%%%%%%%%%%

ccc        call prin2('alpha = *',alpha,1)
        dtwo = 2
        if (alpha .ge. 1/dtwo) then
        do 1000 i=1,nqs(6)
        roots(i) =r6(i)
        weights(i) = w6(i)
 1000 continue

        nquads = nqs(6)
        endif

        if (alpha .ge. dtwo**(-2)) then
        if (alpha .lt. dtwo**(-1)) then

        do 2000 i=1,nqs(5)
        roots(i) =r5(i)
        weights(i) = w5(i)
 2000 continue

        nquads = nqs(5)

        endif
        endif


        if (alpha .ge. dtwo**(-3)) then
        if (alpha .lt. dtwo**(-2)) then

        do 3000 i=1,nqs(4)
        roots(i) =r4(i)
        weights(i) = w4(i)
 3000 continue

        nquads = nqs(4)

        endif
        endif


        if (alpha .ge. dtwo**(-4)) then
        if (alpha .lt. dtwo**(-3)) then

        do 4000 i=1,nqs(3)
        roots(i) =r3(i)
        weights(i) = w3(i)
 4000 continue

        nquads = nqs(3)

        endif
        endif


        if (alpha .ge. dtwo**(-5)) then
        if (alpha .lt. dtwo**(-4)) then

        do 5000 i=1,nqs(2)
        roots(i) =r2(i)
        weights(i) = w2(i)
 5000 continue

        nquads = nqs(2)

        endif
        endif


        if (alpha .lt. dtwo**(-5)) then

        do 6000 i=1,nqs(1)
        roots(i) =r1(i)
        weights(i) = w1(i)
 6000 continue

        nquads = nqs(1)

        endif

        return
        end
c
c
c
c
c
        subroutine lap_load_quads26(alpha,roots,weights,nquads)
        implicit real *8 (a-h,o-z)
        dimension r1(32),r2(31),r3(28),r4(25),r5(22),r6(20),
     1          w1(32),w2(31),w3(28),w4(25),w5(22),w6(20),nqs(6),
     2          roots(1),weights(1)

        data nqs/32,31,28,25,22,20/


        data r1/
     1  0.5837876128823866D-06,0.2372582579915730D-04,
     2  0.2861847811817265D-03,0.1782548538540850D-02,
     3  0.7166016411418857D-02,0.2101696881877009D-01,
     4  0.4863716094525021D-01,0.9379810270891253D-01,
     5  0.1568597060577448D+00,0.2342031124222128D+00,
     6  0.3187515194399580D+00,0.4009321925223221D+00,
     7  0.4706475479122493D+00,0.5221528113384804D+00,
     8  0.5578596562251583D+00,0.5832805358034833D+00,
     9  0.6019464662105250D+00,0.6161194268217042D+00,
     *  0.6275058331723569D+00,0.6373688846740622D+00,
     1  0.6467042249894927D+00,0.6564318643250258D+00,
     2  0.6675884863800610D+00,0.6815800604844175D+00,
     3  0.7005627875834986D+00,0.7278122580448841D+00,
     4  0.7669686851998415D+00,0.8182898554834904D+00,
     5  0.8753389921618713D+00,0.9284643832654573D+00,
     6  0.9697278740651424D+00,0.9941367272145861D+00/

        data w1/
     1  0.2781252537136129D-05,0.6997790423616550D-04,
     2  0.6034375020867443D-03,0.2824245989615989D-02,
     3  0.8727379127773220D-02,0.1989982679046502D-01,
     4  0.3600306084081024D-01,0.5438735970092824D-01,
     5  0.7112178737833893D-01,0.8237014024167698D-01,
     6  0.8510598862041781D-01,0.7748627097320963D-01,
     7  0.6088433677477495D-01,0.4265186876732734D-01,
     8  0.2978230111779781D-01,0.2160933170503439D-01,
     9  0.1609254962508511D-01,0.1253688482809979D-01,
     *  0.1044315553220122D-01,0.9445690767742754D-02,
     1  0.9375279599622999D-02,0.1024877933328119D-01,
     2  0.1229239967616101D-01,0.1604188635983538D-01,
     3  0.2247481418927987D-01,0.3269841563029020D-01,
     4  0.4574320187083685D-01,0.5574228791978200D-01,
     5  0.5665202808780795D-01,0.4826493557488693D-01,
     6  0.3344037025148278D-01,0.1497722606657335D-01/

        data r2/
     1  0.4690061530223821D-06,0.1934721445226495D-04,
     2  0.2364945829713231D-03,0.1490881839245954D-02,
     3  0.6058132516022445D-02,0.1793124922113264D-01,
     4  0.4179899208861163D-01,0.8103049389711506D-01,
     5  0.1359640753500966D+00,0.2034965723868567D+00,
     6  0.2778824782231333D+00,0.3520866054933861D+00,
     7  0.4194456781833978D+00,0.4756362962857711D+00,
     8  0.5198674263887850D+00,0.5540341256692898D+00,
     9  0.5808678260275004D+00,0.6028067400931771D+00,
     *  0.6217724132262006D+00,0.6393199048503349D+00,
     1  0.6568404320717507D+00,0.6757355829009115D+00,
     2  0.6975593180107292D+00,0.7240884481675853D+00,
     3  0.7571353410751644D+00,0.7977350218356586D+00,
     4  0.8446899508348938D+00,0.8937570099213298D+00,
     5  0.9388189268759122D+00,0.9739910836494964D+00,
     *  0.9949470122490412D+00/

        data w2/
     1  0.2241914285263975D-05,0.5733435820783182D-04,
     2  0.5016551360939531D-03,0.2379144946058934D-02,
     3  0.7438688562803416D-02,0.1712636849320462D-01,
     4  0.3120190641013919D-01,0.4732382966527546D-01,
     5  0.6201072940466775D-01,0.7207354838184135D-01,
     6  0.7549892006385631D-01,0.7176320901318477D-01,
     7  0.6220952672284615D-01,0.5005773255222353D-01,
     8  0.3875653867514730D-01,0.3005057886427689D-01,
     9  0.2402528892050023D-01,0.2016959286765832D-01,
     *  0.1801743175835930D-01,0.1730654556796844D-01,
     1  0.1796569011512539D-01,0.2008203887175104D-01,
     2  0.2386365926269232D-01,0.2950951757944387D-01,
     3  0.3678445243563762D-01,0.4423397211610594D-01,
     4  0.4894533244911397D-01,0.4812879836649032D-01,
     5  0.4099267258419089D-01,0.2862536029514790D-01,
     *  0.1289769364570172D-01/

        data r3/
     1  0.6292327509937929D-06,0.2526829758524716D-04,
     2  0.3011837510675821D-03,0.1853583988140776D-02,
     3  0.7358811714943668D-02,0.2128549720326017D-01,
     4  0.4846775311760260D-01,0.9167919584028365D-01,
     5  0.1498718114331119D+00,0.2182566950010716D+00,
     6  0.2900661768600286D+00,0.3589420157836695D+00,
     7  0.4207604389591917D+00,0.4741289865256594D+00,
     8  0.5196925757871882D+00,0.5591077977431364D+00,
     9  0.5943533683870849D+00,0.6274488457172853D+00,
     *  0.6603786802367680D+00,0.6950320189339251D+00,
     1  0.7330159192411424D+00,0.7752291314574433D+00,
     2  0.8212231897552979D+00,0.8687227904867920D+00,
     3  0.9138753336417327D+00,0.9523137664818266D+00,
     4  0.9804405085027045D+00,0.9962978970151610D+00/

        data w3/
     1  0.2990146327356940D-05,0.7422632881927217D-04,
     2  0.6314749215044264D-03,0.2914246576428984D-02,
     3  0.8866298519573573D-02,0.1983798739002391D-01,
     4  0.3502132190220447D-01,0.5122487793421218D-01,
     5  0.6432997622109623D-01,0.7126904271477787D-01,
     6  0.7125771313422402D-01,0.6579874638719537D-01,
     7  0.5762016084293323D-01,0.4924326486386440D-01,
     8  0.4217146157541757D-01,0.3699285453822284D-01,
     9  0.3383570162994728D-01,0.3268667914572798D-01,
     *  0.3349065269407774D-01,0.3609091527236331D-01,
     1  0.4003822144110588D-01,0.4432319686924099D-01,
     2  0.4728188405471589D-01,0.4705565528903310D-01,
     3  0.4249864369406273D-01,0.3376309585740301D-01,
     4  0.2216352998627620D-01,0.9515180069220198D-02/

        data r4/
     1  0.4764656333176531D-06,0.1959443165043396D-04,
     2  0.2386633329755351D-03,0.1498086475970869D-02,
     3  0.6054641603788459D-02,0.1779614840116732D-01,
     4  0.4110743703910756D-01,0.7877459212303952D-01,
     5  0.1303950170267512D+00,0.1924490254747095D+00,
     6  0.2598813979892883D+00,0.3280430720634223D+00,
     7  0.3939035470239751D+00,0.4562574022762871D+00,
     8  0.5153256687246257D+00,0.5722118559192988D+00,
     9  0.6283865394850166D+00,0.6851557943577026D+00,
     *  0.7430598731287203D+00,0.8013183628370309D+00,
     1  0.8576246481829735D+00,0.9085294211559400D+00,
     2  0.9503108244614388D+00,0.9800004498364521D+00,
     *  0.9962759044580585D+00/

        data w4/
     1  0.2276177329771608D-05,0.5800260119538541D-04,
     2  0.5053267660238859D-03,0.2383479259363906D-02,
     3  0.7397480105693692D-02,0.1685548635906449D-01,
     4  0.3026115506993379D-01,0.4500431836611653D-01,
     5  0.5762031259808755D-01,0.6560985957095180D-01,
     6  0.6845927076854674D-01,0.6735682814928669D-01,
     7  0.6417109450492601D-01,0.6058287869285840D-01,
     8  0.5774596461944387D-01,0.5628015865689314D-01,
     9  0.5629716468676489D-01,0.5734147103063520D-01,
     *  0.5834077498506456D-01,0.5778507216836485D-01,
     1  0.5423712687186762D-01,0.4693447778812882D-01,
     2  0.3611764360201318D-01,0.2303539049003177D-01,
     *  0.9616986111413473D-02/

        data r5/
     1  0.4079438304167353D-06,0.1697174726338376D-04,
     2  0.2087972259269805D-03,0.1322141803859463D-02,
     3  0.5386188427857177D-02,0.1595807120778799D-01,
     4  0.3720304474662404D-01,0.7215181138133588D-01,
     5  0.1213950696046413D+00,0.1831157825433372D+00,
     6  0.2542135914352367D+00,0.3316032003670487D+00,
     7  0.4129776743523134D+00,0.4968524970674780D+00,
     8  0.5820516419191586D+00,0.6669748269789489D+00,
     9  0.7490799333889332D+00,0.8248953220737811D+00,
     *  0.8905220145340645D+00,0.9423323674020277D+00,
     1  0.9777034442271580D+00,0.9960394391011756D+00/

        data w5/
     1  0.1953948391311679D-05,0.5042101226490901D-04,
     2  0.4440388002821658D-03,0.2114455081565467D-02,
     3  0.6622707923308921D-02,0.1524798524687594D-01,
     4  0.2777135959635366D-01,0.4222322030200929D-01,
     5  0.5593908791746347D-01,0.6696217969766280D-01,
     6  0.7470618606114605D-01,0.7968825927176325D-01,
     7  0.8282264719882537D-01,0.8474433235917779D-01,
     8  0.8539675591842579D-01,0.8402557844881506D-01,
     9  0.7959677005539454D-01,0.7137201301019582D-01,
     *  0.5926520863004840D-01,0.4390262429271031D-01,
     1  0.2671331444638880D-01,0.1038890078093090D-01/

        data r6/
     1  0.3549045593306105D-06,0.1493905453233523D-04,
     2  0.1860078983917652D-03,0.1193307309943876D-02,
     3  0.4934612660092158D-02,0.1488532083641340D-01,
     4  0.3547885530937722D-01,0.7069012599350122D-01,
     5  0.1227467272912404D+00,0.1916581229310709D+00,
     6  0.2754817956486991D+00,0.3708312807723686D+00,
     7  0.4733143577534458D+00,0.5779275178303295D+00,
     8  0.6794951302370238D+00,0.7731044910511489D+00,
     9  0.8544355459620947D+00,0.9199503322358200D+00,
     *  0.9669693566170180D+00,0.9936831454494705D+00/

        data w6/
     1  0.1704122866194641D-05,0.4455711628617572D-04,
     2  0.3978801866357356D-03,0.1924931537763917D-02,
     3  0.6147734747672850D-02,0.1451892297175960D-01,
     4  0.2735399759756744D-01,0.4344725409054134D-01,
     5  0.6066161418361882D-01,0.7682815691369450D-01,
     6  0.9024769930068867D-01,0.9971296730236038D-01,
     7  0.1044117672634205D+00,0.1039443971293551D+00,
     8  0.9836647840637292D-01,0.8812916503350978D-01,
     9  0.7394211979825537D-01,0.5664085082363000D-01,
     *  0.3709540855035259D-01,0.1618239292364807D-01/
c        
c         %%%%%%%%%%%%%%%%%%%%%

ccc        call prin2('alpha = *',alpha,1)
        dtwo = 2
        if (alpha .ge. 1/dtwo) then
        do 1000 i=1,nqs(6)
        roots(i) =r6(i)
        weights(i) = w6(i)
 1000 continue

        nquads = nqs(6)
        endif

        if (alpha .ge. dtwo**(-2)) then
        if (alpha .lt. dtwo**(-1)) then

        do 2000 i=1,nqs(5)
        roots(i) =r5(i)
        weights(i) = w5(i)
 2000 continue

        nquads = nqs(5)

        endif
        endif


        if (alpha .ge. dtwo**(-3)) then
        if (alpha .lt. dtwo**(-2)) then

        do 3000 i=1,nqs(4)
        roots(i) =r4(i)
        weights(i) = w4(i)
 3000 continue

        nquads = nqs(4)

        endif
        endif


        if (alpha .ge. dtwo**(-4)) then
        if (alpha .lt. dtwo**(-3)) then

        do 4000 i=1,nqs(3)
        roots(i) =r3(i)
        weights(i) = w3(i)
 4000 continue

        nquads = nqs(3)

        endif
        endif


        if (alpha .ge. dtwo**(-5)) then
        if (alpha .lt. dtwo**(-4)) then

        do 5000 i=1,nqs(2)
        roots(i) =r2(i)
        weights(i) = w2(i)
 5000 continue

        nquads = nqs(2)

        endif
        endif


        if (alpha .lt. dtwo**(-5)) then

        do 6000 i=1,nqs(1)
        roots(i) =r1(i)
        weights(i) = w1(i)
 6000 continue

        nquads = nqs(1)

        endif

        return
        end
c
c
c
c
c
        subroutine lap_load_quads27(alpha,roots,weights,nquads)
        implicit real *8 (a-h,o-z)
        dimension r1(32),r2(30),r3(27),r4(24),r5(21),r6(20),
     1          w1(32),w2(30),w3(27),w4(24),w5(21),w6(20),nqs(6),
     2          roots(1),weights(1)

        data nqs/32,30,27,24,21,20/


        data r1/
     1  0.6459130800329935D-06,0.2621179948352761D-04,
     2  0.3156716266593581D-03,0.1962716286712898D-02,
     3  0.7874064460641530D-02,0.2303769587749971D-01,
     4  0.5316090940526970D-01,0.1021782350732437D+00,
     5  0.1702278330153507D+00,0.2531844873522259D+00,
     6  0.3435408074361551D+00,0.4318306131860660D+00,
     7  0.5085313233237165D+00,0.5674004947404752D+00,
     8  0.6085510712261150D+00,0.6366025824953317D+00,
     9  0.6563834531581946D+00,0.6712102829420803D+00,
     *  0.6831381665054657D+00,0.6935066029416544D+00,
     1  0.7033164284167182D+00,0.7134639078683989D+00,
     2  0.7249191421660642D+00,0.7389214839126161D+00,
     3  0.7572309260119017D+00,0.7823075514259365D+00,
     4  0.8166941569167703D+00,0.8603987110258006D+00,
     5  0.9081252278407979D+00,0.9512638152609573D+00,
     6  0.9823434887188413D+00,0.9974638884067871D+00/

        data w1/
     1  0.3076273563009948D-05,0.7727050591254972D-04,
     2  0.6650864584293809D-03,0.3105958554645626D-02,
     3  0.9572156185438551D-02,0.2175251823849768D-01,
     4  0.3918633016442055D-01,0.5887826837045103D-01,
     5  0.7651618209041591D-01,0.8812159077850582D-01,
     6  0.9098546039322832D-01,0.8395005160186381D-01,
     7  0.6833972252939858D-01,0.4945864499709240D-01,
     8  0.3370087772677013D-01,0.2322989426497544D-01,
     9  0.1687618226847789D-01,0.1311025963794310D-01,
     *  0.1096266549588645D-01,0.9937330562849895D-02,
     1  0.9828896943403342D-02,0.1062461530341885D-01,
     2  0.1248809019486845D-01,0.1580549254056505D-01,
     3  0.2122937875279552D-01,0.2939319664817698D-01,
     4  0.3944815389835349D-01,0.4703508407096438D-01,
     5  0.4689243082673177D-01,0.3811047549278377D-01,
     6  0.2335147434185086D-01,0.7363183887321416D-02/

        data r2/
     1  0.6075313346990485D-06,0.2471297069132322D-04,
     2  0.2982158447285461D-03,0.1857189326557453D-02,
     3  0.7458808969098785D-02,0.2182843558616302D-01,
     4  0.5032318491824908D-01,0.9648164830684940D-01,
     5  0.1600423792778514D+00,0.2365867629523543D+00,
     6  0.3187447866090096D+00,0.3982599670325963D+00,
     7  0.4683468586786426D+00,0.5256526278341876D+00,
     8  0.5705593409693365D+00,0.6055786405134335D+00,
     9  0.6335295345179819D+00,0.6567441855180543D+00,
     *  0.6770510517362138D+00,0.6959764599480605D+00,
     1  0.7149299898869658D+00,0.7353402939275740D+00,
     2  0.7587354528903625D+00,0.7866929468628891D+00,
     3  0.8204423616845248D+00,0.8598805870308147D+00,
     4  0.9024003736365376D+00,0.9428980318557213D+00,
     5  0.9754492123176660D+00,0.9952006019814985D+00/

        data w2/
     1  0.2895066706269639D-05,0.7290574161893520D-04,
     2  0.6288521480681053D-03,0.2941549277408733D-02,
     3  0.9072803332779206D-02,0.2060472861530979D-01,
     4  0.3701066652936507D-01,0.5526815080042855D-01,
     5  0.7109467680305966D-01,0.8073011979768921D-01,
     6  0.8216943116413506D-01,0.7569344132048095D-01,
     7  0.6391502888107750D-01,0.5079268265562005D-01,
     8  0.3947304304990426D-01,0.3104723430929444D-01,
     9  0.2524200923463414D-01,0.2148822864743267D-01,
     *  0.1937678918970076D-01,0.1870670982545575D-01,
     1  0.1943697817305571D-01,0.2163862484823974D-01,
     2  0.2542037165307315D-01,0.3071649090992916D-01,
     3  0.3678893043869124D-01,0.4165493031888213D-01,
     4  0.4250097281708548D-01,0.3747245102142377D-01,
     5  0.2680562565084650D-01,0.1223267777860401D-01/

        data r3/
     1  0.6493016629064520D-06,0.2611314539944995D-04,
     2  0.3116799223508985D-03,0.1920410978668235D-02,
     3  0.7631438034450950D-02,0.2209345618146437D-01,
     4  0.5035730221820503D-01,0.9537708928794234D-01,
     5  0.1561918222051925D+00,0.2279904585092754D+00,
     6  0.3038897415063325D+00,0.3773373450817807D+00,
     7  0.4439440321846954D+00,0.5020239174580690D+00,
     8  0.5519824569007307D+00,0.5953315198178416D+00,
     9  0.6339686131228795D+00,0.6698474101865274D+00,
     *  0.7048749426730500D+00,0.7408291608718556D+00,
     1  0.7791395859230250D+00,0.8204474186681401D+00,
     2  0.8640190875708082D+00,0.9073697305478025D+00,
     3  0.9465620165556506D+00,0.9771977244711032D+00,
     *  0.9955572639844275D+00/

        data w3/
     1  0.3086504294508696D-05,0.7674586786529753D-04,
     2  0.6538786836650260D-03,0.3021330099211271D-02,
     3  0.9201566681891352D-02,0.2061085973956261D-01,
     4  0.3644273065871162D-01,0.5343489433536855D-01,
     5  0.6736014697175599D-01,0.7504681875698437D-01,
     6  0.7562723423192536D-01,0.7052782954777632D-01,
     7  0.6241705483380730D-01,0.5382836857235312D-01,
     8  0.4635285420486710D-01,0.4066651820260579D-01,
     9  0.3693447569657861D-01,0.3514195729973813D-01,
     *  0.3521286061363285D-01,0.3694125607201991D-01,
     1  0.3979712249456613D-01,0.4270397667880437D-01,
     2  0.4402037455122556D-01,0.4200648321390212D-01,
     3  0.3562477006969617D-01,0.2501312851365407D-01,
     *  0.1133167690353651D-01/

        data r4/
     1  0.5202052793899751D-06,0.2128174366186101D-04,
     2  0.2578838538911325D-03,0.1610640325577372D-02,
     3  0.6479162898664627D-02,0.1896725910472566D-01,
     4  0.4368058634214315D-01,0.8356632851878043D-01,
     5  0.1383068080394429D+00,0.2043839748602038D+00,
     6  0.2766208794281202D+00,0.3500953049403655D+00,
     7  0.4214060014418314D+00,0.4889876573872297D+00,
     8  0.5527809416396114D+00,0.6136728406727907D+00,
     9  0.6729071334847276D+00,0.7314762272364406D+00,
     *  0.7895089669591294D+00,0.8458216408724084D+00,
     1  0.8978930527959246D+00,0.9423532851097943D+00,
     2  0.9757585268106240D+00,0.9953140146264035D+00/

        data w4/
     1  0.2482377158434247D-05,0.6288703461547606D-04,
     2  0.5447066007946567D-03,0.2554682320921461D-02,
     3  0.7887550879194436D-02,0.1789661949925194D-01,
     4  0.3205108495476348D-01,0.4766567988669566D-01,
     5  0.6120014134716094D-01,0.7005676713849269D-01,
     6  0.7357491862431867D-01,0.7280279105601506D-01,
     7  0.6956434270704529D-01,0.6559885284000705D-01,
     8  0.6214530921586994D-01,0.5985588362236207D-01,
     9  0.5878872228931906D-01,0.5837917996832026D-01,
     *  0.5748972918719079D-01,0.5470661778680038D-01,
     1  0.4886068772334097D-01,0.3947338026025366D-01,
     2  0.2686289467493303D-01,0.1197408800517460D-01/

        data r5/
     1  0.4720330832669798D-06,0.1934112079169802D-04,
     2  0.2349118707468651D-03,0.1472180411587357D-02,
     3  0.5950500067488647D-02,0.1753170948005293D-01,
     4  0.4071756005288150D-01,0.7877002214967363D-01,
     5  0.1322952961595731D+00,0.1992588333802553D+00,
     6  0.2761829023763147D+00,0.3595412188738325D+00,
     7  0.4465725792468900D+00,0.5353131337719505D+00,
     8  0.6240553002602684D+00,0.7106499244649903D+00,
     9  0.7921122547193990D+00,0.8647647848965330D+00,
     *  0.9247660365003220D+00,0.9686948624818672D+00,
     *  0.9939844351980732D+00/

        data w5/
     1  0.2253127977223292D-05,0.5718874989985680D-04,
     2  0.4968465629053955D-03,0.2341027871569283D-02,
     3  0.7276929964894352D-02,0.1667100375705265D-01,
     4  0.3027066847689449D-01,0.4593416469455146D-01,
     5  0.6075451173559494D-01,0.7256771629483709D-01,
     6  0.8067879218869153D-01,0.8557767440301264D-01,
     7  0.8816778974447842D-01,0.8904158074425876D-01,
     8  0.8809177612253275D-01,0.8459716077790999D-01,
     9  0.7770483046351161D-01,0.6694973162546830D-01,
     *  0.5247541684722421D-01,0.3494981731333297D-01,
     *  0.1539311853340208D-01/

        data r6/
     1  0.3020826322838024D-06,0.1318243030243019D-04,
     2  0.1693213829552858D-03,0.1116097963416516D-02,
     3  0.4724290939147464D-02,0.1453321734745424D-01,
     4  0.3519712736629325D-01,0.7100859499219201D-01,
     5  0.1244409848421927D+00,0.1955256732267791D+00,
     6  0.2820802472395234D+00,0.3802956193427083D+00,
     7  0.4853053584681703D+00,0.5916850554660790D+00,
     8  0.6939326233994194D+00,0.7869056392133113D+00,
     9  0.8661745082477063D+00,0.9283220363640400D+00,
     *  0.9712772265749130D+00,0.9946375438244889D+00/

        data w6/
     1  0.1462357921042546D-05,0.3975692241131231D-04,
     2  0.3671024831309960D-03,0.1828512764740910D-02,
     3  0.5987139233340035D-02,0.1443381699235800D-01,
     4  0.2763598393261545D-01,0.4440941193249686D-01,
     5  0.6245530800169025D-01,0.7933957128857146D-01,
     6  0.9312525405487136D-01,0.1024863022868997D+00,
     7  0.1066209965642530D+00,0.1052139521393090D+00,
     8  0.9841830864680857D-01,0.8678472152714197D-01,
     9  0.7117805071592238D-01,0.5276502043343330D-01,
     *  0.3308062190889479D-01,0.1382870581318955D-01/
c        
c         %%%%%%%%%%%%%%%%%%%%%

ccc        call prin2('alpha = *',alpha,1)
        dtwo = 2
        if (alpha .ge. 1/dtwo) then
        do 1000 i=1,nqs(6)
        roots(i) =r6(i)
        weights(i) = w6(i)
 1000 continue

        nquads = nqs(6)
        endif

        if (alpha .ge. dtwo**(-2)) then
        if (alpha .lt. dtwo**(-1)) then

        do 2000 i=1,nqs(5)
        roots(i) =r5(i)
        weights(i) = w5(i)
 2000 continue

        nquads = nqs(5)

        endif
        endif


        if (alpha .ge. dtwo**(-3)) then
        if (alpha .lt. dtwo**(-2)) then

        do 3000 i=1,nqs(4)
        roots(i) =r4(i)
        weights(i) = w4(i)
 3000 continue

        nquads = nqs(4)

        endif
        endif


        if (alpha .ge. dtwo**(-4)) then
        if (alpha .lt. dtwo**(-3)) then

        do 4000 i=1,nqs(3)
        roots(i) =r3(i)
        weights(i) = w3(i)
 4000 continue

        nquads = nqs(3)

        endif
        endif


        if (alpha .ge. dtwo**(-5)) then
        if (alpha .lt. dtwo**(-4)) then

        do 5000 i=1,nqs(2)
        roots(i) =r2(i)
        weights(i) = w2(i)
 5000 continue

        nquads = nqs(2)

        endif
        endif


        if (alpha .lt. dtwo**(-5)) then

        do 6000 i=1,nqs(1)
        roots(i) =r1(i)
        weights(i) = w1(i)
 6000 continue

        nquads = nqs(1)

        endif

        return
        end
c
c
c
c
c
        subroutine lap_load_quads28(alpha,roots,weights,nquads)
        implicit real *8 (a-h,o-z)
        dimension r1(32),r2(30),r3(27),r4(23),r5(21),r6(19),
     1          w1(32),w2(30),w3(27),w4(23),w5(21),w6(19),nqs(6),
     2          roots(1),weights(1)

        data nqs/32,30,27,23,21,19/


        data r1/
     1  0.2911323489919838D-06,0.1325677039515448D-04,
     2  0.1762895548125244D-03,0.1196154477653955D-02,
     3  0.5184427179275974D-02,0.1623880251298205D-01,
     4  0.3978197309584248D-01,0.8056487324844178D-01,
     5  0.1404948992788992D+00,0.2175868789763749D+00,
     6  0.3063516216342664D+00,0.3990032780948535D+00,
     7  0.4868161213317371D+00,0.5617771468903753D+00,
     8  0.6192145028687571D+00,0.6598375705541654D+00,
     9  0.6879716324482770D+00,0.7080821581827965D+00,
     *  0.7233141297410604D+00,0.7356615756678612D+00,
     1  0.7464404803583183D+00,0.7566346271727879D+00,
     2  0.7671145914927965D+00,0.7787994081113557D+00,
     3  0.7928170486058293D+00,0.8106758604374753D+00,
     4  0.8342945694754502D+00,0.8652917366181261D+00,
     5  0.9028359038871204D+00,0.9418059771709013D+00,
     6  0.9746587501135430D+00,0.9950181803787793D+00/

        data w1/
     1  0.1423309331866265D-05,0.4048582202105978D-04,
     2  0.3877897283446128D-03,0.1991139784686856D-02,
     3  0.6678950553350068D-02,0.1637022869546291D-01,
     4  0.3154776846313090D-01,0.5035409487614230D-01,
     5  0.6916672447453218D-01,0.8408026678595056D-01,
     6  0.9212921150301656D-01,0.9169524451096269D-01,
     7  0.8254859975430881D-01,0.6654713940987535D-01,
     8  0.4848356894189465D-01,0.3355361412528167D-01,
     9  0.2347927986986994D-01,0.1725959204215523D-01,
     *  0.1352813040881519D-01,0.1138044766489244D-01,
     1  0.1033785789144725D-01,0.1019297521535156D-01,
     2  0.1091751588325405D-01,0.1263702572250853D-01,
     3  0.1564782011083378D-01,0.2039853268153906D-01,
     4  0.2714754285449896D-01,0.3475299351487149D-01,
     5  0.3944259144154729D-01,0.3716847143756919D-01,
     6  0.2745055188282788D-01,0.1268242063972512D-01/

        data r2/
     1  0.6703879729885028D-06,0.2699551606654534D-04,
     2  0.3226695237450068D-03,0.1991183694452026D-02,
     3  0.7927427181191770D-02,0.2301139657398558D-01,
     4  0.5266388080711246D-01,0.1003555668214737D+00,
     5  0.1657375686118603D+00,0.2444908361264369D+00,
     6  0.3296629361730160D+00,0.4135788255838702D+00,
     7  0.4896571656060988D+00,0.5538595052436261D+00,
     8  0.6053255709864902D+00,0.6456936881467893D+00,
     9  0.6776344345177721D+00,0.7037470337941660D+00,
     *  0.7261939948382177D+00,0.7467552614566225D+00,
     1  0.7669841159057042D+00,0.7883289364141465D+00,
     2  0.8121599804627099D+00,0.8396075625809812D+00,
     3  0.8710697910882472D+00,0.9054346166294106D+00,
     4  0.9396471645618006D+00,0.9693889038195053D+00,
     5  0.9905011031895989D+00,0.9999999999999879D+00/

        data w2/
     1  0.3187513720976413D-05,0.7937598743440715D-04,
     2  0.6774074829878656D-03,0.3135945465279201D-02,
     3  0.9575593317862618D-02,0.2154216900729890D-01,
     4  0.3837534590801358D-01,0.5695146572576183D-01,
     5  0.7306778914926113D-01,0.8324650162397377D-01,
     6  0.8578795954743844D-01,0.8092177710903951D-01,
     7  0.7054926216899135D-01,0.5772484957589055D-01,
     8  0.4550840400202527D-01,0.3569199053007858D-01,
     9  0.2862767828876774D-01,0.2395625430731269D-01,
     *  0.2123146357584042D-01,0.2014782001554650D-01,
     1  0.2055063494590140D-01,0.2236966409897384D-01,
     2  0.2548886852549270D-01,0.2948720631203972D-01,
     3  0.3326180459424083D-01,0.3494626810944431D-01,
     4  0.3273181066918045D-01,0.2604005313033273D-01,
     5  0.1569098659475017D-01,0.2630462717118514D-02/

        data r3/
     1  0.6098385080859034D-06,0.2470735710528536D-04,
     2  0.2969325437689388D-03,0.1841486472033070D-02,
     3  0.7364120026779130D-02,0.2145612248212456D-01,
     4  0.4923896777665781D-01,0.9397076260795838D-01,
     5  0.1552312900177547D+00,0.2288363514872178D+00,
     6  0.3083311701775684D+00,0.3870942805149765D+00,
     7  0.4600899426996261D+00,0.5247013338569561D+00,
     8  0.5805588686374040D+00,0.6287674159614148D+00,
     9  0.6711183241572251D+00,0.7095948506421588D+00,
     *  0.7461242808515587D+00,0.7823987691532795D+00,
     1  0.8196153365375621D+00,0.8580764521509098D+00,
     2  0.8967689905062664D+00,0.9332519302272081D+00,
     3  0.9641504102839573D+00,0.9861884330709820D+00,
     *  0.9976246919520173D+00/

        data w3/
     1  0.2903610835896505D-05,0.7278932992464587D-04,
     2  0.6249192774249149D-03,0.2908893308767637D-02,
     3  0.8925593973467659D-02,0.2015655745239422D-01,
     4  0.3598550833117863D-01,0.5340881778422850D-01,
     5  0.6838728018486444D-01,0.7770956514216872D-01,
     6  0.8014593229934587D-01,0.7652686303699504D-01,
     7  0.6903227981864951D-01,0.6014342901807084D-01,
     8  0.5176966433671333D-01,0.4495399295354430D-01,
     9  0.4008201466794707D-01,0.3719487157145705D-01,
     *  0.3615017054993896D-01,0.3660520508820299D-01,
     1  0.3788624746373016D-01,0.3886992310098357D-01,
     2  0.3809981113139244D-01,0.3428558624502726D-01,
     3  0.2694283772130878D-01,0.1680691521817052D-01,
     *  0.6321427383267046D-02/

        data r4/
     1  0.6121234332699481D-06,0.2472898629895850D-04,
     2  0.2964772454685520D-03,0.1834350809238846D-02,
     3  0.7317115131445419D-02,0.2125898330775940D-01,
     4  0.4863254057044439D-01,0.9250569029281107D-01,
     5  0.1523624395378880D+00,0.2242603863386674D+00,
     6  0.3025438425468967D+00,0.3819001469141181D+00,
     7  0.4586842314925754D+00,0.5312156272192024D+00,
     8  0.5993873182852983D+00,0.6640142379588831D+00,
     9  0.7261201666310330D+00,0.7862118720417174D+00,
     *  0.8436587396962658D+00,0.8964416721214664D+00,
     1  0.9414734261839430D+00,0.9753608554172549D+00,
     *  0.9952335697367850D+00/

        data w4/
     1  0.2912482698206624D-05,0.7278988105798377D-04,
     2  0.6232706340450055D-03,0.2893209572919449D-02,
     3  0.8849167335240291D-02,0.1990754158109471D-01,
     4  0.3538203316872504D-01,0.5227468298116547D-01,
     5  0.6674980696306093D-01,0.7606635766168909D-01,
     6  0.7959242229995796D-01,0.7850881160630961D-01,
     7  0.7478695127884271D-01,0.7026927780324615D-01,
     8  0.6622182096679406D-01,0.6321498238816783D-01,
     9  0.6109240104847119D-01,0.5899072888629224D-01,
     *  0.5556080020711777D-01,0.4947673108124505D-01,
     1  0.4000675588509426D-01,0.2727916016404125D-01,
     *  0.1217738412272375D-01/

        data r5/
     1  0.2813256779251761D-06,0.1228752112781696D-04,
     2  0.1586879092453792D-03,0.1055105955551527D-02,
     3  0.4510301798534743D-02,0.1399375442233932D-01,
     4  0.3405171422332460D-01,0.6863513483778854D-01,
     5  0.1194053076931273D+00,0.1851866378671115D+00,
     6  0.2627610778156185D+00,0.3482476642755748D+00,
     7  0.4381781110107321D+00,0.5298464305590316D+00,
     8  0.6210038500605083D+00,0.7092644922950111D+00,
     9  0.7916779755872670D+00,0.8647525132407212D+00,
     *  0.9248622228290253D+00,0.9687616734267759D+00,
     *  0.9939999288310856D+00/

        data w5/
     1  0.1361398515000742D-05,0.3709889096358788D-04,
     2  0.3452121088173244D-03,0.1738765198978101D-02,
     3  0.5759330848037416D-02,0.1399621567720763D-01,
     4  0.2680975840906854D-01,0.4264349878306337D-01,
     5  0.5868544819130150D-01,0.7231550018131114D-01,
     6  0.8216765407469042D-01,0.8822191015955706D-01,
     7  0.9119400472547766D-01,0.9177894742204895D-01,
     8  0.9014614461820707D-01,0.8588279063795514D-01,
     9  0.7835462148699762D-01,0.6718102746602344D-01,
     *  0.5249058051972637D-01,0.3489510550971759D-01,
     *  0.1535502369233507D-01/

        data r6/
     1  0.5261393246493632D-06,0.2144186850076424D-04,
     2  0.2592648535266039D-03,0.1619080711333190D-02,
     3  0.6531945168839304D-02,0.1926452416301778D-01,
     4  0.4498288000763545D-01,0.8795238057734012D-01,
     5  0.1500528924439346D+00,0.2303572135379158D+00,
     6  0.3255972747633036D+00,0.4309401857152445D+00,
     7  0.5407088806446136D+00,0.6489794233339336D+00,
     8  0.7500608645722172D+00,0.8388414318506019D+00,
     9  0.9110064313422142D+00,0.9631665169832244D+00,
     *  0.9929433775852502D+00/

        data w6/
     1  0.2508211215185674D-05,0.6329562748763872D-04,
     2  0.5473017365202015D-03,0.2569968579277194D-02,
     3  0.7983579291105229D-02,0.1837605619736942D-01,
     4  0.3379828041352072D-01,0.5246033959166871D-01,
     5  0.7157670738789964D-01,0.8847344816953731D-01,
     6  0.1011891534906778D+00,0.1085405806411206D+00,
     7  0.1100014139094821D+00,0.1055854215977028D+00,
     8  0.9572197664868545D-01,0.8111801039879141D-01,
     9  0.6264581568405056D-01,0.4127627466149145D-01,
     *  0.1806986776239656D-01/
c        
c         %%%%%%%%%%%%%%%%%%%%%

ccc        call prin2('alpha = *',alpha,1)
        dtwo = 2
        if (alpha .ge. 1/dtwo) then
        do 1000 i=1,nqs(6)
        roots(i) =r6(i)
        weights(i) = w6(i)
 1000 continue

        nquads = nqs(6)
        endif

        if (alpha .ge. dtwo**(-2)) then
        if (alpha .lt. dtwo**(-1)) then

        do 2000 i=1,nqs(5)
        roots(i) =r5(i)
        weights(i) = w5(i)
 2000 continue

        nquads = nqs(5)

        endif
        endif


        if (alpha .ge. dtwo**(-3)) then
        if (alpha .lt. dtwo**(-2)) then

        do 3000 i=1,nqs(4)
        roots(i) =r4(i)
        weights(i) = w4(i)
 3000 continue

        nquads = nqs(4)

        endif
        endif


        if (alpha .ge. dtwo**(-4)) then
        if (alpha .lt. dtwo**(-3)) then

        do 4000 i=1,nqs(3)
        roots(i) =r3(i)
        weights(i) = w3(i)
 4000 continue

        nquads = nqs(3)

        endif
        endif


        if (alpha .ge. dtwo**(-5)) then
        if (alpha .lt. dtwo**(-4)) then

        do 5000 i=1,nqs(2)
        roots(i) =r2(i)
        weights(i) = w2(i)
 5000 continue

        nquads = nqs(2)

        endif
        endif


        if (alpha .lt. dtwo**(-5)) then

        do 6000 i=1,nqs(1)
        roots(i) =r1(i)
        weights(i) = w1(i)
 6000 continue

        nquads = nqs(1)

        endif

        return
        end
c
c
c
c
c
        subroutine lap_load_quads29(alpha,roots,weights,nquads)
        implicit real *8 (a-h,o-z)
        dimension r1(31),r2(29),r3(26),r4(23),r5(20),r6(19),
     1          w1(31),w2(29),w3(26),w4(23),w5(20),w6(19),nqs(6),
     2          roots(1),weights(1)

        data nqs/31,29,26,23,20,19/


        data r1/
     1  0.5770469552093510D-06,0.2358330345434847D-04,
     2  0.2858837058570146D-03,0.1788565918056037D-02,
     3  0.7218330392507818D-02,0.2124386055982261D-01,
     4  0.4931967586548618D-01,0.9542709652575109D-01,
     5  0.1602356011594014D+00,0.2407376338312094D+00,
     6  0.3311872752120128D+00,0.4244566649631388D+00,
     7  0.5131430895811586D+00,0.5905319912180440D+00,
     8  0.6521069259794256D+00,0.6973154148069430D+00,
     9  0.7292727592834087D+00,0.7521761919399853D+00,
     *  0.7694024848111895D+00,0.7832314965627580D+00,
     1  0.7952021228583537D+00,0.8064671438669625D+00,
     2  0.8180362849553783D+00,0.8309481081670089D+00,
     3  0.8464045923787753D+00,0.8658096198371320D+00,
     4  0.8904085309070769D+00,0.9200176219866295D+00,
     5  0.9513737800662254D+00,0.9785561946598900D+00,
     *  0.9957532065482812D+00/

        data w1/
     1  0.2752580252951610D-05,0.6968156458325384D-04,
     2  0.6041445419240157D-03,0.2841212255364361D-02,
     3  0.8817543532755755D-02,0.2018457296360491D-01,
     4  0.3666411431317095D-01,0.5566769465255787D-01,
     5  0.7342986491533827D-01,0.8659212947816637D-01,
     6  0.9310098984442749D-01,0.9219018976443402D-01,
     7  0.8403895555911598D-01,0.6994864170558210D-01,
     8  0.5311219032981796D-01,0.3787839736358457D-01,
     9  0.2676866339445701D-01,0.1960212840069518D-01,
     *  0.1522316476609640D-01,0.1268532926264618D-01,
     1  0.1144414357432658D-01,0.1125142470899219D-01,
     2  0.1205763116545885D-01,0.1396506511636036D-01,
     3  0.1718623673196312D-01,0.2185262756028182D-01,
     4  0.2735041951715359D-01,0.3131884442094670D-01,
     5  0.3034933647772650D-01,0.2300955887451078D-01,
     *  0.1079235066370393D-01/

        data r2/
     1  0.6047218535290033D-06,0.2462771702422323D-04,
     2  0.2975090616888929D-03,0.1854696475579614D-02,
     3  0.7457099597339026D-02,0.2185517763659671D-01,
     4  0.5049344156735728D-01,0.9712853357661793D-01,
     5  0.1619250445831973D+00,0.2411462975948037D+00,
     6  0.3283096838071551D+00,0.4159572693451606D+00,
     7  0.4973523233089663D+00,0.5678431081745039D+00,
     8  0.6256451615586471D+00,0.6716013469176195D+00,
     9  0.7080227912306802D+00,0.7375005066444126D+00,
     *  0.7623132816867188D+00,0.7843383681454051D+00,
     1  0.8051551122973780D+00,0.8261538269521944D+00,
     2  0.8485670591877207D+00,0.8733509000552039D+00,
     3  0.9008195907702243D+00,0.9300361961540239D+00,
     4  0.9583552577849099D+00,0.9818096084442382D+00,
     *  0.9964082823818435D+00/

        data w2/
     1  0.2882423831391780D-05,0.7268203408735841D-04,
     2  0.6276720059173226D-03,0.2939528045498144D-02,
     3  0.9079909441515947D-02,0.2066790433454603D-01,
     4  0.3726948082647109D-01,0.5603266563464079D-01,
     5  0.7291587302474100D-01,0.8441493158822624D-01,
     6  0.8864019929622387D-01,0.8550905569884561D-01,
     7  0.7649122475799155D-01,0.6419694106811748D-01,
     8  0.5156977870908130D-01,0.4074924400705433D-01,
     9  0.3253479525650156D-01,0.2680106189783569D-01,
     *  0.2313522302723513D-01,0.2117620552042919D-01,
     1  0.2068832823198101D-01,0.2151623569936594D-01,
     2  0.2347551847923068D-01,0.2615558107969089D-01,
     3  0.2863544180395439D-01,0.2935080630205968D-01,
     4  0.2660182041103856D-01,0.1961683878188483D-01,
     *  0.9132170612003001D-02/

        data r3/
     1  0.6050557831616131D-06,0.2458871053490731D-04,
     2  0.2963553690979104D-03,0.1842866492808625D-02,
     3  0.7388511561005695D-02,0.2158120159814583D-01,
     4  0.4965342727470594D-01,0.9502227920826603D-01,
     5  0.1574419188857247D+00,0.2328733216001854D+00,
     6  0.3149397862379358D+00,0.3970027931067651D+00,
     7  0.4738884133603027D+00,0.5427165930220666D+00,
     8  0.6027870881797505D+00,0.6549010998600336D+00,
     9  0.7006248706093851D+00,0.7417779690269049D+00,
     *  0.7801387429639980D+00,0.8172246399346017D+00,
     1  0.8540138469069435D+00,0.8905753209115657D+00,
     2  0.9257462451382724D+00,0.9571394650537937D+00,
     3  0.9816672350168018D+00,0.9964206780869893D+00/

        data w3/
     1  0.2882761836311386D-05,0.7251252180655539D-04,
     2  0.6245281766123087D-03,0.2915894240175417D-02,
     3  0.8973731169662463D-02,0.2032839588404189D-01,
     4  0.3641950975122035D-01,0.5427879242612721D-01,
     5  0.6986036566662798D-01,0.7990260015013957D-01,
     6  0.8309382323035610D-01,0.8015594771381397D-01,
     7  0.7313821921025641D-01,0.6441412369954045D-01,
     8  0.5587336005063914D-01,0.4862412885017251D-01,
     9  0.4313118894526452D-01,0.3947471104998893D-01,
     *  0.3750427144516560D-01,0.3683704679755893D-01,
     1  0.3675901635932357D-01,0.3617059614943282D-01,
     2  0.3376421296486749D-01,0.2849185890614867D-01,
     3  0.2006386521803535D-01,0.9124416661185490D-02/

        data r4/
     1  0.2872555248680219D-06,0.1295203149971113D-04,
     2  0.1709251403779092D-03,0.1152045377225311D-02,
     3  0.4960968734352037D-02,0.1543129754146457D-01,
     4  0.3750455137092974D-01,0.7525396785030660D-01,
     5  0.1298793151033161D+00,0.1990312329342470D+00,
     6  0.2777829297559638D+00,0.3604875552444569D+00,
     7  0.4424518201979423D+00,0.5207730436806290D+00,
     8  0.5943231513692024D+00,0.6632276843809088D+00,
     9  0.7281568291122888D+00,0.7895812228936307D+00,
     *  0.8471173412107315D+00,0.8991849950485430D+00,
     1  0.9431820001183343D+00,0.9761170594547215D+00,
     *  0.9953831623766028D+00/

        data w4/
     1  0.1401018304215768D-05,0.3944535735130373D-04,
     2  0.3747918817187682D-03,0.1910269820816997D-02,
     3  0.6358254299289874D-02,0.1544244756229603D-01,
     4  0.2942015801342017D-01,0.4629062151631903D-01,
     5  0.6255012177247660D-01,0.7489876262425721D-01,
     6  0.8163395472398552D-01,0.8297923065259242D-01,
     7  0.8046680828037478D-01,0.7599711408949056D-01,
     8  0.7113612521078739D-01,0.6679726824548908D-01,
     9  0.6315095244473009D-01,0.5964025043709303D-01,
     *  0.5516557423829811D-01,0.4852456633337568D-01,
     1  0.3895726356497034D-01,0.2646757380169017D-01,
     *  0.1179704411087263D-01/

        data r5/
     1  0.5770605632378326D-06,0.2352013934324844D-04,
     2  0.2841714001089884D-03,0.1770716124602726D-02,
     3  0.7112239272168689D-02,0.2081613598595316D-01,
     4  0.4803180243309093D-01,0.9236306719285350D-01,
     5  0.1542939213727935D+00,0.2312499034944896D+00,
     6  0.3189475085648984D+00,0.4129269076099829D+00,
     7  0.5094196908002150D+00,0.6053695632423746D+00,
     8  0.6979309911440548D+00,0.7839703077918482D+00,
     9  0.8599597027385074D+00,0.9222846883761551D+00,
     *  0.9677133469215410D+00,0.9938009662969447D+00/

        data w5/
     1  0.2751242918137683D-05,0.6942436018763683D-04,
     2  0.5994851026071711D-03,0.2805079843198273D-02,
     3  0.8651226551181051D-02,0.1965581984197020D-01,
     4  0.3540520091126813D-01,0.5333578971688024D-01,
     5  0.7007122960489411D-01,0.8309910831615371D-01,
     6  0.9153800855813067D-01,0.9579327055178783D-01,
     7  0.9669265534467965D-01,0.9474220338035766D-01,
     8  0.8986013737736566D-01,0.8162464957194391D-01,
     9  0.6974330324732884D-01,0.5436120873463431D-01,
     *  0.3608383211826662D-01,0.1586561562424619D-01/

        data r6/
     1  0.4993930594639425D-06,0.2050163806211554D-04,
     2  0.2496919158471735D-03,0.1570489954211665D-02,
     3  0.6380132290668364D-02,0.1893998312956402D-01,
     4  0.4448306120793684D-01,0.8739756891453328D-01,
     5  0.1496584944607259D+00,0.2303271877884174D+00,
     6  0.3260152774379400D+00,0.4317347071124248D+00,
     7  0.5417033019250586D+00,0.6499763860809572D+00,
     8  0.7509068000601279D+00,0.8394532452520243D+00,
     9  0.9113703512974177D+00,0.9633241612447059D+00,
     *  0.9929742957232059D+00/

        data w6/
     1  0.2384441531245178D-05,0.6066866135609926D-04,
     2  0.5288901016123203D-03,0.2503967434114089D-02,
     3  0.7841318381144415D-02,0.1818426526066423D-01,
     4  0.3366159612698156D-01,0.5250441073882110D-01,
     5  0.7185158853746203D-01,0.8890583867720622D-01,
     6  0.1016253688831803D+00,0.1088393166920480D+00,
     7  0.1100994630852908D+00,0.1055004550393570D+00,
     8  0.9551710275819420D-01,0.8086638894589125D-01,
     9  0.6241092592586764D-01,0.4110493064129084D-01,
     *  0.1799111966798666D-01/
c        
c         %%%%%%%%%%%%%%%%%%%%%

ccc        call prin2('alpha = *',alpha,1)
        dtwo = 2
        if (alpha .ge. 1/dtwo) then
        do 1000 i=1,nqs(6)
        roots(i) =r6(i)
        weights(i) = w6(i)
 1000 continue

        nquads = nqs(6)
        endif

        if (alpha .ge. dtwo**(-2)) then
        if (alpha .lt. dtwo**(-1)) then

        do 2000 i=1,nqs(5)
        roots(i) =r5(i)
        weights(i) = w5(i)
 2000 continue

        nquads = nqs(5)

        endif
        endif


        if (alpha .ge. dtwo**(-3)) then
        if (alpha .lt. dtwo**(-2)) then

        do 3000 i=1,nqs(4)
        roots(i) =r4(i)
        weights(i) = w4(i)
 3000 continue

        nquads = nqs(4)

        endif
        endif


        if (alpha .ge. dtwo**(-4)) then
        if (alpha .lt. dtwo**(-3)) then

        do 4000 i=1,nqs(3)
        roots(i) =r3(i)
        weights(i) = w3(i)
 4000 continue

        nquads = nqs(3)

        endif
        endif


        if (alpha .ge. dtwo**(-5)) then
        if (alpha .lt. dtwo**(-4)) then

        do 5000 i=1,nqs(2)
        roots(i) =r2(i)
        weights(i) = w2(i)
 5000 continue

        nquads = nqs(2)

        endif
        endif


        if (alpha .lt. dtwo**(-5)) then

        do 6000 i=1,nqs(1)
        roots(i) =r1(i)
        weights(i) = w1(i)
 6000 continue

        nquads = nqs(1)


        endif

        return
        end
c
c
c
c
c
        subroutine lap_load_quads30(alpha,roots,weights,nquads)
        implicit real *8 (a-h,o-z)
        dimension r1(31),r2(29),r3(26),r4(22),r5(20),r6(19),
     1          w1(31),w2(29),w3(26),w4(22),w5(20),w6(19),nqs(6),
     2          roots(1),weights(1)

        data nqs/31,29,26,22,20,19/


        data r1/
     1  0.6347043439709046D-06,0.2574520389574903D-04,
     2  0.3099166027738289D-03,0.1926184285225288D-02,
     3  0.7725666957265455D-02,0.2260561375511627D-01,
     4  0.5219977472832835D-01,0.1004985703561920D+00,
     5  0.1679736059307111D+00,0.2512787370298330D+00,
     6  0.3443257534867958D+00,0.4398212125335920D+00,
     7  0.5305566297997550D+00,0.6103897502873278D+00,
     8  0.6753006483658965D+00,0.7244688280009890D+00,
     9  0.7601861300287476D+00,0.7861199887709169D+00,
     *  0.8055864726992147D+00,0.8209871517850033D+00,
     1  0.8339755322645063D+00,0.8457575384267297D+00,
     2  0.8573209159309071D+00,0.8695835253978650D+00,
     3  0.8834773228747195D+00,0.8999324556214982D+00,
     4  0.9196058710907990D+00,0.9421331825872929D+00,
     5  0.9652057993250502D+00,0.9848019386145025D+00,
     *  0.9970118831472794D+00/

        data w1/
     1  0.3022591289008485D-05,0.7588352234910224D-04,
     2  0.6528297011162294D-03,0.3047498657429362D-02,
     3  0.9390913021997492D-02,0.2135218356183039D-01,
     4  0.3853535466816003D-01,0.5814501161681592D-01,
     5  0.7622683595615383D-01,0.8934125430859026D-01,
     6  0.9551562770404719D-01,0.9425878580969671D-01,
     7  0.8617120739896962D-01,0.7280337792696338D-01,
     8  0.5688061883593131D-01,0.4187268762930973D-01,
     9  0.3020497796312449D-01,0.2222282266349279D-01,
     *  0.1710611654587888D-01,0.1396545990836117D-01,
     1  0.1220850129250471D-01,0.1151772740931632D-01,
     2  0.1176054977357754D-01,0.1291990789421640D-01,
     3  0.1502649465315413D-01,0.1800592406711487D-01,
     4  0.2129964136590636D-01,0.2336877540497981D-01,
     5  0.2206908869258860D-01,0.1644221340695313D-01,
     *  0.7608706048181227D-02/

        data r2/
     1  0.4301675742577534D-06,0.1813279038993976D-04,
     2  0.2263256467214214D-03,0.1456058342676344D-02,
     3  0.6034139558992707D-02,0.1820133135957393D-01,
     4  0.4320355289531850D-01,0.8520781430332021D-01,
     5  0.1453234862166355D+00,0.2209237147809609D+00,
     6  0.3064398205281922D+00,0.3949667724529043D+00,
     7  0.4799173007111450D+00,0.5562995729400411D+00,
     8  0.6214376812276748D+00,0.6749861117708821D+00,
     9  0.7182728473119619D+00,0.7533852552203523D+00,
     *  0.7824759818228658D+00,0.8074502402940757D+00,
     1  0.8299234371641117D+00,0.8512758176487176D+00,
     2  0.8726775783381328D+00,0.8950066231405831D+00,
     3  0.9186061608046384D+00,0.9428950570040034D+00,
     4  0.9660444996725239D+00,0.9851436765264900D+00,
     *  0.9970608559982486D+00/

        data w2/
     1  0.2065824136582489D-05,0.5411931274756986D-04,
     2  0.4847475436357778D-03,0.2352710004831880D-02,
     3  0.7525720663193385D-02,0.1772141376122077D-01,
     4  0.3301379381186581D-01,0.5119120448039203D-01,
     5  0.6858714995525523D-01,0.8165643139216150D-01,
     6  0.8819643344499534D-01,0.8773800650831560D-01,
     7  0.8131862649333281D-01,0.7099961714338255D-01,
     8  0.5923113893375348D-01,0.4810362309983500D-01,
     9  0.3882963788291274D-01,0.3175879697985385D-01,
     *  0.2674044544506172D-01,0.2347609938708965D-01,
     1  0.2169974757704729D-01,0.2120128064610479D-01,
     2  0.2175256163495730D-01,0.2296968549152986D-01,
     3  0.2413913339167453D-01,0.2413184279108453D-01,
     4  0.2166639148306756D-01,0.1598877638361865D-01,
     *  0.7468798532942217D-02/

        data r3/
     1  0.2857238825099065D-06,0.1274691567091325D-04,
     2  0.1670879929395877D-03,0.1122138211890672D-02,
     3  0.4828040791637550D-02,0.1504138356513247D-01,
     4  0.3669193438466320D-01,0.7402452344919650D-01,
     5  0.1286083796105313D+00,0.1984740209611423D+00,
     6  0.2787483680195813D+00,0.3632285155980007D+00,
     7  0.4460587146424855D+00,0.5229272351405734D+00,
     8  0.5915540647926858D+00,0.6515151686363090D+00,
     9  0.7036684724104984D+00,0.7495283700984124D+00,
     *  0.7907943860818279D+00,0.8290251563525646D+00,
     1  0.8653479271166334D+00,0.9001350885086119D+00,
     2  0.9327146845912911D+00,0.9613238243201492D+00,
     3  0.9834915235783085D+00,0.9967795601049178D+00/

        data w3/
     1  0.1389828124943634D-05,0.3871228341461501D-04,
     2  0.3655277571892601D-03,0.1858242890593763D-02,
     3  0.6189841172579157D-02,0.1509257883934143D-01,
     4  0.2895268798056749D-01,0.4598926756192331D-01,
     5  0.6283794435367779D-01,0.7605271025880855D-01,
     6  0.8343281551002737D-01,0.8453108132331371D-01,
     7  0.8040267923824843D-01,0.7295797072721215D-01,
     8  0.6423552654553297D-01,0.5584061269946737D-01,
     9  0.4872616561139952D-01,0.4328199631355923D-01,
     *  0.3951444948149406D-01,0.3714050260316062D-01,
     1  0.3557249436233587D-01,0.3389118735095948D-01,
     2  0.3096763166767110D-01,0.2582559090994451D-01,
     3  0.1808944670118181D-01,0.8210946028271489D-02/

        data r4/
     1  0.6373572627201451D-06,0.2578815926223228D-04,
     2  0.3094246119736003D-03,0.1915548515021225D-02,
     3  0.7646707070476509D-02,0.2224606187151634D-01,
     4  0.5100812612466716D-01,0.9737203422769555D-01,
     5  0.1611635389571161D+00,0.2386199713696203D+00,
     6  0.3239582269666650D+00,0.4113871272301014D+00,
     7  0.4965501718189550D+00,0.5770104004400132D+00,
     8  0.6519681622954621D+00,0.7215912515501470D+00,
     9  0.7862267909597298D+00,0.8456673735691361D+00,
     *  0.8986864305679547D+00,0.9430639165149655D+00,
     1  0.9761061376884330D+00,0.9953845753493910D+00/

        data w4/
     1  0.3033898710839364D-05,0.7593533617895676D-04,
     2  0.6506793649452663D-03,0.3022410669714423D-02,
     3  0.9255164030482921D-02,0.2087089909328037D-01,
     4  0.3725878498801745D-01,0.5543869828400351D-01,
     5  0.7148801304341262D-01,0.8242482539855107D-01,
     6  0.8726155563378644D-01,0.8686272928808087D-01,
     7  0.8306294527061240D-01,0.7773601799638610D-01,
     8  0.7222093423298281D-01,0.6709585100201182D-01,
     9  0.6214905525150949D-01,0.5652730711447951D-01,
     *  0.4913133874523959D-01,0.3916155222608494D-01,
     1  0.2650685575956718D-01,0.1179541337196141D-01/

        data r5/
     1  0.5257015763512882D-06,0.2159871229671709D-04,
     2  0.2631076846343377D-03,0.1653625372133167D-02,
     3  0.6702005603094575D-02,0.1979910163504920D-01,
     4  0.4611892347198273D-01,0.8951380273220150D-01,
     5  0.1508539426172521D+00,0.2278760215375481D+00,
     6  0.3163540605476136D+00,0.4116054836614612D+00,
     7  0.5094732067410363D+00,0.6065177562790542D+00,
     8  0.6996609900583225D+00,0.7857458766989859D+00,
     9  0.8613883693327729D+00,0.9231891386723091D+00,
     *  0.9681192349619195D+00,0.9938819541201059D+00/

        data w5/
     1  0.2510590490428934D-05,0.6392678092610135D-04,
     2  0.5572574949463175D-03,0.2634157338286593D-02,
     3  0.8213771658537862D-02,0.1888205006311841D-01,
     4  0.3443025704892631D-01,0.5250634443956796D-01,
     5  0.6977439868764557D-01,0.8353810224824491D-01,
     6  0.9262061719888750D-01,0.9718687441741939D-01,
     7  0.9798249766545348D-01,0.9560192077314919D-01,
     8  0.9016097466375183D-01,0.8144069448890593D-01,
     9  0.6927077141657823D-01,0.5381980933075433D-01,
     *  0.3565292427376523D-01,0.1566013942064443D-01/

        data r6/
     1  0.4803909588321807D-06,0.1983010810871589D-04,
     2  0.2428121947910858D-03,0.1535301136680862D-02,
     3  0.6269099378906114D-02,0.1869936185312683D-01,
     4  0.4410467669913465D-01,0.8696139633884347D-01,
     5  0.1493169155116840D+00,0.2302311533694510D+00,
     6  0.3262348819127840D+00,0.4322298264720080D+00,
     7  0.5423552784758903D+00,0.6506467765142488D+00,
     8  0.7514844238118467D+00,0.8398752361625807D+00,
     9  0.9116230510620715D+00,0.9634340922203856D+00,
     *  0.9929959042686297D+00/

        data w6/
     1  0.2296413529518938D-05,0.5878863417442658D-04,
     2  0.5156113144020707D-03,0.2455867568140601D-02,
     3  0.7736005521240065D-02,0.1803812307130093D-01,
     4  0.3354802266957854D-01,0.5251559271647599D-01,
     5  0.7202945493282377D-01,0.8920462118697665D-01,
     6  0.1019383487889765D+00,0.1090638231541256D+00,
     7  0.1101856445080274D+00,0.1054561296675306D+00,
     8  0.9538428867284641D-01,0.8069629247346955D-01,
     9  0.6224918064886505D-01,0.4098579659275122D-01,
     *  0.1793611146476510D-01/
c        
c         %%%%%%%%%%%%%%%%%%%%%

ccc        call prin2('alpha = *',alpha,1)
        dtwo = 2
        if (alpha .ge. 1/dtwo) then
        do 1000 i=1,nqs(6)
        roots(i) =r6(i)
        weights(i) = w6(i)
 1000 continue

        nquads = nqs(6)
        endif

        if (alpha .ge. dtwo**(-2)) then
        if (alpha .lt. dtwo**(-1)) then

        do 2000 i=1,nqs(5)
        roots(i) =r5(i)
        weights(i) = w5(i)
 2000 continue

        nquads = nqs(5)

        endif
        endif


        if (alpha .ge. dtwo**(-3)) then
        if (alpha .lt. dtwo**(-2)) then

        do 3000 i=1,nqs(4)
        roots(i) =r4(i)
        weights(i) = w4(i)
 3000 continue

        nquads = nqs(4)

        endif
        endif


        if (alpha .ge. dtwo**(-4)) then
        if (alpha .lt. dtwo**(-3)) then

        do 4000 i=1,nqs(3)
        roots(i) =r3(i)
        weights(i) = w3(i)
 4000 continue

        nquads = nqs(3)

        endif
        endif


        if (alpha .ge. dtwo**(-5)) then
        if (alpha .lt. dtwo**(-4)) then

        do 5000 i=1,nqs(2)
        roots(i) =r2(i)
        weights(i) = w2(i)
 5000 continue

        nquads = nqs(2)

        endif
        endif


        if (alpha .lt. dtwo**(-5)) then

        do 6000 i=1,nqs(1)
        roots(i) =r1(i)
        weights(i) = w1(i)
 6000 continue

        nquads = nqs(1)

        endif

        return
        end
c
c
c
c
c
        subroutine lap_load_quads31(alpha,roots,weights,nquads)
        implicit real *8 (a-h,o-z)
        dimension r1(31),r2(28),r3(25),r4(22),r5(20),r6(19),
     1          w1(31),w2(28),w3(25),w4(22),w5(20),w6(19),nqs(6),
     2          roots(1),weights(1)

        data nqs/31,28,25,22,20,19/


        data r1/
     1  0.4076416455304407D-06,0.1721950337974619D-04,
     2  0.2154247709585447D-03,0.1389661236601389D-02,
     3  0.5777861450490099D-02,0.1750057084882920D-01,
     4  0.4176261285457318D-01,0.8293486372723580D-01,
     5  0.1426845855219479D+00,0.2192451708691793D+00,
     6  0.3079768040724083D+00,0.4026072582637711D+00,
     7  0.4964663453399308D+00,0.5834192380396842D+00,
     8  0.6586286006354620D+00,0.7193749493027313D+00,
     9  0.7656708741511561D+00,0.7998529462941312D+00,
     *  0.8251527947607973D+00,0.8444540191095777D+00,
     1  0.8598940141352440D+00,0.8729859301471661D+00,
     2  0.8848488200692789D+00,0.8963876227423086D+00,
     3  0.9084012993715866D+00,0.9216154602835482D+00,
     4  0.9365926573697257D+00,0.9534067704251435D+00,
     5  0.9710296226347510D+00,0.9868643880803648D+00,
     *  0.9973457915080912D+00/

        data w1/
     1  0.1958486499245625D-05,0.5143112590379398D-04,
     2  0.4619327902298030D-03,0.2249511950777541D-02,
     3  0.7227200201458268D-02,0.1712104487031926D-01,
     4  0.3216567969938171D-01,0.5046909837987535D-01,
     5  0.6872477664920189D-01,0.8360058184466310D-01,
     6  0.9279603803813888D-01,0.9533805037638175D-01,
     7  0.9134576976694244D-01,0.8173968919573662D-01,
     8  0.6821281362623704D-01,0.5328837745469377D-01,
     9  0.3972084361692820D-01,0.2920423268089686D-01,
     *  0.2188137710801969D-01,0.1707421873296612D-01,
     1  0.1405382983468113D-01,0.1231350447616148D-01,
     2  0.1156130816645852D-01,0.1164881816769713D-01,
     3  0.1250066850618765D-01,0.1402651388007813D-01,
     4  0.1595127200935627D-01,0.1751634295723339D-01,
     5  0.1728707028848859D-01,0.1375552571145391D-01,
     *  0.6710519406952467D-02/

        data r2/
     1  0.5563028753362866D-06,0.2281822823266526D-04,
     2  0.2775489723332310D-03,0.1741968290490248D-02,
     3  0.7050940146746678D-02,0.2080468483345413D-01,
     4  0.4839947120903720D-01,0.9377407324429066D-01,
     5  0.1575408455248458D+00,0.2366056495099069D+00,
     6  0.3252019712846382D+00,0.4165192001244042D+00,
     7  0.5041834538573338D+00,0.5833063329078221D+00,
     8  0.6510983200367029D+00,0.7069902436180903D+00,
     9  0.7521656744022481D+00,0.7887101375691506D+00,
     *  0.8188521632817387D+00,0.8445721620826993D+00,
     1  0.8675103597244186D+00,0.8889861262084672D+00,
     2  0.9099842511160055D+00,0.9310264201354400D+00,
     3  0.9519130085359313D+00,0.9714541581080766D+00,
     4  0.9874961618787801D+00,0.9975221967815566D+00/

        data w2/
     1  0.2655744251026306D-05,0.6750073485722158D-04,
     2  0.5874451948710605D-03,0.2772548797376429D-02,
     3  0.8632647252129533D-02,0.1981604197240367D-01,
     4  0.3606413438600758D-01,0.5479302109862477D-01,
     5  0.7220525061376623D-01,0.8493421652747186D-01,
     6  0.9109806923731124D-01,0.9045676594924847D-01,
     7  0.8404665643761587D-01,0.7373170631735923D-01,
     8  0.6176569030842832D-01,0.5023101482752353D-01,
     9  0.4048078601717982D-01,0.3298577939163087D-01,
     *  0.2762936628034520D-01,0.2408334051851682D-01,
     1  0.2201298835027938D-01,0.2110475366411964D-01,
     2  0.2098298336662735D-01,0.2107337086491441D-01,
     3  0.2050195931339069D-01,0.1821100175705329D-01,
     4  0.1343487534446034D-01,0.6293429732236143D-02/

        data r3/
     1  0.5646099516706670D-06,0.2311146756071291D-04,
     2  0.2804652131252650D-03,0.1755654316756272D-02,
     3  0.7085006514108762D-02,0.2083236947406282D-01,
     4  0.4826632202436619D-01,0.9307339448410833D-01,
     5  0.1555279176459174D+00,0.2322470485167455D+00,
     6  0.3174318419056290D+00,0.4046973696065461D+00,
     7  0.4886231349903918D+00,0.5656033812651928D+00,
     8  0.6339872650257365D+00,0.6937216209798185D+00,
     9  0.7457816883341013D+00,0.7916250032992920D+00,
     *  0.8327624448828556D+00,0.8704117057288457D+00,
     1  0.9051734356866757D+00,0.9367559641481421D+00,
     2  0.9638985027096422D+00,0.9846549119950869D+00,
     *  0.9970126524876983D+00/

        data w3/
     1  0.2694306487937055D-05,0.6831786857971572D-04,
     2  0.5929267641886666D-03,0.2789413600417785D-02,
     3  0.8651479766851180D-02,0.1976386031589069D-01,
     4  0.3575202750412758D-01,0.5391142105598999D-01,
     5  0.7041060403650626D-01,0.8202843510029242D-01,
     6  0.8725020088309007D-01,0.8636333148668492D-01,
     7  0.8088546407717991D-01,0.7280656921601382D-01,
     8  0.6396243103676429D-01,0.5567813927782800D-01,
     9  0.4869193549995641D-01,0.4325209282062574D-01,
     *  0.3923117504449161D-01,0.3617244507836440D-01,
     1  0.3330417595416865D-01,0.2964597874410621D-01,
     2  0.2430374674431501D-01,0.1686099569819848D-01,
     *  0.7620138118880240D-02/

        data r4/
     1  0.5411453594322133D-06,0.2222880199437904D-04,
     2  0.2706115553010992D-03,0.1699067474185927D-02,
     3  0.6876637640528484D-02,0.2027768051165234D-01,
     4  0.4711564612672755D-01,0.9112105914717071D-01,
     5  0.1527474426450709D+00,0.2289364838110622D+00,
     6  0.3143665654907400D+00,0.4032761490271485D+00,
     7  0.4909365430617633D+00,0.5743229824847835D+00,
     8  0.6520339306751371D+00,0.7237549826462226D+00,
     9  0.7895475257190713D+00,0.8491617382731595D+00,
     *  0.9015728600744363D+00,0.9449321568461923D+00,
     1  0.9769576242960940D+00,0.9955559082960139D+00/

        data w4/
     1  0.2584372252880287D-05,0.6578212433192013D-04,
     2  0.5729264175121136D-03,0.2704516050535537D-02,
     3  0.8416854855289785D-02,0.1929546570553039D-01,
     4  0.3503519656688917D-01,0.5305471399036720D-01,
     5  0.6966497679203463D-01,0.8179337961283477D-01,
     6  0.8808091187170000D-01,0.8893919620493240D-01,
     7  0.8587963876808547D-01,0.8066927273133390D-01,
     8  0.7471249922707744D-01,0.6875175364918703D-01,
     9  0.6279996855021778D-01,0.5625858631332315D-01,
     *  0.4824842621901973D-01,0.3807875845500618D-01,
     1  0.2561301807479613D-01,0.1136157344774236D-01/

        data r5/
     1  0.4795057839865142D-06,0.1986000638687454D-04,
     2  0.2438555605724912D-03,0.1544893508093451D-02,
     3  0.6311947404826422D-02,0.1879934755343166D-01,
     4  0.4415149796405072D-01,0.8639883376699688D-01,
     5  0.1467610528147682D+00,0.2233230716990406D+00,
     6  0.3120210471242310D+00,0.4080787246572518D+00,
     7  0.5070591055251580D+00,0.6051986343832748D+00,
     8  0.6991811407115708D+00,0.7857594452156183D+00,
     9  0.8615897401234804D+00,0.9233810909359297D+00,
     *  0.9682212256331953D+00,0.9939038336133024D+00/

        data w5/
     1  0.2293916437455144D-05,0.5893861456376585D-04,
     2  0.5184188933863847D-03,0.2473323106759722D-02,
     3  0.7786855560012618D-02,0.1808262121079973D-01,
     4  0.3332463062397388D-01,0.5137836722504214D-01,
     5  0.6900801058593558D-01,0.8341260009309646D-01,
     6  0.9316902221333234D-01,0.9820174181663188D-01,
     7  0.9913662151245557D-01,0.9659702443381660D-01,
     8  0.9082978166585518D-01,0.8176659781487945D-01,
     9  0.6933958960713260D-01,0.5375215843718920D-01,
     *  0.3555589840271782D-01,0.1560550426598162D-01/

        data r6/
     1  0.3986114598497154D-06,0.1658421370019830D-04,
     2  0.2050906729861476D-03,0.1313118680763639D-02,
     3  0.5444720508405340D-02,0.1653326654239639D-01,
     4  0.3976471792101271D-01,0.7997602276210841D-01,
     5  0.1399170725609546D+00,0.2193092138890591D+00,
     6  0.3149829791310218D+00,0.4217373905492148D+00,
     7  0.5333761630486950D+00,0.6435579245995098D+00,
     8  0.7463505702084536D+00,0.8365375374989927D+00,
     9  0.9097731941485156D+00,0.9626673781693080D+00,
     *  0.9928489515948694D+00/

        data w6/
     1  0.1908305423958106D-05,0.4931135891142304D-04,
     2  0.4377953875095275D-03,0.2118749106017873D-02,
     3  0.6809659464701515D-02,0.1626099896580968D-01,
     4  0.3103957011676431D-01,0.4985768595743345D-01,
     5  0.6997459631915587D-01,0.8826566169130380D-01,
     6  0.1021968935161320D+00,0.1102648328873329D+00,
     7  0.1119467367565429D+00,0.1074232559793713D+00,
     8  0.9729021447117840D-01,0.8235813904410744D-01,
     9  0.6354721136657695D-01,0.4184439522066689D-01,
     *  0.1831238408505978D-01/
c       
c        
c         %%%%%%%%%%%%%%%%%%%%%

ccc        call prin2('alpha = *',alpha,1)
        dtwo = 2
        if (alpha .ge. 1/dtwo) then
        do 1000 i=1,nqs(6)
        roots(i) =r6(i)
        weights(i) = w6(i)
 1000 continue

        nquads = nqs(6)
        endif

        if (alpha .ge. dtwo**(-2)) then
        if (alpha .lt. dtwo**(-1)) then

        do 2000 i=1,nqs(5)
        roots(i) =r5(i)
        weights(i) = w5(i)
 2000 continue

        nquads = nqs(5)

        endif
        endif


        if (alpha .ge. dtwo**(-3)) then
        if (alpha .lt. dtwo**(-2)) then

        do 3000 i=1,nqs(4)
        roots(i) =r4(i)
        weights(i) = w4(i)
 3000 continue

        nquads = nqs(4)

        endif
        endif


        if (alpha .ge. dtwo**(-4)) then
        if (alpha .lt. dtwo**(-3)) then

        do 4000 i=1,nqs(3)
        roots(i) =r3(i)
        weights(i) = w3(i)
 4000 continue

        nquads = nqs(3)

        endif
        endif


        if (alpha .ge. dtwo**(-5)) then
        if (alpha .lt. dtwo**(-4)) then

        do 5000 i=1,nqs(2)
        roots(i) =r2(i)
        weights(i) = w2(i)
 5000 continue

        nquads = nqs(2)

        endif
        endif


        if (alpha .lt. dtwo**(-5)) then

        do 6000 i=1,nqs(1)
        roots(i) =r1(i)
        weights(i) = w1(i)
 6000 continue

        nquads = nqs(1)

        endif

        return
        end
c
c
c
c
c
        subroutine lap_load_quads32(alpha,roots,weights,nquads)
        implicit real *8 (a-h,o-z)
        dimension r1(30),r2(27),r3(24),r4(22),r5(20),r6(19),
     1          w1(30),w2(27),w3(24),w4(22),w5(20),w6(19),nqs(6),
     2          roots(1),weights(1)

        data nqs/30,27,24,22,20,19/


        data r1/
     1  0.6167814504109241D-06,0.2513091538590471D-04,
     2  0.3037844130290187D-03,0.1895448179902020D-02,
     3  0.7630003826551309D-02,0.2240006605781766D-01,
     4  0.5188057466120510D-01,0.1001513233526616D+00,
     5  0.1677940423651845D+00,0.2515771254740608D+00,
     6  0.3455734366980637D+00,0.4428033406791959D+00,
     7  0.5366229948593404D+00,0.6216110506616579D+00,
     8  0.6941041163699947D+00,0.7525552598536620D+00,
     9  0.7975920139064838D+00,0.8314617503048543D+00,
     *  0.8569949979422992D+00,0.8767430713002978D+00,
     1  0.8926640304045233D+00,0.9061880307262907D+00,
     2  0.9183847070006125D+00,0.9300917645267860D+00,
     3  0.9419620927941174D+00,0.9544058944309411D+00,
     4  0.9673986792104856D+00,0.9801898102993850D+00,
     5  0.9911618661114989D+00,0.9982294246925374D+00/

        data w1/
     1  0.2940152993498823D-05,0.7418117347781941D-04,
     2  0.6411366962524341D-03,0.3006008851516150D-02,
     3  0.9301186053853743D-02,0.2122895308289452D-01,
     4  0.3844674991037852D-01,0.5819628286389732D-01,
     5  0.7652934460275401D-01,0.9001956880231222D-01,
     6  0.9678613804441753D-01,0.9655210726325677D-01,
     7  0.9017528713763429D-01,0.7918398643573173D-01,
     8  0.6554224909797730D-01,0.5147407209038518D-01,
     9  0.3899336785001909D-01,0.2923320511321276D-01,
     *  0.2226306699598233D-01,0.1755759143719192D-01,
     1  0.1451959845744799D-01,0.1270446456143821D-01,
     2  0.1182767633880289D-01,0.1169612581954440D-01,
     3  0.1211609873574395D-01,0.1277285726796735D-01,
     4  0.1309117709959745D-01,0.1221498757634064D-01,
     5  0.9364277771112302D-02,0.4485312715865688D-02/

        data r2/
     1  0.5996292173212608D-06,0.2445721349936103D-04,
     2  0.2958781795552633D-03,0.1847167042905562D-02,
     3  0.7437790109639831D-02,0.2183445665094296D-01,
     4  0.5054627823345715D-01,0.9748194696192808D-01,
     5  0.1630832745041082D+00,0.2440449274117304D+00,
     6  0.3344752652300582D+00,0.4275966516193113D+00,
     7  0.5171980990847065D+00,0.5985611356712726D+00,
     8  0.6689232356526535D+00,0.7275419821770731D+00,
     9  0.7753261313358100D+00,0.8141365107838609D+00,
     *  0.8460897326285221D+00,0.8731315292342707D+00,
     1  0.8968745635056141D+00,0.9185463573415911D+00,
     2  0.9389083004196733D+00,0.9580873891514358D+00,
     3  0.9753855469267956D+00,0.9892840361962078D+00,
     *  0.9978823062160877D+00/

        data w2/
     1  0.2859082688324955D-05,0.7221488420755479D-04,
     2  0.6246550286647422D-03,0.2930248019133705D-02,
     3  0.9067809530472771D-02,0.2068726905289894D-01,
     4  0.3742192967985615D-01,0.5652707094887736D-01,
     5  0.7410278669682492D-01,0.8680821240535316D-01,
     6  0.9289807966816103D-01,0.9229442289720135D-01,
     7  0.8611453872927100D-01,0.7614982214606411D-01,
     8  0.6445696358476991D-01,0.5294461381669734D-01,
     9  0.4294616384688865D-01,0.3503437565325410D-01,
     *  0.2919861974831395D-01,0.2515427201361333D-01,
     1  0.2253757982281067D-01,0.2093583909056368D-01,
     2  0.1981294919598218D-01,0.1842954390447654D-01,
     3  0.1590628233310406D-01,0.1155896165163578D-01,
     *  0.5381916568214713D-02/

        data r3/
     1  0.6592070285338512D-06,0.2665280670244157D-04,
     2  0.3197357647618598D-03,0.1979684761346101D-02,
     3  0.7906429681778176D-02,0.2302114620099773D-01,
     4  0.5285496242340828D-01,0.1010802447407185D+00,
     5  0.1676586748675627D+00,0.2487363438801719D+00,
     6  0.3380680557467737D+00,0.4289680947147424D+00,
     7  0.5158921356002563D+00,0.5952540256876615D+00,
     8  0.6654978500895407D+00,0.7266671720754392D+00,
     9  0.7797815860892944D+00,0.8262556830370083D+00,
     *  0.8674345215047937D+00,0.9042083310068104D+00,
     1  0.9366915911153690D+00,0.9640690320511127D+00,
     2  0.9847747155624789D+00,0.9970399605210081D+00/

        data w3/
     1  0.3137186406543054D-05,0.7847014041613969D-04,
     2  0.6723606134125964D-03,0.3124212917279698D-02,
     3  0.9574587582648365D-02,0.2162270040919422D-01,
     4  0.3869089597226118D-01,0.5775417853329507D-01,
     5  0.7473381683184785D-01,0.8634586854026957D-01,
     6  0.9117602845251492D-01,0.8968651818445266D-01,
     7  0.8356329827748267D-01,0.7490914283942563D-01,
     8  0.6559676074249176D-01,0.5692165005074820D-01,
     9  0.4955053376666539D-01,0.4362817634632874D-01,
     *  0.3888348309952265D-01,0.3468627679165746D-01,
     1  0.3014152668259227D-01,0.2434268296200044D-01,
     2  0.1676108180150071D-01,0.7552611275585265D-02/

        data r4/
     1  0.4619559660916242D-06,0.1920444267003952D-04,
     2  0.2365675904819423D-03,0.1502974799817967D-02,
     3  0.6155714963332801D-02,0.1837031963512518D-01,
     4  0.4320074349477810D-01,0.8456364336376991D-01,
     5  0.1434630864935704D+00,0.2175493006561340D+00,
     6  0.3020660899292782D+00,0.3914537641439835D+00,
     7  0.4807853618883916D+00,0.5665686167911344D+00,
     8  0.6468585171525059D+00,0.7208764740406486D+00,
     9  0.7883854713828021D+00,0.8490346167737238D+00,
     *  0.9018792007728019D+00,0.9452675498591404D+00,
     1  0.9771433764655555D+00,0.9955963823426140D+00/

        data w4/
     1  0.2211829712769110D-05,0.5705920130290083D-04,
     2  0.5036400640200850D-03,0.2410188624383694D-02,
     3  0.7607858981507011D-02,0.1770101895136572D-01,
     4  0.3264517827955814D-01,0.5025527776133931D-01,
     5  0.6713410772088713D-01,0.8021677030693661D-01,
     6  0.8786513638379029D-01,0.9007563049799626D-01,
     7  0.8800420249248235D-01,0.8324511711975063D-01,
     8  0.7721581429609145D-01,0.7079200394175515D-01,
     9  0.6417603143895593D-01,0.5696685878610809D-01,
     *  0.4844151545615199D-01,0.3798263216514333D-01,
     1  0.2544090721742595D-01,0.1126083848333521D-01/

        data r5/
     1  0.3948836159428872D-06,0.1667482692585935D-04,
     2  0.2085883589497226D-03,0.1345695228959303D-02,
     3  0.5596752257047088D-02,0.1696070399477369D-01,
     4  0.4050537838030432D-01,0.8053275320883867D-01,
     5  0.1388225468501570D+00,0.2140396182320474D+00,
     6  0.3024430236369007D+00,0.3992172594410658D+00,
     7  0.4996125448697424D+00,0.5994593171523896D+00,
     8  0.6950998688271138D+00,0.7830846328396117D+00,
     9  0.8599987345790982D+00,0.9225608885584452D+00,
     *  0.9678984955321362D+00,0.9938437717352428D+00/

        data w5/
     1  0.1896998175179103D-05,0.4980067950946466D-04,
     2  0.4472732704260754D-03,0.2178611945130227D-02,
     3  0.7002925652485258D-02,0.1660352066356786D-01,
     4  0.3123509560155962D-01,0.4912637273501092D-01,
     5  0.6721692057980802D-01,0.8256939123480991D-01,
     6  0.9341477388741147D-01,0.9933367385881762D-01,
     7  0.1007640591619257D+00,0.9832591574536076D-01,
     8  0.9238349162384548D-01,0.8301491978739675D-01,
     9  0.7026137819754099D-01,0.5437909718777273D-01,
     *  0.3593051933348023D-01,0.1576036185596574D-01/

        data r6/
     1  0.3198327071840780D-06,0.1350585035269355D-04,
     2  0.1703632178963969D-03,0.1117690239439824D-02,
     3  0.4762800401616386D-02,0.1487053944266351D-01,
     4  0.3670079772005431D-01,0.7545800027115629D-01,
     5  0.1343435875027258D+00,0.2133487667798572D+00,
     6  0.3092957969692347D+00,0.4167863741834040D+00,
     7  0.5293841128467109D+00,0.6405590628567254D+00,
     8  0.7442637628475136D+00,0.8352219591142483D+00,
     9  0.9090602697190337D+00,0.9623762880087672D+00,
     *  0.9927936110884602D+00/

        data w6/
     1  0.1535179292365020D-05,0.4039174227612949D-04,
     2  0.3674575862197355D-03,0.1832723302828947D-02,
     3  0.6090230096540609D-02,0.1502997285776353D-01,
     4  0.2953360317720274D-01,0.4853578222133719D-01,
     5  0.6923214795489824D-01,0.8823093654122049D-01,
     6  0.1027425364568360D+00,0.1111499060828131D+00,
     7  0.1129479097984963D+00,0.1083896111156545D+00,
     8  0.9813844218583740D-01,0.8304801563046155D-01,
     9  0.6406073631339929D-01,0.4217371041083859D-01,
     *  0.1845435134608320D-01/
c        
c         %%%%%%%%%%%%%%%%%%%%%

ccc        call prin2('alpha = *',alpha,1)
        dtwo = 2
        if (alpha .ge. 1/dtwo) then
        do 1000 i=1,nqs(6)
        roots(i) =r6(i)
        weights(i) = w6(i)
 1000 continue

        nquads = nqs(6)
        endif

        if (alpha .ge. dtwo**(-2)) then
        if (alpha .lt. dtwo**(-1)) then

        do 2000 i=1,nqs(5)
        roots(i) =r5(i)
        weights(i) = w5(i)
 2000 continue

        nquads = nqs(5)

        endif
        endif


        if (alpha .ge. dtwo**(-3)) then
        if (alpha .lt. dtwo**(-2)) then

        do 3000 i=1,nqs(4)
        roots(i) =r4(i)
        weights(i) = w4(i)
 3000 continue

        nquads = nqs(4)

        endif
        endif


        if (alpha .ge. dtwo**(-4)) then
        if (alpha .lt. dtwo**(-3)) then

        do 4000 i=1,nqs(3)
        roots(i) =r3(i)
        weights(i) = w3(i)
 4000 continue

        nquads = nqs(3)

        endif
        endif


        if (alpha .ge. dtwo**(-5)) then
        if (alpha .lt. dtwo**(-4)) then

        do 5000 i=1,nqs(2)
        roots(i) =r2(i)
        weights(i) = w2(i)
 5000 continue

        nquads = nqs(2)

        endif
        endif


        if (alpha .lt. dtwo**(-5)) then

        do 6000 i=1,nqs(1)
        roots(i) =r1(i)
        weights(i) = w1(i)
 6000 continue

        nquads = nqs(1)


        endif

        return
        end
c
c
c
c
c
        subroutine lap_load_quads33(alpha,roots,weights,nquads)
        implicit real *8 (a-h,o-z)
        dimension r1(29),r2(26),r3(24),r4(22),r5(20),r6(19),
     1          w1(29),w2(26),w3(24),w4(22),w5(20),w6(19),nqs(6), 
     2          roots(1),weights(1)

        data nqs/29,26,24,22,20,19/


        data r1/
     1  0.7032418384953513D-06,0.2829285548501598D-04,
     2  0.3379649985297872D-03,0.2085001120478139D-02,
     3  0.8303465965850814D-02,0.2413291936263401D-01,
     4  0.5537589752195094D-01,0.1059993688412751D+00,
     5  0.1762674276665306D+00,0.2625861423514643D+00,
     6  0.3587825858927442D+00,0.4578193348022750D+00,
     7  0.5531565795309522D+00,0.6395566353572360D+00,
     8  0.7135192936322428D+00,0.7735467580170507D+00,
     9  0.8201663830296315D+00,0.8554555641576993D+00,
     *  0.8821318728828485D+00,0.9027327239187386D+00,
     1  0.9192563434296755D+00,0.9331691834076259D+00,
     2  0.9455318586540960D+00,0.9570886958400773D+00,
     3  0.9682577212593310D+00,0.9790063739763922D+00,
     4  0.9886773558249392D+00,0.9959916054786861D+00,
     *  0.9995664784111441D+00/

        data w1/
     1  0.3342977391116124D-05,0.8316899948027126D-04,
     2  0.7093582995497250D-03,0.3283627357703364D-02,
     3  0.1003566315343365D-01,0.2263688098494124D-01,
     4  0.4054425024833109D-01,0.6074826450025148D-01,
     5  0.7916473253239641D-01,0.9241328790859566D-01,
     6  0.9878361361291027D-01,0.9818712716467377D-01,
     7  0.9160776044180724D-01,0.8060438690392765D-01,
     8  0.6707003512766781D-01,0.5308072446396760D-01,
     9  0.4052156147999625D-01,0.3052442916679944D-01,
     *  0.2325663579042693D-01,0.1827729370538438D-01,
     1  0.1501190811654960D-01,0.1298892165913265D-01,
     2  0.1186060396920357D-01,0.1132491571710651D-01,
     3  0.1101115851154211D-01,0.1037747270290434D-01,
     4  0.8743384011709448D-02,0.5631638862174912D-02,
     *  0.1513851630041533D-02/

        data r2/
     1  0.6795592980159987D-06,0.2741940191714917D-04,
     2  0.3283605001478262D-03,0.2030150499635527D-02,
     3  0.8099228439567282D-02,0.2356877234420199D-01,
     4  0.5411627789643752D-01,0.1035831057147714D+00,
     5  0.1721121980364626D+00,0.2560046783187964D+00,
     6  0.3490458445052719D+00,0.4443076729759751D+00,
     7  0.5356117838448622D+00,0.6183881909041130D+00,
     8  0.6900292796684632D+00,0.7498612108328323D+00,
     9  0.7987596047441867D+00,0.8385049770379222D+00,
     *  0.8711373441938197D+00,0.8985309425476780D+00,
     1  0.9221924394394852D+00,0.9431515081937784D+00,
     2  0.9618384146638339D+00,0.9779586823838805D+00,
     3  0.9905119891661138D+00,0.9981363651277104D+00/

        data w2/
     1  0.3232519272714925D-05,0.8067549418431407D-04,
     2  0.6899726769670822D-03,0.3201256045795155D-02,
     3  0.9801122756436148D-02,0.2213042419335724D-01,
     4  0.3963880330341996D-01,0.5932011924410929D-01,
     5  0.7709668885599239D-01,0.8962036322851644D-01,
     6  0.9528239551440150D-01,0.9419762153867901D-01,
     7  0.8764338243348396D-01,0.7747527034819145D-01,
     8  0.6570116560580692D-01,0.5412011479859139D-01,
     9  0.4398480291721554D-01,0.3585298109289637D-01,
     *  0.2972710274005985D-01,0.2531274808775575D-01,
     1  0.2218390770300022D-01,0.1980841283505714D-01,
     2  0.1751655810070185D-01,0.1455153220616199D-01,
     3  0.1031636864875354D-01,0.4742977111192790D-02/

        data r3/
     1  0.4182814867553561D-06,0.1766500196004854D-04,
     2  0.2208876507062362D-03,0.1423577998444289D-02,
     3  0.5909716379115696D-02,0.1785697456481512D-01,
     4  0.4246438000054591D-01,0.8392806761195634D-01,
     5  0.1435239823739200D+00,0.2189755920792093D+00,
     6  0.3052681603995049D+00,0.3962243737338935D+00,
     7  0.4860310157059487D+00,0.5702569683549698D+00,
     8  0.6462705605366123D+00,0.7131489251439557D+00,
     9  0.7712464905852743D+00,0.8216279534124203D+00,
     *  0.8655299594159105D+00,0.9039091545262719D+00,
     1  0.9370841194688492D+00,0.9645412191183160D+00,
     2  0.9850439350273521D+00,0.9970993526519468D+00/

        data w3/
     1  0.2009550456092772D-05,0.5275575474469191D-04,
     2  0.4734950951660828D-03,0.2302701894029774D-02,
     3  0.7380695050748821D-02,0.1741769811528557D-01,
     4  0.3253166785693303D-01,0.5061915511908395D-01,
     5  0.6817769366092662D-01,0.8186134680814866D-01,
     6  0.8966771758709342D-01,0.9125875152760109D-01,
     7  0.8760931854332237D-01,0.8040588223857050D-01,
     8  0.7147521101351965D-01,0.6235046971722166D-01,
     9  0.5403119005768533D-01,0.4694533467496034D-01,
     *  0.4102424252958071D-01,0.3579089088715464D-01,
     1  0.3047181445775774D-01,0.2422808602333025D-01,
     2  0.1651666859360827D-01,0.7405203243070753D-02/

        data r4/
     1  0.3752586772692679D-06,0.1591152474210759D-04,
     2  0.1997601642414861D-03,0.1292905158354775D-02,
     3  0.5392331047674638D-02,0.1637826579179210D-01,
     4  0.3917325130548542D-01,0.7791787419089607D-01,
     5  0.1341744486630323D+00,0.2062553440863942D+00,
     6  0.2898947537696032D+00,0.3796879887489428D+00,
     7  0.4705253753960139D+00,0.5585088686776791D+00,
     8  0.6412168014372599D+00,0.7174498069606244D+00,
     9  0.7866849563312129D+00,0.8484552821306516D+00,
     *  0.9018626014915889D+00,0.9454140443291322D+00,
     1  0.9772486573048721D+00,0.9956212521809992D+00/

        data w4/
     1  0.1804401689928533D-05,0.4758190254388718D-04,
     2  0.4290250749440842D-03,0.2097068307045790D-02,
     3  0.6760982410542565D-02,0.1606491614611792D-01,
     4  0.3024810688809023D-01,0.4751315431407851D-01,
     5  0.6470586113841752D-01,0.7871894196375400D-01,
     6  0.8763394580620002D-01,0.9108536947757093D-01,
     7  0.8993571359979812D-01,0.8563248249621595D-01,
     8  0.7959309190707323D-01,0.7279659278608054D-01,
     9  0.6560450595375676D-01,0.5778938023021689D-01,
     *  0.4877229372229259D-01,0.3801176574780658D-01,
     1  0.2535745469577838D-01,0.1119996102998560D-01/

        data r5/
     1  0.5358067277720389D-06,0.2204447981328054D-04,
     2  0.2689855754452230D-03,0.1694006455895736D-02,
     3  0.6882821479725457D-02,0.2039561400202403D-01,
     4  0.4768345843918066D-01,0.9293960980150096D-01,
     5  0.1573162694537958D+00,0.2385996293351324D+00,
     6  0.3322488073148688D+00,0.4329031690147107D+00,
     7  0.5355027987509117D+00,0.6356890025363648D+00,
     8  0.7296614216768674D+00,0.8139013067246403D+00,
     9  0.8850845579415663D+00,0.9402676223531549D+00,
     *  0.9773402461426561D+00,0.9960863370017313D+00/

        data w5/
     1  0.2559537368702886D-05,0.6527947001855047D-04,
     2  0.5702000897702765D-03,0.2702287755486031D-02,
     3  0.8454390585184749D-02,0.1951967107092109D-01,
     4  0.3578759982619349D-01,0.5492168906480530D-01,
     5  0.7344281233795624D-01,0.8833974275011612D-01,
     6  0.9804573473275231D-01,0.1024135381741017D+00,
     7  0.1020627710999519D+00,0.9768170090207731D-01,
     8  0.8967826125241760D-01,0.7824419817118794D-01,
     9  0.6362473922898089D-01,0.4637464798411251D-01,
     *  0.2769470075284041D-01,0.1037347521375696D-01/

        data r6/
     1  0.3047512610302608D-06,0.1292234125722028D-04,
     2  0.1638452257411925D-03,0.1081369928569195D-02,
     3  0.4637274076530052D-02,0.1456751267221328D-01,
     4  0.3614977967636419D-01,0.7466320583546892D-01,
     5  0.1333997565746756D+00,0.2123994108001070D+00,
     6  0.3084675083711279D+00,0.4161461401724838D+00,
     7  0.5289377087957155D+00,0.6402747878120616D+00,
     8  0.7440978626373063D+00,0.8351341814614550D+00,
     9  0.9090197460527727D+00,0.9623617355461667D+00,
     *  0.9927910546837736D+00/

        data w6/
     1  0.1463901346185156D-05,0.3870674356633344D-04,
     2  0.3543205138949689D-03,0.1779918803054088D-02,
     3  0.5958994671780512D-02,0.1480893729824568D-01,
     4  0.2927189787264361D-01,0.4832612658287175D-01,
     5  0.6915229144461992D-01,0.8829645156102656D-01,
     6  0.1029084415543303D+00,0.1113494721275843D+00,
     7  0.1131297689219676D+00,0.1085301958606257D+00,
     8  0.9823535873086417D-01,0.8310899930383962D-01,
     9  0.6409589553735395D-01,0.4219172294551090D-01,
     *  0.1846103562487389D-01/
c     
c         %%%%%%%%%%%%%%%%%%%%%

ccc        call prin2('alpha = *',alpha,1)
        dtwo = 2
        if (alpha .ge. 1/dtwo) then
        do 1000 i=1,nqs(6)
        roots(i) =r6(i)
        weights(i) = w6(i)
 1000 continue

        nquads = nqs(6)
        endif

        if (alpha .ge. dtwo**(-2)) then
        if (alpha .lt. dtwo**(-1)) then

        do 2000 i=1,nqs(5)
        roots(i) =r5(i)
        weights(i) = w5(i)
 2000 continue

        nquads = nqs(5)

        endif
        endif


        if (alpha .ge. dtwo**(-3)) then
        if (alpha .lt. dtwo**(-2)) then

        do 3000 i=1,nqs(4)
        roots(i) =r4(i)
        weights(i) = w4(i)
 3000 continue

        nquads = nqs(4)

        endif
        endif


        if (alpha .ge. dtwo**(-4)) then
        if (alpha .lt. dtwo**(-3)) then

        do 4000 i=1,nqs(3)
        roots(i) =r3(i)
        weights(i) = w3(i)
 4000 continue

        nquads = nqs(3)

        endif
        endif


        if (alpha .ge. dtwo**(-5)) then
        if (alpha .lt. dtwo**(-4)) then

        do 5000 i=1,nqs(2)
        roots(i) =r2(i)
        weights(i) = w2(i)
 5000 continue

        nquads = nqs(2)

        endif
        endif


        if (alpha .lt. dtwo**(-5)) then

        do 6000 i=1,nqs(1)
        roots(i) =r1(i)
        weights(i) = w1(i)
 6000 continue

        nquads = nqs(1)

        endif

        return
        end
c
c
c
c
c
        subroutine lap_load_quads34(alpha,roots,weights,nquads)
        implicit real *8 (a-h,o-z)
        dimension r1(28),r2(26),r3(24),r4(21),r5(20),r6(19),
     1          w1(28),w2(26),w3(24),w4(21),w5(20),w6(19),nqs(6),
     2          roots(1),weights(1)

        data nqs/28,26,24,21,20,19/


        data r1/
     1  0.6242164475498522D-06,0.2541659299381960D-04,
     2  0.3070625172062610D-03,0.1915016712958979D-02,
     3  0.7706187779909897D-02,0.2261960276052561D-01,
     4  0.5238923980459758D-01,0.1011561056654545D+00,
     5  0.1695602042970422D+00,0.2544283244538123D+00,
     6  0.3499004111763223D+00,0.4490828160846411D+00,
     7  0.5454416187235970D+00,0.6336710802719170D+00,
     8  0.7101568678083710D+00,0.7732053676854352D+00,
     9  0.8230270722394451D+00,0.8613262898448629D+00,
     *  0.8905378667383814D+00,0.9130766477868764D+00,
     1  0.9309269254899251D+00,0.9455675256461870D+00,
     2  0.9580465331529899D+00,0.9690553709478884D+00,
     3  0.9789363314334616D+00,0.9876216421660761D+00,
     4  0.9945773674977959D+00,0.9989214284435458D+00/

        data w1/
     1  0.2975131901781543D-05,0.7500861822705536D-04,
     2  0.6478960294837863D-03,0.3036271317867641D-02,
     3  0.9392021265922000D-02,0.2143486358789976D-01,
     4  0.3882938934722513D-01,0.5881553007483141D-01,
     5  0.7744243795422123D-01,0.9128754711319394D-01,
     6  0.9848391811196469D-01,0.9877786838902394D-01,
     7  0.9304823997382683D-01,0.8280443380644083D-01,
     8  0.6988462724530438D-01,0.5625428697398206D-01,
     9  0.4368698227038412D-01,0.3333070254802138D-01,
     *  0.2550164566776950D-01,0.1990733633874413D-01,
     1  0.1603808150578158D-01,0.1341642372891444D-01,
     2  0.1165669298696881D-01,0.1041920033364136D-01,
     3  0.9332797249870694D-02,0.7947373329821726D-02,
     4  0.5809402791879744D-02,0.2736046306886066D-02/

        data r2/
     1  0.5028348371873718D-06,0.2081580024330350D-04,
     2  0.2554611181239419D-03,0.1617545165456130D-02,
     3  0.6605008427901176D-02,0.1966006391884804D-01,
     4  0.4613771804022020D-01,0.9017678388598528D-01,
     5  0.1528357851872795D+00,0.2315979227685593D+00,
     6  0.3212581666166853D+00,0.4154446023882412D+00,
     7  0.5080238587371474D+00,0.5940492340117118D+00,
     8  0.6702563442480363D+00,0.7351978601358167D+00,
     9  0.7890435276943019D+00,0.8330833526055253D+00,
     *  0.8691121468300208D+00,0.8989303683573576D+00,
     1  0.9240572764087278D+00,0.9455888663880200D+00,
     2  0.9640968560604642D+00,0.9795396231175709D+00,
     3  0.9912712042180262D+00,0.9982934848049524D+00/

        data w2/
     1  0.2405263333601432D-05,0.6176471791284859D-04,
     2  0.5429577136944046D-03,0.2588737790779501D-02,
     3  0.8144516156221842D-02,0.1889752930126318D-01,
     4  0.3478165568601773D-01,0.5348096603098923D-01,
     5  0.7140063881075437D-01,0.8522882394689634D-01,
     6  0.9300549198318986D-01,0.9432948728544028D-01,
     7  0.8999660417329060D-01,0.8150785384840919D-01,
     8  0.7066848686999000D-01,0.5925158473000244D-01,
     9  0.4866567856881286D-01,0.3972165678134769D-01,
     *  0.3264056008388039D-01,0.2725116922693932D-01,
     1  0.2318691007777523D-01,0.1997326129949541D-01,
     2  0.1703605908476707D-01,0.1373775836051527D-01,
     3  0.9549710064783993D-02,0.4347732143497321D-02/

        data r3/
     1  0.4926268384948025D-06,0.2037497360959796D-04,
     2  0.2497860647223685D-03,0.1579645979939736D-02,
     3  0.6440983730657823D-02,0.1914041715535714D-01,
     4  0.4483578221478671D-01,0.8745726525850785D-01,
     5  0.1479201841322448D+00,0.2237055909163651D+00,
     6  0.3098019207130455D+00,0.4002515204640890D+00,
     7  0.4895428821848296D+00,0.5734836333182704D+00,
     8  0.6495313346960454D+00,0.7166971289408657D+00,
     9  0.7751681516492600D+00,0.8258097394356335D+00,
     *  0.8696829133474737D+00,0.9076381099790301D+00,
     1  0.9400023324246539D+00,0.9664134314342355D+00,
     2  0.9859034345302825D+00,0.9972736948083496D+00/

        data w3/
     1  0.2356029026039860D-05,0.6043685525253762D-04,
     2  0.5306087279677410D-03,0.2525968290845416D-02,
     3  0.7932351135008793D-02,0.1836470512542585D-01,
     4  0.3371377321692695D-01,0.5168850635775588D-01,
     5  0.6879965230831538D-01,0.8190815925860763D-01,
     6  0.8926753451297835D-01,0.9069857801969165D-01,
     7  0.8718021684993778D-01,0.8027949116552344D-01,
     8  0.7165779782669858D-01,0.6271362906234364D-01,
     9  0.5438081049442996D-01,0.4708597116636944D-01,
     *  0.4080744763893999D-01,0.3516050403055029D-01,
     1  0.2950635392904269D-01,0.2314795267678885D-01,
     2  0.1562247248270200D-01,0.6964722838871122D-02/

        data r4/
     1  0.5981525998262264D-06,0.2439920749627709D-04,
     2  0.2952109642497742D-03,0.1843258408721818D-02,
     3  0.7423310556215593D-02,0.2179654107826009D-01,
     4  0.5047237071163602D-01,0.9737930380808145D-01,
     5  0.1630285795690960D+00,0.2442967220710066D+00,
     6  0.3356854583513610D+00,0.4311195544928092D+00,
     7  0.5253987720298246D+00,0.6149017712884013D+00,
     8  0.6975709348055431D+00,0.7724149624622769D+00,
     9  0.8388096953174590D+00,0.8958761142241487D+00,
     *  0.9421895306932650D+00,0.9759331894331705D+00,
     *  0.9953709734345841D+00/

        data w4/
     1  0.2852085162124632D-05,0.7204619915774444D-04,
     2  0.6232840124917209D-03,0.2924329211839750D-02,
     3  0.9051552807093039D-02,0.2065649934686711D-01,
     4  0.3738356322219962D-01,0.5651851147618602D-01,
     5  0.7423582060214543D-01,0.8735185465116804D-01,
     6  0.9438797052461761D-01,0.9560195798018961D-01,
     7  0.9235380700837162D-01,0.8631822263697045D-01,
     8  0.7886720014788588D-01,0.7072785565045405D-01,
     9  0.6192028681607461D-01,0.5197218400727688D-01,
     *  0.4034328733625933D-01,0.2684508241222434D-01,
     *  0.1184183186536503D-01/

        data r5/
     1  0.5091247236209373D-06,0.2102696842631224D-04,
     2  0.2575043425004592D-03,0.1627390401679245D-02,
     3  0.6634605813890598D-02,0.1972498251242953D-01,
     4  0.4626424780759787D-01,0.9045682621301165D-01,
     5  0.1535755644897013D+00,0.2335780124545400D+00,
     6  0.3260571486041043D+00,0.4256943038966468D+00,
     7  0.5273949478449969D+00,0.6267643861409716D+00,
     8  0.7200564262784382D+00,0.8039558576753729D+00,
     9  0.8754964476925046D+00,0.9321918780446502D+00,
     *  0.9722580413628686D+00,0.9947249577529150D+00/

        data w5/
     1  0.2434097691101294D-05,0.6234442602089487D-04,
     2  0.5467799372652296D-03,0.2601624540571841D-02,
     3  0.8171651351726711D-02,0.1894205338636944D-01,
     4  0.3486965741040833D-01,0.5373369260924928D-01,
     5  0.7214516922599570D-01,0.8709887094887543D-01,
     6  0.9695323207733603D-01,0.1014644909334786D+00,
     7  0.1012061884057018D+00,0.9691274013032477D-01,
     8  0.8912033458272093D-01,0.7818162915601604D-01,
     9  0.6448093185582324D-01,0.4861030314920035D-01,
     *  0.3136270135043472D-01,0.1353317042478959D-01/

        data r6/
     1  0.3050479742268840D-06,0.1294819580306244D-04,
     2  0.1642268660696422D-03,0.1083691216119113D-02,
     3  0.4645202070507756D-02,0.1458591781676261D-01,
     4  0.3618422352636602D-01,0.7472304669471874D-01,
     5  0.1335003397076612D+00,0.2125555705236757D+00,
     6  0.3086811884133127D+00,0.4164001114055139D+00,
     7  0.5292013074144849D+00,0.6405158907130246D+00,
     8  0.7442926169799407D+00,0.8352711984985727D+00,
     9  0.9090999748483972D+00,0.9623961878059640D+00,
     *  0.9927977834449693D+00/

        data w6/
     1  0.1465760554214548D-05,0.3879262726672116D-04,
     2  0.3551514747409075D-03,0.1783355376820621D-02,
     3  0.5966974347140952D-02,0.1482193979108523D-01,
     4  0.2929160705762124D-01,0.4835838838228262D-01,
     5  0.6920155670766008D-01,0.8835606898187747D-01,
     6  0.1029604330844681D+00,0.1113757152289880D+00,
     7  0.1131225918341468D+00,0.1084938636234306D+00,
     8  0.9818114431181203D-01,0.8304982873859011D-01,
     9  0.6404317704878737D-01,0.4215406411735346D-01,
     *  0.1844388150537355D-01/
c        

c     
c         %%%%%%%%%%%%%%%%%%%%%

ccc        call prin2('alpha = *',alpha,1)
        dtwo = 2
        if (alpha .ge. 1/dtwo) then
        do 1000 i=1,nqs(6)
        roots(i) =r6(i)
        weights(i) = w6(i)
 1000 continue

        nquads = nqs(6)
        endif

        if (alpha .ge. dtwo**(-2)) then
        if (alpha .lt. dtwo**(-1)) then

        do 2000 i=1,nqs(5)
        roots(i) =r5(i)
        weights(i) = w5(i)
 2000 continue

        nquads = nqs(5)

        endif
        endif


        if (alpha .ge. dtwo**(-3)) then
        if (alpha .lt. dtwo**(-2)) then

        do 3000 i=1,nqs(4)
        roots(i) =r4(i)
        weights(i) = w4(i)
 3000 continue

        nquads = nqs(4)

        endif
        endif


        if (alpha .ge. dtwo**(-4)) then
        if (alpha .lt. dtwo**(-3)) then

        do 4000 i=1,nqs(3)
        roots(i) =r3(i)
        weights(i) = w3(i)
 4000 continue

        nquads = nqs(3)

        endif
        endif


        if (alpha .ge. dtwo**(-5)) then
        if (alpha .lt. dtwo**(-4)) then

        do 5000 i=1,nqs(2)
        roots(i) =r2(i)
        weights(i) = w2(i)
 5000 continue

        nquads = nqs(2)

        endif
        endif


        if (alpha .lt. dtwo**(-5)) then

        do 6000 i=1,nqs(1)
        roots(i) =r1(i)
        weights(i) = w1(i)
 6000 continue

        nquads = nqs(1)

        endif

        return
        end
c
c
c
c
c
        subroutine lap_load_quads35(alpha,roots,weights,nquads)
        implicit real *8 (a-h,o-z)
        dimension r1(27),r2(25),r3(23),r4(21),r5(19),r6(19),
     1          w1(27),w2(25),w3(23),w4(21),w5(19),w6(19),nqs(6),
     2          roots(1),weights(1)

        data nqs/27,25,23,21,19,19/


        data r1/
     1  0.6236706054096086D-06,0.2539575540950875D-04,
     2  0.3068311918457391D-03,0.1913733299030436D-02,
     3  0.7701843070007875D-02,0.2260998973476841D-01,
     4  0.5237658869320579D-01,0.1011563553202354D+00,
     5  0.1696156188742719D+00,0.2546208691065000D+00,
     6  0.3503636167029847D+00,0.4500103015507614D+00,
     7  0.5470887312107913D+00,0.6363410626498169D+00,
     8  0.7141576595874583D+00,0.7787701997433205D+00,
     9  0.8302261990679165D+00,0.8700328625906538D+00,
     *  0.9004790173743373D+00,0.9239201780108630D+00,
     1  0.9423361934921894D+00,0.9572041030669362D+00,
     2  0.9695321361872945D+00,0.9799074828192653D+00,
     3  0.9885000760858900D+00,0.9950597409440913D+00,
     *  0.9990284648152648D+00/

        data w1/
     1  0.2972561589085047D-05,0.7494864373711389D-04,
     2  0.6474306366591981D-03,0.3034415481137171D-02,
     3  0.9387657155176821D-02,0.2142934689133875D-01,
     4  0.3883105549724525D-01,0.5884404183903838D-01,
     5  0.7753090012015530D-01,0.9148199164715779D-01,
     6  0.9884091202202366D-01,0.9935996272235638D-01,
     7  0.9391387902495819D-01,0.8398600824595924D-01,
     8  0.7135249502830195D-01,0.5788635698113880D-01,
     9  0.4528896544761814D-01,0.3471864600440060D-01,
     *  0.2657366007378172D-01,0.2064028621490231D-01,
     1  0.1643697740312194D-01,0.1346595046092205D-01,
     2  0.1128800344077923D-01,0.9491379535799530D-02,
     3  0.7650606355002702D-02,0.5369894827422510D-02,
     *  0.2471255738276175D-02/

        data r2/
     1  0.6134663355741229D-06,0.2499168147535716D-04,
     2  0.3020213891017706D-03,0.1883726951872864D-02,
     3  0.7578878617731325D-02,0.2223503101757844D-01,
     4  0.5145587498278700D-01,0.9923686859533831D-01,
     5  0.1660943613939761D+00,0.2488005272529758D+00,
     6  0.3415638375305889D+00,0.4377326171001018D+00,
     7  0.5312134519302661D+00,0.6173119421913015D+00,
     8  0.6930855419204503D+00,0.7573645762687532D+00,
     9  0.8104954627708507D+00,0.8538403884739547D+00,
     *  0.8891855409571822D+00,0.9182600272049922D+00,
     1  0.9424416728538076D+00,0.9625885783777409D+00,
     2  0.9789341350923977D+00,0.9910881186788929D+00,
     *  0.9982657966608404D+00/

        data w2/
     1  0.2924274527602958D-05,0.7376487524715006D-04,
     2  0.6373191034215725D-03,0.2986575555810059D-02,
     3  0.9234359246874420D-02,0.2105587514928788D-01,
     4  0.3808590667961923D-01,0.5756587070622608D-01,
     5  0.7558885728976018D-01,0.8882774785459368D-01,
     6  0.9556600846210245D-01,0.9574233388986174D-01,
     7  0.9043150284934459D-01,0.8127478259412084D-01,
     8  0.7008174833468976D-01,0.5854279522810660D-01,
     9  0.4795645422602529D-01,0.3903929207996255D-01,
     *  0.3194388409970347D-01,0.2643733791178548D-01,
     1  0.2207082743704643D-01,0.1826223445808236D-01,
     2  0.1436018200409814D-01,0.9808295130248575D-02,
     *  0.4423120559453865D-02/

        data r3/
     1  0.6086046794338888D-06,0.2479489151223875D-04,
     2  0.2996309577982905D-03,0.1868530416551083D-02,
     3  0.7515478594295574D-02,0.2203817754344976D-01,
     4  0.5096304451987315D-01,0.9818784593238070D-01,
     5  0.1641300238691829D+00,0.2454989528359578D+00,
     6  0.3365349862098500D+00,0.4307826330233057D+00,
     7  0.5225582075227375D+00,0.6077682486074126D+00,
     8  0.6841221238680597D+00,0.7509129382004033D+00,
     9  0.8085433142389842D+00,0.8579594007544399D+00,
     *  0.9001167300339106D+00,0.9355478353250314D+00,
     1  0.9641049413161609D+00,0.9849876375678437D+00,
     *  0.9971022134414111D+00/

        data w3/
     1  0.2901157123514683D-05,0.7318421790650061D-04,
     2  0.6322470935745088D-03,0.2962070144785088D-02,
     3  0.9154192495325850D-02,0.2085614402389180D-01,
     4  0.3767709889616064D-01,0.5684519963502043D-01,
     5  0.7446743501424729D-01,0.8727877520574838D-01,
     6  0.9369217412733996D-01,0.9384242656896858D-01,
     7  0.8902339653584099D-01,0.8101897680412503D-01,
     8  0.7157926591942715D-01,0.6207920812561050D-01,
     9  0.5334827173057864D-01,0.4565051143697184D-01,
     *  0.3875750272387946D-01,0.3207979780528277D-01,
     1  0.2489467674152409D-01,0.1667837644777265D-01,
     *  0.7406167148894356D-02/

        data r4/
     1  0.5676719827660570D-06,0.2325534773937706D-04,
     2  0.2825335014630720D-03,0.1771269280525644D-02,
     3  0.7162079254656309D-02,0.2111341472142431D-01,
     4  0.4908429758688361D-01,0.9507377682070905D-01,
     5  0.1597867170381715D+00,0.2403438814958307D+00,
     6  0.3314407636711027D+00,0.4270717292261295D+00,
     7  0.5219666837815173D+00,0.6123346320927983D+00,
     8  0.6959146484590007D+00,0.7715403717056837D+00,
     9  0.8384827059695213D+00,0.8958439539229843D+00,
     *  0.9422534257395817D+00,0.9759838013094874D+00,
     *  0.9953832614410487D+00/

        data w4/
     1  0.2709270017319160D-05,0.6876575510649384D-04,
     2  0.5976708748126353D-03,0.2817260257358292D-02,
     3  0.8761647531599988D-02,0.2009289055589209D-01,
     4  0.3654898392488084D-01,0.5555179383845775D-01,
     5  0.7337240634771652D-01,0.8682691200867186D-01,
     6  0.9434104787797080D-01,0.9603023755740490D-01,
     7  0.9312754086250312D-01,0.8723943703980740D-01,
     8  0.7973780883379265D-01,0.7140286448674382D-01,
     9  0.6233760456312083D-01,0.5215420320739969D-01,
     *  0.4036900952218206D-01,0.2680728205416655D-01,
     *  0.1181192363039439D-01/

        data r5/
     1  0.6452004055402418D-06,0.2618095271679727D-04,
     2  0.3153288850074635D-03,0.1961166880681332D-02,
     3  0.7873016325527964D-02,0.2306465135764351D-01,
     4  0.5335092159528464D-01,0.1029740928255022D+00,
     5  0.1727720944136453D+00,0.2599855375870241D+00,
     6  0.3594718191610763D+00,0.4653123446063190D+00,
     7  0.5719420173674586D+00,0.6745290660608751D+00,
     8  0.7688492317057151D+00,0.8510754211779031D+00,
     9  0.9177497010370273D+00,0.9659415490098674D+00,
     *  0.9934725501233968D+00/

        data w5/
     1  0.3072771667740401D-05,0.7717984842512433D-04,
     2  0.6644205387248934D-03,0.3104370570653860D-02,
     3  0.9578124413523232D-02,0.2181853619209749D-01,
     4  0.3949454365935738D-01,0.5989379749220693D-01,
     5  0.7921922273425042D-01,0.9432152694014570D-01,
     6  0.1036451909719247D+00,0.1071027451799092D+00,
     7  0.1053545550623865D+00,0.9911968763425756D-01,
     8  0.8888319384282950D-01,0.7498985425008483D-01,
     9  0.5786638230299419D-01,0.3815042321440690D-01,
     *  0.1671317238015381D-01/

        data r6/
     1  0.4172242213043255D-06,0.1746524582165125D-04,
     2  0.2168403031651781D-03,0.1390448282276475D-02,
     3  0.5758792310416028D-02,0.1742245202117026D-01,
     4  0.4166392043989175D-01,0.8321827871459177D-01,
     5  0.1445447661837109D+00,0.2250283114139686D+00,
     6  0.3212619438170758D+00,0.4279675223675173D+00,
     7  0.5390176447543820D+00,0.6482318014333793D+00,
     8  0.7498699787599515D+00,0.8388952708573182D+00,
     9  0.9111087155237643D+00,0.9632289792890075D+00,
     *  0.9929574352755507D+00/

        data w6/
     1  0.2000415454611491D-05,0.5201945248159443D-04,
     2  0.4634796123046956D-03,0.2243804796061120D-02,
     3  0.7188247322043841D-02,0.1705131533492570D-01,
     4  0.3225212320824602D-01,0.5127999757883973D-01,
     5  0.7126446694359214D-01,0.8911805071067257D-01,
     6  0.1024517088354028D+00,0.1099248126317872D+00,
     7  0.1111373424778982D+00,0.1063301949127994D+00,
     8  0.9610460070871343D-01,0.8124490444602325D-01,
     9  0.6263412091851150D-01,0.4122169426789763D-01,
     *  0.1803511542634468D-01/
c        

c     
c         %%%%%%%%%%%%%%%%%%%%%

ccc        call prin2('alpha = *',alpha,1)
        dtwo = 2
        if (alpha .ge. 1/dtwo) then
        do 1000 i=1,nqs(6)
        roots(i) =r6(i)
        weights(i) = w6(i)
 1000 continue

        nquads = nqs(6)
        endif

        if (alpha .ge. dtwo**(-2)) then
        if (alpha .lt. dtwo**(-1)) then

        do 2000 i=1,nqs(5)
        roots(i) =r5(i)
        weights(i) = w5(i)
 2000 continue

        nquads = nqs(5)

        endif
        endif


        if (alpha .ge. dtwo**(-3)) then
        if (alpha .lt. dtwo**(-2)) then

        do 3000 i=1,nqs(4)
        roots(i) =r4(i)
        weights(i) = w4(i)
 3000 continue

        nquads = nqs(4)

        endif
        endif


        if (alpha .ge. dtwo**(-4)) then
        if (alpha .lt. dtwo**(-3)) then

        do 4000 i=1,nqs(3)
        roots(i) =r3(i)
        weights(i) = w3(i)
 4000 continue

        nquads = nqs(3)

        endif
        endif


        if (alpha .ge. dtwo**(-5)) then
        if (alpha .lt. dtwo**(-4)) then

        do 5000 i=1,nqs(2)
        roots(i) =r2(i)
        weights(i) = w2(i)
 5000 continue

        nquads = nqs(2)

        endif
        endif


        if (alpha .lt. dtwo**(-5)) then

        do 6000 i=1,nqs(1)
        roots(i) =r1(i)
        weights(i) = w1(i)
 6000 continue

        nquads = nqs(1)

        endif

        return
        end
c
c
c
c
c
        subroutine lap_load_quads36(alpha,roots,weights,nquads)
        implicit real *8 (a-h,o-z)
        dimension r1(27),r2(25),r3(23),r4(21),r5(19),r6(19),
     1          w1(27),w2(25),w3(23),w4(21),w5(19),w6(19),nqs(6),
     2          roots(1),weights(1)

        data nqs/27,25,23,21,19,19/


        data r1/
     1  0.5158702187856276D-06,0.2129117563616083D-04,
     2  0.2605568838687541D-03,0.1645412840814688D-02,
     3  0.6702323314964060D-02,0.1990714813040191D-01,
     4  0.4663960845462937D-01,0.9106411885156981D-01,
     5  0.1543087151545506D+00,0.2340136993687460D+00,
     6  0.3252141619929766D+00,0.4217839872888255D+00,
     7  0.5177072013368467D+00,0.6078819221201593D+00,
     8  0.6885288941481606D+00,0.7573816194613192D+00,
     9  0.8137394312902870D+00,0.8583235698691994D+00,
     *  0.8928582717345634D+00,0.9194705220265045D+00,
     1  0.9401647386129730D+00,0.9565465744219308D+00,
     2  0.9697600615446198D+00,0.9805093324408038D+00,
     3  0.9890834328468541D+00,0.9953945904028977D+00,
     *  0.9991047763382126D+00/

        data w1/
     1  0.2465993852112034D-05,0.6311333763177062D-04,
     2  0.5530708736160492D-03,0.2629098764559175D-02,
     3  0.8248981086513141D-02,0.1909659275025670D-01,
     4  0.3509566308846284D-01,0.5394856898357389D-01,
     5  0.7213042089996119D-01,0.8642608883154698D-01,
     6  0.9493198242132325D-01,0.9719190125471584D-01,
     7  0.9379790985401001D-01,0.8591870001971044D-01,
     8  0.7499641453093607D-01,0.6259436908274577D-01,
     9  0.5024865990723382D-01,0.3921774359535119D-01,
     *  0.3021538860731549D-01,0.2334762577994536D-01,
     1  0.1830909480477888D-01,0.1464468063680352D-01,
     2  0.1189926923141299D-01,0.9647026518724046D-02,
     3  0.7484601323742704D-02,0.5076812996971625D-02,
     *  0.2283754824305153D-02/

        data r2/
     1  0.5516509595002005D-06,0.2265455308226363D-04,
     2  0.2759048707579211D-03,0.1733997067544325D-02,
     3  0.7029237917414139D-02,0.2077672838023507D-01,
     4  0.4843571541307477D-01,0.9408866204836211D-01,
     5  0.1585899282881519D+00,0.2391797357689647D+00,
     6  0.3304898451771029D+00,0.4261170014419296D+00,
     7  0.5199999838282585D+00,0.6072938107055251D+00,
     8  0.6847918174119931D+00,0.7510217117096977D+00,
     9  0.8060571504452502D+00,0.8510663420314899D+00,
     *  0.8877324029232919D+00,0.9177507389920717D+00,
     1  0.9425075762272535D+00,0.9629013080807609D+00,
     2  0.9792433756150984D+00,0.9912613805560927D+00,
     *  0.9983042983916125D+00/

        data w2/
     1  0.2634190657602083D-05,0.6704421588609214D-04,
     2  0.5843239680499963D-03,0.2762320795535029D-02,
     3  0.8617348860423546D-02,0.1982873762037783D-01,
     4  0.3620347349556154D-01,0.5525081456453134D-01,
     5  0.7327264427675512D-01,0.8698406431164390D-01,
     6  0.9454533407429008D-01,0.9568446111636200D-01,
     7  0.9126734435900831D-01,0.8278432519535702D-01,
     8  0.7196941840702894D-01,0.6051051166482002D-01,
     9  0.4976531478076016D-01,0.4054222802617943D-01,
     *  0.3307766425315249D-01,0.2719282168209455D-01,
     1  0.2247229139052974D-01,0.1836733301943000D-01,
     2  0.1426813667558887D-01,0.9651577185452774D-02,
     *  0.4327831870523603D-02/

        data r3/
     1  0.5576466380534680D-06,0.2288083380136712D-04,
     2  0.2783797989085501D-03,0.1747473931691846D-02,
     3  0.7073890272103448D-02,0.2087379065389939D-01,
     4  0.4856657622909505D-01,0.9412983843244230D-01,
     5  0.1582594852221332D+00,0.2380418805287989D+00,
     6  0.3280501496968471D+00,0.4220098448695105D+00,
     7  0.5142313443579502D+00,0.6004604076877469D+00,
     8  0.6781633103519152D+00,0.7463870520084226D+00,
     9  0.8053415892321248D+00,0.8558555613835236D+00,
     *  0.8988391433590539D+00,0.9348383446805828D+00,
     1  0.9637567750101474D+00,0.9848537623615505D+00,
     *  0.9970774076821149D+00/

        data w3/
     1  0.2662364091722838D-05,0.6769253871727234D-04,
     2  0.5892616538515033D-03,0.2781532337803365D-02,
     3  0.8661234563016849D-02,0.1988351492336428D-01,
     4  0.3619865291766404D-01,0.5504884889630864D-01,
     5  0.7270414023321225D-01,0.8592634215372338D-01,
     6  0.9301981859394221D-01,0.9393728366382687D-01,
     7  0.8979230102871704D-01,0.8224620255778311D-01,
     8  0.7300810466670759D-01,0.6348378452523168D-01,
     9  0.5457252581714063D-01,0.4661469290471048D-01,
     *  0.3944840889266713D-01,0.3253417204321755D-01,
     1  0.2517340155096155D-01,0.1683529376720785D-01,
     *  0.7470127406133010D-02/

        data r4/
     1  0.5464918943396281D-06,0.2245423021978383D-04,
     2  0.2735812958203632D-03,0.1719970341645526D-02,
     3  0.6973999519019273D-02,0.2061574186196628D-01,
     4  0.4805928117570460D-01,0.9334454306236810D-01,
     5  0.1573103667477970D+00,0.2372566558736750D+00,
     6  0.3280304199940499D+00,0.4236934199784624D+00,
     7  0.5189431952922299D+00,0.6098813039095951D+00,
     8  0.6941101771331873D+00,0.7703418264555180D+00,
     9  0.8377690707109470D+00,0.8954681613369503D+00,
     *  0.9420840180219585D+00,0.9759248029406740D+00,
     *  0.9953731362928309D+00/

        data w4/
     1  0.2609864477371484D-05,0.6646187846069255D-04,
     2  0.5795072758854350D-03,0.2740469197972421D-02,
     3  0.8550851515378246D-02,0.1967600493237268D-01,
     4  0.3591784923389645D-01,0.5479770943476205D-01,
     5  0.7266343745815044D-01,0.8634038569704884D-01,
     6  0.9419282710673536D-01,0.9623666932797736D-01,
     7  0.9361201742394286D-01,0.8787173750985221D-01,
     8  0.8038246794437257D-01,0.7195690327208244D-01,
     9  0.6274923427783450D-01,0.5242125787834903D-01,
     *  0.4052100872038806D-01,0.2688203462103172D-01,
     *  0.1183855542902926D-01/

        data r5/
     1  0.6225185306558655D-06,0.2533753723444640D-04,
     2  0.3060592081107493D-03,0.1908966160580155D-02,
     3  0.7685327397305201D-02,0.2257932107990898D-01,
     4  0.5237927488710752D-01,0.1013927996816380D+00,
     5  0.1706099166803575D+00,0.2574446999060856D+00,
     6  0.3568646480268007D+00,0.4629486635708267D+00,
     7  0.5700312745097634D+00,0.6731432898952539D+00,
     8  0.7679468975758422D+00,0.8505535570444619D+00,
     9  0.9174914701006673D+00,0.9658440820423050D+00,
     *  0.9934549281679385D+00/

        data w5/
     1  0.2966720584906710D-05,0.7476833320804094D-04,
     2  0.6457686217222776D-03,0.3027181955020276D-02,
     3  0.9371708980770475D-02,0.2142429643522490D-01,
     4  0.3892711074893824D-01,0.5926848914803659D-01,
     5  0.7871323286638323D-01,0.9408944554643086D-01,
     6  0.1037438921683063D+00,0.1074731503375283D+00,
     7  0.1058660211280877D+00,0.9963876792316570D-01,
     8  0.8932088281512298D-01,0.7531116142634585D-01,
     9  0.5807501375357168D-01,0.3826720289710081D-01,
     *  0.1675893819445085D-01/

        data r6/
     1  0.4149377224174119D-06,0.1739066922182721D-04,
     2  0.2161573127493334D-03,0.1387503719300860D-02,
     3  0.5751881652863222D-02,0.1741490175522915D-01,
     4  0.4167013900412724D-01,0.8326262493930829D-01,
     5  0.1446501681027031D+00,0.2252024974083424D+00,
     6  0.3214919573814838D+00,0.4282252791973384D+00,
     7  0.5392709082073906D+00,0.6484536290890723D+00,
     8  0.7500435429998475D+00,0.8390146652356957D+00,
     9  0.9111775689850243D+00,0.9632582642987595D+00,
     *  0.9929631262219172D+00/

        data w6/
     1  0.1989982273200015D-05,0.5181761142230170D-04,
     2  0.4622537277335482D-03,0.2240416470549373D-02,
     3  0.7184507528958765D-02,0.1705575694527823D-01,
     4  0.3227708239626116D-01,0.5133110271477799D-01,
     5  0.7133277726684607D-01,0.8918365361400768D-01,
     6  0.1024951477372621D+00,0.1099362055213037D+00,
     7  0.1111179320076974D+00,0.1062884887618732D+00,
     8  0.9605163510320406D-01,0.8119113109149024D-01,
     9  0.6258803772350004D-01,0.4118947349350227D-01,
     *  0.1802059030205865D-01/
c        

c     
c         %%%%%%%%%%%%%%%%%%%%%

ccc        call prin2('alpha = *',alpha,1)
        dtwo = 2
        if (alpha .ge. 1/dtwo) then
        do 1000 i=1,nqs(6)
        roots(i) =r6(i)
        weights(i) = w6(i)
 1000 continue

        nquads = nqs(6)
        endif

        if (alpha .ge. dtwo**(-2)) then
        if (alpha .lt. dtwo**(-1)) then

        do 2000 i=1,nqs(5)
        roots(i) =r5(i)
        weights(i) = w5(i)
 2000 continue

        nquads = nqs(5)

        endif
        endif


        if (alpha .ge. dtwo**(-3)) then
        if (alpha .lt. dtwo**(-2)) then

        do 3000 i=1,nqs(4)
        roots(i) =r4(i)
        weights(i) = w4(i)
 3000 continue

        nquads = nqs(4)

        endif
        endif


        if (alpha .ge. dtwo**(-4)) then
        if (alpha .lt. dtwo**(-3)) then

        do 4000 i=1,nqs(3)
        roots(i) =r3(i)
        weights(i) = w3(i)
 4000 continue

        nquads = nqs(3)

        endif
        endif


        if (alpha .ge. dtwo**(-5)) then
        if (alpha .lt. dtwo**(-4)) then

        do 5000 i=1,nqs(2)
        roots(i) =r2(i)
        weights(i) = w2(i)
 5000 continue

        nquads = nqs(2)

        endif
        endif


        if (alpha .lt. dtwo**(-5)) then

        do 6000 i=1,nqs(1)
        roots(i) =r1(i)
        weights(i) = w1(i)
 6000 continue

        nquads = nqs(1)

        endif

        return
        end
c
c
c
c
c
        subroutine lapcornmat(alpha,x0s,w0s,u,v,nr0,coefs,
     1     lcoefs,akern)
        implicit real *8 (a-h,o-z)
        real *8, allocatable :: rquads(:),wquads(:),work(:)
        real *8 x0s(nr0),w0s(nr0),u(nr0,nr0),v(nr0,nr0)
        real *8 akern(nr0,nr0),coefs(lcoefs)

        allocate(work(20000),rquads(1000),wquads(1000))

        iw = 1
        iw2 = iw+ nr0*nr0*3+20
        do ipt=1,nr0
          call lap_load_quads(alpha,rquads,wquads,nquads,ipt)
          call lapkern(x0s(ipt),w0s(ipt),nr0,alpha,rquads,
     1      wquads,nquads,v,
     2      work(iw),work(iw2),akern,ipt,coefs,lcoefs)
        enddo

        return
        end
c
c
c
c
c
        subroutine lapkern(x0,w0,nr0,ang,rq,wq,nq,umat,
     1      xmat,work,amat,ipt,coefs,lcoefs)
        implicit real *8 (a-h,o-z)
        real *8 pi
        dimension coefs(lcoefs),rq(*),wq(*),umat(nr0,nr0),
     1      xmat(nq,nr0),work(1),amat(nr0,nr0)


c                           Input parameters
c
c       x0: the target point
c       w0: the weight associated with the target point x0
c       nr0: the number of discretization nodes
c       rq: the quadrature nodes
c       wq: the quadrature weights
c       umat: the val to coef interp matrix for the discretization nodes
c
c                           Output parameters
c
c       amat: the row in the kernel matrix corresponding to the point x0


        done = 1
        pi = atan(done)*4

        do 1000 i=1,nq
    
        top = -x0*sin(ang*pi)
        bottom = (cos(ang*pi)*x0-rq(i))**2+(x0*sin(ang*pi))**2
ccc        bottom = x0**2+rq(i)**2-2*rq(i)*x0*cos(ang)
        work(i) = top/bottom*wq(i)/2.0d0/pi
 1000 continue

        do 2000 i=1,nq

        do 1500 j=1,nr0

        call lapnestev2(ier,coefs,rq(i),j,xmat(i,j))
    
 1500 continue
 2000 continue



        do 3000 i=1,nr0

        val = 0

        do 2500 j=1,nq

        val = val+work(j)*xmat(j,i)

 2500 continue


        amat(ipt,i) = val

 3000 continue

        do 4000 i=1,nr0

        val = 0

        do 3500 j=1,nr0 

        val=val+amat(ipt,j)*umat(j,i)

 3500 continue

        work(i)=val

 4000 continue

        do 4100 i=1,nr0

        amat(ipt,i) = sqrt(w0)*work(i)

 4100 continue

        return
        end

c
c
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c             .  .  . various plotting and debugging
c                     routines
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccci
c
c
c
c
c
        subroutine plotfun(iw,a,b,sc,fun,pars,ifun,npts)
        implicit real *8 (a-h,o-z)
        dimension xs(20 000), ys(20 000)
        external fun


        dten = 10

        call mach_zero(eps)

ccc        hs = (b**(1/sc) - a**(1/sc))/npts
ccc        s = a**(1/sc)

        hs = (log(b)/log(dten) - log(a)/log(dten))/npts
        s = log(a)/log(dten)

        do 1000 i=1,npts

        xs(i) = s
        t = 10**s
ccc        t = s**sc
        call fun(t,pars,ifun,ys(i))
        ys(i) = log(ys(i) )/log(dten)
        s = s+hs

 1000 continue        

        itype = 3
        call quagraph(iw,xs,ys,npts,itype,'fun *')

        return
        end
c
c
c
c
c
        subroutine testinterps(umat,vmat,coefs,rs,ws,n)
        implicit real *8 (a-h,o-z)
        dimension umat(n,n),vmat(n,n),coefs(1),ws(1),rs(1),rhs(1000),
     1      cfs(1000)

c
c         .  .  .  testing matrix inverse
c
        errm = 0
        do 3000 i=1,n
        do 2000 j=1,n

        val = 0
        do 1000 k=1,n
        val = val+umat(i,k)*vmat(k,j)
 1000 continue

        if (i .eq. j) val=val-1
        if (abs(val) .gt. errm) errm = abs(val)

 2000 continue
 3000 continue

        call prin2('max error in inverse interp mat = *',errm,1)

        pow = 0.0d0
        do 4000 i=1,n

        rhs(i) = rs(i)**pow*sqrt(ws(i))

 4000 continue

        call prin2('rhs = *',rhs,n)
        
        do 5000 i=1,n

        val = 0
    
        do 4500 j=1,n

        val = val+vmat(i,j)*rhs(j)

 4500 continue

        cfs(i) = val
        
 5000 continue

        call prin2('coefs=*',cfs,n)

        nxs = 100
        a = log(1.0d-12)
        b = log(0.5d0)
        hv = (b-a)/(nxs-1)

        errm = 0

        do 6000 i=1,nxs

        xval = exp(a+hv*(i-1))
        val = 0

        do 5500 j=1,n

        call lapnestev(ier,coefs,xval,j,vt)
        if (ier .ne. 0) then
            call prinf('Bombing in check interp nestev*',0,0)
            call prinf('ier = *',ier,1)
            stop
        endif

        val = val+vt*cfs(j)

 5500 continue


        vexact = xval**pow

        call prin2('xval = *',xval,1)
        call prin2('error = *',val-vexact,1)

        if (abs(vexact-val) .gt. errm) errm=abs(vexact-val)

 6000 continue

        call prin2('max error in interpolation=*',errm,1)

        return
        end
c
c
c
c
c
        subroutine plotgauspt(gauspt,npts,iplot)
        implicit real *8 (a-h,o-z)
        dimension gauspt(6,1000), xs(10 000),ys(10 000)

        do 1000 i=1,npts

        xs(i) = gauspt(1,i)
        ys(i) = gauspt(2,i)

 1000 continue

        itype = 2 
        call quaplot(iplot,xs,ys,npts,itype,'point *')
        
        return
        end
c
c
c
c
c
        subroutine plotpolydisc(gauspt,ngau,cornpt,ncor,iplot)
        implicit real *8 (a-h,o-z)
        dimension gauspt(6,10 000), cornpt(6,10 000), xgs(10 000),
     1      ygs(10 000),xcs(10 000),ycs(10 000)

        do 1000 i=1,ngau

        xgs(i) = gauspt(1,i)
        ygs(i) = gauspt(2,i)

 1000 continue


        do 2000 i=1,ncor

        xcs(i) = cornpt(1,i)
        ycs(i) = cornpt(2,i)

 2000 continue

        itype = 2 
        call quaplot2(iplot,xgs,ygs,ngau,itype,xcs,ycs,ncor,
     1      itype,'point *')

        return
        end
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c            .  .  .  corners...
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
        subroutine tankcorn(cx,cy,ncs)
        implicit real *8 (a-h,o-z)
        dimension cx(1),cy(1),ctx(40),cty(40)

        ncs = 40

        data ctx/
     *   0,
     *   0.325821,
     *   0.8515777,
     *   1.244044,
     *   1.236639,
     *   1.7179655,
     *   1.7327755,
     *   2.0215715,
     *   2.4584679,
     *   2.8361241,
     *   2.7472638,
     *   3.1397301,
     *   3.3618808,
     *   4.050548,
     *   3.9468776,
     *   3.4803611,
     *   3.3322607,
     *   2.0437865,
     *   -2.6361885,
     *   -3.6580815,
     *   -3.8506125,
     *   -3.7025115,
     *   -3.2656155,
     *   -2.7990985,
     *   -2.0363815,
     *   -2.0067615,
     *   -2.3251775,
     *   -5.3908575,
     *   -5.3908575,
     *   -3.8209925,
     *   -3.6136515,
     *   -3.2878305,
     *   -3.1249205,
     *   -2.1918875,
     *   -2.1104315,
     *   -1.3625245,
     *   -1.3699295,
     *   -0.1332904,
     *   -0.1184804,
     *   -1.5402445/


        data cty/ 
     *   0,
     *   -0.2443657,
     *   -0.2517707,
     *   -0.1406954,
     *   0.0666453,
     *   0.0518353,
     *   -0.037025,
     *   -0.1703155,
     *   -0.0592401,
     *   -0.2813908,
     *   -1.0292982,
     *   -1.0292982,
     *   -0.8811977,
     *   -1.2440439,
     *   -1.6735353,
     *   -1.8660659,
     *   -2.1474568,
     *   -2.7324537,
     *   -2.5695432,
     *   -1.9623312,
     *   -1.5698649,
     *   -1.2292339,
     *   -1.1922087,
     *   -1.1107535,
     *   -1.1699937,
     *   -1.0811334,
     *   -0.9034128,
     *   -0.7627173,
     *   -0.6072119,
     *   -0.659047,
     *   -0.5849968,
     *   -0.5924018,
     *   -0.659047,
     *   -0.651642,
     *   -0.5331616,
     *   -0.4591114,
     *   -0.3036059,
     *   -0.2813908,
     *   -0.1110753,
     *   0.0370251/
        

        do 1000 i=1,ncs

        cx(i) = ctx(ncs-i+1)
        cy(i) = cty(ncs-i+1)
        
 1000 continue

        return
        end
c
c
c
c
c
        subroutine kochcorn(cx,cy,ncs)
        implicit real *8 (a-h,o-z)
        dimension cx(92),cy(92),csx(92),csy(92)

        ncs = 92

        data csx/
     *  0,
     *  -11.4786,
     *  -35.3356,
     *  -59.1926,
     *  -47.0386,
     *  -36.0032,
     *  -47.3179,
     *  -58.6078,
     *  -82.1712,
     *  -105.7345,
     *  -116.8453,
     *  -128.7809,
     *  -141.4037,
     *  -153.2015,
     *  -176.97838,
     *  -200.75528,
     *  -188.92258,
     *  -177.08988,
     *  -199.97858,
     *  -176.58048,
     *  -153.1824,
     *  -141.5612,
     *  -129.9401,
     *  -141.3505,
     *  -152.7609,
     *  -176.36968,
     *  -199.97858,
     *  -177.02998,
     *  -188.85398,
     *  -200.67798,
     *  -177.00058,
     *  -153.3232,
     *  -141.505,
     *  -129.4786,
     *  -117.432,
     *  -105.5936,
     *  -82.1702,
     *  -58.7467,
     *  -48.2872,
     *  -35.7689,
     *  -47.3737,
     *  -58.9786,
     *  -35.2286,
     *  -11.4786,
     *  -0.0477,
     *  12.0214,
     *  24.0905,
     *  35.5214,
     *  59.2714,
     *  83.0214,
     *  71.4165,
     *  59.8117,
     *  72.33,
     *  82.7895,
     *  106.213,
     *  129.6364,
     *  141.4748,
     *  153.5214,
     *  165.5478,
     *  177.366,
     *  201.0434,
     *  224.7208,
     *  212.8968,
     *  202.3863,
     *  224.0214,
     *  200.4125,
     *  176.8037,
     *  165.3933,
     *  153.9829,
     *  165.604,
     *  177.2252,
     *  200.6233,
     *  224.0214,
     *  201.1327,
     *  212.9654,
     *  224.7981,
     *  201.0115,
     *  177.2249,
     *  166.1344,
     *  154.2191,
     *  141.576,
     *  129.7578,
     *  106.2042,
     *  82.6506,
     *  71.3607,
     *  60.046,
     *  71.0814,
     *  82.6884,
     *  59.3784,
     *  35.5214,
     *  24.0428,
     *  12.0214/




        data csy/
     *   0,
     *   -19.9901,
     *   -20.0001,
     *   -20.0101,
     *   -41.9958,
     *   -61.82038,
     *   -82.0101,
     *   -101.5101,
     *   -101.77791,
     *   -102.04572,
     *   -82.77791,
     *   -62.25272,
     *   -81.48519,
     *   -101.97505,
     *   -102.24257,
     *   -102.5101,
     *   -122.9821,
     *   -143.4541,
     *   -183.93604,
     *   -185.00991,
     *   -185.00991,
     *   -205.03769,
     *   -225.06547,
     *   -244.78769,
     *   -264.50991,
     *   -264.77801,
     *   -265.04611,
     *   -306.58044,
     *   -327.04517,
     *   -347.50991,
     *   -347.7773,
     *   -348.04469,
     *   -368.5273,
     *   -389.00991,
     *   -368.49227,
     *   -347.97462,
     *   -348.24227,
     *   -348.50991,
     *   -366.50991,
     *   -388.76566,
     *   -408.8286,
     *   -429.45073,
     *   -430.03791,
     *   -430.06591,
     *   -450.03787,
     *   -470.00983,
     *   -450.03787,
     *   -430.06592,
     *   -430.03782,
     *   -430.00982,
     *   -408.82851,
     *   -388.76557,
     *   -366.50982,
     *   -348.50982,
     *   -348.24218,
     *   -347.97453,
     *   -368.49218,
     *   -389.00982,
     *   -368.52721,
     *   -348.0446,
     *   -347.77721,
     *   -347.50982,
     *   -327.04508,
     *   -304.54508,
     *   -266.13768,
     *   -264.77792,
     *   -264.50982,
     *   -244.7876,
     *   -225.06538,
     *   -205.0376,
     *   -185.00982,
     *   -185.00982,
     *   -185.00982,
     *   -143.45382,
     *   -122.98182,
     *   -102.50982,
     *   -102.24227,
     *   -101.97472,
     *   -82.74227,
     *   -62.25244,
     *   -81.52023,
     *   -102.0454,
     *   -101.77761,
     *   -101.50982,
     *   -82.00982,
     *   -61.8201,
     *   -41.99552,
     *   -21.43524,
     *   -19.99982,
     *   -19.98982,
     *   0.00028,
     *   19.99018/



        do 1000 i=1,ncs
        cx(i) = -(csx(i)+csy(i))/50.0d0
        cy(i) = (csx(i)-csy(i))/50.0d0
 1000 continue
    
        return
        end
c
c
c
c
c
        subroutine tricorn(cx,cy,ncs)
        implicit real *8 (a-h,o-z)
        dimension ctx(3),cty(3),cx(1),cy(1)
        save
        ncs = 3

        data ctx/
     *   -0.5d0,
     *   0.5d0,
     *   0.1d0/


        data cty/
     *   0,
     *   0,
     *   0.5d0/

ccc        data cty/
ccc     *   0,
ccc     *   0,
ccc     *   0.8660254037844386d0/

        do 1000 i=1,ncs
        cx(i) = ctx(i)
        cy(i) = cty(i)
 1000 continue

        return
        end
c
c
c
c
c
        subroutine chevcorn(cx,cy,ncs)
        implicit real *8 (a-h,o-z)
        dimension ctx(4),cty(4),cx(1),cy(1)
        save
        ncs = 4

        data ctx/
     *  -0.80d0,
     *   0,
     *   0.80d0,
     *   0/

        data cty/
     *  -2.80d0,
     *   0,
     *  -2.80d0,
     *   4.95d0/
    
        do 1000 i=1,ncs
        cx(i) = ctx(i)
        cy(i) = cty(i)
 1000 continue

        return
        end
c
c
c
c
c

