      implicit real *8 (a-h,o-z)
      parameter (nmax = 20000)
      parameter (nvmax = 1000)
      real *8, allocatable :: ts(:),wts(:)
      real *8, allocatable :: ts2(:),wts2(:),umat(:),vmat(:)
      real *8 ts3(100),wts3(100),utmp,vtmp,tsloc(100),wtsloc(100)
      real *8 qwtstmp(100)

      real *8, allocatable :: verts(:,:),xverts(:),yverts(:)
      integer, allocatable :: el(:),er(:),iregl(:),iregr(:),imid(:)

      real *8, allocatable :: pl(:),pr(:),rnxe(:),rnye(:)

      real *8, allocatable :: xs(:),ys(:),rnx(:),rny(:),qwts(:),
     1   rkappa(:)
      real *8, allocatable :: xs2(:),ys2(:),rnx2(:),rny2(:),qwts2(:),
     1   rkappa2(:)

      integer, allocatable :: lns(:),rns(:),nepts(:)
      integer, allocatable :: lns2(:),rns2(:),nepts2(:)
      real *8, allocatable :: rlen(:)

      integer, allocatable :: icl(:),icr(:),icsgnl(:),icsgnr(:)
      integer, allocatable :: ixmatc(:)

      real *8, allocatable :: xsgnl(:),xsgnr(:),alpha(:)

      real *8 src(2),trg(2)
      real *8, allocatable :: xsrc(:),ysrc(:),charges(:)
      real *8, allocatable :: ztarg(:,:),xt(:),yt(:)
      real *8, allocatable :: potex(:),pottest(:),pottest2(:)
      real *8, allocatable :: rvals(:),tvals(:)

      real *8 grad(2),gradex(2)

      real *8, allocatable :: rhs(:),soln(:)
      real *8, allocatable :: rhs2(:),soln2(:)

      real *8, allocatable :: rhs_px(:),soln_px(:)
      real *8, allocatable :: rhs2_px(:),soln2_px(:)
      real *8, allocatable :: rhs_py(:),soln_py(:)
      real *8, allocatable :: rhs2_py(:),soln2_py(:)

      real *8, allocatable :: xmatc(:,:,:), xmat(:,:)
      real *8, allocatable :: xmatc2(:,:,:), xmat2(:,:),xmat2copy(:,:)
      real *8, allocatable :: xtmp(:,:)

      real *8 pol_t(2,2),pol_t2(2,2),err_t(2,2)

      real *8, allocatable :: targ(:,:),pottarg(:),pottarg2(:)

      real *8 errs(1000),pars(1000)
      integer, allocatable :: ipiv(:),ipiv2(:),ipiv3(:)

      real *8, allocatable :: xmatcnew(:,:),xmatnew(:,:)
      real *8, allocatable :: xmatnewcopy(:,:)
      real *8, allocatable :: xmatsub(:,:)
      real *8, allocatable :: rhsnew(:),solnnew(:),rmutmp(:)
      real *8, allocatable :: rhstmp(:),rhstmp2(:),rhscoeffs(:)
      real *8, allocatable :: solncomp(:)

      external multaslow

c
c
c      Solve three right hand sides,
c        known solution, x polarization, y polarization
c
c      Metrics to compare,
c
c        1. Error in polarization
c
c        2. potential computed using adaptive integration
c        on a polar grid in the vicinity of the corner at the origin
c
c        use 100x100 grid with exponential spacing in r, and compute
c        solution using adaptive integration for all three cases
c        
c        3. density on the fine mesh
c
c        compare solutions with different 4 different levels of resolve
c          10,20,30,40 for (2) and (3), no refinement needed for (1) 
c
c
c


      call prini(6,13)
      call prin2('Enter n*',n,0)
      read *, n

      done = 1
      pi = atan(done)*4

      call prin2('pi=*',pi,1)

c
c      generate diadically refined mesh
c      with refinement at the origin
c      between [0,1]
c

       itype = 1


       k = 16
       irefinelev = 10
       ncorner = k*irefinelev
       allocate(ts(ncorner),wts(ncorner))
       call getcornerdis(k,irefinelev,ts,wts)


       itype = 1
       lu = 100000
       allocate(umat(lu),vmat(lu))
       ln = 50
       allocate(ts2(ln),wts2(ln))
       call lapdisc(ts2,wts2,umat,vmat,ncorner2,itype)

       do i=1,ncorner2
         ts2(i) = 2*ts2(i)
         wts2(i) = 2*wts2(i)
       enddo


c
cc       generate discretization matrix associated with angle
c        theta
c
cc       call prin2('xmat=*',xmatc,24)
cc       call prinf('ncorner=*',ncorner,1)


c
cc       locations of corner vertices
c

       nverts = 3
       allocate(verts(2,nverts),xverts(nverts),yverts(nverts))
       verts(1,1) = 0
       verts(2,1) = 0

       verts(1,2) = sqrt(done*3)/2 + 0.1d0
       verts(2,2) = -done/2

       verts(1,3) = sqrt(done*3)/2
       verts(2,3) = done/2-0.19d0

c
ccc     declare edges
c
c       el(i)  - denotes the left end vertex on the edge 
c       er(i)  - denotes the right end vertex on the edge 
c       rnxe(i) - x component of normal to the edge
c       rnye(i) - y component of normal to the edge
c       pl(i) - length of panel corresponding to corner discretization
c               on left end vertex of the edge
c       pr(i) - length of panel corresponding to corner discretization
c               on right end vertex of the edge
c
c       Note: the lengths of all panels meeting at the same vertex must
c             be the same, the corner correction matrices are generated
c             based on this assumption
c      
c      lns(i) - denotes the starting location in the discretization 
c               nodes array where the corner discretization nodes
c               corresponding to the left end vertex start
c             
c      rns(i) - denotes the starting location in the discretization 
c               nodes array where the corner discretization nodes
c               corresponding to the right end vertex start
c
c      imid(i) - denotes the number of gauss-legendre panels in
c                the middle of the edge
c      kmid - is the number of nodes on each of the gauss-legendre
c             panels


       nedges = 3

       allocate(el(nedges),er(nedges),iregl(nedges),iregr(nedges))
       allocate(pl(nedges),pr(nedges),rnxe(nedges),rnye(nedges))
       allocate(imid(nedges))
       el(1) = 1
       er(1) = 2
       iregl(1) = 1
       iregr(1) = 2

       el(2) = 2
       er(2) = 3
       iregl(2) = 1
       iregr(2) = 2

       el(3) = 3
       er(3) = 1
       iregl(3) = 1
       iregr(3) = 2

       rtmp = 0.15d0

       n = 0
       n2 = 0
       kmid = 16
       do i=1,nedges

         dxt = verts(1,er(i)) - verts(1,el(i))
         dyt = verts(2,er(i)) - verts(2,el(i))

         dst = sqrt(dxt**2 + dyt**2)
         rnxe(i) = dyt/dst
         rnye(i) = -dxt/dst

         pl(i) = rtmp
         pr(i) = rtmp
         imid(i) = 5
         n = n + imid(i)*kmid + 2*ncorner
         n2 = n2 + imid(i)*kmid + 2*ncorner2
       enddo

       allocate(xs(n),ys(n),rnx(n),rny(n),rkappa(n),qwts(n))
       allocate(xs2(n2),ys2(n2),rnx2(n2),rny2(n2),
     1    rkappa2(n2),qwts2(n2))
       allocate(lns(nedges+1),rns(nedges+1),nepts(nedges+1))
       allocate(lns2(nedges+1),rns2(nedges+1),nepts2(nedges+1))
       allocate(rlen(nedges+1))


c
cc      generate discretization nodes
c
c
       call prinf('setting up geometry*',i,0)
     
       call getgeom(nverts,verts,nedges,el,er,rnxe,rnye,pl,pr,imid,
     1         kmid,ncorner,ts,wts,n,xs,ys,rnx,rny,rkappa,qwts,lns,
     2         rns,nepts,rlen)

       call getgeom(nverts,verts,nedges,el,er,rnxe,rnye,pl,pr,imid,
     1         kmid,ncorner2,ts2,wts2,n2,xs2,ys2,rnx2,rny2,rkappa2,
     2         qwts2,lns2,rns2,nepts2,rlen)

      do i=1,nverts
        xverts(i) = verts(1,i)
        yverts(i) = verts(2,i)
      enddo


      call pyplot2(11,xverts,yverts,nverts,3,xs,ys,n,
     1    1,'a*')

      call pyplot2(12,xverts,yverts,nverts,3,xs2,ys2,n2,
     1    1,'a*')

c
c
c
c       generate targets on an exponential grid
c
      nlat = 10
      ntarg = nlat*nlat
      tmin = atan2(verts(2,2),verts(1,2))
      tmax = atan2(verts(2,3),verts(1,3))
      allocate(ztarg(2,ntarg),potex(ntarg),pottest(ntarg))
      allocate(pottest2(ntarg),xt(ntarg),yt(ntarg))
      allocate(rvals(nlat),tvals(nlat))

      print *, tmin, tmax


      do i=1,nlat
        rr = -8*(i-1)/(nlat-1)
        rvals(i) = (10**rr)*rtmp
        tvals(i) = tmin + (i+0.0d0)/(nlat+1.0d0)*(tmax-tmin)
      enddo
      call prin2('rvals=*',rvals,nlat)
      call prin2('tvals=*',tvals,nlat)


      
      do irr = 1,nlat
        do itt = 1,nlat
          ipt = (irr-1)*nlat + itt
          ztarg(1,ipt) = rvals(irr)*cos(tvals(itt))
          ztarg(2,ipt) = rvals(irr)*sin(tvals(itt))
          xt(ipt) = ztarg(1,ipt)
          yt(ipt) = ztarg(2,ipt)
        enddo
      enddo

      write(33,*) verts(1,1),verts(2,1)      
      write(33,*) verts(1,2),verts(2,2)      
      write(33,*) verts(1,3),verts(2,3)      
      do i=1,ntarg
        write(33,*) xt(i),yt(i)
      enddo




c
cc
c      set up number of corner interactions
c      
c      icl(i) - is the edge number on the left
c      icr(i) - is the edge number on the right
c      icsgnl(i) - if 0, then the left end of the
c                  edge icl(i) is involved in the interaction
c                  if 1, then the right end of the
c                  edge icl(i) is involved in the interaction
c      icsgnr(i) - if 0, then the left end of the
c                  edge icr(i) is involved in the interaction
c                  if 1, then the right end of the
c                  edge icr(i) is involved in the interaction
c      alpha(i) - angle at which these edges intersect
c      ixmatc(i) - corner discretization matrix index
c                  corresponding to angle alpha(i)
c
      ncint = 3
      allocate(icl(ncint),icr(ncint),icsgnl(ncint),icsgnr(ncint))
      allocate(alpha(ncint),xsgnl(ncint),xsgnr(ncint))
      allocate(ixmatc(ncint))
      icl(1) = 1
      icr(1) = 3
      icsgnl(1) = 0
      icsgnr(1) = 1


      icl(2) = 1
      icr(2) = 2
      icsgnl(2) = 1
      icsgnr(2) = 0


      icl(3) = 2
      icr(3) = 3
      icsgnl(3) = 1
      icsgnr(3) = 0

      allocate(xmatc(ncorner,ncorner,ncint))
      allocate(xmatc2(ncorner2,ncorner2,ncint))
      do icint = 1,ncint

        if(icsgnl(icint).eq.0.and.icsgnr(icint).eq.0) then
          ivert1 = er(icl(icint))
          ivert2 = el(icl(icint))
          ivert3 = er(icr(icint))

          xsgnl(icint) = -1
          xsgnr(icint) = 1
        endif

        if(icsgnl(icint).eq.1.and.icsgnr(icint).eq.1) then
          ivert1 = el(icl(icint))
          ivert2 = er(icl(icint))
          ivert3 = el(icr(icint))

          xsgnl(icint) = 1
          xsgnr(icint) = -1
        endif


        if(icsgnl(icint).eq.1.and.icsgnr(icint).eq.0) then
          ivert1 = el(icl(icint))
          ivert2 = er(icl(icint))
          ivert3 = er(icr(icint))

          xsgnl(icint) = -1
          xsgnr(icint) = -1
        endif


        if(icsgnl(icint).eq.0.and.icsgnr(icint).eq.1) then

          ivert1 = er(icl(icint))
          ivert2 = el(icl(icint))
          ivert3 = el(icr(icint))

          xsgnl(icint) = 1
          xsgnr(icint) = 1

        endif

        xvert1 = verts(1,ivert1)
        yvert1 = verts(2,ivert1)


        xvert2 = verts(1,ivert2)
        yvert2 = verts(2,ivert2)

        xvert3 = verts(1,ivert3)
        yvert3 = verts(2,ivert3)

        ixmatc(icint) = icint

        dx1 = xvert1-xvert2
        dy1 = yvert1-yvert2
        dr1 = sqrt(dx1**2 + dy1**2)

        dx2 = xvert3-xvert2
        dy2 = yvert3-yvert2
        dr2 = sqrt(dx2**2 + dy2**2)

        drp = dx1*dx2 + dy1*dy2
        drp = drp/(dr1*dr2)
        alpha(icint) = acos(drp)

        alpha(icint) = atan2(dy2,dx2) - atan2(dy1,dx1)
cc      if(alpha(icint).lt.0) alpha(icint) = alpha(icint)+pi

        rpan = 0.3d0
        thet = alpha(icint)

        call getcornermat(thet,ncorner,rpan,ts,wts,xmatc(1,1,icint))


        thet = alpha(icint)/pi
        if(thet<0) thet = thet+2

        call lapcornmat(thet,xmatc2(1,1,icint),ncorner2)
        do ipt=1,ncorner2
          do jpt=1,ncorner2
            xmatc2(ipt,jpt,icint) = xmatc2(ipt,jpt,icint)/2/pi
          enddo
        enddo

      enddo


      call prin2('alpha=*',alpha,ncint)

      call prinf('n=*',n,1)
      call prinf('n2=*',n2,1)
      allocate(xmat(n,n),xmat2(n2,n2))
      allocate(xmat2copy(n2,n2))


c
cc      generate matrix edge by edge
c


       call prinf('building matrix*',i,0)

      rfac = 1.0d0
      do iedge=1,nedges
        do jedge=1,nedges

          nss = nepts(jedge)
          nts = nepts(iedge)
          allocate(xtmp(nts,nss))

          call getedgemat(iedge,jedge,n,xs,ys,rnx,rny,rkappa,qwts,
     1      lns,rns,nepts,ncint,icl,icr,icsgnl,icsgnr,alpha,ixmatc,
     2      xsgnl,xsgnr,ncorner,xmatc,nts,nss,xtmp)
      
          its = lns(iedge) -1
          iss = lns(jedge) -1

          call xreplmat(nts,nss,n,its,iss,xtmp,xmat,rfac)

          deallocate(xtmp)

          nss = nepts2(jedge)
          nts = nepts2(iedge)
          allocate(xtmp(nts,nss))

          call getedgemat(iedge,jedge,n2,xs2,ys2,rnx2,rny2,rkappa2,
     1      qwts2,lns2,rns2,nepts2,ncint,icl,icr,icsgnl,icsgnr,
     2      alpha,ixmatc,xsgnl,xsgnr,ncorner2,xmatc2,nts,nss,xtmp)
      
          its = lns2(iedge) -1
          iss = lns2(jedge) -1

          call xreplmat(nts,nss,n2,its,iss,xtmp,xmat2,rfac)

          deallocate(xtmp)
          
        enddo
      enddo

      call prinf('done building matrix*',i,0)





c
cc       set up sources for the rhs
c

      ncharges = 2

      allocate(xsrc(ncharges),ysrc(ncharges),charges(ncharges))
 
      do i=1,ncharges
         thet = hkrand(0)*2*pi
         rr = (0.7d0 + hkrand(0))*10.0d0
         xsrc(i) = 0.5d0 + rr*cos(thet)
         ysrc(i) = rr*sin(thet)
         charges(i) = hkrand(0)*10
      enddo

cc      call prin2('xsrc=*',xsrc,ncharges)
cc      call prin2('ysrc=*',ysrc,ncharges)


      allocate(rhs(n),soln(n))
      allocate(rhs_px(n),soln_px(n))
      allocate(rhs_py(n),soln_py(n))
      allocate(rhs2(n2),soln2(n2))
      allocate(rhs2_px(n2),soln2_px(n2))
      allocate(rhs2_py(n2),soln2_py(n2))


c
cc      fix the diagonal and get rhs
c

      call prinf('ncharges=*',ncharges,1)

      ra = 0
      do iedge=1,nedges
        do ipt = 1,nepts(iedge)
          i = lns(iedge) + ipt-1
          xmat(i,i) = 0.5d0

          call getrhs(ncharges,xsrc,ysrc,charges,xs(i),ys(i),pot,grad)

          rhs(i) = sqrt(qwts(i))*(grad(1)*rnxe(iedge)+
     1       grad(2)*rnye(iedge))
          
          rhs_px(i) = sqrt(qwts(i))*rnxe(iedge)
          rhs_py(i) = sqrt(qwts(i))*rnye(iedge)

        enddo

        
        do ipt = 1,nepts2(iedge)
          i = lns2(iedge) + ipt-1
          xmat2(i,i) = 0.5d0

          call getrhs(ncharges,xsrc,ysrc,charges,xs2(i),ys2(i),pot,
     1        grad)

          rhs2(i) = sqrt(qwts2(i))*(grad(1)*rnxe(iedge)+
     1       grad(2)*rnye(iedge))

          rhs2_px(i) = sqrt(qwts2(i))*rnxe(iedge)
          rhs2_py(i) = sqrt(qwts2(i))*rnye(iedge)

        enddo
      enddo

      do i=1,n
        do j=1,n
          xmat(j,i) = xmat(j,i) + sqrt(qwts(j)*qwts(i))
        enddo
      enddo

      do i=1,n2
        do j=1,n2
          xmat2(j,i) = xmat2(j,i) + sqrt(qwts2(j)*qwts2(i))
          xmat2copy(j,i) = xmat2(j,i)
        enddo
      enddo


      do i=1,n
        soln(i) = rhs(i)
        soln_px(i) = rhs_px(i)
        soln_py(i) = rhs_py(i)
      enddo

      do i=1,n2
        soln2(i) = rhs2(i)
        soln2_px(i) = rhs2_px(i)
        soln2_py(i) = rhs2_py(i)
      enddo

       call prinf('end of building matrix*',i,0)

      info = 0
      allocate(ipiv(n))

      call dgetrf(n,n,xmat,n,ipiv,info)

      info = 0
      call dgetrs('t',n,1,xmat,n,ipiv,soln,n,info)

      info = 0
      call dgetrs('t',n,1,xmat,n,ipiv,soln_px,n,info)

      info = 0
      call dgetrs('t',n,1,xmat,n,ipiv,soln_py,n,info)

      info = 0
      allocate(ipiv2(n2))

      call dgetrf(n2,n2,xmat2,n2,ipiv2,info)

      info = 0
      call dgetrs('t',n2,1,xmat2,n2,ipiv2,soln2,n2,info)

      info = 0
      call dgetrs('t',n2,1,xmat2,n2,ipiv2,soln2_px,n2,info)

      info = 0
      call dgetrs('t',n2,1,xmat2,n2,ipiv2,soln2_py,n2,info)



c
cc      test solution 
c

      trg(1) = 0.7d0
      trg(2) = -0.01d0

      call getrhs(ncharges,xsrc,ysrc,charges,trg(1),trg(2),
     1    potex_1,gradex)

      pot1 = 0
      do i=1,n
        rr = (trg(1)-xs(i))**2 + (trg(2)-ys(i))**2
        pot1 = pot1 - log(rr)/4/pi*soln(i)*sqrt(qwts(i))
      enddo

      pot2 = 0
      do i=1,n2
        rr = (trg(1)-xs2(i))**2 + (trg(2)-ys2(i))**2
        pot2 = pot2 - log(rr)/4/pi*soln2(i)*sqrt(qwts2(i))
      enddo

      call prin2('pot=*',pot1,1)
      call prin2('potex=*',potex_1,1)

cc      print *, pot,potex,pot-potex

      rdiff1 = potex_1-pot1
      rdiff2 = potex_1-pot2


      trg(1) = 0.65d0
      trg(2) = -0.01d0

      call getrhs(ncharges,xsrc,ysrc,charges,trg(1),trg(2),
     1    potex_2,gradex)

      pot1_2 = 0
      do i=1,n
        rr = (trg(1)-xs(i))**2 + (trg(2)-ys(i))**2
        pot1_2 = pot1_2 - log(rr)/4/pi*soln(i)*sqrt(qwts(i))
      enddo

      

      pot2_2 = 0
      do i=1,n2
        rr = (trg(1)-xs2(i))**2 + (trg(2)-ys2(i))**2
        pot2_2 = pot2_2 - log(rr)/4/pi*soln2(i)*sqrt(qwts2(i))
      enddo

      call prin2('pot=*',pot1_2,1)
      call prin2('potex=*',potex_2,1)

cc      print *, pot,potex,pot-potex



      rdiff1_2 = potex_2-pot1_2
      rdiff2_2 = potex_2-pot2_2


c
c
c       since rdiff1, and rdiff1_2 agree to machine precision
c        without loss of generality we will treat rdiff1 
c        as the constant differing between the two solutions
c        for the purposes of comparing the potential in
c        the volume
c


      erra = abs(rdiff1-rdiff1_2)/abs(rdiff1)
      erra2 = abs(rdiff2-rdiff2_2)/abs(rdiff2)


      call prin2('error bisection=*',erra,1)
      call prin2('error sp-dis=*',erra2,1)

c
c
c
c        start test 1 - polarization computation
c
c



      pol_t(1,1)=0
      pol_t(1,2)=0
      pol_t(2,1)=0
      pol_t(2,2)=0

      pol_t2(1,1)=0
      pol_t2(1,2)=0
      pol_t2(2,1)=0
      pol_t2(2,2)=0

      do i=1,n
        pol_t(1,1) = pol_t(1,1) + soln_px(i)*xs(i)*sqrt(qwts(i))
        pol_t(2,1) = pol_t(2,1) + soln_px(i)*ys(i)*sqrt(qwts(i))
        
        pol_t(1,2) = pol_t(1,2) + soln_py(i)*xs(i)*sqrt(qwts(i))
        pol_t(2,2) = pol_t(2,2) + soln_py(i)*ys(i)*sqrt(qwts(i))
      enddo

      do i=1,n2
        pol_t2(1,1) = pol_t2(1,1) + soln2_px(i)*xs2(i)*sqrt(qwts2(i))
        pol_t2(2,1) = pol_t2(2,1) + soln2_px(i)*ys2(i)*sqrt(qwts2(i))
        
        pol_t2(1,2) = pol_t2(1,2) + soln2_py(i)*xs2(i)*sqrt(qwts2(i))
        pol_t2(2,2) = pol_t2(2,2) + soln2_py(i)*ys2(i)*sqrt(qwts2(i))
      enddo

      do i=1,2
        do j=1,2
          err_t(i,j)= abs(pol_t(i,j)-pol_t2(i,j))/abs(pol_t(i,j))
        enddo
      enddo

      call prin2('error in polarization tensor=*',err_t,4)

c
c
c
c        start test 2: accuracy in targets in the 
c        vicinity of one of the corners.
c    
c      attempt number 1: try the far-field quadrature
c          - fails due to some of the targets being close
c            to the edge
c
c      attempt number 2: use far-field quadrature except 
c        near the corner, and use adaptive integration
c        use to corner 1
c          
c      
c
c
      do i=1,ntarg
        call getrhs(ncharges,xsrc,ysrc,charges,xt(i),yt(i),potex(i),
     1     grad)
      enddo

      call cpu_time(t1)
C$      t1 = omp_get_wtime()      
      call comppottarg_adap(ntarg,xt,yt,n,xs,ys,soln,qwts,
     1    nedges,nepts,lns,rns,imid,k,irefinelev,pottest)
      

      call cpu_time(t2)
C$      t2 = omp_get_wtime()     

      call prin2('time taken in adaptive integration=*',t2-t1,1)
      
      erra = 0
      ra = 0

      do i=1,ntarg
        pottest(i) = pottest(i)+rdiff1
        ra = ra + potex(i)**2
        erra = erra + (pottest(i)-potex(i))**2
      enddo
      erra = sqrt(erra/ra)
      ra = sqrt(ra)
      call prin2("ra=*",ra,1)
      call prin2('error in targets in volume-bisection=*',erra,1)




c
cc      resolve problem at the corner panel of 
c       vertex 1
c
      
      rpan = rtmp
      nc2 = k + ncorner2
      itype = 1
      call legeexps(itype,k,ts3,utmp,vtmp,wts3)

      xstart = 1.0d0/2
      xend = rpan
      do i=1,k
        tsloc(i) = 0.5d0 + (ts3(i)+1)/4.0d0
        wtsloc(i) = wts3(i)/4.0d0
      enddo

      do i=1,ncorner2
        tsloc(k+i) = ts2(i)/2
        wtsloc(k+i) = wts2(i)/2
      enddo

      do i=1,nc2
        qwtstmp(i) = wtsloc(i)/2
        qwtstmp(i+nc2) = wtsloc(i)/2
      enddo



      icint = 1
      thet = alpha(1)
      nn = 2*nc2

      allocate(xmatcnew(nn,nn))

      call getcornermat(thet,nc2,rpan,tsloc,wtsloc,xmatcnew)

      allocate(xmatnew(nn,nn))
      allocate(xmatnewcopy(nn,nn))


      do i=1,nn
        do j=1,nn
          xmatnew(i,j) = 0
        enddo
      enddo
c
cc      set off-diagonal blocks
c
c      unknowns 1-nc2 are the unknowns
c      on edge 1 at vertex 1,
c      and unknowns nc2+1,nn are the unknowns
c      on edge 3 at vertex 1
c
       call xreplmat(nc2,nc2,nn,0,nc2,xmatcnew,xmatnew,xsgnl(1))
       call xreplmat(nc2,nc2,nn,nc2,0,xmatcnew,xmatnew,xsgnr(1))

       call xreplmat(ncorner2,ncorner2,nn,k,nc2+k,xmatc2(1,1,1),
     1           xmatnew,xsgnl(1))

       call xreplmat(ncorner2,ncorner2,nn,nc2+k,k,xmatc2(1,1,1),
     1           xmatnew,xsgnr(1))
c
cc        the matrix that has been set up are the dirichlet 
c         matrices. Now take the transposes to compute the 
c         neumann matrices
c
      ra = 0
      do i=1,nn
        xmatnew(i,i) = 0.5d0
      enddo

      do i=1,nn
        do j=1,nn
          xmatnew(j,i) = xmatnew(j,i) + sqrt(qwtstmp(i)*qwtstmp(j))
        enddo
      enddo



      allocate(ipiv3(nn))
      call dgetrf(nn,nn,xmatnew,nn,ipiv3,info)

       
      allocate(rhsnew(nn),solnnew(nn))
      allocate(rmutmp(2*ncorner2))
c
cc     now set up the right hand side for the linear system
c

      rint1 = 0
      istart= lns2(1) - 1
      do i=1,ncorner2
        rmutmp(i) = soln2(i+istart)
        rint1 = rint1 + soln2(i+istart)*sqrt(qwts2(istart+i))
      enddo

      istart = rns2(3) -1
      do i=1,ncorner2
        rmutmp(ncorner2+i) = soln2(i+istart)
        rint1 = rint1 + soln2(i+istart)*sqrt(qwts2(istart+i))
      enddo

      call prin2('integral of density on patch=*',rint1,1)



c
cc       initialize matrix for computing right hand side
c
c

      nnn = 2*ncorner2
      allocate(xmatsub(nnn,nnn))
      istart = lns2(1) - 1
      jstart = lns2(1) - 1

      do i=1,ncorner2
        do j=1,ncorner2
          xmatsub(i,j) = xmat2copy(i+istart,j+jstart)
        enddo
      enddo

      istart = lns2(1) - 1
      jstart = rns2(3) - 1
      do i=1,ncorner2
        do j=1,ncorner2
          xmatsub(i,j+ncorner2) = xmat2copy(i+istart,j+jstart)
        enddo
      enddo

      istart = rns2(3) - 1
      jstart = lns2(1) - 1

      do i=1,ncorner2
        do j=1,ncorner2
          xmatsub(i+ncorner2,j) = xmat2copy(i+istart,j+jstart)
        enddo
      enddo

      istart = rns2(3) - 1
      jstart = rns2(3) - 1
      do i=1,ncorner2
        do j=1,ncorner2
          xmatsub(i+ncorner2,j+ncorner2) = xmat2copy(i+istart,j+jstart)
        enddo
      enddo

c
cc        start iterative solve loop
c
      nlev = 1

      allocate(rhstmp(nnn),rhstmp2(ncorner2),rhscoeffs(ncorner2))
      allocate(solncomp(2*nlev*k+ncorner2*2))
      nhalf = nlev*k + ncorner2


      alpha = 1.0d0
      beta = 0
      do ilev = 1,nlev
         call dgemv('t',nnn,nnn,alpha,xmatsub,nnn,rmutmp,1,beta,
     1      rhstmp,1)
cc          call multaslow(nnn,rmutmp,rhstmp,xmatneuloc)

c
cc       extract out the relevant pieces of rhs
c
         do i=1,ncorner2
           rhstmp2(i) = rhstmp(i)
           rhscoeffs(i) = 0
        enddo
       
c
cc       smear this function onto the rest of the grid

cc        call prin2('rhstmp2=*',rhstmp2,ncorner2)
cc        call prin2('vmat=*',vmat,ncorner2*ncorner2)
        call dgemv('n',ncorner2,ncorner2,alpha,vmat,ncorner2,rhstmp2,1,
     1     beta,rhscoeffs,1)
cc         call multaslow(ncorner,rhstmp2,rhscoeffs,vmat)
        call prin2('rhscoeffs edge 1=*',rhscoeffs,ncorner2)

        do ipt = 1,nc2
          x = tsloc(ipt)/2.0d0
          rhsnew(ipt) = 0
          do j=1,ncorner2
            val = 0
            call lapeval(x,j,val)
            rhsnew(ipt) = rhsnew(ipt) + rhscoeffs(j)*val
          enddo
          rhsnew(ipt) = rhsnew(ipt)*sqrt(wtsloc(ipt))
        enddo

        do i=1,ncorner2
          rhstmp2(i) = rhstmp(ncorner2+i)
          rhscoeffs(i) = 0
        enddo
       
c
cc       smear this function onto the rest of the grid
        call dgemv('n',ncorner2,ncorner2,alpha,vmat,ncorner2,rhstmp2,1,
     1     beta,rhscoeffs,1)
cc         call multaslow(ncorner,rhstmp2,rhscoeffs,vmat)
        call prin2('rhscoeffs edge 2=*',rhscoeffs,ncorner2)

        do ipt = 1,nc2
          ii = nc2 + ipt
          x = tsloc(ipt)/2.0d0
          rhsnew(ii) = 0
          do j=1,ncorner2
            val = 0
            call lapeval(x,j,val)
            rhsnew(ii) = rhsnew(ii) + rhscoeffs(j)*val
          enddo
          rhsnew(ii) = rhsnew(ii)*sqrt(wtsloc(ipt))
        enddo
        call prin2('rhsnew=*',rhsnew,nn)

        do i=1,nn
          solnnew(i) = rhsnew(i)
        enddo

cc        do i=1,nn
cc          do j=1,nn
cc            xmatnew(j,i) = xmatnewcopy(j,i) + 
cc     1         sqrt(qwtstmp(i)*qwtstmp(j))/2**(ilev-1)
cc          enddo
cc        enddo
        
        info = 0
        call dgetrs('t',nn,1,xmatnew,nn,ipiv3,solnnew,nn,info)

        rint2 = 0
        do i=1,nn
          rint2 = rint2 + solnnew(i)*sqrt(qwtstmp(i))
        enddo

        call prin2('rint2=*',rint2,1)

cc         call dgmres(ier,nn,multaslow,xmatnew2,rhsnew,eps,numit,
cc     1       solnnew,niter,errs,ngmrec,work)

        istart = (ilev-1)*k
        do i=1,k
          solncomp(istart+i) = solnnew(i)
          solncomp(nhalf+i+istart) = solnnew(nc2+i)
        enddo

        if(ilev.eq.nlev) then
          do i=1,ncorner2
            solncomp(nlev*k+i) = solnnew(k+i)
            solncomp(nhalf+nlev*k+i) = solnnew(nc2+k+i)
          enddo
        endif

c
cc        reinitialize rmutmp
c

        do i=1,ncorner2
          rmutmp(i) = solnnew(i+k)
        enddo
        do i=1,ncorner
          rmutmp(ncorner2+i) = solnnew(i+k+nc2)
        enddo
      enddo
c
c
c
      erra = 0
      ra = 0
      istart = lns(1)-1
      do i=1,nlev*k
        ra = ra + soln(istart+i)**2
        erra = erra + (solncomp(i)-soln(istart+i))**2
      enddo
      call prin2('soln=*',soln(lns(1)),nlev*k)
      call prin2('solncomp=*',solncomp,nlev*k)
      erra = sqrt(erra/ra)

      call prin2("error in density after resolve=*",erra,1)

       stop
       end

c---------------------------------------------

      subroutine comppottarg_adap(ntarg,xt,yt,n,xs,ys,soln,qwts,
     1  nedges,nepts,lns,rns,imid,k,irefinelev,pot)
      implicit real *8 (a-h,o-z)
      integer ntarg,n,nedges
      real *8 xt(ntarg),yt(ntarg),xs(n),ys(n),soln(n),qwts(n)
      integer lns(nedges+1),rns(nedges+1),imid(nedges),nepts(nedges)
      integer k,kmid,irefinelev
      real *8 pot(ntarg)
      real *8, allocatable :: soln_tmp(:)
      real *8, allocatable :: soln_coefs(:)
      real *8, allocatable :: xcoefs(:)
      real *8, allocatable :: ycoefs(:)
      real *8, allocatable :: rpanlen(:)
      real *8 ts(k),umat(k,k),vmat(k,k),wts(k)
      real *8 tsquad(100),wquad(100)
      real *8 stack(2,200),vals(200),par1(4)
      real *8, allocatable :: tpack(:)


      external slp_adap

      done = 1
      pi = atan(done)*4

      itype = 2
      call legeexps(itype,k,ts,umat,vmat,wts)
      allocate(soln_tmp(n),soln_coefs(n),xcoefs(n),ycoefs(n))

      do i=1,n
        soln_tmp(i) = soln(i)/sqrt(qwts(i))
      enddo

      npan = n/k
      allocate(rpanlen(npan))

      ii = 1
      do iedge=1,nedges
        npan_edge = nepts(iedge)/k
        do ipan = 1,npan_edge
          rpanlen(ii) = 0
          do j=1,k
            ipt = lns(iedge) + (ipan-1)*k +j-1
            rpanlen(ii) = rpanlen(ii) + qwts(ipt)
            soln_coefs(ipt) = 0
            xcoefs(ipt) = 0
            ycoefs(ipt) = 0
            do l=1,k
              jpt = lns(iedge) + (ipan-1)*k+l-1
              soln_coefs(ipt) = soln_coefs(ipt) +
     1           umat(j,l)*soln_tmp(jpt)
              xcoefs(ipt) = xcoefs(ipt) + umat(j,l)*xs(jpt)
              ycoefs(ipt) = ycoefs(ipt) + umat(j,l)*ys(jpt)
            enddo
          enddo
          ii = ii+1
        enddo
      enddo

cc      call prin2('rpanlen=*',rpanlen,npan)

      do i=1,n
        write(35,'(5(2x,e11.5))') xs(i),ys(i),xcoefs(i),ycoefs(i),
     1     soln_coefs(i)
      enddo

      a = -1
      b = 1
      m = 20
      eps = 1.0d-13
      maxdepth = 200
      nnmax = 10000

      ifwhts = 1
      call legewhts(m,tsquad,wquad,ifwhts)

      allocate(tpack(3*k))

      call prin2('tsquad=*',tsquad,m)
      call prin2('wquad=*',wquad,m)

      

      par1(1) = k+0.1d0

      


ccC$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(itarg,par1,i,j,ipt)
ccC$OMP$PRIVATE(tpack,pottmp,ier,stack,vals,maxrec,numint)
      do itarg=1,ntarg
        pot(itarg) = 0
        par1(2) = xt(itarg)
        par1(3) = yt(itarg)

        do i=1,npan
          par1(4) = rpanlen(i)/2/4/pi
          do j=1,k
            ipt = (i-1)*k+j
            tpack(j) = xcoefs(ipt)
            tpack(j+k) = ycoefs(ipt)
            tpack(j+2*k) = soln_coefs(ipt)
          enddo

          do j=1,200
            stack(1,j) = 0
            stack(2,j) = 0
            vals(j) = 0
          enddo


          pottmp = 0.0d0
          ier = 0
          maxrec = 0
          numint = 0
          call adinrec(ier,stack,a,b,slp_adap,par1,tpack,tsquad,wquad,
     1      m,vals,nnmax,eps,pottmp,maxdepth,maxrec,numint)
          
          pot(itarg) = pot(itarg) + pottmp
        enddo
      enddo
ccC$OMP END PARALLEL DO      
      
      
      
      return
      end
c
c
c
c
c
      subroutine slp_adap(x,par1,par2,val)
      implicit real *8 (a-h,o-z)
      real *8 par1(3),par2(*)
      real *8 pols(100)

      dd = 0
      xx = 0
      yy = 0
      k = par1(1)
      xt = par1(2)
      yt = par1(3)
      rlen = par1(4)
      call legepols(x,k-1,pols)
      
      do i=1,k
        xx = xx + par2(i)*pols(i)
        yy = yy + par2(i+k)*pols(i)
        dd = dd + par2(i+2*k)*pols(i)
      enddo
      
      r = (xx-xt)**2 + (yy-yt)**2
      val = -log(r)*dd*par1(4)

      

      return
      end




       subroutine getcornermat(thet,n,rlen,ts,wts,xmat)
       implicit real *8 (a-h,o-z)
       real *8 ts(*),wts(*),xmat(n,*)

       done = 1
       pi = atan(done)*4

       sint = sin(thet)
       cost = cos(thet)

       do i=1,n
       do j=1,n

       xmat(i,j) = -ts(i)*sint/(ts(i)**2+ts(j)**2-2*ts(i)*ts(j)*cost)
     1        /(2*pi)*sqrt(wts(i)*wts(j))

       enddo
       enddo

       return
       end
c--------------------------------------------
       
       
      subroutine getrhs(nc,xsrc,ysrc,charges,x,y,pot,grad)
      implicit real *8 (a-h,o-z)
      real *8 xsrc(*),ysrc(*),charges(*),grad(2)

      pot = 0
      grad(1) = 0
      grad(2) = 0


      do i=1,nc

         dx = x-xsrc(i)
         dy = y-ysrc(i)

         rrr = sqrt(dx**2 + dy**2)

         pot = pot + log(rrr)*charges(i)
         grad(1) = grad(1) + dx/rrr**2*charges(i)
         grad(2) = grad(2) + dy/rrr**2*charges(i)

      enddo


      return
      end
c------------------------------------       
      subroutine comppot(trg,n,xs,ys,rnx,rny,qwts,soln,pot)
      implicit real *8 (a-h,o-z)
      real *8 trg(*),xs(*),ys(*),rnx(*),rny(*),qwts(*),soln(*)

      done = 1
      pi = atan(done)*4

      pot = 0
      do i=1,n

      rr = (trg(1)-xs(i))**2 + (trg(2)-ys(i))**2
      rrn = (trg(1)-xs(i))*rnx(i) + (trg(2)-ys(i))*rny(i)

      pot = pot + rrn/rr*sqrt(qwts(i))*soln(i)/(2*pi)

      enddo


      return
      end
c---------------------------------------------      
      subroutine comppot2(trg,n,xs,ys,rnx,rny,qwts,soln,solnneu,
     1  potd,pots)
      implicit real *8 (a-h,o-z)
      real *8 trg(*),xs(*),ys(*),rnx(*),rny(*),qwts(*),soln(*)
      real *8 solnneu(*)


      done = 1
      pi = atan(done)*4

      potd = 0
      pots = 0
      do i=1,n

      rr = (trg(1)-xs(i))**2 + (trg(2)-ys(i))**2
      rrn = (trg(1)-xs(i))*rnx(i) + (trg(2)-ys(i))*rny(i)

      potd = potd + rrn/rr*sqrt(qwts(i))*soln(i)/(2*pi)
      pots = pots - log(rr)*sqrt(qwts(i))*solnneu(i)/(4*pi)

      enddo


      return
      end
c---------------------------------------------      
      subroutine xreplmat(nts,nss,n,its,iss,xmatc,xmat,rfac)
c
cc      this subroutine replaces the subblock of the matrix
c       xmat assoicated with corner interactions at
c       edge ie by the corresponding corner matrix xmatc
c
c       input:
c       ie - edge number
c       nc - size of matrix block to be replaced
c       n - size of full system matrix
c       its - target start point
c       iss - source start point
c       rfac - scaling factor
c


      implicit real *8 (a-h,o-z)
      real *8 xmatc(nts,*),xmat(n,*)

      
      do i=1,nts
      do j=1,nss

      xmat(its+i,iss+j) = rfac*xmatc(i,j)

      enddo
      enddo

      return
      end
c-----------------------------------------
      subroutine getcornerdis(k,nlev,ts,wts)
      implicit real *8 (a-h,o-z)
c
cc       this subroutine computes a diadyically refined
c        mesh on [0,1], with nlev number of panels of
c        k Gauss-legendre nodes each
c
c        input:
c        k - Number of gauss legendre points per panel
c        nlev - number of dyadically refined intervals
c
c        output
c        ts(nlev*k) - location of nodes
c        wts(nlev*k) - quadrature weights at nodes
 

      real *8 ts(*), wts(*), ts0(k),wts0(k)

      itype = 1
      call legeexps(itype,k,ts0,u,v,wts0)

       hcur = 0.5d0

       tstart = 0.5d0
       tend = 1.0d0

       hcur = 1.0d0

       ii = 1
       do i=1,nlev

       if(i.ne.nlev) tstart = hcur/2
       if(i.eq.nlev) tstart = 0
       tend = hcur

       do j=1,k

       ts(ii) = tstart + (tend-tstart)*(ts0(j)+1)/2
       wts(ii) = (tend-tstart)*wts0(j)/2

       ii = ii+1

       enddo
       hcur = hcur/2

       enddo

      

      return
      end
c-------------------------------------------------    
c
       subroutine getgeom(nverts,verts,nedges,el,er,rnxe,rnye,pl,pr,
     1    imid,kmid,ncorner,ts,wts,n,xs,ys,rnx,rny,rkappa,qwts,lns,
     2    rns,nepts,rlen)
       
       implicit real *8 (a-h,o-z)

c        this subroutine computes the discretization for the given
c        geometry
c  
c        input:
c       nverts - number of vertices
c       verts(2,nverts) - locations of vertices
c       nedges - number of edges
c       el(i)  - denotes the left end vertex on the edge 
c       er(i)  - denotes the right end vertex on the edge 
c       rnxe(i) - x component of normal to the edge
c       rnye(i) - y component of normal to the edge
c       pl(i) - length of panel corresponding to corner discretization
c               on left end vertex of the edge
c       pr(i) - length of panel corresponding to corner discretization
c               on right end vertex of the edge
c      imid(i) - denotes the number of gauss-legendre panels in
c                the middle of the edge
c      kmid - is the number of nodes on each of the gauss-legendre
c             panels
c      ncorner - number of discretization nodes at each corner
c      ts - location of corner discretization nodes
c      wts - quadrature weights corresponding to corner discretization
c            nodes
c      
c      output:
c       n  - total number of discretization nodes
c       xs,ys - x and y coordinates of discretization nodes
c       rnx,rny - x and y components of normals at the
c                 discretization nodes
c       rkappa - curvature at the discretization nodes
c       qwts - quadrature weights at the discretization nodes
c      lns(i) - denotes the starting location in the discretization 
c               nodes array where the corner discretization nodes
c               corresponding to the left end vertex start
c             
c      rns(i) - denotes the starting location in the discretization 
c               nodes array where the corner discretization nodes
c               corresponding to the right end vertex start
c      rlen(ie) - length of the edges
c
       real *8 verts(2,*),rnxe(*),rnye(*),pl(*),pr(*),ts(*),wts(*)
       integer el(*),er(*),imid(*)
       real *8 xs(*),ys(*),rnx(*),rny(*),rkappa(*),qwts(*),rlen(*)
       
       integer lns(*),rns(*),nepts(*)

       real *8 ts0(kmid),wts0(kmid)
     
       n = 0
       itype = 1
       call legeexps(itype,kmid,ts0,uk,vk,wts0)

       do ie=1,nedges 

       lns(ie) = n+1

       dxt = verts(1,er(ie))-verts(1,el(ie))
       dyt = verts(2,er(ie))-verts(2,el(ie))

       dst = sqrt(dxt**2+dyt**2)
       dxt = dxt/dst
       dyt = dyt/dst

       rlen(ie) = dst

       do 1100 i=1,ncorner

       xs(n+i) = verts(1,el(ie)) + dxt*pl(ie)*ts(i)
       ys(n+i) = verts(2,el(ie)) + dyt*pl(ie)*ts(i)

       rnx(n+i) = rnxe(ie)
       rny(n+i) = rnye(ie)

       rkappa(n+i) = 0

       qwts(n+i) = pl(ie)*wts(i) 

 1100  continue      

       n = n+ ncorner
       rmid = rlen(ie) - pl(ie)-pr(ie)
       
       rpanel =  rmid/imid(ie)
       rstart = pl(ie)

       ii = 1

       do 1110 i=1,imid(ie)

       do 1105 j=1,kmid

       xs(n+ii) = verts(1,el(ie)) + dxt*(rstart+(ts0(j)+1)/2*rpanel)
       ys(n+ii) = verts(2,el(ie)) + dyt*(rstart+(ts0(j)+1)/2*rpanel)

       rnx(n+ii) = rnxe(ie)
       rny(n+ii) = rnye(ie)

       qwts(n+ii) = wts0(j)*rpanel/2

       ii = ii+1

 1105 continue

      rstart = rstart + rpanel 

 1110  continue       

      n = n + imid(ie)*kmid

      rns(ie) = n+1

       do 1120 i=1,ncorner

       xs(n+i) = verts(1,er(ie)) - dxt*ts(i)*pr(ie)
       ys(n+i) = verts(2,er(ie)) - dyt*ts(i)*pr(ie)

       rnx(n+i) = rnxe(ie)
       rny(n+i) = rnye(ie)

       qwts(n+i) = pr(ie)*wts(i)

 1120  continue    

       nepts(ie) = 2*ncorner + imid(ie)*kmid

       n = n+ ncorner

       enddo

       return
       end
c-----------------------------------------------------------      

      subroutine getedgemat(iedge,jedge,n,xs,ys,rnx,rny,rkappa,qwts,
     1  lns,rns,nepts,ncint,icl,icr,icsgnl,icsgnr,alpha,ixmatc,xsgnl,
     2  xsgnr,ncorner,xmatc,nts,nss,xmat)
       implicit real *8 (a-h,o-z)
c
cc       this subroutine computes the dipole interaction matrix 
c        between edges iedge and jedge excluding the identity
c        term
c
c      input
c      iedge - target edge 
c      jedge - source edge
c      n - number of discretization points
c      xs,ys - x,y locations of discretization points
c      rnx,rny - x,y components of normals at discretization points
c      rkappa - curvature information at discretization points
c      qwts - quadrature weights at discretization points
c      lns - corner panel discretization start at left end
c            of edge
c      rns - corner panel discretization start at right end of edge
c      nepts - number of discretization points on edge
c      ncint - total number of corner interactions
c      icl(i) - is the edge number on the left
c      icr(i) - is the edge number on the right
c      icsgnl(i) - if 0, then the left end of the
c                  edge icl(i) is involved in the interaction
c                  if 1, then the right end of the
c                  edge icl(i) is involved in the interaction
c      icsgnr(i) - if 0, then the left end of the
c                  edge icr(i) is involved in the interaction
c                  if 1, then the right end of the
c                  edge icr(i) is involved in the interaction
c      xsgnl(i) - sign for replacing cornermat when icl is target
c                 edge and icr is source edge
c      xsgnr(i) - sign for replacing cornermat when icr is target
c                 edge and icl is source edge
c                 
c      alpha(i) - angle at which these edges intersect
c      ixmatc(i) - corner discretization matrix index
c                  corresponding to angle alpha(i)
c      ncorner - number of discertization nodes on corner panel
c      xmatc - correction matrices for corner interactions
c      nts - number of discretization points on iedge
c      nss - number of discretization points on jedge
c
c      xmat - interaction matrix
c
       real *8 xs(*),ys(*),rnx(*),rny(*),rkappa(*),qwts(*)
       integer lns(*),rns(*),nepts(*),icl(*),icr(*),icsgnl(*)
       integer icsgnr(*),ixmatc(*)
       real *8 alpha(*),xmatc(ncorner,ncorner,*)
       real *8 xmat(nts,*),xsgnl(*),xsgnr(*)

       done = 1
       pi = atan(done)*4

c
cc      if the edge identity is the same, return 0
c
       if(iedge.eq.jedge) then

       do i=1,nts
       do j=1,nss

       xmat(i,j) = 0

       enddo
       enddo

       return

       endif

       do 1100 ii=1,nts

       i = lns(iedge)+ii-1

       do 1000 jj=1,nss

       j = lns(jedge)+jj-1


       rrr = (xs(i) - xs(j))**2 + (ys(i)-ys(j))**2
       rrn = (xs(i)-xs(j))*rnx(j) + (ys(i)-ys(j))*rny(j)

       xmat(ii,jj) = rrn/rrr/(2.0d0*pi)*sqrt(qwts(j)*qwts(i))

 1000  continue
 1100  continue


c
cc       now check to see if any of the matrices need to be 
c        replaced with appropriate corner matrices
c
       do icint = 1,ncint


       if(iedge.eq.icl(icint).and.jedge.eq.icr(icint)) then
       if(icsgnl(icint).eq.0) its = 0
       if(icsgnl(icint).eq.1) its = nts-ncorner

       if(icsgnr(icint).eq.0) iss = 0
       if(icsgnr(icint).eq.1) iss = nss-ncorner

       call xreplmat(ncorner,ncorner,nts,its,iss,
     1    xmatc(1,1,ixmatc(icint)),xmat,xsgnl(icint))

       endif

       if(iedge.eq.icr(icint).and.jedge.eq.icl(icint)) then

       if(icsgnl(icint).eq.0) iss = 0
       if(icsgnl(icint).eq.1) iss = nss-ncorner

       if(icsgnr(icint).eq.0) its = 0
       if(icsgnr(icint).eq.1) its = nts-ncorner


       call xreplmat(ncorner,ncorner,nts,its,iss,
     1    xmatc(1,1,ixmatc(icint)),xmat,xsgnr(icint))


       endif

 1111 continue       

       enddo

       return
       end
c-------------------------------------------------------

       subroutine multaslow(n,x,y,amat)
       implicit real *8 (a-h,o-z)
c
cc        computes the slow matvec Ax=y
c
       

       real *8 x(*),y(*),amat(n,*)

       do i=1,n
       y(i) = 0

       do j=1,n

       y(i) = y(i) + amat(i,j)*x(j)

       enddo
       enddo

       return
       end
c-----------------------------------
        subroutine genpoterrs(xsamps,ysamps,nx,s1d,s2d,s1n,
     1      s2n,w,cx,cy,ncs,eps,x0,y0,ift,pots)
        implicit real *8 (a-h,o-z)
        dimension xsamps(nx),ysamps(nx),w(1),pots(nx),
     1      cx(1),cy(1)

        call mach_zero(epsz)

        done = 1
        do 2000 i=1,nx

        xt = xsamps(i)
        yt = ysamps(i)

        call prin2('xt = *',xt,1)
        call prin2('yt = *',yt,1)

        call polysign(cx,cy,ncs,xt,yt,vmin)
        call prin2('vmin = *',vmin,1)
c
        if (vmin .gt. eps) then
        call prinf('here = *',0,0)
        ned = 0
        nes = 0
        sfd = s1d
        sfs = s1n
        vte = 0
        if (ift .eq. 1) then
        rx = x0-xt
        ry = y0-yt
        r  = sqrt(rx*rx+ry*ry)
        vte = log(r)
        end if
        call evalfield(xt,yt,ned,nes,sfd,sfs,w,pots(i))
        call prin2('pots(i) = *',pots(i),1)
        pots(i) = (pots(i)+vte)
        end if

        if (vmin .lt. eps) then
        call prinf('there = *',0,0)
        ned = 0
        nes = 0
        sfd = s2d
        sfs = s2n
        call evalfield(xt,yt,ned,nes,sfd,sfs,w,pots(i))
        vte = 0
        if (ift .eq. 0) then
        rx = x0-xt
        ry = y0-yt
        r  = sqrt(rx*rx+ry*ry)
        vte = log(r)
        end if
        pots(i) = pots(i)-vte
        call prin2('pots(i) = *',pots(i),1)
        pots(i) = (pots(i))
        end if

        call prinf('######## *',0,0)

 2000 continue

ccc        call prin2('epsz = *',epsz,1)

        return
        end
c
c
c
c
c
        subroutine genpots(xsamps,ysamps,nx,ny,s1d,s2d,s1n,
     1      s2n,w,cx,cy,ncs,eps,x0,y0,ift,pots)
        implicit real *8 (a-h,o-z)
        dimension xsamps(nx,ny),ysamps(nx,ny),w(1),pots(nx,ny),
     1      cx(1),cy(1)

        call mach_zero(epsz)

        done = 1
        do 2000 i=1,nx
        do 1000 j=1,ny

        xt = xsamps(i,j)
        yt = ysamps(i,j)

ccc        call prin2('xt = *',xt,1)
ccc        call prin2('yt = *',yt,1)

        call polysign(cx,cy,ncs,xt,yt,vmin)
ccc        call prin2('vmin = *',vmin,1)
c
        if (vmin .gt. eps) then
        ned =-1
        nes =-1
        sfd = s1d
        sfs = s1n
        vte = 0
        if (ift .eq. 1) then
        rx = x0-xt
        ry = y0-yt
        r  = sqrt(rx*rx+ry*ry)
        vte = log(r)
        end if
        call evalfield(xt,yt,ned,nes,sfd,sfs,w,pots(i,j))
        pots(i,j) = pots(i,j)-vte
        end if

        if (vmin .lt. eps) then
        ned =-1
        nes =-1
        sfd = s2d
        sfs = s2n
        call evalfield(xt,yt,ned,nes,sfd,sfs,w,pots(i,j))
        pots(i,j) = pots(i,j)+0*vte
        pots(i,j) = pots(i,j)
        end if

 1000 continue
 2000 continue

        call prin2('epsz = *',epsz,1)

        return
        end
c
c
c
c
c
        subroutine genmesh(xmi,xma,ymi,yma,nx,ny,xmas,ymas)
        implicit real *8 (a-h,o-z)
        dimension xmas(nx,ny), ymas(nx,ny)

        done = 1

        do 2000 i=1,nx
        do 1000 j=1,ny

        xt = (i-1)*(xma-xmi)/((nx-1)*done)+xmi
        yt = (j-1)*(yma-ymi)/((ny-1)*done)+ymi
        
        xmas(i,j) = xt
        ymas(i,j) = yt

 1000 continue
 2000 continue

        return
        end
c
c
c
c
c
        subroutine polysign(cx,cy,ncs,xt,yt,vmin)
        implicit real *8 (a-h,o-z)
        dimension cx(1),cy(1)

        call pnpoly(xt,yt,cx,cy,ncs,iot)      
ccc        call prinf('iot = *',iot,1)


        do 1000 i=1,ncs
        x0 = cx(i)
        y0 = cy(i)

        if (i .eq. ncs) then
        x1 = cx(1)
        y1 = cy(1)
        end if

        if (i .lt. ncs) then
        x1 = cx(i+1)
        y1 = cy(i+1)
        end if

        call distpt2line(xt,yt,x0,y0,x1,y1,vtemp)
        if (i .eq. 1) then
            vmin = vtemp
        end if
        if (abs(vtemp) .lt. abs(vmin)) then
            vmin = vtemp
        end if
 1000 continue

        if (iot .eq. -1) then
            vmin = -abs(vmin)
        end if
        if (iot .eq.  1) then
            vmin = abs(vmin)
        end if

        return
        end
c
c
c
c
c
        subroutine distpt2line(xt,yt,x0,y0,x1,y1,val)
        implicit real *8 (a-h,o-z)
       
        rx0 = xt - x0
        ry0 = yt - y0
        rx1 = xt - x1
        ry1 = yt - y1

        r0  = sqrt(rx0*rx0+ry0*ry0)
        r1  = sqrt(rx1*rx1+ry1*ry1)

        rm  = r0
        if (r1 .lt. r0) rm = r1

        sy = x1-x0
        sx = y0-y1
ccc        call prin2('xs = *',sx,1)
ccc        call prin2('ys = *',sy,1)
        sr = sqrt(sx*sx+sy*sy)
        sx = sx/sr
        sy = sy/sr


ccc        call prin2('xt = *',xt,1)
ccc        call prin2('yt = *',yt,1)
ccc        call prin2('x0 = *',x0,1)
ccc        call prin2('y0 = *',y0,1)
ccc        call prin2('x1 = *',x1,1)
ccc        call prin2('y1 = *',y1,1)
ccc        call prin2('sx = *',sx,1)
ccc        call prin2('sy = *',sy,1)


        rdot = rx0*rx1+ry0*ry1
ccc        call prin2('rdot = *',rdot,1)
c
c
        if (rdot .gt. 0) then
        sig =-1
        if ((sx*rx0+sy*ry0) .gt. 0) sig = 1
        val = sig*rm
        return
        end if
c
c

        rv = rx0*sx+ry0*sy
ccc        call prin2('rv = *',rv,1)
        rv = rx1*sx+ry1*sy
ccc        call prin2('rv = *',rv,1)
        val = rv

ccc       stop

        return
        end
c
c
c
c
c
        subroutine evalfield(xt,yt,ned,nes,sfd,sfs,w,val)
        implicit real *8 (a-h,o-z)
        dimension w(1)

        nxs     = w(1)
        npan    = w(2)
        k       = w(3)
        lenp    = w(4)
        isn     = w(7)
        isd     = w(8)
        ipt     = w(9)
        its     = w(10)
        izv     = w(11)

        sf = 1

        call evsgllyracc(w(isn),w(ipt),w(izv),nxs,nes,w(its),npan,
     1      k,lenp,xt,yt,sf,val1)


        call evdblyracc(w(isd),w(ipt),w(izv),nxs,ned,w(its),npan,
     1      k,lenp,xt,yt,sf,val2)

        val = val1*sfs+val2*sfd

        return
        end
c
c
c
c
c
        subroutine wrapevals(sn,sd,pts,izs,nxs,ts,npan,k,lenp,
     1      w)
        implicit real *8 (a-h,o-z)
        dimension sn(1),sd(1),pts(6,1),izs(1),ts(6,1),w(1)

        w(1)    = nxs
        w(2)    = npan
        w(3)    = k
        w(4)    = lenp

        isn     = 50
        lsn     = nxs
        call larrmove(sn,w(isn),lsn)
        w(7)    = isn

        isd     = isn+lsn+50
        lsd     = nxs
        call larrmove(sd,w(isd),lsd)
        w(8)    = isd

        ipt     = isd+lsd+50
        lpt     = nxs*6
        call larrmove(pts,w(ipt),lpt)
        w(9)    = ipt

        its     = ipt+lpt+50
        lts     = (npan+1)*6
        call larrmove(ts,w(its),lts)
        w(10)   = its

        izv     = its+lts+50
        lzv     = nxs
        call larrmove(izs,w(izv),lzv)
        w(11)   = izv

        return
        end
c
c
c
c
c
        subroutine evaledgedd(xt,yt,hx,hy,ts,npan,net,rint)
        implicit real *8 (a-h,o-z)
        dimension ts(6,1),tars(2),vecs(2)
        external fddint

        iflag = 0
        do 2000 i=1,npan
        ip = ts(2,i)
        if (iflag .eq. 0) then
        if (ip .eq. net) then
        iflag = 1
        i0    = i
        end if
        end if
        if (iflag .eq. 1) then
        if (ip .ne. net) then
        i1 = i
        goto 2010
        end if
        end if
 2000 continue
 2010 continue

ccc        call prin2('ts = *',ts,npan*6)

        tars(1) = xt
        tars(2) = yt
        vecs(1) = hx
        vecs(2) = hy

        a = ts(1,i0)
        b = ts(1,i1)
ccc        call prin2('a = *',a,1)
ccc        call prin2('b = *',b,1)
        m =  63
        eps = 1.0d-14
        maxrec = 1000
        call adapgaus_new(ier,a,b,fddint,tars,vecs,m,eps,
     1      rint,maxrec,numint)

ccc        call prinf('numint = *',numint,1)
        
        return
        end
c
c
c
c
c
        subroutine evaledgesd(xt,yt,hx,hy,ts,npan,net,rint)
        implicit real *8 (a-h,o-z)
        dimension ts(6,1),tars(2),vecs(2)
        external fsdint

        iflag = 0
        do 2000 i=1,npan
        ip = ts(2,i)
        if (iflag .eq. 0) then
        if (ip .eq. net) then
        iflag = 1
        i0    = i
        end if
        end if
        if (iflag .eq. 1) then
        if (ip .ne. net) then
        i1 = i
        goto 2010
        end if
        end if
 2000 continue
 2010 continue

ccc        call prin2('ts = *',ts,npan*6)

        tars(1) = xt
        tars(2) = yt
        vecs(1) = hx
        vecs(2) = hy

        a = ts(1,i0)
        b = ts(1,i1)
ccc        call prin2('a = *',a,1)
ccc        call prin2('b = *',b,1)
        m =  63
        eps = 1.0d-14
        maxrec = 1000
        call adapgaus_new(ier,a,b,fsdint,tars,vecs,m,eps,
     1      rint,maxrec,numint)

ccc        call prinf('numint = *',numint,1)
        
        return
        end
c
c
c
c
c
        subroutine evaledged(xt,yt,ts,npan,net,rint)
        implicit real *8 (a-h,o-z)
        dimension ts(6,1)
        external fdint

        iflag = 0
        do 2000 i=1,npan
        ip = ts(2,i)
        if (iflag .eq. 0) then
        if (ip .eq. net) then
        iflag = 1
        i0    = i
        end if
        end if
        if (iflag .eq. 1) then
        if (ip .ne. net) then
        i1 = i
        goto 2010
        end if
        end if
 2000 continue
 2010 continue

ccc        call prin2('ts = *',ts,npan*6)

        a = ts(1,i0)
        b = ts(1,i1)
ccc        call prin2('a = *',a,1)
ccc        call prin2('b = *',b,1)
        m =  63
        eps = 1.0d-14
        maxrec = 1000
        call adapgaus_new(ier,a,b,fdint,xt,yt,m,eps,
     1      rint,maxrec,numint)

ccc        call prinf('numint = *',numint,1)
        
        return
        end
c
c
c
c
c
        subroutine evaledge(xt,yt,ts,npan,net,rint)
        implicit real *8 (a-h,o-z)
        dimension ts(6,1)
        external fsint

        iflag = 0
        do 2000 i=1,npan
        ip = ts(2,i)
        if (iflag .eq. 0) then
        if (ip .eq. net) then
        iflag = 1
        i0    = i
        end if
        end if
        if (iflag .eq. 1) then
        if (ip .ne. net) then
        i1 = i
        goto 2010
        end if
        end if
 2000 continue
 2010 continue

ccc        call prin2('ts = *',ts,npan*6)

        a = ts(1,i0)
        b = ts(1,i1)
ccc        call prin2('a = *',a,1)
ccc        call prin2('b = *',b,1)
        m =  23
        eps = 1.0d-15
        maxrec = 1000
        call adapgaus_new(ier,a,b,fsint,xt,yt,m,eps,
     1      rint,maxrec,numint)

        call prinf('numint = *',numint,1)
        
        return
        end
c
c
c
c
c
        subroutine init_fsint(tdata,npan,k,izi,data,npts,cx,cy,nc,
     1          alens,ifext,ws)
        implicit real *8 (a-h,o-z)
        dimension tdata(6,1),izi(1),data(1),tloc(6,1000),izl(90000),
     1          dloc(100 000),clx(1000),cly(1000),alens(1),aloc(1000),
     2          cx(1),cy(1),ws(1),tars(1),vecs(1)

        save

        ifl = ifext

        do 1500 i=1,(npan+1)
        do 1000 j=1,6
        tloc(j,i) = tdata(j,i)
 1000 continue
 1500 continue

        do 2000 i=1,npts
        dloc(i) = data(i)/sqrt(ws(i))
        izl(i)  = izi(i)
 2000 continue

        do 3000 i=1,nc

        clx(i)  = cx(i)
        cly(i)  = cy(i)
        aloc(i) = alens(i)

 3000 continue        

        nloc  = npts
        ncloc = nc
        nploc = npan
        kloc  = k



        return

        entry fsint(t,x0,y0,val)

        call getcoords(xb,yb,ub,vb,t,aloc,clx,cly,ncloc)
        call interpbndry(t,tloc,nploc,kloc,izl,dloc,vc)

        rx = x0-xb
        ry = y0-yb

        r = sqrt(rx*rx+ry*ry)
        val = log(r) 
        if (ifl .eq. 1) val = (-val+1)*vc
        if (ifl .eq. 0) val = (val)*vc


        return


        entry fsdint(t,tars,vecs,val)

        x0 = tars(1)
        y0 = tars(2)
        hx = vecs(1)
        hy = vecs(2)

        call getcoords(xb,yb,ub,vb,t,aloc,clx,cly,ncloc)
        call interpbndry(t,tloc,nploc,kloc,izl,dloc,vc)

        rx = x0-xb
        ry = y0-yb

        r = sqrt(rx*rx+ry*ry)
        val = (rx*hx+ry*hy)/(r*r)
        val =-(val)*vc

        return


        end
c
c
c
c
c
        subroutine init_fint(tdata,npan,k,izi,data,npts,cx,cy,nc,
     1          alens,ifext,ws)
        implicit real *8 (a-h,o-z)
        dimension tdata(6,1),izi(1),data(1),tloc(6,1000),izl(90000),
     1          dloc(100 000),clx(1000),cly(1000),alens(1),aloc(1000),
     2          cx(1),cy(1),ws(1),tars(1),vecs(1)

        save

        ifl = ifext

        do 1500 i=1,(npan+1)
        do 1000 j=1,6
        tloc(j,i) = tdata(j,i)
 1000 continue
 1500 continue

        do 2000 i=1,npts
        dloc(i) = data(i)/sqrt(ws(i))
        izl(i)  = izi(i)
 2000 continue

        do 3000 i=1,nc

        clx(i)  = cx(i)
        cly(i)  = cy(i)
        aloc(i) = alens(i)

 3000 continue        

        nloc  = npts
        ncloc = nc
        nploc = npan
        kloc  = k

        return

        entry fdint(t,x0,y0,val)


        call getcoords(xb,yb,ub,vb,t,aloc,clx,cly,ncloc)
        call interpbndry(t,tloc,nploc,kloc,izl,dloc,vc)

        rx = x0-xb
        ry = y0-yb

        r = sqrt(rx*rx+ry*ry)
        val = (rx*ub+ry*vb)/(r*r)
        if (ifl .eq. 1) val = (-val+1)*vc
        if (ifl .eq. 0) val = (-val)*vc

        return


        entry fddint(t,tars,vecs,val)

        x0 = tars(1)
        y0 = tars(2)
        hx = vecs(1)
        hy = vecs(2)

        call getcoords(xb,yb,ub,vb,t,aloc,clx,cly,ncloc)
        call interpbndry(t,tloc,nploc,kloc,izl,dloc,vc)

        rx = x0-xb
        ry = y0-yb

        r   = sqrt(rx*rx+ry*ry)
        t1  = (ub*rx+vb*rx)/(r*r)
        t2  = -2*(ub*rx+vb*ry)*(hx*rx+hy*ry)/(r*r*r*r)
        val = t1+t2
        val = (val)*vc

        return

        end
c
c
c
c
c        
        subroutine getcoords(x,y,hx,hy,t,alens,cx,cy,nc)
        implicit real *8 (a-h,o-z)
        dimension alens(1),cx(1),cy(1)

        tlen = 0

        do 3000 i=1,nc
        
        tlen = tlen + alens(i)

        if (t .lt. tlen) then
        cx0 = cx(i)
        cy0 = cy(i)

        if (i .lt. nc) then
        cx1 = cx(i+1)
        cy1 = cy(i+1)
        else
        cx1 = cx(1)
        cy1 = cy(1)
        end if

        toff = tlen - t
        vx = cx0-cx1
        vy = cy0-cy1
        r = sqrt(vx*vx+vy*vy)
        vx = vx/r
        vy = vy/r
        x = cx1+toff*vx
        y = cy1+toff*vy
        hx = vy
        hy =-vx
        goto 3010
        end if

 3000 continue
 3010 continue        

        if (t .lt. 0) then
        cx1 = cx(1)
        cy1 = cy(1)
        cx0 = cx(nc)
        cy0 = cy(nc)
        toff =-t
        vx = cx0 - cx1
        vy = cy0 - cy1
        r = sqrt(vx*vx+vy*vy)
        vx = vx/r
        vy = vy/r
        x = cx1+toff*vx
        y = cy1+toff*vy
        hx =  vy
        hy = -vx
        end if

        return
        end
c
c
c
c
c
        subroutine interpbndry(tval,tdata,npan,k,izi,data,val)
        implicit real *8 (a-h,o-z)
        dimension tdata(6,1),data(1),xvs(1000),wvs(1000),izi(1),
     1          u(100 000), v(100 000),vtemp(1000),coefs(1000),
     2          xc(1000),wc(1000)


        ifirst = 1

        if (ifirst .eq. 1) then
                itype = 2
                call legeexps(itype,k,xvs,u,v,wvs)
                ifirst = 0
        end if

        do 1000 i=1,npan+1
        if (tval .lt. tdata(1,i)) then
                ipan = i-1
                goto 1010
        end if
 1000 continue
 1010 continue

ccc        call prin2('tdata =*',tdata,6*(npan+1))
ccc
ccc        call prin2('tval = *',tval,1)
ccc        call prinf('ipan = *',ipan,1)

        t0 = tdata(1,ipan)
        t1 = tdata(1,ipan+1)

        i0 = izi(tdata(3,ipan))

        if (tdata(2,ipan) .lt. 0) then

        trel = (tval-t0)/(t1-t0)*2 -1
        if (abs(trel) .gt. 1) then
                call prinf('bad interp *',0,0)
                call prin2('trel = *',trel,1)
                call prin2('tval = *',tval,1)
                call prin2('t0 = *',t0,1)
                call prin2('t1 = *',t1,1)
ccc                stop
        end if

        call lapmatvec(u,data(i0),coefs,k)
        call legeexev(trel,val,coefs,k-1)

        end if

        if (tdata(2,ipan) .gt. 0) then

        if ((t1-tval) .le. (tval - t0)) then
        do 4000 i=1,31
        ind = izi(tdata(3,ipan)+31+i-1)
        vtemp(i) = data(ind)
 4000 continue
        trel = (tval-t0)/(t1-t0)-0.5d0
        if (abs(trel-0.25d0) .gt. 0.25d0) then
                call prinf('bad interp *',0,0)
                call prin2('trel = *',trel,1)
                call prin2('tval = *',tval,1)
                call prin2('t0 = *',t0,1)
                call prin2('t1 = *',t1,1)
                stop
        end if
        else
        do 4010 i=1,31
        ind = izi(tdata(3,ipan)+i-1)
        vtemp(i) = data(ind)
 4010 continue        
        trel = (t1-tval)/(t1-t0)-0.5d0
        if (abs(trel-0.25d0) .gt. 0.25d0) then
                call prinf('bad interp *',0,0)
                call prinf('bad interp *',0,0)
                call prin2('trel = *',trel,1)
                call prin2('tval = *',tval,1)
                call prin2('t0 = *',t0,1)
                call prin2('t1 = *',t1,1)
                call prin2('tdata =*',tdata,6*npan)
                stop
        end if

        end if

        call lapinterps(trel,vtemp,val)

        if (trel .lt. 1.0d-15) val = 0

        end if

        return
        end
c
c
c
c
c        
        subroutine labubble(x,n,m,ind)
        implicit real *8 (a-h,o-z)
        dimension x(n,m)

        do 2000 i=1,m-1
        iswap = 0 

        do 1500 j=1,m-1

        if (x(ind,j) .gt. x(ind,j+1)) then
        do 1000 k=1,n
        xt = x(k,j)
        x(k,j) = x(k,j+1)
        x(k,j+1) = xt
 1000 continue      
        iswap = 1
        end if  

 1500 continue

        if (iswap .eq. 0) goto 2500

 2000 continue

 2500 continue        

        return
        end
c
c
c
c
c
        subroutine evdblyracc(sd,pts,izs,npts,ned,ts,npan,
     1      k,lenp,xt,yt,sf,vout)
        implicit real *8 (a-h,o-z)
        dimension pts(6,1),sd(1),izs(1),xs(1),ys(1),ws(1),
     1      vx(1),vy(1),xr(500 000),yr(500 000),wr(500 000),
     1      vr(500 000),ur(500 000),sr(500 000),coefv(500 000),
     1      sdt(500 000),ts(6,1)
        external fdint
c
c
c
c

ccc        call prinf('ned = *',ned,1)

ccc        if (ned .eq. 0) then
ccc
ccc        a = ts(1,1)
ccc        b = ts(1,npan+1)
ccc        m =  62
ccc        eps = 1.0d-14
ccc        maxrec = 1000
ccc        call adapgaus_new(ier,a,b,fdint,xt,yt,m,eps,
ccc     1      rint,maxrec,numint)
ccc
ccc        vout = rint
cccccc        call prinf('numint = *',numint,1)
ccc        return
ccc
ccc        end if


        val = 0

        do 1000 ii=1,npts

        i = izs(ii)
        if (abs(pts(6,i)-ned) .gt. 0.5d0) then
        rx = pts(1,i) - xt
        ry = pts(2,i) - yt
        r  = sqrt(rx*rx+ry*ry)
        vt = rx*pts(4,i)+ry*pts(5,i)
        vt = vt/(r*r)
        vt = vt*sqrt(pts(3,i))*sf*sd(ii)
        val = val + vt
        end if
 1000 continue

        vp = 0

        if (ned .ne. 0) then
        call evaledged(xt,yt,ts,npan,ned,vp)
        end if

        vout = val+ vp

        return
        end
c
c
c
c
c
        subroutine evsgllyracc(sd,pts,izs,npts,ned,ts,npan,
     1      k,lenp,xt,yt,sf,vout)
        implicit real *8 (a-h,o-z)
        dimension pts(6,1),sd(1),izs(1),xs(1),ys(1),ws(1),
     1      vx(1),vy(1),xr(500 000),yr(500 000),wr(500 000),
     1      vr(500 000),ur(500 000),sr(500 000),coefv(500 000),
     1      sdt(500 000),ts(6,1)
c
c
c
c
        val = 0

        do 1000 ii=1,npts

        i = izs(ii)
        if (abs(pts(6,i)-ned) .gt. 0.5d0) then
        rx = pts(1,i) - xt
        ry = pts(2,i) - yt
        r  = sqrt(rx*rx+ry*ry)
        vt = rx*pts(4,i)+ry*pts(5,i)
        vt = log(r)*sqrt(pts(3,i))*sf*sd(ii)
        val = val + vt
        end if
 1000 continue

ccc        call prin2('sd = *',sd,npts)
ccc        stop
ccc        call prinf('ned = *',ned,1)

        vp = 0

        if (ned .ne. 0) then 
        call evaledge(xt,yt,ts,npan,ned,vp)
        end if

        vout = val+ vp

        return
        end
c
c
c
c
c
        subroutine evdbllyr(sd,pts,izs,npts,xt,yt,sf,val)
        implicit real *8 (a-h,o-z)
        dimension pts(6,1),sd(1),izs(1)

        val = 0

        do 1000 ii=1,npts

        i = izs(ii)
        rx = pts(1,i) - xt
        ry = pts(2,i) - yt
        r  = sqrt(rx*rx+ry*ry)
        vt = rx*pts(4,i)+ry*pts(5,i)
        vt = vt/(r*r)*sqrt(pts(3,i))*sf*sd(ii)
        val = val + vt

 1000 continue

        return
        end
c
c
c
c
c
        subroutine evsgllyr(sd,pts,izs,npts,xt,yt,sf,val)
        implicit real *8 (a-h,o-z)
        dimension pts(6,1),sd(1),izs(1)

        val = 0

        do 1000 ii=1,npts

        i  = izs(ii)
        rx = pts(1,i) - xt
        ry = pts(2,i) - yt
        r  = sqrt(rx*rx+ry*ry)
        vt = log(r) 
        vt = vt*sqrt(pts(3,i))*sf*sd(ii)
        val = val + vt

 1000 continue

        return
        end
c
c
c
c
c
        subroutine evdbllyraccd(sd,pts,izs,npts,ned,ts,npan,
     1      k,lenp,xt,yt,hx,hy,sf,vout)
        implicit real *8 (a-h,o-z)
        dimension pts(6,1),sd(1),izs(1),xs(1),ys(1),ws(1),
     1      vx(1),vy(1),xr(500 000),yr(500 000),wr(500 000),
     1      vr(500 000),ur(500 000),sr(500 000),coefv(500 000),
     1      sdt(500 000),ts(6,1)
        external fddint
c
c

ccc        call prinf('ned = *',ned,1)

        if (ned .eq. 0) then

        a = ts(1,1)
        b = ts(1,npan+1)
        m =  62
        eps = 1.0d-14
        maxrec = 1000
        call adapgaus_new(ier,a,b,fddint,xt,yt,m,eps,
     1      rint,maxrec,numint)

        vout = rint
        call prinf('numint = *',numint,1)
        return

        end if


        val = 0

        do 1000 ii=1,npts

        i = izs(ii)
        if (abs(pts(6,i)-ned) .gt. 0.5d0) then
        rx = xt-pts(1,i)
        ry = yt-pts(2,i)
        vxx = pts(4,i)
        vyy = pts(5,i)
        r  = sqrt(rx*rx+ry*ry)
        vt = rx*hx+ry*hy
        va = hx*vxx+hy*vyy
        vb = rx*vxx+ry*vyy
        vt = (hx*vxx+hy*vyy)/(r*r)-2*vt*vb/(r*r*r*r)
        vt = vt*sqrt(pts(3,i))*sf*sd(ii)
        val = val + vt
        end if
 1000 continue

        call evaledgedd(xt,yt,hx,hy,ts,npan,ned,vp)
       
        vout = val+ vp

        return
        end
c
c
c
c
c
        subroutine evsglyraccd(sd,pts,izs,npts,ned,ts,npan,
     1      k,lenp,xt,yt,hx,hy,sf,vout)
        implicit real *8 (a-h,o-z)
        dimension pts(6,1),sd(1),izs(1),xs(1),ys(1),ws(1),
     1      vx(1),vy(1),xr(500 000),yr(500 000),wr(500 000),
     1      vr(500 000),ur(500 000),sr(500 000),coefv(500 000),
     1      sdt(500 000),ts(6,1)
        external fsdint
c
c

ccc        call prinf('ned = *',ned,1)

        if (ned .eq. 0) then

        a = ts(1,1)
        b = ts(1,npan+1)
        m =  62
        eps = 1.0d-14
        maxrec = 1000
        call adapgaus_new(ier,a,b,fsdint,xt,yt,m,eps,
     1      rint,maxrec,numint)

        vout = rint
        call prinf('numint = *',numint,1)
        return

        end if


        val = 0

        do 1000 ii=1,npts

        i = izs(ii)
        if (abs(pts(6,i)-ned) .gt. 0.5d0) then
        rx = xt-pts(1,i)
        ry = yt-pts(2,i)
        r  = sqrt(rx*rx+ry*ry)
        vt = rx*hx+ry*hy
        vt = vt/(r*r)
        vt = vt*sqrt(pts(3,i))*sf*sd(ii)
        val = val + vt
        end if
 1000 continue

        call evaledgesd(xt,yt,hx,hy,ts,npan,ned,vp)
       
        vout = val+ vp

        return
        end
c
c
c
c
c
        subroutine gendirsrcdip(pts,izs,npts,xa,ya,hxa,hya,xb,yb,
     1      hxb,hyb,w1,w2,sf,vals)
        implicit real *8 (a-h,o-z)
        dimension pts(6,1),vals(1),izs(1)

        do 1000 ii=1,npts

        i = izs(ii)
        xi = pts(1,i)
        yi = pts(2,i)
        wi = pts(3,i)
        vx = pts(4,i)
        vy = pts(5,i)

        rx = xi - xa
        ry = yi - ya
        r  = sqrt(rx*rx+ry*ry)
        tp = (rx*hxa+ry*hya)/(r*r)
        tp = log(r)

        rx = xi - xb
        ry = yi - yb
        r  = sqrt(rx*rx+ry*ry)
        tp2 = (rx*hxb+ry*hyb)/(r*r)
        tp2 = log(r)
        vals(ii) = sf*(tp*w1+tp2*w2)*sqrt(wi)

 1000 continue

        return
        end
c
c
c
c
c
        subroutine gendirsrc(pts,izs,npts,x0,y0,sf,vals)
        implicit real *8 (a-h,o-z)
        dimension pts(6,1),vals(1),izs(1)

        do 1000 ii=1,npts

        i = izs(ii)
        xi = pts(1,i)
        yi = pts(2,i)
        wi = pts(3,i)
        vx = pts(4,i)
        vy = pts(5,i)

        rx = xi - x0
        ry = yi - y0
        r  = sqrt(rx*rx+ry*ry)
        vals(ii) = sf*log(r)*sqrt(wi)

 1000 continue

        return
        end
c
c
c
c
c
        subroutine genneusrcdip(pts,izs,npts,xa,ya,hxa,hya,xb,
     1      yb,hxb,hyb,w1,w2,sf,vals)
        implicit real *8 (a-h,o-z)
        dimension pts(6,1),vals(1),izs(1)

        do 1000 ii=1,npts

        i = izs(ii)
        xi = pts(1,i)
        yi = pts(2,i)
        wi = pts(3,i)
        vx = pts(4,i)
        vy = pts(5,i)



        rx = xi - xa
        ry = yi - ya
        r  = sqrt(rx*rx+ry*ry)
        t1  = (vx*rx+vy*ry)/(r*r)
ccc        t2  = -2*(vx*rx+vy*ry)*(hxa*rx+hya*ry)/(r*r*r*r)
        val1 = w1*(t1)

        rx = xi - xb
        ry = yi - yb
        r  = sqrt(rx*rx+ry*ry)
        t1  = (vx*rx+vy*ry)/(r*r)
ccc        t2  = -2*(vx*rx+vy*ry)*(hxb*rx+hyb*ry)/(r*r*r*r)
        val2 = w2*(t1)

        vt = sf*(val1+val2)*sqrt(wi)
        vals(ii) = vt 

 1000 continue

        return
        end
c
c
c
c
c
        subroutine genneusrc(pts,izs,npts,x0,y0,sf,vals)
        implicit real *8 (a-h,o-z)
        dimension pts(6,1),vals(1),izs(1)

        do 1000 ii=1,npts

        i = izs(ii)
        xi = pts(1,i)
        yi = pts(2,i)
        wi = pts(3,i)
        vx = pts(4,i)
        vy = pts(5,i)

        rx = xi - x0
        ry = yi - y0
        r  = sqrt(rx*rx+ry*ry)
        vt = sf*(rx*vx+ry*vy)/(r*r)*sqrt(wi)
        vals(ii) = vt 

 1000 continue

        return
        end
c
c
c
c
c
        subroutine tmatmult(a,x,y,n)
        implicit real *8 (a-h,o-z)
        dimension a(n,n),x(1),y(1)

        do 2000 i=1,n
        y(i) = 0
        do 1000 j=1,n

        y(i) = y(i)+a(i,j)*x(j)

 1000 continue
 2000 continue

        return
        end
c
c
c
c
c
        subroutine getl2(vec,n,val)
        implicit real *8 (a-h,o-z)
        dimension vec(1)
 
        val = 0

        do 1000 i=1,n
        val = val+vec(i)**2
 1000 continue

        return
        end
c
c
c
c
c
        subroutine evpotents(xt,yt,ts,npan,rint)
        implicit real *8 (a-h,o-z)
        dimension ts(6,1)
        external fdint

        a = ts(1,1)
        b = ts(1,npan+1)
        m =  62
        eps = 1.0d-14
        maxrec = 1000
        call adapgaus_new(ier,a,b,fdint,xt,yt,m,eps,
     1      rint,maxrec,numint)

        call prinf('numint = *',numint,1)
        
        return
        end
c
c
c
c
c
        subroutine evalpotential(pts,npts,izs,dens,xt,yt,val)
        implicit real *8 (a-h,o-z)
        dimension pts(6,1),izs(1),dens(1)

        val = 0

        do 2000 i=1,npts
        ind = izs(i)
        xs = pts(1,ind)
        ys = pts(2,ind)
        ws = pts(3,ind)
        vx = pts(4,ind)
        vy = pts(5,ind)

        call evpot(xs,ys,vx,vy,xt,yt,tval)
        val = val + (tval)*sqrt(ws)*dens(i)

 2000 continue

        return
        end
c
c
c
c
c
        subroutine gendipfield(src,npts,pts,izs,x0,y0,hx,hy)
        implicit real *8 (a-h,o-z)
        dimension src(1), pts(6,1), izs(1)

        do 1050 i=1,npts
        
        ii = izs(i)
        xs = pts(1,ii)
        ys = pts(2,ii)
        ws = pts(3,ii)

        call evpot(x0,y0,hx,hy,xs,ys,tval)
        src(i) = -tval*sqrt(ws)
ccc  	   src(i) = sqrt(ws)
 1050 continue

        return
        end
c
c
c
c
c
        subroutine genchrgfield(src,npts,pts,izs,x0,y0)
        implicit real *8 (a-h,o-z)
        dimension src(1),pts(6,1), izs(1)

        do 1050 i=1,npts

        ii = izs(i)
        xs = pts(1,ii)
        ys = pts(2,ii)
        ws = pts(3,ii)
        hx = pts(4,ii)
        hy = pts(5,ii)

        call evpot(x0,y0,hx,hy,xs,ys,tval)
        src(i) = tval*sqrt(ws)

 1050 continue

        return
        end
c
c
c
c
c
ccc        subroutine genfullkern(npts,work,xmat,w,lw)
ccc        implicit real *8 (a-h,o-z)
ccc        dimension work(1),xmat(npts,npts),w(1),s(100 000)
ccc
ccc        do 1010 i=1,npts
ccc        do 1000 j=1,npts
ccc
ccc        call getentry0(i,j,work,inds,xmat(i,j))
ccc
ccc 1000 continue
ccc 1010 continue
ccc
ccc        eps = 1.0d-14
ccc        call svdpivot2(ier,xmat,npts,npts,s,ncols,eps,
ccc     1      work,100 000 000,ltot)
ccc
ccc
ccc        call prin2('s = *',s,npts)
ccc        call prinf('ier = *',ier,1)
ccc        call prinf('npts = *',npts,1)
ccc
ccc        return
ccc        end
c
c
c
c
c
        subroutine ladddiag(xmat,nskel,jo,n)
        implicit real *8 (a-h,o-z)
        dimension xmat(nskel,nskel)

        do 2000 i=1,n
        do 1000 j=1,n

        ii = i+jo-1
        jj = j+jo-1

        if (i .eq. j) xmat(ii,jj) = 1
        if (i .ne. j) xmat(ii,jj) = 0

 1000 continue
 2000 continue

        return
        end
c
c
c
c
c
        subroutine laptransf(tmat,n,m,io,jo,xmat,nskel)
        implicit real *8 (a-h,o-z)
        dimension tmat(n,m), xmat(nskel,nskel)

        do 1010 i=1,n
        do 1000 j=1,m

        ii = i+io-1
        jj = j+jo-1
        xmat(ii,jj) = tmat(i,j)

 1000 continue
 1010 continue        

        return
        end        
c
c
c
c
c
        subroutine getentry0(i,j,a,b,w,val)
        implicit real *8 (a-h,o-z)
        dimension w(1)

        ipts = w(1)
        iamats = w(2)
        nr0 = w(3)
        iizs = w(4)
        ifext = w(5)

ccc        call prinf('ipts = *',ipts,1)
ccc        call prinf('iamats = *',iamats,1)
ccc        call prinf('nr0 = *',nr0,1)
ccc        call prinf('iizs = *',iizs,1)
ccc        call prinf('ifext = *',ifext,1)

ccc        call prin2('iizs = *',w(iizs),550)

        call getentry(i,j,a,b,w(ipts),w(iamats),nr0,w(iizs),
     1      ifext,val)

        return
        end
c
c
c
c
c
        subroutine getentry(i,j,a,b,pts,amats,nr0,izs,ifext,val)
        implicit real *8 (a-h,o-z)
        dimension pts(6,1),amats(nr0,1),izs(1)

        done = 1
        pi = atan(done)*4

ccc        call prinf('i = *',i,1)
ccc        call prinf('j = *',j,1)
       
ccc        call prinf('izs(1) = *',izs(1),1)
        
        wi = pts(3,izs(i))
        wj = pts(3,izs(j))

        val = 0

ccc        if (ifext .eq. 1) val = val +sqrt(wi)*sqrt(wj)

        if (i .eq. j) then
            val = val + a*pi
            return
        endif

ccc        call prinf('ii = *',izs(i),1)
ccc        call prinf('ij = *',izs(j),1)

        ii = pts(6,izs(i))
        ij = pts(6,izs(j))

        if (ii .lt. 0) then
        if (ij .eq. ii) then
ccc        val = 0
        return
        endif
        endif

        if (ii .ne. ij) then

ccc        if (ii .eq. -ij) then
ccc        return
ccc        endif

        xs = pts(1,izs(i))
        ys = pts(2,izs(i))
        ws = pts(3,izs(i))
        vxj = pts(4,izs(i))
        vyj = pts(5,izs(i))
        xt = pts(1,izs(j))
        yt = pts(2,izs(j))
        wt = pts(3,izs(j))
        vx = pts(4,izs(j))
        vy = pts(5,izs(j))

ccc        call prin2('xs = *',xs,1)
ccc        call prin2('ys = *',ys,1)
ccc        call prin2('ws = *',ws,1)
ccc        call prin2('xt = *',xt,1)
ccc        call prin2('yt = *',yt,1)
ccc        call prin2('wt = *',wt,1)
ccc        call prin2('vx = *',vx,1)
ccc        call prin2('vy = *',vy,1)
ccc        call prin2('vxj = *',vxj,1)
ccc        call prin2('vyj = *',vyj,1)
ccc
ccc        call prinf('ii = *',ii,1)
ccc        call prinf('ij = *',ij,1)

        call evpot(xt,yt,vx,vy,xs,ys,tval)
        val = val+b*tval*sqrt(wt)*sqrt(ws)
        return

        endif

        if (ii .eq. ij) then
        if (ii .gt. 0) then
        ncorn = ii
        iloc = mod(izs(i),nr0*2)
        if (iloc .eq. 0) iloc = 2*nr0
        jloc = mod(izs(j),nr0*2)
        if (jloc .eq. 0) jloc = 2*nr0

ccc        call prinf('iloc = *',iloc,1)
ccc        call prinf('jloc = *',jloc,1)

        if (iloc .le. nr0) then
        if (jloc .le. nr0) then
        return
        endif
        endif

        if (iloc .gt. nr0) then
        if (jloc .gt. nr0) then
        return
        endif
        endif

        jloc = mod(jloc,nr0)
        if (jloc .eq. 0) jloc = nr0
        iloc = mod(iloc,nr0)
        if (iloc .eq. 0) iloc = nr0

        nind = nr0*(ncorn-1)+jloc
        tval = amats(iloc,nind)
        val = val+b*tval

        endif
        endif
        
        return
        end
c
c
c
c
c
        subroutine larrmove(a,b,n)
        implicit real *8 (a-h,o-z)
        dimension a(1),b(1)

        do 1000 i=1,n
        b(i)=a(i)
 1000 continue

        return
        end
c
c
c
c
c
        subroutine lapkergen(xmat,nxs,xs,ys,vxs,vys,ws)
        implicit real *8 (a-h,o-z)
        dimension xmat(nxs,nxs),xs(1),ys(1),vxs(1),vys(1),ws(1)

        done = 1
        pi = atan(done)*4

        do 2002 i=1,nxs
        do 2001 j=1,nxs

        xmat(i,j) = 0
        if (i .ne. j) then
        call evpot(xs(j),ys(j),vxs(j),vys(j),xs(i),ys(i),tval)
        tval = tval*sqrt(ws(i))*sqrt(ws(j))
        endif
        if (i .eq. j) then
        tval = pi
        endif
        xmat(i,j) = tval

 2001 continue
 2002 continue

        return
        end
c
c
c
c
c
        subroutine lapkermod(xmat,nxs,ws)
        implicit real *8 (a-h,o-z)
        dimension xmat(nxs,nxs), ws(1)

        do 2000 i=1,nxs
        do 1000 j=1,nxs

        xmat(i,j) = xmat(i,j)+sqrt(ws(i))*sqrt(ws(j))

 1000 continue
 2000 continue

        return
        end
c
c
c
c
c
        subroutine lapmatmult(a,b,c,n,m,k)
        implicit real *8 (a-h,o-z)
        dimension a(n,m),b(m,k),c(n,k)

        do 3000 i=1,n
        do 2000 j=1,k
        val=0
        do 1000 ii=1,m
        val = val+a(i,ii)*b(ii,j)
 1000 continue
        c(i,j) = val
 2000 continue
 3000 continue

        return
        end
c
c
c
c
c
        subroutine lapmatmultt(a,b,c,n,m,k)
        implicit real *8 (a-h,o-z)
        dimension a(n,m),b(k,m),c(n,k)

        do 3000 i=1,n
        do 2000 j=1,k
        val = 0
        do 1000 ii=1,m
        val = val + a(i,ii)*b(j,ii)
 1000 continue
        c(i,j) = val
 2000 continue
 3000 continue

        return
        end
c
c
c
c
c
        subroutinelapmatmult2(a,b,c,n,m,k)
        implicit real *8 (a-h,o-z)
        dimension a(m,n),b(m,k),c(n,k)

        do 3000 i=1,n
        do 2000 j=1,k
        val = 0
        do 1000 ii=1,m
        val = val + a(ii,i)*b(ii,j)
 1000 continue
        c(i,j) = val
 2000 continue
 3000 continue

        return
        end
c
c
c
c
c
        subroutine evpot(x,y,vx,vy,x0,y0,val)
        implicit real *8 (a-h,o-z)
        
        rx = x-x0
        ry = y-y0
        top = (rx*vx+ry*vy)
        bot = rx*rx+ry*ry
        val =top/bot

        return
        end
c
c
c
c
c
        subroutine gendisc(nc,cx,cy,cout,pts,pars,
     1      alens,npts,nr0,izs,izi,ts,npan)
        implicit real *8 (a-h,o-z)
        dimension cx(1),cy(1),pts(6,1),pars(1),
     1      wcs(1),izs(1),alens(1),cout(6,1),xts(1000),wts(1000),
     2      igs(1000),izi(1),ts(6,1),spanl(100)

        cfrac = pars(1)
        lenp = pars(2)
        k = pars(3)
        dmin = pars(4)
ccc        nblev =pars(5)
        kbis = pars(5)

        call prin2('cfrac = *',cfrac,1)
        call prinf('k = *',k,1)
        call prinf('lenp = *',lenp,1)
        call prin2('dmin = *',dmin,1)

        done = 1
        pi = atan(done)*4


ccc              .  .  .   computing lengths

        do 1000 i=1,(nc-1)

        dx = cx(i+1)-cx(i)
        dy = cy(i+1)-cy(i)
        alens(i) = sqrt(dx*dx+dy*dy)

 1000 continue
    
        dx = cx(1) - cx(nc)
        dy = cy(1) - cy(nc)
        alens(nc) = sqrt(dx*dx+dy*dy)


ccc               .  .  .   computing corner panel lengths

        cout(1,1) = min(alens(1),alens(nc))*cfrac

        do 1001 i=2,(nc)

        a1 = alens(i-1)
        a2 = alens(i)
        cout(1,i) = min(a1,a2)*cfrac
ccc        if (i .eq. 2) cout(1,i) = min(a1,a2)*cfrac/6
 1001 continue        

        call prin2('cout = *',cout,6*nc)
ccc        stop

ccc             .  .  .     computing angles

        dx1 = cx(nc) - cx(1)
        dy1 = cy(nc) - cy(1)
        dx2 = cx(2) - cx(1)
        dy2 = cy(2) - cy(1)

        r1 = sqrt(dx1*dx1+dy1*dy1)
        r2 = sqrt(dx2*dx2+dy2*dy2)
        dx1 = dx1/r1
        dy1 = dy1/r1
        dx2 = dx2/r2
        dy2 = dy2/r2

        if ((dx1*dx2+dy1*dy2) .ge.  1) cout(2,1) = 0  
        if ((dx1*dx2+dy1*dy2) .le. -1) cout(2,1) = 1

        if ((dx1*dx2+dy1*dy2) .lt.  1) then
        if ((dx1*dx2+dy1*dy2) .gt. -1) then
            cout(2,1) = acos(dx1*dx2+dy1*dy2)/pi
            if (dy1*dx2-dy2*dx1 .lt. 0) cout(2,1) = 2-cout(2,1)
        endif
        endif

        cout(3,1) = dx2
        cout(4,1) = dy2

        dx1 = cx(nc-1) - cx(nc)
        dy1 = cy(nc-1) - cy(nc)
        dx2 = cx(1) - cx(nc)
        dy2 = cy(1) - cy(nc)

        r1 = sqrt(dx1*dx1+dy1*dy1)
        r2 = sqrt(dx2*dx2+dy2*dy2)
        dx1 = dx1/r1
        dy1 = dy1/r1
        dx2 = dx2/r2
        dy2 = dy2/r2



        if ((dx1*dx2+dy1*dy2) .ge.  1) cout(2,nc) = 0  
        if ((dx1*dx2+dy1*dy2) .le. -1) cout(2,nc) = 1

        if ((dx1*dx2+dy1*dy2) .lt.  1) then
        if ((dx1*dx2+dy1*dy2) .gt. -1) then
            cout(2,nc) = acos(dx1*dx2+dy1*dy2)/pi
            if (dy1*dx2-dy2*dx1 .lt. 0) cout(2,nc) = 2-cout(2,nc)
        endif
        endif

ccc        cout(2,nc) = acos(dx1*dx2+dy1*dy2)/pi
ccc        if (dy1*dx2-dy2*dx1 .lt. 0) cout(2,nc) = 2-cout(2,nc)


        cout(3,nc) = dx2
        cout(4,nc) = dy2


        do 1002 i=2,(nc-1)

        dx1 = cx(i-1) - cx(i)
        dy1 = cy(i-1) - cy(i)
        dx2 = cx(i+1) - cx(i)
        dy2 = cy(i+1) - cy(i)

        r1 = sqrt(dx1*dx1+dy1*dy1)
        r2 = sqrt(dx2*dx2+dy2*dy2)
        dx1 = dx1/r1
        dy1 = dy1/r1
        dx2 = dx2/r2
        dy2 = dy2/r2

ccc        call prin2('cos = *',dx1*dx2+dy1*dy2,1)
c


        if ((dx1*dx2+dy1*dy2) .ge.  1) cout(2,i) = 0  
        if ((dx1*dx2+dy1*dy2) .le. -1) cout(2,i) = 1

        if ((dx1*dx2+dy1*dy2) .lt.  1) then
        if ((dx1*dx2+dy1*dy2) .gt. -1) then
            cout(2,i) = acos(dx1*dx2+dy1*dy2)/pi
            if (dy1*dx2-dy2*dx1 .lt. 0) cout(2,i) = 2-cout(2,i)
        endif
        endif


ccc        cout(2,i) = acos(dx1*dx2+dy1*dy2)/pi
ccc        if (dy1*dx2-dy2*dx1 .lt. 0) cout(2,i) = 2-cout(2,i)


        cout(3,i) = dx2
        cout(4,i) = dy2

 1002 continue

ccc              .  .  .    generating corner panels


        itype = 0
        call lapdisc(xts,wts,uj,vj,nr0,itype)

        toff = 0
        indpt = 1

ccc         .  .  .  deal with the first/last corner separately

        xoff = cx(1)
        yoff = cy(1)
        cplen = cout(1,1)

        do 2910 ii=1,nr0

        pts(1,indpt) =-xts(ii)*cplen*2*cout(3,nc)+xoff 
        pts(2,indpt) =-xts(ii)*cplen*2*cout(4,nc)+yoff 
        pts(3,indpt) = wts(ii)*2*cplen
        pts(4,indpt) = -cout(4,nc)
        pts(5,indpt) =  cout(3,nc)

        ival = 1
        pts(6,indpt) = ival 
        indpt = indpt+1

 2910 continue

ccc         .  .  .  now do the other corners

        indt = 1
        al = 0

        do 3050 i=1,nc

        ts(1,indt) = al-cout(1,i)
        ts(2,indt) = i
        ts(3,indt) = indpt-31
        call prinf('indpt = *',indpt,1)

        al = al +alens(i)
        call prin2('alens = *',alens(i),1)
        call prin2('al = *',al,1)

        indt = indt+1

        xoff = cx(i)
        yoff = cy(i)

        do 3000 ii=1,nr0

        pts(1,indpt) = xts(ii)*cout(1,i)*2*cout(3,i)+xoff 
        pts(2,indpt) = xts(ii)*cout(1,i)*2*cout(4,i)+yoff 
        pts(3,indpt) = wts(ii)*2*cout(1,i)
        pts(4,indpt) =-cout(4,i)
        pts(5,indpt) = cout(3,i)
        pts(6,indpt) = i 

        indpt = indpt+1

 3000 continue

        if (i .lt. nc) then
        xoff = cx(i+1)
        yoff = cy(i+1)
        cplen = cout(1,i+1)
        endif

        if (i .eq. nc) then
        goto 3060
        endif

        do 3010 ii=1,nr0
        
        pts(1,indpt) =-xts(ii)*cplen*2*cout(3,i)+xoff 
        pts(2,indpt) =-xts(ii)*cplen*2*cout(4,i)+yoff 
        pts(3,indpt) = wts(ii)*2*cplen
        pts(4,indpt) = -cout(4,i)
        pts(5,indpt) =  cout(3,i)

        ival = i+1
        if (i+1 .gt. nc) ival = 1
        pts(6,indpt) = ival 

        indpt = indpt+1

 3010 continue

        toff = toff+alens(i)

 3050 continue
 3060 continue


ccc          .  .  .        generating gaussian panels

        call prinf('k=*',k,1)

        itype = 1
        call legeexps(itype,k,xts,uj,vj,wts)

        tlen = 0
        totlen = 0

        do 2010 i=1,nc

        igs(i) = indpt

        toff =cout(1,i)

        xoff =toff*cout(3,i)+cx(i)
        yoff =toff*cout(4,i)+cy(i)

        if (i .lt. nc) then
        d = alens(i)-cout(1,i)-cout(1,i+1)
        endif

        if (i .eq. nc) then
        d= alens(i)-cout(1,i)-cout(1,1)
        endif

        call prin2('d = *',d,1)

        dplen = d/(lenp+0.0d0)
        lenpan = lenp
        if (dplen .gt. dmin) then 
            lenpan = int(d/dmin)
            dplen = d/(lenpan+0.0d0)
        endif
        call prin2('dplen = *',dplen,1)
       
        call calclenbis(d,kbis,lenp,nped,spanl)
        call prin2('spanl= *',spanl,nped)

        tlen = tlen+toff
        if (i .gt. 1) tlen = tlen+toff

        do 2005 ii=1,nped
        dplen = spanl(ii)
        ts(1,indt) = tlen
        ts(2,indt) = -i
        ts(3,indt) = indpt
        call prin2('ts = *',ts(3,indt),1)
        tlen = tlen + dplen
        indt = indt+1

        do 2000 iii =1,k

        pts(1,indpt) = (xts(iii)+1)/2*dplen*cout(3,i)+xoff
        pts(2,indpt) = (xts(iii)+1)/2*dplen*cout(4,i)+yoff
        pts(3,indpt) = wts(iii)/2*dplen
        pts(4,indpt) =-cout(4,i)
        pts(5,indpt) = cout(3,i)
        pts(6,indpt) = -i

        indpt = indpt+1


 2000 continue

        xoff = xoff+cout(3,i)*dplen
        yoff = yoff+cout(4,i)*dplen

 2005 continue

        totlen = totlen+alens(i)

 2010 continue

        npts = indpt-1
ccc        call prinf('ngau = *',ngau,1)
ccc        call plotgauspt(ptgaus,ngau,25)



        iic = 1
        iiz = 1
        iig = 1
        igs(nc+1) = npts+1

        do 4000 i=1,nc
        do 3900 ii =1,nr0

        izs(iiz) = iic+nr0-ii
        iiz = iiz+1

 3900 continue

        iic = iic+nr0
       
        do 3910 ii =1,nr0

        izs(iiz) = iic
        iiz = iiz+1
        iic = iic+1

 3910 continue

        ngs = igs(i+1)-igs(i)

        do 3950 ii=1,ngs

        izs(iiz) = igs(i)+ii-1
        iiz = iiz+1

 3950 continue

 4000 continue

ccc        call prinf('izs = *',izs,iiz-1)
ccc        call prinf('nzs = *',iiz-1,1)

        do 4010 i=1,(iiz-1)

        izi(izs(i)) = i

 4010 continue


        wtot = 0

        do 5000 i=1,npts

        wtot = wtot+pts(3,i)

 5000 continue

        call prin2('total gauss weight=*',wtot,1)

        wtot2 = 0

        do 5001 i=1,nc

        wtot2 = wtot2+alens(i)

 5001 continue

        call prin2('tot. weight = *',wtot2,1)

        effm = wtot-wtot2
        call prin2_long('mass effect = *',effm,1)
        ts(1,indt) = tlen
        npan = indt-1        
           
        return
        end
c
c
c
c
c
        subroutine calclenbis(alen,nlb,nmid,npans,dlens)
        implicit real *8 (a-h,o-z)
        dimension dlens(1)

        atot = 0
        done = 1
        npans = nlb*2+nmid
        do 2000 i=1,nlb
        atot = atot+2/((2*done)**((nlb-i+1)*done))
 2000 continue
        atot = atot+nmid

        call prin2('atot = *',atot,1)
        call prin2('alen = *',alen,1)

        dlpa = alen/atot

        do 3000 i=1,nlb
        dlens(i) = dlpa/(2*done)**((nlb-i+1)*done)
        dlens(npans-i+1) = dlpa/(2*done)**((nlb-i+1)*done)
 3000 continue

        do 4000 i=1,nmid
        dlens(i+nlb) = dlpa
 4000 continue

        wtot = 0
        do 5000 i=1,npans
        wtot = wtot+dlens(i)
 5000 continue
        call prin2('wtot = *',wtot,1)


        return
        end
c
c
c
c
c
        subroutine getangs(cornmat,nc,angs)
        implicit real *8 (a-h,o-z)
        dimension cornmat(6,1),angs(1)


        do 1000 i=1,nc
        angs(i) = cornmat(2,i)

 1000 continue

        return
        end
c
c
c
c
c
c
ccc        call mergedisc(pts,npts,xs,ys,ws,nxs,vxs,vys,izs,izi)


        subroutine mergedisc(pts,npts,xs,ys,ws,
     1      nxs,vxs,vys,izs,izi)
        implicit real *8 (a-h,o-z)
        dimension pts(6,1),xs(1),ys(1),ws(1),vxs(1),
     1      vys(1),izs(1),izi(1)


        nxs = npts

        do 1000 i=1,nxs

        xs(i) = pts(1,izs(i))
        ys(i) = pts(2,izs(i))
        ws(i) = pts(3,izs(i))
        vxs(i) = pts(4,izs(i))
        vys(i) = pts(5,izs(i))

 1000 continue

        return
        end
c
c
c
c
c
        subroutine evalpot(dens,xs,ys,ws,xt,yt,vxs,vys,nxs,val)
        implicit real *8 (a-h,o-z)
        dimension dens(1),xs(1),ys(1),ws(1),vxs(1),vys(1)

        val = 0
        do 1000 i=1,nxs

        rx = xs(i)-xt
        ry = ys(i)-yt
        top = (rx*vxs(i)+ry*vys(i))*dens(i)*sqrt(ws(i))
        bot = rx*rx+ry*ry
        val =val +top/bot

ccc        call prin2('vxs = *',vxs(i),1)
ccc        call prin2('vys = *',vys(i),1)
ccc        call prin2('dens = *', dens(i),1)
ccc        call prin2('wt = *',ws(i),1)
ccc        call prin2('rx = *',rx,1)
ccc        call prin2('top = *',top,1)
ccc        call prin2('bot = *',bot,1)

 1000 continue

        return
        end
c
c
c
c
c
        subroutine lapmatvec(x,v,vn,n)
        implicit real *8 (a-h,o-z)
        dimension x(n,n),v(1),vn(1)

        do 2000 i=1,n
        val = 0
        do 1000 j=1,n
        val = val+x(i,j)*v(j)
 1000 continue
        vn(i) = val
 2000 continue
        
        return
        end
c
c
c
c
c
        subroutine genptsrc(xs,ys,ws,nxs,x0,y0,src)
        implicit real *8 (a-h,o-z)
        dimension xs(1),ys(1),ws(1),src(1)

        do 1000 i=1,nxs

        rx = xs(i)-x0
        ry = ys(i)-y0
        r = sqrt(rx**2+ry**2)
        src(i) = log(r)*sqrt(ws(i))
ccc        src(i) = sqrt(ws(i))
 1000 continue

        return
        end
c
c
c
c
c
        subroutine gendipsrc(xs,ys,ws,nxs,x0,y0,hx,hy,src)
        implicit real *8 (a-h,o-z)
        dimension xs(1),ys(1),ws(1),src(1)

        do 1000 i=1,nxs

        rx = xs(i)-x0
        ry = ys(i)-y0
        r = sqrt(rx**2+ry**2)
        src(i) = log(r)*sqrt(ws(i))
        src(i) = (rx*hx+ry*hy)/(r*r)*sqrt(ws(i))
ccc        src(i) = sqrt(ws(i))
 1000 continue

        return
        end
c
c
c
c
c
        subroutine updatemat(xmat,nx,amat,na,i1,i2,cpanl)
        implicit real *8 (a-h,o-z)
        dimension xmat(nx,nx),amat(na,na)

        dtwo = 2

        do 2000 i=1,na
        do 1000 j=1,na

        xmat(i1-1+i,i2-1+j) = amat(i,j)
        xmat(i2-1+i,i1-1+j) = amat(i,j)

 1000 continue
 2000 continue

        return
        end
c
c
c
c
c
        subroutine addcorkern(xmat,nx,amat,na,inda,alpha)
        implicit real *8 (a-h,o-z)
        dimension xmat(nx,nx),amat(na,na)

        done = 1
        pi = atan(done)*4

        call lapcornmat(alpha,amat,nr0)

        i1 = inda
        i2 = inda+nr0

        do 2000 i=1,na
        do 1000 j=1,na

        xmat(i1-1+i,i2-1+j) = amat(i,j)
        xmat(i2-1+i,i1-1+j) = amat(i,j)
 1000 continue

 2000 continue


        call prinf('na = *',na,1)

        return
        end
c
c
c
c
c
        subroutine getkernmats(amats,cormat,nc)
        implicit real *8 (a-h,o-z)
        dimension amats(1),cormat(6,1)

        ind = 1

        do 3000 in =1,nc

        alpha = cormat(2,in)
        call lapcornmat(alpha,amats(ind),nr0)
        ind = ind+nr0*nr0

 3000 continue

        return
        end

c
c
c
c
c
        subroutine updatekern(xmat,nx,amat,na,cormat,nc,ngau,nxs)
        implicit real *8 (a-h,o-z)
        dimension xmat(nx,nx),amat(na,na),cormat(4,1)

        do 3000 in=1,nc

        alpha = cormat(2,in)
        call lapcornmat(alpha,amat,nr0)

ccc     this might be morally incorrect while still working

        if (in .eq. 1) i2 = (nxs-na+1)
        if (in .gt. 1) i2 = ngau+2*(in-1)*na-na+1

        i1 = ngau+2*(in-1)*na+1
         
        do 2000 i=1,na
        do 1000 j=1,na

        xmat(i1-1+i,i2-1+j) = amat(i,j)
        xmat(i2-1+i,i1-1+j) = amat(i,j)

 1000 continue
 2000 continue

 3000 continue

        call prinf('na = *',na,1)

        return
        end
c
c
c
c
c

        subroutine modkern(xmat,xs,ys,ws,vxs,vys,wt,nxt,nxs,itype)
        implicit real *8 (a-h,o-z)
        dimension xmat(nxt,nxs),xs(1),ys(1),ws(1),vxs(1),vys(1),wt(1)

        do 2000 i=1,nxt
        do 1000 j=1,nxs

        if (itype .eq. 1) then
            xmat(i,j) = xmat(i,j) + sqrt(ws(j))*sqrt(wt(i))
        endif
        if (itype .eq. 2) then
            xmat(i,j) = xmat(i,j) + sqrt(ws(j))
        endif
        if (itype .eq. 3) then
            xmat(i,j) = xmat(i,j) + sqrt(ws(j))
        endif
 1000 continue
 2000 continue

        return
        end
c
c
c
c
c
        subroutine genmat(xmat,xs,ys,vxs,vys,ws,nxs,xt,yt,wt,
     1      nxt,in)
        implicit real *8 (a-h,o-z)
        dimension xmat(nxt,nxs),xs(1),ys(1),ws(1),vxs(1),vys(1),
     1      xt(1),yt(1),wt(1)

        done = 1
        dtwo = 2
        pi = atan(done)*4
        call prinf('nxs = *',nxs,1)

        do 2000 i=1,nxt


        do 1000 j=1,nxs


        rx = xs(j)-xt(i)
        ry = ys(j)-yt(i)
        top = rx*vxs(j)+ry*vys(j)
        bot = rx**2+ry**2
        xmat(i,j) =top/bot*sqrt(ws(j))*sqrt(wt(i))
        if (in .eq. 1) then
        xmat(i,j) = xmat(i,j)+sqrt(ws(j))*sqrt(wt(i))
        endif
 1000 continue

 2000 continue
       
        return
        end
c
c
c
c
c
        subroutine genmatlog(xmat,xs,ys,vxs,vys,ws,nxs,xt,yt,wt,nxt)
        implicit real *8 (a-h,o-z)
        dimension xmat(2*nxt,nxs),xs(1),ys(1),ws(1),vxs(1),vys(1),
     1      xt(1),yt(1),wt(1)

        done = 1
        pi = atan(done)*4

        do 2000 i=1,nxt
        do 1000 j=1,nxs
            rx = xs(j)-xt(i)
            ry = ys(j)-yt(i)
            top = rx*vxs(j)+ry*vys(j)
            bot = rx**2+ry**2
            xmat(i,j) =top/bot*sqrt(ws(j))*sqrt(wt(i))
            xmat(i,j) = xmat(i,j)+sqrt(ws(j))*sqrt(wt(i))
 1000 continue
 2000 continue


        do 4000 i=1,nxt
        do 3000 j=1,nxs
            rx = xs(j)-xt(i)
            ry = ys(j)-yt(i)
            top = rx*vxs(j)+ry*vys(j)
            bot = rx**2+ry**2
            xmat(nxt+i,j) =log(rx*rx+ry*ry)*sqrt(ws(j))*sqrt(wt(i))
 3000 continue
 4000 continue

        return
        end
c
c
c
c
c

        subroutine genmatlnt(xmat,xs,ys,vxs,vys,ws,nxs,xt,yt,wt,nxt)
        implicit real *8 (a-h,o-z)
        dimension xmat(2*nxt,nxs),xs(1),ys(1),ws(1),vxs(1),vys(1),
     1      xt(1),yt(1),wt(1)

        done = 1
        pi = atan(done)*4

        do 2000 i=1,nxt
        do 1000 j=1,nxs
            rx = xs(j)-xt(i)
            ry = ys(j)-yt(i)
            top = rx*vxs(j)+ry*vys(j)
            bot = rx**2+ry**2
            xmat(i,j) =top/bot*sqrt(ws(j))*sqrt(wt(i))
            xmat(i,j) = xmat(i,j)+sqrt(ws(j))*sqrt(wt(i))
 1000 continue
 2000 continue

        return
        end
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
        subroutine lapcornmat(alpha,akern,nr0)
        implicit real *8 (a-h,o-z)
        dimension rquads(100),wquads(100),x0s(31),w0s(31),
     1      u(4000),v(4000),work(19000),akern(1)
    
        save

        data ifirst/1/

        lww = 4000

        if (ifirst .eq. 1) then

        itype = 1
        call lapdisc(x0s,w0s,u,v,nr0,itype)

        ifirst = 0

        endif

        do 2000 ipt=1,nr0

        call lap_load_quads(alpha,rquads,wquads,nquads,ipt)

        iw = 1
        iw2 = iw+ nr0*nr0*3+20

        call lapkern(x0s(ipt),w0s(ipt),nr0,alpha,rquads,
     1      wquads,nquads,v,
     2      work(iw),work(iw2),akern,ipt)

 2000 continue

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

ccc        call prinf('nfuns=*',nfuns,1)
ccc        call prinf('nroots=*',nroots,1)
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
ccc        call prin2('cond = *',cond,1)
c
        return
        end
c
c
c
c
c
c
        subroutine lapkern(x0,w0,nr0,ang,rq,wq,nq,umat,
     1      xmat,work,amat,ipt)
        implicit real *8 (a-h,o-z)
        dimension coefs(1),rq(1),wq(1),umat(nr0,nr0),
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
    
        top = -x0*sin(pi*ang)
        bottom = (cos(pi*ang)*x0-rq(i))**2+(x0*sin(pi*ang))**2
ccc        bottom = x0**2+rq(i)**2-2*rq(i)*x0*cos(pi*ang)
        work(i) = top/bottom*wq(i)
 1000 continue

        do 2000 i=1,nq

        do 1500 j=1,nr0

        call lapeval(rq(i),j,xmat(i,j))
    
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
        subroutine lapeval(x,ifun,val)
        implicit real *8 (a-h,o-z)
        dimension work(9000),coefs(800 000)
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

        data isvd/31/
        data ifirst/1/

        lww = 4000

        if (ifirst .eq. 1) then

        rewind(isvd)
        call lapsvdretr(isvd,coefs)

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
        subroutine lapinterps(x,rhs,val)
        implicit real *8 (a-h,o-z)
        dimension rquads(100),wquads(100),x0s(31),w0s(31),rhs(1),
     1      u(31,31),v(31,31),work(9000),cfs(31)
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

ccc        call prin2('cfs = *',cfs,nr0)

        do 5500 j=1,nr0

        call lapeval(x,j,vt)
ccc        if (ier .ne. 0) then
ccc            call prinf('Bombing in check interp nestev*',0,0)
ccc            call prinf('ier = *',ier,1)
ccc            stop
ccc        endif

        val = val+vt*cfs(j)

 5500 continue

        return
        end
c
c
c
c
c
        subroutine lapsvdretr(ir,w)

        implicit real *8 (a-h,o-z)
        dimension w(1)       
c
c       This subroutine retrieves the singular function information
c       from a file. It is a copy of quaeretr and has no stand-alone
c       purpose.
 
 
c        read various parameters of the problem (as
c       opposed to the parameters of the memory allocation)
 
 1400 format(i8)
        read(ir,1400) ncs 
c
        call laparrread(ir,w(1),ncs)
c
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

        subroutine lapdisc(roots,weights,u,v,npts,itype)
        implicit real *8 (a-h,o-z)
        dimension roots(1), weights(1), u(1), v(1),rts(31),wts(31),
     1      work(8000)

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

        nr0 = 31
        npts = 31

        data rts/
     1  0.4248728091427118D-09,0.2000828012027937D-07,
     2  0.3042641755217164D-06,0.2562427923973320D-05,
     3  0.1462675091746566D-04,0.6284687286972221D-04,
     4  0.2167091099947979D-03,0.6261749223838039D-03,
     5  0.1564147414280510D-02,0.3458334915994540D-02,
     6  0.6894329756217213D-02,0.1257803460957137D-01,
     7  0.2125847666623854D-01,0.3362551681430298D-01,
     8  0.5020500111634076D-01,0.7127384152413946D-01,
     9  0.9681064242170892D-01,0.1264874592356006D+00,
     *  0.1596988736076681D+00,0.1956182421020692D+00,
     1  0.2332684567715073D+00,0.2715952742253952D+00,
     2  0.3095340179759866D+00,0.3460639806479383D+00,
     3  0.3802481937789684D+00,0.4112588528231759D+00,
     4  0.4383904061106130D+00,0.4610632049720723D+00,
     5  0.4788208628728403D+00,0.4913243862378277D+00,
     *  0.4983473798851192D+00/
c
c
        data wts/
     1  0.2069895123644021D-08,0.6281609490556449D-07,
     2  0.7238037620154733D-06,0.4903893574740540D-05,
     3  0.2320916336521615D-04,0.8417196398696858D-04,
     4  0.2477983425514513D-03,0.6160388924715742D-03,
     5  0.1331108275515546D-02,0.2555493978273085D-02,
     6  0.4435469532342073D-02,0.7058540123084172D-02,
     7  0.1041912469025088D-01,0.1440405037351843D-01,
     8  0.1880190605354557D-01,0.2333233573457130D-01,
     9  0.2768609830036848D-01,0.3156566329638601D-01,
     *  0.3471818297719672D-01,0.3695645833905839D-01,
     1  0.3816713868058617D-01,0.3830813817993708D-01,
     2  0.3739856967650554D-01,0.3550477560720021D-01,
     3  0.3272547225643741D-01,0.2917820068112037D-01,
     4  0.2498834903896012D-01,0.2028128933570136D-01,
     5  0.1517763451200640D-01,0.9791567166907443D-02,
     *  0.4237522244825257D-02/



        do 1000 i=1,npts

        roots(i) = rts(i)
        weights(i) = wts(i)

 1000 continue

        if (itype .gt. 0) then

        lww = 8000
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

ccc        if (ipt .eq. 32) then
ccc        call lap_load_quads32(alpha,roots,weights,nqs)
ccc        endif
ccc
ccc        if (ipt .eq. 33) then
ccc        call lap_load_quads33(alpha,roots,weights,nqs)
ccc        endif

        return
        end
c
c
c
c
c
        subroutine lap_load_quads1(alph,roots,weights,nquads)
        implicit real *8 (a-h,o-z)
        dimension r1(24),r2(26),r3(27),r4(28),r5(28),r6(29),
     1          r7(30),r8(32),r9(35),r10(50),
     1          w1(24),w2(26),w3(27),w4(28),w5(28),w6(29),nqs(10),
     1          w7(30),w8(32),w9(35),w10(50),
     2          roots(1),weights(1)
cc
        data nqs/24,26,27,28,28,29,30,32,35,50/

        alpha = 1-alph


        data r1/
     1  0.2011127076322789D-11,0.1969257302595794D-10,
     2  0.6021113130980540D-10,0.1915367138738721D-09,
     3  0.5676312805168901D-09,0.1847475146866003D-08,
     4  0.8839838436213013D-08,0.7968106073916570D-07,
     5  0.1176819293218899D-05,0.1646187592419303D-04,
     6  0.1536916160650335D-03,0.9278903435385204D-03,
     7  0.3878458365472424D-02,0.1207292527109318D-01,
     8  0.2972017693220962D-01,0.6066310530604551D-01,
     9  0.1065496765457009D+00,0.1658108969229099D+00,
     *  0.2339257194017379D+00,0.3046092966087608D+00,
     1  0.3712176608266038D+00,0.4278351943459719D+00,
     2  0.4698631770397662D+00,0.4942001291153860D+00/
     
        data w1/
     1  0.7100345292757701D-11,0.2684569282626586D-10,
     2  0.6800355950243147D-10,0.2145757853355275D-09,
     3  0.6178969071186513D-09,0.2440834965065668D-08,
     4  0.1642537514812518D-07,0.2007341086906871D-06,
     5  0.3246352819887789D-05,0.4049253004564894D-04,
     6  0.3085162526540789D-03,0.1485710485511582D-02,
     7  0.4932307867566604D-02,0.1218551044865509D-01,
     8  0.2378313188881455D-01,0.3839561851254306D-01,
     9  0.5311048020670475D-01,0.6463673592672839D-01,
     *  0.7052069323144228D-01,0.6972062540517694D-01,
     1  0.6250063385389489D-01,0.4996555893511603D-01,
     2  0.3357402606508456D-01,0.1483649150250057D-01/
     
     
        data r2/
     1  0.4439359627634446D-12,0.5657915334484311D-11,
     2  0.3129023378368669D-10,0.1145463968894215D-09,
     3  0.3307173764469297D-09,0.9118998323542441D-09,
     4  0.2897928251272113D-08,0.1303137443315710D-07,
     5  0.9717879499607587D-07,0.1086347562451390D-05,
     6  0.1239900798972195D-04,0.1064706820721105D-03,
     7  0.6359135123362365D-03,0.2720884038553170D-02,
     8  0.8792181919608292D-02,0.2256239159106124D-01,
     9  0.4797894603917560D-01,0.8754666781154806D-01,
     *  0.1410071514976531D+00,0.2050995766892258D+00,
     1  0.2743692335706042D+00,0.3424585746613621D+00,
     2  0.4032900157023882D+00,0.4518293147918628D+00,
     3  0.4844122131871847D+00,0.4988473249456871D+00/
     
        data w2/
     1  0.1558179529875226D-11,0.1116684328395391D-10,
     2  0.4631375925450503D-10,0.1312755857711208D-09,
     3  0.3328649113410076D-09,0.9600136207246004D-09,
     4  0.3751503763328648D-08,0.2271465225489373D-07,
     5  0.2196479155006577D-06,0.2710346642296979D-05,
     6  0.2878797146754644D-04,0.2095897257787942D-03,
     7  0.1025144566587081D-02,0.3548987363219279D-02,
     8  0.9230868737012548D-02,0.1899820658546648D-01,
     9  0.3227632991272587D-01,0.4681981472991026D-01,
     *  0.5954262011156104D-01,0.6771006429427844D-01,
     1  0.6974921339465328D-01,0.6540383458727202D-01,
     2  0.5541497596955883D-01,0.4105270251271147D-01,
     3  0.2373946918289589D-01,0.5246432410994453D-02/
     
     
        data r3/
     1  0.3275782114508195D-12,0.4360719920972993D-11,
     2  0.2397270587267338D-10,0.8414330963289162D-10,
     3  0.2226457796596375D-09,0.5253949775432701D-09,
     4  0.1287345695944449D-08,0.3841171400056475D-08,
     5  0.1667056452709634D-07,0.1190257791968424D-06,
     6  0.1249444329573695D-05,0.1347966279456316D-04,
     7  0.1117669694304118D-03,0.6548647066141768D-03,
     8  0.2772665976521027D-02,0.8905890345362512D-02,
     9  0.2277215433710190D-01,0.4831592934579409D-01,
     *  0.8803180153265061D-01,0.1416473306254799D+00,
     1  0.2058878729362407D+00,0.2752880585674559D+00,
     2  0.3434841080177071D+00,0.4043962418613735D+00,
     3  0.4529902406812061D+00,0.4855983109522895D+00,
     *  0.4998526964048002D+00/
     
        data w3/
     1  0.1165925116310571D-11,0.8696361953051595D-11,
     2  0.3486796861174212D-10,0.9128218638511749D-10,
     3  0.1980591497113071D-09,0.4478448924762366D-09,
     4  0.1243832884847947D-08,0.4804473560319295D-08,
     5  0.2849373476203588D-07,0.2624142062840769D-06,
     6  0.3036584793471055D-05,0.3068149567393105D-04,
     7  0.2171459309602472D-03,0.1046582751883303D-02,
     8  0.3594959059603920D-02,0.9309645346849151D-02,
     9  0.1911102948507574D-01,0.3241617801267170D-01,
     *  0.4697385767997716D-01,0.5969631356718637D-01,
     1  0.6785080749483618D-01,0.6986849656162480D-01,
     2  0.6549763148916156D-01,0.5548262934432210D-01,
     3  0.4109424346665082D-01,0.2374237687978608D-01,
     *  0.4064087110779493D-02/
     
     
        data r4/
     1  0.3006156316382685D-12,0.5189399795309631D-11,
     2  0.2913379731218590D-10,0.9584712288006311D-10,
     3  0.2355548724705449D-09,0.5156108218073636D-09,
     4  0.1172453617560466D-08,0.3266064676458476D-08,
     5  0.7876817166686988D-08,0.1209879678111230D-07,
     6  0.8095572049284754D-07,0.8108429679364577D-06,
     7  0.9108199489462116D-05,0.8231609258965304D-04,
     8  0.5239421804353486D-03,0.2366369352013057D-02,
     9  0.7968186856469824D-02,0.2107411852418792D-01,
     *  0.4579852991978938D-01,0.8487493170831627D-01,
     1  0.1382188508136373D+00,0.2026163387255946D+00,
     2  0.2725348967899044D+00,0.3414648348589110D+00,
     3  0.4031598335630520D+00,0.4524401011108495D+00,
     4  0.4855320694697002D+00,0.4999999999999998D+00/
     
        data w4/
     1  0.1174008412503454D-11,0.1102394548382832D-10,
     2  0.4098039952730729D-10,0.9705405773305240D-10,
     3  0.1919112401308916D-09,0.4007257559087168D-09,
     4  0.1046364814567519D-08,0.3911372860426915D-08,
     5  -.1763221565983404D-08,0.2146683826183176D-07,
     6  0.1719738297525983D-06,0.1968212802758178D-05,
     7  0.2140822672920536D-04,0.1672848505075958D-03,
     8  0.8762320059433728D-03,0.3199010190216963D-02,
     9  0.8646597085841614D-02,0.1828315009487528D-01,
     *  0.3164740723099194D-01,0.4649635806935555D-01,
     1  0.5964063343536000D-01,0.6820826549972456D-01,
     2  0.7052243401607063D-01,0.6628301018335671D-01,
     3  0.5623842243236505D-01,0.4169335007473397D-01,
     4  0.2410049718234526D-01,0.3973773830726014D-02/
     
     
        data r5/
     1  0.2270117694461802D-12,0.2925335321810481D-11,
     2  0.1193122122738662D-10,0.4600999119483418D-10,
     3  0.1247452765925675D-09,0.2620254479735751D-09,
     4  0.4978506791081685D-09,0.9773726564936648D-09,
     5  0.2286966583499732D-08,0.7588128947094793D-08,
     6  0.4153727187952691D-07,0.3688166050885130D-06,
     7  0.4034249413804439D-05,0.3853878897279502D-04,
     8  0.2692371193869361D-03,0.1343397549765375D-02,
     9  0.4964345797044816D-02,0.1425366221166959D-01,
     *  0.3324622715991303D-01,0.6542623016431652D-01,
     1  0.1120778961959518D+00,0.1714322372868422D+00,
     2  0.2389874444408925D+00,0.3086461403523690D+00,
     3  0.3740255552330144D+00,0.4294626392418753D+00,
     4  0.4705563505200123D+00,0.4943348311679374D+00/
     
        data w5/
     1  0.8908741845410956D-12,0.4606537927053741D-11,
     2  0.1721299923292177D-10,0.5421837708012743D-10,
     3  0.1050110663177549D-09,0.1752089008494301D-09,
     4  0.3168981029992326D-09,0.7192451686503957D-09,
     5  0.2273385535257656D-08,0.1087095186554010D-07,
     6  0.8170478161006344D-07,0.8669295509568767D-06,
     7  0.9567688705175203D-05,0.8145471474146255D-04,
     8  0.4771236087759392D-03,0.1947069425325301D-02,
     9  0.5823117624973524D-02,0.1345458591144119D-01,
     *  0.2514146517747909D-01,0.3944773060295872D-01,
     1  0.5355434891473039D-01,0.6438334289292050D-01,
     2  0.6968781632355567D-01,0.6854869893887751D-01,
     3  0.6125751367567537D-01,0.4888090320092781D-01,
     4  0.3281163519654772D-01,0.1449266293040264D-01/
     
     
        data r6/
     1  0.4620544554476023D-12,0.7762405944024359D-11,
     2  0.3758185173562671D-10,0.9912372543208169D-10,
     3  0.1891220570699504D-09,0.3075976000804919D-09,
     4  0.4721052761123797D-09,0.7379610908194151D-09,
     5  0.1277404711060244D-08,0.2765265220644882D-08,
     6  0.8854617026326054D-08,0.4814009936280325D-07,
     7  0.4212903755067134D-06,0.4432541399543053D-05,
     8  0.4070630268382546D-04,0.2772365588457514D-03,
     9  0.1364098307239115D-02,0.5003685574428157D-02,
     *  0.1431076500514463D-01,0.3331123504684821D-01,
     1  0.6548469609298126D-01,0.1121181651383883D+00,
     2  0.1714503064762992D+00,0.2389869862427821D+00,
     3  0.3086351071781445D+00,0.3740121146053990D+00,
     4  0.4294523758995387D+00,0.4705512674415662D+00,
     *  0.4943337677297003D+00/
     
        data w6/
     1  0.1815640705394179D-11,0.1579537749063268D-10,
     2  0.4551506121740744D-10,0.7665338478967442D-10,
     3  0.1031448267482512D-09,0.1364895314454504D-09,
     4  0.2007723679474280D-09,0.3553507947357899D-09,
     5  0.8121693550026701D-09,0.2593729891360207D-08,
     6  0.1253679543260258D-07,0.9450487323254477D-07,
     7  0.9782437134781231D-06,0.1031825008809131D-04,
     8  0.8470100341458758D-04,0.4860870798476958D-03,
     9  0.1963403042171314D-02,0.5842855029108844D-02,
     *  0.1346874857325770D-01,0.2514229359985423D-01,
     1  0.3943433569898342D-01,0.5353271594667627D-01,
     2  0.6436191537236955D-01,0.6967291089246434D-01,
     3  0.6854246535604409D-01,0.6125846969462617D-01,
     4  0.4888568078879299D-01,0.3281667790278071D-01,
     *  0.1449533214270162D-01/
     
     
        data r7/
     1  0.3377183495056858D-12,0.5061173706278048D-11,
     2  0.2747842804573905D-10,0.8301470122638603D-10,
     3  0.1666894989600420D-09,0.2674385562540259D-09,
     4  0.3889828900021385D-09,0.5560023859886725D-09,
     5  0.8391666619559095D-09,0.1476476630611478D-08,
     6  0.3547281632893386D-08,0.1398437218249575D-07,
     7  0.9698048123473151D-07,0.1010833209621538D-05,
     8  0.1112011675348822D-04,0.9576742660432768D-04,
     9  0.5832670889884655D-03,0.2549436143915195D-02,
     *  0.8389187810379170D-02,0.2183506332793756D-01,
     1  0.4692625198220977D-01,0.8629080009835120D-01,
     2  0.1397612477265724D+00,0.2040966200345636D+00,
     3  0.2737936365248348D+00,0.3424065561856170D+00,
     4  0.4037621286111117D+00,0.4527443964495845D+00,
     5  0.4856261036715324D+00,0.4999999999999982D+00/
     
        data w7/
     1  0.1251131687425391D-11,0.1040402711773188D-10,
     2  0.3770508900987969D-10,0.7206745920106545D-10,
     3  0.9307234578302725D-10,0.1090356224889304D-09,
     4  0.1378581503624766D-09,0.2065327365919365D-09,
     5  0.3924259573122326D-09,0.1017258345845801D-08,
     6  0.3893783933615254D-08,0.2320291506307651D-07,
     7  0.2121974557313094D-06,0.2461797821101652D-05,
     8  0.2566755786763914D-04,0.1899406095299511D-03,
     9  0.9530174645131591D-03,0.3376584208051607D-02,
     *  0.8943163645067834D-02,0.1865331144527983D-01,
     1  0.3199223155582472D-01,0.4671340677242683D-01,
     2  0.5967216948686715D-01,0.6805777359805746D-01,
     3  0.7024095523903175D-01,0.6594294134341287D-01,
     4  0.5591040794697538D-01,0.4143308696956146D-01,
     5  0.2394485188766625D-01,0.3947787100279431D-02/
     
     
        data r8/
     1  0.6732202309534888D-12,0.1185225757550403D-10,
     2  0.5662123482562398D-10,0.1346954267758867D-09,
     3  0.2217767823427694D-09,0.3082332257571477D-09,
     4  0.3996366734541743D-09,0.5055229052446436D-09,
     5  0.5129525389156416D-09,0.6814005380679641D-09,
     6  0.1014113554540899D-08,0.1495925348236192D-08,
     7  0.1909746018430992D-08,0.5899838515042498D-08,
     8  0.3470077605428502D-07,0.3475496131949415D-06,
     9  0.4109415012978937D-05,0.4024236449608334D-04,
     *  0.2808420599248469D-03,0.1389714301105937D-02,
     1  0.5091393406554710D-02,0.1451503292807063D-01,
     2  0.3367398290989258D-01,0.6600734749726041D-01,
     3  0.1127538236469479D+00,0.1721198157229152D+00,
     4  0.2396063590431542D+00,0.3091394987010144D+00,
     5  0.3743685914298284D+00,0.4296614212559461D+00,
     6  0.4706410106080489D+00,0.4943512825950656D+00/
     
        data w8/
     1  0.2684801217244245D-11,0.2444014481683718D-10,
     2  0.6511001470040474D-10,0.8600039260946139D-10,
     3  0.8674877910019337D-10,0.8718703858157714D-10,
     4  0.9812056431209263D-10,0.1720235124938439D-10,
     5  0.1148113325090970D-09,0.2217399552494816D-09,
     6  0.5008505253431329D-09,-.4964875999510656D-10,
     7  0.1648705495305984D-08,0.8494119790327150D-08,
     8  0.7231276562058997D-07,0.8543305468191958D-06,
     9  0.9941175254640458D-05,0.8535564251681466D-04,
     *  0.4958494235589468D-03,0.2001655079069901D-02,
     1  0.5931570336154223D-02,0.1361108616861064D-01,
     2  0.2530977775357408D-01,0.3957816297017446D-01,
     3  0.5360933456915790D-01,0.6435227160870305D-01,
     4  0.6958567872717245D-01,0.6840525071344239D-01,
     5  0.6110554573637140D-01,0.4874839354192948D-01,
     6  0.3271851811211082D-01,0.1445067040081394D-01/
     
     
        data r9/
     1  0.9560699576794486D-12,0.1724267122774561D-10,
     2  0.8073283067348429D-10,0.1730026116168509D-09,
     3  0.2521151091269899D-09,0.3136668469394605D-09,
     4  0.3636976578172607D-09,0.3649878012161007D-09,
     5  0.4126394342429477D-09,0.4643636555268132D-09,
     6  0.5307673478643381D-09,0.6340844537987071D-09,
     7  0.7826578755373541D-09,0.8514645608652128D-09,
     8  0.8763980176548085D-09,0.1417713157547040D-08,
     9  0.4148019445337168D-08,0.2600374911238970D-07,
     *  0.2872755124968633D-06,0.3676422115126163D-05,
     1  0.3780384896203315D-04,0.2709773802787324D-03,
     2  0.1360558234004167D-02,0.5025656729606589D-02,
     3  0.1439698577063150D-01,0.3349878253150892D-01,
     4  0.6578589451232781D-01,0.1125098372261191D+00,
     5  0.1718816688450068D+00,0.2393986689216064D+00,
     6  0.3089778812193813D+00,0.3742582205058862D+00,
     7  0.4295982834478152D+00,0.4706143486062518D+00,
     *  0.4943461252766730D+00/
     
        data w9/
     1  0.3834514748066145D-11,0.3584938165020744D-10,
     2  0.8730202871566513D-10,0.8898910032559879D-10,
     3  0.6928775438392730D-10,0.5512673077645739D-10,
     4  0.5392679393011748D-11,0.4297912092032414D-10,
     5  0.4839265050848135D-10,0.5673811441315703D-10,
     6  0.7940883558838117D-10,0.1361986546588512D-09,
     7  0.1329546585935449D-10,0.4315937691367930D-09,
     8  -.1393081301511041D-09,0.1048615983572650D-08,
     9  0.6029491813441373D-08,0.5665880962939382D-07,
     *  0.7328999441886174D-06,0.9127319793209112D-05,
     1  0.8153923737734136D-04,0.4836748014070643D-03,
     2  0.1974225521621125D-02,0.5886049785120631D-02,
     3  0.1355388568942543D-01,0.2525544216171276D-01,
     4  0.3954222295395332D-01,0.5360098212239736D-01,
     5  0.6437161693505256D-01,0.6962567908384806D-01,
     6  0.6845560972870022D-01,0.6115612998974981D-01,
     7  0.4879115176290363D-01,0.3274801664029549D-01,
     *  0.1446384861469970D-01/
     
     
        data r10/
     1  0.2664871539775631D-11,0.4494617887572072D-10,
     2  0.1583243322555276D-09,0.2560402957110188D-09,
     3  0.3145636657177468D-09,0.3504603590199162D-09,
     4  0.3741479658807478D-09,0.3905666991619013D-09,
     5  0.3959433719428723D-09,0.4023159141196272D-09,
     6  0.4109908278457514D-09,0.4177388895130347D-09,
     7  0.4234901250749251D-09,0.4291022046688167D-09,
     8  0.4354584025863766D-09,0.4435855030429317D-09,
     9  0.4445590325378188D-09,0.4445674064664874D-09,
     *  0.4548049085016083D-09,0.4712343653932457D-09,
     1  0.4966324315978329D-09,0.5056918181773838D-09,
     2  0.5319682607791664D-09,0.5372094697102697D-09,
     3  0.5795878723055248D-09,0.6030668191708654D-09,
     4  0.6175101032014307D-09,0.8441605856660029D-09,
     5  0.1766783765524723D-08,0.7930025850753167D-08,
     6  0.6148233478144169D-07,0.6756073606895268D-06,
     7  0.8018463874219462D-05,0.7278305369048915D-04,
     8  0.4520654383684201D-03,0.1952192234423521D-02,
     9  0.5917907788510684D-02,0.1147847224314214D-01,
     *  0.2027887158781514D-01,0.3687113870235724D-01,
     1  0.5799915072492282D-01,0.9085897921896920D-01,
     2  0.1394414512830966D+00,0.1995861070760595D+00,
     3  0.2656162494176793D+00,0.3313727002674152D+00,
     4  0.3910910126583823D+00,0.4401426792791866D+00,
     5  0.4754177332512724D+00,0.4953214428235242D+00/
     
        data w10/
     1  0.1067635018510220D-10,0.8554430427375151D-10,
     2  0.1186164698418405D-09,0.7548995213276111D-10,
     3  0.4478845211821058D-10,0.2865572004321317D-10,
     4  0.1949221812849575D-10,0.1375639980206428D-10,
     5  0.5849280780502558D-14,0.9991692056375765D-11,
     6  0.7544495001740319D-11,0.6104402742514567D-11,
     7  0.5539487838845665D-11,0.5830296137946821D-11,
     8  0.7048514085282304D-11,0.9619001393821893D-11,
     9  -.1498424355208045D-10,0.1477157596293244D-10,
     *  0.1338533805207394D-10,0.2006921567906984D-10,
     1  0.3172970869268323D-10,0.6774748496044721D-12,
     2  -.1797185547264806D-10,0.7216768129792340D-10,
     3  0.1777380765987373D-10,-.8385564547677633D-10,
     4  0.1915188128937285D-09,0.3822762218788367D-09,
     5  0.1948372829807346D-08,0.1456272249496719D-07,
     6  0.1376911277698696D-06,0.1692856588078881D-05,
     7  0.1906757626340283D-04,0.1469070028562622D-03,
     8  0.7408568319250556D-03,0.2520211724846034D-02,
     9  0.5347571221662365D-02,0.5719554963280221D-02,
     *  0.1308486685921028D-01,0.1897252784686465D-01,
     1  0.2508830495933117D-01,0.4115494417661554D-01,
     2  0.5523418600097160D-01,0.6410060501006860D-01,
     3  0.6691685703912846D-01,0.6362257297604202D-01,
     4  0.5503628352325135D-01,0.4254867710056860D-01,
     5  0.2774241651285831D-01,0.1200174053918320D-01/
c       
c
        if (alpha .lt. 0.1d0) then

        do 1000 i=1,nqs(1)
        roots(i) =r1(i)
        weights(i) = w1(i)
 1000 continue

        nquads = nqs(1)

        return

        end if


        if (alpha .lt. 0.2d0) then

        do 2000 i=1,nqs(2)
        roots(i) =r2(i)
        weights(i) = w2(i)
 2000 continue

        nquads = nqs(2)

        return

        end if


        if (alpha .lt. 0.3d0) then

        do 3000 i=1,nqs(3)
        roots(i) =r3(i)
        weights(i) = w3(i)
 3000 continue

        nquads = nqs(3)

        return

        end if


        if (alpha .lt. 0.4d0) then

        do 4000 i=1,nqs(4)
        roots(i) =r4(i)
        weights(i) = w4(i)
 4000 continue

        nquads = nqs(4)
        return

        end if


        if (alpha .lt. 0.5d0) then

        do 5000 i=1,nqs(5)
        roots(i) =r5(i)
        weights(i) = w5(i)
 5000 continue

        nquads = nqs(5)

        return

        end if


        if (alpha .lt. 0.6d0) then

        do 6000 i=1,nqs(6)
        roots(i) =r6(i)
        weights(i) = w6(i)
 6000 continue

        nquads = nqs(6)
        return

        end if


        if (alpha .lt. 0.7d0) then

        do 7000 i=1,nqs(7)
        roots(i) =r7(i)
        weights(i) = w7(i)
 7000 continue

        nquads = nqs(7)

        return

        end if


        if (alpha .lt. 0.8d0) then

        do 8000 i=1,nqs(8)
        roots(i) =r8(i)
        weights(i) = w8(i)
 8000 continue

        nquads = nqs(8)

        end if

        if (alpha .lt. 0.9d0) then

        do 9000 i=1,nqs(9)
        roots(i) =r9(i)
        weights(i) = w9(i)
 9000 continue

        nquads = nqs(9)

        return

        end if



        do 9500 i=1,nqs(10)
        roots(i) =r10(i)
        weights(i) = w10(i)
 9500 continue

        nquads = nqs(10)


        return
        end
c
c
c
c
c
        subroutine lap_load_quads2(alpha,roots,weights,nquads)
        implicit real *8 (a-h,o-z)
        dimension r1(53),
     1          w1(53),nqs(1),
     2          roots(1),weights(1)
cc
        data nqs/53/

        data r1/
     1  0.1721064702155506D-07,0.4221874500580541D-07,
     2  0.4219221511627947D-07,0.1721720421039873D-07,
     3  0.1721199354201211D-07,0.2313609716195958D-01,
     4  0.6892187583661592D-11,0.4659553311399908D-01,
     5  0.9873895410153791D-02,0.8250359085601553D-01,
     6  0.3541389140984136D-02,0.6205481228200723D-10,
     7  0.2896847952313716D-08,0.1310667594278255D+00,
     8  0.1564032637111113D-12,0.1050402128721857D-02,
     9  0.1900766100069715D+00,0.2586296617830319D-03,
     *  0.2552977728327131D+00,0.5537390463642538D-04,
     1  0.3213521275592022D+00,0.1148878171225803D-04,
     2  0.2785686871836292D-09,0.3827018732216904D+00,
     3  0.2676254224794929D-05,0.4344229152051298D+00,
     4  0.8516641780641632D-09,0.7840499528762684D-06,
     5  0.4726495862102693D+00,0.2986990802885736D-06,
     6  0.1440943127509824D-06,0.8395042818529180D-07,
     7  0.3715851700308247D-08,0.5905924598310086D-08,
     8  0.5644253962006173D-07,0.8295614475384170D-08,
     9  0.4947395491971186D+00,0.1980572155581314D-08,
     *  0.1063485256570590D-07,0.3413610513269148D-07,
     1  0.1274939202983254D-07,0.2923691773922947D-07,
     2  0.1455394198565656D-07,0.2610834653029387D-07,
     3  0.1603417862529542D-07,0.2402802585204203D-07,
     4  0.2259385730614130D-07,0.1816364083359098D-07,
     5  0.2156626891337064D-07,0.1892386028702783D-07,
     6  0.2079140283674191D-07,0.1956381559522652D-07,
     *  0.2015672450706506D-07/
     
        data w1/
     1  0.2907664128406606D-08,0.6481778188817093D-08,
     2  0.3960973546732078D-08,0.2405638157085848D-08,
     3  -.4258131867461419D-08,0.1785923490677536D-01,
     4  0.1893273752048088D-10,0.2945330926900063D-01,
     5  0.9231322857916125D-02,0.4240359905306063D-01,
     6  0.3961844790140429D-02,0.1098801083695767D-09,
     7  0.2839210388001240D-11,0.5434747172083565D-01,
     8  0.8416965719084186D-12,0.1377847409752385D-02,
     9  0.6295055328616845D-01,0.3838129149242489D-03,
     *  0.6657826050310981D-01,0.8750139094872912D-04,
     1  0.6459432554802949D-01,0.1768426990659573D-04,
     2  0.3568471981131066D-09,0.5727721659449515D-01,
     3  0.3625890065908139D-05,0.4551823337718371D-01,
     4  0.8258789506703638D-09,0.8589360642781039D-06,
     5  0.3049445825287123D-01,0.2508271868760623D-06,
     6  0.9026836237622779D-07,0.3875990262545341D-07,
     7  0.1998903643986343D-08,0.2335562760401488D-08,
     8  0.1910220310843507D-07,0.2401036730510756D-08,
     9  0.1345839215696265D-01,0.1440504773138114D-08,
     *  0.2248244640621604D-08,0.6164989270461275D-08,
     1  0.1966975855875351D-08,0.3854110421466081D-08,
     2  0.1640207909565596D-08,0.2518087284621762D-08,
     3  0.1325662247441561D-08,0.1707185106778226D-08,
     4  0.1199903312685546D-08,0.8417273273444189D-09,
     5  0.8802152032091548D-09,0.6890781670559930D-09,
     6  0.6877834295710819D-09,0.6028865529933481D-09,
     *  0.5979268620744305D-09/
c        
        do 1000 i=1,nqs(1)
        roots(i) =r1(i)
        weights(i) = w1(i)
 1000 continue

        nquads = nqs(1)

        return
        end
c
c
c
c
c
        subroutine lap_load_quads3(alpha,roots,weights,nquads)
        implicit real *8 (a-h,o-z)
        dimension r1(54),
     1          w1(54),nqs(1),
     2          roots(1),weights(1)
cc
        data nqs/54/

        data r1/
     1  0.3031981229141293D-06,0.3034843014449907D-06,
     2  0.4035447946245645D-06,0.1136109394808162D+00,
     3  0.8170637500341671D-01,0.1505762063177619D+00,
     4  0.5618783116861293D-01,0.3381545074967220D-01,
     5  0.8466260638094636D-08,0.2002648980270258D+00,
     6  0.3100706975835978D-06,0.3945936016747445D-06,
     7  0.2610472053558233D+00,0.2837535437519533D-08,
     8  0.1682527012742313D-01,0.3650567376184804D-11,
     9  0.2944411401020729D-06,0.1824958433841856D-04,
     *  0.5821282831784765D-04,0.3248463379851495D+00,
     1  0.9850337410759095D-10,0.7005596356525398D-02,
     2  0.3191115212409660D-06,0.3848224663139271D+00,
     3  0.2040514557933385D-07,0.6862664254292944D-05,
     4  0.7127034618791439D-09,0.2854314958696790D-06,
     5  0.2078941685359014D-03,0.3146846235425630D-05,
     6  0.2463145708244307D-02,0.4355843950157197D+00,
     7  0.4054213980031248D-07,0.3302723919298636D-06,
     8  0.1726187582185624D-05,0.6838985426059600D-07,
     9  0.1984503069202102D-06,0.1096521589110035D-05,
     *  0.5026897251605327D-06,0.1688516331320449D-06,
     1  0.2746248268538406D-06,0.1013008254905755D-06,
     2  0.6064222738456275D-06,0.2236834226257964D-06,
     3  0.7808190040665221D-06,0.1358022161888669D-06,
     4  0.4731358140496329D+00,0.4374386865703624D-06,
     5  0.3450540632711618D-06,0.2444657017414606D-06,
     6  0.7480905663761844D-03,0.3654701691784059D-06,
     7  0.2612221021704417D-06,0.4948337883804901D+00/
     
        data w1/
     1  0.3327633440548850D-07,-.2598473499122736D-07,
     2  -.4284293412004128D-12,0.3399263400557547D-01,
     3  0.2872712720662715D-01,0.4213893672107261D-01,
     4  0.2353274232280126D-01,0.2043693029761705D-01,
     5  0.8322360954501960D-08,0.5658516858860414D-01,
     6  0.8492459526343308D-08,0.3492201127342938D-07,
     7  0.6356148258523323D-01,0.3438534226270521D-08,
     8  0.1327103591916060D-01,0.1625098557648337D-10,
     9  0.8405054211042835D-08,0.1963504477594415D-04,
     *  0.7173304049450671D-04,0.6290254617235863D-01,
     1  0.2427059484873190D-09,0.6743120024890920D-02,
     2  0.9868114969099135D-08,0.5615108263510127D-01,
     3  0.1589818682259980D-07,0.6024116289551035D-05,
     4  0.1141676116521867D-08,0.9769270737034089D-08,
     5  0.2692464013541854D-03,0.2156121168141820D-05,
     6  0.2770054801558177D-02,0.4470744848868586D-01,
     7  0.2428884206223941D-07,0.1269439392324950D-07,
     8  0.9000934783726207D-06,0.3093965830889290D-07,
     9  0.2748843440119910D-07,0.4293976231492270D-06,
     *  0.8089276217543890D-07,0.3156124676546850D-07,
     1  0.1197562870895306D-07,0.3428131959839436D-07,
     2  0.1316148612006038D-06,0.2296917562736268D-07,
     3  0.2280672211329109D-06,0.3420454598194169D-07,
     4  0.2995626043418303D-01,0.5216004664575758D-07,
     5  0.1719627150851891D-07,0.1867067806531857D-07,
     6  0.9345703932673113D-03,0.2413440827734136D-07,
     7  0.1495687930516568D-07,0.1321782928390900D-01/
        
        do 1000 i=1,nqs(1)
        roots(i) =r1(i)
        weights(i) = w1(i)
 1000 continue

        nquads = nqs(1)

        return
        end
c
c
c
c
c
        subroutine lap_load_quads4(alpha,roots,weights,nquads)
        implicit real *8 (a-h,o-z)
        dimension r1(50),
     1          w1(50),nqs(1),
     2          roots(1),weights(1)
cc
        data nqs/50/

        data r1/
     1  0.4075186792304078D-10,0.1004290431258064D-08,
     2  0.7134609474575359D-08,0.2694611960626939D-07,
     3  0.7238037966216742D-07,0.1592904293572442D-06,
     4  0.3011790135240113D-06,0.5000729774658943D-06,
     5  0.7440956150808451D-06,0.1012127449115350D-05,
     6  0.1281440432010833D-05,0.1533753687284546D-05,
     7  0.1757804259742940D-05,0.1949052601661265D-05,
     8  0.2107975506171448D-05,0.2238111363641073D-05,
     9  0.2344483197482831D-05,0.2432685359333486D-05,
     *  0.2508835685873898D-05,0.2580345685659945D-05,
     1  0.2656329625074474D-05,0.2746730839177187D-05,
     2  0.2862357440390332D-05,0.3016997344295927D-05,
     3  0.3230706639018607D-05,0.3534702060806685D-05,
     4  0.3980051745286484D-05,0.4654791235686868D-05,
     5  0.5719952497255900D-05,0.7490512956284247D-05,
     6  0.1063091020493061D-04,0.1666694127298718D-04,
     7  0.2944824478308081D-04,0.5968915305405296D-04,
     8  0.1398221868334364D-03,0.3717367381018914D-03,
     9  0.1062357842764899D-02,0.3022090143829841D-02,
     *  0.7993888242089610D-02,0.1888698019588791D-01,
     1  0.3932069864969891D-01,0.7237058139285205D-01,
     2  0.1191148109141200D+00,0.1778329146915831D+00,
     3  0.2442477909548539D+00,0.3125431772243340D+00,
     4  0.3765839895192876D+00,0.4308830048468084D+00,
     5  0.4711438971224530D+00,0.4944471905729478D+00/
     
        data w1/
     1  0.1748353124349681D-09,0.2437837642520420D-08,
     2  0.1128456617758954D-07,0.3027320615459072D-07,
     3  0.6336534144827145D-07,0.1128842267302356D-06,
     4  0.1712690202257178D-06,0.2244506929413700D-06,
     5  0.2599349427722403D-06,0.2722355752192134D-06,
     6  0.2633030422770262D-06,0.2394551171539172D-06,
     7  0.2079207217044559D-06,0.1746893899365566D-06,
     8  0.1437729888501023D-06,0.1173515897110639D-06,
     9  0.9632567140657296D-07,0.8108266719644795D-07,
     *  0.7244106835425470D-07,0.7212158706357176D-07,
     1  0.8150565971238234D-07,0.1010585341518666D-06,
     2  0.1324393106740980D-06,0.1801066397066189D-06,
     3  0.2523797342264568D-06,0.3638590148294289D-06,
     4  0.5409654697446045D-06,0.8341416382452460D-06,
     5  0.1345809257617857D-05,0.2298713268066356D-05,
     6  0.4215945202241698D-05,0.8435634391310479D-05,
     7  0.1871047760860164D-04,0.4650823909874288D-04,
     8  0.1287062876769592D-03,0.3813368032724351D-03,
     9  0.1127020478094851D-02,0.3077575460076296D-02,
     *  0.7358051888327522D-02,0.1505391107821526D-01,
     1  0.2635272955538953D-01,0.3994109274313457D-01,
     2  0.5325576054184997D-01,0.6344908474862458D-01,
     3  0.6838980352354944D-01,0.6716197363629001D-01,
     4  0.5999627059028499D-01,0.4788310777178002D-01,
     5  0.3215217086373277D-01,0.1420478598078465D-01/
c        

ccc        call prin2('alpha = *',alpha,1)
        do 1000 i=1,nqs(1)
        roots(i) =r1(i)
        weights(i) = w1(i)
 1000 continue

        nquads = nqs(1)

        return
        end
c
c
c
c
c
        subroutine lap_load_quads5(alpha,roots,weights,nquads)
        implicit real *8 (a-h,o-z)
        dimension r1(51),
     1          w1(51),nqs(1),
     2          roots(1),weights(1)
cc
        data nqs/51/

        data r1/
     1  0.5036332281036991D-01,0.8106811486361801D-01,
     2  0.2883860385806556D-01,0.1243852946227257D+00,
     3  0.1456303033379679D-01,0.6582227858760967D-06,
     4  0.2249069204612260D-03,0.4724419459747517D-03,
     5  0.1803384336693304D+00,0.1313851107779739D-05,
     6  0.2779904600016858D-06,0.7478485051763904D-10,
     7  0.1213272222526594D-03,0.1100665857493214D-02,
     8  0.2449722701278668D+00,0.6523267945200123D-02,
     9  0.2291311970226934D-05,0.9017482930757257D-07,
     *  0.7393089537576561D-04,0.3557987813430031D-05,
     1  0.5010324879193380D-04,0.3123529368164803D+00,
     2  0.5013524183591631D-05,0.3706306234617949D-04,
     3  0.2146124196743610D-08,0.6531168950997837D-05,
     4  0.2940212475666782D-04,0.1920474710319787D-07,
     5  0.2701383861037346D-02,0.7996821620127900D-05,
     6  0.3761047347045070D+00,0.2463722677522962D-04,
     7  0.9330727645938260D-05,0.2153630386593168D-04,
     8  0.1049158477553378D-04,0.1944360833888413D-04,
     9  0.1146992183625065D-04,0.1798759214651191D-04,
     *  0.4304747128980796D+00,0.1227782398596902D-04,
     1  0.1694570953323115D-04,0.1293938031800491D-04,
     2  0.1617793548951332D-04,0.4709357126859764D+00,
     3  0.1348365631380701D-04,0.1559204194214563D-04,
     4  0.1394077696780612D-04,0.4944032892103911D+00,
     5  0.1512358774909918D-04,0.1434178470245712D-04,
     *  0.1472208188214368D-04/
     
        data w1/
     1  0.2553177336265191D-01,0.3659085868260749D-01,
     2  0.1779220778253488D-01,0.5003242231843249D-01,
     3  0.1089156762764697D-01,0.5049757375509308D-06,
     4  0.1531243865429688D-03,0.3780313222469843D-03,
     5  0.6117773440237750D-01,0.8148807178825559D-06,
     6  0.2702539986771738D-06,0.3277736841139696D-09,
     7  0.6725442672349538D-04,0.9713747075329820D-03,
     8  0.6707602716295460D-01,0.5553178639270571D-02,
     9  0.1134313131497263D-05,0.1177786985164449D-06,
     *  0.3246576561110129D-04,0.1381543634674728D-05,
     1  0.1714238212919775D-04,0.6660476965084069D-01,
     2  0.1507718354342124D-05,0.9775524520182279D-05,
     3  0.5649107986912194D-08,0.1508195945868811D-05,
     4  0.5935799958805714D-05,0.3498865257739492D-07,
     5  0.2434834619849551D-02,0.1409767816703560D-05,
     6  0.5993013996999632D-01,0.3789316879611132D-05,
     7  0.1251288744744989D-05,0.2517408630069627D-05,
     8  0.1069028451479205D-05,0.1727835848596710D-05,
     9  0.8898227574271320D-06,0.1220257044557193D-05,
     *  0.4805219589004411D-01,0.7301142709962407D-06,
     1  0.8863034684607012D-06,0.5979095407391597D-06,
     2  0.6643025009740355D-06,0.3235569394011630D-01,
     3  0.4956676863040582D-06,0.5179994430875070D-06,
     4  0.4236776949459906D-06,0.1431482549549890D-01,
     5  0.4271532197142622D-06,0.3841822166793690D-06,
     *  0.3834779436402566D-06/
c        

        do 1000 i=1,nqs(1)
        roots(i) =r1(i)
        weights(i) = w1(i)
 1000 continue

        nquads = nqs(1)
        return
        end
c
c
c
c
c
        subroutine lap_load_quads6(alpha,roots,weights,nquads)
        implicit real *8 (a-h,o-z)
        dimension r1(50),
     1          w1(50),nqs(1),
     2          roots(1),weights(1)
cc
        data nqs/50/

        data r1/
     1  0.8913680384897730D-10,0.3758677627794968D-08,
     2  0.4246343026150192D-07,0.2323934606348723D-06,
     3  0.7723635848790863D-06,0.1665278022140740D-05,
     4  0.2891975741814302D-05,0.5316047287308226D-05,
     5  0.9186411897376972D-05,0.1434924068989690D-04,
     6  0.2041267218662965D-04,0.2685123900609724D-04,
     7  0.3316783863206717D-04,0.3899757633573759D-04,
     8  0.4413567808078195D-04,0.4851569511680622D-04,
     9  0.5216836297108666D-04,0.5518100192875049D-04,
     *  0.5766695836831929D-04,0.5974949406686885D-04,
     1  0.6156389949128864D-04,0.6327596773333924D-04,
     2  0.6509024312825305D-04,0.6722869980208725D-04,
     3  0.6992865950587482D-04,0.7348614962984055D-04,
     4  0.7832351916849669D-04,0.8508626661152368D-04,
     5  0.9480861722767015D-04,0.1092291927523278D-03,
     6  0.1314327880519646D-03,0.1672313938492768D-03,
     7  0.2283108316462041D-03,0.3397936584804771D-03,
     8  0.5591404275903943D-03,0.1023710331534166D-02,
     9  0.2064379911897874D-02,0.4440083182076255D-02,
     *  0.9697780008350328D-02,0.2045349685623962D-01,
     1  0.4014606984228578D-01,0.7197634510457650D-01,
     2  0.1174189976848814D+00,0.1751883302066931D+00,
     3  0.2412569625127708D+00,0.3098034502149253D+00,
     4  0.3745021117412749D+00,0.4296060035570277D+00,
     5  0.4705808436657783D+00,0.4943358130149228D+00/
     
        data w1/
     1  0.4350658926861274D-09,0.1099330607821909D-07,
     2  0.8526115622795454D-07,0.3317311764512352D-06,
     3  0.7598808921048325D-06,0.9631068418010788D-06,
     4  0.1709642730061723D-05,0.3153949056643073D-05,
     5  0.4562814122746929D-05,0.5695270883939543D-05,
     6  0.6341042445370876D-05,0.6451495291562548D-05,
     7  0.6120500947023079D-05,0.5505285871057753D-05,
     8  0.4760575118191132D-05,0.4005278355022624D-05,
     9  0.3315128360517515D-05,0.2729318552628406D-05,
     *  0.2263109830584684D-05,0.1924137460458223D-05,
     1  0.1732055079947308D-05,0.1726726093628383D-05,
     2  0.1938918466619592D-05,0.2376727647798562D-05,
     3  0.3071424010295989D-05,0.4112425857812956D-05,
     4  0.5667311693381996D-05,0.8025664071339248D-05,
     5  0.1169896982959578D-04,0.1763540363691463D-04,
     6  0.2769341342216337D-04,0.4573959693751407D-04,
     7  0.8032810471359067D-04,0.1515076890514650D-03,
     8  0.3083368138821959D-03,0.6721276875968784D-03,
     9  0.1529268736536010D-02,0.3474029766509099D-02,
     *  0.7480493559125862D-02,0.1462276198853265D-01,
     1  0.2532316720389544D-01,0.3860202704214119D-01,
     2  0.5206823478833899D-01,0.6278177160387035D-01,
     3  0.6836234130715143D-01,0.6765775292030058D-01,
     4  0.6077134970891139D-01,0.4867741751645816D-01,
     5  0.3275705406760659D-01,0.1448792190116674D-01/
c        
c        

        do 1000 i=1,nqs(1)
        roots(i) =r1(i)
        weights(i) = w1(i)
 1000 continue

        nquads = nqs(1)
        return
        end
c
c
c
c
c
        subroutine lap_load_quads7(alpha,roots,weights,nquads)
        implicit real *8 (a-h,o-z)
        dimension r1(50),
     1          w1(50),nqs(1),
     2          roots(1),weights(1)
cc
        data nqs/50/

        data r1/
     1  0.9800427825237251D-13,0.2562786643549564D-08,
     2  0.5743279011991618D-07,0.4426791591522213D-06,
     3  0.1849380618729419D-05,0.5058743708441882D-05,
     4  0.1015461463994003D-04,0.1759271911514539D-04,
     5  0.2916565432853425D-04,0.4521376311264107D-04,
     6  0.6470301662530556D-04,0.8601987306609323D-04,
     7  0.1074830009421787D-03,0.1277375081503489D-03,
     8  0.1459266855923845D-03,0.1616732281317443D-03,
     9  0.1749644614665763D-03,0.1860183844441165D-03,
     *  0.1951719487528045D-03,0.2028096963271594D-03,
     1  0.2093445393669845D-03,0.2152588087900483D-03,
     2  0.2211635711087720D-03,0.2277678369096585D-03,
     3  0.2358100851735140D-03,0.2461245029586920D-03,
     4  0.2598195558182565D-03,0.2785193861819734D-03,
     5  0.3047299301059237D-03,0.3424847033727093D-03,
     6  0.3985751684327210D-03,0.4850138761704609D-03,
     7  0.6241830364995569D-03,0.8599986137126943D-03,
     8  0.1282632013437025D-02,0.2082694574270723D-02,
     9  0.3662888419143306D-02,0.6838642413948287D-02,
     *  0.1311660767963404D-01,0.2490334435723039D-01,
     1  0.4531484067437452D-01,0.7730435299300017D-01,
     2  0.1223220080758121D+00,0.1792579216716234D+00,
     3  0.2443289626889013D+00,0.3119136158696771D+00,
     4  0.3758021665132001D+00,0.4302927798377139D+00,
     5  0.4708552254389553D+00,0.4943872752796954D+00/
     
        data w1/
     1  0.3661360085842074D-10,0.1024240757556196D-07,
     2  0.1411176714862381D-06,0.7492263121455093D-06,
     3  0.2213003928727704D-05,0.4209395041465165D-05,
     4  0.6001732351405840D-05,0.9270483358604790D-05,
     5  0.1390658406837058D-04,0.1800698072730429D-04,
     6  0.2069388970745025D-04,0.2165268934113682D-04,
     7  0.2104502626141685D-04,0.1932177243123971D-04,
     8  0.1699579799036191D-04,0.1449749946844336D-04,
     9  0.1212319790076449D-04,0.1004202025010195D-04,
     *  0.8329554045524576D-05,0.7014022597165712D-05,
     1  0.6134664625292422D-05,0.5796323213129145D-05,
     2  0.6133208649643305D-05,0.7197047291055305D-05,
     3  0.9023358466780264D-05,0.1178636970293951D-04,
     4  0.1586867027184345D-04,0.2193854116696026D-04,
     5  0.3113567616641700D-04,0.4546780920883743D-04,
     6  0.6863647276884000D-04,0.1078014957586258D-03,
     7  0.1774797882953124D-03,0.3083448635529795D-03,
     8  0.5669209279544596D-03,0.1096993002605413D-02,
     9  0.2194346345769761D-02,0.4403925271026194D-02,
     *  0.8555357355654279D-02,0.1555243766496481D-01,
     1  0.2578722513565207D-01,0.3845295674668967D-01,
     2  0.5139918510693184D-01,0.6182668764037884D-01,
     3  0.6735614640568536D-01,0.6675929126906305D-01,
     4  0.6005656333399425D-01,0.4816588573596444D-01,
     5  0.3244164095831783D-01,0.1435546853773479D-01/
c        

        do 1000 i=1,nqs(1)
        roots(i) =r1(i)
        weights(i) = w1(i)
 1000 continue

        nquads = nqs(1)


        return
        end
c
c
c
c
c
        subroutine lap_load_quads8(alpha,roots,weights,nquads)
        implicit real *8 (a-h,o-z)
        dimension r1(50),
     1          w1(50),nqs(1),
     2          roots(1),weights(1)

        data nqs/50/

        data r1/
     1  0.1405979088651772D+00,0.3770664007828496D-01,
     2  0.1941793256604013D+00,0.2175583310843908D-01,
     3  0.9654786373850092D-01,0.2552470707435710D+00,
     4  0.1236332680710957D-01,0.4286140777277777D-04,
     5  0.7191666745759328D-02,0.3192274455987489D+00,
     6  0.7083283671594259D-04,0.4420341039543940D-02,
     7  0.2328020152641375D-04,0.6230360257350225D-01,
     8  0.3802439422382859D+00,0.2919717912403516D-02,
     9  0.1607010200161103D-03,0.2077442219690008D-02,
     *  0.2180154724931538D-03,0.1581829331194298D-02,
     1  0.4326250818264827D+00,0.1103982085372605D-03,
     2  0.2778116618466392D-03,0.1275807436226885D-02,
     3  0.3360332933719399D-03,0.1005971879020228D-04,
     4  0.1078377797009089D-02,0.3897853814095172D-03,
     5  0.9460826029546164D-03,0.4374723284515997D-03,
     6  0.8545158221417801D-03,0.4717856598438351D+00,
     7  0.4785905760809708D-03,0.7893271994716224D-03,
     8  0.3114809122228806D-05,0.5133954510652909D-03,
     9  0.6000473828925558D-07,0.6091726713874749D-06,
     *  0.7417082750616985D-03,0.1708461046196341D-08,
     1  0.5425823238575027D-03,0.7060378917995402D-03,
     2  0.5670467905479385D-03,0.4945618116686064D+00,
     3  0.6785914126891615D-03,0.5877503017580200D-03,
     4  0.6567698829810473D-03,0.6057118691627078D-03,
     5  0.6385617059876145D-03,0.6221360809766401D-03/
     
        data w1/
     1  0.4894506663166858D-01,0.1999964001243573D-01,
     2  0.5789317787374240D-01,0.1226810761418808D-01,
     3  0.3914238932800253D-01,0.6345294830302086D-01,
     4  0.6920262378167871D-02,0.2308574149376223D-04,
     5  0.3728440957997423D-02,0.6350164808034770D-01,
     6  0.3353083890218712D-04,0.1998794074925652D-02,
     7  0.1641505848087259D-04,0.2935565845690945D-01,
     8  0.5757190333127644D-01,0.1100686389381919D-02,
     9  0.5454855539930690D-04,0.6334239035577014D-03,
     *  0.5930400654031910D-04,0.3827461347838019D-03,
     1  0.4642132273686547D-01,0.4540414888657398D-04,
     2  0.5960090102036333D-04,0.2422018879306201D-03,
     3  0.5635814593550623D-04,0.9949233023493192D-05,
     4  0.1596154926904850D-03,0.5088410914109870D-04,
     5  0.1089034664469521D-03,0.4441173754295432D-04,
     6  0.7655157218691186D-04,0.3137181715709651D-01,
     7  0.3787276318240238D-04,0.5525687797548835D-04,
     8  0.4285159203908129D-05,0.3185692782864628D-04,
     9  0.1674355052670554D-06,0.1171532032352908D-05,
     *  0.4089750524241331D-04,0.7911337412744574D-08,
     1  0.2666836761244911D-04,0.3105183895855283D-04,
     2  0.2242101179122531D-04,0.1390629982849332D-01,
     3  0.2426561170187988D-04,0.1915353920062464D-04,
     4  0.1970551253736632D-04,0.1696691704685046D-04,
     5  0.1701322667370748D-04,0.1613977368659567D-04/
        
        do 1000 i=1,nqs(1)
        roots(i) =r1(i)
        weights(i) = w1(i)
 1000 continue

        nquads = nqs(1)

        return
        end
c
c
c
c
c
        subroutine lap_load_quads9(alpha,roots,weights,nquads)
        implicit real *8 (a-h,o-z)
        dimension r1(49),
     1          w1(49),nqs(1),
     2          roots(1),weights(1)
cc
        data nqs/49/

        data r1/
     1  0.1727795042124429D-08,0.7055780949695327D-07,
     2  0.8204922812814951D-06,0.4774443272337531D-05,
     3  0.1748749317376857D-04,0.4597534784159532D-04,
     4  0.9522328353437443D-04,0.1676483863597115D-03,
     5  0.2653792851736909D-03,0.3879681616101130D-03,
     6  0.5289757810116877D-03,0.6780042431057753D-03,
     7  0.8246779488161410D-03,0.9611551041015434D-03,
     8  0.1082855889694182D-02,0.1188084213572511D-02,
     9  0.1277216367940568D-02,0.1351877099723081D-02,
     *  0.1414298750462284D-02,0.1466943805852303D-02,
     1  0.1512452108591087D-02,0.1553950395742402D-02,
     2  0.1595433613093840D-02,0.1641530987931195D-02,
     3  0.1697012910644208D-02,0.1767147353776450D-02,
     4  0.1858768630261484D-02,0.1981689124264502D-02,
     5  0.2150734897488111D-02,0.2389170679889077D-02,
     6  0.2734934470727017D-02,0.3252472926402806D-02,
     7  0.4055841476477848D-02,0.5354359800833127D-02,
     8  0.7541341856959388D-02,0.1135442562026047D-01,
     9  0.1811760696065678D-01,0.2997494018925237D-01,
     *  0.4983047699883980D-01,0.8061174234631435D-01,
     1  0.1215152847426034D+00,0.1285320730794659D+00,
     2  0.1804175243850856D+00,0.2444817212351330D+00,
     3  0.3115884805070202D+00,0.3753737009225306D+00,
     4  0.4299695980415889D+00,0.4706984214072730D+00,
     *  0.4943549241597983D+00/
     
        data w1/
     1  0.8286628665545005D-08,0.2068119874799936D-06,
     2  0.1686241629691844D-05,0.7174508338963513D-05,
     3  0.1949510498352906D-04,0.3832515672788459D-04,
     4  0.6050236871873671D-04,0.8475647027939068D-04,
     5  0.1107116673831870D-03,0.1333416518933767D-03,
     6  0.1468666781193831D-03,0.1494193908056084D-03,
     7  0.1426108409418687D-03,0.1295881479575601D-03,
     8  0.1135452402799442D-03,0.9698827138897824D-04,
     9  0.8155765960038281D-04,0.6814228196292226D-04,
     *  0.5711251138971149D-04,0.4861200161905559D-04,
     1  0.4291683051717038D-04,0.4075014694221559D-04,
     2  0.4300046752710092D-04,0.4998187878353866D-04,
     3  0.6184106554099772D-04,0.7954418829593741D-04,
     4  0.1053015760220647D-03,0.1429564558546512D-03,
     5  0.1989147276889356D-03,0.2840949918142403D-03,
     6  0.4178138938733534D-03,0.6355593449063268D-03,
     7  0.1004660832173323D-02,0.1655312546612127D-02,
     8  0.2837399013093792D-02,0.5004839081261214D-02,
     9  0.8878766557478632D-02,0.1533177363117426D-01,
     *  0.2489404595904059D-01,0.3692664893680149D-01,
     1  0.2889967620437316D-01,0.2191710771422538D-01,
     2  0.6047571253813604D-01,0.6664003254584255D-01,
     3  0.6649564536652504D-01,0.6008501095928537D-01,
     4  0.4832299786819307D-01,0.3260008908118094D-01,
     *  0.1443695430420016D-01/
        
        do 1000 i=1,nqs(1)
        roots(i) =r1(i)
        weights(i) = w1(i)
 1000 continue

        nquads = nqs(1)

        return
        end
c
c
c
c
c
        subroutine lap_load_quads10(alpha,roots,weights,nquads)
        implicit real *8 (a-h,o-z)
        dimension r1(48),
     1          w1(48),nqs(1),
     2          roots(1),weights(1)
cc
        data nqs/48/

        data r1/
     1  0.4915940565120498D-09,0.6366396549567589D-07,
     2  0.1052464363706094D-05,0.7300302479610386D-05,
     3  0.2968948216653074D-04,0.8374552603837880D-04,
     4  0.1815960220375512D-03,0.3240508642093744D-03,
     5  0.5033387055626335D-03,0.7239299816886263D-03,
     6  0.9958828101356858D-03,0.1307019177637230D-02,
     7  0.1632399162769291D-02,0.1948693627393997D-02,
     8  0.2239499413984638D-02,0.2496061719284318D-02,
     9  0.2715912339111025D-02,0.2900857855762957D-02,
     *  0.3055130702663833D-02,0.3184064164333204D-02,
     1  0.3293462051057460D-02,0.3389820196767198D-02,
     2  0.3481258876289388D-02,0.3577847991491768D-02,
     3  0.3690443379176457D-02,0.3830397837851448D-02,
     4  0.4011415512432056D-02,0.4252435231767007D-02,
     5  0.4581448239246354D-02,0.5041441959563715D-02,
     6  0.5700623983079858D-02,0.6670539652590317D-02,
     7  0.8137516276064199D-02,0.1041001215283260D-01,
     8  0.1393645222837954D-01,0.1890236811336231D-01,
     9  0.2474106776375934D-01,0.3472070011828413D-01,
     *  0.5277218713175203D-01,0.8143197429207830D-01,
     1  0.1227106581061434D+00,0.1765183085274055D+00,
     2  0.2398322821485779D+00,0.3071733685411125D+00,
     3  0.3719474260180947D+00,0.4278431339790986D+00,
     4  0.4697550523393587D+00,0.4941678022673460D+00/
     
        data w1/
     1  0.3734857781471013D-08,0.2201276804028094D-06,
     2  0.2410400205099120D-05,0.1192610908354360D-04,
     3  0.3559339783298433D-04,0.7469031996163567D-04,
     4  0.1211550132989581D-03,0.1619668043120557D-03,
     5  0.1974159860458423D-03,0.2462054923359959D-03,
     6  0.2953086654337556D-03,0.3225710439284876D-03,
     7  0.3242247069185761D-03,0.3056233678592896D-03,
     8  0.2745551774847368D-03,0.2381934584015340D-03,
     9  0.2018404364394497D-03,0.1687784048088391D-03,
     *  0.1406670912544534D-03,0.1181636626220594D-03,
     1  0.1016959085331816D-03,0.9236317562600105D-04,
     2  0.9222136834937868D-04,0.1027741656345692D-03,
     3  0.1242729141460644D-03,0.1578722212405685D-03,
     4  0.2072589304020128D-03,0.2793515308919256D-03,
     5  0.3856769394109917D-03,0.5453900224748489D-03,
     6  0.7910555761172120D-03,0.1179002688616872D-02,
     7  0.1805530538091293D-02,0.2817469870941902D-02,
     8  0.4297766040163870D-02,0.5404581309340670D-02,
     9  0.6969551766395711D-02,0.1360045625151712D-01,
     *  0.2292630728296836D-01,0.3476033821705765D-01,
     1  0.4779697664210619D-01,0.5930531102053836D-01,
     2  0.6638483045477276D-01,0.6717370658341841D-01,
     3  0.6130707893088744D-01,0.4962768303863527D-01,
     4  0.3360960510372391D-01,0.1491235810723200D-01/
        
        do 1000 i=1,nqs(1)
        roots(i) =r1(i)
        weights(i) = w1(i)
 1000 continue

        nquads = nqs(1)

        return
        end
c
c
c
c
c
        subroutine lap_load_quads11(alpha,roots,weights,nquads)
        implicit real *8 (a-h,o-z)
        dimension r1(48),
     1          w1(48),nqs(1),
     2          roots(1),weights(1)
cc
        data nqs/48/

        data r1/
     1  0.1925446125780938D-02,0.2258394007662382D-02,
     2  0.1491888298549899D-02,0.6986829045448299D-01,
     3  0.1018677651355209D-02,0.1009134366881154D+00,
     4  0.4838882942837815D-01,0.2799798188132364D-02,
     5  0.3419769167436691D-01,0.6246842215566839D-03,
     6  0.1429452996738215D+00,0.4688033344762156D-08,
     7  0.1931103166960575D-06,0.2504458501383201D-01,
     8  0.2304719548796047D-05,0.3353761161528207D-03,
     9  0.1955454707009454D+00,0.1394541001925047D-04,
     *  0.1916266290847425D-01,0.1511103255017971D-03,
     1  0.5385085924183137D-04,0.2559520443428106D+00,
     2  0.1533989134742671D-01,0.1280546650245821D-01,
     3  0.3193959435067395D+00,0.1108556401735481D-01,
     4  0.4006767355971183D-02,0.3800996231015624D+00,
     5  0.4555018642454737D-02,0.9890570693966812D-02,
     6  0.5040057509743092D-02,0.4324138419020599D+00,
     7  0.9041350776173079D-02,0.5457666233386833D-02,
     8  0.3410791504290598D-02,0.8424684720671201D-02,
     9  0.5811095183276898D-02,0.4716563767881768D+00,
     *  0.7967136683548044D-02,0.6107890124473716D-02,
     1  0.7619729791781494D-02,0.6357632015591240D-02,
     2  0.4945323723069244D+00,0.7348635713354441D-02,
     3  0.6570892084945908D-02,0.6939382925596005D-02,
     4  0.6759719516492919D-02,0.7128883480429202D-02/
     
        data w1/
     1  0.3300658706524506D-03,0.4372264257247486D-03,
     2  0.4904941470843468D-03,0.2591918586489932D-01,
     3  0.4409551654215699D-03,0.3642594370995939D-01,
     4  0.1743483841507037D-01,0.5997963247286367D-03,
     5  0.1132987383371117D-01,0.3433781771503504D-03,
     6  0.4757747427425636D-01,0.2246157967200183D-07,
     7  0.5696819810126595D-06,0.7271436445978347D-02,
     8  0.4811500746797552D-05,0.2352058488348148D-03,
     9  0.5715641346267642D-01,0.2159495805776038D-04,
     *  0.4692144358815966D-02,0.1364399282591349D-03,
     1  0.6324791600422571D-04,0.6284377875572660D-01,
     2  0.3079335403131193D-02,0.2066602907991158D-02,
     3  0.6306102036776201D-01,0.1420271131398016D-02,
     4  0.5758461173048052D-03,0.5739057029319380D-01,
     5  0.5181513356334551D-03,0.9988553574941079D-03,
     6  0.4512865550954635D-03,0.4645091802065651D-01,
     7  0.7180366752977967D-03,0.3845263818102377D-03,
     8  0.6109276265474561D-03,0.5272902787934343D-03,
     9  0.3236274551614849D-03,0.3148306855767780D-01,
     *  0.3958202859364866D-03,0.2715775288352290D-03,
     1  0.3045311008650685D-03,0.2296603702305503D-03,
     2  0.1397883947423124D-01,0.2417357348951265D-03,
     3  0.1988326394886317D-03,0.1812008825396122D-03,
     4  0.1813490302665854D-03,0.2012209604432690D-03/
        
        do 1000 i=1,nqs(1)
        roots(i) =r1(i)
        weights(i) = w1(i)
 1000 continue

        nquads = nqs(1)

        return
        end
c
c
c
c
c
        subroutine lap_load_quads12(alpha,roots,weights,nquads)
        implicit real *8 (a-h,o-z)
        dimension r1(47),
     1          w1(47),nqs(1),
     2          roots(1),weights(1)
cc
        data nqs/47/

        data r1/
     1  0.4122667564436619D-02,0.5019913246991335D-02,
     2  0.9302669356486260D-08,0.2094461217665387D-06,
     3  0.2134808297935885D+00,0.2659532100636430D+00,
     4  0.2635441174011053D-02,0.1646989244009733D+00,
     5  0.1755347867378724D-02,0.3225032036154114D+00,
     6  0.2174467172038479D-03,0.1038728230898043D-02,
     7  0.5267372559930134D-03,0.6832958107877562D-04,
     8  0.1218080197241446D+00,0.6133836754766901D-02,
     9  0.3798993774227389D+00,0.1504537702480047D-04,
     *  0.8786745515265810D-01,0.2160953165893883D-05,
     1  0.4314619959225519D+00,0.6335515921279134D-01,
     2  0.3529681952615800D-02,0.4668649077872851D-01,
     3  0.7235268765774187D-02,0.3568076990267976D-01,
     4  0.4710467319454186D+00,0.2844990197713470D-01,
     5  0.8252530602335322D-02,0.2364695642853491D-01,
     6  0.2039558809146438D-01,0.9153337015165129D-02,
     7  0.4943943917704765D+00,0.1814602944473234D-01,
     8  0.9928322850403901D-02,0.1655487789793591D-01,
     9  0.1540476444516984D-01,0.1058305821920502D-01,
     *  0.1266065985595966D-01,0.1455505421665963D-01,
     1  0.1113163319739717D-01,0.1391235730139694D-01,
     2  0.1233111396201624D-01,0.1300843868477772D-01,
     3  0.1341250834319209D-01,0.1159211728775005D-01,
     *  0.1198442133356056D-01/
     
        data w1/
     1  0.6340179932872305D-03,0.1073558706262916D-02,
     2  0.3841878788971486D-07,0.5303201068159114D-06,
     3  0.5064602473010912D-01,0.5455368981012535D-01,
     4  0.9328818929366855D-03,0.4646964508471015D-01,
     5  0.8091473857784379D-03,0.5796747338884250D-01,
     6  0.2193063355550513D-03,0.6175436323249406D-03,
     7  0.4066963319038825D-03,0.9024382035728696D-04,
     8  0.3873910742282586D-01,0.1124784686716045D-02,
     9  0.5563323492649392D-01,0.2587431353812073D-04,
     *  0.2908180858079598D-01,0.4655581117814278D-05,
     1  0.4645246086420830D-01,0.2023413086975272D-01,
     2  0.7551386675848550D-03,0.1348122457315573D-01,
     3  0.1067115685961143D-02,0.8848104251424968D-02,
     4  0.3200685193976228D-01,0.5836352888365041D-02,
     5  0.9622088358706651D-03,0.3913001251665589D-02,
     6  0.2679343673976938D-02,0.8379727874914398D-03,
     7  0.1431981398042963D-01,0.1875696946639954D-02,
     8  0.7130357115237840D-03,0.1342040650884645D-02,
     9  0.9811534550916436D-03,0.5988601553434464D-03,
     *  0.3323659553628552D-03,0.7335445394123155D-03,
     1  0.5013385092936412D-03,0.5623498591054003D-03,
     2  0.3327259840812139D-03,0.3695666516436899D-03,
     3  0.4450341126562244D-03,0.4229352069279540D-03,
     *  0.3653686298079470D-03/
        
        do 1000 i=1,nqs(1)
        roots(i) =r1(i)
        weights(i) = w1(i)
 1000 continue

        nquads = nqs(1)

        return
        end
c
c
c
c
c
        subroutine lap_load_quads13(alpha,roots,weights,nquads)
        implicit real *8 (a-h,o-z)
        dimension r1(46),
     1          w1(46),nqs(1),
     2          roots(1),weights(1)
cc
        data nqs/46/

        data r1/
     1  0.1086450215961553D-07,0.4259798956626428D-06,
     2  0.4914665902307366D-05,0.2939230285012116D-04,
     3  0.1155069251998621D-03,0.3385223522102528D-03,
     4  0.7946698817111337D-03,0.1564119724668502D-02,
     5  0.2673815505030739D-02,0.4089355961242245D-02,
     6  0.5737346664453813D-02,0.7532990712720938D-02,
     7  0.9389775685269570D-02,0.1122013700011049D-01,
     8  0.1294563703975584D-01,0.1451077447311411D-01,
     9  0.1588820952773429D-01,0.1707490649591890D-01,
     *  0.1808441981883518D-01,0.1893943513473061D-01,
     1  0.1966658062849345D-01,0.2029450836487056D-01,
     2  0.2085598524567319D-01,0.2139292650972987D-01,
     3  0.2195754794921503D-01,0.2260564808605361D-01,
     4  0.2339393729866213D-01,0.2438793145253039D-01,
     5  0.2567453149586838D-01,0.2737724628653877D-01,
     6  0.2967725156784978D-01,0.3284575609979534D-01,
     7  0.3729506001254372D-01,0.4365673431668617D-01,
     8  0.5289049654914075D-01,0.6640342021585916D-01,
     9  0.8609569464320100D-01,0.1141513779209429D+00,
     *  0.1523728472828871D+00,0.2011224183780866D+00,
     1  0.2584242368528969D+00,0.3199135191227240D+00,
     2  0.3797499539784828D+00,0.4319344156552361D+00,
     3  0.4713787941489969D+00,0.4944710736167376D+00/
     
        data w1/
     1  0.5143751984777072D-07,0.1238693371455580D-05,
     2  0.1014071464440428D-04,0.4558412283325914D-04,
     3  0.1393766831612646D-03,0.3233816114056010D-03,
     4  0.6029813512429334D-03,0.9405447756352423D-03,
     5  0.1272669686362447D-02,0.1545644005518712D-02,
     6  0.1736077375162409D-02,0.1840848161440747D-02,
     7  0.1857927349249967D-02,0.1789455503909413D-02,
     8  0.1652211264237135D-02,0.1473586434688929D-02,
     9  0.1280916065077064D-02,0.1094865504034811D-02,
     *  0.9280305522722333D-03,0.7864721353490130D-03,
     1  0.6725595178887503D-03,0.5886554320568561D-03,
     2  0.5412413912148163D-03,0.5415023276782697D-03,
     3  0.5971018982977337D-03,0.7083934393547809D-03,
     4  0.8788918974321484D-03,0.1123334989713599D-02,
     5  0.1470111814678849D-02,0.1964991960866214D-02,
     6  0.2679390903219548D-02,0.3725047052303507D-02,
     7  0.5276975962725254D-02,0.7604160875196279D-02,
     8  0.1109666353227845D-01,0.1624876129896401D-01,
     9  0.2351089119474069D-01,0.3291998095691848D-01,
     *  0.4360033497335705D-01,0.5357067033805278D-01,
     1  0.6028862130184469D-01,0.6169284892794251D-01,
     2  0.5696175043706891D-01,0.4655018635018205D-01,
     3  0.3173392287377276D-01,0.1413100492513413D-01/
        
        do 1000 i=1,nqs(1)
        roots(i) =r1(i)
        weights(i) = w1(i)
 1000 continue

        nquads = nqs(1)

        return
        end
c
c
c
c
c
        subroutine lap_load_quads14(alpha,roots,weights,nquads)
        implicit real *8 (a-h,o-z)
        dimension r1(45),
     1          w1(45),nqs(1),
     2          roots(1),weights(1)

        data nqs/45/

        data r1/
     1  0.1209214084757675D-07,0.5386887323091141D-06,
     2  0.6914899222682677D-05,0.4499414480698739D-04,
     3  0.1864859045968508D-03,0.5585284453034897D-03,
     4  0.1311730072318589D-02,0.2559639934799319D-02,
     5  0.4331776310916648D-02,0.6567847664311166D-02,
     6  0.9149367751358062D-02,0.1194314279059931D-01,
     7  0.1482413172776496D-01,0.1767267597147750D-01,
     8  0.2037655463016900D-01,0.2284818153572837D-01,
     9  0.2503787614440809D-01,0.2693327906689178D-01,
     *  0.2854991725194899D-01,0.2992010609382978D-01,
     1  0.3108445312711567D-01,0.3208820616204719D-01,
     2  0.3298403495250103D-01,0.3383980376784996D-01,
     3  0.3474028167265262D-01,0.3577608497955288D-01,
     4  0.3703922141417783D-01,0.3863564092707623D-01,
     5  0.4070484862809234D-01,0.4344289812151313D-01,
     6  0.4713303084400109D-01,0.5219038381222351D-01,
     7  0.5922693301582893D-01,0.6913777309065552D-01,
     8  0.8319088698561906D-01,0.1030556492425053D+00,
     9  0.1306313894378450D+00,0.1675079894224874D+00,
     *  0.2140646158249321D+00,0.2686138443402795D+00,
     1  0.3272192499327028D+00,0.3844256992616677D+00,
     2  0.4344819215686682D+00,0.4724195498045815D+00,
     *  0.4946687856897223D+00/
     
        data w1/
     1  0.5898043746373717D-07,0.1628875359838292D-05,
     2  0.1493858521153569D-04,0.7292732074602620D-04,
     3  0.2318152160063924D-03,0.5385046772668303D-03,
     4  0.9876776321881342D-03,0.1513141793148529D-02,
     5  0.2020492648340146D-02,0.2431083336014109D-02,
     6  0.2709351941936581D-02,0.2857516641756599D-02,
     7  0.2884509165485648D-02,0.2793518667020434D-02,
     8  0.2599398203502705D-02,0.2335604044100154D-02,
     9  0.2041901026331754D-02,0.1751687021185993D-02,
     *  0.1487121451861872D-02,0.1260089199351727D-02,
     1  0.1076100841402537D-02,0.9400253378743847D-03,
     2  0.8628749960576979D-03,0.8630513684273161D-03,
     3  0.9531085764813745D-03,0.1133581585025577D-02,
     4  0.1410017728254213D-02,0.1805701366281027D-02,
     5  0.2364894714111216D-02,0.3157682202674933D-02,
     6  0.4290697911040822D-02,0.5924535040792043D-02,
     7  0.8296111079238045D-02,0.1173608686973799D-01,
     8  0.1665074158199834D-01,0.2340420817726877D-01,
     9  0.3202829907015263D-01,0.4180418033096515D-01,
     *  0.5103159417817041D-01,0.5739736933662250D-01,
     1  0.5888210654137527D-01,0.5455348682162766D-01,
     2  0.4472096530038234D-01,0.3055594191659485D-01,
     *  0.1362367070019044D-01/
        
        do 1000 i=1,nqs(1)
        roots(i) =r1(i)
        weights(i) = w1(i)
 1000 continue

        nquads = nqs(1)

        return
        end
c
c
c
c
c
        subroutine lap_load_quads15(alpha,roots,weights,nquads)
        implicit real *8 (a-h,o-z)
        dimension r1(44),
     1          w1(44),nqs(1),
     2          roots(1),weights(1)

        data nqs/44/

        data r1/
     1  0.3068852684425486D-07,0.1212671383124489D-05,
     2  0.1423002575082112D-04,0.8625096948018828D-04,
     3  0.3376557514025626D-03,0.9660359778062236D-03,
     4  0.2188057673781166D-02,0.4151584897751001D-02,
     5  0.6879088821969035D-02,0.1026885188539555D-01,
     6  0.1413936107215029D-01,0.1828795942074426D-01,
     7  0.2253028034892666D-01,0.2670327429245961D-01,
     8  0.3065954889979968D-01,0.3428103462312016D-01,
     9  0.3749725605230624D-01,0.4028849169535272D-01,
     *  0.4267504267891529D-01,0.4470222758801128D-01,
     1  0.4642820538449870D-01,0.4791858268754365D-01,
     2  0.4925040692840628D-01,0.5052344257513526D-01,
     3  0.5186245329938375D-01,0.5340066977354320D-01,
     4  0.5527282742856303D-01,0.5763296771235707D-01,
     5  0.6068209996557827D-01,0.6469950270163515D-01,
     6  0.7008224541477922D-01,0.7739877121157645D-01,
     7  0.8745809484636896D-01,0.1013820952027924D+00,
     8  0.1206341017298539D+00,0.1468982817776978D+00,
     9  0.1816616192572646D+00,0.2254459454912315D+00,
     *  0.2769597009776597D+00,0.3327412687966288D+00,
     1  0.3876723613337215D+00,0.4361138697687636D+00,
     2  0.4730440339162748D+00,0.4947827402312563D+00/
     
        data w1/
     1  0.1453240456935998D-06,0.3543107568892388D-05,
     2  0.2961464279769124D-04,0.1343246450919961D-03,
     3  0.4024824571819176D-03,0.8920334609848348D-03,
     4  0.1577658932180290D-02,0.2353069222940742D-02,
     5  0.3084177183789030D-02,0.3664034450979075D-02,
     6  0.4042484022095170D-02,0.4224028115335730D-02,
     7  0.4233628246173332D-02,0.4087501196712291D-02,
     8  0.3805008662775037D-02,0.3426176578804992D-02,
     9  0.3003185255639956D-02,0.2582936070381048D-02,
     *  0.2197930259153317D-02,0.1866229306636261D-02,
     1  0.1596595671921814D-02,0.1396769926692256D-02,
     2  0.1283410273771860D-02,0.1283842165291251D-02,
     3  0.1416540280360364D-02,0.1681989746166049D-02,
     4  0.2087528910666046D-02,0.2665746188696454D-02,
     5  0.3478421370320880D-02,0.4621676408500669D-02,
     6  0.6237374461824706D-02,0.8529209052448323D-02,
     7  0.1177569190636902D-01,0.1631780226272767D-01,
     8  0.2247362835698060D-01,0.3031597718232067D-01,
     9  0.3931341657378016D-01,0.4804706581743547D-01,
     *  0.5439421726394086D-01,0.5629516662437995D-01,
     1  0.5260712377297149D-01,0.4342699940553210D-01,
     2  0.2981521170497037D-01,0.1332840353066368D-01/
        
        do 1000 i=1,nqs(1)
        roots(i) =r1(i)
        weights(i) = w1(i)
 1000 continue

        nquads = nqs(1)

        return
        end
c
c
c
c
c
        subroutine lap_load_quads16(alpha,roots,weights,nquads)
        implicit real *8 (a-h,o-z)
        dimension r1(43),
     1          w1(43),nqs(1),
     2          roots(1),weights(1)

        data nqs/43/

        data r1/
     1  0.2906395198270554D-07,0.1226298122149951D-05,
     2  0.1526985926019358D-04,0.9770003808463833D-04,
     3  0.4017971371391552D-03,0.1201761001702416D-02,
     4  0.2830042552815491D-02,0.5547162962591285D-02,
     5  0.9429259816863713D-02,0.1434141114706081D-01,
     6  0.1999572231097785D-01,0.2604741420164713D-01,
     7  0.3217842283643095D-01,0.3813601790280760D-01,
     8  0.4373049696220843D-01,0.4882669797508304D-01,
     9  0.5334681055372111D-01,0.5727192176870551D-01,
     *  0.6063231482172002D-01,0.6349086786166625D-01,
     1  0.6592805580848098D-01,0.6803508625119305D-01,
     2  0.6991963627447038D-01,0.7172156927320920D-01,
     3  0.7361591704222425D-01,0.7578937115756686D-01,
     4  0.7842991059493949D-01,0.8175074565882148D-01,
     5  0.8602720411043915D-01,0.9163652099349241D-01,
     6  0.9910426402106750D-01,0.1091600849510428D+00,
     7  0.1227945083678970D+00,0.1412829329155913D+00,
     8  0.1660981074884974D+00,0.1985896083657034D+00,
     9  0.2393525770054652D+00,0.2874403608213202D+00,
     *  0.3398886646612119D+00,0.3920027202678110D+00,
     1  0.4383500742342197D+00,0.4739180849610615D+00,
     *  0.4949442997845796D+00/
     
        data w1/
     1  0.1396480939055788D-06,0.3659195965436240D-05,
     2  0.3264357000701933D-04,0.1572018804883812D-03,
     3  0.4977844919171622D-03,0.1159969080708331D-02,
     4  0.2142074646301428D-02,0.3305282029462255D-02,
     5  0.4435324841552823D-02,0.5339268921849633D-02,
     6  0.5910410267963488D-02,0.6139325831315434D-02,
     7  0.6080720996821960D-02,0.5803040502952179D-02,
     8  0.5363292240968441D-02,0.4816089191657522D-02,
     9  0.4221162355635414D-02,0.3634406324309255D-02,
     *  0.3097124253510761D-02,0.2633520081551853D-02,
     1  0.2255964219507831D-02,0.1975721704688396D-02,
     2  0.1816558642041303D-02,0.1817026900557232D-02,
     3  0.2003012740602155D-02,0.2374738102786431D-02,
     4  0.2941245172455494D-02,0.3745664484488182D-02,
     5  0.4869360511599342D-02,0.6435924552297926D-02,
     6  0.8620152295196089D-02,0.1165599388614271D-01,
     7  0.1582655970129308D-01,0.2140081490035503D-01,
     8  0.2846528493547270D-01,0.3663128901063857D-01,
     9  0.4474615136827295D-01,0.5093043196102848D-01,
     *  0.5316559404739467D-01,0.5013648448519662D-01,
     1  0.4170881974351262D-01,0.2879368857293260D-01,
     *  0.1291107770850691D-01/
        
        do 1000 i=1,nqs(1)
        roots(i) =r1(i)
        weights(i) = w1(i)
 1000 continue
        
        nquads = nqs(1)
        
        return
        end
c
c
c
c
c
        subroutine lap_load_quads17(alpha,roots,weights,nquads)
        implicit real *8 (a-h,o-z)
        dimension r1(42),
     1          w1(42),nqs(1),
     2          roots(1),weights(1)

        data nqs/42/

        data r1/
     1  0.4741819855613529D-07,0.1909837853055347D-05,
     2  0.2293202596174759D-04,0.1427908317335751D-03,
     3  0.5758816012344594D-03,0.1698726595621155D-02,
     4  0.3960807523393180D-02,0.7708135020217708D-02,
     5  0.1303531010382125D-01,0.1975285668267597D-01,
     6  0.2746458300913045D-01,0.3569401860747725D-01,
     7  0.4399934453148574D-01,0.5203719728679239D-01,
     8  0.5956639418445240D-01,0.6642474710619967D-01,
     9  0.7251816363405521D-01,0.7782119253832387D-01,
     *  0.8237002917822910D-01,0.8624490635414337D-01,
     1  0.8955149185635303D-01,0.9241155324145977D-01,
     2  0.9497021121953998D-01,0.9741676036616651D-01,
     3  0.9998834412824736D-01,0.1029376921933685D+00,
     4  0.1065184006162462D+00,0.1110160748990269D+00,
     5  0.1167954231764001D+00,0.1243476720773044D+00,
     6  0.1343391873907688D+00,0.1476567461295584D+00,
     7  0.1654244549226208D+00,0.1889312477455086D+00,
     8  0.2193691823426130D+00,0.2573016001034534D+00,
     9  0.3019522985703886D+00,0.3506989935461467D+00,
     *  0.3992368331680128D+00,0.4424722558432686D+00,
     1  0.4756673966177014D+00,0.4952847099997786D+00/
     
        data w1/
     1  0.2253729143583074D-06,0.5620006014996524D-05,
     2  0.4829935761391062D-04,0.2264655063177952D-03,
     3  0.7040336061984888D-03,0.1619848451986785D-02,
     4  0.2964682149970369D-02,0.4546248368045184D-02,
     5  0.6075060982810495D-02,0.7291743889330022D-02,
     6  0.8050286429664157D-02,0.8333851071492826D-02,
     7  0.8219607880655869D-02,0.7816410156092943D-02,
     8  0.7215113836217549D-02,0.6485838086485585D-02,
     9  0.5696867985400302D-02,0.4915544521501199D-02,
     *  0.4195835651213109D-02,0.3571704056278007D-02,
     1  0.3061645491674985D-02,0.2682236908180399D-02,
     2  0.2466457296822878D-02,0.2466913179884261D-02,
     3  0.2718696454484927D-02,0.3221625044270126D-02,
     4  0.3986581905099133D-02,0.5068771895712780D-02,
     5  0.6570979467554498D-02,0.8643772200710563D-02,
     6  0.1148663192156741D-01,0.1533678979790718D-01,
     7  0.2041755142890679D-01,0.2680371314883177D-01,
     8  0.3417888302719021D-01,0.4156639748261635D-01,
     9  0.4729545439831023D-01,0.4946416527748495D-01,
     *  0.4674197450884830D-01,0.3892482496193575D-01,
     1  0.2686991023684446D-01,0.1204273659895809D-01/
        
        do 1000 i=1,nqs(1)
        roots(i) =r1(i)
        weights(i) = w1(i)
 1000 continue

        nquads = nqs(1)

        return
        end
c
c
c
c
c
        subroutine lap_load_quads18(alpha,roots,weights,nquads)
        implicit real *8 (a-h,o-z)
        dimension r1(41),
     1          w1(41),nqs(1),
     2          roots(1),weights(1)

        data nqs/41/

        data r1/
     1  0.6089526690320087D-07,0.2509174394039132D-05,
     2  0.3057804898339133D-04,0.1918167432774567D-03,
     3  0.7748340908926795D-03,0.2281075351231155D-02,
     4  0.5300442609208159D-02,0.1027940208137397D-01,
     5  0.1733291798066895D-01,0.2620521042541515D-01,
     6  0.3636959716374256D-01,0.4719020052485871D-01,
     7  0.5807075954293105D-01,0.6854649152832953D-01,
     8  0.7830528220066424D-01,0.8715987788175939D-01,
     9  0.9501586433693553D-01,0.1018569007396839D+00,
     *  0.1077345631826304D+00,0.1127511106308734D+00,
     1  0.1170399605015716D+00,0.1207556614359928D+00,
     2  0.1240837173843032D+00,0.1272671386218090D+00,
     3  0.1306106251490566D+00,0.1344382737216889D+00,
     4  0.1390730172254686D+00,0.1448738178733615D+00,
     5  0.1522909382531084D+00,0.1619160234490404D+00,
     6  0.1745223806031742D+00,0.1910802953849786D+00,
     7  0.2127044266863002D+00,0.2404560738788596D+00,
     8  0.2749205970658421D+00,0.3155879998702825D+00,
     9  0.3603090010616796D+00,0.4052703470305921D+00,
     *  0.4456986950959020D+00,0.4769668939719596D+00,
     *  0.4955296540676741D+00/
     
        data w1/
     1  0.2910626566369943D-06,0.7430737495849879D-05,
     2  0.6473787904730263D-04,0.3050172601115910D-03,
     3  0.9468037026628901D-03,0.2168299978952121D-02,
     4  0.3948257927836743D-02,0.6029424203085728D-02,
     5  0.8033096917530329D-02,0.9620825040693234D-02,
     6  0.1059938814208873D-01,0.1094087353893617D-01,
     7  0.1074247949197136D-01,0.1015756365386686D-01,
     8  0.9329398591114304D-02,0.8364250628090784D-02,
     9  0.7345135421255466D-02,0.6345838593930263D-02,
     *  0.5426885639263515D-02,0.4628639508013114D-02,
     1  0.3974662576973101D-02,0.3487111760378947D-02,
     2  0.3209292002822931D-02,0.3209306980260242D-02,
     3  0.3532121594744482D-02,0.4176218880823104D-02,
     4  0.5152133222556067D-02,0.6523677988341184D-02,
     5  0.8408468912387817D-02,0.1097023309000345D-01,
     6  0.1440488017445641D-01,0.1889978695099859D-01,
     7  0.2453239372103148D-01,0.3107944057705537D-01,
     8  0.3777492940182513D-01,0.4320274931816083D-01,
     9  0.4559105068928297D-01,0.4351600027246052D-01,
     *  0.3655320161144049D-01,0.2538461842300156D-01,
     *  0.1141308393239233D-01/
        
        do 1000 i=1,nqs(1)
        roots(i) =r1(i)
        weights(i) = w1(i)
 1000 continue

        nquads = nqs(1)

        return
        end
c
c
c
c
c
        subroutine lap_load_quads19(alpha,roots,weights,nquads)
        implicit real *8 (a-h,o-z)
        dimension r1(40),
     1          w1(40),nqs(1),
     2          roots(1),weights(1)

        data nqs/40/

        data r1/
     1  0.5869627423615182D-07,0.2619611156978720D-05,
     2  0.3327558357376030D-04,0.2129090539312296D-03,
     3  0.8662474105278688D-03,0.2552249577780637D-02,
     4  0.5925850347171139D-02,0.1150440747296738D-01,
     5  0.1948770641072867D-01,0.2970140860005134D-01,
     6  0.4165200530946879D-01,0.5465617675069200D-01,
     7  0.6799863758148945D-01,0.8107188902390474D-01,
     8  0.9345484610356578D-01,0.1048985691881827D+00,
     9  0.1152540779113703D+00,0.1244332869449220D+00,
     *  0.1324206405081489D+00,0.1392817324914182D+00,
     1  0.1451468865326322D+00,0.1501883005912049D+00,
     2  0.1546143053317489D+00,0.1586922236607729D+00,
     3  0.1627748625105248D+00,0.1672711552701050D+00,
     4  0.1725917555679832D+00,0.1791623469885011D+00,
     5  0.1874864254362437D+00,0.1982055277168548D+00,
     6  0.2121384035253154D+00,0.2302784395921266D+00,
     7  0.2536977759310378D+00,0.2832729179989686D+00,
     8  0.3191645921453985D+00,0.3601415546860665D+00,
     9  0.4031189952941783D+00,0.4433876417809240D+00,
     *  0.4756265528841884D+00,0.4952275129597854D+00/
     
        data w1/
     1  0.2869396106723104D-06,0.7906553654118334D-05,
     2  0.7140465148354611D-04,0.3411754133534956D-03,
     3  0.1061196633366086D-02,0.2424295963757683D-02,
     4  0.4412283070103278D-02,0.6779011661874333D-02,
     5  0.9157526199488777D-02,0.1118535684464750D-01,
     6  0.1259920490224977D-01,0.1328749815101452D-01,
     7  0.1329475524480499D-01,0.1278187508791107D-01,
     8  0.1194406026297122D-01,0.1091923804418038D-01,
     9  0.9776712152258756D-02,0.8578925943530324D-02,
     *  0.7407050138710362D-02,0.6337294303887402D-02,
     1  0.5421844572952366D-02,0.4695285021662129D-02,
     2  0.4200615725416245D-02,0.4014575124129146D-02,
     3  0.4220135032517338D-02,0.4840532229198427D-02,
     4  0.5870428016830244D-02,0.7353646643934701D-02,
     5  0.9400740886038555D-02,0.1217384557139180D-01,
     6  0.1585884982684477D-01,0.2060426551258256D-01,
     7  0.2638929723643502D-01,0.3279998322428375D-01,
     8  0.3878724131883886D-01,0.4264907020603310D-01,
     9  0.4250266343968424D-01,0.3712416084130675D-01,
     *  0.2656724154711847D-01,0.1215851985994324D-01/
        
        do 1000 i=1,nqs(1)
        roots(i) =r1(i)
        weights(i) = w1(i)
 1000 continue

        nquads = nqs(1)

        return
        end
c
c
c
c
c
        subroutine lap_load_quads20(alpha,roots,weights,nquads)
        implicit real *8 (a-h,o-z)
        dimension r1(39),
     1          w1(39),nqs(1),
     2          roots(1),weights(1)

        data nqs/39/

        data r1/
     1  0.2907991640025808D-08,0.7047404547085416D-06,
     2  0.1495406795418753D-04,0.1288158175461680D-03,
     3  0.6430569985794518D-03,0.2202188329699164D-02,
     4  0.5721407513694794D-02,0.1206535391307665D-01,
     5  0.2167605329723243D-01,0.3438954647158920D-01,
     6  0.4951199970062988D-01,0.6605825031793373D-01,
     7  0.8301291030299588D-01,0.9952402089995999D-01,
     8  0.1150009487456899D+00,0.1291178855382610D+00,
     9  0.1417436681392542D+00,0.1528573312116258D+00,
     *  0.1625066175555095D+00,0.1708012913836305D+00,
     1  0.1779057433155208D+00,0.1840257648438098D+00,
     2  0.1894089770005858D+00,0.1943744658498202D+00,
     3  0.1993442401924774D+00,0.2048067709913056D+00,
     4  0.2112491695011393D+00,0.2191678507650295D+00,
     5  0.2291343033647310D+00,0.2418489325326350D+00,
     6  0.2581525888312920D+00,0.2789613232825804D+00,
     7  0.3050615733076412D+00,0.3366999417278026D+00,
     8  0.3729968842795128D+00,0.4114415937860359D+00,
     9  0.4479098931220386D+00,0.4774524047849831D+00,
     *  0.4955704665325898D+00/
     
        data w1/
     1  0.2848709131743317D-07,0.2643008815283230D-05,
     2  0.3762407772971636D-04,0.2374880687113898D-03,
     3  0.8994644222445275D-03,0.2377269153959031D-02,
     4  0.4812149510511885D-02,0.7952626623488830D-02,
     5  0.1123988260698055D-01,0.1406711684135839D-01,
     6  0.1600902848703878D-01,0.1691096461270792D-01,
     7  0.1685474134141653D-01,0.1606901184485692D-01,
     8  0.1483187183691041D-01,0.1338172912386137D-01,
     9  0.1186698483870402D-01,0.1036829497712414D-01,
     *  0.8948692414831327D-02,0.7668331559589407D-02,
     1  0.6574914361419289D-02,0.5705863687114019D-02,
     2  0.5112956743920158D-02,0.4889038788564175D-02,
     3  0.5133491124385627D-02,0.5872410693663000D-02,
     4  0.7093845048344184D-02,0.8837823994856127D-02,
     5  0.1121155013926125D-01,0.1435878077387992D-01,
     6  0.1840442488792620D-01,0.2335095378124508D-01,
     7  0.2890102502062377D-01,0.3423877315696957D-01,
     8  0.3793713286816747D-01,0.3825098199606511D-01,
     9  0.3383783725961545D-01,0.2447649692113642D-01,
     *  0.1127575491491142D-01/
        
        do 1000 i=1,nqs(1)
        roots(i) =r1(i)
        weights(i) = w1(i)
 1000 continue

        nquads = nqs(1)

        return
        end
c
c
c
c
c
        subroutine lap_load_quads21(alpha,roots,weights,nquads)
        implicit real *8 (a-h,o-z)
        dimension r1(38),
     1          w1(38),nqs(1),
     2          roots(1),weights(1)

        data nqs/38/

        data r1/
     1  0.1043960810790154D-06,0.4400261102650955D-05,
     2  0.5443988618966938D-04,0.3451751310600191D-03,
     3  0.1405543111762432D-02,0.4163464341087692D-02,
     4  0.9719677325383295D-02,0.1891290780745890D-01,
     5  0.3196154401200659D-01,0.4838633555944449D-01,
     6  0.6719610994015095D-01,0.8718519644173405D-01,
     7  0.1072002284053786D+00,0.1263079471273086D+00,
     8  0.1438656313659517D+00,0.1595364714405959D+00,
     9  0.1732807968335555D+00,0.1852837648534299D+00,
     *  0.1957898068041391D+00,0.2049805163640499D+00,
     1  0.2129973183003286D+00,0.2200083557156852D+00,
     2  0.2262488045824353D+00,0.2320586627492171D+00,
     3  0.2379087341478420D+00,0.2443479910999556D+00,
     4  0.2519215177521056D+00,0.2611766719034450D+00,
     5  0.2727231683592060D+00,0.2872666694089364D+00,
     6  0.3055733661486898D+00,0.3283209453382436D+00,
     7  0.3557819903092600D+00,0.3873385215880370D+00,
     8  0.4209937894688825D+00,0.4532477173394547D+00,
     9  0.4796579126905007D+00,0.4959906863169810D+00/
     
        data w1/
     1  0.5019163697635540D-06,0.1311113600357899D-04,
     2  0.1159380239688559D-03,0.5520411735578650D-03,
     3  0.1727348182721874D-02,0.3980027762138994D-02,
     4  0.7278705415352023D-02,0.1114534447704862D-01,
     5  0.1486853075359356D-01,0.1781053805455794D-01,
     6  0.1960346196230490D-01,0.2017970168207869D-01,
     7  0.1969376459034200D-01,0.1841440648741393D-01,
     8  0.1664486081571053D-01,0.1469120722565020D-01,
     9  0.1283151714706458D-01,0.1121756191112134D-01,
     *  0.9823979457604733D-02,0.8579727090979646D-02,
     1  0.7481355622867309D-02,0.6579293068788580D-02,
     2  0.5957562162616426D-02,0.5741618378717081D-02,
     3  0.6051852929924815D-02,0.6917199430578531D-02,
     4  0.8319867328411758D-02,0.1029142801833873D-01,
     5  0.1291966223319855D-01,0.1629737016396146D-01,
     6  0.2043419009213892D-01,0.2511485023928850D-01,
     7  0.2971335523794538D-01,0.3306996395399429D-01,
     8  0.3364764633634036D-01,0.3009921881671630D-01,
     9  0.2199344482413428D-01,0.1019784589645511D-01/
        
        do 1000 i=1,nqs(1)
        roots(i) =r1(i)
        weights(i) = w1(i)
 1000 continue

        nquads = nqs(1)

        return
        end
c
c
c
c
c
        subroutine lap_load_quads22(alpha,roots,weights,nquads)
        implicit real *8 (a-h,o-z)
        dimension r1(37),
     1          w1(37),nqs(1),
     2          roots(1),weights(1)

        data nqs/37/

        data r1/
     1  0.6098817800920204D-07,0.2709468410267691D-05,
     2  0.3614885243594329D-04,0.2488399204406959D-03,
     3  0.1097156159320206D-02,0.3493024092547724D-02,
     4  0.8683012275719100D-02,0.1781771083816029D-01,
     5  0.3146421217678794D-01,0.4936708077245642D-01,
     6  0.7055349086602280D-01,0.9365746216825548D-01,
     7  0.1172737850799903D+00,0.1402154914414710D+00,
     8  0.1616392948680944D+00,0.1810602214991919D+00,
     9  0.1982929565501172D+00,0.2133582448322196D+00,
     *  0.2263953442551357D+00,0.2376080101614129D+00,
     1  0.2472406015190726D+00,0.2555709127134459D+00,
     2  0.2629233313855578D+00,0.2697135473882268D+00,
     3  0.2764900578300539D+00,0.2838837988809265D+00,
     4  0.2925087874621776D+00,0.3029534409127380D+00,
     5  0.3158295375535390D+00,0.3317774685143355D+00,
     6  0.3513716841747835D+00,0.3748846129228735D+00,
     7  0.4018961421988175D+00,0.4308508054886382D+00,
     8  0.4588479082504937D+00,0.4820037456475048D+00,
     *  0.4964414144446654D+00/
     
        data w1/
     1  0.2952834137892404D-06,0.8255714398426332D-05,
     2  0.7990889300047605D-04,0.4183730195948806D-03,
     3  0.1432104608690215D-02,0.3575304226037237D-02,
     4  0.7005417754405518D-02,0.1136269887160833D-01,
     5  0.1588560105421553D-01,0.1975208884992558D-01,
     6  0.2238686264897498D-01,0.2358252001364780D-01,
     7  0.2345027613564785D-01,0.2229244041028154D-01,
     8  0.2047540084719624D-01,0.1833723225064372D-01,
     9  0.1613323254465284D-01,0.1402143284053433D-01,
     *  0.1208707731202788D-01,0.1037911293547233D-01,
     1  0.8932198190255991D-02,0.7781694291720199D-02,
     2  0.6991405149214189D-02,0.6681362689550615D-02,
     3  0.6979033059257378D-02,0.7910716055608199D-02,
     4  0.9436238528536875D-02,0.1155484060400544D-01,
     5  0.1430555445298806D-01,0.1768991401649775D-01,
     6  0.2155063016090497D-01,0.2541331058454886D-01,
     7  0.2835498165467791D-01,0.2906508345621387D-01,
     8  0.2626428414437742D-01,0.1937908608413402D-01,
     *  0.9044030663138746D-02/
        
        do 1000 i=1,nqs(1)
        roots(i) =r1(i)
        weights(i) = w1(i)
 1000 continue

        nquads = nqs(1)

        return
        end
c
c
c
c
c
        subroutine lap_load_quads23(alpha,roots,weights,nquads)
        implicit real *8 (a-h,o-z)
        dimension r1(36),
     1          w1(36),nqs(1),
     2          roots(1),weights(1)

        data nqs/36/

        data r1/
     1  0.1176046587610346D-06,0.4917460928880011D-05,
     2  0.6082713995907158D-04,0.3878618337271292D-03,
     3  0.1595328637804210D-02,0.4788364181854442D-02,
     4  0.1134899579322471D-01,0.2243972178062018D-01,
     5  0.3853116407279466D-01,0.5922115319012456D-01,
     6  0.8338602773999432D-01,0.1095245447196912D+00,
     7  0.1361175785410382D+00,0.1618877287847576D+00,
     8  0.1859250931437484D+00,0.2077009700164394D+00,
     9  0.2270101149715495D+00,0.2438794634002787D+00,
     *  0.2584764228243324D+00,0.2710422475686059D+00,
     1  0.2818571444822588D+00,0.2912305716448977D+00,
     2  0.2995179418441881D+00,0.3071708833463489D+00,
     3  0.3147824294753620D+00,0.3230309753609709D+00,
     4  0.3325628984893009D+00,0.3439631228394052D+00,
     5  0.3577794672834437D+00,0.3744850793947704D+00,
     6  0.3943189390140416D+00,0.4169888540611141D+00,
     7  0.4413005430105731D+00,0.4649232314867261D+00,
     8  0.4845995999464879D+00,0.4969465221748939D+00/
     
        data w1/
     1  0.5639049517890573D-06,0.1463405538166760D-04,
     2  0.1297231087197830D-03,0.6234443881181062D-03,
     3  0.1979451823810818D-02,0.4646102319375982D-02,
     4  0.8677505345744786D-02,0.1358546144133847D-01,
     5  0.1852747288917272D-01,0.2265903552800764D-01,
     6  0.2541422993341017D-01,0.2660522284784607D-01,
     7  0.2636591666492649D-01,0.2502261097868373D-01,
     8  0.2296505220716742D-01,0.2055409386249866D-01,
     9  0.1807028257036745D-01,0.1569743938477480D-01,
     *  0.1353724237885937D-01,0.1164127147340691D-01,
     1  0.1003985874286159D-01,0.8765337906080239D-02,
     2  0.7883572774789843D-02,0.7522145064634150D-02,
     3  0.7816574758320405D-02,0.8788519293207530D-02,
     4  0.1037208555363895D-01,0.1252017451570035D-01,
     5  0.1519397562001779D-01,0.1826214115289065D-01,
     6  0.2136356883778827D-01,0.2378197560188908D-01,
     7  0.2444992168438854D-01,0.2223856724640489D-01,
     8  0.1652979166246989D-01,0.7755032478354977D-02/
        
        do 1000 i=1,nqs(1)
        roots(i) =r1(i)
        weights(i) = w1(i)
 1000 continue

        nquads = nqs(1)

        return
        end
c
c
c
c
c
        subroutine lap_load_quads24(alpha,roots,weights,nquads)
        implicit real *8 (a-h,o-z)
        dimension r1(35),
     1          w1(35),nqs(1),
     2          roots(1),weights(1)

        data nqs/35/

        data r1/
     1  0.1529906921185591D-06,0.6302001344919350D-05,
     2  0.7704989661697681D-04,0.4866061801942270D-03,
     3  0.1984655715765149D-02,0.5909897672626635D-02,
     4  0.1389773087270659D-01,0.2726180371819452D-01,
     5  0.4643802894226584D-01,0.7081246787800114D-01,
     6  0.9895275100353557D-01,0.1290461530149840D+00,
     7  0.1593255766082436D+00,0.1883557846549570D+00,
     8  0.2151541573018342D+00,0.2391827981620664D+00,
     9  0.2602678998298497D+00,0.2784951128582915D+00,
     *  0.2941126567779471D+00,0.3074562495700417D+00,
     1  0.3188998917482010D+00,0.3288381900872515D+00,
     2  0.3377111032061088D+00,0.3460656719094051D+00,
     3  0.3545783055055797D+00,0.3639528959909052D+00,
     4  0.3748030239980343D+00,0.3876270423372349D+00,
     5  0.4027797320354713D+00,0.4203349140337457D+00,
     6  0.4398215136620036D+00,0.4599186072794403D+00,
     7  0.4783350465650187D+00,0.4921677052380370D+00,
     *  0.4989887436503226D+00/
     
        data w1/
     1  0.7309363992458877D-06,0.1867212856098463D-04,
     2  0.1635140984539048D-03,0.7778344802392592D-03,
     3  0.2446462735835408D-02,0.5688020973324144D-02,
     4  0.1051672370402065D-01,0.1628693173486180D-01,
     5  0.2195877747487823D-01,0.2654350618709138D-01,
     6  0.2942753586147005D-01,0.3045993927385437D-01,
     7  0.2985779552806163D-01,0.2803912539635027D-01,
     8  0.2546943352353469D-01,0.2256070020533648D-01,
     9  0.1962542335215792D-01,0.1687162756326984D-01,
     *  0.1442024502110343D-01,0.1232960007114905D-01,
     1  0.1062305607745256D-01,0.9325662818174685D-02,
     2  0.8510591278088096D-02,0.8313369092582291D-02,
     3  0.8831577916217150D-02,0.1001988481021451D-01,
     4  0.1176364271121757D-01,0.1394577256668161D-01,
     5  0.1637648404276025D-01,0.1866100639943400D-01,
     6  0.2009486020762288D-01,0.1971190074391007D-01,
     7  0.1662037630261521D-01,0.1059950300156685D-01,
     *  0.3139711781509538D-02/
        
        do 1000 i=1,nqs(1)
        roots(i) =r1(i)
        weights(i) = w1(i)
 1000 continue

        nquads = nqs(1)

        return
        end
c
c
c
c
c
        subroutine lap_load_quads25(alpha,roots,weights,nquads)
        implicit real *8 (a-h,o-z)
        dimension r1(34),
     1          w1(34),nqs(1),
     2          roots(1),weights(1)

        data nqs/34/

        data r1/
     1  0.1225982139267846D-06,0.5131203126864318D-05,
     2  0.6434480125244145D-04,0.4190672778735162D-03,
     3  0.1765843537635791D-02,0.5427932174043136D-02,
     4  0.1314272369570479D-01,0.2645470904672805D-01,
     5  0.4606893724169075D-01,0.7155451610332802D-01,
     6  0.1015074549621373D+00,0.1340016992272090D+00,
     7  0.1670751693733100D+00,0.1990799964403394D+00,
     8  0.2288485215845148D+00,0.2557083493646308D+00,
     9  0.2794072250605848D+00,0.3000010233073754D+00,
     *  0.3177389867465534D+00,0.3329667326613546D+00,
     1  0.3460626473229100D+00,0.3574165367052993D+00,
     2  0.3674514947007251D+00,0.3766842310601269D+00,
     3  0.3857742157799117D+00,0.3954508926870007D+00,
     4  0.4063455101691772D+00,0.4188854844722816D+00,
     5  0.4332186655882935D+00,0.4490520457443817D+00,
     6  0.4654201538526556D+00,0.4805618884649894D+00,
     7  0.4922131951822407D+00,0.4986186788689408D+00/
     
        data w1/
     1  0.5869565622798008D-06,0.1531238311610165D-04,
     2  0.1383987150102813D-03,0.6833796294471367D-03,
     3  0.2234189443370620D-02,0.5390693263563815D-02,
     4  0.1030815809767531D-01,0.1643850761253202D-01,
     5  0.2271491235996275D-01,0.2801405592600432D-01,
     6  0.3156229469091131D-01,0.3309324970262028D-01,
     7  0.3277638201592804D-01,0.3103856991891127D-01,
     8  0.2838779759766542D-01,0.2529176054780605D-01,
     9  0.2211725953403961D-01,0.1911353330096657D-01,
     *  0.1642107966773386D-01,0.1409782864260302D-01,
     1  0.1215860635277912D-01,0.1061860472321877D-01,
     2  0.9536231382247905D-02,0.9040372624660043D-02,
     3  0.9264626679772845D-02,0.1019548097260535D-01,
     4  0.1166588850306174D-01,0.1344143027064669D-01,
     5  0.1518038511701060D-01,0.1632830892849851D-01,
     6  0.1611031305181477D-01,0.1377926095023274D-01,
     7  0.9200296612055107D-02,0.3642243824965730D-02/
        
        do 1000 i=1,nqs(1)
        roots(i) =r1(i)
        weights(i) = w1(i)
 1000 continue

        nquads = nqs(1)
        
        return
        end
c
c
c
c
c
        subroutine lap_load_quads26(alpha,roots,weights,nquads)
        implicit real *8 (a-h,o-z)
        dimension r1(33),
     1          w1(33),nqs(1),
     2          roots(1),weights(1)

        data nqs/33/

        data r1/
     1  0.1204665040217617D-06,0.5101025203294499D-05,
     2  0.6439083133344366D-04,0.4210935853971728D-03,
     3  0.1780510670108872D-02,0.5493864993527599D-02,
     4  0.1336172809850870D-01,0.2703006709241084D-01,
     5  0.4732035964355660D-01,0.7389343136188809D-01,
     6  0.1053795222818997D+00,0.1398181193744556D+00,
     7  0.1751548028147412D+00,0.2096156372767253D+00,
     8  0.2418994861391764D+00,0.2712132146096942D+00,
     9  0.2972076752610182D+00,0.3198700617948134D+00,
     *  0.3394115621145061D+00,0.3561718880631081D+00,
     1  0.3705499330602956D+00,0.3829643519683444D+00,
     2  0.3938489107666381D+00,0.4036870457503434D+00,
     3  0.4130600303900642D+00,0.4226145870476530D+00,
     4  0.4329041272307069D+00,0.4442324861304583D+00,
     5  0.4565600261188530D+00,0.4693991276165626D+00,
     6  0.4817118549431913D+00,0.4919512719332533D+00,
     *  0.4984002683441528D+00/
     
        data w1/
     1  0.5785756171192350D-06,0.1526649829265209D-04,
     2  0.1388190173225910D-03,0.6881739783769908D-03,
     3  0.2259022126916284D-02,0.5478826196268495D-02,
     4  0.1054350125737468D-01,0.1693600403355652D-01,
     5  0.2358541322121302D-01,0.2932443089653659D-01,
     6  0.3331227739780186D-01,0.3521615102398890D-01,
     7  0.3515843681985894D-01,0.3354594116532965D-01,
     8  0.3089077171149392D-01,0.2768038041975135D-01,
     9  0.2430832808595868D-01,0.2105325339867343D-01,
     *  0.1808752759507254D-01,0.1550027422712829D-01,
     1  0.1332568903684578D-01,0.1157483784768076D-01,
     2  0.1027378733447768D-01,0.9499656437570561D-02,
     3  0.9357250752142337D-02,0.9846863063583505D-02,
     4  0.1078410435591925D-01,0.1186778269386669D-01,
     5  0.1270811578259061D-01,0.1279686092322327D-01,
     6  0.1156514100187548D-01,0.8616411790991824D-02,
     *  0.4060121332699757D-02/
        
        do 1000 i=1,nqs(1)
        roots(i) =r1(i)
        weights(i) = w1(i)
 1000 continue

        nquads = nqs(1)

        return
        end
c
c
c
c
c
        subroutine lap_load_quads27(alpha,roots,weights,nquads)
        implicit real *8 (a-h,o-z)
        dimension r1(31),
     1          w1(31),nqs(1),
     2          roots(1),weights(1)

        data nqs/31/

        data r1/
     1  0.1657816593696194D-06,0.6887874631615285D-05,
     2  0.8515070241903721D-04,0.5447592072963921D-03,
     3  0.2253049120428790D-02,0.6804949156412190D-02,
     4  0.1622406469696855D-01,0.3223501310137396D-01,
     5  0.5554298020265824D-01,0.8554222242422277D-01,
     6  0.1205396767818191D+00,0.1582836515155795D+00,
     7  0.1965154057452154D+00,0.2333586969560955D+00,
     8  0.2674979453554888D+00,0.2981857879131508D+00,
     9  0.3251508686846253D+00,0.3484683236423274D+00,
     *  0.3684336618902056D+00,0.3854607276524971D+00,
     1  0.4000115094203033D+00,0.4125606475701516D+00,
     2  0.4235982260675119D+00,0.4336655011352936D+00,
     3  0.4433671349815885D+00,0.4532561836715803D+00,
     4  0.4636152076951332D+00,0.4742793052714006D+00,
     5  0.4845381246632039D+00,0.4931529748901603D+00,
     *  0.4986331574757227D+00/
     
        data w1/
     1  0.7932475183700843D-06,0.2047630777018788D-04,
     2  0.1817378815654529D-03,0.8781865881585523D-03,
     3  0.2809279884355610D-02,0.6645172531954270D-02,
     4  0.1249365441560205D-01,0.1965080014025149D-01,
     5  0.2685984093341972D-01,0.3284691145717391D-01,
     6  0.3676307192715458D-01,0.3834183453818464D-01,
     7  0.3780612080093171D-01,0.3566125628577299D-01,
     8  0.3249436994841713D-01,0.2883834320707313D-01,
     9  0.2510648913466006D-01,0.2157938821686525D-01,
     *  0.1842143588556440D-01,0.1571029542265349D-01,
     1  0.1347033407103790D-01,0.1170891651041003D-01,
     2  0.1045589469738049D-01,0.9781370867709765D-02,
     3  0.9718270418949713D-02,0.1010983946567710D-01,
     4  0.1058238225673951D-01,0.1062758824398102D-01,
     5  0.9679236161051696D-02,0.7291558438833911D-02,
     *  0.3465150113181864D-02/
        
        do 1000 i=1,nqs(1)
        roots(i) =r1(i)
        weights(i) = w1(i)
 1000 continue

        nquads = nqs(1)

        return
        end
c
c
c
c
c
        subroutine lap_load_quads28(alpha,roots,weights,nquads)
        implicit real *8 (a-h,o-z)
        dimension r1(30),
     1          w1(30),nqs(1),
     2          roots(1),weights(1)

        data nqs/30/

        data r1/
     1  0.1375569522112955D-06,0.5816689483332086D-05,
     2  0.7320839865008881D-04,0.4773243501185010D-03,
     3  0.2014140235508927D-02,0.6209707431330917D-02,
     4  0.1510520685759166D-01,0.3057709523636041D-01,
     5  0.5356880150643622D-01,0.8369485694876775D-01,
     6  0.1193822409466576D+00,0.1583752356692833D+00,
     7  0.1983108339904904D+00,0.2371536429564864D+00,
     8  0.2734189291412028D+00,0.3062123725487231D+00,
     9  0.3351543988088275D+00,0.3602541276931786D+00,
     *  0.3817776331227961D+00,0.4001343825471133D+00,
     1  0.4157923789842649D+00,0.4292272844291822D+00,
     2  0.4409080609497751D+00,0.4513118849427533D+00,
     3  0.4609288215172261D+00,0.4701714806231783D+00,
     4  0.4791653658054097D+00,0.4875597858809533D+00,
     5  0.4945095871588293D+00,0.4989054546081549D+00/
     
        data w1/
     1  0.6606564779147015D-06,0.1739412720996103D-04,
     2  0.1575947544408792D-03,0.7788134826662565D-03,
     3  0.2552576317665769D-02,0.6191157709313561D-02,
     4  0.1192667283516361D-01,0.1918133705183158D-01,
     5  0.2673494775439071D-01,0.3324579332223098D-01,
     6  0.3774291209296537D-01,0.3984184663473518D-01,
     7  0.3968649508522210D-01,0.3775118315980115D-01,
     8  0.3463175808993117D-01,0.3089407311238661D-01,
     9  0.2699374840848411D-01,0.2325196902640368D-01,
     *  0.1986471205407857D-01,0.1692764378320932D-01,
     1  0.1446785621697331D-01,0.1247978977886163D-01,
     2  0.1096097191515352D-01,0.9929706498941542D-02,
     3  0.9376885644412217D-02,0.9132883433895076D-02,
     4  0.8796889492024565D-02,0.7847369371409867D-02,
     5  0.5858763849813902D-02,0.2775594339905886D-02/
        
        do 1000 i=1,nqs(1)
        roots(i) =r1(i)
        weights(i) = w1(i)
 1000 continue

        nquads = nqs(1)

        return
        end
c
c
c
c
c
        subroutine lap_load_quads29(alpha,roots,weights,nquads)
        implicit real *8 (a-h,o-z)
        dimension r1(29),
     1          w1(29),nqs(1),
     2          roots(1),weights(1)

        data nqs/29/

        data r1/
     1  0.6885054137900548D-07,0.3468244640079113D-05,
     2  0.4982749514355659D-04,0.3595583106627351D-03,
     3  0.1638894687065347D-02,0.5354859288438784D-02,
     4  0.1360420911151459D-01,0.2845101123658775D-01,
     5  0.5108927727517708D-01,0.8134410363658252D-01,
     6  0.1177380988663090D+00,0.1579789767085625D+00,
     7  0.1995682209656578D+00,0.2402912320164668D+00,
     8  0.2784896150807827D+00,0.3131304178815728D+00,
     9  0.3437402507655530D+00,0.3702759881552716D+00,
     *  0.3929838359707643D+00,0.4122752752137357D+00,
     1  0.4286318920252608D+00,0.4425429214770100D+00,
     2  0.4544745064835552D+00,0.4648600326903417D+00,
     3  0.4740788913060850D+00,0.4823680319747383D+00,
     4  0.4896605975837615D+00,0.4954866143693595D+00,
     *  0.4991051644363292D+00/
     
        data w1/
     1  0.3446989564887765D-06,0.1089469218029237D-04,
     2  0.1130209652948381D-03,0.6184304366377276D-03,
     3  0.2186503629972796D-02,0.5606952369855386D-02,
     4  0.1125115779626816D-01,0.1865032355851628D-01,
     5  0.2659549400424259D-01,0.3365899123973139D-01,
     6  0.3873620949606920D-01,0.4132073738105927D-01,
     7  0.4148417592514879D-01,0.3968447518865935D-01,
     8  0.3654140971713039D-01,0.3266357640054430D-01,
     9  0.2855131524392122D-01,0.2456328626532188D-01,
     *  0.2092258321536701D-01,0.1774160020213789D-01,
     1  0.1505366476807633D-01,0.1284617455608203D-01,
     2  0.1108942153603530D-01,0.9746254526775566D-02,
     3  0.8734024781383638D-02,0.7835268988381664D-02,
     4  0.6669384144366133D-02,0.4852326630179801D-02,
     *  0.2271997641704290D-02/
        
        do 1000 i=1,nqs(1)
        roots(i) =r1(i)
        weights(i) = w1(i)
 1000 continue

        nquads = nqs(1)

        return
        end
c
c
c
c
c
        subroutine lap_load_quads30(alpha,roots,weights,nquads)
        implicit real *8 (a-h,o-z)
        dimension r1(28),
     1          w1(28),nqs(1),
     2          roots(1),weights(1)

        data nqs/28/

        data r1/
     1  0.3257475201557487D-08,0.1386702058029290D-05,
     2  0.3074804103949265D-04,0.2703153580752505D-03,
     3  0.1374821300353866D-02,0.4804272075892332D-02,
     4  0.1275158156157309D-01,0.2747009164166490D-01,
     5  0.5035098482232008D-01,0.8132755564473226D-01,
     6  0.1188981946665100D+00,0.1606385778192019D+00,
     7  0.2038717046542791D+00,0.2462157860454876D+00,
     8  0.2858890421400501D+00,0.3217858999971895D+00,
     9  0.3534023192243696D+00,0.3806921651185125D+00,
     *  0.4039132803438389D+00,0.4234950016012018D+00,
     1  0.4399387380919404D+00,0.4537509429904652D+00,
     2  0.4654004799079226D+00,0.4752852986164209D+00,
     3  0.4836824502802959D+00,0.4906591304413267D+00,
     4  0.4959895736721271D+00,0.4992119800365815D+00/
     
        data w1/
     1  0.4523918430002209D-07,0.5314857506551200D-05,
     2  0.7812794714658741D-04,0.5033200580292610D-03,
     3  0.1949590963469445D-02,0.5287015461982551D-02,
     4  0.1100016784092761D-01,0.1867938382689797D-01,
     5  0.2707047557027510D-01,0.3462295967876754D-01,
     6  0.4010239088099029D-01,0.4292280588833744D-01,
     7  0.4314147076711843D-01,0.4124866821616510D-01,
     8  0.3791503145557391D-01,0.3379721131234715D-01,
     9  0.2943057571233144D-01,0.2519432222111439D-01,
     *  0.2132111775722391D-01,0.1792699070663318D-01,
     1  0.1504556842712415D-01,0.1265747793852465D-01,
     2  0.1070867042105066D-01,0.9108702076124824D-02,
     3  0.7698584054478974D-02,0.6216459446623221D-02,
     4  0.4362573928970773D-02,0.2004977345080585D-02/
        
        do 1000 i=1,nqs(1)
        roots(i) =r1(i)
        weights(i) = w1(i)
 1000 continue

        nquads = nqs(1)

        return
        end
c
c
c
c
c
        subroutine lap_load_quads31(alpha,roots,weights,nquads)
        implicit real *8 (a-h,o-z)
        dimension r1(27),
     1          w1(27),nqs(1),
     2          roots(1),weights(1)

        data nqs/27/

        data r1/
     1  0.1988239476641103D-06,0.8195303461868343D-05,
     2  0.1008910267664202D-03,0.6441149196088527D-03,
     3  0.2660549844415306D-02,0.8024941416908848D-02,
     4  0.1909792865441594D-01,0.3785501727425781D-01,
     5  0.6504299181796167D-01,0.9986360514229324D-01,
     6  0.1402678997506566D+00,0.1835951918645484D+00,
     7  0.2272197631066485D+00,0.2689939102194574D+00,
     8  0.3074384891328340D+00,0.3417358298438240D+00,
     9  0.3716113905919033D+00,0.3971774201329703D+00,
     *  0.4187843985436285D+00,0.4369015074396040D+00,
     1  0.4520308505764572D+00,0.4646498261660985D+00,
     2  0.4751687507800711D+00,0.4838849837122026D+00,
     3  0.4909177143763191D+00,0.4961466405913467D+00,
     *  0.4992480285268402D+00/
     
        data w1/
     1  0.9492617787660564D-06,0.2431589802020993D-04,
     2  0.2150415534613991D-03,0.1037386147842996D-02,
     3  0.3314164423850269D-02,0.7824116241146436D-02,
     4  0.1466653608349520D-01,0.2297752028930293D-01,
     5  0.3125971438093731D-01,0.3802848575590898D-01,
     6  0.4232536666004704D-01,0.4388412735950596D-01,
     7  0.4300405096002137D-01,0.4029806045170908D-01,
     8  0.3645697049673411D-01,0.3209466644866918D-01,
     9  0.2767768610054627D-01,0.2351478721534445D-01,
     *  0.1977878179050915D-01,0.1654004943626825D-01,
     1  0.1379886354385787D-01,0.1150818631223828D-01,
     2  0.9580338763282255D-02,0.7873272644094467D-02,
     3  0.6171807655613819D-02,0.4228359694662475D-02,
     *  0.1916394431151483D-02/

ccc        call prin2('alpha = *',alpha,1)

        do 1000 i=1,nqs(1)
        roots(i) =r1(i)
        weights(i) = w1(i)
 1000 continue

        nquads = nqs(1)

        return
        end
c
c
c
c
c
        subroutine lap_load_quads32(alpha,roots,weights,nquads)
        implicit real *8 (a-h,o-z)
        dimension r1(17),r2(18),r3(19),r4(21),r5(22),r6(24),
     1          w1(17),w2(18),w3(19),w4(21),w5(22),w6(24),nqs(6),
     2          roots(1),weights(1)

        data nqs/17,18,19,21,22,24/

        data r1/
     1  0.8189795324279940D-06,0.2495034098781455D-04,
     2  0.2467080310117823D-03,0.1333518871395668D-02,
     3  0.4844652975941494D-02,0.1324901870546511D-01,
     4  0.2933572415314467D-01,0.5533495342505474D-01,
     5  0.9226113700576302D-01,0.1396923800507081D+00,
     6  0.1958392306854226D+00,0.2576732768559124D+00,
     7  0.3210759492974225D+00,0.3811662127250982D+00,
     8  0.4329039750446157D+00,0.4718123303963446D+00,
     9  0.4945561023623646D+00/
     
        data w1/
     1  0.3614506754719551D-05,0.6770288566370669D-04,
     2  0.4790855770474080D-03,0.1953389905092994D-02,
     3  0.5491147963060689D-02,0.1179855335400285D-01,
     4  0.2075801814326680D-01,0.3141269956712233D-01,
     5  0.4236753123818572D-01,0.5219446737486123D-01,
     6  0.5959253120541354D-01,0.6338400802542333D-01,
     7  0.6259731807222043D-01,0.5673038587109879D-01,
     8  0.4598911785853625D-01,0.3126628088707459D-01,
     9  0.1391414756517462D-01/
     
     
        data r2/
     1  0.3891386226503557D-06,0.1557535483437794D-04,
     2  0.1872518973021824D-03,0.1159023345846999D-02,
     3  0.4610725491031949D-02,0.1334144427977381D-01,
     4  0.3040156332354844D-01,0.5768965114617797D-01,
     5  0.9507754745747577D-01,0.1406606050043246D+00,
     6  0.1917202226807272D+00,0.2455886856377994D+00,
     7  0.2999783126889956D+00,0.3528048754377573D+00,
     8  0.4017625159040585D+00,0.4440013169080214D+00,
     9  0.4762379885079912D+00,0.4953793422882565D+00/
     
        data w2/
     1  0.1840318847947200D-05,0.4586776059335129D-04,
     2  0.3940719974420386D-03,0.1826418489207106D-02,
     3  0.5559087282833656D-02,0.1243724591291399D-01,
     4  0.2201570110217245D-01,0.3253360371194191D-01,
     5  0.4190639339734209D-01,0.4879117908681829D-01,
     6  0.5287815905997144D-01,0.5448107745309342D-01,
     7  0.5395875275329237D-01,0.5131634169102707D-01,
     8  0.4612236267571803D-01,0.3779687494702966D-01,
     9  0.2614529785264657D-01,0.1178972450710861D-01/
     
     
        data r3/
     1  0.1402265463376357D-05,0.4723257355443762D-04,
     2  0.4817095287345028D-03,0.2586095776872020D-02,
     3  0.9117927038684004D-02,0.2378635399632536D-01,
     4  0.4947919558699672D-01,0.8647353673512377D-01,
     5  0.1321282879746890D+00,0.1822315954015305D+00,
     6  0.2328006293642555D+00,0.2811373127291608D+00,
     7  0.3259114097107341D+00,0.3667105889117933D+00,
     8  0.4034764996698214D+00,0.4359900423985806D+00,
     9  0.4634684300273203D+00,0.4844185163494968D+00,
     *  0.4969555227402398D+00/
     
        data w3/
     1  0.6402718345050133D-05,0.1309178334042694D-03,
     2  0.9395839028963587D-03,0.3736720249636860D-02,
     3  0.9978202441824002D-02,0.1987417592979565D-01,
     4  0.3156322813049448D-01,0.4195921284094597D-01,
     5  0.4861501608694098D-01,0.5090987627166046D-01,
     6  0.4977852023639329D-01,0.4668170950116090D-01,
     7  0.4280745677400670D-01,0.3878981784390334D-01,
     8  0.3471254596484818D-01,0.3018714181496717D-01,
     9  0.2451563204862972D-01,0.1705599003893589D-01,
     *  0.7757849371210723D-02/
     
     
        data r4/
     1  0.1647399925891663D-05,0.5497925761449626D-04,
     2  0.5536785988376915D-03,0.2935964163110965D-02,
     3  0.1024742403771503D-01,0.2654051747863092D-01,
     4  0.5494083767717407D-01,0.9563633149178491D-01,
     5  0.1453587072063393D+00,0.1987507614695764D+00,
     6  0.2505869343638130D+00,0.2974079453459422D+00,
     7  0.3378704029125830D+00,0.3721457728724760D+00,
     8  0.4011091317814074D+00,0.4257617771689390D+00,
     9  0.4469200424935372D+00,0.4650400869946909D+00,
     *  0.4800897973727621D+00,0.4915000220144611D+00,
     *  0.4983363836252091D+00/
     
        data w4/
     1  0.7512412149773471D-05,0.1517463324294542D-03,
     2  0.1072790464677998D-02,0.4208674788542379D-02,
     3  0.1112544787605294D-01,0.2201966380334552D-01,
     4  0.3482197004419813D-01,0.4599854301215081D-01,
     5  0.5251237309690235D-01,0.5337177301716301D-01,
     6  0.4971932796715488D-01,0.4371091731428094D-01,
     7  0.3725716407252070D-01,0.3144898617585484D-01,
     8  0.2664941071552552D-01,0.2279461847721739D-01,
     9  0.1960014783346610D-01,0.1663579799360251D-01,
     *  0.1336625580639031D-01,0.9289885650337377D-02,
     *  0.4236993146037066D-02/
    

        data r5/
     1  0.4740612523461186D-05,0.1329222647844369D-03,
     2  0.1162080336859950D-02,0.5487506243677310D-02,
     3  0.1739491087857232D-01,0.4154351107698845D-01,
     4  0.8026391425578440D-01,0.1316942883142226D+00,
     5  0.1901995941218938D+00,0.2487081453283713D+00,
     6  0.3014001141424566D+00,0.3452765559570492D+00,
     7  0.3800862176144144D+00,0.4071934806958686D+00,
     8  0.4284095237589926D+00,0.4453427642998798D+00,
     9  0.4592080856198088D+00,0.4708356757310370D+00,
     *  0.4807022055556212D+00,0.4889220387238194D+00,
     1  0.4952281915680793D+00,0.4990597734039171D+00/
     
        data w5/
     1  0.2070260085722365D-04,0.3460782929966694D-03,
     2  0.2103324041832648D-02,0.7288692220295358D-02,
     3  0.1735354863351635D-01,0.3136717067328795D-01,
     4  0.4575509164756819D-01,0.5612320876085601D-01,
     5  0.5966014788001369D-01,0.5635949328820882D-01,
     6  0.4852753464635365D-01,0.3920525254006406D-01,
     7  0.3066229318973997D-01,0.2386454022132943D-01,
     8  0.1883807086042689D-01,0.1523032764377324D-01,
     9  0.1263802582171818D-01,0.1069681076140714D-01,
     *  0.9055922827425304D-02,0.7338198229589671D-02,
     1  0.5175167226955832D-02,0.2390397991783699D-02/
     
     
        data r6/
     1  0.5725664834327429D-05,0.1585549341215552D-03,
     2  0.1363964638033536D-02,0.6336843117420739D-02,
     3  0.1980116779574722D-01,0.4673999177354833D-01,
     4  0.8947732527321974D-01,0.1457555766638765D+00,
     5  0.2092335875678320D+00,0.2719423084501454D+00,
     6  0.3271553187337524D+00,0.3712987137983211D+00,
     7  0.4041984406176724D+00,0.4278611764187079D+00,
     8  0.4448700565991547D+00,0.4574047294976355D+00,
     9  0.4670019589134379D+00,0.4746689006336106D+00,
     *  0.4810484796044318D+00,0.4865291491924725D+00,
     1  0.4912848207251137D+00,0.4952626385524154D+00,
     2  0.4981804583148322D+00,0.4997043203141944D+00/
     
        data w6/
     1  0.2495230786177181D-04,0.4104347223951812D-03,
     2  0.2445966593398190D-02,0.8318420559109963D-02,
     3  0.1949561945361371D-01,0.3480534752096664D-01,
     4  0.5028315510076251D-01,0.6117128311154185D-01,
     5  0.6440532650365916D-01,0.5985747950128832D-01,
     6  0.4996910212563039D-01,0.3831208237255233D-01,
     7  0.2785142133314935D-01,0.1992087439883614D-01,
     8  0.1446294939683498D-01,0.1085991998481181D-01,
     9  0.8499439776359525D-02,0.6939759136495816D-02,
     *  0.5883511384420381D-02,0.5107134862096627D-02,
     1  0.4395813981164109D-02,0.3510663350828704D-02,
     2  0.2258857958157812D-02,0.8104845640647460D-03/
c        
c         %%%%%%%%%%%%%%%%%%%%%



ccc        call prin2('alpha = *',alpha,1)
        dtwo = 2
        if (alpha .ge. 1/dtwo) then
        do 1000 i=1,nqs(1)
        roots(i) =r1(i)
        weights(i) = w1(i)
 1000 continue

        nquads = nqs(1)
        endif

        if (alpha .ge. dtwo**(-2)) then
        if (alpha .lt. dtwo**(-1)) then

        do 2000 i=1,nqs(2)
        roots(i) =r2(i)
        weights(i) = w2(i)
 2000 continue

        nquads = nqs(2)

        endif
        endif


        if (alpha .ge. dtwo**(-3)) then
        if (alpha .lt. dtwo**(-2)) then

        do 3000 i=1,nqs(3)
        roots(i) =r3(i)
        weights(i) = w3(i)
 3000 continue

        nquads = nqs(3)

        endif
        endif


        if (alpha .ge. dtwo**(-4)) then
        if (alpha .lt. dtwo**(-3)) then

        do 4000 i=1,nqs(4)
        roots(i) =r4(i)
        weights(i) = w4(i)
 4000 continue

        nquads = nqs(4)

        endif
        endif


        if (alpha .ge. dtwo**(-5)) then
        if (alpha .lt. dtwo**(-4)) then

        do 5000 i=1,nqs(5)
        roots(i) =r5(i)
        weights(i) = w5(i)
 5000 continue

        nquads = nqs(5)

        endif
        endif


        if (alpha .ge. dtwo**(-6)) then
        if (alpha .lt. dtwo**(-5)) then

        do 6000 i=1,nqs(6)
        roots(i) =r6(i)
        weights(i) = w6(i)
 6000 continue

        nquads = nqs(6)

        endif
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
        dimension r1(16),r2(18),r3(19),r4(20),r5(22),r6(23),
     1          w1(16),w2(18),w3(19),w4(20),w5(22),w6(23),nqs(6),
     2          roots(1),weights(1)

        data nqs/16,18,19,20,22,23/

        data r1/
     1  0.1544485658250285D-05,0.4488235390273785D-04,
     2  0.4218614492226822D-03,0.2164442736857530D-02,
     3  0.7469229243835309D-02,0.1945506446434533D-01,
     4  0.4119394758440851D-01,0.7464527220086736D-01,
     5  0.1200814476798851D+00,0.1760445612922790D+00,
     6  0.2395167327604341D+00,0.3061018657018772D+00,
     7  0.3703489653397081D+00,0.4264179120192769D+00,
     8  0.4689792191181763D+00,0.4939975820077109D+00/
     
        data w1/
     1  0.6747583969802286D-05,0.1195755816494533D-03,
     2  0.7973866194660992D-03,0.3057701018050175D-02,
     3  0.8089748139129156D-02,0.1641331451492855D-01,
     4  0.2740056359233415D-01,0.3954722556595329D-01,
     5  0.5107870966828019D-01,0.6034195189554069D-01,
     6  0.6586501094429211D-01,0.6638536150941619D-01,
     7  0.6111943027490623D-01,0.5011665151950217D-01,
     8  0.3432569866936899D-01,0.1533492290321277D-01/
     
     
        data r2/
     1  0.1495778769726298D-06,0.9001664438166460D-05,
     2  0.1306330757818964D-03,0.9003379748746093D-03,
     3  0.3840323375117847D-02,0.1166927154285213D-01,
     4  0.2756066912009599D-01,0.5371380297547435D-01,
     5  0.9030773481716716D-01,0.1355948918472505D+00,
     6  0.1868361720982411D+00,0.2412382813695678D+00,
     7  0.2963688812762562D+00,0.3500291101616886D+00,
     8  0.3998328792851005D+00,0.4428592037860798D+00,
     9  0.4757383755834960D+00,0.4952803617163517D+00/
     
        data w2/
     1  0.8031332392793145D-06,0.2883569723602466D-04,
     2  0.2927936751724413D-03,0.1494876110134413D-02,
     3  0.4851214503602691D-02,0.1136254816840876D-01,
     4  0.2080641716134139D-01,0.3152711425476171D-01,
     5  0.4135146406473976D-01,0.4874943307813905D-01,
     6  0.5326075271681575D-01,0.5514084048094027D-01,
     7  0.5476179279162170D-01,0.5216776102347603D-01,
     8  0.4695184853680810D-01,0.3852795069016708D-01,
     9  0.2668242428562237D-01,0.1204112962777317D-01/
     
     
        data r3/
     1  0.1251921498728987D-05,0.4312212914653661D-04,
     2  0.4475063920576115D-03,0.2435842158299089D-02,
     3  0.8685091401409833D-02,0.2287274901772383D-01,
     4  0.4797340527878985D-01,0.8445644805355412D-01,
     5  0.1298704437330058D+00,0.1800784604795876D+00,
     6  0.2310377659232045D+00,0.2799191459064549D+00,
     7  0.3252595346931680D+00,0.3665433051552861D+00,
     8  0.4036426017041649D+00,0.4363048183269062D+00,
     9  0.4637598153927526D+00,0.4845817671696417D+00,
     *  0.4969921300259089D+00/
     
        data w3/
     1  0.5748151168743391D-05,0.1203699247764588D-03,
     2  0.8797645415258223D-03,0.3549702879495336D-02,
     3  0.9593376079613957D-02,0.1931208676820788D-01,
     4  0.3097506950920358D-01,0.4155854084776281D-01,
     5  0.4854695592223939D-01,0.5117606925689601D-01,
     6  0.5026974568528850D-01,0.4725727003989769D-01,
     7  0.4334726278510207D-01,0.3920759741876214D-01,
     8  0.3495576021046162D-01,0.3024317769303875D-01,
     9  0.2442425312522484D-01,0.1690973981495444D-01,
     *  0.7667509346379970D-02/
     
     
        data r4/
     1  0.2396808744155500D-05,0.7586554340524748D-04,
     2  0.7350567746801663D-03,0.3776675204815248D-02,
     3  0.1281473064832015D-01,0.3230925786612758D-01,
     4  0.6514868815560120D-01,0.1105579412315238D+00,
     5  0.1640890725609744D+00,0.2196643950691566D+00,
     6  0.2720745925497708D+00,0.3183712075398600D+00,
     7  0.3577915303509737D+00,0.3908992045280614D+00,
     8  0.4187343421837728D+00,0.4422809825582075D+00,
     9  0.4621767825464026D+00,0.4785318387234534D+00,
     *  0.4908525009692616D+00,0.4982111994667882D+00/
     
        data w4/
     1  0.1077689118922597D-04,0.2060430234197731D-03,
     2  0.1399036657280767D-02,0.5302605169145572D-02,
     3  0.1356225763927704D-01,0.2594815272211722D-01,
     4  0.3960144193286306D-01,0.5043707937578688D-01,
     5  0.5556635968482009D-01,0.5469210435408690D-01,
     6  0.4964024044916539D-01,0.4284234668948034D-01,
     7  0.3610703958130801D-01,0.3029020288551736D-01,
     8  0.2554853204235907D-01,0.2165402422471887D-01,
     9  0.1815759963223517D-01,0.1446781896039877D-01,
     *  0.1000967823749725D-01,0.4556659847333244D-02/
     
     
        data r5/
     1  0.4600505718186255D-05,0.1291227641858519D-03,
     2  0.1130324124370463D-02,0.5345808356168057D-02,
     3  0.1697691032794878D-01,0.4063575514281529D-01,
     4  0.7872476790964139D-01,0.1295959599930984D+00,
     5  0.1878934834164793D+00,0.2467456195762233D+00,
     6  0.3003269527160944D+00,0.3454207385679524D+00,
     7  0.3814823637175387D+00,0.4096637098840626D+00,
     8  0.4316908983542703D+00,0.4491669963088722D+00,
     9  0.4633190533532988D+00,0.4749668492463271D+00,
     *  0.4845370983025195D+00,0.4920668908379259D+00,
     1  0.4972491532436520D+00,0.4996780429423127D+00/
     
        data w5/
     1  0.2009456035419539D-04,0.3363353692242274D-03,
     2  0.2047452490960842D-02,0.7109457065880726D-02,
     3  0.1697074806517478D-01,0.3078192983747638D-01,
     4  0.4511456087892152D-01,0.5569453074018246D-01,
     5  0.5970897444273454D-01,0.5699712968714038D-01,
     6  0.4962967642865127D-01,0.4048589733130080D-01,
     7  0.3184948405818057D-01,0.2481237727377888D-01,
     8  0.1951274162127875D-01,0.1564405077765671D-01,
     9  0.1279633127102319D-01,0.1057045469861079D-01,
     *  0.8575812205998415D-02,0.6429313471062607D-02,
     1  0.3854772751375509D-02,0.1057874973032457D-02/
     
     
        data r6/
     1  0.5966944885794545D-05,0.1644374942885832D-03,
     2  0.1407773444153891D-02,0.6510331548321358D-02,
     3  0.2025801200946931D-01,0.4764431719577157D-01,
     4  0.9093486740548463D-01,0.1477872139547291D+00,
     5  0.2118115480712318D+00,0.2750513557195279D+00,
     6  0.3308170863517420D+00,0.3755319773278415D+00,
     7  0.4089617755104104D+00,0.4330467461101423D+00,
     8  0.4503435913425476D+00,0.4630380947233654D+00,
     9  0.4726771577431689D+00,0.4802632372490062D+00,
     *  0.4864096948345098D+00,0.4914417964906121D+00,
     1  0.4954357528635326D+00,0.4982524981148072D+00,
     *  0.4997080124013129D+00/
     
        data w6/
     1  0.2597514358307715D-04,0.4248633246001313D-03,
     2  0.2517851541088801D-02,0.8517304057896361D-02,
     3  0.1986554267183284D-01,0.3531983151928869D-01,
     4  0.5085975562964311D-01,0.6173428149366038D-01,
     5  0.6493776274561703D-01,0.6039464794798102D-01,
     6  0.5053728638368303D-01,0.3887574684954097D-01,
     7  0.2833591684497053D-01,0.2027670371206434D-01,
     8  0.1468462714429811D-01,0.1095988944601046D-01,
     9  0.8482497352792509D-02,0.6790203807696801D-02,
     *  0.5556393967082444D-02,0.4520268303961989D-02,
     1  0.3442245716383337D-02,0.2152141929906659D-02,
     *  0.7882624664173875D-03/
c        
c         %%%%%%%%%%%%%%%%%%%%%

ccc        call prin2('alpha = *',alpha,1)
        dtwo = 2
        if (alpha .ge. 1/dtwo) then
        do 1000 i=1,nqs(1)
        roots(i) =r1(i)
        weights(i) = w1(i)
 1000 continue

        nquads = nqs(1)
        endif

        if (alpha .ge. dtwo**(-2)) then
        if (alpha .lt. dtwo**(-1)) then

        do 2000 i=1,nqs(2)
        roots(i) =r2(i)
        weights(i) = w2(i)
 2000 continue

        nquads = nqs(2)

        endif
        endif


        if (alpha .ge. dtwo**(-3)) then
        if (alpha .lt. dtwo**(-2)) then

        do 3000 i=1,nqs(3)
        roots(i) =r3(i)
        weights(i) = w3(i)
 3000 continue

        nquads = nqs(3)

        endif
        endif


        if (alpha .ge. dtwo**(-4)) then
        if (alpha .lt. dtwo**(-3)) then

        do 4000 i=1,nqs(4)
        roots(i) =r4(i)
        weights(i) = w4(i)
 4000 continue

        nquads = nqs(4)

        endif
        endif


        if (alpha .ge. dtwo**(-5)) then
        if (alpha .lt. dtwo**(-4)) then

        do 5000 i=1,nqs(5)
        roots(i) =r5(i)
        weights(i) = w5(i)
 5000 continue

        nquads = nqs(5)

        endif
        endif


        if (alpha .ge. dtwo**(-6)) then
        if (alpha .lt. dtwo**(-5)) then

        do 6000 i=1,nqs(6)
        roots(i) =r6(i)
        weights(i) = w6(i)
 6000 continue

        nquads = nqs(6)

        endif
        endif

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
        dimension cx(1),cy(1),csx(92),csy(92)

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

