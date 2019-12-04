      subroutine freddethelm(rzk,izk,coefs,par2,fdet)
      implicit real *8 (a-h,o-z)
      real *8 rzk,coefs(*),par2(*)
      complex *16 fdet,zk
c
c
c      temporary variables
c
      real *8, allocatable :: xys(:,:),dxys(:,:),pl(:),pr(:),qwts(:)
      real *8, allocatable :: xs(:),ys(:)
      integer, allocatable :: el(:),er(:),imid(:)
      integer, allocatable :: iscorn(:),ixys(:),ks(:),nepts(:)
      integer, allocatable :: lns(:),rns(:)
      integer, allocatable :: ichl(:),ichr(:)

      real *8, allocatable :: vtmp(:,:),rpan(:),rmid(:),angs(:)

      integer, allocatable :: icl(:),icr(:)
      real *8, allocatable :: xsgnl(:),xsgnr(:),verts(:,:)
      integer, allocatable :: icsgnl(:),icsgnr(:)
      
      complex *16, allocatable :: xmatc(:,:,:),xmat(:,:)
      integer, allocatable :: ipiv(:)

      real *8, allocatable :: tsc(:),wtsc(:),umatc(:,:),vmatc(:,:)
      complex *16 ima,imainv4,h0,h1,z
      data imainv4/(0.0d0,0.25d0)/
      data ima/(0.0d0,1.0d0)/

      done = 1
      pi = atan(done)*4

      nr0 = 36


      ncorner = 36
      allocate(tsc(ncorner),wtsc(ncorner),umatc(ncorner,ncorner),
     1   vmatc(ncorner,ncorner))


      call getcornstruct(ncorner,tsc,wtsc,umatc,vmatc)


      zk = rzk
      ifexpon = 1

      nverts = par2(1)
      ifdn = par2(2)
      lcoefs = par2(3)
      a0 = par2(4)
      b0 = par2(5)

      allocate(verts(2,nverts))
      do i=1,nverts
        verts(1,i) = par2(5+i)
        verts(2,i) = par2(5+nverts+i)
      enddo



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
        rr = pi/b0
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
      kmid = 24



      do i=1,nverts
        rmid(i) = rpan(i) -pl(i)-pr(i)
        hmid = 2.0d0*min(pl(i),pr(i))
        hmida = 2.0d0*min(pl(i)*abs(tan(angs(i))),
     1      pr(i)*abs(tan(angs(i+1))))
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

cc      call prinf('iscorn=*',iscorn,nch)


c
c
c       set up corner interactions
c
      ncint = nverts
      allocate(icl(ncint),icr(ncint),icsgnl(ncint),icsgnr(ncint))
      allocate(xsgnl(ncint),xsgnr(ncint),ichl(ncint),ichr(ncint))

      do i=1,ncint
        icl(i) = i-1
        icr(i) = i
        icsgnl(i) = 1
        icsgnr(i) = 0
        xsgnl(i) = -1.0d0
        xsgnr(i) = -1.0d0
      enddo
      icl(1) = nedges

      do ich=1,nch
        if(iscorn(ich).gt.0) then
          ichr(iscorn(ich)) = ich
        endif

        if(iscorn(ich).lt.0) then
          ichl(-iscorn(ich)) = ich
        endif
      enddo

c
c
c
c
      allocate(xmatc(ncorner,ncorner,ncint))


      allocate(xmat(npts,npts))
      allocate(ipiv(npts))


      do i=1,ncint
        call helmcornmat(angs(i),zk,pl(i),tsc,wtsc,umatc,vmatc,
     1     ncorner,coefs,lcoefs,xmatc(1,1,i))
      enddo



c
c        generate the matrix
c
      do i=1,npts
        do j=1,npts
          xmat(j,i) = 0
        enddo
      enddo




      istart = 0
      do iedge=1,nedges
        do ii=1,nepts(iedge)
          i = istart+ii
          dst = sqrt(dxys(1,i)**2 + dxys(2,i)**2)
          rnx = dxys(2,i)
          rny = -dxys(1,i)

          jstart = 0
          do jedge=1,nedges
            if(jedge.eq.iedge) goto 1000
            do jj=1,nepts(jedge)
              j = jj+jstart
              rrr = sqrt((xys(1,i)-xys(1,j))**2 + 
     1          (xys(2,i)-xys(2,j))**2)
              z = rrr*zk
              call hank103(z,h0,h1,ifexpon)
              rrn = ((xys(1,j)-xys(1,i))*rnx+
     1            (xys(2,j)-xys(2,i))*rny)/dst

              xmat(j,i) = imainv4*h1*rrn/rrr*
     1            zk*sqrt(qwts(j)*qwts(i))
            enddo
 1000       continue              
            jstart = jstart + nepts(jedge)
          enddo
        enddo
        istart = istart + nepts(iedge)
      enddo

c
c       fix correction matrices
c
      do icint=1,ncint
        if(icsgnl(icint).eq.0) ipts1 = lns(icl(icint))
        if(icsgnl(icint).eq.1) ipts1 = rns(icl(icint))
        
        if(icsgnr(icint).eq.0) ipts2 = lns(icr(icint))
        if(icsgnr(icint).eq.1) ipts2 = rns(icr(icint))

        do i=1,ncorner
          do j=1,ncorner
            xmat(j+ipts1-1,i+ipts2-1) = xmatc(j,i,icint)*xsgnl(icint)
            xmat(j+ipts2-1,i+ipts1-1) = xmatc(j,i,icint)*xsgnr(icint)
          enddo
        enddo
      enddo

      do i=1,npts
        do j=1,npts
          xmat(j,i) = -2.0d0*(-1)**(ifdn)*xmat(j,i)
        enddo
      enddo

      do i=1,npts
        xmat(i,i) = 1.0d0
      enddo



      info = 0

      call zgetrf(npts,npts,xmat,npts,ipiv,info)

      fdet = 1.0d0


      do i=1,npts
        fdet = fdet*xmat(i,i)
        if(ipiv(i).ne.i) fdet = -fdet
      enddo

      return
      end
