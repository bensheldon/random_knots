      program knotter
ccc 
c     for random walks and knotting
c     revised number length of walk+1 =Number of points
c
c     parameter (repete = 10000)
      parameter (sts = 4)
      parameter (pi = 3.14159265358979)
ccc
c
      real*8
     :       X,             Y,              Z,
     1       rdomang,       diff,           cdiam,
     2       newcurv,       nlcurv,         x1,
     3       y1,            z1,             zj,
     4       t1,            tj,             sig1,
     5       sgg1,          qq1,            qq2,
     6       pp1,           pp2,            raport,
     7       subangle,      xxc,            sdiam,
     8       diam,          elength,        avcross,
     9       totalcross,    tdiam,          avdiam,
     :       crgm,          trgm,           avrgm,
     :       xxy,           xxz,            cx,
     :       cy,            cz,             xxx,
     :       ndiam,         xj,             yj,
     :       radius,        sradius,        rope,
     :       ssss,          rrss,           length,
     :       tt,            rr,             
     :       omega,         theta,
     :       getx,          gety,           getz,
     :       balldiam
c
      real*8
     :       points(501,3),  ronecords(3),  movepts(501,3),
     2       setpts(501,3),  rrrr,          matricks(3,3),
     3       rotpts(501,3),  newpts(501,3),  newvect(501,3),
     4       newc(501),                     pts(501,3),
     5       d_pts(501,3),   r_pts(501,3),   r_delta(501,3),
     6       t(501),         rotp(501,3),    rotd(501,3),
     7       p(501,6),       sols(2),
     8       solutions(2001,2), coords(3),   coordinates(2001,3),
     :       orso(2001,2),    ordsols(2001,2,2001), sutions(20001,2)
c
      real*8
     :       tvalues(2001),    svalues(2001),  point(3),    
     3       xx(1001),         yy(1001),         zz(1001), 
     4       reverse(10001,2), p_app(10001,6),   tintix(10001),
     5       p_end(10001,6),   inter(10001),     th(10001,6),
     6       center(1,3),      ss(6),            
c
      character*1
     :       ilet(10001),      olet(10001),
     :       cdwr(10001,9,4),  lett(10001,24)
c
      character*3
     : ch1, ch2, ch4, ch5, ch6
c
      character*4
     :       cd(10001,9), minus
c
      integer  rone,         rtwo,           r,
     :         counter,      tag(501),       iseed,
     :         r_tag(501),   stp,            failure,
     :         rotag(501),   sn,             snn,
     :         ith(10001,4),  inlet(10001),    onlet(10001),
     :         num(10001),    walk,           number,
     :         irstrt,       istef1,         maxcross,
     :         n,            type
c
c      Set Maximum Number of Crossings
c
      maxcross = 999
c
c
      open(8,file='walker.out',status='unknown',form='formatted')
c
c
c     Read number of points in walk
      read(8,*) number
      n = number
      ndiam = number - 1
c     
c     Read  Number of closures
c$$$      read(8,*) repete
c$$$      write(6,*) repete
c
c     Read in points
      do 1982, idm = 1,number
         read(8,533) getx, gety, getz
         points(idm,1)=getx
         points(idm,2)=gety
         Points(idm,3)=getz
 1982    continue
c
         idm = 0
c
c     Random Seed
c     Uses UNIX time as a random seed
c     Adds one in case walker.f runs in less than 1 second
      istef1 = irand(time()+1)
c
c
c     File to write knot data to
      open(7,file='knotter.out', status='unknown',form ='formatted')
c     File to write closure points to
      open(10,file='closures.out', status='unknown',form ='formatted')
c     
c      
c     
c$$$      open(11,file = 'center.type', status='unknown', form ='formatted') 
c$$$      read(11,*) type
c      

c     Read in Center
      open(12,file = 'center.out', status='unknown', form='formatted')
c$$$      if(type.eq.0) then
c$$$c        Origin Centered
c$$$         center(1,1) = 0
c$$$         center(1,2) = 0
c$$$         center(1,3) = 0
c$$$         else
c$$$            if (type.eq.1) then
c$$$c              Endpoint average
c$$$               center(1,1) = points(number,1)*.5
c$$$               center(1,2) = points(number,2)*.5
c$$$               center(1,3) = points(number,3)*.5
c$$$               else
c$$$c              Smallball 
c$$$                  read(12,*) center(1,1)
c$$$                  read(12,*) center(1,2)
c$$$                  read(12,*) center(1,3)
c$$$             endif
c$$$       endif
c$$$       close(11)

       read(12,*) center(1,1)
       read(12,*) center(1,2)
       read(12,*) center(1,3)
       read(12,*) balldiam

       close(12)
c
c
c
c
c       Calculate  Diameter
      write(6,*) 'calculate diameter'
      cdiam = 0
       do 560 idm = 1, number
        do 562 jdm = idm + 1, number
        diam = 0
         do 563 kdm = 1, 3
         xxc = points(idm,kdm) - points(jdm,kdm)
         diam = diam + xxc*xxc
  563    continue
       diam = dsqrt(diam)
c      write(6,*) 'diam = ',diam
        if(diam.gt.cdiam) then
        cdiam = diam
        endif
  562   continue
  560  continue
c      write(6,*) ' Diameter = ',cdiam
c        End of Diameter Calculation
c
c       End Point Distance
c 
      diam = dsqrt(points(n,1)*points(n,1) + points(n,2)*
     :points(n,2) + points(n,3)*points(n,3))
c      write(6,*) 'End Distance = ',diam
c
c
      counter = 0
      tdiam = 0
      avdiam = 0
      totalcross = 0
      avcross = 0
      ircont = 1
      intix = 0 
      jrcont = 1
      walk = 1
      trgm = 0
      avrgm = 0
      newcurv = 100
c
      do 1001 kcm = 1, walk
c
      do 1011 ish = 1, repete
      intix = intix + 1
c
 192  continue
c 
c -   set to zero .........
c   
      lp = 1
      lrev = 1
      do 561 i = 1, 3001
      sutions(i,1) = 0.0
      reverse(i,1) = 0.0
      sutions(i,2) = 0.0
      reverse(i,2) = 0.0
c	
      p_app(i,1) = 0.0
      p_app(i,2) = 0.0
      p_app(i,3) = 0.0
      p_app(i,4) = 0.0
      p_app(i,5) = 0.0
      p_app(i,6) = 0.0
c
	
      p_end(i,1) = 0.0
      p_end(i,2) = 0.0
      p_end(i,3) = 0.0
      p_end(i,4) = 0.0
      p_end(i,5) = 0.0
      p_end(i,6) = 0.0
c
      inter(i) = 0.0
      num(i) = 0
      tintix(i) = 0.0
      inlet(i)  = 0
      onlet(i) = 0
      ilet(i) = ' '
      olet(i) = ' '
      ith(i,1) = 0
      th(i,1) = 0.0
      ith(i,2) = 0
      th(i,2) = 0.0
      ith(i,3) = 0
      th(i,3) = 0.0
      ith(i,4) = 0
      th(i,4) = 0.0
  561 continue
      do 567 j = 1, 9 
      do 569 i = 1, 3001 
      cd(i,j) = ' '
      cdwr(i,j,1) = ' '
      cdwr(i,j,2) = ' '
      cdwr(i,j,3) = ' '
      cdwr(i,j,4) = ' '
  569 continue
  567 continue
      do 572 j = 1, 24 
      do 574 i = 1, 3001 
      lett(i,j) = ' '
  574 continue
  572 continue
c
      X = 0.0
      Y = 0.0
      Z = 0.0
      do 1501 i = 1, number+1
      do 1203 j = 1,3
      movepts(i,j) = 0.0
 1203 continue
 1501 continue
      do 1207 i = 1,3
      ronecords(i) = 0.0
 1207 continue
      do 1209 i = 1, number+1
      xx(i) = 0.0
      yy(i) = 0.0
      zz(i) = 0.0
      r_tag(i) = 0
      tag(i) =  0
      rotag(i) = 0
      setpts(i,1) = 0.0
      setpts(i,2) = 0.0
      setpts(i,3) = 0.0
 1209 continue
c
c
c
c
      seed = iseed
c
c
c     Random Closure of Walk
c     Millett's fixed distribution
c
c
 1832  rr = rand(0)
       tt = rand(0) 
        if(0.lt.tt.and.tt.lt.dcos(-pi/2 +rr*pi)) then
        tt = tt/dcos(-pi/2 +rr*pi)
        points(number+1,1) = center(1,1) + 
     :  .5*balldiam*dcos(-pi/2 + rr*pi)*dcos(tt*2*pi)
        points(number+1,2) = center(1,2) +
     :  .5*balldiam*dcos(-pi/2 + rr*pi)*dsin(tt*2*pi)
        points(number+1,3) = center(1,3) + 
     :  .5*balldiam*dsin(-pi/2 + rr*pi)
       else
        go to 1832
       endif
c     Output the closures to file
       write(10,533) points(number+1,1), points(number+1,2), 
     :  points(number+1,3)
c
c
c
c	initialization
c
      do 1215 i = 1,3
      do 1217 j = 1, number+1
      newvect(j,i) = 0.0
      newpts(j,i) = points(j,i)
      pts(j,i) = 0.0
      d_pts(j,i) = 0.0
      r_pts(j,i) = 0.0
      r_delta(j,i) = 0.0
      rotp(j,i) = 0.0
      rotd(j,i) = 0.0
 1217 continue
 1215 continue
c
c
c
C	we now have the new knot
C
c
c        Calculate Radius of Gryation
       cx = 0
       cy = 0
       cz = 0
        do 570 irm = 1, number + 1
	    cx = cx + newpts(irm,1)
            cy = cy + newpts(irm,2)
            cz = cz + newpts(irm,3)
  570   continue
	 cx = cx/(number + 1)
         cy = cy/(number + 1)
         cz = cz/(number + 1)
	 crgm = 0
       do 580 idm = 1, number + 1
 	   xxx = newpts(idm,1) - cx
 	   xxy = newpts(idm,2) - cy
	   xxz = newpts(idm,3) - cz
	   crgm = crgm + dsqrt(xxx*xxx + xxy*xxy + xxz*xxz)
  580    continue
c      write(6,*) 'diam = ',diam
	 rgm = crgm/(number + 1)
       trgm = trgm + rgm
       counter = counter + 1
c	   End of Radius of Gyration Calculation
c
c
c533   format(f9.6,x,f9.6,x,f10.6)
 533   format(3(2x,f21.16))
 534   format(a10)
 535   format(3(2x,f21.17))
 536   format((2x,f1.6))
c535   format(f9.6,x,f9.6,x,f9.6)
c
c1035 continue
c
 3779 format(' # ',i6,'  rgm =  ',f21.17,'  diam = ',f21.17)
c
c
c
c --- first while -------------------------------------------------------|
c                                                                        |
c     if(newcurv.le.pi/10000) go to 999 
c
      do 35 i = 1, number+1
      pts(i,1) = newpts(i,1)
      pts(i,2) = newpts(i,2)
      pts(i,3) = newpts(i,3)
  35  continue
c
      nn = number+1
c
c     write(6,*) counter, '"next case = "',kcm,'"."',ish,
c    :              '"pts = "'
c     write(6,*) 'number = ',number+1
c      do 973 i = 1, number+1
c      write(6,533) pts(i,1), pts(i,2), pts(i,3)
c973  continue
c     write(6,*) ' d_pts = '
c
      do 37 j=1,(nn-1)
      d_pts(j,1) = pts(j+1,1) - pts(j,1)
      d_pts(j,2) = pts(j+1,2) - pts(j,2)
      d_pts(j,3) = pts(j+1,3) - pts(j,3)
c     write(6,*) d_pts(j,1)
c     write(6,*) d_pts(j,2)
c     write(6,*) d_pts(j,3)
  37  continue
c
      d_pts(nn,1) = pts(1,1) - pts(nn,1)
      d_pts(nn,2) = pts(1,2) - pts(nn,2)
      d_pts(nn,3) = pts(1,3) - pts(nn,3)
c     write(6,*) d_pts(nn,1)
c     write(6,*) d_pts(nn,2)
c     write(6,*) d_pts(nn,3)
 
c
c -tag------------------
      do 36 j=1,nn
      tag(j) = j
  36  continue
c - end tag
c
c - rotate ---------------------
c
      do 739 j=1,(nn-1)
      r_pts(j,1) = pts(j+1,1)
      r_pts(j,2) = pts(j+1,2)
      r_pts(j,3) = pts(j+1,3)
      r_delta(j,1) = d_pts(j+1,1)
      r_delta(j,2) = d_pts(j+1,2)
      r_delta(j,3) = d_pts(j+1,3)
 739  continue
      r_pts(nn,1) = pts(1,1)
      r_pts(nn,2) = pts(1,2)
      r_pts(nn,3) = pts(1,3)
      r_delta(nn,1) = d_pts(1,1)
      r_delta(nn,2) = d_pts(1,2)
      r_delta(nn,3) = d_pts(1,3)
c
c - rotate tag
      do 738 i = 1, (nn-1)
      r_tag(i) = tag(i+1)
c      write(6,*) r_tag(i)
 738  continue
      r_tag(nn) = tag(1)
c
c - end rotate ------------------
c
c ------ second while ------------------------------------------|
c                                                               |
      ind = 0
      stp = nn
c
c - set to zero ---------------------
      do 517 i = 1, number+1 
      do 519 k = 1, number+1
      ordsols(i,1,k) = 0.0
      ordsols(i,2,k) = 0.0
  519 continue
  517 continue
      snn  = 0
      sn = 0
c
c - set to zero ---------------------
c
  704 continue
c
      if( ind.ge.stp) go to 703
      ind = ind + 1
c - rotate ---------------------
      do 39 j=1,nn
      rotp(j,1) = r_pts(j,1)
      rotp(j,2) = r_pts(j,2)
      rotp(j,3) = r_pts(j,3)
      rotd(j,1) = r_delta(j,1)
      rotd(j,2) = r_delta(j,2)
      rotd(j,3) = r_delta(j,3)
  39  continue
c
c - rotate tag
      do 38 i = 1, nn
      rotag(i) = r_tag(i)
  38  continue
c
c - end rotate ------------------
c
c - initialize p, solutions, coordinates
c
      npts = 0
      do 7501 i = 1,number+1
      p(i,1) = 0.0
      p(i,2) = 0.0
      p(i,3) = 0.0
      p(i,4) = 0.0
      p(i,5) = 0.0
      p(i,6) = 0.0
 7501 continue
      do 503 i = 1,number+1
      solutions(i,1) = 0.0
      solutions(i,2) = 0.0 
      orso(i,1) = 0.0
      orso(i,2) = 0.0
 503  continue
      do 505 i = 1,number+1
      coordinates(i,1) = 0.0
      coordinates(i,2) = 0.0
      coordinates(i,3) = 0.0
 505  continue
      do 507 i =1,number+1
      tvalues(i) = 0.0
      svalues(i) = 0.0
      t(i) = 0.0
 507  continue
      t1 = 0.0
      tj = 0.0
      x1 = 0.0
      y1 = 0.0
      z1 = 0.0
      zj = 0.0
      qq1 = 0.0
      pp1 = 0.0
      pp2 = 0.0
      do 511 i = 1,6
      ss(i) = 0.0
  511 continue
      do 513 i = 1,2
      sols(i) = 0.0
  513 continue
      do 515 i = 1,3
      coords(i) = 0.0
  515 continue
c
c - end of initialize: p, solutions, coordinates
c
      do 41 k = 3, (nn-1)
      diff = r_delta(k,1)*r_delta(1,2) - r_delta(1,1)*  
     1       r_delta(k,2)
c     write(6,*) 'diff = ',diff
c    
      if(diff.eq.0.0) go to 40 
      t(r_tag(1)) = ((r_pts(k,2) - r_pts(1,2))*r_delta(k,1) +
     :             (r_pts(1,1) - r_pts(k,1))*r_delta(k,2) )/
     :             diff
      t(r_tag(k)) = ((r_pts(k,2) - r_pts(1,2))*r_delta(1,1) +
     :             (r_pts(1,1) - r_pts(k,1))*r_delta(1,2) )/
     :             diff
c
      t1 = t(r_tag(1))
      tj = t(r_tag(k))
      z1 = r_pts(1,3) + r_delta(1,3) * t1
      x1 = r_pts(1,1) + r_delta(1,1) * t1
      y1 = r_pts(1,2) + r_delta(1,2) * t1
      zj = r_pts(k,3) + r_delta(k,3) * tj
      xj = r_pts(k,1) + r_delta(k,1) * tj
      yj = r_pts(k,2) + r_delta(k,2) * tj
c     write(6,*) 't1 = ',t1,' tj = ',tj
c     write(6,*) 'x1 = ',x1,' xj = ',xj
c     write(6,*) 'y1 = ',y1,' yj = ',yj
      sig1 = z1 - zj
c     write(6,*) ' sig1 = ', sig1
c     write(6,*) ' at 801 lp = ',lp
c
      if( sig1.lt.0) then 
        sn = -1
        go to 801
        endif
      if( sig1.eq.0) then 
        sn = 0
        go to 801
        endif
      if( sig1.gt.0) then 
        sn = 1 
        go to 801
        endif
c
 801  continue
c
      qq1 = r_delta(1,1)
      qq2 = r_delta(1,2)
      pp1 = r_delta(k,1)
      pp2 = r_delta(k,2)
c
      sgg1 = qq1*pp2 - qq2*pp1
      if( sgg1.lt.0) then 
        snn = -1
        go to 802
        endif
      if( sgg1.eq.0) then 
        snn = 0
        go to 802
        endif
      if( sgg1.gt.0) then  
        snn = 1
        go to 802
        endif
c
 802  continue
c
       if (0.0.lt.t1.and.1.0.gt.t1) then
c
          if(0.0.lt.tj.and.1.0.gt.tj) then
c
              if(sn.eq.0) go to 40
              failure = 0
c
              npts = npts + 1
c
              ss(1) = 0.0
              ss(2) = snn
              ss(3) = sn
c
              ss(4) = x1
              ss(5) = y1
              ss(6) = t1
c
              do 51 i = 1,6
              p(npts,i) = ss(i)
  51          continue
c
c     write(6,*) 'ind =',ind,' npts =',npts,' p =',(p(npts,i), i =1,6) 
c
c
              sols(1) = t1
              sols(2) = tj
c
              do 53 i =1,2
              solutions(npts,i) = sols(i)
  53          continue
c
           do 63 i = 1, npts
           orso(i,1) = solutions(i,1)
           orso(i,2) = solutions(i,2)
  63       continue
c
              coords(1) = x1
              coords(2) = y1
              coords(3) = z1
c
              do 55 i=1,3
              coordinates(npts,i) = coords(i)
  55          continue
c
          else
              junk = 0
          endif
        else
          junk = 0
        endif
  40  continue
c
      sn = 1
      failure = 2
c
      if(t1.eq.1.0) go to 42 
         failure = 0
c
      if(t1.eq.0.0) go to 42 
         failure = 0
c
      if(tj.eq.1.0) go to 42 
         failure = 0
c
      if(tj.eq.0.0) go to 42 
         failure = 0
c
      if(sn.eq.0) go to 42 
         failure = 0
c
 41   continue
c
 42   continue
c
       if(npts.ge.199) then
         write(6,*) ' ........ attention .....'
         write(6,*) 'npts = ', npts
         write(6,*) ' ........ attention .....'
       endif
c
c      write(6,*) '551 npts = ',npts
c      write(6,*) '551 lp = ',lp
c
              do 553 j = 1, npts
              do 551 i = 1,6
              p_app(lp,i) = p(j,i)
 551    continue
c       write(6,*) 'lp =',lp,' npts=',npts,' p_app =',(p(j,i), i =1,6)
              lp = lp + 1
 553    continue
c      
c       write(6,*) 'control of lp =',lp
         if(npts.gt.0) then
c             write(6,*) 'reduction of lp'
              lp = lp - 1
         endif
        if(npts.gt.1) then      
              do 1553 j = 1, npts      
              inter(j) = p_app(lp-npts+j,6)
 1553         continue
c       write(6,*) ' lp = ',lp
c
c       write(6,*) ' npts = ',npts,'  inter = ',(inter(j), j=1,npts)
c
              call ssort(npts,inter)
c
              do 1555  j = 1, npts 
                 do 1557 k = 1, npts 
                    if(inter(j).eq.p_app(lp-npts+k,6)) then
                       do 1559 m = 1, 6
                          p_end(lp-npts+j,m) = p_app(lp-npts+k,m)
 1559                  continue
                    endif
 1557             continue
 1555          continue
c
        elseif (npts.eq.1) then
              do 1563 i = 1,6
              p_end(lp,i) = p_app(lp,i)
 1563         continue
c       write(6,*) ' lp = ',lp
c
c       write(6,*) ' npts = ',npts,'p_end  = ',(p_end(lp,j), j=1,6)

        else
        endif
      continue
c       
        if(npts.gt.0) then
         lp = lp + 1
        endif
c 
      do 61 i=1,npts
      tvalues(i) = p(i,6)
  61  continue
c
      do 62 i=1,npts
      svalues(i) = tvalues(i)
  62  continue 
      if(npts.gt.1) then
c
      call ssort(npts,svalues)
c
      call ssort(npts, orso)
c
      endif
      do 65 i = 1, npts
         do 67 j = 1,npts
           if (orso(i,1).eq.solutions(j,1))  then
               orso(i,2) = solutions(j,2)
           endif
  67     continue
  65  continue
c
      do 69 k = 1, npts
         ordsols(ind,1,k) = orso(k,1)
         ordsols(ind,2,k) = orso(k,2)
 69   continue
c
      if(ind.lt.stp) then
c
c - rotate ---------------------
c
      do 839 j=1,(nn-1)
      r_pts(j,1) = rotp(j+1,1)
      r_pts(j,2) = rotp(j+1,2)
      r_pts(j,3) = rotp(j+1,3)
      r_delta(j,1) = rotd(j+1,1)
      r_delta(j,2) = rotd(j+1,2)
      r_delta(j,3) = rotd(j+1,3)
 839  continue
      r_pts(nn,1) = rotp(1,1)
      r_pts(nn,2) = rotp(1,2)
      r_pts(nn,3) = rotp(1,3)
      r_delta(nn,1) = rotd(1,1)
      r_delta(nn,2) = rotd(1,2)
      r_delta(nn,3) = rotd(1,3)
c - rotate tag
      do 838 i = 1, (nn-1)
      r_tag(i) = rotag(i+1)
c      write(6,*) r_tag(i)
 838  continue
      r_tag(nn) = rotag(1)
c - end rotate ------------------
      go to 704 
      endif
c
 703  continue
c
      if(failure.ne.0) go to 999
      junk = 0
c
c     means not in Gen Posn
c
      l=1
      do 71 i = 1, number+1 
         do 73 k = 1, number  + 1
            if( ordsols(i,1,k).eq.0.0) then
                 go to 71
             else 
                sutions(l,1) = ordsols(i,1,k)
                sutions(l,2) = ordsols(i,2,k)
                l = l + 1
            endif
 73      continue
 71   continue
      lmax = l - 1
c     write(6,*) ' lmax = ', lmax
c
      if (lmax.ge.2*maxcross) then
       totalcross = totalcross + lmax/2
       go to 999
c        write(6,*) ' ATTENTION '
c        write(6,*) ' lmax greater that 1001 - increase memory - ', lmax
c        write(6,*) ' ATTENTION '
      endif
c
      failure = 1
      if(lmax.eq.0) then
       write(7,3779) counter, rgm, cdiam
       write(7,534) '1+1b1a1d1c'
       write(7,1795)
        if(cdiam.lt.ndiam) then
         ndiam = cdiam
        else
        go to 999
        endif
c
       go to 999
      endif
      failure = 0
      if(lmax.le.4) then
       write(7,3779) counter, rgm, cdiam 
       write(7,534) '1+2c2b1d1c'
       write(7,534) '2-2d1b1a2a'
       write(7,1795)
        if(cdiam.lt.ndiam) then
         ndiam = cdiam
        else
        go to 999
        endif
       go to 999
      endif
c
c
c     means has 0 crossings ------------
c      write(6,*) ' .......... sutions ............. '
c      do 75 l = 1, lmax
c      write(6,711)  sutions(l,1), sutions(l,2) 
c 75   continue 
  711  format(2x,2(f10.6,2x))
c
c
c
      do 77 k = 1, lmax
      reverse(k,1) = sutions(k,2)
      reverse(k,2) = sutions(k,1)
 77   continue
c      write(6,*) '.......... reverse ...............'
c      do 79 l = 1, lmax
c      write(6,711) reverse(l,1), reverse(l,2)
c 79   continue
c
c      write(6,*) ' .......p_end ................'
c      do 81 l  = 1,  lmax
c      write(6,811) ( p_end(l,i), i = 1,6 )
c 81   continue
  811  format(6(2x,f10.6))
c
      do 83 k = 1, lmax
         do 85 j = 1, lmax
            if(sutions(j,1).eq.reverse(k,1)) then
              tintix(lrev) = j 
              lrev =  lrev +  1
c             write(6,*) 'j= ',j,' k= ',k,' lrev= ',lrev
            endif
  85     continue
  83  continue
c
c      write(6,*) ' .........  t(lrev) .....'
c      write(6,*) (tintix(i), i = 1, lmax)
c
      do 91 i = 1, lmax
      th(i,1) = p_end(i,2)
      th(i,2) = p_end(i,3)
      th(i,3) = i
      th(i,4) = tintix(i)
  91  continue
      do 95 i = 1, lmax
      ith(i,1) = idint(th(i,1))
      ith(i,2) = idint(th(i,2))
      ith(i,3) = idint(th(i,3))
      ith(i,4) = idint(th(i,4))
  95  continue
c     write(6,*) ' ... Thistlethwaite code ...'
c     do 97 i = 1, lmax
c     write(6,923) (ith(i,j), j = 1, 4)
c 97  continue
 923  format(4(2x,i6))
c     write(6,*) ' '
c
c -  calculate the standard (Millett/Ewing convention)
c    matrix of the link, starting with the Thistlethwaite-
c    Representation in ith() -
c
      do 101 k =1, lmax
      if(mod(ith(k,3),2).eq.0) then
         num(k) = ith(k,3)/2
      else
         num(k) = ith(k,4)/2
      endif
 101  continue
      do 103 k =1, lmax
      if(ith(k,2).eq.1) then
         ilet(k) = 'c'
         inlet(k) = 5
         olet(k) = 'a'
         onlet(k) = 3
      else
         if(ith(k,1).eq.1) then
            ilet(k) = 'd'
            inlet(k) = 6
            olet(k) = 'b'
            onlet(k) = 4
         else
            ilet(k) = 'b'
            inlet(k) = 4
            olet(k) = 'd'
            onlet(k) = 6
         endif
      endif
 103  continue
c - start of the while5 loop --------------|
c
         write(ch2, '(i3)' ) num(1)
         write(ch4, '(i3)' ) num(lmax)
c
      do 501 i = 1, lmax
c
      if(ith(i,2).eq.1) then
         j = num(i)
c - stef
         write(ch1, '(i3)' ) j
c
       if(i.lt.lmax) then
         write(ch5, '(i3)' ) num(i+1)
       else
         write(ch5, '(i3)' ) num(1)
       endif
c
       if(i.gt.1) then
         write(ch6, '(i3)' ) num(i-1)
       else
         write(ch6, '(i3)' ) num(lmax)
       endif 
c - stef
c
         cd(j,1) = ch1 
         iminus = ith(i,1)*ith(i,2)
         write(minus, '(i3)') iminus
         cd(j,2) = minus 
c
           if(cd(j,2)(1:1).eq.'-') then
              go to 1071
              elseif (cd(j,2)(2:2).eq.'-') then
              go to 1071
              elseif (cd(j,2)(3:3).eq.'-') then
              go to 1071
              elseif (cd(j,2)(4:4).eq.'-') then
              go to 1071
           else
              cd(j,2)='+'
              go to 1072
           endif
c
 1071           cd(j,2)='-'
c
 1072     continue
c
         if(i.eq.lmax) then
            cd(j,3) = ch2 // ilet(1)
         else 
            cd(j,3) = ch5 // ilet(i+1)
c
         endif
         if(i.eq.1) then
            cd(j,5) = ch4 // olet(lmax)
         else
            cd(j,5) = ch6 // olet(i-1)
         endif
         if(i.eq.lmax) then
            cd(num(1),inlet(1)) = ch1 // 'a'
         else
            cd(num(i+1),inlet(i+1)) = ch1 // 'a'
         endif
         if(i.eq.1) then
            cd(num(lmax),onlet(lmax)) = ch1 // 'c'
         else
            cd(num(i-1),onlet(i-1)) = ch1 // 'c'
         endif
       else 
         if((ith(i,1)*ith(i,2)).eq.1) then
            j = num(i)
c - stef
         write(ch1, '(i3)' ) j
c
       if(i.lt.lmax) then
         write(ch5, '(i3)' ) num(i+1)
       else
         write(ch5, '(i3)' ) num(1)
       endif
c
       if(i.gt.1) then
         write(ch6, '(i3)' ) num(i-1)
       else
         write(ch6, '(i3)' ) num(lmax)
       endif
c - stef
            if(i.eq.lmax) then
              cd(j,6) = ch2 // ilet(1)
            else
              cd(j,6) = ch5 // ilet(i+1)
            endif
            if(i.eq.1) then
              cd(j,4) = ch4 // olet(lmax)
            else
              cd(j,4) = ch6 // olet(i-1)
            endif
            if(i.eq.lmax) then
              cd(num(1), inlet(1)) = ch1 // 'd'
            else
              cd(num(i+1), inlet(i+1)) = ch1 // 'd'
            endif
            if(i.eq.1) then
              cd(num(lmax),onlet(lmax)) = ch1 // 'b'
            else
              cd(num(i-1),onlet(i-1)) = ch1 // 'b'
            endif
         else 
            j = num(i)
c - stef
         write(ch1, '(i3)' ) j
c
       if(i.lt.lmax) then
         write(ch5, '(i3)' ) num(i+1)
       else
         write(ch5, '(i3)' ) num(1)
       endif
c
       if(i.gt.1) then
         write(ch6, '(i3)' ) num(i-1)
       else
         write(ch6, '(i3)' ) num(lmax)
       endif
c - stef
            if(i.eq.lmax) then
              cd(j,4) = ch2 // ilet(1)
            else
              cd(j,4) = ch5 // ilet(i+1)
            endif
            if(i.eq.1) then
              cd(j,6) = ch4 // olet(lmax)
            else
              cd(j,6) = ch6 // olet(i-1)
            endif
            if(i.eq.lmax) then
              cd(num(1), inlet(1)) = ch1 // 'b'
            else
              cd(num(i+1), inlet(i+1)) = ch1 // 'b'
            endif
            if(i.eq.1) then
              cd(num(lmax),onlet(lmax)) = ch1 // 'd'
            else
              cd(num(i-1),onlet(i-1)) = ch1 // 'd'
            endif
         endif 
       endif
 501   continue
c  
c - end of while#5 loop
c
      if(failure.ne.0) then
        go to 999
      else
        junk = 0
      endif
c
c - the next section prints out the matrix and calculates
c - the writhe
c - start of while#6 loop
c
       do 773 ist = 1,lmax/2
       do 777 jst = 1,6
       cdwr(ist,jst,1) = cd(ist,jst)(1:1)
       cdwr(ist,jst,2) = cd(ist,jst)(2:2) 
       cdwr(ist,jst,3) = cd(ist,jst)(3:3)
       cdwr(ist,jst,4) = cd(ist,jst)(4:4)
 777   continue
 773   continue
c
        if(cdiam.lt.ndiam) then
         ndiam = cdiam
        else
        go to 575
        endif
c
 3778 format('"#"',i6,'"case = "',i6,'"-"',i6,'"diam =  "',f21.17)
 575  write(7,3779) counter, rgm, cdiam 
c
      totalcross = totalcross + lmax/2
       do 775 ist = 1, lmax/2
          ia = 1
       do 779 jst = 1, 6
       do 781 kst = 1, 4 
       if (cdwr(ist,jst,kst).ne.' ') then
           lett(ist,ia) = cdwr(ist,jst,kst)
           ia = ia + 1
       endif
 781   continue
 779   continue
       ib = ia - 1
       write(7,1793) (lett(ist,lst), lst=1, ib)
 1793  format(24a1)
 775   continue	
c
       write(7,1795)
 1795  format(a1)
c
       if(failure.eq.0) then
          failure = 5 
       else
          junk = 0
       endif
c
c - next section is for when the knot has no crossings,
c - or is not in general position
c
      if(failure.lt.5) then
c
         if(failure.eq.1) then
            write(6,*) 'NO CROSSINGS'
         else
            junk = 0
         endif
c
         if(failure.eq.2) then
            write(6,*) 'Not in General Position'
          else
            junk = 0
         endif
c
         if(failure.eq.3) then
            write(6,*) 'The same point was chosen twice.
     : Execution stopped'
         else
            junk = 0
         endif
c
      endif
c
 999  continue
 1011 continue
 1001 continue
      avcross = totalcross/counter
      avdiam = tdiam/counter
      avrgm = trgm/counter
c      write(6,*) 'average number of crossings is  ',avcross
c      write(6,*) 'average diameter is ',avdiam
c      write(6,*) 'average radius of gyration is ',avrgm
c      write(6,*) 'closing minimal diameter is ',ndiam
      close(7)
      close(8)
      close(9)
      close(10)
      stop
      end
      subroutine ssort(n,ra)
      real*8 ra(n), rra 
      l=n/2+1
      ir=n
10    continue
        if(l.gt.1)then
          l=l-1
          rra=ra(l)
        else
          rra=ra(ir)
          ra(ir)=ra(1)
          ir=ir-1
          if(ir.eq.1)then
            ra(1)=rra
            return
          endif
        endif
        i=l
        j=l+l
20      if(j.le.ir)then
          if(j.lt.ir)then
            if(ra(j).lt.ra(j+1))j=j+1
          endif
          if(rra.lt.ra(j))then
            ra(i)=ra(j)
            i=j
            j=j+j
          else
            j=ir+1
          endif
        go to 20
        endif
        ra(i)=rra
      go to 10
      end
