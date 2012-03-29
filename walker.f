      program walker
ccc 
c     A split of rwalk.f: This creates a random walk
c     Inputs:
c      1   # of edges (# of points = # of edges + 1
c      2   # of closures (This gets carried over into knotter.f)
c     Outputs
c      walker.out  Contains Inputs and list of points in walk 
c
c     Random numbers are done using time() as a seed
c
c     Initialized variables are kept
c
c     Many unused statements removed for readability
c
c
c      parameter (repete = 10000)
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
     :       omega,         theta
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
     6       center(1,3),      ss(6) 
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
     :         n,            repete
      character(len=15)
     : filename
c
c
c
c     write(6,*) 'number of edges - integer number'
      read(5,*) number
c
c     Update the number as most statements use the number of points
c     not the number of edges
c
      number = number+1
c
c     write(6,*) 'number of points(# edges+1)= ', number
c
c
c     write(6,*) 'how many closures?'
      read(5,*) repete
c
c     Uses UNIX time as a random seed
c     Initialize random number sequence
      istef1 = irand(time())
c
c     write(6,*) 'random seed = ',istef1
c
c
      open(9,file='walker.out', status='unknown',form ='formatted')
c
      write(9,*) number
      write(9,*) repete
c
c - generation of random walk  -
c
c     Begin at origin
c
       points(1,1) = 0
       points(1,2) = 0
       points(1,3) = 0
c
c     Ben's distribution of points
c
       do 1995 i = 2, number
       rr = rand(0)
       tt = rand(0)
       theta = 2*pi*rr
       omega = 2*tt - 1
c       
       points(i,1) = points(i-1,1) + sqrt(1-omega*omega)*dcos(theta)
       points(i,2) = points(i-1,2) + sqrt(1-omega*omega)*dsin(theta)
       points(i,3) = points(i-1,3) + omega
c
 1995  continue
c
      do 1997 i = 1, number
      write(9,533) points(i,1), points(i,2), points(i,3)
 1997 continue
 533  format(3(2x,f21.16))
c
      close(9)
c

      stop
      end
