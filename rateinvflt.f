c---- rateinvflt.f for public distribution
c copyright Art Frankel  2025
c--- finds rates of multiple segment ruptures
c--- given slip rates and segment jumping probabilities
c--- for single fault
c--- written By Art Frankel April-Sept 2025.
c  solving for rate for nucleation on each segment j that ruptures through segment i
c--- data are slip rates on segment i
c--- uses slip tapers
c  complile with gfortran -ffixed-line-length-128 -o rateinv rateinv.f nnls.f
c--- make sure dimension statements in nnls.f are consitent with rateinvflt
c---- uses non-negative least squares subroutine nnls.f from Hanson and Lawson
c-----(1974). Solving Least Squares Problems, Prentice Hall
c---  c(mm)= jump probability of segment boundary mm
c---- find probs of all possible scenarios
c--- prob(i,j) = prob of scenario j that starts at segment i 
c--- oiutput files:   rateinvflt.out, fort.2, fort.3
      dimension a(40,20),atmp(40,20),a2(40,20),c(20),xlen0(20),pred(20)
      dimension xlen(20,20),prob(20,500),sliprate(20),slip(500)
      dimension b(40),x(40),w(40),zz(40),index(40),xmag(500)
      dimension indsc(20,20,500),rate(20),rate2(500)
      dimension x2(20),pred2(20),sumbin(40),xmag2(40)
      dimension istart(500),iend(500)
      open(unit=10,file='rateinvflt.out',status='old')
c---- has 20 segments;  user can change
      mm= 20
      nn=20
      mmall=39
c test values for inputs
      facsm= 7000.
c--- facsm is smoothing factor 
      write(6,*) "enter jump prob, jump prob at iseg 10,facsm"
c--- user can specify different jumping probs for each segment boundary
      read(5,*) c0,c10,facsm
c  sample inputs
c   0.5 0.5 1000
c   0.7 0.7 2000
c   0.9 0.9 3000
c   0.95 9.95 5000
      xlentot= 20.*mm
      xlentemp= 0.
      totmo= 0.
c      totmo= 400.*15*1000.*1000.*xmu*sr
      xmu= 3.e10
c--------- segment length set here to 20 km
c--------  slip rate set here to 25 mm/yr
c--------   user can reset  theses values
      do 20 i=1,mm
      c(i)= c0
      xlen0(i)= 20.
      sliprate(i)= 25.
c      f1= 0.5*xlen0(1) +xlen0(1)*abs(float(itemp-iitemp))
      f1= 0.5*xlen0(1) + xlentemp
      f2= xlentot
      if(f2.ne.0.) frac= f1/f2
      if(f2.eq.0.) frac= 0.
      sliptap= 1.0
      if(frac.lt.0.15) sliptap= 0.75
      if(frac.lt.0.1) sliptap= 0.5
      if(frac.lt.0.05) sliptap=.25
      if(frac.gt.0.85) sliptap= 0.75
      if(frac.gt.0.9) sliptap= 0.5
      if(frac.gt.0.95) sliptap=.25
c---- no factor to adjust for fault
c      sliptap= sliptap*1.18
      sliprate(i)= sliprate(i)*sliptap
      totmo = totmo + sliprate(i)*xlen0(i)*15.*xmu*1000.
      write(6,*) f1,f2,frac,sliptap
c      if((i.eq.1).and.(j.eq.1)) sliptap=1.
c      write(11,*) ii,jj,i,sliptap,xlen(ii,jj),prob(ii,jj)
      b(i)= sliprate(i)
      xlentemp= xlentemp+ xlen0(i)
 20   continue
c--- append 0's in data vector for smoothing
      do 22 i=mm+1,mmall
      b(i)= 0.
  22  continue
ccccc
      c(10)= c10
c      c(9)= c10*2.
c      c(11)= c10*2.
      do 421 i=1,mm
      do 422 j=1,mm
      do 423 k=1,500
      indsc(i,j,k) =0
  423 continue
  422 continue
  421 continue
      do 425 i=1,mm
      do 426 iscen=1,500
      prob(i,iscen)= 0
  426 continue
 425  continue
      do 111 i=1,mmall
      do 112 j=1,mm
      a(i,j)=0.
      a2(i,j)= 0.
  112  continue
  111  continue
ccc
      write(6,*) "build scenarios"
      iscen= 0
      do 11 i=2,mm-2
      xlen(i,i+1)= 0.
      do 12 j=i+1,mm-1
      term1= 1.
      xlen(i,j)= xlen0(j-1)
      do 14 k=i,j-1
      term1= term1*c(k)
c---- check xlen indices
      xlen(i,j)= xlen(i,j)+ xlen0(j)
 14   continue
      iscen= iscen+1
      prob2= term1 - term1*sumprob(c(i-1),c(j))
      area= xlen(i,j)*15.
c---- assumes 15 km fault width; user can change
      xmag(iscen)= alog10(area) + 4.2
      xmo= 1.5*xmag(iscen)+ 9.05
      xmu= 3.0e10
      area= area*1000.*1000.
      slipt= xmo- alog10(area)-alog10(xmu)
      slipt= 10.**slipt
      slip(iscen) = 1000.*slipt
      do 194 ii=i,j
      prob(ii,iscen)= prob2
      do 195 jj=i,j     
       indsc(ii,jj,iscen)= 1
  195 continue
  194 continue
      istart(iscen)= i
      iend(iscen)= j
  12  continue
  11  continue
cccc   scenarios starting at segment 1
      i=1
      xlen(i,i+1)=0.
      do 116 j=i+1,mm-1
      term1= 1.
      xlen(i,j)= xlen0(j-1)
      do 114 k=i,j-1
      term1= term1*c(k)
      xlen(i,j)= xlen(i,j)+ xlen0(j)
 114  continue
      iscen= iscen+1
      prob2= term1 - term1*c(j) 
c      xlen(i,j)= xlen(i,j-1)+ xlen0(j-1)
      area= 15.*xlen(i,j)
      xmag(iscen)= alog10(area) + 4.2
      xmo= 1.5*xmag(iscen)+ 9.05
      xmu= 3.0e10
      area= area*1000.*1000.
      slipt= xmo- alog10(area)-alog10(xmu)
      slipt= 10.**slipt
      slip(iscen) = 1000.*slipt
      do 294 ii=i,j
      prob(ii,iscen)= prob2
      do 295 jj=i,j     
       indsc(ii,jj,iscen)= 1
  295 continue
  294 continue
      istart(iscen)= i
      iend(iscen)= j
 116  continue
c      for scenario starting at 1 and ending at mm
      i=1
      j=mm
      prob2= 1.
      xlen(i,j)= 0.
      do 440 ii=1,mm-1
      prob2= prob2*c(ii)
      xlen(i,j)= xlen(i,j) + xlen0(ii)
  440 continue
      iscen= iscen+1
      xlen(i,j)= xlen(i,j)+ xlen0(mm)
      area= 15.*xlen(i,j)
      xmag(iscen)= alog10(area) + 4.2
      xmo= 1.5*xmag(iscen)+ 9.05
      xmu= 3.0e10
      area= area*1000.*1000.
      slipt= xmo- alog10(area)-alog10(xmu)
      slipt= 10.**slipt
      slip(iscen) = 1000.*slipt
      do 394 ii=i,j
      prob(ii,iscen)= prob2
      do 395 jj=i,j     
       indsc(ii,jj,iscen)= 1
  395 continue
  394 continue
      istart(iscen)= i
      iend(iscen)= j
c---- scenarios ending at segment mm
      j=mm
      do 211 i=2,mm-1
      xlen(i,j)=xlen0(i-1)
      term1= 1.
      do 214 k=i,j-1
      term1= term1*c(k)
      xlen(i,j)= xlen(i,j)+ xlen0(i)
 214  continue
      iscen= iscen+1
      prob2= term1 - term1*c(i-1)
      area= 15.*xlen(i,j)
      xmag(iscen)= alog10(area) + 4.2
      xmo= 1.5*xmag(iscen)+ 9.05
      xmu= 3.0e10
      area= area*1000.*1000.
      slipt= xmo- alog10(area)-alog10(xmu)
      slipt= 10.**slipt
      slip(iscen) = 1000.*slipt
      do 494 ii=i,j
      prob(ii,iscen)= prob2
      do 495 jj=i,j     
       indsc(ii,jj,iscen)= 1
  495 continue
  494 continue
      istart(iscen)= i
      iend(iscen)= j
 211  continue     
c-----  for single segment ruptures i=j
      do 311 i=2,mm-1
      iscen= iscen+1
      prob(i,iscen)= 1.- c(i-1) - c(i) + c(i-1)*c(i)
      xlen(i,i)= xlen0(i-1)
      area= 15.*xlen(i,i)
      xmag(iscen)= alog10(area) +4.2
      xmo= 1.5*xmag(iscen) + 9.05
      xmu= 3.0e10
      area= area*1000.*1000.
      slipt= xmo- alog10(area)-alog10(xmu)
      slipt= 10.**slipt
      slip(iscen) = 1000.*slipt
      indsc(i,i,iscen)= 1
      istart(iscen)=i
      iend(iscen)= i
  311 continue
c---- for first segment rupture
      iscen= iscen+1
      prob(1,iscen)= 1.-c(1)
      xlen(1,1)= xlen0(1)
      area= 15.*xlen(1,1)
      xmag(iscen)= alog10(area) + 4.2
      xmo= 1.5*xmag(iscen)+ 9.05
      xmu= 3.0e10
      area= area*1000.*1000.
      slipt= xmo- alog10(area)-alog10(xmu)
      slipt= 10.**slipt
      slip(iscen) = 1000.*slipt
      indsc(1,1,iscen)= 1
      istart(iscen)= 1
      iend(iscen)= 1
c    for mm segment rupture
      iscen= iscen+1
      prob(mm,iscen)= 1.-c(mm-1)
      xlen(mm,mm)= xlen0(mm)
      area= 15.*xlen(mm,mm)
      xmag(iscen)= alog10(area) + 4.2
      xmo= 1.5*xmag(iscen)+ 9.05
      xmu= 3.0e10
      area= area*1000.*1000.
      slipt= xmo- alog10(area)-alog10(xmu)
      slipt= 10.**slipt
      slip(iscen) = 1000.*slipt
      indsc(mm,mm,iscen)= 1
      istart(iscen) =mm
      iend(iscen)= mm
      nscen= iscen
      do 40 i=1,nscen
      write(2,*) i,istart(i),iend(i),slip(i),xmag(i)
  40  continue
cccoc--- find matrix elements a(i,j)
      do 1 i=1,mm
      do 2 j=1,mm
cc  loop through all scenarios
      do 3 iscen= 1,nscen
      if(indsc(i,j,iscen).eq.1) then
c---- old code; works;  test if segments i,j are both between ii and jj
c------for slip taper at ends of rupture
      iitemp = istart(iscen)
      itemp =i
c      f1= 0.5*xlen0(1) +xlen0(1)*abs(float(itemp-iitemp))
      f1= 0.5*xlen0(1) + xlen0(1)*abs(float(itemp-iitemp))
      ii= istart(iscen)
      jj= iend(iscen)
      f2= xlen(ii,jj)
      if(f2.ne.0.) frac= f1/f2
      if(f2.eq.0.) frac= 0.
      sliptap= 1.0
      if(frac.lt.0.15) sliptap= 0.75
      if(frac.lt.0.1) sliptap= 0.5
      if(frac.lt.0.05) sliptap=.25
      if(frac.gt.0.85) sliptap= 0.75
      if(frac.gt.0.9) sliptap= 0.5
      if(frac.gt.0.95) sliptap=.25
c---- factor to get proper average slip
      sliptap= sliptap*1.18
      a(i,j)= a(i,j) +sliptap*slip(iscen)*prob(j,iscen)
      a2(i,j)= a2(i,j) + prob(j,iscen)
      endif
    4 continue
    3 continue
    2 continue
    1 continue
ccccccccccccc
      do 15 i=1,mm
      do 16 j=1,mm
      atmp(i,j)= a(i,j)
      if(a(i,j).eq.0.) write(6,*) "a(i,j)=0.",i,j
      if(a(i,j).lt.0.) then
      write(6,*) "a(i,j) less than 0",i,j
      read(5,*) idum
      endif
  16  continue
 15   continue
ccc   smoothing
      j=1
      do 244 i=mm+1,mmall
      a(i,j)= facsm
      a(i,j+1)= -facsm
      atmp(i,j)= a(i,j)
      atmp(i,j+1)= a(i,j+1)
      j=j+1
 244  continue
      m= mmall
      n= mm
      ncol= mm
ccc----------inversion stuff here
      ma= mmall
      write(6,*) "before nnls"
      call nnls(a,ma,m,n,b,x,rnorm,w,zz,index,mode)
      do 71 i=1,mm
      write(6,*) "x=",x(i)
      write(10,*) (i-1)*xlen0(i),x(i)
      write(10,*) (i-1)*xlen0(i)+xlen0(i),x(i)
 71   continue
      write(6,*) "mode",mode
      do 100 i=1,mm
      pred(i)=0.
      pred2(i)= 0.
      do 101 j=1,mm
      pred(i)= pred(i)+ x(j)*atmp(i,j)
 101  continue
 100  continue
       do 102 i=1,mm
      write(6,*) "a2(i,i)=", i,a2(i,i)
 102  continue
      sumres= 0.
      sumres2 = 0.
      do 200 i=1,mm
      write(1,*) sliprate(i),pred(i)
      write(6,*) "sliprate,pred=",i,sliprate(i),pred(i)
      write(10,*) (i-1)*xlen0(i),sliprate(i),pred(I)
      write(10,*) (i-1)*xlen0(i)+xlen0(i),sliprate(i),pred(I)
      sumres= sumres+ (sliprate(i)-pred(i))**2
 200  continue
      sumres= sqrt(sumres/float(mm))
      write(6,*) "sumres",sumres
      write(6,*) "c0,c10,facsm",c0,c10,facsm
      write(10,*) "c0,c10,facsm",c0,c10,facsm
      write(10,*) "sumres",sumres
c-- find rates of all rupture through segment i
      do 220 i=1,mm
      rate(i)= 0.
      do 221 j=1,mm
      rate(i)= rate(i)+x(j)*a2(i,j)
 221  continue
 220  continue
      write(3,*) "sumres",sumres
      write(3,*) "c0,c10,wl,facwl",c0,c10,wl,facwl
      do 202 i=1,mm
      write(3,*) i,rate(i)
      write(10,*) (i-1)*xlen0(i),rate(i)
      write(10,*) (i-1)*xlen0(i)+xlen0(i),rate(i)
 202  continue
      do 454 iscen=1,nscen
      rate2(iscen)=0.
 454  continue
      magmax=0
      do 401 iscen=1,nscen
      do 403 i=1,mm
      rate2(iscen)= rate2(iscen) + x(i)*prob(i,iscen)
  403 continue
      mag= iend(iscen)-istart(iscen) +1
      sumbin(mag)= sumbin(mag) +rate2(iscen)
      xmag2(mag)= xmag(iscen) 
      if(mag.gt.magmax) magmax=mag
  401 continue
      sumscen=0.
      do 410 iscen=1,nscen
      write(3,*) istart(iscen),iend(iscen),xmag(iscen),rate2(iscen)
      sumscen= sumscen+ rate2(iscen)
  410 continue
      write(6,*) "MFD",magmax
      write(3,*) "MFD",magmax
      sum=0.
      summor= 0.
      do 314 i=1,magmax
      write(3,*) xmag2(i),sumbin(i)
      write(6,*) xmag2(i),sumbin(i)
      write(10,*) xmag2(i),sumbin(i)
      sum= sum +sumbin(i)
      xmo= 1.5*xmag2(i) + 9.05
      xmo= 10.**xmo
      summor= summor + xmo*sumbin(i)
 314  continue
      sr= 25./1000.
      write(6,*) "total sum,sumscen,nscen",sum,sumscen,nscen
      write(6,*) "model and actual mo rates",summor,totmo
      end
      function sumprob(p1,p2) result(psum)
      psum= p1 + p2 -p1*p2
      end function
