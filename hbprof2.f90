module tableP2
  real :: zKuSJ(300)
  real :: dmJ(300)
  real :: rJ(300)
  integer :: nJ
end module tableP2
subroutine initP2
  use tableP2
  use tables2
  implicit none
  integer :: i
  do i=1,nbins
     zKuSJ(i)=zmin+(i-1)*dzbin
     dmJ(i)=d013Table(i,1)
     rJ(i)=pr13Table(i,1)
  end do
  nJ=nbins
  !print*, d013Table(200:240,1),nbins,nbinS2,nbinH
end subroutine initP2

subroutine iter_profcv(btop,bzd,bcf,bsfc,zKuL,zKaL,dr,n1d,eps,imu,&
     dn1d,dm1d,rrate1d,zKuC,zKaSim,epst,piaKu,piaKa)
  use tables2
  use tableP2
  implicit none
  integer :: btop, bzd, bcf, bsfc, n1d, imu
  real :: zKuL(n1d), dr, zKaL(n1d)
  real :: dmCoeff(2)=(/0.027773772993636318,-0.6624076086959071/)
  real,intent(out) :: dm1d(n1d)
  real,intent(out) :: epst,piaKu,piaKa
  real,intent(out) :: dn1d(n1d), rrate1d(n1d), zKuC(n1d), zKaSim(n1d)
  real :: dzKa(n1d), dns(n1d)
  real :: eps, dm, pia
  integer :: it, ik, k, n1, n1H
  real :: ztrueS,ztrueH,ztrue,f, attKuH, attKaH, attKuS, attKaS
  real :: ftran, zetaS, zeta1d(n1d), q
  real :: piaKuS, piaKaS, beta, piamax, attKu, attKa, dn
  real :: dnCoeff(2)=(/-0.00570364,  0.13319214/)
  piaKu=0
  !dmCoeff=array([ 0.02893781, -0.69481455])

  zKuC=zKuL
  dns=0
  do it=1,2
     piaKu=0.0
     piaKa=0.0
     do k=btop,bzd-1
        if (zKuC(k+1).gt.10) then
           if (k.lt.bzd) then
              ztrue=zKuC(k+1)
           end if
           if(ztrue.lt.40) then
              ztrueS=ztrue
              ztrueH=0.0
           else
              ztrueS=40
              ztrueH=10*log10(10**(0.1*ztrue)-10**(0.1*39.999))
           endif
           !n1,n2=bisection(zKuSJ,ztrueS-10*dns[k+1])
           call bisection2(zKuSJ,nbinS2,ztrueS-10*dns(k+1), n1)
           attKuS=att13TableS2(n1,imu)*10**dns(k+1)
           attKaS=att35TableS2(n1,imu)*10**dns(k+1)
           !n1H,n2H=bisection(zKuHJ,ztrueH)
           call bisection2(zKuSJ,nbinS2,ztrueH, n1H)
           attKuH=att13TableH(n1H,imu)
           attKaH=att35TableH(n1H,imu)
           piaKu=piaKu+(attKuS+attKuH)*dr
           piaKa=piaKa+(attKaS+attKaH)*dr
           zKuC(k+1)=zKuL(k+1)+piaKu
           zKaSim(k+1)=10*log10(10**(0.1*z35TableS2(n1,imu)+dns(k+1))+&
                10**(0.1*z35TableH(n1H,imu)))-piaKa
           if(zKaL(k+1).gt.10) then
              dZKa(k+1)=zKaSim(k+1)-zKaL(k+1)
           end if
           piaKu=piaKu+(attKuS+attKuH)*dr
           piaKa=piaKa+(attKaS+attKaH)*dr
           dm1d(k+1)=0.5*(d013TableS2(n1,imu)+d013TableH(n1H,imu))
           dn1d(k+1)=dns(k+1)
           rrate1d(k+1)=pr13TableS2(n1,imu)*10**dns(k+1)+pr13TableH(n1H,imu)
        end if
     end do
     if(it.eq.1) then
        do k=btop,bzd
           dns(k+1)=dns(k+1)-0.6*dZKa(k+1)
           if(dns(k+1).gt.0.975) dns(k+1)=0.975
           if(dns(k+1).lt.-1.0) dns(k+1)=-1.0
        end do
     end if
  end do
  
  piaKuS=piaKu+0.0
  piaKaS=piaKa+0.0

  zeta1d=0
  beta=0.76
  
  piamax=52-zKuL(bcf+1)
  q=0.2*log(10.0)
  zKuC(bzd+1:bcf+1)=zKuL(bzd+1:bcf+1)+piaKuS
  epst=eps
  attKu=0.0
  attKa=0.0
  f=1.0
  do it=1,2
     zetaS=0.0
     piaKa=piaKaS+0.0
     do k=bzd,bcf
        if(zKuC(k+1)>0) then
           ztrue=zKuC(k+1)
           ftran=1
           dn=1.5*dnCoeff(1)*ztrue+dnCoeff(2)+0.2*log(epst)
           if(10*dn+ztrue.gt.zKuSJ(nbins)) dn=(zKuSJ(nbins)-ztrue)/10.0
           if(10*dn+ztrue.lt.zKuSJ(1)) dn=(zKuSJ(1)-ztrue)/10.0
           n1=((ztrue-10*dn-zmin)/dzbin)+1
           if(n1.lt.1) n1=1
           if(n1.gt.nbins) n1=nbins
           dm=d013Table(n1,imu)
           attKu=att13Table(n1,imu)*10**dn
           attKa=att35Table(n1,imu)*10**dn
           zetaS=zetaS+attKu/10**(zKuC(k+1)*0.1*beta)*&
                10**(zKuL(k+1)*0.1*beta)*dr
           rrate1d(k+1)=pr13Table(n1,imu)*10**dn
           dn1d(k+1)=dn
           dm1d(k+1)=dm
           piaKa=piaKa+attKa*2*dr
        end if
        zeta1d(k+1)=zetaS
     end do
     eps=min(1.0,(1-10**(-0.1*beta*f*piamax))/(q*beta*zetaS))
     epst=epst*eps
     do k=bzd,bcf
        zKuC(k+1)=zKuL(k+1)+piaKuS-10/beta*log10(1-eps*q*beta*zeta1d(k+1))
     end do
     piaKu=(-10/beta*log10(1-eps*q*beta*zeta1d(bcf+1)))+attKu*(bsfc-bcf)*2*dr
  end do
end subroutine iter_profcv


subroutine iter_profst(btop,bzd,bb,bbt,bbb,bcf,bsfc,zKuL,zKaL,dr,n1d,eps,imu,&
     dn1d,dm1d,rrate1d,zKuC,zKaSim,epst,piaKu,piaKa)
  use tables2
  use tableP2
  implicit none
  integer :: btop, bzd, bcf, bsfc, n1d, imu, bb,bbt,bbb
  real :: zKuL(n1d), dr, zKaL(n1d)
  real :: dmCoeff(2)=(/0.027773772993636318,-0.6624076086959071/)
  real,intent(out) :: dm1d(n1d)
  real,intent(out) :: epst,piaKu,piaKa
  real,intent(out) :: dn1d(n1d), rrate1d(n1d), zKuC(n1d), zKaSim(n1d)
  real :: dzKa(n1d), dns(n1d)
  real :: eps, dm, pia
  integer :: it, ik, k, n1, n1H
  real :: ztrueS,ztrueH,ztrue,f, attKuH, attKaH, attKuS, attKaS
  real :: ftran, zetaS, zeta1d(n1d), q
  real :: piaKuS, piaKaS, beta, piamax, attKu, attKa, dn
  real :: dnCoeff(2)=(/-0.00570364,  0.13319214/)
  integer :: n1S,n1BB,n1b
  real :: attKuBB, attKaBB, Z1, Z2
  piaKu=0
  !dmCoeff=array([ 0.02893781, -0.69481455])

  zKuC=zKuL
  dns=0
  do it=1,2
     piaKu=0.0
     piaKa=0.0
     do k=btop,bbt-1
        if (zKuC(k+1).gt.10) then
           if (k.lt.bzd) then
              ztrue=zKuC(k+1)
           end if
           if(ztrue.lt.40) then
              ztrueS=ztrue
              ztrueH=0.0
           else
              ztrueS=40
              ztrueH=10*log10(10**(0.1*ztrue)-10**(0.1*39.999))
           endif
           !n1,n2=bisection(zKuSJ,ztrueS-10*dns[k+1])
           call bisection2(zKuSJ,nbinS2,ztrueS-10*dns(k+1), n1)
           attKuS=att13TableS2(n1,imu)*10**dns(k+1)
           attKaS=att35TableS2(n1,imu)*10**dns(k+1)
           !n1H,n2H=bisection(zKuHJ,ztrueH)
           call bisection2(zKuSJ,nbinS2,ztrueH, n1H)
           attKuH=att13TableH(n1H,imu)
           attKaH=att35TableH(n1H,imu)
           piaKu=piaKu+(attKuS+attKuH)*dr
           piaKa=piaKa+(attKaS+attKaH)*dr
           zKuC(k+1)=zKuL(k+1)+piaKu
           zKaSim(k+1)=10*log10(10**(0.1*z35TableS2(n1,imu)+dns(k+1))+&
                10**(0.1*z35TableH(n1H,imu)))-piaKa
           if(zKaL(k+1).gt.10) then
              dZKa(k+1)=zKaSim(k+1)-zKaL(k+1)
           end if
           piaKu=piaKu+(attKuS+attKuH)*dr
           piaKa=piaKa+(attKaS+attKaH)*dr
           dm1d(k+1)=0.5*(d013TableS2(n1,imu)+d013TableH(n1H,imu))
           dn1d(k+1)=dns(k+1)
           rrate1d(k+1)=pr13TableS2(n1,imu)*10**dns(k+1)+pr13TableH(n1H,imu)
        end if
     end do
     if(it.eq.1) then
        do k=btop,bbt-1
           dns(k+1)=dns(k+1)-0.6*dZKa(k+1)
           if(dns(k+1).gt.0.975) dns(k+1)=0.975
           if(dns(k+1).lt.-1.0) dns(k+1)=-1.0
        end do
     end if
  end do
  print*, bbt, bb
  do k=bbt,bb
     if(zKuC(k+1)>10) then
        f=(k-bbt+0.0)/(bb-bbt)
        if(k.eq.bbt) then
           f=0.875
        else
           f=1
        end if
        Z1=log10((1-f)*10**(0.1*zKuC(k+1))+1e-9)*10.
        Z2=log10(f*10**(0.1*zKuC(k+1))+1e-9)*10.
        n1S=(Z1-zmin-10*dns(bbt))/dzbin+1
        n1BB=(Z2-zmin-10*dns(bbt))/dzbin+1
        if(n1S.gt.nbinS2) n1S=nbinS2
        if(n1BB.gt.nbins) n1BB=nbins
        if(n1S.lt.1) n1S=1
        if(n1BB.lt.1) n1BB=1
        print*, n1S,n1BB,f, Z1,Z2,zKuC(k+1)
        attKuS=att13TableS2(n1S,imu)*10**dns(bbt)
        attKaS=att35TableS2(n1S,imu)*10**dns(bbt)
        attKuBB=att13TableBB(n1S,imu)*10**dns(bbt)
        attKaBB=att35TableBB(n1S,imu)*10**dns(bbt)
        piaKu=piaKu+(attKuS+attKuBB)*dr
        piaKa=piaKa+(attKaS+attKaBB)*dr
        zKuC(k+1)=zKuL(k+1)+piaKu
        zKaSim(k+1)=10*log10(10**(0.1*z35TableS2(n1S,imu)+dns(bbt))+&
             10**(0.1*z35TableBB(n1BB,imu)+dns(bbt)))-piaKa
        piaKu=piaKu+(attKuS+attKuBB)*dr
        piaKa=piaKa+(attKaS+attKaBB)*dr
        dm1d(k+1)=(1-f)*d013TableS2(n1S,imu)+f*d013TableBB(n1BB,imu)
        dn1d(k+1)=dns(bbt)
        rrate1d(k+1)=pr13TableS2(n1S,imu)*10**dns(bbt)+pr13TableBB(n1BB,imu)*10**dns(bbt)
     end if
  end do

  do k=bb+1,bbb-1
     if(zKuC(k)>100) then
        f=(k-bb+0.0)/(bbb-bb)
        f=0.99
        Z1=(1-f)*zKuC(k+1)
        Z2=(f)*zKuC(k)
        ztrue=Z2
        n1=(Z1-zmin-10*dns(bbt))/dzbin+1
        dn=1.5*dnCoeff(1)*Z2+dnCoeff(2)+0.2*log(epst)
        if(10*dn+ztrue.gt.zKuSJ(nbins)) dn=(zKuSJ(nbins)-ztrue)/10.0
        if(10*dn+ztrue.lt.zKuSJ(1)) dn=(zKuSJ(1)-ztrue)/10.0
        n1b=(Z2-zmin-10*dn)/dzbin+1
        if(n1.gt.nbins) n1=nbins
        if(n1b.gt.nbins) n1b=nbins
        if(n1.lt.1) n1=1
        if(n1b.lt.1) n1b=1
        attKu=att13Table(n1b,imu)*10**dn
        attKa=att35Table(n1b,imu)*10**dn
        attKuBB=att13TableBB(n1,imu)*10**dns(bbt)
        attKaBB=att35TableBB(n1,imu)*10**dns(bbt)
        piaKu=piaKu+(attKu+attKuBB)*dr
        piaKa=piaKa+(attKa+attKaBB)*dr
        zKuC(k+1)=zKuL(k+1)+piaKu
        !zKaSim(k+1)=10*log10(10**(0.1*z35TableS2(n1,imu)+dns(bbt))+&
        !     10**(0.1*z35TableBB(n1H,imu)+dns(bbt)))-piaKa
        piaKu=piaKu+(attKu+attKuBB)*dr
        piaKa=piaKa+(attKa+attKaBB)*dr
        dm1d(k+1)=(1-f)*d013TableBB(n1,imu)+f*d013Table(n1b,imu)
        dn1d(k+1)=(1-f)*dns(bbt)+f*dn
        rrate1d(k+1)=pr13TableBB(n1,imu)*10**dns(bbt)+pr13Table(n1b,imu)*10**dn
     end if
  end do
  
  piaKuS=piaKu+0.0
  piaKaS=piaKa+0.0

  zeta1d=0
  beta=0.76
  
  piamax=52-zKuL(bcf+1)
  q=0.2*log(10.0)
  zKuC(bzd+1:bcf+1)=zKuL(bzd+1:bcf+1)+piaKuS
  epst=eps
  attKu=0.0
  attKa=0.0
  f=1.0
  do it=1,2
     zetaS=0.0
     piaKa=piaKaS+0.0
     do k=bb+1,bcf
        if(zKuC(k+1)>0) then
           ztrue=zKuC(k+1)
           ftran=1
           dn=1.5*(dnCoeff(1)*ztrue+dnCoeff(2))+0.2*log(epst)
           if(10*dn+ztrue.gt.zKuSJ(nbins)) dn=(zKuSJ(nbins)-ztrue)/10.0
           if(10*dn+ztrue.lt.zKuSJ(1)) dn=(zKuSJ(1)-ztrue)/10.0
           n1=((ztrue-10*dn-zmin)/dzbin)+1
           if(n1.lt.1) n1=1
           if(n1.gt.nbins) n1=nbins
           dm=d013Table(n1,imu)
           attKu=att13Table(n1,imu)*10**dn
           attKa=att35Table(n1,imu)*10**dn
           zetaS=zetaS+attKu/10**(zKuC(k+1)*0.1*beta)*&
                10**(zKuL(k+1)*0.1*beta)*dr
           rrate1d(k+1)=pr13Table(n1,imu)*10**dn
           dn1d(k+1)=dn
           dm1d(k+1)=dm
           piaKa=piaKa+attKa*2*dr
        end if
        zeta1d(k+1)=zetaS
     end do
     eps=min(1.0,(1-10**(-0.1*beta*f*piamax))/(q*beta*zetaS))
     epst=epst*eps
     do k=bzd,bcf
        zKuC(k+1)=zKuL(k+1)+piaKuS-10/beta*log10(1-eps*q*beta*zeta1d(k+1))
     end do
     piaKu=(-10/beta*log10(1-eps*q*beta*zeta1d(bcf+1)))+attKu*(bsfc-bcf)*2*dr
  end do
end subroutine iter_profst

