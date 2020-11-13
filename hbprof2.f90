module tableP2
  real :: zKuSJ(300)
end module tableP2
subroutine initP2
  use tableP2
  use tables2
  implicit none
  integer :: i
  do i=1,nbins
     zKuSJ(i)=zmin+(i-1)*dzbin
  end do
  print*, d013Table(200:240,1),nbins,nbinS2,nbinH
end subroutine initP2

subroutine iter_prof(btop,bzd,bcf,bsfc,zKuL,zKaL,dr,n1d,eps,imu,&
     dn1d,dm1d,rrate1d,zKuC,zKaSim,epst)
  use tables2
  use tableP2
  implicit none
  integer :: btop, bzd, bcf, bsfc, n1d, imu
  real :: zKuL(n1d), dr, zKaL(n1d)
  real :: dmCoeff(2)=(/0.027773772993636318,-0.6624076086959071/)
  real,intent(out) :: dm1d(n1d)
  real,intent(out) :: epst
  real,intent(out) :: dn1d(n1d), rrate1d(n1d), zKuC(n1d), zKaSim(n1d)
  real :: piaKu, piaKa, dzKa(n1d), dns(n1d)
  real :: eps, dm, pia
  integer :: it, k, n1, n1H
  real :: ztrueS,ztrueH,ztrue,f, attKuH, attKaH, attKuS, attKaS
  real :: ftran, zetaS, zeta1d(n1d), q
  real :: piaKuS, piaKaS, beta, piamax, attKu, attKa, dn
  piaKu=0
  
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
           !else
           !f=(bzd-k)/6.0
           !ztrue=10*log10(10**(0.1*zKuC(k+1))*f+0.1)
           !end if
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
           !if(k<bzd+6) then
           !   ftran=0.85+1/6.0*(k-bzd)
           !else
           !   ftran=1.0
           !end if
           ftran=1
           dm=ftran*1.05*exp(dmCoeff(1)*ztrue+dmCoeff(2))/epst**0.25
           !print*,ztrue,dm,k,it,zKuL(k+1),epst
           if(dm.gt.3.0) dm=3.0
           if(dm.lt.0.31) dm=0.31
           call bisection2(d013Table(:,1),nbins,dm,n1)
           dn=(ztrue-zKuSJ(n1))/10.0
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
     !print*, piamax
     eps=min(1.0,(1-10**(-0.1*beta*f*piamax))/(q*beta*zetaS))
     !if(zetaS.gt.1e-3) then
     !   f=1.05*(-10/beta*log10(1-eps*q*beta*zeta1d(bcf+1))/&
     !   ((-10/beta*log10(1-eps*q*beta*zeta1d(bcf+1)))+attKu*(bsfc-bcf)*2*dr)
     !end if
     !if(eps.gt.10*f) eps=10*f
     !if(eps.lt.0.01) eps=0.01
     epst=epst*eps
     do k=bzd,bcf
        zKuC(k+1)=zKuL(k+1)+piaKuS-10/beta*log10(1-eps*q*beta*zeta1d(k+1))
     end do
     piaKu=(-10/beta*log10(1-eps*q*beta*zeta1d(bcf+1)))+attKu*(bsfc-bcf)*2*dr
  end do
end subroutine iter_prof

