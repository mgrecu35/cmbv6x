module tableP2
  real :: zKuSJ(300)
  real :: dmJ(300)
  real :: rJ(300)
  real :: attKaJ(300)
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
     attKaJ(i)=att35Table(i,1)
     rJ(i)=pr13Table(i,1)
  end do
  nJ=nbins
  !print*, d013Table(200:240,1),nbins,nbinS2,nbinH
end subroutine initP2

subroutine prof1d(btop,bzd,bcf,bsfc,binBB,binBBT,zKuL,zKaL,pType,dr,n1d,eps,imu,&
     dn1d,dm1d,rrate1d,zKuC,zKaSim,epst,piaKu,piaKa,dnCoeff_new,dn,dsrtPIA,relSRT)
  implicit none
  real :: dnp(n1d), dn
  integer :: btop, bzd, bcf, bsfc, n1d, imu, bb,bbt,bbb, pType
  real :: zKuL(n1d), dr, zKaL(n1d)
  real :: dnCoeff_new(2)
  integer :: binBB,binBBT
  real,intent(out) :: dm1d(n1d)
  real,intent(out) :: epst,piaKu,piaKa
  real,intent(out) :: dn1d(n1d), rrate1d(n1d), zKuC(n1d), zKaSim(n1d)
  real :: eps
  real :: dZKa(n1d), dPIA
  real :: dn1d1(n1d),dm1d1(n1d),rrate1d1(n1d),zKuC1(n1d),zKaSim1(n1d),epst1,piaKu1,piaKa1
  real :: ddn, gradZS, gradZe, dzdn(n1d,n1d)
  integer :: i, nbz, relSRT
  real :: sigma, dnbz(n1d), dsrtPIA,   dpiadn
  !print*, bbb, binBB

  dnp=0
  if (pType.eq.2) then
     !print*, pType, dnCoeff_new,dn
     call iter_profcv(btop,bzd,bcf,bsfc,zKuL,zKaL,dr,n1d,eps,imu,&
          dn1d,dm1d,rrate1d,zKuC,zKaSim,epst,piaKu,piaKa,ptype,dnCoeff_new,dn,dnp,dzdn)
     !print*, rrate1d(bcf+1), epst
  else
     !dnCoeff_new(1)=-0.01257341
     !dnCoeff_new(2)=-0.00933038
     if (binBB.lt.1) then
        call iter_profcv(btop,bzd,bcf,bsfc,zKuL,zKaL,dr,n1d,eps,imu,&
             dn1d,dm1d,rrate1d,zKuC,zKaSim,epst,piaKu,piaKa,ptype,dnCoeff_new,dn,dnp,dzdn)
        call iter_profcv(btop,bzd,bcf,bsfc,zKuL,zKaL,dr,n1d,eps,imu,&
             dn1d1,dm1d1,rrate1d1,zKuC1,zKaSim1,epst1,piaKu1,piaKa1,ptype,dnCoeff_new,dn-0.1,dnp,dzdn)
        bbb=bzd+3
     else
        bbb=binBB+2
        bbt=binBBT
        bb=binBB        
        call iter_profst(btop,bzd,bb,bbt,bbb,bcf,bsfc,zKuL,zKaL,dr,n1d,eps,imu,&
             dn1d,dm1d,rrate1d,zKuC,zKaSim,epst,piaKu,piaKa,dn,dnCoeff_new,dnp,dzdn)
        call iter_profst(btop,bzd,bb,bbt,bbb,bcf,bsfc,zKuL,zKaL,dr,n1d,eps,imu,&
             dn1d1,dm1d1,rrate1d1,zKuC1,zKaSim1,epst1,piaKu1,piaKa1,dn-0.1,dnCoeff_new,dnp,dzdn)
     end if
     gradZS=1
     gradZe=0.
     do i=bbb+3,bcf
        if(zKuL(i+1).gt.10.0.and.zKaL(i+1).gt.10) then
           gradZS=gradZS+(-(zKaSim1(i+1)-zKaSim(i+1))/0.1)**2
           gradZe=gradZe+(-(zKaSim1(i+1)-zKaSim(i+1))/0.1)*0.1*(zKaL(i+1)-zKaSim(i+1))
        endif
     enddo
     if(relSRT.eq.1) then
        dpiadn=(piaKa1-piaKu1-piaKa+piaKu)/(-0.1)
        gradZS=gradZS+dpiadn**2
        gradZe=gradZe-dpiadn*(piaKa-piaKu-dsrtPIA)
        !print*, 'dpiadn',piaKa-piaKu, dsrtPIA, dpiadn, piaKu, piaKa, binBBT
     end if

     if(bbb.lt.bcf-3) then
        sigma=1.5
        nbz=bcf-bbb
        call gauss_newton(dzdn(bbb+2:bcf+1,bbb+2:bcf+1), zKasim(bbb+2:bcf+1), zKaL(bbb+2:bcf+1), nbz, sigma,&
             dnbz)
        dnp(bbb+2:bcf+1)=dnbz(1:nbz)
     end if
     do i=bbb+2,bcf
        if(dnp(i+1).gt.1) dnp(i+1)=1
        if(dnp(i+1).lt.-1) dnp(i+1)=-1
     enddo
     !dy=transpose(dZdn)*(dY)-0.1*invCov*(xsol-xMean);
     !A=transpose(dZdn)*dZdn+0.1*invCov;
     !dn=A\dy;
     
     ddn=0.975*gradZe/gradZS
     !ddn=0.0
     if(ddn.lt.-1.5) ddn=-1.5
     if(ddn.gt.1.5) ddn=1.5
     !dnp=0.0
     !print*, ddn, bbb, bcf
     if (binBB.lt.1) then
        call iter_profcv(btop,bzd,bcf,bsfc,zKuL,zKaL,dr,n1d,eps,imu,&
             dn1d,dm1d,rrate1d,zKuC,zKaSim,epst,piaKu,piaKa,ptype,dnCoeff_new,dn+ddn,dnp,dzdn)
     else
        bbb=binBB+2
        bbt=binBBT
        bb=binBB
        !dnp=0
        call iter_profst(btop,bzd,bb,bbt,bbb,bcf,bsfc,zKuL,zKaL,dr,n1d,eps,imu,&
             dn1d,dm1d,rrate1d,zKuC,zKaSim,epst,piaKu,piaKa,dn+ddn,dnCoeff_new,dnp,dzdn)
     end if
     
  end if
end subroutine prof1d
subroutine iter_profcv(btop,bzd,bcf,bsfc,zKuL,zKaL,dr,n1d,eps,imu,&
     dn1d,dm1d,rrate1d,zKuC,zKaSim,epst,piaKu,piaKa,itype,dnCoeff_new,dncv,dnp,dzdn)
  use tables2
  use tableP2
  use ran_mod
  implicit none
  integer :: btop, bzd, bcf, bsfc, n1d, imu
  real :: zKuL(n1d), dr, zKaL(n1d), dnp(n1d)
  !real :: dmCoeff(2)=(/0.027773772993636318,-0.6624076086959071/)
  real,intent(out) :: dm1d(n1d)
  real,intent(out) :: epst,piaKu,piaKa
  real,intent(out) :: dn1d(n1d), rrate1d(n1d), zKuC(n1d), zKaSim(n1d), dzdn(n1d,n1d)
  real :: dnCoeff_new(2)
  real :: dzKa(n1d), dns(n1d)
  real :: eps, dm, pia, dncv
  integer :: it, ik, k, n1, n1H, itype
  real :: ztrueS,ztrueH,ztrue,f, attKuH, attKaH, attKuS, attKaS
  real :: ftran, zetaS, zeta1d(n1d), q
  real :: piaKuS, piaKaS, beta, piamax, attKu, attKa, dn, dni
  real :: dnCoeff(2)=(/-0.00570364,  0.13319214/)
  real :: dmCoeff(3)= (/-2.56885571e-04,  7.18909743e-02, -1.60889523e+00/)
  !(/-4.57599678e-04,  7.852247958e-02, -1.78316499e+00/)
  !(/-4.57599678e-04,  8.52247958e-02, -1.78316499e+00/)
  real :: rn, dm_old
  real :: zka1,zka2,attka2
  integer :: kk
  piaKu=0
  !dmCoeff=array([ 0.02893781, -0.69481455])
  dzdn=0.0
  !if(itype.eq.1) then
  dnCoeff=dnCoeff_new
  !end if
  zKuC=zKuL
  dns=0
  dZKa=0
  do it=1,2
     piaKu=0.0
     piaKa=0.0
     do k=btop,min(bcf,bzd-1)
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
        do k=btop,min(bcf,bzd-1)
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
  attKa=0.0
  do it=1,4
     zetaS=0.0
     piaKa=piaKaS+0.0
     do k=bzd,bcf
        if(zKuC(k+1)>0) then
           ztrue=zKuC(k+1)
           ftran=1
           if(itype.eq.2) then
              dn=1.5*(dnCoeff(1)*ztrue+dnCoeff(2))+0.2*log(epst)+dnp(k+1)+dncv
           else
              dn=1.0*(dnCoeff(1)*ztrue+dnCoeff(2))+0.2*log(epst)+dnp(k+1)+dncv
           end if
           n1=(ztrue-1*dn-zmin)/dzbin
           dm=d013Table(n1,imu)
           !-----------------------------------------------------
           !dm=0.9*(exp(dmCoeff(1)*ztrue**2+dmCoeff(2)*ztrue+dmCoeff(3)))
           !dm_old=dm
           !if(dm.gt.0.5.and.dm.lt.2.8) then
           !if(it.eq.2) then
           !   rn=normal2(0.0,1.0)
           !   dm=0.9*(0.875*dm+0.125*(1.3+0.25*normal2(0.0,1.0)))
           !end if
           !end if
           !if(dm<0.31) dm=0.31
           !if(dm>3) dm=3
           !call bisection2(d013Table(:,imu),nbins,dm, n1)
           !dn=(ztrue-zKuSJ(n1))/10
           !dn=0.5*(dn+dn)
           !-------------------------------------------------------
           if(10*dn+ztrue.gt.zKuSJ(nbins)) dn=(zKuSJ(nbins)-ztrue)/10.0
           if(10*dn+ztrue.lt.zKuSJ(1)) dn=(zKuSJ(1)-ztrue)/10.0
           n1=((ztrue-10*dn-zmin)/dzbin)+1
           if(n1.lt.1) n1=1
           if(n1.gt.nbins) n1=nbins
           dm=d013Table(n1,imu)
           attKu=att13Table(n1,imu)*10**dn
           attKa=att35Table(n1,imu)*10**dn
           zKa1=z35Table(n1,imu)+10*dn
           if(it.eq.4) then
              if(n1.lt.nbins-4) then
                 zKa2=z35Table(n1+4,imu)+10*(dn-0.1)
                 attKa2=att35Table(n1+4,imu)*10**(dn-0.1)
                 dzdn(k+1,k+1)=dzdn(k+1,k+1)+(zka2-zka1)/(-0.1)-(attKa2-attKa)/(-0.1)*dr
                 do kk=k+1,bcf
                    dzdn(kk+1,k+1)=dzdn(kk+1,k+1)-(attKa2-attKa)/(-0.1)*dr
                 end do
              end if
           end if
           zetaS=zetaS+attKu/10**(zKuC(k+1)*0.1*beta)*&
                10**(zKuL(k+1)*0.1*beta)*dr
           rrate1d(k+1)=pr13Table(n1,imu)*10**dn
           dn1d(k+1)=dn
           dm1d(k+1)=dm
           piaKa=piaKa+attKa*dr
           zKaSim(k+1)=z35Table(n1,imu)+10*dn-piaKa
           piaKa=piaKa+attKa*dr
        else
           attKu=0.0
           attKa=0.0
        end if
        zeta1d(k+1)=zetaS
     end do
     eps=min(1.0,(1-10**(-0.1*beta*f*piamax))/(q*beta*zetaS))
     epst=epst*eps
     do k=bzd,bcf
        zKuC(k+1)=zKuL(k+1)+piaKuS-10/beta*log10(1-eps*q*beta*zeta1d(k+1))
     end do
     piaKu=(-10/beta*log10(1-eps*q*beta*zeta1d(bcf+1)))+attKu*(bsfc-bcf)*2*dr
     if(isnan(piaKu).eqv..true.) then
        print*, eps, zetaS, piamax, zeta1d(bcf+1)
        print*,'conv'
        stop
     end if
     if(isnan(piaKa).eqv..true.) then
        !print*, eps, zetaS, piamax, zeta1d(bcf+1), piaKa, attKa
        print*,dns(btop:min(bcf,bzd-1))
        print*,zKuC(btop:min(bcf,bzd-1))
        print*, 'cv',eps, zetaS, piamax, zeta1d(bcf+1), piaKa, piaKas,attKa, bsfc
        stop
     end if
  end do
end subroutine iter_profcv


subroutine iter_profst(btop,bzd,bb,bbt,bbb,bcf,bsfc,zKuL,zKaL,dr,n1d,eps,imu,&
     dn1d,dm1d,rrate1d,zKuC,zKaSim,epst,piaKu,piaKa,dnst,dnCoeff_new,dnp,dzdn)
  use tables2
  use tableP2
  use ran_mod
  implicit none
  integer :: btop, bzd, bcf, bsfc, n1d, imu, bb,bbt,bbb
  real :: zKuL(n1d), dr, zKaL(n1d), dzdn(n1d,n1d)
  real :: dnCoeff_new(2)
  !real :: dmCoeff(2)=(/0.027773772993636318,-0.6624076086959071/)
  real,intent(out) :: dm1d(n1d)
  real,intent(out) :: epst,piaKu,piaKa
  real,intent(out) :: dn1d(n1d), rrate1d(n1d), zKuC(n1d), zKaSim(n1d)
  real :: dzKa(n1d), dns(n1d), dnp(n1d)
  real :: eps, dm, pia, dnst
  integer :: it, ik, k, n1, n1H
  real :: ztrueS,ztrueH,ztrue,f, attKuH, attKaH, attKuS, attKaS
  real :: ftran, zetaS, zeta1d(n1d), q
  real :: piaKuS, piaKaS, beta, piamax, attKu, attKa, dn, dni
  real :: dnCoeff(2)=(/-0.00570364,  0.13319214/)
  integer :: n1S,n1BB,n1b
  real :: attKuBB, attKaBB, Z1, Z2
  real :: dmCoeff(3)=(/-2.56885571e-04,  7.18909743e-02, -1.60889523e+00/)
  !=(/-4.57599678e-04,  7.752247958e-02, -1.78316499e+00/)
!  array([-2.56885571e-04,  7.18909743e-02, -1.60889523e+00])

  real :: rn, r1, zka1, zka2, attKa2
  integer :: kk
  piaKu=0
  !dmCoeff=array([ 0.02893781, -0.69481455])
  dnCoeff=dnCoeff_new
  zKuC=zKuL
  dns=0
  dzdn=0
  dZKa=0
  do it=1,2
     piaKu=0.0
     piaKa=0.0
     do k=btop,min(bcf,bbt-1)
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
  !print*, bbt, bb
  do k=bbt,min(bb,bcf)
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
        !print*, n1S,n1BB,f, Z1,Z2,zKuC(k+1)
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
  epst=eps
  do k=bb+1,min(bbb-1,bcf)
     if(zKuC(k)>0) then
        f=(k-bb+0.0)/(bbb-bb)
        f=0.99
        Z1=(1-f)*zKuC(k+1)
        Z2=(f)*zKuC(k+1)
        ztrue=Z2
        n1=(Z1-zmin-10*dns(bbt))/dzbin+1
        !dn=1.0*(dnCoeff(1)*Z2+dnCoeff(2))+0.2*log(epst)
        dn=1.0*(dnCoeff(1)*ztrue+dnCoeff(2))+dnp(k+1)+0.2*log(epst)+dnst
        if(10*dn+ztrue.gt.zKuSJ(nbins)) dn=(zKuSJ(nbins)-ztrue)/10.0
        if(10*dn+ztrue.lt.zKuSJ(1)) dn=(zKuSJ(1)-ztrue)/10.0
        n1b=((Z2-zmin-10*dn)/dzbin)+1
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
        zKaSim(k+1)=10*log10(10**(0.1*z35TableBB(n1,imu)+dns(bbt))+&
             10**(0.1*z35Table(n1b,imu)+dn))-piaKa
        piaKu=piaKu+(attKu+attKuBB)*dr
        piaKa=piaKa+(attKa+attKaBB)*dr
        dm1d(k+1)=(1-f)*d013TableBB(n1,imu)+f*d013Table(n1b,imu)
        dn1d(k+1)=(1-f)*dns(bbt)+f*dn
        rrate1d(k+1)=pr13TableBB(n1,imu)*10**dns(bbt)+pr13Table(n1b,imu)*10**dn
        !print*, zKuC(k+1),Z1,Z2,dn,n1b,n1
     end if
  end do
  
  piaKuS=piaKu+0.0
  piaKaS=piaKa+0.0

  zeta1d=0
  beta=0.76
  
  piamax=52-zKuL(bcf+1)
  q=0.2*log(10.0)
  zKuC(bzd+1:bcf+1)=zKuL(bzd+1:bcf+1)+piaKuS
 
  attKu=0.0
  attKa=0.0
  f=1.0
  do it=1,3
     zetaS=0.0
     piaKa=piaKaS+0.0
     do k=bbb,bcf
        if(zKuC(k+1)>0) then
           ztrue=zKuC(k+1)
           ftran=1
           !dni=1.5*(dnCoeff(1)*ztrue+dnCoeff(2))+0.2*log(epst)
           dn=1.0*(dnCoeff(1)*ztrue+dnCoeff(2))+dnp(k+1)+0.2*log(epst)+dnst
           n1=(ztrue-1*dn-zmin)/dzbin
           dm=d013Table(n1,imu)
           !------------------------------------------------
           !dm=exp(dmCoeff(1)*ztrue**2+dmCoeff(2)*ztrue+dmCoeff(3))
           !if(dm.gt.0.5.and.dm.lt.2.8) then
           !if (it.eq.2) then
           !   rn=normal2(0.0,1.0)
           !   dm=0.875*dm+0.125*(1.6+0.25*normal2(0.0,1.0))
           !end if
           !end if
           !if(dm<0.31) dm=0.31
           !if(dm>3) dm=3
           !call bisection2(d013Table(:,imu),nbins,dm, n1)
           !dn=(ztrue-zKuSJ(n1))/10
           !dn=0.5*(dn+dn)
           !---------------------------------------------------
           if(10*dn+ztrue.gt.zKuSJ(nbins)) dn=(zKuSJ(nbins)-ztrue)/10.0
           if(10*dn+ztrue.lt.zKuSJ(1)) dn=(zKuSJ(1)-ztrue)/10.0
           n1=((ztrue-10*dn-zmin)/dzbin)+1
           if(n1.lt.1) n1=1
           if(n1.gt.nbins) n1=nbins
           dm=d013Table(n1,imu)
           attKu=att13Table(n1,imu)*10**dn
           attKa=att35Table(n1,imu)*10**dn
           zKa1=z35Table(n1,imu)+10*dn
           if(it.eq.3) then
              if(n1.lt.nbins-4) then
                 zKa2=z35Table(n1+4,imu)+10*(dn-0.1)
                 attKa2=att35Table(n1+4,imu)*10**(dn-0.1)
                 dzdn(k+1,k+1)=dzdn(k+1,k+1)+(zka2-zka1)/(-0.1)-(attKa2-attKa)/(-0.1)*dr
                 do kk=k+1,bcf
                    dzdn(kk+1,k+1)=dzdn(kk+1,k+1)-(attKa2-attKa)/(-0.1)*dr
                 end do
              end if
           end if
           zetaS=zetaS+attKu/10**(zKuC(k+1)*0.1*beta)*&
                10**(zKuL(k+1)*0.1*beta)*dr
           rrate1d(k+1)=pr13Table(n1,imu)*10**dn
           dn1d(k+1)=dn
           dm1d(k+1)=dm
           piaKa=piaKa+attKa*dr
           zKaSim(k+1)=z35Table(n1,imu)+10*dn-piaKa
           piaKa=piaKa+attKa*dr
        else
           attKu=0.0
           attKa=0.0
        end if
        zeta1d(k+1)=zetaS
     end do
     eps=min(1.0,(1-10**(-0.1*beta*f*piamax))/(q*beta*zetaS))
     epst=epst*eps
     do k=bzd,bcf
        zKuC(k+1)=zKuL(k+1)+piaKuS-10/beta*log10(1-eps*q*beta*zeta1d(k+1))
     end do
     piaKu=(-10/beta*log10(1-eps*q*beta*zeta1d(bcf+1)))+attKu*(bsfc-bcf)*2*dr
     piaKa=piaKa+attKa*(bsfc-bcf)*2*dr
     if(isnan(piaKu).eqv..true.) then
        print*, eps, zetaS, piamax, zeta1d(bcf+1)
        stop
     end if
     if(isnan(piaKa).eqv..true.) then
        print*, 'st',eps, zetaS, piamax, zeta1d(bcf+1), piaKa, piaKas,attKa, bsfc
        stop
     end if
  end do
end subroutine iter_profst

