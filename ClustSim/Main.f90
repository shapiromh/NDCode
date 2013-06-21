module mstype

    use nrtype
    ! Program-specific parameters
    real(DP), parameter :: x = 500 ! Production value
    real(DP), parameter :: lambda = .3_dp ! Density benefit, must be between 0,1
    integer, parameter :: J = 2 ! Number of firms
    integer, parameter :: N = 50 ! Sq. root of number of sections
    ! TODO Make more heterogenous so there are many larger and smaller firms
    integer, parameter :: parmin = 1 ! Minimum support of Pareto distribution, >1 for finite mean
    real(DP), parameter :: alpha = 2.5_dp ! Shape parameter of Pareto distribution, >2 for finite var

end module mstype

module firms
    
    use mstype
    use nrtype
    use nrutil, ONLY : cumsum
	implicit none
    
    ! Declare variables for both subroutines
    real(DP) :: temp,pr,prcl,prcv
    integer :: i,k,ND,firm
    integer :: NN = N*N
    real(DP), dimension(J) :: c
    
contains 
    subroutine baseline(cost,CVbl)

        real(DP), dimension(J), intent(in) :: cost
        real(DP), intent(out) :: CVbl
        real(DP) :: blah
        integer, dimension(N,N) :: grid
        real(DP), dimension(N,N) :: clus
        integer, dimension(J) :: dem, sections
        real(DP) :: aa=0.0_dp, bb=3.0_dp, cc=1000.0_dp
        PROCEDURE(template_function), POINTER :: func
        
        ! Find equilibrium price and demand in baseline model 
        c=cost
        func => exdemand
        blah = brent(aa,bb,cc,func,1D-10,pr)
        dem = NINT(demandbl(pr))
    
        ! Allocate firms on a grid 
        sections = dem
        do i = 1,N
            do k = 1,N
                firm = firmsec(sections)
                grid(i,k) = firm
                sections(firm) = sections(firm)-1 ! Subtract one for firm that just won a slot 
            end do 
        end do

    
        ! Calculate the clustering value for the baseline model
        do i = 1,N
            do k = 1,N
                clus(i,k)=cluster(grid,i,k)
            end do
        end do
    
        temp=sum(clus)
        CVbl=temp/sum(dem)
    
    end subroutine baseline
    
    subroutine density(cost,CVden)
        
        real(DP), dimension(J), intent(in) :: cost
        real(DP), intent(out) :: CVden
        real(DP), dimension(J,2) :: dem
        integer :: iter
        real(DP) :: toler
        real(DP), dimension(3,2) :: startprice
        real(DP), dimension(3) :: funceval
        PROCEDURE(template_function3), POINTER :: func2
        
        ! Find equilibrium objects in density model
        c=cost
        func2 => clexdemand
        toler = 1.0D-10
        startprice(1,:) = (/ 600.0_dp , 1000.0_dp /)
        startprice(2,:) = (/ 620.0_dp , 1020.0_dp /)
        startprice(3,:) = (/ 580.0_dp , 960.0_dp /)
        funceval(1) = clexdemand(startprice(1,:))
        funceval(2) = clexdemand(startprice(2,:))
        funceval(3) = clexdemand(startprice(3,:))
        
        call amoeba(startprice,funceval,toler,func2,iter)
        ! print*,startprice, iter
        prcl = startprice(1,1)
        prcv = startprice(1,2)
        dem = demandcl(prcl,prcv)
!         print* , dem(:,1), " ", sum(dem(:,1))-NN
!         print*, dem(:,2)
                    
        ! Calculate the expected clustering value for the density model
        CVden=sum(dem(:,2))*(1/real(J))
        
        ! Allocate firms on a grid
        
       !  sections = dem
!         do i = 1,N
!             do k = 1,N
!                 firm = firmsec(sections)
!                 grid(i,k) = firm
!                 sections(firm) = sections(firm)-1 ! Subtract one for firm that just won a slot 
!             end do 
!         end do
!         
!         ! 
        
        
    end subroutine density
    
    SUBROUTINE amoeba(p,y,ftol,func,iter)
        USE nrtype; USE nrutil, ONLY : assert_eq,imaxloc,iminloc,nrerror,swap
        IMPLICIT NONE
        INTEGER(I4B), INTENT(OUT) :: iter
        REAL(dp), INTENT(IN) :: ftol
        REAL(dp), DIMENSION(:), INTENT(INOUT) :: y   ! "func" evaluated at the n vertices provided in "p"
        REAL(dp), DIMENSION(:,:), INTENT(INOUT) :: p ! vertices. If we have n vertices, then we must be
                                                         ! in n-1 dimensional space (we need one extra vertex
                                                         ! than dimensions. For each row, the n-1 vector
                                                         ! specifies the vertex
        PROCEDURE(template_function3), POINTER, INTENT(in) :: func
        INTEGER(I4B), PARAMETER :: ITMAX=5000
        REAL(dp), PARAMETER :: TINY=1.0e-10
        INTEGER(I4B) :: ihi,ndim
        REAL(dp), DIMENSION(size(p,2)) :: psum
        call amoeba_private
    CONTAINS
        !BL
        SUBROUTINE amoeba_private
            IMPLICIT NONE
            INTEGER(I4B) :: i,ilo,inhi
            REAL(dp) :: rtol,ysave,ytry,ytmp

            ndim=assert_eq(size(p,2),size(p,1)-1,size(y)-1,'amoeba')
            iter=0
            psum(:)=sum(p(:,:),dim=1)
            do
                ilo=iminloc(y(:))
                ihi=imaxloc(y(:))
                ytmp=y(ihi)
                y(ihi)=y(ilo)
                inhi=imaxloc(y(:))
                y(ihi)=ytmp
                rtol=2.0_dp*abs(y(ihi)-y(ilo))/(abs(y(ihi))+abs(y(ilo))+TINY)
                if (rtol < ftol) then
                    call swap(y(1),y(ilo))
                    call swap(p(1,:),p(ilo,:))
                    RETURN
                end if
                if (iter >= ITMAX) call nrerror('ITMAX exceeded in amoeba')
                ytry=amotry(-1.0_dp)
                iter=iter+1
                if (ytry <= y(ilo)) then
                    ytry=amotry(2.0_dp)
                    iter=iter+1
                else if (ytry >= y(inhi)) then
                    ysave=y(ihi)
                    ytry=amotry(0.5_dp)
                    iter=iter+1
                    if (ytry >= ysave) then
                        p(:,:)=0.5_dp*(p(:,:)+spread(p(ilo,:),1,size(p,1)))
                        do i=1,ndim+1
                            if (i /= ilo) y(i)=func(p(i,:))
                        end do
                        iter=iter+ndim
                        psum(:)=sum(p(:,:),dim=1)
                    end if
                end if
            end do
        END SUBROUTINE amoeba_private
        !BL
        FUNCTION amotry(fac)
            IMPLICIT NONE
            REAL(dp), INTENT(IN) :: fac
            REAL(dp) :: amotry
            REAL(dp) :: fac1,fac2,ytry
            REAL(dp), DIMENSION(size(p,2)) :: ptry
            fac1=(1.0_dp-fac)/ndim
            fac2=fac1-fac
            ptry(:)=psum(:)*fac1-p(ihi,:)*fac2
            ytry=func(ptry)
            if (ytry < y(ihi)) then
                y(ihi)=ytry
                psum(:)=psum(:)-p(ihi,:)+ptry(:)
                p(ihi,:)=ptry(:)
            end if
            amotry=ytry
        END FUNCTION amotry
    END SUBROUTINE amoeba
   
    FUNCTION exdemand(price) RESULT(Y)
        ! Calculates excess demand in the baseline model
        
        real(DP), intent(in) :: price
        real(DP), dimension(J) :: Z
        real(DP) :: Y
    
        Z = demandbl(price)
        Y = abs(sum(Z)-NN)
        
    END FUNCTION exdemand
    
    function demandbl(price) RESULT(Z)
        ! Calculates demand
        
        real(DP), intent(in) :: price
        real(DP), dimension(J) :: Z
        
        Z = (x-price)/(2*c)
        
    end function demandbl
    
    FUNCTION firmsec(secs) RESULT(Y)
        ! Returns the firm that gets a section
        
        integer, dimension(J), intent(IN) :: secs
        real(DP), dimension(J) :: prob, cumprob
        real(DP) :: ran
        integer, dimension(1) :: Z
        integer :: Y
        
        if (sum(secs)==0) then
            Y = 0
            return
        end if 
        prob = real(secs) / sum(secs)
        cumprob = cumsum(prob)
        ran = rand()
        cumprob = cumprob-ran
        Z = minloc(cumprob,mask=(cumprob>0))
        Y = Z(1) ! Because fortran is stupid and minloc returns array
        
    END FUNCTION firmsec
    
    function cluster(gr,l,m) RESULT(Y)
        ! Creates cluster value for one section
        
        integer, dimension(N,N), intent(IN) :: gr
        real(DP) :: Y
        integer, intent(in) :: l,m
        integer :: tot, count,ii,kk,sr,er,sc,ec
    
        tot = (min(l-1,1)+1+min(N-l,1))*(min(m-1,1)+1+min(N-m,1))
        sr = l-min(l-1,1)
        er = l+min(N-l,1)
        sc = m-min(m-1,1)
        ec = m+min(N-m,1)
        
        count=0
        do ii= sr,er
            do kk = sc,ec 
                if (gr(l,m)==gr(ii,kk) .and. (gr(l,m) /= 0)) then
                    count=count+1
                end if 
            end do
        end do
        
        Y = real(count) / tot
        
    end function cluster
    
    FUNCTION brent(ax,bx,cx,func,tol,xmin)
            USE nrtype; USE nrutil, ONLY : nrerror
            IMPLICIT NONE
            REAL(dp), INTENT(IN) :: ax,bx,cx,tol
            REAL(dp), INTENT(OUT) :: xmin
            REAL(dp) :: brent
            PROCEDURE(template_function), POINTER :: func
            INTEGER(I4B), PARAMETER :: ITMAX=100
            REAL(dp), PARAMETER :: CGOLD=0.3819660_dp,ZEPS=1.0e-3_dp*epsilon(ax)
            INTEGER(I4B) :: iter
            REAL(dp) :: a,b,d,e,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm
            a=min(ax,cx)
            b=max(ax,cx)
            v=bx
            w=v
            x=v
            e=0.0
            fx=func(x)
            fv=fx
            fw=fx
            do iter=1,ITMAX
                xm=0.5_dp*(a+b)
                tol1=tol*abs(x)+ZEPS
                tol2=2.0_dp*tol1
                if (abs(x-xm) <= (tol2-0.5_dp*(b-a))) then
                    xmin=x
                    brent=fx
                    RETURN
                end if
                if (abs(e) > tol1) then
                    r=(x-w)*(fx-fv)
                    q=(x-v)*(fx-fw)
                    p=(x-v)*q-(x-w)*r
                    q=2.0_dp*(q-r)
                    if (q > 0.0) p=-p
                    q=abs(q)
                    etemp=e
                    e=d
                    if (abs(p) >= abs(0.5_dp*q*etemp) .or. &
                        p <= q*(a-x) .or. p >= q*(b-x)) then
                        e=merge(a-x,b-x, x >= xm )
                        d=CGOLD*e
                    else
                        d=p/q
                        u=x+d
                        if (u-a < tol2 .or. b-u < tol2) d=sign(tol1,xm-x)
                    end if
                else
                    e=merge(a-x,b-x, x >= xm )
                    d=CGOLD*e
                end if
                u=merge(x+d,x+sign(tol1,d), abs(d) >= tol1 )
                fu=func(u)
                if (fu <= fx) then
                    if (u >= x) then
                        a=x
                    else
                        b=x
                    end if
                    call shft(v,w,x,u)
                    call shft(fv,fw,fx,fu)
                else
                    if (u < x) then
                        a=u
                    else
                        b=u
                    end if
                    if (fu <= fw .or. w == x) then
                        v=w
                        fv=fw
                        w=u
                        fw=fu
                    else if (fu <= fv .or. v == x .or. v == w) then
                        v=u
                        fv=fu
                    end if
                end if
            end do
            call nrerror('brent: exceed maximum iterations')
        CONTAINS
            !BL
            SUBROUTINE shft(a,b,c,d)
                REAL(dp), INTENT(OUT) :: a
                REAL(dp), INTENT(INOUT) :: b,c
                REAL(dp), INTENT(IN) :: d
                a=b
                b=c
                c=d
            END SUBROUTINE shft
    END FUNCTION brent
    
    function demandcl(price,pricecv) result(z)
        ! Calculates demand for sections and clustering value in density model
        
        real(DP), intent(in) :: price
        real(DP), intent(in) :: pricecv
        real(DP), dimension(J) :: clus, demand
        real(DP), dimension(J,2) :: Z
        integer :: iii
        
        clus = 1.0/lambda - (x-price)/(2.0*sqrt(c*lambda*pricecv))
        
        do iii=1,J
            if (clus(iii)>1) then
                Z(iii,2)=1.0_dp
                Z(iii,1)=(x-price)/(2.0*c(iii)*(1.0-lambda))
            else if (clus(iii)<0) then
                Z(iii,2)=0.0_dp
                Z(iii,1)=(x-price)/(2.0*c(iii)*lambda)
            else 
                Z(iii,2)=clus(iii)
                Z(iii,1)=sqrt(pricecv/(c(iii)*lambda))
            end if
        end do
        
    end function demandcl
    
    function clexdemand(pricev) result(z)
        ! Calculates excess demand in density model
        
        real(DP), dimension(2), intent(in) :: pricev
        real(DP), dimension(J,2) :: Y
        real(DP), dimension(J) :: int
        real(DP) :: z, demz, cvz, price, pricecv
    
        price=pricev(1)
        pricecv=pricev(2)
        Y = demandcl(price,pricecv)
        int = Y(:,1)*Y(:,2)*((3.0_dp-(2.0_dp/N))**2) ! TODO This condition doesn't work
        demz = abs(sum(Y(:,1))-NN) 
        cvz = abs(sum(int)-(3*N-2)**2)
        
        z = demz+cvz
        
    end function clexdemand
    
    ! function clustercl
!         
!     end function clustercl
    
end module firms

program main
    
    ! TODO feed random seed
    
    use nrtype
    use mstype
    use firms
    implicit none
    
    real(DP), dimension(J) :: cost
    real(DP) :: temp2, CVden, CVbl
    integer :: jj,argcount,S ! Number of simulations 
    integer :: sim = 1
    character(LEN=15) :: arg1
    real(DP), dimension(:,:), allocatable :: output ! To save CVden and CVbl results for output

    
    argcount = Command_argument_count()
    if (argcount /= 1) then
        print *, "Error"
        return
    end if
    call Get_command_argument(1, arg1, S)
    read (arg1,*) S                         ! argument is of type S

    
    ! Write(*,*) "How many simulations?"
    ! Read(*,*) S Check this
    
    allocate(output(S,2)) ! determine size of allocate
    
    ! Firms draw their cost from an inverted Pareto distribution 
    do while (sim<=S)
        do jj = 1, J
            temp2=rand()
            cost(jj)=1/(parmin*(1-temp2)**(-1/alpha)-.7) ! Last real is shift in Pareto draw left
        end do
        call baseline(cost,CVbl)
        output(sim,1)=CVbl
        call density(cost,CVden)
        output(sim,2)=CVden
        sim=sim+1
    end do
    
!     open(1, file="Simulation")
!     do jj = 1,S
!         write(1,*) output(jj,:)
!     end do
!     close(1)
    do jj = 1,S
        print*, output(jj,:)
    end do 
    
    deallocate(output) ! to prevent increased memory usage
    
    
end program main

