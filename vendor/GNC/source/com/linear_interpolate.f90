subroutine linear_int(x,y,n, xev,yev)
	implicit none
	integer n,i
	real(8) x(n),y(n),xev,yev
    if(n.le.1)then
        print*, "linear_int, n<=1", n
        return
    end if
    if(xev<x(1)) then
!		print*, "ev out of range"
!		pause
        yev=y(1)-(y(2)-y(1))/(x(2)-x(1))*(x(1)-xev)	
        return
    end if
    if(xev>x(n))then
        yev=y(n)+(y(n)-y(n-1))/(x(n)-x(n-1))*(xev-x(n))	
        !print*, "n=",n
        !print*, "x=",x
        !print*, "y=",y
        return
    end if
    if(xev.eq.x(n))then
		yev=y(n)
		return
	end if
	do i=1, n-1
		if(x(i).ge.x(i+1)) then
			print*,"bad input in linear_int:x(i)>x(i+1)"
			print*, i, x(i), x(i+1)
			pause
		end if
		if(x(i).lt.xev.and.x(i+1).gt.xev)then
            yev=y(i+1)-(y(i+1)-y(i))/(x(i+1)-x(i))*(x(i+1)-xev)	
            return
		end if
		if(x(i).eq.xev)then
			yev=y(i)
			return
		end if
	end do
    print*, "error?, xev=",xev
    print*, "x=",x
    print*, "y=",y
    stop
end subroutine
subroutine linear_int_2d(xmin,ymin,nx,ny,xstep, ystep, z, vx,vy, vz)
	! must be grid type
	!   x=|=x=|=x=|=x  
	implicit none
	integer nx, ny
	real(8) xmin,ymin,z(nx,ny)
	real(8) vx, vy, vz
	real(8) rdx,rdy
	real(8) xstep, ystep
	integer idx,idy,idxn,idyn,idxm,idym
	real(8) y1,y2,y3,y4,t,u
	!print*, "xmin,ymin,nx,ny,xstep, ystep, vx,vy, vz=", xmin,ymin,nx,ny,xstep, ystep, vx,vy, vz
	rdx=(vx-xmin)/xstep
	rdy=(vy-ymin)/ystep
	idx=int(rdx)+1
	idy=int(rdy)+1
	!print*, "rdx,rdy,idx,idy=",rdx,rdy,idx,idy

	if(idx<0.or.idx>nx.or.idy<0.or.idy>ny)then
		vz=0
		return
	end if

	if(idx.eq.nx)then
		idxn=idx-1
		idxm=idx
	elseif(idx.eq.1)then
		idxn=1
		idxm=2
	else
		idxn=idx
		idxm=idx+1
	end if
	if(idy.eq.ny)then
		idyn=idy-1
		idym=idy
	elseif(idy.eq.1)then
		idyn=1
		idym=2
	else
		idyn=idy
		idym=idy+1
	end if

	y1=z(idxn,idyn)
	y2=z(idxm,idyn)
	y3=z(idxm,idym)
	y4=z(idxn,idym)
	!print*, "y1,y2,y3,y4=",y1,y2,y3,y4
	t=rdx-idxn+1
	u=rdy-idyn+1
	!print*, "t,u=",t,u

	vz=(1-t)*(1-u)*y1+t*(1-u)*y2+t*u*y3+(1-t)*u*y4
	!print*, "vz=",vz
	!read(*,*)
end subroutine

subroutine linear_int_arb(x,y,n, xev,yev)
	implicit none
	integer n,i
	real(8) x(n),y(n),xev,yev
	real(8) x1,xn, y1, yn 
	logical increase
    if(n.le.1)then
        print*, "linear_int, n<=1", n
        return
    end if
	if(x(1)>x(n))then
		x1=x(n); xn=x(1)
		y1=y(n); yn=y(1)
		increase=.false.
	else
		x1=x(1); xn=x(n)
		y1=y(1); yn=y(n)
		increase=.true.
	end if

	if(xev<x1.or.xev>xn) then
		print*, "xev out of range ", xev, x1, xn
		stop
		return
	end if
	
    
	if(increase)then
		if(xev.eq.xn)then
			yev=yn
			return
		end if
		do i=1, n-1
			if(x(i).lt.xev.and.x(i+1).gt.xev)then
				yev=y(i+1)-(y(i+1)-y(i))/(x(i+1)-x(i))*(x(i+1)-xev)	
			!	print*, "yev,i=",yev,i
				return
			end if
			if(x(i).eq.xev)then
				yev=y(i)
				return
			end if
		end do
	else
		if(xev.eq.xn)then
			yev=yn
			return
		end if
		
		do i=n, 2, -1
			!print*, "i=",i,x(i),xev, x(i-1)
			if(x(i).lt.xev.and.x(i-1).gt.xev)then
				yev=y(i-1)-(y(i-1)-y(i))/(x(i-1)-x(i))*(x(i-1)-xev)	
			!	print*, "yev,i=",yev,i
				return
			end if
			if(x(i).eq.xev)then
				yev=y(i)
				return
			end if
		end do
	end if

    print*, "error?, xev=",xev
    print*, "x=",x
    print*, "y=",y
    stop
end subroutine