       parameter (nc=10)
       dimension xc(0:nc)
       read(5,*) x0,R,h
c logistic dx/dt=-x+R*x^2 with Euler fwd marching
c carle0
       t0 = 0.0
       nt = 10./h
       xc(0) = x0
       write(10,*) t0,x0
       do it=1,nt
        t = it*h
        exa= x0*exp(-t)/(1.0-R*x0*(1.0-exp(-t)))
        xc(0) = (1.0-h)*xc(0)
        write(10,*) t,xc(0),exa
       end do
c carle1
       xc(1)=x0
       x2 = x0*x0
       write(11,*) t0,x0
       do it=1,nt
        t = it*h
        exa= x0*exp(-t)/(1.0-R*x0*(1.0-exp(-t)))
        x2 = (1.0-2.0*h)*x2
        xc(1) = (1.0-h)*xc(1)+h*R*x2
        write(11,*) h*it,xc(1),exa
       end do
c carle2
       xc(1)= x0
       xc(2)= x0*x0
       x3   = x0*x0*x0
       write(12,*) t0,x0
       do it=1,nt
        t = it*h
        exa= x0*exp(-t)/(1.0-R*x0*(1.0-exp(-t)))
        x3 = (1.0-3.0*h)*x3
        xc(2) = (1.0-2.0*h)*xc(2)+2.0*h*R*x3
        xc(1) = (1.0-h)*xc(1)    +h*R*xc(2)
        write(12,*) t,xc(1),exa
       end do
c carle3
       xc(1)= x0
       xc(2)= x0*x0
       xc(3)= x0*x0*x0
       x4   = x0*x0*x0*x0
       write(13,*) t0,x0
       do it=1,nt
        t = it*h
        exa= x0*exp(-t)/(1.0-R*x0*(1.0-exp(-t)))

        x4 = (1.0-4.0*h)*x4
        xc(3) = (1.0-3.0*h)*xc(3)+3.0*h*R*x4
        xc(2) = (1.0-2.0*h)*xc(2)+2.0*h*R*xc(3)
        xc(1) = (1.0-h)*xc(1)    +h*R*xc(2)
        write(13,*) t,xc(1),exa
       end do

       stop
       end
