      subroutine runge2(nn,tmin,tmax,xini,vini)
c written by Yoshitaro Takaesu -Aug/11/2011 @KEK
      implicit none

      integer i

      integer nn
      real*8 tmin,tmax,xini,vini

      real*8 x,v,t,h,rk1,rl1,rk2,rl2,rk3,rl3,rk4,rl4
      real*8 df,ddf
      external ddf,df

      h = (tmax -tmin) / nn
      t = tmin
      x = xini
      v = vini

      open(1, file='runge2.dat') 
      write(1, *) t, x
      do i=1, nn
         rk1 = h * df(x,v,t)
         rl1 = h * ddf(x, v, t)
         rk2 = h * df(x + rk1/2., v + rl1/2., t + h/2.)
         rl2 = h * ddf(x + rk1/2., v + rl1/2., t + h/2.)
         rk3 = h * df(x + rk2/2., v + rl2/2., t + h/2.)
         rl3 = h * ddf(x + rk2/2., v + rl2/2., t + h/2.)
         rk4 = h * df(x + rk3, v + rl3, t + h)
         rl4 = h * ddf(x + rk3, v + rl3, t + h)

         x = x + (rk1 + 2.*rk2 + 2.*rk3 + rk4)/6. 
         v = v + (rl1 + 2.*rl2 + 2.*rl3 + rl4)/6. 
         t = t + h

         write(1,*) t, x
      enddo

      close(1)
      
      return
      end


      real*8 function df(x,v,t)
      implicit none

      real*8 x,v,t

      df = v
      
      return
      end
