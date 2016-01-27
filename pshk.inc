	DOUBLE PRECISION P(5,50),M(3,50)
        COMMON /PARTCL/P,M

        DOUBLE PRECISION KPT2(20),KETA(20),KPhi(20),KPT(20)
        COMMON /PARTCL2/ KPT2, KETA, KPhi,KPT

        REAL*8 PTcut,AYEcut,rcut,QCut,Pi,Qmax
        COMMON /cut/ PTcut,AYEcut,rcut,QCut,Pi,Qmax

        REAL*8 PTcut2,AYEcut2,AYEcut3,AYEcut4,ptmin_cut(20),ptmax_cut(20)
        COMMON /cut2/ ptmin_cut,ptmax_cut,PTcut2,AYEcut2,AYEcut3,AYEcut4

	real*8 qcut2,minmij,mij,lptcut,QQ,ptijcut
	common /cut3/ qcut2,minmij,mij,lptcut,QQ,ptijcut

	real*8 ss,sqrts,x1,x2
        common /kinvar/ ss,sqrts,x1,x2

        real*8 wpsn1,wpsn2,jacob1,jacob2,wpsn,jacob		
	common /jacob/ wpsn1,wpsn2,jacob1,jacob2,wpsn,jacob 