#include <Rcpp.h>
#include <Python.h>
	#####################  Librerias
from numpy import isnan,savetxt,genfromtxt, array, log, unique, exp, append,concatenate,zeros, repeat,linspace,matrix
import pytwalk
import cProfile
from scipy.stats import uniform
from matplotlib.pyplot import plot, close, show, savefig,hist, xlabel, ylabel, title,axis,subplot, figure, setp
from numpy.random import seed


def runmod(dirt,plomo,carbon,Dircc,T_mod,T_mod_C,S_year,num_sup,det_lim,iterations, by,shape1_m,mean_m,shape_acc,mean_acc,fi_mean,fi_acc,As_mean,As_acc,cc,ccpb,resolution,seeds):
	seed(int(seeds))
    ##################### Data
	shape2_m= (shape1_m*(1-mean_m) )/mean_m
	scale_acc=mean_acc/shape_acc
	Data=genfromtxt (dirt+plomo, delimiter = ',')
	data=genfromtxt (dirt+carbon, delimiter = ',')
	#Radiocarbon data should be orded as
	#Radiocarbon age, sd, depth
	if data.ndim ==1:
		data=array([data,data])
	fimean=fi_mean
	shapefi=fi_acc
	ASmaean=As_mean
	shapeAS=As_acc
	shape2_m= (shape1_m*(1-mean_m) )/mean_m
	scale_acc=mean_acc/shape_acc
	scale_fi=fimean/shapefi
	scale_As=ASmaean/shapeAS

	##################### Data definition 210Pb
	num_sup=len(Data[:,0])-num_sup
	density=Data[:num_sup,1] * 10.
	activity=Data[:num_sup,2]
	sd_act=Data[:num_sup,3]
	thic=Data[:num_sup,4]
	depth=Data[:num_sup,0]
	if num_sup == 0:
		supp=Data[num_sup:,5]
		sd_supp=Data[num_sup:,6]
	if num_sup != 0:
		supp=Data[num_sup:,2]
		sd_supp=Data[num_sup:,3]

	activity=activity*density
	sd_act=sd_act*density

	Sample_year=1950.-S_year

	##################### Constant definition
	lam=0.03114
	##################### Read calibration curve

	if cc==1:
		intcal = genfromtxt(Dircc+'IntCal13.14C', delimiter = '\t')
		#intcal = intcal[::-1,...]
		ic = intcal[:,0:3]
		print("IntCal13 will be use")
	if cc==2:
		intcal = genfromtxt(Dircc+'Marine13.14C', delimiter = '\t')
    #intcal = intcal[::-1,...]
		ic = intcal[:,0:3]
		print("Marine13 will be use")
	if cc==3:
		intcal = genfromtxt(Dircc+'SHCal13.14C', delimiter = ',')
		#intcal = intcal[::-1,...]
		ic = intcal#[:,0:3]
		print("SHCal13 will be use")
	if ccpb==0:
		npost=0
		ic=ic
	if ccpb==1:
		intcalpost = genfromtxt(Dircc+'postbomb_NH1.14C', delimiter =  "\t"  )
		ic=concatenate((intcalpost,ic), axis=0)
		npost=len(intcalpost[:,1])
	if ccpb==2:
		intcalpost = genfromtxt(Dircc+'postbomb_NH2.14C', delimiter =  "\t"  )
		ic=concatenate((intcalpost,ic), axis=0)
		npost=len(intcalpost[:,1])
	if ccpb==3:
		intcalpost = genfromtxt(Dircc+'postbomb_NH3.14C', delimiter =  "\t"  )
		ic=concatenate((intcalpost,ic), axis=0)
		npost=len(intcalpost[:,1])
	if ccpb==4:
		intcalpost = genfromtxt(Dircc+'postbomb_SH3.14C', delimiter =  "\t"  )
		ic=concatenate((intcalpost,ic), axis=0)
		npost=len(intcalpost[:,1])
	if ccpb==5:
		intcalpost = genfromtxt(Dircc+'postbomb_SH1-2.14C', delimiter =  "\t"  )
		ic=concatenate((intcalpost,ic), axis=0)
		npost=len(intcalpost[:,1])
	##################### Depths calculation
	#print(ic)
	#print(npost)


	hard_lim=invlookup(max(data[:,0]),max(data[:,1]),cc,ccpb,Dircc)[-1]+150




	m=1
	breaks=array(m*by)
	maxd=max(max(depth),max(data[:,2]))

	while m*by< maxd:
	    m += 1
	    breaks=append(breaks,m*by)

	m_last_Pb=0
	while  breaks[m_last_Pb]< depth[-1]:
		       m_last_Pb += 1

	dep_time_data=append(depth-thic,depth)
	dep_time_data=list(set(dep_time_data))
	X1, X0= [], []
	for i1 in range(len(depth)):
		for k1 in range(len(dep_time_data)):
			if depth[i1]== dep_time_data[k1]:
				X1=append(X1,int(k1))
			if (depth-thic)[i1]== dep_time_data[k1]:
				X0=append(X0,int(k1))


	#################### Functions

	def support(param):
	    tmp3=True
	    for i in param:
             if i <= 0.:
		       tmp3=False
	    if param[2]>=1.:
		   tmp3=False
	    if times([depth[-1]],param)[-1]>last_t(param[0]):
		   tmp3=False
	    if times([breaks[-1]],param)>hard_lim:
		   tmp3=False
	    return tmp3

	def last_t(fi):
	    return ( (1/lam)*log(fi/det_lim) )


	def times(x,param):
	    w=param[2]
	    a=param[3:]
	    t1=m-1
	    ms1=array([a[m-1]])
	    while t1 >0 :
		   ms1= append(ms1, w*ms1[-1]+(1-w)*a[t1-1])
		   t1 -= 1

	    ms=ms1[::-1]
	    ages=array([])
	    y_last=append([0],array([sum(ms[:i+1]*by) for i in range(len(ms))] ) )
	    for i in range(len(x)):
		   k1=0
		   while  breaks[k1]< x[i]:
		       k1 += 1
		   ages=append(ages,y_last[k1]+(ms[k1]*(by-(breaks[k1]-x[i]))))
	    return ages

	def pendi(param):
		w=param[2]
		a=param[3:]
		t1=m-1
		ms1=array([a[m-1]])
		while t1 >0 :
			ms1= append(ms1, w*ms1[-1]+(1-w)*a[t1-1])
			t1 -= 1
			ms=ms1[::-1]
		return ms

	def incallookup(points):
	    result =[]
	    for i in points:
		   if i<0:
		       down=1
		       if i > -61.2:
		           while i > ic[down][0]:
		               down += 1
		       else:
		           down=1
		       up = down-1
		   elif i < 14000:
		       down = int(i/5.+1)+npost
		       up = down-1
		   elif i < 25000:
		       down = (2791+int((i-14000.)/10))+npost
		       up = down -1
		   else:
		       down = +(3891+int((i-25000)/20)) +npost
		       up=down -1
		   icdown0=ic[down,0]  #CalBp
		   icdown1=ic[down,1]  #CalRC
		   icdown2=ic[down,2]  #SD
		   prop = (i - icdown0)/(ic[up,0]-icdown0)
		   mean = prop*(ic[up,1]-icdown1)+icdown1
		   var = prop*(ic[up,2]-icdown2)+icdown2
		   result.append ([mean,var])
	    return result

#######lead likelihood

	def ln_like_data(param):
		Asup=param[1]*density
		loglike = 0.
		tmp2=param[0]/lam
		ts=times(dep_time_data,param)
		for i in range(len(activity)):
			A_i= Asup[i] + tmp2 *(exp(-lam*ts[int(X0[i])] ) - exp(-lam*ts[int(X1[i])]) )
			Tau=.5*(sd_act[i]**(-2.))
			loglike = loglike + Tau*((A_i-activity[i])**2.)
		return loglike

	def ln_like_T(param):
		loglike = 0.
		Asup=param[1]*density
		tmp2=param[0]/lam
		ts=times(dep_time_data,param)
		for i in range(len(activity)):
			A_i= Asup[i] + tmp2 *(exp(-lam*ts[int(X0[i])] ) - exp(-lam*ts[int(X1[i])]) )
			Tau=.5*(sd_act[i]**(-2.))
			loglike = loglike + 3.5*log(4. + Tau*((A_i-activity[i])**2.) )
		return loglike


	def ln_like_supp(param):
	    logsupp=0.
	    for i in range(len(supp)):
		   logsupp = logsupp + ((param[1]-supp[i])**2.)/(2.*(sd_supp[i])**2.)
	    return logsupp


	def ln_prior_supp(param):
		prior=0.
		prior= prior -  (  (shapefi-1.)*log(param[0])-(param[0]/scale_fi) )# prior for fi
		prior= prior -  (  ( (shapeAS-1.)-1.)*log(param[1])-(param[1]/scale_As) )# prior for supp
		prior= prior -  ( ((1./by)-1.)*log(param[2]) +  ((1./by)*(shape1_m-1.))*log(param[2]) + (shape2_m - 1.)*log(1.-param[2]**(1./by) ) )
		for ms in range(m):
		   prior= prior -  (  (shape_acc-1.)*log(param[ms+3])-(param[ms+3]/scale_acc) )
		return prior

	if T_mod:
	    log_data=ln_like_T
	else:
	    log_data=ln_like_data


	def Ux(param):
	    u=0
	    dat=times(data[:,2],param) + Sample_year
	    inc=incallookup(dat)
	    for i in range(len(inc)):
		   mu=inc[i][0]
		   sigm=(inc[i][1]**2+data[i][1]**2)
		   u=u+(7./2.)*log(4. + ((data[i,0]-mu)**2.)/((2.*sigm)) ) + .5*log(sigm)
		   #log(lamda/alpha1)/2-((alpha1+1)/2)*log(1+(lamda/alpha1)*(data[i,0]-mu)**2)
	    return u

	def UxN(param):
	    u=0
	    dat=times(data[:,2],param) + Sample_year
	    inc=incallookup(dat)
	    for i in range(len(inc)):
		   mu=inc[i][0]
		   sigm=(inc[i][1]**2+data[i][1]**2)
		   u=u+( ((data[i,0]-mu)**2.)/((2.*sigm)) ) + .5*log(sigm)
		   #log(lamda/alpha1)/2-((alpha1+1)/2)*log(1+(lamda/alpha1)*(data[i,0]-mu)**2)
	    return u

	if T_mod_C:
	    log_dataC=Ux
	else:
	    log_dataC=UxN


	def obj(param):
	    objval= ln_like_supp(param) + ln_prior_supp(param) + log_dataC(param)  + log_data(param)
	    return objval

	def Calib(numfi,maxe):
	    figure(numfi)
	    for c in data:
		   fecha=c[:2]
		   if fecha[0]<=0.:
			x=linspace(Sample_year,200,1500)
		   else:
			lowa=max((fecha[0]-400.),Sample_year)
			x=linspace(lowa,maxe,2500)
		   u , newx=[], []
		   dat=array(x)
		   inc=incallookup(dat)
		   for i in range(len(inc)):
		       mu=inc[i][0]
		       sigm=(inc[i][1]**2 + fecha[1]**2)
		       utest=((4. + ((fecha[0]-mu)**2.)/((2.*sigm)) )**(-7./2)  )
		       u.append([utest]  )
		   u=array(u)
		   maxu=max(u)
		   u=(u/maxu)
		   y=[]
		   for i in range(len(u)):
		       if u[i]>.08:
		           y.append(u[i])
		           newx.append([x[i]])
		   y=array(y)*5.
		   newx=array(newx)- Sample_year
		   plot((-y+c[2]),newx,color="blue",alpha=.6)
		   plot( (y+c[2]),newx,color="blue", alpha=.6)

	#################### Initial valules
	print("Seaching initial values")
	fi_ini_1= uniform.rvs(size=1,loc=50, scale=200)  #200.
	fi_ini_2= uniform.rvs(size=1,loc=200, scale=150) #100.
	supp_ini_1= uniform.rvs(size=1,loc=5, scale=20) #5.
	supp_ini_2= uniform.rvs(size=1,loc=0, scale=20) #20.
	w_ini = uniform.rvs(size=1) #.3
	w_ini0 = uniform.rvs(size=1)  #.7
	m_ini_1=uniform.rvs(size=m,loc=0, scale=5)  #  repeat(array(3.1),m,axis=0)
	m_ini_2=uniform.rvs(size=m,loc=0, scale=3)  # repeat(array(.5),m,axis=0)

	x=append(append(append(fi_ini_1,supp_ini_1),w_ini), m_ini_1)
	xp=append(append(append(fi_ini_2,supp_ini_2),w_ini0), m_ini_2)

	while not support(x):
	    m_ini_1=uniform.rvs(size=m,loc=0, scale=1)
	    x=append(append(append(fi_ini_1,supp_ini_1),w_ini), m_ini_1)


	while not support(xp):
	    m_ini_2=uniform.rvs(size=m,loc=0, scale=1)
	    xp=append(append(append(fi_ini_2,supp_ini_2),w_ini0), m_ini_2)

	print("initial values were obtained")



	################## New MCMC test
	print("the number of itrations,")
	print(iterations)
	thi = int((len(x)))*50 #100
	print("Thining,")
	print(thi)
	burnin=10000*len(xp) #20000
	print("Burnin,")
	print(burnin)
	print("Total iterations,")
	print(burnin + iterations*thi)


	leadchrono = pytwalk.pytwalk(n=len(x),U=obj,Supp=support)
	i, k ,k0, n=0 , 0, 0, len(x)
	U , Up = obj(x), obj(xp)
	por=int(iterations/10.)
	Output = zeros((iterations+1, n+1))
	Output[ 0, 0:n] = x.copy()
	Output[ 0, n] = U
	por2=int(burnin/5.)
	while i< iterations:
		onemove=leadchrono.onemove(x, U, xp, Up)
		k+= 1
		if (all([k<burnin,k % por2==0]) ):
			print("burn in progress")
			print int(100*(k+.0)/burnin)
		if (uniform.rvs() < onemove[3] ):
			x, xp, ke, A, U, Up =onemove
			k0+=1
			if all([k % thi ==0 , k>int(burnin)]):
				Output[i+1,0:n] = x.copy()
				Output[i+1,n] = U
				if any([i % por==0, i==0]) :
					print int(100*(i+.0)/iterations),"%"
				   #print((time.clock()-tiempomedir)/60)
				i+= 1
		else:
			if all([k % thi ==0 , k>int(burnin)]):
				Output[i+1,0:n] = x.copy()
				Output[i+1,n] = U
				if any([i % por==0, i==0]) :
					print int(100*(i+.0)/iterations),"% Done"
				   #print((time.clock()-tiempomedir)/60)
				i+= 1



	#Output=array(Output)
	print("Acceptance rate")
	print(k0/(i+.0))
	print("The twalk did", k, "iterations")


	savetxt(dirt+'Results/Results_output.csv', Output,delimiter=',')
	estim=[]
	for i in range(iterations):
	    estim.append(times(breaks,Output[i+1,:-1])  )
	estim=array(estim)
	savetxt(dirt+'Results/dates.csv', estim  )
	intervals=[]
	for i in range(len(estim[1,])):
	    sort=sorted(estim[:,(i)])
	    mean=sum(sort)/len(sort)
	    disc=int(len(sort)*.025)
	    sort=sort[disc:]
	    sort=sort[:-disc]
	    intervals.append([breaks[i],sort[0],mean,sort[-1]])
	savetxt(dirt+'Results/intervals.csv', intervals ,delimiter=',')
	depths=array([append([0.0],breaks)])
	savetxt(dirt+"Results/depths.csv", depths,delimiter=',')

	grafdepts=linspace(0,breaks[-1],resolution)
	grafdepts2=grafdepts+(grafdepts[1]-grafdepts[0])/2
	grafdepts2=grafdepts2[0:(len(grafdepts2)-1)]

	grafage=linspace(0,(max(estim[:,-1])+.10),resolution)
	y=[]
	for i in range(len(depths[0,:])-1):
		logvect=array(grafdepts2>depths[0,i])*array(grafdepts2<=depths[0,i+1])
		for k in range(len(logvect)):
			if logvect[k]:
				if i!=0:
					y1=estim[:,i-1]+((estim[:,i]-estim[:,i-1])/(depths[0,i+1]-depths[0,i]))*(grafdepts[k]-depths[0,i])
					porc=[]
					for posi in range(len(grafage)-1):
						porc.append(sum(array(y1>=grafage[posi])*array(y1<grafage[posi+1]) ))
					y.append(porc/(max(porc)+0.0 ))
				else:
					y1=((estim[:,i]/depths[0,i+1])*(grafdepts[k]) )
					porc=[]
					for posi in range(len(grafage)-1):
						porc.append(sum(array(y1>=grafage[posi])*array(y1<grafage[posi+1]) ))
					y.append(porc/(max(porc)+0.0 ))

	savetxt(dirt+"Results/Graphs.csv", array(y),delimiter=',')
	slopes=[]
	for i in range(iterations-1):
		slopes.append(pendi(Output[(i+1),:-1])  )
	savetxt(dirt+"Results/Slopes.csv", array(slopes),delimiter=',')




#single lookup





def invlookup(date,sd,cc,ccpb,dirt):

    Xlow=date-5*sd
    Xhigh=date+5*sd
    if cc==1:
        intcal = genfromtxt(dirt+'IntCal13.14C', delimiter = '\t')
        intcal = intcal[::-1,...]
        ic = intcal[:,0:3]
    if cc==2:
        intcal = genfromtxt(dirt+'Marine13.14C', delimiter = '\t')
        intcal = intcal[::-1,...]
        ic = intcal[:,0:3]
    if cc==3:
        intcal = genfromtxt(dirt+'SHCal13.14C', delimiter = ',')
        intcal = intcal[::-1,...]
        ic = intcal[:,0:3]
    if ccpb==0:
        ic2=ic
        if(Xlow<min(ic2[:,1])):
            Xlow=min(ic2[:,1])+.1
    if ccpb==1:
        intcalpost = genfromtxt(dirt+'postbomb_NH1.14C', delimiter =  "\t"  )
        #intcalpost = intcalpost[::-1,...]
        ic2=concatenate((intcalpost,ic), axis=0)

    if ccpb==2:
        intcalpost = genfromtxt(dirt+'postbomb_NH2.14C', delimiter =  "\t"  )
        #intcalpost = intcalpost[::-1,...]
        ic2=concatenate((intcalpost,ic), axis=0)

    if ccpb==3:
        intcalpost = genfromtxt(dirt+'postbomb_NH3.14C', delimiter =  "\t"  )
        #intcalpost = intcalpost[::-1,...]
        ic2=concatenate((intcalpost,ic), axis=0)

    if ccpb==4:
        intcalpost = genfromtxt(dirt+'postbomb_SH3.14C', delimiter =  "\t"  )
        #intcalpost = intcalpost[::-1,...]
        ic2=concatenate((intcalpost,ic), axis=0)

    if ccpb==5:
        intcalpost = genfromtxt(dirt+'postbomb_SH1-2.14C', delimiter =  "\t"  )
        #intcalpost = intcalpost[::-1,...]
        ic2=concatenate((ic,intcalpost), axis=0)



    placeslow=[]
    placeshigh=[]



    for i in range(len(ic2)-1):
        if all([ic2[i,1]>Xhigh,ic2[i+1,1]<=Xhigh]):
            placeshigh.append(i)
#            print("PLACESHIGH")
        if all([ic2[i,1]>Xlow,ic2[i+1,1]<=Xlow]):
            placeslow.append(i)
#            print("PLACES LOW")
    places=[ ic2[placeslow[-1],0], ic2[placeshigh[0],0] ]
    return places



def incallookup2(i,cc,ccpb,dirt):
    if cc==1:
        intcal = genfromtxt(dirt+'Calibration Curves/IntCal13.14C', delimiter = '\t')
        #intcal = intcal[::-1,...]
        ic = intcal[:,0:3]
    if cc==2:
        intcal = genfromtxt(dirt+'Calibration Curves/Marine13.14C', delimiter = '\t')
        #intcal = intcal[::-1,...]
        ic = intcal[:,0:3]
    if cc==3:
        intcal = genfromtxt(dirt+'Calibration Curves/SHCal13.14C', delimiter = ',')
        #intcal = intcal[::-1,...]
        ic = intcal[:,0:3]
    if ccpb==0:
        npost=0
        ic2=ic
    if ccpb==1:
        intcalpost = genfromtxt(dirt+'Calibration Curves/postbomb_NH1.14C', delimiter =  "\t"  )
        intcalpost = intcalpost[::-1,...]
        ic2=concatenate((ic,intcalpost), axis=0)
        npost=len(intcalpost[:,1])
    if ccpb==2:
        intcalpost = genfromtxt(dirt+'Calibration Curves/postbomb_NH2.14C', delimiter =  "\t"  )
        intcalpost = intcalpost[::-1,...]
        ic2=concatenate((ic,intcalpost), axis=0)
        npost=len(intcalpost[:,1])
    if ccpb==3:
        intcalpost = genfromtxt(dirt+'Calibration Curves/postbomb_NH3.14C', delimiter =  "\t"  )
        intcalpost = intcalpost[::-1,...]
        ic2=concatenate((ic,intcalpost), axis=0)
        npost=len(intcalpost[:,1])
    if ccpb==4:
        intcalpost = genfromtxt(dirt+'Calibration Curves/postbomb_SH3.14C', delimiter =  "\t"  )
        intcalpost = intcalpost[::-1,...]
        ic2=concatenate((ic,intcalpost), axis=0)
        npost=len(intcalpost[:,1])
    if ccpb==5:
        intcalpost = genfromtxt(dirt+'Calibration Curves/postbomb_SH1-2.14C', delimiter =  "\t"  )
        intcalpost = intcalpost[::-1,...]
        ic2=concatenate((ic,intcalpost), axis=0)
        npost=len(intcalpost[:,1])

    if i<0:
		down=1
		if i > -61.2:
			while i > ic2[down][0]:
				down += 1
		else:
			down=1
		up = down-1
    elif i < 14000:
		down = int(i/5.+1.)+npost
		up = down-1
    elif i < 25000:
		down = (2791+int((i-14000)/10))+npost
		up = down -1
    elif i < 49980:
		down = (3891+int((i-25000)/20))+npost
		up=down -1
    else:
		down=5789
		up=5790
    icdown0=ic2[down,0]
    icdown1=ic2[down,1]
    icdown2=ic2[down,2]
    prop = (i - icdown0)/(ic2[up,0]-icdown0)
    mean = prop*(ic2[up,1]-icdown1) + icdown1
    var = prop*(ic2[up,2]-icdown2)  + icdown2
    result=([mean,var])
    return result





#dirt="/home/endymion/Documents/PLUM-14C/Test1/"
#lead="Data-14C.csv"
#carbon="Data-14C.csv"
#genfromtxt (dir1+lead, delimiter = ',')
#dircc="/home/aquinom/R/x86_64-pc-linux-gnu-library/3.4/CPb/"
#incallookup([5])



#runmod("/home/aquinom/Documents/PLUM-14C/Test/","Data-210.csv","Data-14C.csv","/home/aquinom/R/x86_64-pc-linux-gnu-library/3.4/CPb/Calibration Curves/",False,True,2011,4,.01,100, 20,2,.7,2,20,1,1,200)

#dirt,plomo,carbon,Dircc,T_mod,T_mod_C,S_year,num_sup,det_lim,iterations, by,shape1_m,mean_m,shape_acc,mean_acc,cc,ccpb,resolution = "/home/aquinom/Documents/PLUM-14C/Test/","Data-210.csv","Data-14C.csv","/home/aquinom/R/x86_64-pc-linux-gnu-library/3.4/CPb/Calibration Curves/",False,True,2011,4,.01,100, 20,2,.7,2,20,1,1,200


