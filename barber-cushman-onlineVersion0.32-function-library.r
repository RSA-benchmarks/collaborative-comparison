#!R script
#=======================================================
# author: J.A.Postma
# contact: digigram at gmail dot com
# date: summer 2014
# license: gpl
# purpose: Simulating nutrient uptake by roots
# reference: Itoh S, Barber SA. 1983. A numerical solution of whole plant nutrient uptake for soil-root systems with root hairs. Plant and Soil 70: 403-413.
# content: functions use by the model. 
# usage: in R type
#     source("barber-cushman-onlineVersion0.1-function-library.r")
#     param=paramBarberBookP
#     m<-barber(param,msg="my test run")
#     plotbarber(param,m)
# disclaimer: This software comes without any warrenty, and may contain errors.
# changes: Out of cortesy, please submit any changes to sourceforge project page (https://sourceforge.net/projects/barber-cushman/) or e-mail to above address
#
#=======================================================


#parameter set used in roothair paper
paramRoothairPublication1983=list(
   Imax=1e-6,       #1e-6 umol/cm2/s
   Imaxh=1e-6,
   Km=3e-3,            #3e3 umol/cm3
   Cmin=3e-4,         #3e-4
   b= 167,             #unit less 167
   Cli=1.55e-2,        #1.55 umol/cm3
   De = 2.e-9,        #2e-9 cm2/s 
   r0=1e-2,           #root radius in cm 1e-2
   rh=5e-4,                #radius of the root hairs 
   r1=3e-1,             #radius of the depletion zone cm 3e-1
   lh=3e-2,             #length of the root hairs 3e-2
   Nh=1000,               #Nh which is number of hairs / cm root
   v0 = 1E-6,          #1e-6 cm/s
   E=0, #not used

   dr=0.005,               #delta r (which will be used to determine number of spatial steps)
   dt=0.05,             #timestep day
   case=2 ,            #outerboundary condition 
   roothairs=1,        #include roothairs
   endtime=24,         #end time in days
   tcfactor=86400     #time conversion factor s to day
)

paramNST3_KNST=list(
  #note adjust a=1 and c=0 in nst and k= 0.0001
  # theta=1,f=1
  #linear growth, start with 1 m, with small k
  Imax=2.4e-6,       
  Imaxh=2.4e-6,
  Km=4.0e-1,            
  Cmin=1e-4,        #nst has 0 but better to have small Cmin
  b= 15,            
  Cli=2e-1,        
  De = 1.5e-6/15, #nst expresses De differently, divide by b        
  r0=0.019,           
  rh=0.0005,               
  r1=0.33,            
  lh=0.0,            
  Nh=1000,              
  v0 = 2.7E-7,          
  E=0, #not used
  
  dr=0.01,               
  dt=0.1,
  case=2 ,            
  roothairs=0,        
  endtime=15,         
  tcfactor=86400,    
  unit="mmol"
)


paramNST3spinach=list(
#p l rh P spinache limear b wit root hairs  
#note adjust a=1 and c=0 in nst and k= 0.0001
# theta=1,f=1
#units are in pmols  
#results seem quite different, so not sure if this is right, but think NST might not do a good job.   
   Imax=5.85e-1,       
   Imaxh=5.85e-1,
   Km=4.0e2,            
   Cmin=1e2,        
   b= 1500,            
   Cli=1.4e3,        
   De = 8.9e-6/1500,        
   r0=0.0107,           
   rh=0.0005,               
   r1=0.2107,            
   lh=0.02, # roothairlength is diffuse in NST00
   Nh=1000,              
   v0 = 2.7E-7,          
   E=0, #not used
   dr=0.005,              
   dt=0.1,
   case=2,            
   roothairs=1,        
   endtime=30,         
   tcfactor=86400,    
   unit="pmol"
)


paramBarberBookP=list(
   v0 = 5E-7/5,          #5*10-7 cm/s = cm/day0.0432
   De = 2.3e-9,        #cm2/s 1.9872e-4
   b= 163,             #unit less 163
   Cli=13.6e-3,        #umol/mL 13.6e-3
   Imax=3.21E-7,       #6.43 nmol/m2/s=6.43 10-7 umol/cm2/s
   Km=5.45e-3,         #umol/mL 5.45e-3
   E=0, #1.626881e-08  #very important for stability. efflux cmin=0.2e-2 μmol/ml so E=(Imax*0.2e-2/(Km+0.2e-2))=1.726e-7
   Cmin=2.e-4,
   #or corrected for circumference E=1.626881e-08
   r0=0.05,           #root radius in cm 0.015
   r1=0.55,             #radius of the depletion zone cm
   k=60,               #number of spatial steps
   small=1e-2,         #accuracy constant for iterative method
   dt=0.1,             #timestep day
   case=2,             #outerboundary condition 
   endtime=30,         #end time in days
   tcfactor=86400,     #time conversion factor s to day
   #for including roothairs
   roothairs=1,        #include roothairs
   rh=5e-4,                #radius of the root hairs 
   Nh=1000,               #Nh which is number of hairs / cm root
   lh=0.2            #length of the root hairs
)


#solver for tridiagonal matrix
tridiagonalMatrixSolver <- function(a,b,c,f){
   n=length(a)
   w=a;v=a;z=a;y=a;
   #avoid zero's in w
   w[w==0]<-1e-6
   v[1]=c[1]/w[1]; z[1]=f[1]/w[1]
   for(i in (2:n)){
      w[i]=a[i]-(b[i]*v[i-1])
      v[i]=c[i]/w[i]
      z[i]=(f[i]-(b[i]*z[i-1]))/w[i]
   }
   y[n]=z[n]
   for (i in (n-1):1){
      y[i]=z[i]-(v[i]*y[i+1])
   }
   return(y)
}

#setting all constants
setConstants<-function(p){
   #do some calculations that are needed for other calculations down
   Imax=p$Imax*p$dt*86400
   if(length(names(p)[names(p)=="Imaxh"])==1) {
      Imaxh=p$Imaxh*p$dt*86400
   }else{
      Imaxh=Imax
   }       #1e-6 umol/cm2/s
   De=p$De*p$dt*86400
   v0=p$v0*p$dt*86400
   if(length(names(p)[names(p)=="small"])==1){
     small=p$small
   }else{
     small=p$Cli/10000
   }
   if(length(names(p)[names(p)=="dr"])==1){
      dr=p$dr
      k=floor((p$r1-p$r0)/p$dr)
      r1=p$r0+k*dr
      r=((0:(k-1))*dr)+p$r0
      dr = array(0,k)
      dr[1:(k-1)]=diff(r)
      dr[k]=dr[k-1]
   }else if(length(names(p)[names(p)=="k"])==1){   
      k=p$k
      dr=(p$r1-p$r0)/(p$k-1)
      r1=p$r1
      r=((0:(k-1))*dr)+p$r0
      dr = array(0,k)
      dr[1:(k-1)]=diff(r)
      dr[k]=dr[k-1]
   }else{
     #variable dr and k
     print("Warning: you are using experimental code which is known to produce bogus results, specify k or dr to have fixed spatial steps")
     lh=p$lh
     if(is.null(lh)) lh=0
     offset=lh/2 #at what distance dr starts to increase
     lhd=lh-offset
     if(lhd<0) lhd=0
     const_dr=p$initial_dr #probably this is better controled by De/V0
     scale=0.1 #relative increase of dr
     N1 = ceiling(lhd/const_dr);
     N2 = floor(log(1+scale*(p$r1-p$r0-lhd)/const_dr )/log(scale+1) )
     k = N1+N2;
     dr = array(const_dr,k)
     dr[(N1+1):k]= const_dr*cumprod(array(scale+1,k-N1))
     r=p$r0+c(0,cumsum(dr[2:k]))
     r[k]=p$r1;
     dr[1:(k-1)]=diff(r)
     dr[k]=dr[k-1] 
     r1=p$r1
   }
   S=(dr/2)*(1+((v0*p$r0)/(De*p$b)))
   rh1= sqrt(r*pi/(2*p$Nh)) #this assuming squared organization, random would be better, but the the model is not very sensitive for rh1, so it does not matter so much I think
   vol=(pi*(r+dr)^2)-(pi*(r)^2) #const$vol=pi*((r+0.5*dr)^2-(r-0.5*dr)^2)
   #vol=(pi*(r+dr/2)^2)-(pi*(r-dr/2)^2)
   Ah=(p$Nh*pmin(dr,pmax(0,(p$lh+p$r0)-r))*2*pi*p$rh)/vol   #   #roothairs only, surface area of the root hairs in each compartment cm2/cm3
   Ah[r>(p$lh+p$r0)]<-0
   #simply for testing if Ah is not in the boundary compartment
   #Ah[2]=Ah[1]*vol[1]/vol[2]
   #Ah[1]=0
   #unit
   if(is.character(p$unit)){
     unit=p$unit
   }else{
     unit="μmol"
   }
   #construct list of parameters
   D3=dr^2/(De*p$b)
   D3[1]=D3[1]*2
   results=list(
         Imax=Imax,
         Imaxh=Imaxh,
         Km=p$Km,
         Cmin=p$Cmin,
         b=p$b,
         Cli=p$Cli,
         De=De,
         r0=p$r0,
         rh=p$rh,
         r1=r1,             
         dr=dr, 
         r=r,
         vol=vol,
         lh=p$lh,
         Nh=p$Nh,
         v0=v0,
         E=p$E, 
         k=k,               #number of spatial steps
         small=small,     #accuracy constant for iterative method
         dt=p$dt,             #timestep day
         endtime=p$endtime,         #end time in days
         time=1:(p$endtime/p$dt),

         case=p$case ,            #outerboundary condition 
         roothairs=p$roothairs,        #include roothairs
         D1=(2*(dr^2)/De)+2,
         D2=(2*(dr^2)/De)-2,
         D3=D3, #for roothairs only
         S1=(S/r)-1,
         S2=(S/r)+1,
         S3=((2*dr[1])/(De*p$b)),
         A1=(2*dr[k]*p$r0*v0)/(De*p$b*r1),
         rh1= rh1, #this assuming squared organization, random would be better, but the the model is not very sensitive for rh1, so it does not matter so much I think
         S4=((Imaxh*p$rh)/(De*p$b))*log(rh1/(1.6487*p$rh)),#roothairs only, log computes natural logarithm by default (ln)
         Ah=Ah,      #roothairs only, surface area of the root hairs in each compartment
         unit=unit
      )
   return (results)
}

#########roothair functions################
calcCh<-function(C,Km,Cmin,S4){
    return ( ( (   (C-Km+Cmin-S4)+sqrt(  (-C+Km-Cmin+S4)^2-4*C*Cmin+4*S4*Cmin+4*C*Km) )/2) - Cmin)
}
calcdCh<-function(C,Km,Cmin,S4){
   return ( 0.5 + 0.5*(C + Km - Cmin - S4)/sqrt( (-C+Km-Cmin+S4)^2-4*C*Cmin+4*S4*Cmin+4*C*Km )   )
}

calcIh<-function(Ch,Imax,Km,Cmin,Ah){
    return  ((Ah*Imax*Ch)/(Km+Ch))
}
calcdIh<-function(Ch,dCh,Imax,Km,Cmin,Ah){
   return ( Ah*Imax*( (dCh/(Km+Ch)) - (Ch*dCh/((Km+Ch)^2)) ) )
}

MMuptake<-function(const,cr){ 
   return ((const$r0*pi*2*((const$Imax*(cr-const$Cmin)/(const$Km+cr-const$Cmin))-const$E))/const$dt)
}

##########main model function#######
barber<-function(param,msg){
   #print message
   print(msg)
   #calculate constants
   const<-setConstants(param)
   k=const$k

   #build tridiagonal jacobian matrix
   Jad=array(const$D1,k)
   Jad[k]=1
   if(const$case==2) Jad[k]=const$D1[k]+const$S2[k]*const$A1 # case II only
   Jb=array(const$S1,k) 
   Jb[k]=0
   if(const$case==2) Jb[k]=-2 # case II only
   Jc=array(-const$S2,k)
   Jc[1]=-2

   #result matrix
   m=matrix(0,k,length(const$time)+1)
   #mIh=matrix(0,k,length(const$time)+1)

   #initial concentration
   Cn=array(const$Cli,k)
   Cn1=Cn
   Cnp=Cn
   Cnp[1]=Cn[1]*2
   #try a better initial estimate for Cn1[1]
   #Cn1[1]=0.7*Cli

   #right hand side: f
   f=array(0,k)

   #warning messages on
   warn=TRUE
   warn1=TRUE
   warn2=TRUE
   warn3=TRUE

   #time loop
   for (t in const$time){
      m[,t]=Cn # store data

      #concentration at the root surface of the root hairs (assuming steady state diffusion)
      Ch=calcCh(Cn,const$Km,const$Cmin,const$S4)
      dCh=calcdCh(Cn,const$Km,const$Cmin,const$S4)
      
      #uptake by root hairs
      Ih=calcIh(Ch,const$Imaxh,const$Km,Cmin,const$Ah)
      dIh=calcdIh(Ch,dCh,const$Imaxh,const$Km,const$Cmin,const$Ah)

      #Hypothetical Concentrations for inner and outer boundary conditions
      Co=0
      if(const$case==2) Co=Cn[k-1]-const$A1*Cn[k]
      Cea=max(0,Cn[1]-const$Cmin)
      Ci=Cn[2]-const$S3*(((const$Imax*Cea)/(const$Km+Cea))-(const$v0*Cn[1]))

      #mIh[,t]=Ih # store data

      #Solve the tridiagonal matrix using an iterative method
      err=1e8   #error
      count=0
      maxcount=20
      perr=array(-1,maxcount)
      countreset=TRUE

      #initial estimate for itterative method
      Cn1=Cn  #-0.9*(Cnp-Cn)
      pCd=0
      
      #newton-ralphson itterative loop
      while(err>const$small){
         #loop counter
         count = count+1




         #concentration at the root surface of the root hairs (assuming steady state diffusion), and uptake flux by root hairs
         Ch1=calcCh(Cn1,const$Km,const$Cmin,const$S4)
         dCh1=calcdCh(Cn1,const$Km,const$Cmin,const$S4)
         Ih1=calcIh(Ch1,const$Imaxh,const$Km,const$Cmin,const$Ah)
         dIh1=calcdIh(Ch1,dCh1,const$Imaxh,const$Km,const$Cmin,const$Ah)
         
      #inner root concentrations
         Cea=max(0,Cn1[1]-const$Cmin)
         Ci1=Cn1[2]-const$S3*(((const$Imax*Cea)/(const$Km+Cea))-(const$v0*Cn1[1])) 
   
         #f based on estimates (f should become 0)
         f[1]=const$S1[1]*Ci1 + const$D1[1]*Cn1[1] - const$S2[1]*Cn1[2] +  const$S1[1]*Ci - const$D2[1]*Cn[1] - const$S2[1]*Cn[2]
         f[2:(k-1)]=   const$S1[2:(k-1)]*Cn1[1:(k-2)] + const$D1[2:(k-1)]*Cn1[2:(k-1)] - const$S2[2:(k-1)]*Cn1[3:k] +  const$S1[2:(k-1)]*Cn[1:(k-2)] - const$D2[2:(k-1)]*Cn[2:(k-1)] - const$S2[2:(k-1)]*Cn[3:k]
         f[k]= Cn1[k]-const$Cli
         if(const$case==2) f[k]=-2*Cn1[k-1] + (const$D1[k]+const$S2[k]*const$A1)*Cn1[k] +  const$S1[k]*Cn[k-1] - const$D2[k]*Cn[k] - const$S2[k]*Co
         if(const$roothairs){
            ft=const$D3*(Ih1+Ih)
            #ft[1]=2*ft[1] #todo now the two is in D3. Check what is right
            #print(paste(count,err,max(Ih),max(abs(MMuptake(const,Cn1[1])))))
            f=f+ft
         }
         
         #time -1 for -1*f(x1)=(x2-x1)/f'(x1). 
         f=f*-1

         #Jacobian matrix based on estimate for inner boundary condition
         Ja=Jad #reset
         Ja[1]= -1 * const$S1[1]*const$S3*(((-const$Imax*Cea)/((const$Km+Cea)^2))+(const$Imax/(const$Km+Cea))-const$v0)+const$D1[1] 
         #add root hairs
         if(const$roothairs){
            Ja=Ja+const$D3*dIh1   
         }

         #solve  J*Cd=-f
         Cd=tridiagonalMatrixSolver(Ja,Jb,Jc,f)

         #error estimate
         perr[count]=err
         err=max(abs(Cd))
         if(is.nan(err)) {
            if(warn1) {
              print(paste("NaN values",t,"further warnings suppressed. Extras info: ",msg))
              warn1=FALSE
            }
            err=0 #terminate the loop
         }else if(err-perr[count]>err*1e-3){
            if(warn2) {
               print(paste("Not converging at time:",t,"further warnings suppressed. Extras info: ",msg))
               warn2=FALSE
            }
            err=0 #terminate the loop
         }else if(count>maxcount && err>100*const$small) {
            if(warn3) {
               print(paste("Can't resolve at time:",t,"Error",err,"Extras info: ",msg))
               warn3=FALSE
            }
            err=0  #terminate the loop
         }

         #new estimates for the next itteration
         Cn1=Cn1+Cd 
         Cn1=Cn1+((Cn1-const$Cmin)<0)*const$Cmin
         pCd=Cd
         
      }#end of itteration
      
      #use new concentrations for next timestep
      if(warn){
         if(min(Cn1)<0) {
            print(paste("Negative concentrations. ",msg));
            warn=FALSE
         }
         if(max(Cn1)>1.1*const$Cli){
            print(paste("Concentrations higher than initial concentrations. ",msg))
            warn=FALSE
         }
      }
      #print(paste("Solved in ",count," itterations"))
      Cnp=Cn
      Cn=Cn1
   }
   #store the last result
   m[,(length(const$time)+1)]=Cn1
   #return (list(m=m,Ih=mIh))
   return (m)
}

###### roothair function plot to proof that the derivative is right
#for low C, the dIh plot is a little off. I think it may just be a numerical error as Ih is very small?
plotCn<-function(param){
   #calculate constants
   const<-setConstants(param)
   k=const$k

   #build array with concentrations
   concentration=1.5*((1:k)/k)*const$Cli

   #calculate concentrations at root surface
   Ch=calcCh(concentration,const$Km,const$Cmin,const$S4)
   #analytical differentiation
   dCh=calcdCh(concentration,const$Km,const$Cmin,const$S4)
   #numerical differentiation
   dChc=dCh
   dChc[1]=NA
   dChc[k]=NA
   for(i in 2:(k-1)){
         dChc[i]=(Ch[(i+1)]-Ch[(i-1)])/(concentration[(i+1)]-concentration[(i-1)])
   }
   
   #calculate uptake by root hairs
   Ih=calcIh(Ch,const$Imaxh,const$Km,const$Cmin,const$Ah[1])
   #analytical differentiation
   dIh=calcdIh(Ch,dCh,const$Imaxh,const$Km,const$Cmin,const$Ah[1])
   #numerical differentiation
   dIhc=dIh
   dIhc[1]=NA
   dIhc[k]=NA
   for(i in 2:(k-1)){
         dIhc[i]=(Ih[(i+1)]-Ih[(i-1)])/(concentration[(i+1)]-concentration[(i-1)])
   }
   
   #plot 4 plots
   par(mfcol=c(2,2) )   
   plot(concentration,Ch,type="l")
   plot(concentration,dCh,type="l")
   lines(concentration,dChc,col=2,lty=2)
   legend("bottomright",legend=c("analytical solution (used in model)","numerical solution"),col=c(1,2),lty=c(1,2),bty="n")
   plot(Ch,Ih,type="l")
   plot(Ch,dIhc,type="l")
   lines(Ch,dIh,col=2,lty=2)
   legend("topright",legend=c("analytical solution (used in model)","numerical solution"),col=c(1,2),lty=c(1,2),bty="n")
}

##################nutrient uptake##########################
barberTotN<-function(param,m){
   #calculate constants
   const<-setConstants(param)
   k=const$k
   dt=const$dt
   l=length(const$time)
   jm=l+1
   time=c(0,const$time*const$dt)
   
   #calculate nutrient uptake rate per cm root
   Nupt=MMuptake(const,m[1,])

   #calculate total (cumulative) nutrient uptake 
   TNupt=0
   TotNupt<-array(0,jm)
   for (i in 2:jm){
      #TNupt=TNupt+(dt*(Nupt[i-1]+Nupt[i])/2)
      TNupt=TNupt+(dt*(Nupt[i-1]))
      TotNupt[i]=TNupt
   }

   #uptake by roothairs
   TNupth=0
   TotNupth<-array(0,jm)
   Nupth=array(0,jm) #total value of Ih for each time, Ih in uMol/timestep/cylinder
   if(const$roothairs){
      #Ah=matrix(const$Ah*const$vol,k,jm) #convert Ah from cm2/cm3 to cm2/cm by multiplying but the volume/cm root
     Ah=matrix(const$Ah*const$vol,k,jm)
      Ch=calcCh(m,const$Km,const$Cmin,const$S4)
      Ih=calcIh(Ch,const$Imaxh,const$Km,const$Cmin,Ah)
      #print(paste("Ah=",max(Ah)))
      #print(paste("Ch=",max(Ch)))
      #print(paste("C=",max(m)))
      #print(paste("Ih=",max(Ih)))
      Nupth=colSums(Ih) #total value of Ih for each time step, Ih in μmol/timestep/cylinder
      Nupth=Nupth/dt #per day
      for (i in 2:jm){
         #TNupth=TNupth+(dt*(Nupth[i-1]+Nupth[i])/2)
         TNupth=TNupth+(dt*(Nupth[i-1]))
         TotNupth[i]=TNupth
      }
   };
   
   res<-list(
      Nupt=Nupt, #array of uptake rates over time
      Nupth=Nupth, #array of roothair uptake rates over time
      TNupt=TNupt, #cumulative uptake at end of simulation 
      TNupth=TNupth, #cumulative uptake at end of simulation for root hairs
      TotNupt= TotNupt, #array of cummulative uptake over time
      TotNupth= TotNupth
   )
   return(res)
}
##################plotting functions##########################
plotbarber<-function(param,m,filename="none"){
   #if(filename!="none") png(filename=paste(filename,".png"),width=1200,height=1200)
   if(filename!="none") svg(filename=paste(filename,".svg"),width=10,height=7)
   #calculate constants
   const<-setConstants(param)
   k=const$k
   dt=const$dt
   l=length(const$time)
   jm=l+1
   time=c(0,const$time*const$dt)
   fr=8 #number of lines per plot
   
   #unit
   unit=const$unit

   #calculate nutrient uptake
   data<-barberTotN(param,m)   

   #4 graphs
   oldpar=par(no.readonly=TRUE)
   par(mfcol=c(2,3),ps=11,lwd=2 ) #one graph with large text

   #plot the m matrix
   ylim=range(pretty(c(0,1.1*max(m))))
   plot(const$r,m[,l+1], main="Concentration/Space",xlab="Distance from the root",ylab="Concentration", type='l',ylim=ylim)
   a=as.integer(l/fr)
   for(i in c(1,(1:(l/a))*a)){
      lines(const$r,m[,i])
   }

   #plot the depletion 
   plot (time,m[1,], main="Concentration/Time",xlab="Time",ylab="Concentration", type='l', ylim=ylim) #plot
   a=as.integer(k/fr)
   for(i in c(1,(1:(k/a))*a)){
      lines(time,m[i,]) #plot
   }

   #plot the uptake rate by time
   sumNupt=100*(data$Nupt+data$Nupth)
   plot (time,sumNupt, main="Uptake rate",xlab="Time (day)",ylab=paste("Nutrient uptake rate (",unit,"/m/day)",sep=""), type='l', ylim=range(pretty(c(0,max(sumNupt))))) #plot
   lines(time,100*data$Nupt, lty=2)
   lines(time,100*data$Nupth, lty=3)
   legend("topright", legend = c("total","root","roothairs"),lty=c(1,2,3),bty='n')

   #plot the michaelis mente uptake curve
   crm=max(10*const$Km, 1.1*max(m))
   cr=seq(const$Cmin,crm,(crm-const$Cmin)/100) #concentration range
   vm=100*const$r0*pi*2*const$Imax/const$dt
   
   ylab=paste("Nutrient uptake (",unit,"/m/day)",sep="")
   plot(cr,100*MMuptake(const,cr), main="Michaelis Menten",xlab="Concentration uM",ylab=ylab, type='l',ylim=range(pretty(c(0,vm))), xlim=range(pretty(c(0,crm))) ) #plot
   lines(m[1,],100*data$Nupt,col=2,lty=3,lwd=4)
   lines(range(cr),c(vm,vm),col=2,lty=1,lwd=2)
   
   #plot uptake curves
   sumNupt=100*(data$TotNupt+data$TotNupth)
   plot (time,sumNupt, main="Cum. uptake",xlab="Time (day)",ylab=ylab, type='l', ylim=range(pretty(c(0,max(sumNupt))))) #plot
   lines(time,100*data$TotNupt, lty=2)
   lines(time,100*data$TotNupth, lty=3)
   legend("topleft", legend = c("total","root","roothairs"),lty=c(1,2,3),bty='n')

   #total change in mass ()
   #sumNupt=100*(data$TotNupt+data$TotNupth)
   dc=const$Cli - m #concentration change
   mv=sweep(dc,MARGIN=1,(const$b)*const$vol,`*`) #mulitply concentration with volume and buffer to get amounts
   dcs=100*colSums(mv)
   #dcs=100*colSums((mv[1:(k-1),]+mv[2:k,])/2) #spatial integration
   plot (time,sumNupt, main="Cum. uptake",xlab="Time (day)",ylab=ylab, type='l', ylim=range(pretty(c(0,max(sumNupt,dcs))))) #plot
   lines(time,dcs, lty=2)
   legend("topleft", legend = c("total uptake","total solute change"),lty=c(1,2),bty='n')
   
   
   par(oldpar)
   if(filename!="none") dev.off()

}

######sensitivity analysis########
runRange<-function(param,range,pname){
   spl=unlist(strsplit(pname,'~'))
   p<-param
   unit=setConstants(param)$unit
   ov<-as.numeric(param[spl])
   Ntot<-array(0,length(range))
   Nroot<-array(0,length(range))
   Nhairs<-array(0,length(range))
   #todo this could run in parallel with the for each package
   count=0
   for(i in range){
     count=count+1
     p[spl]<-i*ov
     m<-barber(p,msg=paste("Sensitivity analysis for",pname,"=",p[pname]))
     res<-barberTotN(p,m)
     Ntot[count]<-res$TNupt+res$TNupth
     Nroot[count]<-res$TNupt
     Nhairs[count]<-res$TNupth
     s=''
     for(trp in i*ov) s=paste(s,trp)
     print(paste("Result for",pname,"=",s,"is",Ntot[count],paste(" (",unit,"/m)",sep="")))#,"root only =",res$TNupt,"root hairs =",res$TNupth))
   }
   return (list(Ntot=Ntot,Nroot=Nroot,Nhairs=Nhairs))
}


sensitivityBarber<-function(param,sensitivity_names=c("Cli","r0","r1","Imax","Km","De","b","v0","lh","rh","Nh","Nh~lh")){
   #range=(1:10)/5
   #range=(10*((1:15)^3)+500)/5620 #0.09-6 in 15 exponentiele stappen
   #range=(10*((1:10)^2)+200)/450  #.07-2
   #range=c(0.5,1.5)
   range=c(0.5,0.6,0.7,0.8,0.9,1.0,1.2,1.4,1.6,1.8,2.0,2.2)
   
   #names=c("r0","Imax")
   data<-list(range=range)
   for (pn in sensitivity_names){
      new=runRange(param,range,pn)
      data=c(data,list(new=new$Ntot))
      names(data)[names(data)=="new"]=pn
   }
   return(data)
} 
  
plotSensitivityBarber<-function(data,filename="none",exclude=NULL,include=NULL,maintitle="Sensitivity Analysis Barber-Cushman Model"){
   range=data$range
   
   #exclude data (todo: don't know how to do this in an R ellegant way without the loop)
   drm=NULL
   if(length(exclude)>0){
      for(i in 2:length(names(data))){
         if(sum(exclude==names(data)[i])!=0) drm=c(drm,-i)
      }
   }
   if(length(include)>0){
      for(i in 2:length(names(data))){
         if(sum(include==names(data)[i])==0) drm=c(drm,-i)
      }
   }
   if(length(drm)) data=data[drm]

   #4 graphs combined
   kn=length(names(data))
   if(filename!="none")  svg(filename=paste(filename,".svg"))#,width=800,height=800)
   oldpar=par(no.readonly=TRUE)
   par(mfcol=c(1,1),ps=11,lwd=4 ) #one graph with large text
   ymax=0  
   for(p in 2:kn){
      ymax=max(ymax,unlist(data[p]))   
   }
   ymax=100*ymax
   yrange=c(0,ymax)#/max(range) 
   ylab=paste("Cumulative nutrient uptake (amount/m)",sep="")
   plot(c(min(range),max(range)),yrange, main=maintitle,xlab="Relative change parameter",ylab=ylab, type='l',lty=2,col=0)
   legend(mean(range),min(yrange),names(data[2:kn]), col = 2:kn,lty=1, bty='n',xjust=1,yjust=0,y.intersp=1.2)
   #unity=round(Nv0[(range-1)**2==min((range-1)**2)],3)
   #text(1.1, unity,unity,pos=1)
   for (p in 2:kn){
      lines(range,100*unlist(data[p]),col=p)
   }
   par(oldpar)
   if(filename!="none") dev.off()

}



