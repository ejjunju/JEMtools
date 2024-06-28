#run devtools::document() to create documentation

#HYDROLOGICAL MODELLING (with HBV)###############################################################################

#' This functions computes the day number since january first of the repsective year
#' @param date
#' @return juliaday
#' @examples
fx.julian.day<-function(x=as.Date("2005-01-31")){
  # Julian days from date **********************************************************************************************
  yr<-as.numeric(format(x,'%Y'))
  origin<-as.Date(paste0(yr,"-01-01"))
  jdate<-as.numeric(x-origin+1) #julian(x,origin = origin)+1 #days since 1970-01-01
  return(jdate)
}

#' This functions computes the hydrological year for a date given the month number in which the hydrological year starts
#' The water year is designated by the calendar year in which it ends,
#' so ifthe 2024 water year started on september 1, 2023, it will end on August 30, 2024
#' @param date
#' @return juliaday
#' @examples
fx.hydrological.year<-function(x=as.Date("2005-01-31"),month=9){
  # Hydrological years from date given starting month*******************************************************************
  yr<-as.numeric(format(x,'%Y'))
  mo<-as.numeric(format(x,'%m'))
  hyear<-apply(data.frame(mo,yr),1,function(x) {if(x[1]>=month){return(x[2]+1)} else {return(x[2])}})
  return(hyear)
}

#' volumetric efficiency
#' The VE represents the fractional volumetric difference between the simulated
#' and observed streamflows. It ranges from 0 to 1.
#' A perfect value of 1 indicates that the volume of water
#' predicted by the model matches the observed volume.
#' This criterion gives the same weight to any flow range
#' (e.g., slow recession and rapid rising flow),
#' relaxing the constraint of model residuals heteroscedasticity
#' references
#' https://hess.copernicus.org/articles/26/197/2022/hess-26-197-2022-supplement.pdf
#' https://www.scielo.br/j/rbrh/a/RhjqdXxzfkPH79GpQS9VjWd/?lang=en&format=pdf
#' @import zoo
#' @param Qsim
#' @param Qobs
#' @return KGE
#' @examples
fx.ve<-function(Qsim,Qobs){  #calc
  # Volumetric Efficiency (VE)  criteria********************************************************************************
  qna<-data.frame(Qsim,Qobs)
  q<-qna[!apply(qna,1,function(x) any(is.na(x))),]#remone nas
  VE<-1 - (sum(abs(q$Qsim-q$Qobs))/sum(q$Qobs))
  return(VE)}


#' Nash-Sutcliffe Efficiency (NSE).
#' This functions calculates the Nash-Sutcliffe Efficiency (NSE).
#' NSE = 1 indicates perfect correspondence between simulations and observations
#' NSE = 0 indicates that the model simulations have the same explanatory power as the mean of the observations
#' NSE < 0 indicates that the model is a worse predictor than the mean of the observations
#' #' references
#' https://hess.copernicus.org/articles/26/197/2022/hess-26-197-2022-supplement.pdf
#' https://www.scielo.br/j/rbrh/a/RhjqdXxzfkPH79GpQS9VjWd/?lang=en&format=pdf
#' @import zoo
#' @param Qsim
#' @param Qobs
#' @param with.VE Relative Error  logical. if TRUE calculates NSE*VE
#' @return NSE
#' @examples
fx.nse<-function(Qsim,Qobs,with.VE=FALSE){
  # Nash Sutcliffe Efficiency (NSE criteria ****************************************************************************
  #calc
  qna<-data.frame(Qsim,Qobs)
  q<-qna[!apply(qna,1,function(x) any(is.na(x))),]#remone nas
  nse<-1-sum((q$Qsim-q$Qobs)^2)/ sum((q$Qobs-mean(q$Qobs))^2)

  if(with.VE ==TRUE){
    VE<-fx.ve(Qsim,Qobs)
    nse<-nse*VE
  }
  return(nse)
}

#' Kling-Gupta efficiency scores
#' This functions calculates the Kling-Gupta efficiency scores (KGE).
#' KGE = 1 indicates perfect agreement between simulations and observations. Analogous to NSE = 0,
#' certain authors state that KGE < 0 indicates that the mean of observations provides better estimates than simulations
#' #' references
#' https:#hess.copernicus.org/preprints/hess-2019-327/hess-2019-327.pdf
#' https://hess.copernicus.org/articles/26/197/2022/hess-26-197-2022-supplement.pdf
#' https://www.scielo.br/j/rbrh/a/RhjqdXxzfkPH79GpQS9VjWd/?lang=en&format=pdf
#'
#' @import zoo
#' @param Qsim
#' @param Qobs
#' @return KGE
#' @examples
fx.kge<-function(Qsim,Qobs){
  # Kling-Gupta Efficiency (KGE) scores ********************************************************************************
  #calc
  qna<-data.frame(Qsim,Qobs)
  q<-qna[!apply(qna,1,function(x) any(is.na(x))),]#remone nas
  r<-cor(q$Qsim,q$Qobs)#linear correlation between observations and simulation
  AlphA<-sd(q$Qsim)/sd(q$Qobs) #𝛼 a measure of the flow variability error,
  BetA<- mean(q$Qsim)/sd(q$Qobs)  #𝛽 a bias term
  kge <- 1- ((r-1)^2 + (AlphA-1)^2 + (BetA-1)^2)^0.5
  return(kge)
}

#' Plot Timesereies Qobs and Qsim with objective criteria
#' @import zoo
#' @param Qsim
#' @param Qobs
#' @examples
fx.plot.hbv.calib<-function(Qsim=NULL,Qobs=NULL,Time=NULL,las=2){
  # Simple plots of HBV calibration results Qobs,Qsim and Time *********************************************************
  if(is.null(Time)){
    Time<-1:length(Qobs)
    xp<-1:length(Time)
    xl<-xp
  } else{
    xp<-seq(Time[1],Time[length(Time)]+31,by="month")
    xl<-paste0(as.numeric(format(xp,'%Y')),
               "-",
               sprintf("%02d",as.numeric(format(xp,'%m'))),
               "-",
               sprintf("%02d",as.numeric(format(xp,'%d')))

               )
  }
  qna<-data.frame(Qsim,Qobs)
  yp<-pretty(range(qna,na.rm=TRUE))
  #layout(matrix(1:2,nrow=1,ncol=2),widths=c(5,2))
  #par(mar=c(4,4,1,1),oma=c(1,1,1,1))
  plot(Time,Qobs,type="l",col="blue",
       ylim=c(min(qna,na.rm=TRUE),1.2*max(qna,na.rm=TRUE)),
       ylab="Flow (m³/s)",xlab="",las=2,yaxs="i",xaxt="n",yaxt="n")
  VE<-fx.ve(Qsim,Qobs)
  NSE<-fx.nse(Qsim,Qobs,with.VE = F)
  NSE.VE<-fx.nse(Qsim,Qobs,with.VE = T)
  KGE<-fx.kge(Qsim,Qobs)
  obj<-data.frame(NSE=NSE,VE=VE,NSE.VE=NSE.VE,KGE=KGE)
  lines(Time,Qsim,col="red")
  #plot(0, type = 'n', axes = FALSE, ann = FALSE)
  legend("top",paste(names(obj),"=",round(obj,2)),cex=0.6,horiz=TRUE)
  #par(mar=c(5,4,1,1),oma=c(1,1,1,1))
  #box(which = "outer",col="blue")
  axis(1,xp,xl,las=2,cex.axis=0.8)
  axis(2,yp,yp,las=2,cex.axis=0.8)
  abline(v=xp,lw=0.1,lty=2,col="grey80")
  abline(h=yp,lw=0.1,lty=2,col="grey80")
}


#' This functions reads data from the input files for the HBV model
#' sample files are provided
#' @import zoo
#' @import readxl
#' @import tidyverse
#' @import gdata
#' @param dir
#' @param Timeseries Timeseries csv with Time,Prec,PET,Qobs. all must be there. Qobs can be 0 if not run
#' @param Free.txt Free parameters csv Parname, Value,low,high,
#' @param Fixed.txt Fixed param,eteres csv csv file; tempalte provided
#' @param serie.start.date
#' @return hbv data object
#' @examples
fx.hbv.daily.getdata<-function(dir,
                               Timeseries.txt,
                               Free_Par.txt,
                               Fixed_Par.txt,
                               serie.start.date=NULL){
  ## Data input for HBV Model ******************************************************************************************

  #dir<-"D:/WORK/Tools/R/Packages and Apps/ejjhydrotools/data/hbv/"
  Timeseries.txt=paste0(dir,"/",Timeseries.txt)
  Free.txt=paste0(dir,"/",Free_Par.txt)
  Fixed.txt=paste0(dir,"/",Fixed_Par.txt)


  fri<-read.table(file = Free.txt)
  fxd<-read.table(file = Fixed.txt,nrow=4)
  zone<-read.table(file = Fixed.txt,skip=4,nrow=4)
  timeseries<-read.table(file=Timeseries.txt,header = TRUE)
  if(is.null(serie.start.date)) {serie.start.date=timeseries[1,1]}
  timeseries$Date<-seq(as.Date(serie.start.date),as.Date(serie.start.date)+ nrow(timeseries)-1,1)
  timeseries$jdate<-fx.julian.day(timeseries$Date)
  timeseries$hyear<-fx.hydrological.year(timeseries$Date,month=9)

  input<-list(files=list(dir = dir,
                        Fixed.txt = Fixed.txt,
                        Free.txt = Free.txt,
                        Timeseries.txt = Timeseries.txt),
             data=list(fri = fri,
                       fxd = fxd,
                       zones = zone,
                       serie.start.date = serie.start.date,
                       timeseries = timeseries,
                       hyear=unique(timeseries$hyear))
  )
  return(input)

}

#' This is the hbv model for the given parameters and entire dataset
#' @import zoo
#' @import dplyr
#' @import readxl
#' @import tidyverse
#' @import gdata
#' @param prm free-parameters
#' @param input inputdata from fx.hbv.daily.getdata()
#' @param Snow2GlacierOption TRUE, convert snow to glacier at end of hydrological year
#' @param Snow2GlacierJulianDate  as.Date("2000-09-01"), Date for calc of Snow2GlacierJulianday.
#' @param calibrate Makes the function return obj.crit if TRUE to enable optimisation
#' @param obj.criteria Choose an optmising objective criteria from ("NSE","VE","NSE.VE","KGE")
#' @param icall controls printing of visual obj.criteria during calibration if provided by external calib function
#' @return either objective critria or results object
#' @examples
fx.hbv.model.serial<-function(prm = NULL,
                              input, #data
                              Snow2GlacierOption = TRUE,#convert snow tro glacier at end of hydrological year
                              Snow2GlacierJulianDate= as.Date("2000-09-01"), #For calc of Snow2GlacierJulianday
                              calibrate=TRUE, # makes the function return obj.crit if TRUE to enable optimisation
                              obj.criteria="NSE", # chhose an objective criteria
                              icall=1,
                              calib.plots=FALSE){
  # HBV Model (Serial simulation for entire series) *******************************************************************
  #comment after testing
  # prm<-NULL
  # input<-fx.hbv.daily.getdata(dir="D:/WORK/Tools/R/Packages and Apps/ejjhydrotools/data/hbv",
  #                             Timeseries ="Timeseries.txt",
  #                             Free_Par.txt="Free_Par.txt",
  #                             Fixed_Par.txt="Fixed_Par.txt",
  #                             serie.start.date="2020-09-01")
  # Snow2GlacierOption = TRUE
  # Snow2GlacierJulianDate= as.Date("2000-09-01")
  # calibrate=TRUE
  # obj.criteria="NSE"
  # icall=1

  #inputs
  fri<-input$data$fri
  fxd<-input$data$fxd
  zones<-input$data$zones
  timeseries<-input$data$timeseries

  #Free parameters
  if(!is.null(prm)){fri[,2] =prm}#read from prm and update fri
  prm<-fri[,2];names(prm)<-fri[,1]
  PKORR= prm[1];SKORR= prm[2];HPKORR= prm[3];TX= prm[4];#snow/rain
  TS= prm[5];CX= prm[6];CPRO= prm[7];PRO = prm[8]; #snowmelt
  FC= prm[9];BETA= prm[10];LP= prm[11];#Soil mositure
  KUZ2= prm[12];KUZ1= prm[13];UZ1= prm[14];PERC= prm[15];KLZ= prm[16];
  GCX= prm[17];#galcier melt
  Tlp= prm[18];Tlo= prm[19];
  SM0= prm[20]*FC; ##initial SM SMO is SMO/FC so we have converted it to mm here
  UZ0= prm[21];LZ0= prm[22]#initial upper and lower zone


  #zonal data
  zone.names<-zones[,1]
  zone<-as.data.frame(t(zones[,-1]))
  names(zone)<-zone.names
  nZ<-nrow(zone)
  z<-zone$masl
  SN0<-zone$SNO #initial snow by zone
  SW0<-zone$SWO #initial snow water by zone
  gIce<-zone$gIce #initial glacier proportion per zone in %


  #Precipitation lapse rates
  hpkorr<-rep( HPKORR,nZ)

  #Mountain Melt Increase factor MiMt
  #not implemented so set to one (forests)
  #can be larger for bare mountains
  MiMt = 1;


  #Fixed parameters
  Zp=fxd$V2[1];#elevation of precipitation station
  Zt=fxd$V2[2];#elevation of temperature station
  pLKs=fxd$V2[3]; #"lake percentage"
  Akm2 = fxd$V2[4]; #VCatchment area

  #data
  DY = nrow(timeseries)#length of dataset
  tp<-timeseries$Temp
  pr<-timeseries$Prec
  PE<-timeseries$PET
  Qobs<-timeseries$Qobs
  Time<-timeseries$Date
  JULIANS<-timeseries$jdate
  Snow2GlacierJulianDay<-fx.julian.day(Snow2GlacierJulianDate)


  #model result containers
  ##zonal
  P<- array(0,dim=c(DY,nZ))#Precipitation for each zone
  Ta<- array(0,dim=c(DY,nZ)) #Temperature for each zone
  R<- array(0,dim=c(DY,nZ)) #Rainfsll for each zone
  S<- array(0,dim=c(DY,nZ)) #Snowfall for each zone
  SN<- array(0,dim=c(DY,nZ)) # Dry Snow for each zone
  ATM<-array(0,dim=c(DY,nZ)) # Snow that is available to melt
  M<-array(0,dim=c(DY,nZ)) # Snow melt
  SW<-array(0,dim=c(DY,nZ)) # Snow water
  FR<-array(0,dim=c(DY,nZ)) # Snow freeze
  ST<-array(0,dim=c(DY,nZ)) # Snow water threshold
  #Portion of rainwater that is converted to snow water when there is a snow water defict
  RainToSW<-array(0,dim=c(DY,nZ))
  SWToSoil<-array(0,dim=c(DY,nZ)) #Excess sow water that goes to soil
  InS<-array(0,dim=c(DY,nZ)) #remaining rainwater that goes to INSOIL
  SNOW = array(0,dim=c(DY,nZ)) #Total Snow
  GM = array(0,dim=c(DY,nZ)) #Glacier melt

  #model result containers
  #Not zonal
  Pav = array(0,dim=c(DY,1)) #Average Areal Prec
  Tav =  array(0,dim=c(DY,1)) #Average Temperature
  SNOWav = array(0,dim=c(DY,1)) #Average Snow
  SnowCover = array(0,dim=c(DY,1)) #Snow cover
  Meltav =  array(0,dim=c(DY,1)) #Average Melt
  INSOIL =  array(0,dim=c(DY,1)) #Average INSOIL
  GMLT =  array(0,dim=c(DY,1)) #Average Glacier Melt
  SM =  array(0,dim=c(DY,1)) #Soil mositure
  EA =  array(0,dim=c(DY,1)) #Actual Evapotranspiration
  dUZ =  array(0,dim=c(DY,1)) # direct runff
  UZ =  array(0,dim=c(DY,1)) # Upper Zone
  LZ =  array(0,dim=c(DY,1)) # Lower Zone
  perc =  array(0,dim=c(DY,1)) # Percolation to Lower Zone
  QUZ1 =  array(0,dim=c(DY,1)) # Flow upper zone 1;mm
  QUZ2 =  array(0,dim=c(DY,1)) # Flow upper zone 2;mm
  QLZ =  array(0,dim=c(DY,1)) # Flow lower zone;mm
  Qsim  =  array(0,dim=c(DY,1)) # Flow m3/s


  for (tstep in 1:DY){
    for (zn in 1:nZ){ #for each elevation zone
      ## Snow Routine**************************************************************************************************
      #calculate rainfall and snowfall
      #temperature at precipitation gauge calculated from lapse rate and temp gauge elevation
      TZp =tp[tstep] + (Zt -Zp) / 100 * Tlp;
      #corrected rainfall
      if (TZp > TX) {prc = PKORR *pr[tstep]}else { prc = SKORR *pr[tstep] }

      #Lapse precipitation
      P[tstep,zn] = max((prc * (1 + ((z[zn] -Zp) / 100) * hpkorr[zn] / 100)), 0);

      #lapse Temperature
      TLR = Tlo #temperature lapse rate is either wet (Tlp) or dry (Tlo)
      if (P[tstep,zn] > 0) { TLR = Tlp; } else { TLR = Tlo; }
      Ta[tstep,zn] =tp[tstep] + TLR * (z[zn] -Zt) / 100;

      #Rainfall & snow fall for eath zone
      #R[tstep,zn] = (Ta[zn] > TX) ?P[tstep,zn] : 0;
      if(Ta[tstep,zn] > TX){R[tstep,zn]<-P[tstep,zn]}else{R[tstep,zn]<-0}
      S[tstep,zn] =P[tstep,zn] -R[tstep,zn];

      #Snow-routine
      #Begin snow storage % cumulating per tstep
      if (zn == 1) {SnowCover[tstep] = 0}

      #Melt M
      CX = MiMt*CX #in future cosider separating melt in mountains and forests. Mimt=1 forests and >1 mts
      Mmax<-0
      Mmax1 = CX * (Ta[tstep,zn] - TS); #maximum melt capacity
      if(0 > Mmax1){Mmax = 0 } else {Mmax = Mmax1}#real melt Capacity shouldnt be less than 0

      #Existing snow
      SNE = 0;#initialize #existing dry snow temporary variable
      if (tstep == 1){ SNE =SN0[zn]}else {SNE =SN[tstep - 1,zn]};#existing dry snow
      #Snow2Glacier
      if ((JULIANS[tstep] ==Snow2GlacierJulianDay) & (Snow2GlacierOption == TRUE)) { SNE = 0; }

      #Snow available to melt
      ATM[tstep,zn] = SNE +S[tstep,zn];
      #M[tstep,zn] = (ATMS[tstep,zn] > Mmax) ? Mmax : ATMS[tstep,zn];
      if(ATM[tstep,zn] > Mmax){M[tstep,zn] =Mmax}else{M[tstep,zn] = ATM[tstep,zn]}

      #dry snow state
      SN[tstep,zn] = max(SNE +S[tstep,zn] -M[tstep,zn], 0);
      SN[tstep,zn] = SNE +S[tstep,zn] -M[tstep,zn];

      #Refreeze snow water
      Fmax=0
      Fmax1 = CX * PRO* (TS -Ta[tstep,zn]);
      if(0 > Fmax1) {Fmax =0 } else {Fmax =Fmax1} #freeze capacity

      SWE=0;#initialize #existing snow water temporary variable
      if (tstep == 1) {SWE = SW0[zn]} else {SWE =SW[tstep - 1,zn]} #exisiting snow water
      #Snow2Glacier
      if ((JULIANS[tstep] ==Snow2GlacierJulianDay) & (Snow2GlacierOption == TRUE)) { SWE = 0}
      ATF = SWE;#available to fz
      if(ATF > Fmax) {FR[tstep,zn] = Fmax} else {FR[tstep,zn] = ATF}

      #Snow water State
      ST[tstep,zn] = CPRO * 0.01*SN[tstep,zn];#max free water
      SW1 = SWE +M[tstep,zn] -FR[tstep,zn];
      SN[tstep,zn] =SN[tstep,zn] +FR[tstep,zn];# Add re-freeze to dry-snow 22-12-2011
      swd =ST[tstep,zn] - SW1;# snow water deficit after melt and freeze
      #update SW1 if greater than ST
      ##check effect in result
      ##possibley ipdate cpp model to show eefetcive melt(SWToSoil)
      if(swd<0){
        SWToSoil[tstep,zn]=SW1-ST[tstep,zn] #effective snow melt
        SW1=ST[tstep,zn]
        swd=0} else {
          SWToSoil[tstep,zn]=0;
          SW1=SW1
          swd=swd
          }
      if(swd >R[tstep,zn]) {RainToSW[tstep,zn] =R[tstep,zn]} else { RainToSW[tstep,zn] =swd};
      SW[tstep,zn] = SW1 +RainToSW[tstep,zn];

      #Total SNOW
      SNOW[tstep,zn] =SN[tstep,zn] +SW[tstep,zn];
      #track snow storage cover in each zone 22-12-201
      if (SNOW[tstep,zn] > 0) {SnowCover[tstep] = SnowCover[tstep] + 10; }

      #To soil moisture
      InS[tstep,zn] =R[tstep,zn] - RainToSW[tstep,zn] + SWToSoil[tstep,zn];

      #Glacier Melt
      #Goes directlt to upper zone
      if (SN[tstep,zn] < 0.01){
        GM[tstep,zn] = GCX * CX * max(Ta[zn] - TS, 0) * (gIce[zn] / 100);
      }   else {GM[tstep,zn] = 0}

    } #End elevation zones

      ## Snow Routine #####calculate Averages of zonal parameters for timestep
    Pav[tstep] = mean(P[tstep,],na.rm=TRUE);#Average Areal Prec
    SNOWav[tstep] = mean(SNOW[tstep,],na.rm=TRUE);#Average Snow
    Meltav[tstep] =  mean(SWToSoil[tstep,],na.rm=TRUE);#Average Melt
    INSOIL[tstep] =  mean(InS[tstep,],na.rm=TRUE);#Average INSOIL
    GMLT[tstep] =  mean(GM[tstep,],na.rm=TRUE);#Average Glacier Melt
    Tav[tstep] =  mean(Ta[tstep,],na.rm=TRUE);#Average Glacier Melt

    ## Soil Routine*****************************************************************************************************
    SMi=0;#starting Soil Mosture SMO (can be as a % of FC in par-file)
    if (tstep == 1) { SMi = SM0 } else { SMi = (SM[tstep - 1]); } #SM0

    #Actual Evaporation
    evap1=0;evap2=0; #22-12-2011
    evap1 = PE[tstep] * (1 -SnowCover[tstep] / 100); # SM>LP
    evap2 = PE[tstep] * SMi / (LP / 100 * FC) * (1 -SnowCover[tstep] / 100);# SM<=LP
    EA[tstep] = min(evap1, evap2);#22-12-2011

    #Soil Moisture Water Balance
    SMdef = FC - SMi;
    dUZ[tstep] = max((INSOIL[tstep] - SMdef -EA[tstep]),INSOIL[tstep] * ((SMi / FC)^BETA));
    dSM =INSOIL[tstep] - dUZ[tstep];
    SM[tstep] = max(0,SMi + dSM -EA[tstep]); #check this in c++

    ## Dynamic Routine**************************************************************************************************
    ###UPPER ZONE
    UZi=0 #temp UZ
    if (tstep == 1){ UZi = UZ0;} else {UZi = UZ[tstep - 1];} #initital conditions
    UZ[tstep] = UZi + dUZ[tstep] +GMLT[tstep];
    if(UZ[tstep] < PERC){ perc[tstep] = UZ[tstep]}else{ perc[tstep] = PERC}
    UZ[tstep] =UZ[tstep] - perc[tstep];
    QUZ1[tstep] = min(UZ[tstep], UZ1) * KUZ1;
    QUZ2[tstep] = max(UZ[tstep] - UZ1, 0) * KUZ2;
    UZ[tstep] = UZ[tstep] - (QUZ1[tstep] + QUZ2[tstep]);#Update State

    #LOWER ZONE
    if (tstep == 1){LZi = LZ0} else {LZi =LZ[tstep - 1]}#initital conditions
    LZ[tstep] = LZi + perc[tstep] + (Pav[tstep] -PE[tstep]) * (pLKs / 100);
    LZ[tstep] = max(LZ[tstep], 0);#To prevent negative LZ
    QLZ[tstep] =LZ[tstep] * KLZ;
    LZ[tstep] =LZ[tstep] - QLZ[tstep];#Update State

    #QSIM
    ##hard-coded daily flow m3/s with 24*3.6.
    ###modified added the multiplier *(1-pLKs/100) 22-03-2011
    Qsim[tstep] = ((QUZ1[tstep] + QUZ2[tstep]) * (1 -pLKs / 100) + QLZ[tstep]) * (Akm2 / (86.4));


  }#end  t-steps

  ## Objective Criteria ************************************************************************************************
  NSE<-fx.nse(Qsim,Qobs)
  VE<-fx.kge(Qsim,Qobs)
  NSE.VE<-fx.nse(Qsim,Qobs,with.VE=TRUE)
  KGE<-fx.kge(Qsim,Qobs)
  objective.criteria<-data.frame(NSE,VE,NSE.VE,KGE)
  i.return<-which(names(objective.criteria) ==obj.criteria)
  OBJ<-objective.criteria[,i.return]

  ## List of results ####
  Zonal<-list(Precipitation=P,Temperature=Ta,Rainfall=R,Snowfall=S,
              SnowDry=SN,AvailableToMelt=ATM,Melt=M,SnowWater=SW,Refreeze=FR,
              SnowWaterThreshold=ST,RainToSW=RainToSW,SWToSoil=SWToSoil,
              SnowTotal = SNOW,GlacierMelt = GM,
              InSoil=InS)

  Lumped<- data.frame(Time=Time,
                      Precipitation= Pav,Temperature= Tav,
                      SnowStorage = SNOWav,
                      SnowCoverPercentage =SnowCover,
                      SnowMeltWater = Meltav,
                      GlacierMelt= GMLT,
                      INSOIL = INSOIL,
                      SoilMositure= SM,ActualEvaporation= EA,DirectToUZ = dUZ,
                      UpperZone =UZ,LowerZone = LZ,
                      Percolation = perc,
                      QUZ1 = QUZ1, QUZ2 = QUZ2, QLZ = QLZ,
                      Qsim = Qsim,
                      Qobs = Qobs)

  AvQobs<-sprintf("%5.2f",mean(Qobs,na.rm=TRUE))
  AvQsim<-sprintf("%5.2f",mean(Qsim,na.rm=TRUE))
  #return
  if(calibrate==TRUE){
    objdf<-round(objective.criteria,2)
    for(k in 1:length(objdf)){objdf[k]<-sprintf("%5.2f",objdf[k]) }

    ishow<-cbind(icall,obj.criteria,OBJ=round(OBJ,2),objdf,AvQobs,AvQsim)
    names(ishow)[1]<-c("HBV_RUN:")
    if (icall == 1||icall %% 500 == 0) {
      if (icall > 0) {print(ishow)}
    }
    if(calib.plots==TRUE){
      fx.plot.hbv.calib(Qsim = Qsim,Qobs = Qobs,Time=Time)
      title(main=paste("\n Calibration [", obj.criteria," = ",sprintf("%5.2f",OBJ),"]"),outer=FALSE)
      }

    return(OBJ)
  } else {
    out<-list(Zonal=Zonal,Lumped=Lumped,objective.criteria=objective.criteria,AvQobs=AvQobs,AvQsim=AvQsim,criteria=obj.criteria)
    fx.plot.hbv.calib(Qsim = Qsim,Qobs = Qobs,Time=Time)
    title(main=paste("\n Simulation"),outer=FALSE)
    #title(main=paste("\n Simulation [", obj.criteria," = ",sprintf("%5.2f",OBJ),"]"),outer=FALSE)
    return(out)
  }


} #end function



#' This is the hbv model for the given parametes and entire dataset
#' But runs each hydrological year independently and
#' reuses the initial conditions for each hydrological year
#' useful for calibration for year by ear
#' @import zoo
#' @import dplyr
#' @import readxl
#' @import tidyverse
#' @import gdata
#' @param prm free-parameters
#' @param input inputdata from fx.hbv.daily.getdata()
#' @param Snow2GlacierOption TRUE, convert snow to glacier at end of hydrological year
#' @param Snow2GlacierJulianDate  as.Date("2000-09-01"), Date for calc of Snow2GlacierJulianday.
#' @param calibrate Makes the function return obj.crit if TRUE to enable optimisation
#' @param obj.criteria Choose an optmising objective criteria from ("NSE","VE","NSE.VE","KGE"),
#' @param hyears  Provide hydrological years to run. if NULL all years in dataset are used
#' @param icall controls printing of visual obj.criteria during calibration if provided by external calib function
#' @return either objective critria or results object
#' @examples
fx.hbv.model.parallel<-function(prm=NULL,
                                input, #data
                                Snow2GlacierOption = TRUE,#convert snow tro glacier at end of hydrological year
                                Snow2GlacierJulianDate= as.Date("2000-09-01"), #For calc of Snow2GlacierJulianday
                                calibrate=TRUE, # makes te dunction return obj.crit if TRUE to enable optimisation
                                obj.criteria="NSE", # chhose an objective criteria),
                                icall=1,
                                hyears=NULL, #provide hydrological years to run. if NULL all are used
                                calib.plots=FALSE
){
  # hbv.model.parallel simulations (each year simulated alone) ####
  library(dplyr)
  if(is.null(hyears)){hyears<-input$data$hyear} else {hyears=hyears}
  nhyr<-length(hyears)
  out<-vector("list",length = nhyr)
  names(out)<-hyears


 par(mfrow=c(ceiling(nhyr),1))

  for(i in 1:nhyr){
    hyr.input<-input
    hyr.input$data$timeseries<-input$data$timeseries %>% filter(hyear==hyears[i]) #susbset data
    hyr.input$data$hyear=hyears[i]
    hyr.input$data$serie.start.date = hyr.input$data$timeseries$Date[1]
    out[[i]]<-fx.hbv.model.serial(prm=prm,
                                  input=hyr.input,
                                  Snow2GlacierOption = Snow2GlacierOption,
                                  Snow2GlacierJulianDate= Snow2GlacierJulianDate,
                                  calibrate=calibrate,
                                  obj.criteria = obj.criteria,
                                  icall = 1,
                                  calib.plots = calib.plots)
  }
  par(mfrow=c(1,1))
  obj1<-lapply(out,function(x) x["objective.criteria"])
  objective.criteria<-do.call(rbind,lapply(obj1,function(x)x[[1]]))
  objdf<-apply(objective.criteria,2,function(x) sprintf("%5.2f",x))
  i.return<-which(names(objective.criteria) ==obj.criteria)
  OBJ<-objective.criteria[,i.return]
  OBJ<-mean(OBJ)


  if(calibrate==TRUE){

    return(OBJ)

  } else {
    # to be developed
    title(main=paste("\n[ Average ", obj.criteria," = ",sprintf("%5.2f",OBJ),"]"),outer=TRUE)
    return(out)

  }

}

#' Plot Timesereies Qobs and Qsim and Prec Using ggplot
#' @import zoo
#' @param Qsim Simulated flow - same length as Qobs, Prec, Time
#' @param Qobs Observed flow
#' @param Prec Precipitaion
#' @param Time Time as.Date("yyyy-dd-mm")
#' @param las x-axosd label orientation (unused)
#' @examples
fx.ggplot.Qsim.Qobs.Time<-function(Qsim=NULL,Qobs=NULL,Time=NULL,Prec=NULL,las=2){
  # ggplot of Qsim and Qobs with Time***********************************************************************************
  # Load the required library
  library(ggplot2)

  #check nulls
  if(is.null(Time)){Time<-1:length(Qobs)}
  if(is.null(Prec)){Prec<-rep(NA,length(Qobs))}

  # Create a data frame for timeseries
  data_timeseries <- data.frame(Time,Qobs,Qsim,Prec)

  # Melt the data frame to long format
  data2plot <- reshape2::melt(data_timeseries, id.vars = "Time")

  # Create the timeseries plot
  p <- ggplot(data = data2plot[data2plot$variable != "Prec", ], aes(x = Time, y = value, color = variable)) +
    geom_line() +
    scale_color_manual(values = c("Qobs" = "blue", "Qsim" = "red","Prec" ="#0000FF1A")) +
    ylab("Flow (m³/s)") +
    theme_bw()

  # Add the inverted bar graph
  p <- p + geom_bar(data = data2plot[data2plot$variable == "Prec", ], aes(x = Time, y = -value, fill = variable), stat = "identity", show.legend = FALSE) +
    scale_fill_manual(values = c("Qobs" = "transparent", "Qsim" = "transparent","Prec" = "#0000FF1A"))  # Set the fill color for Prec

  # Add secondary axis
  max_value <- max(data_timeseries$Qobs, data_timeseries$Qsim)
  sec.axis.trans <- ~ . * -1
  p <- p + scale_y_continuous(sec.axis = sec_axis(trans = sec.axis.trans, name = "Precipitation (mm)"))
  p<-p+ theme(legend.position="bottom")

  # Print the plot
  print(p)
}

#' Compute Potential evaporation based on Temperature, latitide (degree) and day of the year
#' Suggested that the first column is an index or a date n the form yyy-mm-dd H:M
#' handles multiple plots in one graph
#' enables use of secondary y-axis
#' you can define which series go to the second or first y axis by wrting the clumn number
#' @import zoo
#' @import RSEIS
#' @param Ta air-temperature numerical vector<7sere
#' @param latdeg latitude in degrees
#' @param strdates dates in the format as.Date("YYYY-MM-DD")
#' @return  PETts PE sereis in mm/d.
fx.OudinPE<-function(Ta=c(22,	23,	22,	22,	22	,21	,21	,21	,21	,21,	21,	21), #kampala
                     latdeg=0.3476,
                     strdates=seq(as.Date("1990-01-01"),
                                  as.Date("1990-12-01"),
                                  by="months")){
  # Potential Evaporation (PET) from Temperature, Latitude and Date ***************************************************
  # #lat<-1 #latitude
  #Ta<-(sin(1:12)+19) #temperature
  Ra<-Ta;Ra[]<-NA;Ret<-Ra
  da<-substr(strdates,9,10)
  mons<-substr(strdates,6,7)
  yrs<-substr(strdates,1,4)
  if(!class(strdates)=="Date"){day<-strdates}else {day<-as.Date(strdates)}
  # deg to radians for latitude
  j=pi/180*latdeg
  #MJ/m2 min Solar constant
  Gsc=0.0820
  #Julian Day or number of the day in the year between 1 (1 January) and 365 or 366 (31 December)
  library(RSEIS)
  J1<-tojul(as.numeric(yrs),1,1) #january 1st
  J<-tojul(as.numeric(yrs),as.numeric(mons),
           as.numeric(da))-J1
  #inverse relative distance Earth-Sun, dr,
  dr<-1+0.033*cos(2*pi/365*J)
  #solar declination, d
  decl<-0.409*sin(2*pi/365*J-1.39)
  #sunset hour angle, w s
  #ws = arccos(-tan (j)*tan (decl))
  X<-(1-(tan(j))^2 *(tan(decl))^2)
  X[X<0.00001]<-0.00001
  ws<-pi/2-atan(-tan(j)*tan(decl)/(X^0.5))
  #MJ /m2 day Extra Terrestrial Radiation
  Ra=24*60/pi*Gsc*dr**(sin(j)*sin(decl)+cos(j)*cos(decl)*sin(ws)) #[MJ m-2 d-1]
  Ret<-Ra*1e6/(24*3600) #J/m2 s-1 or W/m2
  #Convert to mm/d
  alpha<-2.45e6
  PE<-1000*24*3600 *(Ret/(alpha*1000))*(Ta+5)/(100)# "mm/d"
  fun<-function(x){mean(x,na.rm=T)}
  library(zoo)
  PETts<-zoo(PE,day) #Timeseries
  plot(PETts,main="Oudin Temperature based PE",ylab="mm/d",xlab="Date",las=2) #plot
}

#' color ramp with transparency alpha
#' @param ... a string with colors or colornames
#' @param n number of colors
#' @param alpha transparency of colors
colorRampAlpha <- function(..., n, alpha) {
  # Color Ramp with transparency(alpha) ####
  colors <- colorRampPalette(...)(n)
  paste(colors, sprintf("%x", ceiling(255 * alpha)),
        sep = "")
}

#' Terrain color ramp with transparency alpha for n bins
#' Uses colorRampAlpha function
#' @param col a string with colors or colornames
#' @param n number of breaks/colors in terrain bins (reults in n-1 colors)
#' @param alpha transparency of colors
fx.terrain.colors<-function(col=c("blue", "cadetblue1",
                                  "yellowgreen", "darkgreen", "gold3",
                                  "firebrick", "darkred", "chocolate4",
                                  "gray81", "white"), n=32,alpha=0.35){
# Terrain color ramp with transparency alpha ***************************************************************************
terrain.kol <- colorRampAlpha(col, n = n-1, alpha = 0.35)
return(terrain.kol )
}

#' Delineate catchment(s) using whitebox package given outlet point shapefile or coordinates and a DEM (tif)
#' Depends on the whitebox package whgoch must be installed
#' Path to white box must be provided
#' coordinates for ctahcment outlet must be given as a data.frame with x (longitude or xUTM), y (latitude or yUTM) and epsg (the projection) OR as an already projected point file
#' @param x outlet points (path to shapefile or to csv with a dataframe havig x,y,epsg )
#' @param dem path to tif file with DEM. Should preferably be same epsg as x
#' @param wbt path to whitebox tools
#' @import whitebox
#' @import raster
#' @import terra
#' @import sf
fx.delineate.catchment<-function(x){
   #Delineate a catchment**********************************************************************************************

}

#READ AND PLOT GIS FEATURES FROM SHAPEFILES OR GEODATABASES############################################################

#' This function reads and plots shapefiles and geodatabases
#' @param folder Folder containing preferable only the esri-shapefiles and/or geodatabase.
#' @param use.gdb Format to read. TRUE for gdb and FALSE for shapefiles
#' @param cols Colors for painting shapefiles
#' @param borders colors for shapefile borders
#' @param custom.epsg epsg to  use for unprojected files. Default utm32(epsg=25832).More information see https://epsg.io/
#' @param use.custom.epsg TRUE Use custom.epsg for unprojected shapes, FALSE use custom.epsg only if all shapes are unprojected. But if not, then the projection is taken from the projected files.
#' @param  lokale.nb TRUE if we want to set the lokale to Norwegian for Pcs struggling with reading special characters in paths to files
#' @param return.default TRUE for "default" ; FALSE for "map"
#' @return A List with the sf objects and the map (default) or just a mapview map (map)
fx.read.mapview.shp.gdb<-function(folder,
                                  use.gdb=FALSE,
                                  borders=c("red","pink","blue","lightblue","black","yellow","orange","chartreuse1","darkolivegreen1","green","darkgreen","palegreen"),
                                  cols=c("red","pink","blue","lightblue","black","yellow","orange","chartreuse1","darkolivegreen1","green","darkgreen","palegreen"),
                                  custom.epsg=25832,
                                  use.custom.epsg=FALSE,
                                  lokale.nb=TRUE,
                                  return.default=TRUE
){

  if(lokale.nb==TRUE){
    #LOKALE*************************************************************************************************************
    #attempt to solve the issue with not being able to read the paths conatiing norwegian characters
    Sys.setlocale("LC_ALL", "nb-NO.UTF-8")
    Sys.getlocale() #check lokale

    #To make the locale permanent in RStudio, you can add the following line to your ~/.Rprofile file:
    #candidates <- c(
    # Sys.getenv("R_PROFILE"),
    # file.path(Sys.getenv("R_HOME"), "etc", "Rprofile.site"),
    # Sys.getenv("R_PROFILE_USER"),
    # file.path(getwd(), ".Rprofile"),
    # file.path(Sys.getenv("HOME"), ".Rprofile")
    # )
    #Filter(file.exists, candidates)
    Sys.setlocale("LC_ALL", "nb-NO.UTF-8")
    Sys.getlocale() #check lokale
  }

  #PACKAGES*************************************************************************************************************
  library(sf)
  library(kableExtra)
  library(mapview)
  mapviewOptions(viewer.suppress = TRUE)
  library(leaflet)
  library(leaflet.extras)
  library(tcltk2)
  library(crayon)


  #BACKGROUND MAPS******************************************************************************************************
  urlg<-c("http://services.geodataonline.no/arcgis/services/Geocache_UTM33_WGS84/GeocacheBasis/MapServer/WMSServer",
          "http://services.geodataonline.no/arcgis/services/Geocache_UTM33_WGS84/GeocacheBilder/MapServer/WMSServer",
          "http://services.geodataonline.no/arcgis/services/Geocache_UTM33_WGS84/GeocacheGraatone/MapServer/WMSServer",
          "http://services.geodataonline.no/arcgis/services/Geocache_UTM33_WGS84/GeocacheLandskap/MapServer/WMSServer",
          "http://services.geodataonline.no/arcgis/services/Geocache_UTM33_WGS84/GeocacheHybrid/MapServer/WMSServer")

  #PROJECTIONS**********************************************************************************************************
  #possible projections to use when data doesnt have a projection
  utm32<-st_crs(25832)
  utm33<-st_crs(25833)

  # SHAPEFILED#PACKAGES*************************************************************************************************************
  if(use.gdb==FALSE){
    #paths to shapefiles
    shps<-list.files(folder, pattern = ".shp",full.names = TRUE)
    shps<-shps[!grepl(".xml",shps)]
    #print(shps)

    #Names of shapefiles
    layer_names<-unlist(lapply(strsplit(shps,"/"),function(x) x[length(x)]))
    #Remove special characters from the names
    #layer_names<-iconv(layer_names, from = 'latin1', to = 'ASCII//TRANSLIT')
    #print(layer_names)

    #shapefile list
    shps.list<-vector("list",length = length(layer_names))
    names(shps.list)<-layer_names
    #coordinate system list
    shps.crs<-shps.list
    #color.list
    col.list<-shps.list
    borders.list<-rep(NA,length(layer_names))

    #loop through the files
    for( i in 1:length(layer_names)){
      shps.list[[i]]<-sf::st_read(shps[i],quiet = TRUE) #read shapefile
      shps.list[[i]]<-st_zm( shps.list[[i]])
      #print(paste("Coordinate system:", layer_names[i]))
      shps.crs[[i]]<-st_crs(shps.list[[i]]) #save coordinate system
      num.shapes<-nrow(shps.list[[i]])
      col.list[[i]]<-cols[i]
      borders.list[i]<-borders[i]
    }

    #provide projection to unprojected assuming projection is the same for all
    #assumes at least one of the shapes are projected
    isNA.proj<-which(is.na(lapply(shps.crs,function(x)x$'Coordinate Reference System')))
    #stopifnot("The shapefiles have no projection"=length(isNA.proj)<length(layer_names))

    if(length(isNA.proj)<length(layer_names)){
      if(length(isNA.proj)>0){
        i.proj<-which(!is.na(lapply(shps.crs,function(x)x$'Coordinate Reference System')))[1]
        if(use.custom.epsg==TRUE){pjn<-custom.epsg} else {pjn<-shps.crs[[i.proj]]}
        for( k in 1:length(isNA.proj)){
          st_crs(shps.list[[isNA.proj[k]]])<-pjn
          st_write(obj = shps.list[[isNA.proj[k]]],
                   dsn = paste0(folder,"Reprojected_",layer_names[isNA.proj[k]]),
                   delete_layer=TRUE)
        }

      }
    } else {
      cat(red("The shapefiles have no projection!!!"))
      #prj<-custom.epsg#c("utm32","utm33")
      #w1 <- tktoplevel()
      #use.prj<-tk_select.list(choices=prj)
      #tkwait.window(w1)
      #tkdestroy(w1)
      #cat(blue("Reprojecting to: ",use.prj))
      #if(use.prj=="utm32"){pjn=utm32}
      #if(use.prj=="utm33"){pjn=utm33}
      pjn<-custom.epsg

      for( k in 1:length(isNA.proj)){
        st_crs(shps.list[[isNA.proj[k]]])<-pjn
        st_write(obj = shps.list[[isNA.proj[k]]],
                 dsn = paste0(folder,"Reprojected_",layer_names[isNA.proj[k]]),
                 delete_layer=TRUE)
      }

    }

  }

  #GEODATABASE**********************************************************************************************************
  if(use.gdb==TRUE){
    #paths to gdb
    gdb_path<-list.files(folder, pattern = ".gdb",full.names = TRUE)
    #print(gdb_path)

    #list the layres
    layer_info <- st_layers(dsn = gdb_path)
    layer_names <- layer_info$name
    #print(layer_names)

    #shapefile list
    shps.list<-vector("list",length = length(layer_names))
    names(shps.list)<-layer_names
    #coordinate system list
    shps.crs<-shps.list
    #color.list
    col.list<-shps.list
    borders.list<-rep(NA,length(layer_names))

    #loop through the files
    for( i in 1:length(layer_names)){
      shps.list[[i]]<-st_read(gdb_path, layer = layer_names[i],quiet = TRUE) #read shapefile
      shps.list[[i]]<-st_zm( shps.list[[i]])
      #print(paste("Coordinate system:", layer_names[i]))
      shps.crs[[i]]<-st_crs(shps.list[[i]]) #save coordinate system
      num.shapes<-nrow(shps.list[[i]])
      col.list[[i]]<-cols[i]
      borders.list[i]<-borders[i]
    }

    #remove empty layers
    n.f<-lapply(shps.list, function(x) nrow(x)) # number of shapes in each file
    #print(t(as.data.frame(n.f)))
    i.keep<-which(n.f>0)
    shps.list<-shps.list[i.keep]
    layer_names<-layer_names[i.keep]
    col.list<-col.list[i.keep]
    borders.list<-borders.list[i.keep]

  }


  #MAPVIEW**************************************************************************************************************
  #mapview
  m<-mapview(shps.list,col.regions=col.list,color=borders.list,alpha.regions = 0.5)
  m@map<-m@map %>%
    addFullscreenControl()  %>%
    addWMSTiles(
      baseUrl=urlg[3],
      layers = "0",
      options = WMSTileOptions(format = "image/png", transparent = TRUE),
      group="GeocacheBasisGraatone"
    )%>%
    addWMSTiles(
      baseUrl=urlg[1],
      layers = "0",
      options = WMSTileOptions(format = "image/png", transparent = TRUE),
      group="GeocacheBasis"
    )%>%
    addFullscreenControl() %>%
    #addTiles %>%
    addLayersControl( #"OpenStreetMap", "Stamen.Toner","Stamen.Terrain", "Esri.WorldStreetMap", "Wikimedia", "CartoDB.Positron", "Esri.WorldImagery"
      baseGroups = c("GeocacheBasis","GeocacheBasisGraatone"),
      overlayGroups = c(layer_names))%>%
    addLegend("bottomleft",
              colors = col.list,
              title = "Tegnforklaring",
              labels = layer_names,
              opacity = c(0.2,1,0.1,0.2,0.2,1)
    )

  m@map #show the map

  if(return.default==TRUE){
    out<-list(map=m,data=shps.list)
    return(out)

  } else
  {
    return(m)
  }

}


# FLOW CHARTS ##########################################################################################################
# fx.draw.system
#
# This functions draws a flow chart based on data in an xls file
# Nodes are in a sheet call Nodes. case sensitive
# Edges are in a sheet called edges. case sensitive

#' This functions draws a flow chart based on data in an xls file
#' @import DiagrammeR
#' @import readxl
#' @import tidyverse
#' @param nodes_edges_xls path to xls file with niods and edges. temlate is provided
#' @param rankdir rankdirection. can be "LR" or "TB".
#' @param bgcolor background color
#' @return a Grphiz graph obkect.
#' @examples
#' fx.draw.system(system.file(system.file("data/Flowchart_Nodes_and_edges.xlsx",package = "ejjhydrotools")))
fx.draw.system<-function(nodes_edges_xls,
                         rankdir="LR",
                         bgcolor="#8080800D"){
  library(DiagrammeR)
  #https://rich-iannone.github.io/DiagrammeR/articles/graphviz-mermaid.html
  #attributes
  #https://graphviz.org/docs/attrs/dir/
  #graph https://graphviz.org/docs/graph/
  #
  library(readxl) #https://readxl.tidyverse.org/
  library(tidyverse)
  Nodes<-read_excel(nodes_edges_xls, sheet = "Nodes")
  Edges<-read_excel(nodes_edges_xls, sheet = "Edges")
  #creating a node data frame
  # N=29 #29
  # Nodes<-Nodes[1:N,] ; print(Nodes,n=N);
  # E=27 #27
  # Edges<-Edges[1:E,] ; print(Edges,n=E)
  nodes<- create_node_df(n= nrow(Nodes),
                         style= Nodes$style,
                         label= Nodes$label,
                         color= Nodes$color,
                         shape= Nodes$shape,
                         fillcolor=Nodes$fillcolor,
                         width=Nodes$width,
                         height=Nodes$height,
                         fixedsize=Nodes$fixedsize,#stretch to fit the words
                         fontcolor=Nodes$fontcolor,
                         fontsize= Nodes$fontsize)
  # data frame of edges
  edges<-create_edge_df(from = Edges$from,
                        to=Edges$to,
                        label = Edges$label,
                        rel= Edges$rel,
                        color= Edges$color,
                        fontcolor=Edges$fontcolor,
                        style=Edges$style,
                        arrowhead=Edges$arrowhead,
                        dir=Edges$dir,
                        fontsize=Edges$fontsize
  )
  graph<-create_graph(nodes, edges) %>%
    add_global_graph_attrs(
      attr = c("layout", "rankdir", "splines","bgcolor"),
      value = c("dot", rankdir, "false",bgcolor),
      attr_type = c("graph", "graph", "graph","graph"))
  graf<-render_graph(graph)
  return(graf)
}

#' fx.draw.system.clusters
#'
#' This functions draws a flow chart based on data in an xls file
#' Includes option for draawing clusters
#' cluster data is also in the xls file

#' This functions draws a flow chart based on data in an xls file
#' @import DiagrammeR
#' @import readxl
#' @import tidyverse
#' @param nodes_edges_xls path to xls file with niods and edges. temlate is provided
#' @param rankdir rankdirection. can be "LR" or "TB".
#' @param bgcolor background color
#' @param splines Determines how the nodes are connected. can be "true","ortho","line","none","curved","spline","polyline"
#' @return a Grphiz graph obkect.
#' @examples
#' fx.draw.system(system.file(system.file("data/Flowchart_Nodes_and_edges.xlsx",package = "ejjhydrotools")))
fx.draw.system.clusters<-function(nodes_edges_xls,
                                  rankdir="LR",
                                  bgcolor="#8080800D",
                                  splines="true"
                                  #ortho,line,none,curved,spline,polyline
){
  #library(DiagrammeR)
  #https://rich-iannone.github.io/DiagrammeR/articles/graphviz-mermaid.html
  #attributes
  #https://graphviz.org/docs/attrs/dir/
  #graph https://graphviz.org/docs/graph/
  #
  library(readxl) #https://readxl.tidyverse.org/
  library(tidyverse)
  Nodes<-read_excel(nodes_edges_xls, sheet = "Nodes")
  Edges<-read_excel(nodes_edges_xls, sheet = "Edges")
  #creating a node data frame
  # N=29 #29
  # Nodes<-Nodes[1:N,] ; print(Nodes,n=N);
  # E=27 #27
  # Edges<-Edges[1:E,] ; print(Edges,n=E)
  nodes<- create_node_df(n= nrow(Nodes),
                         style= Nodes$style,
                         label= Nodes$label,
                         color= Nodes$color,
                         shape= Nodes$shape,
                         fillcolor=Nodes$fillcolor,
                         width=Nodes$width,
                         height=Nodes$height,
                         fixedsize=Nodes$fixedsize,#stretch to fit the words
                         fontcolor=Nodes$fontcolor,
                         fontsize= Nodes$fontsize,
                         cluster= Nodes$cluster)
  # data frame of edges
  edges<-create_edge_df(from = Edges$from,
                        to=Edges$to,
                        label = Edges$label,
                        rel= Edges$rel,
                        color= Edges$color,
                        fontcolor=Edges$fontcolor,
                        style=Edges$style,
                        arrowhead=Edges$arrowhead,
                        dir=Edges$dir,
                        fontsize=Edges$fontsize,
                        penwidth=Edges$penwidth,
                        splines=Edges$splines
  )
  graph<-create_graph(nodes, edges) %>%
    add_global_graph_attrs(
      attr = c("layout", "rankdir", "splines","bgcolor","splines"),
      value = c("dot", rankdir, "false",bgcolor,splines),
      attr_type = c("graph", "graph", "graph","graph","graph"))
  graf<-render_graph(graph)
  return(graf)
}

#INTERACTIVE PLOTS WITH HIGHCHART######################################################################################

#' Plot Timesereies Qobs and Qsim and Prec Using highcharter interactive plot
#' @import zoo
#' @param Qsim Simulated flow - same length as Qobs, Prec, Time
#' @param Qobs Observed flow
#' @param Prec Precipitaion
#' @param Time Time as.Date("yyyy-dd-mm")
#' @param las x-axosd label orientation (unused)
#' @param xts logical converts data to xts te enable some zooming possibilities
#' @param theme a number from 1 to 3 . chooses between three themes
#' @examples
fx.highchart.hbv<-function(Qsim=NULL,Qobs=NULL,Time=NULL,Prec=NULL,las=2,xts=TRUE,theme=1){
  # Interactive of Qsim and Qobs with Time*****************************************************************************
  ## #https://www.r-bloggers.com/2017/02/hyetographs-hydrographs-and-highcharter/
  #https://jkunst.com/highcharter/reference/hc_add_yAxis.html
  #best help
  #https://rkabacoff.github.io/datavis/Interactive.html
  library(highcharter)

  library(dplyr)
  library(reshape2)
  library(xts)

  # Set highcharter options and themes
  #'a theme for hc
  thm <- hc_theme(
    colors = c("red", "green", "blue"),
    chart = list(
      backgroundColor =  "#E6E6E61A"
    ),
    title = list(
      style = list(
        color = "#333333",
        fontFamily = "Erica One"
      )
    ),
    subtitle = list(
      style = list(
        color = "#666666",
        fontFamily = "Shadows Into Light"
      )
    ),
    legend = list(
      itemStyle = list(
        fontFamily = "Tangerine",
        color = "black"
      ),
      itemHoverStyle = list(
        color = "gray"
      )
    ),
    tooltip = list(valueDecimals = 2)
  )
  if(theme==1){options(highcharter.theme = thm )}
  if(theme==2){options(highcharter.theme = hc_theme_ggplot2(tooltip = list(valueDecimals = 2)))}
  if(theme==3){options(highcharter.theme = hc_theme_handdrawn(tooltip = list(valueDecimals = 2)))}


  #check nulls
  if(is.null(Time)){Time<-1:length(Qobs)}
  if(is.null(Prec)){Prec<-rep(NA,length(Qobs))}

  # Create a data frame for timeseries
  data_timeseries <- data.frame(Time,Qobs,Qsim,Prec)

  #round to 2 dp
  #data_timeseries[,-1]<-apply(data_timeseries[,-1],2,function(x) round(x,2))

  #ylim
  qmaxlim<-1.5*max(c(data_timeseries$Qobs,data_timeseries$Qsim),na.rm=TRUE)
  rmaxlim<-1.5*max(c(data_timeseries$Qobs,data_timeseries$Prec),na.rm=TRUE)

  #xts - needed for some types of zooming
  if(xts==TRUE) {
    data2plot <- as.xts( data_timeseries,order_by=data_timeseries$Time )
  } else {
    data2plot<-data_timeseries
  }

  if(xts==TRUE) {hc <- highchart(type = "stock")} else {hc <- highchart()}
  hc <- hc %>%
    hc_title(text = "HBV") %>%
    #hc_plotOptions(series = list(marker = list(enabled = FALSE))) %>%
    hc_exporting(enabled = TRUE) %>%
    hc_chart(zoomType = "xy") %>%
    # hc_yAxis_multiples(
    #   list(lineWidth = 3,
    #        title = list(text = "Flow (m³/s)"),
    #        max= qmaxlim),
    #   list(title = list(text = "Precipitation (mm)"),
    #        showLastLabel = FALSE,
    #        opposite = TRUE,
    #        reversed = TRUE,
    #        max= rmaxlim)
    # ) %>%
    hc_yAxis_multiples(list(title = list(text = "Flow (m³/s)"),
                            opposite = FALSE,
                            max= qmaxlim),
                       list(showLastLabel = FALSE,
                            opposite = TRUE,
                            reversed = TRUE,
                            max= rmaxlim,
                            title = list(text = "Precipitation (mm)"))
    )%>%
    hc_add_series(data = data2plot$Qobs,color="blue", type = "spline",name="Qobs") %>%
    hc_add_series(data = data2plot$Qsim,color="red", type = "spline",name="Qsim")%>%
    hc_add_series(data = data2plot$Prec, type = "column", yAxis = 1,color="#0000FF1A",name="Prec")%>%
    hc_xAxis(categories = Time, title = list(text = "date")) %>%
    hc_tooltip(shared = TRUE) %>%
    hc_rangeSelector(selected = 4)
  return(hc)
}



# Load necessary libraries
library(highcharter)
#' Plot Timesereies/series in a dataframe sing highcharter interactive plot
#' Suggested that the first column is an index or a date n the form yyy-mm-dd H:M
#' handles multiple plots in one graph
#' enables use of secondary y-axis
#' you can define which series go to the second or first y axis by wrting the clumn number
#' @import zoo
#' @import xts
#' @import highcharter
#' @import dplyr
#' @import reshape2
#' @import lubridate
#' @param data_timeseries dataframe with timeseries data
#' @param x  column for x-axis data
#' @param y0 columns for y-axis data (defaukt 2)
#' @param y1  columns for secondary y-axis data (default 34)
#' @param cols0 serie colorsfor  y0 c("lightblue","blue"),
#' @param cols1 serie colors for y1 c("orange","red") ,
#' @param type0 types for y0 c("spline","column"),
#' @param type1 types for y1 c("column","spline"),
#' @param opposite  opposite axes - c(FALSE,ALSE),
#' @param reversed  reversed c(FALSE, FALSE),
#' @param las axis text orientation (not impelmented)
#' @param xts logical for date polotting default is TRUE
#' @param theme appearance theme selection default=1
#' @paramy0lab   xaxis text default="Flow (m³/s)",
#' @paramy1lab yaxis text default="WSE (m)",#axis lbels
#' @param main Title- default is "Timeseries plot",
#' @param zoomType "zoomtype xy is defasult"
#' @param y0max factor for max value of ylim on primary y axis y0
#' @param y1max factor for max value of ylim on secondary y axis y1
#' @examples
#' fx.highchart.many.series.and.yaxes(data_timeseries=data_timeseries)
#'fx.highchart.many.series.and.yaxes(data_timeseries=data_timeseries ,#dataframe
#'                                 x=1 ,#column for x-axis data
#'                                 y0=2 ,#columns for y0-axis data
#'                                 y1=3:4 ,#columns for y1-axis data
#'                                 cols0=c("blue","cadetblue"),#serie colors y0
#'                                 cols1=c("red","magenta") ,#serie colors y1
#'                                 type0=c("spline","column") ,#types for y0
#'                                 type1=c("spline","spline") ,##types for y1
#'                                 opposite =c(FALSE,TRUE),
#'                                 reversed= c(FALSE, TRUE),
#'                                 las=2, #axis text orientation
#'                                 xts=TRUE, #for date polotting
#'                                 theme=1, #theme selction
#'                                y0lab = "Flow (m³/s)",
#'                                y1lab="WSE (m)",#axis lbels
#'                                 main="Timeseries plot",
#'                                 zoomType="xy")

# Universal function for plotting two series with dual y-axes
fx.highchart.many.series.and.yaxes<-function(data_timeseries,
                                             x=1,#column for x-axis data
                                             y0=2,y1=NA,#columns for y-axis data
                                             cols0=c("lightblue","blue"),#serie colors y0
                                             cols1=c("orange","red") ,#serie colors y1
                                             type0=c("spline","column"), #types for y0
                                             type1=c("spline","spline"), ##types for y1
                                             opposite =c(FALSE,TRUE),
                                             reversed= c(FALSE, TRUE),
                                             las=2, #axis text orientation
                                             xts=TRUE, #for date polotting
                                             theme=1, #theme selction
                                             y0lab = "y-Axis 0",
                                             y1lab=  "y-Axis 1)",#axis lbels
                                             main="Timeseries plot",
                                             zoomType="xy",
                                             y0max=1.1,
                                             y1max=1
){
  # Interactive plot of dataframe with multiple series and primary and secondary y-axes********************************
  #https://www.r-bloggers.com/2017/02/hyetographs-hydrographs-and-highcharter/
  #https://jkunst.com/highcharter/reference/hc_add_yAxis.html
  #best help
  #https://rkabacoff.github.io/datavis/Interactive.html
  library(highcharter)
  library(dplyr)
  library(reshape2)
  library(xts)
  library(lubridate)

  # # #debugging inputs
  # # data_timeseries=data_timeseries #dataframe
  # x=1 #column for x-axis data
  # y0=2:3 #columns for y0-axis data
  # y1=4 #columns for y1-axis data
  # cols0=c("blue","red")#serie colors y0
  # cols1=c("cyan","magenta") #serie colors y1
  # type0=c("spline","spline") #types for y0
  # type1=c("column","column") ##types for y1
  # opposite =c(FALSE,TRUE)
  # reversed= c(FALSE, TRUE)
  # las=2 #axis text orientation
  # xts=TRUE #for date polotting
  # theme=1 #theme selction
  #y0lab = "Flow (m³/s)"
  #y1lab="WSE (m)"#axis lbels
  # main="Timeseries plot"
  # zoomType="xy"
  #  y0max=1.0
  #  y1max=1.1


  # Set highcharter options and themes
  #'a theme for hc
  thm <- hc_theme(
    colors = c("red", "green", "blue"),
    chart = list(
      backgroundColor =  "#E6E6E61A"
    ),
    title = list(
      style = list(
        color = "#333333",
        fontFamily = "Calibri"
      )
    ),
    subtitle = list(
      style = list(
        color = "#666666",
        fontFamily = "Shadows Into Light"
      )
    ),
    legend = list(
      itemStyle = list(
        fontFamily = "Calibri",
        color = "black"
      ),
      itemHoverStyle = list(
        color = "gray"
      )
    ),
    tooltip = list(valueDecimals = 2)
  )
  if(theme==1){options(highcharter.theme = thm )}
  if(theme==2){options(highcharter.theme = hc_theme_ggplot2(tooltip = list(valueDecimals = 2)))}
  if(theme==3){options(highcharter.theme = hc_theme_handdrawn(tooltip = list(valueDecimals = 2)))}




  #ylim
  ylim0<-y0max*max(data_timeseries[,y0],na.rm=TRUE)
  if(is.na(y1)){ylim1<-ylim0} else {
    ylim1<-y1max*max(data_timeseries[,y1],na.rm=TRUE)
  }

  #xts - needed for some types of zooming
  if(xts==TRUE) {
    if(!is.Date(data_timeseries[,x])){
      data_timeseries[,x]<-as.POSIXct(data_timeseries[,x],tz="UTC")
    }
    data2plot <- zoo( data_timeseries[,-x],data_timeseries[,x])
    data2plot<-as.xts(data2plot)
    Time<-data_timeseries[,x]
    y0<-pmax(1,y0-1);y1<-pmax(1,y1-1)
  } else {
    data2plot<-data_timeseries
    Time<-data2plot[,x]
  }

  if(xts==TRUE) {hc <- highchart(type = "stock")} else {hc <- highchart()}
  hc <- hc %>%
    hc_title(text = main) %>%
    #hc_plotOptions(series = list(marker = list(enabled = FALSE))) %>%
    hc_exporting(enabled = TRUE) %>%
    hc_chart(zoomType = zoomType) %>%
    hc_yAxis_multiples(list(showLastLabel = TRUE,
                            title = list(text =y0lab),
                            opposite = opposite[1],
                            reversed = reversed[1],
                            max= ylim0,
                            lineColor="lightblue"
    ),
    list(showLastLabel = FALSE,
         title = list(text =y1lab),
         opposite = opposite[2],
         reversed = reversed[2],
         max= ylim1,
         lineColor="pink"
    )
    )
  for(i in 1:length(y0)){
    if(type0[i]=="column"){
      marker0= list(symbol = "square",enabled = TRUE, fillColor = cols0[i])
    } else{
      marker0 =list(enabled = TRUE, fillColor = cols0[i])
    }
    hc <- hc %>%hc_add_series(data = data2plot[,y0[i]],
                              color=cols0[i],
                              type = type0[i],
                              name=names(data2plot)[y0][i],
                              yAxis = 0,
                              marker = marker0, #list(enabled = TRUE, fillColor = cols0[i])
                              showInLegend = TRUE
    )
  }

  #add secondary axis data
  if(!is.na(y1)){
    for(i in 1:length(y1)){
      if(type1[i]=="column"){
        marker1= list(symbol = "square",enabled = TRUE, fillColor = cols1[i])} else{
          marker1 =list(enabled = TRUE, fillColor = cols1[i])
        }
      hc <- hc %>%hc_add_series(data = data2plot[,y1[i]],
                                color=cols1[i],
                                type = type1[i],
                                name=names(data2plot)[y1][i],
                                yAxis = 1,
                                marker = marker1,
                                showInLegend = TRUE
      )
    }
  }

  hc<-hc %>%
    hc_xAxis(categories = Time, title = list(text = "date"),dateTimeLabelFormats = list(day = '%Y-%m-%d %H:%S'), type = "datetime") %>%
    hc_tooltip(
      shared = TRUE,
      split = FALSE,
      valueDecimals = 2,
      sort = FALSE,
      table = FALSE#TRUE colors the names but takes away the markers for the default
    ) %>%
    hc_rangeSelector(selected = 5)


  # Symbols
  # 1."circle" ●
  # 2."diamond" ♦
  # 3."square" ■
  # 4."triangle" ▲
  # 5."triangle-down" ▼


  hc<-hc %>% hc_tooltip(
    shared = TRUE,
    split = FALSE,
    valueDecimals = 2,
    sort = FALSE,
    table = FALSE,#TRUE colors the names but takes away the markers
    pointFormat = "<span style=\"color:{series.color}\">{series.name}</span>:
             <b>{point.percentage:.1f}%</b> ({point.y:,.0f} millions)<br/>",
    formatter = JS("
                   function () {
                   var tooltip = '<b>' + Highcharts.dateFormat('%Y-%m-%d %H:%M', this.x) + '</b><br/>';
                   $.each(this.points, function (i, point) {
                   let symbolMap = new Map();
                   symbolMap.set('circle',        '&#9679');
                   symbolMap.set('diamond',       '&#9670');
                   symbolMap.set('square',        '&#9632');
                   symbolMap.set('triangle',      '&#9650');
                   symbolMap.set('triangle-down', '&#9660');
                   var color = point.series.color;
                   var symbolName = point.series.symbol;
                   if(typeof point.series.symbol === 'undefined'){
                   tooltip += '<span style=\"color:' + color + ';\">' + symbolMap.get('square') + '</span><span style=\"color:' + color + ';\">: '+ point.series.name + '</b>: ' + point.y.toFixed(2) + '</span><br/>';
                   } else {
                   tooltip += '<span style=\"color:' + color + ';\">' + symbolMap.get(symbolName) + '</span><span style=\"color:' + color + ';\">: '+ point.series.name + '</b>: ' + point.y.toFixed(2) + '</span><br/>';
                   }

    });
    return tooltip;
}
"
    )
  )

  hc<-hc%>%
    hc_legend(enabled = TRUE) |>
    hc_legend(
      align = "left",
      verticalAlign = "middle",
      layout = "vertical",
      floating=TRUE, #FALSE puts it outside plot
      x = 90,
      y = 45,
      title = list(
        text = "Legend",
        style = list(
          textDecoration = "underline"
        )
      )
    )
  return(hc)
}
