
#install.packages("minpack.lm")
library(nls2)
library(minpack.lm)

#Selecting the data
OILW13<-read.csv(file.choose(),header=T)

#Creating pdf to store graphs
pdf("Leon County.pdf")

count = 0
par(oma=c(0, 0, 0, 5))

#Creting table to store the Np and model errors
tab=data.frame(WEll_Number=0,Np_Actual_cuft=0,Np_Arps_cuft=0,Np_LGA_cuft=0,Np_Duong_cuft=0,Np_PLE_cuft=0,Error_ARPS=0,Error_LGA=0,Error_Duong=0,Error_PLE=0)


# Selecting unique wells from the county data 
for (api in unique(OILW13$API.UWI)){
 
  well = OILW13[OILW13$API.UWI == api,]
  rows = nrow(well)
  
  # Selecting only wells with 100 and over months of production
  if (nrow(well) > 100){
    if(count<50){
    count = count  + 1
    time_abs<-1:length(well$Monthly.Gas)
    plot(time_abs[1:100],well$Monthly.Gas[1:100], xlab="Time (month)", ylab="Production (cuft)",  main="Model Fit", sub=(toString(api)))
    time_train = time_abs[1:50]
    time_test = time_abs[51:100]
    gas_train = well$Monthly.Gas[1:50]
    
    #ARPS
    b_lower = 0.001
    b_upper = 20
    qi_lower  = 0
    qi_upper  = 1000000
    di_lower = 0.001
    di_upper = 100
    
    lower_list = c(qi_lower,b_lower,di_lower)
    upper_list = c(qi_upper,b_upper,di_upper)
    time = time_train
    trend1.mod<-nlsLM(gas_train~(q_i/(1+b*d_i * time)^(1/b)), start=list(q_i=40000,b=0.5, d_i=1),algorithm = "port",trace=T, control=list(maxiter=1024),lower=lower_list, upper=upper_list)
    summary(trend1.mod)
    fit1<-fitted.values(trend1.mod)
    time <- data.frame(time = time_abs[1:100])
    pred  <- predict(trend1.mod,newdata=time)
    lines(time_abs[1:100],pred,col="red")
    
    #LGA
    K_lower = 0.001
    K_upper = 1000000000000
    n_lower  = 0
    n_upper  = 1000000
    a_lower = 0.001
    a_upper = 10000
    
    lower_list = c(K_lower,n_lower,a_lower)
    upper_list = c(K_upper,n_upper,a_upper)
    time = time_train
    trend2.mod<-nlsLM(gas_train~(K*n*a*time^(n-1))/((a+time^n)^2), start=list(K=1000000,n=0.8, a=25),algorithm = "port",trace=T, control=list(maxiter=1024),lower=lower_list, upper=upper_list)
    summary(trend2.mod)
    fit2<-fitted.values(trend2.mod)
    time <- data.frame(time = time_abs[1:100])
    pred2  <- predict(trend2.mod,newdata=time)
    lines(time_abs[1:100],pred2,col="green")
    
    #Duong
    q1_lower = 0.001
    q1_upper = 1000000000000
    m_lower  = 0
    m_upper  = 1000000
    a_lower = 0.001
    a_upper = 10000
    
    lower_list = c(q1_lower,m_lower,a_lower)
    upper_list = c(q1_upper,m_upper,a_upper)
    time = time_train
    trend3.mod<-nlsLM(gas_train~(q1*time^(-m)*exp((a/(1-m))^(time^(1-m)-1))), start=list(q1=100000,m=0.8, a=1),algorithm = "port",trace=T, control=list(maxiter=1024),lower=lower_list, upper=upper_list)
    summary(trend3.mod)
    fit3<-fitted.values(trend3.mod)
    time <- data.frame(time = time_abs[1:100])
    pred3  <- predict(trend3.mod,newdata=time)
    lines(time_abs[1:100],pred3,col="blue")
    
    #Power Law Decline
    qi_lower = 10
    qi_upper = 1000000000000
    Dinf_lower  = 0.0000001
    Dinf_upper  = 100
    Di_lower = 0.001
    Di_upper = 100
    n_lower = 0.0001
    n_upper = 10
    
    lower_list = c(qi_lower,Dinf_lower,Di_lower,n_lower)
    upper_list = c(qi_upper,Dinf_upper,Di_upper,n_upper)
    time = time_train
    trend4.mod<-nlsLM(gas_train~(qi*exp(-Dinf*time-Di*time^n)), start=list(qi=100000,Dinf=0.001, Di=5, n=0.01),algorithm = "port",trace=T, control=list(maxiter=1024),lower=lower_list, upper=upper_list)
    summary(trend4.mod)
    fit4<-fitted.values(trend4.mod)
    time <- data.frame(time = time_abs[1:100])
    pred4  <- predict(trend4.mod,newdata=time)
    lines(time_abs[1:100],pred4,col="purple")
    
    legend(par('usr')[2], par('usr')[4],bty='n', xpd=NA, legend=c("Arps","LGA","Duong","PLE"), col=c("red","green","blue","purple"), lty=1)
    
    #cumm production, creating tables of Np for each well
    dataACT=data.frame( time_MNTHS =time_abs[1:100], Np_CuFt = cumsum(well$Monthly.Gas[1:100]*5.615))
    data_ARPS=data.frame( time_MNTHS =time_abs[1:100], Np_CuFt = cumsum(pred*5.615))
    data_LGA=data.frame( time_MNTHS =time_abs[1:100], Np_CuFt = cumsum(pred2*5.615))
    data_Duong=data.frame( time_MNTHS =time_abs[1:100], Np_CuFt = cumsum(pred3*5.615))
    data_PLE=data.frame( time_MNTHS =time_abs[1:100], Np_CuFt = cumsum(pred4*5.615))

    # Ploting Np of each model vs actual NP
    plot(dataACT, type="l", col="black")
    title(main= "Cummulative Production", sub=(toString(api)))
    lines(data_ARPS, col="red")
    lines(data_LGA, col="green")
    lines(data_Duong, col="blue")
    lines(data_PLE, col="purple")
    legend(par('usr')[2], par('usr')[4],bty='n', xpd=NA, legend=c("Actual", "Arps","LGA","Duong","PLE"), col=c("black", "red","green","blue","purple"), lty=1)
    
    
    #Calculating Errors of each model
    Error1=(abs(dataACT$Np_CuFt[100]-data_ARPS$Np_CuFt[100])/dataACT$Np_CuFt[100])*100
    Error2=(abs(dataACT$Np_CuFt[100]-data_LGA$Np_CuFt[100])/dataACT$Np_CuFt[100])*100
    Error3=(abs(dataACT$Np_CuFt[100]-data_Duong$Np_CuFt[100])/dataACT$Np_CuFt[100])*100
    Error4=(abs(dataACT$Np_CuFt[100]-data_PLE$Np_CuFt[100])/dataACT$Np_CuFt[100])*100
    
    #Creating tables of Np and error for each well
    data_Np=data.frame(WEll_Number=api,Np_Actual_cuft=dataACT$Np_CuFt[100],Np_Arps_cuft=data_ARPS$Np_CuFt[100],Np_LGA_cuft=data_LGA$Np_CuFt[100],Np_Duong_cuft=data_Duong$Np_CuFt[100],Np_PLE_cuft=data_PLE$Np_CuFt[100],Error_ARPS=Error1,Error_LGA=Error2,Error_Duong=Error3,Error_PLE=Error4)
    tab[count,]=data_Np

        }
  }
}

# Save the data in your computer
write.csv(tab,'C:\\Users\\ucegbe\\Desktop\\NAT GAS\\Leon County.csv', row.names = FALSE)
dev.off()

