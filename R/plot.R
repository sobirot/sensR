plotTemplate<-function(data,xtitle,ytitle,title,annotation){
  
  min.range=apply(matrix(sapply(data,min),2,dim(data)[2]),1,min)
  max.range=apply(matrix(sapply(data,max),2,dim(data)[2]),1,max)
  
  p <- plot_ly()%>%
    layout(title = title,
           xaxis = list(title = xtitle,range = c(0.95*min.range[1],1.05*max.range[1])),
           yaxis = list(title = ytitle,range = c(0.95*min.range[2],1.05*max.range[2]))
           ,hovermode = 'y',annotations=annotation
    )
  
  for (i in 1: dim(data)[2]){
    p <- add_trace(p,x=data[,i][[1]],y=data[,i][[2]],type = "scatter", mode = 'lines+markers',name =colnames(data)[i])
    
  }
  return(p)
}

power.SS=function(pdA, pd0, alpha, target.power,pGuess,sample.size,statistic,test){
  power <- seq(0.51,0.99,0.01) #power>0.5 for discrimSS
  ss <- sapply(power, function(s) {discrimSS(pdA=pdA, pd0 = pd0,alpha = alpha, target.power=s,pGuess=pGuess,statistic = statistic,test =test)})
  return(list(ss=ss,power=power))
}
power.ES=function(pdA, pd0, alpha, target.power,pGuess,sample.size,statistic,test){
  ES <- seq(0.01,0.99,0.01)
  power <- sapply(ES, function(s) {discrimPwr(pdA=s, pd0 = pd0,alpha = alpha, sample.size=sample.size,pGuess=pGuess,statistic = statistic,test =test)})
  return(list(ES=ES,power=power))
}
SS.ES=function(pdA, pd0, alpha, target.power,pGuess,sample.size,statistic,test){
  ES <- seq(0.51,0.99,0.01) #power>0.5 for discrimSS
  ss <- sapply(ES, function(s) {discrimSS(pdA=s, pd0 = pd0,alpha = alpha, target.power=target.power,pGuess=pGuess,statistic = statistic,test =test)})
  return(list(ES=ES,ss=ss))
}

discrimCalc <-  function(type=c("power-SS","power-ES","SS-ES"),display=c("protocol","alpha","ES","power","SS")
                         ,pdA=c(0.25,0.375,0.5), pd0 = 0, target.power = c(0.9,0.85,0.8), alpha = c(0.05,0.10,0.01), pGuess = c(1/3,1/2,1/9,1/4), sample.size=c(20,50,100), 
                         test = c("difference", "similarity"), 
                         statistic = c("exact","stable.exact", "normal", "cont.normal"))
{
  
  ## power-SS curve
  if (type =="power-SS"){
    if (display=="protocol"){
      #vec=pGuess
      data <- sapply(pGuess, function(s) power.SS(pdA=pdA, pd0=pd0, alpha = alpha, target.power=target.power,pGuess=s,sample.size=sample.size,statistic = statistic,test=test))
      colnames(data)=paste("pGuess = ",fractions(pGuess),sep="")
      text.anot=paste("pdA = ", round(pdA*100,1),"%<br>alpha = ",round(alpha*100,0), "%<br>statistic = ",statistic,ifelse(test=="similarity",paste("<br>pd0 = ", round(pd0*100,1), "%"),""))
      xtitle="Sample size"
      ytitle="Power"
      title=paste ("Power of",ifelse(test=='difference', 'difference',ifelse(test=='similarity','similarity')), "test")
    }else if (display=="alpha"){
      #vec=alpha
      data <- sapply(alpha, function(s) power.SS(pdA=pdA, pd0=pd0,alpha =s, target.power=target.power,pGuess=pGuess,sample.size=sample.size,statistic = statistic,test=test))
      colnames(data)=paste("alpha = ",round(alpha*100,0),"%",sep="")
      text.anot=paste("pdA = ", round(pdA*100,1),"%<br>pGuess = ",fractions(pGuess),"<br>statistic = ",statistic,ifelse(test=="similarity",paste("<br>pd0 = ", round(pd0*100,1), "%"),"") )
      xtitle="Sample size"
      ytitle="Power"
      title=paste ("Power of",ifelse(test=='difference', 'difference',ifelse(test=='similarity','similarity')), "test")
    }else if (display=="ES"){
      #vec=pdA
      data <- sapply(pdA, function(s) power.SS(pdA=s, pd0=pd0,alpha =alpha, target.power=target.power,pGuess=pGuess,sample.size=sample.size,statistic = statistic,test=test))
      colnames(data)=paste("pdA = ",round(pdA*100,1),"%",sep="")
      text.anot=paste("alpha = ", round(alpha*100,0)," %<br>pGuess = ",fractions(pGuess),"<br>statistic = ",statistic,ifelse(test=="similarity",paste("<br>pd0 = ", round(pd0*100,1), "%"),"") )
      xtitle="Sample size"
      ytitle="Power"
      title=paste ("Power of",ifelse(test=='difference', 'difference',ifelse(test=='similarity','similarity')), "test")
    }
    
    anot=list(
      text=text.anot,
      align = "left",
      bordercolor = "black",
      xref= "paper",
      yref="paper",
      x=1,
      y=0,
      xanchor = "right",
      yanchor = "bottom",
      showarrow=F
    )
  }
  else if (type =="power-ES"){
    if (display=="protocol"){
      #vec=pGuess
      data <- sapply(pGuess, function(s) power.ES(pdA=pdA, pd0=pd0, alpha = alpha, target.power=target.power,pGuess=s,sample.size=sample.size,statistic = statistic,test=test))
      colnames(data)=paste("pGuess = ",fractions(pGuess),sep="")
      text.anot=paste("SS = ", round(sample.size,0),"<br>alpha = ",round(alpha*100,0), "%<br>statistic = ",statistic,ifelse(test=="similarity",paste("<br>pd0 = ", round(pd0*100,1), "%"),""))
      xtitle="Effect size (pdA)"
      ytitle="Power"
      title=paste ("Power of",ifelse(test=='difference', 'difference',ifelse(test=='similarity','similarity')), "test")
    }else if (display=="alpha"){
      #vec=alpha
      data <- sapply(alpha, function(s) power.ES(pdA=pdA, pd0=pd0,alpha =s, target.power=target.power,pGuess=pGuess,sample.size=sample.size,statistic = statistic,test=test))
      colnames(data)=paste("alpha = ",round(alpha*100,0),"%",sep="")
      text.anot=paste("SS = ", round(sample.size,0),"<br>pGuess = ",fractions(pGuess),"<br>statistic = ",statistic,ifelse(test=="similarity",paste("<br>pd0 = ", round(pd0*100,1), "%"),"") )
      xtitle="Effect size (pdA)"
      ytitle="Power"
      title=paste ("Power of",ifelse(test=='difference', 'difference',ifelse(test=='similarity','similarity')), "test")
    }else if (display=="SS"){
      #vec=sample.size
      data <- sapply(sample.size, function(s) power.ES(pdA=pdA, pd0=pd0,alpha =alpha, target.power=target.power,pGuess=pGuess,sample.size=s,statistic = statistic,test=test))
      colnames(data)=paste("SS = ",round(sample.size,0),sep="")
      text.anot=paste("alpha = ", round(alpha*100,0),"%<br>pGuess = ",fractions(pGuess),"<br>statistic = ",statistic,ifelse(test=="similarity",paste("<br>pd0 = ", round(pd0*100,1), "%"),"") )
      xtitle="Effect size (pdA)"
      ytitle="Power"
      title=paste ("Power of",ifelse(test=='difference', 'difference',ifelse(test=='similarity','similarity')), "test")
    }
    anot=list(
      text=text.anot,
      align = "left",
      bordercolor = "black",
      xref= "paper",
      yref="paper",
      x=1,
      y=0,
      xanchor = "right",
      yanchor = "bottom",
      showarrow=F
    )
  }  else if (type =="SS-ES"){
    if (display=="protocol"){
      #vec=pGuess
      data <- sapply(pGuess, function(s) SS.ES(pdA=pdA, pd0=pd0, alpha = alpha, target.power=target.power,pGuess=s,sample.size=sample.size,statistic = statistic,test=test))
      colnames(data)=paste("pGuess = ",fractions(pGuess),sep="")
      text.anot=paste("Power = ", round(sample.size,0),"%<br>alpha = ",round(alpha*100,0), "%<br>statistic = ",statistic,ifelse(test=="similarity",paste("<br>pd0 = ", round(pd0*100,1), "%"),""))
      xtitle="Effect size (pdA)"
      ytitle="Sample Size"
      title=paste ("Sample size of",ifelse(test=='difference', 'difference',ifelse(test=='similarity','similarity')), "test")
    }else if (display=="alpha"){
      #vec=alpha
      data <- sapply(alpha, function(s) SS.ES(pdA=pdA, pd0=pd0,alpha =s, target.power=target.power,pGuess=pGuess,sample.size=sample.size,statistic = statistic,test=test))
      colnames(data)=paste("alpha = ",round(alpha*100,0),"%",sep="")
      text.anot=paste("Power = ", round(sample.size,0),"%<br>pGuess = ",fractions(pGuess),"<br>statistic = ",statistic,ifelse(test=="similarity",paste("<br>pd0 = ", round(pd0*100,1), "%"),"") )
      xtitle="Effect size (pdA)"
      ytitle="Sample Size"
      title=paste ("Sample size of",ifelse(test=='difference', 'difference',ifelse(test=='similarity','similarity')), "test")
    }else if (display=="power"){
      #vec=alpha
      data <- sapply(target.power, function(s) SS.ES(pdA=pdA, pd0=pd0,alpha =alpha, target.power=s,pGuess=pGuess,sample.size=sample.size,statistic = statistic,test=test))
      colnames(data)=paste("Power = ",round(target.power*100,0),"%",sep="")
      text.anot=paste("alpha = ", round(alpha*100,0),"%<br>pGuess = ",fractions(pGuess),"<br>statistic = ",statistic,ifelse(test=="similarity",paste("<br>pd0 = ", round(pd0*100,1), "%"),"") )
      xtitle="Effect size (pdA)"
      ytitle="Sample Size"
      title=paste ("Sample size of",ifelse(test=='difference', 'difference',ifelse(test=='similarity','similarity')), "test")
    }
    anot=list(
      text=text.anot,
      align = "left",
      bordercolor = "black",
      xref= "paper",
      yref="paper",
      x=1,
      y=1,
      xanchor = "right",
      yanchor = "top",
      showarrow=F
    )
  }
  return(list(data=data,xtitle=xtitle,ytitle=ytitle,title=title, anot=anot))
  
}

discrimPlot=function(type=c("power-SS","power-ES","SS-ES"),display=c("protocol","alpha","ES","power","SS")
                     ,pdA=c(0.25,0.375,0.5), pd0 = 0, target.power = c(0.9,0.85,0.8), alpha = c(0.05,0.10,0.01), pGuess = c(1/3,1/2,1/9,1/4), sample.size=c(20,50,100), 
                     test = c("difference", "similarity"), 
                     statistic = c("exact","stable.exact", "normal", "cont.normal")){
  test <- match.arg(test)
  stat <- match.arg(statistic)
  stopifnot(is.numeric(pdA)&& pdA >= 0 && pdA <= 1)
  stopifnot(is.numeric(pd0) && length(pd0) == 1 && pd0 >= 0 && pd0 <= 1)
  stopifnot(is.numeric(alpha) && alpha > 0 && alpha < 1)
  stopifnot(is.numeric(target.power) && target.power > 0 && target.power < 1)
  stopifnot(is.numeric(pGuess) && pGuess >= 0 && pGuess < 1)
  stopifnot(is.numeric(sample.size) && isTRUE(all.equal(round(sample.size), sample.size)) && sample.size > 0)
  sample.size <- as.integer(round(sample.size))
  
  if (type=="power-SS" && !display%in%c("protocol","alpha","ES") )
    stop("type = power-SS: display has be chosen among protocol, alpha and ES")
  if (type=="power-ES" && !display%in%c("protocol","alpha","SS") )
    stop("type = power-ES: display has be chosen among protocol, alpha and SS")
  if (type=="SS-ES" && !display%in%c("protocol","alpha","power") )
    stop("type = SS-ES: display has be chosen among protocol, alpha and power")
  
  if(display=="ES") pdA else pdA=pdA[1]
  if(display=="protocol") pGuess else pGuess=pGuess[1]
  if(display=="power") target.power else target.power=target.power[1]
  if(display=="alpha") alpha=alpha[order(alpha,decreasing = T)] else alpha=alpha[1]
  if(display=="SS") sample.size else sample.size=sample.size[1]
  statistic = ifelse(length(statistic)>1,statistic[1],statistic)
  test = ifelse(length(test)>1,test[1],test)
  
  if(type=="power-SS"&&test == "difference" && any(pdA <= pd0))
    stop("type = power-SS: pdA has to be larger than pd0 for difference tests")
  if(type=="power-SS"&&test == "similarity" && any(pdA >= pd0))
    stop("type = power-SS: pdA has to be less than pd0 for similarity tests")
  
  if(type=="power-ES"&&test == "difference" && any(pd0 > 0))
    stop("type = power-ES: 0 < pdA < 1,  pd0 has to be less than pdA for difference tests")
  if(type=="power-ES"&&test == "similarity" && any(pd0 < 1))
    stop("type = power-ES: 0 < pdA < 1, pd0 has to be larger than pdA for similarity tests")
  
  if(type=="SS-ES"&&test == "difference" && any(pd0 > 0.5))
    stop("type = SS-ES: 0.5 < pdA < 1,  pd0 has to be less than pdA for difference tests")
  if(type=="SS-ES"&&test == "similarity" && any(pd0 < 1))
    stop("type = SS-ES: 0.5 < pdA < 1, pd0 has to be larger than pdA for similarity tests")
  
  input=discrimCalc(type = type, display =display, pdA=pdA,pd0=pd0,target.power = target.power, alpha = alpha,pGuess = pGuess, sample.size = sample.size, test = test, statistic = statistic)
  plotTemplate(data=input$data,xtitle=input$xtitle,ytitle=input$ytitle,title=input$title,annotation=input$anot)
}