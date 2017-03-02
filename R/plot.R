plotTemplate<-function(data,data.opt=NA, trace,xrange,yrange,xtitle,ytitle){
  p <-plot_ly()%>%
    add_trace(x = eval(parse(text=paste0("data$",as.name(trace$first$x)))),
              y = eval(parse(text=paste0("data$",as.name(trace$first$y)))),
              name = trace$first$name,text=~ifelse(is.na(trace$first$text),"", trace$first$text),
              type = 'scatter', mode = 'lines+markers'
    )%>%
    add_trace(x=eval(parse(text=paste0("data.opt$",as.name(trace$opt$x)))),
              y=eval(parse(text=paste0("data.opt$",as.name(trace$opt$y)))),
              type = 'scatter', mode = 'lines+markers',name=trace$opt$name,
              hoverinfo='skip y + x + text',text = ifelse(is.na(trace$opt$text),"", trace$opt$text)
    )%>%
    layout(
      xaxis = list(title = xtitle,range = ytitle),
      yaxis = list(title = ytitle,range = yrange)
      ,hovermode = 'x'
    )
  return(p)
}

discrimPwrCurve <-  function(pdA, pd0 = 0, alpha = 0.05, pGuess = 1/2, 
                             test = c("difference", "similarity"), 
                             statistic = c("exact", "normal", "cont.normal"),
                             SS=c("first","stable"))
{
  
  ss <- c(5:150)
  vecPower <- sapply(ss, function(s) discrimPwr(pdA=pdA, pd0 = pd0, sample.size=s,pGuess=pGuess))
  lag.vecPower <- c(tail(vecPower,-1),T)
  
  #Checking the growth of the zig-zag curve
  subsetPower=vecPower[vecPower>lag.vecPower]
  subsetSS=ss[vecPower>lag.vecPower]
  
  #Stable sample size
  stablePower=subsetPower[c(diff(subsetPower)>0,T)]
  stableSS=subsetSS[c(diff(subsetPower)>0,T)]
  
  #First exact sample size
  firstPower=subsetPower[c(T,diff(subsetPower)>0)]
  firstSS=subsetSS[c(T,diff(subsetPower)>0)]
  
  #Number of discriminator
  nb.discrim=pdA*subsetSS
  nb.discrimStable=nb.discrim[c(diff(subsetPower)>0,T)]
  nb.discrimFirst=nb.discrim[c(T,diff(subsetPower)>0)]

  #Number of correct answer
  nb.correct=pd2pc(pdA,pGuess)*subsetSS
  nb.correctStable=nb.correct[c(diff(subsetPower)>0,T)]
  nb.correctFirst=nb.correct[c(T,diff(subsetPower)>0)]
  
  #Data frame for plotting with plot_ly
  data.stable=cbind.data.frame(SS=stableSS,power=stablePower,nb.correct=nb.correctStable,nb.discrim=nb.discrimStable)
  data.first=cbind.data.frame(SS=firstSS,power=firstPower,nb.correct=nb.correctFirst,nb.discrim=nb.discrimFirst)
  data=cbind.data.frame(SS=ss,power=vecPower,nb.correct=pd2pc(pdA,pGuess)*ss,nb.discrim=pdA*ss)
  
  #Data to input to the plot_ly template                  ))
  input=NULL
  
  input$first$x=names(data)[1]
  input$first$y=names(data)[2]
  input$first$name='Sample size'
  input$first$text=NA
  
  data.opt=suppressWarnings(if (SS=="first") data.first else data.stable)
  input$opt$x=names(data.opt)[1]
  input$opt$y=names(data.opt)[2]
  input$opt$name=suppressWarnings(if (SS=="first") "First exact" else "Stable exact")
  input$opt$text=with(data.opt,paste('Power: ', round(power*100,1),'%','</br> Nb. discriminator:',nb.discrim,'</br> Nb. correct answer:',nb.correct))
  p=plotTemplate(data=data,data.opt=data.opt, trace=input,xrange=ss,yrange=c(0,1),xtitle='Sample size',ytitle='Power')

  return(plot=p)
}
