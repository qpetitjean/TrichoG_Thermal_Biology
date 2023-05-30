##########################
#       plot data        #  
##########################

##### visualization activity, nb detected ind, temperature

strains_plot <-
  function(
    finalDatList = finalDatList, # the cleaned dataset, needed to compute the average number of individuals over the video
    smoothDf = SmoothedResDf, # the smoothed dataset
    TempAxis = "y", # position of the temperature axis
    TempCol = "Temp", # name of the temperature column
    TempScale = NULL, # a vector containing 2 values (same unit as TempCol) indicating lower and upper limit of the temperature axis (either x or y)
    TimeCol = "runTimelinef", # name of the time column,
    TimeLimit = NULL, # a vector containing 2 values (same unit as TimeCol) indicating lower and upper limit of the x axis
    colTemp = NULL, # color of the temperature ramp
    nbrIndcounted = nbrIndcounted, # the number\r of individual manually counted
    IC = SplittedIC, # the dataset containing 95%Ci
    variables = NULL, # a vector containing the names of variables to plot (names of columns present in both smoothDf and IC)
    colvariable = NULL, # A vector containing the color of each variable on the plot
    scaleVariable = NULL, # a list containing vectors which indicates the upper and lower value to specify the scale for each variable
    Lab = NULL) # a list of label to add on the plot for each variable
{
    # in case TimeLimnit is specified, trunc the dataset to remove values below and above lower and upper limits respectively
    if(!is.null(TimeLimit)){ 
      # for finalDatList
      if(is.list(finalDatList)){
        finalDatList <- as.data.frame(finalDatList)
      }
      finalDatList <- finalDatList[which(finalDatList[[TimeCol]] >= TimeLimit[1] & finalDatList[[TimeCol]] <= TimeLimit[2]),]
      # for smoothDf
      smoothDf <- smoothDf[which(smoothDf[[TimeCol]] >= TimeLimit[1] & smoothDf[[TimeCol]] <= TimeLimit[2]),]
      # for IC
      IC <- lapply(IC, function(x) 
        x[which(x[[TimeCol]] >= TimeLimit[1] & x[[TimeCol]] <= TimeLimit[2]),]
        )
       }
    
    if(TempAxis == "x") {
      # build temperature x axis
      TempRange <- range(as.numeric(smoothDf[[TempCol]]))
      
      if(TempRange[1] > 0) {
        # extract first temp treshold (the first frame with min temp rounded to the 10)
        a <- round(TempRange[1], digits = -1) 
        a2 <-
          smoothDf[[TimeCol]][which(as.numeric(smoothDf[[TempCol]]) == a)[1]]
        
        # extract second temp treshold (the first frame with min temp rounded to the 10 + half the the first frame with min temp rounded to the 10)
        b <- round(TempRange[1], digits = -1) /
          2 + round(TempRange[1], digits = -1)
        b2 <-
          smoothDf[[TimeCol]][which(
            as.numeric(smoothDf[[TempCol]]) == b)[1]]
        
        
        # extract third temp treshold (the first frame with min temp rounded to the 10 + the first frame with the min temp rounded to the 10)
        c <- round(TempRange[1], digits = -1) + round(TempRange[1], digits = -1)
        c2 <-
          smoothDf[[TimeCol]][which(
            as.numeric(smoothDf[[TempCol]]) == c)[1]]
        
        # extract fourth temp treshold (the first frame with max temp)
        d <- TempRange[2]
        d2 <-
          smoothDf[[TimeCol]][which(as.numeric(smoothDf[[TempCol]]) == d)[1]]
        
        
        # extract fifth temp treshold (the last frame with min temp rounded to the 10 + the first frame with the min temp rounded to the 10)
        e2 <-
          smoothDf[[TimeCol]][which(
            as.numeric(smoothDf[[TempCol]]) == c)[length(which(
              as.numeric(smoothDf[[TempCol]]) == round(TempRange[1], digits = -1) + round(TempRange[1], digits = -1)
            ))]]
        
        # extract sixth temp treshold (the last frame with min temp rounded to the 10 + half the the first frame with min temp rounded to the 10)
        f2 <-
          smoothDf[[TimeCol]][which(
            as.numeric(smoothDf[[TempCol]]) == b)[length(which(
              as.numeric(smoothDf[[TempCol]]) == round(TempRange[1], digits = -1) /
                2 + round(TempRange[1], digits = -1)
            ))]]
        
        # extract seventh temp treshold (the last frame with min temp rounded to the 10)
        g2 <-
          smoothDf[[TimeCol]][which(as.numeric(smoothDf[[TempCol]]) == a)[length(which(
            as.numeric(smoothDf[[TempCol]]) == round(TempRange[1], digits = -1)
          ))]]
      } else { 
        # extract first temp threshold 
        a <-  round(TempRange[2], digits = -1) - 4
        a2 <-
          smoothDf[[TimeCol]][which(as.numeric(smoothDf[[TempCol]]) == a)[1]]
        
        # extract second temp threshold (the first frame with min temp / 2)
        b <-  a / 2
        b2 <-
          smoothDf[[TimeCol]][which(
            as.numeric(smoothDf[[TempCol]]) == b)[1]]
        
        # extract third temp threshold (the first frame with second threshold subtracted by second threshold)
        c <- b - b
        c2 <-
          smoothDf[[TimeCol]][which(
            as.numeric(smoothDf[[TempCol]]) == c)[1]]
        
        # extract fourth temp threshold (the first frame with the third threshold + the second threshold)
        d <- TempRange[1]
        d2 <-
          smoothDf[[TimeCol]][which(as.numeric(smoothDf[[TempCol]]) == d)[1]]
        
        # extract fifth temp treshold (the last frame with min temp rounded to the 10 + the first frame with the min temp rounded to the 10)
        e2 <-
          smoothDf[[TimeCol]][which(
            as.numeric(smoothDf[[TempCol]]) == c)[length(which(
              as.numeric(smoothDf[[TempCol]]) == c
            ))]]
        
        # extract sixth temp treshold (the last frame with min temp rounded to the 10 + half the the first frame with min temp rounded to the 10)
        f2 <-
          smoothDf[[TimeCol]][which(
            as.numeric(smoothDf[[TempCol]]) == b)[length(which(
              as.numeric(smoothDf[[TempCol]]) == b
            ))]]
        
        # extract seventh temp treshold (the last frame with min temp rounded to the 10)
        g2 <-
          smoothDf[[TimeCol]][which(as.numeric(smoothDf[[TempCol]]) == a)[length(which(
            as.numeric(smoothDf[[TempCol]]) == a
          ))]]
      }
    } else if(TempAxis == "y") {
      # build temperature x axis
      TempRange <- range(as.numeric(smoothDf[[TempCol]]), na.rm= T)
      
      if(TempRange[1] > 0) {
        # extract first temp treshold (the first frame with min temp rounded to the 10)
        a <- round(TempRange[1], digits = -1) 
        a2 <-
          smoothDf[[TimeCol]][which(as.numeric(smoothDf[[TempCol]]) == a)[1]]
        
        # extract second temp treshold (the first frame with min temp rounded to the 10 + half the the first frame with min temp rounded to the 10)
        b <- round(TempRange[1], digits = -1) /
          2 + round(TempRange[1], digits = -1)
        b2 <-
          smoothDf[[TimeCol]][which(
            as.numeric(smoothDf[[TempCol]]) == b)[1]]
        
        
        # extract third temp treshold (the first frame with min temp rounded to the 10 + the first frame with the min temp rounded to the 10)
        c <- round(TempRange[1], digits = -1) + round(TempRange[1], digits = -1)
        c2 <-
          smoothDf[[TimeCol]][which(
            as.numeric(smoothDf[[TempCol]]) == c)[1]]
        
        # extract fourth temp treshold (the first frame with max temp)
        d <- TempRange[2]
        d2 <-
          smoothDf[[TimeCol]][which(as.numeric(smoothDf[[TempCol]]) == d)[1]]
        
        
        # extract fifth temp treshold (the last frame with min temp rounded to the 10 + the first frame with the min temp rounded to the 10)
        e2 <-
          smoothDf[[TimeCol]][which(
            as.numeric(smoothDf[[TempCol]]) == c)[length(which(
              as.numeric(smoothDf[[TempCol]]) == round(TempRange[1], digits = -1) + round(TempRange[1], digits = -1)
            ))]]
        
        # extract sixth temp treshold (the last frame with min temp rounded to the 10 + half the the first frame with min temp rounded to the 10)
        f2 <-
          smoothDf[[TimeCol]][which(
            as.numeric(smoothDf[[TempCol]]) == b)[length(which(
              as.numeric(smoothDf[[TempCol]]) == round(TempRange[1], digits = -1) /
                2 + round(TempRange[1], digits = -1)
            ))]]
        
        # extract seventh temp treshold (the last frame with min temp rounded to the 10)
        g2 <-
          smoothDf[[TimeCol]][which(as.numeric(smoothDf[[TempCol]]) == a)[length(which(
            as.numeric(smoothDf[[TempCol]]) == round(TempRange[1], digits = -1)
          ))]]
        
      } else { 
        # extract first temp threshold 
        a <-  round(TempRange[2], digits = -1) - 4
        a2 <-
          smoothDf[[TimeCol]][which(as.numeric(smoothDf[[TempCol]]) == a)[1]]
        
        # extract second temp threshold (the first frame with min temp / 2)
        b <-  a / 2
        b2 <-
          smoothDf[[TimeCol]][which(
            as.numeric(smoothDf[[TempCol]]) == b)[1]]
        
        # extract third temp threshold (the first frame with second threshold subtracted by second threshold)
        c <- b - b
        c2 <-
          smoothDf[[TimeCol]][which(
            as.numeric(smoothDf[[TempCol]]) == c)[1]]
        
        # extract fourth temp threshold (the first frame with the third threshold + the second threshold)
        d <- TempRange[1]
        d2 <-
          smoothDf[[TimeCol]][which(as.numeric(smoothDf[[TempCol]]) == d)[1]]
        
        # extract fifth temp treshold (the last frame with min temp rounded to the 10 + the first frame with the min temp rounded to the 10)
        e2 <-
          smoothDf[[TimeCol]][which(
            as.numeric(smoothDf[[TempCol]]) == c)[length(which(
              as.numeric(smoothDf[[TempCol]]) == c
            ))]]
        
        # extract sixth temp treshold (the last frame with min temp rounded to the 10 + half the the first frame with min temp rounded to the 10)
        f2 <-
          smoothDf[[TimeCol]][which(
            as.numeric(smoothDf[[TempCol]]) == b)[length(which(
              as.numeric(smoothDf[[TempCol]]) == b
            ))]]
        
        # extract seventh temp treshold (the last frame with min temp rounded to the 10)
        g2 <-
          smoothDf[[TimeCol]][which(as.numeric(smoothDf[[TempCol]]) == a)[length(which(
            as.numeric(smoothDf[[TempCol]]) == a
          ))]]
      }
    }
    
    # concatenate then within vector and then df
    tempAxisTemp<- c(a,b,c,d,c,b,a)
    if(TempAxis == "x"){
      tempAxislab <- c(a2,b2,c2,d2,e2,f2,g2)
      tempdf <- as.data.frame(cbind(tempAxisTemp, tempAxislab))
    } else if(TempAxis == "y"){
      tempAxislab <- c(a2,b2,c2,d2,e2,f2,g2)
      tempdf <-
        as.data.frame(cbind(tempAxisTemp[!duplicated(tempAxisTemp)], tempAxislab[!duplicated(tempAxisTemp)]))
    }
    
    ## add temperature curve
    # compute the mean temperature per frame 
    meanTempf <- tapply(as.numeric(smoothDf[[TempCol]]), smoothDf[[TimeCol]], mean)
    meanTempfRange <- range(meanTempf, na.rm = TRUE)
    if(meanTempfRange[1] < 0){
      minT <- round(meanTempfRange[1], -1)
      maxT <- round(meanTempfRange[2], -1)
    } else {
      minT <- meanTempfRange[1]
      maxT <- round(meanTempfRange[2], 1)
    }
    Tdf <-as.data.frame(cbind(Temp = meanTempf, frame = unique(smoothDf[[TimeCol]])))
    if(min(Tdf[["Temp"]], na.rm = T) > 0) { 
      Pal <- colorRampPalette(c("yellow", "red"))
    }else{ 
      Pal <- colorRampPalette(c("blue", "yellow"))
    }
    coloration <- Pal(meanTempfRange[2] - meanTempfRange[1] + 1)
    dfCol <- as.data.frame(cbind(Temp = seq(meanTempfRange[1], meanTempfRange[2], by = 1), coloration))
    Tdf$color <-
      coloration[match(round(Tdf[["Temp"]]), dfCol[[TempCol]])]
    
    # X axis ranges
    minTime <- ifelse(is.null(TimeLimit), 0, TimeLimit[1]/25/60)
    maxTime <- ifelse(
      is.null(TimeLimit),
      round(max(smoothDf[[TimeCol]]/25/60, na.rm = T), 0),
      TimeLimit[2]/25/60
    )
    xscale <- pretty(c(minTime, maxTime), n = 5, bounds = TRUE)
    xscale2 <- xscale
    xscale2[which(xscale2==0)] <- ""
    
    opar <- par()
    graphT <- function(){
      # specify plot window parameters
      opar <- par()
      par(mar = c(5, 4, 4, 4) + 0.3)
      corners = par("usr") # Gets the four corners of plot area (x1, x2, y1, y2)
      par(xpd = TRUE) #Draw outside plot area
      
     if(!is.null(TempScale)){
        yscale <- pretty(c(TempScale[1], TempScale[2]), n = 5, bounds = TRUE)
      }else{
      if(minT < 0 ){ 
        #if(minT < -8){ 
        #  yscale <- c(seq(floor(minT), maxT+5, by = 5))
        yscale <- c(-8, seq(-5, 25, by = 5))
        #}else{
        #  yscale <- c(-8, seq(-5, maxT+5, by = 5))
       }
      #}else{ 
        if(minT > 0 | minT > 18){ 
          yscale <- c(minT, seq(minT, 50, by = 5))
       # }else{
        #  yscale <- c(18, seq(round(floor(minT),-1), maxT+5, by = 5))
        }
      }
        #}
      # Create an empty plot 
      plot(
        NULL,
        yaxt = "n",
        ylab = "",
        xlab = "",
        axes = FALSE,
        xlim = c(
          xscale[1]*25*60,
          xscale[length(xscale)]*25*60
        ),
        ylim = c(yscale[1], maxT + 5),
      )
      
       # add temperature curve
      if(is.null(colTemp)){ 
        with(Tdf,
             segments(
               head(frame, -1),
               head(Temp, -1),
               frame[-1],
               Temp[-1],
               color
             ))
      }else{
        lines(Tdf[["Temp"]] ~ Tdf[["frame"]] , col = colTemp)
      }
      # add plot axis title and display x axis in minute instead of frames
      par(new=TRUE)
      plot(
        NULL,
        yaxt = "n",
        ylab = "",
        xlab = "",
        axes = FALSE,
        xlim = c(
          xscale[1],
          xscale[length(xscale)]
        ),
        ylim = c(yscale[1], maxT+5),
      )
      
      title(ylab="Temperature (\u00B0C)", line=1.5, cex.lab=1.2, family="sans")
      title(xlab = "Time (min)", line=2.5, cex.lab=1.2, family="sans")
      
      # draw y axis
      segments(
        x0 = xscale[1],
        y0 = yscale[1],
        x1 = xscale[1],
        y1 = maxT+5
      )
      
      text(
        rep(xscale[1], length(yscale)),
        yscale+0.15,
        yscale,
        xpd = TRUE,
        srt = 0,
        adj = 1.5,
        pos = 2
      )
      text(
        rep(xscale[1], length(yscale)),
        yscale+0.15,
        "-",
        xpd = TRUE,
        srt = 0,
        adj = 0.9
      )
      
      # draw x axis
      segments(
        x0 = xscale[1],
        y0 = yscale[1],
        x1 = xscale[length(xscale)],
        y1 = yscale[1]
      )
      
      text(
        xscale+2,
        rep(minT-0.5, length(xscale)),
        xscale,
        xpd = TRUE,
        srt = 45,
        adj = 2,
        pos = 2
      )
      
      text(
        xscale,
        rep(yscale[1], length(xscale)),
        "-",
        xpd = TRUE,
        srt = 90,
        adj = 1.2)
      
      # add right y axis
      par(new=TRUE)
      plot(
        NULL,
        yaxt = "n",
        ylab = "",
        xlab = "",
        axes = FALSE,
        xlim = c(
          xscale[1],
          xscale[length(xscale)]
        ),
        ylim = c(0, 1),
      )
      segments(
        x0 = xscale[length(xscale)],
        y0 = 0,
        x1 = xscale[length(xscale)],
        y1 = 1,
        col = "black"
      )
      
      # add variable  
      linepar <- c(rep(1.8, 3), rep(2.8, 3))
      #IC <- base::split(IC, IC[["customfunc"]])
      for (v in seq(length(variables))){ 
        smoothDf[[variables[[v]]]][which(is.infinite(smoothDf[[variables[[v]]]]))] <- NA
        smoothDfNoNA <- list()
        if(length(which(is.na(smoothDf[[variables[[v]]]])== TRUE)) > 0){
          smoothDfNoNA[[variables[[v]]]] <- smoothDf[[variables[[v]]]][-c(which(is.na(smoothDf[[variables[[v]]]]) == TRUE))]
          smoothDfNoNA[[TimeCol]] <- smoothDf[[TimeCol]][-c(which(is.na(smoothDf[[variables[[v]]]]) == TRUE))]
          } else{ 
          smoothDfNoNA = smoothDf
          }
        if(is.null(scaleVariable)){
        if (!is.null(IC)) {
          Scale <-
            pretty(
              c(IC[[variables[v]]]$"2.5%", IC[[variables[[v]]]]$"97.5%"),
              min.n = 0,
              n = 2,
              bounds = TRUE
            )
          
        }else{
          Scale <- pretty(smoothDfNoNA[[variables[[v]]]], n = 5, bounds = TRUE)
        }}else{
          if(length(scaleVariable[[v]]) > 2) {
            Scale <- seq(from = scaleVariable[[v]][1], to = scaleVariable[[v]][2], by = scaleVariable[[v]][3])
          }else{
          Scale <- pretty(scaleVariable[[v]], n = 5, bounds = TRUE)
          }
        }
        # add right y axis name labels
        step <- Scale[2] - Scale[1]
        if(Scale[1] < 0) {
          Scale <- Scale[-1]}
        if(Scale[1] > 0){
          Scale[length(Scale)+1] <- 0
          Scale <- sort(Scale)
        }
        if(length(Scale) < 6){
          for(i in seq((6 - length(Scale)))){
          Scale[length(Scale) + 1] <- Scale[length(Scale)] + step
          }
        }else if(length(Scale) > 6){
          Scale <- Scale[-length(Scale)]
        }
        
        par(new=TRUE)
        plot(
          NULL,
          yaxt = "n",
          ylab = "",
          xlab = "",
          axes = FALSE, 
          xlim = c(
            xscale[1]*25*60,
            xscale[length(xscale)]*25*60
          ),
          ylim = c(0, max(Scale)),
        )
        
        # add variable smoothed line(s)
        lines(smoothDfNoNA[[variables[[v]]]] ~ smoothDfNoNA[[TimeCol]], col =  colvariable[[v]])
        if(!is.null(IC)){
          IC[[variables[[v]]]]$`2.5%.smoothed` <- MoveR::slidWindow(IC[[variables[[v]]]]$"2.5%", Tstep = 4, customFunc = function(x) mean(x, na.rm = T))
          IC[[variables[[v]]]]$`97.5%.smoothed` <- MoveR::slidWindow(IC[[variables[[v]]]]$"97.5%", Tstep = 4, customFunc = function(x) mean(x, na.rm = T))
          IC$temporaryNoNA <- na.omit(IC[[variables[[v]]]])
          polygon(x = c(IC$temporaryNoNA[[TimeCol]], rev(IC$temporaryNoNA[[TimeCol]])),
                  y = c(IC$temporaryNoNA$`2.5%.smoothed`, rev(IC$temporaryNoNA$`97.5%.smoothed`)), 
                  col=adjustcolor(colvariable[[v]],alpha=0.15), border=NA
                  , density = NA)}
        axislab <- c(0, step/4, -step/4, step/2, -step/2, step/1.5)
        if(v == 1){
          text(rep(ifelse(
            is.null(TimeLimit),
            max(smoothDfNoNA[[TimeCol]], na.rm = T),
            TimeLimit[2]
          ), 6),
          Scale + step/100,
          "-",
          xpd = TRUE,
          adj = 0)
          
          text(rep(ifelse(
            is.null(TimeLimit),
            max(smoothDfNoNA[[TimeCol]], na.rm = T),
            TimeLimit[2]
          ) + ifelse(
            is.null(TimeLimit),
            max(smoothDfNoNA[[TimeCol]], na.rm = T),
            TimeLimit[2]
          ) / 1520, 6),
          Scale + axislab[[v]],
          Scale,
          xpd = TRUE,
          col = colvariable[[v]], 
          pos = 4,
          cex = 1 - length(variables)/20)
        }
        if(v > 1){
        text(rep(ifelse(
          is.null(TimeLimit),
          max(smoothDfNoNA[[TimeCol]], na.rm = T),
          TimeLimit[2]
        ) + ifelse(
          is.null(TimeLimit),
          max(smoothDfNoNA[[TimeCol]], na.rm = T),
          TimeLimit[2]
        ) / 1520, 5),
        Scale[-1] + axislab[[v]],
        Scale[-1],
        xpd = TRUE,
        col = colvariable[[v]], 
        pos = 4,
        cex = 1 - length(variables)/20)
        }
        usrpar <- rep(c(Scale[3] + step/2,  Scale[5] + step/2, Scale[1] + step/3), 2)
        if(is.null(Lab)){
          mtext(variables[[v]], side = 4, line = linepar[[v]], at = usrpar[[v]], col = colvariable[[v]], cex=1.2, family="sans")
        } else {
          mtext(Lab[[v]], side = 4, line = linepar[[v]],  at = usrpar[[v]], col = colvariable[[v]], cex=1.2, family="sans")
        }
      }
      }
      
    ## add relative detected indiv number in the upper part of the graph
    # compute the number of ID detected per frame and divide it by the true number of Id (manually counted)
    nbrInd <- tapply(finalDatList[["identity"]], finalDatList[[TimeCol]], length)
    nbrInd[is.na(nbrInd)] <- 0
    nbrInd[nbrInd == 1] <- 0
    nbrIndRelative <- nbrInd /nbrIndcounted
    # Compute the mean relative number of ID per frame using sliding window
    nbrIndRelativeMean <- MoveR::slidWindow(nbrIndRelative,
                                  Tstep = 100 * 25, function (x)
                                    mean(x, na.rm = T))
    
    # add the smoothed relative number of ID per frame to the previous plot 
    # compute range of the axis and round them
    nbrIndRange <- range(nbrIndRelativeMean, na.rm = TRUE)
    minRID <- 0
    if(minT < 0) { 
      maxRID <- 1
    } else { 
      maxRID <- round(nbrIndRange[2], 0)
    }

    graphRID <- function() {
      plot(
        NULL,
        yaxt = "n",
        ylab = "",
        xlab = "",
        axes = FALSE,
        xlim = c(
          xscale[1],
          xscale[length(xscale)]
        ),
        ylim = c(minRID, maxRID),
      )
      # draw y axis
      segments(
        x0 = xscale[1],
        y0 = minRID,
        x1 = xscale[1],
        y1 = maxRID
      )
      yscaleRID <- seq(minRID, maxRID, by = 0.2)
      text(
        rep(xscale[1], length(yscaleRID)),
        yscaleRID+0.02,
        yscaleRID,
        xpd = TRUE,
        srt = 0,
        adj = 1.6, 
        cex = 0.7, 
        pos = 2
      )
      text(
        rep(xscale[1], length(yscaleRID)),
        yscaleRID+0.03,
        "-",
        xpd = TRUE,
        srt = 0,
        adj = 0.9
      )
      # draw x axis
      segments(
        x0 = xscale[1],
        y0 = 0,
        x1 = xscale[length(xscale)],
        y1 = 0
      )
      text(mean(xscale), 1.15, "Detected individuals", cex = 0.8)
      par(new=TRUE)
      plot(
        NULL,
        yaxt = "n",
        ylab = "",
        xlab = "",
        axes = FALSE,
        xlim = c(xscale[1]*60*25, xscale[length(xscale)]*60*25),
        ylim = c(minRID, maxRID),
      )
      lines(nbrIndRelativeMean ~ unique(finalDatList[[TimeCol]])[order(unique(finalDatList[[TimeCol]]))], col = "#FF99CC")
    }
    
    par(fig=c(0,1,0,0.88))
    graphT()
    par(fig=c(0,0.99,0.55,1), new=TRUE)
    graphRID()
    par <- opar
    
  }

