library(shiny)
library(textClust)
library(stream)
library(DT)
library(mongolite)
library(Rtsne)
library(plotly)
library(shinyWidgets)

options(shiny.maxRequestSize=10000*1024^2) ## increase file size limit to 10GB


## Download buttons require orca command line tool: https://github.com/plotly/orca
Sys.setenv('MAPBOX_TOKEN' = 'secret token') ## bug: orca requires a mapbox token even if unused

colormap = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628", "#f781bf", "#999999")


server <- shinyServer(function(input, output, session) {


  # Online Version -------------------------------------------------------------------

  ## query result from mongo
  mongo_get = function(dbconnection, time, isTime =T){
    if(isTime == T){
      request = paste('{"time" :"',time,'"}', sep="")
    } else{
      request = paste('{"observation" :',time,'}', sep="")
    }
    iter <- dbconnection$iterate(request)

    data = iter$one()

    data = transformMongo(data)

    return(data)
  }

  ## transform mongo result to correct format
  transformMongo = function(data){
    # some postprocessing
    data$microDistance = matrix(unlist(data$microDistance), ncol=length(data$microRepresentatives))
    data$macroDistance = matrix(unlist(data$macroDistance), ncol=length(data$macroRepresentatives))

    # right format for microcluster
    data$microRepresentatives = lapply(1:length(data$microRepresentativesWeights),function(x){
      names(data$microRepresentativesWeights[[x]])=data$microRepresentatives[[x]]
      unlist(data$microRepresentativesWeights[[x]])
    })

    #right format for macrocluster
    data$macroRepresentatives = lapply(1:length(data$macroRepresentativesWeights),function(x){
      names(data$macroRepresentativesWeights[[x]])=data$macroRepresentatives[[x]]
      unlist(data$macroRepresentativesWeights[[x]])
    })

    data$microWeights = unlist(data$microWeights)
    data$macroWeights = unlist(data$macroWeights)

    data$dendrogram$merge = matrix(unlist(data$dendrogram$merge), ncol=2, byrow=T)
    data$dendrogram$height = unlist(data$dendrogram$height)
    data$dendrogram$order = unlist(data$dendrogram$order)

    #clean up
    data$microRepresentativesWeights= NULL
    data$macroRepresentativesWeights= NULL

    data$clusterAssignment = sapply(data$clusterAssignment, function(x){
      ifelse(x=="NA", NA_integer_, as.integer(x))
    })
    data$microToMacroAssignment = sapply(data$microToMacroAssignment, function(x){
      ifelse(x=="NA", NA_integer_, as.integer(x))
    })


    return(data)
  }


  ## connection changes when user changes stream
  collection <- reactive({
    req(input$streamSelect != "")
    tryCatch({
      mongo(url= paste("mongodb://", isolate(input$user), ":", isolate(input$password), "@", isolate(input$url), ":", isolate(input$port), "/admin", sep=""), db="streams", collection = input$streamSelect)
    }, error=function(x){
      sendSweetAlert(session, title="Database Error", text="I could not connect to the MongoDB database.", type="error")
      req(F)
    })
  })




  # Local Version ---------------------------------------------------------

  ## initialize local algorithm
  algorithm <<- DSC_textClust()
  clusteringProgress = NULL
  rows <<- 0

  ## create the source stream
  stream <- reactive({
    input$stream_stop ## to force reaction when stream is resetted
    if(input$streamType=="File"){
      req(!is.null(react$file))

      ## get number of lines without loading into memory
      testcon <- file(react$file,open="r")
      readsizeof <- 20000
      nooflines <- 0
      ( while((linesread <- length(readLines(testcon,readsizeof))) > 0 )
        nooflines <- nooflines+linesread )
      close(testcon)
      rows <<- nooflines
      sep = gsub("\\\\t", "\t", input$file_separator) ## this is not ideal but allows to support \t input
      numCols = count.fields(textConnection(readLines(react$file, n=1)), sep=sep) ## count number of columns and force to character
      DSD_ReadCSV(react$file, quote = "", comment.char = "", sep = sep, encoding="UTF-8", colClasses = rep("character",numCols))
    } else if(input$streamType=="Text"){
      sep = gsub("\\\\t", "\t", input$file_separator) ## this is not ideal but allows to support \t input
      data =  lapply(strsplit(input$streamInput, "\n")[[1]], function(x){
        strsplit(x, sep)[[1]]
      })
      data = data.frame(do.call(rbind,data), stringsAsFactors = F)
      rows <<- nrow(data)
      DSD_Memory(data)
    }
  })


  # run the algorithm locally
  observe({
    if(input$streamType != "Server"){
      req(react$run) ## require that processing should start
      dsd = stream() ## get stream
      update(algorithm, dsd, n = min(c(rows - algorithm$RObj$C$n, 1)))

      ## cluster in increments of 1 to keep UI responsive, however get result only after time interval
      if(algorithm$RObj$C$n %% isolate(input$stream_horizon) == 0) getClusteringResult()
      invalidateLater(0, session) ## invalidate immediately to trigger next run

      ## stop when we reached the end of the stream
      if(rows == algorithm$RObj$C$n){
        getClusteringResult()
        react$run = F
        if(!is.null(clusteringProgress)) clusteringProgress$close()
        clusteringProgress <<- NULL
      }
    }
  })


  # store the local result set
  getClusteringResult = function(){
    observation = as.character(algorithm$RObj$C$n)
    isolate({
      result[[observation]]$observation <<- algorithm$RObj$C$n
      result[[observation]]$time <<- algorithm$RObj$currenttime
      result[[observation]]$microDistance <<- as.matrix(algorithm$RObj$dist(get_centers(algorithm, type="micro")))
      result[[observation]]$macroDistance <<- as.matrix(algorithm$RObj$dist(get_centers(algorithm, type="macro")))
      result[[observation]]$microRepresentatives <<- algorithm$RObj$getRepresentatives(Inf, "micro")
      result[[observation]]$macroRepresentatives <<- algorithm$RObj$getRepresentatives(Inf, "macro")
      result[[observation]]$microWeights <<- get_weights(algorithm, type="micro")
      result[[observation]]$macroWeights <<- get_weights(algorithm, type="macro")
      result[[observation]]$clusterAssignment <<- algorithm$RObj$get_assignment()
      result[[observation]]$microToMacroAssignment <<- algorithm$RObj$microToMacro()
      result[[observation]]$dendrogram <<- algorithm$RObj$hc[1:3] ## store stripped hclust object (merge, height and order is sufficient information)
    })
    if(input$time_obs==T){
      updateObservationSlider(sapply(result, function(x){
        x$time
      }))
    } else{
      updateObservationSlider(names(result))
      ## or
      # updateObservationSlider(sapply(result, function(x){
      #   x$observation
      # }))
    }
  }


  output$export <- downloadHandler(
    filename = function(){
      paste("Export.RData", sep="")
    },
    content = function(file){
      inputfile = isolate(react$file)
      inputTextColumn=input$file_textColumn
      inputTimeColumn=input$file_timeColumn
      inputSeparator=input$file_separator
      inputGroupByColumn=input$file_groupByColumn
      inputParentTextColumn=input$file_parentTextColumn
      inputParentTimeColumn=input$file_parentTimeColumn

      if(input$streamType == "Server"){
        it <- collection()$iterate(query='{}')

        res = list()
        while(!is.null(x <- it$one())){
          res[[as.character(x$observation)]] = x
        }
        res = lapply(res, transformMongo)
        result = res

        save(result, file=file)
      } else if(input$streamType == "File"){
        save(result,inputfile,inputTextColumn,inputTimeColumn,inputSeparator,inputGroupByColumn,inputParentTextColumn,inputParentTimeColumn,file=file)
      }


    }
  )

  observe({
    req(input$import$datapath)
    load(input$import$datapath, envir = .GlobalEnv)
    react$file <- inputfile
    updateNumericInput(session, "file_textColumn", value = inputTextColumn)
    updateNumericInput(session, "file_timeColumn", value = inputTimeColumn)
    updateNumericInput(session, "file_groupByColumn", value = inputGroupByColumn)
    updateNumericInput(session, "file_parentTextColumn", value = inputParentTextColumn)
    updateNumericInput(session, "file_parentTimeColumn", value = inputParentTimeColumn)

    updateTextInput(session, "file_separator", value = inputSeparator)
    if(input$time_obs==T){
      updateObservationSlider(sapply(result, function(x){
        x$time
      }))
    } else{
      updateObservationSlider(names(result))
    }
  })

  observe({
    react$file <<- input$file$datapath
  })
  ## reinitialize the local algorithm
  reset <- function(){
    result <<- list()
    react$run <<- NULL
    react$microRepresentatives <<- NULL
    if(!is.null(input$stream_setting_stopwords_sources)){
        langs = strsplit(input$stream_setting_stopwords_sources, split = " - ")
        words = lapply(langs, function(x){
          stopwords(language=x[2], source=x[1])
        })
        stpwords = unique(unlist(words))
    } else{
      stpwords = character()
    }

    k <- ifelse(is.na(input$stream_setting_k), 5L, input$stream_setting_k)
    r <- ifelse(is.na(input$stream_setting_r), 0.6 , input$stream_setting_r)
    tgap <- ifelse(is.na(input$stream_setting_tgap), 100, input$stream_setting_tgap)
    lambda <- ifelse(is.na(input$stream_setting_lambda), 5L, input$stream_setting_lambda)
    minWeight <- ifelse(is.na(input$stream_setting_minWeight), 3, input$stream_setting_minWeight)
    algorithm <<- DSC_textClust(r = r, lambda = lambda, tgap = tgap, nmin = input$stream_setting_ngrams[1], nmax = input$stream_setting_ngrams[2], k = k, termFading = input$stream_setting_termFading, stopword = stpwords, linkage = input$stream_setting_linkage, weightedReclustering = input$stream_setting_weightedReclustering, textCol = input$file_textColumn, timeCol = input$file_timeColumn, minWeight = minWeight, verbose=F, groupByCol = input$file_groupByColumn, parentTextCol = input$file_parentTextColumn, parentTimeCol = input$file_timeColumn, timeFormat=input$file_timeFormat, timePrecision=input$file_timePrecision)
  }






  # Update logic ---------------------------------------------------------
  react = reactiveValues(file = NULL, run=FALSE, microDistance=NULL, macroDistance=NULL, microRepresentatives=NULL, macroRepresentatives=NULL, microWeights=NULL, macroWeights=NULL, clusterAssignment=NULL, microToMacroAssignment=NULL, dendrogram=NULL)

  updateObservationSlider <- function(choices, selected=NULL){
    choices=as.vector(choices)
    if(input$live){
      updateSliderTextInput(session, "observation", choices=c(0,choices), selected=tail(choices,1))
    } else{
      if(is.null(selected)){
        updateSliderTextInput(session, "observation", choices=c(0,choices), selected=selected)
      } else{
        updateSliderTextInput(session, "observation", choices=c(0,choices))
      }
    }
  }


  ## query the result whenever time selection changes or to poll new data
  observe({
    if(input$streamType != "Server"){
      ## either locally
      observation = as.character(input$observation)
      if(isolate(input$time_obs)==T && observation != "0"){
        time = sapply(result, function(x){
          x$time
        })
        select = names(time)[time == observation] ## TODO this can get slow
        observation = as.character(result[[select]]$observation)
      }
      react$microDistance <- result[[observation]]$microDistance
      react$microRepresentatives = lapply(result[[observation]]$microRepresentatives, function(x){
        head(x, input$tokens)
      })
      react$microWeights <- result[[observation]]$microWeights

      react$macroDistance <- result[[observation]]$macroDistance
      react$macroRepresentatives = lapply(result[[observation]]$macroRepresentatives, function(x){
        head(x, input$tokens)
      })
      react$macroWeights <- result[[observation]]$macroWeights

      react$clusterAssignment = result[[observation]]$clusterAssignment
      react$microToMacroAssignment = result[[observation]]$microToMacroAssignment

      react$dendrogram = result[[observation]]$dendrogram
    } else{
      withProgress(message = 'Polling...', value = 1, {
        ## or from database
        if(input$time_obs == T){
          choices = collection()$find(query='{}', fields='{"time": 1}')$time
        }
        else{
          choices = collection()$find(query='{}', fields='{"observation": 1}')$observation
        }
        updateObservationSlider(choices)

        tryCatch({
          data = mongo_get(collection(), input$observation, input$time_obs)
          react$microRepresentatives = lapply(data$microRepresentatives, function(x){
            head(x, input$tokens)
          })
          react$microDistance = data$microDistance
          react$microWeights = data$microWeights

          react$macroRepresentatives = lapply(data$macroRepresentatives, function(x){
            head(x, input$tokens)
          })
          react$macroDistance = data$macroDistance
          react$macroWeights = data$macroWeights


          react$clusterAssignment = data$clusterAssignment
          react$microToMacroAssignment = data$microToMacroAssignment

          react$dendrogram = data$dendrogram

        }, error=function(x){
          react$microRepresentatives = NULL
        })
        if(input$refreshRate!="None") invalidateLater(1000*as.integer(input$refreshRate), session) ## polling rate when from server
      })
    }
  })




  # Buttons ----------------------------------------------------------------
  stopAndReset <- function(){
    react$run <- FALSE
    if(!is.null(clusteringProgress)) clusteringProgress$close()
    clusteringProgress <<- NULL
    reset()
    updateObservationSlider(0, 0)
  }

  ## start button
  observeEvent(input$stream_start, {
    react$run <- TRUE
    if(is.null(clusteringProgress)) clusteringProgress <<- shiny::Progress$new()
    clusteringProgress$set(message = "In progress...", value = 1)
  })

  ## pause button
  observeEvent(input$stream_pause, {
    react$run <- FALSE
    if(!is.null(clusteringProgress)) clusteringProgress$set(message = "Paused...", value = 1)
  })

  ## stop button
  observeEvent(input$stream_stop, {
    stopAndReset()
  })

  ## stream type selection
  observeEvent(input$streamType,{
    stopAndReset()
  })

  ## connect button
  observeEvent(input$connect,{
    tryCatch({
      m <<- mongo(url= paste("mongodb://", input$user, ":", input$password, "@", input$url, ":", input$port, "/admin", sep=""), db="streams", collection = "mapping")
      streams = m$find(query='{}', fields='{"_id":0, "name":1}')
      choices = sort(streams$name)
      updateSelectInput(session, "streamSelect", choices = sort(choices), selected = choices[1])
    }, error=function(x){
      updateSelectInput(session, "streamSelect", choices = "", selected = "")
      sendSweetAlert(session, title="Database Error", text="I could not connect to the MongoDB database.", type="error")
      req(F)
    })
  })

  ## disconnect button
  observeEvent(input$disconnect,{
    m$disconnect()
    updateSelectInput(session, "streamSelect", choices = "", selected = "")
    stopAndReset()
  })

  ## parent text column selection
  observeEvent(input$file_parentTextColumn,{
    stopAndReset()
  })

  ## parent time column selection
  observeEvent(input$file_parentTimeColumn,{
    stopAndReset()
  })

  # groupBy column selection
  observeEvent(input$file_groupByColumn,{
    stopAndReset()
  })

  ## text column selection
  observeEvent(input$file_textColumn,{
    stopAndReset()
  })

  ## time column selection
  observeEvent(input$file_timeColumn,{
    stopAndReset()
  })

  ## file separator selection
  observeEvent(input$file_separator,{
    stopAndReset()
  })

  ## real time vs t switch
  observeEvent(input$time_obs,{
    if(input$streamType != "Server"){
      if(input$time_obs==T){
        updateObservationSlider(sapply(result, function(x){
          x$time
        }))
      } else{
        updateObservationSlider(names(result))
      }
    }
  })



  # Observe settings --------------------------------------------------------

  observe({algorithm$RObj$k <- ifelse(is.na(input$stream_setting_k), 5L, input$stream_setting_k)})
  observe({algorithm$RObj$linkage <- input$stream_setting_linkage})
  observe({algorithm$RObj$C$r <- ifelse(is.na(input$stream_setting_r), 0.6 , input$stream_setting_r)})
  observe({algorithm$RObj$C$tgap <- ifelse(is.na(input$stream_setting_tgap), 100, input$stream_setting_tgap)})
  observe({algorithm$RObj$C$lambda <- ifelse(is.na(input$stream_setting_lambda), 5L, input$stream_setting_lambda)})
  observe({algorithm$RObj$nmin <- input$stream_setting_ngrams[1]})
  observe({algorithm$RObj$nmax <- input$stream_setting_ngrams[2]})
  observe({algorithm$RObj$C$termFading <- input$stream_setting_termFading})
  observe({algorithm$RObj$weightedReclustering <- input$stream_setting_weightedReclustering})
  observe({algorithm$RObj$minWeight <- ifelse(is.na(input$stream_setting_minWeight), 3, input$stream_setting_minWeight)})
  observe({
    if(!is.null(input$stream_setting_stopwords_sources)){
      langs = strsplit(input$stream_setting_stopwords_sources, split = " - ")
      words = lapply(langs, function(x){
        stopwords(language=x[2], source=x[1])
      })
      algorithm$RObj$stopword <- unique(unlist(words))
    } else{
      algorithm$RObj$stopword <- character()
    }
  })



  # N-Grams -------------------------------------------------------------------
  getNgramData = function(){

    if(input$type=="micro"){
      representatives = react$microRepresentatives
      weights = react$microWeights
      if(!input$ngrams_showUnassigned){
        select = !is.na(react$microToMacroAssignment)
        representatives = representatives[select]
        weights = weights[select]
      }
    } else{
      representatives = react$macroRepresentatives
      weights = react$macroWeights
    }

    req(length(representatives)>0)

    data = data.frame(Index=seq_along(representatives), weight=round(weights,4), Tokens=sapply(representatives, function(x){
      paste(paste(names(x), " (",round(x,2),")", sep=""), collapse = ", ")
    }))

    return(data)
  }

  output$ngrams <- DT::renderDataTable({
    withProgress(message = 'Generating Table', value = 1, {
      getNgramData()
    })
  }, options=list(pageLength=20, order = list(list(1, 'desc'))), rownames= FALSE)


  output$downloadNgramData <- downloadHandler(
    filename = function(){
      paste(input$observation, "_NGram.RData", sep="")
    },
    content = function(file){
      data = getNgramData()
      save(data, file=file)
    }
  )

  # Barplots --------------------------------------------------------------
  getBarplotData = function(){

    if(input$type=="micro"){
      representatives = react$microRepresentatives
      weights = react$microWeights

      if(!input$barplot_showUnassigned){
        select = !is.na(react$microToMacroAssignment)
        representatives = representatives[select]
        weights = weights[select]
      }
      updateSliderInput(session, "barplot_numPlots", max = length(representatives))
    } else{
      representatives = react$macroRepresentatives
      weights = react$macroWeights
      updateSliderInput(session, "barplot_numPlots", max = length(representatives))
    }

    req(length(representatives)>0)


    ## order by cluster weight
    ord = order(weights, decreasing = T)
    representatives = representatives[ord]
    weights = weights[ord]


    representatives = representatives[input$barplot_numPlots[1]:input$barplot_numPlots[2]]
    weights = weights[input$barplot_numPlots[1]:input$barplot_numPlots[2]]

    return(list(representatives=representatives, weights=weights))
  }

  drawBarplot = function(){

    list2env(getBarplotData(), .GlobalEnv)

    if(input$barplot_plotType=="Individual"){
      
      ymax = max(sapply(representatives, function(x){
        x[1]
      }))
      ## create plots
      plots=lapply(seq_along(representatives),function(i){

        name = names(representatives[[i]])

        if(input$barplot_order=="Alphabetical"){
          name = factor(name, level=name[order(as.character(name))])
        } else if(input$barplot_order=="Weight"){
          name = factor(name, level=name)
        }

        weight = representatives[[i]]

        plot_ly(x=name, y=weight, type="bar", height=400*ceiling(length(representatives)/4), color=factor(i, levels = seq_along(representatives)), colors = colormap) %>% layout(showlegend = FALSE, yaxis = list(range = c(0, ymax)))
      })
      subplot(plots, nrows = ceiling(length(representatives)/4))

    } else if(input$barplot_plotType=="Grouped"){

      ## transform list into data frame
      reps = lapply(representatives, function(x){
        data.frame(token=names(x), weight=x)
      })
      suppressWarnings(reps <- Reduce(function(...) merge(..., by="token", all=TRUE), reps)) ## and merge by common tokens
      reps[is.na(reps)]=0

      if(input$barplot_order=="Alphabetical"){
        reps$token = factor(reps$token, level=reps$token[order(as.character(reps$token))])
      } else if(input$barplot_order=="Weight"){
        reps$token = factor(reps$token, level=reps$token[order(apply(reps[,-1], 1, max), decreasing = T)])
      }
      ## plot one by one
      lvls = seq_along(representatives)
      p <- plot_ly(x=reps$token, y=reps[,2], type="bar", name = "Cluster 1", width=max(c(1500,nrow(reps)*20)), height=800, color=factor(1, levels = lvls), colors = colormap)
      for(i in (1+seq_along(representatives))[-1]){
        p <- p %>% add_trace(y = reps[,i], name=paste("Cluster", i-1), color = factor(i-1, levels = lvls))
      }
      p <- p %>% layout(barmode = 'group')
    }

  }

  output$barplot <- renderPlotly({
    withProgress(message = 'Drawing Barplots...', value = 1, {
      drawBarplot()
    })
  })

  output$downloadBarplotData <- downloadHandler(
    filename = function(){
      paste(input$observation, "_Barplot.RData", sep="")
    },
    content = function(file){
      data = getBarplotData()
      save(data, file=file)
    }
  )

  output$downloadBarplot <- downloadHandler(
    filename = function(){
      paste(input$observation, "_Barplot.pdf", sep="")
    },
    content = function(file){
      p <- drawBarplot()
      orca(p, file = "temp.pdf") ## bug: orca cannot handle paths
      file.rename("temp.pdf", file) ## therefore move to directory
    }
  )


  # MDS ---------------------------------------------------------------------

  getMDSData = function(){

    if(input$type=="micro"){
      d = react$microDistance
      representatives = react$microRepresentatives
      weights = react$microWeights
      microToMacroAssignment = react$microToMacroAssignment

      if(!input$mds_showUnassigned){
        select = !is.na(microToMacroAssignment)
        d = d[select, select]
        representatives = representatives[select]
        weights = weights[select]
        microToMacroAssignment = microToMacroAssignment[select]
      }
      col = as.factor(paste("Cluster", microToMacroAssignment))

    } else{
      d = react$macroDistance
      representatives = react$macroRepresentatives
      weights = react$macroWeights
      col = as.factor(paste("Cluster", seq_along(representatives)))
    }

    req(length(representatives)>0)


    dims = as.integer(input$mds_dimensionality)

    mds <- cmdscale(d, k=dims)
    ## get word with highest weight for each mc
    tokens = sapply(representatives, function(x){
      paste(names(x), collapse=", ")
    })

    return(list(mds=mds, tokens=tokens, weights=weights, col=col))
  }

  drawMDS = function(){

    list2env(getMDSData(), .GlobalEnv)

    dims = as.integer(input$mds_dimensionality)

    if(dims==2){
      if(input$mds_plotType=="Bubble"){
        plot_ly(x = mds[,1], y = mds[,2], text = tokens, type = 'scatter', mode = 'markers', size=weights, sizes = c(10, 500), color=col, colors=colormap)
      } else if(input$mds_plotType=="Text"){
        plot_ly(x = mds[,1], y = mds[,2], text =  tokens, type = 'scatter', mode = 'markers+text', size=weights, sizes = c(10, 500), color=col, colors=colormap, textposition = 'middle right', textfont = list(size = 16))
      }
    } else if(dims==3){
      if(input$mds_plotType=="Bubble"){
        plot_ly(x = mds[,1], y = mds[,2], z = mds[,3], text = tokens, type = 'scatter3d', mode = 'markers', size=weights, sizes = c(200, 1500), color=col, colors=colormap)
      } else if(input$mds_plotType=="Text"){
        plot_ly(x = mds[,1], y = mds[,2], z = mds[,3], text = tokens, type = 'scatter3d', mode = 'markers+text', size=weights, sizes = c(200, 1500), color=col, colors=colormap, textposition = 'middle right', textfont = list(size = 16))
      }
    }

  }

  output$MDS <- renderPlotly({
    withProgress(message = 'Computing Multidimensional Scaling', value = 1, {
      drawMDS()
    })
  })

  output$downloadMDSData <- downloadHandler(
    filename = function(){
      paste(input$observation, "_MDS.RData", sep="")
    },
    content = function(file){
      data = getMDSData()
      save(data, file=file)
    }
  )

  output$downloadMDS <- downloadHandler(
    filename = function(){
      paste(input$observation, "_MDS.pdf", sep="")
    },
    content = function(file){
      p <- drawMDS()
      orca(p, file = "temp.pdf") ## bug: orca cannot handle paths
      file.rename("temp.pdf", file) ## therefore move to directory
    }
  )


  # T-sne -------------------------------------------------------------------
  getTsneData = function(){

    if(input$type=="micro"){
      d = react$microDistance
      representatives = react$microRepresentatives
      weights = react$microWeights
      microToMacroAssignment = react$microToMacroAssignment

      if(!input$tsne_showUnassigned){
        select = !is.na(microToMacroAssignment)
        d = d[select, select]
        representatives = representatives[select]
        weights = weights[select]
        microToMacroAssignment = microToMacroAssignment[select]
      }
      col = as.factor(paste("Cluster", microToMacroAssignment))

    } else{
      d = react$macroDistance
      representatives = react$macroRepresentatives
      weights = react$macroWeights
      col = as.factor(paste("Cluster", seq_along(representatives)))
    }

    req(length(representatives)>0)


    dims = as.integer(input$tsne_dimensionality)

    updateSliderInput(session, "tsne_perplexity", max = floor((nrow(d)-1)/3))

    tsne_out <- Rtsne(d, dims=dims, is_distance = T, perplexity = input$tsne_perplexity, pca=input$tsne_pca, theta=input$tsne_theta, max_iter=input$tsne_iterations)

    ## get word with highest weight for each mc
    tokens = sapply(representatives, function(x){
      paste(names(x), collapse=", ")
    })

    return(list(tsne_out=tsne_out, tokens=tokens, weights=weights, col=col))
  }


  drawTsne = function(){

    list2env(getTsneData(), .GlobalEnv)

    dims = as.integer(input$tsne_dimensionality)

    if(dims==2){
      if(input$tsne_plotType=="Bubble"){
        plot_ly(x = tsne_out$Y[,1], y = tsne_out$Y[,2], text = tokens, type = 'scatter', mode = 'markers', size=weights, sizes = c(10, 500), color=col, colors=colormap)
      } else if(input$tsne_plotType=="Text"){
        plot_ly(x = tsne_out$Y[,1], y = tsne_out$Y[,2], text = tokens, type = 'scatter', mode = 'markers+text', size=weights, sizes = c(10, 500), color=col, colors=colormap, textposition = 'middle right', textfont = list(size = 16)) ## TODO i think the marker size is not correctly drawn
      }
    } else if(dims==3){
      if(input$tsne_plotType=="Bubble"){
        plot_ly(x = tsne_out$Y[,1], y = tsne_out$Y[,2], z = tsne_out$Y[,3], text = tokens, type = 'scatter3d', mode = 'markers', size=weights, sizes = c(200, 1500), color=col, colors=colormap)
      } else if(input$tsne_plotType=="Text"){
        plot_ly(x = tsne_out$Y[,1], y = tsne_out$Y[,2], z = tsne_out$Y[,3], text = tokens, type = 'scatter3d', mode = 'markers+text', size=weights, sizes = c(200, 1500), color=col, colors=colormap, textposition = 'middle right', textfont = list(size = 16))
      }
    }
  }


  output$tsne <- renderPlotly({
    withProgress(message = 'Computing t-SNE', value = 1, {
      drawTsne()
    })
  })

  output$downloadTsneData <- downloadHandler(
    filename = function(){
      paste(input$observation, "_tsne.RData", sep="")
    },
    content = function(file){
      data = getTsneData()
      save(data, file=file)
    }
  )



  output$downloadTsne <- downloadHandler(
    filename = function(){
      paste(input$observation, "_tsne.pdf", sep="")
    },
    content = function(file){
      p <- drawTsne()
      orca(p, file = "temp.pdf") ## bug: orca cannot handle paths
      file.rename("temp.pdf", file) ## therefore move to directory
    }
  )




  # Assignment -------------------------------------------------------------------

  getData <- function(){
    ## get data
    sep = gsub("\\\\t", "\t", input$file_separator) ## this is not ideal but allows to support \t input
    if(input$streamType=="File"){
      data = read.table(react$file, sep=sep, quote = "", comment.char = "", stringsAsFactors = F, nrows=input$observation)
    } else if(input$streamType=="Text"){
      data =  lapply(strsplit(input$streamInput, "\n")[[1]], function(x){
        strsplit(x, sep)[[1]]
      })
      data = data.frame(do.call(rbind,data), stringsAsFactors = F)[seq_len(input$observation),, drop=FALSE]
    }
    else if(input$streamType=="Server"){
      name = strsplit(input$streamSelect,"--")[[1]][1]
      data = read.table(paste("/data/",name,sep=""), sep=sep, quote = "", comment.char = "", stringsAsFactors = F, nrows=input$observation)
    }
  }

  getAssignmentData = function(){
    req(length(react$microRepresentatives)>0)

    if(input$type=="micro"){
      representatives = react$microRepresentatives
      assignment = react$clusterAssignment
    } else{
      representatives = react$macroRepresentatives
      assignment = react$microToMacroAssignment[react$clusterAssignment]
    }

    updateNumericInput(session, "selectCluster", max=length(representatives))


    data = getData()
    data = data.frame(Index=seq_len(input$observation), data)
    select = which(assignment==input$selectCluster)
    select = select[!is.na(select)] ## filter faded texts

    ## filter for selected cluster
    data = data[select,]
    return(data)
  }


  output$assignment <- DT::renderDataTable({
    withProgress(message = 'Gathering Texts', value = 1, {
      getAssignmentData()
    })
  }, options=list(pageLength=100), rownames=FALSE)


  output$downloadAssignmentData <- downloadHandler(
    filename = function(){
      paste(input$observation, "_sankey.RData", sep="")
    },
    content = function(file){
      data = getAssignmentData()
      save(data, file=file)
    }
  )




  # Sankey ------------------------------------------------------------------
  getSankeyData = function(){
    req(length(react$microRepresentatives)>0)

    micro = react$microRepresentatives
    micro = sapply(micro, function(x){
      paste(names(x), collapse = ", ")
    })
    macro = react$macroRepresentatives
    macro = sapply(macro, function(x){
      paste(names(x), collapse = ", ")
    })
    microWeights = react$microWeights
    microToMacroAssignment = react$microToMacroAssignment
    clusterAssignment = react$clusterAssignment


    if(!input$sankey_showUnassigned){
      select = !is.na(microToMacroAssignment)
      selectData = !is.na(microToMacroAssignment[clusterAssignment])

      microToMacroAssignment = microToMacroAssignment[select]
      micro = micro[select]
      microWeights = microWeights[select]
    }


    if(input$sankey_microAndMacro || input$type=="micro"){
      ## get data
      data = getData()[,algorithm$RObj$textCol]

      if(!input$sankey_showUnassigned){
        clusterAssignment = clusterAssignment[selectData]
        clusterAssignment = as.numeric(as.factor(clusterAssignment)) ## re-encode e.g. 1,2,5 to 1,2,3
        data = data[selectData]
      }

      if(!is.na(input$sankey_trimTexts)){
        data = sapply(data, function(x){
          strtrim(x, input$sankey_trimTexts)
        })
      }
    }


    if(input$sankey_microAndMacro){
      labels = c(data, micro, macro)
      color = c(colormap[((seq_along(data)-1)%%length(colormap))+1], colormap[((seq_along(micro)-1)%%length(colormap))+1], colormap[((seq_along(macro)-1)%%length(colormap))+1])

      source = seq_len(length(data)+length(micro))-1
      target = c(length(data)+clusterAssignment-1, length(data)+length(micro)+microToMacroAssignment-1)
      value = c(rep(1, length(data)), microWeights)

      height=max(c(sum(!is.na(clusterAssignment))*20, 500))
    } else if(input$type=="micro"){
      labels = c(data, micro)
      color = c(colormap[((seq_along(data)-1)%%length(colormap))+1], colormap[((seq_along(micro)-1)%%length(colormap))+1])

      source = seq_len(length(data))-1
      target = c(length(data)+clusterAssignment-1)
      value = c(rep(1, length(data)))

      height=max(c(sum(!is.na(clusterAssignment))*20, 500))
    } else if(input$type=="macro"){
      labels = c(micro, macro)
      color = c(colormap[((seq_along(micro)-1)%%length(colormap))+1], colormap[((seq_along(macro)-1)%%length(colormap))+1])

      source = seq_len(length(micro))-1
      target = c(length(micro)+microToMacroAssignment-1)
      value = c(microWeights)

      height=max(c(sum(!is.na(microToMacroAssignment))*30, 500))
    }
    return(list(labels=labels, color=color, source=source, target=target, value=value, height=height))
  }

  drawSankey = function(){

    list2env(getSankeyData(), .GlobalEnv)

    plot_ly(
      type = "sankey",
      orientation = "h",
      arrangement ="freeform",

      node = list(
        label = labels,
        color = color,
        pad = 15,
        thickness = 20,
        line = list(
          color = color,
          width = 0.5
        )
      ),

      link = list(
        source = source,
        target = target,
        value =  value
      )
      , height = height)
  }

  output$sankey <- renderPlotly({
    withProgress(message = 'Drawing Sankey', value = 1, {
      drawSankey()
    })
  })

  output$downloadSankeyData <- downloadHandler(
    filename = function(){
      paste(input$observation, "_sankey.RData", sep="")
    },
    content = function(file){
      data = getSankeyData()
      save(data, file=file)
    }
  )


  output$downloadSankey <- downloadHandler(
    filename = function(){
      paste(input$observation, "_sankey.pdf", sep="")
    },
    content = function(file){
      p <- drawSankey()
      orca(p, file = "temp.pdf") ## bug: orca cannot handle paths
      file.rename("temp.pdf", file) ## therefore move to directory
    }
  )



  # Dendrogram ---------------------------------------------------------------
  getDendrogramData <- function(){
    req(length(react$macroRepresentatives)>0 && input$type == "macro")

    micro = react$microRepresentatives
    microToMacroAssignment = react$microToMacroAssignment

    micro = micro[!is.na(microToMacroAssignment)]
    micro = sapply(micro, function(x){
      paste(names(x), collapse = ", ")
    })

    hc <- react$dendrogram
    return(list(micro=micro, hc=hc))
  }

  drawDendrogram <- function(){

    list2env(getDendrogramData(), .GlobalEnv)
    class(hc) <- "hclust"
    plot(hc, labels=micro)
  }

  output$dendrogram <- renderPlot({
    withProgress(message = 'Drawing Dendrogram', value = 1, {
      drawDendrogram()
    })
  })

  output$downloadDendrogram <- downloadHandler(
    filename = function(){
      paste(input$observation, "_dendrogram.pdf", sep="")
    },
    content = function(file){
      pdf(file=file)
      drawDendrogram()
      dev.off()
    }
  )

  output$downloadDendrogramData <- downloadHandler(
    filename = function(){
      paste(input$observation, "_dendrogram.RData", sep="")
    },
    content = function(file){
      data = getDendrogramData()
      save(data, file = file)
    }
  )


})
