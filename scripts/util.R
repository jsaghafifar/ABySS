# install.packages(c("tidyverse","devtools","phytools","Matrix","TreeTools","phytools","seqinr","TeachingDemos"))
# library("devtools")
# devtools::install_github("walterxie/TraceR")

# validation pipeline for WCSS, stores coverage in CSV
coverage = function(patterns, path.in = ".", path.out = ".", states = 4) {
  require(tidyverse)
  require(phytools)
  require(TraceR)
  require(Matrix)
  require(TeachingDemos)
  if (substr(path.in,nchar(path.in),nchar(path.in)) != "/") {
    path.in = paste0(path.in,"/")
  }
  if (substr(path.out,nchar(path.out),nchar(path.out)) != "/") {
    path.out = paste0(path.out,"/")
  }
  if ("psi" %in% patterns && patterns[length(patterns)]!="psi") {
    cat("Tree parameters (psi) must be listed last in patterns.\n")
    next
  }
  
  # check what files are present in dir (after taking out any failed runs)
  name = gsub("-.*$", "", list.files(path=path.in, pattern="-([0-9]+).log")[1])
  last = sort(as.numeric(
    gsub("\\D", "", list.files(path = path.in, pattern = "-([0-9]+).log"))),
              decreasing = TRUE)[1]
  len = length(list.files(path = path.in, pattern = "-([0-9]+).log"))
  if(len!=last+1){ cat(len, "files are present out of", last+1, "\n") }
  
  allparams = colnames(read.delim(paste0(path.in,list.files(path = path.in, 
                                             pattern = "-([0-9]+).log")[1])))
  allparams = allparams[! allparams %in% c('Sample', 'posterior', 'likelihood', 
                                           'prior', 'D.treeLikelihood')]
  params = c()
  for (j in 1:length(patterns)) {
    params = c(params, allparams[grepl(paste0("^",patterns[j]), allparams)])
  }
  
  results = tibble()
  for (i in 0:last) {
    beastlog = file.path(paste0(path.in,name,"-", i,".log"))
    cat("\n\nStarting analysis", i+1, "out of",last+1,"...(",beastlog,")\n")
    # skip if file absent
    if(!file.exists(beastlog)) {
      cat("\nFile", i, "absent.\n")
      next
    }
    mcmc.log = readMCMCLog(beastlog)
    traces = TraceR::getTraces(mcmc.log, burn.in=0.1)
    stats = TraceR::analyseTraces(traces) %>% column_to_rownames(var = "trace")
    
    # extract HPD intervals from beast log
    median = stats %>% select(!!params) %>% 
      filter(row.names(stats) == "median")
    hpd.lo = stats %>% select(!!params) %>% 
      filter(row.names(stats) == "HPD95.lower")
    hpd.up = stats %>% select(!!params) %>% 
      filter(row.names(stats) == "HPD95.upper")
    ess = stats %>% select(!!params) %>% 
      filter(row.names(stats) == "ESS") %>% mutate_if(is.character, as.numeric)
    
    # combine all posterior stats here
    df = tibble(Parameter = names(hpd.lo), HPD95.lower = hpd.lo %>% unlist)
    df$HPD95.upper = hpd.up %>% unlist
    df$Median = median %>% unlist
    df$ESS = as.numeric(ess) %>% unlist
    
    # true value
    tru.file = file.path(paste0(path.in,name,"-", i,"_true.log"))
    tru.fil = read_delim(tru.file, delim="\t", comment = "#")
    tru = tru.fil %>% select(starts_with(patterns)) %>% unlist()
    
    if (!is.null(tru)) {
      names(tru)[grepl("^freq",names(tru))] = params[grepl("^freq",params)]
      names(tru)[grepl("^rates",names(tru))] = params[grepl("^rates",params)]
      names(tru)[grepl("^nqRates",names(tru))]=params[grepl("^nqRates",params)]
      names(tru)[grepl("^qRates",names(tru))] = params[grepl("^qRates",params)]
      names(tru)[grepl("^qFreq",names(tru))] = params[grepl("^qFreq",params)]
      
      if (any(grepl("equilibrium.freq",patterns))) {
        tru.nq = tru.fil %>% select(starts_with(c('Q_','nonrevQ_')))%>%unlist()
        Qmat = matrix(tru.nq,nrow = states, ncol = states)
        P = expm(Qmat*10000)
        tru.nf = c()
        for (j in 1:states) {
          tru.nf[j] = P[j,1]
        }
        names(tru.nf) = params[grepl("equilibrium.freq",params)]
        if (patterns[1]!="freq") {tru=c(tru.nf, tru)} else {tru=c(tru,tru.nf)}
      }
      
    }
    
    if (any(grepl("psi",params))) {
      # true tree
      tru.tre.f <- file.path(paste0(path.in,name,"-", i,"_true_psi.trees"))
      tru.tre <- read.nexus(tru.tre.f)
      tru["psi.height"] <- max(nodeHeights(tru.tre))
      tru["psi.treeLength"] <- sum(tru.tre$edge.length)
    }
    
    # combine all truth here
    tru.df <- tibble(Parameter = params, True.Value = tru)
    
    # merge posterior dataframe with true
    df <- merge(df, tru.df)
    stopifnot(nrow(df) == length(params))
    df_sorted <- df %>%
      mutate(Parameter = factor(Parameter, levels = params)) %>%
      arrange(Parameter) %>%
      mutate(i = i)
    
    # add ith processed result into the final dataframe
    results <- bind_rows(results, df_sorted)
  }
  
  # Convert HPD95.lower and HPD95.upper from character to numeric
  results <- results %>%  mutate(
    i = i,
    HPD95.lower = as.numeric(HPD95.lower),
    HPD95.upper = as.numeric(HPD95.upper),
    Median = as.numeric(Median),
    is.in = True.Value >= HPD95.lower & True.Value <= HPD95.upper, 
    ESS = as.numeric(ESS)
  )
  
  file.pre = ifelse(path.in!="./",basename(path.in),basename(getwd()))
  for (j in 1:length(patterns)) {
    write.csv(results[grepl(paste0("^",patterns[j]),results$Parameter),], 
             paste0(path.out,file.pre,"-",patterns[j],"-coverage.csv"))
  }
  
  # coverage
  cov = results %>%
    group_by(Parameter) %>%        # Group by Parameter
    summarize(Coverage = sum(is.in))
  cov[,"Coverage"] = cov[,"Coverage"]*100/len # in %
  write.csv(cov,paste0(path.out,file.pre,"-coverage.csv"))
  
  # check if parameters are mixed (>=200 effective sample size)
  ESS = results %>%
        group_by(Parameter) %>% 
        summarize(is.mixed = sum(ESS>=200)*100/len)
  write.csv(ESS,paste0(path.out,file.pre,"-ess.csv"))
  
}

# extract quantified time-nonreversibility from logs
nonrev = function(path.in = ".", path.out = ".") {
  require(TraceR)
  require(tidyverse)
  if (substr(path.in,nchar(path.in),nchar(path.in)) != "/") {
    path.in = paste0(path.in,"/")
  }
  if (substr(path.out,nchar(path.out),nchar(path.out)) != "/") {
    path.out = paste0(path.out,"/")
  }
  file.pre = ifelse(path.in!=".",basename(path.in),basename(getwd()))
  files = paste0(path.in,list.files(pattern = ".+_[Nn][Qq].log", path=path.in))
  nonreversibility = tibble()
  for (f.nq in files) {
    name = gsub(".+/","",gsub("_.+","",f.nq))
    allparams = colnames(read.delim(f.nq, comment.char = "#")) 
    params = c(allparams[grepl("balance", allparams)],
               allparams[grepl("flux", allparams)])
    
    log = readMCMCLog(f.nq)
    traces = TraceR::getTraces(log, burn.in=0.1)
    stats = TraceR::analyseTraces(traces) %>% 
      column_to_rownames(var = "trace")
    median = stats %>% select(!!params) %>%
      filter(row.names(stats) == "median") %>% 
      mutate_if(is.character, as.numeric) %>% unlist()
    hpd.lo = stats %>% select(!!params) %>% 
      filter(row.names(stats) == "HPD95.lower") %>% 
      mutate_if(is.character, as.numeric) %>% unlist()
    hpd.up = stats %>% select(!!params) %>% 
      filter(row.names(stats) == "HPD95.upper") %>% 
      mutate_if(is.character, as.numeric) %>% unlist()
    
    df = tibble(Parameter = names(hpd.lo),
                   HPD95.lower = hpd.lo %>% unlist)
    df$HPD95.upper = hpd.up %>% unlist
    df$Median = median %>% unlist
    
    df = df %>%
      mutate(Parameter = factor(Parameter, levels = params)) %>%
      arrange(Parameter) %>%
      mutate(Name = name) %>% 
      mutate(Method = ifelse(grepl("balance",Parameter),"Detailed Balance",
                             "Net Flux")) %>%
      mutate(Directory = file.pre)
    nonreversibility = bind_rows(nonreversibility, df)
    cat("File",name,"in",file.pre,"finished time-nonreversibility analysis.\n")
  }
  nonreversibility = nonreversibility %>% arrange(Name)
  write.csv(nonreversibility,paste0(path.out,file.pre,"-nonreversibility.csv"))
}

# extract RMS of deviation between all NQs and Qs in a folder
allDeviationsQ = function(path.in = ".", path.out = ".", states = 4) {
  require(TraceR)
  require(Matrix)
  require(tidyverse)
  if (substr(path.in,nchar(path.in),nchar(path.in)) != "/") {
    path.in = paste0(path.in,"/")
  }
  if (substr(path.out,nchar(path.out),nchar(path.out)) != "/") {
    path.out = paste0(path.out,"/")
  }
  file.pre = ifelse(path.in!=".",basename(path.in),basename(getwd()))
  files.nq = paste0(path.in,list.files(pattern = ".+_NQ.log", path=path.in))
  files.q = paste0(path.in,list.files(pattern = ".+_Q.log", path=path.in))
  deviations = tibble()
  for (i in 1:length(files.nq)) {
    name = gsub(".+/","",gsub("_.+","",files.nq[i]))
    nq = constructQ(files.nq[i],states=states,is.rev=F)
    q = constructQ(files.q[i],states=states,is.rev=T)
    sum.diffs = 0
    x = 0
    for (i2 in 1:states) {
      for (j in 1:states) {
        if (i2==j) {
          next
        }
        sum.diffs = sum.diffs + (nq[i2,j]-q[i2,j])^2
        x = x + 1
      }
    }
    dev = sqrt(sum.diffs/x)
    df = tibble(Name = name,
                Deviation = dev,
                Directory = file.pre)
    deviations = bind_rows(deviations,df)
    cat("File",i,"in",file.pre,"finished Q deviation analysis.\n")
  }
  deviations = deviations %>% arrange(Name)
  write.csv(deviations,paste0(path.out,file.pre,"-deviations.csv"))
}

# extract root clade splits from NQ and Q to compare with truth
roots = function(path.in = ".", path.out = ".") {
  library(phytools)
  library(TreeTools)
  library(tidyverse) 
  if (substr(path.in,nchar(path.in),nchar(path.in)) != "/") {
    path.in = paste0(path.in,"/")
  }
  if (substr(path.out,nchar(path.out),nchar(path.out)) != "/") {
    path.out = paste0(path.out,"/")
  }
  file.pre = ifelse(path.in!=".",basename(path.in),basename(getwd()))
  
  files.tru = paste0(path.in,list.files(pattern = ".+_true_psi.trees", path=path.in))
  files.nq = paste0(path.in,list.files(pattern = ".+_NQ.trees", path=path.in))
  files.q = paste0(path.in,list.files(pattern = ".+_Q.trees", path=path.in))
  
  results = tibble()
  for (i in 1:length(files.tru)) {
    tree.tru = readNexus(files.tru[i])
    trees.nq = readNexus(files.nq[i])
    trees.q = readNexus(files.q[i])
    
    # get what tips are in one clade at the root split
    true.split = extract.clade(tree.tru,RootNode(tree.tru)+1)$tip.label
    
    df = tibble()
    # skip initial random tree and 10% burnin
    for (j in ((length(trees.nq)-1)/10):length(trees.nq)) { 
      tree.nq = trees.nq[[j]]
      tree.q = trees.q[[j]]
      
      nq.split = extract.clade(tree.nq,RootNode(tree.nq)+1)$tip.label
      q.split = extract.clade(tree.q,RootNode(tree.q)+1)$tip.label
      
      # if tips match exactly OR if tips don't match exactly
      nq.root = ifelse(all(nq.split %in% true.split) && all(true.split %in% nq.split) || 
                         all(!(nq.split %in% true.split)) && all(!(true.split %in% nq.split)), T, F)
      q.root = ifelse(all(q.split %in% true.split) && all(true.split %in% q.split) || 
                        all(!(q.split %in% true.split)) && all(!(true.split %in% q.split)), T, F) 
      
      root = tibble(Simulation = paste0(file.pre,"-",i), i = i, Tree = j, 
                    CorrectNQ = nq.root, CorrectQ = q.root)
      root$i = root$i %>% as.numeric()
      root$Tree = root$Tree %>% as.numeric()
      df = bind_rows(df, root)
      
    }
    
    results = bind_rows(results, df)
    cat("File",i,"in",file.pre,"finished root analysis.\n")
  }
  results = results %>% arrange(Simulation)
  write.csv(results,paste0(path.out,file.pre,"-roots_raw.csv"))
}

# sub-sampling sequences
sample_month = function(n, file) {
  require(seqinr)
  sample = c()
  total = 0
  fasta = read.fasta(file)
  names = names(fasta)
  for (year in 2013:2022) {
    for (month in 1:12) {
      month_data = if (month>9)
        grep(paste0(year,"-",month),names) else
          grep(paste0(year,"-0",month),names)
    if (length(month_data) != 0 & length(month_data) >= n) {
        month_sample = sample(fasta[month_data], size = n)
        sample = c(sample, month_sample)
        total = total + n
      }
    if (length(month_data) != 0 & length(month_data) < n) {
      sample = c(sample, fasta[month_data])
      total = total + length(month_data)
    }
    }
  }
  write.fasta(sample, names(sample), 
              file.out = gsub(".fasta",paste0("_",total,".fasta"),file))
  cat(total,"sequences.\n")
}

# extract any parameter's HPD and median from logs
extractTrace = function(patterns, path.in = ".", path.out = ".") {
  require(TraceR)
  require(tidyverse)
  if (substr(path.in,nchar(path.in),nchar(path.in)) != "/") {
    path.in = paste0(path.in,"/")
  }
  if (substr(path.out,nchar(path.out),nchar(path.out)) != "/") {
    path.out = paste0(path.out,"/")
  }
  file.pre = ifelse(path.in!=".",basename(path.in),basename(getwd()))
  files = paste0(path.in,list.files(pattern = ".log", path=path.in))
  results = tibble()
  for (f in files) {
    name = gsub(".+/","",gsub("\\.log","",f))
    allparams = colnames(read.delim(f, comment.char = "#")) 
    params = c()
    for (j in 1:length(patterns)) {
      params = c(params, allparams[grepl(paste0("^",patterns[j]), allparams)])
    }
    if (length(params)==0) {
      cat("Parameters not found for file",name,"in",file.pre,"\n")
      next
    }
    
    log = readMCMCLog(f)
    traces = TraceR::getTraces(log, burn.in=0.1)
    stats = TraceR::analyseTraces(traces) %>% 
      column_to_rownames(var = "trace")
    median = stats %>% select(!!params) %>%
      filter(row.names(stats) == "median") %>% 
      mutate_if(is.character, as.numeric) %>% unlist()
    mean = stats %>% select(!!params) %>%
      filter(row.names(stats) == "mean") %>% 
      mutate_if(is.character, as.numeric) %>% unlist()
    hpd.lo = stats %>% select(!!params) %>% 
      filter(row.names(stats) == "HPD95.lower") %>% 
      mutate_if(is.character, as.numeric) %>% unlist()
    hpd.up = stats %>% select(!!params) %>% 
      filter(row.names(stats) == "HPD95.upper") %>% 
      mutate_if(is.character, as.numeric) %>% unlist()
    
    df = tibble(Parameter = names(hpd.lo),
                HPD95.lower = hpd.lo %>% unlist)
    df$HPD95.upper = hpd.up %>% unlist
    df$Median = median %>% unlist
    df$mean = mean %>% unlist
    
    df = df %>%
      mutate(Parameter = factor(Parameter, levels = params)) %>%
      arrange(Parameter) %>%
      mutate(Name = name) %>% 
      mutate(Directory = file.pre)
    results = bind_rows(results, df)
    cat("File",name,"in",file.pre,"finished output analysis.\n")
  }
  results = results %>% arrange(Name)
  write.csv(results,paste0(path.out,file.pre,"-results.csv"))
}

# extract site-mixture compositions from a log
siteMixture = function(file, partitionIndex = NULL, partitionName = "partition",
                       path.out = ".") {
  require(TraceR)
  if (substr(path.out,nchar(path.out),nchar(path.out)) != "/") {
    path.out = paste0(path.out,"/")
  }
  mcmc.log = readMCMCLog(file)
  traces = TraceR::getTraces(mcmc.log, burn.in=0.1)
  stats = TraceR::analyseTraces(traces) %>% column_to_rownames(var = "trace")
  params = colnames(stats)[grepl("^siteLikelihood", colnames(stats))]
  
  mean = stats %>% select(!!params) %>% 
    filter(row.names(stats) == "mean") %>% 
    mutate_if(is.character, as.numeric) %>% unlist
  hpd.lo = stats %>% select(!!params) %>% 
    filter(row.names(stats) == "HPD95.lower") %>% 
    mutate_if(is.character, as.numeric) %>% unlist
  hpd.up = stats %>% select(!!params) %>% 
    filter(row.names(stats) == "HPD95.upper") %>% 
    mutate_if(is.character, as.numeric) %>% unlist
  ess = stats %>% select(!!params) %>% 
    filter(row.names(stats) == "ESS") %>% 
    mutate_if(is.character, as.numeric) %>% unlist
  
  if (is.null(partitionIndex)) {
    partitionIndex = rep("Site",length(params))
  }
  
  df = tibble(Site = c(1:length(params)),
              Partition = partitionIndex,
              Mean = mean,
              HPD.Lower = hpd.lo,
              HPD.Upper = hpd.up,
              is.mixed = ess>=200)
  name = gsub(".+/","",gsub("_.+","",file))
  write.csv(df,paste0(path.out,name,"-",partitionName,".csv"))
  
}

# construct Q-matrices from a log
constructQ = function(file, states = 4, models = 1, is.rev = T) {
  require(TraceR)
  require(Matrix)
  if (length(is.rev) != models) {
    cat("Warning: assuming all models are",ifelse(is.rev,"time-reversible.","time-nonreversible."),
        "Specify is.rev as a vector if not.\n")
    is.rev = rep(is.rev,models)
  }
  
  name = gsub(".+/","",gsub("_.+","",file))
  
  if (grepl('.log',file)) {
    mcmc.log = readMCMCLog(file)
    traces = TraceR::getTraces(mcmc.log, burn.in=0.1)
    stats = TraceR::analyseTraces(traces) %>% column_to_rownames(var = "trace")
  } else {
    stats = read.csv(file) %>% column_to_rownames(var = "Parameter") %>% t() %>% as_tibble(rownames = NA)
  }
  params = colnames(stats)
    
  state.names = if (states == 4) {
    c("A", "C", "G", "T") 
  } else if (states == 20) { 
    c('A','C','D','E','F','G','H','I','K','L',
      'M','N','P','Q','R','S','T','V','W','Y') 
    } else {
      as.character(1:states)
      }
  Qr.param = params[grepl("^rates",params) |
                      grepl("^Q.+ates",params) |
                      grepl("^N.+rates",params)]
  Qf.param = params[grepl("^freq",params) |
                      grepl("^Q.+req",params)]
  x.r = 1
  x.f = 1
  y.f = 1
  Q.list = list()
  
  for (mod in 1:models) {
    n.r = states*states-states
    if (is.rev[mod]) n.r = n.r/2
    y.r = x.r+n.r
    Qr.p = Qr.param[x.r:(y.r-1)]
    
    if (is.rev[mod]) {
      y.f = x.f+states
      Qf.p = Qf.param[x.f:(y.f-1)]
      Qf = stats %>% select(!!Qf.param) %>% 
        filter(row.names(stats) == "mean") %>% unlist() %>% as.numeric()
      names(Qf) = state.names
      Qf = Qf / sum(Qf)
    } else {
      Qf = NA
    }
    
    Qm = matrix(NA, nrow = states, ncol = states,
                dimnames = list(state.names, state.names))
    for (i in 1:states) {
      for (j in 1:states) {
        if (i == j){ next }
        # Construct the column name for pair i,j
        rate.name = paste0(state.names[i],".",state.names[j])
        if (any(grepl(rate.name,Qr.p))) {
          Qr.colname = Qr.p[grepl(rate.name,Qr.p)]
        } else {
          if (!is.rev[mod]) {
            cat("bad")
            next
          }
          rate.name2 = paste0(state.names[j],".",state.names[i])
          if (any(grepl(rate.name2,Qr.p))) {
            Qr.colname = Qr.p[grepl(rate.name2,Qr.p)]
          } else {
            cat("super bad")
            next
          }
        }
        # Calculate the mean 
        Qr.mean = stats %>% select(!!Qr.colname) %>% 
          filter(row.names(stats) == "mean") %>% unlist() %>% as.numeric()
        
        if (is.rev[mod]) {
          Qm[state.names[i], state.names[j]] = Qr.mean*Qf[j]
        } else {
          # if nonreversible, do not include freq in Qmat
          Qm[state.names[i], state.names[j]] = Qr.mean
        }
      }
    }
    
    # Construct diagonals
    for (i in 1:states) {
      Qm[i,i] = 0
      for (j in 1:states) {
        if (i==j) { next }
        Qm[i,i] = Qm[i,i] - Qm[i,j]
      }
    }
    
    # Get NQ equilibrium frequencies
    if (!is.rev[mod]) {
      P = expm(Qm*100000)
      Qf = P[1,1:states]
    }
    
    # Normalise matrix
    subst = 0
    for (i in 1:states) {
      subst = subst - Qm[i,i]*Qf[i]
    }
    Qm = Qm / subst
    
    Q.list[[mod]] = Qm
    names(Q.list)[mod] = paste0(name,"-",mod,"_",ifelse(is.rev[mod],"Q","NQ"))
    
    x.r = y.r
    x.f = y.f
  }
  if (models==1) {
    Q.list[[1]]
  } else {
    Q.list
  }
}

constructQdiffs = function(Q1, Q2) {
  n = length(Q1[1,])
  Qdiffs = Q1
  for (i in 1:n) {
    for (j in 1:n) {
      Qdiffs[i,j] = Q1[i,j]-Q2[i,j]
    }
  }
  Qdiffs
}

# extract RMS of deviation between two Q-matrices
deviationQ = function(Q1, Q2) {
  sum.diffs = 0
  x = 0
  n = length(Q1[1,])
  for (i in 1:n) {
    for (j in 1:n) {
      if (i==j) {
        next
      }
      sum.diffs = sum.diffs + (Q1[i,j]-Q2[i,j])^2
      x = x + 1
    }
  }
  dev = sqrt(sum.diffs/x)
  dev
}

# plot a Q-matrix
plotQ = function(Q, name = 'Q-matrix', title = "Q-matrix", path.out = ".", monochrome = T) {
  if (substr(path.out,nchar(path.out),nchar(path.out)) != "/") {
    path.out = paste0(path.out,"/")
  }
  n = length(Q[1,])
  png.name = paste0(path.out,name,".png")
  if (n==4) {
    states = c('A','C','G','T')
  } else {
    states = c('W','Y','F','I','V','L','M','D','E','Q',
               'N','H','K','R','S','T','A','C','G','P') # group by AA chem props
  }
  png(png.name, res=300, width=2400, height=2400)
  par(mar=c(4,4,2,2),cex=ifelse(n==4,2,0.8),family="Noto Serif")
  
  if (monochrome) {
    COLS = rev(c("#6f0000", "darkred", "#a23333", "#b96666", 
                 "#d19999", "#e8cccc", "gray97"))
    FONT.COLS = rev(c("gray100", "gray100", "gray100", "gray100", 
                      "gray100", "black", "black"))
    CUTOFFS = if (n==4) {
      c(0, 1/20, 1/10, 1/5, 3/4, 4/3, 2)
    } else {
      c(0, 1/100, 1/20, 1/10, 1/4, 1/2, 1)
    }
  } else {
    COLS = rev(c("#6f0000", "#a23333", "#d19999", "gray97",
                 "#99D1D1", "#33A2A2", "#006F6F"))
    FONT.COLS = rev(c("gray100", "gray100", "black", "black", 
                      "black", "gray100", "gray100"))
    CUTOFFS = c(-1/10, -1/20, -1/100, 0, 1/100, 1/20, 1/10)
  }
  
  # Plot the matrix
  plot(0,0, type="n", xlim=c(0.5, (n+1)), ylim = c(0.5,(n+0.5)), axes=F, 
       xlab="To", ylab="From", xaxs="i", yaxs="i", 
       main = title)
  
  for (i in 1:n) {
    for (j in 1:n) { 
      x = j
      y = n-i+1
      if (i == j){
        rect(x-0.5, y-0.5, x+0.5, y+0.5, col="black")
        text(x, y, states[i], cex=2, adj=c(0.5, 0.5), col="white")
      }else{
        rate = Q[states[i], states[j]]
        
        index = which(rate >= CUTOFFS)
        index = index[length(index)]
        col = COLS[index]
        font.col = FONT.COLS[index]
        
        rect(x-0.5, y-0.5, x+0.5, y+0.5, col = col)
        text(x, y, signif(rate, 2), cex=ifelse(states==4,1,0.8), adj=c(0.5, 0.5), col=font.col)
      }
    }
  }
  
  axis(1, at = 1:n, labels = states)
  axis(2, at = 1:n, labels = rev(states), las=2)
  
  axis(1, at = 0:(n+1), labels = rep("", (n+2)), lwd.ticks=0)
  axis(2, at = 0:(n+1), labels = rep("", (n+2)), lwd.ticks=0)
  
  dev.off()
  cat(name,"plotted.\n")
}
