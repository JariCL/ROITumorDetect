
#################################################################################
# CODE PART 2: modeling, selection of variables, prediction visualization, etc. #
#################################################################################

list.of.names <- paste0("sample_",
                        LETTERS[1:8])


# Variable selection -----

# Function that does forward selection, choosing the best fitting variable in each step according to the deviance of the model:
variable.selection.f.bin <- function(df1){
  
  # t: # tumor cells, nt: # non-tumor cells
  t <- df1$t; nt <- df1$nt
  data <- df1[,1:(ncol(df1)-2)]
  
  cat('\n Subtype 1 -- Variable 1')
  
  # loop over all variables, fit them univariately, save there deviance, select best:
  coef <- sapply(1:ncol(data), function(var){
    cat(var, "---")
    x <- data.matrix(data[,var])
    mod <- glm(cbind(t,nt) ~ x, family = "binomial")
    return(mod$deviance)
  })
  df <- data.frame(var = 1:ncol(data), coef = abs(coef)); df <- na.omit(df); df <- df %>% arrange(coef)
  var.incl <- df[1,1]
  colnames(b.df)[var.incl]
  
  
  cat('\n Subtype 1 -- Variable 2')
  # remove previously already chosen variable from set of possible variables to loop over, and start process again:
  coef <- sapply((1:ncol(data))[-var.incl], function(var){
    cat(var, "---")
    x <- data.matrix(data[,c(var, var.incl)])
    mod <- glm(cbind(t,nt) ~ x, family = "binomial")
    return(mod$deviance)
  })
  df <- data.frame(var = (1:ncol(data))[-var.incl], coef = abs(coef)); df <- na.omit(df); df <- df %>% arrange(coef)
  K <- T
  i <- 0
  while(K){
    i <- i + 1
    cor.df <- cor(b.df[,c(var.incl, df[i,1])])
    K <- ifelse(sum(abs(cor.df[nrow(cor.df),1:(ncol(cor.df)-1)])>.5)>0, T, F)
  }
  cor(b.df[,c(var.incl, df[i,1])])
  var.incl <- c(var.incl, df[i,1])
  
  cat('\n Subtype 1 -- Variable 3')
  coef <- sapply((1:ncol(data))[-var.incl], function(var){
    cat(var, "---")
    x <- data.matrix(data[,c(var, var.incl)])
    mod <- glm(cbind(t,nt) ~ x, family = "binomial")
    return(mod$deviance)
  })
  df <- data.frame(var = (1:ncol(data))[-var.incl], coef = abs(coef)); df <- na.omit(df); df <- df %>% arrange(coef)
  K <- T
  i <- 0
  while(K){
    i <- i + 1
    cor.df <- cor(b.df[,c(var.incl, df[i,1])])
    K <- ifelse(sum(abs(cor.df[nrow(cor.df),1:(ncol(cor.df)-1)])>.5)>0, T, F)
  }
  cor(b.df[,c(var.incl, df[i,1])])
  var.incl <- c(var.incl, df[i,1])
  
  
  cat('\n Subtype 1 -- Variable 4')
  coef <- sapply((1:ncol(data))[-var.incl], function(var){
    cat(var, "---")
    x <- data.matrix(data[,c(var, var.incl)])
    mod <- glm(cbind(t,nt) ~ x, family = "binomial")
    return(mod$deviance)
  })
  df <- data.frame(var = (1:ncol(data))[-var.incl], coef = abs(coef)); df <- na.omit(df); df <- df %>% arrange(coef)
  K <- T
  i <- 0
  while(K){
    i <- i + 1
    cor.df <- cor(b.df[,c(var.incl, df[i,1])])
    K <- ifelse(sum(abs(cor.df[nrow(cor.df),1:(ncol(cor.df)-1)])>.5)>0, T, F)
  }
  cor(b.df[,c(var.incl, df[i,1])])
  var.incl <- c(var.incl, df[i,1])
  
  cat('\n Subtype 1 -- Variable 5')
  coef <- sapply((1:ncol(data))[-var.incl], function(var){
    cat(var, "---")
    x <- data.matrix(data[,c(var, var.incl)])
    mod <- glm(cbind(t,nt) ~ x, family = "binomial")
    return(mod$deviance)
  })
  df <- data.frame(var = (1:ncol(data))[-var.incl], coef = abs(coef)); df <- na.omit(df); df <- df %>% arrange(coef)
  K <- T
  i <- 0
  while(K){
    i <- i + 1
    cor.df <- cor(b.df[,c(var.incl, df[i,1])])
    K <- ifelse(sum(abs(cor.df[nrow(cor.df),1:(ncol(cor.df)-1)])>.5)>0, T, F)
  }
  cor(b.df[,c(var.incl, df[i,1])])
  var.incl <- c(var.incl, df[i,1])
  var.incl1 <- var.incl
  
  
  
  # Following code starts not from best, but from second best (look at 11th line from this line)
  t <- df1$t; nt <- df1$nt
  data <- df1[,1:(ncol(df1)-2)]
  cat('\n Subtype 2 -- Variable 1')
  coef <- sapply(1:ncol(data), function(var){
    cat(var, "---")
    x <- data.matrix(data[,var])
    mod <- glm(cbind(t,nt) ~ x, family = "binomial")
    return(mod$deviance)
  })
  df <- data.frame(var = 1:ncol(data), coef = abs(coef)); df <- na.omit(df); df <- df %>% arrange(coef)
  var.incl <- df[2,1] # Here, the second best variable is chosen as starting point for further selection process
  colnames(b.df)[var.incl]
  
  cat('\n Subtype 2 -- Variable 2')
  coef <- sapply((1:ncol(data))[-var.incl], function(var){
    cat(var, "---")
    x <- data.matrix(data[,c(var, var.incl)])
    mod <- glm(cbind(t,nt) ~ x, family = "binomial")
    return(mod$deviance)
  })
  df <- data.frame(var = (1:ncol(data))[-var.incl], coef = abs(coef)); df <- na.omit(df); df <- df %>% arrange(coef)
  K <- T
  i <- 0
  while(K){
    i <- i + 1
    cor.df <- cor(b.df[,c(var.incl, df[i,1])])
    K <- ifelse(sum(abs(cor.df[nrow(cor.df),1:(ncol(cor.df)-1)])>.5)>0, T, F)
  }
  cor(b.df[,c(var.incl, df[i,1])])
  var.incl <- c(var.incl, df[i,1])
  
  cat('\n Subtype 2 -- Variable 3')
  # var.incl <- var.incl[1:2]
  coef <- sapply((1:ncol(data))[-var.incl], function(var){
    cat(var, "---")
    x <- data.matrix(data[,c(var, var.incl)])
    mod <- glm(cbind(t,nt) ~ x, family = "binomial")
    return(mod$deviance)
  })
  df <- data.frame(var = (1:ncol(data))[-var.incl], coef = abs(coef)); df <- na.omit(df); df <- df %>% arrange(coef)
  K <- T
  i <- 0
  while(K){
    i <- i + 1
    cor.df <- cor(b.df[,c(var.incl, df[i,1])])
    K <- ifelse(sum(abs(cor.df[nrow(cor.df),1:(ncol(cor.df)-1)])>.5)>0, T, F)
  }
  cor(b.df[,c(var.incl, df[i,1])])
  var.incl <- c(var.incl, df[i,1])
  
  
  cat('\n Subtype 2 -- Variable 4')
  coef <- sapply((1:ncol(data))[-var.incl], function(var){
    cat(var, "---")
    x <- data.matrix(data[,c(var, var.incl)])
    mod <- glm(cbind(t,nt) ~ x, family = "binomial")
    return(mod$deviance)
  })
  df <- data.frame(var = (1:ncol(data))[-var.incl], coef = abs(coef)); df <- na.omit(df); df <- df %>% arrange(coef)
  K <- T
  i <- 0
  while(K){
    i <- i + 1
    cor.df <- cor(b.df[,c(var.incl, df[i,1])])
    K <- ifelse(sum(abs(cor.df[nrow(cor.df),1:(ncol(cor.df)-1)])>.5)>0, T, F)
  }
  cor(b.df[,c(var.incl, df[i,1])])
  var.incl <- c(var.incl, df[i,1])
  
  cat('\n Subtype 2 -- Variable 5')
  coef <- sapply((1:ncol(data))[-var.incl], function(var){
    cat(var, "---")
    x <- data.matrix(data[,c(var, var.incl)])
    mod <- glm(cbind(t,nt) ~ x, family = "binomial")
    return(mod$deviance)
  })
  df <- data.frame(var = (1:ncol(data))[-var.incl], coef = abs(coef)); df <- na.omit(df); df <- df %>% arrange(coef)
  K <- T
  i <- 0
  while(K){
    i <- i + 1
    cor.df <- cor(b.df[,c(var.incl, df[i,1])])
    K <- ifelse(sum(abs(cor.df[nrow(cor.df),1:(ncol(cor.df)-1)])>.5)>0, T, F)
  }
  cor(b.df[,c(var.incl, df[i,1])])
  var.incl <- c(var.incl, df[i,1])
  var.incl2 <- var.incl
  
  
  
  
  
  # Following code starts not from best, but from third best (look at 11th line from this line)
  t <- df1$t; nt <- df1$nt
  data <- df1[,1:(ncol(df1)-2)]
  cat('\n Subtype 3 -- Variable 1')
  coef <- sapply(1:ncol(data), function(var){
    cat(var, "---")
    x <- data.matrix(data[,var])
    mod <- glm(cbind(t,nt) ~ x, family = "binomial")
    return(mod$deviance)
  })
  df <- data.frame(var = 1:ncol(data), coef = abs(coef)); df <- na.omit(df); df <- df %>% arrange(coef)
  var.incl <- df[3,1]  # Here, the third best variable is chosen as starting point for further selection process
  colnames(b.df)[var.incl]
  
  cat('\n Subtype 3 -- Variable 2')
  coef <- sapply((1:ncol(data))[-var.incl], function(var){
    cat(var, "---")
    x <- data.matrix(data[,c(var, var.incl)])
    mod <- glm(cbind(t,nt) ~ x, family = "binomial")
    return(mod$deviance)
  })
  df <- data.frame(var = (1:ncol(data))[-var.incl], coef = abs(coef)); df <- na.omit(df); df <- df %>% arrange(coef)
  K <- T
  i <- 0
  while(K){
    i <- i + 1
    cor.df <- cor(b.df[,c(var.incl, df[i,1])])
    K <- ifelse(sum(abs(cor.df[nrow(cor.df),1:(ncol(cor.df)-1)])>.5)>0, T, F)
  }
  cor(b.df[,c(var.incl, df[i,1])])
  var.incl <- c(var.incl, df[i,1])
  
  cat('\n Subtype 3 -- Variable 3')
  # var.incl <- var.incl[1:2]
  coef <- sapply((1:ncol(data))[-var.incl], function(var){
    cat(var, "---")
    x <- data.matrix(data[,c(var, var.incl)])
    mod <- glm(cbind(t,nt) ~ x, family = "binomial")
    return(mod$deviance)
  })
  df <- data.frame(var = (1:ncol(data))[-var.incl], coef = abs(coef)); df <- na.omit(df); df <- df %>% arrange(coef)
  K <- T
  i <- 0
  while(K){
    i <- i + 1
    cor.df <- cor(b.df[,c(var.incl, df[i,1])])
    K <- ifelse(sum(abs(cor.df[nrow(cor.df),1:(ncol(cor.df)-1)])>.5)>0, T, F)
  }
  cor(b.df[,c(var.incl, df[i,1])])
  var.incl <- c(var.incl, df[i,1])
  
  
  cat('\n Subtype 3 -- Variable 4')
  coef <- sapply((1:ncol(data))[-var.incl], function(var){
    cat(var, "---")
    x <- data.matrix(data[,c(var, var.incl)])
    mod <- glm(cbind(t,nt) ~ x, family = "binomial")
    return(mod$deviance)
  })
  df <- data.frame(var = (1:ncol(data))[-var.incl], coef = abs(coef)); df <- na.omit(df); df <- df %>% arrange(coef)
  K <- T
  i <- 0
  while(K){
    i <- i + 1
    cor.df <- cor(b.df[,c(var.incl, df[i,1])])
    K <- ifelse(sum(abs(cor.df[nrow(cor.df),1:(ncol(cor.df)-1)])>.5)>0, T, F)
  }
  cor(b.df[,c(var.incl, df[i,1])])
  var.incl <- c(var.incl, df[i,1])
  
  cat('\n Subtype 3 -- Variable 5')
  coef <- sapply((1:ncol(data))[-var.incl], function(var){
    cat(var, "---")
    x <- data.matrix(data[,c(var, var.incl)])
    mod <- glm(cbind(t,nt) ~ x, family = "binomial")
    return(mod$deviance)
  })
  df <- data.frame(var = (1:ncol(data))[-var.incl], coef = abs(coef)); df <- na.omit(df); df <- df %>% arrange(coef)
  K <- T
  i <- 0
  while(K){
    i <- i + 1
    cor.df <- cor(b.df[,c(var.incl, df[i,1])])
    K <- ifelse(sum(abs(cor.df[nrow(cor.df),1:(ncol(cor.df)-1)])>.5)>0, T, F)
  }
  cor(b.df[,c(var.incl, df[i,1])])
  var.incl <- c(var.incl, df[i,1])
  var.incl3 <- var.incl
  
  return(list(var.incl1, var.incl2, var.incl3)) 
}


# Function that collects all variables per grid cell + the response variable (tumor counts vs. non-tumor counts):
make.df.bin <- function(p.name, nx = 25){
  name1 <- paste0("https://raw.githubusercontent.com/JariCL/ROITumorDetect/refs/heads/main/nx", nx,
                  "/complete_variables_nx_", nx, "_", p.name, ".txt")
  
  read.table(file = name1)
}

# Function that executes the two previous functions for all specified samples:
make.all.sel.var <- function(nx = 25){
  
  
  df.big <- do.call(rbind, lapply(list.of.names, function(pattern.name){
    cat(pattern.name, " - ")
    make.df.bin(pattern.name, nx = nx)
  }))
  
  
  df <- df.big
  
  which.na.per.sample <- lapply(list.of.names, function(x){
    sub <- df %>% filter(sample == x)
    which(is.na(sub$region))
  })
  
  df <- na.omit(df)
  df <- df[,-c(1,2)]
  
  b.df <- df
  
  df <- b.df[,1:(ncol(b.df)-1)]
  
  lapply(list.of.names, function(s){
    cat(s, "---\n \n")
    st <- Sys.time()
    train <- b.df %>% filter(!sample %in% s)
    train <- train[,!colnames(train) %in% "sample"]
    var.incl <- variable.selection.f.bin(train)
    
    en <- Sys.time()
    passed <- difftime(as.POSIXct(en), as.POSIXct(st), unit="secs")
    cat("\n --- PASSED: , ", passed, " seconds \n")
    return(var.incl)
  })-> sel.var.list
  
  names(sel.var.list) <- list.of.names
  
  # save(sel.var.list, file = paste0("sel.5var.list.tvntbin.", nx, ".Rda"))
  
}

# Do this process for nx = 25 ...:
make.all.sel.var(25)

# ... and also do it for nx = 10 and 20:
lapply(c(10, 20), function(x){
  make.all.sel.var(x)
})




















# Do all predictions for all model subtypes and save them into a result list:-----
make.df.bin <- function(p.name, nx = 25){
  name1 <- paste0("https://raw.githubusercontent.com/JariCL/ROITumorDetect/refs/heads/main/nx", nx,
                  "/complete_variables_nx_", nx, "_", p.name, ".txt")
  
  read.table(file = name1)
}


make.all.pred <- function(nx = 25, substyle = 1, varcount = 5){
  
  df.big <- lapply(list.of.names, function(pattern.name){
    cat(pattern.name, " - ")

    make.df.bin(pattern.name, nx = nx)
  })
  
  
  df <- df.big
  
  df <- na.omit(df)
  df <- df[,-c(1,2)]
  
  b.df <- df
  
  df <- b.df[,1:(ncol(b.df)-1)]
  
  url1 <- paste0("https://github.com/JariCL/ROITumorDetect/blob/main/nx",nx,
                 "/sel.5var.list.tvntbin.", nx, ".Rda?raw=true")
  load(url(url1))
  
  
  tot.list.of.names <- paste0("sample_", LETTERS[1:8])
  
  which.of.tot.names <- which(tot.list.of.names %in% list.of.names)
  
  sel.var.list <- lapply(which.of.tot.names, function(x)sel.var.list[[x]])
  
  names(sel.var.list) <- list.of.names
  
  sub <- b.df
  
  
  lapply(list.of.names, function(s){
    cat(s, "---\n \n")
    train <- b.df %>% filter(!sample %in% s)
    train <- train[,!colnames(train) %in% "sample"]
    
    var.incl <-  sel.var.list[[s]][[substyle]][1:varcount]
    
    train <- cbind(train[,var.incl], train[,c("t", "nt")])%>%as.data.frame
    
    validation <- b.df %>% filter(sample %in% s)
    validation.x <- validation[,!colnames(validation) %in% c("region", "sample")]
    validation.x <- validation.x[,var.incl]
    x <- as.data.frame(train[,1:varcount])
    t <- train$t; nt <- train$nt
    mod <- glm(cbind(t,nt) ~ ., data = x, family = "binomial")
    
    predicted <- predict(mod, newdata = validation.x, type = "response")
    return(list(predicted, observed = validation$t/(validation$t+validation$nt)))
  })-> list1
  
  lapply(1:length(list1), function(x){list1[[x]][[1]]}) %>% unlist -> pred
  lapply(1:length(list1), function(x){list1[[x]][[2]]}) %>% unlist -> obs
  
  
  
  indices <- lapply(list.of.names, function(s){
    which(sub$sample == s)
  }); names(indices) <- list.of.names
  
  
  
  result.list <- list(pvo = cbind(pred, obs),
                      indices = indices)

  # save(result.list, file = paste0("result.list.tvntbin.", nx, ".", substyle, ".", varcount, ".Rda"))
  
  
  return(result.list)
}


# Do these predictions for nx = 25:
lapply(c(1,2,3), function(substyle1){
  lapply(c(2:5), function(varcount1){
    cat("\n \n ", 25, "___", substyle1, "___", varcount1, "___", "\n\n")
    make.all.pred(type = type1, nx = 25, substyle = substyle1, varcount = varcount1)
  })
})

#... and also do it for nx = 10 and 20
lapply(c(10,20), function(nx1){
  lapply(c(1,2,3), function(substyle1){
    lapply(c(2:5), function(varcount1){
      cat("\n \n ", 25, "___", substyle1, "___", varcount1, "___", "\n\n")
      make.all.pred(type = type1, nx = nx1, substyle = substyle1, varcount = varcount1)
    })
  })
})













#SAVE ROC dfs -----

# Make a data frame that for each threshold (defining the percentage of tumor cells a region needs to have to be deemed a tumor region)
# gives the amount of true positives, false positives, false negatives and true negatives:

count <- 0
do.call(rbind, lapply(c(10,20,25), function(nx1){
  # do.call(rbind, lapply(c(25), function(nx1){
  do.call(rbind, lapply(c(1,2,3), function(substyle1){
    do.call(rbind, lapply(c(2:5), function(varcount1){
      st <- Sys.time()
      count <- count + 1 
      cat("\n \n ", "Count: ", count, " of 72........ ", type1, "___", nx1, "___", substyle1, "___", varcount1, "___", "\n\n")
      url1 <- paste0("https://github.com/JariCL/ROITumorDetect/blob/main/nx",nx1,
                     "/results/result.list.tvntbin.", nx1, ".", substyle1, ".", varcount1, ".Rda?raw=true")
      load(url(url1))
      
      ROCdf <- do.call(rbind, lapply(1:1000/1000, function(threshold){
        p <- ifelse(result.list[[1]][,1]>threshold, 1,0)
        o <- ifelse(result.list[[1]][,2]>threshold, 1,0)
        TP <- sum(p==1 & o==1); TN <- sum(p==0 & o==0); FP <- sum(p==1 & o==0); FN <- sum(p==0 & o==1)
        cbind(TP = TP, FP = FP, TN = TN, FN = FN, type = type1, nx = nx1, substyle = substyle1, varcount = varcount1, threshold = threshold)
      }))
      # save(ROCdf, file = paste0("ROCdftvntbin", nx1, substyle1, varcount1, ".Rda"))
      en <- Sys.time()
      passed <- difftime(as.POSIXct(en), as.POSIXct(st), unit="secs")
      cat("\n --- PASSED: , ", passed, " seconds \n")
    }))
    
  }))
})) %>% as.data.frame 



















# Set up dataframes: -----
count <- 0
do.call(rbind, lapply(c(10,20,25), function(nx1){
  do.call(rbind, lapply(c(1,2,3), function(substyle1){
    do.call(rbind, lapply(c(2:5), function(varcount1){
      st <- Sys.time()
      count <<- count + 1 
      cat("\n \n ", "Count: ", count, " of 72........ ", type1, "___", nx1, "___", substyle1, "___", varcount1, "___", "\n\n")

      url1 <- paste0("https://github.com/JariCL/ROITumorDetect/blob/main/nx",nx1,
                     "/results/ROCdftvntbin", nx1, substyle1, varcount1, ".Rda?raw=true")
      load(url(url1))
      ROCdf <- ROCdf[1:(nrow(ROCdf)-1),]
      ROCdfc <- ROCdf
      ROCdf <- ROCdf[,1:4] %>% as.data.frame %>% apply(2,as.numeric)
      Sens <- Rec <- ROCdf[,1]/(ROCdf[,1] + ROCdf[,4])
      Spec <- ROCdf[,3]/(ROCdf[,2] + ROCdf[,3])
      Prec <- ROCdf[,1]/(ROCdf[,1] + ROCdf[,2])
      Jstat <- Sens + Spec - 1
      Fstat <- 2*(Prec*Rec)/(Prec+Rec)
      Dstat <- -sqrt((1-Sens)^2+(1-Spec)^2)
      Astat <- (ROCdf[,1] + ROCdf[,3])/(ROCdf[,1] + ROCdf[,3]+ ROCdf[,2] + ROCdf[,4])
      return(c((lapply(list(Jstat,Fstat,Dstat,Astat), function(x){
        nna <- which(!is.na(x))
        xnna <- x[!is.na(x)]
        ROCdfc[mean(nna[which(xnna == max(xnna))]),9]
      }) %>% unlist %>% as.numeric), type = type1, nx = nx1, substyle = substyle1, varcount = varcount1) )
      en <- Sys.time()
      passed <- difftime(as.POSIXct(en), as.POSIXct(st), unit="secs")
      cat("\n --- PASSED: , ", passed, " seconds \n")
    })) %>% as.data.frame %>% rename(J = 1, F = 2, D = 3, A = 4, type = 5, nx = 6, substyle = 7, varcount = 8)
  }))
})) %>% as.data.frame %>% rename(J = 1, F = 2, D = 3, A = 4) -> df.thresh


do.call(rbind, lapply(c(10,20,25), function(nx1){
  do.call(rbind, lapply(c(1,2,3), function(substyle1){
    do.call(rbind, lapply(c(2:5), function(varcount1){
      st <- Sys.time()
      count <<- count + 1 
      cat("\n \n ", "Count: ", count, " of 72........ ", type1, "___", nx1, "___", substyle1, "___", varcount1, "___", "\n\n")
      url1 <- paste0("https://github.com/JariCL/ROITumorDetect/blob/main/nx",nx1,
                     "/results/ROCdftvntbin", nx1, substyle1, varcount1, ".Rda?raw=true")
      load(url(url1))
      ROCdf <- ROCdf[1:(nrow(ROCdf)-1),]
      ROCdfc <- ROCdf
      ROCdf <- ROCdf[,1:4] %>% as.data.frame %>% apply(2,as.numeric)
      Sens <- Rec <- ROCdf[,1]/(ROCdf[,1] + ROCdf[,4])
      Spec <- ROCdf[,3]/(ROCdf[,2] + ROCdf[,3])
      Prec <- ROCdf[,1]/(ROCdf[,1] + ROCdf[,2])
      Jstat <- Sens + Spec - 1
      Fstat <- 2*(Prec*Rec)/(Prec+Rec)
      Dstat <- -sqrt((1-Sens)^2+(1-Spec)^2)
      Astat <- (ROCdf[,1] + ROCdf[,3])/(ROCdf[,1] + ROCdf[,3]+ ROCdf[,2] + ROCdf[,4])
      return(c((lapply(list(Jstat,Fstat,Dstat,Astat), function(x){
        nna <- which(!is.na(x))
        xnna <- x[!is.na(x)]
        max(xnna)
      }) %>% unlist %>% as.numeric), type = type1, nx = nx1, substyle = substyle1, varcount = varcount1) )
      en <- Sys.time()
      passed <- difftime(as.POSIXct(en), as.POSIXct(st), unit="secs")
      cat("\n --- PASSED: , ", passed, " seconds \n")
    })) %>% as.data.frame %>% rename(J = 1, F = 2, D = 3, A = 4, type = 5, nx = 6, substyle = 7, varcount = 8)
  }))
})) %>% as.data.frame %>% rename(J = 1, F = 2, D = 3, A = 4) -> df.max






# Plotting visual results observed/predicted for best model: -----

# Specifically load the best model (nx = 25, best selected first, variable count 5):
nx1 = 25; substyle1 = 1; varcount1 = 5
url1 <- paste0("https://github.com/JariCL/ROITumorDetect/blob/main/nx",nx1,
               "/results/result.list.tvntbin.", nx1, ".", substyle1, ".", varcount1, ".Rda?raw=true")

# Function that plots the observed next to the predicted:
plot.res <- function(valid){
  ind <- result.list$indices[[valid]]
  pr <- result.list[[1]][ind,1]
  ob <- result.list[[1]][ind,2]

  hull <- geojson_sf(paste0("https://raw.githubusercontent.com/JariCL/ROITumorDetect/refs/heads/main/GEOJSON/",
                            sample, ".geojson"))$geometry[[1]] %>%
    as.matrix %>% as.data.frame %>% rename(X = 1, Y = 2)
  hull <- Polygon(hull)
  hull = Polygons(list(hull),1)
  hull = SpatialPolygons(list(hull))
  hull <- st_as_sf(hull, coords = c("X", "Y"))
  tess <- st_make_grid(hull, n = 100)
  tess <- st_as_sf(tess)
  grid <- tess
  
  url1 <- paste0("https://github.com/JariCL/ROITumorDetect/blob/main/detailed%20grids/which.na.",
                 valid, ".Rda?raw=true")
  load(url(url1))

  grid <- grid[!(1:nrow(grid) %in% which.na),]
  
  grid2 <- rbind(cbind(grid,
                       value = pr, type = "Predicted"),
                 cbind(grid,
                       value = ob,  type = "Observed"))
  plot <- ggplot()+
    geom_sf(data = grid2, aes(fill = value), col = NA)+
    scale_fill_distiller(palette = "Reds", direction = 1, na.value = "black")+
    facet_wrap(~type, nrow = 1)+
    labs(fill = "")+
    theme_bw()+
    theme(text = element_text(size = 30), strip.background = element_rect(fill = "#fdd9b4"),
          axis.text = element_blank()
    )+
    scale_x_continuous(expand=c(0,0))+ 
    scale_y_continuous(expand=c(0,0))+
    guides(fill=guide_colorbar(ticks.colour = NA))
  return(plot)
}

# Run this plot for one sample to be able to extract the legend:
plot.for.legend <- plot.res(list.of.names[1])

# extract the legend:
tmp <- ggplot_gtable(ggplot_build(plot.for.legend))
leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
legend <- tmp$grobs[[leg]]

# Function that plots the observed next to the predicted, now without including the individual legends:
plot.res <- function(valid){
  ind <- result.list$indices[[valid]]
  pr <- result.list[[1]][ind,1]
  ob <- result.list[[1]][ind,2]

  hull <- geojson_sf(paste0("https://raw.githubusercontent.com/JariCL/ROITumorDetect/refs/heads/main/GEOJSON/",
                            sample, ".geojson"))$geometry[[1]] %>%
    as.matrix %>% as.data.frame %>% rename(X = 1, Y = 2)
  hull <- Polygon(hull)
  hull = Polygons(list(hull),1)
  hull = SpatialPolygons(list(hull))
  hull <- st_as_sf(hull, coords = c("X", "Y"))
  tess <- st_make_grid(hull, n = 100)
  tess <- st_as_sf(tess)
  grid <- tess
  
  url1 <- paste0("https://github.com/JariCL/ROITumorDetect/blob/main/detailed%20grids/which.na.",
                 valid, ".Rda?raw=true")
  load(url(url1))
  grid <- grid[!(1:nrow(grid) %in% which.na),]
  
  grid2 <- rbind(cbind(grid,
                       value = pr, type = "Predicted"),
                 cbind(grid,
                       value = ob,  type = "Observed"))
  plot <- ggplot()+
    geom_sf(data = grid2, aes(fill = value), col = NA)+
    scale_fill_distiller(palette = "Reds", direction = 1, na.value = "black")+
    facet_wrap(~type, nrow = 1)+
    labs(fill = "")+
    theme_bw()+
    theme(text = element_text(size = 30), strip.background = element_rect(fill = "#fdd9b4"),
          axis.text = element_blank()
    )+
    scale_x_continuous(expand=c(0,0))+ 
    scale_y_continuous(expand=c(0,0))+
    guides(fill=guide_colorbar(ticks.colour = NA))+
    theme(legend.position = "none")
  return(plot)
}

# Run this plot for each sample and then collect into one plot:
plot1 <- cowplot::plot_grid(cowplot::plot_grid(plot.res(list.of.names[1]),
                                               plot.res(list.of.names[2]),
                                               plot.res(list.of.names[3]),
                                               plot.res(list.of.names[4]), 
                                               plot.res(list.of.names[5]),
                                               plot.res(list.of.names[6]),
                                               plot.res(list.of.names[7]),
                                               plot.res(list.of.names[8]),
                                               nrow = 4, ncol = 2, label_size = 18,
                                               labels = c(letters[1:8])),
                            legend,
                            ncol = 2, rel_widths=c(10, 1))

library(grid)

jpeg(filename = "plotpredictions_all_tvnt.jpeg", width= 4200, height = 4200, res = 300)
grid.draw(plot1)
dev.off()

















# Performance/threshold plot for nx 25: -----


nx1 <- 25; substyle1 <- 1; varcount1 <- 5
url1 <- paste0("https://github.com/JariCL/ROITumorDetect/blob/main/nx",nx1,
               "/results/ROCdftvntbin", nx1, substyle1, varcount1, ".Rda?raw=true")
load(url(url1))
ROCdf <- ROCdf[1:(nrow(ROCdf)-1),]
ROCdfc <- ROCdf
ROCdf <- ROCdf[,1:4] %>% as.data.frame %>% apply(2,as.numeric)
Sens <- Rec <- ROCdf[,1]/(ROCdf[,1] + ROCdf[,4])
Spec <- ROCdf[,3]/(ROCdf[,2] + ROCdf[,3])
Prec <- ROCdf[,1]/(ROCdf[,1] + ROCdf[,2])
Jstat <- Sens + Spec - 1
Fstat <- 2*(Prec*Rec)/(Prec+Rec)
Dstat <- -sqrt((1-Sens)^2+(1-Spec)^2)
Astat <- (ROCdf[,1] + ROCdf[,3])/(ROCdf[,1] + ROCdf[,3]+ ROCdf[,2] + ROCdf[,4])


url1 <- paste0("https://github.com/JariCL/ROITumorDetect/blob/main/nx",nx1,
               "/results/ROCdftvntbin", nx1, substyle1, varcount1, ".Rda?raw=true")
load(url(url1))
thresh <- as.numeric(ROCdf[1:(nrow(ROCdf)-1),"threshold"])

ggdf <- data.frame(thresh, J = Jstat, `F` = Fstat, D = Dstat, A = Astat)




type1 = "tvnt"
url1 <- paste0("https://github.com/JariCL/ROITumorDetect/blob/main/nx",nx1,
               "/results/ROCdftvntbin", nx1, substyle1, varcount1, ".Rda?raw=true")
load(url(url1))
ROCdf <- ROCdf[1:(nrow(ROCdf)-1),]
ROCdfc <- ROCdf
ROCdf <- ROCdf[,1:4] %>% as.data.frame %>% apply(2,as.numeric)
Sens <- Rec <- ROCdf[,1]/(ROCdf[,1] + ROCdf[,4])
Spec <- ROCdf[,3]/(ROCdf[,2] + ROCdf[,3])
Prec <- ROCdf[,1]/(ROCdf[,1] + ROCdf[,2])
Jstat <- Sens + Spec - 1
Fstat <- 2*(Prec*Rec)/(Prec+Rec)
Dstat <- -sqrt((1-Sens)^2+(1-Spec)^2)
Astat <- (ROCdf[,1] + ROCdf[,3])/(ROCdf[,1] + ROCdf[,3]+ ROCdf[,2] + ROCdf[,4])
do.call(rbind, lapply(list(Jstat,Fstat,Dstat,Astat), function(x){
  nna <- which(!is.na(x))
  xnna <- x[!is.na(x)]
  c(ROCdfc[mean(nna[which(xnna == max(xnna))]),9], max(xnna))
})) %>% cbind(c("J", "F", "D", "A")) %>% as.data.frame() %>%
  rename(x=1, y = 2, name=3) %>% mutate(x=as.numeric(x), y = as.numeric(y)) -> opt.thresh




y.extdf <- data.frame(min = c(0.38,-1.05,-.03,-.03), max = c(.84,-0.22,.78,.63)) %>% mutate(name = c("A", "D", "F", "J"))

y.extdf <- y.extdf %>% arrange(factor(name, levels = c("J", "F", "D", "A")))

y.extdf <- y.extdf[,1:2]

y.ext <- data.frame(min = c(0.38,-1.05,-.03,-.03), name = c("A", "D", "F", "J"))
color.df <- data.frame(name = c("A", "D", "F", "J"), col1 = c("#fcba03",
                                                              "#fc9595",
                                                              "#0abcf2",
                                                              "#1ae888"), col2 = c("black", "black", "black", "black"))




plot.thresh1 <- function(type1 = "tvntbin", Ttype1 = "A"){
  df <- df.thresh %>%
    pivot_longer(cols = c(J,F,D,A), names_to = "Ttype", values_to = "threshold") %>%
    mutate(threshold = as.numeric(threshold)) %>%
    filter(type == type1, Ttype == Ttype1)
  
  r <- range(df$threshold)
  
  df <- df %>% mutate(indicator = as.factor(ifelse(threshold > r[1] + .7*(r[2]-r[1]),1,0 )))
  
  df <- df[,2:ncol(df)]
  
  return(ggplot(df,
                aes(x = factor(varcount), y = factor(substyle), fill = threshold)) +
           geom_point() +
           geom_tile(aes(group = nx)) +
           geom_tile(data = data.frame(varcount = as.factor(5), substyle = as.factor(1)), fill="transparent",
                     colour="black", size=1.5)+
           geom_text(aes(label = round(threshold, digits = 3), col = indicator), size = 4)+
           scale_color_manual(values = c("0" = "black", "1" = "white"))+
           labs(
             x = "",
             y = "",
             fill = "Threshold"
           ) +
           scale_fill_gradient(high = "#1944bd", low = "white", na.value = NA)+
           
           theme_bw()+
           guides(col = "none")+
           theme(panel.grid = element_blank(),
                 plot.background = element_rect(fill = "white"),
                 panel.background = element_rect(fill = "white"),
                 strip.background = element_rect(fill = color.df$col1[color.df$name == Ttype1]),
                 legend.background = element_rect(fill = "white"),
                 legend.title = element_blank(),
                 text = element_text(size = 10),
                 strip.text = element_text(size = 14, face = "bold")))
}

plot.thresh2 <- function(type1 = "tvntbin", Ttype1 = "A"){
  df <- df.max %>%
    pivot_longer(cols = c(J,F,D,A), names_to = "Ttype", values_to = "threshold") %>%
    mutate(threshold = as.numeric(threshold)) %>%
    filter(type == type1, Ttype == Ttype1)
  
  r <- range(df$threshold)
  
  df <- df %>% mutate(indicator = as.factor(ifelse(threshold > r[1] + .7*(r[2]-r[1]),1,0 )))
  
  return(ggplot(df,
                aes(x = factor(varcount), y = factor(substyle), fill = threshold)) +
           geom_point() +
           geom_tile(aes(group = nx)) +
           geom_tile(data = data.frame(varcount = as.factor(5), substyle = as.factor(1)), fill="transparent",
                     colour="white", size=1.5)+
           geom_text(aes(label = round(threshold, digits = 3), col = indicator), size = 4)+
           scale_color_manual(values = c("0" = "black", "1" = "white"))+
           
           labs(
             
             x = "",
             y = "",
             fill = "Max"
             
           ) +
           
           scale_fill_gradient(high = "#99000D", low = "white", na.value = NA)+
           theme_bw()+
           guides(col = "none")+
           theme(panel.grid = element_blank(),
                 plot.background = element_rect(fill = "white"),
                 panel.background = element_rect(fill = "white"),
                 strip.background = element_rect(fill = color.df$col1[color.df$name == Ttype1]),
                 legend.background = element_rect(fill = "white"),
                 legend.title = element_blank(),
                 text = element_text(size = 10),
                 strip.text = element_text(size = 14, face = "bold")))
}


facets <- lapply(c("A", "D", "F", "J"), function(NAME){
  check <- plot.thresh1("tvntbin", NAME)
  tmp <- ggplot_gtable(ggplot_build(check))
  leg <- which(sapply(tmp$grobs, function(x) x$name) ==  "strip")
  legend <- tmp$grobs[[leg]] 
  return(legend)
})

make.Xplot <- function(xx){
  ggdf %>% pivot_longer(cols = -thresh) %>% filter(name == xx) %>%
    ggplot()+
    geom_line(aes(x=thresh*100, y = value, col = name), alpha = .8, linewidth = 2.5, linejoin = "bevel", lineend = "butt")+
    # geom_point(data = opt.thresh[opt.thresh$name == xx,], aes(x=x*100, y = y), col = "grey30", size = 2, alpha = .8)+
    
    geom_segment(x = opt.thresh[opt.thresh$name == xx, "x"]*100,
                 xend = opt.thresh[opt.thresh$name == xx, "x"]*100,
                 y = y.ext$min[y.ext$name == xx],
                 yend = opt.thresh[opt.thresh$name == xx, "y"],
                 col = "grey30", linetype = "1111", linewidth = 1.5, alpha = .5)+
    geom_segment(x = -5,
                 xend = opt.thresh[opt.thresh$name == xx, "x"]*100,
                 y = opt.thresh[opt.thresh$name == xx, "y"],
                 yend = opt.thresh[opt.thresh$name == xx, "y"],
                 col = "grey30", linetype = "1111", linewidth = 1.5, alpha = .5)+
    geom_point(data = opt.thresh[opt.thresh$name == xx,], aes(x=x*100, y = y), fill = "grey80", col = "grey30",
               size = 3, alpha = .8, pch = 22)+
    scale_color_manual(values= c("A"= "#fcba03",
                                 "D" = "#fc9595",
                                 "F" = "#0abcf2",
                                 "J" = "#1ae888"))+
    ylab("")+
    xlab("")+
    theme_bw()+
    theme(strip.background = element_rect(fill = color.df$col1[color.df$name == xx]),
          legend.position = "none",
          text = element_text(size = 14))
}


Aplot <- make.Xplot("A")
Dplot <- make.Xplot("D")
Fplot <- make.Xplot("F")
Jplot <- make.Xplot("J")

plot.f <- plot_grid(facets[[1]],facets[[2]],facets[[3]],facets[[4]],
                    plot.thresh1("tvntbin", "A"), plot.thresh1("tvntbin", "D"),
                    plot.thresh1("tvntbin", "F"), plot.thresh1("tvntbin", "J"),
                    plot.thresh2("tvntbin", "A"), plot.thresh2("tvntbin", "D"),
                    plot.thresh2("tvntbin", "F"), plot.thresh2("tvntbin", "J"),
                    Aplot, Dplot, Fplot, Jplot, nrow = 4, ncol = 4, align = "hv", axis="tblr",
                    rel_heights = c(1,4,4,5))

plot.f <- plot.f +
  draw_label("a", x = 0.02, y = .92, size = 20, fontface = "bold") +
  draw_label("b", x = 0.02, y = .64, size = 20, fontface = "bold")+
  draw_label("c", x = 0.02, y = .35, size = 20, fontface = "bold")+
  draw_label("Variable count", x=  0.525, y=0.655, hjust= 1, angle=0)+
  draw_label("Variable count", x=  0.525, y=0.37, hjust= 1, angle=0)+
  draw_label("Cutoff (tumor %)", x=  0.525, y=0.01, hjust= 1, angle=0)+
  draw_label("Model subtype", x=  0.02, y=0.58, hjust= 1, angle=90)+
  draw_label("Model subtype", x=  0.02, y=0.85, hjust= 1, angle=90)


jpeg(filename = "dfthreshandmax.tvnt.maincheck.jpeg", width= 5600, height = 4200, res = 300)
grid.draw(plot.f)
dev.off()

















# Variable importance -----


variable.selection.f.bin <- function(df1){
  y <- df1$region
  data <- df1[,1:(ncol(df1)-1)]
  
  coef <- sapply(1:ncol(data), function(var){
    cat(var, "---")
    x <- data.matrix(data[,var])
    mod <- glm(y ~ x, family = "binomial")
    return(mod$deviance)
  })
  df <- data.frame(var = 1:ncol(data), coef = abs(coef)); df <- na.omit(df); df <- df %>% arrange(coef)
  
  return(df) 
}

make.all.sel.var <- function(nx = 20){
  
  df.big <- do.call(rbind, lapply(list.of.names, function(pattern.name){
    make.df.bin(pattern.name, nx = nx)
  }))
  
  df <- df.big
  
  which.na.per.sample <- lapply(list.of.names, function(x){
    sub <- df %>% filter(sample == x)
    which(is.na(sub$region))
  })
  
  df <- na.omit(df)
  df <- df[,-c(1,2)]
  
  b.df <- df
  
  df <- b.df[,1:(ncol(b.df)-1)]
  
  train <- df
  
  var.incl <- variable.selection.f.bin(train)
  
  
  var.incl[,1] <- colnames(train)[var.incl[,1]]
  var.ranking <- var.incl
  # save(var.ranking, file = paste0("var.ranking.tvntbin.", nx, ".Rda"))
}

lapply(c(10,20,25), function(nx1){
  make.all.sel.var(type1, nx1)
})







rankbin <- do.call(rbind, lapply(c("tvntbin"), function(type){
  do.call(rbind, lapply(c(25), function(nx){
    url1 <- paste0("https://github.com/JariCL/ROITumorDetect/blob/main/nx",nx,
                   "/var.ranking.", type, nx, ".Rda?raw=true")
    load(url(url1))
    cbind(var.ranking[,1], value = nrow(var.ranking):1, type = type, nx = nx)
  }))
})) %>% as.data.frame %>% rename(variable = 1, value = 2) %>% mutate(value= as.numeric(value))


copy <- rankbin$variable
copy <- sub("variance_red", "red_v", copy)
copy <- sub("variance_blue", "blue_v", copy)
copy <- sub("variance_green", "green_v", copy)
copy <- sub("red3", "env_red", copy)
copy <- sub("green3", "env_green", copy)
copy <- sub("blue3", "env_blue", copy)
copy <- sub("red7var", "env_red_v", copy)
copy <- sub("green7var", "env_green_v", copy)
copy <- sub("blue7var", "env_blue_v", copy)
copy <- sub("variance_green", "green_v", copy)
copy <- sub("variance_green", "green_v", copy)
copy <- sub("obs_surf", "sm", copy)
copy <- sub("obs_real", "rm", copy)
copy <- sub("obs_surf_var", "sv", copy)
copy <- sub("grad_est", "gm", copy)
copy <- sub("grad_var", "gv", copy)


rankbin$variable <- copy
un <- unique(rankbin$variable)
rank.df <- cbind(lapply(un, function(x)sum(rankbin[rankbin$variable == x,2])) %>% unlist, un) %>%
  as.data.frame %>% rename(value = 1, variable = 2) %>% mutate(value = as.numeric(value)) %>% arrange(desc(value))

keep <- rank.df$variable[1:31] 
rankbin <- do.call(rbind, lapply(c("tvntbin"), function(type){
  do.call(rbind, lapply(c(25), function(nx){
    url1 <- paste0("https://github.com/JariCL/ROITumorDetect/blob/main/nx",nx,
                   "/var.ranking.", type, nx, ".Rda?raw=true")
    load(url(url1))
    cbind(var.ranking, type = type, nx = nx)
  }))
})) %>% as.data.frame %>% rename(variable = 1, value = 2) %>% mutate(value= as.numeric(value))

copy <- rankbin$variable
copy <- sub("variance_red", "red_v", copy)
copy <- sub("variance_blue", "blue_v", copy)
copy <- sub("variance_green", "green_v", copy)
copy <- sub("red3", "env_red", copy)
copy <- sub("green3", "env_green", copy)
copy <- sub("blue3", "env_blue", copy)
copy <- sub("red7var", "env_red_v", copy)
copy <- sub("green7var", "env_green_v", copy)
copy <- sub("blue7var", "env_blue_v", copy)
copy <- sub("variance_green", "green_v", copy)
copy <- sub("variance_green", "green_v", copy)
copy <- sub("obs_surf", "sm", copy)
copy <- sub("obs_real", "rm", copy)
copy <- sub("obs_surf_var", "sv", copy)
copy <- sub("grad_est", "gm", copy)
copy <- sub("grad_var", "gv", copy)

rankbin$variable <- copy
rankbin <- rankbin %>% filter(variable %in% keep)
rankbin$variable = factor(rankbin$variable, levels = keep)

plot2 <- ggplot(rbind(rankbin)) +
  geom_tile(aes(x = variable, y = 1, fill = value), col = alpha("grey", .8))+
  scale_fill_gradient(low = "#99000D", high = "white", na.value = NA)+
  ylab("")+xlab("")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = -45, vjust = 0.5, hjust=0),
        plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white"),
        strip.background = element_rect(fill = "white"),
        legend.background = element_rect(fill = "white"),
        legend.title = element_blank(),
        axis.ticks = element_line(linewidth = .1),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        # axis.text.x = element_blank(),
        plot.margin = unit(c(1, .05, .2, .5), "cm"),
        text = element_text(size = 24))




library(cowplot)
plot.a <- plot_grid(plot2, plot1, ncol  = 1, align = "hv", axis = "tblr", labels = c("a", "b"), label_size = 30)
jpeg(filename = "varimpandrank.main.jpeg", width= 4000, height = 2800, res = 300)
plot.a
dev.off()























# Error contribution per cell type plot -----

# Reasoning: 
#Redenering hier: dus adhv de vorige test gebruiken we hier de andere, finale score: het gewicht is steeds de ratio v de counts van dat celtype tov 
# het totaal van counts VAN DAT CELTYPE (ipv het totaal van counts IN DIE GRID)!!!


nx1 <- 25
substyle1 = 1
varcount1 = 5

url1 <- paste0("https://github.com/JariCL/ROITumorDetect/blob/main/nx",nx1,
               "/results/result.list.tvntbin.", nx1, ".", substyle1, ".", varcount1, ".Rda?raw=true")
load(url(url1))
p <- result.list[[1]][,1]
o <- result.list[[1]][,2]

# Data frame that per grid (rows) shows the count for each cell type as well the predicted and observed percentage of tumor cells
cellcountstyped <- do.call(rbind, lapply(list.of.names, function(sample){
  ps <- p[result.list[[2]][[sample]]]
  os <- o[result.list[[2]][[sample]]]
  url1 <- paste0("https://github.com/JariCL/ROITumorDetect/blob/main/detailed%20grids/which.na.",
                 sample, ".Rda?raw=true")
  load(url(url1))
  # load(file =paste0("cellcounts_extrasplit_100.", sample, ".Rda"))
  return(cellcounts[(1:10000)[!(1:10000 %in% which.na)],] %>%
           cbind(O = os, P = ps) %>%
           cbind(sample = sample, type = type1, nx = nx1, substyle = substyle1, varcount = varcount1))
}))

# same data frame as above, just with per grid cell not the count of each cell type, but the
# percentage of how much they make up in the total grid cell:
cellpropstyped <- cbind(t(apply(cellcountstyped[,1:13], 1, function(x)x/sum(x))), cellcountstyped[,14:ncol(cellcountstyped)])


# Using the real data, making a metric that per non-tumor cell type makes a product of its weight in that grid cell
# as defined in 'cellpropstyped' and the difference between predicted - observed score in that grid cell 
# (essentially indicating an association metric of how strongly it was "involved" in erroneous areas):

savetrue <- lapply(list.of.names, function(sample1){
  truecountlist <- cellcountstyped[cellcountstyped$sample == sample1,c(1:6,8:13)]
  NT <- apply(truecountlist[1:ncol(truecountlist)], 1, function(x)as.numeric(x) %>% sum) %>% as.vector
  PTY <- apply(truecountlist[1:ncol(truecountlist)], 2, function(x)sum(as.numeric(x)))%>% as.vector
  D <-  cellcountstyped[cellcountstyped$sample == sample1, "P"] - cellcountstyped[cellcountstyped$sample == sample1, "O"]
  c <- do.call(cbind, lapply(c(1:12), function(x)cellcountstyped[cellcountstyped$sample == sample1,c(1:6,8:13)[x]]/PTY[x]))
  
  lapply(c(1:12), function(x){
    if (sum(c[,x]!=0) & !is.na(sum(c[,x]))){
      return(sum(c[,x]*(D)))
    } else {
      return(NA)
    }
  }) %>% unlist %>% cbind(colnames(cellcountstyped)[c(1:6,8:13)]) %>% cbind(sample = sample1, boot = "real")
})

# Per sample, calculate previously defined metric, but using permutated data 
save <- do.call(rbind, lapply(list.of.names, function(sample1){
  cat("\n", sample1, " ---- ")
  truecountlist <- cellcountstyped[cellcountstyped$sample == sample1,c(1:6,8:13)]
  NT <- apply(truecountlist[1:ncol(truecountlist)], 1, function(x)as.numeric(x) %>% sum) %>% as.vector
  PTY <- apply(truecountlist[1:ncol(truecountlist)], 2, function(x)sum(as.numeric(x)))%>% as.vector
  D <-  cellcountstyped[cellcountstyped$sample == sample1, "P"] - cellcountstyped[cellcountstyped$sample == sample1, "O"]
  c <- do.call(cbind, lapply(c(1:12), function(x)cellcountstyped[cellcountstyped$sample == sample1,c(1:6,8:13)[x]]/PTY[x]))
  
  NTcs <- c(0,cumsum(NT))
  
  lapply(1:12, function(x){
    rep(colnames(truecountlist[1:ncol(truecountlist)])[x], PTY[x])
  }) %>% unlist -> total
  
  
  lapply(1:length(NT), function(x)rep(x,NT[x])) %>% unlist -> totales
  
  st <- Sys.time()
  boot <- do.call(rbind, lapply(1:500, function(ii){
    
    cat(ii, '-')
    samp <- sample(total)
    grouped <- split(samp, totales)
    cellcounts.b <- do.call(rbind, lapply(1:length(NT), function(x){
      gcell <- grouped[[as.character(x)]]
      lapply(colnames(cellcountstyped)[c(1:6,8:13)], function(y){
        sum(gcell== y)
      }) %>% unlist
    })) %>% as.data.frame 
    
    c.b <- do.call(cbind, lapply(c(1:12), function(x)cellcounts.b[,x]/PTY[x]))
    
    res <- lapply(c(1:12), function(x){
      if (sum(c.b[,x]!=0) & !is.na(sum(c.b[,x]))){
        return(sum(c.b[,x]*(D)))
      } else {
        return(NA)
      }
    }) %>% unlist %>% cbind(colnames(cellcountstyped)[c(1:6,8:13)]) %>% cbind(sample = sample1, boot = "boot")
    
    return(res)
    
  }))
  en <- Sys.time()
  difftime(as.POSIXct(en), as.POSIXct(st), units = "secs")
  return(boot)
}))

# save(save, file = "save.under.over.error.extrasplit.tvnt.final.2.Rda")
# load(file = "save.under.over.error.extrasplit.tvnt.final.2.Rda")

data_all <- rbind(do.call(rbind, savetrue),save)%>%
  as.data.frame %>% rename(value = 1, celltype = 2, sample = 3, boot = 4) %>% mutate(value = as.numeric(value))


repres <- do.call(rbind, lapply(list.of.names, function(x){
  
  do.call(rbind, lapply(unique(data_all$celltype), function(y){
    cat(y)
    if (!is.na(data_all$value[data_all$sample == x & data_all$celltype == y & data_all$boot == "real"])){
      v <- ifelse(quantile(data_all$value[data_all$sample == x & data_all$celltype == y & data_all$boot == "boot"], .005)>
                    data_all$value[data_all$sample == x & data_all$celltype == y & data_all$boot == "real"], -1, 0)+
        ifelse(quantile(data_all$value[data_all$sample == x & data_all$celltype == y & data_all$boot == "boot"], .995)<
                 data_all$value[data_all$sample == x & data_all$celltype == y & data_all$boot == "real"], 1, 0)
      return(c(v, x, y))
    } else {
      return(c(NA, x, y))
    }
    
  }))
})) %>% as.data.frame %>% rename(value = 1, sample = 2, celltype = 3)

represv <- lapply(1:nrow(data_all %>% filter(boot == "real")), function(x){
  repres$value[repres$sample == (data_all %>% filter(boot == "real"))$sample[x] & repres$celltype == (data_all %>% filter(boot == "real"))$celltype[x]]
}) %>% unlist


indices <- do.call(rbind, lapply(1:8, function(x)cbind(ind = which(data_all$sample == list.of.names[x]), new = letters[x]))) %>%
  as.data.frame
indices[,1] <- as.numeric(indices[,1])
data_all$sample[indices[,1]] <- indices[,2]


ploterror <- ggplot(data = data_all, aes(y = value, x = sample))+
  geom_boxplot(data = data_all %>% filter(boot == "boot"), alpha = .8, col= "grey")+
  geom_point(data = data_all %>% filter(boot == "real") %>% mutate(represv), pch = 22, size = 2, aes(fill = as.factor(represv), col = as.factor(represv)))+
  
  scale_color_manual(values = c("-1" = "#1483de", "1"= "#de1414", "0" = "black"))+
  # scale_color_manual(values = c("-1" = "black", "1"= "black", "0" = "black"))+
  scale_fill_manual(values = c("-1" = "#1483de", "1"= "#de1414", "0" = "white"))+
  
  facet_wrap(vars(celltype),ncol = 2, scales = "free_y")+
  theme_bw()+
  theme(
    # axis.text.x = element_text(angle = -45, hjust = 0, vjust = .5),
    strip.background = element_rect(fill = "white"),
    legend.position = "none")

library(grid)
jpeg(filename = "error_under_over_tvnt_splitextra_final.jpeg", height=2100, width = 1400, res = 300)
print(ploterror)
dev.off()



















# Performance/thresholds plots for nx 10 and 20 (appendix): -----


allnames <- c(lapply(var_n_tr[,1], function(x)paste0(x, c("_os", "_ov"))) %>% unlist,
              lapply(var_n_tr[,1], function(x)paste0(x, c("_gs", "_gv"))) %>% unlist,
              paste0(var_n_tr[,1], "_or"))



count <- 0

type1 <- "bin"; nx1 <- 25; substyle1 <- 1; varcount1 <- 5

do.call(rbind, lapply(c("tvntbin"), function(type1){
  do.call(rbind, lapply(c(10,20), function(nx1){
    do.call(rbind, lapply(c(1,2,3), function(substyle1){
      do.call(rbind, lapply(c(2:5), function(varcount1){
        st <- Sys.time()
        count <<- count + 1 
        cat("\n \n ", "Count: ", count, " of 72........ ", type1, "___", nx1, "___", substyle1, "___", varcount1, "___", "\n\n")
        url1 <- paste0("https://github.com/JariCL/ROITumorDetect/blob/main/nx",nx1,
                       "/results/ROCdf", type1, nx1, substyle1, varcount1, ".Rda?raw=true")
        load(url(url1))
        ROCdf <- ROCdf[1:(nrow(ROCdf)-1),]
        ROCdfc <- ROCdf
        ROCdf <- ROCdf[,1:4] %>% as.data.frame %>% apply(2,as.numeric)
        Sens <- Rec <- ROCdf[,1]/(ROCdf[,1] + ROCdf[,4])
        Spec <- ROCdf[,3]/(ROCdf[,2] + ROCdf[,3])
        Prec <- ROCdf[,1]/(ROCdf[,1] + ROCdf[,2])
        Jstat <- Sens + Spec - 1
        Fstat <- 2*(Prec*Rec)/(Prec+Rec)
        Dstat <- -sqrt((1-Sens)^2+(1-Spec)^2)
        Astat <- (ROCdf[,1] + ROCdf[,3])/(ROCdf[,1] + ROCdf[,3]+ ROCdf[,2] + ROCdf[,4])
        return(c((lapply(list(Jstat,Fstat,Dstat,Astat), function(x){
          nna <- which(!is.na(x))
          xnna <- x[!is.na(x)]
          ROCdfc[mean(nna[which(xnna == max(xnna))]),9]
        }) %>% unlist %>% as.numeric), type = type1, nx = nx1, substyle = substyle1, varcount = varcount1) )
        en <- Sys.time()
        passed <- difftime(as.POSIXct(en), as.POSIXct(st), unit="secs")
        cat("\n --- PASSED: , ", passed, " seconds \n")
      })) %>% as.data.frame %>% rename(J = 1, F = 2, D = 3, A = 4, type = 5, nx = 6, substyle = 7, varcount = 8)
    }))
  }))
})) %>% as.data.frame %>% rename(J = 1, F = 2, D = 3, A = 4) -> df.thresh


do.call(rbind, lapply(c("tvntbin"), function(type1){
  do.call(rbind, lapply(c(10,20), function(nx1){
    do.call(rbind, lapply(c(1,2,3), function(substyle1){
      do.call(rbind, lapply(c(2:5), function(varcount1){
        st <- Sys.time()
        count <<- count + 1 
        cat("\n \n ", "Count: ", count, " of 72........ ", type1, "___", nx1, "___", substyle1, "___", varcount1, "___", "\n\n")
        url1 <- paste0("https://github.com/JariCL/ROITumorDetect/blob/main/nx",nx1,
                       "/results/ROCdf", type1, nx1, substyle1, varcount1, ".Rda?raw=true")
        load(url(url1))        ROCdf <- ROCdf[1:(nrow(ROCdf)-1),]
        ROCdfc <- ROCdf
        ROCdf <- ROCdf[,1:4] %>% as.data.frame %>% apply(2,as.numeric)
        Sens <- Rec <- ROCdf[,1]/(ROCdf[,1] + ROCdf[,4])
        Spec <- ROCdf[,3]/(ROCdf[,2] + ROCdf[,3])
        Prec <- ROCdf[,1]/(ROCdf[,1] + ROCdf[,2])
        Jstat <- Sens + Spec - 1
        Fstat <- 2*(Prec*Rec)/(Prec+Rec)
        Dstat <- -sqrt((1-Sens)^2+(1-Spec)^2)
        Astat <- (ROCdf[,1] + ROCdf[,3])/(ROCdf[,1] + ROCdf[,3]+ ROCdf[,2] + ROCdf[,4])
        return(c((lapply(list(Jstat,Fstat,Dstat,Astat), function(x){
          nna <- which(!is.na(x))
          xnna <- x[!is.na(x)]
          max(xnna)
        }) %>% unlist %>% as.numeric), type = type1, nx = nx1, substyle = substyle1, varcount = varcount1) )
        en <- Sys.time()
        passed <- difftime(as.POSIXct(en), as.POSIXct(st), unit="secs")
        cat("\n --- PASSED: , ", passed, " seconds \n")
      })) %>% as.data.frame %>% rename(J = 1, F = 2, D = 3, A = 4, type = 5, nx = 6, substyle = 7, varcount = 8)
    }))
  }))
})) %>% as.data.frame %>% rename(J = 1, F = 2, D = 3, A = 4) -> df.max


color.df <- data.frame(name = c("A", "D", "F", "J"), col1 = c("#fcba03",
                                                              "#fc9595",
                                                              "#0abcf2",
                                                              "#1ae888"), col2 = c("black", "black", "black", "black"))

plot.thresh1 <- function(type1 = "tvntbin", Ttype1 = "A"){
  df <- df.thresh %>%
    pivot_longer(cols = c(J,F,D,A), names_to = "Ttype", values_to = "threshold") %>%
    mutate(threshold = as.numeric(threshold)) %>%
    filter(type == type1, Ttype == Ttype1)
  
  r <- range(df$threshold)
  
  df <- df %>% mutate(indicator = as.factor(ifelse(threshold > r[1] + .7*(r[2]-r[1]),1,0 )))
  
  df <- df[,2:ncol(df)]
  
  
  df$Ttype[df$Ttype == "A"] = "Accuracy"
  
  df$Ttype[df$Ttype == "D"] = "(-)Distance"
  
  df$Ttype[df$Ttype == "F"] = "F1"
  
  df$Ttype[df$Ttype == "J"] = "Youden's J"
  
  return(ggplot(df,
                aes(x = factor(varcount), y = factor(substyle), fill = threshold)) +
           geom_point() +
           geom_tile(aes(group = nx)) +
           geom_text(aes(label = round(threshold, digits = 3), col = indicator), size = 4)+
           scale_color_manual(values = c("0" = "black", "1" = "white"))+
           facet_grid(cols = vars(Ttype), rows = vars(nx)) +
           labs(
             x = "",
             y = "",
             fill = "Threshold"
           ) +
           scale_fill_gradient(high = "#1944bd", low = "white", na.value = NA)+
           theme_bw()+
           guides(col = "none")+
           theme(panel.grid = element_blank(),
                 plot.background = element_rect(fill = "white"),
                 panel.background = element_rect(fill = "white"),
                 strip.background = element_rect(fill = color.df$col1[color.df$name == Ttype1]),
                 legend.background = element_rect(fill = "white"),
                 legend.title = element_blank(),
                 text = element_text(size = 10),
                 strip.text = element_text(size = 14, face = "bold")))
}

plot.thresh2 <- function(type1 = "tvntbin", Ttype1 = "A"){
  df <- df.max %>%
    pivot_longer(cols = c(J,F,D,A), names_to = "Ttype", values_to = "threshold") %>%
    mutate(threshold = as.numeric(threshold)) %>%
    filter(type == type1, Ttype == Ttype1)
  
  r <- range(df$threshold)
  
  df <- df %>% mutate(indicator = as.factor(ifelse(threshold > r[1] + .7*(r[2]-r[1]),1,0 )))
  
  
  df$Ttype[df$Ttype == "A"] = "Accuracy"
  
  df$Ttype[df$Ttype == "D"] = "(-)Distance"
  
  df$Ttype[df$Ttype == "F"] = "F1"
  
  df$Ttype[df$Ttype == "J"] = "Youden's J"
  
  return(ggplot(df,
                aes(x = factor(varcount), y = factor(substyle), fill = threshold)) +
           geom_point() +
           geom_tile(aes(group = nx)) +
           geom_text(aes(label = round(threshold, digits = 3), col = indicator), size = 4)+
           scale_color_manual(values = c("0" = "black", "1" = "white"))+
           facet_grid(cols = vars(Ttype), rows = vars(nx)) +
           labs(
             
             x = "",
             y = "",
             fill = "Max"
           ) +
           scale_fill_gradient(high = "#99000D", low = "white", na.value = NA)+
           
           theme_bw()+
           guides(col = "none")+
           theme(panel.grid = element_blank(),
                 plot.background = element_rect(fill = "white"),
                 panel.background = element_rect(fill = "white"),
                 strip.background = element_rect(fill = color.df$col1[color.df$name == Ttype1]),
                 legend.background = element_rect(fill = "white"),
                 legend.title = element_blank(),
                 text = element_text(size = 10),
                 strip.text = element_text(size = 14, face = "bold")))
}

library(grid)


bottom = richtext_grob("Variable count")

yleft = richtext_grob("Model subtype", rot=90)

bottom = richtext_grob("Variable count")



plot1 <- arrangeGrob(plot.thresh1("tvntbin", "A"), plot.thresh1("tvntbin", "D"),
                     plot.thresh1("tvntbin", "F"), plot.thresh1("tvntbin", "J"),
                     nrow = 2, left = yleft, bottom = bottom)

plot2 <- arrangeGrob(plot.thresh2("tvntbin", "A"), plot.thresh2("tvntbin", "D"),
                     plot.thresh2("tvntbin", "F"), plot.thresh2("tvntbin", "J"),
                     nrow = 2, left = yleft, bottom = bottom)


jpeg(filename = "dfthreshandmax.tvnt.appendix1.jpeg", width= 2800, height = 2000, res = 300)
grid.draw(plot1)
dev.off()


jpeg(filename = "dfthreshandmax.tvnt.appendix2.jpeg", width= 2800, height = 2000, res = 300)
grid.draw(plot2)
dev.off()














# Variable importance for nx 10 and 20 (appendix): -----


rankbin <- do.call(rbind, lapply(c("tvntbin"), function(type){
  do.call(rbind, lapply(c(10,20), function(nx){
    url1 <- paste0("https://github.com/JariCL/ROITumorDetect/blob/main/nx",nx,
                   "/var.ranking.", type, nx, ".Rda?raw=true")
    load(url(url1))
    cbind(var.ranking[,1], value = nrow(var.ranking):1, type = type, nx = nx)
  }))
})) %>% as.data.frame %>% rename(variable = 1, value = 2) %>% mutate(value= as.numeric(value))


copy <- rankbin$variable
copy <- sub("variance_red", "red_v", copy)
copy <- sub("variance_blue", "blue_v", copy)
copy <- sub("variance_green", "green_v", copy)
copy <- sub("red3", "env_red", copy)
copy <- sub("green3", "env_green", copy)
copy <- sub("blue3", "env_blue", copy)
copy <- sub("red7var", "env_red_v", copy)
copy <- sub("green7var", "env_green_v", copy)
copy <- sub("blue7var", "env_blue_v", copy)
copy <- sub("variance_green", "green_v", copy)
copy <- sub("variance_green", "green_v", copy)
copy <- sub("obs_surf", "sm", copy)
copy <- sub("obs_real", "rm", copy)
copy <- sub("obs_surf_var", "sv", copy)
copy <- sub("grad_est", "gm", copy)
copy <- sub("grad_var", "gv", copy)


rankbin$variable <- copy
un <- unique(rankbin$variable)

rank.df <- cbind(lapply(un, function(x)sum(rankbin[rankbin$variable == x,2])) %>% unlist, un) %>%
  as.data.frame %>% rename(value = 1, variable = 2) %>% mutate(value = as.numeric(value)) %>% arrange(desc(value))


keep <- rank.df$variable[1:31] 


rankbin <- do.call(rbind, lapply(c("tvntbin"), function(type){
  do.call(rbind, lapply(c(10,20), function(nx){
    url1 <- paste0("https://github.com/JariCL/ROITumorDetect/blob/main/nx",nx,
                   "/var.ranking.", type, nx, ".Rda?raw=true")
    load(url(url1))
    cbind(var.ranking, type = type, nx = nx)
  }))
})) %>% as.data.frame %>% rename(variable = 1, value = 2) %>% mutate(value= as.numeric(value))

copy <- rankbin$variable
copy <- sub("variance_red", "red_v", copy)
copy <- sub("variance_blue", "blue_v", copy)
copy <- sub("variance_green", "green_v", copy)
copy <- sub("red3", "env_red", copy)
copy <- sub("green3", "env_green", copy)
copy <- sub("blue3", "env_blue", copy)
copy <- sub("red7var", "env_red_v", copy)
copy <- sub("green7var", "env_green_v", copy)
copy <- sub("blue7var", "env_blue_v", copy)
copy <- sub("variance_green", "green_v", copy)
copy <- sub("variance_green", "green_v", copy)
copy <- sub("obs_surf", "sm", copy)
copy <- sub("obs_real", "rm", copy)
copy <- sub("obs_surf_var", "sv", copy)
copy <- sub("grad_est", "gm", copy)
copy <- sub("grad_var", "gv", copy)

rankbin$variable <- copy
rankbin <- rankbin %>% filter(variable %in% keep)
rankbin$variable = factor(rankbin$variable, levels = keep)


rankbin$variable = factor(rankbin$variable, levels = keep)
plot2 <- ggplot(rbind(rankbin)) +
  geom_tile(aes(x = variable, y = 1, fill = value), col = alpha("grey", .8))+
  
  scale_fill_gradient(low = "#99000D", high = "white", na.value = NA)+
  ylab("")+xlab("")+
  facet_grid(rows = vars(nx))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = -45, vjust = 0.5, hjust=0),
        plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white"),
        strip.background = element_rect(fill = "white"),
        legend.background = element_rect(fill = "white"),
        legend.title = element_blank(),
        axis.ticks = element_line(linewidth = .1),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        # axis.text.x = element_blank(),
        plot.margin = unit(c(1, .05, .2, .5), "cm"),
        text = element_text(size = 24))



library(cowplot)
plot.a <- plot_grid(plot2, plot1, ncol  = 1, align = "hv", axis = "tblr", labels = c("a", "b"), label_size = 30)

jpeg(filename = "varimpandrank.appendix.jpeg", width= 4000, height = 4000, res = 300)
plot.a
dev.off()

























# Morisita - Horn appendix plot -----

nx1 <- 25
substyle1 = 1
varcount1 = 5
standardize = T

url1 <- paste0("https://github.com/JariCL/ROITumorDetect/blob/main/nx",nx1,
               "/results/result.list.tvntbin.", nx1, ".", substyle1, ".", varcount1, ".Rda?raw=true")
load(url(url1))

valid <- sample <- list.of.names[5]
p <- result.list[[1]][result.list[[2]][[sample]],1]
o <- result.list[[1]][result.list[[2]][[sample]],2]

all.samples.MH <- do.call(rbind, lapply(list.of.names, function(sample){
  url1 <- paste0("https://github.com/JariCL/ROITumorDetect/blob/main/cellcounts/",
                 "cellcounts_extrasplit_10.", sample, ".Rda?raw=true")
  load(url(url1))
  which.na <- which(as.vector(apply(cellcounts %>% as.data.frame, 1, sum)) == 0)

  cellcounts <- cellcounts[(1:nrow(cellcounts))[!((1:nrow(cellcounts)) %in% which.na)],]
  
  # saves totals of all cell subtypes over the entire sample:
  sums <- apply(cellcounts, 2, sum)
  
  # percentage of that particular grid cell subtype count relative to totals of that cell subtype
  cc <- do.call(cbind, lapply(1:ncol(cellcounts), function(x){
    cellcounts[,x]/sums[x]
  }))
  
  
  # saves totals of per grid cells:
  PG <- cellcounts[,c(1:6,8:13)] %>% apply(1, sum) %>% as.vector
  
  # saves totals of per cell subtype over the entire sample:
  PT <- cellcounts[,c(1:6,8:13)] %>% apply(2, sum) %>% as.vector
  
  # saves entire list of named cells (later used to resample)
  all <- rep(colnames(cellcounts[,c(1:6,8:13)]), times = as.numeric(PT))
  
  # converts totals per grid cell into long string that indicates how many cell are in each grid by letting that grid cell index appear
  # that many times (e.g. if grid cell 1 has 5 cells, then '1' appears 5 times in this total string):
  times <- rep(1:length(PG), times = PG)
  
  # Resampling chunk (named cell vector gets randomly permutated) and then remake the permutated version of initial cellcount matrix
  MH.b <- do.call(rbind, lapply(1:100, function(ii){
    cat(ii, '-')
    samp <- sample(all)
    grouped <- split(samp, times)
    cellcounts.b <- do.call(rbind, lapply(1:length(PG), function(x){
      gcell <- grouped[[as.character(x)]]
      lapply(colnames(cellcounts)[c(1:6,8:13)], function(y){
        sum(gcell== y)
      }) %>% unlist
    })) %>% as.data.frame 
    colnames(cellcounts.b) <- colnames(cellcounts)[c(1:6,8:13)]
    
    # use this permutated version of the cellcount matrix to calculate per non-tumor cell type its Morisita-Horn index relating to tumor cells.
    lapply(1:12, function(x){
      X <- cellcounts.b[,x]
      Y <- cellcounts$Tumor
      XT <- sum(X)
      YT <- sum(Y)
      real <- 2*sum((X/XT)*(Y/YT))/(sum((X/XT)^2)+sum((Y/YT)^2))
    }) %>% unlist
    
  }))
  
  # now do this Morisita-Horn index relating each non-tumor cell type to tumor for the real data:
  MH <- lapply(c(1:6, 8:13), function(x){
    X <- cellcounts[,x]
    Y <- cellcounts$Tumor
    XT <- sum(X)
    YT <- sum(Y)
    real <- 2*sum((X/XT)*(Y/YT))/(sum((X/XT)^2)+sum((Y/YT)^2))
  }) %>% unlist
  
  # Save both the distribution of bootstrap/permutated MHs versus the real MH:
  cbind(rbind(cbind(MH.b, type = "b"), c(MH, type = "r")), sample = sample)
  
}))
# save(all.samples.MH, file = "all.samples.MH.10extrasplit.Rda")
# load(file = "all.samples.MH.10extrasplit.Rda")

mh <- all.samples.MH %>% as.data.frame
colnames(mh)[1:12] <- colnames(cellcounts)[c(1:6,8:13)]
mh <- mh %>% pivot_longer(cols = -c(type, sample)) %>% mutate(value = as.numeric(value))


repres <- do.call(rbind, lapply(list.of.names, function(x){
  do.call(rbind, lapply(unique(mh$name), function(y){
    if (!is.na(mh$value[mh$sample == x & mh$name == y & mh$type == "r"])){
      v <- ifelse(quantile(mh$value[mh$sample == x & mh$name == y & mh$type == "b"], .005)>mh$value[mh$sample == x & mh$name == y & mh$type == "r"], -1, 0)+
        ifelse(quantile(mh$value[mh$sample == x & mh$name == y & mh$type == "b"], .995)<mh$value[mh$sample == x & mh$name == y & mh$type == "r"], 1, 0)
      return(c(v, x, y))
    } else {
      return(c(NA, x, y))
    }
    
  }))
})) %>% as.data.frame %>% rename(value = 1, sample = 2, name = 3)

represv <- lapply(1:nrow(mh %>% filter(type == "r")), function(x){
  repres$value[repres$sample == (mh %>% filter(type == "r"))$sample[x] & repres$name == (mh %>% filter(type == "r"))$name[x]]
}) %>% unlist


indices <- do.call(rbind, lapply(1:8, function(x)cbind(ind = which(mh$sample == list.of.names[x]), new = letters[x]))) %>%
  as.data.frame
indices[,1] <- as.numeric(indices[,1])
mh$sample[indices[,1]] <- indices[,2]

plotmh <- ggplot()+
  geom_boxplot(data = mh %>% filter(type == "b"), aes(x= sample, y = value),
               alpha = .8, col= "grey")+
  geom_point(data = mh %>% filter(type == "r") %>% mutate(represv),
             pch = 22, size = 2, aes(x= sample, y = value,  fill = as.factor(represv), col = as.factor(represv)))+
  scale_color_manual(values = c("-1" = "#1483de", "1"= "#de1414", "0" = "black"))+
  scale_fill_manual(values = c("-1" = "#1483de", "1"= "#de1414", "0" = "white"))+
  facet_wrap(vars(name),ncol = 2, scales = "free_y")+
  theme_bw()+
  theme(
    strip.background = element_rect(fill = "white"),
    legend.position = "none")

library(grid)

jpeg(filename = "MH_under_over10_extrasplit.jpeg", height=2100, width = 1400, res = 300)
print(plotmh)
dev.off()










