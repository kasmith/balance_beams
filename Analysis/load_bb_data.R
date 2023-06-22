## Helper functions for mdAnalysis that loads in data in consistent forms for analysis
# (Allows quick switching between different model types)

library(jsonlite)

TYPE_LEVELS = c('Bal','CB','Weight','Dist','CW','CD')
GEOMAT_TYPE_LEVELS = c('CBO','CWCB','CWCD','WB','WW')

loadShapesAgg = function(modfl) {
  
  dat = read.csv(modfl)
  
  # Assign trial types & identifiers
  dat$BaseName = sapply(as.character(dat$Trial), function(tnm) {
    spl = strsplit(tnm,'_')[[1]]
    return(paste(spl[1],spl[2],sep='_'))
  })
  dat$Type = factor(sapply(as.character(dat$Trial), function(tnm) {
    sp = strsplit(tnm,'_')
    return(sp[[1]][1])
  }), levels=TYPE_LEVELS)
  dat$Shape = sapply(as.character(dat$Trial), function(tnm) {
    sp = strsplit(tnm,'_')
    return(sp[[1]][3])
  })
  dat$Shape = with(dat, ifelse(Shape == 'CB', 'BC',Shape))
  dat$Shape = factor(dat$Shape)
  
  # Get proportion of empirical choices
  names(dat)[3:5] = c('NLeft','NBal','NRight')
  ns = with(dat, NLeft+NBal+NRight)
  dat$EmpLeft = dat$NLeft / ns
  dat$EmpBal = dat$NBal / ns
  dat$EmpRight = dat$NRight / ns
  dat$EmpAcc = with(dat, ifelse(Type %in% c('Bal','CB'), EmpBal, EmpLeft))
  dat$ModAcc = with(dat, ifelse(Type %in% c('Bal','CB'), ModBal, ModLeft))
  dat$Falls = with(dat, ifelse(Type %in% c('Bal','CB'), "B", "L"))
  
  # Reorder to look nice
  dat = dat %>% select(Trial, BaseName, Type, Shape, EmpLeft, EmpBal, EmpRight, EmpAcc, 
               ModLeft, ModBal, ModRight, ModAcc, LLH, WasFit, NLeft, NBal, NRight, Falls) %>% 
    arrange(BaseName, Shape)
  
  return(dat)
}

loadShapesInd = function(modfl, return_all = F, use_rules = F) {
  
  dat = read.csv(modfl)
  
  # Assign trial types & identifiers
  dat$BaseName = sapply(as.character(dat$Trial), function(tnm) {
    spl = strsplit(tnm,'_')[[1]]
    return(paste(spl[1],spl[2],sep='_'))
  })
  dat$Type = factor(sapply(as.character(dat$Trial), function(tnm) {
    sp = strsplit(tnm,'_')
    return(sp[[1]][1])
  }), levels=TYPE_LEVELS)
  dat$Shape = sapply(as.character(dat$Trial), function(tnm) {
    sp = strsplit(tnm,'_')
    return(sp[[1]][3])
  })
  dat$Shape = with(dat, ifelse(Shape == 'CB', 'BC',Shape))
  dat$Shape = factor(dat$Shape)
  
  if(!use_rules) {
    names(dat)[4:6] = c('EmpLeft','EmpBal','EmpRight')
  }
  
  dat$EmpAcc = with(dat, ifelse(Type %in% c('Bal','CB'), EmpBal, EmpLeft))
  dat$ModAcc = with(dat, ifelse(Type %in% c('Bal','CB'), ModBal, ModLeft))
  dat$Falls = with(dat, ifelse(Type %in% c('Bal','CB'), "B", "L"))
  
  # Reorder to look nice
  if (use_rules) {
    dat = dat %>% select(WID, Trial, BaseName, Type, Shape, EmpLeft, EmpBal, EmpRight, EmpAcc, 
                         ModLeft, ModBal, ModRight, ModAcc, LLH, Rule, Falls) %>% 
      arrange(WID, BaseName, Shape)
  } else {
    dat = dat %>% select(WID, Trial, BaseName, Type, Shape, EmpLeft, EmpBal, EmpRight, EmpAcc, 
                         ModLeft, ModBal, ModRight, ModAcc, LLH, WasFit, Falls) %>% 
      arrange(WID, BaseName, Shape)
  }
  
  if(!return_all) {
    
    if(use_rules) {
      dat = dat %>% group_by(Trial,BaseName,Type,Shape,Falls) %>%
        summarize(NLeft = sum(EmpLeft), NBal = sum(EmpBal), NRight = sum(EmpRight),
                  EmpLeft = mean(EmpLeft), EmpBal = mean(EmpBal), EmpRight = mean(EmpRight), 
                  ModLeft = mean(ModLeft), ModBal = mean(ModBal), ModRight = mean(ModRight),
                  LLH = sum(LLH)) %>%
        mutate(EmpAcc = mean(ifelse(Type %in% c('Bal','CB'), EmpBal, EmpLeft)),
               ModAcc = mean(ifelse(Type %in% c('Bal','CB'), ModBal, ModLeft))) %>%
        select(Trial, BaseName, Type, Shape, EmpLeft, EmpBal, EmpRight, EmpAcc, 
               ModLeft, ModBal, ModRight, ModAcc, LLH, NLeft, NBal, NRight) %>%
        arrange(BaseName, Shape) %>% ungroup()
    } else {
      dat = dat %>% group_by(Trial,BaseName,Type,Shape,WasFit,Falls) %>%
        summarize(NLeft = sum(EmpLeft), NBal = sum(EmpBal), NRight = sum(EmpRight),
                  EmpLeft = mean(EmpLeft), EmpBal = mean(EmpBal), EmpRight = mean(EmpRight), 
                  ModLeft = mean(ModLeft), ModBal = mean(ModBal), ModRight = mean(ModRight),
                  LLH = sum(LLH)) %>%
        mutate(EmpAcc = mean(ifelse(Type %in% c('Bal','CB'), EmpBal, EmpLeft)),
               ModAcc = mean(ifelse(Type %in% c('Bal','CB'), ModBal, ModLeft))) %>%
        select(Trial, BaseName, Type, Shape, EmpLeft, EmpBal, EmpRight, EmpAcc, 
               ModLeft, ModBal, ModRight, ModAcc, LLH, WasFit, NLeft, NBal, NRight) %>%
        arrange(BaseName, Shape) %>% ungroup()
    }
  }
  return(dat)
}

loadMatAgg = function(modfl, use_geomat = F, geomat_info_nm = '../Scene_Creation/StimuliInfo_GeomMat/GeomMatSum.csv') {
  dat = read.csv(modfl)
  
  if(use_geomat) {
    types = GEOMAT_TYPE_LEVELS
  } else {
    types = TYPE_LEVELS
  }
  
  # Assign trial types & identifiers
  if(use_geomat) {
    dat$BaseName = factor(sapply(as.character(dat$Trial), function(tnm) {
      sp = strsplit(tnm,'_')[[1]]
      return(paste(sp[1],sp[3],sep='_'))
    }))
  } else {
    dat$BaseName = dat$Trial
  }
  dat$Type = factor(sapply(as.character(dat$Trial), function(tnm) {
    sp = strsplit(tnm,'_')
    return(sp[[1]][1])
  }), levels=types)
  dat$Material = factor(sapply(as.character(dat$Trial), function(tnm) {
    sp = strsplit(tnm,'_')
    return(sp[[1]][2])
  }), levels = c('Pure','Mixed'))
  
  # Get proportion of empirical choices
  names(dat)[3:5] = c('NLeft','NBal','NRight')
  ns = with(dat, NLeft+NBal+NRight)
  dat$EmpLeft = dat$NLeft / ns
  dat$EmpBal = dat$NBal / ns
  dat$EmpRight = dat$NRight / ns
  if(use_geomat) {
    gmi = read.csv(geomat_info_nm) %>% select(Trial, Falls)
    dat = merge(dat, gmi)
  } else {
    dat$Falls = with(dat, ifelse(Type %in% c('Bal','CB'), 'B', 'L'))
  }
  dat$EmpAcc = with(dat, ifelse(Falls=='B', EmpBal, ifelse(Falls=='L', EmpLeft, EmpRight)))
  dat$ModAcc = with(dat, ifelse(Falls=='B', ModBal, ifelse(Falls=='L', ModLeft, ModRight)))
  
  # Reorder to look nice
  if(use_geomat) {
    dat = dat %>% select(Trial, BaseName, Type, Material, Falls, EmpLeft, EmpBal, EmpRight, EmpAcc, 
                         ModLeft, ModBal, ModRight, ModAcc, LLH, WasFit, NLeft, NBal, NRight) %>% 
      arrange(BaseName, Material)
  } else {
    dat = dat %>% select(Trial, BaseName, Type, Material, Falls, EmpLeft, EmpBal, EmpRight, EmpAcc, 
                         ModLeft, ModBal, ModRight, ModAcc, LLH, WasFit, NLeft, NBal, NRight) %>% 
      arrange(BaseName, Material)
  }
  
  
  return(dat)
}

loadMatInd = function(modfl, return_all = F, use_geomat = F, geomat_info_nm = '../Scene_Creation/StimuliInfo_GeomMat/GeomMatSum.csv', use_rules = F) {
  
  dat = read.csv(modfl)
  
  if(use_geomat) {
    types = GEOMAT_TYPE_LEVELS
  } else {
    types = TYPE_LEVELS
  }
  
  # Assign trial types & identifiers
  if(use_geomat) {
    dat$BaseName = factor(sapply(as.character(dat$Trial), function(tnm) {
      sp = strsplit(tnm,'_')[[1]]
      return(paste(sp[1],sp[3],sep='_'))
    }))
  } else {
    dat$BaseName = dat$Trial
  }
  dat$Type = factor(sapply(as.character(dat$Trial), function(tnm) {
    sp = strsplit(tnm,'_')
    return(sp[[1]][1])
  }), levels= types)
  dat$Material = factor(sapply(as.character(dat$Trial), function(tnm) {
    sp = strsplit(tnm,'_')
    return(sp[[1]][2])
  }), levels = c('Pure','Mixed'))
  
  if(!use_rules) {
    names(dat)[4:6] = c('EmpLeft','EmpBal','EmpRight')
  }
  
  if(use_geomat) {
    gmi = read.csv(geomat_info_nm) %>% select(Trial, Falls)
    dat = merge(dat, gmi)
  } else {
    dat$Falls = with(dat, ifelse(Type %in% c('Bal','CB'), 'B', 'L'))
  }
  dat$EmpAcc = with(dat, ifelse(Falls=='B', EmpBal, ifelse(Falls=='L', EmpLeft, EmpRight)))
  dat$ModAcc = with(dat, ifelse(Falls=='B', ModBal, ifelse(Falls=='L', ModLeft, ModRight)))
  
  # Reorder to look nice
  if (use_rules) {
    dat = dat %>% select(WID, Trial, BaseName, Type, Material, Falls, EmpLeft, EmpBal, EmpRight, EmpAcc, 
                         ModLeft, ModBal, ModRight, ModAcc, LLH, Rule) %>% 
      arrange(WID, BaseName, Material)
  } else {
    dat = dat %>% select(WID, Trial, BaseName, Type, Material, Falls, EmpLeft, EmpBal, EmpRight, EmpAcc, 
                         ModLeft, ModBal, ModRight, ModAcc, LLH, WasFit) %>% 
      arrange(WID, BaseName, Material)
  }
  
  
  
  if(!return_all) {
    
    if (use_rules) {
      dat = dat %>% group_by(Trial,BaseName,Type,Material,Falls) %>%
        summarize(NLeft = sum(EmpLeft), NBal = sum(EmpBal), NRight = sum(EmpRight),
                  EmpLeft = mean(EmpLeft), EmpBal = mean(EmpBal), EmpRight = mean(EmpRight), 
                  ModLeft = mean(ModLeft), ModBal = mean(ModBal), ModRight = mean(ModRight),
                  LLH = sum(LLH)) %>%
        mutate(EmpAcc = mean(ifelse(Type %in% c('Bal','CB'), EmpBal, EmpLeft)),
               ModAcc = mean(ifelse(Type %in% c('Bal','CB'), ModBal, ModLeft))) %>%
        select(Trial, BaseName, Type, Material, Falls, EmpLeft, EmpBal, EmpRight, EmpAcc, 
               ModLeft, ModBal, ModRight, ModAcc, LLH, NLeft, NBal, NRight) %>%
        arrange(BaseName, Material) %>% ungroup()
    } else {
      dat = dat %>% group_by(Trial,BaseName,Type,Material,Falls,WasFit) %>%
        summarize(NLeft = sum(EmpLeft), NBal = sum(EmpBal), NRight = sum(EmpRight),
                  EmpLeft = mean(EmpLeft), EmpBal = mean(EmpBal), EmpRight = mean(EmpRight), 
                  ModLeft = mean(ModLeft), ModBal = mean(ModBal), ModRight = mean(ModRight),
                  LLH = sum(LLH)) %>%
        mutate(EmpAcc = mean(ifelse(Type %in% c('Bal','CB'), EmpBal, EmpLeft)),
               ModAcc = mean(ifelse(Type %in% c('Bal','CB'), ModBal, ModLeft))) %>%
        select(Trial, BaseName, Type, Material, Falls, EmpLeft, EmpBal, EmpRight, EmpAcc, 
               ModLeft, ModBal, ModRight, ModAcc, LLH, WasFit, NLeft, NBal, NRight) %>%
        arrange(BaseName, Material) %>% ungroup()
    }
  }
  return(dat)
}



loadBalAgg = function(modfl, trinfonm = "../Scene_Creation/StimuliInfo_Balance/BalMaterialSum.csv") {
  
  dat = read.csv(modfl)
  trinfo = read.csv(trinfonm)
  
  trialinfo = lapply(as.character(dat$Trial), function(tnm) {
    spl = strsplit(tnm,'_')[[1]]
    if (spl[1] == 'BalAdjust') {
      return(c(spl[1:3],spl[5],spl[4],'NA'))
    } else{
      return(c('SizeAdjust',spl[1:2],spl[4],spl[5],spl[3]))
    }
  })
  
  # Assign trial types & identifiers
  dat$BaseName = sapply(trialinfo, function(ti) {
    if (ti[1] == 'BalAdjust') {
      return(paste(ti[1],ti[3],ti[4],ti[5],sep='_'))
    } else {
      return(paste(ti[1],ti[2],ti[3],ti[4],ti[6],sep='_'))
    }
  })
  dat$Type = factor(sapply(trialinfo, function(ti) {
    return(ifelse(ti[3] %in% c('CD','CW'),ti[3],
           ifelse(ti[3]=='D','Dist','Weight')))
  }), levels=TYPE_LEVELS)
  dat$Class = sapply(trialinfo, function(ti) {
    return(ti[1])
  })
  dat$Centering = sapply(trialinfo, function(ti) {
    return(ifelse(ti[1]=='BalAdjust',ifelse(ti[2]=='Raw','Centered','Uncentered'),ti[2]))
  })
  dat$StrutWidth = sapply(trialinfo, function(ti) {
    return(factor(ti[5], levels = c('0.25', '0.5', '1.0', '2.0')))
  })
  dat$SubClass = with(dat, interaction(StrutWidth,Centering,sep='_'))
  
  # Get proportion of empirical choices
  names(dat)[3:5] = c('NLeft','NBal','NRight')
  ns = with(dat, NLeft+NBal+NRight)
  dat$EmpLeft = dat$NLeft / ns
  dat$EmpBal = dat$NBal / ns
  dat$EmpRight = dat$NRight / ns
  
  dat = merge(dat, trinfo[c('Trial','Falls')])
  
  dat$EmpAcc = with(dat, ifelse(Falls=='R',EmpRight,ifelse(Falls=='L',EmpLeft,EmpBal)))
  dat$ModAcc = with(dat, ifelse(Falls=='R',ModRight,ifelse(Falls=='L',ModLeft,ModBal)))
  
  # Reorder to look nice
  dat = dat %>% select(Trial, BaseName, Type, Class, Centering, StrutWidth, SubClass,
                       Falls, EmpLeft, EmpBal, EmpRight, EmpAcc, 
                       ModLeft, ModBal, ModRight, ModAcc, LLH, WasFit, NLeft, NBal, NRight) %>% 
    arrange(BaseName, Class, SubClass)
  
  return(dat)
}

loadBalInd = function(modfl, return_all = F, trinfonm = "../Scene_Creation/StimuliInfo_Balance/BalMaterialSum.csv", use_rules = F) {
  
  dat = read.csv(modfl)
  trinfo = read.csv(trinfonm)
  
  trialinfo = lapply(as.character(dat$Trial), function(tnm) {
    spl = strsplit(tnm,'_')[[1]]
    if (spl[1] == 'BalAdjust') {
      return(c(spl[1:3],spl[5],spl[4],'NA'))
    } else{
      return(c('SizeAdjust',spl[1:2],spl[4],spl[5],spl[3]))
    }
  })
  
  # Assign trial types & identifiers
  dat$BaseName = sapply(trialinfo, function(ti) {
    if (ti[1] == 'BalAdjust') {
      return(paste(ti[1],ti[3],ti[4],ti[5],sep='_'))
    } else {
      return(paste(ti[1],ti[2],ti[3],ti[4],ti[6],sep='_'))
    }
  })
  dat$Type = factor(sapply(trialinfo, function(ti) {
    return(ifelse(ti[3] %in% c('CD','CW'),ti[3],
                  ifelse(ti[3]=='D','Dist','Weight')))
  }), levels=TYPE_LEVELS)
  dat$Class = sapply(trialinfo, function(ti) {
    return(ti[1])
  })
  dat$Centering = sapply(trialinfo, function(ti) {
    return(ifelse(ti[1]=='BalAdjust',ifelse(ti[2]=='Raw','Centered','Uncentered'),ti[2]))
  })
  dat$StrutWidth = sapply(trialinfo, function(ti) {
    return(factor(ti[5], levels = c('0.25', '0.5', '1.0', '2.0')))
  })
  dat$SubClass = with(dat, interaction(StrutWidth,Centering,sep='_'))
  
  if(!use_rules) {
    names(dat)[4:6] = c('EmpLeft','EmpBal','EmpRight')
  }
  dat = merge(dat, trinfo[c('Trial','Falls')])
  
  dat$EmpAcc = with(dat, ifelse(Falls=='R',EmpRight,ifelse(Falls=='L',EmpLeft,EmpBal)))
  dat$ModAcc = with(dat, ifelse(Falls=='R',ModRight,ifelse(Falls=='L',ModLeft,ModBal)))
  
  # Reorder to look nice
  if (use_rules) {
    dat = dat %>% select(WID, Trial, BaseName, Type, Class, Centering, StrutWidth, SubClass,
                         Falls, EmpLeft, EmpBal, EmpRight, EmpAcc, 
                         ModLeft, ModBal, ModRight, ModAcc, LLH, Rule) %>% 
      arrange(WID, BaseName, Class, SubClass)
  } else {
    dat = dat %>% select(WID, Trial, BaseName, Type, Class, Centering, StrutWidth, SubClass,
                         Falls, EmpLeft, EmpBal, EmpRight, EmpAcc, 
                         ModLeft, ModBal, ModRight, ModAcc, LLH, WasFit) %>% 
      arrange(WID, BaseName, Class, SubClass)
  }
  
  if(!return_all) {
    if (use_rules) {
      dat = dat %>% group_by(Trial,BaseName,Type,Class, Centering,StrutWidth,SubClass, Falls,WasFit) %>%
        summarize(NLeft = sum(EmpLeft), NBal = sum(EmpBal), NRight = sum(EmpRight),
                  EmpLeft = mean(EmpLeft), EmpBal = mean(EmpBal), EmpRight = mean(EmpRight), 
                  ModLeft = mean(ModLeft), ModBal = mean(ModBal), ModRight = mean(ModRight),
                  LLH = sum(LLH)) %>%
        mutate(EmpAcc = mean(ifelse(Falls=='R',EmpRight,ifelse(Falls=='L',EmpLeft,EmpBal))),
               ModAcc = mean(ifelse(Falls=='R',ModRight,ifelse(Falls=='L',ModLeft,ModBal)))) %>%
        select(Trial, BaseName, Type, Class, Centering, StrutWidth, SubClass, Falls, EmpLeft, EmpBal, EmpRight, EmpAcc, 
               ModLeft, ModBal, ModRight, ModAcc, LLH, NLeft, NBal, NRight) %>%
        arrange(BaseName, Class, SubClass) %>% ungroup()
    } else {
      dat = dat %>% group_by(Trial,BaseName,Type,Class, Centering,StrutWidth,SubClass, Falls,WasFit) %>%
        summarize(NLeft = sum(EmpLeft), NBal = sum(EmpBal), NRight = sum(EmpRight),
                  EmpLeft = mean(EmpLeft), EmpBal = mean(EmpBal), EmpRight = mean(EmpRight), 
                  ModLeft = mean(ModLeft), ModBal = mean(ModBal), ModRight = mean(ModRight),
                  LLH = sum(LLH)) %>%
        mutate(EmpAcc = mean(ifelse(Falls=='R',EmpRight,ifelse(Falls=='L',EmpLeft,EmpBal))),
               ModAcc = mean(ifelse(Falls=='R',ModRight,ifelse(Falls=='L',ModLeft,ModBal)))) %>%
        select(Trial, BaseName, Type, Class, Centering, StrutWidth, SubClass, Falls, EmpLeft, EmpBal, EmpRight, EmpAcc, 
               ModLeft, ModBal, ModRight, ModAcc, LLH, WasFit, NLeft, NBal, NRight) %>%
        arrange(BaseName, Class, SubClass) %>% ungroup()
    }
  }
  return(dat)
}

loadCombAgg = function(modfl) {
  
  dat = read.csv(modfl)
  
  trialinfo = lapply(as.character(dat$Trial), function(tnm) {
    return(strsplit(tnm,'_')[[1]])
  })
  
  # Assign trial types & identifiers
  dat$BaseName = sapply(trialinfo, function(ti) {
    return(paste(ti[1],ti[2],ti[3],ti[4],ti[5],sep='_'))
  })
  dat$Type = factor(sapply(trialinfo, function(ti) {
    return(ifelse(ti[1] %in% c('CD','CW','CB'),ti[1],
                  ifelse(ti[1]=='D','Dist',ifelse(ti[1]=='W','Weight','Bal'))))
  }), levels=TYPE_LEVELS)
  dat$Materials = sapply(trialinfo, function(ti) {
    return(ti[2])
  })
  dat$Centering = sapply(trialinfo, function(ti) {
    return(ti[3])
  })
  dat$StrutWidth = sapply(trialinfo, function(ti) {
    return(factor(ti[4], levels = c('0.25', '0.5', '1.0', '2.0')))
  })
  dat$ShapeType = sapply(trialinfo, function(ti) {
    return(ti[6])
  })
  
  # Get proportion of empirical choices
  names(dat)[3:5] = c('NLeft','NBal','NRight')
  ns = with(dat, NLeft+NBal+NRight)
  dat$EmpLeft = dat$NLeft / ns
  dat$EmpBal = dat$NBal / ns
  dat$EmpRight = dat$NRight / ns
  
  dat$Falls = with(dat, ifelse(Type %in% c('Bal','CB'), 'B', 'L'))
  
  dat$EmpAcc = with(dat, ifelse(Falls=='R',EmpRight,ifelse(Falls=='L',EmpLeft,EmpBal)))
  dat$ModAcc = with(dat, ifelse(Falls=='R',ModRight,ifelse(Falls=='L',ModLeft,ModBal)))
  
  dat$NConds = factor(with(dat, (Materials == 'mixed') + (Centering == 'uncentered') + (ShapeType == 'mixed') + (StrutWidth == "1.0")))
  
  # Reorder to look nice
  dat = dat %>% select(Trial, BaseName, Type, Materials, Centering, StrutWidth, ShapeType, NConds,
                       Falls, EmpLeft, EmpBal, EmpRight, EmpAcc, 
                       ModLeft, ModBal, ModRight, ModAcc, LLH, WasFit, NLeft, NBal, NRight) %>% 
    arrange(BaseName, Materials, Centering, StrutWidth, ShapeType)
  
  return(dat)
}

loadCombInd = function(modfl, return_all = F) {
  
  dat = read.csv(modfl)
  
  trialinfo = lapply(as.character(dat$Trial), function(tnm) {
    return(strsplit(tnm,'_')[[1]])
  })
  
  # Assign trial types & identifiers
  dat$BaseName = sapply(trialinfo, function(ti) {
    return(paste(ti[1],ti[2],ti[3],ti[4],ti[5],sep='_'))
  })
  dat$Type = factor(sapply(trialinfo, function(ti) {
    return(ifelse(ti[1] %in% c('CD','CW','CB'),ti[1],
                  ifelse(ti[1]=='D','Dist',ifelse(ti[1]=='W','Weight','Bal'))))
  }), levels=TYPE_LEVELS)
  dat$Materials = sapply(trialinfo, function(ti) {
    return(ti[2])
  })
  dat$Centering = sapply(trialinfo, function(ti) {
    return(ti[3])
  })
  dat$StrutWidth = sapply(trialinfo, function(ti) {
    return(factor(ti[4], levels = c('0.25', '0.5', '1.0', '2.0')))
  })
  dat$ShapeType = sapply(trialinfo, function(ti) {
    return(ti[6])
  })
  
  # Get proportion of empirical choices
  names(dat)[4:6] = c('NLeft','NBal','NRight')
  ns = with(dat, NLeft+NBal+NRight)
  dat$EmpLeft = dat$NLeft / ns
  dat$EmpBal = dat$NBal / ns
  dat$EmpRight = dat$NRight / ns
  
  dat$Falls = with(dat, ifelse(Type %in% c('Bal','CB'), 'B', 'L'))
  
  dat$EmpAcc = with(dat, ifelse(Falls=='R',EmpRight,ifelse(Falls=='L',EmpLeft,EmpBal)))
  dat$ModAcc = with(dat, ifelse(Falls=='R',ModRight,ifelse(Falls=='L',ModLeft,ModBal)))
  
  dat$NConds = factor(with(dat, (Materials == 'mixed') + (Centering == 'uncentered') + (ShapeType == 'mixed') + (StrutWidth == "1.0")))
  
  # Reorder to look nice
  dat = dat %>% select(WID, Trial, BaseName, Type, Materials, Centering, StrutWidth, ShapeType, NConds, 
                       Falls, EmpLeft, EmpBal, EmpRight, EmpAcc, 
                       ModLeft, ModBal, ModRight, ModAcc, LLH, WasFit, NLeft, NBal, NRight) %>% 
    arrange(WID, BaseName, Materials, Centering, StrutWidth, ShapeType)
  
  
  
  if(!return_all) {
    dat = dat %>% group_by(Trial,BaseName,Type,Materials, Centering,StrutWidth,ShapeType,NConds, Falls, WasFit) %>%
      summarize(NLeft = sum(EmpLeft), NBal = sum(EmpBal), NRight = sum(EmpRight),
                EmpLeft = mean(EmpLeft), EmpBal = mean(EmpBal), EmpRight = mean(EmpRight), 
                ModLeft = mean(ModLeft), ModBal = mean(ModBal), ModRight = mean(ModRight),
                LLH = sum(LLH)) %>%
      mutate(EmpAcc = mean(ifelse(Falls=='R',EmpRight,ifelse(Falls=='L',EmpLeft,EmpBal))),
             ModAcc = mean(ifelse(Falls=='R',ModRight,ifelse(Falls=='L',ModLeft,ModBal)))) %>%
      select(Trial, BaseName, Type, Materials, Centering,StrutWidth,ShapeType,NConds, Falls, EmpLeft, EmpBal, EmpRight, EmpAcc, 
             ModLeft, ModBal, ModRight, ModAcc, LLH, WasFit, NLeft, NBal, NRight) %>%
      arrange(BaseName, Class, SubClass) %>% ungroup()
  }
  
  
  return(dat)
}


kludge_cv_modtype = function(prefix = '../output/complete/cv_modtype_', nparts = 3, nsims = 1) {
  
  maxn = 0
  
  for (i in 1:nparts) {
    if (i == 1) {
      dat = read.csv(paste(prefix,nsims,'sim/cv_llh_grid.csv',sep=''))
      maxn = max(dat$Order)
    } else {
      newdat = read.csv(paste(prefix,nsims,'sim_pt',i,'/cv_llh_grid.csv',sep=''))
      newdat$Order = newdat$Order + maxn
      dat = rbind(dat, newdat)
      maxn = max(dat$Order)
    }
  }
  return(dat)
}

loadFerrettiAgg = function(modfl) {
  dat = read.csv(modfl)
  types = TYPE_LEVELS
  
  dat$BaseName = dat$Trial
  dat$Type = factor(sapply(as.character(dat$Trial), function(tnm) {
    sp = strsplit(tnm,'_')
    return(sp[[1]][1])
  }), levels=types)
  dat$DiffType = factor(sapply(as.character(dat$Trial), function(tnm) {
    sp = strsplit(tnm,'_')
    ttype = sp[[1]][1]
    if (ttype %in% c('Bal','CB')) {
      return(0)
    } else {
      return (sp[[1]][2])
    }
  }))
  
  # Get proportion of empirical choices
  names(dat)[3:5] = c('NLeft','NBal','NRight')
  ns = with(dat, NLeft+NBal+NRight)
  dat$EmpLeft = dat$NLeft / ns
  dat$EmpBal = dat$NBal / ns
  dat$EmpRight = dat$NRight / ns
  dat$Falls = with(dat, ifelse(Type %in% c('Bal','CB'), 'B', 'L'))

  dat$EmpAcc = with(dat, ifelse(Falls=='B', EmpBal, ifelse(Falls=='L', EmpLeft, EmpRight)))
  dat$ModAcc = with(dat, ifelse(Falls=='B', ModBal, ifelse(Falls=='L', ModLeft, ModRight)))
  
  # Reorder to look nice
  dat = dat %>% select(Trial, BaseName, Type, DiffType, Falls, EmpLeft, EmpBal, EmpRight, EmpAcc, 
                       ModLeft, ModBal, ModRight, ModAcc, LLH, WasFit, NLeft, NBal, NRight) %>% 
    arrange(BaseName, DiffType)
  
  return(dat)
}


read_json_strat_params = function(jsonfile) {
  dat = read_json(jsonfile)
  inds = dat$ind_strategies
  nwid = length(inds)
  ret = NULL
  for (i in 1:nwid) {
    wid = names(inds)[i]
    idat = inds[[i]]
    ret = rbind(ret,
                data.frame(WID=wid, Model=idat$model, SP=idat$strat_params$sp,
                           SWP=idat$strat_params$smp, Guess=idat$strat_params$guess))
  }
  levels(ret$Model) = c("Shapes", "Pivot", "Materials")
  return(ret)
}
