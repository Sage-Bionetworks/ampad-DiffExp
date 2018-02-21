### R function to combine all differential expression results

## Load libraries
library(synapseClient)
library(knitr)
library(githubr)

library(CovariateAnalysis)
library(data.table)
library(plyr)
library(tidyverse)

synapseLogin()

diff.exp = lapply(c(MAYO = 'syn8468023',
                    MSSM = 'syn10526259',
                    ROSMAP = 'syn8456721'),
                  function(id){
                    fread(synGet(id)@filePath, data.table = FALSE, header = T)
                  })

############################################
## Combine all differential expression results
#### Fix format issues with MSSM data
diff.exp$MSSM = diff.exp$MSSM %>%
  plyr::ddply(.(Model), .fun = function(x){
    if (any(x$Model %in% c('Diagnosis.SEX', 'Diagnosis.AOD'))){
      x %>%
        dplyr::select(-Tissue) %>%
        tidyr::separate(Comparison, c('c1','c2'), sep = '\\-') %>%
        tidyr::separate(c1, c('Tissue', 'c11','SEX'), sep = '\\.') %>%
        tidyr::separate(c2, c('Tissue1', 'c21','SEX1'), sep = '\\.') %>%
        tidyr::unite(Comparison, c11, c21, sep = '-')
    } else {
     return(x) 
    }
  })
diff.exp$MSSM$SEX[diff.exp$MSSM$SEX == 'AOD'] = 'ALL'
diff.exp$MSSM$SEX[is.na(diff.exp$MSSM$SEX)] = 'ALL'

#### Organise model names
diff.exp$MSSM$Model = factor(diff.exp$MSSM$Model)
levels(diff.exp$MSSM$Model) = c("APOE4", "CDR", "Diagnosis", "Diagnosis.AOD", "Diagnosis.Sex")

diff.exp$MAYO$Model = factor(diff.exp$MAYO$Model)
levels(diff.exp$MAYO$Model) = c("APOE4", "Diagnosis", "Diagnosis.AOD", "Diagnosis.Sex", "SourceDiagnosis")

diff.exp$ROSMAP$Model = factor(diff.exp$ROSMAP$Model)
levels(diff.exp$ROSMAP$Model) = c("APOE4", "cogdx", "Diagnosis", "Diagnosis.AOD", "Diagnosis.Sex")

#### Organise tissue names
diff.exp$MSSM$Tissue = factor(diff.exp$MSSM$Tissue)
levels(diff.exp$MSSM$Tissue) = c("FP", "FP", "IFG", "IFG", "PHG", "PHG", "STG", "STG")

setnames(diff.exp$MAYO, 'Tissue.ref', 'Tissue')
diff.exp$MAYO$Tissue.ag = NULL
diff.exp$MAYO$Tissue = factor(diff.exp$MAYO$Tissue)

setnames(diff.exp$ROSMAP, 'Region', 'Tissue')
diff.exp$ROSMAP$Tissue = factor(diff.exp$ROSMAP$Tissue)

#### Organise Sex and fix comparison
diff.exp$ROSMAP$Sex = 'ALL'
diff.exp$ROSMAP$Sex[diff.exp$ROSMAP$Comparison == "AD-CONTROL.IN.FEMALE"] = "FEMALE"
diff.exp$ROSMAP$Sex[diff.exp$ROSMAP$Comparison == "AD-CONTROL.IN.MALE"] = "MALE"

diff.exp$ROSMAP$Comparison = factor(diff.exp$ROSMAP$Comparison)
levels(diff.exp$ROSMAP$Comparison) = c("AD-CONTROL", "AD-CONTROL", "AD-CONTROL.IN.FEMALE-MALE", "AD-CONTROL",
                                       "AOD", "AD-CONTROL", "1-0", "2-0", "2-1", "2-1", "4-1", "4-2", 
                                       "FEMALE-MALE")

setnames(diff.exp$MSSM, 'SEX','Sex')

diff.exp$MSSM$Comparison = factor(diff.exp$MSSM$Comparison)
levels(diff.exp$MSSM$Comparison) = c("0-1", "0-2", "1-2", "AD-CONTROL", "AD-OTHER", "ALL-ALL", "ALL-NA", "CDR", "OTHER-CONTROL")     

setnames(diff.exp$MAYO, 'Sex.ag', 'Sex')
diff.exp$MAYO$Sex[is.na(diff.exp$MAYO$Sex)] = 'ALL'

diff.exp$MAYO$Comparison = factor(diff.exp$MAYO$Comparison)
levels(diff.exp$MAYO$Comparison) = c("AD-CONTROL", "AD-OTHER", "AD-PATH_AGE", "1-0",  
                                     "2-0", "2-1", "OTHER-CONTROL", "PATH_AGE-CONTROL", "PSP-CONTROL")   

# Fix study name
diff.exp$ROSMAP$Study = NULL

colNames = sapply(diff.exp, colnames) %>% Reduce(intersect,.)
diff.exp = data.table::rbindlist(diff.exp, use.names = T, fill = T, idcol = 'Study') %>%
  data.frame() %>%
  dplyr::select(one_of(colNames), Study)

# Assign direction
diff.exp$Direction = 'NONE'
diff.exp$Direction[diff.exp$adj.P.Val <= 0.05 & diff.exp$logFC >= log2(1.2)] = 'UP'
diff.exp$Direction[diff.exp$adj.P.Val <= 0.05 & diff.exp$logFC <= -log2(1.2)] = 'DOWN'

# Get gene sets
amp.ad.de.geneSets = diff.exp %>%
  dplyr::filter(Direction != 'NONE') %>%
  plyr::dlply(.(Model, Tissue,Comparison, Sex, Direction), .fun = function(x){
    unique(x$ensembl_gene_id)
  })
  
# Get github commit
thisFileName <- 'collateDiffExp.R'
thisRepo <- getRepo(repository = "th1vairam/ampad-DiffExp", ref="branch", refName='geneLevelAnalysis')
thisFile <- getPermlink(repository = thisRepo, repositoryPath=paste0('utilityFunctions/',thisFileName))

# Store results in synapse
save(list='amp.ad.de.geneSets', file = 'all.diff.exp.gs.RData')
obj = File('all.diff.exp.gs.RData', name = 'All Differential Expression GeneSets (RData format)', parentId = 'syn8672415')
obj = synStore(obj, used = c('syn8468023', 'syn10157628', 'syn8456721'), 
               executed = thisFile, activityName = 'Collate differential expression results')

# Merged differential expression
write.table(diff.exp, 'differentialExpressionSummary.tsv', row.names = F, quote = F, sep = '\t')
obj = File('differentialExpressionSummary.tsv', 
           name = 'All Differential Expression (Merged)', 
           parentId = 'syn8672415')
obj = synStore(obj, used = c('syn8468023', 'syn10526259', 'syn8456721'), 
               executed = thisFile, activityName = 'Collate differential expression results')