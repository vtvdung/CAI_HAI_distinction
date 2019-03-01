library(dplyr)
library(caret)
library(lubridate)
library(naivebayes)
library(purrr)
"%ni%" <- Negate("%in%")


# (1) Training data: MALDITOF
load("baseline.rda")
baseline.df <- baseline.df[,c("Sex","Age","AdmissionDate","AdmissionTime","ward","Specimen","Specimen.Origin","InfectionReason","GDetectedDate","GDetectedTime","BacName")]
baseline.df$adm.growth <- (baseline.df$GDetectedDate - baseline.df$AdmissionDate)*24+ 
  ifelse(is.na(baseline.df$GDetectedTime - baseline.df$AdmissionTime),0,(baseline.df$GDetectedTime - baseline.df$AdmissionTime) )
baseline.df$AgeCat <- cut(baseline.df$Age, breaks=c(-1,10, 20, 40, 60,100), right = TRUE)

baseline.df$BacGroup <- as.character(baseline.df$BacName)
baseline.df$BacGroup[baseline.df$BacName %ni% 
                                 c("Escherichia coli","Streptococcus suis","Cryptococcus neoformans",
                                   "Staphylococcus aureus","Klebsiella pneumoniae","Streptococcus suis type 2",
                                   "Penicillium marneffei","Streptococcus pneumoniae","Salmonella sp.","Pseudomonas pseudomallei",
                                   "Pseudomonas aeruginosa","Acinetobacter baumannii","Pseudomonas pseudomallei")] <- "Others"
baseline.df$BacGroup[baseline.df$BacGroup=="Streptococcus suis type 2"] <- "Streptococcus suis" # correct S.suis name

baseline.df$Specimen <- ifelse(baseline.df$Specimen.Origin %ni% c("BLOOD","CSF"),"Other",baseline.df$Specimen.Origin)
baseline.df$HAICAI <- ifelse(baseline.df$adm.growth>48,"HAI","CAI")

baseline.df[,c("Sex","AgeCat","ward","Specimen","BacGroup")] <- apply(baseline.df[,c("Sex","AgeCat","ward","Specimen","BacGroup")],2,as.character)


# (2) generate models from predictors
generate_formula <- function(input_vector){
  # input_vector: a vector of all variable names
  l <- list()
  for (i in 1:length(input_vector)){
    l[[i]] <-c(TRUE,FALSE) 
  }

  l %>% expand.grid() %>% apply(1,function(x){paste("HAICAI ~ ",paste(input_vector[unlist(x)],collapse = " + "))})->all.models
  
  all.models[1:(2^length(input_vector)-1)]
  
}


# (3) prepare training and testing data
set.seed(111)
inTrain <- createDataPartition(
  y = baseline.df$HAICAI,
  ## the outcome data are needed
  p = .80,
  ## The percentage of data in the training set
  list = FALSE
)

training <- baseline.df[ inTrain,c("Sex","AgeCat","ward","Specimen","BacGroup","HAICAI")]
#training <- baseline.df[ inTrain,]
#testing  <- baseline.df[-inTrain,c("Sex","AgeCat","ward","Specimen","BacGroup")]
testing  <- baseline.df[-inTrain,c("Sex","AgeCat","ward","Specimen","BacGroup")]
#trainOutcome <- baseline.df$HAICAI[inTrain]
testOutcome <- baseline.df$HAICAI[-inTrain]
#prop.table(table(baseline.df$HAICAI))
#prop.table(table(trainOutcome))

#bootControl <- trainControl(number = 100)
list.method <- c("nb","svmRadial","lda","knn","multinom")
#lda: cannot predict

# (4) train all models

all.models.methods <- expand.grid(generate_formula(c("Sex","AgeCat","ward","Specimen","BacGroup")),list.method)
all.models.methods.for.testing <- expand.grid(generate_formula(c("Sex","AgeCat","ward","Specimen","BacGroup"))[c(1,28)],list.method)

#call.train <- function(formula,method){
#  train(as.formula(as.character(formula)),data = training,method = method)#, trControl=bootControl)
#} 

#trained.models <- map2(all.models.methods$Var1,all.models.methods$Var2,call.train)

# (5) predict and cross-validation using testing data
#trained.models %>% lapply(function(x){predict(x$finalModel,newdata=testing)})

predValues <- extractPrediction(trained.models.for.testing[1:3], testX = testing, testY = testOutcome)

nbmodel <- naive_bayes(HAICAI ~  Sex + AgeCat + ward + Specimen + BacGroup,data = training)

caretcontrol <- trainControl(number = 200)
caretmodel <- train(HAICAI ~ AgeCat,data = training, 
                    method="gbm", trcontrol=caretcontrol)

table(predict(nbmodel,newdata = testing))
table(predict(caretmodel,newdata = testing))

trainingValues <- subset(predValues,dataType=="Training")
testValues <- subset(predValues,dataType=="Test")
#by(testValues,testValues$object,function(x){confusionMatrix(x$obs,x$pred)})

confusionMatrix(testValues$obs,testValues$pred)





# (6) importance of predictors
varImp(model.ward.bac,scale = F)



# load VINARES data
load("C:/Users/vtvdu/Box Sync/VTVDUNG/PhD/Works/04 HAI CAI distinction/VINARES data/vinares.rda")
#result.groupDM <- result.groupDM[,c("SEX","AGE_GROUP","AGE","DEPARTMENT","SPEC","Bacgroup")]
result.groupDM$Sex <- ifelse(result.groupDM$SEX=="f","FEMALE",ifelse(result.groupDM$SEX=="m","MALE",NA))
result.groupDM$AgeCat <- plyr::mapvalues(result.groupDM$AGE_GROUP,
                                         from = as.character(sort(unique(result.groupDM$AGE_GROUP))),to = as.character(sort(unique(baseline.df$AgeCat))))
result.groupDM$ward <- ifelse(result.groupDM$DEPARTMENT %in% 
                                c("C?p C?","cc","CC","ccdk","cchs","cctm","HS","hscc","HSCC","hsn","HSN","hsng","icu","NICU"),
                              "Critical Care","Other")

result.groupDM$Specimen <- ifelse(result.groupDM$SPEC =='bl',"BLOOD",
                                  ifelse(result.groupDM$SPEC =='sf',"CSF","Other"))

result.groupDM$BacGroup <- plyr::mapvalues(result.groupDM$Bacgroup,
                                         from = as.character(sort(unique(result.groupDM$Bacgroup))),
                                         to = c("Acinetobacter baumannii","Escherichia coli",
                                                "Others","Others","Others","Klebsiella pneumoniae",
                                                "Others","Pseudomonas aeruginosa","Staphylococcus aureus",
                                                "Streptococcus pneumoniae"))

result.groupDM[,c("Sex","AgeCat","ward","Specimen","BacGroup")] <- apply(result.groupDM[,c("Sex","AgeCat","ward","Specimen","BacGroup")],2,as.character)

# test scenarios

