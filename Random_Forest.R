
data <- t(voomed_counts) #colnames should be genes and rownames should be samples; TRANSPOSED
#colnames(data) <- ev_markers$gene_symbol[match(ev_markers$entrez_id,colnames(data))]
data <- data.frame(data)
#add label
data$phenotype <- phenotype$condition[match(rownames(data), phenotype$sample)] #make sure the order is right
set.seed(100)
#set proportion
i <- sample(nrow(data), 0.8*nrow(data), replace = FALSE)
train <- data[i,]
table(train$phenotype)
#         Metastatic       Primary Tumor Solid Tissue Normal 
#                  1                 396                  43 
test <- data[-i,]
table(test$phenotype)
#      Primary Tumor Solid Tissue Normal 
#                101                   9 

#control
control <- trainControl(method="repeatedcv", number=10, repeats=10)
#training
set.seed(100)
#mtry <- sqrt(ncol(train)) #number of variables randomly sampled as attribute at each split
#tunegrid <- expand.grid(.mtry=mtry)
rf_untuned <- train(phenotype~., 
	data=train, 
	method="rf", 
	metric="Accuracy", 
	trControl=control)
print(rf_untuned)
#Random Forest 

#100 samples
# 15 predictor
#  2 classes: 'Primary Tumor', 'Solid Tissue Normal' 

#No pre-processing
#Resampling: Cross-Validated (10 fold, repeated 10 times) 
#Summary of sample sizes: 395, 397, 396, 396, 395, 395, ... 
#Resampling results across tuning parameters:

#  mtry  Accuracy  Kappa
#   2    0.832     0.664
#   8    0.839     0.678
#  15    0.827     0.654

#Accuracy was used to select the optimal model using the largest value.
#The final value used for the model was mtry = 7.

#identify best mtry
set.seed(100)
tuneGrid <- expand.grid(.mtry = c(1: 10))
rf_mtry <- train(phenotype~.,
    data = train,
    method = "rf",
    metric = "Accuracy",
    tuneGrid = tuneGrid,
    trControl = control,
    importance = TRUE,
    nodesize = 14,
    ntree = 300)
print(rf_mtry)
max(rf_mtry$results$Accuracy)
#store best mtry
best_mtry <- rf_mtry$bestTune$mtry

####identify best maxnodes; tree depth
#kappa = compares observed accuracy with expected accuracy; 
	#not only evaluate a single classifier but amongst classifier themselves
store_maxnode <- list()
tuneGrid <- expand.grid(.mtry = best_mtry)
for (maxnodes in c(3: 25)) {
    set.seed(100)
    rf_maxnode <- train(phenotype~.,
        data = train,
        method = "rf",
        metric = "Accuracy",
        tuneGrid = tuneGrid,
        trControl = control,
        importance = TRUE,
        nodesize = 3,
        maxnodes = maxnodes,
        ntree = 300)
    current_iteration <- toString(maxnodes)
    store_maxnode[[current_iteration]] <- rf_maxnode
}
results_mtry <- resamples(store_maxnode)
summary(results_mtry)
max_node <- 5

#identify best num of trees
store_maxtrees <- list()
for (ntree in c(20, 30, 50, 100, 200, 250, 300, 350, 400, 450, 500, 550, 600, 800, 1000, 2000)) {
    set.seed(100)
    rf_maxtrees <- train(phenotype~.,
        data = train,
        method = "rf",
        metric = "Accuracy",
        tuneGrid = tuneGrid,
        trControl = control,
        importance = TRUE,
        nodesize = 3,
        maxnodes = max_node,
        ntree = ntree)
    key <- toString(ntree)
    store_maxtrees[[key]] <- rf_maxtrees
}
results_tree <- resamples(store_maxtrees)
summary(results_tree)
best_ntree <- 250

fit_rf <- train(phenotype~.,
    data = train,
    method = "rf",
    metric = "Accuracy",
    tuneGrid = tuneGrid,
    trControl = control,
    importance = TRUE,
    nodesize = 3,
    ntree = best_ntree,
    maxnodes = max_node)

#evaluate model
prediction <- predict(fit_rf, test)
confusionMatrix(prediction, as.factor(test$phenotype))
#Confusion Matrix and Statistics
#                     Reference
#Prediction            Primary Tumor Solid Tissue Normal
#  Primary Tumor                 101                   5
#  Solid Tissue Normal             2                   2                                        
#               Accuracy : 0.9364         
#                 95% CI : (0.8733, 0.974)
#    No Information Rate : 0.9364         
#    P-Value [Acc > NIR] : 0.5988                                                 
#                  Kappa : 0.3328                                                 
# Mcnemar's Test P-Value : 0.4497                                                 
#            Sensitivity : 0.9806         
#            Specificity : 0.2857         
#         Pos Pred Value : 0.9528         
#         Neg Pred Value : 0.5000         
#             Prevalence : 0.9364         
#         Detection Rate : 0.9182         
#   Detection Prevalence : 0.9636         
#      Balanced Accuracy : 0.6331         
#                                         
#       'Positive' Class : Primary Tumor 


varImp(fit_rf)
#rf variable importance
#        Importance
#NKX3.1     100.000
#SMYD3       99.733
#TET3        97.114
#CYLD        73.502
#SNORA54     66.577
#HAS2        48.142
#MXD4        24.572
#UPK1B       21.923
#FOXO3       19.612
#ESR1        19.475
#IRF1         3.842
#BRCA1        0.000

#ROC
library(pROC)
pred <- predict(fit_rf,test, type = "prob")
roc.multi <- multiclass.roc(test$phenotype, pred)
rs <- list(bacterial = roc.multi$rocs$`bacterial/covid`[[2]], healthy = roc.multi$rocs$`covid/healthy`[[2]], 
    others = roc.multi$rocs$`covid/others`[[1]], viral = roc.multi$rocs$`covid/viral`[[1]])

roc_all = data.frame()
for (i in 1:length(rs)) {
    obj <- data.frame(x = 1-rs[[i]]$specificities, y = rs[[i]]$sensitivities, Category = rep(names(rs[i]), length(rs[[i]]$sensitivities)))
    roc_all <- rbind(roc_all, obj)
}
roc_all <- roc_all[order(roc_all$Category,roc_all$y),]

roc_all_val = data.frame()
for (i in 1:length(rs)) {
    auc <- data.frame(auc = auc(rs[[i]]))
    roc_all_val <- rbind(roc_all_val, auc)
}
roc_all_val <- cbind(roc_all_val, roc_all[match(c("90","227","329","416"), rownames(roc_all)),])


pdf("roc_all.pdf", height = 8, width = 8) 
ggplot <- ggplot(roc_all, aes(x, y, color = Category),alpha = 0.6) + geom_line(lwd = 1.5) + 
    labs(x = "False Positive Rate", y = "True Positive Rate") + 
    geom_abline(intercept = 0, slope = 1, color = "red", 
                linetype = "dashed", lwd=0.8) + theme_classic(base_size = rel(5)) +
    theme(text = element_text(size=rel(5)), 
          #strip.text = element_text(size = rel(4)), 
          legend.text=element_text(size= rel(4)), legend.position="bottom") + 
    geom_label_repel(data=roc_all_val, aes(label = paste("AUC: ", format(round(auc, 4), nsmall = 2))), 
                     xlim = c(0.25,NA), ylim = c(0.3, NA), size = 6, show.legend = F)
ggplot + theme(legend.text=element_text(size= rel(4)))
dev.off()




