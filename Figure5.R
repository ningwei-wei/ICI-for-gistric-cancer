library(data.table)
library(CIBERSORT)
library(dplyr)
library(limma)
library(stringr)
########################################
FPKM_data<-fread("d:/gastric cancer/data/TCGA-STAD.htseq_fpkm-uq.tsv.gz",data.table = F)
id_data<-fread("d:/gastric cancer/data/gencode.v22.annotation.gene.probeMap",data.table = F)
id_data<-id_data[,1:2]
colnames(id_data)[1]<-colnames(FPKM_data)[1]

FPKM_data<-inner_join(id_data,FPKM_data)
FPKM_data[,-c(1:2)]<-2^FPKM_data[,-c(1:2)]-1

fpkmToTpm <- function(fpkm)
{
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}
FPKM_data[,-c(1:2)]<-apply(FPKM_data[,-c(1:2)], 2, fpkmToTpm)
data_clean<-avereps(FPKM_data[,-c(1:2)],ID=FPKM_data$gene)

a<-grep("01A",colnames(data_clean))
data_clean<-data_clean[,a]
data_scaled<-log2(data_clean+0.1)
save(data_scaled,file = "G:/gastric cancer/TCGA_scale.Rdata")

####################################################
load("g:/gastric cancer/Sig.RData")
write.table(Sig,"g:/gastric cancer/Allcell_sig.txt",sep = "\t",col.names = T,row.names = T)    

results1 <- cibersort(sig_matrix ="g:/gastric cancer/Allcell_sig.txt", mixture_file = data_scaled,
                      perm = 100,QN = T)

write.csv(results,file = "d:/gastric cancer/result/cibersort/Allcell.csv")
results<-read.csv("d:/gastric cancer/result/cibersort/Allcell.csv")
rownames(results)<-results$X
results<-results[,-1]

###################### Tcell

results<-t(results)
data_clean1<-rbind(results,data_clean)

data_clean1<-data_clean1[-c(2:14),]

scaling_factor <- data_clean1[1, ]  # 第一行，长度等于列数

# 使用sweep：每一行乘以第一行（按列缩放）,提取Tcell表达谱
data_scaled <- sweep(data_clean1[-1,], 2, scaling_factor, FUN = "*")

Tcell_results <- cibersort(sig_matrix ="g:/gastric cancer/Tcell_sig.txt", 
                           mixture_file = data_scaled,perm = 100,QN = T)
Tcell_results<-as.data.frame(Tcell_results)
Tcell_results$ID<-rownames(Tcell_results)



a<-t(results)
a<-as.data.frame(a)
a$ID<-rownames(a)

data<-inner_join(Tcell_results,a,by = "ID")

#################### core
data$ratio <- (data$`PD-L1_sen`/data$`PD-L1_res`)
data$ratio<-data$ratio*data$T_cells_NK_cells
##################

data$group<-ifelse(data$ratio>median(data$ratio),"High","Low")

pheno_data<-read.table("d:/gastric cancer/data/TCGA-STAD.survival.tsv",sep = "\t",header = T)
pheno_data<-pheno_data[,c(1,2,4)]
colnames(pheno_data)[1]<-"ID"

data<-inner_join(pheno_data,data)
############################
library("survival")
library("survminer")
library(data.table)
library(dplyr)

data$OS.time<-data$OS.time/30
write.csv(data,file = "G:/gastric cancer/TCGA_group.csv")

######################################
value <- surv_cutpoint(data, time = "OS.time", event = "OS", variables = "ratio")
cut_off <- as.numeric(value[["cutpoint"]][1, 1])


data <- data %>% 
  dplyr::mutate(new_group = ifelse(data[, "ratio"] >= cut_off,"High","Low"))

#########################  Figure 5A
fit <- survfit(Surv(OS.time, OS) ~ new_group, data = data)

# 手动计算Wilcoxon检验的p值
survdiff_obj <- survdiff(Surv(OS.time, OS) ~ new_group, data = data)
p_value_wilcoxon <- 1 - pchisq(survdiff_obj$chisq, length(survdiff_obj$n) - 1)

library(gridExtra)
p<-ggsurvplot(fit, data = data,
              pval=F, #显示p值
              pval.coord = c(0, 0.2),#p值位置坐标
              pval.size =5,#p值字体大小
              pval.method = T,
              #conf.int=T, #是否显示生存率的95%CI
              risk.table=TRUE,#显示风险表
              risk.table.height = 0.25, #风险表的高度
              palette="jco", #柳叶刀配色
              #legend.labs=c("Cold Tumor", "Hot Tumor"), 
              legend.title="",  #改图例名称 
              surv.median.line = "hv" #中位生存期。
)
p
# 合并主图和风险表
combined_plot <- grid.arrange(p$plot, p$table, ncol = 1, heights = c(3, 1))

ggsave("d:/gastric cancer/12_15/plot/TCGA_KM.pdf",combined_plot,width = 8,height = 6)    

########################
library(DALEX)
library(DALEXtra)
library(dplyr)
library(data.table)

load("G:/gastric cancer/TCGA_scale.Rdata")
load("g:/gastric cancer/Tcell_sig.Rdata")


TF_gene<-read.csv("g:/gastric cancer/TF_gene.csv")
deg_gene<-unlist(sig)

a<-intersect(rownames(data_scaled),c(TF_gene$tf,deg_gene))

data<-data_scaled[a,]
data<-t(data)
data<-as.data.frame(data)
data$ID<-rownames(data)

group<-read.csv("g:/gastric cancer/Tcga_subgroup.csv")
group<-group[,c(5,31)]

pheno<-fread("d:/gastric cancer/data/TCGA-STAD.GDC_phenotype.tsv.gz")
pheno<-pheno[,c(1,6,71)]
colnames(pheno)<-c("ID","Age","Sex")



data<-inner_join(group,data)
data<-inner_join(data,pheno)
data$Sex<-ifelse(data$Sex=="male",1,0)


#############################
# 加载必要的包
library(caret)        # 用于数据分割和模型训练
library(randomForest) # 随机森林
library(e1071)        # SVM
library(DALEX)        # 模型解释
library(pROC)         # ROC曲线
library(ROSE)         # 处理不平衡数据
library(dplyr)        # 数据处理
library(xgboost)      # XGBoost

# 设置随机种子确保可重复性
set.seed(123)

data<-na.omit(data)
# 检查类别分布
table(data$new_group)
prop.table(table(data$new_group))
data<-data[,-1]

colnames(data)<-gsub("-","_",colnames(data))

# 方法1：使用createDataPartition进行分层抽样（保持类别比例）
train_index <- createDataPartition(data$new_group, 
                                   p = 0.7,           # 70%训练集
                                   list = FALSE,
                                   times = 1)
train_data <- data[train_index, ]
test_data <- data[-train_index, ]

# # 方法2：使用SMOTE处理训练集的不平衡（可选但推荐）
# train_data_balanced <- ROSE(new_group ~ ., data = train_data, seed = 123)$data
# table(train_data_balanced$new_group)


# 准备特征矩阵和标签
prepare_data <- function(data) {
  x <- data %>% select(-new_group) %>% as.data.frame()
  y <- as.numeric(as.factor(data$new_group))
  return(list(x = x, y = y))
}


train_orig <- prepare_data(train_data)
test <- prepare_data(test_data)

# 为XGBoost准备DMatrix格式
dtrain <- xgb.DMatrix(data = as.matrix(train_orig$x), 
                      label = as.integer(train_orig$y) - 1 )
dtest <- xgb.DMatrix(data = as.matrix(test$x), 
                     label = as.integer(test$y)-1 )

################################################
# 定义特征和标签
x_train <- train_data %>% select(-new_group)
y_train <- as.factor(train_data$new_group)
x_test <- test_data %>% select(-new_group)
y_test <- as.factor(test_data$new_group)


fitControl <- trainControl(
  method = "cv",           
  number = 5,              
  classProbs = TRUE,       
  summaryFunction = twoClassSummary,  
  savePredictions = "final",
  sampling = "up"         
)

# 2.1 
rf_model <- train(
  x = x_train,
  y = y_train,
  method = "rf",
  trControl = fitControl,
  tuneLength = 5,          
  metric = "ROC",          
  importance = TRUE        
)

# 2.2 
glm_model <- train(
  x = x_train,
  y = y_train,
  method = "glm",
  family = binomial(),
  trControl = fitControl,
  metric = "ROC"
)

# 2.3 
svm_model <- train(
  x = x_train,
  y = y_train,
  method = "svmRadial",    
  trControl = fitControl,
  tuneLength = 5,          
  metric = "ROC",
  prob.model = TRUE        
)
# 2.4 xgboost 

class_table <- table(train_orig$y)
scale_pos_weight <- class_table[1] / class_table[2] 


xgb_params <- list(
  objective = "binary:logistic",
  eval_metric = "auc",
  max_depth = 6,
  eta = 0.1,
  gamma = 0,
  colsample_bytree = 0.8,
  min_child_weight = 1,
  subsample = 0.8,
  scale_pos_weight = scale_pos_weight,
  tree_method = "hist"
)

xgb_cv <- xgb.cv(
  params = xgb_params,
  data = dtrain,
  nrounds = 100,
  nfold = 5,
  showsd = TRUE,
  stratified = TRUE,
  print_every_n = 10,
  early_stopping_rounds = 50,
  maximize = FALSE
)
print(xgb_cv$evaluation_log)

best_nrounds <- 35

xgb_model <- xgboost(
  params = xgb_params,
  data = dtrain,
  nrounds = best_nrounds,
  verbose = 0
)

########################
evaluate_model_enhanced <- function(model, x_test, y_test, model_name, 
                                    model_type = "auto") {
  
  if (model_type == "auto") {
    if (inherits(model, "xgb.Booster")) {
      model_type <- "xgb_direct"
    } else if (inherits(model, "train")) {
      model_type <- "caret"
    } else if (inherits(model, "svm")) {
      model_type <- "svm_direct"
    } else if (inherits(model, "randomForest")) {
      model_type <- "rf_direct"
    } else if (inherits(model, "glm")) {
      model_type <- "glm_direct"
    } else {
      stop("no model")
    }
  }
  
  
  if (model_type == "xgb_direct") {
    # xgb.Booster
    if (!is.matrix(x_test)) {
      x_test_matrix <- as.matrix(x_test)
    } else {
      x_test_matrix <- x_test
    }
    dtest <- xgb.DMatrix(data = x_test_matrix)
    prob_vector <- predict(model, dtest)
    
  } else if (model_type == "caret") {
    # caret
    pred_prob <- predict(model, newdata = x_test, type = "prob")
    prob_vector <- pred_prob[, 2]  # 第二列是正类概率
    
  } else if (model_type == "svm_direct") {
    # SVM
    prob_vector <- predict(model, newdata = x_test, probability = TRUE)
    prob_vector <- attr(prob_vector, "probabilities")[, 2]
    
  } else if (model_type == "rf_direct") {
    # RF
    prob_vector <- predict(model, newdata = x_test, type = "prob")[, 2]
    
  } else if (model_type == "glm_direct") {
    # GLM
    prob_vector <- predict(model, newdata = x_test, type = "response")
    
  } else {
    stop("no model: ", model_type)
  }
  
  if (!is.factor(y_test)) {
    y_test <- factor(y_test)
  }
  
  roc_obj <- roc(y_test, prob_vector, quiet = TRUE)
  best_thresh <- coords(roc_obj, "best", ret = "threshold", 
                        best.method = "youden")$threshold
  
  
  pred_class <- ifelse(prob_vector > best_thresh, 
                       levels(y_test)[2], 
                       levels(y_test)[1])
  pred_class <- factor(pred_class, levels = levels(y_test))
  
  
  cm <- confusionMatrix(pred_class, y_test, positive = levels(y_test)[2])
  
  auc_val <- auc(roc_obj)
  precision <- cm$byClass["Precision"]
  recall <- cm$byClass["Recall"]
  f1 <- 2 * (precision * recall) / (precision + recall)
  specificity <- cm$byClass["Specificity"]
  
  # 绘制ROC曲线
  plot(roc_obj, main = paste(model_name, "ROC"),
       print.auc = TRUE, auc.polygon = TRUE, grid = TRUE)
  points(1 - specificity, recall, col = "red", pch = 19, cex = 1.2)
  
  return(list(
    model_name = model_name,
    model_type = model_type,
    confusion = cm,
    auc = auc_val,
    roc = roc_obj,
    predictions = pred_class,
    probabilities = prob_vector,
    best_threshold = best_thresh,
    metrics = data.frame(
      AUC = auc_val,
      Accuracy = cm$overall["Accuracy"],
      Sensitivity = recall,
      Specificity = specificity,
      Precision = precision,
      F1 = f1
    )
  ))
}

results_enhanced <- list()

results_enhanced$RF <- evaluate_model_enhanced(rf_model, test$x, test$y, 
                                               "RF", model_type = "caret")
results_enhanced$GLM <- evaluate_model_enhanced(glm_model, test$x, test$y, 
                                                "GLM", model_type = "caret")
results_enhanced$SVM <- evaluate_model_enhanced(svm_model, test$x, test$y, 
                                                "SVM", model_type = "caret")
results_enhanced$XGBoost_direct <- evaluate_model_enhanced(xgb_model, test$x, test$y, 
                                                           "XGBoost", 
                                                           model_type = "xgb_direct")
########################## Figure 5B
library(plotROC)


roc_data <- data.frame(
  predictions = c(
    results_enhanced$RF$probabilities,
    results_enhanced$GLM$probabilities,
    results_enhanced$SVM$probabilities,
    results_enhanced$XGBoost_direct$probabilities
  ),
  labels = rep((test$y-1), 4),
  model = rep(c("Random Forest", "GLM", "SVM", "XGBoost"), 
              each = length(test$y))
)

# 
p<-ggplot(roc_data, aes(d = labels, m = predictions, color = model)) +
  geom_roc(n.cuts = 0, size = 1) +  
  style_roc(theme = theme_gray) +
  labs(title = "ROC Curves Comparison") +
  scale_color_manual(values = color_use) +
  theme(legend.position = "bottom",
        legend.title = element_blank()) +
  annotate("text", x = 0.7, y = 0.25, 
           label = paste("Random Forest AUC = 0.896\n",
                         "GLM AUC = 0.738\n",
                         "SVM AUC = 0.905\n",
                         "XGBoost AUC = 0.898"),
           size = 3, hjust = 0)+theme_bw()
ggsave("g:/gastric cancer/plot/ROC.pdf",p,width = 8,height = 6)

############################### Figure 5C
library(ggsci)

color_use<-pal_frontiers()(10)

#########################
xgb_predict_function <- function(model, newdata) {
  if (is.data.frame(newdata)) {
    newdata_matrix <- data.matrix(newdata)
  } else {
    newdata_matrix <- as.matrix(newdata)
  }
  pred_prob <- predict(model, newdata_matrix)
  return(pred_prob)
}

# 1.2 
explainer_XGBoost <- DALEX::explain(
  model = xgb_model,
  data = x_test,
  y = as.numeric(y_test) - 1,  # 
  predict_function = xgb_predict_function,
  label = "XGBoost",
  verbose = TRUE
)


per_xgb <- model_performance(explainer_XGBoost)
plot(per_xgb)

explainer_rf <- DALEX::explain(
  model = rf_model,
  data = x_test,
  y = as.numeric(y_test) - 1,
  predict_function = caret_predict_for_dalex,
  label = "RF",
  verbose = TRUE
)

explainer_glm<-DALEX::explain(glm_model,label = "GLM",
                              data = x_test,
                              y = as.numeric(y_test)-1)

explainer_svm<-DALEX::explain(svm_model,label = "SVM",
                              data = x_test,
                              y = as.numeric(y_test)-1)

per_rf<-model_performance(explainer_rf)
per_glm<-model_performance(explainer_glm)
per_svm<-model_performance(explainer_svm)
per_xgb<-model_performance(explainer_XGBoost)


pdf("g:/gastric cancer/plot/redidual.pdf",width = 8,height = 6)
plot(per_rf,per_glm,per_svm,per_xgb)
dev.off()
plot(per_rf,per_glm,per_svm,per_xgb,geom = "boxplot")

########################## Figure 5D

###############################
library(shapviz)

shp <- shapviz(xgb_model, X_pred = data.matrix(x_train), X = x_train)
sv_importance(shp, kind = "bar")

pdf("g:/gastric cancer/plot/shape_1.pdf",width = 8,height = 6)
sv_importance(shp, kind = "both", show_numbers = TRUE, 
              color_bar_title = "Feature value")
dev.off()

sv_importance(shp, kind = "bee")
sv_waterfall(shp, row_id = 1) +
  ggtitle("Waterfall plot for second prediction")

shap.dependence_plot("CCL5", shap_values, X)


shapviz_res_list <- list()

for (i in seq_along(shap_values_list)) {
  shapviz_res_list[[i]] <- shapviz(shap_values_list[[i]])
  
  cat("Processed shapviz for sample:", i, "\n")
}

print(shapviz_res_list)

pdf("g:/gastric cancer/plot/shape.pdf",width = 8,height = 6)
sv_waterfall(shapviz_res_list[[2]])
dev.off()



############################### Figure 5E
library(limma)

gbm_data<-readRDS("g:/gastric cancer/Melanoma-GSE91061.Response.Rds")
gbm_data<-as.data.frame(gbm_data)
gbm_data<-na.omit(gbm_data)

fpkmToTpm <- function(fpkm)
{
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}
a<-apply(gbm_data[,-c(1)], 2, fpkmToTpm)
data_clean<-avereps(a,ID=gbm_data$GENE_SYMBOL)
gbm_data<-as.data.frame(t(data_clean))

gbm_pheno<-read.table("g:/gastric cancer/Melanoma-GSE91061.Response.tsv",sep = "\t",fill = T)
gbm_pheno<-gbm_pheno[,c(1,6,7)]

a<-intersect(colnames(x_train),colnames(gbm_data))
#a<-c("CCL5","IRF1","GZMA","LAG3","CCR7")

gbm_data<-gbm_data[,a]
gbm_data$sample_id<-rownames(gbm_data)

gbm_data<-inner_join(gbm_pheno,gbm_data)

gbm_data<-gbm_data[-which(gbm_data$response_NR=="UNK"),]

x_stad<-gbm_data[,-c(1:3)]
x_stad<-log2(x_stad+0.1)
y_stad<-ifelse(gbm_data$response_NR=="R","High","Low")
y_stad<-factor(y_stad)


x_train<-x_train[,a]
x_test<-x_test[,a]
####################################
svm_model <- train(
  x = x_train,
  y = y_train,
  method = "svmRadial",
  trControl = fitControl,
  tuneGrid = expand.grid(
    sigma = seq(0.01, 0.1, length.out = 5), 
    C = 2^seq(-5, 10, length.out = 8)        
  ),
  metric = "ROC",
  prob.model = TRUE
)

evaluate_model_enhanced(svm_model, x_test,y_test, 
                        "SVM", model_type = "caret")

a<-evaluate_model_enhanced(svm_model, x_stad,y_stad, 
                           "SVM", model_type = "caret")


#####################
# install.packages("plotROC")
library(plotROC)


roc_data <- data.frame(
  predictions = c(
    a$probabilities
  ),
  labels = rep((y_stad), 1),
  model = rep(c("STAD-PRJEB25780"), 
              each = length(y_stad))
)


p<-ggplot(roc_data, aes(d = labels, m = predictions, color = model)) +
  geom_roc(n.cuts = 0, size = 1) +  # n.cuts=0去掉切点标记
  style_roc(theme = theme_gray) +
  labs(title = "ROC Curves Comparison") +
  scale_color_manual(values = color_use) +
  theme(legend.position = "bottom",
        legend.title = element_blank()) +
  annotate("text", x = 0.7, y = 0.25, 
           label = paste(" AUC = 0.73\n"
           ),
           size = 3, hjust = 0)+theme_bw()
ggsave("g:/gastric cancer/plot/PRJEB25780_ROC.pdf",p,width = 8,height = 6)

#################################### Figure 5F
library(tidyr)

result_matrix <- as.data.frame(table(a$predictions, gbm_data$response))
print(result_matrix)
result_matrix<-result_matrix%>%group_by(Var2)%>%mutate(pro= Freq/sum(Freq))
result_matrix<-result_matrix[,-3]

chji_width <- spread(result_matrix, key = "Var2",
                     value = "pro")
chji_width<-as.data.frame(chji_width)
rownames(chji_width)<-chji_width[,1]
chji_width<-chji_width[,-1]

pdf("g:/gastric cancer/plot/confus.pdf",width = 8,height = 4)
pheatmap::pheatmap(chji_width,scale = "none",display_numbers = T,angle_col = 0)
dev.off()




