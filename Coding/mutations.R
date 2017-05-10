## library
library(e1071)
library(ggplot2)

# mutation 
# preprocessing the protein info file
aa_info <- read.delim("/Users/jerome/Projects/deleteriousness_estimation/aa_info2.txt",header = T)

## scale the matrix
aa_matrix <- aa_info[,c(2,3,4,5,7,8,10:13)]
for (i in 3:ncol(aa_matrix)) {
    tmp_scale <- scale(aa_matrix[,i])
    aa_matrix[,i] <- tmp_scale
}

sum(aa_matrix[1,2:9] - aa_matrix[2,2:9])^2

score_list <- c()
for (i in 1:20) {
    for (j in 1:20) {
        tmp <- sum(abs(aa_matrix[i,2:9] - aa_matrix[j,2:9]))
        score_list <- c(score_list, tmp)
    }
    
}

write.table(aa_matrix,"scaled_aa_info.txt",quote = F, row.names = F)

## read the polyphen results
polyphen2 <- read.delim("/Users/jerome/Projects/deleteriousness_estimation/polyphen2.txt",header = T, )

## read the deviance data 
data <- read.delim("/Users/jerome/Projects/deleteriousness_estimation/training_data.txt", header = F)
names(data) <- c(names(aa_matrix)[2:10],"outcome")

set.seed(42)
train_index <- sample(1:2831,2000,replace = F) 
train_data <- data[train_index,]
test_data <- data[-train_index,]

## using linear model to see which predictors are significant 
## linear model
model1 <- lm(outcome ~.,data = train_data)
pred1 <- predict(model1, newdata = test_data)
mse1 <- mean((pred1 - test_data$outcome)^2)

## Support Vector Regression
model2 <- svm(outcome ~.,data = train_data)
pred2 <- predict(model2, newdata = test_data)
mse2 <- mean((pred2 - test_data$outcome)^2)

## re weight the matrix based on the previous analysis
aa_matrix2 <- aa_matrix
aa_matrix2[,c(3,6,7,8)] <- aa_matrix2[,c(3,6,7,8)]*2
write.table(aa_matrix2,"scaled_aa_info2.txt",quote = F, row.names = F)

## Compare the results
ph_ds_df <- read.delim("/Users/jerome/Projects/deleteriousness_estimation/ph_ds_list.txt",header = F)
ph_top_index <- order(ph_ds_df$V1,decreasing = T)[1:200]
ds_top_index <- order(ph_ds_df$V2,decreasing = T)[1:200]
sum(!is.na(match(ph_top_index,ds_top_index)))

ph_top_index <- order(ph_ds_df$V1,decreasing = F)[1:200]
ds_top_index <- order(ph_ds_df$V2,decreasing = F)[1:200]
sum(!is.na(match(ph_top_index,ds_top_index)))

pred_ds <- read.delim("/Users/jerome/Projects/deleteriousness_estimation/deleterious_score.txt",header = F)

## compare the overall results
ggplot(pred_ds, aes(x=V4)) +geom_histogram(position="identity", bins = 60) + ggtitle("V2 Estimation Score", subtitle = NULL) + xlab("Score")
ggplot(ph_ds_df, aes(x=V1)) +geom_histogram(position="identity", bins = 60) + ggtitle("Polyphen2 Estimation Score",subtitle = NULL) + xlab("Score")

