install.packages('googleAuthR')
install.packages('googleCloudStorageR')
library(googleAuthR)
library(googleCloudStorageR)
library(dplyr)

#Sys.time() ->""2020-02-16 01:06:31 PST""
sm1 <- read.csv('smartwatch.csv')
sm2 <- sm1[sm1$source == "heart_rate",]
sm2$source <- NULL
sm2$timediff <- Sys.time() - as.POSIXct(sm2$timestamp)
sm2$index <- NULL

#heart_rate data processing
df3 <- tidyr::extract(sm2,values, into = c('values', 'D'), "(\\d+\\.\\d+).*(\\d)") %>%dplyr::filter(D == 1) %>%dplyr::select(-D)
df3$timediff <- round(df3$timediff,4)
df3 <- chk
#chk <- read.csv('chk4.csv')
chk2 <- round(chk,4)
chk2$t2 <- chk2$timediff+0.4
chk2$t3 <- chk2$timediff+0.9
chk2$t4 <- chk2$timediff+1.2
chk2$t5 <- chk2$timediff+1.7
colnames(chk2)[2] <- "t"
chk2$t5 = ifelse(round(chk2$t5,0) == round(chk2$t4,0), chk2$heart_rate, NA)
chk2$t4 = ifelse(round(chk2$t4,0) == round(chk2$t3,0), chk2$heart_rate, NA)
chk2$t3 = ifelse(round(chk2$t3,0) == round(chk2$t2,0), chk2$heart_rate, NA)
chk2$t2 = ifelse(round(chk2$t2,0) == round(chk2$t,0), chk2$heart_rate, NA)
chk3 <- data.frame(na.approx(chk2))
colnames(chk3)[1] <- "t"
colnames(chk3)[2] <- "time1"
chk3 = chk3[!duplicated(chk3$time1),]
chk4 <- chk3[complete.cases(chk3), ]

set.seed(123)
smp_size <- floor(0.75 * nrow(chk4))
train_ind <- sample(seq_len(nrow(chk4)), size = smp_size)

train <- chk4[train_ind, ]
test <- chk4[-train_ind, ]

train <- chk4[1:15857,]
test <- chk4[15858:16858,]
train$time2<- round(train$time1,2)
train <- train[complete.cases(train), ]
test <- test[complete.cases(test), ]
min_train <-min(training_data$time1)
min_test <- min(test_data$time1)
max_train <- max(training_data$time1)
max_test <- max(test_data$time1)


linRidgeMod <- linearRidge(time1 ~ t + t2 + t3 + t4 + t5, data = training_data)
predicted <- data.frame(predict(linRidgeMod, test_data))
colnames(predicted)[1] <- "time1"
predicted1 <- predicted
predicted1$time1 <- round(predicted1$time1,3)
test_data$time1 <- round(test_data$time1,3)

predicted_heart_rate <- merge(predicted1,test_data[1:2],by = "time1")[1:2]


sohit_test <- test
sohit_train <- train
write.csv(test_data,'test_data_heartrate_1individual.csv')
write.csv(training_data,'training_data_heartrate_1individual.csv')
predicted_heart_rate <- merge(training_data,predicted1,by = "time1")[1:2]
predicted_heart_rate1$time1 <- round(predicted_heart_rate1$time1,2)
predicted_heart_rate2 <- predicted_heart_rate1[!duplicated(predicted_heart_rate1$time1),]

plot(predicted_heart_rate2$time1,predicted_heart_rate2$t,type = 'l',xlab = "time(in days)",ylab = "heart_rate")


#BMI Data processing


sohit_train <- training_data
sohit_test <- test
sohit_train <- train
sohit_train$time1 <- round(sohit_train$time1,2)
sohit_test$time1 <- round(sohit_test$time1,2)
View(sohit_train)
chk4 <- chk3[complete.cases(chk3), ]
View(chk4)
train <- chk4[1:15857,]
test <- chk4[15858:16858,]
train$time2<- round(train$time1,2)
train <- train[complete.cases(train), ]
test <- test[complete.cases(test), ]
sohit_test <- test
sohit_train <- train
sohit_train$time1 <- round(sohit_train$time1,2)
sohit_test$time1 <- round(sohit_test$time1,2)
View(sohit_test)
View(sohit_train)
sohit_train$time2 <- NULL
View(sohit_train)
sohit_train$time2 <- round(sohit_train$time1,1)
View(sohit_train)
sohit_train$time1 <- NULL
sohit_train %>%
arrange
library(plyr)
ddply(sohit_train, .(time2), summarize,  Rate1=mean(t))
library(plyr)
sohit_train2 <- ddply(sohit_train, .(time2), summarize,  Rate1=mean(t))
View(sohit_train2)
library(plyr)
sohit_train2 <- ddply(sohit_train, .(time2), summarize,  heart_rate=mean(t))
sohit_test2 <- ddply(sohit_test, .(time1), summarize,  heart_rate=mean(t))
View(sohit_test2)
write.csv(sohit_test2,'heart_rate_test_plot.csv')
write.csv(sohit_train2,'heart_rate_train_plot.csv')



#BMI processing





