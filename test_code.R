setwd('/Users/sayaneshome/20newsgroups/')
df <- read.csv('train_data.csv', header = T)
df1 <- read.csv('train_label.csv',header = T)
names(df) <- c("d","w","w_d")
names(df1) <- c("n")
df1$d <- seq.int(nrow(df1))


# for training data
# # setwd('/Users/sayaneshome/20newsgroups/')
# # 
# # df <- read.csv('train_data.csv', header = T)
# # df1 <- read.csv('train_label.csv',header = T)
# # names(df) <- c("d", "w","w_d")
# # names(df1) <- c("n")
# # df1$d <- seq.int(nrow(df1))
# 
# # setwd('/Users/sayaneshome/20newsgroups/')
# # df <- read.table(text =
# #                    "d w w_d
# #                  1 2 2
# #                  2 3 1
# #                  3 4 10
# #                  4 5 11
# #                  5 3 12
# #                  6 5 7
# #                  7 2 9
# #                  8 5 8 ", header = T)
# # 
# # df1 <- read.table(text =
# #                     "n d
# #                   1 1
# #                   1 2
# #                   3 8
# #                   2 7
# #                   2 5
# #                   4 6
# #                   3 4
# #                   1 3 ", header = T)



v1 <- merge(df,df1,by="d",all = TRUE)
library(tidyverse)
c1 <- v1 %>%
  group_by(d) %>%
  mutate(tw_d = sum(w_d))

#classifier_2 <- merge(classifier_1,train_data_label,by = 'docIdx',all=TRUE)
#classifier_2[ , colSums(is.na(classifier_2)) == 0]

c2 <- c1 %>%
  group_by(w,n) %>%
  mutate(w_n = sum(w_d)) # %>%
#summarise(newsgroup_ID,w_n)
#till here

c3 <- c2 %>%
  group_by(n) %>%
  mutate(tw_n = sum(w_n))

docs <- length(unique(c3$d))
p1 <- c3 %>%
  group_by(n) %>%
  mutate(nd = length(unique(d))) 
p1$p_wj <- p1$nd/docs
p1 <- select(p1,n,p_wj)
p1 <- unique(p1)
p1 <- p1[!p1$n == "NA", ]
view(p1)
write.table(p1,'priors.txt')
vc <- length(vocab)
n1 <- c3 %>%
  group_by(n) %>%
  filter(n,tw_n)
n1 <- unique(n1)

#merger <- merge(n1,c2,by = "n",all=TRUE)

w1 <- c3 %>%
  ungroup() %>%
  expand(w, n)
w1 <- w1[!w1$n == "NA", ]
w2 <- c3 %>%
  select(w,n,w_n,tw_n)

#c3$MLE <- c3$w_n/c3$tw_n
#c3$Bayesianestimator <- (c3$w_n+1)/(vc+c3$tw_n)
merger <- merge(w1,w2,by=c('w','n'),all = TRUE)
view(merger)
merger[is.na(merger)] <- 0
merger <- merger %>%
  group_by(n) %>%
  mutate(tw_n = sum(w_n))
#merger <- merger[!merger$d == "0", ]
merger <- unique(merger)
docs <- merge(df,merger,by=c('w'))
docs_1 <- docs %>%
  ungroup() %>%
  expand(d, n)
view(docs_1)
merger <- data.matrix(merge(docs_1,merger,all = TRUE))
view(merger)

#merger <- merger[!merger$d == "0", ]

merger$MLE <- log(merger$w_n/merger$tw_n)
merger$Bayesianestimator <- log((merger$w_n+1)/(4+merger$tw_n))
merger <- merge(merger,p1,by = 'n')
merger$p_wj_log <- log(merger$p_wj)

#MLE calculation
merger_MLE <- subset(merger, MLE!="-Inf")
merger_MLE$Bayesianestimator <- NULL
merger_MLE <- merger_MLE %>%  group_by(d,n) %>% 
  summarise(sumMLE = sum((MLE)))
dd <- merge(merger_MLE,p1,by=c('n'))
dd$p_wj_log <- log(dd$p_wj)
dd <- dd[is.finite(rowSums(dd)),]
#dd$Bayesianestimator <- NULL
#dd$Omega_NB<-dd$sumB+dd$p_wj_log
dd$Omega_MLE<-dd$sumMLE+dd$p_wj_log
maxMLE <- dd %>% group_by(d) %>% filter(Omega_MLE == max(Omega_MLE))

results_MLE <- maxMLE %>%
  group_by(d,n) %>%
  select(Omega_MLE)
results_MLE <- unique(results_MLE)

#empty <- array(numeric(),c(2,3,0)) 

#BE calculation
merger_BE <- merger
merger_BE$MLE <- NULL  
merger_BE <- merger_BE %>%  group_by(d,n) %>% 
  summarise(sumBE = sum((Bayesianestimator)))
dd <- merge(merger_BE,p1,by=c('n'))
dd$p_wj_log <- log(dd$p_wj)
dd$Omega_BE<-dd$sumBE+dd$p_wj_log
#dd$Omega_MLE<-dd$sumMLE+dd$p_wj_log
maxBE <- dd %>% group_by(d) %>% filter(Omega_BE == max(Omega_BE))

results_BE <- maxBE %>%
  group_by(d,n) %>%
  select(Omega_BE)
results_BE <- unique(results_BE)

#docs <- merge(df,merger,by='w')

# library(tidyverse)
# #for MLE
# #sumMLE<-c3 %>%  group_by(d) %>% summarise(sumMLE = sum(log(MLE)))
# docs <- docs %>%  group_by(d) %>% 
#   summarise(sumMLE = sum((MLE)),sumB=sum((Bayesianestimator)))
# dd <- merge(merger,sumMLE,by='d')
# dd$Omega_NB<-dd$sumB+dd$p_wj_log
# dd$Omega_MLE<-dd$sumMLE+dd$p_wj_log
# maxMLE <- dd %>% group_by(d) %>% filter(Omega_MLE == max(Omega_MLE))
# maxBayesianb <- dd %>% group_by(d) %>% filter(Omega_NB == max(Omega_NB))
# write.csv(max,'results.csv')
# 
# results_BE <- maxBayesianb %>%
#   group_by(d,n) %>%
#   select(Omega_NB)
# results_BE <- unique(results_NB)
# 
# results_MLE <- maxMLE %>%
#   group_by(d,n) %>%
#   select(Omega_MLE)
# results_MLE <- unique(results_MLE)

d_count <- length(unique(df1$d))
n_count <- length(unique(df1$n))

predict_MLE_n <- results_MLE %>%
  group_by(d) %>%
  select(n)

predict_BE_n <- results_BE %>%
  group_by(d) %>%
  select(n)

count_MLE_match <- nrow(merge(df1, predict_MLE_n))
accuracy_MLE <- count_MLE_match/d_count
count_BE_match <- nrow(merge(df1, predict_BE_n))
accuracy_BE <- count_BE_match/d_count
print(accuracy_BE)
print(accuracy_MLE)

#confusion matrix

d_count <- length(unique(df1$d))
n_count <- length(unique(df1$n))

M1<-matrix(0L, nrow=n_count, ncol=n_count)
for (i in 1:d_count){
  true=df1$n[df1$d==i]
  predict=results_MLE$n[results_MLE$d==i]
  M1[true,predict]=M1[true,predict]+1
}     
print(M1)
write.table(M1,'confusionmatrix_BE_trainingset.txt')

M<-matrix(0L, nrow=n_count, ncol=n_count)

for (i in 1:d_count){
  true=df1$n[df1$d==i]
  predict=results_BE$n[results_BE$d==i]
  M[true,predict]=M[true,predict]+1
} 
print(M)




