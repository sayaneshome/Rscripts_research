setwd('/Users/sayaneshome/20newsgroups')
#processing test data and train data 
test_data <- read.csv('test_data.csv',header = TRUE)
test_data_label <- read.csv('test_label.csv',header = TRUE)
train_data <- read.csv('train_data.csv',header = TRUE)
train_data_label <- read.csv('test_label.csv',header = TRUE)
names(train_data) <- c("docIdx", "wordIDx","count")
names(test_data) <- c("docIdx", "wordIDx","count")
train_data_label$docIdx <- seq.int(nrow(train_data_label))
test_data_label$docIdx <- seq.int(nrow(train_data_label))
colnames(test_data_label)[1] <- "newsgroup_ID"
colnames(train_data_label)[1] <- "newsgroup_ID"
test_data_v1 <- merge(test_data, test_data_label, by = "docIdx")
train_data_v1 <- merge(train_data, train_data_label, by = "docIdx")
colnames(test_data_v1)[3] <- "w_d"
colnames(train_data_v1)[3] <- "w_d"

#processing vocabulary text,map.csv
vocab <- read.table('vocabulary.txt',header = FALSE)
vocab$ID <- seq.int(nrow(vocab))
map <- read.csv('map.csv',header = FALSE)

library(plyr)
classifier <- count(train_data_v1, vars=c("docIdx","newsgroup_ID"))
classifier[, c("sum.freq", "length.freq")] <- with(classifier, sapply(c(sum, length), function(x)
ave(freq, newsgroup_ID, FUN = x)))
colnames(classifier)[3] <- "totalnumberofwords/document"
colnames(classifier)[4] <- "totalnumberofwords/newsgroup"
colnames(classifier)[5] <- "numberofdocs/newsgroup"
classifier_1 <- merge(train_data_v1, classifier, by = c("docIdx","newsgroup_ID"))

colnames(vocab)[2] <- "wordIDx"
classifier_2 <- merge(classifier_1, vocab, by = c("wordIDx"))
n_k <- count(classifier_2, vars=c("wordIDx","newsgroup_ID"))
classifier_3 <- merge(classifier_2, n_k, by = c("wordIDx","newsgroup_ID"))

colnames(classifier_3)[9] <- "occurenceofwordID/newsgroup_n_k"
colnames(classifier_3)[6] <- "totalnumberofwords/newsgroup_n"

classifier_3$W_j <- classifier_3$`numberofdocs/newsgroup`/11269
colnames(classifier_3)[10] <- "P(W_j)"

classifier_3$MLE <- classifier_3$`occurenceofwordID/newsgroup_n_k`/classifier_3$`totalnumberofwords/newsgroup_n`
vocab_count <- nrow(vocab)
classifier_3$Bayesianestimator <- (classifier_3$`occurenceofwordID/newsgroup_n_k`+1)/(vocab_count+classifier_3$`totalnumberofwords/newsgroup_n`)


print(unique(classifier_3$newsgroup_ID)) 
print(unique(classifier_3$`P(W_j)`))


#Naive bayes classifier for training data

#Omega_NB_ <- arg.max(log(classifier_3$`P(W_j)`+ (sum(log(classifier_3$Bayesianestimator))))#needs completion

library(dplyr)
classifier_3 %>%
  group_by(docIdx) %>%
  mutate(Omega_NB = pmax(log(P(W_j)+sum(`occurenceofwordID/newsgroup_n_k`))))

#classifier_3$Omega_NB<- pmax(log(classifier_3$`P(W_j)`)+sum(log(clas))


#colnames(classifier_3)[10] <- "Prior"
#priors <- select(classifier_3,newsgroup_ID,Prior)
#priors <- unique(select(classifier_3,newsgroup_ID,Prior))
#View(priors)
#colnames(classifier_3)[11] <- "Prior_l"
#priors <- unique(select(classifier_3,newsgroup_ID,Prior_l))

#Naive bayes classifier for test data 


library(plyr)
classifier_t <- count(test_data_v1, vars=c("docIdx","newsgroup_ID"))
classifier_t[, c("sum.freq", "length.freq")] <- with(classifier_t, sapply(c(sum, length), function(x)
  ave(freq, newsgroup_ID, FUN = x)))
colnames(classifier_t)[3] <- "totalnumberofwords/document"
colnames(classifier_t)[4] <- "totalnumberofwords/newsgroup"
colnames(classifier_t)[5] <- "numberofdocs/newsgroup"
classifier_1_t <- merge(test_data_v1, classifier_t, by = c("docIdx","newsgroup_ID"))

colnames(vocab)[2] <- "wordIDx"
classifier_2_t <- merge(classifier_1_t, vocab, by = c("wordIDx"))
n_k_t <- count(classifier_2_t, vars=c("wordIDx","newsgroup_ID"))
classifier_3_t <- merge(classifier_2_t, n_k_t, by = c("wordIDx","newsgroup_ID"))

colnames(classifier_3_t)[9] <- "occurenceofwordID/newsgroup_n_k"
colnames(classifier_3_t)[6] <- "totalnumberofwords/newsgroup_n"

classifier_3_t$W_j <- classifier_3_t$`numberofdocs/newsgroup`/11269
colnames(classifier_3_t)[10] <- "P(W_j)"

classifier_3_t$MLE <- classifier_3_t$`occurenceofwordID/newsgroup_n_k`/classifier_3_t$`totalnumberofwords/newsgroup_n`
vocab_count <- nrow(vocab)
classifier_3_t$Bayesianestimator <- (classifier_3_t$`occurenceofwordID/newsgroup_n_k_t`+1)/(vocab_count+classifier_3_t$`totalnumberofwords/newsgroup_n`)

#print the values of priors
# subset the dataset into values 
priors <- priors[order(priors$newsgroup_ID),]
View(priors)

classifier_3$Omega_NB<- pmax(log(classifier_3$`P(W_j)`)+
                               (log(classifier_3$`occurenceofwordID/newsgroup_n_k`/classifier_3$`totalnumberofwords/newsgroup_n`) +
                                  log(classifier_3$`occurenceofwordID/newsgroup_n_k`/classifier_3$`totalnumberofwords/newsgroup_n`)),
                             log(classifier_3$`occurenceofwordID/newsgroup_n_k`/classifier_3$`totalnumberofwords/newsgroup_n`))


# how to add the additional column in classifier_3 Doubts : 
#if classifier_3$`P(W_j)`== priors$Prior_l{
 # classifier_3$newsgroup <- rbindpriors[1,] 
#}

classify_compare <- select(classifier_3,classifier_3$docIdx)

#Overall accuracy

#Accuracy of training dataset
#Accuracy_trainingdb <- #Accuracy <- length(intersect(A,train_data_label))/11269
#Accuracy_testdb <- #Accuracy <- length(intersect(A1,test_data_label))/7505

#Confusion Matrix 

#pos_or_neg <- ifelse(probability_prediction > threshold, positive_class, negative_class)
#p_class <- factor(pos_or_neg, levels = levels(test_values))

#confusionMatrix(p_class, test_values)

#to extract the values 

library(dplyr)

classifier_3 %>%
  select(docIdx,MLE)