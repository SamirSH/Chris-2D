########################################
###### Download Required Packages ######
########################################

library("HIBAG")
library("dplyr")
library("ggplot2")
library("RColorBrewer")
library("tidyr")
library("caret")

########################################

# TRue hla types for dimelli
dimelli.hlatypes = read.csv(file = 'genotypes_metadata_dimelli.csv', sep = ';')
Labels = read.csv(file = 'Dimelli_Predicions_Labels.csv', sep = ';')

# filter for valid HLA types
dimelli.hlatypes = dimelli.hlatypes[dimelli.hlatypes$dataset == "dimelli" & dimelli.hlatypes$hla != "Incomplete" & dimelli.hlatypes$hla != "Unknown", ]

# Download estimated parameters of European Ancestery population for HIBAG
params = get(load(file = "European-HLA4-hg19.RData"))
hla.id = c("A","B","C","DRB1","DQA1","DQB1","DPB1")

# Define trained HLA models
model.A = hlaModelFromObj(params[[hla.id[1]]])
model.B = hlaModelFromObj(params[[hla.id[2]]])
model.C = hlaModelFromObj(params[[hla.id[3]]])
model.DRB1 = hlaModelFromObj(params[[hla.id[4]]])
model.DQA1 = hlaModelFromObj(params[[hla.id[5]]])
model.DQB1 = hlaModelFromObj(params[[hla.id[6]]])
model.DPB1 = hlaModelFromObj(params[[hla.id[7]]])

summary(model.A)
summary(model.B)
summary(model.C)
summary(model.DRB1)
summary(model.DQA1)
summary(model.DQB1)
summary(model.DPB1)



# Convert dimelli test data into HIBAG readible file format
# Convert GDSC data into HIBAG readible file format

dimelli.genom = hlaBED2Geno(bed.fn="plink.bed", fam.fn="plink.fam", bim.fn="plink.bim")
summary(dimelli.genom)

# Predict hla types for Dimelli
pred.A = predict(model.A, dimelli.genom, type="response+prob")
pred.B = predict(model.B, dimelli.genom, type="response+prob")
pred.C = predict(model.C, dimelli.genom, type="response+prob")
pred.DRB1 = predict(model.DRB1, dimelli.genom, type="response+prob")
pred.DQA1 = predict(model.DQA1, dimelli.genom, type="response+prob")
pred.DQB1 = predict(model.DQB1, dimelli.genom, type="response+prob")
pred.DPB1 = predict(model.DPB1, dimelli.genom, type="response+prob")

# Save predictions
save(pred.A, pred.B, pred.C, pred.DRB1, pred.DQA1, pred.DQB1, pred.DPB1, file = "predictions.RData")
# load("predictions.RData")

summary(pred.A)
summary(pred.B)
summary(pred.C)
summary(pred.DRB1)
summary(pred.DQA1)
summary(pred.DQB1)
summary(pred.DPB1)


# Data preprocessing: Define probabilities and filter dimelli test data, remove unnecessary columns
A = pred.A$value
B = pred.B$value
C = pred.C$value
DPB1 = pred.DPB1$value
DQA1 = pred.DQA1$value
DQB1 = pred.DQB1$value
DRB1 = pred.DRB1$value


dimelli.hlatypes$dataset = NULL
dimelli.hlatypes$sex = NULL
dimelli.hlatypes$t1d = NULL
dimelli.hlatypes$t1d_onset_age = NULL
dimelli.hlatypes$fdr = NULL
dimelli.hlatypes$risk_score = NULL
names(dimelli.hlatypes)[1] = names(A)[1]

# check 
names(A)[1] == names(dimelli.hlatypes)[1]


# Cross-map hla type predictions by sample ID
lists = list(DQA1,DQB1,DRB1)

multi_inner = Reduce(
  function(x, y, ...) merge(x, y, by = "sample.id", ...), 
  lists
)

# remove unnecessary columns and assign proper names
multi_inner[c(5,9,13)] = NULL

names(multi_inner) = c("sample.id", 
                       "HLA.DQA1_allele1", "HLA.DQA1_allele2", "DQA1", 
                       "HLA.DQB1_allele1", "HLA.DQB1_allele2", "DQB1",
                       "HLA.DRB1_allele1", "HLA.DRB1_allele2", "DRB1")


# retrieve PredictLabels.R function to convert predicted HLA alleles to Dr-DQ serotypes
source("PredictLabels.R")
preds = PredictLabels(multi_inner)

# add converted columns to data frame
cross = merge(multi_inner, preds, by = "sample.id")

# add also true labels to assess performance in later use case
cross = merge(cross, dimelli.hlatypes, by = "sample.id")

# convert HLA 'DR3/3' to DR3/DR3' for clear reading
cross$hla = gsub('DR3/3', 'DR3/DR3', cross$hla)

# find wrong predictions
bool = cross$prediction == cross$hla
index = which(bool %in% "FALSE")
wrong.pred = cross[index, ]



##########################################################################################
################ Statistics and Confusion Matrix #########################################
##########################################################################################

predicted = cross$prediction
actual = cross$hla
ConfusionMatrix = table(actual, predicted)



# Plot confusion matrix as heatmap


heat = as.data.frame(ConfusionMatrix)

ggplot(data =  heat, mapping = aes(x = actual, y = predicted)) +
  geom_tile(aes(fill = Freq), colour = "white") +
  geom_text(aes(label = sprintf("%1.0f", Freq)), vjust = 1) +
  scale_fill_gradient(low = "grey80", high = "green") +
  theme_bw() + theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  # ggtitle("Confusion Matrix") + theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Actual Classes") +
  ylab("Predicted Classes") +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"), plot.title = element_text(size = 14, face = "bold"))


# Performance assessment

N = sum(ConfusionMatrix) # number of instances
NumClasses = ncol(ConfusionMatrix) # number of classes
diag = diag(ConfusionMatrix) # number of correctly classified instances per class 
rowsums = apply(ConfusionMatrix, 1, sum) # number of instances per class
colsums = apply(ConfusionMatrix, 2, sum) # number of predictions per class


# accuracy metrics 
accuracy = round(sum(diag) / N, 2)
precision = diag / colsums 
recall = diag / rowsums 
f1 = 2 * precision * recall / (precision + recall) 
stats = data.frame(precision, recall, f1) 
stats = tibble::rownames_to_column(stats, "classes")
stats = stats %>% gather(scores, values, -classes)

# Plot Accuracy metrics
ggplot(data =  stats, mapping = aes(x = scores, y = classes)) +
  geom_tile(aes(fill = values), colour = "white") +
  geom_text(aes(label = sprintf("%1.2f", values)), vjust = 1) +
  scale_fill_gradient(low = "grey80", high = "lightgreen") +
  theme_bw() + theme(legend.position = "none") +
  xlab("Scores") +
  ylab("Predicted Classes") +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"), plot.title = element_text(size = 14, face = "bold"))




#plot frequency for predicted classes
freq = table(cross$prediction)
freq = as.data.frame(freq)
names(freq) = c("Classes", "Frequency")


ggplot(data=freq, aes(x=reorder(Classes, -Frequency), y=Frequency)) +
  geom_bar(stat="identity", fill="steelblue")+
  geom_text(aes(label=Frequency), vjust=-0.5, color="black", size=3.5)+
  theme_minimal() +
  #ggtitle("Class frequencies") + theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Predicted classes") +
  ylab("Frequencies") +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"), plot.title = element_text(size = 14, face = "bold")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


