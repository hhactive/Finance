library(caret,quietly=FALSE)
library(dplyr,quietly = FALSE)
#load data
train=read.csv("train.csv",stringsAsFactors =FALSE )
test=read.csv("test.csv",stringsAsFactors = FALSE)
#deal with the empty data with more explanatory 
train$Alley[is.na(train$Alley)] = "No alley access"
test$Alley[is.na(test$Alley)] = "No alley access"
train$BsmtQual[is.na(train$BsmtQual)] = "No basement"
test$BsmtQual[is.na(test$BsmtQual)] = "No basement"
train$BsmtCond[is.na(train$BsmtCond)] = "No basement"
test$BsmtCond[is.na(test$BsmtCond)] = "No basement"
train$BsmtExposure[is.na(train$BsmtExposure)] = "No basement"
test$BsmtExposure[is.na(test$BsmtExposure)] = "No basement"
train$BsmtFinType1[is.na(train$BsmtFinType1)] = "No basement"
test$BsmtFinType1[is.na(test$BsmtFinType1)] = "No basement"
train$BsmtFinType2[is.na(train$BsmtFinType2)] = "No basement"
test$BsmtFinType2[is.na(test$BsmtFinType2)] = "No basement"
train$FireplaceQu[is.na(train$FireplaceQu)] = "No fireplace"
test$FireplaceQu[is.na(test$FireplaceQu)] = "No fireplace"
train$GarageType[is.na(train$GarageType)] = "No garage"
test$GarageType[is.na(test$GarageType)] = "No garage"
train$GarageFinish[is.na(train$GarageFinish)] = "No garage"
test$GarageFinish[is.na(test$GarageFinish)] = "No garage"
train$GarageQual[is.na(train$GarageQual)] = "No garage"
test$GarageQual[is.na(test$GarageQual)] = "No garage"
train$GarageCond[is.na(train$GarageCond)] = "No garage"
test$GarageCond[is.na(test$GarageCond)] = "No garage"
train$PoolQC[is.na(train$PoolQC)] = "No pool"
test$PoolQC[is.na(test$PoolQC)] = "No pool"
train$Fence[is.na(train$Fence)] = "No fence"
test$Fence[is.na(test$Fence)] = "No fence"
train$MiscFeature[is.na(train$MiscFeature)] = "None"
test$MiscFeature[is.na(test$MiscFeature)] = "None"

#deal with the missing value

#replace missing values in MasVnrType column in training and test sets with "None"
#, which is most commonly occuring type
train$MasVnrType[is.na(train$MasVnrType)] = "None"
test$MasVnrType[is.na(test$MasVnrType)] = "None"
#replace missing values in Electrical column in training set with _
#"SBrkr, which is most frequently occuring electrical set-up
train$Electrical[is.na(train$Electrical)] = "SBrkr"
#replace missing values in MSZoning column in test set with RL (most frequently occuring)
test$MSZoning[is.na(test$MSZoning)] = "RL"
#replace missing values in Utilities column in test set with AllPub (most frequently occuring)
test$Utilities[is.na(test$Utilities)] = "AllPub"
#replace missing values in Exterior1st column in test set with VinylSd (most frequently occuring)
test$Exterior1st[is.na(test$Exterior1st)] = "VinylSd"
#replace missing values in Exterior2nd column in test set with VinylSd (most frequently occuring)
test$Exterior2nd[is.na(test$Exterior2nd)] = "VinylSd"
#replace missing value in KitchenQual column in test set with TA (most common)
test$KitchenQual[is.na(test$KitchenQual)] = "TA"
#replace missing values in Functional column in test set with Min2 (most common)
test$Functional[is.na(test$Functional)] = "Min2"
#replace missing value in SaleType column in test set with WD (most common)
test$SaleType[is.na(test$SaleType)] = "WD"

#deal with some numeric variable 

#replace missing value in BsmtFinSF1 column in test set with 0 (has no basement)
test$BsmtFinSF1[is.na(test$BsmtFinSF1)] = 0
#replace missing value in BsmtFinSF2 column in test set with 0 (has no basement)
test$BsmtFinSF2[is.na(test$BsmtFinSF2)] = 0
#replace missing value in BsmtUnfSF column in test set with 0 (has no basement)
test$BsmtUnfSF[is.na(test$BsmtUnfSF)] = 0
#replace missing value in TotalBsmtSF column in test set with 0 (has no basement)
test$TotalBsmtSF[is.na(test$TotalBsmtSF)] = 0
#replace missing values in BsmtFullBath column in test set with 0 (has no basement)
test$BsmtFullBath[is.na(test$BsmtFullBath)] = 0
#replace missing values in BsmtHalfBath column in test set with 0 (has no basement)
test$BsmtHalfBath[is.na(test$BsmtHalfBath)] = 0

#replace missing value in GarageCars column in test set with 0 (no garage)
test$GarageCars[is.na(test$GarageCars)] = 0
#replace missing value in GarageArea column in test set with 0 (no garage)
test$GarageArea[is.na(test$GarageArea)] = 0

#using median imputation 
#replace missing values in LotFrontage column in training and test sets with median
train$LotFrontage[is.na(train$LotFrontage)] = median(train$LotFrontage, na.rm = TRUE)
test$LotFrontage[is.na(test$LotFrontage)] = median(test$LotFrontage, na.rm = TRUE)

#replace missing values in GarageYrBlt column in training and test sets with -999
train$GarageYrBlt[is.na(train$GarageYrBlt)] = -999
test$GarageYrBlt[is.na(test$GarageYrBlt)] = -999
#replace missing values in MasVnrArea column in training and test sets with -999
train$MasVnrArea[is.na(train$MasVnrArea)] = -999
test$MasVnrArea[is.na(test$MasVnrArea)] = -999

#start run!

#linear regression 
Contr = trainControl(method = "cv", number = 5, verboseIter = FALSE)
set.seed(77)
model_lm = train(SalePrice ~ ., 
                 data = train,
                 method = "lm",
                 trControl = Contr)
#check 
model_lm

prediction_lm = predict(model_lm, test)
#make a submission file
submit <- data.frame(Id = test$Id, SalePrice = prediction_lm)

write.csv(submit, file = "submission.csv", row.names = FALSE)
