# install the packages if not installed on your computer 
# install.packages("glmnet", repos='http://cran.us.r-project.org')
# install.packages("performance",repos = "https://cloud.r-project.org")

library(glmnet, quietly = TRUE)
library(faraway)
library(broom)
library(glmnet)
library(tidyverse)
library(performance)
library(insight)
library(magrittr)
library(ISLR2)

# Reusable function for reading csv file 
# File must be saved on your computer in downloads. 

# Outer wrapped function 

Final <- function(){

if(!exists('d')){
  d <- read.csv("WineQT.csv", header = T)
}


#d<-read.csv("/Users/mariahdiaz/Desktop/WineQT.csv")
  
# Downloaded wine data set used as example in the project below
# Aquired from Kaggle: https://www.kaggle.com/datasets/yasserh/wine-quality-dataset 


head(d)
x <- c(1:5)
y<- c(6:13)
m<- matrix(c(x,y), nrow = 1143, ncol = 13)
class(m)

#For the specified data 
datamat <- as.matrix(d)


#Pre processing <-- Change the function pic in here. 

preproc <- function(d,quality){
  temp <- na.omit(d)
  y <- temp$quality
  indx <- which(names(d) == d$quality)
  x <- temp[ ,-indx]
  return(list(x,y))
  
}

#Step 1 function 
print("Step 1: LASSO, RIDGE, OLS")
  d <- na.omit(d)
  x <- model.matrix(quality ~., data = d)[, -1]
  y <- d$quality
  
  grid.lambda <- 10^seq(10, -2, length = 100)
  ridge.model <- glmnet(x, y, alpha = 0, lambda = grid.lambda)
  
  #Plot the L1 norm of the coefficients
  print("plot l1 norm of the coefficients")
  plot(ridge.model)
  
  ell2.norm <- numeric()
  for(i in 1:length(grid.lambda)){
    ell2.norm[i] <- sqrt(sum(coef(ridge.model)[-1, i]^2))
  }
  
  plot(x = grid.lambda, y = ell2.norm, xlab = expression(lambda), ylab = "L2 Norm", xlim = c(10,10000))
  
  set.seed(1)
  train <- sample(1:nrow(x), nrow(x) / 2)
  test <- (-train)
  y.train <- y[train]
  y.test <- y[test]
  
  ridge.model.train <- glmnet(x[train, ], y.train, alpha = 0, lambda = grid.lambda)
  
  set.seed(1)
  cv.out <- cv.glmnet(x[train, ], y.train, alpha = 0)
  plot(cv.out)
  
  best.lambda <- cv.out$lambda.min
  plot(cv.out)
  abline(v = log(best.lambda), col = "blue", lwd = 2)
  
  ridge.pred <- predict(ridge.model.train, s = best.lambda, newx = x[test, ])
  mspe.ridge <- mean((ridge.pred - y.test)^2)
  
  final.model <- glmnet(x, y, alpha = 0, lambda = best.lambda)
  Coef.Ridge <- coef(final.model)[1:13, ]
  
  #First, let's look at the shrinkage effects of Lasso on the entire data set
  lasso.model <- glmnet(x, y, alpha = 1, lambda = grid.lambda)
  plot(lasso.model)
  
  #Now, fit a Lasso regression model to the training data
  lasso.model.train <- glmnet(x[train, ], y.train, alpha = 1, lambda = grid.lambda)
  
  
  #Perform cross validation on the training set to select the best lambda
  set.seed(1) #For reproducibility
  cv.out <- cv.glmnet(x[train, ], y.train, alpha = 1)
  plot(cv.out)
  
  #Find the best lambda value
  best.lambda <- cv.out$lambda.min
  plot(cv.out)
  abline(v = log(best.lambda), col = "blue", lwd = 2)
  
  #Calculate the MSPE of the model on the test set
  lasso.pred <- predict(lasso.model.train, s = best.lambda, newx = x[test,])
  mspe.lasso <- mean((lasso.pred - y.test)^2)
  
  #Fit the final model to the entire data set using the chosen lambda
  final.model <- glmnet(x, y, alpha = 1, lambda = best.lambda)
  Coef.Lasso <- coef(final.model)[1:13,]
  
  #Shrinkage effects of Lasso 
  EN.model <- glmnet(x, y, alpha = 0.5, lambda = grid.lambda)
  plot(EN.model)
  
  #Lasso regression model to the training data
  EN.model.train <- glmnet(x[train, ], y.train, alpha = 0.5, lambda = grid.lambda)
  
  #Cross validation on the training set 
  set.seed(1)
  cv.out <- cv.glmnet(x[train, ], y.train, alpha = 0.5)
  plot(cv.out)
  
  #Best lambda value
  best.lambda <- cv.out$lambda.min
  plot(cv.out)
  abline(v = log(best.lambda), col = "blue", lwd = 2)
  
  #MSPE of the model test set
  EN.pred <- predict(EN.model.train, s = best.lambda, newx = x[test,])
  mspe.EN <- mean((EN.pred - y.test)^2)
  
  #Final model to the entire data set
  final.model <- glmnet(x, y, alpha = 0.5, lambda = best.lambda)
  Coef.EN <- coef(final.model)[1:13,]
  Coefficients <- data.frame(Ridge = Coef.Ridge, Lasso = Coef.Lasso, Elastic.Net = Coef.EN)
  MSPE <- data.frame(Ridge = mspe.ridge, Lasso = mspe.lasso, Elastic.Net = mspe.EN)

#Step 2 Function 
print("Step 2: Dealing with Outliers")

# Every time this chunk is run, more outliers relative to the current data 
# (the data set updates at the end of this chunk)
# are removed from the data frame. Please be mindful of how many times this 
# chunk is run before moving to test other chunks reading the same data frame. 

  model <- lm(quality ~ chlorides+total.sulfur.dioxide+volatile.acidity+sulphates+alcohol, data=d)
  summary(model)
  
  cooksd <- cooks.distance(model)
  
  plot(cooksd, pch="*", cex=2, main="Influential Obs by Cooks distance")  # plot cook's distance
  abline(h = 4*mean(cooksd, na.rm=T), col="red")  # cutoff line
  text(x=1:length(cooksd)+1, y=cooksd, labels=ifelse(cooksd>4*mean(cooksd, na.rm=T),names(cooksd),""), col="red")  #labels
  influential <- as.numeric(names(cooksd)[(cooksd > 4*mean(cooksd, na.rm=T))])

print("Outliers/influential points can be observed in the plot highlighted with red labels.")
print("In our data set we have decided to removes all points for which Cook’s distance") 
print("for that point is at least 4 times the mean Cook’s distance for all the points.")

#Step 3 function 
print("Step 3:Performing Model Selection")
# Make sure to load islr and leaps if you do not already have it. 

install.packages("ISLR2",repos = "https://cloud.r-project.org")
library(ISLR2)
library(leaps)


  regfit.full <- regsubsets(quality ~ ., d)
  summary(regfit.full)
  
  regfit.full <- regsubsets(quality ~., data = d ,
                            nvmax = 100)
  
  reg.summary <- summary(regfit.full)
  names (reg.summary)
  reg.summary$rsq
  
  plot (reg.summary$rss , xlab = " Number of Variables ",
        ylab = " RSS ", type = "l")
  print("Looking to minimize the RSS")
  plot (reg.summary$adjr2 , xlab = " Number of Variables ",
        ylab = " Adjusted RSq ", type = "l")
  best<-which.max(reg.summary$adjr2)
  
  points(best,reg.summary$adjr2[best],col = " red ",cex = 2,pch = 20)
  print("Looking for highest R squared")
  
  
  plot(reg.summary$cp,xlab = " Number of Variables ",ylab = "Cp",type = "l")
  
  bestcp<-which.min(reg.summary$cp)
  
  points(bestcp, reg.summary$cp[bestcp], col = " red ", cex = 2,pch = 20)
  
  print("Looking for lowest CP")
  

  
  plot(reg.summary$bic , xlab = " Number of Variables ",ylab = " BIC ", type = "l")
  bestBic<-which.min(reg.summary$bic)
  points(bestBic, reg.summary$bic[bestBic], col = " red ", cex = 2,pch = 20)
  
  print("Want smallest BIC value")
  
  
  print("BIC Plot")
  plot(regfit.full , scale = "r2")
  
 print("Adjusted R squared plot")
  plot(regfit.full , scale = "adjr2")
  
  print("Cp plot")
  plot(regfit.full , scale = "Cp")
  
  print("Bic plot")
  plot(regfit.full , scale = "bic")
  
 
  model<-lm(d$quality~., data=d)
  print(summary(model))
  print("We can identify which variables are significant by looking at the") 
  print("codes in the summary. For our wine set this would include, volatile-")
  print("acidity, chlorides, total sulfer, sulphates, and alcohol")
  
  regfit.full <- regsubsets(d$quality~d$volatile.acidity+d$chlorides+d$total.sulfur.dioxide+d$sulphates+d$alcohol, data=d)
  reg.summary <- summary(regfit.full)
  print("These are the variables we have decided to include in our model.")
  print("5 of the most signifcant variables were chosen based on our")
  print("exploratory plots and summaries.")



#Step 4: F-test
print("F-Test:")
  Full.Model<-lm(d$quality ~ ., data=d)
  Reduced.Model<-lm(d$quality~d$volatile.acidity+d$chlorides+d$total.sulfur.dioxide+d$sulphates+d$alcohol, data=d)
 
#Run on chosen model vs. full original model 
 print(anova(Reduced.Model,Full.Model)) 
  
print("H0: No significant difference in these two models: P-value>0.05")
print("Reduced model is significantly better: p-value<0.05")
print("When we run the f-test on our chosen Wine data set, we are comparing the") 
print("full model with the reduced model with a select few variables. From the") 
print("output we can see our p-value is 0.03247 which is less than our")
print("significance level of 0.05 and we reject the null hypothesis that there")
print("is no significant difference between the two models. From this output we")
print("have more evidence to go with the reduced chosen model. ")


# Step 5
print("Step 5:")
model<-lm(d$quality~d$volatile.acidity+d$chlorides+d$total.sulfur.dioxide+d$sulphates+d$alcohol, data=d)

#Linearity 
print("Observing linearity")
print("We are using linearity to see if there is any noticeable trend between")
print("the residuals.")
fitted.values <- model$fitted.values
residuals <- model$residuals

plot(y = fitted.values, x = residuals, xlab = "Residuals", ylab = "Fitted Values") + abline(h = mean(fitted.values), col = "red")

cor(fitted.values, residuals)
print("Homoscedasticity:")

#Homoscedasticity
print("If the data is homoscedastic, the points on the plot should be evenly")
print("distributed #around the zero line.")
print("We use this to check whether the spread of the residuals change as a")
  print("function of the predictors and/or fitted values.") 

plot(y = fitted.values, x = residuals, xlab = "Residuals", ylab = "Fitted Values") + abline(h = mean(fitted.values), col = "blue")

#normality
print("Testing for Normality")
print("We use normality to check if the error terms are iid normal.")
par(mfrow = c(1, 2))
qqnorm(residuals)
hist(residuals, n = 40)
#testing the relationship between residuals and normal distribution using the K-S test
mean.resid <- mean(residuals)
sd.resid <- sd(residuals)

print("Plots appear to show normality.")

normal.sample <- rnorm(1000, mean = mean.resid, sd = sd.resid)
print(ks.test(residuals, normal.sample)) #p-value = 0.04086 
print("with p-value<0.05 we fail to reject the Null Hypothesis")
  
  
  
  # Step 6: Determining with transformations to apply 
print("Step 6: Transformations")
  library(MASS)
  
  #create data
  y=c(fitted.values)
  x=c(residuals)
  
  #fit linear regression model
  model <- lm(y~x)
  
  #find optimal lambda for Box-Cox transformation 
  print("Find optimal lambda for Box-Cox transformation")
  bc <- boxcox(y ~ x)
  (lambda <- bc$x[which.max(bc$y)])
  
  
  
  # optimal output of lambda that was recommended by the transformation on the 
  # residuals and fitted values was -1.
  
  y_box <- (y^lambda - 1)/lambda
  
 # Re-fitting the model with y_box as the new response, and the same predictors.
  
  model_box <- lm(y_box ~ x)
  
  
  # Compare this with the MSE of the transformed model. 
  
  new_fitted_values <- (lambda * model_box$fitted.values + 1)^(1/lambda)
  
  # These new_fitted_values are directly comparable to the original 
  # model$fitted.values (without the transformation). 
  # So the MSE of the transformed model is

  mean((y - new_fitted_values)^2)
  
  # and we displayed this alongside the original MSE to see which one is lower 
  
  mean((y - model$fitted.values)^2)
  
  # Now we will ignore the new_fitted_values in the previous step 
  # and run the usual diagnostic tests on the model_box with its output,
  
  m1<-model_box$fitted.values 
  
  m2<-model_box$residuals
  
  print("Q-Q plot for orignial Model:Right")
  op<-par(pty="s", mfrow=c(1,2))
  qqnorm(model$residuals)
  qqline(model$residuals, col="red")
  
  print("Q-Q plot for Box-Cox transformed model:Left")
  qqnorm(model_box$residuals)
  qqline(model_box$residuals, col="red")
  
  par(op)
  print("Data appears to be linear.")
}







