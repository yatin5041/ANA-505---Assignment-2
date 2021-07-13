# Traditional Conjoint Analysis (R)

# R preliminaries to get the user-defined function for spine chart: 
# place the spine chart code file <R_utility_program_1.R>
# in your working directory and execute it by
#     source("R_utility_program_1.R")
# Or if you have the R binary file in your working directory, use
#     load(file="mtpa_spine_chart.Rdata")

# spine chart accommodates up to 45 part-worths on one page
# |part-worth| <= 40 can be plotted directly on the spine chart
# |part-worths| > 40 can be accommodated through standardization

print.digits <- 2  # set number of digits on print and spine chart

library(support.CEs)  # package for survey construction ## The package support.CEs provides basic functions that support an implementation of choice experiments

# generate a balanced set of product profiles for survey ## Line 21-31 organizing data using Lma.design function having format Lma.design(candidate.array = NULL, attribute.names,
## nalternatives, nblocks, row.renames = TRUE, seed = NULL)

provider.survey <- Lma.design(attribute.names = 
  list(brand = c("AT&T","T-Mobile","US Cellular","Verizon"), ## attribute.names created using list() using variables brand,startup,monthly,service,reatil,apple,samsung,google
  startup = c("$100","$200","$300","$400"), ## c() combines values in the list
  monthly = c("$100","$200","$300","$400"),
  service = c("4G NO","4G YES"), 
  retail = c("Retail NO","Retail YES"),
  apple = c("Apple NO","Apple YES"), 
  samsung = c("Samsung NO","Samsung YES"), 
  google = c("Nexus NO","Nexus YES")), nalternatives = 1, nblocks=1, seed=9999) ##nalternatives is set as 1 for a binary choice experiment with an opt-out or common base option
## nblocks is an integer value describing the number of blocks into which a choice experiment design is divided, in this case 1, Seed is random number generator
print(questionnaire(provider.survey))  # print survey design for review ## function questionnaire () convertes CE design into CE questions used in a questionnaire survey.

sink("questions_for_survey.txt")  # send survey to external text file ## seurvey is saved in external text file
questionnaire(provider.survey)
sink() # send output back to the screen ##output displayed on screen

# user-defined function for plotting descriptive attribute names ## created function effect.name.map and populated based on function effect.name
effect.name.map <- function(effect.name) { 
  if(effect.name=="brand") return("Mobile Service Provider")
  if(effect.name=="startup") return("Start-up Cost")
  if(effect.name=="monthly") return("Monthly Cost")
  if(effect.name=="service") return("Offers 4G Service")
  if(effect.name=="retail") return("Has Nearby Retail Store")
  if(effect.name=="apple") return("Sells Apple Products")
  if(effect.name=="samsung") return("Sells Samsung Products")
  if(effect.name=="google") return("Sells Google/Nexus Products")
  } 

# read in conjoint survey profiles with respondent ranks ## created conjoint data frame by using mobile_services_ranking.csv file
conjoint.data.frame <- read.csv("mobile_services_ranking.csv")

# set up sum contrasts for effects coding as needed for conjoint analysis ##options() function examines a variety of global options which affect the way in which R computes and displays its results
options(contrasts=c("contr.sum","contr.poly"))  ## a logical indicating whether contrasts should be computed. Here sum contarsts is used.

# main effects model specification ## ranking is column in the dataset populated by author after going through given portfolio options
main.effects.model <- {ranking ~ brand + startup + monthly + service + 
  retail + apple + samsung + google}

# fit linear regression model using main effects only (no interaction terms) ## im is linear model used as lm(formula, data)
main.effects.model.fit <- lm(main.effects.model, data=conjoint.data.frame)
print(summary(main.effects.model.fit)) ## printed summary of the model

# save key list elements of the fitted model as needed for conjoint measures
conjoint.results <- 
  main.effects.model.fit[c("contrasts","xlevels","coefficients")] ## combined values in of contrasts, xlevels and coefficients, as we are interested in those

conjoint.results$attributes <- names(conjoint.results$contrasts) ## $ operator is used to extract list or subset a specific part of a data object, here contrasts

# compute and store part-worths in the conjoint.results list structure ## this section of the code is computing and storing most important aspect part worths
part.worths <- conjoint.results$xlevels  # list of same structure as xlevels
end.index.for.coefficient <- 1  # intitialize skipping the intercept
part.worth.vector <- NULL # used for accumulation of part worths ## line 73-85 used "for loop", seq() creates sequence for elements in vector
for(index.for.attribute in seq(along=conjoint.results$contrasts)) {
  nlevels <- length(unlist(conjoint.results$xlevels[index.for.attribute])) ## unlist converts list to vector
  begin.index.for.coefficient <- end.index.for.coefficient + 1
  end.index.for.coefficient <- begin.index.for.coefficient + nlevels -2
  last.part.worth <- -sum(conjoint.results$coefficients[
    begin.index.for.coefficient:end.index.for.coefficient])
  part.worths[index.for.attribute] <- 
    list(as.numeric(c(conjoint.results$coefficients[
      begin.index.for.coefficient:end.index.for.coefficient],
      last.part.worth)))
  part.worth.vector <- 
    c(part.worth.vector,unlist(part.worths[index.for.attribute]))    
  } 
conjoint.results$part.worths <- part.worths

# compute standardized part-worths ## these values are used to plot spine chart
standardize <- function(x) {(x - mean(x)) / sd(x)}
conjoint.results$standardized.part.worths <- 
  lapply(conjoint.results$part.worths,standardize) ## lapply acts on a list and returns a list of the same length as input list object
 
# compute and store part-worth ranges for each attribute 
part.worth.ranges <- conjoint.results$contrasts
for(index.for.attribute in seq(along=conjoint.results$contrasts)) ## used for loop
  part.worth.ranges[index.for.attribute] <- 
  dist(range(conjoint.results$part.worths[index.for.attribute])) ## Range() function returns the maximum and minimum value of the vector and column of the dataframe
conjoint.results$part.worth.ranges <- part.worth.ranges          ## dist () computes and returns the distance matrix

sum.part.worth.ranges <- sum(as.numeric(conjoint.results$part.worth.ranges))

# compute and store importance values for each attribute 
attribute.importance <- conjoint.results$contrasts
for(index.for.attribute in seq(along=conjoint.results$contrasts)) 
  attribute.importance[index.for.attribute] <- 
  (dist(range(conjoint.results$part.worths[index.for.attribute]))/
  sum.part.worth.ranges) * 100
conjoint.results$attribute.importance <- attribute.importance
 
# data frame for ordering attribute names ## names() is used to get or set the name of an Object
attribute.name <- names(conjoint.results$contrasts)
attribute.importance <- as.numeric(attribute.importance) ## as.numeric() method takes an object that needs to be coerced and returns the converted numeric value
temp.frame <- data.frame(attribute.name,attribute.importance)
conjoint.results$ordered.attributes <- 
  as.character(temp.frame[sort.list(                      ## as.character() attempts to coerce its argument to character type
  temp.frame$attribute.importance,decreasing = TRUE),"attribute.name"])

# respondent internal consistency added to list structure ## $ operator is used to extract list or subset a specific part of a data object, here r.squared
conjoint.results$internal.consistency <- summary(main.effects.model.fit)$r.squared 
 
# user-defined function for printing conjoint measures ## sprintf("%1.2f") dsiplayes all digits before decimal point and 2 digits after decimal point
if (print.digits == 2) 
  pretty.print <- function(x) {sprintf("%1.2f",round(x,digits = 2))} 
if (print.digits == 3) 
  pretty.print <- function(x) {sprintf("%1.3f",round(x,digits = 3))} 
 
# report conjoint measures to console 
# use pretty.print to provide nicely formated output ## "\n" to print on new line and cat() combines character values and print them to the screen or directly save in a file
for(k in seq(along=conjoint.results$ordered.attributes)) {
  cat("\n","\n")
  cat(conjoint.results$ordered.attributes[k],"Levels: ",
  unlist(conjoint.results$xlevels[conjoint.results$ordered.attributes[k]]))
  
  cat("\n"," Part-Worths:  ")
  cat(pretty.print(unlist(conjoint.results$part.worths
    [conjoint.results$ordered.attributes[k]])))
    
  cat("\n"," Standardized Part-Worths:  ")
  cat(pretty.print(unlist(conjoint.results$standardized.part.worths
    [conjoint.results$ordered.attributes[k]])))  
    
  cat("\n"," Attribute Importance:  ")
  cat(pretty.print(unlist(conjoint.results$attribute.importance
    [conjoint.results$ordered.attributes[k]])))
  }

# plotting of spine chart begins here
# all graphical output is routed to external pdf file ## conjoint.results is used to plot spine chart
pdf(file = "fig_preference_mobile_services_results.pdf", width=8.5, height=11)
spine.chart(conjoint.results)
dev.off()  # close the graphics output device

# Suggestions for the student:
# Enter your own rankings for the product profiles and generate
# conjoint measures of attribute importance and level part-worths.
# Note that the model fit to the data is a linear main-effects model.
# See if you can build a model with interaction effects for service
# provider attributes.

