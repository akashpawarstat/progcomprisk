#' SurviCompet
#'
#' 'SurviCompet' function from 'progcomprisk' enables the user to obtain DEGs
#' (Differentially expressed genes) which leads to death given it leads to
#' progression of cancer as well as accounts for competing risks.

#' @details SurviCompet function first adds a new column 'Status' to the input
#' dataset and assigns status 0, 1, 2 to each cell. It assigns 0 when progression is
#' 1, death is 0 or progression 0 death is 0, assigns 1 when progression is 1 and
#' death is 1, whereas assigns 2 when progression is 0 and death is 1.
#'
#' Further, it subsets 2 data sets, one data set named 'deathdata' which includes Ids/
#' Subjects with status 0 and 1 and applies Cox PH on it. Another subset data named as
#' 'compdata' includes Ids/Subjects with status 0 and 2, then applies Cox PH after
#' substituting 2 by 1. Then it filters out DEGs having P-value > 0.05 from both
#' subset data. At last the it merges the commom DEGs from both data and creates a
#' new dataframe which consists 'Id','OS','Death','PFS','Prog','Status' and observations
#' of common significant DEGs(leading to death given they leads to progression of cancer
#' as well as accounts for competing risks)
#' function it lists out the common DEGs from
#'
#'
#' SurviCompet(m,n,data), m - starting
#' column number of data from where DEG observations are to be read,    n - ending
#' column number of data till which DEG observations are to be read,    data - name of
#' input data consisting survival and progression outcomes with DEG observations.
#'
#' Concern should be given to rename the columns of the input data before using 'progcomprisk' package, as follows.
#'
#' 1) Subject Id column should be named as 'Id'.
#'
#' 2) OS column must be named as 'OS' .
#'
#' 3) Death status/event column should be named as 'Death'.
#'
#' 4) PFS column should be named as 'PFS'.
#'
#' 5) Progression event column should be named as 'Prog'.
#'

#' @import survival
#' @export
#' @param m numeric staring column of variables
#' @param n numeric ending column of variables
#' @param data Data frame by user
#' @param Status Column states 0, 1, 2 on the basis of Prog column and Death column

#' @examples
#' SurviCompet(6,104,user_data)
#' @author Akash Pawar

SurviCompet <- function(m,n,data){

  data$Status[data$Prog == 0 & data$Death == 1] <- 2
  data$Status[data$Prog == 0 & data$Death == 0] <- 0
  data$Status[data$Prog == 1 & data$Death == 1] <- 1
  data$Status[data$Prog == 1 & data$Death == 0] <- 0

  data1 <- data
  data2 <- data

  data2$Status[data2$Status == 1] <- 3 #death1 = 1 in mydata2 is assigned "3"
  data2$Status[data2$Status == 2] <- 1 #death1 = 2 in mydata2 is assigned "1"

  deathdata <- subset(data1, Status != 2) #deleting "2" row from mydata1
  compdata <- subset(data2, Status != 3) #deleting "3" row from mydata2

  da<-matrix(nrow = 0, ncol = 5) #creating dummy matrix
  dc<-matrix(nrow = 0, ncol = 5)

  for(i in m:n){
    b = deathdata[,i]
    c = compdata[,i]
    coxdeath <- coxph(Surv(OS, Status) ~ b, data = deathdata)
    coxcomp <- coxph(Surv(OS, Status) ~ c, data = compdata)
    dsum <- summary(coxdeath)
    csum <- summary(coxcomp)

    da <- rbind(da,dsum$coefficients)
    dc <- rbind(dc,csum$coefficients)
  }
  coefficientsA <- data.frame(da)
  coefficientsB <- data.frame(dc)
  namevect <- names(data)
  coefficientsA$variable = namevect[m:n]
  coefficientsB$variable = namevect[m:n]
  significantpvalueA <- subset(coefficientsA, Pr...z.. < 0.05)
  significantpvalueB <- subset(coefficientsB, Pr...z.. < 0.05)
  print("Results for data with death status 0 and 1, coefficientsA")
  print(coefficientsA)
  print("Results for data with death status 0 and 2, coefficientsB")
  print(coefficientsB)
  print("Significant genes from the data consisting progressions, significantpvalueA")
  print(significantpvalueA)
  print("Significant genes from the data which does not consist progressions, significantpvalueB")
  print(significantpvalueB)
  commongenes <- merge(significantpvalueA, significantpvalueB, by ="variable")
  print("Genes/variables playing role in competing risk as well due to progression, commongenes")
  print(commongenes)

  cvar <- commongenes$variable
  print("Following are the variables/genes found to be common")
  print(cvar)
  commondata <- subset(data, select = c("Id", "OS", "Death", "PFS", "Prog", "Status", cvar))
  print("Following is the data which is subset of the original data showing the ID, OS, Death, PFS, Prog, Status & Significant common variables")
  print(commondata)
}
