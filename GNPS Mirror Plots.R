list.of.packages <- c("officer", "tidyverse", "RSelenium","netstat","magrittr","flextable")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(officer)
library(RSelenium)
library(tidyverse)
library(netstat)
library(magrittr)
library(flextable)

#GNPS URL
url <- 'https://gnps.ucsd.edu/ProteoSAFe/status.jsp?task=a41e05f0699e45328675aac49f667753'

#set working directory
setwd("~/Downloads/")

#Load Matches
Matches <- read_csv("~/Downloads/Liver_Pos_Matches(Done).csv")

#Output presentation
output <- "Liver_Pos_Matches(Done).pptx"

############
#If there are any errors, try incrementing by 1
loadingtime <- 1
############

Matches %<>% drop_na(Annotation)
Matches %<>% mutate('Shared Peaks'= Matches$'No..of.Shared.Peaks')
Matches %<>% mutate('Cosine Score'= Matches$'Cosine.Score')
MatchesAnalog <- Matches[(Matches$`Analog?`=="Yes"),]
MatchesAnalog %<>% select(Annotation,'M/Z',RT,'Cosine Score','Shared Peaks',`Cluster Index`)
MatchesAnnotation <- Matches[(Matches$`Analog?`=="No"),]
MatchesAnnotation %<>% select(Annotation,'M/Z',RT,'Cosine Score','Shared Peaks',`Cluster Index`)

#Make Images folder
unlink("~/Downloads/GNPS_Mirror_Matches_Temp",recursive = TRUE)
dir.create("~/Downloads/GNPS_Mirror_Matches_Temp")
file.create("~/Downloads/GNPS_Mirror_Matches_Temp/null.png")

#Make Powerpoint
pres <- read_pptx()
pres <- add_slide(pres,layout="Title Slide",master = "Office Theme")
pres <- ph_with(pres, value = "Mirror Matches",location=ph_location_type(type="ctrTitle"))
pres <- ph_with(pres, value = paste(nrow(MatchesAnnotation)," Annotations + ",nrow(MatchesAnalog)," Analogs"),location=ph_location_type(type="subTitle"))

pres <- add_slide(pres,layout="Title Slide",master = "Office Theme")
pres <- ph_with(pres, value = "Annotations",location=ph_location_type(type="ctrTitle"))

#Set Firefox profile
fprof <- makeFirefoxProfile(list(browser.download.manager.showWhenStarting=FALSE,
                                 browser.download.dir = "~/Downloads/GNPS_Mirror_Matches_Temp", 
                                 browser.helperApps.alwaysAsk.force= FALSE,
                                 browser.helperApps.neverAsk.saveToDisk="image/png",
                                 browser.helperApps.neverAsk.openFile = "image/png",
                                 browser.download.folderList = 2L))

#Create Remote Driver
rD <- rsDriver(port=free_port(),browser = "firefox",extraCapabilities=fprof)

remDr <- rD[["client"]]

remDr$navigate(url)
Sys.sleep(loadingtime)

remDr$findElement(value='/html/body/div/div[2]/table/tbody/tr[2]/td/table/tbody/tr[2]/td/div[4]/a[1]')$clickElement()
Sys.sleep(loadingtime+1)

x <- 1
delay <- 0

#Get image from GNPS function
GetImage <- function() {
  remDr$findElement(using="id",value='main.#Scan#_lowerinput')$clearElement()
  remDr$findElement(using="id",value='main.#Scan#_upperinput')$clearElement()
  
  remDr$findElement(using="id",value='main.#Scan#_lowerinput')$sendKeysToElement(list(clusterIndex))
  remDr$findElement(using="id",value='main.#Scan#_upperinput')$sendKeysToElement(list(clusterIndex))
  
  Sys.sleep(loadingtime)
  
  remDr$findElement(using="id",value="main.filter")$clickElement()
  Sys.sleep(loadingtime)
  remDr$findElement(value='/html/body/div[3]/div[1]/table/tbody/tr[1]/td[1]/button[1]')$clickElement()
  Sys.sleep(loadingtime)
  while(remDr$findElement(value='/html/body/div[3]/div[1]/table/tbody/tr[3]/td[2]/div/div/table/tbody/tr[3]/td/div[2]/nobr[1]/input[4]')$isElementDisplayed()==FALSE){
    Sys.sleep(0.1)
    delay = delay + 0.1
    
    if(delay>5){
      remDr$findElement(using="id",value="main.filter")$clickElement()
      Sys.sleep(loadingtime)
      delay <- 0
    }
  }
  Sys.sleep(loadingtime)
  
  remDr$findElement(value='/html/body/div[3]/div[1]/table/tbody/tr[3]/td[2]/div/div/table/tbody/tr[3]/td/div[2]/nobr[1]/input[4]')$getElementLocationInView()
  Sys.sleep(loadingtime)
  
  remDr$findElement(value='/html/body/div[3]/div[1]/table/tbody/tr[3]/td[2]/div/div/table/tbody/tr[3]/td/div[2]/nobr[1]/input[4]')$clickElement()
  Sys.sleep(loadingtime)
}
#
#webElemtest <-NULL
#while(is.null(webElemtest)){
#  webElemtest <- tryCatch({remDr$findElement(value='/html/body/div[3]/div[1]/table/tbody/tr[3]/td[2]/div/div/table/tbody/tr[3]/td/div[2]/nobr[1]/input[4]')},
#                          error = function(e){NULL})
#}


#Add slide to powerpoint function
AddSlide <- function() {
  pres <- add_slide(pres,layout="Title and Content",master = "Office Theme")
  
  img_loc <- ph_location(left=0.5,top=1.5,width=9,height=5.7)
  pres <- ph_with(pres, value = external_img(paste("~/Downloads/GNPS_Mirror_Matches_Temp/null(",x,").png",sep="")),location=img_loc)
  
  subset_x <- MatchesAnnotation[x,]
  
  table_loc <- ph_location(left=1,top=0.4,width=10,height=1)
  
  a <- width(flextable(subset_x),j=1,width=4.63)
  a %<>% bold(i=1)
  a %<>% theme_zebra(odd_header = "#5A80B8",odd_body = "#D1D8E6")
  a %<>% color(part="header",color="#FFFFFF")
  a %<>% hline_top(border = fp_border(color="white",width=2))
  a %<>% vline(border = fp_border(color="white",width=1))
  a %<>% hrule(rule="atleast",part="all")
  a %<>% height(height=0.5,part="body")
  a %<>% colformat_num(j=6,big.mark="")
  a %<>% colformat_double(j=4,digits=4)
  
  pres <- ph_with(pres, value = a,location= table_loc)
}


#Loop for Annotations
repeat{
  
  if(dim(MatchesAnnotation)[1]==0){
    break
  }
  
  clusterIndex <- toString(MatchesAnnotation$`Cluster Index`[x])
  
  GetImage()
  
  AddSlide()
  
  x <- x+1
  
  if(x>nrow(MatchesAnnotation)){
    x <- 1
    break
  }
}


########
#Analogs
########

#Make Images folder
unlink("~/Downloads/GNPS_Mirror_Matches_Temp",recursive = TRUE)
dir.create("~/Downloads/GNPS_Mirror_Matches_Temp")
file.create("~/Downloads/GNPS_Mirror_Matches_Temp/null.png")

x <- 1
delay <- 1

pres <- add_slide(pres,layout="Title Slide",master = "Office Theme")
pres <- ph_with(pres, value = "Analogs",location=ph_location_type(type="ctrTitle"))

#Get image from GNPS function
GetImageAnalog <- function() {
  remDr$findElement(using="id",value='main.#Scan#_lowerinput')$clearElement()
  remDr$findElement(using="id",value='main.#Scan#_upperinput')$clearElement()
  
  remDr$findElement(using="id",value='main.#Scan#_lowerinput')$sendKeysToElement(list(clusterIndex))
  remDr$findElement(using="id",value='main.#Scan#_upperinput')$sendKeysToElement(list(clusterIndex))
  
  Sys.sleep(loadingtime)
  
  remDr$findElement(using="id",value="main.filter")$clickElement()
  Sys.sleep(loadingtime)
  remDr$findElement(value='/html/body/div[3]/div[1]/table/tbody/tr[1]/td[1]/img')$clickElement()
  Sys.sleep(loadingtime)
  while(remDr$findElement(value='/html/body/div[3]/div[1]/table/tbody/tr[2]/td[2]/div/div/table/tbody/tr[3]/td/div[2]/nobr[1]/input[4]')$isElementDisplayed()==FALSE){
    Sys.sleep(0.1)
    delay = delay + 0.1
    
    if(delay>5){
      remDr$findElement(using="id",value="main.filter")$clickElement()
      Sys.sleep(loadingtime)
      delay <- 0
    }
  }
  Sys.sleep(loadingtime)
  
  remDr$findElement(value='/html/body/div[3]/div[1]/table/tbody/tr[2]/td[2]/div/div/table/tbody/tr[3]/td/div[2]/nobr[1]/input[4]')$getElementLocationInView()
  Sys.sleep(loadingtime)
  remDr$findElement(value='/html/body/div[3]/div[1]/table/tbody/tr[2]/td[2]/div/div/table/tbody/tr[3]/td/div[2]/nobr[1]/input[4]')$clickElement()
  Sys.sleep(loadingtime)
} 

#Add slide to powerpoint function
AddSlideAnalog <- function() {
  pres <- add_slide(pres,layout="Title and Content",master = "Office Theme")
  
  img_loc <- ph_location(left=0.5,top=1.5,width=9,height=5.7)
  pres <- ph_with(pres, value = external_img(paste("~/Downloads/GNPS_Mirror_Matches_Temp/null(",x,").png",sep="")),location=img_loc)
  
  subset_x <- MatchesAnalog[x,]
  
  table_loc <- ph_location(left=1,top=0.4,width=10,height=1)
  
  a <- width(flextable(subset_x),j=1,width=4.63)
  a %<>% bold(i=1)
  a %<>% theme_zebra(odd_header = "#5A80B8",odd_body = "#D1D8E6")
  a %<>% color(part="header",color="#FFFFFF")
  a %<>% hline_top(border = fp_border(color="white",width=2))
  a %<>% vline(border = fp_border(color="white",width=1))
  a %<>% hrule(rule="atleast",part="all")
  a %<>% height(height=0.5,part="body")
  a %<>% colformat_num(j=6,big.mark="")
  a %<>% colformat_double(j=4,digits=4)
  
  pres <- ph_with(pres, value = a,location= table_loc)
}

remDr$navigate(url)
Sys.sleep(loadingtime)

remDr$findElement(value='/html/body/div/div[2]/table/tbody/tr[2]/td/table/tbody/tr[2]/td/div[4]/a[3]')$clickElement()
Sys.sleep(loadingtime+1)

#Loop for Analogs
repeat{
  
  if(dim(MatchesAnalog)[1]==0){
    break
  }
  
  clusterIndex <- toString(MatchesAnalog$`Cluster Index`[x])
  
  GetImageAnalog()
  
  AddSlideAnalog()
  
  x <- x+1
  
  if(x>nrow(MatchesAnalog)){
    x <- 1
    break
  }
}

print(pres,target=output)

unlink("~/Downloads/GNPS_Mirror_Matches_Temp",recursive = TRUE)

remDr$close()
rD$server$stop()
rm(rD)
gc()

