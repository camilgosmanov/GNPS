list.of.packages <- c("dplyr","tidyverse","magrittr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(tidyverse)
library(dplyr)
library(magrittr)

setwd("~/Downloads/")
Nodes <- read_csv("~/Downloads/Nodes.csv") #Node table from Cytoscape
Edge1 <- read_csv("~/Downloads/Edge.csv") #Edge table from Cytoscape
Masses <- read_csv("~/Downloads/Matches.csv") #Fold Change/Features table

output <- "~/Downloads/Matches(Completed).csv" #Set path for final file


#Splits column into M/Z and RT

Masses$features2 <- Masses$features
Masses %<>% separate(features2,into= c('M/Z','RT'),sep="_",convert = TRUE)
Masses %<>% mutate(`M/Z`=gsub("X","",`M/Z`))
Masses$`M/Z` <- as.numeric(Masses$`M/Z`)
Masses$features2 <- NULL

#Add other columns 

Masses %<>% add_column('Annotation'=NA,'Cosine.Score'=NA,'No..of.Shared.Peaks'=NA,
                       'Mass.Diff..to.Library.Reference'=NA,'ppm.error'=NA,'Cluster Index'=NA,
                       'Analog?'=NA,'Nearest Neighbors'=NA,'Proximity'=NA)


Masses %<>% mutate("Mass" = Masses$`M/Z`)
Nodes %<>% mutate("Mass" = Nodes$`precursor mass`)

###########
#Thresholds

Masses %<>% mutate("MassHigh"= Mass+ 0.0011) #Thresholds for Mass 
Masses %<>% mutate("MassLow"= Mass-0.0001)

Masses %<>% mutate("RTHigh"= RT+0.011) #Thresholds for RT 
Masses %<>% mutate("RTLow"= RT-0.001)

Masses %<>% mutate("Threshold_Increased" = "No")
numThresholds <- 0


Edge1 %<>% select("node1","node2","EdgeAnnotation","property1")

Edge2 <- Edge1
Edge2 %<>% mutate(store=0)
Edge2$store <- Edge2$node1 
Edge2$node1 <- Edge2$node2
Edge2$node2 <- Edge2$store
Edge2$store <- NULL

Edge <- rbind(Edge1,Edge2)

x <- 1

repeat {#Match masses with node cluster index
  indicies <- which(((Nodes$Mass <= Masses$MassHigh[x]) & (Nodes$Mass >= Masses$MassLow[x])) & 
  ((Nodes$RTMean <= Masses$RTHigh[x]) & (Nodes$RTMean >= Masses$RTLow[x])))
  subset_Nodes <- Nodes[indicies,] 
  
  if(dim(subset_Nodes)[1]==0){#If no nodes match the mass
    
    Masses$MassHigh[x] <- Masses$MassHigh[x] + 0.0001 #increase thresholds by amount
    Masses$MassLow[x] <- Masses$MassLow[x] - 0.0001
    Masses$RTHigh[x] <- Masses$RTHigh[x] + 0.001
    Masses$RTLow[x] <- Masses$RTLow[x] - 0.001
    
    numThresholds <- numThresholds + 1 #number of times thresholds were increased
    
    if(numThresholds==1){ #New column shows number of threshold increases
      Masses$Threshold_Increased[x] <- paste("Yes ", numThresholds," time")
      } else{Masses$Threshold_Increased[x] <- paste("Yes ", numThresholds," times")}
    
  } else{
  
  Masses$Annotation[x] <- subset_Nodes$Compound_Name[1] #set annotation as compound name
  Masses$`Cluster Index`[x] <- subset_Nodes$'cluster index'[1] #set cluster index equal to found node
  
  if(is.na(Masses$Annotation[x])){ #if no compound name
    Masses$Annotation[x] <- subset_Nodes$`Analog:Compound_Name`[1] #set to analog name
    if(is.na(Masses$Annotation[x])){ #if no analog name
      Masses$`Analog?`[x] = "No" #no annotation
    } else{ #annotation is analog
      
      Masses$`Analog?`[x] = "Yes"
      
    #copy over info from found node
    Masses$'Mass.Diff..to.Library.Reference'[x] <- subset_Nodes$'Analog:MassDiff'[1]
    Masses$'No..of.Shared.Peaks'[x] <- subset_Nodes$`Analog:SharedPeaks`[1]
    Masses$'ppm.error'[x] <- subset_Nodes$`Analog:MZErrorPPM`[1]
    Masses$'Cosine.Score'[x] <- subset_Nodes$'Analog:MQScore'[1]
    }
  } else { #annotation is compound name
    Masses$`Analog?`[x] = "No"
    
    Masses$'Mass.Diff..to.Library.Reference'[x] <- subset_Nodes$'MassDiff'[1]
    Masses$'No..of.Shared.Peaks'[x] <- subset_Nodes$`SharedPeaks`[1]
    Masses$'ppm.error'[x] <- subset_Nodes$`MZErrorPPM`[1]
    Masses$'Cosine.Score'[x] <- subset_Nodes$'MQScore'[1]
    
  }
  
  x <- x+1} #node match was completed
  
  if(x>nrow(Masses))
  {
    break
  }
}
  

########
#Nearest Neighbors
########

x <- 1
i <- 1

Masses$Proximity <- "None"

repeat{
  indicies <- which(Edge$node1 == Masses$`Cluster Index`[x]) #find all first neighbors
  subset_Edge <- Edge[indicies,] #create subset_Edge
  
  if(subset_Edge$property1[1] == 1){ #find all Self-Loops
    Masses$`Nearest Neighbors`[x] <- "Self-Loop"
    Masses$Proximity[x] <- "Self-Loop"
    
  } else { #if not Self-Loop
    
    repeat{ #find first neighbor names and analog names
      node_match <- which(Nodes$`cluster index` == subset_Edge$node2[i])
      
      if(!is.na(Nodes$Compound_Name[node_match])){
        Masses$`Nearest Neighbors`[x] <- Nodes$Compound_Name[node_match]
        Masses$Proximity[x] <- "First"
        
        i <- 1
        break
        
      } else if(!is.na(Nodes$`Analog:Compound_Name`[node_match])){
        
        Masses$`Nearest Neighbors`[x] <- Nodes$`Analog:Compound_Name`[node_match]
        Masses$Proximity[x] <- "First Analog"
        
      } 
      
      i <- i + 1
      
      if(i>nrow(subset_Edge)){
        i <- 1
        break}
    }
  }
  
  x <- x + 1
  
  if(x>nrow(Masses)){
    x <- 1
    break
  }
}


x <- 1
i <- 1
i2 <- 1
z <- 1


repeat{#Second Neighbors
  if(!(Masses$Proximity[x] == "First") & !(Masses$Proximity[x] == "Self-Loop")){     
    indicies <- which(Edge$node1 == Masses$`Cluster Index`[x]) #find all first neighbors
    subset_Edge <- Edge[indicies,] #create subset_Edge 
    
    repeat{
      indicies2 <- which(Edge$node1 == subset_Edge$node2[i])
      subset_Edge2 <- Edge[indicies2,] #find all first neighbors of first neighbor
      
      
      repeat{ #find second neighbors and analogs
        if(Masses$`Cluster Index`[x] != subset_Edge2$node2[i2]){
          node_match <- which(Nodes$`cluster index` == subset_Edge2$node2[i2])
          
          if(!is.na(Nodes$Compound_Name[node_match])){
            Masses$`Nearest Neighbors`[x] <- Nodes$Compound_Name[node_match]
            Masses$Proximity[x] <- "Second"
            
            i2 <- 1
            z <- 2
            break
            
          } else if(!is.na(Nodes$`Analog:Compound_Name`[node_match])){
            
            if(Masses$Proximity[x] != "First Analog"){
              Masses$`Nearest Neighbors`[x] <- Nodes$`Analog:Compound_Name`[node_match]
              Masses$Proximity[x] <- "Second Analog" }
            
          } 
          
        } #end if
        
        i2 <- i2 + 1
        
        if(i2>nrow(subset_Edge2)){
          i2 <- 1
          break}
      }
      
      if(z == 2){
        z <- 1
        i <- 1
        break
      } else {i <- i + 1}
      
      if(i>nrow(subset_Edge)){
        i <- 1
        break
      }
      
    }
    
  } 
  
  x <- x + 1
  
  if(x>nrow(Masses)){
    x <- 1
    break
  }
}



x <- 1
i <- 1
i2 <- 1
i3 <- 1
z <- 1


repeat{#Third Neighbors
  if((Masses$Proximity[x] == "None") | (Masses$Proximity[x] == "Second Analog")){     
    indicies <- which(Edge$node1 == Masses$`Cluster Index`[x]) #find all first neighbors
    subset_Edge <- Edge[indicies,] #create subset_Edge 
    
    repeat{
      indicies2 <- which(Edge$node1 == subset_Edge$node2[i])
      subset_Edge2 <- Edge[indicies2,] #find all first neighbors of first neighbor (second neighbors)
      
      repeat{#Find all first neighbors of second neighbors (third neighbors)
        indicies3 <- which(Edge$node1 == subset_Edge2$node2[i2])
        subset_Edge3 <- Edge[indicies3,]
        
        ################
        #find third neighbors names
        
        repeat{ #find third neighbors and analogs
          if(Masses$`Cluster Index`[x] != subset_Edge3$node2[i3]){
            node_match <- which(Nodes$`cluster index` == subset_Edge3$node2[i3])
            
            if(!is.na(Nodes$Compound_Name[node_match])){
              Masses$`Nearest Neighbors`[x] <- Nodes$Compound_Name[node_match]
              Masses$Proximity[x] <- "Third"
              
              i3 <- 1
              z <- 2
              break
              
            } else if(!is.na(Nodes$`Analog:Compound_Name`[node_match])){
              
              if(Masses$Proximity[x] != "Second Analog"){
                Masses$`Nearest Neighbors`[x] <- Nodes$`Analog:Compound_Name`[node_match]
                Masses$Proximity[x] <- "Third Analog" }
              
            } 
            
          } #end if
          
          i3 <- i3 + 1
          
          if(i3>nrow(subset_Edge3)){
            i3 <- 1
            break}
        }
        
        #find third neighbors names
        ####################
        
        if(z == 2){
          i2 <- 1
          break
        } else{i2 <- i2 + 1}
        
        if(i2>nrow(subset_Edge2)){
          i2 <- 1
          break
        }
        
      }
      
      
      if(z == 2){
        z <- 1
        i <- 1
        break
      } else{i <- i + 1}
      
      if(i>nrow(subset_Edge)){
        i <- 1
        break
      }
      
    }
    
  } 
  
  x <- x + 1
  
  if(x>nrow(Masses)){
    x <- 1
    break
  }
}


x <- 1
i <- 1
i2 <- 1
i3 <- 1
i4 <- 1
z <- 1


repeat{#Fourth Neighbors
  if((Masses$Proximity[x] == "None") | (Masses$Proximity[x] == "Third Analog")){     
    indicies <- which(Edge$node1 == Masses$`Cluster Index`[x]) #find all first neighbors
    subset_Edge <- Edge[indicies,] #create subset_Edge 
    
    repeat{
      indicies2 <- which(Edge$node1 == subset_Edge$node2[i])
      subset_Edge2 <- Edge[indicies2,] #find all first neighbors of first neighbor (second neighbors)
      
      repeat{#Find all first neighbors of second neighbors (third neighbors)
        indicies3 <- which(Edge$node1 == subset_Edge2$node2[i2])
        subset_Edge3 <- Edge[indicies3,]
        
        repeat{#Find all fourth neighbors
          indicies4 <- which(Edge$node1 == subset_Edge3$node2[i3])
          subset_Edge4 <- Edge[indicies4,]
          
          ################
          #find fourth neighbors names
          
          repeat{ #find fourth neighbors and analogs
            if(Masses$`Cluster Index`[x] != subset_Edge4$node2[i4]){
              node_match <- which(Nodes$`cluster index` == subset_Edge4$node2[i4])
              
              if(!is.na(Nodes$Compound_Name[node_match])){
                Masses$`Nearest Neighbors`[x] <- Nodes$Compound_Name[node_match]
                Masses$Proximity[x] <- "Fourth"
                
                i4 <- 1
                z <- 2
                break
                
              } else if(!is.na(Nodes$`Analog:Compound_Name`[node_match])){
                
                if(Masses$Proximity[x] != "Third Analog"){
                  Masses$`Nearest Neighbors`[x] <- Nodes$`Analog:Compound_Name`[node_match]
                  Masses$Proximity[x] <- "Fourth Analog" }
                
              } 
              
            } #end if
            
            i4 <- i4 + 1
            
            if(i4>nrow(subset_Edge4)){
              i4 <- 1
              break}
          }
          
          #find fourth neighbors names
          ####################
          
          if(z == 2){
            i3 <- 1
            break
          } else{i3 <- i3 + 1}
          
          if(i3>nrow(subset_Edge3)){
            i3 <- 1
            break
          }
          
        }
        
        if(z == 2){
          i2 <- 1
          break
        } else{i2 <- i2 + 1}
        
        if(i2>nrow(subset_Edge2)){
          i2 <- 1
          break
        }
        
      }
      
      
      if(z == 2){
        z <- 1
        i <- 1
        break
      } else{i <- i + 1}
      
      if(i>nrow(subset_Edge)){
        i <- 1
        break
      }
      
    }
    
  } 
  
  x <- x + 1
  
  if(x>nrow(Masses)){
    x <- 1
    break
  }
}


#remove unnecessary columns
Masses$RTHigh <- NULL
Masses$RTLow <- NULL
Masses$MassHigh <- NULL
Masses$MassLow <- NULL
Masses$Mass <- NULL

if(length(unique(Masses$Threshold_Increased)) == 1){
  Masses$Threshold_Increased <- NULL
}

Masses %<>% arrange(Masses$`M/Z`)

write.csv(Masses,output,row.names=FALSE)                     


                

                      