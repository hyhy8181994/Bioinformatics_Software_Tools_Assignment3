#The question I proposed for Callyspongiidae is finding how the species richness of Callyspongiidae in Northern part of earth is compare to Southern part of earth. The Callyspongiidae is widely spread in world from North to South. 
#I assume latitude from 0 to 30 or -30 to 0 is Southern part, and latitude 30 to 90 or -90 to -30 as Northern part.  
#

library(tidyverse)
library(iNEXT)
library(vegan)

Callyspongiidae <- read_tsv("Callyspongiidae.txt")   #load file
View(Callyspongiidae)

hist(Callyspongiidae$lat)     #histogram of Callyspongiidae with latitude
Callyspongiidae$lat[which.max(Callyspongiidae$lat)] # showing the most Northern sample is collected 62.6628
Callyspongiidae$lat[which.min(Callyspongiidae$lat)] # 

Callys_South <- Callyspongiidae %>%   #create the subdata which were collected in southern part of each hemisphere
  filter(lat > -30) %>%   # latitude between -30 to 30 
  filter(lat < 30)        #


#funciton m_s to find the most southern sample collected place, 
m_s<-function(x){
  a<-NULL
  b<-NULL
  for (i in x) { #for loop  each latitude from beginning to end 
  if(i>0){
    a<-c(a,i)  #filter the latitude is over zero 
    }else{
    b<-c(b,i)  # latitude below zero
    }
  }
  print(a[which.min(a)]) #print the min latitude over zero
  print(b[which.max(b)]) #print the max latitude below zero
}

m_s(Callys_South$lat) #load the m_s function first#the most southern sample is collected in lat 7.555

Callys_North <- Callyspongiidae %>%  #create the subdata which were collected in northern part of each hemisphere
  filter(lat < -30 | lat > 30)  # latitude from -90 to -30 & 30 to 90

Callys_North$lat <- round(Callys_North$lat) #round the latitude to integer, narrow down the range
Callys_South$lat <- round(Callys_South$lat)

Callys_N_Bin <- Callys_North %>%  #turn north subdata in to community object with bin as proxy aganist latitude
  group_by(lat, bin_uri) %>%
  count(bin_uri)
  

Callys_S_Bin <- Callys_South %>%  #turn south subdata in to community object with bin as proxy aganist latitude
  group_by(lat, bin_uri) %>%
  count(bin_uri)


Callys_N_Bin <- spread(Callys_N_Bin, bin_uri, n)  #spread values across multiple columns
Callys_S_Bin <- spread(Callys_S_Bin, bin_uri, n)  #same above
Callys_N_Bin[is.na(Callys_N_Bin)]<-0  # switch na data to 0
Callys_S_Bin[is.na(Callys_S_Bin)]<-0  # same above

Callys_N_Bin <- Callys_N_Bin %>%   # remove row name
  remove_rownames %>%             
  column_to_rownames(var = "lat")#  set row name 

Callys_S_Bin <- Callys_S_Bin %>%  # # remove row name
  remove_rownames %>%          #
  column_to_rownames(var = "lat")  # set row name 


Callys_N_Bin<-as.data.frame(Callys_N_Bin) # change data type to data frame
Callys_S_Bin<-as.data.frame(Callys_S_Bin) 

accul_N <- specaccum(Callys_N_Bin)   # species acculumation curve of North part
accul_S <- specaccum(Callys_S_Bin)  # species acculumation curve of South part
accul_N
accul_S
plot(accul_N)   #plot the curve 
plot(accul_S)   #

modN<-iNEXT(Callys_N_Bin,q= 0, datatype = "abundance")   # Interpolation and extrapolation of Hill number with order q with North
modS<-iNEXT(Callys_S_Bin,q= 0, datatype = "abundance")  #  Interpolation and extrapolation of Hill number with order q with South

plot(modN) 
plot(modS)

specpool(Callys_N_Bin)  # estimate extrapolated species richness of north
specpool(Callys_S_Bin)  # estimate extrapolated species richness of sorth


Callyspongiidae_lat_s <- Callyspongiidae %>%    #arrange data with descending order
  group_by(lat) %>%
  arrange(desc(lat))

sum(is.na(Callyspongiidae_lat_s$lat)) #number of latitude with NA data


Callyspongiidae_lat_N<-Callyspongiidae_lat_s %>% #remove the sample without latitude data
  filter(!is.na(lat))   


#function to determine whether the sample is north or south, it create a list of S and N
search1<-function(x){
  n<-NULL
  for(i in x){    # loop for each latitude
    if (i > -30 & i< 30){   #if the latitude is between -30 and 30, label as S
      b<-c("S")
    } else {
      b<-c("N")  #otherwise label as N
    }
    n<-c(n,b)  #compose the vector contain N and S
    }
  return(n)
}

list_of_SN<-search1(Callyspongiidae_lat_N$lat) #load the search1 function first#perform the funciton
View(list_of_SN)


Callyspongiidae_lat_N$S_or_N<-list_of_SN   #add North or South information into data
View(Callyspongiidae_lat_N)


Callyspongiidae_NS_Bin <- Callyspongiidae_lat_N %>%   #turn data in to community object with bin as proxy aganist south or north
  group_by(S_or_N,bin_uri) %>%
  count(bin_uri)

View(Callyspongiidae_NS_Bin)

Callyspongiidae_NS_Bin$bin_uri[duplicated(Callyspongiidae_NS_Bin$bin_uri)] #show the sample appear in both North and South




#The Callyspongiidae species are widely spread from North to South. The most northern sample is collected in latitude 62.6628, and the most southern sample is collect in latitude 7.555. In samples with latitude data, there are 10 species appeared in both south and north ("BOLD:AAA4795" "BOLD:AAA4797" "BOLD:AAA6745" "BOLD:AAB0641" "BOLD:AAC5002" "BOLD:AAD0321" "BOLD:AAD4297" "BOLD:AAE7849" "BOLD:ACE9210" "BOLD:ACH8757).
#According to the species richness, the species richness in North is Higher than South. Maybe the species perfer cold than warm. 