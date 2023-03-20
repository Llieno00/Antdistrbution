
install.packages('dismo')  ##for my installation, get a question on install for 'terra' package and need to click "no" instead of "yes"
install.packages('rworldmap')
install.packages('sf')
install.packages('geodata')
install.packages('jsonlite')
install.packages("dplyr")
##every time you run it, need to activate the packages
library(dismo)
library(rworldmap)
library(sf)
library(geodata)
library(jsonlite)
library(raster)
library(ggplot2) # for data visualization
library(dplyr) # for data manipulation

wrld_simpl<-getMap(resolution = "coarse")
plot(wrld_simpl)
Fusca.gbif <- gbif("Formica","fusca", geo=TRUE)

##pulls out the lat and lon columns => there is more info if you want it, but we don't need it here
Fusca.coords<-cbind(Fusca.gbif$lon,Fusca.gbif$lat)
##delete rows with missing data = NA
Fusca.coords<-na.omit(Fusca.coords)
Fusca.coords<-data.frame(Fusca.coords)
colnames(Fusca.coords)<-c("lon","lat")


trim.coords<-function (x,latmin,latmax,lonmin,lonmax) {
  if (sum(x$lon < lonmin)>0) {
    tmp<-which(x$lon < lonmin)
    x<-x[-tmp,]}
  if (sum(x$lon > lonmax)>0) {
    tmp<-which(x$lon > lonmax)
    x<-x[-tmp,]}
  if (sum(x$lat < latmin)>0) {
    tmp<-which(x$lat < latmin)
    x<-x[-tmp,]}
  if (sum(x$lat > latmax)>0) {
    tmp<-which(x$lat > latmax)
    x<-x[-tmp,]}
  return(x) }

##Then use the function to make a new table of coordinates, just within the ranges specified
Fusca.coords.trim<-trim.coords(Fusca.coords,latmin=20,latmax=90,lonmin=-180,lonmax=180)


Fusca.coords<-Fusca.coords.trim




##download bioclimatic data from the worldclim database and convert to Raster format
bio.data<-worldclim_global(var="bio",res=10,path=getwd())
names(bio.data)<-paste0("bio",1:19)


plot(bio.data,1)
plot(bio.data,2)
plot(bio.data,15)
plot(bio.data,18)
plot(bio.data,19)

##you can extract bioclimatic data for the focal localities in Africa where your Species is found
bio.values <- extract(bio.data, Fusca.coords)[,-1]
rownames(bio.values)<-rownames(Fusca.coords)


par(mar = c(5, 5, 4, 2) + 0.1)
plot(bio.values[,1],bio.values[,2], xlab = "Annual Mean Temperature", ylab = "Mean Diurnal Range (Mean of monthly (max temp - min temp))")
plot(bio.values[,15],bio.values[,2], xlab = "Precipitation Seasonality (Coefficient of Variation)", ylab = "Mean Diurnal Range (Mean of monthly (max temp - min temp))")
plot(bio.values[,18],bio.values[,19], xlab = "Precipitation of Warmest Quarter", ylab = "Precipitation of Coldest Quarter")

##or look at how different variables correlate - here focusing on the first five
pairs(bio.values[,1:2])
pairs(bio.values[,c(15,2)])
pairs(bio.values[,c(18,19)])


##append to lat long and save to file
Fusca.data<-cbind(Fusca.coords,bio.values)
##remove any rows with missing data
Fusca.data<-na.omit(Fusca.data)
##write to file
write.csv(Fusca.data,file="Fusca.data.csv")

##extract mean values, min and max
rbind(mean=colMeans(Fusca.data),
      min=apply(Fusca.data, 2, min),
      max=apply(Fusca.data, 2, max))





##This bit makes a raster from the world map, which is a shaded region
ext <- extent(wrld_simpl)
xy <- abs(apply(as.matrix(bbox(ext)), 1, diff))
n <- 5
r <- raster(ext, ncol=xy[1]*n, nrow=xy[2]*n)
mask <-rasterize(wrld_simpl, r)


##I'm setting a plotting window using the same lat and lon as above
plot(mask, axes = TRUE, col = "grey",
     xlim = c(-180, 180),
     ylim = c(20,90)); box()
points(Fusca.coords,col = "red", pch = 20, cex = 1)



##Now I can set an area of extent for considering background using this new range
e <- extent(-180, 180, 20,90)
plot(e, add=TRUE, col='black')
##this generates 500 random points within that area
bg <- randomPoints(mask, 500, ext=e,extf=1)
points(bg,col="black", pch=20)
colnames(bg)<-c("lon","lat")


##It isn't essential but speeds up some later steps
bio.data<-crop(bio.data,e)


##b) Next step, combine the presence data and the background data in one table
train <- rbind(Fusca.coords, bg)
##and make a vector record 1=presence, 0=background
pb_train <- c(rep(1, nrow(Fusca.coords)), rep(0, nrow(bg))) 
##extract bioclimatic data for these points
envtrain <- extract(bio.data, train)
envtrain <- data.frame(cbind(pa=pb_train, envtrain) )
##and for each set separately
testpres <- data.frame( extract(bio.data, Fusca.coords) )
testbackg <- data.frame( extract(bio.data, bg) )



##	 To start with I am including the first 5 bioclim variables as predictors, but you could include 1,2 or more.I used 1,2,15,18,19 as before
gm1 <- glm(pa ~ bio1 + bio2 + bio15 + bio18 + bio19,
           family = binomial(link = "logit"), data=envtrain)

##look at a summary of the results
summary(gm1)


pg <- predict(bio.data, gm1, ext=e,type="response")
pg<-crop(pg,e)



##plot this
plot(pg, main='GLM probability of occurrence')
##add country boundaries
plot(wrld_simpl, add=TRUE, border='light blue') 
##add our observed locality data
points(Fusca.coords,col = "black", pch = 4, cex = 0.5)


##first we evaluate how well the model predicts presence/absence at each point
ge <- evaluate(testpres, testbackg, gm1)
ge


##then we use this evaluation to pick a threshold probability for defining presence/absence
##using the model that gives the most accurate match to observed presence/absence
tr <- threshold(ge, 'prevalence')
plot(pg > tr, main='presence/absence')
plot(wrld_simpl, add=TRUE, border='light blue') 
points(Fusca.coords,col = "black", pch = 4, cex = 0.5)



##This model uses next 5 climatic variables instead (I used 5 climatic values BO3,4,5,6,7)
gm2 <- glm(pa ~ bio3 + bio4 + bio5 + bio6 + bio7,family = binomial(link = "logit"), data=envtrain)
summary(gm2)

##You can compare two models, even with different predictor variables, using the Akaike Information Criterion
AIC(gm1,gm2)

##You can also compare the metrics of how well the model predicts the data
evaluate(testpres, testbackg, gm1)
evaluate(testpres, testbackg, gm2)


future.bio.data<-cmip6_world(model="CanESM5",var="bio",ssp="245",res=10,time="2061-2080",path=getwd())
names(future.bio.data)<-names(bio.data)


future.bio.data<-crop(future.bio.data,e)


par(mfrow=c(1,2))


#As before i changed the BO things to 1,2,15,18,19
gm1 <- glm(pa ~ bio1 + bio2 + bio15 + bio18 + bio19,
           family = binomial(link = "logit"), data=envtrain)

##and plot the predicted probability of occurrence at present
pg <- predict(bio.data, gm1, ext=e,type="response")
pg<-crop(pg,e)
plot(pg, main='A) GLM present')
plot(wrld_simpl, add=TRUE, border='light blue') 
points(Fusca.coords,col = "black", pch = 4, cex = 0.5)

## Now predict future distribution, i.e. feeding future bio data to model and plot this
pg.future <- predict(future.bio.data, gm1, ext=e,type="response")
pg.future<-crop(pg.future,e)
plot(pg.future, main="b) GLM, 2060-2081")
plot(wrld_simpl, add=TRUE, border='light blue') 
points(Fusca.coords,col = "black", pch = 4, cex = 0.5)

##Question 7: How is the distribution of the Fusca expected to change in the future?

##We can extract some numbers about how the range will change. For example
predict.localities.now <- extract(pg>=tr, Fusca.coords)[,-1]
predict.localities.future <- extract(pg.future>=tr, Fusca.coords)[,-1]
##number of localities where currently absent but present in the future, i.e. range expansion
matches.range.expansion <- predict.localities.now == FALSE & predict.localities.future == TRUE
sum(matches.range.expansion, na.rm = TRUE)

##number of localities where currently present but absent in the future, i.e. range contraction
matches.range.contraction <- predict.localities.now == TRUE & predict.localities.future == FALSE
sum(matches.range.contraction, na.rm = TRUE)

##Among all of the localities and background points, how many is it predicted to be in presently? 

predict.localities.now<-na.omit(predict.localities.now)
sum(predict.localities.now)


predict.localities.future<-na.omit(predict.localities.future)
sum(predict.localities.future)
#returned 14379

##Question 8: Use these numbers to calculate the percentage contraction or expansion in the range.

##A final plot that can be useful to understand causes of changes is to plot how climate will change

##work out the change in bioclim variables from now to the future
change.bio.data<-future.bio.data-bio.data




# Q2
par(mfrow=c(1,1))
Polyergus.gbif <- gbif("Polyergus", geo=TRUE)


##pulls out the lat and lon columns => there is more info if you want it, but we don't need it here
Polyergus.coords<-cbind(Polyergus.gbif$lon,Polyergus.gbif$lat)
##delete rows with missing data = NA
Polyergus.coords<-na.omit(Polyergus.coords)
Polyergus.coords<-data.frame(Polyergus.coords)
colnames(Polyergus.coords)<-c("lon","lat")

##Then use the function to make a new table of coordinates, just within the ranges specified
Polyergus.coords.trim<-trim.coords(Polyergus.coords,latmin=20,latmax=90,lonmin=-180,lonmax=180)

##plot world map again now using the trimmed data
plot(wrld_simpl, xlim=range(Polyergus.coords.trim$lon), ylim=range(Polyergus.coords.trim$lat), axes=TRUE, col="light yellow")
##add points for this Species
points(Polyergus.coords.trim, col='red', cex=0.75)

##if it looks OK, let's use the trimmed version as our new coordinates
Polyergus.coords<-Polyergus.coords.trim

##download bioclimatic data from the worldclim database and convert to Raster format
bio.data<-worldclim_global(var="bio",res=10,path=getwd())
names(bio.data)<-paste0("bio",1:19)



##you can extract bioclimatic data for the focal localities in Africa where your Species is found
bio.values <- extract(bio.data, Polyergus.coords)[,-1]
rownames(bio.values)<-rownames(Polyergus.coords)

##and plot them to see the range of environmental values the Species lives in
#Interesting combos: 1&2, 15&2, 18&19
par(mar = c(5, 5, 4, 2) + 0.1)
plot(bio.values[,1],bio.values[,2], xlab = "Annual Mean Temperature", ylab = "Mean Diurnal Range (Mean of monthly (max temp - min temp))")
plot(bio.values[,15],bio.values[,2], xlab = "Precipitation Seasonality (Coefficient of Variation)", ylab = "Mean Diurnal Range (Mean of monthly (max temp - min temp))")
plot(bio.values[,18],bio.values[,19], xlab = "Precipitation of Warmest Quarter", ylab = "Precipitation of Coldest Quarter")


##or look at how different variables correlate - here focusing on the first five
pairs(bio.values[,1:2])
pairs(bio.values[,c(15,2)])
pairs(bio.values[,c(18,19)])


##append to lat long and save to file
Polyergus.data<-cbind(Polyergus.coords,bio.values)
##remove any rows with missing data
Polyergus.data<-na.omit(Polyergus.data)
##write to file
write.csv(Polyergus.data,file="Polyergus.data.csv")

##extract mean values, min and max
rbind(mean=colMeans(Polyergus.data),
      min=apply(Polyergus.data, 2, min),
      max=apply(Polyergus.data, 2, max))


ext <- extent(wrld_simpl)
xy <- abs(apply(as.matrix(bbox(ext)), 1, diff))
n <- 5
r <- raster(ext, ncol=xy[1]*n, nrow=xy[2]*n)
mask <-rasterize(wrld_simpl, r)

##this plots it and shows your points
##I'm setting a plotting window using the same lat and lon as above
plot(mask, axes = TRUE, col = "grey",
     xlim = c(-180, 180),
     ylim = c(20,90)); box()
points(Polyergus.coords,col = "red", pch = 20, cex = 1)



##Now I can set an area of extent for considering background using this new range
e <- extent(-180, 180, 20,90)
plot(e, add=TRUE, col='black')
##this generates 500 random points within that area
bg <- randomPoints(mask, 500, ext=e,extf=1)
points(bg,col="black", pch=20)
colnames(bg)<-c("lon","lat")


bio.data<-crop(bio.data,e)


##b) Next step, combine the presence data and the background data in one table
train <- rbind(Polyergus.coords, bg)
##and make a vector record 1=presence, 0=background
pb_train <- c(rep(1, nrow(Polyergus.coords)), rep(0, nrow(bg))) 
##extract bioclimatic data for these points
envtrain <- extract(bio.data, train)
envtrain <- data.frame(cbind(pa=pb_train, envtrain) )
##and for each set separately
testpres <- data.frame( extract(bio.data, Polyergus.coords) )

testbackg <- data.frame( extract(bio.data, bg) )



##	 To start with I am including the first 5 bioclim variables as predictors, but you could include 1,2 or more.I used 1,2,15,18,19 as before
gm1 <- glm(pa ~ bio1 + bio2 + bio15 + bio18 + bio19,
           family = binomial(link = "logit"), data=envtrain)

##look at a summary of the results
summary(gm1)


pg <- predict(bio.data, gm1, ext=e,type="response")
pg<-crop(pg,e)


##plot this
plot(pg, main='GLM probability of occurrence')
##add country boundaries
plot(wrld_simpl, add=TRUE, border='light blue') 

##add our observed locality data
points(Polyergus.coords,col = "black", pch = 4, cex = 0.5)


ge <- evaluate(testpres, testbackg, gm1)
ge


##then we use this evaluation to pick a threshold probability for defining presence/absence
##using the model that gives the most accurate match to observed presence/absence
tr <- threshold(ge, 'prevalence')
plot(pg > tr, main='presence/absence')
plot(wrld_simpl, add=TRUE, border='light blue') 
points(Polyergus.coords,col = "black", pch = 4, cex = 0.5)



##This model uses next 5 climatic variables instead (I used 5 climatic values BO3,4,5,6,7)
gm2 <- glm(pa ~ bio3 + bio4 + bio5 + bio6 + bio7,family = binomial(link = "logit"), data=envtrain)
summary(gm2)

##You can compare two models, even with different predictor variables, using the Akaike Information Criterion
AIC(gm1,gm2)


evaluate(testpres, testbackg, gm1)
evaluate(testpres, testbackg, gm2)




future.bio.data<-cmip6_world(model="CanESM5",var="bio",ssp="245",res=10,time="2061-2080",path=getwd())
names(future.bio.data)<-names(bio.data)

future.bio.data<-crop(future.bio.data,e)


par(mfrow=c(1,2))


#As before i changed the BO things to 1,2,15,18,19
gm1 <- glm(pa ~ bio1 + bio2 + bio15 + bio18 + bio19,
           family = binomial(link = "logit"), data=envtrain)

##and plot the predicted probability of occurrence at present
pg <- predict(bio.data, gm1, ext=e,type="response")
pg<-crop(pg,e)
plot(pg, main='A) GLM present')
plot(wrld_simpl, add=TRUE, border='light blue') 
points(Polyergus.coords,col = "black", pch = 4, cex = 0.5)

## Now predict future distribution, i.e. feeding future bio data to model and plot this
pg.future <- predict(future.bio.data, gm1, ext=e,type="response")
pg.future<-crop(pg.future,e)
plot(pg.future, main="b) GLM, 2060-2081")
plot(wrld_simpl, add=TRUE, border='light blue') 
points(Polyergus.coords,col = "black", pch = 4, cex = 0.5)


predict.localities.now <- extract(pg>=tr, Polyergus.coords)[,-1]
predict.localities.future <- extract(pg.future>=tr, Polyergus.coords)[,-1]
##number of localities where currently absent but present in the future, i.e. range expansion
matches.range.expansion <- predict.localities.now == FALSE & predict.localities.future == TRUE
sum(matches.range.expansion, na.rm = TRUE)
##number of localities where currently present but absent in the future, i.e. range contraction
matches.range.contraction <- predict.localities.now == TRUE & predict.localities.future == FALSE
sum(matches.range.contraction, na.rm = TRUE)



predict.localities.now<-na.omit(predict.localities.now)
sum(predict.localities.now)


predict.localities.future<-na.omit(predict.localities.future)
sum(predict.localities.future)

change.bio.data<-future.bio.data-bio.data


# Q3
# combine occurrence records for both species
# combine occurrence records for both species
occurrences <- rbind(
  data.frame(species = "Polyergus", Polyergus.coords),
  data.frame(species = "Fusca", Fusca.coords)
)

# create a scatter plot of the occurrence records
par(mfrow=c(1,1))
ggplot(occurrences, aes(x = lon, y = lat, color = species)) +
  geom_point() +
  geom_point(data = inner_join(data.frame(Polyergus.coords), data.frame(Fusca.coords),
                               by = c("lon", "lat"), multiple = "all"),  
             aes(color = "Overlap"), size = 3, shape = 21) +
  labs(title = "Distribution overlap of Polyergus and Fusca") +
  scale_color_manual(values = c("Polyergus" = "darkred", "Fusca" = "darkblue", "Overlap" = "green"))

# Calculate the number of overlapping points
overlap_count <- nrow(inner_join(data.frame(Polyergus.coords), data.frame(Fusca.coords), by = c("lon", "lat"), multiple = "all"))


# Calculate the total number of unique points for both species
total_unique_points <- nrow(distinct(occurrences, lon, lat))

# Calculate the Jaccard index
jaccard_index <- overlap_count / total_unique_points

# Print the Jaccard index
cat("Jaccard index:", jaccard_index)


###Q4###

##### Q5####
## PREDICT FUTURE SPECIES RANGES
# combine occurrence records for both species
occurrences_poly <- rbind(
  data.frame(species = "Polyergus", Polyergus.coords),
  data.frame(species = "Fusca", Fusca.coords)
)

# create a raster layer covering the extent of both distributions
r_poly <- raster(extent(min(Polyergus.coords$lon, Fusca.coords$lon),
                        max(Polyergus.coords$lon, Fusca.coords$lon),
                        min(Polyergus.coords$lat, Fusca.coords$lat),
                        max(Polyergus.coords$lat, Fusca.coords$lat)),
                 nrow = 100, ncol = 100)

library(raster)
# extract bioclimatic variables from the raster layer
bio.data_poly <- extract(biovars, r_poly)

# create a data frame with the bioclimatic variables and presence/absence data
pa_poly <- extract(pa,r_poly) 
pa_poly <- ifelse(is.na(pa_poly),0,pa_poly) 
envtrain_poly <- data.frame(pa=pa_poly,bio.data_poly)

# fit the model with the present data in envtrain_poly
gm_poly <- glm(pa ~ bio1 + bio2 + bio15 + bio18 + bio19,
               family = binomial(link = "logit"), data=envtrain_poly)
## Download one set of future climate data for period 2061-2080
## There are multiple future climate models, and versions of each, which are stored in CMIP6
## This code just extracts one model and set of parameters for you to try out the methods
future.bio.data <- cmip6_world(model = "CanESM5", var = "bio", ssp = "245", res = 10, time = "2061-2080", path = getwd())
names(future.bio.data) <- names(bio.data)

## Crop to just keep the region of interest to speed up some steps
future.bio.data <- crop(future.bio.data, e)

## Plot present and future distributions next to each other
par(mfrow = c(1, 2))

## Fit the model with the present data in envtrain_poly for Polyergus
gm_poly <- glm(pa ~ bio1 + bio2 + bio15 + bio18 + bio19, family = binomial(link = "logit"), data = envtrain_poly)

## Predict the probability of occurrence at present for Polyergus
pg_poly <- predict(bio.data, gm_poly, ext = e, type = "response")
pg_poly <- crop(pg_poly, e)
plot(pg_poly, main = "A) GLM present - Polyergus")
plot(wrld_simpl, add = TRUE, border = "light blue") 
points(Polyergus.coords, col = "black", pch = 4, cex = 0.5)

## Predict future distribution for Polyergus
pg_future_poly <- predict(future.bio.data, gm_poly, ext = e, type = "response")
pg_future_poly <- crop(pg_future_poly, e)
plot(pg_future_poly, main = "B) GLM, 2061-2081 - Polyergus")
plot(wrld_simpl, add = TRUE, border = "light blue") 
points(Polyergus.coords, col = "black", pch = 4, cex = 0.5)

## Fit the model with the present data in envtrain_fusca for Fusca
gm_fusca <- glm(pa ~ bio1 + bio2 + bio3 + bio4 + bio5 + bio6 + bio7 + bio8 + bio9 + bio10 + bio11 + bio12 + bio13 + bio14 + bio15 + bio16 + bio17 + bio18 + bio19, family = binomial(link = "logit"), data = envtrain_fusca)

## Predict the probability of occurrence at present for Fusca
pg_fusca <- predict(bio.data, gm_fusca, ext = e, type = "response")
pg_fusca <- crop(pg_fusca, e)
plot(pg_fusca, main = "A) GLM present - Fusca")
plot(wrld_simpl, add = TRUE, border = "light blue") 
points(Fusca.coords, col = "black", pch = 4, cex = 0.5)

## Predict future distribution for Fusca
pg_future_fusca <- predict(future.bio.data, gm_fusca, ext = e, type = "response")
pg_future_fusca <- crop(pg_future_fusca, e)
plot(pg_future_fusca, main = "B) GLM, 2061
