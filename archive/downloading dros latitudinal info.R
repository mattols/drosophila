library("spocc")

gs<-read.csv("drosophila2_17_22.csv", as.is=T)
str(gs)
head(gs)

i<-1
latinfo<-as.data.frame(matrix(,153,11))
colnames(latinfo)<-c("Species", "GbifMinlat", "GbifMaxlat", "BisonMaxLat", "BisonMinLat", "inatMinLat", "inatMaxLat", "alaMinLat", "alaMaxLat", "idigbioMinLat", "idigbioMaxlat")

for(i in 1:length(gs$Name)){
  test<-occ(query=(gs[i,1]), from='gbif', limit=200)
  foo<-test$gbif$data[[1]]$latitude
  foo2<-na.omit(foo)
  minlat<-min(foo2)
  maxlat<-max(foo2)
  latinfo[i,1]<-gs[i,1]
  latinfo[i,2]<-minlat
  latinfo[i,3]<-maxlat
  
  test<-occ(query=(gs[i,1]), from='bison', limit=200)
  foo<-test$bison$data[[1]]$latitude
  foo2<-na.omit(foo)
  minlat<-min(foo2)
  maxlat<-max(foo2)
  latinfo[i,4]<-minlat
  latinfo[i,5]<-maxlat
  
  
  test<-occ(query=(gs[i,1]), from='inat', limit=200)
  foo<-test$inat$data[[1]]$latitude
  foo2<-na.omit(foo)
  minlat<-min(foo2)
  maxlat<-max(foo2)
  latinfo[i,6]<-minlat
  latinfo[i,7]<-maxlat
  
  test<-occ(query=(gs[i,1]), from='ala', limit=200)
  foo<-test$ala$data[[1]]$latitude
  foo2<-na.omit(foo)
  minlat<-min(foo2)
  maxlat<-max(foo2)
  latinfo[i,8]<-minlat
  latinfo[i,9]<-maxlat
  
  test<-occ(query=(gs[i,1]), from='idigbio', limit=200)
  foo<-test$idigbio$data[[1]]$latitude
  foo2<-na.omit(foo)
  minlat<-min(foo2)
  maxlat<-max(foo2)
  latinfo[i,10]<-minlat
  latinfo[i,11]<-maxlat
  
  
}

head(latinfo)

write.csv(latinfo, "droslatitudeinfogs.csv")


test<-occ(query=("Drosophila melanogaster"), from='idigbio', limit=200)
test$idigbio$data$Drosophila_melanogaster$latitude

latinfo
test<-occ(query=("Drosophila"), from='idigbio', limit=200)
test$idigbio$data$Drosophila$latitude

##absolute value when fitting model
##cluniomax<-max(abs(foo2))
##cluniomin<-min(abs(foo2))


?paste



