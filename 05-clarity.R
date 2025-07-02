#!/usr/local/bin/R


##remotes::install_github("danjlawson/CLARITY/Clarity")
library("readxl")
library("writexl")

library("Clarity")
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggrepel)

read_excel_as=function(...,type=as.matrix){
    df=as.data.frame(read_excel(...))
    rownames(df)<-df[[1]]
    df[[1]]<-NULL
    df[] <- lapply(df, function(x) as.numeric(as.character(x)))
    type(df)
}
get_order<-function(x){
    hcr <- hclust(dist(x))
    ddr <- as.dendrogram(hcr)
    Rowv <- rowMeans(x)
    ddr <- reorder(ddr, Rowv)
    labels(ddr)
}
really_as_matrix<-function(A){
    ret=matrix(as.numeric(A),dim(A)[1],dim(A)[2]) # Coerce to standard matrix
    rownames(ret)=rownames(A)
    colnames(ret)=colnames(A)
    ret
}

## Read data
gene=read_excel_as("Pipeline results/matched_genetic_similarity.xlsx",type=as.matrix)
lang=read_excel_as("Pipeline results/matched_linguistic_similarity.xlsx",type=as.matrix)

## Read metadata

genref=as.data.frame(read_excel("data-manual/genetic_mapping_v1.xlsx"))
genref=genref[genref[,"LinguisticName"]%in%rownames(lang),]
rownames(genref)<-genref[,"LinguisticName"]

## Make groupings
grouping<-genref[rownames(lang),"manual_cluster"]
names(grouping) <-rownames(lang)
#sort(grouping)

## Add to metadata
metadata=as.data.frame(read_excel("Metadata/merged_metadata.xlsx"))
rownames(metadata)=metadata[,"Reference"]
metadata=metadata[rownames(lang),]
metadata[,"grouping"]<-factor(paste("Cluster",grouping))
write_xlsx(metadata,"Metadata/merged_metadata_withclusters.xlsx")
## 
                              
## Plot the data on a map for reference

world <- ne_countries(scale = "medium", returnclass = "sf")
points_df=data.frame(lon=metadata[,"Coordinate 2"],
                     lat=metadata[,"Coordinate 1"],
                     label=rownames(metadata),
                     group=metadata[,"grouping"]),
                     negrito=metadata[,"Cluster 1"]=="Negrito")
rownames(points_df)=points_df[,"label"]
grouping_map=c("darkgrey","red","blue","orange",
               "darkgreen","purple")
#cluster_cols=get_colors(length(unique(metadata[,"grouping"])))

pdf("Pipeline results/mergedData_Map.pdf",height=8,width=8)
ggplot(data = world) +
  geom_sf(fill = "antiquewhite") +
    geom_point(data = points_df,
               aes(x = lon, y = lat, color = group,shape = negrito), size = 2) +
    coord_sf(xlim = range(metadata[,"Coordinate 2"])+c(-5,5),
             ylim = range(metadata[,"Coordinate 1"])+c(-5,5), expand = FALSE) +
    scale_color_manual(values = grouping_map) +
    scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16)) +
    geom_text_repel(data = points_df,
              aes(x = lon, y = lat, label = label, color = group), 
              size = 3,max.overlaps = 20) +
  theme_minimal() +
  labs(title = "",
       x = "Longitude", y = "Latitude")
dev.off()

## Use this to make a partial ordering
tdata=scale(gene) + scale(lang)
tdata=tdata+t(tdata)
ugroups=levels(points_df[,"group"])
myorder=do.call("c",lapply(ugroups,function(x)get_order(tdata[points_df[,"group"]==x,])))

metadata=metadata[myorder,]
gene=gene[myorder,myorder]
lang=lang[myorder,myorder]
points_df=points_df[myorder,]

names(grouping_map)=ugroups
grouping_colors<-grouping_map[points_df[,"group"]]


kmax=25
lang=scale(lang)
lang = (lang+t(lang))/2
gene=scale(gene)
gene = (gene + t(gene))/2
langscan<-Clarity_Scan(lang,kmax=kmax)
genescan<-Clarity_Scan(gene,kmax=kmax)
langpredictsgene<-Clarity_Predict(gene,langscan)
genepredictslang<-Clarity_Predict(lang,genescan)

langpredictsgene_resid=Clarity_Extract(langpredictsgene,summary=I,k=10)
genepredictslang_resid=Clarity_Extract(genepredictslang,summary=I,k=10)
langpredictsgene_resid=really_as_matrix(langpredictsgene_resid)
genepredictslang_resid=really_as_matrix(genepredictslang_resid)

tv=sort(as.numeric(langpredictsgene_resid))
gene_thresh=quantile(abs(tv),0.9)
tv=sort(as.numeric(genepredictslang_resid))
lang_thresh=quantile(abs(tv),0.9)

par(mfrow=c(1,2))
plot(tv,col=c("black","red")[1+(abs(tv)>gene_thresh)])
plot(tv,col=c("black","red")[1+(abs(tv)>lang_thresh)])


my_palette <- colorRampPalette(c("blue", "white", "red"))(100)

mypar=c(5,5,2,1)
pdf("Pipeline results/OCSEAN_Clarity_Genes_Language.pdf",height=6,width=12)
par(mfrow=c(2,3))
Clarity_Chart(gene,scalefun=I,cex.axis=0.4,las=2,mar=mypar,
              col.axis.Y=grouping_colors,col.axis.X=grouping_colors)
mtext("Genes",side=3)
Clarity_Chart(langpredictsgene_resid,scalefun=I,cex.axis=0.4,las=2,mar=mypar,
              zlim=c(-1,1)*max(abs(range(langpredictsgene_resid))),
              col.axis.Y=grouping_colors,col.axis.X=grouping_colors,
              cols=my_palette,
              signif=abs(langpredictsgene_resid)>gene_thresh)
mtext("Language predicting Genes",side=3)
plot(langpredictsgene,rotate=TRUE,cex.axis=0.4,las=2,mar=mypar,
              col.axis.Y=grouping_colors)
abline(v=10)
## Second row
Clarity_Chart(lang,scalefun=I,cex.axis=0.4,las=2,mar=mypar,
              col.axis.Y=grouping_colors,col.axis.X=grouping_colors)
mtext("Language",side=3)
Clarity_Chart(genepredictslang_resid,scalefun=I,cex.axis=0.4,las=2,mar=mypar,
              zlim=c(-1,1)*max(abs(range(genepredictslang_resid))),
              col.axis.Y=grouping_colors,col.axis.X=grouping_colors,
              cols=my_palette,
              signif=abs(genepredictslang_resid)>lang_thresh)
mtext("Genes predicting Language",side=3)
plot(genepredictslang,rotate=TRUE,cex.axis=0.4,las=2,mar=mypar,
              col.axis.Y=grouping_colors)
abline(v=10)
dev.off()

#################################


library(igraph)
library(umap)

clarity_network<-function(A,main="",legtext=c("smaller than predicted","larger than predicted")){
    diag(A)=0 # remove self-loops
    Alab=substr(rownames(A),1,3)
    
    set.seed(1)

    g <- graph_from_adjacency_matrix(A, mode = "undirected",
                                     weighted = TRUE, diag = FALSE)
    layout_coords <- layout_with_fr(g, weights = abs(E(g)$weight))
                                        #layout_coords=umap(langpredictsgene_resid)$layout
    par(mar=c(2,2,2,1))
    plot(g,
         layout = layout_coords,
         vertex.size = 11,
         vertex.label = Alab,
         vertex.label.cex = 0.7,
         vertex.color=grouping_colors,
         edge.color = ifelse(E(g)$weight > 0, "blue", "red"),
         edge.lty = ifelse(E(g)$weight > 0, 2, 1),
         xlim=c(-1,1.8))
    mtext(main,3)
    legend("right",legend=paste(Alab,rownames(A),sep=": "),cex=0.5,text.col=grouping_colors)
    legend("bottomright",
           legend=legtext,
           lty=1,col=c("blue","red"),text.col=c("blue","red"))
    
}

Agene <- (abs(langpredictsgene_resid)>gene_thresh)*sign(langpredictsgene_resid)
Alang <- (abs(genepredictslang_resid)>lang_thresh)*sign(genepredictslang_resid)
pdf("results/Residuals_Network.pdf",height=6,width=12)
par(mfrow=c(1,2))
clarity_network(Agene,"Genetic residuals network",
                legtext=c("Genetically less similar than linguistics",
                          "Genetically more similar than linguistics"))
clarity_network(Alang,"Language residuals network",
                legtext=c("Linguistically less similar than genetics",
                          "Linguistically more similar than genetics"))
dev.off()
