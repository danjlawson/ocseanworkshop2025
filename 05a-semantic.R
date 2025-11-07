library(RColorBrewer)
library(colorspace)
#install.packages(c("sf", "dplyr", "tmap", "tmaptools"))
library(sf)
require(concaveman)
library(dplyr)
library(tmap)
library(tmaptools)

# Darken all colors by a factor
darken_colors <- function(cols, amount = 0.3) {
  darken(cols, amount = amount)
}
get_colors<-function(n,darken_by=0.2,seed=2){
    set.seed(seed)
    exclude_palettes <- c("Pastel1", "Pastel2", "Set3","Oranges")
    qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual' & 
                                     !(rownames(brewer.pal.info) %in% exclude_palettes), ]
    col_vector = unlist(mapply(brewer.pal,
                               qual_col_pals$maxcolors,
                               rownames(qual_col_pals)))
    myorder=sample(1:length(col_vector),n)
    darken_colors(col_vector[myorder], amount = darken_by)
}
get_dists<-function(X,refnames,refvec){
    factors=unique(refvec)
    unames=unique(refnames)
    faclist=lapply(factors,function(fac){
        tw=which(refvec==fac)
        X[tw,]
    })
    names(faclist)=factors
    distlists=lapply(factors,function(f1){
        r=lapply(factors,function(f2){
            r=sqrt(rowSums((faclist[[f1]] - faclist[[f2]])^2))
            names(r)=unames
            r
        })
        names(r)<-factors
        r
    })
    names(distlists)=factors
    distlists
}
get_mean_matrix<-function(x){
    sapply(x,function(xx){
        sapply(xx,mean)
    })
}
get_edge_list<-function(d_lists,thresh){
    dnames=names(d_lists)
    ret=data.frame("d1"=numeric(),"d2"=numeric(),"dist"=numeric(),lang=character())
    for(i1 in 1:(length(dnames)-1)){
        d1=dnames[i1]
        for(i2 in (i1+1):length(dnames)){
            d2=dnames[i2]
            tw=d_lists[[d1]][[d2]]>thresh
            for(i in which(tw)){
                ret<-rbind(ret,
                           data.frame("d1"=d1,"d2"=d2,
                                      dist=d_lists[[d1]][[d2]][i],
                                      lang=names(d_lists[[d1]][[d2]])[i]))
            }
        }        
    }
    ret
}
get_edges<-function(edge_list,X,refdf){
    ret=data.frame()
    for(i in 1:dim(edge_list)[1]){
        tw1=which(refdf[,"field"]==edge_table[i,"d1"] & refdf[,"lang"]==edge_table[i,"lang"])
        tw2=which(refdf[,"field"]==edge_table[i,"d2"] & refdf[,"lang"]==edge_table[i,"lang"])
        ret=rbind(ret,data.frame(i=tw1,
                                 j=tw2,
                                 lang=edge_table[i,"lang"],
                                 s1=edge_table[i,"d1"],
                                 s2=edge_table[i,"d2"]))
    }
    ret
}
draw_edges<-function(edges,Xi,Xj,col=rep("grey",dim(edges)[1]),...){
    startXi = Xi[edges[,"i"]]
    endXi = Xj[edges[,"i"]]
    startXj = Xi[edges[,"j"]]
    endXj = Xj[edges[,"j"]]
    for(i in 1:dim(edges)[1]){
        lines(c(startXi[i],startXj[i]),
              c(endXi[i],endXj[i]),col=col[i],...)
    }
}
get_means<-function(X,refdf,order,what="lang"){
    ulang=order#levels(refdf[,"lang"])
    ret=sapply(ulang,function(x){
        colMeans(X[refdf[,what]==x,])
    })
    t(ret)
}
polygon_containing<-function(points){
    colnames(points)[1:2]<- c("lon", "lat")
    dat_sf <- st_as_sf(as.data.frame(points), coords = c("lon", "lat"), crs = 4326)
    polygons <- concaveman(dat_sf)
    st_coordinates(st_as_sf(polygons))[,1:2]
}
                                                                      
## Currently requires functions from 05-clarity.R!!

## Read the distances of all semantic fields (and the list of all semantic fields to get that)
semantic_df=as.data.frame(read_excel("semantic_field_list.xlsx"))
semantic_fields=semantic_df[,"Semantic_Field"]

semantic_list_full=lapply(semantic_fields,function(x){
    read_excel_as(paste0("Pipeline\ results/Semantic/matched_linguistic_semantic_field_",x,".xlsx"))
})
names(semantic_list_full)=semantic_fields

##semantic_list_complete=c(lapply(semantic_list_full,function(x)log(1+x)),
semantic_list_complete=c(semantic_list_full,
                         list("gene"=gene,"lang"=lang))
pdf("Pipeline results/OCSEAN_SemanticDistancesSeparate.pdf",height=6,width=6)
for(i in 1:length(semantic_list_complete)){
    Clarity_Chart(semantic_list_complete[[i]][rownames(gene),colnames(gene)],
                  scalefun=I,cex.axis=0.4,las=2,mar=mypar,
                  col.axis.Y=grouping_colors,col.axis.X=grouping_colors,
                  main=names(semantic_list_complete)[i])
}
dev.off()

## Keep only chosen semantic fields and languages 
keep_semantic_fields=names(semantic_list_full)[!names(semantic_list_full)%in%c("Social and political relations","Law","Religion and belief","Miscellaneous function words","Animals","Warfare and hunting")]
keep_languages=rownames(semantic_list_full[[1]])[!rownames(semantic_list_full[[1]])%in%c("Sangil","Ati")]
# Warfare and hunting

## Read the raw data for count information
library("data.table")
wordlist<-data.frame(fread("Pipeline\ results/OCSEAN_processed_joineddata_semantic.tsv",sep="\t"))
semantic_full_table=table(wordlist[,c("Language","Semantic_Field")])
ttab=data.frame("Concepts"=rowSums(semantic_full_table),"Matched"=rownames(semantic_full_table)%in%rownames(gene))
write.csv(ttab,"reporting_table.csv")
semantic_table=semantic_full_table[keep_languages,keep_semantic_fields]
concepts_in_fields=lapply(semantic_df[,2],function(x)
    unique(wordlist[wordlist[,"Semantic_Field"]==x,"Concept"]))
names(concepts_in_fields)=semantic_df[,2]

## Sequentially remove row/columns that to ensure that all semantic fields/languages have a minimum number of concepts
## We used this to suggest rows/columns to remove, which have mow been specidied above.
print(dim(semantic_table))
min_val=10
while(any(semantic_table<min_val)){
    worst_cols=sort(colSums(semantic_table<min_val),decreasing=T)
    worst_rows=sort(rowSums(semantic_table<min_val),decreasing=T)
    if(max(worst_cols)>=max(worst_rows)){
        worst_cols=names(worst_cols)[worst_cols==max(worst_cols)]
        print(paste("Dropping cols:",worst_cols))
        semantic_table=semantic_table[,!colnames(semantic_table)%in%worst_cols]        
    }else{
        worst_rows=worst_rows[worst_rows==max(worst_rows)]
        print(paste("Dropping rows",rownames(semantic_table)[worst_rows]))
        semantic_table=semantic_table[!rownames(semantic_table)%in%worst_rows,]
    }
    print(dim(semantic_table))
}

## Summarise the semantic structure
thm=heatmap(semantic_table,keep.dendro=T)
##    labels(thm$Rowv)
semantic_lang_order=rownames(lang)[rownames(lang)%in%rownames(semantic_table)]
semantic_colors=grouping_colors[semantic_lang_order]
semantic_order=labels(thm$Colv)

pdf("Pipeline results/SemanticConceptCountTable.pdf",height=6,width=15)
Clarity_Chart((semantic_table[semantic_lang_order,semantic_order]),scalefun=I,text=T,las=2,
              main="Semantic concept counts",mar=c(10,12,3,1),
              col.axis.X=semantic_colors)
dev.off()


keep_semantic_fields=colnames(semantic_table)
keep_languages=rownames(semantic_table)
keep_languages%in%rownames(gene)
keep_languages=rownames(gene)[rownames(gene) %in% keep_languages]

semantic_list=semantic_list_full[keep_semantic_fields]
semantic_list=lapply(semantic_list,function(x)x[keep_languages,keep_languages])
semantic_list=c(list("lang"=lang[keep_languages,keep_languages]),semantic_list)
semantic_list=c(list("gene"=gene[keep_languages,keep_languages]),semantic_list)
semantic_list=lapply(semantic_list,scale)

## Make the unfolded matrix
ndata=length(semantic_list)
nperdata=dim(semantic_list[[1]])[1]
semantic_data=do.call("cbind",semantic_list)
colnames(semantic_data)=paste(rep(names(semantic_list),each=nperdata),
                              rep(rownames(semantic_list[[1]]),
                                  by=ndata))
refdf=data.frame(name=colnames(semantic_data),
                 field=as.factor(rep(names(semantic_list),each=nperdata)),
                 lang=as.factor(rep(rownames(semantic_list[[1]]),
                                    by=ndata)))
refdf$group=metadata[refdf[,"lang"],"grouping"]
ref_order=as.character(refdf[1:nperdata,"lang"])
plot_labels=paste(rep(names(semantic_list),each=nperdata),
                  rep(rownames(semantic_list[[1]]),by=ndata),
                  sep="\n")

## Make a version that is a vector for each semantic field
flat_data=as.data.frame(sapply(semantic_list,function(x)as.numeric(x[upper.tri(x)])))
tcor=cor(flat_data)
tthresh=tcor["gene","lang"]
mean(tcor[upper.tri(tcor)]<tthresh)
#plot(flat_data)
#heatmap(cor(flat_data),symm=TRUE,scale="none")

thm=heatmap(cor(semantic_table),symm=TRUE,scale="none",keep.dendro=TRUE)
semantic_order=labels(thm$Rowv)

## UASE
##scaled_semantic_data<-scale(semantic_data)
mysvd<-svd(semantic_data)
d<-5 # Or 12. Assessed from plot(mysvd$d)
X<-mysvd$v[,1:d] %*% sqrt(diag(mysvd$d[1:d]))
rownames(X)<-colnames(semantic_data)
meanX=get_means(X,refdf,ref_order)
meanXsemantic=get_means(X,refdf,semantic_order,what="field")



## Distances between semantic categories
semantic_d_lists=get_dists(X,refdf[,"lang"],refdf[,"field"])
thresh=1
plot(sort(unlist(semantic_d_lists))) ## threshold chosen from empirical elbow
abline(h=thresh)
edge_table<-get_edge_list(semantic_d_lists,thresh)
edges<-get_edges(edge_table,X,refdf)

edge_col_idx=sapply(edges[,"lang"],function(x)which(levels(refdf[,"lang"])==x))
highlight=unique(names(edge_col_idx))

leglocs=c("topleft","bottom","topleft","topleft","topleft","topleft","topleft","topleft","topleft","topleft")
##Most useful visualisation
pdf("Pipeline results/GeneLangSemanticChanges.pdf",height=6,width=8)
ii=1
for(i in 1:4) for (j in (i+1):5)
{
    par(mar=c(4,4,2,1))
plot(X[,i],X[,j],type="p",
     pch=as.numeric(refdf[,"field"]),cex=c(0.5,0.7)[1+rownames(meanX)%in%highlight],
     col=get_colors(nperdata)[refdf[,"lang"]],xlab=paste("PC",i),ylab=paste("PC",j))
legend(leglocs[ii],pch=1:length(levels(refdf[,"field"])),
       legend=as.character(levels(refdf[,"field"])),cex=0.7)
draw_edges(edges,X[,i],X[,j],col=get_colors(nperdata)[edge_col_idx])
text(meanX[,i],meanX[,j],labels=rownames(meanX),
     col=get_colors(nperdata)[refdf[1:nperdata,"lang"]],
     cex=c(0.4,0.6)[1+rownames(meanX)%in%highlight])
    ii=ii+1
}
dev.off()

## X[refdf[,"lang"]=="bali",]


###########################
## Geographic linkage
geo=points_df[rownames(points_df)%in%semantic_lang_order,1:2]
geo2=data.frame("lon"=rep(geo[,1],times=ndata),
                "lat"=rep(geo[,2],times=ndata))
rownames(geo2)=rownames(X)
points_df2=points_df[rownames(geo),]
points_df2["Balinese",c("lon","lat")]<-c(120,6)


## UMAP
set.seed(1)
umap_params=umap.defaults
umap_params$n_neighbors=40
Xumap=umap(X,umap_params)$layout
meanXumap=get_means(Xumap,refdf,ref_order)
meanXumapsemantic=get_means(Xumap,refdf,semantic_order,what="field")

df_umap=data.frame("x"=Xumap[,1],
                   "y"=Xumap[,2],
                   lang=as.factor(refdf[,"lang"]),
                   field=as.factor(refdf[,"field"]),
                   group=as.factor(refdf[,"group"]))
df_umap <- df_umap %>%
  mutate(class = case_when(
    field == "gene" ~ "gene",
    field == "lang" ~ "lang",
    TRUE            ~ "semantic"
    ))
df_meanumap=data.frame("x"=meanXumap[,1],
                       "y"=meanXumap[,2],
                       "lang"=rownames(meanXumap),
                       group=points_df2$group)
df_meanumapsemantic=data.frame("x"=meanXumapsemantic[,1],
                       "y"=meanXumapsemantic[,2],
                       "field"=rownames(meanXumapsemantic))
bali_poly_umap <- polygon_containing(df_umap[df_umap[,"lang"]=="Balinese",1:2])

pdf("Pipeline results/GeneLangSemanticChangesUmap.pdf",height=6,width=6)
##ggplot(df_umap, aes(x = x, y = y)) +
##    geom_point(aes(color = lang,shape=class, alpha = class), size = 2) +
##    scale_alpha_manual(values = c("gene" = 1,"lang"=1, "semantic" = 0.2)) +
##    geom_text_repel(data = df_meanumap,max.overlaps=20,
##                    aes(x = x, y = y, label = lang,color=lang)) +
##    guides(color = "none") +
##    theme_minimal()
ggplot(df_umap, aes(x = x, y = y)) +
    geom_point(aes(shape=class, alpha = class), size = 2) +
    scale_color_manual(values = grouping_map) +
    scale_alpha_manual(values = c("gene" = 1,"lang"=1, "semantic" = 0.2)) +
    geom_text_repel(data = df_meanumap,max.overlaps=20,
                    aes(x = x, y = y, label = lang,color=group),show.legend = FALSE) +
    guides(class = "none") +
    theme_minimal()
dev.off()

plot_base<-function(Xall,meanX,highlight=character(),highlightcex=c(0.5,1),
                    xlab="",ylab="",legpos="topleft",...){
    plot(Xall[,1],Xall[,2],type="p",
         pch=as.numeric(refdf[,"field"]),
         cex=highlightcex[1+rownames(meanX)%in%highlight],
         col=get_colors(nperdata)[refdf[,"lang"]],
         xlab="",ylab="",...)
    legend(legpos,pch=1:length(levels(refdf[,"field"])),
           legend=as.character(levels(refdf[,"field"])),cex=0.7)
    text(meanX[,1],meanX[,2],labels=rownames(meanX),
         col=get_colors(nperdata)[refdf[1:nperdata,"lang"]],
         cex=highlightcex[1+rownames(meanX)%in%highlight])
}
make_edges<-function(refdf,
                     language=as.character(unique(refdf[,"lang"])),
                     field=as.character(unique(refdf[,"field"]))){
    tw=which((refdf[,"field"]%in%field)&(refdf[,"lang"]%in%language))
    tw2=expand.grid(tw,tw)
    tw2=tw2[tw2[,2]>tw2[,1],]
    data.frame(i=tw2[,1],j=tw2[,2])
}

##pdf("Pipeline results/GeneLangSemanticChangesUmap.pdf",height=6,width=8)
##plot(Xumap[,1],Xumap[,2],type="p",
##     pch=as.numeric(refdf[,"field"]),cex=c(0.5,1)[1+(refdf[,"field"]=="gene")],
##     col=get_colors(nperdata)[refdf[,"lang"]],xlab=paste("UMAP",1),ylab=paste("UMAP",2))
##legend("topleft",pch=1:length(levels(refdf[,"field"])),
##       legend=as.character(levels(refdf[,"field"])),cex=0.7)
#draw_edges(edges,Xumap[,1],Xumap[,2],col=get_colors(nperdata)[edge_col_idx])
##text(meanXumap[,1],meanXumap[,2],labels=rownames(meanXumap),
##     col=get_colors(nperdata)[refdf[1:nperdata,"lang"]],
##     cex=c(0.4,0.6)[1+rownames(meanXumap)%in%highlight])
##dev.off()

## Plotting the semantic correlations
semantic_d=get_mean_matrix(semantic_d_lists)
myorder=c("lang",semantic_order,"gene")
myorder2=semantic_order#c("The body", "Agriculture and vegetation", "Motion", "Food and drink", "Basic actions and technology")
#thm=heatmap(semantic_d,symm=TRUE,scale="none",keep.dendro=TRUE)
#myorder=labels(thm$Rowv)

pdf("Pipeline results/Semantic_Structure.pdf",height=10,width=10)
Clarity_Chart(semantic_d[myorder,myorder],scalefun=I,text=T,las=2,
              main="Average distance moved between...",mar=c(12,12,3,1),
              signif=semantic_d[myorder,myorder]>semantic_d["gene","lang"],
              signiffade="ff")
dev.off()

genelangdists=data.frame(semantic_d[-(1:2),c("gene","lang")])
genelangdists[,"Semantic field"]=rownames(genelangdists)

pdf("Pipeline results/Semantic_DistanceFromGeneLang.pdf",height=5,width=5)
ggplot(genelangdists, aes(x = gene, y = lang)) +
  geom_point(color = "blue", size = 2) +
    geom_text_repel(aes(label = `Semantic field`), size = 3,max.overlaps = 20) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
    labs(x = "Distance from Genetics", y = "Distance from typical Language") +
    annotate("text",
           x = 0.2,0.6,             # left & top
           label = "Closer to Genes",
           hjust = 0, vjust = 1,             # flush left, flush top
           size = 4, fontface = "italic") +
    annotate("text",
           x = 0.7, y = 0.2,             # right & bottom
           label = "Closer to Language",
           hjust = 1, vjust = 0,             # flush right, flush bottom
           size = 4, fontface = "italic") +
    coord_sf(xlim = c(0.1,0.85),ylim =c(0.1,0.85), expand = FALSE) +
    theme_minimal()
dev.off()


pdf("Pipeline results/Semantic_Table.pdf",height=10,width=10)
Clarity_Chart(t(semantic_table[semantic_lang_order,myorder2]),scalefun=I,text=T,las=2,
              main="Semantic concept counts",mar=c(10,12,3,1),
              col.axis.Y=semantic_colors)
dev.off()


## Procrust the estimate onto geography
library("vegan")
Xpro_vegan=procrustes(geo2,X, scale=TRUE)## S3 method for class 'procrustes':
Xpro=predict(Xpro_vegan)
my_scale=function(x1,x2){
    mx1=colMeans(x1)
    mx2=colMeans(x2)
    for(i in 1:dim(x1)[2]){
        x1[,i]=(x1[,i]-mx1[i])/sd(x1[,i])*sd(x2[,i])+mx2[i]
    }
    x1
}
Xpro=my_scale(Xpro[,1:2],geo)
meanXpro=data.frame(get_means(Xpro,refdf,ref_order))

## Construct a data frame containing this information
points_df2[,"lonhat"]<-meanXpro[,1]
points_df2[,"lathat"]<-meanXpro[,2]
library(dplyr)
df_actual <- points_df2 %>%
  mutate(x = lon, y = lat, type = "Geography")
df_estimate <- points_df2 %>%
  mutate(x = lonhat, y = lathat, type = "Estimate")
df_combined <- bind_rows(df_actual, df_estimate)

gene_poly_geo <- polygon_containing(Xpro[refdf[,"field"]=="gene",1:2])
bali_poly_geo <- polygon_containing(Xpro[refdf[,"lang"]=="Balinese",1:2])
palawano_poly_geo <- polygon_containing(Xpro[refdf[,"lang"]=="Palawano",1:2])
meranaw_poly_geo <- polygon_containing(Xpro[refdf[,"lang"]=="Meranaw",1:2])
gaddang_poly_geo <- polygon_containing(Xpro[refdf[,"lang"]=="Gaddang",1:2])
ata_poly_geo <- polygon_containing(Xpro[refdf[,"lang"]=="Ata",1:2])
# Extract coordinates

## Plot it
pdf("Pipeline results/GeneLangGeo.pdf",height=6,width=8)
ggplot(data = world) +
  geom_sf(fill = "antiquewhite") +
  # Segment from estimate to actual
    geom_segment(data = points_df2,aes(x = lonhat, y = lathat,
                                       xend = lon, yend = lat, color = group), alpha = 0.5) +
    geom_point(data=df_combined,aes(x = x, y = y, color = group, shape = type), size = 1) +
    geom_polygon(data = gaddang_poly_geo,
               aes(x = X, y = Y),fill = "#0000ff33", color = NA) +
#    geom_polygon(data = bali_poly_geo,
#               aes(x = X, y = Y),fill = "#33333333", color = NA) +
    geom_polygon(data = ata_poly_geo,
               aes(x = X, y = Y),fill = "#ffff0033", color = NA) +
    geom_polygon(data = palawano_poly_geo,
               aes(x = X, y = Y),fill = "#ff000033", color = NA) +
    scale_shape_manual(values = c("Geography" = 1, "Estimate" = 16)) +
    labs(shape = "Point Type", color = "Group") +
    coord_sf(xlim = range(metadata[,"Coordinate 2"])+c(-5,5),
                          ylim =c(2,20), expand = FALSE) +
    scale_color_manual(values = grouping_map) +
    geom_text_repel(data = points_df2,
              aes(x = lon, y = lat, label = label, color = group), 
              size = 3,max.overlaps = 40,show.legend = FALSE) +
    annotate("text",
             x = 117,13,             # left & top
             label = "Gaddang\nFeatures",
             hjust = 0, vjust = 1, col="blue",            # flush left, flush top
             size = 4, fontface = "italic") +
#    annotate("text",
#             x = 121,4,             # left & top
#             label = "Balinese\nFeatures",
#             hjust = 0, vjust = 1, col="darkgrey",            # flush left, flush top
#             size = 4, fontface = "italic") +
    annotate("text",
             x = 125,5,             # left & top
             label = "Ata\nFeatures",
             hjust = 0, vjust = 1, col="orange",            # flush left, flush top
             size = 4, fontface = "italic") +
    annotate("text",
             x = 114,10,             # left & top
             label = "Palawano\nFeatures",
             hjust = 0, vjust = 1, col="red",            # flush left, flush top
             size = 4, fontface = "italic") +
    theme_minimal() +
  labs(title = "",
       x = "Longitude", y = "Latitude")
dev.off()

## Plot the semantic components
pdf("Pipeline results/GeneLangSemanticChangesGeo.pdf",height=6,width=8)
plot(Xpro[,1],Xpro[,2],type="p",
     pch=as.numeric(refdf[,"field"]),
     cex=c(0.5,0.7)[1+rownames(meanXpro)%in%highlight],
     col=get_colors(nperdata)[refdf[,"lang"]],
     xlab="Long",ylab="Lat")
legend("topright",pch=1:length(levels(refdf[,"field"])),
       legend=as.character(levels(refdf[,"field"])),cex=0.7)
draw_edges(edges,Xpro[,1],Xpro[,2],
           col=get_colors(nperdata)[edge_col_idx])
text(meanXpro[,1],meanXpro[,2],labels=rownames(meanXpro),
     col=get_colors(nperdata)[refdf[1:nperdata,"lang"]],
     cex=c(0.4,0.6)[1+rownames(meanXpro)%in%highlight])
dev.off()

###############################
## Visualisation of results for one language group at a time
kmax=25
lang=scale(lang)
lang = (lang+t(lang))/2
gene=scale(gene)
gene = (gene + t(gene))/2

slangscan<-Clarity_Scan(semantic_list[["lang"]],kmax=kmax)
sgenescan<-Clarity_Scan(semantic_list[["gene"]],kmax=kmax)
slangpredictsx<-lapply(semantic_list,Clarity_Predict,clist=slangscan,verbose=FALSE)
sgenepredictsx<-lapply(semantic_list,Clarity_Predict,clist=sgenescan,verbose=FALSE)

slangpredictsx_resid=lapply(slangpredictsx,Clarity_Extract,summary=I,k=15)
sgenepredictsx_resid=lapply(sgenepredictsx,Clarity_Extract,summary=I,k=15)
slangpredictsx_resid=lapply(slangpredictsx_resid,really_as_matrix)
sgenepredictsx_resid=lapply(sgenepredictsx_resid,really_as_matrix)


get_thresh<-function(r){
    diag(r)<-0
    max(abs(r))
}
as_network<-function(r,thresh){
    ret=(r>thresh)-((r< -thresh))
    diag(ret)=0
    ret=((ret+t(ret))/2)
    ret[ret<0]=-1
    ret[ret>0]=1
    ret
}
s_lang_thresh=get_thresh(slangpredictsx_resid[["lang"]])
s_gene_thresh=get_thresh(sgenepredictsx_resid[["gene"]])
slang_networks<-lapply(slangpredictsx_resid,as_network,thresh=s_lang_thresh)
sgene_networks<-lapply(sgenepredictsx_resid,as_network,thresh=s_gene_thresh)
sapply(sgene_networks,mean)
sapply(slang_networks,mean)

image(slangpredictsx_resid[["lang"]])
sapply(slangpredictsx_resid,range)[,"lang"]
sapply(sgenepredictsx_resid,range)

for(pop in ref_order) {
    print(pop)
    popsignifmat<-sapply(slang_networks,function(x)x[,pop])
    popdf=cbind(points_df2,popsignifmat)
    ## Plot it
    pdf(paste0("Pipeline results/bypop/SimplerVisualisation",pop,".pdf"),height=6,width=8)
    for(field in names(semantic_list)) {
        field_sym <- sym(field) 
        r=ggplot(data = world) +
            geom_sf(fill = "antiquewhite") +
            geom_point(data=popdf,aes(x = lon, y = lat, 
                color = group,
                alpha=factor(!!field_sym), 
                shape = factor(!!field_sym)),
                size = 1) +
            coord_sf(xlim = range(metadata[,"Coordinate 2"])+c(-5,5),
                     ylim =c(2,20), expand = FALSE) +
            scale_color_manual(values = grouping_map) +
            scale_size_manual(values = c("0" = 2, "1" = 2.5)) + # smaller for faded labels
            scale_alpha_manual(values = c("0" = 0.2, "1" = 1)) +
            scale_shape_manual(values = c("0" = 16, "1" = 17, "-1" = 4)) +  # example shapes
            geom_text_repel(data = popdf,
                            aes(x = lon, y = lat, label = label, color = group,
                                alpha = factor(!!field_sym),
                                fontface = ifelse(!!field_sym == -1, "italic", "plain")
                        ), size=2.5,
                            max.overlaps=40, seed = 42,show.legend = FALSE) +
            theme_minimal() +
            guides(alpha = "none", size = "none",shape="none") +
            labs(title = field,
                 x = "Longitude", y = "Latitude")
        print(r)
    }
    dev.off()
}


###
## Sort the function words by overall borrowing rate according to an external reference
## The last 8 (from Quantity) are 20.5% and below borrowing rate
# Tadmor, Uri. 2009. Loanwords in the world’s languages: Findings and results. In M. Haspelmath & U. Tadmor (Eds.), Loanwords in the world’s languages: A comparative handbook. (pp. 55–75). Berlin/New York: De Gruyter Mouton.

semantic_fields_sorted=c("Miscellaneous function words",
"Religion and belief","Clothing and grooming","The house","Law",
"Social and political relations","Agriculture and vegetation","Food and drink",
"Warfare and hunting","Possession","Animals","Cognition","Basic actions and technology",
"Time","Speech and language","Quantity","Emotions and values","The physical world","Motion",
"Kinship","The body","Spatial relations","Sense perception")

weightedMean<-function(xl,w){
    res=xl[[1]]*w[1]
    for(i in 2:length(w)){
        res=res + xl[[2]]*w[2]
    }
    return(res/sum(w))
}

sf1<-tail(semantic_fields_sorted,8)
sf2<-tail(head(semantic_fields_sorted,9),8)
weightings<-colMeans(semantic_table)
sf2=sf2[!is.na(weightings[sf2])]
sfmatrix1<-weightedMean(semantic_list_full[sf1],weightings[sf1])
sfmatrix2<-weightedMean(semantic_list_full[sf2],weightings[sf2])

pdf("Pipeline results/OCSEAN_SemanticDistancesGrouped.pdf",height=6,width=6)
    Clarity_Chart(sfmatrix1[rownames(gene),colnames(gene)],
                  scalefun=I,cex.axis=0.4,las=2,mar=mypar,
                  col.axis.Y=grouping_colors,col.axis.X=grouping_colors,
                  main="Highly Conserved Semantic Fields")
    Clarity_Chart(sfmatrix2[rownames(gene),colnames(gene)],
                  scalefun=I,cex.axis=0.4,las=2,mar=mypar,
                  col.axis.Y=grouping_colors,col.axis.X=grouping_colors,
                  main="Poorly Conserved Semantic Fields")
dev.off()
