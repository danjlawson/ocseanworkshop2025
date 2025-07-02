library(RColorBrewer)
library(colorspace)

# Darken all colors by a factor
darken_colors <- function(cols, amount = 0.3) {
  darken(cols, amount = amount)
}
get_colors<-function(n,darken_by=0.2){
    exclude_palettes <- c("Pastel1", "Pastel2", "Set3")
    qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual' & 
                                 !(rownames(brewer.pal.info) %in% exclude_palettes), ]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))[1:n]
    darken_colors(col_vector, amount = darken_by)
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
get_means<-function(X,refdf){
    ulang=levels(refdf[,"lang"])
    ret=sapply(ulang,function(x){
        colMeans(X[refdf[,"lang"]==x,])
    })
    t(ret)
}


    
## Currently requires functions from 05-clarity.R!!

## Read the distances of all semantic fields (and the list of all semantic fields to get that)
semantic_df=as.data.frame(read_excel("semantic_field_list.xlsx"))
semantic_fields=semantic_df[,"Semantic_Field"]

semantic_list_full=lapply(semantic_fields,function(x){
    read_excel_as(paste0("Pipeline\ results/Semantic/matched_linguistic_semantic_field_",x,".xlsx"))
})
names(semantic_list_full)=semantic_fields

## Keep only chosen semantic fields and languages 
keep_semantic_fields=names(semantic_list_full)[!names(semantic_list_full)%in%c("Law","Religion and belief")]
keep_languages=rownames(semantic_list_full[[1]])[!rownames(semantic_list_full[[1]])%in%c("Sangil","Ati")]

## Read the raw data for count information
wordlist<-read.table("Pipeline\ results/OCSEAN_processed_joineddata_semantic.tsv",sep="\t",header=T)
semantic_full_table=table(wordlist[,c("Language","Semantic_Field")])
semantic_table=semantic_full_table[keep_languages,keep_semantic_fields]

## Sequentially remove row/columns that to ensure that all semantic fields/languages have a minimum number of concepts
print(dim(semantic_table))
min_val=5
while(any(semantic_table<min_val)){
    worst_cols=sort(colSums(semantic_table<min_val),decreasing=T)
    worst_rows=sort(rowSums(semantic_table<min_val),decreasing=T)
    if(max(worst_cols)>max(worst_rows)){
        worst_cols=worst_cols[worst_cols==max(worst_cols)]
        semantic_table=semantic_table[,!colnames(semantic_table)%in%names(worst_cols)]
    }else{
        worst_rows=worst_rows[worst_rows==max(worst_rows)]
        semantic_table=semantic_table[!rownames(semantic_table)%in%names(worst_rows),]
    }
    print(dim(semantic_table))
}
semantic_table

thm=heatmap(semantic_table,keep.dendro=T)
ling_order=labels(thm$Rowv)
semantic_order=labels(thm$Colv)
pdf("Pipeline results/SemanticConceptCountTable.pdf",height=12,width=6)
Clarity_Chart(t(semantic_table[ling_order,semantic_order]),scalefun=I,text=T,las=2,
              main="Semantic concept counts",mar=c(8,8,3,1))
dev.off()


keep_semantic_fields=colnames(semantic_table)
keep_languages=rownames(semantic_table)
keep_languages%in%rownames(gene)
keep_languages=rownames(gene)[rownames(gene) %in% keep_languages]

semantic_list=semantic_list_full[keep_semantic_fields]
semantic_list=lapply(semantic_list,function(x)x[keep_languages,keep_languages])
semantic_list=c(list("gene"=gene[keep_languages,keep_languages]),semantic_list)
semantic_list=c(list("lang"=lang[keep_languages,keep_languages]),semantic_list)
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
plot_labels=paste(rep(names(semantic_list),each=nperdata),
                  rep(rownames(semantic_list[[1]]),by=ndata),
                  sep="\n")

## Make a version that is a vector for each semantic field
flat_data=as.data.frame(sapply(semantic_list,function(x)as.numeric(x)))
plot(flat_data)
heatmap(cor(flat_data),symm=TRUE,scale="none")

heatmap(cor(semantic_table))

##scaled_semantic_data<-scale(semantic_data)
mysvd<-svd(semantic_data)
d<-5 #Assessed from plot(mysvd$d)

## Simple exploratory plots
X<-mysvd$v[,1:d] %*% sqrt(diag(mysvd$d[1:d]))
rownames(X)<-colnames(semantic_data)
i=1;j=2
plot(X[,i],X[,j],type="n")
#text(X[,i],X[,j],labels=refdf[,"lang"],col=get_colors(5)[refdf[,"field"]],cex=0.5)
text(X[,i],X[,j],labels=refdf[,"lang"],col=get_colors(nperdata)[refdf[,"lang"]],cex=0.5)


meanX=get_means(X,refdf)

semantic_d_lists=get_dists(X,refdf[,"lang"],refdf[,"field"])
thresh=0.8
plot(sort(unlist(semantic_d_lists))) ## threshold chosen from empirical elbow
abline(h=thresh)
edge_table<-get_edge_list(semantic_d_lists,thresh)
edges<-get_edges(edge_table,X,refdf)

edge_col_idx=sapply(edges[,"lang"],function(x)which(levels(refdf[,"lang"])==x))
highlight=unique(names(edge_col_idx))

leglocs=c("bottomright","bottom","topleft","topleft","topleft","topleft","topleft","topleft","topleft","topleft")
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
     col=get_colors(nperdata)[refdf[1:nperdata,"lang"]],cex=c(0.4,0.6)[1+rownames(meanX)%in%highlight])
    ii=ii+1
}
dev.off()

semantic_d=get_mean_matrix(semantic_d_lists)
myorder=c("The body", "Agriculture and vegetation", "lang", "Motion", "Food and drink","gene", "Basic actions and technology")
#thm=heatmap(semantic_d,symm=TRUE,scale="none",keep.dendro=TRUE)
#myorder=labels(thm$Rowv)

Clarity_Chart(semantic_d[myorder,myorder],scalefun=I,text=T,las=2,main="Average distance moved between...",mar=c(8,8,3,1))

