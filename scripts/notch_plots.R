library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(showtext)
library(ape)
library(ggrepel)
library('Biostrings')
library(grid)
library(gridExtra)
library(cowplot)
showtext_auto()

plotNormGeneFeatures <- function(sub_gff) {
  
  GFFbyL2 <- sub_gff %>% 
    dplyr::arrange(desc(spp),L2) %>%
    dplyr::group_by(spp,L2) %>%
    dplyr::mutate(ngroup=cur_group_id()) %>%
    dplyr::ungroup()
  
  neg_str <- GFFbyL2 %>% dplyr::filter(strand=="-")
  neg_L3 <- neg_str %>% 
    dplyr::filter(type=="exon" | type =="intron" | type =="CDS") %>%
    dplyr::group_by(L2) %>%
    dplyr::mutate(start=abs(start-(max(end)))) %>%
    dplyr::mutate(end=abs(end-max(end))) %>%
    dplyr::mutate(temp=start) %>%
    dplyr::mutate(start=end) %>%
    dplyr::mutate(end=temp) %>%
    dplyr::select(-temp) %>%
    dplyr::ungroup()
  neg_L1L2 <- neg_str %>% dplyr::filter(!(type=="exon") & !(type =="intron") & !(type =="CDS"))
  neg_corr <- rbind(neg_L1L2,neg_L3) %>% dplyr::arrange(start)
  
  pos_str <- GFFbyL2 %>% dplyr::filter(strand=="+")
  
  grouped_gff <- rbind(neg_corr,pos_str)
  
  labels <- grouped_gff %>%
    dplyr::filter(type=="gene") %>%
    dplyr::select(tag,ngroup,strand,end)
  
  CDS_start <- grouped_gff %>% 
    dplyr::group_by(L2) %>%
    dplyr::filter(type=="CDS") %>% 
    dplyr::filter(start==(min(start))) %>%
    dplyr::ungroup() %>%
    dplyr::filter(start==max(start)) %>%
    dplyr::distinct(start,.keep_all = T)
  
  maxCDS <- CDS_start %>% dplyr::select(start,ngroup)
  
  reference_tran <- grouped_gff %>% dplyr::filter(ngroup==maxCDS$ngroup)
  sub_tran <- grouped_gff %>% 
    dplyr::filter(!(ngroup==maxCDS$ngroup)) %>%
    dplyr::group_by(L2) %>%
    dplyr::mutate(CDSstr=ifelse(type=="CDS",start,1e9)) %>%
    dplyr::mutate(hjust_cds=min(CDSstr)) %>%
    dplyr::mutate(shift=maxCDS$start-hjust_cds) %>%
    dplyr::mutate(start=start+(maxCDS$start-hjust_cds)) %>%
    dplyr::mutate(end=end+(maxCDS$start-hjust_cds)) %>%
    dplyr::ungroup() %>%
    dplyr::select(-CDSstr,-hjust_cds,-shift)
  
  plottable_df <- rbind(reference_tran,sub_tran)
  
  geneL <- plottable_df %>%
    dplyr::filter(type=="gene") %>%
    dplyr::mutate(len=end-start)
  
  span <- max(geneL$len)
  
  tipE <- plottable_df %>% 
    dplyr::group_by(L2) %>%
    dplyr::filter(type=="exon") %>%
    dplyr::mutate(tip=ifelse(strand=="+",max(start),min(start))) %>%
    dplyr::filter(start==tip) %>%
    dplyr::ungroup()
  
  poly<- tipE %>% 
    group_by(L2) %>%
    dplyr::mutate(xpol=list(c(end/1e3+0.015*(span/1e3),end/1e3+0.015*(span/1e3),end/1e3+0.03*(span/1e3)))) %>%
    dplyr::mutate(ypol=list(c(ngroup+0.05,ngroup-0.05,ngroup))) %>%
    dplyr::mutate(xmin=end/1e3,xmax=end/1e3+0.015*(span/1e3),ymin=ngroup-0.025,ymax=ngroup+0.025)
  
  restE <- plottable_df %>% 
    dplyr::filter(type=="exon")
  
  cdsE <- plottable_df %>% 
    dplyr::filter(type=="CDS")
  
  ids <- factor(seq(1,nrow(tipE),1))
  t2 <- data.frame(x=unlist(dplyr::pull(poly,xpol)),y=unlist(dplyr::pull(poly,ypol)),z=rep(ids,each=3))
  
  ggplot() + geom_rect(data = restE, aes(xmin = start/1e3 ,xmax = end/1e3 ,ymin=ngroup-0.25 , ymax=ngroup+0.25),fill="grey") +
    geom_rect(data = cdsE, aes(xmin = start/1e3 ,xmax = end/1e3 ,ymin=ngroup-0.25 , ymax=ngroup+0.25, fill=spp)) +
    geom_segment(data = plottable_df %>% dplyr::filter(type=="intron"), aes(x=start/1e3,xend=(start+((end-start)/2))/1e3,y=ngroup,yend=ngroup+0.25),color="darkgrey") +
    geom_segment(data = plottable_df %>% dplyr::filter(type=="intron"), aes(x=(start+((end-start)/2))/1e3,xend=end/1e3,y=ngroup+0.25,yend=ngroup),color="darkgrey") +
    #geom_rect(data = poly, aes(xmin = xmin ,xmax = xmax ,ymin=ymin , ymax=ymax),fill="black") +
    #geom_polygon(data = t2,aes(x=x,y=y,group=z),fill="black") +
    annotate('text',x=0,y=labels$ngroup+0.45,label=paste0((labels$tag)),parse=F,fontface='italic',hjust=0) +
    annotate('text',x=labels$end/1e3+0.015*(span/1e3),y=labels$ngroup,label=paste0("(",labels$strand,")"),parse=F,hjust=0) +
    theme(axis.title.y  = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          axis.line.x = element_line(),
          panel.background = element_blank(),
          legend.position="none",
          # legend.title = element_blank(),
          # legend.position = c(.95, .25),
          # legend.justification = c("right", "top"),
          # legend.box.just = "right",
          plot.title = element_text(hjust = 0,size = 12)) + xlab("Transcript length (kb)") + scale_fill_manual(values=c("pink","lightblue"))
}

getGeneFeature <- function(path,gene) {
  
  gff <- ape::read.gff(file = path) %>% 
    dplyr::mutate(source="AndersenLab")
  
  geneF <- gff %>% 
    dplyr::filter(type =="gene") %>%
    tidyr::separate(attributes, into = c("ID", "post"),sep = ";Name=") %>%
    dplyr::mutate(ID=gsub("ID=","",ID)) %>%
    dplyr::mutate(ID=gsub("gene:","",ID)) %>%
    dplyr::mutate(ID=gsub("Gene:","",ID)) %>%
    tidyr::separate(post,into=c("pre","aliases"), sep=";Alias=") %>%
    tidyr::separate(aliases,into=c("mAlias","rest"), sep=",",extra = "merge") %>%
    dplyr::mutate(L1=ID)
  
  L1 <- geneF %>% 
    dplyr::select(seqid,start,end,type,strand,ID,L1)
  
  geneAlias <- geneF %>% 
    dplyr::select(ID,mAlias) 
  
  L2 <- gff %>% 
    dplyr::filter(type=="mRNA") %>%
    tidyr::separate(attributes,into = c("ID","post"),sep=";Parent=") %>%
    dplyr::mutate(ID=gsub("ID=","",ID)) %>%
    dplyr::mutate(ID=gsub("transcript:","",ID)) %>%
    dplyr::mutate(ID=gsub("Transcript:","",ID)) %>%
    tidyr::separate(post,into = c("L1","rest"),sep=";Name=") %>%
    dplyr::mutate(L1=gsub("gene:","",L1)) %>%
    dplyr::mutate(L1=gsub("Gene:","",L1)) %>%
    dplyr::select(seqid,start,end,type,strand,ID,L1)
  
  parentL2 <- L2 %>% 
    dplyr::select(ID,L1)
  
  children <- gff %>% 
    dplyr::filter(!(type=="mRNA") & !(type=="gene")) %>%
    dplyr::filter(!(type=="three_prime_UTR") & !(type=="five_prime_UTR")) %>%
    tidyr::separate(attributes,into=c("L3","L2"),sep=";") %>%
    dplyr::mutate(L2=ifelse(is.na(L2) | grepl("Note=",L2),L3,L2)) %>%
    dplyr::mutate(L3=gsub("ID=","",L3)) %>%
    dplyr::mutate(L3=gsub("CDS:","",L3)) %>%
    dplyr::mutate(L2=gsub("Parent=","",L2)) %>%
    dplyr::mutate(L2=gsub("transcript:","",L2)) %>%
    dplyr::mutate(L2=gsub("Transcript:","",L2)) %>%
    dplyr::mutate(L3=gsub("Parent=","",L3)) %>%
    dplyr::mutate(L3=gsub("transcript:","",L3)) %>%
    dplyr::mutate(L3=gsub("Transcript:","",L3)) 
  
  if(any(grepl(",",children$L2))) {
    children <- children %>%
      tidyr::separate_rows(L2,sep = ",") %>%
      dplyr::mutate(L3=L2) %>%
      dplyr::left_join(parentL2,by = c("L2"="ID")) %>%
      dplyr::mutate(L3=paste0(L3,"_",type)) %>%
      dplyr::group_by(L3) %>%
      dplyr::mutate(L3=paste0(L3,row_number()))%>%
      dplyr::ungroup()
  } else {
    children <- children %>%
      dplyr::left_join(parentL2,by = c("L2"="ID")) 
  }
  
  parentL3 <- children %>% 
    dplyr::select(L2,L3)
  
  L3 <- children %>% 
    dplyr::rename(ID=L3) %>%
    dplyr::select(seqid,start,end,type,strand,ID,L1)
  
  tabGFF <- BiocGenerics::rbind(L1,L2,L3) %>%
    dplyr::arrange(seqid,start) %>%
    dplyr::filter(L1==gene)
  
  geneMain <- tabGFF %>% 
    dplyr::filter(type=="gene")
  
  subHead <- tabGFF %>% 
    dplyr::filter(type=="mRNA")
  
  geneHead <- geneMain %>%
    dplyr::mutate(L2=paste(subHead$ID,collapse = ",")) %>%
    tidyr::separate_rows(L2,sep = ",")
  
  remFeatures <- tabGFF %>% dplyr::filter(!(type=="gene")) %>%
    dplyr::left_join(parentL3,by=c("ID"="L3")) %>%
    dplyr::mutate(L2=ifelse(type=="mRNA",ID,L2)) 
  
  if (nrow(remFeatures %>% dplyr::filter(type=="intron")) == 0) {
    introns <- remFeatures %>% 
      dplyr::filter(type=="exon") %>% 
      dplyr::group_by(L2) %>%
      dplyr::mutate(istart=end+1,iend=lead(start)-1) %>%
      dplyr::select(-start,-end) %>%
      dplyr::rename(start=istart,end=iend) %>%
      dplyr::filter(!is.na(end)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(type="intron") %>%
      dplyr::select(seqid,start,end,type,strand,ID,L1,L2) %>%
      dplyr::mutate(ID=gsub("exon","intron",ID))
    
    finalFeatures <- rbind(geneHead,remFeatures,introns) %>%
      dplyr::group_by(L2) %>%
      dplyr::arrange(start) %>%
      dplyr::left_join(geneAlias,by = c("L1"="ID"))%>%
      dplyr::mutate(startPos=min(start)) %>%
      dplyr::mutate(end=end-min(start)) %>%
      dplyr::mutate(start=start-min(start)) 
    
  } else {
    finalFeatures <- rbind(geneHead,remFeatures) %>%
      dplyr::group_by(L2) %>%
      dplyr::arrange(start) %>%
      dplyr::left_join(geneAlias,by = c("L1"="ID")) %>%
      dplyr::mutate(startPos=min(start)) %>%
      dplyr::mutate(end=end-min(start)) %>%
      dplyr::mutate(start=start-min(start))
    
  }
  
  
  features <- dplyr::group_split(finalFeatures)
  
  return(features)
}

#dev
# plotGeneFeature <- function(sub_gff,species) {
#   if (sub_gff$strand[1]=="+") {
#     tipE <- sub_gff %>% dplyr::filter(type=="exon" & start==max(start))
#     restE <- sub_gff %>% dplyr::filter(type=="exon")
#     poly<- tipE %>% dplyr::mutate(xpol=list(c((end),(end),end+(0.1*(end-start))))) %>%
#       dplyr::mutate(ypol=list(c(1.5,0.5,1))) 
#     t2 <- data.frame(x=unlist(dplyr::pull(poly,xpol)),y=unlist(dplyr::pull(poly,ypol)))
#     p1 <- ggplot() + geom_rect(data = restE, aes(xmin = start ,xmax = end ,ymin=0.5 , ymax=1.5),fill="lightblue") +
#       #geom_rect(data = tipE, aes(xmin = start ,xmax = start+(end-start)/2 ,ymin=0.5 , ymax=1.5),fill="lightblue") +
#       geom_segment(data = sub_gff %>% dplyr::filter(type=="intron"), aes(x=start,xend=start+((end-start)/2),y=1,yend=1.5),color="darkgrey") +
#       geom_segment(data = sub_gff %>% dplyr::filter(type=="intron"), aes(x=start+((end-start)/2),xend=end,y=1.5,yend=1),color="darkgrey") +
#       geom_polygon(data = t2,aes(x=x,y=y),fill="grey") +
#       theme(axis.title.y  = element_blank(),
#             axis.text.y = element_blank(),
#             axis.ticks.y = element_blank(),
#             axis.line.y = element_blank(),
#             axis.line.x = element_line(),
#             panel.background = element_blank(),
#             plot.title = element_text(hjust = 0,size = 12)) + xlab(paste0(species," pysical position"))
#   } else {
#     tipE <- sub_gff %>% dplyr::filter(type=="exon" & start==min(start))
#     restE <- sub_gff %>% dplyr::filter(type=="exon" & !(start==min(start)))
#     poly<- tipE %>% dplyr::mutate(xpol=list(c((end-(end-start)/2),(end-(end-start)/2),start))) %>%
#       dplyr::mutate(ypol=list(c(1.5,0.5,1))) 
#     t2 <- data.frame(x=unlist(dplyr::pull(poly,xpol)),y=unlist(dplyr::pull(poly,ypol)))
#     p1 <- ggplot() + geom_rect(data = restE, aes(xmin = start ,xmax = end ,ymin=0.5 , ymax=1.5),fill="pink") +
#       geom_rect(data = tipE, aes(xmin = end-((end-start)/2) ,xmax = end ,ymin=0.5 , ymax=1.5),fill="pink") +
#       geom_segment(data = sub_gff %>% dplyr::filter(type=="intron"), aes(x=start,xend=start+((end-start)/2),y=1,yend=1.5),color="darkgrey") +
#       geom_segment(data = sub_gff %>% dplyr::filter(type=="intron"), aes(x=start+((end-start)/2),xend=end,y=1.5,yend=1),color="darkgrey") +
#       geom_polygon(data = t2,aes(x=x,y=y),fill="pink") +
#       theme(axis.title.y  = element_blank(),
#             axis.text.y = element_blank(),
#             axis.ticks.y = element_blank(),
#             axis.line.y = element_blank(),
#             axis.line.x = element_line(),
#             panel.background = element_blank(),
#             plot.title = element_text(hjust = 0,size = 12)) + xlab(paste0(species," pysical position"))
#   }
#   return(p1)
# }


#read QX1410 gff and get transcript list
curated_gff <- ape::read.gff("/gff/all_pc_genes/Curation-VF-230612.PC.clean.renamed.WB.gff3") %>%
  dplyr::filter(type=="mRNA") %>%
  tidyr::separate(attributes,into=c("preID","Other"),sep=";Parent") %>%
  dplyr::mutate(ID=gsub("ID=","",preID)) %>%
  dplyr::select(ID,seqid,start,end)
colnames(curated_gff) <- c("ID","QX1410_chr","QX1410_start","QX1410_end")

#read AF16 gff and get transcript list
af16_gff <- ape::read.gff("/gff/all_pc_genes/c_briggsae.PRJNA10731.WS280.protein_coding.gff") %>%
  dplyr::filter(type=="mRNA") %>%
  tidyr::separate(attributes,into=c("preID","Other"),sep=";Parent") %>%
  tidyr::separate(Other,into=c("Parent","Other"),sep=";Name") %>%
  dplyr::mutate(ID=gsub("ID=","",preID)) %>%
  dplyr::mutate(Gene=gsub("=Gene:","",Parent)) %>%
  dplyr::select(Gene,ID,seqid,start,end) 

#read N2 gff and get transcript list
N2_gff <- ape::read.gff("/gff/all_pc_genes/c_elegans.PRJNA13758.WS283.protein_coding.gff3") %>%
  dplyr::filter(type=="mRNA") %>%
  tidyr::separate(attributes,into=c("preID","Other"),sep=";Parent") %>%
  tidyr::separate(Other,into=c("Parent","Other"),sep=";Name") %>%
  dplyr::mutate(ID=gsub("ID=","",preID)) %>%
  dplyr::mutate(Gene=gsub("=Gene:","",Parent)) %>%
  dplyr::select(Gene,ID,seqid,start,end)

#get QX1410 protein sequence lengths
pred <- readAAStringSet("/prot/Curation-VF-230612.PC.clean.renamed.WB.prot.fa")
prednames <- names(pred)
predseq <- paste(pred)
QXdf <- data.frame(prednames,predseq) %>%
  dplyr::mutate(first=substr(predseq,1,1)) %>%
  dplyr::mutate(len=(nchar(predseq))) %>%
  dplyr::mutate(transcript=gsub(":","_",prednames)) %>%
  dplyr::select(transcript,first,len) 

#get N2 protein sequence lengths
elegansprot <- readAAStringSet("/prot/c_elegans.PRJNA13758.WS283.protein_coding.prot.fa")
N2names <- names(elegansprot)
N2seq <- paste(elegansprot)
N2df <- data.frame(N2names,N2seq) %>%
  dplyr::mutate(first=substr(N2seq,1,1)) %>%
  dplyr::mutate(len=(nchar(N2seq))) %>%
  dplyr::mutate(transcript=gsub(":","_",N2names)) %>%
  dplyr::select(transcript,first,len) 

#get AF16 protein sequence lengths
AFprot <-readAAStringSet("/prot/c_briggsae.PRJNA10731.WS280.protein_coding.prot.fa")
AFnames <- names(AFprot)
AFseq <- paste(AFprot)
AFdf <- data.frame(AFnames,AFseq) %>%
  dplyr::mutate(first=substr(AFseq,1,1)) %>%
  dplyr::mutate(len=(nchar(AFseq))) %>%
  dplyr::mutate(transcript=gsub(":","_",AFnames)) %>%
  dplyr::select(transcript,first,len) 

#read AF16 gff and get gene list
AF16_genelist <- ape::read.gff("/gff/all_pc_genes/c_briggsae.PRJNA10731.WS280.protein_coding.gff") %>%
  dplyr::filter(type=="gene") %>%
  dplyr::select(seqid,start,end, attributes) %>%
  dplyr::filter(!grepl("cb",seqid)) %>%
  tidyr::separate(attributes, into = c("pre","post"), sep = ";sequence_name=") %>%
  tidyr::separate(post, into = c("name","post2"), sep = ";biotype") %>%
  tidyr::separate(post2, into = c("pre2","aliases"),sep = ";Alias=") %>%
  tidyr::separate(aliases, into= c("alias","alt"), sep=",") %>%
  dplyr::select(-pre,-pre2,-alt)

# get AF16vQX1410 nucmer alignments
transformed_coords <- readr::read_tsv("/alignments/AF16vQX1410/AF16_transformed.tsv") %>% 
  tidyr::separate(`[TAGS]`,into = c("QX1410_chr","AF16_chr"),sep="\t") %>%
  dplyr::mutate(misplaced=ifelse(QX1410_chr==AF16_chr,"CP","MISPLACED"))

#read orthofinder results
orthogroups <- readr::read_tsv("/orthology/Orthogroups.tsv")
colnames(orthogroups) <- c("Orthogroup","AF16","QX1410","N2","JU1422")

#subset to only matching chromosome alignments
cp_coords <- transformed_coords %>% dplyr::filter(misplaced=="CP")

#write vector with N2 lin/notch gene lists
goi_class <- c("lin-12","glp-1","lag-2","apx-1","dsl-1","arg-1","epn-1","sup-17","adm-4","sel-12","hop-1","aph-1","aph-2","pen-2","lag-1","sel-8")
goi_sec <- c("R107.8","F02A9.6","Y73C8B.4","K08D9.3","W09G12.4","F31A9.3","T04C10.2","DY3.7","ZK154.7","F35H12.3","C18E3.8","VF36H2L.1","ZC434.6","T28D6.9","K08B4.1","C32A3.1") 
goi_N2 <- as.data.frame(cbind(goi_class,goi_sec)) 

#subset to lin/notch orthogroups
orthogroups_lin <- orthogroups %>%
  dplyr::filter(grepl(paste(goi_sec,collapse = "|"),N2)) %>%
  dplyr::mutate(N2=gsub("Transcript_","",N2)) %>%
  dplyr::mutate(QX1410=gsub('QX1410.g714.1','QX1410.g714.2',QX1410))

orthoTable <- orthogroups_lin %>%
  dplyr::mutate(N2 = strsplit(as.character(N2), ",")) %>%
  tidyr::unnest(N2) %>%
  dplyr::mutate(N2=trimws(N2)) %>%
  dplyr::mutate(N2_gene=substr(N2,1,nchar(N2)-2)) %>%
  dplyr::left_join(goi_N2,by=c("N2_gene"="goi_sec")) %>%
  dplyr::mutate(N2_gene=ifelse(is.na(goi_class),substr(N2_gene,1,nchar(N2_gene)-1),N2_gene)) %>%
  dplyr::select(-goi_class) %>%
  dplyr::left_join(goi_N2,by=c("N2_gene"="goi_sec")) %>%
  dplyr::mutate(QX1410_count=stringr::str_count(pattern = "QX1410", QX1410)) %>%
  dplyr::mutate(AF16_count=stringr::str_count(pattern = "Transcript_",AF16)) %>%
  dplyr::select(goi_class,QX1410_count,AF16_count) %>%
  dplyr::rename(N2=goi_class,QX1410=QX1410_count,AF16=AF16_count)

write.table(orthoTable,"orthoCount_table.tsv",quote = F,sep = '\t',row.names = F)

#get 1:1 lin/notch orthologs and calculate protein-length accuracy
one2one_lin <- orthogroups_lin %>%
  dplyr::filter(!grepl(",",AF16)&!grepl(",",QX1410)&!grepl(",",N2)) %>%
  dplyr::mutate(N2_gene=substr(N2,1,nchar(N2)-2)) %>%
  dplyr::left_join(goi_N2,by=c("N2_gene"="goi_sec")) %>%
  dplyr::mutate(N2_gene=ifelse(is.na(goi_class),substr(N2_gene,1,nchar(N2_gene)-1),N2_gene)) %>%
  dplyr::select(-goi_class) %>%
  dplyr::left_join(goi_N2,by=c("N2_gene"="goi_sec")) %>%
  dplyr::mutate(N2=paste0("Transcript_",N2)) %>%
  dplyr::left_join(AFdf,by=c("AF16"="transcript")) %>%
  dplyr::rename("AF16_length"=len,"AF16_first"=first) %>%
  dplyr::left_join(N2df,by=c("N2"="transcript")) %>%
  dplyr::rename("N2_length"=len,"N2_first"=first) %>%
  dplyr::left_join(QXdf,by=c("QX1410"="transcript")) %>%
  dplyr::rename("QX1410_length"=len,"QX1410_first"=first) %>%
  dplyr::mutate("AF_N2_ratio"=AF16_length/N2_length) %>%
  dplyr::mutate("QX_N2_ratio"=QX1410_length/N2_length) %>%
  dplyr::mutate(absdiff=abs(1-QX_N2_ratio)-abs(1-AF_N2_ratio)) %>%
  dplyr::mutate(fill=ifelse(abs(absdiff) < 0.01,"concordant",ifelse(abs(1-QX_N2_ratio) > abs(1-AF_N2_ratio),"not improved","improved"))) 

#get X:1 lin/notch orthologs
many2one_lin <- orthogroups_lin %>%
  dplyr::filter(grepl(",",AF16)|grepl(",",QX1410)|grepl(",",N2)) %>%
  separate_longer_delim(c(N2,AF16,QX1410),delim=", ") %>%
  dplyr::mutate(N2_gene=substr(N2,1,nchar(N2)-2)) %>%
  dplyr::left_join(goi_N2,by=c("N2_gene"="goi_sec")) %>%
  dplyr::mutate(N2_gene=ifelse(is.na(goi_class),substr(N2_gene,1,nchar(N2_gene)-1),N2_gene)) %>%
  dplyr::select(-goi_class) %>%
  dplyr::left_join(goi_N2,by=c("N2_gene"="goi_sec")) %>%
  dplyr::mutate(N2=paste0("Transcript_",N2)) %>%
  dplyr::left_join(AFdf,by=c("AF16"="transcript")) %>%
  dplyr::rename("AF16_length"=len,"AF16_first"=first) %>%
  dplyr::left_join(N2df,by=c("N2"="transcript")) %>%
  dplyr::rename("N2_length"=len,"N2_first"=first) %>%
  dplyr::left_join(QXdf,by=c("QX1410"="transcript")) %>%
  dplyr::rename("QX1410_length"=len,"QX1410_first"=first) %>%
  dplyr::mutate("AF_N2_ratio"=AF16_length/N2_length) %>%
  dplyr::mutate("QX_N2_ratio"=QX1410_length/N2_length) %>%
  dplyr::mutate(absdiff=abs(1-QX_N2_ratio)-abs(1-AF_N2_ratio)) %>%
  dplyr::mutate(fill=ifelse(abs(absdiff) < 0.01,"concordant",ifelse(abs(1-QX_N2_ratio) > abs(1-AF_N2_ratio),"not improved","improved"))) %>%
  dplyr::group_by(goi_class) %>%
  dplyr::arrange(desc(QX_N2_ratio)) %>%
  dplyr::distinct(goi_class,.keep_all = T) %>%
  dplyr::ungroup()

#merge one and many orthologs
ortho_lin <- rbind(one2one_lin,many2one_lin) %>%
  dplyr::filter(!(goi_class=="arg-1")) %>%
  dplyr::mutate(exprGene=c(paste0("italic(",goi_class,")")))

#plot protein-length accuracy of QX1410 lin/notch genes as a function of AF16 lin/notch protein-length accuracy values
x_limits <- c(1.56, NA)
y_limits <- c(0.1,1.6)
colors <- c("concordant"="blue","not improved"="red","improved"="green")
ggplot(ortho_lin) + 
  geom_point(aes(x=AF_N2_ratio,y=QX_N2_ratio,color=fill)) + 
  geom_abline(slope=1,intercept = 0,linetype="dashed") + 
  xlim(0.2,1.8) + 
  ylim(0.2,1.8) + 
  theme(panel.background = element_blank(),
        legend.title = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(fill=NA),
        legend.position = c(.17, .97),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.key = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size=12),
        legend.text = element_text(size=12)) + 
  geom_label_repel(aes(x=AF_N2_ratio,y=QX_N2_ratio,label=exprGene),
                   parse=T,
                   box.padding   = 0.5, 
                   point.padding = 0.5,
                   force=80,
                   xlim=x_limits,
                   ylim=y_limits,
                   max.overlaps = 40,
                   segment.color = 'grey50',
                   family='helvetica') +
  xlab("AF16 protein-length accuracy") + 
  ylab("QX1410 protein-length accuracy") + 
  scale_color_manual(values=colors)

#plot APH-1 duplication event
aph_ortho1 <- c(920817,924399)
aph_ortho2 <- c(874503,874944)
aph_qx <- c(995580,998848)
cp_coords <- transformed_coords %>% dplyr::filter(misplaced=="CP")
misplaced <- transformed_coords %>% dplyr::filter(misplaced=="MISPLACED")
I_coords <- cp_coords %>% dplyr::filter(QX1410_chr == "I")
p1 <- ggplot(I_coords) + geom_rect(aes(xmin=aph_ortho1[1],xmax=aph_ortho1[2],ymin=0,ymax=15000000),fill="lightpink") +
  geom_rect(aes(xmin=aph_ortho2[1],xmax=aph_ortho2[2],ymin=0,ymax=15000000),fill="lightpink") +
  geom_rect(aes(xmin=0,xmax=15000000,ymin=aph_qx[1],ymax=aph_qx[2]),fill="lightblue") +
  geom_rect(aes(xmin=aph_ortho1[1],xmax=aph_ortho1[2],ymin=aph_qx[1],ymax=aph_qx[2]),fill="#8077D5") +
  geom_rect(aes(xmin=aph_ortho2[1],xmax=aph_ortho2[2],ymin=aph_qx[1],ymax=aph_qx[2]),fill="#8077D5") +
  geom_segment(aes(x=`[S2]`,xend=`[E2]`,y=`[S1]`,yend=`[E1]`)) + coord_cartesian(xlim=c(870000,930000),ylim=c(9.65e5,1.025e6)) +
  annotate("text",x=aph_ortho1[1]-1.5e3,y=1.01e6,size=4,label="paste('AF16 ',italic('Cbr-aph-1'),' locus')", color="red", family="helvetica",fontface="bold", parse=T,angle=90) +
  annotate("text",x=aph_ortho2[1]-1e3,y=1.01e6,size=4,label="paste('AF16 ',italic('CBG19473'),' locus')", color="red", family="helvetica",fontface="bold", parse=T,angle=90) +
  annotate("text",x=9e5,y=1e6,size=4,label="paste('QX1410 ',italic('aph-1'),' locus')", color="blue", family="helvetica",fontface="bold", parse=T) + theme_classic() +
  xlab("AF16 Physical position") + ylab("QX1410 Physical position")

onegene <- ape::read.gff("/gff/specific_genes/AF16_aph1.gff") %>%
  dplyr::filter(type=="gene") %>%
  dplyr::select(start,end)

one <- ape::read.gff("/gff/specific_genes/AF16_aph1.gff") %>%
  dplyr::select(type,start,end,attributes) %>%
  dplyr::mutate(tair=ifelse(type=="CDS","coding_region",
                            ifelse(type=="five_prime_UTR","5' utr",
                                   ifelse(type=="three_prime_UTR","3' utr",as.character(type))))) %>%
  dplyr::filter(!(tair=="mRNA" | tair =="gene" | tair =="stop_codon" | tair =="start_codon")) 

cds <- one %>% dplyr::filter(tair=="coding_region") %>%
  tidyr::separate(attributes,into = c("ID","Parent","Other"), sep = ";",extra="merge") %>%
  dplyr::select(-ID,-type,-Other) 

other_l2l3 <-one %>% dplyr::filter(!(tair=="coding_region")) %>%
  tidyr::separate(attributes,into = c("Parent","Note"), sep = ";",extra = "merge") %>%
  dplyr::select(-Note,-type) 

one_clean <- rbind(other_l2l3,cds) %>%
  dplyr::group_by(Parent) %>%
  dplyr::mutate(intron=ifelse(any(tair=="intron"),"wi","ni")) %>%
  dplyr::ungroup()

# if (nrow(one_clean %>% dplyr::filter(intron=="ni"))>0) {
# intronless_one <- one_clean %>%
#   dplyr::filter(intron=="ni" & tair == "exon") %>% 
#   dplyr::mutate(intstart=end+1) %>%
#   dplyr::mutate(intend=lead(start)-1) %>%
#   dplyr::filter(!(is.na(intend))) %>%
#   dplyr::mutate(tair="intron") %>%
#   dplyr::select(intstart,intend,Parent,tair) %>%
#   dplyr::rename(start=intstart) %>%
#   dplyr::rename(end=intend)
# }

cds_one <- one_clean %>%
  dplyr::filter(tair=="coding_region") %>%
  dplyr::group_by(Parent) %>%
  dplyr::mutate(cdsstart=min(start)) %>%
  dplyr::mutate(cdsend=max(end)) %>%
  dplyr::distinct(cdsstart,cdsend) %>%
  dplyr::mutate(tair="ORF") %>%
  dplyr::select(cdsstart,cdsend,Parent,tair) %>%
  dplyr::ungroup() %>%
  dplyr::rename(start=cdsstart) %>%
  dplyr::rename(end=cdsend)

tair_one <- rbind(cds_one,one_clean %>% dplyr::select(-intron))  %>%
  dplyr::arrange(Parent) %>%
  dplyr::select(tair,start,end,Parent) %>%
  dplyr::rename(type=tair) %>%
  dplyr::arrange(start) %>%
  dplyr::mutate(start=as.numeric(start)) %>%
  dplyr::mutate(end=as.numeric(end)) %>%
  tidyr::unite("coordinates",start,end,sep = "-")

tair_one_clean <- tair_one%>%
  tidyr::separate(coordinates, into = c("start","end"), sep="-") %>%
  dplyr::filter(!(type=="ORF")) %>%
  dplyr::filter(!(type=="exon")) %>%
  dplyr::mutate(start=as.numeric(start)) %>%
  dplyr::mutate(end=as.numeric(end)) %>%  
  dplyr::group_by(Parent) %>%
  dplyr::mutate(ngroup=cur_group_id()*1.5) %>%
  dplyr::ungroup()

tips <-tair_one_clean %>% 
  dplyr::filter(!(type=="intron")) %>%
  dplyr::group_by(Parent) %>%
  dplyr::arrange(desc(start)) %>%
  dplyr::filter(row_number()==n()) %>%
  dplyr::mutate(xpol=list(c((end-(end-start)/2),(end-(end-start)/2),start))) %>%
  dplyr::mutate(ngroup=ifelse(ngroup==4.5,3,ngroup)) %>%
  dplyr::mutate(ypol=list(c(ngroup+0.5,ngroup-0.5,ngroup))) %>%
  dplyr::ungroup() 

ids <- factor(c("t1"))
t2 <- data.frame(x=unlist(dplyr::pull(tips %>% dplyr::filter(Parent=="Parent=Transcript:CBG25153.1"),xpol)),y=unlist(dplyr::pull(tips %>% dplyr::filter(Parent=="Parent=Transcript:CBG25153.1"),ypol)),z=rep(ids,each=3))

p2<-ggplot() + geom_rect(data = tair_one_clean %>% 
                           dplyr::filter(Parent=="Parent=Transcript:CBG25153.1") %>% 
                           dplyr::filter(!(type=="intron")) %>%
                           dplyr::group_by(Parent) %>%
                           dplyr::arrange(desc(start)) %>%
                           dplyr::filter(!(row_number()==n())) %>%
                           dplyr::ungroup(), aes(xmin = start ,xmax = end+1 ,ymin =ngroup-0.5 , ymax= ngroup+0.5),fill="pink") +
  geom_rect(data=tair_one_clean %>% 
              dplyr::filter(Parent=="Parent=Transcript:CBG25153.1") %>%
              dplyr::filter(!(type=="intron")) %>%
              dplyr::group_by(Parent) %>%
              dplyr::arrange(desc(start)) %>%
              dplyr::filter(row_number()==n()) %>%
              dplyr::ungroup(), aes(xmin=end+1,xmax=end-(end-start)/2,ymin = ngroup-0.5 , ymax=ngroup+0.5),fill="pink") +
  geom_polygon(data = t2,aes(x=x,y=y,group=z),fill="pink") +
  geom_segment(data = tair_one_clean %>% dplyr::filter(Parent=="Parent=Transcript:CBG25153.1") %>% dplyr::filter(type=="intron"), aes(x=start,xend=start+((end-start)/2),y=ngroup,yend=ngroup+0.5),color="darkgrey") +
  geom_segment(data = tair_one_clean %>% dplyr::filter(Parent=="Parent=Transcript:CBG25153.1") %>% dplyr::filter(type=="intron"), aes(x=start+((end-start)/2),xend=end,y=ngroup+0.5,yend=ngroup),color="darkgrey") + 
  #annotate("text",x=onegene$start[1]+1.8e3,y=3.55,label="paste('AF16 ',italic('Cbr-aph-1'),' locus')",parse=T,color="red",family="helvetica") +         
  theme(axis.title  = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_line(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0,size = 12)) 

ids <- factor(c("t1"))
t2 <- data.frame(x=unlist(dplyr::pull(tips %>% dplyr::filter(Parent=="Parent=Transcript:CBG19473.1"),xpol)),y=unlist(dplyr::pull(tips %>% dplyr::filter(Parent=="Parent=Transcript:CBG19473.1"),ypol)),z=rep(ids,each=3))

p3<-ggplot() + geom_rect(data = tair_one_clean %>% 
                           dplyr::filter(Parent=="Parent=Transcript:CBG19473.1") %>% 
                           dplyr::filter(!(type=="intron")) %>%
                           dplyr::group_by(Parent) %>%
                           dplyr::arrange(desc(start)) %>%
                           dplyr::filter(!(row_number()==n())) %>%
                           dplyr::ungroup(), aes(xmin = start ,xmax = end ,ymin =ngroup-0.5 , ymax= ngroup+0.5),fill="pink") +
  geom_rect(data=tair_one_clean %>% 
              dplyr::filter(Parent=="Parent=Transcript:CBG19473.1") %>%
              dplyr::filter(!(type=="intron")) %>%
              dplyr::group_by(Parent) %>%
              dplyr::arrange(desc(start)) %>%
              dplyr::filter(row_number()==n()) %>%
              dplyr::ungroup(), aes(xmin=end+1,xmax=end-(end-start)/2,ymin = ngroup-0.5 , ymax=ngroup+0.5),fill="pink") +
  geom_polygon(data = t2,aes(x=x,y=y,group=z),fill="pink") +
  geom_segment(data = tair_one_clean %>% dplyr::filter(Parent=="Parent=Transcript:CBG19473.1") %>% dplyr::filter(type=="intron"), aes(x=start,xend=start+((end-start)/2),y=ngroup,yend=ngroup+0.5),color="darkgrey") +
  geom_segment(data = tair_one_clean %>% dplyr::filter(Parent=="Parent=Transcript:CBG19473.1") %>% dplyr::filter(type=="intron"), aes(x=start+((end-start)/2),xend=end,y=ngroup+0.5,yend=ngroup),color="darkgrey") + 
  #annotate("text",x=onegene$end[2]-1.95e3,y=2.05,label="paste('AF16 ',italic('CBG19473'),' locus')",parse=T,color="red",family="helvetica") +
  geom_rect(aes(xmin=874500,xmax=874506,ymin=1,ymax=2),fill="pink")+
  theme(axis.title  = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_line(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0,size = 12))  + xlim(onegene$start[2]-(onegene$end[1]-onegene$start[1]),onegene$end[2])#+ xlab("AF16 Physical position")

plot_grid(p2,p3,nrow=2,ncol=1)


onegene <- ape::read.gff("/gff/specific_genes/QX1410_aph1.gff") %>%
  dplyr::filter(type=="gene") %>%
  dplyr::select(start,end)

one <- ape::read.gff("/gff/specific_genes/QX1410_aph1.gff") %>%
  dplyr::select(type,start,end,attributes) %>%
  dplyr::mutate(tair=ifelse(type=="CDS","coding_region",
                            ifelse(type=="five_prime_UTR","5' utr",
                                   ifelse(type=="three_prime_UTR","3' utr",as.character(type))))) %>%
  dplyr::filter(!(tair=="mRNA" | tair =="gene" | tair =="stop_codon" | tair =="start_codon")) %>%
  tidyr::separate(attributes,into = c("ID","Parent","Other"), sep = ";",extra="merge") %>%
  dplyr::select(-ID,-type,-Other) 


one_clean <- one %>%
  dplyr::group_by(Parent) %>%
  dplyr::mutate(intron=ifelse(any(tair=="intron"),"wi","ni")) %>%
  dplyr::ungroup()

if (nrow(one_clean %>% dplyr::filter(intron=="ni"))>0) {
  intronless_one <- one_clean %>%
    dplyr::filter(intron=="ni" & tair == "exon") %>%
    dplyr::mutate(intstart=end+1) %>%
    dplyr::mutate(intend=lead(start)-1) %>%
    dplyr::filter(!(is.na(intend))) %>%
    dplyr::mutate(tair="intron") %>%
    dplyr::select(intstart,intend,Parent,tair) %>%
    dplyr::rename(start=intstart) %>%
    dplyr::rename(end=intend)
}

cds_one <- one_clean %>%
  dplyr::filter(tair=="coding_region") %>%
  dplyr::group_by(Parent) %>%
  dplyr::mutate(cdsstart=min(start)) %>%
  dplyr::mutate(cdsend=max(end)) %>%
  dplyr::distinct(cdsstart,cdsend) %>%
  dplyr::mutate(tair="ORF") %>%
  dplyr::select(cdsstart,cdsend,Parent,tair) %>%
  dplyr::ungroup() %>%
  dplyr::rename(start=cdsstart) %>%
  dplyr::rename(end=cdsend)

tair_one <- rbind(cds_one,one_clean %>% dplyr::select(-intron),intronless_one)  %>%
  dplyr::arrange(Parent) %>%
  dplyr::select(tair,start,end,Parent) %>%
  dplyr::rename(type=tair) %>%
  dplyr::arrange(start) %>%
  dplyr::mutate(start=start) %>%
  dplyr::mutate(end=end) %>%
  tidyr::unite("coordinates",start,end,sep = "-")

tair_one_clean <- tair_one%>%
  tidyr::separate(coordinates, into = c("start","end"), sep="-") %>%
  dplyr::filter(!(type=="ORF")) %>%
  dplyr::filter(!(type=="exon")) %>%
  dplyr::mutate(start=as.numeric(start)) %>%
  dplyr::mutate(end=as.numeric(end)) %>%  
  dplyr::group_by(Parent) %>%
  dplyr::mutate(ngroup=cur_group_id()*1.5) %>%
  dplyr::ungroup()

tips <-tair_one_clean %>% 
  dplyr::filter(!(type=="intron")) %>%
  dplyr::group_by(Parent) %>%
  dplyr::arrange(desc(start)) %>%
  dplyr::filter(row_number()==n()) %>%
  dplyr::mutate(xpol=list(c((end-(end-start)/2),(end-(end-start)/2),start))) %>%
  dplyr::mutate(ngroup=ifelse(ngroup==4.5,3,ngroup)) %>%
  dplyr::mutate(ypol=list(c(ngroup+0.5,ngroup-0.5,ngroup))) %>%
  dplyr::ungroup()

ids <- factor(c("t1"))
t2 <- data.frame(x=unlist(dplyr::pull(tips,xpol)),y=unlist(dplyr::pull(tips,ypol)),z=rep(ids,each=3))

p4 <- ggplot() + geom_rect(data = tair_one_clean %>% 
                             dplyr::filter(!(type=="intron")) %>%
                             dplyr::group_by(Parent) %>%
                             dplyr::arrange(desc(start)) %>%
                             dplyr::filter(!(row_number()==n())) %>%
                             dplyr::ungroup(), aes(xmin = start ,xmax = end+1 ,ymin =ngroup-0.5 , ymax= ngroup+0.5),fill="lightblue") +
  geom_rect(data=tair_one_clean %>% 
              dplyr::filter(!(type=="intron")) %>%
              dplyr::group_by(Parent) %>%
              dplyr::arrange(desc(start)) %>%
              dplyr::filter(row_number()==n()) %>%
              dplyr::ungroup(), aes(xmin=end+1,xmax=end-(end-start)/2,ymin = ngroup-0.5 , ymax=ngroup+0.5),fill="lightblue") +
  geom_polygon(data = t2,aes(x=x,y=y,group=z),fill="lightblue") +
  geom_segment(data = tair_one_clean %>% dplyr::filter(type=="intron"), aes(x=start,xend=start+((end-start)/2),y=ngroup,yend=ngroup+0.5),color="darkgrey") +
  geom_segment(data = tair_one_clean %>% dplyr::filter(type=="intron"), aes(x=start+((end-start)/2),xend=end,y=ngroup+0.5,yend=ngroup),color="darkgrey") + 
  #annotate("text",x=onegene$start+1.68e3,y=2.05,label="paste('QX1410 ',italic('aph-1'),' locus')",parse=T,color="blue",family="helvetica") +         
  theme(axis.title.y  = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_line(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0,size = 12)) + xlab("Physical position") + xlim(NA,998835)

blank <- grid.rect(gp=gpar(col="white"))
plot_grid(p1, arrangeGrob(arrangeGrob(p2,
                                      arrangeGrob(p3,top = textGrob(expression(italic('CBG19473')),vjust = 1,hjust=5.75,gp=gpar(col="red"))),
                                      top=textGrob(expression(paste0("AF16",italic('AF16 Cbr-aph-1'))), vjust = 1, hjust=4, gp=gpar(col="red")),
                                      nrow=2,ncol=1,heights = c(1.38,2)),
                          arrangeGrob(p4,top = textGrob(expression(italic('QX1410 Cbr-aph-1')),vjust = 1,hjust=3.35,gp=gpar(col="blue"))),
                          ncol = 1,nrow=2,heights = c(3.2,2)),
          ncol = 1,nrow = 2,rel_heights = c(1.5,0.6),labels = c("A","B"),align = "hv")




#plot LIN-12 duplication event                                                               
lin12_1 <- c(9950920,9959304)
lin12_2 <- c(10474323,10481764)
lin12_3 <- c(10458881,10466015)
lin12a <- c(10607219,10617555)
lin12b <- c(10588085,10598994)
III_coords <- cp_coords %>% dplyr::filter(QX1410_chr == "III")
colors <- c("high"="black","low"="red")
plindup<-ggplot(III_coords %>% dplyr::mutate(Identity=ifelse(`[% IDY]`>99,"high","low")) %>% dplyr::filter(Identity=="high")) +
  geom_rect(aes(xmin=0,xmax=15000000,ymin=lin12a[1],ymax=lin12a[2]),fill="lightblue") +
  geom_rect(aes(xmin=0,xmax=15000000,ymin=lin12b[1],ymax=lin12b[2]),fill="lightblue") +
  geom_rect(aes(xmin=lin12_1[1],xmax=lin12_1[2],ymin=0,ymax=15000000),fill="lightpink") +
  geom_rect(aes(xmin=lin12_2[1],xmax=lin12_2[2],ymin=0,ymax=15000000),fill="lightpink") +
  geom_rect(aes(xmin=lin12_3[1],xmax=lin12_3[2],ymin=0,ymax=15000000),fill="lightpink") +
  geom_rect(aes(xmin=lin12_1[1],xmax=lin12_1[2],ymin=lin12a[1],ymax=lin12a[2]),fill="#8077D5") +
  geom_rect(aes(xmin=lin12_2[1],xmax=lin12_2[2],ymin=lin12b[1],ymax=lin12b[2]),fill="#8077D5") +
  geom_rect(aes(xmin=lin12_1[1],xmax=lin12_1[2],ymin=lin12a[1],ymax=lin12a[2]),fill="#8077D5") +
  geom_segment(aes(x=`[S2]`,xend=`[E2]`,y=`[S1]`,yend=`[E1]`)) + coord_cartesian(xlim=c(9945000,10490000),ylim=c(10520000,10670000)) + scale_color_manual(values=colors) +
  annotate("text",x=lin12_1[1]-6.5e3,y=10.55e6,size=4,label="paste('AF16 ',italic('Cbr-lin-12.1'),' locus')", color="red", family="helvetica",fontface="bold", parse=T,angle=90) +
  annotate("text",x=lin12_2[1]+1.3e4,y=10.55e6,size=4,label="paste('AF16 ',italic('Cbr-lin-12.2'),' locus')", color="red", family="helvetica",fontface="bold", parse=T,angle=90)+
  annotate("text",x=lin12_3[1]-6.5e3,y=10.55e6,size=4,label="paste('AF16 ',italic('Cbr-lin-12.3'),' locus')", color="red", family="helvetica",fontface="bold", parse=T,angle=90) +
  annotate("text",x=10220000,y=10601500,size=4,label="paste('QX1410 ',italic('Cbr-lin-12.1'),' locus')", color="blue", family="helvetica",fontface="bold", parse=T) +
  annotate("text",x=10220000,y=10620000,size=4,label="paste('QX1410 ',italic('Cbr-lin-12.2'),' locus')", color="blue", family="helvetica",fontface="bold", parse=T) +
  theme_classic() +
  xlab("AF16 Physical position") + ylab("QX1410 Physical position")


path <- "/gff/all_pc_genes/Curation-VF-230612.PC.clean.renamed.WB.gff3"
gene <- "QX1410.13860"
features <- getGeneFeature(path,gene)

path <- "/gff/all_pc_genes/Curation-VF-230612.PC.clean.renamed.WB.gff3"
gene <- "QX1410.13858"
features2 <- getGeneFeature(path,gene)

path <- "/gff/all_pc_genes/c_briggsae.PRJNA10731.WS280.protein_coding.gff"
gene <- "WBGene00029134"
features3 <- getGeneFeature(path,gene)

path <- "/gff/all_pc_genes/c_briggsae.PRJNA10731.WS280.protein_coding.gff"
gene <- "WBGene00029035"
features4 <- getGeneFeature(path,gene)

path <- "/gff/all_pc_genes/c_briggsae.PRJNA10731.WS280.protein_coding.gff"
gene <- "WBGene00029037"
features5 <- getGeneFeature(path,gene)




sub_gff <- rbind(features[[1]] %>% dplyr::mutate(tag="QX1410 Cbr-lin-12.1") %>% dplyr::mutate(spp="QX1410"),
                 features2[[1]]%>% dplyr::mutate(tag="QX1410 Cbr-lin-12.2") %>% dplyr::mutate(spp="QX1410"),
                 features3[[1]]%>% dplyr::mutate(tag="AF16 Cbr-lin-12.1") %>% dplyr::mutate(spp="AF16"),
                 features4[[1]]%>% dplyr::mutate(tag="AF16 Cbr-lin-12.2") %>% dplyr::mutate(spp="AF16"),
                 features5[[1]]%>% dplyr::mutate(tag="AF16 Cbr-lin-12.3") %>% dplyr::mutate(spp="AF16"))

plingene <- plotNormGeneFeatures(sub_gff)


plot_grid(plindup,plingene,nrow=2,ncol=1,rel_heights = c(1,0.45),labels=c("A","B"))


#plot apx-1 orthologs
path <- "/gff/all_pc_genes/Curation-VF-230612.PC.clean.renamed.WB.gff3"
gene <- "QX1410.714"
apxqx_features <- getGeneFeature(path,gene)

path <- "/gff/all_pc_genes/c_briggsae.PRJNA10731.WS280.protein_coding.gff"
gene <- "WBGene00040504"
apxaf_features <- getGeneFeature(path,gene)

sub_gff <- rbind(apxqx_features[[1]] %>% dplyr::mutate(tag="QX1410 Cbr-apx-1.1") %>% dplyr::mutate(spp="QX1410"),
                 apxqx_features[[2]] %>% dplyr::mutate(tag="QX1410 Cbr-apx-1.2") %>% dplyr::mutate(spp="QX1410"),
                 apxaf_features[[1]]%>% dplyr::mutate(tag="AF16 Cbr-apx-1") %>% dplyr::mutate(spp="AF16"))

papxgene <- plotNormGeneFeatures(sub_gff)

#plot DSL-1 paralogs
dsl1_1 <- c(5052340,5053049)
dsl1_2 <- c(5059716,5060425)
dsl1_3 <- c(5067095,5068119)
dsl1_4 <- c(15344089,15344507)

af_dsl1 <- c(5556658,5558766)
IV_coords <- cp_coords %>% dplyr::filter(QX1410_chr == "IV")
pdsl13<-ggplot(IV_coords) + 
  geom_segment(aes(x=`[S2]`,xend=`[E2]`,y=`[S1]`,yend=`[E1]`)) + coord_cartesian(xlim=c(5100000,5700000),ylim=c(5000000,5100000)) +
  geom_rect(aes(xmin=0,xmax=15000000,ymin=dsl1_1[1],ymax=dsl1_1[2]),fill="lightblue") + 
  geom_rect(aes(xmin=0,xmax=15000000,ymin=dsl1_2[1],ymax=dsl1_2[2]),fill="lightblue")+
  geom_rect(aes(xmin=0,xmax=15000000,ymin=dsl1_3[1],ymax=dsl1_3[2]),fill="lightblue") +
  geom_rect(aes(xmin=af_dsl1[1],xmax=af_dsl1[2],ymin=0,ymax=15000000),fill="lightpink") + 
  theme_classic() +xlab("AF16 Physical position") + ylab("QX1410 Physical position") +
  geom_rect(aes(xmin=af_dsl1[1],xmax=af_dsl1[2],ymin=dsl1_3[1],ymax=dsl1_3[2]),fill="#8077D5") +
  annotate("text",x=af_dsl1[1]-8e3,y=5027500,size=4,label="paste('AF16 ',italic('Cbr-dsl-1'),' locus')", color="red", family="helvetica",fontface="bold", parse=T,angle=90) +
  annotate("text",x=af_dsl1[1]-1e5,y=dsl1_1[2]+1e3,size=4,label="paste('QX1410 ',italic('Cbr-dsl-1.3'),' locus')", color="blue", family="helvetica",fontface="bold", parse=T) +
  annotate("text",x=af_dsl1[1]-1e5,y=dsl1_2[2]+1e3,size=4,label="paste('QX1410 ',italic('Cbr-dsl-1.2'),' locus')", color="blue", family="helvetica",fontface="bold", parse=T)+
  annotate("text",x=af_dsl1[1]-1e5,y=dsl1_3[2]+1e3,size=4,label="paste('QX1410 ',italic('Cbr-dsl-1.1'),' locus')", color="blue", family="helvetica",fontface="bold", parse=T)

#plot DSL-1d 
pseudo <- c(15883646,15884013)
pdsl4<-ggplot(IV_coords) + 
  geom_rect(aes(xmin=pseudo[1],xmax=pseudo[2],ymin=0,ymax=20000000),fill="lightpink") + 
  geom_rect(aes(xmin=0,xmax=20000000,ymin=dsl1_4[1],ymax=dsl1_4[2]),fill="lightblue") +
  geom_rect(aes(xmin=pseudo[1],xmax=pseudo[2],ymin=dsl1_4[1],ymax=dsl1_4[2]),fill="#8077D5")+
  geom_segment(aes(x=`[S2]`,xend=`[E2]`,y=`[S1]`,yend=`[E1]`)) + coord_cartesian(xlim=c(15880000,15885000),ylim=c(15340000,15346000)) +
  annotate("text",x=pseudo[1]-1e2,y=dsl1_4[1]-2e3,size=4,label="paste('AF16 ',italic('CBG13949'),' locus (pseudogene)')", color="red", family="helvetica",fontface="bold", parse=T,angle=90) +
  annotate("text",x=pseudo[1]-1.2e3,y=dsl1_4[2]+1e2,size=4,label="paste('QX1410 ',italic('Cbr-dsl-1.4'),' locus')", color="blue", family="helvetica",fontface="bold", parse=T) +
  theme_classic() +xlab("AF16 Physical position") 

plot_grid(pdsl13,pdsl4,nrow=1,ncol=2,labels=c("A","B"))
