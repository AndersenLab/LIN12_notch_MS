library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(ape)
library('Biostrings')
library(grid)
library(gridExtra)
library(cowplot)
library(showtext)
library(ggtree)
showtext_auto()
showtext_opts(dpi=600)

plotNormGeneFeatures <- function(sub_gff,cust) {
  
  GFFbyL2 <- sub_gff %>% 
    dplyr::arrange(desc(spp),L2) %>%
    dplyr::group_by(spp,L2) %>%
    dplyr::mutate(ngroup=cur_group_id()) %>%
    dplyr::ungroup()
  
  if (cust == T) {
    groups <- data.frame(tag=c("JU1422 010.g27410","QX1410 Cbr-lin-12.1"),ngroup=c(1,2))
    #groups <- data.frame(tag=c("JU1422 010.g27411","QX1410 Cbr-lin-12.2"),ngroup=c(1,2))
    #groups <- data.frame(tag=c("JU1422 010.g27411","JU1422 010.g27410","QX1410 Cbr-lin-12.2","QX1410 Cbr-lin-12.1","N2 lin-12"),ngroup=c(1,2,3,4,5))
    
    GFFbyL2 <- GFFbyL2 %>%
      dplyr::select(-ngroup) %>%
      dplyr::left_join(groups,by="tag")
  }
  
  neg_str <- GFFbyL2 %>% dplyr::filter(strand=="-")
  neg_L3 <- neg_str %>% 
    dplyr::filter(type=="exon" | type =="intron" | type =="CDS" | type=="DELETION") %>%
    dplyr::group_by(L2) %>%
    dplyr::mutate(start=abs(start-(max(end)))) %>%
    dplyr::mutate(end=abs(end-max(end))) %>%
    dplyr::mutate(temp=start) %>%
    dplyr::mutate(start=end) %>%
    dplyr::mutate(end=temp) %>%
    dplyr::select(-temp) %>%
    dplyr::ungroup()
  neg_L1L2 <- neg_str %>% dplyr::filter(!(type=="exon") & !(type =="intron") & !(type =="CDS") & !(type=="DELETION"))
  neg_corr <- rbind(neg_L1L2,neg_L3) %>% dplyr::arrange(start)
  
  pos_str <- GFFbyL2 %>% dplyr::filter(strand=="+")
  
  grouped_gff <- rbind(neg_corr,pos_str)
  
  labels <- grouped_gff %>%
    dplyr::filter(type=="gene") %>%
    dplyr::select(tag,ngroup,strand,end)
  
  # CDS_start <- grouped_gff %>% 
  #   dplyr::group_by(L2) %>%
  #   dplyr::filter(type=="CDS") %>% 
  #   dplyr::filter(start==(min(start))) %>%
  #   dplyr::ungroup() %>%
  #   dplyr::filter(start==max(start)) %>%
  #   dplyr::distinct(start,.keep_all = T)
  # 
  # maxCDS <- CDS_start %>% dplyr::select(start,ngroup)
  # 
  # reference_tran <- grouped_gff %>% dplyr::filter(ngroup==maxCDS$ngroup)
  # sub_tran <- grouped_gff %>% 
  #   dplyr::filter(!(ngroup==maxCDS$ngroup)) %>%
  #   dplyr::group_by(L2) %>%
  #   dplyr::mutate(CDSstr=ifelse(type=="CDS",start,1e9)) %>%
  #   dplyr::mutate(hjust_cds=min(CDSstr)) %>%
  #   dplyr::mutate(shift=maxCDS$start-hjust_cds) %>%
  #   dplyr::mutate(start=start+(maxCDS$start-hjust_cds)) %>%
  #   dplyr::mutate(end=end+(maxCDS$start-hjust_cds)) %>%
  #   dplyr::ungroup() %>%
  #   dplyr::select(-CDSstr,-hjust_cds,-shift)
  
  # plottable_df <- rbind(reference_tran,sub_tran)
  plottable_df <- grouped_gff
  
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
  
  delE <- plottable_df %>%
    dplyr::filter(type=="DELETION")
  
  ids <- factor(seq(1,nrow(tipE),1))
  t2 <- data.frame(x=unlist(dplyr::pull(poly,xpol)),y=unlist(dplyr::pull(poly,ypol)),z=rep(ids,each=3))
  
  ggplot() +
    geom_segment(data = plottable_df %>% dplyr::filter(type=="intron"), aes(x=(start/1e3)-(20/1e3),xend=((start+((end-start)/2))/1e3)+(1/1e3),y=ngroup,yend=ngroup+0.25),color="darkgrey",linewidth=0.3) +
    geom_segment(data = plottable_df %>% dplyr::filter(type=="intron"), aes(x=((start+((end-start)/2))/1e3)-(1/1e3),xend=(end/1e3)+(20/1e3),y=ngroup+0.25,yend=ngroup),color="darkgrey",linewidth=0.3) +
    geom_rect(data = restE, aes(xmin = start/1e3 ,xmax = end/1e3 ,ymin=ngroup-0.25 , ymax=ngroup+0.25),fill="grey") +
    geom_rect(data = cdsE, aes(xmin = start/1e3 ,xmax = end/1e3 ,ymin=ngroup-0.25 , ymax=ngroup+0.25, fill=spp)) +
    geom_rect(data = delE, aes(xmin = start/1e3 ,xmax = end/1e3 ,ymin=ngroup-0.35 , ymax=ngroup-0.40),fill="black") +
    #geom_rect(data = poly, aes(xmin = xmin ,xmax = xmax ,ymin=ymin , ymax=ymax),fill="black") +
    #geom_polygon(data = t2,aes(x=x,y=y,group=z),fill="black") +
    annotate('text',x=0,y=labels$ngroup+0.35,label=paste0((labels$tag)),parse=F,fontface='italic',hjust=0,size=8/.pt) +
    #geom_text(data = labels, aes(label=paste0("(",strand,")"),x=(end/1e3)+0.015*(span/1e3),y=ngroup,hjust=0),size=3)+
    theme(axis.title.y  = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          axis.line.x = element_line(),
          panel.background = element_blank(),
          legend.position="none",
          axis.title = element_text(size = 7.5),
          axis.text = element_text(size = 6)) + 
    xlab("Transcript length (kb)") + scale_fill_manual(values=c("lightgreen","lightblue"))
}

getGeneFeature <- function(path,gene,dels) {
  
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
    dplyr::select(L2,L3) %>%
    dplyr::distinct(L3,.keep_all = T)
  
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
  # 
  # delFeatures <- dels %>%
  #   dplyr::left_join(geneAlias,by = 'mAlias') %>%
  #   dplyr::filter(ID==gene) %>%
  #   dplyr::rename(L1=ID) %>%
  #   dplyr::mutate(type="DELETION",strand=unique(gFeatures$strand)) %>%
  #   dplyr::left_join(parentL2,by="L1") %>%
  #   dplyr::rename(L2=ID) %>%
  #   dplyr::mutate(ID=paste(L2,type,sep="_")) %>%
  #   dplyr::select(seqid,start,end,type,strand,ID,L1,L2)
  
  #   remFeatures <- rbind(gFeatures,delFeatures)
  
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

###### lin-12 orthogroup

#ele
path <- "/gff/all_pc_genes/c_elegans.PRJNA13758.WS283.protein_coding.gff3"
gene <- "WBGene00003001"
features <- getGeneFeature(path,gene)

#bri - qx.2
path <- "/gff/all_pc_genes/Curation-VF-230612.PC.clean.renamed.WB.gff3"
gene <- "QX1410.13858"
features2 <- getGeneFeature(path,gene)

#bri - qx.1
path <- "/gff/all_pc_genes/Curation-VF-230612.PC.clean.renamed.WB.gff3"
gene <- "QX1410.13860"
features3 <- getGeneFeature(path,gene)

# #bri - af.1
# path <- "/projects/b1059/projects/Nicolas/c.briggsae/wormbase/WS280/c_briggsae.PRJNA10731.WS280.protein_coding.gff"
# gene <- "WBGene00029134"
# features4 <- getGeneFeature(path,gene)
# 
# #bri - af.2
# path <- "/projects/b1059/projects/Nicolas/c.briggsae/wormbase/WS280/c_briggsae.PRJNA10731.WS280.protein_coding.gff"
# gene <- "WBGene00029035"
# features5 <- getGeneFeature(path,gene)
# 
# #bri - af.3
# path <- "/projects/b1059/projects/Nicolas/c.briggsae/wormbase/WS280/c_briggsae.PRJNA10731.WS280.protein_coding.gff"
# gene <- "WBGene00029037"
# features6 <- getGeneFeature(path,gene)

#nig - ju.1
path <- "/gff/all_pc_genes/c_nigoni.PRJNA384657.WBPS19.protein_coding.gff3"
gene <- "nigoni.pc_2016.07.11_010.g27410"
features7 <- getGeneFeature(path,gene)

#nig - ju.2
path <- "/gff/all_pc_genes/c_nigoni.PRJNA384657.WBPS19.protein_coding.gff3"
gene <- "nigoni.pc_2016.07.11_010.g27411"
features8 <- getGeneFeature(path,gene)



sub_gff <- rbind(#features[[1]] %>% dplyr::mutate(tag="N2 lin-12") %>% dplyr::mutate(spp="N2"),
                 #features2[[1]] %>% dplyr::mutate(tag="QX1410 Cbr-lin-12.2") %>% dplyr::mutate(spp="QX1410"),
                 features3[[1]] %>% dplyr::mutate(tag="QX1410 Cbr-lin-12.1") %>% dplyr::mutate(spp="QX1410"),
                 #features4[[1]] %>% dplyr::mutate(tag="AF16 Cbr-lin-12.1") %>% dplyr::mutate(spp="AF16"),
                 #features5[[1]] %>% dplyr::mutate(tag="AF16 Cbr-lin-12.2") %>% dplyr::mutate(spp="AF16"),
                 #features6[[1]] %>% dplyr::mutate(tag="AF16 Cbr-lin-12.3") %>% dplyr::mutate(spp="AF16"),
                 features7[[1]] %>% dplyr::mutate(tag="JU1422 010.g27410") %>% dplyr::mutate(spp="JU1422"))
                 #features8[[1]] %>% dplyr::mutate(tag="JU1422 010.g27411") %>% dplyr::mutate(spp="JU1422"))

p1 <- plotNormGeneFeatures(sub_gff,F) +
  theme(axis.title.y  = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_line(),
        panel.background = element_blank(),
        legend.position="none",
        axis.title.x = element_blank(),
        axis.text = element_text(size = 8))


sub_gff <- rbind(#features[[1]] %>% dplyr::mutate(tag="N2 lin-12") %>% dplyr::mutate(spp="N2"),
  features2[[1]] %>% dplyr::mutate(tag="QX1410 Cbr-lin-12.2") %>% dplyr::mutate(spp="QX1410"),
  #features3[[1]] %>% dplyr::mutate(tag="QX1410 Cbr-lin-12.1") %>% dplyr::mutate(spp="QX1410"),
  #features4[[1]] %>% dplyr::mutate(tag="AF16 Cbr-lin-12.1") %>% dplyr::mutate(spp="AF16"),
  #features5[[1]] %>% dplyr::mutate(tag="AF16 Cbr-lin-12.2") %>% dplyr::mutate(spp="AF16"),
  #features6[[1]] %>% dplyr::mutate(tag="AF16 Cbr-lin-12.3") %>% dplyr::mutate(spp="AF16"),
  #features7[[1]] %>% dplyr::mutate(tag="JU1422 010.g27410") %>% dplyr::mutate(spp="JU1422"),
  features8[[1]] %>% dplyr::mutate(tag="JU1422 010.g27411") %>% dplyr::mutate(spp="JU1422"))

p11 <-  plotNormGeneFeatures(sub_gff,F) +
  theme(axis.title.y  = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_line(),
        panel.background = element_blank(),
        legend.position="none",
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 8)) 

#ggsave(plot = p1,"/projects/b1059/projects/Nicolas/collabs/forHelen/sp3_lin12_altgr.png",device = 'png',dpi = 600,width = 7,height = 4.5 ,units = 'in')

p111 <- plot_grid(p1,p11,nrow = 2,ncol = 1,rel_heights = c(1,1.15))
p111

align <- readr::read_delim("/alignments/JU1422vQX1410/JUtQX_transformed.tsv",col_names = T) %>%
  tidyr::separate(`[TAGS]`,into = c("CB","CN"),sep = "\t")
colnames(align) <- c("S1","E1","S2","E2","L1","L2","IDY","LENR","LENQ","CB","CN")

nigo_1 <- ape::read.gff("/gff/all_pc_genes/c_nigoni.PRJNA384657.WBPS19.protein_coding.gff3") %>%
  dplyr::filter(grepl("010.g27410",attributes)) %>%
  dplyr::filter(type=="CDS") %>%
  dplyr::select(seqid,type,start,end)

nigo_2 <- ape::read.gff("/gff/all_pc_genes/c_nigoni.PRJNA384657.WBPS19.protein_coding.gff3") %>%
  dplyr::filter(grepl("010.g27411",attributes)) %>%
  dplyr::filter(type=="CDS") %>%
  dplyr::select(seqid,type,start,end)

brigg_1 <- ape::read.gff("/gff/all_pc_genes/Curation-VF-230612.PC.clean.renamed.WB.gff3") %>%
  dplyr::filter(grepl("QX1410.13858",attributes)) %>%
  dplyr::filter(type=="CDS") %>%
  dplyr::select(seqid,type,start,end)

brigg_2 <- ape::read.gff("/gff/all_pc_genes/Curation-VF-230612.PC.clean.renamed.WB.gff3") %>%
  dplyr::filter(grepl("QX1410.13860",attributes)) %>%
  dplyr::filter(type=="CDS") %>%
  dplyr::select(seqid,type,start,end)

p2 <- ggplot() +
  #geom_rect(aes(xmin=0,xmax=Inf,ymin=10607219/1e6,ymax=10617555/1e6),alpha=0.5,fill="lightgrey") +
  #geom_rect(aes(xmin=0,xmax=Inf,ymin=10588085/1e6,ymax=10598994/1e6),alpha=0.5,fill="lightgrey") + 
  #geom_rect(aes(xmin=120007/1e6,xmax=130413/1e6,ymin=0,ymax=Inf),alpha=0.5,fill="lightgrey") + 
  #geom_rect(aes(xmin=143512/1e6,xmax=147489/1e6,ymin=0,ymax=Inf),alpha=0.5,fill="lightgrey") + 
  geom_rect(data=nigo_1,aes(xmin=start/1e6,xmax=end/1e6,ymin=0,ymax=Inf),fill='lightgreen',alpha=0.6)+
  #geom_rect(data=nigo_2,aes(xmin=start/1e6,xmax=end/1e6,ymin=0,ymax=Inf),fill='lightgreen')+
  #geom_rect(data=brigg_1,aes(ymin=start/1e6,ymax=end/1e6,xmin=0,xmax=Inf),fill='lightblue')+
  geom_rect(data=brigg_2,aes(ymin=start/1e6,ymax=end/1e6,xmin=0,xmax=Inf),fill='lightblue',alpha=0.6)+
  geom_segment(data=align %>% dplyr::filter(CB=="III" & CN=="PDUG01000010.1" & E2<S2) ,aes(x=S2/1e6,xend=E2/1e6,y=S1/1e6,yend=E1/1e6)) + 
  coord_cartesian(ylim=c(10.605,10.62),xlim=c(0.115,0.135)) +
  #ggtitle(expression(paste("Alignment of ",italic('C. briggsae '),"and ", italic('C. nigoni lin-12.1'), " orthologs"))) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill=NA),
        panel.grid = element_blank()) +
  xlab("JU1422 Physical position (Mb)") +
  ylab("QX1410 Physical position (Mb)") 

p3 <- ggplot() +
  #geom_rect(aes(xmin=0,xmax=Inf,ymin=10607219/1e6,ymax=10617555/1e6),alpha=0.5,fill="lightgrey") +
  #geom_rect(aes(xmin=0,xmax=Inf,ymin=10588085/1e6,ymax=10598994/1e6),alpha=0.5,fill="lightgrey") + 
  #geom_rect(aes(xmin=120007/1e6,xmax=130413/1e6,ymin=0,ymax=Inf),alpha=0.5,fill="lightgrey") + 
  #geom_rect(aes(xmin=143512/1e6,xmax=147489/1e6,ymin=0,ymax=Inf),alpha=0.5,fill="lightgrey") + 
  #geom_rect(data=nigo_1,aes(xmin=start/1e6,xmax=end/1e6,ymin=0,ymax=Inf),fill='lightgreen',alpha=0.6)+
  geom_rect(data=nigo_2,aes(xmin=start/1e6,xmax=end/1e6,ymin=0,ymax=Inf),fill='lightgreen',alpha=0.6)+
  geom_rect(data=brigg_1,aes(ymin=start/1e6,ymax=end/1e6,xmin=0,xmax=Inf),fill='lightblue',alpha=0.6)+
  #geom_rect(data=brigg_2,aes(ymin=start/1e6,ymax=end/1e6,xmin=0,xmax=Inf),fill='lightblue',alpha=0.6)+
  geom_segment(data=align %>% dplyr::filter(CB=="III" & CN=="PDUG01000010.1" & E2<S2) ,aes(x=S2/1e6,xend=E2/1e6,y=S1/1e6,yend=E1/1e6)) + 
  coord_cartesian(ylim=c(10.585,10.60),xlim=c(0.135,0.15)) +
  #ggtitle(expression(paste("Alignment of ",italic('C. briggsae '),"and ", italic('C. nigoni lin-12.1'), " orthologs"))) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill=NA),
        panel.grid = element_blank()) +
  xlab("JU1422 Physical position (Mb)") +
  ylab("QX1410 Physical position (Mb)") 


tree <- ggtree::read.tree("/tree/OG0001886.mafft.treefile")
groups <- data.frame(seq=c("Cbr-lin-12.1","Cbr-lin-12.2","lin-12","010.g27410","010.g27411"),group=c(1,1,2,3,3))
#tree <- ggtree::groupOTU(tree,branches)

pt<-ggtree(tree,layout="daylight") + coord_cartesian(clip="off",ylim = c(-0.25,0.38),xlim = c(-0.55,0.1)) 
pt2 <- pt %<+% groups +  
  geom_tiplab2(aes(color=factor(group)),offset=-0.003,fontface="italic",size=2.7) +
  theme(legend.position = 'none') +
  scale_color_manual(values=c("lightblue","orange","lightgreen"))
pt2

p231 <- cowplot::plot_grid(p2,p3,p111,pt2,labels = c("A","B","C","D"),rel_heights = c(1.4,1))

ggsave(plot=p231,"sp2_CBCN_lin12.png",device = 'png',dpi = 600,width = 7,height = 6.5 ,units = 'in')
