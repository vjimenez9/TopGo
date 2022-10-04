# Description:
#  este programa realiza el analisis de enriquecimiento de GOs
#  Recibe de entrada, el path de trabajo y 2 archivos:
#   MyDirWork:  Directorio de trabajo
#   GeneId_GOFile  : universo de anotaciones. es un archivo con 2 columnas, en la primera los geneIds y la segunda todos lo GOsIds correspondientes
#   InterestGeneFile :  lista de genes de interes.  Uno por renglon. 
# ejemplo de Ejecución:
# RScript /Dropbox/bin/R/MyTopGo.R "~/Dropbox/Projects/Project_FrecAlim/OnlyActinbopterygii" "Actinopterygii_GT75_GeneId_GOs.txt" "UPFA12_DOWNFA8_Ids.txt"
#if (!requireNamespace("BiocManager", quietly=TRUE))
#  install.packages("BiocManager")
#BiocManager::install()

#BiocManager::install("topGO")

library(topGO)
args <- commandArgs(trailingOnly = TRUE)
MyDirWork = args[1]
MyGeneId_GOfile = args[2]
InterestGenesFile = args[3]



#MyDirWork="~/Dropbox/Projects/Project_FrecAlim/OnlyActinbopterygii"
#MyGeneId_GOfile="Actinopterygii_GT75_GeneId_GOs.txt"
#InterestGenesFile="UPFA12_DOWNFA8_Ids.txt"

buildingObjet<-function(Outputfile,geneList, geneID2GO,OntologyType,myInterestingGenes){
  # Construir el objeto GoData:
  #OntologyType="CC"
  #Outputfile=outputfile1
  print (paste("Ontology Type:",OntologyType))
  # Outputfile="UPFA12_DOWNFA8_Ids.MolecularFunction.out"
  myGOdata <- new("topGOdata", description="My project", ontology=OntologyType, allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO)
  # Lista de genes de interes que pueden ser accesados:
  sg <- sigGenes(myGOdata)
  str(sg)
  numSigGenes(myGOdata) 
  print(paste("Numero de Genes significativos:",numSigGenes(myGOdata)) )
  # Realizamos el enriquecimiento funcional:
  resultClassic <- runTest(myGOdata, algorithm="classic", statistic="fisher")
  resultElim <- runTest(myGOdata, algorithm="elim", statistic="fisher")
  resultTopgo <- runTest(myGOdata, algorithm="weight01", statistic="fisher")
  resultParentchild <- runTest(myGOdata, algorithm="parentchild", statistic="fisher") 
  
  # see how many results we get where weight01 gives a P-value <= 0.001:
  mysummary <- summary(attributes(resultTopgo)$score <= 0.009)
  numsignif <- as.integer(mysummary[[3]]) # how many terms is it true that P <= 0.001
  
  # Solicitando la lista de los 10 mas revelantes:
  allRes <- GenTable(myGOdata, classicFisher = resultClassic, elimFisher = resultElim, topgoFisher = resultTopgo, parentchildFisher = resultParentchild, orderBy = "topgoFisher", ranksOf = "classicFisher", topNodes = numsignif)
  write.table(allRes,Outputfile)
  
  # tratemos de visualizar:
  #showSigOfNodes(myGOdata, score(resultFisher), firstSigNodes = 5, useInfo ='all')
  #printGraph(myGOdata, resultFisher, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)
  
  # print out the genes that are annotated with the significantly enriched GO terms:
  myterms <- allRes$GO.ID
  mygenes <- genesInTerm(myGOdata, myterms)
  df=data.frame()
  for (i in 1:length(myterms))
  {
    myterm <- myterms[i]
    mygenesforterm <- mygenes[myterm][[1]] 
    myfactor <- mygenesforterm %in% myInterestingGenes # find the genes that are in the list of genes of interest
    mygenesforterm2 <- mygenesforterm[myfactor == TRUE] 
    mygenesforterm2 <- paste(mygenesforterm2, collapse=',')
    line=c(myterm,mygenesforterm2)
    df<-rbind(df,line)
  }
  names(df)<-c("Term","Genes")
  write.table(df,file=paste("Genes_",Outputfile,sep=""),quote = FALSE)
}

Enrichment<-function(MyDirWork,GeneId_GOFile,InterestGeneFile){
  setwd(MyDirWork)
  # Leemos el universo de los genes:
  geneID2GO<-readMappings(MyGeneId_GOfile)
  geneNames <- names(geneID2GO)
  # tomamos genes de interes:
  myInterestingGenes<-read.table(InterestGenesFile,header=FALSE)
  outputfile=tools::file_path_sans_ext(InterestGenesFile)
  myInterestingGenes<-as.character(myInterestingGenes$V1)
  geneList <- factor(as.integer(geneNames %in% myInterestingGenes))
  names(geneList) <- geneNames
  outputfile1=paste(outputfile,".MolecularFunction.out",sep="")
  buildingObjet(outputfile1,geneList, geneID2GO,"MF",myInterestingGenes)
  outputfile1=paste(outputfile,".BiologicalProcess.out",sep="")
  buildingObjet(outputfile1,geneList, geneID2GO,"BP",myInterestingGenes)
  outputfile1=paste(outputfile,".CelullarComponents.out",sep="")
  buildingObjet(outputfile1,geneList, geneID2GO,"CC",myInterestingGenes)
}

Enrichment(MyDirWork, MyGeneId_GOfile, InterestGeneFile)