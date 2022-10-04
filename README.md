# TopGo
ESte programa calcula el enriquecimiento de los GOs para un conjunto de genes de interés y el universo de genes con anotación

#Parametros de entrada:
MyDirWork = args[1]:  Directorio de Trabajo.  ejemplo:  Dropbox/Projects/Project_FA/  

MyGeneId_GOfile = args[2]. : Archivo de texto con 2 columnas:  la primer columna con los ids de todo el universo de los genes o transcritos, 
      la segunda columna los GO_id separados con , correspondientes al gene_id.  Ejemplo:
      
       more  ~/Dropbox/Projects/Project_FA/GeneId_GOs.txt
       
      TRINITY_DN28190_c0_g3   GO:0005829,GO:0002144,GO:0016779,GO:0000049,GO:0032447,GO:0034227,GO:0002143,GO:0002098
      
      TRINITY_DN27280_c0_g1   GO:0016021
      
      TRINITY_DN28490_c0_g2   GO:0030054,GO:0031410,GO:0030425,GO:0005768,GO:0035838,GO:0044231,GO:0005887,GO:0031224,GO:0045121,GO:0005739
      
      TRINITY_DN27682_c0_g2   GO:0097361,GO:0051539,GO:0046872,GO:0016226,GO:0002040
      
      TRINITY_DN27332_c1_g1   GO:0030122,GO:0005905,GO:0030136
      
      ...
      
InterestGenesFile = args[3] :  una lista de texto, con lo gene_id o transcritos_id de interes (porque fueron UP o Down en un análisis de expresion
      diferencial,  por ejemplo:
      
      more UP_cond1_DOWN_cond2_Ids.txt
      
      TRINITY_DN24527_c0_g1
      TRINITY_DN27947_c2_g1
      TRINITY_DN30921_c1_g5
      TRINITY_DN27452_c0_g1
      TRINITY_DN23466_c0_g1
      TRINITY_DN27278_c0_g1
      TRINITY_DN23340_c0_g2
      TRINITY_DN5046_c0_g1
      TRINITY_DN31282_c0_g1
      TRINITY_DN28626_c0_g2
      TRINITY_DN25611_c1_g1
      TRINITY_DN25849_c0_g5
      TRINITY_DN24582_c0_g2
      TRINITY_DN27406_c2_g2
      ...
# Ejemplo de ejecución:
Desde la consola de RStudio:

RScript MyTopGo.R "~/Dropbox/Projects/Project_FA/" "GeneId_GOs.txt" "UP_Cond1_DOWN_cond2_Ids.txt"
