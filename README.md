# Venn_Count_Matrix
The venn diagram is generated from the counts matrix.
The bars represent the number of genes that are present in samples.
The dots represent which samples those numbers are present.

### Run from command line
`python venn_count_matrix.py -h`  
`usage: venn_count_matrix.py [-h] -i <inputfile> -m <metadatafile> -s <sample> -p <pattern> -o <outfiel>`  

`Using the count matrix generates venn diagram with genes that have mean normalixzed counts > 10 for each group.`  

`optional arguments:                     
  -h, --help         show this help message and exit                        
  -i <inputfile>     input count file                
  -m <metadatafile>  metafile with description of samples                 
  -s <sample>        name given to sample file               
  -p <pattern>       name given to pattern file              
  -o <outfiel>       figure name`              

![Figure](https://github.com/sue02/Venn_Count_Matrix/blob/master/figure.png)
