## crmSubpathway: Identify cancer-related metabolic subpathways

> Description: The crmSubpathway package is a systematic biological tool to identify cancer-related metabolic subpathways. The main capabilities of this tool are as follows: 
1. This function uses the k-clique algorithm to split the metabolic pathways in the KEGG database into metabolic subpathways. 
2. A stable metabolic subpathway activity matrix is constructed by GSVA or ssGSEA methods. 
3. Cancer-related metabolic subpathways are identified through differential analysis. 
4. Visualization


## Installation: 
```R
library(devtools)
install_github("hanjunwei-lab/crmSubpathway")
use：library(crmSubpathway)
```

The `crmSubpathway` is published in Frontiers in Oncology. Please cite the following article when using `crmSubpathway`:  
Xudong Han, Donghua Wang, Qingfei Kong and Junwei Han. Inference of Subpathway Activity Profiles Reveals Metabolism Abnormal Subpathway Regions in Glioblastoma Multiforme. Frontiers in Oncology. https://doi.org/10.3389/fonc.2020.01549. 
