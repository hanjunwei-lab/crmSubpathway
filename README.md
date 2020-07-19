crmSubpathway: Identify cancer-related metabolic subpathways

Description: This function uses the k clique algorithm to split the metabolic pathways in the KEGG database into metabolic subpathways. Subsequently, a stable metabolic subpathway activity matrix is constructed by GSVA or ssGSEA methods. Eventually, cancer-related metabolic subpathways are identified through differential analysis.

Installation: library(devtools); install_github("hanjunwei-lab/crmSubpathway")

use：library(crmSubpathway)
