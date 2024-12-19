# Multicenter Next-Generation Sequencing Studies between Theory and Practice

<p align="center">
    <img height="500" src="https://uni-muenster.sciebo.de/s/Yqv8hcXYvcHryZD/download">
</p>

Sandmann S, de Graaf AO, Tobiasson M, Kosmider O, Abáigar M, Clappier E, Gallì A, van der Reijden BA, Malcovati L, Fenaux P, Díez-Campelo M, Fontenay M, Hellström-Lindberg E, Jansen JH, Dugas M. Multicenter Next-Generation Sequencing Studies between Theory and Practice: Harmonization of Data Analysis Using Real-World Myelodysplastic Syndrome Data. J Mol Diagn. 2021 Mar;23(3):347-357. doi: 10.1016/j.jmoldx.2020.12.001.

## Background

In the age of personalized medicine, genetic testing by means of targeted sequencing has taken a key role. However, when comparing different sets of targeted sequencing data, these are often characterized by a considerable lack of harmonization. Laboratories follow their own best practices, analyzing their own target regions. The question on how to best integrate data from different sites remains unanswered. Studying the example of myelodysplastic syndrome (MDS), we analyzed 11 targeted sequencing sets, collected from six different centers (n = 831). An intersecting target region of 43,076 bp (30 genes) was identified; whereas, the original target regions covered up to 499,097 bp (117 genes). Considering a region of interest in the context of MDS, a target region of 55,969 bp (31 genes) was identified. For each gene, coverage and sequencing data quality was evaluated, calculating a sequencing score. Analyses revealed huge differences between different data sets as well as between different genes. Analysis of the relation between sequencing score and mutation frequency in MDS revealed that most genes with high frequency in MDS could be sequenced without expecting low coverage or quality. Still, no gene appeared consistently unproblematic for all data sets. To allow for comparable results in a multicenter setting analyzing MDS, we propose to use a predefined target region of interest and to perform centralized data analysis using harmonized criteria. 

To determine our sequencing score, we defined several thresholds, identifying low coverage, mapping quality (MQ), base quality (BQ), and background noise (BN). Any change in the scoring as well as the considered data sets may lead to a change in sequencing score. To test the influence of different scoring (absolute versus relative, different thresholds, and different weighting of quality parameters), selection of data sets, and selection of target region (intersecting versus of interest), an R Shiny interface has been developed. Every parameter influencing the sequencing score and mutation frequency can easily be adjusted and an updated plot can be generated.


## Requirements
To run the Shiny app, you need R Version 4.1.0 or higher.

Download the RData-files ofInterest.RData (https://uni-muenster.sciebo.de/s/xNl1hIXmhbaMQip) and intersecting.RData (https://uni-muenster.sciebo.de/s/OQKXPSF0Gi02VpT). Make sure they lie in the same folder as server.R and ui.R when you execute the Shiny app.

## Detailed information
For detailed information on the harmonization of NGS studies, please consider our publication (https://doi.org/10.1016/j.jmoldx.2020.12.001).

In case of errors or feature requests, do not hesitate to open an issue or contact Sarah Sandmann (sarah.sandmann@uni-muenster.de).

## Citation
Sandmann S, de Graaf AO, Tobiasson M, Kosmider O, Abáigar M, Clappier E, Gallì A, van der Reijden BA, Malcovati L, Fenaux P, Díez-Campelo M, Fontenay M, Hellström-Lindberg E, Jansen JH, Dugas M. Multicenter Next-Generation Sequencing Studies between Theory and Practice: Harmonization of Data Analysis Using Real-World Myelodysplastic Syndrome Data. J Mol Diagn. 2021 Mar;23(3):347-357. doi: 10.1016/j.jmoldx.2020.12.001.
