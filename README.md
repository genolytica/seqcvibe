# SeqCVIBE
A shiny application for visualization of Big Next Generation Sequencing datasets

----

# Overview

Analyzing and visualizing an RNA-Seq analysis entirely is quite a demanding task. There is a wide range of tools required in the process of depicting all the significant information produced by an experiment. In almost all cases they have to be combined in a well-orchestrated way in order for a complete picture to be drawn and a full image of the analysis to be delivered, as there is no such thing as a 'swiss army knife' for this process. SeqCVIBE does exactly that!

## SeqCVIBE Overview

SeqCVIBE is an R-shiny tool that aggregates all the analyses involved in an RNA-Seq experiment in a well-organized and clear maner. The user can select an analysis from a variety of experimental designs derived from two of the most enriched dataset repositories, Gene Expression Omnibus (GEO) and European Nucleotide Archive (ENA). These experiments cover areas like profiling of genes involved in cancer or other neurological diseases, organ evolution, stem cell expression signature profiling etc. 

These analyses are fully configurable allowing the user to control every parameter in very simple steps, especially compared to the often confusing way a command line tool works. In addition to that, all analyses are reproducible and can be easily shared across a group of users who can pick-up from any given point and move on anyway they prefer.

# Data Selector

The 'Data Selector' tab is the root of the analysis. The user is prompted to select one of the available experiments to analyze. The datasets are sourced from two of the largest databases, GEO and ENA. After selection, the dataset is loaded and some basic information about the dataset are being displayed. At this point, SeqCVIBE is ready to use.

## Input data selection

The first step of any analysis is the dataset selection. This panel allows the user to select an input from one of the currently available sources. Each source has its own bunch of datasets and each dataset has a reference species genome which is always automatically selected based on the dataset. The option to select a dataset based on the reference genome is not provided at the moment.

Once the selection is made, the relevant short summary of the experiment will appear right below. Clicking on the title will redirect the user to the relevant GEO or ENA page were the analysis design and details are more thoroughly described. 

Loading the dataset allows for exclusion of specific samples per condition or exclusion of all samples belonging to a single condition of the experiment. This way, the user has the ability to create a fully customizeable dataset.

After a dataset is loaded, the option to inspect its quality is enabled with MultiQC. By clicking on 'Dataset MultiQC' a pop-up presenting an aggregated analysis of all quality metrics for all samples of the dataset will appear. This inspection can prove essential in the selection or exclusion of samples to be included in the final dataset based on specific properties the user might find significant.

When the dataset selections are made the user can hit the 'Create dataset' button to finalize the dataset. This will "lock" the sample options and all specific inclusions applied and will set the base for all follow-up analyses. The included samples of the final dataset will be displayed in the tab on the right named 'Current dataset'. After this point, anytime the user wants to run a fresh analysis with a new dataset they have to 'Clear dataset' and start over again.

## Bookmark application

Bookmark's purpose is quite straight-forward. It places a bookmark in every tab and every visualization of the current status of SeqCVIBE, and allows for easy restore to that checkpoint at anytime. In other words it generates a save screenshot which can be loaded by any user of the group. At the end of the analysis the user has to enter a name for the bookmark and click the 'Add Bookmark' button. Then the bookmark will be added to the list of available bookmarks and any user will have the ability to restore it to its current state by clicking on the relevant URL. After the bookmarked state is fully loaded, users can make additional changes and save the analysis on top of the previous bookmark.

# Signal viewer

The Signal viewer allows for visual representation of one or more genomic regions of interest with referece to normalized signal. There are two options available. Gene Signal and Area Signal.

## Gene signal

Gene signal plots the combined normalized signal of the selected gene(s) features/conditions. The user can select plotting input between specific genes by their name, or by a custom input genomic region. Further down they can optionally change the upstream and downstream regions to be plotted, select the different colours for each feature/condition of the experiment and finally, select one method for the gene profile averaging among Mean, Median and Trimmed mean. It the last case, the trim's fraction can be manually defined.

By hitting the 'Engage' button the signal of all the selected genes and custom areas will be calculated and averaged. Once calculations are complete the plot will appear in the center of the tab. If there's more than one selection of gene or custom region (or both) to be displayed, an equivalent number of plots will be displayed in the same tab.

### Gene plotting

Giving a gene input is very simple. By clicking on the space below 'Plot data' a drop-down menu of all known gene names of the concerning experiment species will appear. The search box supports smart-input method which allows for fastest searching of a specific gene. Alternatively, the user can always input the gene name they are considering, however this is not adviced due to the fact that gene name is case sensitive. If the gene is selected successfully it will be displayed as a new entry in the below table. There, its genomic coordinates will be also displayed.

### Custom region plotting

Other than plotting a specific gene region there is an option to plot a custom genomic area of interest, by selecting the 'Custom' tab in the 'Plot data' panel. The user will be prompted to fill in the requested area coordinates, like chromosome number, start-end position of the area and strand. The given area should be alo given a unique name to be separated from any other regions that may be entered. The 'Add' and 'Remove' buttons below allow for the addition of a new region to the table of regions to be plotted or the removal of the selected ones.

## Area signal

The Area signal works in a similar manner to the Gene signal utility. It's difference lies in the display style of the results which is of a wider area than what's just covered by the gene and resembles that of a genome browser. It can also plot multiple genes in a region without separating the plots for every region, like it happens with the Gene signal option.

### Custom genomic area

Custom genomic area will visualize a whole chromosomal region of interest input by the user. By default Known genes overlaping the given area will be automatically detected. However, if the user wants any specific regions of interest to be also included to the automatic detection, along with the known genes, they just have to click on the option "Include custom transcribed areas" and provide the desired regions either by chosing to include the gene explorer entries, or by manually inserting the custom regions.

### Around gene area

This options resembles the 'Gene signal' to a great extend as it plots the requested genomic area that is centered around a specific gene of interest. The gene can be selected from the drop-down menu right below, while the upstream and downstream limits of the plot can be given from the 'Flanking' options bars

## Plot options

Plot options are to a large extend identical between Gene signal and Area signal plotting. First option is the flanking region limits. Presets are 2000 base-pairs upstream and downstream the region of interest. Other plotting options include selecting a color for every condition of the experiment and finally, selecting a method to be used in the profile  averaging of the area to-be-plotted. These methods include Mean signal, Signal median and the trimmed mean. In the latter case the user can select the fraction to be trimmed from the dataset manually. Here, it's worth noting that these options are relevant to their own tab selection and are therefore separate between area and gene signal plotting.

## Exporting plot

There are currently three ways a user can export a plot. The 'Export ggplot2' method exports a .Rda file. This file can be loaded in R (*load(file='file.rda')*) and allow a user with relevant experience to further modify the exported visualization. Alternative exporting options are in the png and pdf file formats. All exports are stored by default in the system's dedicated 'Downloads' directory. These exporting options are common across all other plot panels.

# Expression viewer

The expression viewer tab is a presentation of the features count table, which is the backbone of a differential expression analysis in an RNA-Seq experiment. The table demonstrates the expression measurements for the input selection of genes or genomic regions. Two tables will be displayed after the selections are made. The top one shows expression measurements for the first ordered condition of the experiment (usually the wild-type) and the bottom one shows the equivalent measurements for the 2nd condition. 

The first column of the table is the name of the gene or the selected region, the second is the mean (by default) averaged expression and the third is the expression deviation between samples. Default deviation measure is the standard deviation. If the experiment contains more than two conditions, an equivalent number of additional tables will be displayed.

## Known genes expression

First option of the Expression viewer is the 'Known genes'. This will allow the user to obtain a table measuring expression of a single gene, a custom set of genes, or even all genes of the studied organism. The input genes can be selected either from the given list of known gene names of the organism, or as a custom gene name input, where every entry is separated by a new line. Finally, there is also the choice of viewing expression for all genes, if the third option is clicked.

## Expression calculator

The expression calculator works the exact same way as the 'Known genes', the only difference being that the user needs to give a pre-defined genomic region as input. There are two ways this can be done. Either with the gene explorer or the custom regions. Genomic regions given to the gene explorer will be used to produce the expression tables in the first scenario. Alternatively, if the user selects 'Custom regions' an input area will appear right below where they can input the genomic coordinates they wish to inspect and click the 'Add' button. All added entries will appear in the mini-table below. Hitting the 'Engage' button will generate the requested tables.

## Expression settings

There are a number of configurable parameters the user can play with to create a desired output expression table. First one has to do with how the RNA-Seq expression will be measured. There are three available options here: 
Raw counts is the raw number of gene/genomic area observations in each conditions. RPKM and RPGM refer to the Reads Per Kilobase/Gigabase of transcript, per Million mapped reads. It is a normalized unit of transcript expression which scales by transcript length to compensate for the fact that most RNA-seq protocols will generate more sequencing reads from longer RNA molecules.
Second option is a scaling setting for the expression measure. The user gets to select between Natural and log2 scale. Log2 scale is usually preferable as it models proportional changes rather than additive changes. This is typically biologically more relevant.
Next is the selection of the expression measure averaging method. Available options are 'Mean' and 'Median'
Lastly, there is the expression deviation measure selection. Here the user can select among three methods by which the deviation of expression between the samples will be measured: Standard deviation, Mean absolute deviation and Interquartile range.

## Exporting table

The resulting table can be exported in different ways. Cells of the resulting tables can be selected by clicking on them. Buttons below each table control what you can do with the selection. 'Clear selection' removes all selected rows from the final selection. 'Invert selection' will include all non-selected rows to the final selection. When selection set has been finalized the user can click on 'Export selection' to obtain a csv file with all the configured information. At any case, the option 'Export all' allows the user to export the whole dataset of gene expression. Exported tables and their settings are unique per condition.

# Analysis

Analysis tab is the most essential as it is a toolkit for complete visualization of an RNA-Seq data experiment. The four visualizations provided include Differential expression, Clustering analysis, Correlation analysis and MDS/PCA analysis.

## Differential Expression analysis

Differential expression analysis means taking the normalised read count data and performing statistical analysis to discover quantitative changes in expression levels between experimental groups. For example, in a typical experiment statistical testing would bee used to decide whether, for a given gene, an observed difference in read counts is significant, that is, whether it is greater than what would be expected just due to natural random variation.

### Pipeline settings

General settings of the analysis include the following: First, the analysis pipeline that will be used. This only changes the tool that runs the analysis, which means there is a slight differentiation in the algorithm. However the results will not vary significantly. The three widely used tools that are included in the available options are DESeq, edgeR and voom. Further down, the user has to select which condition of the experiment will be used as control for the analysis, the multiple testing correction method and when will the selected filters (next tab) be applied to the results (before or after the normalization process). Last setting is the inclusion of the custom regions from the expression calculator to the analysis results. If the tick box is clicked the results will include expression measurements for the custom regions. Filtering settings, on the tab on the right will be applied to the analysis after it is complete, or beforehand, depending on the user's choice in the General settings tab.
 Basic gene filter options include Mean/Median expression, qualtile and known genes. If quantile gene filter is selected user will be prompted to fill in the desired percentile which will be used to filter genes. In the case of 'known genes' a drop-down list will appear with all the species specific genes. Users selections will restrict the results to only those genes selected. Additional filtering options include filtering a gene by its size in base-pairs and by specific biotype.

### Results table

When the analysis is complete an MA plot will appear in the central panel of SeqCVIBE along with the results table right below. The table contains four tabs, first of which is the 'Summary'. This contains all the basic information about the outcome of the analysis. It consists of the following columns: gene_name, the p-value of the findings for that gene, the false discovery rate, the expression ratio between different conditions of the experiment and lastly the mean and standard deviation counts for each condition. Second tab is the annotation tab. It provides a table with the genomic coordinates of each differentially expressed gene and its biotype description. Flags tab #################################################### . The last tab provides the option of combining all the aforementioned tables in one. Results table can be further configured with the 'Results table settings' panel. The user can control which significant gene hits will show up by choosing a threshold for one of the following measures: p-value, FDR or fold change. Post analysis filtering is vailable on the 'Filters' tab which is next. A mix of gene_name, chromosome and biotype filters can be applied to further narrow down the results to something more meaningful to the user. Scale tab will allow for configuration of the measures used in the results table. By default, summary value components will be the gene counts, however RPKM and RPGM are also available. Similarly, the user might switch between natural and log2 count for the differential expression ratio, mean or median for averaging and sd, median absolute deviation or interquartile range for intra-sample deviation. The table can be fully or selectively exported, just like with the Expression viewer table.

### Plot options

The resulting MA plot will appear after the analysis is complete. If real time table update is ticked, selecting an area of the plot will immediately reduce the results table only to the genes of interest contained in that area. By default gray spots indicate a gene that is not differentially expressed, red ones indicate an overregulated gene while greens indicate a downregulated one. Colouring options are available for all three categories in the 'Plot colors' panel.

## Clustering analysis

Clustering analysis is a classification algorithm that aims to group the set of samples included in the final dataset in order to create a set of sample obects where each are similar within the same cluster and dissimilar to samples in other clusters. It helps visualizing the intra-sample diversity of the datasets and helps determine if each of the experiment conditions are separated in well-differentiated categories. Clusteringsettings allow for the customization of the resulting heatmap.

### Gene settings

Gene settings tab allows the user to control which genes will be used to produce the heatmap. They can either be custom genes/regions or a specified gene list or just the genes that were found to be differentially expressed (if a DE analysis has already ran. Next filter is the Clustering variable selection. The user has to choose one of the available variables that will be used to define the clusters of the analysis. An expression measurment filter is also available here, allowing for normalized counts vs RPKM vs RPGM measure or log2 expression scaling vs Natural coount.

### Heatmap settings

Heatmap settings determine the mathematical parameters of the visualization. The first choice defines the metric that will be used to define the clusters distance. Default option is the Euclidean distance. Right below the user gets to select the linkage function that will be used to provide the relevant dendrogram of the resulting clusters. The default here is Average linkage. The clustering dimensions can also be modified right below along with the plot colors. 

## Correlation analysis

Correlation analysis will give an insight on how strong the relation is between samples included in the dataset in terms of their specific gene expression measurements or between gene expression levels among samples. 

### Settings

Before starting the analysis the gene settings will have to be configured first. First option has to do with what set of expressed genes will be used to calculate the correlation metric. The user can select all genes with non-zero counts, which is usually not feasible due to the large number of genes that will delay the calculations significantly, all expressed genes, which requires a differential expression analysis to have already run, a custom gene list, or a custom selection of genes from the available drop-down box. Expression measurement settings will affect how the counts will be measured in the 'Data matrix' table below the correlation plot (see additional matrices). Correlation settings will determine the type of the expression correlation. 'Sample-wise' will correlate all samples based on their overall gene expression levels, 'Gene-wise' will correlate all genes given in the  input list provided and 'Reference gene-wise' will calculate the correlation based on a given reference gene of interest. Lastly, the correlation method can be selected next and can either be Pearson's or Spearman's. Finally the panel to the left is the color selection panel. There correlation intensity colors can be selected manually and the resulting plots can be exported in PNG or PDF format.

### Additional matrices

Right under the main correlation plot three additional supporting visualization will be displayed after a successfully complete correlation analysis. From left to right there's the MDS plot giving a quick look on how well differentiated the samples are based on the genes selected for correlation. The middle table is the correlation matrix which shows the relevant correlation measures for the selected samples. The last table is more of a supporting data matrix which demonstrates the gene expression measures for the selected genes. These threee visualizations combined can give a well-concentrated amount of information for the selected gene set.

## MDS/PCA Analysis

The goal of MDS Analysis analysis is to detect meaningful underlying dimensions that allow the user to explain observed similarities or dissimilarities (distances) between the investigated samples. Principal component analysis (PCA) is one of the commonly used statistical technique for finding patterns in data of high dimensions and expressing the data in such a way as to highlight their similarities and differences. 

### Settings

Gene settings here, as in previous analysis tabs, are related to gene set selection and what expression measure and scale will be used for determining the intra-sample distances and patterns. 'MDS/PCA' tab sets the variance projection method that will be used between the two. Based on this selection, the relevant settings will appear right below. MDS settings allow for the selection of the (dis)similarity metric and the number of dimensions to be included. PCA settings allow selecting whether to scale the values or center them. The resulting plots can be exported as R objects as well, appart from PNG or PDF.

### Coordinate plotting

Coordinate plotting sets which two principal components or coordinates will be used in the plot and in which axis. Additional settings further down enable point names in the resulting plot or highlighting of any set of genes that might be of interest.

# Genome browser

pending...