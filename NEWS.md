# mapa 2.1.0

* 1. Allow data visualization for expression data and functional module bipartite graph.
* 2. Allow using SiliconFlow for text embedding and LLM interpretation.

# mapa 2.0.0

* 1. Enable to use both ChatGPT and Gemini models to perform llm interpretation.
* 2. Enable to help user determine the optimal clustering parameters and assess clustering quality.
* 3. More clustering algorithms are provided, including: binary cut, hierarchical clustering, and eight graph-based clustering methods.

# mapa 1.6.2

* Combined `merge_modules()` and `merge_pathways_bioembedsim()` into one generia function `get_functional_modules()`.

# mapa 1.6.0

* Enable clustering evaluation.

# mapa 1.5.1

* Fixed bugs in `merge_pathways()`.

# mapa 1.5.0

* Now for non-model organism, users can provide AnnotationHub ID to perform the following analysis. 

# mapa 1.4.1

* For metabolite analysis, if a pathway do not have any metabolite, no metabolite information will be available for PubMed Search.

# mapa 1.4.0

* ID conversion is now available for mapa package.

# mapa 1.3.0

* Two more clustering methods "binary cut" and "hierarchical" are now available for clustering overlapped based similarity calculation results.

# mapa 1.2.0

* Fixed LLM interpretation bugs during RAG.

# mapa 0.1.1

* Fix bugs.

# mapa 0.1.2

* Fix bugs.

# mapa 0.1.3

* Fix bugs.

# mapa 0.1.4

* Fix bugs.

# mapa 0.1.5

* Fix bugs.

# mapa 0.1.6

* Fix bugs.

# mapa 0.1.7

* Improve the plot_pathway_bar function.

# mapa 0.1.8

* Support ChatGPT and Gemini.

# mapa 0.1.9

* Use the ChatGPT to interpretate the pathways.

# mapa 0.1.10

* Add the launch_mapa_shiny function.

# mapa 0.1.11

* Improve launch_mapa_shiny function and add filter_functional_module function.

# mapa 0.1.12

* Translate language functions can be used for all the enriched functional results.

# mapa 0.1.13

* Update shiny app.

# mapa 0.1.14

* Update shiny app.

# mapa 0.1.15

* Fix bugs.

# mapa 0.1.16

* Fix bugs.

# mapa 0.1.17

* Fix bugs.

# mapa 0.1.18

* Fix bugs for shiny app.

# mapa 0.1.19

* Added GSEA methods.
* Fixed one bug in the show function.

# mapa 0.1.20

* Fixed one typo.

# mapa 0.1.21

* Modify term similarity calculation and plot_pathway_bar.

# mapa 0.1.22

* Fixed bug in term semantic similarity calculation.

# mapa 0.1.23

* Fixed bug in plot_pathway_bar for "functional_module". Added title to the plot.

# mapa 0.1.24

* 1.Fixed bug in `do_gsea()` and `plot_pathway_bar()` for GSEA. 
* 2.Added database color to the barplot for `plot_module_info()` when level equals to "functional_module". 
* 3.Added visualization for GSEA result.

# mapa 0.1.25

* Modified `plot_relationship_network()` for the visualization for GSEA result.

# mapa 0.1.26

* Added metabolite enrichment analysis based on HMDB and KEGG.

# mapa 0.1.27

* New `get_bioembedsim()` supports pathway similarity calculation via embedding.

# mapa 0.1.28

* Added clustering based on biotext embedding similarity calculation.
* Added cellular component sub-ontology.

# mapa 0.1.29

* New `llm_interpret_module()` utilizes LLMs integrated with RAG to annotate functional modules.

# mapa 0.1.30

* Allow analysis of the association between functional modules and specified phenotype using LLM.

