# PDXplorer App

## To run locally

### Install Dependencies

```bash
R -q -s -e 'if(!requireNamespace("pak",quietly=TRUE)){install.packages("pak")};'

R -q -s -e 'pak::pkg_install(pkg=c(
  "shiny",
  "shinydashboard",
  "shinyjs",
  "shinythemes",
  "shinymanager",
  "bslib",
  "rhino",
  "cachem",
  "data.table",
  "reshape2",
  "fs",
  "dplyr",
  "tidyr",
  "readr",
  "ggdendro",
  "ggplot2",
  "grid",
  "gridExtra",
  "DT",
  "scales",
  "pheatmap",
  "RColorBrewer",
  "bioc::maftools",
  "bioc::limma",
  "bioc::DESeq2",
  "bioc::scran",
  "bioc::EnhancedVolcano",
  "bioc::chimeraviz",
  "github::statistikat/codeModules",
  "github::marchionniLab/ABCutilities"
  ),upgrade=TRUE,ask=FALSE);'

```

### Run Application

```bash
R -q -s -e "shiny::runApp(appDir='./app',host='0.0.0.0',port=8383,launch.browser=FALSE)"
```

## Run app locally

```r
shiny::runApp(
  host = "0.0.0.0",
  port = 8080,
  launch.browser = FALSE
)
```

## To build Docker container

```bash
docker buildx build -t pdxplorer-app:latest .
```

To run Docker container:

```bash
docker run -d -p <EXTERNAL_PORT>:80 --name pdxplorer-app pdxplorer-app:latest;
```

## Aknowledges and Citation

This web application is part of the following publication, and the software is derived from work and data of the following contributors.

```
A patient-derived T-cell lymphoma biorepository uncovers new pathogenetic mechanisms and host-related therapeutic vulnerabilities.
Danilo Fiore1,2,3,41, Luca Vincenzo Cappelli1,4,41 Liu Zhaoqi5,6,7, Nikita Kotlov8, Maria Sorokina8, Jude Phillip9,10, Paul Zumbo11,12, Liron Yoffe1,13, Paola Ghione14, Anqi Wang5, Xueshuai Han5,6,7, William Chiu1, Valentina Fragliasso1,15, Fabrizio Tabbo1,16, Nahuel Zamponi9, Clarisse Kayembe1, Marcello Gaudiano1,, Rodolfo Machiorlatti17, Giuseppina Astone1, Maria Teresa Cacciapuoti1, Sanjay Patel1, Francesca Zammarchi18, Claudio Zanettini1, Lucio Queiroz1, Anastasia Nikitina8, Olga Kudryasho-va8, Anton Karelin8, Daniil Nikitin8, Dmitry Tychinin8, Ekaterina Postovalova8, Alexander Bagaev8, Viktor Svekolkin8, Ekaterina Belova8, Katerina Tikhonova8, Sandrine Degryse8, Domenico Novero19, Maurilio Ponzoni20, Enrico Tiacci21, Brunangelo Falini21, Joo Song22, Inna Khodos23, Elisa De Stanchina23 Gabriele Macari24, Luciana Cafforio24, Simone Gardini24, Roberto Piva25,26, Enzo Medi-co27, Samuel Y Ng28,29, Allison Moskowitz14, Zachary Epstein14, Andrew Intlekofer14, Dogan Ahmed30, Peter Martin31, Jia Ruan31, Francesco Bertoni32,33, Robin Foà4, Joshua D. Brody34,35,36,37, Jaspreet Osan38, Laura Santambrogio38, Doron Betel9,11,12, Wing C. Chan22, Wayne Tam1,39, David M. Wein-stock28,40, Leandro Cerchietti9, Raul Rabadan5, Steven Horwitz14, Giorgio Inghirami1*

1Pathology and Laboratory Medicine, New York Presbyterian Hospital, Weill Cornell Medicine, New York, New York, US.
2Department of Molecular Medicine and Medical Biotechnology, University of Naples Federico II, Na-ples, Italy. 
3Institute for Experimental Endocrinology and Oncology, "G.Salvatore" IEOS, Consiglio Nazionale delle Ricerche (CNR), 80131 Naples, Italy.
4Department of Translational and Precision Medicine, Sapienza University of Rome, Rome, Italy.
5Program for Mathematical Genomics, Department of Systems Biology, Department of Biomedical Informatics, Columbia University, New York, NY 10027 US. 
6Chinese Academy of Sciences Key Laboratory of Genomic and Precision Medicine, Beijing Institute of Genomics, Chinese Academy of Sciences and China National Center for Bioinformation, Beijing 100101, China.
7University of Chinese Academy of Sciences, Beijing 100049, China. 
8BostonGene Corporation, Waltham, Massachusetts, US.
9Division of Hematology and Medical Oncology, Weill Cornell Medicine, New York, NY, US.
10Chemical and Biomolecular Engineering, Oncology, Sidney Kimmel Comprehensive Cancer Cen-ter, Core Member, Institute for Nanobiotechnology (INBT), Whiting School of Engineering, Johns Hopkins University, US.
11Applied Bioinformatics Core, Weill Cornell Medicine, New York, NY, US.
12Department of Physiology and Biophysics, Weill Cornell Medicine, New York, NY, US.
13Englander Institute for Precision Medicine, Institute for Computational Biomedicine, Weill Cornell Medicine, US.
14Department of Medicine, Lymphoma Service, Memorial Sloan Kettering Cancer Center, New York, NY, US.
15Laboratory of translational research, Azienda USL – IRCCS di Reggio Emilia, Italy.
16SC Oncologia ASL CN2 Alba Bra Ospedale Michele e Pietro Ferrero - Verduno (CN), Italy.
17Department of Pathology, Center for Experimental Research and Medical Studies, University of Torino, Torino, Italy.
18ADC Therapeutics (UK) Limited, London, United Kingdom.
19Division of Pathological Anatomy, Quality and Safety of Diagnosis and Treatment, Città della Salute e della Scienza, Turin, Italy.
20Pathology Unit, San Raffaele Scientific Institute, Milan, Italy; Unit of Lymphoid Malignancies, San Raffaele Scientific Institute, Milan, Italy.
21Institute of Hematology, University of Perugia, Ospedale S. Maria della Misericordia, S. Andrea del-le Fratte, Perugia, 06156 Italy.
22Department of Pathology, City of Hope Medical Center, Duarte CA, 91010, US. 
23Antitumor Assessment Core Facility, Memorial Sloan Kettering Cancer Center, New York, NY, US.
24GenomeUp S.r.l., Rome, Italy. 
25Department of Molecular Biotechnology and Health Sciences, University of Turin, 10126 Turin, Ita-ly.
26Medical Genetics Unit, Città della Salute e della Scienza University Hospital, 10126 Turin, Italy.
27Department of Oncology, University of Torino, Candiolo, TO, Italy; Candiolo Cancer Institute, FPO-IRCCS, Candiolo, TO, Italy.
28Department of Medical Oncology, Dana-Farber Cancer Institute, 450 Brookline Ave, Boston, MA, 02215, US.
29National Cancer Institute, Bethesda, MD, US.
30Department of Pathology, Memorial Sloan-Kettering Cancer Center, 1275 York Avenue, New York, NY 10065, US.
31Lymphoma Service, Weill Cornell Medical Center, New York, NY, US.
32Lymphoma Genomics, Institute of Oncology Research, Faculty of Biomedical Sciences, USI, 6500 Bellinzona, Switzerland.
33Oncology Institute of Southern Switzerland, EOC,6500 Bellinzona, Switzerland.
34Department of Medicine, Hematology and Medical Oncology, Icahn School of Medicine at Mount Sinai, New York, NY, US.
35Tisch Cancer Institute, Icahn School of Medicine at Mount Sinai, New York, NY, US.
36Department of Oncological Sciences, Icahn School of Medicine at Mount Sinai, New York, NY, US.
37Precision Immunology Institute, Icahn School of Medicine at Mount Sinai, New York, NY, US.
38Department of Radiation Oncology, Weill Cornell Medicine, New York, NY, US. 
39Division of Hematopathology, Northwell Health, New York, NY, US.
40Merck Research Laboratories, Boston, MA 02115, US.
41D. Fiore and L.V. Cappelli contributed equally to this article.
```
