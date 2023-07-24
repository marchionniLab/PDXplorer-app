# PDXplorer App

## To build Docker container

```bash
docker buildx build -t pdxplorer-app:latest .
```

To run Docker container:

```bash
docker run -d -p <EXTERNAL_PORT>:80 --name pdxplorer-app pdxplorer-app:latest;
```

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


## TODO

- [ ] Move data to a separate repo or data package.
- [ ] Remove dependency on `reshape2`, `data.table`, and `magrittr`.
- [ ] Fix warnings to replace `aes_string()` on `ggplot2` calls.
- [ ] Fix `size` for `linewidth` on line-based `geom_*` in `ggplot2` calls.
- [ ] Create Icon, explorer mouse hugging a cell or RNA molecule.
- [ ] Move `app/app.R` to `app/main`.
- [ ] Replace `callModules()` for `moduleServer()` to allow the use of `testServer()`.
- [ ] There is a panel without the image generation layer behind the plot.

## Run app locally

```r
shiny::runApp(
  host = "0.0.0.0",
  port = 8080,
  launch.browser = FALSE
)
```
