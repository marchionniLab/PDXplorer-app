# PDXplorer

```r
shiny::runApp()
```

## TODO

- [ ] Move data to a separate repo or data package.
- [ ] Remove dependency on `reshape2`, `data.table`, and `magrittr`.
- [ ] Fix warnings to replace `aes_string()` on `ggplot2` calls.
- [ ] Fix `size` for `linewidth` on line-based `geom_*` in `ggplot2` calls.
- [x] Fix `scran` call to C++ `cxx_fit_linear_model` function on [`findMarkers.R`](./app/R/findMarkers.R).
  - Currently using an older version of `scran` from GitHub commit "e8c73c25bade05483d60025578955c2b357f08f9".
  - Moved C++ code ABCutilities:::fit_linear_model.
- [ ] Create Icon, explorer mouse hugging a cell or RNA molecule.
- [ ] Create Login screen with [`shinymanager`](https://datastorm-open.github.io/shinymanager/).
- [ ] Move `app/app.R` to `app/main`
- [ ] Replace CallModules for
- [x] Crash on tables panel for gene explorer, further investigation
- [x] `ggjoy` warning, complaining to move to `ggridges`
- [ ] Intersections panel is crashing, under _Fusion explorer_ tab, with the error below:

```r
Warning: Error in systemPipeR::overLapper: Unexpected input.
             The input 'setlist' needs to be of class 'list' where each list component stores a
             label set as 'vector' and the name of each label set is provided under the
             name slot of each list component.
```

## Run app locally

```r
shiny::runApp(
  host = "0.0.0.0",
  port = 8080,
  launch.browser = FALSE
)
```
