# https://raw.githubusercontent.com/road2stat/ggsci/master/data-raw/data-generator.R

# from ggsci and Rcolorbrewer

ggsci_db <- vector("list")

# Discrete Color Palettes
ggsci_db$"auto"$"default" <- c(scales::hue_pal()(1))
ggsci_db$"ggplot"$"default" <- c(scales::hue_pal()(15))


# Color palette inspired by plots in Nature Reviews Cancer
ggsci_db$"npg"$"default" <- c(
  "Cinnabar"     = "#E64B35", "Shakespeare"    = "#4DBBD5",
  "PersianGreen" = "#00A087", "Chambray"       = "#3C5488",
  "Apricot"      = "#F39B7F", "WildBlueYonder" = "#8491B4",
  "MonteCarlo"   = "#91D1C2", "Monza"          = "#DC0000",
  "RomanCoffee"  = "#7E6148", "Sandrift"       = "#B09C85"
)

# Color palette inspired by plots in Science from AAAS
ggsci_db$"aaas"$"default" <- c(
  "Chambray"      = "#3B4992", "Red"           = "#EE0000",
  "FunGreen"      = "#008B45", "HoneyFlower"   = "#631879",
  "Teal"          = "#008280", "Monza"         = "#BB0021",
  "ButterflyBush" = "#5F559B", "FreshEggplant" = "#A20056",
  "Stack"         = "#808180", "CodGray"       = "#1B1919"
)

# Color palette inspired by plots in The New England Journal of Medicine
ggsci_db$"nejm"$"default" <- c(
  "TallPoppy"      = "#BC3C29", "DeepCerulean" = "#0072B5",
  "Zest"           = "#E18727", "Eucalyptus"   = "#20854E",
  "WildBlueYonder" = "#7876B1", "Gothic"       = "#6F99AD",
  "Salomie"        = "#FFDC91", "FrenchRose"   = "#EE4C97"
)

# Color palette inspired by plots in Lancet Oncology
ggsci_db$"lancet"$"default" <- c(
  "CongressBlue" = "#00468B", "Red"       = "#ED0000",
  "Apple"        = "#42B540", "BondiBlue" = "#0099B4",
  "TrendyPink"   = "#925E9F", "MonaLisa"  = "#FDAF91",
  "Carmine"      = "#AD002A", "Edward"    = "#ADB6B6",
  "CodGray"      = "#1B1919"
)

# Color palette inspired by plots in The Journal of the American Medical Association
ggsci_db$"jama"$"default" <- c(
  "Limed Spruce" = "#374E55", "Anzac"         = "#DF8F44",
  "Cerulean"     = "#00A1D5", "Apple Blossom" = "#B24745",
  "Acapulco"     = "#79AF97", "Kimberly"      = "#6A6599",
  "Makara"       = "#80796B"
)

# Color palette inspired by plots in Journal of Clinical Oncology
ggsci_db$"jco"$"default" <- c(
  "Lochmara" = "#0073C2", "Corn"         = "#EFC000",
  "Gray"     = "#868686", "ChestnutRose" = "#CD534C",
  "Danube"   = "#7AA6DC", "RegalBlue"    = "#003C67",
  "Olive"    = "#8F7700", "MineShaft"    = "#3B3B3B",
  "WellRead" = "#A73030", "KashmirBlue"  = "#4A6990"
)

# Color palette inspired by University of Chicago Color Palette
ggsci_db$"uchicago"$"default" <- c(
  "Maroon" = "#800000", "DarkGray"   = "#767676",
  "Yellow" = "#FFA319", "LightGreen" = "#8A9045",
  "Blue"   = "#155F83", "Orange"     = "#C16622",
  "Red"    = "#8F3931", "DarkGreen"  = "#58593F",
  "Violet" = "#350E20"
)



# Color palette inspired by The Simpsons
ggsci_db$"simpsons"$"default" <- c(
  "HomerYellow" = "#FED439", "HomerBlue"  = "#709AE1",
  "HomerGrey"   = "#8A9197", "HomerBrown" = "#D2AF81",
  "BartOrange"  = "#FD7446", "MargeGreen" = "#D5E4A2",
  "MargeBlue"   = "#197EC0", "LisaOrange" = "#F05C3B",
  "NedGreen"    = "#46732E", "MaggieBlue" = "#71D0F5",
  "BurnsPurple" = "#370335", "BurnsGreen" = "#075149",
  "DuffRed"     = "#C80813", "KentRed"    = "#91331F",
  "BobGreen"    = "#1A9993", "FrinkPink"  = "#FD8CC1"
)

# Color palette inspired by Futurama
ggsci_db$"futurama"$"default" <- c(
  "FryOrange"    = "#FF6F00", "FryRed"        = "#C71000",
  "FryBlue"      = "#008EA0", "LeelaPurple"   = "#8A4198",
  "BenderIron"   = "#5A9599", "ZoidbergRed"   = "#FF6348",
  "ZoidbergBlue" = "#84D7E1", "AmyPink"       = "#FF95A8",
  "HermesGreen"  = "#3D3B25", "ProfessorBlue" = "#ADE2D0",
  "ScruffyGreen" = "#1A5354", "LeelaGrey"     = "#3F4041"
)

# Color palette inspired by Rick and Morty
ggsci_db$"rickandmorty"$"default" <- c(
  "MortyYellow" = "#FAFD7C", "MortyBrown"   = "#82491E",
  "MortyBlue"   = "#24325F", "RickBlue"     = "#B7E4F9",
  "BethRed"     = "#FB6467", "JerryGreen"   = "#526E2D",
  "SummerPink"  = "#E762D7", "SummerOrange" = "#E89242",
  "BethYellow"  = "#FAE48B", "RickGreen"    = "#A6EEE6",
  "RickBrown"   = "#917C5D", "MeeseeksBlue" = "#69C8EC"
)

# Color palette inspired by Star Trek
ggsci_db$"startrek"$"default" <- c(
  "Engineering" = "#CC0C00", "Sciences" = "#5C88DA",
  "Senior"      = "#84BD00", "Command"  = "#FFCD00",
  "Teal"        = "#7C878E", "Cerulean" = "#00B5E2",
  "Jade"        = "#00AF66"
)


# R color_brewer set3
ggsci_db$"set3"$"default" <- c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")

# R color_brewer set2
ggsci_db$"set2"$"default" <- c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3")

# R color_brewer set1
ggsci_db$"set1"$"default" <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999")

# R color_brewer pastel2
ggsci_db$"pastel2"$"default" <- c("#B3E2CD", "#FDCDAC", "#CBD5E8", "#F4CAE4", "#E6F5C9", "#FFF2AE", "#F1E2CC", "#CCCCCC")

# R color_brewer pastel1
ggsci_db$"pastel1"$"default" <- c("#FBB4AE", "#B3CDE3", "#CCEBC5", "#DECBE4", "#FED9A6", "#FFFFCC", "#E5D8BD", "#FDDAEC", "#F2F2F2")

# R color_brewer dark
ggsci_db$"dark"$"default" <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666")

# R color_brewer accent
ggsci_db$"accent"$"default" <- c("#7FC97F", "#BEAED4", "#FDC086", "#FFFF99", "#386CB0", "#F0027F", "#BF5B17", "#666666")
