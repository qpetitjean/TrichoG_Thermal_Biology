



########################################
# Install needed package               #
########################################

if(!require(data.table)){
  install.packages("data.table")
}
if(!require(knitr)){
  install.packages("knitr")
}
if(!require(kableExtra)){
  install.packages("kableExtra")
}
if(!require(ggplot2)){
  install.packages("ggplot2")
}
if(!require(ggspatial)){
  install.packages("ggspatial")
}
if(!require(ggrepel)){
  install.packages("ggrepel")
}

########################################
# Specify some path                    #
########################################

ImportDir <- "W:/Postdoc_INRAE_SAM/Data_Bidime"
  
SavingDir <- "W:/Postdoc_INRAE_SAM/Draft/BIDIME/FigsAndTables"

########################################
# Import final dataset                 #
########################################

dat = as.data.frame(data.table::fread(file.path(ImportDir, "FullDatasetBidime.csv"),
                                      sep = ";",
                                      dec = ".",
                                      na = "NA"
))

####################################################################################################
# Generate the table 1 - Summary of Trichogramma Strains and Capture Locations (WGS 84)            #
####################################################################################################

# Select relevant columns
datSubs <- dat[order(dat$Sp, dat$Strain), c("Sp", "Strain", "Haplotype", "Annee.capture", "Lon", "Lat")]
rownames(datSubs) <- NULL

# Generate HTML table
StrainSummaryTable <- kableExtra::kable_classic_2(
  kableExtra::kbl(
    datSubs,
    align = "c",
    col.names = c("Species", "Strain", "Haplotype", "Year of Capture", "Longitude", "Latitude"),
    caption = "<span style='font-size:12px; font-weight: bold; font-style: italic'>Summary of Trichogramma Strains and Capture Locations (WGS84)</span>"
  ),
  bootstrap_options = c("striped", "hover", "condensed", "responsive"),
  full_width = FALSE,
  html_font = "arial",
  font_size = 10
)

# save it
kableExtra::save_kable(StrainSummaryTable, file.path(SavingDir, "StrainSummaryTable.HTML"))

####################################################################################################
# Generate the fig 1 - The map of the locations of Trichogramma Strains by Species                 #
####################################################################################################

# Convert dataframe to sf object
trichoSf <- sf::st_as_sf(dat, coords = c("Lon", "Lat"), crs = 4326)

# Define main map region
xlim = c(-17, 11)
ylim = c(26, 51.5)

# Define inset region
inset_xlim <- c(3.5, 7.5)
inset_ylim <- c(43, 45.2)

# Separate data for labels: inset vs main
inInset <- dat$Lon >= inset_xlim[1] & dat$Lon <= inset_xlim[2] &
  dat$Lat >= inset_ylim[1] & dat$Lat <= inset_ylim[2]
dat_inset <- dat[inInset, ]
dat_main <- dat[!inInset, ]

# Get world map
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

# Main map
mainMap <- ggplot2::ggplot(data = world) +
  ggplot2::geom_sf(fill = "gray95") +
  ggplot2::geom_sf(data = trichoSf, ggplot2::aes(color = Sp), size = 3, shape = 16) +
  ggrepel::geom_text_repel(data = dat_main,
                            ggplot2::aes(x = Lon, y = Lat, label = Strain),
                            size = 3.5, max.overlaps = 15) +
  ggplot2::geom_rect(ggplot2::aes(xmin = inset_xlim[1], xmax = inset_xlim[2],
                                  ymin = inset_ylim[1], ymax = inset_ylim[2]),
                     color = "black", fill = NA, linewidth = 0.8, linetype = "solid") +
  ggplot2::geom_rect(ggplot2::aes(xmin = -17, xmax = -4.5,
                                  ymin = 26, ymax = 27.25),
                     color = "black", fill = "white", linewidth = 0.0, linetype = "solid") +
  ggspatial::annotation_scale(location = "bl", width_hint = 0.25, text_cex = 1) +
  ggspatial::annotation_north_arrow(location = "tl", style = ggspatial::north_arrow_fancy_orienteering()) +
  ggplot2::scale_color_viridis_d(option = "A", begin = 0.1, end = 0.9) +
  ggplot2::coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
  ggplot2::theme_minimal() +
  ggplot2::theme(legend.position = c(1.15, 0.8),
                 legend.background = ggplot2::element_rect(fill = "white", color = NA),
                 legend.title = ggplot2::element_text(size = 11),
                 legend.text = ggplot2::element_text(size = 9)) +
  ggplot2::labs(title = "Locations of Trichogramma Strains by Species",
                x = "Longitude", y = "Latitude",
                color = "Species")

# Inset map (zoomed in southeastern France)
insetMap <- ggplot2::ggplot(data = world) +
  ggplot2::geom_sf(fill = "gray95") +
  ggplot2::geom_sf(data = trichoSf, ggplot2::aes(color = Sp), size = 3, shape = 16) +
  ggrepel::geom_text_repel(data = dat_inset,
                           ggplot2::aes(x = Lon, y = Lat, label = Strain),
                           size = 3.5, max.overlaps = 50) +
  ggplot2::coord_sf(xlim = inset_xlim, ylim = inset_ylim, expand = FALSE) +
  ggplot2::scale_color_viridis_d(option = "A", begin = 0.1, end = 0.9) +
  ggplot2::theme_void() +
  ggplot2::theme(legend.position = "none",
                 panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.5))

# Combine main map with inset
## specify inset location
insetX <- 0.58
insetY <- 0.07
insetWidth <- 0.38
insetHeight <- 0.38

# Combine main map with inset
finalMap <- cowplot::ggdraw() +
  cowplot::draw_plot(mainMap) +
  cowplot::draw_plot(
    insetMap,
    x = insetX,
    y = insetY,
    width = insetWidth,
    height = insetHeight
  )

finalMap

# Save final map (some adjustments has been manually made using Inkskape meaning that the map in the article could slightly differ from this one)
ggsave(file.path(SavingDir, "TrichogrammaStrainMap.svg"), finalMap, width = 12, height = 8, dpi = 900)


####################################################################################################
# Test .....                #
####################################################################################################



