# ===================================================================
# GEOSPECTRA JOMBANG: WEBGIS INDEKS SPEKTRAL & TUTUPAN LAHAN
# Departemen Teknik Geomatika - ITS
# ===================================================================

# 1. LOAD LIBRARY WAJIB
library(shiny)
library(shinydashboard)
library(leaflet)
library(leaflet.extras)
library(sf)
library(raster)
library(tidyverse)
library(DT)

# 2. DEFINISIKAN BASE DIRECTORY UTAMA
base_dir <- "D:/KULIAH S1 TEKNIK GEOMATIKA/SEMESTER 5/PEMSPA/TUBES PEMSPA/DATA"
raster_output_dir <- file.path(base_dir, "FILE RASTER/R10m")

if (!dir.exists(raster_output_dir)) {
  dir.create(raster_output_dir, recursive = TRUE)
}

# 3. LOAD VEKTOR KECAMATAN JOMBANG (EPSG:4326 untuk Leaflet)
path_shp <- file.path(base_dir, "FILE SHP JOMBANG/shp_Jombang/KAB. JOMBANG/ADMINISTRASIDESA_AR_25K.shp")
kec_jombang <- read_sf(path_shp) %>% 
  group_by(WADMKC) %>% 
  summarise(geometry = st_union(geometry)) %>% 
  ungroup() %>% 
  st_transform(crs = 4326)


# ===================================================================
# 4. PREPROCESSING AUTOMATION (CROP, MASK, & KALKULASI)
# ===================================================================

path_ndvi <- file.path(raster_output_dir, "NDVI_Jombang.tif")
path_ndbi <- file.path(raster_output_dir, "NDBI_Jombang.tif")
path_ndwi <- file.path(raster_output_dir, "NDWI_Jombang.tif")
path_ndci <- file.path(raster_output_dir, "NDCI_Jombang.tif")
path_csv  <- file.path(raster_output_dir, "Tabel_Zonal_Statistics_Jombang.csv")

if (!file.exists(path_ndvi) || !file.exists(path_ndbi) || !file.exists(path_ndwi) || !file.exists(path_ndci) || !file.exists(path_csv)) {
  message("--> Memulai proses pengolahan, pemotongan (masking) wilayah, & statistik zonal...")
  
  # Cari folder .SAFE secara otomatis
  safe_folders <- list.files(base_dir, pattern = "\\.SAFE$", full.names = TRUE)
  if (length(safe_folders) == 0) {
    stop("Folder .SAFE tidak ditemukan di direktori: ", base_dir)
  }
  safe_folder <- safe_folders[1]
  
  # Path folder IMG_DATA
  granule_dir <- file.path(safe_folder, "GRANULE")
  sub_granule <- list.files(granule_dir)[1]
  
  r10m_path <- file.path(granule_dir, sub_granule, "IMG_DATA", "R10m")
  r20m_path <- file.path(granule_dir, sub_granule, "IMG_DATA", "R20m")
  
  # Cari file band (B02, B03, B04, B08 di R10m, B11 di R20m)
  b2_file  <- list.files(r10m_path, pattern = "_B02_10m\\.(jp2|tif)$", full.names = TRUE)
  b3_file  <- list.files(r10m_path, pattern = "_B03_10m\\.(jp2|tif)$", full.names = TRUE)
  b4_file  <- list.files(r10m_path, pattern = "_B04_10m\\.(jp2|tif)$", full.names = TRUE)
  b8_file  <- list.files(r10m_path, pattern = "_B08_10m\\.(jp2|tif)$", full.names = TRUE)
  b11_file <- list.files(r20m_path, pattern = "_B11_20m\\.(jp2|tif)$", full.names = TRUE)
  
  # 1. Baca raster acuan untuk mengambil CRS asli citra
  r_temp <- raster(b2_file)
  
  # 2. Transformasi Vektor Kecamatan ke CRS Raster
  kec_jombang_proj <- st_transform(kec_jombang, crs = crs(r_temp))
  kec_sp <- as(kec_jombang_proj, "Spatial")
  
  # 3. Load, Crop, & Masking (Persis mengikuti bentuk batas wilayah poligon)
  r_b2  <- r_temp %>% crop(kec_sp) %>% mask(kec_sp)
  r_b3  <- raster(b3_file) %>% crop(kec_sp) %>% mask(kec_sp)
  r_b4  <- raster(b4_file) %>% crop(kec_sp) %>% mask(kec_sp)
  r_b8  <- raster(b8_file) %>% crop(kec_sp) %>% mask(kec_sp)
  r_b11 <- raster(b11_file) %>% crop(kec_sp) %>% mask(kec_sp) %>% resample(r_b4, method = "bilinear")
  
  # Skala Reflektansi Level-2A (dibagi 10000)
  r_b2  <- r_b2 / 10000
  r_b3  <- r_b3 / 10000
  r_b4  <- r_b4 / 10000
  r_b8  <- r_b8 / 10000
  r_b11 <- r_b11 / 10000
  
  # Kalkulasi Indeks Spektral
  ndvi_calc <- (r_b8 - r_b4) / (r_b8 + r_b4)
  ndbi_calc <- (r_b11 - r_b8) / (r_b11 + r_b8)
  ndwi_calc <- (r_b3 - r_b8) / (r_b3 + r_b8)
  ndci_calc <- (r_b3 - r_b4) / (r_b3 + r_b4)
  
  # Simpan Hasil ke Format .tif
  writeRaster(ndvi_calc, path_ndvi, overwrite = TRUE)
  writeRaster(ndbi_calc, path_ndbi, overwrite = TRUE)
  writeRaster(ndwi_calc, path_ndwi, overwrite = TRUE)
  writeRaster(ndci_calc, path_ndci, overwrite = TRUE)
  
  # Hitung Zonal Statistics per Kecamatan
  z_ndvi <- raster::extract(ndvi_calc, kec_sp, fun = mean, na.rm = TRUE)
  z_ndbi <- raster::extract(ndbi_calc, kec_sp, fun = mean, na.rm = TRUE)
  z_ndwi <- raster::extract(ndwi_calc, kec_sp, fun = mean, na.rm = TRUE)
  z_ndci <- raster::extract(ndci_calc, kec_sp, fun = mean, na.rm = TRUE)
  
  tabel_zonal_calc <- data.frame(
    Kecamatan = kec_jombang$WADMKC,
    Rata_NDVI = round(as.vector(z_ndvi), 4),
    Rata_NDBI = round(as.vector(z_ndbi), 4),
    Rata_NDWI = round(as.vector(z_ndwi), 4),
    Rata_NDCI = round(as.vector(z_ndci), 4)
  )
  write.csv(tabel_zonal_calc, path_csv, row.names = FALSE)
  
  message("--> Pengolahan selesai! File raster & tabel zonal baru berhasil dibuat.")
}


# ===================================================================
# 5. LOAD RASTER TER-DOWNSAMPLE (UNTUK VISUALISASI CEPAT)
# ===================================================================

load_fast_raster <- function(file_path) {
  r <- raster(file_path)
  if (ncell(r) > 300000) {
    r <- aggregate(r, fact = 3, fun = mean)
  }
  return(r)
}

r_ndvi <- load_fast_raster(path_ndvi)
r_ndbi <- load_fast_raster(path_ndbi)
r_ndwi <- load_fast_raster(path_ndwi)
r_ndci <- load_fast_raster(path_ndci)

v_ndvi <- getValues(r_ndvi)
v_ndbi <- getValues(r_ndbi)
v_ndwi <- getValues(r_ndwi)
v_ndci <- getValues(r_ndci)

pal_ndvi <- colorNumeric("YlGn", domain = c(-1, 1), na.color = "transparent")
pal_ndbi <- colorNumeric("YlOrRd", domain = c(-1, 1), na.color = "transparent")
pal_ndwi <- colorNumeric("Blues", domain = c(-1, 1), na.color = "transparent")
pal_ndci <- colorNumeric("PuBuGn", domain = c(-1, 1), na.color = "transparent")

tabel_zonal <- read.csv(path_csv)


# ===================================================================
# 6. USER INTERFACE (UI)
# ===================================================================

ui <- dashboardPage(
  skin = "green",
  dashboardHeader(
    title = "GeoSPECTRA Jombang",
    titleWidth = 250,
    tags$li(
      class = "dropdown",
      style = "padding: 5px 20px; color: white; text-align: right;",
      div(style = "font-weight: bold; font-size: 14px; line-height: 1.2;", 
          "Monitoring Indeks Spektral Citra Sentinel-2 - Kabupaten Jombang"),
      div(style = "font-size: 11px; opacity: 0.85; margin-top: 2px;", 
          "Departemen Teknik Geomatika-ITS")
    )
  ),
  
  dashboardSidebar(
    width = 250,
    sidebarMenu(
      menuItem("Peta Interaktif", tabName = "map_tab", icon = icon("map")),
      menuItem("Tabel Zonal Statistics", tabName = "table_tab", icon = icon("table"))
    )
  ),
  
  dashboardBody(
    tabItems(
      # TAB 1: PETA INTERAKTIF
      tabItem(
        tabName = "map_tab",
        fluidRow(
          box(
            width = 12,
            title = "Peta Interaktif Indeks Spektral Jombang",
            status = "primary",
            solidHeader = TRUE,
            leafletOutput("map", height = "78vh")
          )
        )
      ),
      
      # TAB 2: TABEL DATA
      tabItem(
        tabName = "table_tab",
        fluidRow(
          box(
            width = 12,
            title = "Tabel Statistik Indeks Spektral per Kecamatan",
            status = "success",
            solidHeader = TRUE,
            DT::dataTableOutput("zonal_table")
          )
        )
      )
    )
  )
)


# ===================================================================
# 7. SERVER LOGIC
# ===================================================================

server <- function(input, output, session) {
  
  output$map <- renderLeaflet({
    leaflet(kec_jombang) %>%
      # Basemaps
      addTiles(group = "OpenStreetMap") %>%
      addProviderTiles(providers$CartoDB.Positron, group = "CartoDB Light") %>%
      addProviderTiles(providers$Esri.WorldImagery, group = "Satelit Esri") %>%
      
      # Layer Raster Indeks Spektral
      addRasterImage(r_ndvi, colors = pal_ndvi, opacity = 0.8, group = "NDVI", maxBytes = 12 * 1024 * 1024) %>%
      addRasterImage(r_ndbi, colors = pal_ndbi, opacity = 0.8, group = "NDBI", maxBytes = 12 * 1024 * 1024) %>%
      addRasterImage(r_ndwi, colors = pal_ndwi, opacity = 0.8, group = "NDWI", maxBytes = 12 * 1024 * 1024) %>%
      addRasterImage(r_ndci, colors = pal_ndci, opacity = 0.8, group = "NDCI", maxBytes = 12 * 1024 * 1024) %>%
      
      # Legenda
      addLegend(pal = pal_ndvi, values = v_ndvi, title = "NDVI", group = "NDVI", position = "bottomright") %>%
      addLegend(pal = pal_ndbi, values = v_ndbi, title = "NDBI", group = "NDBI", position = "bottomright") %>%
      addLegend(pal = pal_ndwi, values = v_ndwi, title = "NDWI", group = "NDWI", position = "bottomright") %>%
      addLegend(pal = pal_ndci, values = v_ndci, title = "NDCI", group = "NDCI", position = "bottomright") %>%
      
      # Batas Vektor Kecamatan
      addPolygons(
        fillColor = "transparent",
        weight = 2,
        opacity = 1,
        color = "white",
        dashArray = "3",
        highlightOptions = highlightOptions(
          weight = 3,
          color = "#ffff00",
          dashArray = "",
          bringToFront = TRUE
        ),
        label = ~WADMKC,
        popup = ~paste0("<b>Kecamatan: </b>", WADMKC),
        group = "Batas Kecamatan"
      ) %>%
      
      # Kontrol Layer
      addLayersControl(baseGroups = c("Satelit Esri", "CartoDB Light", "OpenStreetMap"),
                       overlayGroups = c("Batas Kecamatan", "NDVI", "NDBI", "NDWI", "NDCI"),
                       options = layersControlOptions(collapsed = FALSE)) %>%
      
      hideGroup(c("NDBI", "NDWI", "NDCI")) %>%
      
      addScaleBar(position = "bottomleft") %>%
      addFullscreenControl()
  })
  
  output$zonal_table <- DT::renderDataTable({
    DT::datatable(
      tabel_zonal,
      options = list(pageLength = 10, scrollX = TRUE),
      rownames = FALSE
    )
  })
}


# ===================================================================
# 8. RUN APPLICATION
# ===================================================================

shiny::shinyApp(ui = ui, server = server) 

