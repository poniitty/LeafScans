library(tidyverse)
library(sf)
library(terra)
library(scales)
# install.packages("benchmarkme", lib = "/projappl/project_2003061/Rpackages/")
# devtools::install_github("atsyplenkov/centerline", lib = "/projappl/project_2003061/Rpackages/")
library(smoothr, lib.loc = "/projappl/project_2003061/Rpackages/")
library(centerline, lib.loc = "/projappl/project_2003061/Rpackages/")
library(geos, lib.loc = "/projappl/project_2003061/Rpackages/")
library(Gmedian)
library(centerline, lib.loc = "/projappl/project_2003061/Rpackages/")
library(zoo)
library(nngeo, lib.loc = "/projappl/project_2003061/Rpackages/")
library(moments)

terraOptions(memfrac=0.1)

my_median = function(x) {
  if(length(x) %% 2 == 0L) { return(median(x, na.rm = TRUE)) }
  if(length(x) %% 2 == 1L) { return(median(x[-ceiling(0.5*length(x))], na.rm = TRUE)) }
}

scan_path <- "/scratch/project_2007415/leaf_scans/2024"

imgs <- list.files(scan_path, pattern = ".jpg$", full.names = TRUE, recursive = TRUE)

for(i in sample(imgs)){
  # i <- imgs[1]
  # i <- imgs[30]
  # i <- imgs[3]
  # i <- imgs[14]
  # i <- imgs[1324]
  print(i)
  
  if(!file.exists(gsub(".jpg",".csv",i))){
    # i <- imgs[5]
    r <- rast(i)
    # plot(r)
    names(r) <- c("R","G","B")
    
    # Get resolution from meta data
    rmeta <- terra::meta(r) %>% as.data.frame() %>% pull(X2) %>% as.list
    names(rmeta) <- terra::meta(r) %>% as.data.frame() %>% pull(X1)
    
    xdim <- rmeta$EXIF_PixelXDimension %>% parse_number
    xres <- rmeta$EXIF_XResolution %>% parse_number
    
    reso <- (xdim/xres*2.53999863)/xdim
    
    if(!file.exists(gsub(".jpg","_mask.tif",i))){
      # Needed indices
      r$RCC <- r$R / (r$R + r$G + r$B)
      r$NDYI <- (r$G - r$B) / (r$G + r$B)
      
      m <- r[[3]]
      names(m) <- "mask"
      
      ma <- st_bbox(m) %>% 
        st_as_sfc %>% 
        st_as_sf %>% 
        st_cast("LINESTRING") %>% 
        st_buffer(50)
      
      mar <- st_bbox(m) %>% 
        st_as_sfc %>% 
        st_as_sf %>% 
        st_cast("LINESTRING") %>% 
        st_buffer(20) %>% 
        mutate(value = 1) %>% 
        rasterize(., m, "value")
      
      # Simple mask
      m <- ifel(m < 180, 1, NA)
      m[r$R < 30] <- NA
      m[r$NDYI < 0.13 & mar == 1] <- NA
      
      pol1 <- terra::as.polygons(m) %>% 
        st_as_sf %>% 
        setNames(c("value", "geometry")) %>% 
        filter(value == 1) %>% 
        st_cast("POLYGON") %>% 
        rowid_to_column("id") %>% 
        mutate(a = st_area(.)) %>% 
        filter(a > 200)
      exl <- st_contains(ma, pol1) %>% unlist
      if(length(exl) > 0){
        pol1 <- pol1[-exl,]
      }
      pol1 <- fill_holes(pol1, threshold = 1000)
      
      gc()
      #
      polf <- pol1 %>% slice(0)
      for(ii in seq_len(nrow(pol1))){
        # ii <- 2
        pt <- pol1 %>% slice(ii)
        # plot(pt[,"id"])
        
        mar <- pt %>% st_cast("LINESTRING") %>% st_buffer(10)
        rt <- crop(r, mar)
        # plot(rt)
        # plot(pt[,"id"] %>% st_geometry(), add = TRUE)
        rt <- mask(rt, rasterize(pt %>% st_buffer(10), rt, "value"))
        mart <- rasterize(mar, rt, "value")
        
        points <- st_coordinates(pt) %>% as.data.frame()
        shp <- diff(st_bbox(pt)[c("xmin","xmax")])/diff(st_bbox(pt)[c("ymin","ymax")])
        if(shp > 1){
          epts <- bind_rows(points %>% slice_max(X) %>% slice(round(nrow(.)/2)),
                            points %>% slice_min(X) %>% slice(round(nrow(.)/2))) %>% 
            st_as_sf(coords = c("X", "Y")) %>% 
            geos::as_geos_geometry()
        } else {
          epts <- bind_rows(points %>% slice_max(Y) %>% slice(round(nrow(.)/2)),
                            points %>% slice_min(Y) %>% slice(round(nrow(.)/2))) %>% 
            st_as_sf(coords = c("X", "Y")) %>% 
            geos::as_geos_geometry()
        }
        
        cl_buff <- 0
        success <- FALSE
        while(!success) {
          e <- try({
            cl <- cnt_path(cnt_skeleton(pt %>% 
                                          st_buffer(cl_buff) %>% 
                                          geos::as_geos_geometry()), 
                           start_point = epts[1], 
                           end_point = epts[2]) %>% 
              st_as_sf
          }, silent = TRUE)
          if(class(e)[1] == "try-error"){
            cl_buff <- cl_buff + 10
          } else {
            success <- TRUE
          }
        }
        
        # plot(st_geometry(pt))
        # plot(epts, add = TRUE)
        # plot(cl, add = TRUE, col = "red")
        
        dr <- rasterize(cl %>% mutate(var = 2) %>% vect,
                        rt, "var")
        rm(cl)
        dr[is.na(dr)] <- 1
        dr <- terra::costDist(dr, target = 2)
        sp <- cos((terrain((dr)*(-1), "aspect") - 225) * pi / 180)
        filter <- matrix(1, nrow=15, ncol=15)
        filter[lower.tri(filter, diag = TRUE)] <- NA
        filter[lower.tri(filter, diag = TRUE)[,15:1]] <- NA
        fcor7 <- focalPairs(c(rt$B, dr), w = filter, fun = "pearson", na.rm = TRUE)
        fcor7[is.na(rt$B)] <- NA
        
        rt2 <- rt
        rt2[fcor7 > 0.5 & sp < 0 & rt2$NDYI < 0.15 & rt2$RCC > 0.2 & mart == 1] <- NA
        rt2[rt$B > 200] <- NA
        rt2 <- mask(rt2, rasterize(pt, rt, "value"))
        
        # plot(rt2)
        # plot(pt %>% st_geometry(), add = TRUE)
        
        rt3 <- ifel(is.na(rt2$R), 0, 1)
        rt4 <- focal(rt3, 5, mean)
        rt4 <- ifel(rt4 >= 0.5, 1, 0)
        
        polt <- terra::as.polygons(rt4) %>% 
          st_as_sf %>% 
          setNames(c("value", "geometry")) %>% 
          filter(value == 1) %>% 
          st_cast("POLYGON") %>% 
          rowid_to_column("id") %>% 
          mutate(a = st_area(.)) %>% 
          filter(a > 200)
        
        polt <- fill_holes(polt, threshold = 100)
        
        # plot(polt %>% st_simplify(dTolerance = 0) %>% st_geometry())
        # plot(pt %>% st_simplify(dTolerance = 0) %>% st_geometry(), add = TRUE)
        if(nrow(polt) > 0){
          polf <- bind_rows(polf, polt)
        }
        gc()
      }
      
      # plot(st_geometry(pol1))
      # plot(st_geometry(polf), col = "red", add = TRUE)
      
      polf <- polf %>% 
        mutate(id = seq_len(nrow(.)))
      
      r1 <- rasterize(polf, r, "value")
      # plot(r1)
      writeRaster(r1, gsub(".jpg","_mask.tif",i), datatype = "INT1U")
      unlink(tmpFiles())
      unlink(list.files(tempdir(), full.names = TRUE, recursive = TRUE))
      
    } else {
      r1 <- rast(gsub(".jpg","_mask.tif",i))
      
      polf <- terra::as.polygons(r1) %>% 
        st_as_sf %>% 
        setNames(c("value", "geometry")) %>% 
        filter(value == 1) %>% 
        st_cast("POLYGON") %>% 
        mutate(id = seq_len(nrow(.))) %>% 
        relocate(id) %>% 
        mutate(a = st_area(.)) %>% 
        filter(a > 200)
    }
    
    # Minimum rectangles ----
    
    rects <- polf %>% 
      st_minimum_rotated_rectangle(pch) %>% 
      st_cast("MULTIPOINT") %>%
      st_as_sf()
    
    rres <- lapply(seq_len(nrow(polf)), function(x){
      suppressWarnings({
        xy <- rects %>% 
          slice(x)%>% 
          st_cast("MULTIPOINT") %>%
          st_as_sf() %>% 
          st_cast("POINT")
      })
      
      seg_lengths <- st_distance(xy[-1,], xy[-nrow(xy),], by_element = TRUE)
      
      return(tibble(id = xy$id[[1]],
                    major_axis_length = max(seg_lengths)*reso,
                    minor_axis_length = min(seg_lengths)*reso,
                    lw_ratio = major_axis_length/minor_axis_length))
      
    }) %>% bind_rows()
    
    # Convex hull roundness ---------
    
    chs <- polf %>% 
      group_by(id) %>% 
      st_convex_hull()
    
    chs$convex_hull_roundness = (2*pi*sqrt((st_area(chs)*(reso^2))/pi))/(st_perimeter(chs)*reso)
    chs <- chs%>% 
      select(id, convex_hull_roundness) %>% 
      st_drop_geometry()
    
    # Basic measurements ------------
    
    bmeas <- polf %>% st_simplify(dTolerance = 1) %>% 
      mutate(perimeter = st_perimeter(.)*reso) %>% 
      mutate(area = st_area(.)*(reso^2)) %>% 
      mutate(peri2 = (2*pi*sqrt(area/pi))) %>% 
      mutate(poly_complexity = (perimeter/peri2)-1,
             dissection_index = perimeter/sqrt(area),
             interior_edge_ratio = (perimeter / area),
             shape_index = perimeter / (2 * sqrt(pi * area))) %>% 
      select(-value,-a,-peri2) %>% 
      st_drop_geometry()
    
    # Centerline & widths -----------
    
    shps <- lapply(seq_len(nrow(polf)), function(ii){
      # ii <- 3
      pt <- polf %>% slice(ii)
      
      points <- st_coordinates(pt) %>% as.data.frame()
      shp <- diff(st_bbox(pt)[c("xmin","xmax")])/diff(st_bbox(pt)[c("ymin","ymax")])
      if(shp > 1){
        epts <- bind_rows(points %>% slice_max(X) %>% slice(round(nrow(.)/2)),
                          points %>% slice_min(X) %>% slice(round(nrow(.)/2))) %>% 
          st_as_sf(coords = c("X", "Y")) %>% 
          geos::as_geos_geometry()
      } else {
        epts <- bind_rows(points %>% slice_max(Y) %>% slice(round(nrow(.)/2)),
                          points %>% slice_min(Y) %>% slice(round(nrow(.)/2))) %>% 
          st_as_sf(coords = c("X", "Y")) %>% 
          geos::as_geos_geometry()
      }
      
      cl_buff <- 0
      success <- FALSE
      while(!success) {
        e <- try({
          cl <- cnt_path(cnt_skeleton(pt %>% 
                                        st_buffer(cl_buff) %>% 
                                        geos::as_geos_geometry()), 
                         start_point = epts[1], 
                         end_point = epts[2]) %>% 
            st_as_sf
        }, silent = TRUE)
        if(class(e)[1] == "try-error"){
          cl_buff <- cl_buff + 10
        } else {
          success <- TRUE
        }
      }
      
      cl <- cl %>% 
        st_simplify(dTolerance = 10)
      
      sample_points <- sf::st_line_sample(cl, density = 0.8/(st_length(cl)*reso), type = "regular") %>% 
        sf::st_cast(to = "POINT")
      
      points_from_boundary <- cl %>% 
        sf::st_coordinates() %>% 
        as.data.frame() %>% 
        sf::st_as_sf(coords = c("X", "Y"))
      
      pls <- lapply(seq_len(length(sample_points)), function(x){
        # x <- 2
        p1 <- sample_points[x]
        p2 <- points_from_boundary$geometry[sf::st_nearest_feature(p1, points_from_boundary)]
        
        theta = base::atan2(y = st_coordinates(p1)[2] - st_coordinates(p2)[2], x = st_coordinates(p1)[1] - st_coordinates(p2)[1])  # Angle between points on the line in radians
        
        tlen <- 1000 # distance in m
        thetaT = theta+pi/2         # Get the angle
        dx_poi <- tlen*cos(thetaT) # coordinates of point of interest as defined by position length (sep)
        dy_poi <- tlen*sin(thetaT)
        
        p3 <- p1 + c(dx_poi, dy_poi)
        p4 <- p1 - c(dx_poi, dy_poi)
        
        l1 <- st_cast(st_union(p3, p4), "LINESTRING") %>% 
          st_as_sf
        
      }) %>% bind_rows()
      
      wdths <- st_intersection(pls, pt) %>% st_length
      wdths <- wdths*reso
      med_wdths <- rollapply(wdths, width = 21, FUN=my_median, fill = NA, partial = FALSE, align = "center")
      
      # plot(wdths)
      # plot(med_wdths, col = "red")
      
      plss <- na.omit(med_wdths) %>% as.numeric
      plss <- sample(seq_len(length(plss)), size = length(plss)*100, prob = plss, replace = TRUE)
      
      ress <- tibble(id = ii,
                     length = st_length(cl)*reso,
                     median_width = median(med_wdths, na.rm = TRUE),
                     mean_width = mean(med_wdths, na.rm = TRUE),
                     max_width = max(med_wdths, na.rm = TRUE),
                     widest_loc = abs(diff(c(0.5, which.max(med_wdths)/length(med_wdths)))),
                     shape_skew = abs(skewness(plss)))
      
      rm(cl)
      return(ress)
      
    }) %>% bind_rows()
    
    # Leaf color & spectral stuff ----------
    rmask <- rasterize(polf, r$R, "id")
    
    va <- values(c(rmask, r[[c("R","G","B")]]), data_frame = TRUE, na.rm = TRUE) %>% 
      as.data.frame() %>% 
      mutate(across(R:B, ~rescale(.x, to = c(0, 1), from = c(0, 255)))) %>% 
      mutate(GCC = G / (R + G + B),
             BCC = B / (R + G + B),
             RCC = R / (R + G + B),
             BITM = (((B**2.0)+(G**2.0)+(R**2.0))/3.0)**0.5,
             CRI550 = (1.0 / B) - (1.0 / G),
             DSWI4 = G/R,
             ExG = 2 * G - R - B,
             ExR = 1.3 * R - G,
             ExGR = ExG - ExR,
             GLI = (2.0 * G - R - B) / (2.0 * G + R + B),
             IKAW = (R - B)/(R + B),
             MGRVI = (G ** 2.0 - R ** 2.0) / (G ** 2.0 + R ** 2.0),
             NDTI = (R-G)/(R+G),
             NDYI = (G - B) / (G + B),
             NGRDI = (G - R) / (G + R),
             OSI = (G + R)/B,
             RGBVI = (G ** 2.0 - B * R)/(G ** 2.0 + B * R),
             RGRI = R/G,
             RI = (R - G)/(R + G),
             TGI = ((-0.5) * (190 * (R - G) - 120 * (R - B))),
             VARI = (G - R) / (G + R - B),
             VIG = (G - R) / (G + R)) %>% 
      mutate(senesced = ifelse(RCC > 0.4, 1, 0)) %>% 
      relocate(senesced, .after = id) %>% 
      mutate(across(GCC:VIG, ~ifelse(is.infinite(.x), NA, .x)))
    
    gc()
    # summary(va)
    
    cls <- lapply(unique(va$id), function(ii){
      vat <- va %>% filter(id == ii) %>% 
        select(R:B)
      
      cls <- Gmedian(vat)
      
      cls2 <- vat %>% summarise(across(R:B, median))
      cls3 <- vat %>% summarise(across(R:B, mean))
      
      rtr <- cls %>% as.data.frame() %>% tibble() %>% 
        setNames(c("R","G","B")) %>% 
        mutate(clr = rgb(cls[,1], cls[,2], cls[,3]),
               clr2 = rgb(cls2[,1], cls2[,2], cls2[,3]),
               clr3 = rgb(cls3[,1], cls3[,2], cls3[,3])) %>% 
        mutate(id = ii) %>% 
        relocate(id)
      return(rtr)
      
    }) %>% bind_rows()
    
    inds <- va %>% 
      select(-c(R:B)) %>% 
      group_by(id) %>% 
      summarise(across(everything(), ~mean(.x, na.rm = TRUE)))
    
    # Combine everything ----------
    all <- left_join(bmeas, shps) %>% 
      left_join(rres) %>% 
      left_join(chs) %>% 
      left_join(cls) %>% 
      left_join(inds) %>% 
      rename(leaf_id = id) %>% 
      mutate(scan_path = i) %>% 
      relocate(scan_path)
    
    st_write(polf %>% select(-a), gsub(".jpg",".gpkg",i))
    write_csv(all, gsub(".jpg",".csv",i))
    
    rm(va);rm(r);rm(rmask)
    unlink(list.files(tempdir(), full.names = TRUE))
    
  }
}
