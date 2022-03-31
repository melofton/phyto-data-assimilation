#Function to load and re-format NOAA forecasts
#Courtesy of the fabulous T.N. Moore, developed for Macrosystems EDDIE module use

#' load NOAA GEFS forecast
#'
#' @param siteID name of lake site
#' @param start_date start date for forecast
#' @param my_directory filepath for local directory where forecasts are stored;
#' the assumption is that within this directory, you can then retrieve forecasts
#' using the key for each forecast from the S3 bucket (aka your local directory
#' file structure is the same as the folder structure on the S3 bucket)

format_noaa_forecast <- function(siteID, start_date, my_directory){
  
  fpath <- file.path(my_directory, "noaa/NOAAGEFS_1hr", siteID)
  fpath2 <- file.path(fpath, start_date, "00")
  fils <<- list.files(fpath2)
  fils <<- fils[-c(grep("ens00", fils))]
  fid <- nc_open(file.path(fpath2, fils[1]))
  vars <- fid$var # Extract variable names for selection
  fc_vars <<- names(vars)
  membs <<- length(fils)
  
  fc <- lapply(start_date, function(dat) {
    idx <- which(start_date == dat)
    
    fpath2 <- file.path(fpath, dat, "00")
    fils <- list.files(fpath2)
    fils <- fils[-c(grep("ens00", fils))]
    
    for( i in seq_len(length(fils))) {
      
      fid <- ncdf4::nc_open(file.path(my_directory, "noaa/NOAAGEFS_1hr", siteID, dat,
                                      "00", fils[i]))
      tim = ncvar_get(fid, "time")
      tunits = ncatt_get(fid, "time")
      lnam = tunits$long_name
      tustr <- strsplit(tunits$units, " ")
      step = tustr[[1]][1]
      tdstr <- strsplit(unlist(tustr)[3], "-")
      tmonth <- as.integer(unlist(tdstr)[2])
      tday <- as.integer(unlist(tdstr)[3])
      tyear <- as.integer(unlist(tdstr)[1])
      tdstr <- strsplit(unlist(tustr)[4], ":")
      thour <- as.integer(unlist(tdstr)[1])
      tmin <- as.integer(unlist(tdstr)[2])
      origin <- as.POSIXct(paste0(tyear, "-", tmonth, 
                                  "-", tday, " ", thour, ":", tmin), 
                           format = "%Y-%m-%d %H:%M", tz = "UTC")
      if (step == "hours") {
        tim <- tim * 60 * 60
      }
      if (step == "minutes") {
        tim <- tim * 60
      }
      time = as.POSIXct(tim, origin = origin, tz = "UTC")
      var_list <- lapply(fc_vars, function(x) {
        data.frame(time = time, value = ncdf4::ncvar_get(fid, x))
      }) 
      
      
      ncdf4::nc_close(fid)
      names(var_list) <- fc_vars
      
      mlt1 <- reshape::melt(var_list, id.vars = "time")
      mlt1 <- mlt1[, c("time", "L1", "value")]
      
      # df <- get_vari(file.path("data", fils[i]), input$fc_var, print = F)
      cnam <- paste0("ens", formatC(i, width = 2, format = "d", flag = "0"))
      if(i == 1) {
        df2 <- mlt1
        colnames(df2)[3] <- cnam
      } else {
        df2 <- merge(df2, mlt1, by = c(1,2))
        colnames(df2)[ncol(df2)] <- cnam
      }
      
    }
    return(df2)
  })
  
  names(fc) <- start_date
  fc_idx <- fc[[start_date]]
  
  fc_conv_list <- lapply(1:30, function(x) {
    df <- fc_idx
    sub <- df[(df[, 2] %in% fc_vars), c(1, 2, 2 + x)]
    df2 <- tidyr::pivot_wider(data = sub, id_cols = time, names_from = L1, values_from = 3)
    df2$air_temperature <- df2$air_temperature - 273.15
    df2$date <- as.Date(df2$time)
    df2$time <- NULL
    df3 <- plyr::ddply(df2, "date", function(x){
      colMeans(x[, 1:9], na.rm = TRUE)
    })
    # df3 <- df3[2:16, ]
    fc_out_dates <<- df3$date
    
    #df3 <- df3[, fc_vars]
    df3$fc_date <- start_date
    return(df3)
  })
  
  l1 <- fc_conv_list
  idvars <- colnames(l1[[1]])
  mlt1 <- reshape::melt(l1, id.vars = idvars)
  mlt2 <- mlt1 %>%
    rename(issue_date = fc_date,
           fc_date = date,
           ensemble_member = L1)
  
  
  return(mlt2)
  
  
}