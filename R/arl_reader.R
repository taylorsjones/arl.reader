#R version of the ARL reader
require('proj4')   #needed for CRS definitions
require('raster')  #creates raster objects
require('rts')     #creates raster time series (RTS) objects


ARL_read_main_header <- function(f){
	h <- list()
	fh <- readChar(f,50)
	h$met_type <- readChar(f,4)
	h$forecast_hr <- readChar(f,3)
	h$forecast_mn <- readChar(f,2)
	h$pole_lat    <- as.numeric( readChar(f,7) )
	h$pole_lon    <- as.numeric( readChar(f,7) )
	h$ref_lat     <- as.numeric( readChar(f,7) )
	h$ref_lon     <- as.numeric( readChar(f,7) )
	h$dx          <- as.numeric( readChar(f,7) )
	h$orientation <- readChar(f,7)
	h$cone_angle  <- as.numeric( readChar(f,7) )
	h$x_sync      <- as.numeric( readChar(f,7) )
	h$y_sync      <- as.numeric( readChar(f,7) )
	h$sync_lat    <- as.numeric( readChar(f,7) )
	h$sync_lon    <- as.numeric( readChar(f,7) )
	h$reserved    <- readChar(f,7)
	h$nx          <- as.numeric( readChar(f,3) )
	h$ny          <- as.numeric( readChar(f,3) )
	if( h$met_type == "HRRR"){
		h$nx <- h$nx + 1000
		h$ny <- h$ny + 1000
	}
	h$nz          <- as.numeric( readChar(f,3) )
	h$v_coord     <- readChar(f,2)
	h$index_length<- as.numeric( readChar(f,4) )
	d2_vars <- NULL
	z <- NULL
	for( j  in 1:h$nz){
		d3_vars <- NULL
		z <- c(z,readChar(f,6))
		num_of_vars <- as.numeric( readChar(f,2) )
		for( i in 1:num_of_vars){
			if( j == 1){
				d2_vars <- c(d2_vars, readChar(f,4) )
			}else{
				d3_vars <- c(d3_vars, readChar(f,4) )
			}
			readChar(f,3) #check sum...ignore it.
			readChar(f,1) #reserved
		}
	}
	h$z <- z
	h$d2_vars <- d2_vars
	h$d3_vars <- d3_vars

	#build the projection:
	lat_1 = acos( h$cone_angle /360.0 )*(180/pi)
	h$proj <- paste0("+proj=lcc +lat_1=",h$cone_angle," +lat_2=",h$cone_angle," +lat_0=",h$ref_lat," +lon_0=",h$ref_lon,' +units=m +e=0 +a=6371229')
	h$xx <- seq(-1*((h$nx*h$dx*1000)-1)/2,((h$nx*h$dx*1000)-1)/2,h$dx*1000) - 520.143 #offset for HRRR
	h$yy <- seq(-1*((h$ny*h$dx*1000)-1)/2,((h$ny*h$dx*1000)-1)/2,h$dx*1000) - 306.153 #offset for HRRR

	return(h)
}

ARL_read_variable <- function(f,h,x1,x2,y1,y2){
	output <- list()
	while( readChar(f,1) == ''){}
	seek(f,-1,"current")
	datetime <- readChar(f,10)
	dont_know   <- readChar(f,4)
	var_name    <- readChar(f,4)
	scale       <- as.numeric( readChar(f,4) )
	precision   <- as.numeric( readChar(f,14) )
	last_value  <- as.numeric( readChar(f,14) )
	var_map     <- matrix(0,1+y2-y1,1+x2-x1)
	for( i in 1:h$ny){
		if(i < y1){ #skip rows early rows
			last_value <- (readBin(f,"int",n=1,size=1,signed=F)-127.0)/(2^(7-scale)) + last_value
			seek(f,h$nx-1,"current")
		}else if( i > y2){ #leave after the rows we care about
			seek(f,h$nx*(1+h$ny-i),"current")
			break
		}else{
			first_col_last_value <- last_value
			for( j in 1:h$nx){
				value <- (readBin(f,"int",n=1,size=1,signed=F)-127.0)/(2^(7-scale)) + last_value
				if(j > x1){
        	var_map[i-y1,j-x1] <- value
				}
		    if(j > x2){
		    	seek(f,h$nx-j,"current")
		    	break
		    }
		    last_value <- value
			}
		last_value <- first_col_last_value
		}
	}
	var_map <- apply(var_map,2,rev)
	output$name <- var_name
	output$map <- var_map
	return(output)
}

#' Read an ARL file
#'
#' @param file_to_read the name of the file to read in
#' @param var_i_want the name of the variable to extract
#' @param ll the lat and lon of the lower-left corner of the area to extract
#' @param ur the lat and lon of the upper-right corner of the area to extract
#' @param latlon whether to convert the raster to lat-lon grid. Otherwise it stays in the native projection
#' @param verbose whether to output additional messages during extraction'
#' @return for 2d variables, a raster time series object is returned. For 3d variables, a list of RTS objects is returned, one for each altitude.
#' @examples
#' file_to_read <- 'hysplit.20161020.00z.hrrra'
#' lower_left <- c(-73,40)
#' upper_left <- c(-70,43)
#' pbl <- ARL_read(file_to_read,"PBLH",ll=lower_left,ur=upper_right)
#' plot(pbl)
ARL_read <- function(file_to_read,var_i_want,ll,ur,latlon=TRUE, verbose=FALSE){
	f = file(file_to_read,'rb')
	raster_list <- list()
	levels <- NULL
	datetimes <- NULL
	for(step in seq(1,2000,1)){
		date <- readChar(f,10)
		level <- readChar(f,2)
		AA <- readChar(f,2)
		block_ID <- readChar(f,4)
		seek(f,-18,"current")
		if( length(block_ID) < 1 ){
			if(verbose){ message("...end of file") }
			break
		}
		if(block_ID == "INDX"){ #header detected...load in all the header info.
			h <- ARL_read_main_header(f)
			seek(f,h$nx*h$ny - h$index_length , "current")
			ll_xy <- project( ll, h$proj )
			x1 <- findInterval(ll_xy[1],h$xx)
			y1 <- findInterval(ll_xy[2],h$yy)
			ur_xy <- project( ur, h$proj )
			x2 <- findInterval(ur_xy[1],h$xx)
			y2 <- findInterval(ur_xy[2],h$yy)
			if(verbose){
				message("Array Dimensions: ")
				message(x2-x1)
				message(y2-y1)
			}
			if(var_i_want %in% h$d2_vars){
				if(verbose){message("2-d var detected")}
			}else if(var_i_want %in% h$d3_vars){
				if(verbose){message("3-d var detected")}
			}else{
				message("var not found! These vars are available:")
				message("2-dimensional:")
				print(h$d2_vars)
				message("3-dimensional:")
				print(h$d3_vars)
				break
			}

		}else if(block_ID == var_i_want){
			out <- ARL_read_variable(f,h,x1,x2,y1,y2)
			datetime <- paste0("20",gsub(" ","0",date))
			ras <- raster(out$map ,crs=h$proj, xmn=ll_xy[1], xmx=ur_xy[1], ymn=ll_xy[2], ymx=ur_xy[2] )
			if(latlon){
				ras <- projectRaster(ras,crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
			}
			raster_list <- c(raster_list,ras)
			datetimes <- c(datetimes,datetime)
			levels <- c(levels,level)
		}else{
			seek(f,50 + h$nx*h$ny,"current")
		}
	}

	close(f)

	unique_levels <- unique(levels)
	unique_datetimes <- unique(datetimes)

	if( length( unique_levels ) > 1){
		rts_list <- list()
  	output <- list()

		for( level in unique_levels){
			sub_list_of_rasters <- raster_list[ level == levels ]
			sub_list_of_datetimes <- datetimes[ level == levels ]
			rstack <- stack(sub_list_of_rasters)
			rts_object <- rts(rstack, as.POSIXct(sub_list_of_datetimes,format="%Y%m%d%H%M",tz="GMT"))
			rts_list <- c(rts_list,rts_object)
		}

		output$levels <- h$z
		output$rasters <- rts_list
		return(output)
	}else{
		rstack <- stack(raster_list)
		rts_object <- rts(rstack, as.POSIXct(datetimes,format="%Y%m%d%H%M",tz="GMT"))
		return(rts_object)
	}
}
