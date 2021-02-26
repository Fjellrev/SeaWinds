 
  # This function converts angles given in Radians **relative to the X-Axis**
  # into values in degrees relative to the North (Y-Axis)
  # or the reverse (if using 'deg2rad' with conv argument)
  # 
  # For use with adehabitatLT

  Conv.Radians <- function (MyAngles, conv = c("deg2rad", "rad2deg")) {
  
    # require (CircStats)
    
    if( conv == "rad2deg"){
      tmp          <- MyAngles[!is.na(MyAngles)]
      tmp[tmp>0]   <- tmp[tmp>0] - 2*pi
      tmp          <- tmp - pi/2
      tmp          <- abs(deg(tmp))
      tmp[tmp>360] <- tmp[tmp>360] -360
      
      MyAngles[!is.na(MyAngles)] <- tmp      
    } else {     
      if( conv == "deg2rad"){
      
        tmp               <- MyAngles[!is.na(MyAngles)]
        tmp               <- 90 - tmp
        tmp               <- rad(tmp)
        tmp[tmp<= (-pi)] <- tmp[tmp<= (-pi)] + 2*pi
      
        MyAngles[!is.na(MyAngles)] <- tmp
      }
    }
            
    return (MyAngles)
    
  }
  
  

  


