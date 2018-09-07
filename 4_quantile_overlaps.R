
# -----------------------------------------------------------------------
# Ecosystem service capacity-pressure-demand overlap
# -----------------------------------------------------------------------


# read in files -----------------------------------------------------------
pred_lyr <- raster('hi_otp_all_osds_nitrogen.tif') # download from: http://www.pacioos.hawaii.edu/projects/oceantippingpoints/#data
resp_lyr <- raster('recreation_raster.tif')
# pred_lyr <- projectRaster(OTP_MHI_Sedimentation,crs=projection(resp_lyr)) #need to transform the predictor layer, but this is not currently working, exported to ArcGIS and reimported after transformation

# extract values ----------------------------------------------------------
pred_out <- getValues(pred_lyr)
resp_out <- getValues(resp_lyr)

# plot --------------------------------------------------------------------
col.vec <- data.frame('col'=rep("grey",length(pred_out)),"resp"=resp_out,"pred"=pred_out,stringsAsFactors =F)
col.vec <- col.vec[!is.na(col.vec$pred),]
col.vec <- col.vec[!is.na(col.vec$resp),]
brks_pred <- quantile(LBSP_NutrientsOSDS_Nflux_prj_noZero,c(0,0.25,0.75,1))
brks_resp <- quantile(col.vec$resp,c(0,0.25,0.75,1))

col.vec$col[c(col.vec$pred < brks_pred[2] & col.vec$resp < brks_resp[2])] <- "black"
col.vec$col[c(col.vec$pred > brks_pred[3] & col.vec$resp < brks_resp[2])] <- "red"
col.vec$col[c(col.vec$pred > brks_pred[3] & col.vec$resp > brks_resp[3])] <- "blue"
col.vec$col[c(col.vec$pred < brks_pred[2] & col.vec$resp > brks_resp[3])] <- "green"
col.vec$resp2 <- col.vec$resp
col.vec$resp2[col.vec$col=="grey"] <- NA
col.vec$pred2 <- col.vec$pred
col.vec$pred2[col.vec$col=="grey"] <- NA

png(file="scatter_Nflux_v_RecAdj_061518.png",height=2000,width=2700,res=300)
par(mar=c(5,5.2,0.5,0.5))
plot(log(col.vec$pred+1),log(col.vec$resp+1),col="grey",pch=19,xaxt="n",yaxt="n",ylab="Recreation ES",xlab="Nitrogen Pollution (OSDS)",cex.axis=2.5,cex.lab=2.5,mgp=c(3.3,1.5,0))
points(log(col.vec$pred2+1),log(col.vec$resp2+1),col=col.vec$col,pch=19)
axis(1,at=c(0,log(100+1),log(500+1),log(1000+1),log(5000+1),log(8000+1)),labels=c(0,100,500,1000,5000,8000),cex.axis=2.5,mgp=c(3.1,1.5,0))
axis(2,at=c(0,log(0.5+1),log(1+1),log(1.5+1),log(2+1),log(2.5+1)),labels=c(0,0.5,1,1.5,2,2.5),cex.axis=2.5,mgp=c(3,0.9,0))
dev.off()


