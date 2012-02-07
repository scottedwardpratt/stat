## ccs, 06.02.2012
## common filepaths that are used in every script here, source this file first

lumOutputFile  <- paste("./",mwString, "/output/lum_fun_outps_", mwString.Up, "_5par.dat", sep="")
metOutputFile <- paste("./", mwString, "/output/metallicity_MV_outputs_", mwString.Up, "_5par.dat", sep="")
designFile <- paste("./", mwString, "/design/design_", mwString.Up, "_5par_sorted.dat", sep="")
## this is a 5 param system, so we have more names than before..
desNames <- c("Zr", "Fescp", "Fbary", "sfe", "yfe2")

if(exists("mwStringCompare")){
  expDataFileLum <- paste("./", mwStringCompare, "/lum_fun_observed_", mwStringCompare.Up, "_5par.dat", sep="")
  expDataFileMet <- paste("./", mwStringCompare, "/metallicity_MV_observed_", mwStringCompare.Up, "_5par.dat", sep="")
}
