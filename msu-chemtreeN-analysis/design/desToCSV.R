# convert a design file to a csv table to be used in paraview
#
# ccs, cec24@phy.duke.edu

design.in <- as.matrix(read.table("design_sorted_wave_1.dat"))
des.names <-  c("Zr", "Fescp", "Fbary")

fname <- "design_sorted_wave_1.csv"

# we'll scale the design so that paraveiw can fit it all in a cube nicely
write.table(scale(design.in), file=fname, row.names=FALSE, col.names=des.names, sep=",", qmethod="double", dec=".")
