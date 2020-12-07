 MODULE_TOPDIR = ../..
 
 # Nom de l'exécutable
 PGM = r.waterbalance
 
 LIBES = $(GISLIB) $(RASTERLIB) $(SEGMENTLIB)
 DEPENDENCIES = $(GISDEP) $(RASTERDEP) $(SEGMENTDEP)
  
 include $(MODULE_TOPDIR)/include/Make/Module.make
 
 default: cmd 