# This makefile was originally intended to be used on EXPANSE but has since been modified to work on my home computer

FC           = gfortran 
FLAGS        = -fdefault-real-8 -fdefault-double-8
LIBS 	     = -L/usr/lib
INCLUDES     = -I/usr/include

srms: 
	$(FC) $(FLAGS) $(INCLUDES) $(LIBS) -lnetcdff x-avg-justin.f90 -o srms

autocorel:  
	$(FC) $(FLAGS) $(INCLUDES) $(LIBS) -lnetcdff autocorel-PADDI-fromsimdat.f90 -o autocorel -O3

autocorel-saniya:  
	$(FC) $(FLAGS) $(INCLUDES) $(LIBS) -lnetcdff autocorel-PADDI-fromsimdat-saniya-fast.f90 -o autocorel-saniya -O3

extract:
	$(FC) $(FLAGS) $(INCLUDES) $(LIBS) -lnetcdff netcdf-extract.f90 -o extract -O3

zaverage:
	$(FC) $(FLAGS) $(INCLUDES) $(LIBS) -lnetcdff vertical_average_pert.f90 -o zaverage -O3

pdf:
	$(FC) $(FLAGS) $(INCLUDES) $(LIBS) -lnetcdff pdf_eta.f90 -o pdf_eta -O3

wrms:
	$(FC) $(FLAGS) $(INCLUDES) $(LIBS) -lnetcdff wrms_vorticity_weighted.f90 -o wrms -O2

vort_plot:
	$(FC) $(FLAGS) $(INCLUDES) $(LIBS) -lnetcdff -lnetcdf data_extraction.f90 -o vort_plot -O2

