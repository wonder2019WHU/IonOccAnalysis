all: ip2rep ip2info ipfromdig ipmatch 

#compiler
CPP = g++
	
#install path
INSTDIR = ../bin
	
#include path
INCLUDEDIR = ./include

#libray path

#library name
LIB1 = netcdf_c++4
LIB2 = netcdf
	
#development options
VERSION = -std=c++11
	
#compile rules
#bin
ip2rep: ip2rep.o IonPrf.o Coordinate.o
	$(CPP) $(VERSION) ip2rep.o IonPrf.o Coordinate.o -o ip2rep -l$(LIB1) -l$(LIB2) 
	
ip2info: ip2info.o IonPrf.o Coordinate.o
	$(CPP) $(VERSION) ip2info.o IonPrf.o Coordinate.o -o ip2info -l$(LIB1) -l$(LIB2) 
	
ipfromdig: ipfromdig.o IonPrf.o Coordinate.o
	$(CPP) $(VERSION) ipfromdig.o IonPrf.o Coordinate.o -o ipfromdig -l$(LIB1) -l$(LIB2) 
	
ipmatch: ipmatch.o IonPrf.o Coordinate.o
	$(CPP) $(VERSION) ipmatch.o IonPrf.o Coordinate.o -o ipmatch -l$(LIB1) -l$(LIB2) 
#main.o	
ip2rep.o: ip2rep.cpp
	$(CPP) $(VERSION) ip2rep.cpp -c -I$(INCLUDEDIR)
	
ip2info.o: ip2info.cpp
	$(CPP) $(VERSION) ip2info.cpp -c -I$(INCLUDEDIR)
	
ipfromdig.o: ipfromdig.cpp
	$(CPP) $(VERSION) ipfromdig.cpp -c -I$(INCLUDEDIR)
	
ipmatch.o: ipmatch.cpp
	$(CPP) $(VERSION) ipmatch.cpp -c -I$(INCLUDEDIR)
	
#help.o	
IonPrf.o: IonPrf.cpp
	$(CPP) $(VERSION) IonPrf.cpp -c -I$(INCLUDEDIR)
	
Coordinate.o: Coordinate.cpp
	$(CPP) $(VERSION) Coordinate.cpp -c -I$(INCLUDEDIR)

install: ip2rep ip2info ipfromdig ipmatch 
	@if [ ! -d $(INSTDIR) ]; then \
		@mkdir $(INSTDIR); \
	fi
	@mv ip2rep $(INSTDIR)
	@mv ip2info $(INSTDIR)
	@mv ipfromdig $(INSTDIR)
	@mv ipmatch $(INSTDIR)
	@echo "Install Finished"
	
clean:
	-rm ip2rep.o ip2info.o ipfromdig.o ipmatch.o IonPrf.o Coordinate.o

	

	