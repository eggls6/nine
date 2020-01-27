#Makefile for p-type hierarchical triple simualtion in 2d 


#compiler
cf=gfortran
#compilder flags
cflags=-O  
#program output name
pg=nine.exe

modk=global_m.o
mod0=transform_m.o
mod1=yark_m.o
mod2=input_m.o
mod3=output_m.o
mod4=planetoid_m.o
mod5=cutoffmerge_m.o
mod6=radau_m.o
mod7=symp_m.o
mod8=radeih_m.o
mod9=lie_m.o
mod10=bs_m.o
mod11=main.o


vpath  %.f src
vpath  %.f90 src

obj=$(modk) $(mod0) $(mod1) $(mod2) $(mod3) $(mod4) $(mod5) $(mod6) $(mod7) \
    $(mod8) $(mod9) $(mod10) $(mod11) 


targ:$(obj)
	$(cf) -o $(pg) $(obj)
	
%.o: %.f
	$(cf) ${cflags} -c   $<    
    
%.o: %.f90
	$(cf) ${cflags} -c $<


    
# Utility targets
.PHONY: clean veryclean

clean:
	rm -f *.o *.mod *.MOD

veryclean:       
	rm -f *.o *.mod *.MOD ./src/*.o ./src/*.mod



