#INC_DIRS+=./
#LIB_DIRS+=./
#LIBS+=./

mod_usr.o: mod_xns.o
gmunu: mod_xns.o

%.o: %.f
	$(F90) $(F90FLAGS) -c $< -o $@ $(addprefix -I,$(INC_DIRS))
