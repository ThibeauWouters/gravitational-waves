# Extra dependencies for this case
mod_usr.o: mod_xns.o
gmunu: mod_xns.o

# How to generate object files
%.o: %.f90
	$(F90) $(F90FLAGS) -c $< -o $@ $(addprefix -I,$(INC_DIRS))
