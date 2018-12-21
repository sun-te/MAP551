# Makefile

FC = gfortran
LD = gfortran

FCFLAGS = -J ./obj
LDFLAGS = 

EXE = run/bz_1d

OBJDIR = ./obj

OBJS = obj/mod_precision.o \
       obj/mod_cartesian_grid.o \
       obj/mod_utils.o \
       obj/mod_bz_2eq_1d.o \
       obj/mod_bz_3eq_1d.o \
       obj/radau5.o \
       obj/decsol.o \
       obj/dc_decsol.o \
       obj/mod_radau.o \
       obj/rock4.o \
       obj/rho.o \
       obj/mod_rock.o \
       obj/mod_strang.o \
       obj/mod_strang_drd.o \
       obj/mod_lie_rd.o \
       obj/mod_lie_dr.o \
       obj/mod_integration.o \
       obj/bz_1d_main.o

default: $(OBJDIR) $(EXE)

$(OBJDIR):
	@echo "Creation du repertoire $@"
	@mkdir -p $@
	@echo

$(EXE): $(OBJS)
	$(LD) $(LDFLAGS) $^ -o $@ $(LIBS)

obj/%.o : src/%.f90
	$(FC) $(FCFLAGS) -c $< -o $@

obj/%.o : src/%.f
	$(FC) $(FCFLAGS) -c $< -o $@


#dependencies

obj/mod_utils.o: obj/mod_precision.o obj/mod_cartesian_grid.o

obj/mo1_bz_2eq_1d.o: obj/mod_precision.o obj/mod_cartesian_grid.o obj/mod_csr_matrix.o

obj/mod_integration : obj/mod_precision.o obj/mod_radau.o obj/mod_strang.o obj/mod_bz_2eq_1d.o

obj/bz_2eq_1d_main.o: obj/mod_precision.o obj/mod_cartesian_grid.o obj/mod_utils.o obj/mod_bz_2eq_1d.o obj/mod_integration.o

clean:
	rm -f $(OBJS) obj/*.mod $(EXE)
