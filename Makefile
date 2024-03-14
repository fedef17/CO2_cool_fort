# Compiler
FC = ifort
#FC = gfortran

# Flags
FFLAGS = -O3 -m64

# Directories
SRC_DIR = source
OBJ_DIR = obj
BIN_DIR = .

# Module sources
MODULE_SOURCES = precision.f90 constants.f90 coedat.f90 varsub.f90 co2cool.f90
MODULE_OBJECTS = $(patsubst %.f90,$(OBJ_DIR)/%.o,$(MODULE_SOURCES))

# Main source
MAIN_SOURCE = main.f90
MAIN_OBJECT = $(OBJ_DIR)/main.o

# Output file
OUTPUT = $(BIN_DIR)/run_cool

# Targets
all: $(OUTPUT)

$(OUTPUT): $(MODULE_OBJECTS) $(MAIN_OBJECT)
	$(FC) $(FFLAGS) $^ -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/modules/%.f90 | $(OBJ_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/main.o: $(SRC_DIR)/$(MAIN_SOURCE) $(MODULE_OBJECTS) | $(OBJ_DIR)
	$(FC) $(FFLAGS) -c $< -o $@ -I$(SRC_DIR)/modules

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

clean:
	rm -rf $(OBJ_DIR) $(OUTPUT)

.PHONY: all clean
