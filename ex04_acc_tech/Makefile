# configuration

TEST_CASE = molecular_dynamics
BUILD_DIR = build
OUTPUT_DIR = visualization

# all

TIME_STEP    = 0.003
NUM_STEPS = 5000
SIGMA   = 1
EPSILON = 2
BOX_SIZE = 3
INPUT_FILE = points.txt

PARAMETERS = ${TIME_STEP} ${NUM_STEPS} ${SIGMA} ${EPSILON} $(INPUT_FILE) ${BOX_SIZE}


targets = $(BUILD_DIR)/$(TEST_CASE)


.PHONY: all
all: mk-target-dir $(targets)


mk-target-dir:
	mkdir -p $(BUILD_DIR)

# build rules


$(BUILD_DIR)/$(TEST_CASE): $(TEST_CASE).cu Force.h Molecule.h util.h
	rm -r $(OUTPUT_DIR) && mkdir $(OUTPUT_DIR)
	cp $(INPUT_FILE) $(BUILD_DIR)
#	g++ -O3 createfile.cpp -o do_file
#	./do_file
	nvcc -O3 -arch=sm_89 $(TEST_CASE).cu -o $(BUILD_DIR)/$(TEST_CASE)
	./$(BUILD_DIR)/$(TEST_CASE) $(PARAMETERS)

# aliases without build directory

.PHONY: $(TEST_CASE)
$(TEST_CASE): $(BUILD_DIR)/$(TEST_CASE)

# clean target

.PHONY: clean
clean:
	rm $(targets)