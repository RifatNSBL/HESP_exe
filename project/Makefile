# Configuration
TEST_CASE = molecular_dynamics
BUILD_DIR = build
OUTPUT_DIR = visualization

# Parameters
TIME_STEP = 0.003
NUM_STEPS = 1
SPRING_CONSTANT = 50
DAMPING_COEFFICIENT = 0.4
BOX_SIZE = 25
CUTOFF_RAD = 1.5
INPUT_FILE = points.txt

PARAMETERS = ${TIME_STEP} ${NUM_STEPS} ${SPRING_CONSTANT} ${DAMPING_COEFFICIENT} ${INPUT_FILE} ${CUTOFF_RAD} ${BOX_SIZE}

targets = $(BUILD_DIR)/$(TEST_CASE)

.PHONY: all
all: mk-target-dir $(targets)

mk-target-dir:
	@mkdir -p $(BUILD_DIR)

# Build rules
$(BUILD_DIR)/$(TEST_CASE): $(TEST_CASE)_submission.cu Force.h Molecule.h util.h
	rm -r $(OUTPUT_DIR) && mkdir $(OUTPUT_DIR)
	cp $(INPUT_FILE) $(BUILD_DIR)
	g++ -O3 createfile.cpp -o do_file
#	./do_file $(BOX_SIZE) $(CUTOFF_RAD)
	nvcc -O3 -arch=sm_89 $(TEST_CASE)_submission.cu -o $(BUILD_DIR)/$(TEST_CASE)
	./$(BUILD_DIR)/$(TEST_CASE) $(PARAMETERS)

# Aliases without build directory
.PHONY: $(TEST_CASE)
$(TEST_CASE): $(BUILD_DIR)/$(TEST_CASE)

# Clean target
.PHONY: clean
clean:
	@rm -f $(targets)
	@rm -rf $(OUTPUT_DIR)
