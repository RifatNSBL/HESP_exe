# Configuration
TEST_CASE = molecular_dynamics
BUILD_DIR = build
OUTPUT_DIR = visualization

# Parameters
TIME_STEP = 0.025
NUM_STEPS = 10000
BOX_SIZE = 100
INPUT_FILE = points.txt
POISSON_RATIO = 0.49
YOUNG_MODULUS = 0.05
BALL_SIZE = 1

PARAMETERS = ${TIME_STEP} ${NUM_STEPS} ${INPUT_FILE} ${BOX_SIZE} ${POISSON_RATIO} ${YOUNG_MODULUS}

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
	./do_file $(BOX_SIZE) $(BALL_SIZE)
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
