# configuration

TEST_CASE = molecular_dynamics
BUILD_DIR = build

# all

TIME_STEP    = 0.000001
NUM_STEPS = 1000000
DELTA   = 1
EPSILON = 1

PARAMETERS = ${TIME_STEP} ${NUM_STEPS} ${DELTA} ${EPSILON}


targets = $(BUILD_DIR)/$(TEST_CASE)


.PHONY: all
all: mk-target-dir $(targets)


mk-target-dir:
	mkdir -p $(BUILD_DIR)

# build rules

$(BUILD_DIR)/$(TEST_CASE): $(TEST_CASE).cu
	nvcc -O3 -arch=sm_89 $(TEST_CASE).cu -o $(BUILD_DIR)/$(TEST_CASE)
	./$(BUILD_DIR)/$(TEST_CASE) $(PARAMETERS)

# aliases without build directory

.PHONY: $(TEST_CASE)
$(TEST_CASE): $(BUILD_DIR)/$(TEST_CASE)

# clean target

.PHONY: clean
clean:
	rm $(targets)