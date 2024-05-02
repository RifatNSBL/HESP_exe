# configuration

TEST_CASE = stream
BUILD_DIR = ../../build/$(TEST_CASE)


# default benchmark parameters

NX    = 67108864
NWARM = 2
NIT   = 10

PARAMETERS = ${NX} ${NWARM} ${NIT}


# all

targets = $(BUILD_DIR)/$(TEST_CASE)-base $(BUILD_DIR)/$(TEST_CASE)-omp-host $(BUILD_DIR)/$(TEST_CASE)-cuda


.PHONY: all
all: mk-target-dir $(targets)


mk-target-dir:
	mkdir -p $(BUILD_DIR)


# build rules

$(BUILD_DIR)/$(TEST_CASE)-base: $(TEST_CASE)-base.cpp $(TEST_CASE)-util.h ../util.h
	g++ -O3 -march=native -std=c++17 -o $(BUILD_DIR)/$(TEST_CASE)-base $(TEST_CASE)-base.cpp

$(BUILD_DIR)/$(TEST_CASE)-omp-host: $(TEST_CASE)-omp-host.cpp $(TEST_CASE)-util.h ../util.h
	g++ -O3 -march=native -std=c++17 -fopenmp -o $(BUILD_DIR)/$(TEST_CASE)-omp-host $(TEST_CASE)-omp-host.cpp

$(BUILD_DIR)/$(TEST_CASE)-cuda: $(TEST_CASE)-cuda.cu $(TEST_CASE)-util.h ../util.h
	nvcc -O3 -arch=sm_89 -o $(BUILD_DIR)/$(TEST_CASE)-cuda $(TEST_CASE)-cuda.cu

# aliases without build directory

.PHONY: $(TEST_CASE)-base
$(TEST_CASE)-base: $(BUILD_DIR)/$(TEST_CASE)-base

.PHONY: $(TEST_CASE)-omp-host
$(TEST_CASE)-omp-host: $(BUILD_DIR)/$(TEST_CASE)-omp-host

.PHONY: $(TEST_CASE)-cuda
$(TEST_CASE)-cuda: $(BUILD_DIR)/$(TEST_CASE)-cuda



# automated benchmark target

.PHONY: bench
bench: all
	@echo "Base:"
	$(BUILD_DIR)/$(TEST_CASE)-base $(PARAMETERS)
	@echo ""
	@echo "OpenMP Host:"
	$(BUILD_DIR)/$(TEST_CASE)-omp-host $(PARAMETERS)
	@echo ""
	@echo "CUDA Host:"
	$(BUILD_DIR)/$(TEST_CASE)-cuda $(PARAMETERS)
	@echo ""


# clean target

.PHONY: clean
clean:
	rm $(targets)
