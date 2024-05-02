# configuration

TEST_CASE = julia_set
BUILD_DIR = build
LODEPNG_OBJ = lodepng.o
SHOW_IMAGE = plot_image

# all

targets = $(BUILD_DIR)/$(TEST_CASE)


.PHONY: all
all: mk-target-dir $(targets)


mk-target-dir:
	mkdir -p $(BUILD_DIR)


# build rules

$(BUILD_DIR)/$(TEST_CASE): $(TEST_CASE).cu
	g++ -O3 -c  lodepng.cpp -o $(BUILD_DIR)/$(LODEPNG_OBJ)
	nvcc -O3 -arch=sm_89 $(TEST_CASE).cu $(BUILD_DIR)/$(LODEPNG_OBJ) -o $(BUILD_DIR)/$(TEST_CASE)
	./$(BUILD_DIR)/$(TEST_CASE) 
	g++ -O3 plot_image.cpp -o $(BUILD_DIR)/$(SHOW_IMAGE) -Dcimg_use_png -lpng -lX11 -lz -lpthread
	./$(BUILD_DIR)/$(SHOW_IMAGE)

# aliases without build directory

.PHONY: $(TEST_CASE)
$(TEST_CASE): $(BUILD_DIR)/$(TEST_CASE)

# clean target

.PHONY: clean
clean:
	rm $(targets)
