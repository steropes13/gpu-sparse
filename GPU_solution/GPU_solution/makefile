CC = nvcc
CUDA_TOOLKIT := $(shell dirname $$(command -v nvcc))/..
LIBS =
INCLUDES = -I../../ -I$(CUDA_TOOLKIT)/include
LIB_FLAGS = -lm -lcusparse

BIN_FOLDER := bin
OBJ_FOLDER := obj
SRC_FOLDER := src
BATCH_OUT_FOLDER := outputs

MAIN_NAME = spmv
MAIN_BIN = $(MAIN_NAME)
MAIN_SRC = $(MAIN_NAME).cu

OBJECTS = $(OBJ_FOLDER)/my_time_lib.o \
          $(OBJ_FOLDER)/mmio.o \
          $(OBJ_FOLDER)/spmvCPU.o 

OBJECTS_CU = $(OBJ_FOLDER)/cuSparseComputation.o

all: $(BIN_FOLDER)/$(MAIN_BIN)

$(OBJ_FOLDER)/%.o: $(SRC_FOLDER)/%.c
	@mkdir -p $(BIN_FOLDER) $(OBJ_FOLDER) $(BATCH_OUT_FOLDER)
	$(CC) -c $< -o $@ $(LIB_FLAGS)

$(OBJ_FOLDER)/%.o: $(SRC_FOLDER)/%.cu
	@mkdir -p $(BIN_FOLDER) $(OBJ_FOLDER) $(BATCH_OUT_FOLDER)
	$(CC) -c $< -o $@ $(LIB_FLAGS) $(INCLUDES) 


$(BIN_FOLDER)/$(MAIN_BIN): $(MAIN_SRC) $(OBJECTS) $(OBJECTS_CU)
	@mkdir -p $(BIN_FOLDER)
	$(CC) -arch=sm_80 $^ -g -o $@ $(LIBS) $(INCLUDES) $(LIB_FLAGS)

clean:
	rm -rf $(BIN_FOLDER) $(OBJ_FOLDER) 

#we decided here to add the -arch=sm_80 to specify the architecture of the GPU card (NVIDIA A30, compute capability 8.0)). it allows the activation of specific functionnalities of the material such as atomicAdd with (double*, double, see the function computeSpmvCOOGPU) that are not available in older architecture. 
