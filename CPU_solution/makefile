CC=gcc

LIBS=
INCLUDES=-I../../
LIB_FLAGS=-lm


BIN_FOLDER := bin
OBJ_FOLDER := obj
SRC_FOLDER := src
BATCH_OUT_FOLDER := outputs


MAIN_NAME=sequentialSpMV
MAIN_BIN=$(MAIN_NAME)
MAIN_SRC=$(MAIN_NAME).c

OBJECTS = $(OBJ_FOLDER)/my_time_lib.o
OBJECTS_MMIO = $(OBJ_FOLDER)/mmio.o


all: $(BIN_FOLDER)/$(MAIN_BIN)

$(OBJ_FOLDER)/mmio.o: $(SRC_FOLDER)/mmio.c
	@mkdir -p $(BIN_FOLDER) $(OBJ_FOLDER) $(BATCH_OUT_FOLDER)
	$(CC) -c $(SRC_FOLDER)/mmio.c -o $@ $(LIB_FLAGS)

$(OBJ_FOLDER)/my_time_lib.o: $(SRC_FOLDER)/my_time_lib.c
	@mkdir -p $(BIN_FOLDER) $(OBJ_FOLDER) $(BATCH_OUT_FOLDER)
	$(CC) -c $(SRC_FOLDER)/my_time_lib.c -o $@ $(LIB_FLAGS)

$(BIN_FOLDER)/$(MAIN_BIN): $(MAIN_SRC) $(OBJECTS) $(OBJECTS_MMIO)
	mkdir -p $(BIN_FOLDER)
	$(CC) $^ -g -o $@ $(LIBS) $(INCLUDES) $(LIB_FLAGS)

clean:
	rm -rf $(BIN_FOLDER) $(OBJ_FOLDER)
