CC=clang

LIBS=
INCLUDES=
LIB_FLAGS=-lm


SRC_FOLDER := src
BIN_FOLDER := bin
OBJ_FOLDER := obj
BATCH_OUT_FOLDER := outputs


MAIN_NAME=example_read
MAIN_BIN=$(MAIN_NAME)
MAIN_SRC=$(MAIN_NAME).c

OBJECTS = $(OBJ_FOLDER)/mmio.o

all: $(BIN_FOLDER)/$(MAIN_BIN)

$(OBJ_FOLDER)/mmio.o: $(SRC_FOLDER)/mmio.c
	@mkdir -p $(BIN_FOLDER) $(OBJ_FOLDER) $(BATCH_OUT_FOLDER)
	$(CC) -c $(SRC_FOLDER)/mmio.c -o $@ $(LIB_FLAGS)

$(BIN_FOLDER)/$(MAIN_BIN): $(MAIN_SRC) $(OBJECTS)
	mkdir -p $(BIN_FOLDER)
	$(CC) $^ -o $@ $(LIBS) $(INCLUDES) $(LIB_FLAGS)


clean:
	rm -rf $(BIN_FOLDER) $(OBJ_FOLDER)
