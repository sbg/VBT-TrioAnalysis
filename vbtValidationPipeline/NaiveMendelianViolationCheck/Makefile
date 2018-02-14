CC := g++
# CC := clang --analyze # and comment out the linker last line for sanity

BUILDDIR := build
TARGET := naive-mendelian


COMMON := common
VBTNAIVE := naive


CFLAGS := -std=c++11 -Wall -O2 -g
LIB := -lz -pthread -lhts
INC := -I htslib -I $(COMMON) -I $(VBTNAIVE) -I $(shell pwd)
 
SOURCESCOMMON := $(shell find $(COMMON) -type f -name '*.cpp')
SOURCESVBTNAIVE := $(shell find $(VBTNAIVE) -type f -name '*.cpp')

OBJECTSCOMMON := $(subst $(COMMON), $(BUILDDIR), $(SOURCESCOMMON:.cpp=.o))
OBJECTSVBTNAIVE := $(subst $(VBTNAIVE), $(BUILDDIR), $(SOURCESVBTNAIVE:.cpp=.o))

OBJECTS := $(OBJECTSCOMMON) $(OBJECTSVBTNAIVE) $(BUILDDIR)/main.o

all: $(TARGET)
	@echo "SUCCESSFULLY COMPILED!!"
	

$(TARGET): $(OBJECTS)
	@echo " Linking..."
	@echo " OBJECTS COMMON: $(OBJECTSCOMMON)"
	@echo " OBJECTS VBT NAIVE: $(OBJECTSVBTNAIVE)"
	@echo " $(CC) $^ -o $(TARGET) $(LIB)"; $(CC) $^ -o $(TARGET) $(LIB) $(INC)

$(BUILDDIR)/%.o: $(COMMON)/%.cpp
	@mkdir -p $(BUILDDIR)
	@echo " COMMON: $(CC) $(CFLAGS) $(INC) -c -o $@ $<"; $(CC) $(CFLAGS) $(INC) -c $< -o $@

$(BUILDDIR)/%.o: $(VBTNAIVE)/%.cpp
	@mkdir -p $(BUILDDIR)
	@echo " NAIVE: $(CC) $(CFLAGS) $(INC) -c -o $@ $<"; $(CC) $(CFLAGS) $(INC) -c $< -o $@

$(BUILDDIR)/main.o: main.cpp
	@mkdir -p $(BUILDDIR)
	@echo " MAIN: $(CC) $(CFLAGS) $(INC) -c -o $@ $<"; $(CC) $(CFLAGS) $(INC) -c $< -o $@
 
clean:
	@echo " Cleaning..."; 
	@echo " $(RM) -r $(BUILDDIR) $(TARGET)"; $(RM) -r $(BUILDDIR) $(TARGET)


.PHONY: clean
