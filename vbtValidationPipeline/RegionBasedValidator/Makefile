CC := g++
# CC := clang --analyze # and comment out the linker last line for sanity

BUILDDIR := build
TARGET := vbtvalidator


COMMON := common
VBTVALIDATE := VBTvalidate


CFLAGS := -std=c++11 -Wall -O2 -g
LIB := -lz -pthread -lhts
INC := -I htslib -I $(COMMON) -I $(VBTVALIDATE)
 
SOURCESCOMMON := $(shell find $(COMMON) -type f -name '*.cpp')
SOURCESVBTVALIDATE := $(shell find $(VBTVALIDATE) -type f -name '*.cpp')

OBJECTSCOMMON := $(subst $(COMMON), $(BUILDDIR), $(SOURCESCOMMON:.cpp=.o))
OBJECTSVBTVALIDATE := $(subst $(VBTVALIDATE), $(BUILDDIR), $(SOURCESVBTVALIDATE:.cpp=.o))

OBJECTS := $(OBJECTSCOMMON) $(OBJECTSVBTVALIDATE) $(BUILDDIR)/main.o

all: $(TARGET)
	@echo "SUCCESSFULLY COMPILED!!"
	

$(TARGET): $(OBJECTS)
	@echo " Linking..."
	@echo " OBJECTS COMMON: $(OBJECTSCOMMON)"
	@echo " OBJECTS VBTVALIDATE: $(OBJECTSVBTVALIDATE)"
	@echo " $(CC) $^ -o $(TARGET) $(LIB)"; $(CC) $^ -o $(TARGET) $(LIB) $(INC)

$(BUILDDIR)/%.o: $(COMMON)/%.cpp
	@mkdir -p $(BUILDDIR)
	@echo " CORE: $(CC) $(CFLAGS) $(INC) -c -o $@ $<"; $(CC) $(CFLAGS) $(INC) -c $< -o $@

$(BUILDDIR)/%.o: $(VBTVALIDATE)/%.cpp
	@mkdir -p $(BUILDDIR)
	@echo " DUO: $(CC) $(CFLAGS) $(INC) -c -o $@ $<"; $(CC) $(CFLAGS) $(INC) -c $< -o $@

$(BUILDDIR)/main.o: main.cpp
	@mkdir -p $(BUILDDIR)
	@echo " MAIN: $(CC) $(CFLAGS) $(INC) -c -o $@ $<"; $(CC) $(CFLAGS) $(INC) -c $< -o $@
 
clean:
	@echo " Cleaning..."; 
	@echo " $(RM) -r $(BUILDDIR) $(TARGET)"; $(RM) -r $(BUILDDIR) $(TARGET)


.PHONY: clean
