CC := g++
# CC := clang --analyze # and comment out the linker last line for sanity

BUILDDIR := build
TARGET := vbt


INCCORE := Core/include
INCDUO := DuoComparison/include
INCTRIO := MendelianViolation/include
INCVCFIO := VcfIO/include
INCGRAPH := GraphComparison/include

CFLAGS := -std=c++11 -Wall -O2 -g
LIB := -lz -pthread -lhts
INC := -I $(INCCORE) -I htslib -I $(INCDUO) -I $(INCTRIO) -I $(INCVCFIO) -I $(INCGRAPH) -I $(shell pwd)

SRCCORE := Core/src
SRCDUO := DuoComparison/src
SRCTRIO := MendelianViolation/src
SRCVCFIO := VcfIO/src
SRCGRAPH := GraphComparison/src

 
SOURCESCORE := $(shell find $(SRCCORE) -type f -name '*.cpp')
SOURCESDUO := $(shell find $(SRCDUO) -type f -name '*.cpp')
SOURCESTRIO := $(shell find $(SRCTRIO) -type f -name '*.cpp')
SOURCESVCFIO := $(shell find $(SRCVCFIO) -type f -name '*.cpp')
SOURCESGRAPH := $(shell find $(SRCGRAPH) -type f -name '*.cpp')

OBJECTSCORE := $(subst $(SRCCORE), $(BUILDDIR), $(SOURCESCORE:.cpp=.o))
OBJECTSDUO := $(subst $(SRCDUO), $(BUILDDIR), $(SOURCESDUO:.cpp=.o))
OBJECTSTRIO := $(subst $(SRCTRIO), $(BUILDDIR), $(SOURCESTRIO:.cpp=.o))
OBJECTSVCFIO := $(subst $(SRCVCFIO), $(BUILDDIR), $(SOURCESVCFIO:.cpp=.o))
OBJECTSGRAPH := $(subst $(SRCGRAPH), $(BUILDDIR), $(SOURCESGRAPH:.cpp=.o))

OBJECTS := $(OBJECTSCORE) $(OBJECTSDUO) $(OBJECTSTRIO) $(OBJECTSVCFIO) $(OBJECTSGRAPH) $(BUILDDIR)/main.o
SOURCES := $(SOURCESCORE) $(SOURCESDUO) $(SOURCESTRIO) $(SOURCESVCFIO) $(SOURCESGRAPH) main.cpp

all: $(TARGET)
	@echo "SUCCESSFULLY COMPILED!!"
	

$(TARGET): $(OBJECTS)
	@echo " Linking..."
	@echo " OBJECTS CORE: $(OBJECTSCORE)"
	@echo " OBJECTS DUO : $(OBJECTSDUO)"
	@echo " OBJECTS TRIO: $(OBJECTSTRIO)"
	@echo " $(CC) $^ -o $(TARGET) $(LIB)"; $(CC) $^ -o $(TARGET) $(LIB)

$(BUILDDIR)/%.o: $(SRCCORE)/%.cpp Constants.h
	@mkdir -p $(BUILDDIR)
	@echo " CORE: $(CC) $(CFLAGS) $(INC) -c -o $@ $<"; $(CC) $(CFLAGS) $(INC) -c $< -o $@

$(BUILDDIR)/%.o: $(SRCDUO)/%.cpp Constants.h
	@mkdir -p $(BUILDDIR)
	@echo " DUO: $(CC) $(CFLAGS) $(INC) -c -o $@ $<"; $(CC) $(CFLAGS) $(INC) -c $< -o $@

$(BUILDDIR)/%.o: $(SRCTRIO)/%.cpp Constants.h
	@mkdir -p $(BUILDDIR)
	@echo " TRIO: $(CC) $(CFLAGS) $(INC) -c -o $@ $<"; $(CC) $(CFLAGS) $(INC) -c $< -o $@

$(BUILDDIR)/%.o: $(SRCVCFIO)/%.cpp Constants.h
	@mkdir -p $(BUILDDIR)
	@echo " VCFIO: $(CC) $(CFLAGS) $(INC) -c -o $@ $<"; $(CC) $(CFLAGS) $(INC) -c $< -o $@

$(BUILDDIR)/%.o: $(SRCGRAPH)/%.cpp Constants.h
    @mkdir -p $(BUILDDIR)
    @echo " VCFIO: $(CC) $(CFLAGS) $(INC) -c -o $@ $<"; $(CC) $(CFLAGS) $(INC) -c $< -o $@

$(BUILDDIR)/main.o: main.cpp
	@mkdir -p $(BUILDDIR)
	@echo " MAIN: $(CC) $(CFLAGS) $(INC) -c -o $@ $<"; $(CC) $(CFLAGS) $(INC) -c $< -o $@
 
clean:
	@echo " Cleaning..."; 
	@echo " $(RM) -r $(BUILDDIR) $(TARGET)"; $(RM) -r $(BUILDDIR) $(TARGET)


.PHONY: clean
