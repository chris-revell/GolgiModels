# Compile model with clang++ compiler

# REMEMBER:
# export CPATH=/opt/homebrew/include
# export LIBRARY_PATH=/opt/homebrew/lib

CC := /opt/homebrew/Cellar/llvm/13.0.0_1/bin/clang++
SRCDIR := src
SRCEXT := cpp
SOURCES := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
BUILDDIR := build
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))
TARGET := program
CFLAGS := -c -Wall -std=c++11 -g -O0 #Uncomment this line to compile for debugging
#CFLAGS := -c -O3 -std=c++11	#Uncomment this line to compile for fast runs
LIB :=
INC := -I include

$(TARGET): $(OBJECTS)
	@echo " Linking..."
	@echo " $(CC) $^ -o $(TARGET) $(LIB)"; $(CC) $^ -o $(TARGET) $(LIB)

$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(BUILDDIR)
	@echo " $(CC) $(CFLAGS) $(INC) -c -o $@ $<"; $(CC) $(CFLAGS) $(INC) -c -o $@ $<

clean:
	@echo " Cleaning...";
	@echo " $(RM) -r $(BUILDDIR) $(TARGET)"; $(RM) -r $(BUILDDIR) $(TARGET)
