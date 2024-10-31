TARGET = WLat
# compilers
CXX = g++

#compiling flags
CPPFLAGS := $(STDINCDIR) -I/N/u/smithwya/BigRed200/install/eigen
CPPFLAGS += -O3
CPPFLAGS += -mtune=native
CPPFLAGS += -march=native

SRCDIR = src
OBJDIR = obj
BINDIR = bin

SOURCES := $(wildcard $(SRCDIR)/*.cpp)
INCLUDES := $(wildcard $(SRCDIR)/*.h)
OBJECTS := $(SOURCES:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)
rm = rm -f

all: $(BINDIR)/$(TARGET)

$(BINDIR)/$(TARGET): $(OBJECTS)
	@if [ ! -d $(BINDIR) ]; then mkdir $(BINDIR); fi
	$(CXX) $(CPPFLAGS) -o $(BINDIR)/$(TARGET) $(OBJECTS)
	@echo "Linking complete"

$(OBJECTS): $(OBJDIR)/%.o : $(SRCDIR)/%.cpp
	@if [ ! -d $(OBJDIR) ]; then mkdir $(OBJDIR); fi
	$(CXX) $(CPPFLAGS) -c $< -o $@
	@echo "Compiled "$<" successfully"
clean:
	$(rm) $(OBJECTS)
	@echo "Cleanup complete"
