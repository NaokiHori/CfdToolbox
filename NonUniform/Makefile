CC     := cc
CFLAG  := -std=c99 -Wall -Wextra -Wpedantic -O3#-fopenmp
INC    := -Iinclude
LIB    := -lm
SRCDIR := src
OBJDIR := obj
SRCS   := $(shell find $(SRCDIR) -type f -name "*.c")
OBJS   := $(patsubst %.c, $(OBJDIR)/%.o, $(SRCS))
DEPS   := $(patsubst %.c, $(OBJDIR)/%.d, $(SRCS))
TARGET := a.out
OUTDIR := output

help:
	@echo "all     : create \"$(TARGET)\""
	@echo "clean   : remove \"$(TARGET)\" and object files under \"$(OBJDIR)\""
	@echo "output  : create \"$(OUTDIR)\""
	@echo "datadel : clean-up \"$(OUTDIR)\""
	@echo "help    : show this message"

all: $(TARGET)

clean:
	$(RM) -r $(OBJDIR) $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAG) -o $@ $^ $(LIB)

$(OBJDIR)/%.o: %.c
	@if [ ! -e $(dir $@) ]; then \
		mkdir -p $(dir $@); \
	fi
	$(CC) $(CFLAG) -MMD $(INC) -c $< -o $@

output:
	@if [ ! -e $(OUTDIR)/log ]; then \
		mkdir -p $(OUTDIR)/log; \
	fi
	@if [ ! -e $(OUTDIR)/save ]; then \
		mkdir -p $(OUTDIR)/save; \
	fi

datadel:
	@rm -rf $(OUTDIR)
	@make output

-include $(DEPS)

.PHONY : help all clean output

