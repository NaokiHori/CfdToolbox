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

help:
	@echo "all   : create \"$(TARGET)\""
	@echo "clean : remove \"$(TARGET)\" and object files under \"$(OBJDIR)\""
	@echo "help  : show this message"

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

test_derivative:
	$(CC) -DTEST_DERIVATIVE $(CFLAG) $(INC) \
		$(SRCDIR)/constant.c \
		$(SRCDIR)/derivative.c \
		$(SRCDIR)/domain.c \
		-o $(TARGET)

test_domain:
	$(CC) -DTEST_DOMAIN $(CFLAG) $(INC) \
		$(SRCDIR)/constant.c \
		$(SRCDIR)/domain.c \
		-o $(TARGET)

test_runge_kutta:
	$(CC) -DTEST_RUNGE_KUTTA $(CFLAG) $(INC) \
		$(SRCDIR)/runge_kutta.c \
		-o $(TARGET)

test_transform:
	$(CC) -DTEST_TRANSFORM $(CFLAG) $(INC) \
		$(SRCDIR)/constant.c \
		$(SRCDIR)/transform.c \
		$(SRCDIR)/transform/dft.c \
		$(SRCDIR)/transform/rdft.c \
		-o $(TARGET)

-include $(DEPS)

.PHONY : help all clean

