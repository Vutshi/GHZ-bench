CFLAGS := -march=native -Ofast
LIBS := -lm

ifeq ($(shell arch), e2k)
CFLAGS += -ffast
endif

NAME := ghz-bench
SOURCES := ghz-bench.c
OBJECTS := $(SOURCES:.c=.o)

.PHONY : all clean

all: $(NAME)

clean:
	$(RM) $(NAME) $(OBJECTS)

$(NAME): $(OBJECTS)
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)
