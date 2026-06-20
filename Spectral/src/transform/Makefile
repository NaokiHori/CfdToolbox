CC     := cc
CFLAG  := -std=c99 -Wall -Wextra -Wpedantic -O3
TARGET := a.out

test_dft:
	$(CC) -DTEST_DFT $(CFLAG) \
		dft.c \
		test_dft.c \
		-o $(TARGET)

test_rdft:
	$(CC) -DTEST_RDFT $(CFLAG) \
		rdft.c \
		test_rdft.c \
		-o $(TARGET)
