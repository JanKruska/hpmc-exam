CC         := gcc
LINKER     := $(CC)

# ---------------------

CFLAGS     := -g -O0 -m64 -I${MKLROOT}/include -m64 -mavx2 -march=native -fopenmp
FFLAGS     := $(CFLAGS) 
LDFLAGS    := -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl

all: main.out dgemm.out

# ---------------------

main.out: main.o utils.o
	$(LINKER) $(CFLAGS) -o $@ $^ $(LDFLAGS)


dgemm.out: dgemm.o utils.o
	$(LINKER) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# ---------------------

clean:
	rm -f *.o *~ core *.out *.pdf
