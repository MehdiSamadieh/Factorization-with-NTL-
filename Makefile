# Use: ./sup <Degree>
CILKPP  = g++
INCLUDE_ARGS=-I$(NTL_HOME)/include -L $(NTL_HOME)/src
LIBARG  = -lbpas -lm -lntl -lgmp -lgmpxx  -lcilkrts -lmodpnLINUXINTEL64  -lmps -lpthread 
COMPILEARG = -c -O2 -g -fcilkplus -DLINUXINTEL64=1 -march=native
TARGET  = sup

all: $(TARGET)

$(TARGET): test.o
	$(CILKPP) -o $@ $^ $(LIBARG)

test.o: test.cpp $(BPAS_HOME)/include/Polynomial/upolynomial.h
	 $(CILKPP) $< $(LIBARG) $(COMPILEARG) $(INCLUDE_ARGS)

debug: COMPILEARG += -DBPASDEBUG=1
debug: $(TARGET)

serial: COMPILEARG += -DSERIAL=1
serial: $(TARGET)

test:
	./$(TARGET) 3	# SUP operations

clean:
	rm -rf $(TARGET) *.out* *~ *.log *.o *.dat
