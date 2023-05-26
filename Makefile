CC = g++

CFLAGS = 

TARGET = fft_1d

SRCS = fft_1d.cpp

OBJS = ${SRCS:.cpp=.o}

INCDIR  =

LIBDIR  = 

LIBS    = 

$(TARGET): $(OBJS)
	$(CC) -o $@ $^ $(LIBDIR) $(LIBS) 
	./fft_1d

$(OBJS): $(SRCS) 
	$(CC) $(CFLAGS) $(INCDIR) -c $(SRCS) 

all: clean $(OBJS) $(TARGET)
	
clean: 
	-rm -f $(OBJS) $(TARGET) *.d 
