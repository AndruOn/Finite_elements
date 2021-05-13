CC=gcc
CFLAGS=-Werror -g -lm 
INCLUDE_HEADERS_DIRECTORY=-Iheaders
LIBS=-lm

proj: src/main.c src/motor.o src/fem.o
	$(CC) $(INCLUDE_HEADERS_DIRECTORY) $(CFLAGS) -o $@ $^ $(LIBS) 

src/%.o: src/%.c                  
	$(CC) $(INCLUDE_HEADERS_DIRECTORY) $(CFLAGS) -o $@ -c $< 

clean :
	rm -f src/*.o
	rm -f *.o
.PHONY: clean 