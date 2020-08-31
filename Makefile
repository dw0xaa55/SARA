OBJS=src/sara.c
CC=gcc
LFLAGS=-Wno-format -Wno-deprecated-declarations -Wno-format-security -lm `pkg-config --cflags --libs gtk+-3.0` -export-dynamic
CFLAGS=-Wall -g
TARGET=bin/sara

$(TARGET):$(OBJS)
	$(CC) $(OBJS) $(LFLAGS) $(CFLAGS) -o $(TARGET)
	cp src/sara.glade bin/
clean:
	$(RM) $(TARGET)
	$(RM) bin/sara.glade

install:
	cp bin/sara /usr/bin/
	mkdir /usr/share/sara/
	cp bin/sara.glade /usr/share/sara/sara.glade

uninstall:
	rm -f /usr/bin/sara
	rm -rf /usr/share/sara/
