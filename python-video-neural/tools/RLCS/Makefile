SRCDIR = ./src
OBJDIR = ./obj
LIBDIR = ./lib
INCDIR = ./include

PREFIX = /usr/local

SHELL  =  /bin/bash
CC     =  gcc
CFLAGS  =  -c -Wall -Werror -fPIC
LDFLAGS = -shared

all: $(LIBDIR)/librlcs.so

$(LIBDIR)/librlcs.so: $(OBJDIR)/types.o $(OBJDIR)/rlcs.o $(OBJDIR)/maxmin.o $(OBJDIR)/sp.o $(OBJDIR)/yamlwrapper.o
	$(CC) $(LDFLAGS) -o $(LIBDIR)/librlcs.so $(OBJDIR)/types.o $(OBJDIR)/rlcs.o $(OBJDIR)/maxmin.o $(OBJDIR)/sp.o $(OBJDIR)/yamlwrapper.o


#$(LIBDIR)/librlcs.so: $(OBJDIR)/types.o $(OBJDIR)/maxmin.o $(OBJDIR)/sp.o
#	$(CC) $(LDFLAGS) -o $(LIBDIR)/librlcs.so $(OBJDIR)/types.o $(OBJDIR)/maxmin.o $(OBJDIR)/sp.o


$(OBJDIR)/yamlwrapper.o: $(SRCDIR)/yamlwrapper.c
	$(CC) $(CFLAGS) -I $(INCDIR) -o $(OBJDIR)/yamlwrapper.o $(SRCDIR)/yamlwrapper.c 


$(OBJDIR)/types.o: $(SRCDIR)/types.c
	$(CC) $(CFLAGS) -I $(INCDIR) -o $(OBJDIR)/types.o $(SRCDIR)/types.c 

$(OBJDIR)/maxmin.o: $(SRCDIR)/maxmin.c
	$(CC) $(CFLAGS) -I $(INCDIR) -o $(OBJDIR)/maxmin.o $(SRCDIR)/maxmin.c 

$(OBJDIR)/rlcs.o: $(SRCDIR)/rlcs.c
	$(CC) $(CFLAGS) -I $(INCDIR) -o $(OBJDIR)/rlcs.o $(SRCDIR)/rlcs.c 

$(OBJDIR)/sp.o: $(SRCDIR)/sp.c
	$(CC) $(CFLAGS) -I $(INCDIR) -o $(OBJDIR)/sp.o $(SRCDIR)/sp.c

install: 
	cp -v $(LIBDIR)/librlcs.so $(PREFIX)/lib/
	cp -v $(INCDIR)/types.h $(INCDIR)/maxmin.h $(INCDIR)/rlcs.h   $(INCDIR)/sp.h $(PREFIX)/include/

remove:
	rm -rf $(PREFIX)/lib/librlcs.so
	rm -rf $(PREFIX)/include/types.h  $(PREFIX)/include/maxmin.h  $(PREFIX)/include/rlcs.h $(PREFIX)/include/sp.h 

clean:
	rm -rf  $(OBJDIR)/*.o  $(LIBDIR)/*.so
