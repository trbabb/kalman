CC = clang++
LD = clang++
AR = ar

# compile
INCLUDES = include 
CFLAGS   = -g -std=c++11 -Os -Wall -c -fmessage-length=0 -Wno-unused -Wno-unused-local-typedefs
IFLAGS   = $(addprefix -I, $(INCLUDES))

# link
LIBDIRS  =  
LIBS     = geomc 
LDFLAGS  = $(addprefix -l, $(LIBS)) \
           $(addprefix -L, $(LIBDIRS)) \

# sources
MODULES  = imu
SRC      = $(wildcard src/*.cpp) \
           $(foreach m, $(MODULES), $(wildcard src/$(m)/*.cpp))
OBJ      = $(patsubst src/%.cpp, build/%.o, $(SRC))


all: test

clean:
	rm -rf ./build/*

## binaries


test: build/test.o bin
	$(CC) $(LDFLAGS) build/test.o -o bin/test

build/%.o : src/%.cpp build
	mkdir -p $(patsubst src/%, build/%, $(dir $<))
	$(CC) $(CFLAGS) $(IFLAGS) -o $@ $<

build :
	mkdir build

bin :
	mkdir bin

