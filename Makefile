SRC_DIR = src/
INC_DIR = include/
EXE_BIN_DIR = bin/release/
OBJ_BIN_DIR = bin/objs/
SRC_MOTIONDETECT_DIR := $(SRC_DIR)libmotiondetection/
EXE_NAME = autraco
JPEG_INC_DIR 	= -Idependencies/jpeg/include	
PNG_INC_DIR		=	-Idependencies/png/include
ZLIB_INC_DIR 	= -Idependencies/zlib/include
JPZ_INC_DIR		=  $(JPEG_INC_DIR) $(PNG_INC_DIR) $(ZLIB_INC_DIR)

JPEG_LIB_DIR 	= -Ldependencies/jpeg/lib 
PNG_LIB_DIR		=	-Ldependencies/png/lib
ZLIB_LIB_DIR	= -Ldependencies/zlib/lib	
JPZ_LIB_DIR 	= $(JPEG_LIB_DIR) $(PNG_LIB_DIR) $(ZLIB_LIB_DIR)
JPZ_LIBS 			= -ljpeg -lpng -lz

export PKG_CONFIG_PATH = dependencies/ffmpeg/lib/pkgconfig
#$(info $(PKG_CONFIG_PATH))
FFMPEG_LIBS=    libavdevice                        \
                libavformat                        \
                libavfilter                        \
                libavcodec                         \
                libswresample                      \
                libswscale                         \
                libavutil                          \
#$(info $(SRC_MOTIONDETECT_DIR))
CC = gcc 
CFLAGS = -Wall -O3 -c -Iinclude $(JPZ_INC_DIR) -I$(SRC_MOTIONDETECT_DIR)
CFLAGS += $(shell pkg-config --cflags $(FFMPEG_LIBS)) 
LFLAGS := $(JPZ_LIB_DIR) $(shell pkg-config --libs $(FFMPEG_LIBS)) $(JPZ_LIBS)

OBJS = $(OBJ_BIN_DIR)main.o $(OBJ_BIN_DIR)image_reader_writer.o \
 				$(OBJ_BIN_DIR)jpeghandler.o $(OBJ_BIN_DIR)pnghandler.o \
				$(OBJ_BIN_DIR)ffmpeghandler.o $(OBJ_BIN_DIR)affinetransforms.o\
				$(OBJ_BIN_DIR)detect.o

all: $(EXE_BIN_DIR)$(EXE_NAME)

$(EXE_BIN_DIR)$(EXE_NAME): $(OBJS)
	$(CC)  -o $(EXE_BIN_DIR)$(EXE_NAME) $(OBJS) $(LFLAGS)

$(OBJ_BIN_DIR)main.o: $(SRC_DIR)main.c $(INC_DIR)main.h
	$(CC) $(CFLAGS)  -o $(OBJ_BIN_DIR)main.o  $(SRC_DIR)main.c

$(OBJ_BIN_DIR)image_reader_writer.o: $(SRC_DIR)image_reader_writer.c $(INC_DIR)image_reader_writer.h
	$(CC)  $(CFLAGS) -o $(OBJ_BIN_DIR)image_reader_writer.o $(SRC_DIR)image_reader_writer.c

$(OBJ_BIN_DIR)jpeghandler.o: $(SRC_DIR)jpeghandler.c $(INC_DIR)jpeghandler.h
	$(CC) $(CFLAGS) -o $(OBJ_BIN_DIR)jpeghandler.o $(SRC_DIR)jpeghandler.c

$(OBJ_BIN_DIR)pnghandler.o: $(SRC_DIR)pnghandler.c $(INC_DIR)pnghandler.h
	$(CC)  $(CFLAGS) -o $(OBJ_BIN_DIR)pnghandler.o $(SRC_DIR)pnghandler.c

$(OBJ_BIN_DIR)ffmpeghandler.o: $(SRC_DIR)ffmpeghandler.c $(INC_DIR)ffmpeghandler.h
	$(CC) $(CFLAGS) -o $(OBJ_BIN_DIR)ffmpeghandler.o $(SRC_DIR)ffmpeghandler.c

$(OBJ_BIN_DIR)affinetransforms.o: $(SRC_DIR)affinetransforms.c $(INC_DIR)affinetransforms.h
	$(CC) $(CFLAGS) -o	$(OBJ_BIN_DIR)affinetransforms.o $(SRC_DIR)affinetransforms.c

$(OBJ_BIN_DIR)detect.o: $(SRC_MOTIONDETECT_DIR)detect.c $(SRC_MOTIONDETECT_DIR)detect.h
	$(CC) $(CFLAGS) -o $(OBJ_BIN_DIR)detect.o $(SRC_MOTIONDETECT_DIR)detect.c

clean_objs:
	rm -f $(OBJ_BIN_DIR)*.o

clean_exes:
	rm -f $(EXE_BIN_DIR)*.exe
