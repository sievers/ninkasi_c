CC=gcc
CPP=g++
BINDIR= ./
SRCDIR= ./
LIBDIR= ./

INCLUDES= -I./
#INCLUDES= -I./include/ -I../  -I/cita/d/raid-sievers/sievers/util/fftw-3.2/include/ -I /cita/d/raid-nolta/local/include/

CFLAGS= -fPIC -O3 -std=c99 ${INCLUDES}

#COBJ=dirfile.o getdata.o iers_bulletin_a.o  astro.o

#CPP_OBJ =  bit_constants.o bitstream.o raw_section.o slim_channel.o  slim_codec_constant.o  slim_codec_default.o  slim_codec_factory.o  slim_codec_runlength.o  slim_codec_reduced_binary.o  slim_control.o  slim_exceptions.o  slim_file.o  slimlib.o


#CPP_OBJ=act_spec.o bit_constants.o bitstream.o compute_entropy.o crc.o raw_section.o slim_acthk.o slim_channel.o slim_codec_constant.o slim_codec_default.o slim_codec_factory.o slim_codec_reduced_binary.o slim_codec_runlength.o slim_control_acthk.o slim_control.o  slim_exceptions.o slim_file.o slimlib.o test_bitstream.o test_building_blocks.o

#CPP_OBJ=act_spec.o bit_constants.o bitstream.o build_bit_constants_h.o compute_entropy.o crc.o raw_section.o slim_acthk.o slim_channel.o slim_codec_constant.o slim_codec_default.o slim_codec_factory.o slim_codec_reduced_binary.o slim_codec_runlength.o slim_control_acthk.o slim_control.o slim.o slim_dump.o slim_exceptions.o slim_file.o slimlib.o test_bitstream.o test_building_blocks.o


COBJ=dirfile.o getdata.o ${CPP_OBJ}
LINK= -shared  -lc 
VERSION= 1.0.1
#gcc -shared -Wl,-soname,libmystuff.so.1 \
#    -o libmystuff.so.1.0.1 a.o b.o -lc


%.o:%.c
	$(CC) $(CFLAGS)  -c $< -o $@
#	$(CC) $(CFLAGS)  -E $< -o ${<}.cpp


CPPFLAGS = -O3 -Wall -fPIC -I.
%.o:%.cpp
	$(CPP) $(CPPFLAGS)  -c $< -o $@


all: libreadtod


libreadtod: ${COBJ}
	${CC} ${LINK} -Wl,-soname,${@}.so.1 ${COBJ} -o ${@}.so.${VERSION}	
	rm -f ${@}.so
	ln -s ${@}.so.${VERSION} ${@}.so


distclean:
	rm -f lib*.so*
	rm -f *.o
clean:
	rm -f *.o 


#jgunzip: jgunzip.o
#        $(CC) $(FFLAGS) jgunzip.o -o $(BINDIR)$@ $(ZLIB) $(CLIBS) $(FLIBS)
