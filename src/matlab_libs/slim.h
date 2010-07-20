// -*- mode: c++; -*-

/// \file slim.h
/// Prototypes and class declarations for all slim objects.
/// Includes files, channels, and all encoders/decoders.

#ifndef SLIM_H
#define SLIM_H

#include <iostream>
#include <cstdlib>
#include <cassert>
#include "bitstream.h"
#include <string.h>


//#define DEBUG_ENCODING

using namespace std;

class slim_compressor_t;
class slim_expander_t;
class slim_channel;
class slim_channel_array;
class slim_channel_encode;
class slim_channel_decode;

class encoder;
class decoder;
class encoder_reduced_binary;
class decoder_reduced_binary;
class encoder_runlength;
class decoder_runlength;
class raw_section;

/// Allowed coder/decoder methods.
enum code_t {
  SLIM_ENCODER_DEFAULT,   ///< The default non-compressing encoder.
  SLIM_ENCODER_REDUCED_BINARY,///< The original reduced_binary encoder.
  SLIM_ENCODER_CODE_A,    ///< Variant A of REDUCED_BINARY (better params).
  SLIM_ENCODER_CODE_B,    ///< Variant B of REDUCED_BINARY (shorter overflows).
  SLIM_ENCODER_HUFFMAN,   ///< Uses Huffman coding for upper bits.
  SLIM_ENCODER_RUNLENGTH, ///< Uses (value, repeats) pairs.
  SLIM_ENCODER_CONSTANT,  ///< For strictly constant values.
};

/// Allowed data types.
enum data_t {
  SLIM_TYPE_UNDETERMINED, ///< This channel's type is not known yet.
  SLIM_TYPE_U32,          ///< Type uint32_t
  SLIM_TYPE_I32,          ///< Type int32_t
  SLIM_TYPE_U16,          ///< Type uint16_t
  SLIM_TYPE_I16,          ///< Type int16_t
  SLIM_TYPE_FLOAT,        ///< Type IEEE-754 float
  SLIM_TYPE_DOUBLE,       ///< Type IEEE-754 double-precition float
  SLIM_TYPE_U8,           ///< Type uint8_t
  SLIM_TYPE_I8,           ///< Type int8_t
};

static const size_t slim_type_size[] = {
  0, 4, 4, 2, 2, 4, 8, 1, 1}; ///< Sizes of the types in the enum data_t list.


/// Process state: are we encoding raw data or decoding a slim file?
enum slim_mode_t {
  SLIM_ENCODE,         ///< Encoding raw data.
  SLIM_DECODE,         ///< Decoding slim data.
  SLIM_MODE_UNKNOWN,   ///< Encode/decode mode not yet determined.
};

/// Sizes (in bits) of some header data.
enum {
  BITS_SLIM_SECT_SIZE = 32,  ///< Bits for size of raw section (in bytes).
  BITS_SLIM_NUM_CHAN = 24,   ///< Bits for number of channels in section.
  BITS_SLIM_REPETITIONS = 24,///< Bits for number of reps/channel.
  BITS_SLIM_NBITS = 5,
};

/// Ghost bytes are those padded onto the end of the last raw data section
/// in a file.  Their purpose is to ensure that we can encode the last word
/// as usual, even if the raw data happen to end in the middle of that word.
const size_t MAX_GHOST_BYTES=sizeof(double)-1;


/// File header flags.
enum file_header_flags_t {
  FLAG_SIZE =    0x01, ///< x01 Raw file size is present in file hdr.
  FLAG_NAME =    0x02, ///< x02 Raw file name is present in file hdr.
  FLAG_XTRA =    0x04, ///< x04 Unspecified extra data is present in file hdr.
  FLAG_TOC =     0x08, ///< x08 Section pointers are present in file hdr.
  FLAG_ONECHAN = 0x10, ///< x10 All sections have only one channel.
  FLAG_NOREPS =  0x20, ///< x20 No channels repeat within a frame.
  FLAG_CRC =     0x40, ///< x40 Raw data CRC-32 is present at end of sections.
};

#define FILE_MAGIC "SL"     ///< ASCII string at byte 0 of all slim files.
#define SLIM_SUFFIX "slm"   ///< Slim file suffix.


//---------------------------------------------------------------------------
// Exceptions
//---------------------------------------------------------------------------

class bad_file
{
public:
  bad_file(const char *filename, const char *message);
  bad_file(const bad_file &bf);
  ~bad_file();
  void mesg() const;

protected:
  char *str;
};


class bad_output_file : public bad_file
{
public:
  bad_output_file(const char *filename, const char *read_write_gerund);
  ~bad_output_file();
};


//---------------------------------------------------------------------------
// Inline 
//---------------------------------------------------------------------------

/// Integer division, rounding up any remainder
/// \param x     Number to be divided.
/// \param denom Denominator (divisor).
inline int divide_round_up(int x, int denom) {
  return 1+(x-1)/denom;
}


//---------------------------------------------------------------------------
// Command-line options.
//---------------------------------------------------------------------------

class slim_control {
public:
  slim_control();
  virtual ~slim_control() {;}   ///< Do-nothing destructor.
  void process_args(int argc, char * const argv[]);
  char flags() const;
  char flags(char in);
  void unslim();
  void slimcat();
  void handle_one_file(const char *fname);
  virtual void usage() const;
  virtual void version() const;
  void usage_printoptions() const;

protected:
  // Private methods
  void set_defaults();
  enum slim_mode_t detect_file_mode(const char *fname) const;
  virtual void compress_one_file(const char *fname);
  virtual void expand_one_file(const char *fname);
  void debug_compress_from_memory(const char *rawname);
  void debug_expand_from_memory(const char *rawname);

protected:
  // Private attributes
  bool deltas;           ///< Will all channels encode deltas?
  bool force_clobber;    ///< Willing to clobber or to slim files with multiple links?
  bool preserve_input;   ///< Do we save the input file (or delete it)?
  bool practice;         ///< Do we discard the compressed file?
  bool permit_bitrotation; ///< Do we have channels try to rotate up low bits?
  bool use_stdout;       ///< Write to stdout instead of to a file.
  code_t code_method;    ///< Code method for all channels.
  data_t data_type;      ///< Data type for all channels.
  int  nchan;            ///< How many channels are in the raw file?
  int  nframes;          ///< How many frames allowed per section?
  int  repeats;          ///< How many repeats per channel in a frame?
  int  sample_pct;       ///< What pct of data per channel to use when sampling?
  size_t debug_buf_size; ///< Buffer size to use in debugging write()/read()?
  enum slim_mode_t mode; ///< Is the SLIM_ENCODE or SLIM_DECODE mode?

  bool save_filename;    ///< Compressed file saves the filename.
  bool save_rawsize;     ///< Compressed file saves the raw file size.
  bool have_xtra;        ///< Compressed file has XTRA header data.
  bool have_toc;         ///< Compressed file has a Section Table of Contents.
  bool onechan;          ///< Compressed file has only 1 channel.
  bool noreps;           ///< Compressed file has only 1 datum per frame.
  bool crc;              ///< CRC-32 appears at the end of each section.
  bool ignore_crc;       ///< Do not test the CRC value on expansion.
  bool reserved0;        ///< Not used.
};


// A container class for an array of slim_channel objects.

class slim_channel_array {
public:
  slim_channel_array(int n_alloc=20);
  ~slim_channel_array();
  int offset(int i) const;
  slim_channel *operator[](int i) const;
  void push(slim_channel *c, size_t frame_offset);
  int size() const {return num_chan;} ///< Return the # of channels known.
  void clear();

private:
  int num_chan;              ///< Number of channels stored in object.
  int *offsets_in_frame;     ///< Array of frame starting position per chan.
  slim_channel **chan_array; ///< Array of ptrs to all channels.
  int num_allocated;         ///< Allocated size of the arrays (>= num_chan).

  void resize_arrays(int n);
};


//---------------------------------------------------------------------------
// Files
//---------------------------------------------------------------------------

class slim_compressor_t {
public:
  slim_compressor_t(const char *out_name, 
		    char flags_in, bool deltas=false, int samplepct_in=50);
  virtual ~slim_compressor_t();
  void close_output();

  slim_channel_encode * add_channel(slim_channel_encode *c);
  slim_channel_encode * add_channel(int reps, enum code_t code, 
				    enum data_t data_type, 
				    bool deltas, bool rotate);
  int num_channels() const;
  void set_section_frames(unsigned int nf); 

  size_t write(const unsigned char *buf, size_t max);
  size_t write_onesection(const unsigned char *buf, size_t max);
  int compress_from_file(const char *raw_file_name);

public:
  // Public attributes--so that they can be faked by caller.
  time_t mtime;          ///< Last modification time of raw file.
  size_t raw_size;       ///< Raw file (uncompressed) size (bytes).

public:
  // Inline read-only access to private data.
  size_t get_raw_size() const {return raw_size;} ///< Read raw file size.
  size_t get_frame_size() const {return frame_size;} ///< Read raw frame size.

private: 
  // Private methods
  bool no_reps() const;
  int num_data(int chan_num, int frames_used=-1) const;
  int write_file_header(const char *in_filename);
  int write_section_header();
  size_t encode_write_section(size_t length);
  int compute_section_params(size_t length);
  long data_offset(int i_data, int chan_num);
  void clear_channel_history();
  void confirm_flags();
  void write_last_section_foot();

private:
  // Private attributes.
  char *out_filename;    ///< Path of the output (compressed) data file.

  char flags;            ///< File header flags.

  slim_channel_array channels; ///< Easily-resizable vector of slim_channel_encode.
  size_t frame_size;     ///< Size (bytes) of a single frame.
  ///////  size_t this_sect_size; ///< Size (bytes) of the current section.
  size_t total_bytes_compressed; ///< Total compressed so far.
  unsigned int num_frames;///< Number of frames in the current section.
  int max_frames_per_section; ///< Limit on # of frames in each data section.
  int sections_written;  ///< Number of sections written to disk.
  int sample_pct;        ///< What percent of data to sample.
  raw_section *section;  ///< Buffer for holding entire section in memory.
  unsigned char *curptr; ///< Points beyond currently filled part of section buffer.
  size_t sec_bytes_stored;///< Bytes stored into current section.
  obitstream *ob;        ///< The bitstream for compressed output.
  bool encode_deltas;    ///< Should all channels encode deltas.
};



class slim_expander_t {
public:
  slim_expander_t(const char *in_name);
  ~slim_expander_t();

  int num_channels() const;
  int expand_to_file(const char *raw_file_name);
  int expand_to_stdout();
  size_t read(unsigned char *buf, size_t max);
  size_t read_onesection(const unsigned char **bufptr);
  int dump_sliminfo(); 

public:
  // Inline methods for read-only access to attribues.
  time_t get_mtime() const {return mtime;} ///< Read slim file's mtime.
  size_t get_rawsize() const {return raw_size;} ///< Read slim file's raw size.
  size_t get_slimsize() const {return slim_size;} ///< Read slim file size.
  void   set_ignore_crc(bool ic=true) {ignore_crc=ic;} ///< Ignore CRCs
  bool   is_open() {return ib->is_open();}  ///< Is low-level file open?


private:
  // Private methods
  int read_file_header();
  int read_section_header();
  size_t load_decode_section();
  slim_channel_decode * add_channel(slim_channel_decode *c, int bit_rotat);
  slim_channel_decode * add_channel(int reps, enum code_t code, 
				    enum data_t data_type,
				    bool deltas, int bit_rotat);

private:
  // Private attributes
  char *in_filename;     ///< Path of the input (compressed) data file.
  time_t mtime;          ///< Last modification time of raw file.
  size_t raw_size;       ///< Raw file (uncompressed) size (bytes).
  size_t slim_size;      ///< Slim file (compressed) size (bytes).
  char flags;            ///< File header flags.

  size_t bytes_read;     ///< Bytes read (raw) from file.
  size_t sec_bytes_read; ///< Bytes read (raw) from this section.
  size_t current_section_size; ///< Size of the currently open section.
  bool eof_tag_found;    ///< Have we read the End-of-File tag?

  raw_section *section;  ///< Buffer for holding entire section in memory.
  unsigned char *curptr; ///< Points to decoded, unconsumed data in section.
  bool used_read;        ///< User has called ::read()
  bool used_r_onesection;///< User has called ::read_onesection()
  bool ignore_crc;       ///< Do we ignore the CRC-32 checking.

  slim_channel_array channels;///< Array of all channels in this data file.
  unsigned int num_frames;///< Number of frames in the current section.
  ibitstream *ib;        ///< The bitstream for compressed input
};





//---------------------------------------------------------------------------
// Channels: single fields of data to be encoded/decoded.
//---------------------------------------------------------------------------

// A single data channel, for either encoding or decoding.

class slim_channel {
public:
  slim_channel(unsigned int reps, int size, bool deltas);
  virtual ~slim_channel();
  int get_deltas() const {return encode_deltas;} ///< Whether deltas are used.
  int get_bit_rotation() const {return bit_rotation;} ///< Read bit_rotation
  int get_repetitions() const {return repetitions;} ///< Number of reps / frame.
  size_t get_raw_size() const {return raw_size;} ///< Raw word size.
  size_t get_frame_size() const {return frame_size;} ///< Size of a frame.
  virtual void reset_previous() {;}  ///< Clear the "history" value.

public:
  slim_channel *next_chan;       ///< Next channel in a linked list

protected:
  const unsigned int repetitions;///< Number of repetitions of data per frame.
  const size_t raw_size;         ///< The size (bytes) of a single data word.
  const size_t frame_size;       ///< The size (bytes) of all data words per frame.
  int bit_rotation;              ///< Cyclically rotate raw data by X bits.
  int bit_unrotation;            ///< Cyclically rotate data by Y bits->raw.
  bool encode_deltas;            ///< Whether to encode the deltas
};


// Specialization of slim_channel for encoding only.

class slim_channel_encode : public slim_channel {
public:
  slim_channel_encode(int reps, int size, 
		      bool deltas, bool permit=false);
  ~slim_channel_encode();

  bool set_output(char * const out_name);
  bool set_output(obitstream *ob);
  void set_encoder(encoder *e);
  template <typename T> int compute_params(T *data, int ndata);
  int write_params() const;
  size_t encode_frame(void *buf);
  size_t encode_partial_frame(void *buf, size_t size);
  size_t encode_frame_singlevalue(void *buf);
  bool expect_zero_compression() const;
  encoder *replace_encoder();
  encoder *replace_constant(int d0);
  encoder *restore_encoder();
  virtual void reset_previous();

private:
  encoder *enc;              ///< Encoder object for this channel.
  obitstream *ob;            ///< Output-bitstream for writing encoded data.
  encoder *usual_encoder;    ///< Encoder object to use on most sections.
  bool permit_rotation;      ///< Is bit-rotation permitted?
  bool usual_deltas;             ///< Whether the usual encoder does deltas.
  int ndata_sampled;             ///< How many data values were used in sampling.
  const static int MIN_SAMPLES=5;///< How many data are required for a valid sample.

private:
  template <typename T>
  int constant_low_bits(const T *data, int ndata) const;
  uint32_t rotate(uint32_t u) const;
  uint32_t rotate(int i) const;
};


// Specialization of slim_channel for decoding only.

class slim_channel_decode : public slim_channel {
public:
  slim_channel_decode(int reps, int size, bool deltas);
  ~slim_channel_decode();

  bool set_input(char * const in_name);
  bool set_input(ibitstream *in_bs);
  void set_decoder(decoder *d);
  int read_params(int bit_rotation_in);
  size_t decode_frame(void *buf, size_t size);
  size_t decode_frame_singlevalue(void *buf);
  void dump_info(ostream &fout=cout) const;

private:
  decoder *dec;             ///< Decoder object for this channel.
  ibitstream *ib;           ///< Input-bitstream for reading encoded data.

public:
  uint32_t rotate(uint32_t u) const;
  uint32_t rotate(int i) const;
};



//---------------------------------------------------------------------------
// raw_section is the buffer to hold a raw data _section_ for improved
// I/O efficiency.
//---------------------------------------------------------------------------

enum section_mode_t {
  SECTION_COMPRESS_MODE,     ///< This is a section buffer for compressing.
  SECTION_EXPAND_MODE,       ///< This is a section buffer for expanding.
};

enum {
  MAX_SECTION_LENGTH=0x1000000, ///< Maximum raw size of a section (bytes). 
};


class raw_section {
public:
  typedef unsigned char Buffer_t;  ///< Alias for fundamental buffer data type.

  raw_section(enum section_mode_t mode_in);
  ~raw_section();

  void add_channel(int reps, int chan_size);
  void reset_channels();
  size_t fill(FILE *fp, size_t size);
  size_t flush(FILE *fp, size_t size);
  size_t resize(size_t size);
  unsigned long crc(size_t length=0) const;
  int set_num_frames(int i);

  void use_external_buffer(const unsigned char *buffer, size_t length);
  void use_internal_buffer();

  int8_t& cval(int ichan, int i);
  int8_t& cval(int ichan, int iframe, int i);
  int16_t& sval(int ichan, int i);
  int16_t& sval(int ichan, int iframe, int i);
  int32_t& ival(int ichan, int i);
  int32_t& ival(int ichan, int iframe, int i);
  uint32_t& uval(int ichan, int i);
  uint32_t& uval(int ichan, int iframe, int i);
  unsigned char *ptr(int ichan, int iframe) const;

public:
  // Public inline methods for read-only access to private data.
  size_t get_size() const {return buf_size;} ///< Read the raw sect size.
  size_t get_framesize() const {return frame_size;}

private:
  Buffer_t *private_buf;     ///< The private raw data buffer.
  Buffer_t *buf;             ///< The section raw data buffer (=private, usually).
  size_t private_buf_size;   ///< Private buffer allocated size
  size_t buf_size;           ///< Buffer allocated size
  size_t frame_size;         ///< Size (bytes) of each data frame.
  int num_frames;            ///< Number of whole frames stored in this buffer.
  const enum section_mode_t 
    mode;                    ///< Is this a write or read buffer?

  unsigned int *chan_reps;   ///< List of number of repetitions per chan.
  int *chan_offset;          ///< List of channel offsets within a frame.
  int n_chan_alloc;          ///< Number of channels allocated.
  int n_chan_used;           ///< Number of channels used.
};



//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
// Encoders/decoders: Objects containing the compression/decompression 
// algorithm.
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
// encoder / decoder: Base class for all compression/decompression.  Can
// be used, but they contain only a trivial non-compressing algorithm.
//---------------------------------------------------------------------------

class encoder {
public:
  encoder(enum data_t dt, bool deltas, obitstream *ob=NULL);
  virtual ~encoder();

  bool set_output(obitstream *ob);
  bool set_data_type(enum data_t dt_in);
  void use_signed_data_type();
  
  virtual void encode(uint32_t datum) const;
  virtual void encode(uint16_t datum) const;
  virtual void encode(uint8_t datum) const;
  void encode_scalar(const uint32_t *data);
  void encode_scalar(const uint16_t *data);
  void encode_scalar(const uint8_t *data);
  virtual void encode_vector(const uint32_t *data, int ndata);
  virtual void encode_vector(const uint16_t *data, int ndata);
  virtual void encode_vector(const uint8_t *data, int ndata);
  virtual int compute_params(const uint32_t *data, const int ndata);
  virtual int compute_params(const uint16_t *data, const int ndata);
  virtual int compute_params(const uint8_t *data, const int ndata);
  virtual int write_params() const;
  virtual bool expect_zero_compression() const;
  bool is_signed() const;
  virtual encoder *replacement_encoder();
  virtual encoder *constant_encoder(int d0);
  enum data_t get_data_type() const {return data_type;} ///< Read data_type.
  void reset_previous() {prev_datum = 0u; prev_sdatum = 0u;}   ///< Clear delta history.

protected:
  const bool use_deltas; ///< Whether to encode successive difference values.
  obitstream *out_bs;    ///< The bitstream for writing encoded data.
  enum data_t data_type; ///< The type to be encoded.
  unsigned int data_size_bytes;   ///< Raw data word size in bytes
  unsigned int data_size_bits;    ///< Raw data word size in bits
  uint32_t prev_datum;        ///< Previous value for deltas (uint32_t).
  uint16_t prev_sdatum; ///< Previous value for deltas (uint16_t).
  uint8_t prev_cdatum; ///< Previous value for deltas (uint8_t).

  template <typename T>
  void compute_mean(double& mean,
		    const T *data, int ndata) const;
  template <typename T>
  void compute_mean_and_mad(double& mean, double& meanabsdev,
			    const T *data, int ndata) const;

private:
  // Base default encoder has no parameters to the algorithm.
  const static enum code_t 
  ALGORITHM_CODE = SLIM_ENCODER_DEFAULT;///< Encoder algorithm
};


class decoder {
public:
  decoder(enum data_t dt, bool deltas, ibitstream *ib=NULL);
  virtual ~decoder();

  virtual int read_params();
  virtual void dump_info(ostream &fout=cout) const;

  bool set_input(ibitstream *in_bs);
  bool set_data_type(enum data_t dt_in);

  void decode_scalar(uint32_t *data);
  void decode_scalar(uint16_t *data);
  void decode_scalar(uint8_t *data);
  void decode_vector(uint32_t *data, int ndata=1);
  void decode_vector(uint16_t *data, int ndata=1);
  void decode_vector(uint8_t *data, int ndata=1);

protected:
  virtual uint32_t decode_u32();
  virtual uint16_t decode_u16();
  virtual uint8_t decode_u8();

protected:
  const bool use_deltas; ///< Whether to encode successive difference values.
  ibitstream *in_bs;     ///< The bitstream for reading encoded data.
  enum data_t data_type; ///< The type to be encoded.
  int data_size_bytes;   ///< Raw data word size in bytes
  int data_size_bits;    ///< Raw data word size in bits
  uint32_t prev_datum;   ///< Previous value for deltas (uint32_t).
  uint16_t prev_sdatum;  ///< Previous value for deltas (uint16_t).
  uint8_t prev_cdatum;   ///< Previous value for deltas (uint8_t).

private:
  // Base default encoder has no parameters to the algorithm.
};

//---------------------------------------------------------------------------
// encoder_reduced_binary / decoder_reduced_binary: 
// A first stab at compressing by limiting the number of bits used for (most)
// values in the data.
//---------------------------------------------------------------------------

// Derived class for encoding simple full reduced_binary of data.
class encoder_reduced_binary : public encoder {
public:
  encoder_reduced_binary(enum data_t dt, bool deltas, obitstream *ob=NULL);
  virtual ~encoder_reduced_binary();

  void encode(uint32_t datum) const;
  void encode(uint16_t datum) const;
  void encode(uint8_t datum) const;
  virtual int compute_params(const uint32_t *data, const int ndata);
  virtual int compute_params(const uint16_t *data, const int ndata);
  virtual int write_params() const;
  virtual bool expect_zero_compression() const;

protected:
  virtual int overflow_waste(const int histogram[33], unsigned int n);
  int best_code_length(const int histogram[33], int ndata);

protected:
  unsigned int nbits;     ///< Number of bits per encoded symbol
  uint32_t max;       ///< Maximum codable value (after offset removed)
  uint32_t offset;    ///< Offset (subtract to encode; add back to decode)
  uint32_t Overflow;  ///< Special overflow (range failure) symbol.

private:
  const static enum code_t ALGORITHM_CODE = SLIM_ENCODER_REDUCED_BINARY;///< ID code #
};


// Derived class for decoding simple full range of data.
class decoder_reduced_binary : public decoder {
public:
  decoder_reduced_binary(enum data_t dt, bool deltas, ibitstream *ib=NULL);
  virtual ~decoder_reduced_binary();

  virtual int read_params();
  virtual void dump_info(ostream &fout=cout) const;

protected:
  virtual uint32_t decode_u32();
  virtual uint16_t decode_u16();

  unsigned int nbits; ///< Number of bits per encoded symbol
  uint32_t max;       ///< Maximum codable value
  uint32_t offset;    ///< Offset (subtract to encode; add back to decode)
  uint32_t Overflow;  ///< Special overflow (range failure) symbol.

private:
  const static enum code_t ALGORITHM_CODE = SLIM_ENCODER_REDUCED_BINARY;///< ID code #
};


//---------------------------------------------------------------------------
// encoder_runlength / decoder_runlength:
// Class for encoding highly repetitive data streams.
//---------------------------------------------------------------------------

class encoder_runlength : public encoder {
public:
  encoder_runlength(enum data_t dt, bool deltas, obitstream *ob=NULL);
  virtual ~encoder_runlength();

  virtual void encode(uint32_t datum) const;
  virtual void encode_vector(const uint32_t *data, int ndata);
  virtual int compute_params(const uint32_t *data, const int ndata);
  virtual int compute_params(const uint16_t *data, const int ndata);
  virtual int write_params() const;
  virtual bool expect_zero_compression() const;
  virtual encoder *replacement_encoder();

protected:
  int ndata_checked;    ///< Number of data tested in compute_params
  int nchanges_data;    ///< Number of data changes detected in compute_params().

private:
  const static enum code_t ALGORITHM_CODE = SLIM_ENCODER_RUNLENGTH;///< ID code #
};


// Derived class for decoding simple full range of data.
class decoder_runlength : public decoder {
public:
  decoder_runlength(enum data_t dt, bool deltas, ibitstream *ib=NULL);
  virtual ~decoder_runlength();

  virtual int read_params();
  virtual void dump_info(ostream &fout=cout) const;

protected:
  virtual uint32_t decode_u32();

  uint32_t repeated_value;  ///< Value to be used in future.
  uint32_t uses_remaining;  ///< Undecoded consecutive uses of value.

private:
  const static enum code_t ALGORITHM_CODE = SLIM_ENCODER_RUNLENGTH;///< ID code #
};




//---------------------------------------------------------------------------
// encoder_constant / decoder_constant:
// Class for encoding strictly constant data streams.
//---------------------------------------------------------------------------

class encoder_constant : public encoder {
public:
  encoder_constant(int value, enum data_t dt, bool deltas, obitstream *ob=NULL);
  virtual ~encoder_constant();

  virtual void encode(uint32_t datum) const;
  virtual void encode(uint16_t datum) const;
  virtual void encode(uint8_t datum) const;
  virtual int write_params() const;

protected:
  uint32_t fixed_data;    ///< The fixed value for this channel.
  uint16_t fixed_sdata;   ///< The fixed value for this channel.
  uint8_t fixed_cdata;    ///< The fixed value for this channel.

private:
  const static enum code_t ALGORITHM_CODE = SLIM_ENCODER_CONSTANT;///< ID code #
};



class decoder_constant : public decoder {
public:
  decoder_constant(enum data_t dt, bool deltas, ibitstream *ib=NULL);
  virtual ~decoder_constant();

  virtual int read_params();
  virtual void dump_info(ostream &fout=cout) const;

protected:
  virtual uint32_t decode_u32();
  virtual uint16_t decode_u16();
  virtual uint8_t decode_u8();

  uint32_t fixed_data;    ///< The fixed value for this channel.
  uint16_t fixed_sdata;   ///< The fixed value for this channel.
  uint8_t fixed_cdata;    ///< The fixed value for this channel.

private:
  const static enum code_t ALGORITHM_CODE = SLIM_ENCODER_CONSTANT;///< ID code #
};




encoder * encoder_generator(enum code_t code, enum data_t data_type, bool deltas);
decoder * decoder_generator(enum code_t code, enum data_t data_type, bool deltas);

#endif  // #ifndef SLIM_H
