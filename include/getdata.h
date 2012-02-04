
#ifndef GETDATA_H
#define GETDATA_H

#define GD_E_OK                0
#define GD_E_OPEN_FORMAT       1
#define GD_E_FORMAT            2
#define GD_E_FIELD             4
#define GD_E_BAD_CODE          5
#define GD_E_BAD_RETURN_TYPE   6
#define GD_E_OPEN_RAWFIELD     7
#define GD_E_OPEN_INCLUDE      8
#define GD_E_OPEN_RAWFIELD_    10
#define GD_E_SIZE_MISMATCH    12
#define GD_E_OPEN_LINFILE     13
#define GD_E_RECURSE_LEVEL    14

#ifdef __cplusplus
extern "C" {
#endif

extern const char *GD_ERROR_CODES[15];
struct FormatType;

/***************************************************************************/
/*                                                                         */
/*  GetData: read BLAST format files.                                      */
/*    filename_in: the name of the file directory (raw files are in here)  */
/*    field_code: the name of the field you want to read                   */
/*    first_frame, first_samp: the first sample read is                    */
/*              first_samp + samples_per_frame*first_frame                 */
/*    num_frames, num_samps: the number of samples read is                 */
/*              num_samps + samples_per_frame*num_frames                   */
/*    return_type: data type of *data_out.  's': 16 bit signed             */
/*              'u' 16bit unsiged. 'S' 32bit signed 'U' 32bit unsigned     */
/*              'c' 8 bit unsigned                                         */
/*    void *data_out: array to put the data                                */
/*    *error_code: error code is returned here. If error_code==null,       */
/*               GetData prints the error message and exits                */
/*                                                                         */
/*    return value: returns number of samples actually read into data_out  */
/*                                                                         */
/***************************************************************************/
int GetData(const struct FormatType *F, const char *field_code,
             int first_sframe, int first_samp,
             int num_sframes, int num_samp,
             char return_type, void *data_out,
            int *error_code);

/***************************************************************************/
/*                                                                         */
/*    Get the number of samples for each frame for the given field         */
/*                                                                         */
/***************************************************************************/
int GetSamplesPerFrame(const struct FormatType *F, const char *field_name, int *error_code);

/***************************************************************************/
/*                                                                         */
/*    Get the number of frames available                                   */
/*                                                                         */
/***************************************************************************/
int GetNFrames(const struct FormatType *F, int *error_code, const char *field);

/***************************************************************************/
/*                                                                         */
/*    Get the number format structure                                      */
/*                                                                         */
/***************************************************************************/
struct FormatType *GetFormat(const char *filedir, const char *linterp_prefix, int *error_code);

/***************************************************************************/
/*                                                                         */
/*    Close all open file descriptors for a given dirfile                  */
/*                                                                         */
/***************************************************************************/
void GetDataClose(struct FormatType *F);


#ifdef __cplusplus
} // extern "C"
#endif

#endif // GETDATA_H
