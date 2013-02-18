#ifndef NOISE_TYPES_H
#define NOISE_TYPES_H

#define MB_NOISE_MAX_PARAM 20  //put this in for now so that we may not have to change
                               //the depot for future, more parameter-endowed noise models.


///Define noise model types
typedef enum {  MBNOISE_UNDEFINED,     ///< Unknown noise model
                MBNOISE_LINEAR_POWLAW,  ///< model as linear combination of power laws 
		MBNOISE_CONSTANT,      ///  Noise is a constant
		MBNOISE_INTERP,        ///  Noise is fit as in constant, but linearly interpolated when applied
		MBNOISE_FULL           ///  Save a whole timestream
		

} mbNoiseType;


/************************************************************************************************************/
/*!
 *  Structure to hold the parameters describing the noise for a single pixel.
 */

typedef struct {
  //int noise_type;
  mbNoiseType noise_type;
  int initialized;  
  actData knee;
  actData old_knee;
  actData white;     //The white level in cts^2/hz 
  
  actData params[MB_NOISE_MAX_PARAM]; //Two parameters is all we need for simple 1/f.
  int nparam;
  //actData powlaw[MB_NOISE_MAX_PARAM]; //Be able to specify an index other than 1/f
  actData powlaw; //Be able to specify an index other than 1/f
  actData minfreq;  //If >0, lower limit of frequency for 1/f fitting, in case constant term is messed up.
  actData maxfreq;  //If >0, upper limit of frequency for 1/f fitting. 
  actData centerfreq; //Frequency over which we expect the noise to be flat.
  int converged;          //set to true if the noise-fitting found a solution with which it's happy.
  
} NoiseParams1Pix;



/************************************************************************************************************/
/*!
 *  Structure to hold the parameters describing the noise for an array
 */


typedef struct {
  int ndet;
  NoiseParams1Pix *noises;

} mbNoiseVectorStruct;










/*--------------------------------------------------------------------------------*/
typedef struct {
  mbNoiseType noise_type;
  actData dt;
  actData nu_low;
  actData nu_high;
  
  int i_low;
  int i_high;

  actData white;
  actData knee;



  actData *noise_data;

} mbNoiseParams1PixBand;


/************************************************************************************************************/
/*!
 *  Structure to hold the parameters describing noise as combination of detector and correlated
 *  noise in bands
 */

typedef struct {
  int ndet;
  int nband;
  int *band_edges;
  int *nvecs;
  actData **noises;
  actData ***vecs;
  
}  mbNoiseStructBandsVecs;


/************************************************************************************************************/
/*!
 *  Structure to hold the parameters describing the noise for an array in frquency bands
 */


typedef struct {
  int ndet;
  actData deltat;
  int nband;
  actData *bands;
  int *ibands;

  
  bool *do_rotations;


  actData ***rot_mats;
  actData ***inv_rot_mats_transpose;
  
  mbNoiseParams1PixBand **noise_params;
  
  
} mbNoiseVectorStructBands;



/*--------------------------------------------------------------------------------*/
/*structure to hold hwp demodulated things*/
typedef struct {
  actComplex **data;  
  actData hwp_freq;
  actData lowpass_freq; //in units of the hwp
  actData lowpass_taper; 
  actData highpass_freq; //in units of the hwp
  actData highpass_taper;
  actData *freqs;  //if you want multiple demodulated timestreams, this are the frequencies
  int nfreq;
  int nmode;  //number of frequencies to keep
} DemodData;


#endif
