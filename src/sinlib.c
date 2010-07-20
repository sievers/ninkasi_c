#include <stdio.h>
#include <math.h>
#include <stdlib.h>

typedef double mytype;
//const mytype sin4_params[3]={0.999771340672793,-0.16582687904326,0.00757417450801057};                                                                                                                           
const mytype sin4_params[3]={   9.99694901254465e-01,  -1.65670041043745e-01,   7.51339148607193e-03};
//const mytype sin4_params_flip[3]={     -1.65670041043745e-01,    7.51339148607193e-03,  9.99694901254465e-01};                                                                                                   

//const mytype sin6_params[4]={0.999997486206621 ,-0.166651676590364, 0.00830951228208954,-0.000184470857616342};                                                                                                  
const mytype sin6_params[4]={   9.99996601050174e-01,  -1.66648235616735e-01,   8.30628614181131e-03,  -1.83627485772683e-04};

//const mytype sin8_params[5]={   9.99999982781442e-01,  -1.66666515196234e-01,   8.33296400727452e-03,  -1.98047545833302e-04,   2.59810891065399e-06};                                                           
const mytype sin8_params[5]={     9.99999976513754e-01,  -1.66666475934886e-01,   8.33289922280274e-03,  -1.98008653054027e-04,   2.59043003065907e-06};


const mytype cos3_params[3]={   9.99579502755762e-01,  -4.96392260254040e-01,   3.72092848991399e-02};
const mytype cos5_params[4]={   9.99995282503834e-01,  -4.99930917750131e-01,   4.15117334679269e-02,  -1.27871281238398e-03};
const mytype cos7_params[5]={     9.99999967268498e-01,  -4.99999268956354e-01,   4.16640912969955e-02,  -1.38574213281038e-03,   2.32376335545823e-05};


const mytype asin3_params1[7]={   2.29287730500315e+00,  -4.84236566639965e+00,   4.12272163981868e+00,  -1.47181251782873e+00,   3.05738342281475e-01,   9.77273703943210e-01,   3.98168115947950e-04};
const mytype asin3_params2[3]={  -8.19252532048417e-02,  -1.39802729122732e+00,   1.56995464465493e+00};
const mytype asin_split=0.85;

static inline mytype asin3(mytype x)
{
  if (x<asin_split)
    return asin3_params1[6]+x*(asin3_params1[5]+x*(asin3_params1[4]+x*(asin3_params1[3]+x*(asin3_params1[2]+x*(asin3_params1[1]+x*asin3_params1[0])))));
  else {
    mytype xx=sqrt(1-x);
    return asin3_params2[2]+xx*(asin3_params2[1]+xx*asin3_params2[0]);
  }

}


const mytype asin4_params1[5]={   6.32559537178112e-05,   9.97002719101181e-01,   3.23729856176963e-02,   3.89287300071597e-02,   1.93549238398372e-01};
const mytype asin4_params2[7]={   2.09625797161885e+01,  -1.74835553411477e+02,   6.13575281494908e+02,  -1.14033116228467e+03,   1.19159992307311e+03,  -6.63957441058529e+02,   1.54421991537526e+02};
const mytype asin4_params3[4]={   1.57080010233116e+00,  -1.41437401362252e+00,   1.84777752400778e-03,  -1.24625163381900e-01};
const mytype asin4_split1=0.6;
const mytype asin4_split2=0.925;
static inline mytype asin4(mytype x)
{
  if (x<asin4_split1)
    return asin4_params1[0]+x*(asin4_params1[1]+x*(asin4_params1[2]+x*(asin4_params1[3]+x*(asin4_params1[4]))));
  if (x<asin4_split2)
    return asin4_params2[0]+x*(asin4_params2[1]+x*(asin4_params2[2]+x*(asin4_params2[3]+x*(asin4_params2[4]+x*(asin4_params2[5]+x*asin4_params2[6])))));
  mytype xx=sqrt(1-x);
  return asin4_params3[0]+xx*(asin4_params3[1]+xx*(asin4_params3[2]+xx*asin4_params3[3]));
}


static inline mytype cos3(mytype x)
{
  mytype xx=x*x;
  return (cos3_params[0]+xx*(cos3_params[1]+xx*cos3_params[2]));
}



static inline mytype sin42(mytype x)
{
  const mytype xx=x*x;
  return x*(  9.99694901254465e-01+xx*( -1.65670041043745e-01 + xx* 7.51339148607193e-03));
}
static inline  mytype sin4(mytype x)  //good to 1.6*10^-4 between -pi/2 and pi/2                                                                                                                                   
{
  const mytype xx=x*x;
  return x*(sin4_params[0]+xx*(sin4_params[1]+xx*sin4_params[2]));

}
static inline mytype cos5(mytype x)
{
  mytype xx=x*x;
  return (cos5_params[0]+xx*(cos5_params[1]+xx*(cos5_params[2]+xx*cos5_params[3])));
}


static inline mytype sin6(mytype x)  //good to 6*10^-7 between -pi/2 and pi/2                                                                                                              
{
  mytype xx=x*x;
  return x*(sin6_params[0]+xx*(sin6_params[1]+xx*(sin6_params[2]+xx*sin6_params[3])));
}


static inline mytype cos7(mytype x)
{
  mytype xx=x*x;
  return (cos7_params[0]+xx*(cos7_params[1]+xx*(cos7_params[2]+xx*(cos7_params[3]+xx*cos7_params[4]))));
}


static inline mytype sin8(mytype x)  //good to <10^-8 between -pi/2 and pi/2                                                                                                                                       
{
  mytype xx=x*x;
  return x*(sin8_params[0]+xx*(sin8_params[1]+xx*(sin8_params[2]+xx*(sin8_params[3]+xx*sin8_params[4]))));
}


/*--------------------------------------------------------------------------------*/
mytype get_max_err(mytype x1,mytype x2,mytype dx,mytype (*f1)(mytype),mytype (*f2)(mytype))
{
  mytype x,maxerr,dy;
  maxerr=0;
  for (x=x1;x<=x2;x+=dx) {
    dy=fabs((f1)(x)-(f2)(x));
    if (dy>maxerr)
      maxerr=dy;
  }

  return dy;

}

/*================================================================================*/
int main(int argc, char *argv[])
{
  mytype n=400*60*15;
  mytype  m=1000;

#if 0
  printf("maxsin4 err is %14.4e\n",get_max_err(-M_PI/2.0,M_PI/2.0,0.0001,sin,sin4));
  printf("maxsin6 err is %14.4e\n",get_max_err(-M_PI/2.0,M_PI/2.0,0.0001,sin,sin6));
  printf("maxsin8 err is %14.4e\n",get_max_err(-M_PI/2.0,M_PI/2.0,0.0001,sin,sin8));
  printf("maxasin4 err is %14.4e\n",get_max_err(0,1,0.0001,asin,asin4));
  return 0;
#endif
#if 0

  mytype x=atof(argv[1]);
  printf("sin is %14.8f, fit is %14.8f, err is %14.4e\n",sin(x),sin4(x),(sin4(x)-sin(x)));
  printf("sin is %14.8f, fit is %14.8f, err is %14.4e\n",sin(x),sin6(x),(sin6(x)-sin(x)));
  printf("sin is %14.8f, fit is %14.8f, err is %14.4e\n",sin(x),sin8(x),(sin8(x)-sin(x)));
  printf("cos is %14.8f, fit is %14.8f, err is %14.4e\n",cos(x),cos3(x),(cos3(x)-cos(x)));
  printf("cos is %14.8f, fit is %14.8f, err is %14.4e\n",cos(x),cos5(x),(cos5(x)-cos(x)));
  printf("cos is %14.8f, fit is %14.8f, err is %14.4e\n",cos(x),cos7(x),(cos7(x)-cos(x)));

  printf("asin is %14.8f, fit is %14.8f, err is %14.4e\n",asin(x),asin3(x),asin3(x)-asin(x));
#endif
  mytype i,j;
  mytype sum=0;
  for (i=0;i<n;i++)
    for (j=0;j<m;j++)
      {
        double val=0.001*j;
        //sum+=fabs(asin4(val)-asin(val));                                                                                                                                                                         
        sum +=sin4(val);
        //sum+=cos5(1.57-val);                                                                                                                                                                                     
      }
  printf("sum is %18.8e\n",sum/n/m);

  return 0;
}


