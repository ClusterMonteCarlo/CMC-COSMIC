/* software by kabel: last modified on Feb 18, 2002 */
/* header for single.c and binary.c */


#define Vicky 1                     /* 0--evol continus, 1--evol stops in cases pointed out by Vicky */

#define num_tested 1000000                                               /* number of tested systems */
#define NOUT 100000                        /* frequency of output into info.dat (every NOUT systems) */

#define Random 1                   /* 1-- uses ran3.c function to generate random numbers (standard) */
                                              /* 2-- uses ran2.c function to generate random numbers */

#define Pi (2.0*asin(1.0))
#define Rsun 6.9599e+10                                             /* [cm], sun radius in cgs units */
#define Msun 1.989e+33                                                 /* [g], sun mass in cgs units */
#define MCh 1.44                                      /* Chandrasekhar mass [Msun] (WD/NS formation) */

#define ZZ METALLICITY                                         /* metallicity (solar value assumed for now) */
#define Mmaxns 3.0                                         /* maximum mass of NS, over it we have BH */

#define cflag 3          /* 1-Jarrod prescription of compact object mass, 2-my for single stars only */
                                            /* 3 - STANDARD my for binaries, fb calculated with wind */
                                  /* 4-my for binaries, fb calculated with no wind (used in singl.c) */


/*#define Sal (-2.7) */                       /* exponent of IMF function for get_M(), Salpeter is -2.35 */
#define Mhmxb 5.0                /* minimal current mass of donor in High Mass X-ray Binaries [Msun] */
#define Mlmxb1 1.4           /* mass of donor (K=0,1)  M_donor< Mlmxb1 x mass of acceptor to be LMXB */
#define Mlmxb2 0.666667    /* mass of donor(2,3,4,5,6) M_donor< Mlmxb2 x mass of acceptor to be LMXB */
#define Mhecon 4.5              /* below this mass Helium giant has convective env., otherwise it is */ 
                                                          /* radiative, so no common envelope issues */ 

#define NEW1 1                               /* 1--new MB,GR tidial breaking treatment, 0--old (PhD) */

#define MB 1                                     /* 0-- no magnetic breaking, 1-- new MB, 2-- old MB */
#define gamMB 2                                      /* only for MB=2, exponent in the old law [0-4] */
#define minMB 0.35                               /* minimum mass of MS star over which MB is applied */
#define maxMB 1.5                              /* maximum mass of MS star below which MB is applied */
                                   /* MS star: below minMB fully convective, over maxMB no conv env. */ 

#define WIND WIND_FACTOR           /* calculated wind mass loss (Hurley prescription) is multiplied by WIND */
                            /* 1.0 -- standard */

#define HCE 0                                 /* 0- NS in any CE do accrete with hyper critical rate */
                                           /* 1-NS in CE with H-rich stars do not accrete ANY matter */
                                                      /* 2-NS/BH do not accrete ANY matter in ANY CE */
                                      
#define Rat (0.0)
#define MarkRat 1                                               /* 1 -- flat q distribution, Rat=0.0 */
                                                          /* 2 -- q^(Rat) distribution, Rat=positive */                                       
                               /* 3 -- q=const + q^(Rat) distribution (limit at q=0.2), Rat=negative */          


#define kick 2                          /* 1-Cordes&Chernoff kick distribution: NS and BH same kicks */
                      /* 2-BH kicks decreased, no direct collapse kicks otherwise like 1: (STANDARD) */
                                               /* 3-no BH kicks, otherwise like 2, 4-no kicks at all */
                                                 /* 5-Paczynski kicks distribution, otherwise like 2 */

#define Divide 0.8                                    /* 0.8 for Cordes&Chernoff kick distribution */
#define Sigma1 175.0                                  /* 175.0 for Cordes&Chernoff kick distribution */
#define Sigma2 700.0                                  /* 700.0 for Cordes&Chernoff kick distribution */
 
#define Alfa 1.0        /* was 0.8, efficiency of orbital energy loss for the ejection of common envelope */
#define Beta 1.0             /* Podsiadlowski et al. 1992 -specific angular momentum of matter [2Pia^2/P] */
#define Lambda 1.0             /* was 0.5, parameter describing strength of giant envelope binding energy */
#define Fa 0.5                 /* Meurs,Heuvel 1989 -fraction of mass lost by donor attached to companion */
#define Cd 6.0                        /* drag coefficent, see Bethe and Brown 1998 ApJ 506, 780 eq. (5.7) */
#define Cr1 1.0                 /* if M_donor<cr*M_acceptor then we have quasi-dynamic mass transfer/loss */
#define Cr2 2.5                                  /* if M_donor>=cr*M_acceptor then common envelope begins */
                                      /* cr2 -- donor K=2, acceptor K=0,1 (Vicky), cr1 -- all other cases */


#define acc 1.0e-10                                             /* accuracy for variable comparisons */

                                      /* time steps for calculation at different evolutionary stages */
#define delms 0.01                                                /* Main Sequence, HS Main Sequence */
#define delhg 0.01                                                                 /* Hartzprung Gap */
#define delrg 0.01                           /* Red Giant Branch, HS Hartzprung Gap, HS Giant Branch */
#define delhe 0.01                                                            /* Core Heluim Burning */
#define delag 0.01                                                        /* Asymptotic Giant Branch */


#define min(a,b) ((a)<(b)?(a):(b))
#define max(a,b) ((a)>(b)?(a):(b))


                                                                                 /* Global variables */ 
#define Hub 1                      /* Hub=0: outburst of star formation 10Gyrs ago, t_hubble=10000.0 */
                                 /* Hub=1: continous star form. through 10Gyrs, t_hubble=[0:10000.0] */
double t_hubble;        /* Hubble time [10^6 yrs] - for SN and mergers counters, defined in binary.c */
#define hub_val 10000.0                               /* value of t_hubble -- Hubble time [10^6 yrs] */ 
                            /* value for out Galatic disk is hub_val=10000.0, and maximum is 15000.0 */
            /* the age of univerese, but artificially Hubble can be supplied here as large as needed */
               /* evolution of any system stops at time = t_hubble which has maximum value of Hubble */
                          /* whether system or any particular star has finished its evolution or not */


#define MM 1            /* 1 --singl.c fills in tend, and may return at diffrent time than requested */
                                          /* by external driver (binary.c), 0--always return at tend */

int DNStype;               /* DNStype=1: non-recycled, DNStype=0: recycled, DNStype=-1: do not apply */
int DWDtype;               /* DWDtype=1: new path(for DNS would be non-recycled), DWDtype=1 old path */  
             /* initialized inside for{} loop, housekeeping of DNStype and DWDtype done in systype() */

double xPacz[2000],yPacz[2000];                          /* for Paczynski kicks filled in get_Vkick5 */    

FILE *fp0;                         /* file with errors and warnings, common for singl.c and binary.c */

double M_hook,M_HeF,M_FGB;                                           /* critical evolutionary masses */

double aa1,aa2,aa3,aa4,aa5,aa6,aa7,aa8,aa9,aa10;                          /* ZZ dependent parameters */
double aa11,aa12,aa13,aa14,aa15,aa16,aa17,aa18,aa19,aa20;
double aa21,aa22,aa23,aa24,aa25,aa26,aa27,aa28,aa29,aa30;
double aa31,aa32,aa33,aa34,aa35,aa36,aa37,aa38,aa39,aa40;
double aa41,aa42,aa43,aa44,aa45,aa46,aa47,aa48,aa49,aa50;
double aa51,aa52,aa53,aa54,aa55,aa56,aa57,aa58,aa59,aa60;
double aa61,aa62,aa63,aa64,aa65,aa66,aa67,aa68,aa69,aa70;
double aa71,aa72,aa73,aa74,aa75,aa76,aa77,aa78,aa79,aa80,aa81;

double bb1,bb2,bb3,bb4,bb5,bb6,bb7,bb9,bb10;
double bb11,bb12,bb13,bb14,bb15,bb16,bb17,bb18,bb19,bb20;
double bb21,bb22,bb23,bb24,bb25,bb26,bb27,bb28,bb29,bb30;
double bb31,bb32,bb33,bb34,bb36,bb37,bb38,bb39,bb40;
double bb41,bb42,bb43,bb44,bb45,bb46,bb47,bb48,bb49;
double bb50,bb51,bb52,bb53,bb54,bb55,bb56,bb57;   
/* bb50 set in Ragbf() called by single() not in coeff_bb() */


/* counters for different SN types and for fraction of NS formed in different SNs */
long int ss0,ss1a,ss1b,ss1c,ss2a,ss2b,ss2c,ss3a,ss3b,ss3c,ss4a,ss4b,ss4c;
long int ns1a,ns1b,ns1c,ns2a,ns2b,ns2c,ns3a,ns3b,ns3c,ns4a,ns4b,ns4c;

long int st0,st1a,st1b,st1c,st2a,st2b,st2c,st3a,st3b,st3c,st4a,st4b,st4c;
long int nt1a,nt1b,nt1c,nt2a,nt2b,nt2c,nt3a,nt3b,nt3c,nt4a,nt4b,nt4c;



long int dce1;                        /* number of Double Common Envelope events of H-rich stars */
long int dce2;               /* number of Double CE events of H-rich stars that result in merger */
long int dce3;                /* number of Double Common Envelope events of H-rich+He-rich stars */
long int dce4;       /* number of Double CE events of H-rich+He-rich stars that result in merger */
long int dce5;                       /* number of Double Common Envelope events of He-rich stars */
long int dce6;              /* number of Double CE events of He-rich stars that result in merger */
long int hce10;                                   /* number of hyper-crit Common Envelope events */


/* count formation of different eccentric and circular compact object binaries through double CE */
long int dce10,dce11,dce12,dce13,dce14,dce15,dce16,dce17,dce18,dce19,dce20,dce21;         

long int dce30,dce31,dce32,dce33,dce34,dce35,dce36,dce37;    /* mergeing and not merging compact */
long int dce40,dce41,dce42,dce43,dce44,dce45,dce46,dce47;  /* object binaries which went through */
long int dce50,dce51,dce52,dce53,dce54,dce55,dce56,dce57;        /* different types of double CE */
long int dce60,dce61,dce62,dce63,dce64,dce65,dce66,dce67;
long int dce70,dce71,dce72,dce73,dce74,dce75,dce76,dce77;
long int dce80,dce81,dce82,dce83,dce84,dce85,dce86,dce87;




long int hce1,hce2,hce3,hce4;                           /* number of Hyper Common Envelope events */ 
long int Hmxb1,Hmxb2,Hmxb3,Hmxb4;                              /* number of all (wind+RLOF) HMXBs */
long int hmxb1,hmxb2,hmxb3,hmxb4;                                         /* number of RLOF HMXBs */
long int lmxb1,lmxb2,lmxb3,lmxb4,lmxb5,lmxb6,lmxb7,lmxb8;                      /* number of LMXBs */
long int msns,msbh;                                             /* number of MS-NS, MS-BH systems */
long int hehe1;           /* number of MT of two He giants when both are more massive then Mhecon */
long int hehe2;           /* number of MT of two He giants when both are less massive then Mhecon */     
long int hehe3;   /* number of MT of two He giants when donor is not massive and acceptor massive */
long int hehe4;   /* number of MT of two He giants when donor is massive and acceptor not massive */




void singl(double *MzamsR, double *M0R, double *MR, int *KR, double *TBR, double *TVIRR, double *TER,
           double *LR, double *RR, double *MCR, double *MHER, double *MC0R, int *FLAGR, double *DTR, 
           double *MPRER, int *KPR, double *TSTART, double *FRACR);



void coef_aa(void);                                              /* Initializing functions */
void coef_bb(void);
double M_hookf(void);                                           
double M_HeFf(void);
double M_FGBf(void);


void ZAMSf(double *M, double *t, double *Mhe, double *Mco, int *stop);          /* Evolutionary functions */
int MSf(double *M, double *t, double *Mhe, double *Mco, int *stop);
int HGf(double *M, double *t, double *tstart, double *Mhe, double *Mco, int *Kp, int *stop);
int RGf(double *M, double *t, double *tstart, double *Mhe, double *Mco, int *Kp, int *stop);
int HEf(double *M, double *t, double *tstart, double *Mhe, double *Mco, int *stop);
int AGBf(double *M, double *t, double *tstart, double *Mhe, double *Mco, int *Kp, int *stop, int Klast);
int HSMSf(double *M, double *t, double *tstart, double *Mhe, double *Mco, int *stop);
int HSGBf(double *M, double *t, double tstart, double *Mhe, double *Mco, int *Kp, int *stop);
int REMNANTf(double *M, double *t, double *Mhe, double *Mco, int Kp);
void sn_type(double M, double Mhe, double Mco, double t, int K, int fb);
double windf(double M, double L, double R, double Mc, int K);
int perturb(double M, double Mc, double *L, double *R, double t, int K);                             
double tnf(double M, int K);

double Lzamsf(double M);                                         /* Main Sequance functions */
double Rzamsf(double M);
double thookf(double M);
double tmsf(double M);
double Lmsf(double M, double t);
double Rmsf(double M, double t);
double Ltmsf(double M);
double Rtmsf(double M);

double thgf(double M);                                           /* Hartzprung Gap functions */
double Mchgf(double M, double t);
double Lhgf(double M, double t);
double Rhgf(double M, double Mtmp, double t);
double Mcehgf(double M);
double Lehgf(double M);
double Rehgf(double M, double Mtmp);

double tbgbf(double M);                                          /* Red Giant Branch functions */
double Mcbgbf1a(double M);                                       /* real core mass on BRGB */
double Mcbgbf1b(double M);                                       /* dummy core mass on BRGB */
double Lbgbf(double M);
double Ahpf(double M);
double Bf(double M);
double Df(double M);
double pf(double M);
double qf(double M);
double tinf1f(double M);
double tinf2f(double M);
double txf(double M);
double Mxf(double M);
double Lxf(double M);
double Mcgbf1a(double M, double t);                              /* real core mass on RGB */   
double Mcgbf1b(double M, double t);                              /* dummy core mass on RGB */
double Mcgbf2(double M, double L);
double Lgbf1(double M, double t);
double Lgbf2(double M, double Mc);
double Rgbf(double M, double L);

double Rhe1f(double M);                                          /* Core Helium Burning Phase functions */
double Mche1f(double M);
double Lhe1f(double M);
double the1f(double M);
double thef(double M);
double Mchef(double M, double t);
double Lminhef(double M);
double Lzahbf(double M, double Mc);
double Rzahbf(double M, double Mc);
double Rmhef(double M);
double Rmhefa(double M, double Lzahb);
double fblf(double M);
double taublf(double M);
void CHeBf(double M, double t, double *R, double *L, double *Mc);

double tbagbf(double M);                                         /* Asymptotic Giant Branch functions */
double Mcbagbf(double M);
double Lbagbf(double M);
double Rbagbf(double M);
double Ahef(void);
double tinf1ef(double M);
double tinf2ef(double M);
double txef(double M);
double Mceagbf1(double M, double t);
double Mceagbf2(double M, double L);
double Leagbf1(double M, double t);
double Leagbf2(double M, double Mcco);
double Mcduf(double M);
double Lduf(double M);
double tduf(double M);
double Ahhef(void);
double tinf1tf(double M);
double tinf2tf(double M);
double txtf(double M);
double Mctagbf1a(double M, double t);                            /* real core mass on TP AGB */
double Mctagbf1b(double M, double t);                            /* dummy core mass on TP AGB */
double Mctagbf2(double M, double L);
double Ltagbf1(double M, double t);
double Ltagbf2(double M, double Mc);
double lambdaf(double M);
double Ragbf(double M, double L);

double Lzhsf(double M);                                          /* Naked Helium Stars functions */
double Rzhsf(double M);
double thsmsf(double M);
double Lhsmsf(double M, double t);
double Rhsmsf(double M, double t);
double Lthsf(double M);
double Bhsf(void);
double Dhsf(double M);
double phsf(void);
double qhsf(void);
double tinf1hsf(double M);
double tinf2hsf(double M);
double txhsf(double M);
double Lxhsf(double M);
double Mxhsf(double M);
double thsevolf(double M, double Mc);
double Mchsgbf1(double M, double t);
double Mchsgbf2(double M, double L);
double Lhsgbf(double M, double Mc);
double Rhsgbf(double M, double L, double Lths, double Rzhs, int *K);
double Mcmaxf(double M);

double compactf1(double Mco);                                    /* Compact Objects functions */
double compactf2(double M, double Mco, int *fb);
double compactf3(double M, double Mco, int *fb);
double compactf4(double M, double Mco, int *fb);
double Lwdf(double M, double t, double K);                      
double Rwdf(double M);
double Rnsf(void);
double Rbhf(double M);
double Mcsnf(double Mcbagb, double Mcco0);                       
double Mupf(void);
double Mecf(void);

double get_T(double L, double R);                                /* Other functions */
double inter_line(double x1, double y1, double x2, double y2, double x);
          
double tnfI(double M, double m0, int K);                         /* functions for binary.c */

