#define c1 2.0    //variable cost


#define LAMDA 0.7425     //revenue sharing coefficient *****
#define Y0  40//40
#define K   2
#define KR  2//0.3//0.2// 8

#define END 100 //100//0~END iteration
#define L0  2
#define LM  70 // maximum life-time
#define PRI 20.0
double a1 = 60.0;
double b1 = 1.0;


#define PINTY 702 //maximum Point for ordering
#define PINTL 200 //maximum Point for life-span
#define PRICEP 50 //maximum Point for price
#define PINTHETA 50
double l[PINTL+1];

double p[PRICEP];

int POINTS;
int POINTL;
int PRICE;
int THETA;

double Q[PINTY];  //*****
double theta[PINTHETA];//****

#define MEAN0 2.0      //additive demand d(l)+epsilon, epsilon~ N(Mean0, SIGMA0)
#define SIGMA0 2.0       //
#define TRUNCT  0.977218197  //Phi(4)-Phi(-mu/sigma);

#define MEAN1 2.0      //multiplicative demand d(l)epsilon, epsilon~N(Mean1, SIGMA1)
#define SIGMA1 0.4    //

//#define Alpha1 5.0
#define Alpha_r 1// 0.1 // m(l)=Alpha_r(l-l0)^2

#define STEP 2 //2//0.5//0.2//follow step search
#define STEPIntg 0.05//theta's integral precision
//#define TT  2   //theta
#define NEGATIVE -99999.0

//double Vs[PINTY][PINTL+1][PRICEP];
double Vs[PINTY];
double Vr[PINTHETA][PINTL+1][PRICEP];
double VR_opt[PINTY][PINTHETA];
double z[PINTY][PINTHETA];
int Nl[PINTY][PINTHETA], Np[PINTY][PINTHETA];
int Nq;

//Another Revenue Sharing
double Vs_ER[PINTY];
double Vr_ER[PINTHETA][PINTL+1][PRICEP];
double Vc_ER[PINTY]; // Centralized Total Profit
double VR_ER_opt[PINTY][PINTHETA];
double z_ER[PINTY][PINTHETA];
int Nl_ER[PINTY][PINTHETA], Np_ER[PINTY][PINTHETA];
int Nq_ER;
int Nq_c_ER; // Centralized Total Profit
//Supplier Leads lifespan investment
double Vs_SL[PINTY][PINTL+1];
double Vr_SL[PINTHETA][PINTL+1][PRICEP];
double VR_SL_opt[PINTY][PINTL+1][PINTHETA];
double z_SL[PINTY][PINTHETA];
int  Np_SL[PINTY][PINTL+1][PINTHETA];
int Nl_SL, Nq_SL;

int Ny_c,Ny_cb,Nl_cb,Nm_db,Nw_dr,Nm_theta;

double VRtheta[PINTL+1][PINTY];
/*********************function declaration ******************************/
double max(double x,double y);
double min(double x,double y);
double N(double x,int flg);
double MaintFee(double l);//maintainance cost function in life-span l
double LifeSpan(double maintainance);

double Demand(double price,double l);//d(l)=y_0*l*p^{-K}
double Profit(double l, double price, double q, double epsilon);//计算R(Q,l,theta) //min(D(l,p),theta*Q)


//void opt_cb(double v[PINTY][PINTL+1], int &I, int &J);
void opt_R_pl(double v[PINTL+1][PRICEP], int &L, int &P);
void opt_RSL_pl(double v[PRICEP], int &P);
void opt_S_q(double v[PINTY], int &nq);
void opt_S_ql(double v[PINTY][PINTL+1], int &L, int &q);
//void opt_r_db(int k, double vv[PINTL][PINTY], int &Ny, int &Nl, double &v_opt);
void opt_r(double v[PINTY], int &N_y,double &v_opt);//
double Romberg0(double q,double l, double price); //Romberg0(Q[i],l[j], p[k])
void output_db();
void save_json();
