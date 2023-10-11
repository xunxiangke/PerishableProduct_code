// Fish_db.cpp : Defines the entry point for the console application.
//Newsvendor Model for perishable products based on price-only contracts
//consider both parties invest for given wholesale price w: supplier--M, retailer--m(l)

#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include "string.h"
#include "stdlib.h"
#include "constants.h"
#include "json/json.h"
#include <ctime>

#include <unistd.h>

using namespace std;


double a1 = 60.0;
double b1 = 1.0;

int main(int argc, char *argv[]) {
    POINTS = (int) (END / STEP + .1); // ordering quantity
    POINTL = (int) ((LM - L0) / STEP + .1); //lifespan
    PRICE = (int) ((PRI - c1) / STEP + .1);
    THETA = (int) ((1.0 - 0.5) / STEPIntg + .1);


    for (int i = 0; i <= POINTS; i++) {
        if (i == 490)
            i = i;

        Q[i] = i * STEP; //initial value

    }
    for (int j = 0; j <= POINTL; j++) {
        l[j] = L0 + j * STEP;
    }
    for (int k = 0; k <= PRICE; k++) {
        p[k] = c1 + k * STEP;
    }

    for (int k = 0; k <= THETA; k++) {
        theta[k] = 0.5 + k * STEPIntg;
    }

    double tmp_integral = 0.0;

    //**Two Models for Basic and Revenue sharing scheme
    for (int i = 0; i <= POINTS; i++)//Q
    {
        int tmp_theta = 0;
        double tmp_vs = 0.0, tmp_vs_ER = 0.0;
        do {

            for (int j = 0; j <= POINTL; j++)//l
                for (int k = 0; k <= PRICE; k++)//p price
                {
                    tmp_integral = Romberg0(theta[tmp_theta] * Q[i], l[j], p[k]);//integral on D
                    //	tmp_integral = Romberg0(tmp_theta*STEPIntg*Q[i],l[j], p[k]);//integral on D

                    //Vs[i]  = -c*Q[i]-MaintFee(l[j]) + tmp_integral;
                    //Vs[i]  =(1-LAMDA )*Vs[i];

                    Vr[tmp_theta][j][k] = -MaintFee(l[j]) + tmp_integral;
                    Vr_ER[tmp_theta][j][k] = LAMDA * tmp_integral - MaintFee(l[j]);

                }//end j k - l p

            opt_R_pl(Vr[tmp_theta], Nl[tmp_theta], Np[tmp_theta]);
            opt_R_pl(Vr_ER[tmp_theta], Nl_ER[tmp_theta], Np_ER[tmp_theta]);

            z[i][tmp_theta] = Y0 * l[Nl[tmp_theta]] * pow(p[Np[tmp_theta]], -K) / (theta[tmp_theta] * Q[i]);
            z_ER[i][tmp_theta] = Y0 * l[Nl_ER[tmp_theta]] * pow(p[Np_ER[tmp_theta]], -K) / (theta[tmp_theta] * Q[i]);

            VR_opt[i][tmp_theta] = Vr[tmp_theta][Nl[tmp_theta]][Np[tmp_theta]];
            VR_ER_opt[i][tmp_theta] = Vr_ER[tmp_theta][Nl_ER[tmp_theta]][Np_ER[tmp_theta]];

            tmp_vs += VR_opt[i][tmp_theta] * STEP *
                      2.0;  //at present, VR_opt[i][tmp_theta]=E_D[R(Q,l,theta)-Kr-ar(l-l0)^2], assume uniform of theta, U[0,1],density is 1
            tmp_vs_ER += (VR_ER_opt[i][tmp_theta] + MaintFee(l[Nl_ER[tmp_theta]])) * STEP * 2.0 / LAMDA;

            VR_opt[i][tmp_theta] = LAMDA * VR_opt[i][tmp_theta];

            cout << "\nQ=" << Q[i] << "\ttheta=" << tmp_theta * STEP << "\topt_l=" << Nl[tmp_theta] << "\topt_p="
                 << Np[tmp_theta] << "\tVr=" << VR_opt[i][tmp_theta] << "\topt_l_ER=" << Nl_ER[tmp_theta]
                 << "\topt_p_ER=" << Np_ER[tmp_theta] << "\tVr_ER=" << VR_ER_opt[i][tmp_theta] << endl;
            tmp_theta++;
        } while (tmp_theta <= THETA);

        Vs[i] = (1 - LAMDA) * tmp_vs - c1 * Q[i];

        Vs_ER[i] = (1 - LAMDA) * tmp_vs_ER - c1 * Q[i];

    }//end i -Q

    opt_S_q(Vs, Nq);
    opt_S_q(Vs_ER, Nq_ER);

    cout << "\n\nopt_Q=" << Q[Nq] << "\tVs=" << Vs[Nq] << "opt_Q_ER=" << Q[Nq_ER] << "\tVs_ER=" << Vs_ER[Nq_ER] << endl;


    //Model-Supplier Leads lifespan investment
    for (int i = 0; i <= POINTS; i++)//Q
        for (int j = 0; j <= POINTL; j++)//l
        {
            int tmp_theta = 0;
            double tmp_vs_SL = 0.0;
            do {        //given tmp_theta, retailer has random yield quantity--theta*Q[i]
                for (int k = 0; k <= PRICE; k++)//p price
                {
                    tmp_integral = Romberg0(theta[tmp_theta] * Q[i], l[j], p[k]);//integral on D
                    Vr_SL[tmp_theta][j][k] = -MaintFee(l[j]) + tmp_integral;

                }//end k-p
                opt_RSL_pl(Vr_SL[tmp_theta][j], Np_SL[i][j][tmp_theta]);

                VR_SL_opt[i][j][tmp_theta] = Vr_SL[tmp_theta][j][Np_SL[i][j][tmp_theta]]; //VR_opt[i][tmp_theta] = LAMDA * VR_opt[i][tmp_theta];

                tmp_vs_SL += VR_SL_opt[i][j][tmp_theta] * STEP * 2.0;

                VR_SL_opt[i][j][tmp_theta] = LAMDA * VR_SL_opt[i][j][tmp_theta];
                tmp_theta++;
            } while (tmp_theta <= THETA);

            Vs_SL[i][j] = (1 - LAMDA) * tmp_vs_SL - c1 * Q[i];
        }//end i -Q j-l

    opt_S_ql(Vs_SL, Nl_SL, Nq_SL);
    cout << "\n\n" << "Opt_L_SL=" << l[Nl_SL] << "Opt_Q_SL=" << Q[Nq_SL] << "\tVs_SL=" << Vs_SL[Nq_SL][Nl_SL] << endl;

    output_db();
//    save_json();
    return 0;
}


//int main(){
//    cout<<Romberg0(2,1,3);
//}

//////////////////////// Subfunctions ///////////////////////////////////
double max(double x, double y)  //maximum of two double numbers.
{
    double maximum;
    maximum = x > y ? x : y;
    return (maximum);
}

double min(double x, double y)  //minimum of two double numbers.
{
    double minimum;
    minimum = x < y ? x : y;
    return (minimum);
}

double N(double x, int flg)  //normal density function
{
    double y, tmp;
    if (flg == 0) {
        tmp = SIGMA0 * pow(6.28, 0.5);
        y = (1 / tmp) * exp(-(x - MEAN0) * (x - MEAN0) / (2 * SIGMA0 * SIGMA0));
        y = y / TRUNCT;             //truncated normal
    } else {
        tmp = SIGMA1 * pow(6.28, 0.5);
        y = (1 / tmp) * exp(-(x - MEAN1) * (x - MEAN1) / (2 * SIGMA1 * SIGMA1));
    }
    return (y);
}

double MaintFee(double l) //maintainance cost function in life-span l
{
    if (l == L0) return (Alpha_r * (l - L0) * (l - L0));
    else return (KR + Alpha_r * (l - L0) * (l - L0));
}


double Demand(double price, double l)//d(l)=y_0*l*p^{-K}
{
    double y;
    y = Y0 * l * pow(price, -K);
    return (y);
}

double Profit(double l, double price, double q, double epsilon)//计算R(Q,l,theta) //min(D(l,p),theta*Q)
{
    double tmp;

    tmp = Demand(price, l);//*epsilon;

    tmp = price * min(tmp, q);
    return (tmp * N(epsilon, 1));
}

void opt_R_pl(double v[PINTL + 1][PRICEP], int &L, int &P) { //find  global v_cb(x)=min_{y,l}v[][]
    double value;
    int nl, np;

    value = v[0][0];
    nl = 0;
    np = 0;
    for (int i = 0; i <= POINTL; i++) //i - l
        for (int j = 0; j <= PRICE; j++)//j - price
        {
            //tmp=value;
            if (v[i][j] > value) {
                nl = i;
                np = j;
                value = v[i][j];
            }
        }
    P = np;
    L = nl;
}

//opt_R_pl(Vr_SL[tmp_theta][j], Np_SL[tmp_theta]);
void opt_RSL_pl(double v[PRICEP], int &P) { //find  global v_sl(x)=min_{p}v[][]
    double value;
    int nl, np;

    value = v[0];
    nl = 0;
    np = 0;
    for (int j = 0; j <= PRICE; j++)//j - price
    {
        //tmp=value;
        if (v[j] > value) {
            np = j;
            value = v[j];
        }
    }
    P = np;
}


void opt_S_q(double v[PINTY], int &nq) {//find  global v_s=min_{q}v[]

    double value = v[0];
    int q = 0;
    for (int i = 0; i <= POINTS; i++) {
        if (v[i] > value) {
            q = i;
            value = v[i];
        }
    }
    nq = q;
}


void opt_S_ql(double v[PINTY][PINTL + 1], int &L, int &q) { //find  global v_cb(x)=min_{q,l}v[][]
    double value;
    int nl, nq;

    value = v[0][0];
    nl = 0;
    nq = 0;
    for (int i = 0; i <= POINTS; i++) //i -q
        for (int j = 0; j <= POINTL; j++)//j - l
        {
            //tmp=value;
            if (v[i][j] > value) {
                nq = i;
                nl = j;
                value = v[i][j];
            }
        }
    q = nq;
    L = nl;
}

void opt_r(double v[PINTY], int &N_y, double &v_opt) {
    double value = v[0];
    int ny = 0;
    for (int i = 0; i <= POINTS; i++) {
        if (v[i] > value) //+STEP/100.0   +STEP/5.0
        {
            ny = i;
            value = v[i];
        }
    }
    N_y = ny;
    v_opt = value;
}


/********************************************************************
* 用龙贝格法计算函数f(x)从a到b的积分值
* 输入: f--函数f(x)的指针
*       a--积分下限
*       b--积分上限
*       eps--计算精度
*       max_it--最大迭代次数
* 输出: 返回值为f(x)从a?          //a = MEAN0-4*SIGMA0; b = MEAN0+4*SIGMA0;
* 迭代上限为64次?
*******************************************************************/

//double Romberg0(double q,double l, double price)  //Romberg0(Q[i],l[j], p[k])
//{
//    double eps = 0.0000001;
//    double	a,  b;
//    ///////////////////////////////////////////
//    a = MEAN1-4*SIGMA1; b = MEAN1+4*SIGMA1;
//    //////////////////////////////////////////H(double x,double d, int flg)
//    double T[64];
//    double hh=b-a;//, a=dt, b=ut;//a=x-ut, b= x-dt;
//    int n=1;
//    int tmp; // ++储存循环次数
//    //Profit(double l, double price, double q, double epsilon)
//
//    T[0]=hh*(Profit(l,price,q,a)/4.0+Profit(l,price,q,b)/4.0+Profit(l,price,q,a+hh/2.0)/2.0);//??????
//    T[1]=hh*(Profit(l,price,q,b)/6.0+Profit(l,price,q,a)/6.0+Profit(l,price,q,a+hh/2.0)/1.5);//?????
//
//    for(int i=2;fabs(T[i-2]-T[i-1])>eps;++i)
//    {
//        //计算递推梯形值,base
//        hh/=2.0;
//        n<<=1;//分数点
//        int j;
//        double base=T[0]/hh;
//        double y=a+hh/2.0;//(ut+dt)/2.0;//x - (ut+dt)/2.0,
//        for(j=0;j<n;++j)
//        {
//            //if(flg==0)
//            base+=Profit(l,price,q,y);
//            //else base+=Demand(x,flg)*N(y,flg);
//            y+=hh;
//        }
//        base=base*hh/2.0;
//        //计算外推加速值,T[0]->T[i]
//        double k1=4.0/3.0,k2=1.0/3.0;
//        for(j=0;j<i;++j)
//        {
//            double hand=k1*base-k2*T[j];
//            T[j]=base;
//            base=hand;
//            k2=k2/(4.0*k1-k2);
//            k1=k2+1.0;
//        }
//        //tmp = i;
//        T[i]=base;
//    }
//    return T[3];
//}


double Romberg0(double q,double l, double price)
{
    size_t max_steps = 100000;
    double acc = 0.0000001;
    double a,b;
    a = MEAN1-4*SIGMA1; b = MEAN1+4*SIGMA1;
    double R1[max_steps], R2[max_steps]; // buffers
    double *Rp = &R1[0], *Rc = &R2[0]; // Rp is previous row, Rc is current row
    double h = b-a; //step size
    Rp[0] = (Profit(l,price,q,a) + Profit(l,price,q,b)) * h * 0.5; // first trapezoidal step

    // print_row(0, Rp);

    for (size_t i = 1; i < max_steps; ++i) {
        h /= 2.;
        double c = 0;
        size_t ep = 1 << (i-1); //2^(n-1)
        for (size_t j = 1; j <= ep; ++j) {
            c += Profit(l,price,q,a + (2 * j - 1) * h);
        }
        Rc[0] = h*c + .5*Rp[0]; // R(i,0)

        for (size_t j = 1; j <= i; ++j) {
            double n_k = pow(4, j);
            Rc[j] = (n_k*Rc[j-1] - Rp[j-1]) / (n_k-1); // compute R(i,j)
        }

        // Print ith row of R, R[i,i] is the best estimate so far
        // print_row(i, Rc);

        if (i > 1 && fabs(Rp[i-1]-Rc[i]) < acc) {
            return Rc[i];
        }

        // swap Rn and Rc as we only need the last row
        double *rt = Rp;
        Rp = Rc;
        Rc = rt;
    }
    return Rp[max_steps-1]; // return our best guess
}

//////////////////////////////
void output_db() {
    //double value0;
    int i = 0;

    int timestamp  = std::time(0);  // Using timestamp to name the json file
    ostringstream tmp_str;
    tmp_str<<timestamp<<".ttmm";
    string file_name = tmp_str.str();

    ofstream g(file_name);

    //value0 = Alpha_r*(l[Nm_db]*l[Nm_db]-L0*L0);
    g << "Q=\t" << Q[Nq] << "\tVs=" << Vs[Nq] << endl << endl;
    g << "Q_ER=\t" << Q[Nq_ER] << "\tVs_ER=" << Vs_ER[Nq_ER] << endl << endl;
    g
            << "theta_i\t\tOpt_l\tOpt_p\tOpt_z\t\tVR_opt\t\tOpt_l_ER\tOpt_p_ER\tOpt_z_ER\t\tVR_ER_opt\t\tOpt_p_SL\t\tVR_SL_opt"
            << endl;
    do//for(int i =0;i<=THETA;i++)
    {
        g << theta[i] << "\t\t" << l[Nl[i]] << "\t" << p[Np[i]] << "\t" << z[Nq][i] << "\t\t" << VR_opt[Nq][i]
          << "\t\t" << l[Nl_ER[i]] << "\t" << p[Np_ER[i]] << "\t" << z[Nq_ER][i] << "\t\t" << VR_ER_opt[Nq_ER][i]
          << "\t\t" << p[Np_SL[Nq_SL][Nl_SL][i]] << "\t\t" << VR_SL_opt[Nq_SL][Nl_SL][i] << endl;
        i++;
    } while (i <= THETA);

    g << "\n\nOpt_l_SL=" << l[Nl_SL] << "\tOpt_Q_SL=" << Q[Nq_SL] << "\t\tVS_SL=" << Vs_SL[Nq_SL][Nl_SL] << endl<< endl;

    g << "Lamda\t" << LAMDA << "\n" << "l0\t" << L0 << "\n" << "Kr\t" << KR << "\n" << "Y0\t" << Y0 << "\n" << "K\t" << K
      << "\n" << "c1\t" << c1 << "\n" << "MEAN0\t" << MEAN0 << "\n" << "SIGMA0\t" << SIGMA0 << "\n" << "MEAN1\t" << MEAN1
      << "\n" << "SIGMA1\t" << SIGMA1 << "\n" << "Alpha_r\t" << Alpha_r << "\n" << "Alpha1\t" << a1 << "\n" << "Beta1\t"
      << b1 << "\n" << "STEP\t" << STEP << "\n" << "STEPIntg\t" << STEPIntg << "\n" << "TT\t" << THETA << "\n" << "END\t"
      << END << "\n" << "LM\t" << LM << "\n" << "PRI\t" << PRI << "\n" << "POINTS\t" << POINTS << "\n" << "POINTL\t"
      << POINTL << "\n" << "PRICE\t" << PRICE << "\n" << "PINTY\t" << PINTY << "\n" << "PINTL\t" << PINTL << "\n"
      << "PRICEP\t" << PRICEP << "\n" << "PINTHETA\t" << PINTHETA << "\n" << "LAMDA\t" << LAMDA << "\n" << "L0\t" << L0
      << "\n" << "KR\t" << KR << "\n" << "Y0\t" << Y0 << "\n" << "K\t" << K << "\n" << "c1\t" << c1 << "\n" << "MEAN0\t"
      << MEAN0 << "\n" << "SIGMA0\t" << SIGMA0 << "\n" << "MEAN1\t" << MEAN1 << "\n" << "SIGMA1\t" << SIGMA1 << "\n"
      << "Alpha_r\t" << Alpha_r << "\n" << "Alpha1\t" << a1 << "\n" << "Beta1\t" << b1 << "\n" << "STEP\t" << STEP
      << endl;

    g.close();
}

void save_json() {
    int i = 0;
    stringstream B_str;  // Basic model detail string
    B_str<<"theta_i\tOpt_l\tOpt_p\tOpt_z\tVR_opt"<<endl;

    stringstream ER_str; // ER model detail string
    ER_str<<"theta_i\tOpt_l_ER\tOpt_p_ER\tOpt_z_ER\tVR_ER_opt"<<endl;

    stringstream SL_str; // Supplier lead model detail string
    SL_str<<"theta_i\tOpt_p_SL\tVR_SL_opt"<<endl;

    while(i <= THETA){
        B_str << theta[i] << "\t" << l[Nl[i]] << "\t" << p[Np[i]] << "\t" << z[Nq][i] << "\t" <<VR_opt[Nq][i] << endl;
        ER_str << theta[i] << l[Nl_ER[i]] << "\t" << p[Np_ER[i]] << "\t" << z[Nq_ER][i] << "\t" << VR_ER_opt[Nq_ER][i]
               << endl;
        SL_str << theta[i] << "\t" << p[Np_SL[Nq_SL][Nl_SL][i]] << "\t" << VR_SL_opt[Nq_SL][Nl_SL][i] << endl;
        i++;
    }

    string B_str_res(B_str.str());
    string ER_str_res(ER_str.str());
    string SL_str_res(SL_str.str());


    Json::Value root;

    // Baisc model
    Json::Value basic;
    basic["opt_q"] = Q[Nq];
    basic["opt_Vs"] = Vs[Nq];
    basic["retailer"] = B_str_res;
    root["basic"] = Json::Value(basic);

    // Revenue sharing model(ER)
    Json::Value ER;
    ER["opt_q"] = Q[Nq_ER];
    ER["opt_Vs"] = Vs_ER[Nq_ER];
    ER["retailer"] = ER_str_res;
    root["ER"] = Json::Value(ER);

    // Supplier lead lifespan investment model(SL)
    Json::Value SL;
    SL["opt_q"] = Q[Nq_SL];
    SL["opt_l"] = l[Nl_SL];
    SL["opt_Vs"] = Vs_SL[Nq_SL][Nl_SL];
    SL["retailer"] = SL_str_res;
    root["SL"] = Json::Value(SL);

    // Ootput without indent
    cout << "StyledWriter:" << endl;
    Json::StyledWriter sw;
    cout << sw.write(root) << endl << endl;

    // Output json to file
    int timestamp  = std::time(0);  // Using timestamp to name the json file
    ostringstream tmp_str;
    tmp_str<<timestamp<<".json";
    string file_name = tmp_str.str();

    fstream os;
    os.open(file_name, std::ios::out | std::ios::app);
    if (!os.is_open())
        cout << "error：can not find or create the file which named \" demo.json\"." << endl;
    os << sw.write(root);
    os.close();
}