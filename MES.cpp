#include<iostream>
#include<iomanip>
#include<cmath>
#include<math.h>
#include<vector>

#define SCH 3

using namespace std;

struct node {
    double x, y;
    short BC;
};

struct Surface {
    double N[SCH][4] = { 0.0 };
};

struct elem2d_4 {
#if SCH==2
    double pc_ksi[2] = { -1 / sqrt(3), 1 / sqrt(3) };
    double pc_eta[4] = { -1 / sqrt(3), 1 / sqrt(3), 1 / sqrt(3), -1 / sqrt(3) };
    double w[2] = { 1.0, 1.0 };
    double wsq[4] = { 1.0, 1.0, 1.0, 1.0 };
#elif SCH==3
    double pc_ksi[3] = { -sqrt(0.6), 0, sqrt(0.6) };
    double pc_eta[9] = { -sqrt(0.6), 0, sqrt(0.6), sqrt(0.6), 0, -sqrt(0.6), -sqrt(0.6), 0, sqrt(0.6) };
    double w[3] = { (double)5 / (double)9, (double)8 / (double)9, (double)5 / (double)9 };
    double w3[3] = { w[0] * w[0], w[0] * w[1], w[1] * w[1] };
    double wsq[9] = { w3[0], w3[1], w3[0], w3[1], w3[2], w3[1], w3[0], w3[1], w3[0] };
#endif // SCH
    double pelne[2] = { -1.0, 1.0 };
    double dNdKsi[SCH * SCH][4];
    double dNdEta[SCH * SCH][4];
    Surface sur[4];
    void values() {
        for (int i = 0; i < (SCH * SCH); i++) {
            dNdKsi[i][0] = -0.25 * (1 - pc_ksi[i / SCH]);
            dNdKsi[i][1] = 0.25 * (1 - pc_ksi[i / SCH]);
            dNdKsi[i][2] = 0.25 * (1 + pc_ksi[i / SCH]);
            dNdKsi[i][3] = -0.25 * (1 + pc_ksi[i / SCH]);
        }

        for (int i = 0; i < (SCH * SCH); i++) {
            dNdEta[i][0] = -0.25 * (1 - pc_eta[i]);
            dNdEta[i][1] = -0.25 * (1 + pc_eta[i]);
            dNdEta[i][2] = 0.25 * (1 + pc_eta[i]);
            dNdEta[i][3] = 0.25 * (1 - pc_eta[i]);
        }

        cout << "dN/dKsi" << endl;
        for (int i = 0; i < (SCH * SCH); i++) {
            for (int j = 0; j < 4; j++) {
                cout << dNdKsi[i][j] << " ";
            }
            cout << endl;
        }
        cout << "dN/dEta" << endl;
        for (int i = 0; i < (SCH * SCH); i++) {
            for (int j = 0; j < 4; j++) {
                cout << dNdEta[i][j] << " ";
            }
            cout << endl;
        }
        cout << endl;

    }

};

struct element {
    int ID[4];
    double H[4][4] = { 0.0 };
    double Hbc[4][4] = { 0.0 };
    double C[4][4] = { 0.0 };
    double P[4] = { 0.0 };
};

struct grid {
    double H = 0.1;
    double B = 0.1;
    int nH = 31;
    int nB = 31;
    int nN = nH * nB;
    double deltaX = B / (nB - 1);
    double deltaY = H / (nH - 1);
    int nE = (nH - 1) * (nB - 1);

    void print() {
        cout << "Height of the grid: " << H << endl;
        cout << "Width of the grid: " << B << endl;
        cout << "Nodes up: " << nH << endl;
        cout << "Nodes across: " << nB << endl;
        cout << "Number of elements: " << nE << endl;
        cout << "Delta X = " << deltaX << endl;
        cout << "Delta Y = " << deltaY << endl;
    }

    void calcNodes(node* wezel) {
        wezel[0].x = 0.0;
        wezel[0].y = 0.0;
        for (int i = 1; i <= nN; i++) {
            if (i % nH != 0) {
                wezel[i].y = wezel[i - 1].y + deltaY;
                wezel[i].x = wezel[i - (i % nH)].x;
            }
            else {
                wezel[i].x = wezel[i - nH].x + deltaX;
                wezel[i].y = 0.0;
            }
        }

        /*for(int j=0; j<nN; ++j){
            cout<<"Node No."<<j+1<<": "<<wezel[j].x<<", "<<wezel[j].y<<endl;
        }
        cout<<endl;*/
    }

    void calcElements(element* el) {
        el[0].ID[0] = 1;
        el[0].ID[1] = el[0].ID[0] + nH;
        el[0].ID[2] = el[0].ID[1] + 1;
        el[0].ID[3] = el[0].ID[0] + 1;
        for (int i = 1; i < nE; i++) {
            if (i % (nH - 1) != 0) {
                el[i].ID[0] = el[i - 1].ID[0] + 1;
                el[i].ID[1] = el[i].ID[0] + nH;
                el[i].ID[2] = el[i].ID[1] + 1;
                el[i].ID[3] = el[i].ID[0] + 1;
            }
            else {
                el[i].ID[0] = el[i - 1].ID[3] + 1;
                el[i].ID[1] = el[i - 1].ID[2] + 1;
                el[i].ID[2] = el[i].ID[1] + 1;
                el[i].ID[3] = el[i].ID[0] + 1;
            }
        }
        /*cout<<"Nodes' IDs of each element"<<endl;
        for(int i=0; i<nE; ++i){
            cout<<"Element "<<i+1<<": "<<el[i].ID[0]<<" "<<el[i].ID[1]<<" "<<el[i].ID[2]<<" "<<el[i].ID[3]<<endl;
        }
        cout<<endl;*/
    }

    void setBC(node* tab) {
        for (int j = 0; j < nH; ++j) {
            tab[j].BC = 1;
        }
        for (int j = nN-1; j >= nN - nH; --j) {
            tab[j].BC = 1;
        }
        for (int j = 0; j <= nN - nH; j += nB) {
            tab[j].BC = 1;
        }
        for (int j = nH-1; j < nN; j += nB) {
            tab[j].BC = 1;
        }
    }

    void checkBC(node* tab) {
        for (int i = 0; i < nN; i++) {
            if (tab[i].BC == 1) {
                cout << tab[i].x << " " << tab[i].y << endl;
            }
        }
        cout << endl;
    }

};

void zeroing(double tab[4][4]) {
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            tab[i][j] = 0;
        }
    }
}

vector<double> rozwUkladu(int nN, double** macierz, double* wektor) {

    vector<double> n(nN), x1(nN), x2(nN);
    vector<vector<double>> M(nN);

    for (int i = 0; i < nN; i++)
        M[i].resize(nN);

    for (int i = 0; i < nN; i++)
        n[i] = 1 / macierz[i][i];

    for (int i = 0; i < nN; i++)
        for (int j = 0; j < nN; j++)
            if (i == j)
                M[i][j] = 0;
            else
                M[i][j] = -(macierz[i][j] * n[i]);

    for (int k = 0; k < 100; k++) {
        for (int i = 0; i < nN; i++) {
            x2[i] = n[i] * wektor[i];
            for (int j = 0; j < nN; j++)
                x2[i] += M[i][j] * x1[j];
        }
        for (int i = 0; i < nN; i++)
            x1[i] = x2[i];
    }
    return x1;
}


int main() {
    grid g1;
    g1.print();
    node* wezel = new node[g1.nN];
    element* el = new element[g1.nE];
    g1.calcNodes(wezel);
    g1.calcElements(el);
    g1.setBC(wezel);
    g1.checkBC(wezel);
    elem2d_4 el2d;
    el2d.values();

    double cond = 25.0;
    double Tot = 1200.0;
    double c = 700.0;
    double ro = 7800.0;
    int time_step = 1;
    int time = 100;
    double t0 = 100.0;
    double alfa = 300.0;

    double J[2][2];
    double Jinv[2][2];
    double dX[4] = { 0.0 };
    double dY[4] = { 0.0 };
    double HX[4][4] = { 0.0 };
    double HY[4][4] = { 0.0 };
    double Hpc[4][4] = { 0.0 };
    double Cpc[SCH * SCH][4] = { 0.0 };
    double C[4][4] = { 0.0 };
    double detJh = (g1.deltaY) / 2.0;
    double detJb = (g1.deltaX) / 2.0;
    double Hbcpc[4][4]={0.0};
    double Hbcsur[4][4][4] = { 0.0 };
    double Ppc[4] = { 0.0 };
    double Psur[4][4] = { 0.0 };
    double* init_temp = new double[g1.nN];
    for (int i = 0; i < g1.nN; ++i) {
        init_temp[i] = t0;
    }
    for (int i = 0; i < g1.nE; ++i) {
        for (int j = 0; j < 4; ++j) {
            el[i].P[j] = 0.0;
        }
    }

    double** H_global = new double* [g1.nN];
    double** Hbc_global = new double* [g1.nN];
    double** C_global = new double* [g1.nN];
    double* P_global = new double[g1.nN];
    for (int i = 0; i < g1.nN; ++i) {
        H_global[i] = new double[g1.nN];
        Hbc_global[i] = new double[g1.nN];
        C_global[i] = new double[g1.nN];
    }

    for (int i = 0; i < g1.nN; ++i) {
        P_global[i] = 0.0;
        for (int j = 0; j < g1.nN; ++j) {
            H_global[i][j] = 0.0;
            Hbc_global[i][j] = 0.0;
            C_global[i][j] = 0.0;
        }
    }


    for (int i = 0; i < g1.nE; ++i) {
        for (int j = 0; j < (SCH * SCH); ++j) {
            J[0][0] = 0;
            J[0][1] = 0;
            J[1][0] = 0;
            J[1][1] = 0;
            for (int k = 0; k < 4; ++k) {
                J[0][0] += el2d.dNdKsi[j][k] * wezel[el[i].ID[k] - 1].x;
                J[0][1] += el2d.dNdKsi[j][k] * wezel[el[i].ID[k] - 1].y;
                J[1][0] += el2d.dNdEta[j][k] * wezel[el[i].ID[k] - 1].x;
                J[1][1] += el2d.dNdEta[j][k] * wezel[el[i].ID[k] - 1].y;
                Jinv[0][0] = J[1][1];
                Jinv[0][1] = -1 * J[1][0];
                Jinv[1][0] = -1 * J[0][1];
                Jinv[1][1] = J[0][0];
            }

            //cout << "Punkt calkowania nr " << j + 1 << " elementu " << i + 1 << endl;
            /*cout << fixed << setprecision(4) << J[0][0] << "\t" << J[0][1] << endl;
            cout<<fixed<<setprecision(4)<<J[1][0]<<"\t"<<J[1][1]<<endl;
            cout<<"Jakobian odwrotny"<<endl;
            cout<<fixed<<setprecision(4)<<Jinv[0][0]<<"\t"<<Jinv[0][1]<<endl;
            cout<<fixed<<setprecision(4)<<Jinv[1][0]<<"\t"<<Jinv[1][1]<<endl;*/
            double detJ = J[0][0] * J[1][1] - J[0][1] * J[1][0];
            //cout<<"det[J] = "<<detJ<<endl;
            double detInv;
            if (detJ == 0)
                detInv = 0;
            else
                detInv = 1 / detJ;
            //cout<<"1/det[J] = "<<detInv<<endl;
            Jinv[0][0] *= detInv;
            Jinv[0][1] *= detInv;
            Jinv[1][0] *= detInv;
            Jinv[1][1] *= detInv;
            /*cout << "Ostateczny Jakobian" << endl;
            cout<<Jinv[0][0]<<"\t"<<Jinv[0][1]<<endl;
            cout<<Jinv[1][0]<<"\t"<<Jinv[1][1]<<endl;*/

            for (int k = 0; k < 4; k++) {
                dX[k] = Jinv[0][0] * el2d.dNdKsi[j][k] + Jinv[0][1] * el2d.dNdEta[j][k];
                dY[k] = Jinv[1][0] * el2d.dNdKsi[j][k] + Jinv[1][1] * el2d.dNdEta[j][k];
            }
            /*cout << "dN/dX" << endl;
            for(int k = 0; k<4; k++)
                cout<<dX[k]<<" ";
            cout<<endl;

            cout<<"dN/dY"<<endl;
            for(int k = 0; k<4; k++)
                cout<<dY[k]<<" ";
            cout<<endl<<endl;*/

            for (int k = 0; k < 4; k++) {
                for (int l = 0; l < 4; ++l) {
                    HX[k][l] = dX[k] * dX[l];
                    HY[k][l] = dY[k] * dY[l];
                }
            }
            for (int k = 0; k < 4; ++k) {
                for (int l = 0; l < 4; ++l) {
                    Hpc[k][l] = el2d.wsq[j] * (cond * (HX[k][l] + HY[k][l])) * detJ;
                }
            }
            for (int k = 0; k < 4; ++k) {
                for (int l = 0; l < 4; ++l) {
                    el[i].H[k][l] += Hpc[k][l];
                }
            }

            Cpc[j][0] = 0.25 * (1 - el2d.pc_eta[j]) * (1 - el2d.pc_ksi[j / SCH]);
            Cpc[j][1] = 0.25 * (1 + el2d.pc_eta[j]) * (1 - el2d.pc_ksi[j / SCH]);
            Cpc[j][2] = 0.25 * (1 + el2d.pc_eta[j]) * (1 + el2d.pc_ksi[j / SCH]);
            Cpc[j][3] = 0.25 * (1 - el2d.pc_eta[j]) * (1 + el2d.pc_ksi[j / SCH]);

            for (int k = 0; k < 4; ++k) {
                for (int l = 0; l < 4; ++l) {
                    C[k][l] = c * ro * el2d.wsq[j] * Cpc[j][k] * Cpc[j][l] * detJ;
                }
            }
            for (int k = 0; k < 4; ++k) {
                for (int l = 0; l < 4; ++l) {
                    el[i].C[k][l] += C[k][l];
                }
            }


            /*cout << "HX" << endl;
            for(int k=0; k<4; ++k){
                for(int l=0; l<4; ++l){
                    cout<<HX[k][l]<<" ";
                }
                cout<<endl;
            }
            cout<<endl;
            cout<<"HY"<<endl;
            for(int k=0; k<4; ++k){
                for(int l=0; l<4; ++l){
                    cout<<HY[k][l]<<" ";
                }
                cout<<endl;
            }
            cout<<endl;
            cout<<"Hpc"<<endl;
            for(int k=0; k<4; ++k){
                for(int l=0; l<4; ++l){
                    cout<<Hpc[k][l]<<" ";
                }
                cout<<endl;
            }
            cout<<endl;
            cout<<"C"<<endl;
            for(int k=0; k<4; ++k){
                cout<<Cpc[j][k]<<" ";
                //cout<<endl;
            }
            cout<<endl<<endl;*/
        }

        /*cout << "Macierz H:" << endl;
        for(int j=0; j<4; j++){
            for(int k=0; k<4; k++){
                cout<<el[i].H[j][k]<<" ";
            }
            cout<<endl;
        }
        cout<<"Macierz C:"<<endl;
        for(int j=0; j<4; j++){
            for(int k=0; k<4; k++){
                cout<<el[i].C[j][k]<<" ";
            }
            cout<<endl;
        }*/
    }

    cout << "DetJh = " << detJh << endl;
    cout << "DetJb = " << detJb << endl;
    for(int i=0; i<g1.nE; ++i){
        //cout << "Element " << i + 1 << endl;
        for (int j = 0; j < 4; ++j) {
            for (int k = 0; k < 4; ++k) {
                for (int l = 0; l < 4; ++l) {
                    Hbcsur[j][k][l] = 0.0;
                }
            }
        }
        zeroing(el[i].Hbc);
        if (wezel[el[i].ID[0] - 1].BC == 1 && wezel[el[i].ID[1] - 1].BC == 1) {
            //cout << "Sciana dolna" << endl;
            for (int j = 0; j < SCH; ++j) {
                zeroing(Hbcpc);
                el2d.sur[0].N[j][0] = 0.25 * ((double)1 - el2d.pc_ksi[j]) * ((double)1 - el2d.pelne[0]);
                el2d.sur[0].N[j][1] = 0.25 * ((double)1 + el2d.pc_ksi[j]) * ((double)1 - el2d.pelne[0]);
                for (int k = 0; k < 4; ++k) {
                    Ppc[k] = el2d.w[j] * el2d.sur[0].N[j][k] * Tot;
                    for (int l = 0; l < 4; ++l) {
                        Hbcpc[k][l] = el2d.w[j] * el2d.sur[0].N[j][k] * el2d.sur[0].N[j][l];
                        Psur[0][k] += cond * Ppc[k] * detJb;
                    }
                }
                for (int k = 0; k < 4; ++k) {
                    for (int l = 0; l < 4; ++l) {
                        Hbcsur[0][k][l] += alfa * Hbcpc[k][l] * detJb;
                    }
                }
            }

            for (int j = 0; j < 4; ++j) {
                el[i].P[j] += Psur[0][j];
                for (int k = 0; k < 4; ++k) {
                    el[i].Hbc[j][k] += Hbcsur[0][j][k];
                }
            }
            /*for (int j = 0; j < 4; ++j) {
                for (int k = 0; k < 4; ++k) {
                    cout << Hbcsur[0][j][k] << " ";
                }
                cout << endl;
            }
            cout << endl;
            for (int j = 0; j < 4; ++j) {
                cout << Psur[0][j] << " ";
            }
            cout << endl;*/
        }

        zeroing(Psur);
        for (int j = 0; j < 4; ++j) {
            for (int k = 0; k < 4; ++k) {
                for (int l = 0; l < 4; ++l) {
                    Hbcsur[j][k][l] = 0.0;
                }
            }
        }
        for (int j = 0; j < 4; ++j) {
            Ppc[j] = 0.0;
        }

        if (wezel[el[i].ID[1] - 1].BC == 1 && wezel[el[i].ID[2] - 1].BC == 1) {
            //cout << "Sciana prawa" << endl;
            for (int j = 0; j < SCH; ++j) {
                zeroing(Hbcpc);
                el2d.sur[1].N[j][1] = 0.25 * (1 + el2d.pelne[1]) * (1 - el2d.pc_ksi[j]);
                el2d.sur[1].N[j][2] = 0.25 * (1 + el2d.pelne[1]) * (1 + el2d.pc_ksi[j]);
                for (int k = 0; k < 4; ++k) {
                    Ppc[k] = el2d.w[j] * el2d.sur[1].N[j][k] * Tot;
                    for (int l = 0; l < 4; ++l) {
                        Hbcpc[k][l] = el2d.w[j] * el2d.sur[1].N[j][k] * el2d.sur[1].N[j][l];
                        Psur[1][k] += cond * Ppc[k] * detJh;
                    }
                }
                for (int k = 0; k < 4; ++k) {
                    for (int l = 0; l < 4; ++l) {
                        Hbcsur[1][k][l] += alfa * Hbcpc[k][l] * detJh;
                    }
                }
            }
            
            for (int j = 0; j < 4; ++j) {
                el[i].P[j] += Psur[1][j];
                for (int k = 0; k < 4; ++k) {
                    el[i].Hbc[j][k] += Hbcsur[1][j][k];
                }
            }
            /*for (int j = 0; j < 4; ++j) {
                for (int k = 0; k < 4; ++k) {
                    cout << Hbcsur[1][j][k] << " ";
                }
                cout << endl;
            }
            for (int j = 0; j < 4; ++j) {
                cout << Psur[1][j] << " ";
            }
            cout << endl;*/
        }

        for (int j = 0; j < 4; ++j) {
            for (int k = 0; k < 4; ++k) {
                for (int l = 0; l < 4; ++l) {
                    Hbcsur[j][k][l] = 0.0;
                }
            }
        }
        for (int k = 0; k < 4; ++k) {
            Ppc[k] = 0.0;
        }

        zeroing(Psur);

        if (wezel[el[i].ID[2] - 1].BC == 1 && wezel[el[i].ID[3] - 1].BC == 1) {
            //cout << "Sciana gorna" << endl;
            for (int j = 0; j < SCH; ++j) {
                zeroing(Hbcpc);
                el2d.sur[2].N[j][2] = 0.25 * (1 + el2d.pc_ksi[j]) * (1 + el2d.pelne[1]);
                el2d.sur[2].N[j][3] = 0.25 * (1 - el2d.pc_ksi[j]) * (1 + el2d.pelne[1]);
                for (int k = 0; k < 4; ++k) {
                    Ppc[k] = el2d.w[j] * el2d.sur[2].N[j][k] * Tot;
                    for (int l = 0; l < 4; ++l) {
                        Hbcpc[k][l] = el2d.w[j] * el2d.sur[2].N[j][k] * el2d.sur[2].N[j][l];
                        Psur[2][k] += cond * Ppc[k] * detJb;
                    }
                   
                }
                for (int k = 0; k < 4; ++k) {
                    for (int l = 0; l < 4; ++l) {
                        Hbcsur[2][k][l] += alfa * Hbcpc[k][l] * detJb;
                    }
                }
            }
            for (int j = 0; j < 4; ++j) {
                el[i].P[j] += Psur[2][j];
                for (int k = 0; k < 4; ++k) {
                    el[i].Hbc[j][k] += Hbcsur[2][j][k];
                }
            }
            /*for (int j = 0; j < 4; ++j) {
                for (int k = 0; k < 4; ++k) {
                    cout << Hbcsur[2][j][k] << " ";
                }
                cout << endl;
            }
            for (int j = 0; j < 4; ++j) {
                cout << Psur[2][j] << " ";
            }
            cout << endl;*/
        }

        for (int j = 0; j < 4; ++j) {
            for (int k = 0; k < 4; ++k) {
                for (int l = 0; l < 4; ++l) {
                    Hbcsur[j][k][l] = 0.0;
                }
            }
        }
        for (int k = 0; k < 4; ++k) {
            Ppc[k] = 0.0;
        }

        zeroing(Psur);

        if (wezel[el[i].ID[3] - 1].BC == 1 && wezel[el[i].ID[0] - 1].BC == 1) {
            //cout << "Sciana lewa" << endl;
            for (int j = 0; j < SCH; ++j) {
                zeroing(Hbcpc);
                el2d.sur[3].N[j][3] = 0.25 * (1 - el2d.pelne[0]) * (1 + el2d.pc_ksi[j]);
                el2d.sur[3].N[j][0] = 0.25 * (1 - el2d.pelne[0]) * (1 - el2d.pc_ksi[j]);
                for (int k = 0; k < 4; ++k) {
                    Ppc[k] = el2d.w[j] * el2d.sur[3].N[j][k] * Tot;
                    for (int l = 0; l < 4; ++l) {
                        Hbcpc[k][l] = el2d.w[j] * el2d.sur[3].N[j][k] * el2d.sur[3].N[j][l];
                        Psur[3][k] += cond * Ppc[k] * detJh;
                    }
                    
                }
                for (int k = 0; k < 4; ++k) {
                    for (int l = 0; l < 4; ++l) {
                        Hbcsur[3][k][l] += alfa * Hbcpc[k][l] * detJh;
                    }
                }
            }
            
            
            for (int j = 0; j < 4; ++j) {
                el[i].P[j] += Psur[3][j];
                for (int k = 0; k < 4; ++k) {
                    el[i].Hbc[j][k] += Hbcsur[3][j][k];
                }
            }
            /*for (int j = 0; j < 4; ++j) {
                for (int k = 0; k < 4; ++k) {
                    cout << Hbcsur[3][j][k] << " ";
                }
                cout << endl;
            }*/
            /*cout << "Macierz Hbc elementu " << i + 1 << endl;
            for (int j = 0; j < 4; ++j) {
                for (int k = 0; k < 4; ++k) {
                    cout << el[i].Hbc[j][k] << " ";
                }
                cout << endl;
            }*/
            /*for (int j = 0; j < 4; ++j) {
                cout << Psur[3][j] << " ";
            }
            cout << endl;*/
        }
        //cout<<endl;
    }

    /*for (int i = 0; i < g1.nE; ++i) {
        cout<<"Macierz Hbc elementu "<<i+1<<endl;
        for(int j=0; j<4; ++j){
            for(int k=0; k<4; ++k){
                cout<<el[i].Hbc[j][k]<<" ";
            }
            cout<<endl;
        }
        cout<<endl;
        cout<<"Wektor P"<<endl;
        for(int j=0; j<4; ++j)
            cout<<el[i].P[j]<<" ";
        cout<<endl;
    }*/

    for (int i = 0; i < g1.nE; ++i) {
        for (int j = 0; j < 4; ++j) {
            for (int k = 0; k < 4; ++k) {
                H_global[el[i].ID[j] - 1][el[i].ID[k] - 1] = H_global[el[i].ID[j] - 1][el[i].ID[k] - 1] + el[i].H[j][k];
                Hbc_global[el[i].ID[j] - 1][el[i].ID[k] - 1] = Hbc_global[el[i].ID[j] - 1][el[i].ID[k] - 1] + el[i].Hbc[j][k];
                C_global[el[i].ID[j] - 1][el[i].ID[k] - 1] = C_global[el[i].ID[j] - 1][el[i].ID[k] - 1] + el[i].C[j][k];
            }
            P_global[el[i].ID[j] - 1] = P_global[el[i].ID[j] - 1] + el[i].P[j]*3;
        }
    }
    /*cout << "H global" << endl;
    for(int i=0; i<g1.nN; ++i){
        for(int j=0; j<g1.nN; ++j){
            cout<<H_global[i][j]<<" ";
        }
        cout<<endl;
    }*/
    /*cout << "Hbc global" << endl;
    for(int i=0; i<g1.nN; ++i){
        for(int j=0; j<g1.nN; ++j){
            cout<<Hbc_global[i][j]<<" ";
        }
        cout<<endl;
    }*/


    for (int j = 0; j < g1.nN; ++j) {
        for (int k = 0; k < g1.nN; ++k) {
            H_global[j][k] = H_global[j][k] + Hbc_global[j][k];
        }
    }

    /*for (int j = 0; j<g1.nN; ++j) {
        for(int k=0; k<g1.nN; ++k){
            HC_global[j][k] = H_global[j][k]+(C_global[j][k]/time_step);
        }
    }
    for(int j=0; j<g1.nN; ++j){
        for(int k=0; k<g1.nN; ++k){
            vecP[j] += (C_global[j][k]/time_step)*vecT[k];
        }
        vecP[j] += P_global[j];
    }*/



    /*cout << "H+Hbc global" << endl;
    for(int i=0; i<g1.nN; ++i){
        for(int j=0; j<g1.nN; ++j){
            cout<<H_global[i][j]<<" ";
        }
        cout<<endl;
    }
    cout<<endl;

    cout<<endl;
    cout<<"C global"<<endl;
    for(int i=0; i<g1.nN; ++i){
        for(int j=0; j<g1.nN; ++j){
            cout<<C_global[i][j]<<" ";
        }
        cout<<endl;
    }
    cout<<endl;*/
    /*cout << "P global" << endl;
    for(int i=0; i<g1.nN; ++i){
        cout<<P_global[i]<<" ";
    }
    cout<<endl;*/


    double** HC = new double* [g1.nN];
    double* vecP = new double[g1.nN];
    for (int i = 0; i < g1.nN; ++i) {
        HC[i] = new double[g1.nN];
    }
    for (int i = 0; i < g1.nN; ++i) {
        vecP[i] = 0.0;
    }
    for (int j = 0; j < g1.nN; ++j) {
        for (int k = 0; k < g1.nN; ++k) {
            HC[j][k] = H_global[j][k] + (C_global[j][k] / (double)time_step);
        }
    }

    for (int j = 0; j < g1.nN; ++j) {
        for (int k = 0; k < g1.nN; ++k) {
            vecP[j] += (C_global[j][k] / (double)time_step) * init_temp[k];
        }
        vecP[j] += P_global[j];
    }


    /*cout << "H=H+C/dt" << endl;
    for(int i=0; i<g1.nN; ++i){
        for(int j=0; j<g1.nN; ++j){
            cout<<HC[i][j]<<" ";
        }
        cout<<endl;
    }
    cout<<endl;
    cout<<"P=P+(C/dT)*temp"<<endl;
    for(int i=0; i<g1.nN; ++i){
        cout<<vecP[i]<<" ";
    }
    cout<<endl;*/

    double maksimum, minimum;

    vector<double> vecT(g1.nN);
    cout << "Czas: " << time_step << "s" << endl;
    vecT = rozwUkladu(g1.nN, HC, vecP);
    minimum = vecT[0];
    maksimum = vecT[0];
    for (int j = 1; j < g1.nN; ++j) {
        if (vecT[j] > maksimum)
            maksimum = vecT[j];
        else if (vecT[j] < minimum)
            minimum = vecT[j];
    }
    cout << "Tmin = " << minimum << "\tTmax = " << maksimum << endl;
    /*for(int j=0; j<g1.nN; ++j){
        cout<<"T["<<j+1<<"]: "<<vec[j]<<endl;
    }*/


    for (int i = 2 * time_step; i <= time; i += time_step) {
        vecT = rozwUkladu(g1.nN, HC, vecP);
        for (int j = 0; j < g1.nN; ++j) {
            vecP[j] = 0.0;
        }

        for (int j = 0; j < g1.nN; ++j) {
            for (int k = 0; k < g1.nN; ++k) {
                vecP[j] += (C_global[j][k] / (double)time_step) * vecT[k];
            }
            vecP[j] += P_global[j];
        }
        /*for(int i=0; i<g1.nN; ++i){
            cout<<vecP[i]<<" ";
        }
        cout<<endl;*/

        vecT = rozwUkladu(g1.nN, HC, vecP);
        minimum = vecT[0];
        maksimum = vecT[0];
        for (int j = 1; j < g1.nN; ++j) {
            if (vecT[j] > maksimum)
                maksimum = vecT[j];
            else if (vecT[j] < minimum)
                minimum = vecT[j];
        }

        cout << "Czas: " << i << "s" << endl;
        cout << "Tmin = " << minimum << "\tTmax = " << maksimum << endl;
        /*for(int j=0; j<g1.nN; ++j){
            cout<<"T["<<j+1<<"]: "<<vec[j]<<endl;
        }*/


    }

    delete[] P_global;
    delete[] vecP;
    for (int i = 0; i < g1.nN; ++i)
        delete[] HC[i];
    delete[] HC;

    for (int i = 0; i < g1.nN; ++i)
        delete[] H_global[i];
    delete[] H_global;

    for (int i = 0; i < g1.nN; ++i)
        delete[] Hbc_global[i];
    delete[] Hbc_global;

    for (int i = 0; i < g1.nN; ++i)
        delete[] C_global[i];
    delete[] C_global;

    return 0;
}
