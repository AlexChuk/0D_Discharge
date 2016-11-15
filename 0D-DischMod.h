# include <stdio.h>
# include <math.h>
# include <string.h>
# include <stdlib.h>

#ifndef MODEL_EEDF_H
#define MODEL_EEDF_H

# define NEmax 1000  //dots per energy range
# define Nmax 70  //components
# define CSmax 100  //max number of CS
# define NRmax 3000  //max number of chem reactions

const double
	pi = 3.141592653589,
	c = 2.997924562e+10,
	ma = 1.67e-24,//[�]
	me = 9.1e-28,//[�]
	e = 4.8e-10,//[���]
	E0 = 300,//E[�/��]=300E[���.���]
	Na = 6.022e+23,
	kb = 1.38e-16,//[erg/K]
	eV = 11605,//1 �� � ��������� �
	p0 = 1333.22, // ����-� �������� ���� --> ��� [���/��^3]
	cm_eV = 8065.5447;//����-� �������� cm-1 --> ��

extern int N,Nt,Nte,Nchem,Nedf,Ndots,NR;

extern double Ne[NEmax],Ni[Nmax],Mi[Nmax],LJi[Nmax][2],Roi[Nmax],Pgas,Tgas,Ngas,Hgas,Rogas;
extern double E,E_N,Nel,Ee,Te,Tv,Vdr,Muel,Jel,Qel,QE;
extern double dTgas,dTe,dNel;
extern double Rad;
extern double tau,dt,dte;
extern double Emax,dE,dEev;
extern double CXi[Nmax][2][8];
extern double Kel[CSmax];
extern double HCpSi[3][Nmax];

extern char Rtype[NRmax][10];
extern char Spec[Nmax][10],Spec_R[Nmax][10];

#endif
