
#ifndef __SUPERCUBE_H__
#define __SUPERCUBE_H__
#include <stdlib.h>
#include <stdio.h>

#define NDEM    3
#define MAXSEG  50
//#define NSEG    17
#define NCOLOR  4
typedef float  Colorsf[NCOLOR];
typedef float  SPos[NDEM];
typedef struct _SuperUniverse{
	float     width;
	SPos	  center;
	Colorsf   cls[15];
	int       nseg;
	int      * life,*prelife;
	float    * power, *prepow;
}SuperUniverse;
typedef struct _SuperInit{
	float width;
	int   nseg;
	float nu, sigma;
}SuperInit;
typedef struct _SuperControl{
	float pow_start, pow_end;
	int effect_start, effect_end;
	float move_limit;
	float life_period, life_sigma;
	float(*Effect)(SuperUniverse*pUni, int idx, struct _SuperControl*pCtr);
	void(*Move)(SuperUniverse*pUni, int idx, struct _SuperControl*pCtr, SPos move);
}SuperControl;
SuperUniverse* CreateSuperUniverse(SuperInit init, SuperControl*pCtr);
void ShowSuperUniverse(SuperUniverse*pUni, SuperControl*pCtr);
void NextSuperUniverse(SuperUniverse*pUni, SuperControl*pCtr);
void FreeSuperUniverse(SuperUniverse*pUni);
float Effect(SuperUniverse*pUni, int idx, SuperControl*pCtr);
void Move(SuperUniverse*pUni, int idx, SuperControl*pCtr, SPos move);
void DrawLine(Colorsf cls, SPos move);
void DrawCube(Colorsf cls, float hw);
#endif