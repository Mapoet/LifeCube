
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
typedef struct _Cells{
	int		idx;
	SPos	loc;
	int		iface;
	int     nface;
}Cells;
typedef struct _SuperCubes{
	float     width;
	SPos	  center;
	int       nseg;
	Colorsf	  facecolor[NDEM * 2];
	Cells	  cell[MAXSEG*MAXSEG*NDEM * 2];
}SuperCubes;
int  CompareCellLoc(const void*ce1, const void*ce2);
int  CompareCellIdx(const void*ce1, const void*ce2);
int  IsCellOnLayer(SuperCubes*cube, Cells*cell, int idem, int ilayer);
void UpdateCellFace(SuperCubes*cube, Cells* cell);
void FillCell(SuperCubes* cube);
void InitSuperCube(SuperCubes* cube, SPos center, double width);
void BrokenSuperCube(SuperCubes* cube);
void ShowSuperCube(SuperCubes* cube);
void Operating(SuperCubes* cube, int idem, int ilayer, int times);
void Roting(SuperCubes* cube, int idem, float phi);
SuperCubes* CreateSuperCubes();
#endif