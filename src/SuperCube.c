
#include "SuperCube.h"
#include "base.h"
#include <time.h>
int  CompareCellLoc(const void*ce1, const void*ce2){
	Cells*c1 = (Cells*)ce1;
	Cells*c2 = (Cells*)ce2;
	for (int idem = 0; idem < NDEM; idem++)
		if (fabs(c1->loc[idem] - c2->loc[idem])>1e-6){
			if (c1->loc[idem] > c2->loc[idem])
				return 1;
			else
				return -1;
		}
	return 0;
}
int  CompareCellIdx(const void*ce1, const void*ce2){
	Cells*c1 = (Cells*)ce1;
	Cells*c2 = (Cells*)ce2;
	if (c1->idx > c2->idx)return 1;
	if (c1->idx < c2->idx)return -1;
	return 0;
}
int  IsCellOnLayer(SuperCubes*cube,Cells*cell,int idem,int ilayer){
	if (ilayer == 0)
		return 0;
	if (ilayer < 0)
		ilayer = NSEG + ilayer + 1;
	double dw = cube-> width / NSEG;
	if (ilayer == 1)
	{
		if (fabs(cell->loc[idem] + dw*NSEG / 2.) < dw)
			return 1;
		else
			return 0;
	}
	if (ilayer == NSEG)
	{
		if (fabs(cell->loc[idem] - dw*NSEG / 2.) < dw)
			return 1;
		else
			return 0;
	}
	if (fabs(cell->loc[idem] - dw*(ilayer - NSEG / 2. - 0.5)) < dw/4)
		return 1;
	else
		return 0;
	return 0;

}
void UpdateCellFace(SuperCubes*cube, Cells* cell){
	double dw = cube->width / NSEG;
	for (int idem = 0; idem < NDEM; idem++)
	{
		if (fabs(cell->loc[idem] + dw*NSEG / 2.) < dw / 2)
		{
			cell->iface = idem;
		}
		if (fabs(cell->loc[idem] - dw*NSEG / 2.) < dw / 2)
		{
			cell->iface = idem+NDEM;
		}
	}
}
void FillCell(SuperCubes* cube)
{
	int icell = 0;
	double dw = cube->width /NSEG;
	for (int iface = 0; iface < NDEM * 2; iface++)
	for (int idx1 = 0; idx1 < NSEG; idx1++)
	for (int idx2 = 0; idx2 < NSEG; idx2++)
	{
		if (iface >= NDEM)
			cube->cell[icell].loc[(iface%NDEM) % NDEM] = cube->center[(iface%NDEM) % NDEM] + dw*NSEG / 2.;
		else
			cube->cell[icell].loc[(iface%NDEM) % NDEM] = cube->center[(iface%NDEM) % NDEM] - dw*NSEG / 2.;
		cube->cell[icell].loc[(iface%NDEM + 1) % NDEM] = cube->center[(iface%NDEM + 1) % NDEM] + dw*(idx1 - NSEG / 2. + 0.5);
		cube->cell[icell].loc[(iface%NDEM + NDEM - 1) % NDEM] = cube->center[(iface%NDEM + NDEM - 1) % NDEM] + dw*(idx2 - NSEG / 2. + 0.5);
		cube->cell[icell].iface = iface;
		cube->cell[icell].nface = iface;
		cube->cell[icell].idx = icell;
		icell++;
	}
}
void InitSuperCube(SuperCubes* cube, SPos center, double width){
	srand(time(NULL));
	for (int i = 0; i < NDEM * 2;i++)
	for (int j = 0; j < 4; j++)
		cube->facecolor[i][j] = rand()%255 /255.;
	//cube->facecolor[0][0] = 1;//White
	//cube->facecolor[0][1] = 1;
	//cube->facecolor[0][2] = 1;
	//cube->facecolor[3][0] = 1;//Yellow
	//cube->facecolor[3][1] = 1;
	//cube->facecolor[3][2] = 0;
	//cube->facecolor[1][0] = 0;//Blue
	//cube->facecolor[1][1] = 0;
	//cube->facecolor[1][2] = 1;
	//cube->facecolor[4][0] = 0;//Green
	//cube->facecolor[4][1] = 1;
	//cube->facecolor[4][2] = 0;
	//cube->facecolor[2][0] = 1;//Red
	//cube->facecolor[2][1] = 0;
	//cube->facecolor[2][2] = 0;
	//cube->facecolor[5][0] = 1;//Or
	//cube->facecolor[5][1] = 0.32;
	//cube->facecolor[5][2] = 0;
	for (int i = 0; i < NDEM;i++)cube->center[i] = center[i];
	cube->width = width;
	FillCell(cube);
	BrokenSuperCube(cube);
}
void ShowSuperCube(SuperCubes* cube){
	qsort(cube->cell, NSEG*NSEG*NDEM * 2, sizeof(Cells), CompareCellLoc);
/*	for (int icell = 0; icell < NSEG*NSEG*NDEM * 2; icell++)
		fprintf(stdout, "Cell %3.3d:loc[%4.1f,%4.1f,%4.1f],color(%4.2f,%4.2f,%4.2f,%4.2f),idx %2.2d,face %2.2d.\n",
		icell, cube->cell[icell].loc[0], cube->cell[icell].loc[1], cube->cell[icell].loc[2],
		cube->facecolor[cube->cell[icell].iface][0], cube->facecolor[cube->cell[icell].iface][1], 
		cube->facecolor[cube->cell[icell].iface][2], cube->facecolor[cube->cell[icell].iface][3], 
		cube->cell[icell].idx,cube->cell[icell].iface);
		*/
	//int c = 0;
	double s[NDEM],e[NDEM],d;
	double dw = cube->width / NSEG;
	float v3f[NDEM];
	int idx[4][2] = { {0,0}, {0,1}, {1,1}, {1,0} };
	for (int iface = 0; iface < NDEM * 2; iface++)
	{
		////c = 0;
		//fprintf(stdout, "the %2.2d face with color(%4.2f,%4.2f,%4.2f,%4.2f):\n", iface,
		//	cube->facecolor[iface][0], cube->facecolor[iface][1],
		//	cube->facecolor[iface][2], cube->facecolor[iface][3]);
		//glDepthMask(GL_FALSE);
		for (int icell = 0; icell < NSEG*NSEG*NDEM * 2; icell++)
		if (cube->cell[icell].iface == iface)
		{
			d = cube->cell[icell].loc[iface%NDEM]; 
			s[0] = cube->cell[icell].loc[(iface%NDEM + 1) % NDEM] - dw / 2;
			s[1] = cube->cell[icell].loc[(iface%NDEM + 1) % NDEM] + dw / 2;
			e[0] = cube->cell[icell].loc[(iface%NDEM + NDEM-1) % NDEM] - dw / 2;
			e[1] = cube->cell[icell].loc[(iface%NDEM + NDEM - 1) % NDEM] + dw / 2;
			glBegin(GL_POLYGON);
			glColor4f(cube->facecolor[cube->cell[icell].nface][0], cube->facecolor[cube->cell[icell].nface][1],
				cube->facecolor[cube->cell[icell].nface][2], cube->facecolor[cube->cell[icell].nface][3]);
			for (int ii = 0; ii < 4; ii++)
			{
				v3f[iface%NDEM] = d;
				v3f[(iface%NDEM + 1) % NDEM] = s[idx[ii][0]];
				v3f[(iface%NDEM + NDEM - 1) % NDEM] = e[idx[ii][1]];
				glVertex3fv(v3f);
			}
			glEnd();
			glEnable(GL_BLEND);
			glEnable(GL_LINE_SMOOTH);
			glEnable(GL_LINE_STIPPLE);
			glLineWidth(8.f);
			glBegin(GL_LINE_LOOP);
			glColor3f(0.5f, 0.5f, 0.5f);
			for (int ii = 0; ii < 4; ii++)
			{
				v3f[iface%NDEM] = d;
				v3f[(iface%NDEM + 1) % NDEM] = s[idx[ii][0]];
				v3f[(iface%NDEM + NDEM - 1) % NDEM] = e[idx[ii][1]];
				glVertex3fv(v3f);
			}
			glEnd();
		}
	}
}
void Operating(SuperCubes* cube, int idem, int ilayer, int times)
{
	int rot[NDEM][NDEM] = { 0 };
	float phi = atan2(1.,0.)*times,cp=cos(phi),sp=sin(phi);
	Cells cell;
	for (int iidem = 0; iidem < NDEM; iidem++)
		rot[iidem][iidem] = 1.;
	rot[(idem + 1) % NDEM][(idem + 1) % NDEM] = cp;
	rot[(idem + 1) % NDEM][(idem - 1 + NDEM) % NDEM] = sp;
	rot[(idem - 1 + NDEM) % NDEM][(idem + 1) % NDEM] = -sp;
	rot[(idem - 1 + NDEM) % NDEM][(idem - 1 + NDEM) % NDEM] = cp;
	for (int icell = 0; icell < NSEG*NSEG*NDEM * 2;icell++)
	if (IsCellOnLayer(cube, cube->cell + icell, idem, ilayer))
	{
		for (int iidem = 0; iidem < NDEM; iidem++)
			cell.loc[iidem] = cube->cell[icell].loc[iidem] - cube->center[iidem];
		for (int iidem = 0; iidem < NDEM; iidem++)
		{
			cube->cell[icell].loc[iidem] = cube->center[iidem];
			for (int jjdem = 0; jjdem < NDEM; jjdem++)
				cube->cell[icell].loc[iidem] += rot[iidem][jjdem]*cell.loc[jjdem];
		}
		UpdateCellFace(cube, cube->cell + icell);
	}
}
void Roting(SuperCubes* cube, int idem, float phi)
{
	int rot[NDEM][NDEM] = { 0 };
	float cp = cos(phi), sp = sin(phi);
	Cells cell;
	for (int iidem = 0; iidem < NDEM; iidem++)
		rot[iidem][iidem] = 1.;
	rot[(idem + 1) % NDEM][(idem + 1) % NDEM] = cp;
	rot[(idem + 1) % NDEM][(idem - 1 + NDEM) % NDEM] = sp;
	rot[(idem - 1 + NDEM) % NDEM][(idem + 1) % NDEM] = -sp;
	rot[(idem - 1 + NDEM) % NDEM][(idem - 1 + NDEM) % NDEM] = cp;
	for (int icell = 0; icell < NSEG*NSEG*NDEM * 2; icell++)
	{
		for (int iidem = 0; iidem < NDEM; iidem++)
			cell.loc[iidem] = cube->cell[icell].loc[iidem] - cube->center[iidem];
		for (int iidem = 0; iidem < NDEM; iidem++)
		{
			cube->cell[icell].loc[iidem] = cube->center[iidem];
			for (int jjdem = 0; jjdem < NDEM; jjdem++)
				cube->cell[icell].loc[iidem] += rot[iidem][jjdem] * cell.loc[jjdem];
		}
	}
}
SuperCubes* CreateSuperCubes()
{
	SuperCubes*cube = (SuperCubes*)malloc(sizeof(SuperCubes));
	SPos center = { 0, 0, 0 };
	InitSuperCube(cube, center, 1); 
	return cube;
}
void BrokenSuperCube(SuperCubes* cube){
	FILE*fid = fopen("Opera.txt", "w");
	INT op = 0;
	INT nlayer = 0;
	int nrand = 5;
	char ops[] = { 'L', 'R', 'U', 'D', 'F', 'B' };
	srand(time(NULL));
	for (int i = 0; i < nrand; i++){
		op = rand() % 6;
		nlayer = rand() % NSEG+1;
		if(op%2==0)
			Operating(cube, op / 2, nlayer, 1);
		else
			Operating(cube, op / 2, -nlayer, 1);
		fprintf(fid, "%6d %c %3d \n", i, ops[op], nlayer);
	}
	fclose(fid);
}