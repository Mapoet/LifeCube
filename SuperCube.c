
#include "SuperCube.h"
#include "base.h"
#include "math.h"
#include <time.h>
SuperUniverse* CreateSuperUniverse(SuperInit init, SuperControl*pCtr){
	SuperUniverse*pUni = (SuperUniverse*)malloc(sizeof(SuperUniverse));
	memset(pUni, 0, sizeof(SuperUniverse));
	int inum = 0;
	srand(time(NULL));
	pUni->width = init.width;
	pUni->center[0] = 0;
	pUni->center[1] = 0;
	pUni->center[2] = 0;
	pUni->nseg = init.nseg;
	pUni->prelife = (int*)malloc(sizeof(int)*pUni->nseg*pUni->nseg*pUni->nseg);
	memset(pUni->prelife, 0, sizeof(int)*pUni->nseg*pUni->nseg*pUni->nseg);
	pUni->life = (int*)malloc(sizeof(int)*pUni->nseg*pUni->nseg*pUni->nseg);
	pUni->prepow = (float*)malloc(sizeof(float)*pUni->nseg*pUni->nseg*pUni->nseg);
	memset(pUni->prepow, 0, sizeof(float)*pUni->nseg*pUni->nseg*pUni->nseg);
	pUni->power = (float*)malloc(sizeof(float)*pUni->nseg*pUni->nseg*pUni->nseg);
	for (int i = 0; i < pUni->nseg*pUni->nseg*pUni->nseg; i++){
		pUni->power[i] = gasdev(&inum)*init.sigma + init.nu ;
		pUni->life[i] = ((pUni->power[i] >= pCtr->pow_start&&pUni->power[i] <= pCtr->pow_end) ? 0 : -2);
		pUni->prelife[i] = -2;
	}
	for (int i = 0; i < 15;i++)
	for (int j = 0; j < 4; j++)
		pUni->cls[i][j] = ran1(&inum);
	return pUni;
}
void NextSuperUniverse(SuperUniverse*pUni, SuperControl*pCtr){
	float*power = (float*)malloc(sizeof(float)*pUni->nseg*pUni->nseg*pUni->nseg);
	int*life = (int*)malloc(sizeof(int)*pUni->nseg*pUni->nseg*pUni->nseg);
	double dw = pUni->width / pUni->nseg;
	double hw = dw / 2;
	double  effect = 0;
	float enu[3][3] = { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
	float phi[3] = { 0, 0, 0 };
	SPos move;
	SPos  c;
	float m;
	int inum;
	srand(time(NULL));
	for (int i = 0; i < pUni->nseg*pUni->nseg*pUni->nseg; i++)
	{
		c[0] = i / (pUni->nseg*pUni->nseg);
		c[1] = (i / (pUni->nseg)) % pUni->nseg;
		c[2] = i%pUni->nseg;
		for (int ii = 0; ii < 3; ii++)c[ii] = (c[ii] - pUni->nseg / 2.)*dw + pUni->center[ii];
		pCtr->Move(pUni, i, pCtr, move);
		power[i] = pCtr->Effect(pUni, i, pCtr);
		m = sqrt(sqrt(move[0] * move[0] + move[1] * move[1] + move[2] * move[2])*hw);
		if (m>pCtr->move_limit&&pUni->life[i] >= 0&&power[i] >= pCtr->pow_start&&power[i] <= pCtr->pow_end)
			life[i] = -4;
		if (pUni->life[i] <= -4)
			power[i] = (pCtr->pow_start + pCtr->pow_end) / 2;
		if (power[i] >= pCtr->pow_start&&power[i] <= pCtr->pow_end){
			if (pUni->life[i] <= -4)
				life[i] = pUni->life[i] - 1;
			else if (pUni->life[i] == -2 || pUni->life[i] == -1 || pUni->life[i] == -3)
				life[i] = 0;
			else
				life[i] = pUni->life[i] + 1;
			if (pUni->life[i]<-4&&abs(pUni->life[i])>10*pCtr->life_period + gasdev(&inum)*pCtr->life_sigma)
				life[i] = -3;
			if (pUni->life[i]>=4&&abs(pUni->life[i])>pCtr->life_period + gasdev(&inum)*pCtr->life_sigma)
				life[i] = -3;
		}
		else{
			if (pUni->power[i] >= pCtr->pow_start&&pUni->power[i] <= pCtr->pow_end)
				life[i] = -1;
			else
				life[i] = -2;
		}
		if (pUni->life[i] == -2)
		{
			if (power[i] <= pCtr->pow_start)
				power[i] = pCtr->pow_start - 5 * (pCtr->pow_end - pCtr->pow_start);
			else
				power[i] = pCtr->pow_end + 5 * (pCtr->pow_end - pCtr->pow_start);
		}
	}
	memcpy(pUni->prepow, pUni->power, sizeof(float)*pUni->nseg*pUni->nseg*pUni->nseg);
	memcpy(pUni->power, power, sizeof(float)*pUni->nseg*pUni->nseg*pUni->nseg);
	memcpy(pUni->prelife, pUni->life, sizeof(int)*pUni->nseg*pUni->nseg*pUni->nseg);
	memcpy(pUni->life, life, sizeof(int)*pUni->nseg*pUni->nseg*pUni->nseg);
	free(power);
}
void DrawLine(Colorsf cls, SPos move){
	glBegin(GL_LINE_LOOP);
	glColor4f(cls[0], cls[1], cls[2], cls[3]);
	glVertex3f(0.f, 0.f, 0.f);
	glVertex3f(move[0], move[1], move[2]);
	glEnd();
}
void DrawCube(Colorsf cls, float hw){
	hw *= 0.6;
	////glutSolidSphere(hw, 20, 5);
	//glColor4f(cls[0], cls[1], cls[2], cls[3]);
	//glutWireSphere(hw, 20, 5);
	int pos[8][3] = { { 1, 1, 1 }, { 1, 1, -1 }, { 1, -1, -1 }, { 1, -1, 1 }, { -1, 1, 1 }, { -1, 1, -1 }, { -1, -1, -1 }, { -1, -1, 1 } };
	for (int i = 0; i < 3; i++){
		glBegin(GL_POLYGON);
		glColor4f(cls[0], cls[1], cls[2], cls[3]);
		for (int j = 0; j < 4; j++)
			glVertex3f(pos[j][(i + 0) % NDEM] * hw, pos[j][(i + 1) % NDEM] * hw, pos[j][(i + 2) % NDEM] * hw);
		glEnd();
		glBegin(GL_POLYGON);
		glColor4f(cls[0], cls[1], cls[2], cls[3]);
		for (int j = 0; j < 4; j++)
			glVertex3f(pos[j + 4][(i + 0) % NDEM] * hw, pos[j + 4][(i + 1) % NDEM] * hw, pos[j + 4][(i + 2) % NDEM] * hw);
		glEnd();
		glBegin(GL_LINE_LOOP);
		glColor4f(0, 0, 0, 0);
		for (int j = 0; j < 4; j++)
			glVertex3f(pos[j][(i + 0) % NDEM] * hw, pos[j][(i + 1) % NDEM] * hw, pos[j][(i + 2) % NDEM] * hw);
		glEnd();
		glBegin(GL_LINE_LOOP);
		glColor4f(0, 0, 0, 0);
		for (int j = 0; j < 4; j++)
			glVertex3f(pos[j + 4][(i + 0) % NDEM] * hw, pos[j + 4][(i + 1) % NDEM] * hw, pos[j + 4][(i + 2) % NDEM] * hw);
		glEnd();
	}
}
void ShowSuperUniverse(SuperUniverse*pUni, SuperControl*pCtr){
	double dw = pUni->width / pUni->nseg;
	double hw = dw / 2;
	float enu[3][3] = { {1,0,0}, {0,1,0}, {0,0,1} };
	float phi[3] = { 0, 0, 0 };
	int idxs[NDEM], idxsp[NDEM], idxp = 0;
	SPos  c,cp;
	for (int ik = 0; ik < pUni->nseg*pUni->nseg*pUni->nseg; ik++){
		if (pUni->life[ik] == -2)
			continue;
		idxs[0] = ik / (pUni->nseg*pUni->nseg);
		idxs[1] = (ik / (pUni->nseg)) % pUni->nseg;
		idxs[2] = ik%pUni->nseg;
		for (int ii = 0; ii < 3; ii++)c[ii] = (idxs[ii] - pUni->nseg / 2.)*dw + pUni->center[ii];
		glPushMatrix();
		glTranslatef(c[0], c[1], c[2]);
		for (int i = 0; i < 3; i++)glRotatef(phi[i], enu[i][0], enu[i][1], enu[i][2]);
		if (pUni->life[ik] <= -4)DrawCube(pUni->cls[0], hw);
		if (pUni->life[ik] == -3)DrawCube(pUni->cls[1], hw);
		if (pUni->life[ik] == -1)DrawCube(pUni->cls[2], hw);
		if (pUni->life[ik] == 0)DrawCube(pUni->cls[3], hw);
		if (pUni->life[ik] == 1)DrawCube(pUni->cls[4], hw);
		if (pUni->life[ik] >= 2)DrawCube(pUni->cls[5], hw);
		for (int i = -pCtr->effect_end; i < pCtr->effect_end; i++)if (i<-pCtr->effect_start || i>pCtr->effect_start)
		for (int j = -pCtr->effect_end; j < pCtr->effect_end; j++)if (j<-pCtr->effect_start || j>pCtr->effect_start)
		for (int k = -pCtr->effect_end; k < pCtr->effect_end; k++)if (k<-pCtr->effect_start || k>pCtr->effect_start){
			if (idxs[0] + i < 0 || idxs[0] + i >= pUni->nseg)continue;
			if (idxs[1] + j < 0 || idxs[1] + j >= pUni->nseg)continue;
			if (idxs[2] + k < 0 || idxs[2] + k >= pUni->nseg)continue;
			idxsp[0] = (idxs[0] + pUni->nseg + i) % pUni->nseg;
			idxsp[1] = (idxs[1] + pUni->nseg + j) % pUni->nseg;
			idxsp[2] = (idxs[2] + pUni->nseg + k) % pUni->nseg;
			idxp = (idxsp[0] * pUni->nseg + idxsp[1]) * pUni->nseg + idxsp[2];
			cp[0] = i*dw;
			cp[1] = j*dw;
			cp[2] = k*dw;
			if (pUni->prelife[idxp] <= -4)
			{
				DrawLine(pUni->cls[6], cp);
				if (pUni->life[idxp] == -2)
				{
					glPushMatrix();
					glTranslatef(cp[0], cp[1], cp[2]);
					DrawCube(pUni->cls[9], hw);
					glPopMatrix();
				}
			}
			if (pUni->prelife[idxp] >= 0)
			{
				DrawLine(pUni->cls[7], cp);
				if (pUni->life[idxp] == -2)
				{
					glPushMatrix();
					glTranslatef(cp[0], cp[1], cp[2]);
					DrawCube(pUni->cls[10], hw);
					glPopMatrix();
				}
			}
			if (pUni->prelife[idxp] == -1)
			{
				DrawLine(pUni->cls[8], cp);
				if (pUni->life[idxp] == -2)
				{
					glPushMatrix();
					glTranslatef(cp[0], cp[1], cp[2]);
					DrawCube(pUni->cls[11], hw);
					glPopMatrix();
				}
			}
		}
		glPopMatrix();
	}
}
void FreeSuperUniverse(SuperUniverse*pUni)
{
	if (pUni){
		if (pUni->life)free(pUni->life);
		if (pUni->power)free(pUni->power);
		if (pUni->prepow)free(pUni->prepow);
		free(pUni);
	}
}
float Effect(SuperUniverse*pUni, int idx,SuperControl*pCtr){
	int idxs[NDEM], idxsp[NDEM], idxp = 0;
	float effect = 0;// pUni->power[idx];
	float avepow = (pCtr->effect_start + pCtr->effect_end) / 2., dist;
	float rate = 0;
	float dw = pUni->width / pUni->nseg;
	idxs[0] = idx / pUni->nseg / pUni->nseg ;
	idxs[1] = idx / pUni->nseg % pUni->nseg;
	idxs[2] = idx % pUni->nseg;
	for (int i = -pCtr->effect_end; i < pCtr->effect_end; i++)if (i<-pCtr->effect_start || i>pCtr->effect_start)
	for (int j = -pCtr->effect_end; j < pCtr->effect_end; j++)if (j<-pCtr->effect_start || j>pCtr->effect_start)
	for (int k = -pCtr->effect_end; k < pCtr->effect_end; k++)if (k<-pCtr->effect_start || k>pCtr->effect_start){
		idxsp[0] = (idxs[0] + pUni->nseg + i) % pUni->nseg;
		idxsp[1] = (idxs[1] + pUni->nseg + j) % pUni->nseg;
		idxsp[2] = (idxs[2] + pUni->nseg + k) % pUni->nseg;
		idxp = (idxsp[0] * pUni->nseg + idxsp[1]) * pUni->nseg + idxsp[2];
		dist = pow(sqrt(i*i + j*j + k*k), 1);
		if (pUni->life[idxp] == -3)effect += (pUni->power[idxp] ) * 1 / dist;
		if (pUni->life[idxp] == -1)effect += (pUni->power[idxp] ) * 0.01 / dist;
		if (pUni->life[idxp] >= 0){ 
			rate = (pCtr->life_period - abs(pUni->life[idxp] - pCtr->life_period / 2)) / pCtr->life_period;
			effect += (pUni->power[idxp] ) * rate / dist;
		}
	}
	if (effect > 1e-6){
		int  fff = effect;
	}
	return effect;
}
void Move(SuperUniverse*pUni, int idx, SuperControl*pCtr, SPos move){
	int idxs[NDEM], idxsp[NDEM], idxp = 0;
	idxs[0] = idx / pUni->nseg / pUni->nseg;
	idxs[1] = idx / pUni->nseg % pUni->nseg;
	idxs[2] = idx % pUni->nseg;
	move[0] = move[1] = move[2] = 0;
	float rate = 0;
	float dw = pUni->width / pUni->nseg;
	float dist;
	for (int i = -pCtr->effect_end; i < pCtr->effect_end; i++)if (i<-pCtr->effect_start || i>pCtr->effect_start)
	for (int j = -pCtr->effect_end; j < pCtr->effect_end; j++)if (j<-pCtr->effect_start || j>pCtr->effect_start)
	for (int k = -pCtr->effect_end; k < pCtr->effect_end; k++)if (k<-pCtr->effect_start || k>pCtr->effect_start){
		idxsp[0] = (idxs[0] + pUni->nseg + i) % pUni->nseg;
		idxsp[1] = (idxs[1] + pUni->nseg + j) % pUni->nseg;
		idxsp[2] = (idxs[2] + pUni->nseg + k) % pUni->nseg;
		idxp = (idxsp[0] * pUni->nseg + idxsp[1]) * pUni->nseg + idxsp[2];
		dist = pow(sqrt(i*i + j*j + k*k)*dw, 3);
		if (pUni->life[idxp] <= -4){ move[0] += i / dist*dw; move[1] += j / dist*dw; move[2] += k / dist*dw; }
		//if (pUni->life[idxp] == -3){ move[i] += i; move[1] += j; move[2] += k; }
		if (pUni->life[idxp] == -1){ move[0] += i*0.01 / dist*dw; move[1] += j*0.01 / dist*dw; move[2] += k*0.01 / dist*dw; }
		if (pUni->life[idxp] >= 0){ rate = (pCtr->life_period - abs(pUni->life[idxp] - pCtr->life_period / 2)) / pCtr->life_period;
		move[0] += i*rate / dist*dw; move[1] += j*rate / dist*dw; move[2] += k*rate / dist*dw;
		}
	}
}