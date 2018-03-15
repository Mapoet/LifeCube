#pragma once

#ifndef __BASE__
#define __BASE__
//common includes
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<string.h>
#include<Windows.h>
#ifndef __GL_H__
#include "GL.h"
#endif
#ifndef __GLU_H__
#include "GLU.h"
#endif
#ifndef __glut_h__
#include "GLUT.h"
#endif
#ifndef __GLAUX_H__
#include "GLAUX.h"
#endif
#pragma comment(lib,"OpenGL32.Lib")
#pragma comment(lib, "glu32.lib")
#pragma comment(lib, "glut32.lib")
#pragma comment(lib, "GLAUX.lib")
//#define __CPlusPlus__

typedef unsigned int	uint;
typedef unsigned char	uchar;
typedef long double		real;

#define NDEM	3
#define NROW	19
#define NCOL	19

extern uint		nDem;
extern uint		pidx[2];
extern uint		gidx[3];
extern uint		nRow;
extern uint		nCol;
extern POINT	nPos;
extern uint		nWidth;
extern uint		nHight;

typedef struct{  
	real pos[NDEM];
}Pos;
#ifdef __CPlusPlus__
real abs(Pos p);
Pos operator +(Pos p1, Pos p2);
Pos operator -(Pos p1, Pos p2);
Pos operator *(Pos p1, Pos p2);
Pos operator /(Pos p1, Pos p2);
Pos operator +(real p1, Pos p2);
Pos operator -(real p1, Pos p2);
Pos operator *(real p1, Pos p2);
Pos operator /(real p1, Pos p2);
Pos operator +(Pos p1, real p2);
Pos operator -(Pos p1, real p2);
Pos operator *(Pos p1, real p2);
Pos operator /(Pos p1, real p2);
#endif
real sprod(Pos p1, Pos p2);
Pos  vprod(Pos p1, Pos p2);
Pos  norms(Pos p);
real pabs(Pos p);
Pos ppp(Pos p1, Pos p2);
Pos psp(Pos p1, Pos p2);
Pos pmp(Pos p1, Pos p2);
Pos pdp(Pos p1, Pos p2);
Pos vpp(real p1, Pos p2);
Pos vsp(real p1, Pos p2);
Pos vmp(real p1, Pos p2);
Pos vdp(real p1, Pos p2);
Pos ppv(Pos p1, real p2);
Pos psv(Pos p1, real p2);
Pos pmv(Pos p1, real p2);
Pos pdv(Pos p1, real p2);
POINT np2p2(Pos p);

typedef struct _Wnd{
	HDC			hDC;		// Private GDI Device Context
	HGLRC		hRC;		// Permanent Rendering Context
	HWND		hWnd;		// Holds Our Window Handle
	HINSTANCE	hInstance;
	HINSTANCE	hPrevInstance;
	RECT		hRect;
	uint		nRow;
	uint		nCol;
	uint		nBits;
	BOOL		isFullScreen;
	BOOL        isGL;
	char		nTitle[30];
	uint        nPixelFormat;
	WNDPROC		nProc;
	INT(*Draw)(struct _Wnd*);
	BOOL		isDone;
}*pWnd,Wnd;
typedef struct{
	BOOL		Keys[256];
	uint		mKeys[3];
	uint        Event;
	POINT		M, Mp,Mv,nSize;
	BOOL		isActive;
	BOOL		isResize;
}*pCtrl,Ctrl;
typedef struct{
	Pos Center;
	Pos Eye;
	Pos Up;
	Pos Rot;
}*pView,View;
extern  pWnd pWindow;
extern  pCtrl pControl;
extern  pView pViewer;
INT initGlobelView(HINSTANCE	hInstance,			// Instance
	HINSTANCE	hPrevInstance,		// Previous Instance
	LPSTR		lpCmdLine,			// Command Line Parameters
	INT			nCmdShow);
VOID destroyGlobelView();
INT createWindow(pWnd pWin);
VOID destroyWindow(pWnd pWin);
LRESULT CALLBACK WndProc(HWND	hWnd,			// Handle For This Window
	UINT	uMsg,			// Message For This Window
	WPARAM	wParam,			// Additional Message Information
	LPARAM	lParam);			// Additional Message Information
VOID KillWindow(pWnd pWin);
BOOL InitGL();
VOID ReSizeWnd(pWnd pWin);
INT UpdateKeys(pCtrl pCon, pWnd pWin);
void drawString(const char* str);
void GetString(char*str);
// func
//BOOL isInRect(POINT p, RECT rect, BOOL isWithEdge);
//BOOL isInRange(POINT p, POINT sp, POINT sz, BOOL isWithEdge);
//VOID drawShapeN(HDC hDC, uint x, uint y, uint N, uint Edge, real s);
//VOID drawSharpShapeN(HDC hDC, uint x, uint y, uint N, uint Edge, real rate, real s);
//
INT WinDraw(pWnd pWin);
VOID PrepareGL(pWnd pWin);
INT GLDraw(pWnd pWin);


#endif