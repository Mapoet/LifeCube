

#include "base.h"
#include "SuperCube.h"
uint	nDem = 3;
uint	pidx[2] = { 0, 1 };
uint	gidx[3] = { 0, 1, 2 };
uint	nRow = 19;
uint	nCol = 19;
POINT	nPos = {320,240};
uint	nWidth = 640;
uint	nHight = 480;
pWnd	pWindow = NULL;
pCtrl	pControl = NULL;
pView	pViewer = NULL;
SuperUniverse*pUni = NULL;
/*0.45~0.55*/
SuperControl pCtr = {2.0,2.9,2,4,1.0,50,10,Effect,Move};
SuperInit init = {1,40,2.0,1.5};
#ifdef __CPlusPlus__
real abs(Pos p){ return pabs(p); }
Pos operator +(Pos p1, Pos p2){ return ppp(p1, p2); }
Pos operator -(Pos p1, Pos p2){ return psp(p1, p2); }
Pos operator *(Pos p1, Pos p2){ return psp(p1, p2); }
Pos operator /(Pos p1, Pos p2){ return psp(p1, p2); }
Pos operator +(real p1, Pos p2){ return vsp(p1, p2); }
Pos operator -(real p1, Pos p2){ return vsp(p1, p2); }
Pos operator *(real p1, Pos p2){ return vsp(p1, p2); }
Pos operator /(real p1, Pos p2){ return vsp(p1, p2); }
Pos operator +(Pos p1, real p2){ return psv(p1, p2); }
Pos operator -(Pos p1, real p2){ return psv(p1, p2); }
Pos operator *(Pos p1, real p2){ return psv(p1, p2); }
Pos operator /(Pos p1, real p2){ return psv(p1, p2); }
#endif
real sprod(Pos p1, Pos p2)
{
	real result = 0.;
	for (uint i = 0; i < nDem; i++)
		result += p1.pos[i] * p2.pos[i];
	return result;
}
Pos vprod(Pos p1, Pos p2)
{
	Pos result;
	for (uint i = 0; i < nDem; i++)
		result.pos[i] = p1.pos[(i + 1) % nDem] * p2.pos[(i - 1 + nDem) % nDem] - p2.pos[(i + 1) % nDem] * p1.pos[(i - 1 + nDem) % nDem];
	return result;
}
Pos norms(Pos p)
{
	Pos result;
	real aa = sqrt(sprod(p, p));
	if (fabs(aa) < 1e-10)
		return p;
	for (uint i = 0; i < nDem; i++)result.pos[i] = p.pos[i] / aa;
	return result;
}
real pabs(Pos p)
{
	return sqrt(sprod(p, p));
}
Pos ppp(Pos p1, Pos p2)
{
	Pos result;
	for (uint i = 0; i < nDem; i++)result.pos[i] = p1.pos[i] + p2.pos[i];
	return result;
}
Pos psp(Pos p1, Pos p2)
{
	Pos result;
	for (uint i = 0; i < nDem; i++)result.pos[i] = p1.pos[i] - p2.pos[i];
	return result;
}
Pos pmp(Pos p1, Pos p2)
{
	Pos result;
	for (uint i = 0; i < nDem; i++)result.pos[i] = p1.pos[i] * p2.pos[i];
	return result;
}
Pos pdp(Pos p1, Pos p2)
{
	Pos result;
	for (uint i = 0; i < nDem; i++)result.pos[i] = p1.pos[i] / p2.pos[i];
	return result;
}
Pos ppv(Pos p1, real p2)
{
	Pos result;
	for (uint i = 0; i < nDem; i++)result.pos[i] = p1.pos[i] + p2;
	return result;
}
Pos psv(Pos p1, real p2)
{
	Pos result;
	for (uint i = 0; i < nDem; i++)result.pos[i] = p1.pos[i] - p2;
	return result;
}
Pos pmv(Pos p1, real p2)
{
	Pos result;
	for (uint i = 0; i < nDem; i++)result.pos[i] = p1.pos[i] * p2;
	return result;
}
Pos pdv(Pos p1, real p2)
{
	Pos result;
	for (uint i = 0; i < nDem; i++)result.pos[i] = p1.pos[i] / p2;
	return result;
}
Pos vpp(real p1, Pos p2)
{
	Pos result;
	for (uint i = 0; i < nDem; i++)result.pos[i] = p1 + p2.pos[i];
	return result;
}
Pos vsp(real p1, Pos p2)
{
	Pos result;
	for (uint i = 0; i < nDem; i++)result.pos[i] = p1 - p2.pos[i];
	return result;
}
Pos vmp(real p1, Pos p2)
{
	Pos result;
	for (uint i = 0; i < nDem; i++)result.pos[i] = p1 * p2.pos[i];
	return result;
}
Pos vdp(real p1, Pos p2)
{
	Pos result;
	for (uint i = 0; i < nDem; i++)result.pos[i] = p1 / p2.pos[i];
	return result;
}
POINT np2p2(Pos p){
	POINT P = { (LONG)p.pos[pidx[0]], (LONG)p.pos[pidx[1]] };
	return P;
}
INT initGlobelView(HINSTANCE	hInstance, HINSTANCE	hPrevInstance, LPSTR lpCmdLine, INT	nCmdShow){
	pWindow = (pWnd)malloc(sizeof(Wnd));
	pUni = CreateSuperUniverse(init,&pCtr);
	memset(pWindow, 0, sizeof(Wnd));
	pViewer = (pView)malloc(sizeof(View));
	memset(pViewer, 0, sizeof(View));
	pControl = (pCtrl)malloc(sizeof(Ctrl));
	memset(pControl, 0, sizeof(Ctrl));
	pWindow->hInstance = hInstance;
	pWindow->hPrevInstance = hPrevInstance;
	pWindow->hRect.left		= nPos.x - nWidth / 2;
	pWindow->hRect.right	= nPos.x + nWidth / 2;
	pWindow->hRect.top		= nPos.y - nHight / 2;
	pWindow->hRect.bottom	= nPos.y + nHight / 2;
	strcpy(pWindow->nTitle,"Magic Cube");
	pWindow->isGL =  TRUE;// FALSE;
	if (pWindow->isGL)
		pWindow->Draw = GLDraw;
	else
		pWindow->Draw = WinDraw;
	pWindow->nProc = WndProc;
	pWindow->nBits = 24;
	pWindow->nRow = NROW;
	pWindow->nCol = NCOL;
	pWindow->isFullScreen = FALSE;
	pWindow->isDone = FALSE;
	pViewer->Eye.pos[0] = 1;
	pViewer->Up.pos[2] = 1;
	for (uint i = 0; i < 256; i++)pControl->Keys[i] = 0;
	pControl->mKeys[0] = 0;
	pControl->mKeys[1] = 0;
	pControl->mKeys[2] = 0;
	return 0;
}
VOID destroyGlobelView(){
	free(pWindow);
	free(pViewer);
	free(pControl);
	FreeSuperUniverse(pUni);
}
INT createWindow(pWnd pWin){
	WNDCLASS	wc;						// Windows Class Structure
	DWORD		dwExStyle;				// Window Extended Style
	DWORD		dwStyle;				// Window Style

	pWin->hInstance = GetModuleHandle(NULL);				// Grab An Instance For Our Window
	wc.style = CS_HREDRAW | CS_VREDRAW | CS_OWNDC;	// Redraw On Size, And Own DC For Window.
	wc.lpfnWndProc =pWin->nProc;					// WndProc Handles Messages
	wc.cbClsExtra = 0;									// No Extra Window Data
	wc.cbWndExtra = 0;									// No Extra Window Data
	wc.hInstance = pWin->hInstance;							// Set The Instance
	wc.hIcon = LoadIcon(NULL, IDI_WINLOGO);			// Load The Default Icon
	wc.hCursor = LoadCursor(NULL, IDC_ARROW);			// Load The Arrow Pointer
	wc.hbrBackground = NULL;									// No Background Required For GL
	wc.lpszMenuName = NULL;									// We Don't Want A Menu
	wc.lpszClassName = "OpenGL";								// Set The Class Name

	if (!RegisterClass(&wc))									// Attempt To Register The Window Class
	{
		MessageBox(NULL, "Failed To Register The Window Class.", "ERROR", MB_OK | MB_ICONEXCLAMATION);
		return FALSE;											// Return FALSE
	}

	if (pWin->isFullScreen)												// Attempt Fullscreen Mode?
	{
		DEVMODE dmScreenSettings;								// Device Mode
		memset(&dmScreenSettings, 0, sizeof(dmScreenSettings));	// Makes Sure Memory's Cleared
		dmScreenSettings.dmSize = sizeof(dmScreenSettings);		// Size Of The Devmode Structure
		dmScreenSettings.dmPelsWidth = pWin->hRect.right - pWin->hRect.left;				// Selected Screen Width
		dmScreenSettings.dmPelsHeight = pWin->hRect.bottom - pWin->hRect.top;				// Selected Screen Height
		dmScreenSettings.dmBitsPerPel = pWin->nBits;					// Selected Bits Per Pixel
		dmScreenSettings.dmFields = DM_BITSPERPEL | DM_PELSWIDTH | DM_PELSHEIGHT;

		// Try To Set Selected Mode And Get Results.  NOTE: CDS_FULLSCREEN Gets Rid Of Start Bar.
		if (ChangeDisplaySettings(&dmScreenSettings, CDS_FULLSCREEN) != DISP_CHANGE_SUCCESSFUL)
		{
			// If The Mode Fails, Offer Two Options.  Quit Or Use Windowed Mode.
			if (MessageBox(NULL, "The Requested Fullscreen Mode Is Not Supported By\nYour Video Card. Use Windowed Mode Instead?", "NeHe GL", MB_YESNO | MB_ICONEXCLAMATION) == IDYES)
			{
				pWin->isFullScreen = FALSE;		// Windowed Mode Selected.  Fullscreen = FALSE
			}
			else
			{
				// Pop Up A Message Box Letting User Know The Program Is Closing.
				MessageBox(NULL, "Program Will Now Close.", "ERROR", MB_OK | MB_ICONSTOP);
				return FALSE;									// Return FALSE
			}
		}
	}

	if (pWin->isFullScreen)												// Are We Still In Fullscreen Mode?
	{
		dwExStyle = WS_EX_APPWINDOW;								// Window Extended Style
		dwStyle = WS_POPUP;										// Windows Style
		ShowCursor(FALSE);										// Hide Mouse Pointer
	}
	else
	{
		dwExStyle = WS_EX_APPWINDOW | WS_EX_WINDOWEDGE;			// Window Extended Style
		dwStyle = WS_OVERLAPPEDWINDOW;							// Windows Style
	}

	AdjustWindowRectEx(&pWindow->hRect, dwStyle, FALSE, dwExStyle);		// Adjust Window To True Requested Size

	// Create The Window
	if (!(pWin->hWnd = CreateWindowEx(dwExStyle,							// Extended Style For The Window
		"OpenGL",							// Class Name
		pWin->nTitle,								// Window Title
		dwStyle |							// Defined Window Style
		WS_CLIPSIBLINGS |					// Required Window Style
		WS_CLIPCHILDREN,					// Required Window Style
		0, 0,								// Window Position
		pWin->hRect.right - pWin->hRect.left,	// Calculate Window Width
		pWin->hRect.bottom - pWin->hRect.top,	// Calculate Window Height
		NULL,								// No Parent Window
		NULL,								// No Menu
		pWin->hInstance,							// Instance
		NULL)))								// Dont Pass Anything To WM_CREATE
	{
		KillWindow(pWin);								// Reset The Display
		MessageBox(NULL, "Window Creation Error.", "ERROR", MB_OK | MB_ICONEXCLAMATION);
		return FALSE;								// Return FALSE
	}

	PIXELFORMATDESCRIPTOR pfd =				// pfd Tells Windows How We Want Things To Be
	{
		sizeof(PIXELFORMATDESCRIPTOR),				// Size Of This Pixel Format Descriptor
		1,											// Version Number
		PFD_DRAW_TO_WINDOW |						// Format Must Support Window
		PFD_SUPPORT_OPENGL |						// Format Must Support OpenGL
		PFD_DOUBLEBUFFER,							// Must Support Double Buffering
		PFD_TYPE_RGBA,								// Request An RGBA Format
		pWin->nBits,										// Select Our Color Depth
		0, 0, 0, 0, 0, 0,							// Color Bits Ignored
		0,											// No Alpha Buffer
		0,											// Shift Bit Ignored
		0,											// No Accumulation Buffer
		0, 0, 0, 0,									// Accumulation Bits Ignored
		16,											// 16Bit Z-Buffer (Depth Buffer)  
		0,											// No Stencil Buffer
		0,											// No Auxiliary Buffer
		PFD_MAIN_PLANE,								// Main Drawing Layer
		0,											// Reserved
		0, 0, 0										// Layer Masks Ignored
	};

	if (!(pWin->hDC = GetDC(pWin->hWnd)))							// Did We Get A Device Context?
	{
		KillWindow(pWin);								// Reset The Display
		MessageBox(NULL, "Can't Create A GL Device Context.", "ERROR", MB_OK | MB_ICONEXCLAMATION);
		return FALSE;								// Return FALSE
	}

	if (!(pWin->nPixelFormat = ChoosePixelFormat(pWin->hDC, &pfd)))	// Did Windows Find A Matching Pixel Format?
	{
		KillWindow(pWin);								// Reset The Display
		MessageBox(NULL, "Can't Find A Suitable PixelFormat.", "ERROR", MB_OK | MB_ICONEXCLAMATION);
		return FALSE;								// Return FALSE
	}

	if (!SetPixelFormat(pWin->hDC, pWin->nPixelFormat, &pfd))		// Are We Able To Set The Pixel Format?
	{
		KillWindow(pWin);								// Reset The Display
		MessageBox(NULL, "Can't Set The PixelFormat.", "ERROR", MB_OK | MB_ICONEXCLAMATION);
		return FALSE;								// Return FALSE
	}
	if (pWin->isGL){
		if (!(pWin->hRC = wglCreateContext(pWin->hDC)))				// Are We Able To Get A Rendering Context?
		{
			KillWindow(pWin);								// Reset The Display
			MessageBox(NULL, "Can't Create A GL Rendering Context.", "ERROR", MB_OK | MB_ICONEXCLAMATION);
			return FALSE;								// Return FALSE
		}

		if (!wglMakeCurrent(pWin->hDC, pWin->hRC))					// Try To Activate The Rendering Context
		{
			KillWindow(pWin);								// Reset The Display
			MessageBox(NULL, "Can't Activate The GL Rendering Context.", "ERROR", MB_OK | MB_ICONEXCLAMATION);
			return FALSE;								// Return FALSE
		}
	}
	else{
		pWin->hRC = NULL;
	}
	ShowWindow(pWin->hWnd, SW_SHOW);						// Show The Window
	SetForegroundWindow(pWin->hWnd);						// Slightly Higher Priority
	SetFocus(pWin->hWnd);									// Sets Keyboard Focus To The Window
	// Set Up Our Perspective GL Screen
	if (pWin->isGL){ ReSizeWnd(pWin); pControl->isResize = FALSE; }

	if (pWin->isGL&&!InitGL())									// Initialize Our Newly Created GL Window
	{
		KillWindow(pWin);								// Reset The Display
		MessageBox(NULL, "Initialization Failed.", "ERROR", MB_OK | MB_ICONEXCLAMATION);
		return FALSE;								// Return FALSE
	}

	return TRUE;									// Success

}
BOOL InitGL()										// All Setup For OpenGL Goes Here
{
	glShadeModel(GL_SMOOTH);							// Enable Smooth Shading
	glClearDepth(1.0f);									// Depth Buffer Setup
	glEnable(GL_DEPTH_TEST);							// Enables Depth Testing
	glEnable(GL_TEXTURE_2D);
	glEnable(GL_POINT_SMOOTH);
	glDepthFunc(GL_LEQUAL);								// The Type Of Depth Testing To Do
	glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);	// Really Nice Perspective Calculations
	glMatrixMode(GL_PROJECTION);							// Select The Modelview Matrix
	glLoadIdentity();									// Reset The Current Modelview Matrix
	/***************/
	glMatrixMode(GL_MODELVIEW);							// Select The Modelview Matrix
	glClearColor(0.1f, 0.1f, 0.15f, 0.5f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	// Clear Screen And Depth Buffer
	return TRUE;										// Initialization Went OK
}
VOID ReSizeWnd(pWnd pWin)		// Resize And Initialize The GL Window
{
	if (pWin->isGL){
		if (pWin->hRect.bottom == 0)										// Prevent A Divide By Zero By
		{
			pWin->hRect.bottom = 1;										// Making Height Equal One
		}
		pWin->Draw(pWin);
	}
	else{
		pWin->Draw(pWin);
	}
}
VOID destroyWindow(pWnd pWin){

}
LRESULT CALLBACK WndProc(HWND	hWnd,			// Handle For This Window
	UINT	uMsg,			// Message For This Window
	WPARAM	wParam,			// Additional Message Information
	LPARAM	lParam)			// Additional Message Information
{
	pControl->Event = uMsg;
	switch (uMsg)									// Check For Windows Messages
	{
		case WM_ACTIVATE:							// Watch For Window Activate Message
		{
			if (!HIWORD(wParam))					// Check Minimization State
			{
				pControl->isActive = TRUE;						// Program Is Active
			}
			else
			{
				pControl->isActive = FALSE;						// Program Is No Longer Active
			}
			return 0;								// Return To The Message Loop
		}
		case WM_SYSCOMMAND:
		{
			switch (wParam)
			{
				case SC_SCREENSAVE:
				case SC_MONITORPOWER:
					return 0;
			}
			break;
		}
		case WM_CLOSE:								// Did We Receive A Close Message?
		{
			PostQuitMessage(0);						// Send A Quit Message
			return 0;								// Jump Back
		}
		case WM_KEYDOWN:							// Is A Key Being Held Down?
		{
			pControl->Keys[wParam] = TRUE;					// If So, Mark It As TRUE
			return 0;								// Jump Back
		}
		case WM_KEYUP:								// Has A Key Been Released?
		{
			pControl->Keys[wParam] = FALSE;					// If So, Mark It As FALSE
			return 0;								// Jump Back
		}
		case WM_SIZE:								// Resize The OpenGL Window
		{
			pControl->nSize.x = LOWORD(lParam);
			pControl->nSize.y = HIWORD(lParam);
			pControl->isResize = TRUE;
			return 0;								// Jump Back
		}
		case WM_MOUSEMOVE:{
			pControl->Mv.x = LOWORD(lParam);
			pControl->Mv.y = HIWORD(lParam);
			return 0;
		}
		case WM_LBUTTONDOWN:
		{
			pControl->mKeys[0] = TRUE;
			pControl->Mp = pControl->M;
			pControl->M.x = LOWORD(lParam);
			pControl->M.y = HIWORD(lParam);
			return 0;								// Jump Back
		}
		case WM_LBUTTONDBLCLK:
		{
			pControl->mKeys[0] = 2;
			pControl->Mp = pControl->M;
			pControl->M.x = LOWORD(lParam);
			pControl->M.y = HIWORD(lParam);
			return 0;								// Jump Back
		}
		case WM_LBUTTONUP:
		{
			pControl->mKeys[0] = FALSE;
			pControl->Mp = pControl->M;
			pControl->M.x = LOWORD(lParam);
			pControl->M.y = HIWORD(lParam);
			return 0;								// Jump Back
		}
		case WM_MBUTTONDOWN:
		{
			pControl->mKeys[1] = TRUE;
			pControl->Mp = pControl->M;
			pControl->M.x = LOWORD(lParam);
			pControl->M.y = HIWORD(lParam);
			return 0;								// Jump Back
		}
		case WM_MBUTTONDBLCLK:
		{
			pControl->mKeys[1] = 2;
			pControl->Mp = pControl->M;
			pControl->M.x = LOWORD(lParam);
			pControl->M.y = HIWORD(lParam);
			return 0;								// Jump Back
		}
		case WM_MBUTTONUP:
		{
			pControl->mKeys[1] = FALSE;
			pControl->Mp = pControl->M;
			pControl->M.x = LOWORD(lParam);
			pControl->M.y = HIWORD(lParam);
			return 0;								// Jump Back
		}
		case WM_RBUTTONDOWN:
		{
			pControl->mKeys[2] = FALSE;
			pControl->Mp = pControl->M;
			pControl->M.x = LOWORD(lParam);
			pControl->M.y = HIWORD(lParam);
			return 0;								// Jump Back
		}
		case WM_RBUTTONDBLCLK:
		{
			pControl->mKeys[2] = 2;
			pControl->Mp = pControl->M;
			pControl->M.x = LOWORD(lParam);
			pControl->M.y = HIWORD(lParam);
			return 0;								// Jump Back
		}
		case WM_RBUTTONUP:
		{
			pControl->mKeys[2] = FALSE;
			pControl->Mp = pControl->M;
			pControl->M.x = LOWORD(lParam);
			pControl->M.y = HIWORD(lParam);
			return 0;								// Jump Back
		}
	}
	// Pass All Unhandled Messages To DefWindowProc
	return DefWindowProc(hWnd, uMsg, wParam, lParam);
}
VOID KillWindow(pWnd pWin)								// Properly Kill The Window
	{
	if (pWin->isFullScreen)										// Are We In Fullscreen Mode?
		{
			ChangeDisplaySettings(NULL, 0);					// If So Switch Back To The Desktop
			ShowCursor(TRUE);								// Show Mouse Pointer
		}

		if (pWin->isGL&&pWin->hRC)											// Do We Have A Rendering Context?
		{
			if (!wglMakeCurrent(NULL, NULL))					// Are We Able To Release The DC And RC Contexts?
			{
				MessageBox(NULL, "Release Of DC And RC Failed.", "SHUTDOWN ERROR", MB_OK | MB_ICONINFORMATION);
			}

			if (!wglDeleteContext(pWin->hRC))						// Are We Able To Delete The RC?
			{
				MessageBox(NULL, "Release Rendering Context Failed.", "SHUTDOWN ERROR", MB_OK | MB_ICONINFORMATION);
			}
			pWin->hRC = NULL;										// Set RC To NULL
		}

		if (pWin->hDC && !ReleaseDC(pWin->hWnd, pWin->hDC))					// Are We Able To Release The DC
		{
			MessageBox(NULL, "Release Device Context Failed.", "SHUTDOWN ERROR", MB_OK | MB_ICONINFORMATION);
			pWin->hDC = NULL;										// Set DC To NULL
		}

		if (pWin->hWnd && !DestroyWindow(pWin->hWnd))					// Are We Able To Destroy The Window?
		{
			MessageBox(NULL, "Could Not Release hWnd.", "SHUTDOWN ERROR", MB_OK | MB_ICONINFORMATION);
			pWin->hWnd = NULL;										// Set hWnd To NULL
		}

		if (!UnregisterClass("OpenGL", pWin->hInstance))			// Are We Able To Unregister Class
		{
			MessageBox(NULL, "Could Not Unregister Class.", "SHUTDOWN ERROR", MB_OK | MB_ICONINFORMATION);
			pWin->hInstance = NULL;									// Set hInstance To NULL
		}
	}
INT UpdateKeys(pCtrl pCon, pWnd pWin){
	static BOOL isReading=FALSE;
	static INT  nPos = 0;
	static INT  nType = 0;
	static char inStr[300];
	if (pCon->isResize){
		//GetWindowRect(pWin->hWnd, &pWin->hRect);
		pWin->hRect.right = pWin->hRect.left + pCon->nSize.x+10;
		pWin->hRect.bottom = pWin->hRect.top + pCon->nSize.y+30;
		ReSizeWnd(pWin);
		pCon->isResize = FALSE;
	}
	if (isReading){
		if (pCon->Keys[VK_RETURN] == 1){
			pCon->Keys[VK_RETURN] = 0;
			isReading = 0;
			FILE*fid=NULL;
			switch (nType){
			case 5:
				 fid = fopen(inStr, "w");
				 fprintf(fid, "%6d%15.6f\n", pUni->nseg,pUni->width);
				 for (int i = 0; i < NDEM; i++)
					 fprintf(fid, "%15.6f", pUni->center[i]);
				 fprintf(fid, "\n");
				 for (int i = 0; i < 6; i++)
				 {
					 for (int j = 0; j < NCOLOR; j++)
						 fprintf(fid, "%15.6f", pUni->cls[i][j]);
					fprintf(fid, "\n");
				 }
				 for (int i = 0; i < pUni->nseg*pUni->nseg*pUni->nseg; i++)
				 {
					 fprintf(fid, "%6d%6d%15.6f%15.6f", pUni->prelife[i], pUni->life[i], pUni->prepow[i], pUni->power[i]);
					 if (i % 10 == 9)
						 fprintf(fid, "\n");
				 }
				 fprintf(fid, "%15.6f%15.6f\n", init.nu, init.sigma);
				 fprintf(fid, "%15.6f%15.6f%15.6f%15.6f%15.6f%15.6f%15.6f", pCtr.pow_start, pCtr.pow_end, pCtr.effect_start, pCtr.effect_end,
					 pCtr.move_limit,pCtr.life_period,pCtr.life_sigma);
				 fprintf(fid, "FIGURE:\n");
				 fprintf(fid, "****************************************************************\n");
				 for (int i = 0; i < pUni->nseg*pUni->nseg*pUni->nseg; i++){
					 if (pUni->life[i] >= 0)
						 fprintf(fid, "%2.2d", pUni->life[i]);
					 else if (pUni->life[i] <= -4)
						 fprintf(fid, "%2.2d", 50);
					 else
						 fprintf(fid, "  ");
					 if (i % pUni->nseg == pUni->nseg-1)
						 fprintf(fid, "\n");
					 if (i % (pUni->nseg * pUni->nseg) == pUni->nseg * pUni->nseg - 1)
						 fprintf(fid, "****************************************************************\n");
				 }
				 fclose(fid);
				break;
			case 6:	
				fid = fopen(inStr, "r");
				if (!fid){
					MessageBox(NULL, "Warning", "Open File Err", IDOK);
					nPos = 0;
					nType = 0;
					memset(inStr, 0, sizeof(char)* 300);
					return 0; 
				}
				fscanf(fid, "%d%f", &pUni->nseg, &pUni->width);
				if (pUni->life){
					free(pUni->life);
					pUni->life = (int*)malloc(sizeof(int)*pUni->nseg*pUni->nseg*pUni->nseg);
				}
				if (pUni->power){
					free(pUni->power);
					pUni->power = (float*)malloc(sizeof(float)*pUni->nseg*pUni->nseg*pUni->nseg);
				}
				if (pUni->prepow){
					free(pUni->prepow);
					pUni->prepow = (float*)malloc(sizeof(float)*pUni->nseg*pUni->nseg*pUni->nseg);
				}
				for (int i = 0; i < NDEM; i++)
					fscanf(fid, "%f", &pUni->center[i]);
				for (int i = 0; i < 6; i++)
					for (int j = 0; j < NCOLOR; j++)
						fscanf(fid, "%f", &pUni->cls[i][j]);
				for (int i = 0; i < pUni->nseg*pUni->nseg*pUni->nseg; i++)
					fscanf(fid, "%d%d%f%f", &pUni->prelife[i], &pUni->life[i], &pUni->prepow[i], &pUni->power[i]);
				fscanf(fid, "%f%f", &init.nu, &init.sigma);
				fscanf(fid, "%f%f%f%f%f%f%f", &pCtr.pow_start, &pCtr.pow_end, &pCtr.effect_start,& pCtr.effect_end,
					&pCtr.move_limit,& pCtr.life_period,& pCtr.life_sigma);
				fclose(fid);
				break;
			}
			nPos = 0;
			nType = 0;
			memset(inStr, 0, sizeof(char)* 300);
			return 0;
		}
		if (pCon->Keys[VK_ESCAPE] == 1){
			pCon->Keys[VK_ESCAPE] = 0;
			nPos = 0;
			nType = 0;
			memset(inStr, 0, sizeof(char)* 300);
			return 0;
		}
		for (int i = 0; i < 256;i++)
		if (((i >= '0'&&i <= '9') || (i >= 'A'&&i <= 'Z')||i=='_') && pCon->Keys[i] != 0)
		{
			inStr[nPos++] = i;
			pCon->Keys[i] = 0;
		}
		return 0;
	}
	if (pCon->Keys['R']==1)
	{
		//pCon->Keys['R'] = FALSE;
		if (pCon->Keys['X'] == 1){
			//pCon->Keys['X'] = FALSE;
			if (pCon->Keys[VK_SHIFT])
			{
				pViewer->Rot.pos[0] += 5;
				if (pViewer->Rot.pos[0] > 360.0)
					pViewer->Rot.pos[0] -= 360;
				//pCon->Keys[VK_SHIFT] = FALSE;
			}
			else{
				pViewer->Rot.pos[0] -= 5;
				if (pViewer->Rot.pos[0] < -360.0)
					pViewer->Rot.pos[0] += 360;
			}
		}
		if (pCon->Keys['Y'] == 1){
			//pCon->Keys['Y'] = FALSE;
			if (pCon->Keys[VK_SHIFT])
			{
				pViewer->Rot.pos[1] += 5;
				if (pViewer->Rot.pos[1] > 360.0)
					pViewer->Rot.pos[1] -= 360;
				//pCon->Keys[VK_SHIFT] = FALSE;
			}
			else{
				pViewer->Rot.pos[1] -= 5;
				if (pViewer->Rot.pos[1] < -360.0)
					pViewer->Rot.pos[1] += 360;
			}
		}
		if (pCon->Keys['Z'] == 1){
			//pCon->Keys['Z'] = FALSE;
			if (pCon->Keys[VK_SHIFT])
			{
				pViewer->Rot.pos[2] += 5;
				if (pViewer->Rot.pos[2] > 360.0)
					pViewer->Rot.pos[2] -= 360;
				//pCon->Keys[VK_SHIFT] = FALSE;
			}
			else{
				pViewer->Rot.pos[2] -= 5;
				if (pViewer->Rot.pos[2] < -360.0)
					pViewer->Rot.pos[2] += 360;
			}
		}
	}
	if (pCon->Keys['D'] == 1){
		//pCon->Keys['D'] = FALSE;
		if (pCon->Keys['X'] == 1)
		{
			if (pCon->Keys[VK_SHIFT])
				pViewer->Eye.pos[0] -= 0.1;
			else
				pViewer->Eye.pos[0] += 0.1;
		}
		if (pCon->Keys['Y'] == 1)
		{
			if (pCon->Keys[VK_SHIFT])
				pViewer->Eye.pos[1] -= 0.1;
			else
				pViewer->Eye.pos[1] += 0.1;
		}
		if (pCon->Keys['Z'] == 1)
		{
			if (pCon->Keys[VK_SHIFT])
				pViewer->Eye.pos[2] -= 0.1;
			else
				pViewer->Eye.pos[2] += 0.1;
		}
	}
	if (pCon->Keys['O'] == 1){
		pCon->Keys['O'] = 0;
		isReading = 1;
		nType = 5;
	}
	if (pCon->Keys['I'] == 1){
		pCon->Keys['I'] = 0;
		isReading = 1;
		nType = 6;
	}
	if (pCon->Keys[' '] == 1){
		pCon->Keys[' '] = 0;
		NextSuperUniverse(pUni, &pCtr);
	}
	if (pCon->mKeys[0] == 1&&pCon->Event==WM_MOUSEMOVE){
		POINT dM = { pCon->Mv.x - pCon->M.x, pCon->Mv.y - pCon->M.y };

	}
	return 0;
}
BOOL isInRect(POINT p, RECT rect,BOOL isWithEdge){
	if (isWithEdge){
		if (p.x >= rect.left&&p.x <= rect.right&&p.y >= rect.top&&p.y <= rect.bottom)return TRUE;
		else return FALSE;
	}
	else{
		if (p.x > rect.left&&p.x < rect.right&&p.y > rect.top&&p.y < rect.bottom)return TRUE;
		else return FALSE;
	}
}
BOOL isInRange(POINT p, POINT sp, POINT sz, BOOL isWithEdge){
	RECT rect = {sp.x,sp.y,sp.x+sz.x,sp.y+sz.y};
	return isInRect(p, rect, isWithEdge);
}
VOID drawShapeN(HDC hDC, uint x,uint y, uint N, uint Edge, real s){
	real pa=atan2(0,-1)/N*2;
	MoveToEx(hDC, x + Edge*cos(s), y + Edge*sin(s), NULL);
	for (uint i = 0; i < N; i++)
		LineTo(hDC, x + Edge*cos(s + pa*(i + 1)), y + Edge*sin(s + pa*(i + 1)));
}
VOID drawSharpShapeN(HDC hDC, uint x, uint y, uint N, uint Edge,real rate, real s){
	real pa = atan2(0, -1) / N;
	MoveToEx(hDC, x + Edge*cos(s), y + Edge*sin(s), NULL);
	for (uint i = 0; i < 2*N; i++)
		if(i%2)LineTo(hDC, x + Edge*cos(s + pa*(i + 1)), y + Edge*sin(s + pa*(i + 1)));
		else LineTo(hDC, x + rate*Edge*cos(s + pa*(i + 1)), y + rate*Edge*sin(s + pa*(i + 1)));

}
INT WinDraw(pWnd pWin){
	HBRUSH hBrush = CreateSolidBrush(RGB(255, 255, 255));
	FillRect(pWin->hDC, &pWin->hRect, hBrush);
	const char *str = "How old are you?";
	RECT rect = { 0, 0, pWin->hRect.right - pWin->hRect.left, 50 };
	for (int i = 0; i < 5; i++){
		rect.top = i * 20;
		rect.bottom = 20 + i * 20;
		DrawText(pWin->hDC, str, strlen(str), &rect, 1);
		DrawFocusRect(pWin->hDC, &rect);
	}
	drawShapeN(pWin->hDC, 200, 200, 3, 50, 0);
	drawShapeN(pWin->hDC, 200, 300, 4, 50, 0);
	drawShapeN(pWin->hDC, 300, 300, 5, 50, 0);
	drawShapeN(pWin->hDC, 300, 200, 6, 50, 0);
	drawShapeN(pWin->hDC, 400, 300, 7, 50, 0);
	drawShapeN(pWin->hDC, 400, 200, 8, 50, 0);
	drawShapeN(pWin->hDC, 500, 300, 9, 50, 0);
	drawShapeN(pWin->hDC, 500, 200, 10, 50, 0);
	drawSharpShapeN(pWin->hDC, 600, 200, 3, 50, 0.2, 0);
	drawSharpShapeN(pWin->hDC, 600, 300, 4, 50, 0.2, 0);
	drawSharpShapeN(pWin->hDC, 700, 300, 5, 50, 0.2, 0);
	drawSharpShapeN(pWin->hDC, 700, 200, 6, 50, 0.2, 0);
	drawSharpShapeN(pWin->hDC, 800, 300, 7, 50, 0.2, 0);
	drawSharpShapeN(pWin->hDC, 800, 200, 8, 50, 0.2, 0);
	drawSharpShapeN(pWin->hDC, 900, 300, 9, 50, 0.2, 0);
	drawSharpShapeN(pWin->hDC, 900, 200, 10, 50, 0.2, 0);
	return 0;
}

VOID PrepareGL(pWnd pWin){
	glViewport(0, 0, pWin->hRect.right, pWin->hRect.bottom);						// Reset The Current Viewport
	glMatrixMode(GL_PROJECTION);						// Select The Projection Matrix
	glLoadIdentity();									// Reset The Projection Matrix
	// Calculate The Aspect Ratio Of The Window
	//gluLookAt(pViewer->Eye.pos[0], pViewer->Eye.pos[1], pViewer->Eye.pos[2], pViewer->Center.pos[0], pViewer->Center.pos[1], pViewer->Center.pos[2], pViewer->Up.pos[0], pViewer->Up.pos[1], pViewer->Up.pos[2]);
	glMatrixMode(GL_MODELVIEW);							// Select The Modelview Matrix
	glLoadIdentity();									// Reset The Modelview Matrix
	glClearColor(0.5f, 0.5f, 0.5f, 1.f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	// Clear Screen And Depth Buffer
}
INT GLDraw(pWnd pWin){
	PrepareGL(pWin);
	//glRectf(-0.5f, -0.5f, 0.5f, 0.5f);
	//glFlush();
	glRotatef(pViewer->Rot.pos[0], 1.f, 0.f, 0.f);
	glRotatef(pViewer->Rot.pos[1], 0.f, 1.f, 0.f);
	glRotatef(pViewer->Rot.pos[2], 0.f, 0.f, 1.f);
	ShowSuperUniverse(pUni,&pCtr);
	return 0;
}
#define FONT_SIZE       64
#define TEXTURE_SIZE    FONT_SIZE
void selectFont(int size, int charset, const char* face) {
	HFONT hFont = CreateFontA(size, 0, 0, 0, FW_MEDIUM, 0, 0, 0,
		charset, OUT_DEFAULT_PRECIS, CLIP_DEFAULT_PRECIS,
		DEFAULT_QUALITY, DEFAULT_PITCH | FF_SWISS, face);
	HFONT hOldFont = (HFONT)SelectObject(wglGetCurrentDC(), hFont);
	DeleteObject(hOldFont);
}
GLuint drawChar_To_Texture(const char* s) {
	wchar_t w;
	HDC hDC = wglGetCurrentDC();

	// 选择字体字号、颜色
	// 不指定字体名字，操作系统提供默认字体
	// 设置颜色为白色
	selectFont(FONT_SIZE, DEFAULT_CHARSET, "");
	glColor3f(1.0f, 1.0f, 1.0f);

	// 转化为宽字符
	MultiByteToWideChar(CP_ACP, MB_PRECOMPOSED, s, 2, &w, 1);

	// 计算绘制的位置
	{
		int width, x, y;
		GetCharWidth32W(hDC, w, w, &width);    // 取得字符的宽度
		x = (TEXTURE_SIZE - width) / 2;
		y = FONT_SIZE / 8;
		//glWindowPos2iARB(x, y); // 一个扩展函数
	}

	// 绘制字符
	// 绘制前应该将各种可能影响字符颜色的效果关闭
	// 以保证能够绘制出白色的字符
	{
		GLuint list = glGenLists(1);

		glDisable(GL_DEPTH_TEST);
		glDisable(GL_LIGHTING);
		glDisable(GL_FOG);
		glDisable(GL_TEXTURE_2D);

		wglUseFontBitmaps(hDC, w, 1, list);
		glCallList(list);
		glDeleteLists(list, 1);
	}

	// 复制字符像素到纹理
	// 注意纹理的格式
	// 不使用通常的GL_RGBA，而使用GL_LUMINANCE4
	// 因为字符本来只有一种颜色，使用GL_RGBA浪费了存储空间
	// GL_RGBA可能占16位或者32位，而GL_LUMINANCE4只占4位
	{
		GLuint texID;
		glGenTextures(1, &texID);
		glBindTexture(GL_TEXTURE_2D, texID);
		glCopyTexImage2D(GL_TEXTURE_2D, 0, GL_LUMINANCE4,
			0, 0, TEXTURE_SIZE, TEXTURE_SIZE, 0);
		glTexParameteri(GL_TEXTURE_2D,
			GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D,
			GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		return texID;
	}
}
void drawString(const char* str) {
	static int isFirstCall = 1;
	static GLuint lists;

	if (isFirstCall) { // 如果是第一次调用，执行初始化
		// 为每一个ASCII字符产生一个显示列表
		isFirstCall = 0;

		// 申请MAX_CHAR个连续的显示列表编号
		lists = glGenLists(128);

		// 把每个字符的绘制命令都装到对应的显示列表中
		wglUseFontBitmaps(wglGetCurrentDC(), 0, 128, lists);
	}
	// 调用每个字符对应的显示列表，绘制每个字符
	for (; *str != '\0'; ++str)
		glCallList(lists + *str);
}
HWND CreateEdit(
	const HWND hParent,      //父窗口类
	const HINSTANCE hInst,   //应用程序实例
	DWORD dwStyle,                //窗口样式
	const RECT rc,                  //相对于父窗口的位置矩形
	const int id,                  //文本框ID号
	const char* caption)  //文本框标题（初始文字）
{
	dwStyle |= WS_CHILD | WS_VISIBLE;
	return CreateWindowEx(
		0,
		TEXT("EDIT"),
		caption,
		dwStyle,
		rc.left,
		rc.top,
		rc.right - rc.left,
		rc.bottom - rc.top,
		hParent,
		(HMENU)((INT_PTR)(id)),
		hInst,
		NULL);
}
void GetString(char*str){
	RECT editrc_1 = { 120, 30, 220, 100 };
	HWND hEdit_1 = CreateEdit(
		NULL,
		(HINSTANCE)GetWindowLong(pWindow->hWnd, GWL_HINSTANCE),
		WS_POPUP | ES_LEFT ,
		editrc_1,
		201, ""   //文本框初始文字
		);
	GetWindowText(hEdit_1, str, 5);
}