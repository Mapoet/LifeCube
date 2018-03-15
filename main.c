
#include"base.h"
int WINAPI WinMain(HINSTANCE	hInstance,			// Instance
	HINSTANCE	hPrevInstance,		// Previous Instance
	LPSTR		lpCmdLine,			// Command Line Parameters
	int			nCmdShow)			// Window Show State
{
	MSG		msg;									// Windows Message Structure
	initGlobelView(hInstance, hPrevInstance, lpCmdLine, nCmdShow);
	createWindow(pWindow);
	while (!pWindow->isDone)									// Loop That Runs While done=FALSE
	{
		if (PeekMessage(&msg, NULL, 0, 0, PM_REMOVE))	// Is There A Message Waiting?
		{
			if (msg.message == WM_QUIT)				// Have We Received A Quit Message?
			{
				pWindow->isDone = TRUE;							// If So done=TRUE
			}
			else									// If Not, Deal With Window Messages
			{
				TranslateMessage(&msg);				// Translate The Message
				DispatchMessage(&msg);				// Dispatch The Message
			}
		}
		else										// If There Are No Messages
		{
				//// Draw The Scene.  Watch For ESC Key And Quit Messages From DrawGLScene()
			if (UpdateKeys(pControl, pWindow) == 1)
					return 0;
				//if (!isPause)
				//	UpdateData();
			if (pControl->isActive && pWindow->Draw(pWindow) || pControl->Keys[VK_ESCAPE])	// Active?  Was There A Quit Received?
			{
					pWindow->isDone = TRUE;							// ESC or DrawGLScene Signalled A Quit
			}
			else									// Not Time To Quit, Update Screen
			{
				SwapBuffers(pWindow->hDC);					// Swap Buffers (Double Buffering)
			}
			Sleep(50);
		}
	}
	destroyGlobelView();
	return 0;
}