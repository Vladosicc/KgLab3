#include "Render.h"

#include <future>
#include <sstream>
#include <iostream>
#include <fstream>

#include <windows.h>
#include <GL\GL.h>
#include <GL\GLU.h>

#include "MyOGL.h"

#include "Camera.h"
#include "Light.h"
#include "Primitives.h"

#include <array>

#include "GUItextRectangle.h"

#define PI 3.14159265

#define TOP_RIGHT 1.0f,1.0f
#define TOP_LEFT 0.0f,1.0f
#define BOTTOM_RIGHT 1.0f,0.0f
#define BOTTOM_LEFT 0.0f,0.0f

bool textureMode = true;
bool lightMode = true;
bool changeTexture = false;
bool alpha = false;
bool CTRLisDown = false;
bool EnableWriteInFile = true;
bool EnableReadInFile = true;

// Края поверхности
std::vector<std::vector<Vector3>> points(
	{
		{
			{0, 0, 3},
			{0, -1, 2},
			{0, -2, 2},
			{0, -3, 1}
		},
		{
			{1, 0, 1},
			{1, -1, 2},
			{1, -2, -2},
			{1, -3, 1}
		},
		{
			{2, 0, 1},
			{2, -1, 2},
			{2, -2, 2},
			{2, -3, 1}
		},
		{
			{3, 0, 1},
			{3, -1, 2},
			{3, -2, 2},
			{3, -3, 1}
		}
	});

//класс для настройки камеры
class CustomCamera : public Camera
{
public:
	//дистанция камеры
	double camDist;
	//углы поворота камеры
	double fi1, fi2;

	
	//значния масеры по умолчанию
	CustomCamera()
	{
		camDist = 15;
		fi1 = 1;
		fi2 = 1;
	}

	
	//считает позицию камеры, исходя из углов поворота, вызывается движком
	void SetUpCamera()
	{
		//отвечает за поворот камеры мышкой
		lookPoint.setCoords(0, 0, 0);

		pos.setCoords(camDist*cos(fi2)*cos(fi1),
			camDist*cos(fi2)*sin(fi1),
			camDist*sin(fi2));

		if (cos(fi2) <= 0)
			normal.setCoords(0, 0, -1);
		else
			normal.setCoords(0, 0, 1);

		LookAt();
	}

	void CustomCamera::LookAt()
	{
		//функция настройки камеры
		gluLookAt(pos.X(), pos.Y(), pos.Z(), lookPoint.X(), lookPoint.Y(), lookPoint.Z(), normal.X(), normal.Y(), normal.Z());
	}



}  camera;   //создаем объект камеры


//Класс для настройки света
class CustomLight : public Light
{
public:
	CustomLight()
	{
		//начальная позиция света
		pos = Vector3(1, 1, 3);
	}

	
	//рисует сферу и линии под источником света, вызывается движком
	void  DrawLightGhismo()
	{
		glDisable(GL_LIGHTING);

		
		glColor3d(0.9, 0.8, 0);
		Sphere s;
		s.pos = pos;
		s.scale = s.scale*0.08;
		s.Show();
		
		if (OpenGL::isKeyPressed('G'))
		{
			glColor3d(0, 0, 0);
			//линия от источника света до окружности
			glBegin(GL_LINES);
			glVertex3d(pos.X(), pos.Y(), pos.Z());
			glVertex3d(pos.X(), pos.Y(), 0);
			glEnd();

			//рисуем окруность
			Circle c;
			c.pos.setCoords(pos.X(), pos.Y(), 0);
			c.scale = c.scale*1.5;
			c.Show();
		}

	}

	void SetUpLight()
	{
		GLfloat amb[] = { 0.2, 0.2, 0.2, 0 };
		GLfloat dif[] = { 1.0, 1.0, 1.0, 0 };
		GLfloat spec[] = { .7, .7, .7, 0 };
		GLfloat position[] = { pos.X(), pos.Y(), pos.Z(), 1. };

		// параметры источника света
		glLightfv(GL_LIGHT0, GL_POSITION, position);
		// характеристики излучаемого света
		// фоновое освещение (рассеянный свет)
		glLightfv(GL_LIGHT0, GL_AMBIENT, amb);
		// диффузная составляющая света
		glLightfv(GL_LIGHT0, GL_DIFFUSE, dif);
		// зеркально отражаемая составляющая света
		glLightfv(GL_LIGHT0, GL_SPECULAR, spec);

		glEnable(GL_LIGHT0);
	}


} light;  //создаем источник света

void CheckClickPoint(LPPOINT poi, double dy, double dx) {
	GLint    viewport[4];    // параметры viewport-a.
	GLdouble projection[16]; // матрица проекции.
	GLdouble modelview[16];  // видовая матрица.

	glGetIntegerv(GL_VIEWPORT, viewport);           // узнаём параметры viewport-a.
	glGetDoublev(GL_PROJECTION_MATRIX, projection); // узнаём матрицу проекции.
	glGetDoublev(GL_MODELVIEW_MATRIX, modelview);   // узнаём видовую матрицу.
	double delta = 10;
	for (auto& v : points) 
	{
		for (auto& elem : v) 
		{
			double tempPoint[3];
			gluProject(elem.X(), elem.Y(), elem.Z(), modelview, projection, viewport, &tempPoint[0], &tempPoint[1], &tempPoint[2]); // Получаем экранные координаты обьекта
			if (tempPoint[0] > poi->x - delta && tempPoint[0] < poi->x + delta &&
				tempPoint[1] > poi->y - delta && tempPoint[1] < poi->y + delta) 
			{
				tempPoint[0] -= dx;
				tempPoint[1] += dy;
				gluUnProject(tempPoint[0], tempPoint[1], tempPoint[2], modelview, projection, viewport, elem.GetLinkX(), elem.GetLinkY(), elem.GetLinkZ()); // Переводим экранные координаты обратно в объектные
			}
		}
	}
}

//старые координаты мыши
int mouseX = 0, mouseY = 0;

void mouseEvent(OpenGL *ogl, int mX, int mY)
{
	int dx = mouseX - mX;
	int dy = mouseY - mY;
	mouseX = mX;
	mouseY = mY;

	//меняем углы камеры при нажатой левой кнопке мыши
	if (OpenGL::isKeyPressed(VK_RBUTTON))
	{
		camera.fi1 += 0.01*dx;
		camera.fi2 += -0.01*dy;
	}

	
	//двигаем свет по плоскости, в точку где мышь
	if (OpenGL::isKeyPressed('G') && !OpenGL::isKeyPressed(VK_LBUTTON))
	{
		LPPOINT POINT = new tagPOINT();
		GetCursorPos(POINT);
		ScreenToClient(ogl->getHwnd(), POINT);

		POINT->y = ogl->getHeight() - POINT->y;

		Ray r = camera.getLookRay(POINT->x, POINT->y);

		double z = light.pos.Z();

		double k = 0, x = 0, y = 0;
		if (r.direction.Z() == 0)
			k = 0;
		else
			k = (z - r.origin.Z()) / r.direction.Z();

		x = k*r.direction.X() + r.origin.X();
		y = k*r.direction.Y() + r.origin.Y();

		light.pos = Vector3(x, y, z);
	}

	if (OpenGL::isKeyPressed(VK_LBUTTON))
	{
		LPPOINT POINT = new tagPOINT();
		GetCursorPos(POINT);
		ScreenToClient(ogl->getHwnd(), POINT);
		POINT->y = ogl->getHeight() - POINT->y;
		CheckClickPoint(POINT, dy, dx);
	}

	if (OpenGL::isKeyPressed('G') && OpenGL::isKeyPressed(VK_LBUTTON))
	{
		light.pos = light.pos + Vector3(0, 0, 0.02*dy);
	}
	
}

void mouseWheelEvent(OpenGL *ogl, int delta)
{

	if (delta < 0 && camera.camDist <= 1)
		return;
	if (delta > 0 && camera.camDist >= 100)
		return;

	camera.camDist += 0.01*delta;

}

void SavePointsInTxt()
{
	std::ofstream inTXT("ConfigPointForPlane.txt");

	for (auto i : points)
	{
		for (auto item : i)
		{
			inTXT << item.ToString();
		}
	}

	inTXT.close();
}

void UploadPointsFromTxt()
{
	std::ifstream inTXT("ConfigPointForPlane.txt");
	std::string str;

	for (int i = 0; i < points.size(); i++)
	{
		for (int j = 0; j < points[0].size(); j++)
		{
			double cordtemp[3];
			for (int i = 0; i < 3; i++)
			{
				inTXT >> str;
				cordtemp[i] = atof(str.c_str());
			}
			points[i][j].setCoords(cordtemp[0], cordtemp[1], cordtemp[2]);
		}
	}
	inTXT.close();
}

void keyDownEvent(OpenGL *ogl, int key)
{
	if (key == 'L')
	{
		lightMode = !lightMode;
	}

	if (key == 'T')
	{
		textureMode = !textureMode;
	}

	if (key == 'Q')
	{
		changeTexture = !changeTexture;
	}

	if (key == 'B')
	{
		alpha = !alpha;
	}

	if (key == 'R')
	{
		camera.fi1 = 1;
		camera.fi2 = 1;
		camera.camDist = 15;

		light.pos = Vector3(1, 1, 3);
	}

	if (key == 'F')
	{
		light.pos = camera.pos;
	}

	if (key == 17)
	{
		CTRLisDown = true;
	}

	if (key == 'S' && CTRLisDown && EnableWriteInFile)
	{
		auto future = std::async(std::launch::async,SavePointsInTxt);
		EnableWriteInFile = false;
	}

	if (key == 'U' && CTRLisDown && EnableReadInFile)
	{
		UploadPointsFromTxt();
		EnableReadInFile = false;
	}
}

void keyUpEvent(OpenGL *ogl, int key)
{
	if (key == 17)
	{
		CTRLisDown = false;
	}

	if (key == 'S')
	{
		EnableWriteInFile = true;
	}

	if (key == 'U')
	{
		EnableReadInFile = true;
	}
}

GLuint texId[2];

void UploadTextureInTexId(const char* name, const int NumberOfTexId)
{
	//массив трехбайтных элементов  (R G B)
	RGBTRIPLE* texarray;
	
	//массив символов, (высота*ширина*4      4, потомучто   выше, мы указали использовать по 4 байта на пиксель текстуры - R G B A)
	char* texCharArray;
	int texW, texH;
	OpenGL::LoadBMP(name, &texW, &texH, &texarray);
	OpenGL::RGBtoChar(texarray, texW, texH, &texCharArray);


	//генерируем ИД для текстуры
	glGenTextures(1, &texId[NumberOfTexId]);
	//биндим айдишник, все что будет происходить с текстурой, будте происходить по этому ИД
	glBindTexture(GL_TEXTURE_2D, texId[NumberOfTexId]);

	//загружаем текстуру в видеопямять, в оперативке нам больше  она не нужна
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, texW, texH, 0, GL_RGBA, GL_UNSIGNED_BYTE, texCharArray);

	//отчистка памяти
	free(texCharArray);
	free(texarray);

	//наводим шмон
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
}

//выполняется перед первым рендером
void initRender(OpenGL *ogl)
{
	//настройка текстур

	//4 байта на хранение пикселя
	glPixelStorei(GL_UNPACK_ALIGNMENT, 4);

	//настройка режима наложения текстур
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

	//включаем текстуры
	glEnable(GL_TEXTURE_2D);
	
	UploadTextureInTexId("texture.bmp", 0);
	UploadTextureInTexId("texture2.bmp", 1);

	//камеру и свет привязываем к "движку"
	ogl->mainCamera = &camera;
	ogl->mainLight = &light;

	// нормализация нормалей : их длины будет равна 1
	glEnable(GL_NORMALIZE);

	// устранение ступенчатости для линий
	glEnable(GL_LINE_SMOOTH); 


	//   задать параметры освещения
	//  параметр GL_LIGHT_MODEL_TWO_SIDE - 
	//                0 -  лицевые и изнаночные рисуются одинаково(по умолчанию), 
	//                1 - лицевые и изнаночные обрабатываются разными режимами       
	//                соответственно лицевым и изнаночным свойствам материалов.    
	//  параметр GL_LIGHT_MODEL_AMBIENT - задать фоновое освещение, 
	//                не зависящее от сточников
	// по умолчанию (0.2, 0.2, 0.2, 1.0)

	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 0);

	camera.fi1 = -1.3;
	camera.fi2 = 0.8;
}

double* ConvertInPoint(double x, double y, double z)
{
	double Point[] = { x,y,z };
	return Point;
}

void NormalizeVector(double* vec)
{
	double modVector = -sqrt(pow(vec[0], 2) + pow(vec[1], 2) + pow(vec[2], 2));

	for (int i = 0; i < 3; ++i)
	{
		vec[i] /= modVector;
	}

}


double* FindNormal(double x, double y, double z, double x1, double y1, double z1, double x2, double y2, double z2, int FlagSwap = 0) // A - "общая"
{
	double vectorA[3], vectorB[3];
	double a[3] = { x,y,z }, b[3] = { x1,y1,z1 }, c[3] = { x2,y2,z2 };

	for (int i = 0; i < 3; ++i) // Получаем вектор A и B
	{
		vectorA[i] = a[i] - c[i];
		vectorB[i] = b[i] - c[i];
	}

	double VectorNormal[3];

	//VectorNormal[0] = a[1] * b[2] - a[2] * b[1];
	//VectorNormal[1] = a[2] * b[0] - a[0] * b[2];
	//VectorNormal[2] = a[0] * b[1] - a[1] * b[0];

	VectorNormal[0] = vectorA[1] * vectorB[2] - vectorB[1] * vectorA[2];
	VectorNormal[1] = -vectorA[0] * vectorB[2] + vectorB[0] * vectorA[2];
	VectorNormal[2] = vectorA[0] * vectorB[1] - vectorB[0] * vectorA[1];

	NormalizeVector(VectorNormal);

	if (FlagSwap != 0)
	{
		for (int i = 0; i < 3; ++i) // Получаем вектор A и B
		{
			VectorNormal[i] *= -1;
		}
	}

	return VectorNormal;
}

double* FindNormal(const double* a,const double* b,const double* c, int FlagSwap = 0) // A - "общая"
{
	double vectorA[3], vectorB[3];

	for (int i = 0; i < 3; ++i) // Получаем вектор A и B
	{
		vectorA[i] = a[i] - c[i];
		vectorB[i] = b[i] - c[i];
	}

	double VectorNormal[3];

	//VectorNormal[0] = a[1] * b[2] - a[2] * b[1];
	//VectorNormal[1] = a[2] * b[0] - a[0] * b[2];
	//VectorNormal[2] = a[0] * b[1] - a[1] * b[0];

	VectorNormal[0] = vectorA[1] * vectorB[2] - vectorB[1] * vectorA[2];
	VectorNormal[1] = -vectorA[0] * vectorB[2] + vectorB[0] * vectorA[2];
	VectorNormal[2] = vectorA[0] * vectorB[1] - vectorB[0] * vectorA[1];

	NormalizeVector(VectorNormal);

	if (FlagSwap != 0)
	{
		for (int i = 0; i < 3; ++i) // Получаем вектор A и B
		{
			VectorNormal[i] *= -1;
		}
	}

	return VectorNormal;
}


double NewCoordForTexture(double oldCoord, double maxCoord)
{
	return oldCoord / maxCoord;
}

double f(double p1, double p2, double p3, double t)
{
	return p1 * (1 - t) * (1 - t) * (1 - t) + 2 * p2 * t * (1 - t) + p3 * t * t; //посчитанная формула
}

double f(double a, double b, double t)
{
	return a * (1 - t) + b * t;
}

double f(double p1, double p2, double p3, double p4, double t)
{
	return p1 * (1 - t) * (1 - t) * (1 - t) + 3 * p2 * t * (1 - t) * (1 - t) + 3 * p3 * t * t * (1 - t) + p4 * t * t * t; //посчитанная формула
}

double fE(double p1, double p4, double vec1, double vec4, double t)
{
	return p1 * (2 * t * t * t - 3 * t * t + 1) + p4 * (-2 * t * t * t + 3 * t * t) + vec1 * (t * t * t - 2 * t * t + t) + vec4 * (t * t * t - t * t); //посчитанная формула
}

void bize(double P1[3], double P2[3], double P3[3], double P4[3])
{
	glDisable(GL_LIGHTING);
	glColor3d(1, 0.5, 0);
	glLineWidth(3); //ширина линии
	glBegin(GL_LINE_STRIP);
	for (double t = 0; t <= 1.0001; t += 0.01)
	{
		double P[3];
		P[0] = f(P1[0], P2[0], P3[0], P4[0], t);
		P[1] = f(P1[1], P2[1], P3[1], P4[1], t);
		P[2] = f(P1[2], P2[2], P3[2], P4[2], t);
		glVertex3dv(P); //Рисуем точку P
	}
	glEnd();
	glColor3d(0, 0, 1);
	glLineWidth(1); //возвращаем ширину линии = 1

	glBegin(GL_LINES);
	glVertex3dv(P1);
	glVertex3dv(P2);
	glVertex3dv(P2);
	glVertex3dv(P3);
	glVertex3dv(P3);
	glVertex3dv(P4);
	glEnd();
	glEnable(GL_LIGHTING);
}

Vector3 bizeWithoutDraw(double P1[3], double P2[3], double P3[3], double P4[3], double t)
{
	Vector3 Vec;
	Vec.setCoords(f(P1[0], P2[0], P3[0], P4[0], t), f(P1[1], P2[1], P3[1], P4[1], t), f(P1[2], P2[2], P3[2], P4[2], t));
	return Vec;
}

void Draw_Cube(double P1[3])
{
	//glBegin(GL_QUADS);

	////низ
	//glColor3d(0.7, 0.7, 0.7);
	//glVertex3d(P1[0], P1[1], P1[2]);
	//glVertex3d(P1[0] + 1, P1[1], P1[2]);
	//glVertex3d(P1[0] + 1, P1[1] + 1, P1[2]);
	//glVertex3d(P1[0], P1[1] + 1, P1[2]);

	////верх
	//glColor3d(0.2, 0.2, 0.2);
	//glVertex3d(P1[0], P1[1], P1[2] + 1);
	//glVertex3d(P1[0] + 1, P1[1], P1[2] + 1);
	//glVertex3d(P1[0] + 1, P1[1] + 1, P1[2] + 1);
	//glVertex3d(P1[0], P1[1] + 1, P1[2] + 1);


	//glColor3d(0.3, 0.3, 0.3);
	//glVertex3d(P1[0], P1[1], P1[2]);
	//glVertex3d(P1[0], P1[1] + 1, P1[2]);
	//glVertex3d(P1[0], P1[1] + 1, P1[2] + 1);
	//glVertex3d(P1[0], P1[1], P1[2] + 1);

	//glColor3d(0.4, 0.4, 0.4);
	//glVertex3d(P1[0] + 1, P1[1], P1[2]);
	//glVertex3d(P1[0] + 1, P1[1] + 1, P1[2]);
	//glVertex3d(P1[0] + 1, P1[1] + 1, P1[2] + 1);
	//glVertex3d(P1[0] + 1, P1[1], P1[2] + 1);

	//glColor3d(0.5, 0.5, 0.5);
	//glVertex3d(P1[0], P1[1], P1[2]);
	//glVertex3d(P1[0] + 1, P1[1], P1[2]);
	//glVertex3d(P1[0] + 1, P1[1], P1[2] + 1);
	//glVertex3d(P1[0], P1[1], P1[2] + 1);

	//glColor3d(0.6, 0.6, 0.6);
	//glVertex3d(P1[0], P1[1] + 1, P1[2]);
	//glVertex3d(P1[0] + 1, P1[1] + 1, P1[2]);
	//glVertex3d(P1[0] + 1, P1[1] + 1, P1[2] + 1);
	//glVertex3d(P1[0], P1[1] + 1, P1[2] + 1);
	//glEnd();

	//glBegin(GL_TRIANGLES);
	//glVertex3d(P1[0] + 0.5, P1[1] + 0.5, P1[2] + 2);
	//glVertex3d(P1[0], P1[1], P1[2] + 1);
	//glVertex3d(P1[0] + 1, P1[1], P1[2] + 1);

	//glVertex3d(P1[0] + 0.5, P1[1] + 0.5, P1[2] + 2);
	//glVertex3d(P1[0] + 1, P1[1], P1[2] + 1);
	//glVertex3d(P1[0] + 1, P1[1] + 1, P1[2] + 1);

	//glVertex3d(P1[0] + 0.5, P1[1] + 0.5, P1[2] + 2);
	//glVertex3d(P1[0] + 1, P1[1] + 1, P1[2] + 1);
	//glVertex3d(P1[0], P1[1] + 1, P1[2] + 1);

	//glVertex3d(P1[0] + 0.5, P1[1] + 0.5, P1[2] + 2);
	//glVertex3d(P1[0], P1[1] + 1, P1[2] + 1);
	//glVertex3d(P1[0], P1[1], P1[2] + 1);
	//glEnd();

	//glBegin(GL_LINES);
	//glColor3d(1, 0, 0);
	//glVertex3d(0, 0, 0);
	//glVertex3d(1, 0, 0);

	//glColor3d(0, 1, 0);
	//glVertex3d(0, 0, 0);
	//glVertex3d(0, 1, 0);

	//glColor3d(0, 0, 1);
	//glVertex3d(0, 0, 0);
	//glVertex3d(0, 0, 1);
	//glEnd();
	glPushMatrix();
	glRotated(90, 0, 0, 1);
	glRotated(90, 1, 0, 0);


	glBegin(GL_TRIANGLES);
	// Конус носа /////////////////////////////
	  // Белый
	glColor3ub(255, 255, 255);
	glNormal3dv(FindNormal(0.0f, 0.0f, 0.600f, -0.150f, 0.0f, 0.300f, 0.150f, 0.0f, 0.300f, 1));
	glTexCoord2d(0, 0); glVertex3f(0.0f, 0.0f, 0.600f);
	glTexCoord2d(0.5, 1); glVertex3f(-0.150f, 0.0f, 0.300f);
	glTexCoord2d(1, 0); glVertex3f(0.150f, 0.0f, 0.300f);

	// Черный
	glColor3ub(0, 0, 0);
	glNormal3dv(FindNormal(0.150f, 0.0f, 0.300f, 0.0f, 0.150f, 0.300f, 0.0f, 0.0f, 0.600f, 1));
	glTexCoord2d(0, 0); glVertex3f(0.150f, 0.0f, 0.300f);
	glTexCoord2d(0.5, 1); glVertex3f(0.0f, 0.150f, 0.300f);
	glTexCoord2d(1, 0); glVertex3f(0.0f, 0.0f, 0.600f);

	// Красный
	glColor3ub(255, 0, 0);
	glNormal3dv(FindNormal(0.0f, 0.0f, 0.600f, 0.0f, 0.150f, 0.300f, -0.150f, 0.0f, 0.300f, 1));
	glTexCoord2d(0, 0); glVertex3f(0.0f, 0.0f, 0.600f);
	glTexCoord2d(0.5, 1); glVertex3f(0.0f, 0.150f, 0.300f);
	glTexCoord2d(1, 0); glVertex3f(-0.150f, 0.0f, 0.300f);

	// Тело самолета //////////////////////////
	  // Зеленый
	glColor3ub(0, 255, 0);
	glNormal3dv(FindNormal(-0.150f, 0.0f, 0.300f, 0.0f, 0.150f, 0.300f, 0.0f, 0.0f, -0.560f, 1));
	glTexCoord2d(0, 0); glVertex3f(-0.150f, 0.0f, 0.300f);
	glTexCoord2d(0.5, 1); glVertex3f(0.0f, 0.150f, 0.300f);
	glTexCoord2d(1, 0); glVertex3f(0.0f, 0.0f, -0.560f);

	// Желтый
	glColor3ub(255, 255, 0);
	glNormal3dv(FindNormal(0.0f, 0.0f, -0.560f, 0.0f, 0.15f, 0.30f, 0.15f, 0.0f, 0.30f, 1));
	glTexCoord2d(0, 0); glVertex3f(0.0f, 0.0f, -0.560f);
	glTexCoord2d(0.5, 1); glVertex3f(0.0f, 0.15f, 0.30f);
	glTexCoord2d(1, 0); glVertex3f(0.15f, 0.0f, 0.30f);

	// Голубой
	glColor3ub(0, 255, 255);
	glNormal3dv(FindNormal(0.15f, 0.0f, 0.30f, -0.15f, 0.0f, 0.30f, 0.0f, 0.0f, -0.560f, 1));
	glTexCoord2d(0, 0); glVertex3f(0.15f, 0.0f, 0.30f);
	glTexCoord2d(0.5, 1); glVertex3f(-0.15f, 0.0f, 0.30f);
	glTexCoord2d(1, 0); glVertex3f(0.0f, 0.0f, -0.560f);

	///////////////////////////////////////////
	// Левое крыло
	// Большой треугольник для основания крыла
	  // Серый
	glColor3ub(128, 128, 128);
	glNormal3dv(FindNormal(0.0f, 0.02f, 0.270f, -0.60f, 0.02f, -0.08f, 0.60f, 0.02f, -0.08f, 1));
	glTexCoord2d(0, 0); glVertex3f(0.0f, 0.02f, 0.270f);
	glTexCoord2d(0.5, 1); glVertex3f(-0.60f, 0.02f, -0.08f);
	glTexCoord2d(1, 0); glVertex3f(0.60f, 0.02f, -0.08f);

	// Темно - серый
	glColor3ub(64, 64, 64);
	glNormal3dv(FindNormal(0.60f, 0.02f, -0.08f, 0.0f, 0.07f, -0.08f, 0.0f, 0.02f, 0.270f, 1));
	glTexCoord2d(0, 0); glVertex3f(0.60f, 0.02f, -0.08f);
	glTexCoord2d(0.5, 1); glVertex3f(0.0f, 0.07f, -0.08f);
	glTexCoord2d(1, 0); glVertex3f(0.0f, 0.02f, 0.270f);

	// Светло - серый
	glColor3ub(192, 192, 192);
	glNormal3dv(FindNormal(0.60f, 0.02f, -0.08f, -0.60f, 0.02f, -0.08f, 0.0f, 0.07f, -0.08f, 1));
	glTexCoord2d(0, 0); glVertex3f(0.60f, 0.02f, -0.08f);
	glTexCoord2d(0.5, 1); glVertex3f(-0.60f, 0.02f, -0.08f);
	glTexCoord2d(1, 0); glVertex3f(0.0f, 0.07f, -0.08f);

	// Другое крыло верхней секции
	  // Темно - серый
	glColor3ub(64, 64, 64);
	glNormal3dv(FindNormal(0.0f, 0.02f, 0.270f, 0.0f, 0.07f, -0.08f, -0.60f, 0.02f, -0.08f, 1));
	glTexCoord2d(0, 0); glVertex3f(0.0f, 0.02f, 0.270f);
	glTexCoord2d(0.5, 1); glVertex3f(0.0f, 0.07f, -0.08f);
	glTexCoord2d(1, 0); glVertex3f(-0.60f, 0.02f, -0.08f);

	///////////////////////////////////////////
	// Хвост 
	  // Розовый
	glColor3ub(255, 128, 255);
	glNormal3dv(FindNormal(-0.30f, 0.00f, -0.76f, 0.30f, 0.00f, -0.76f, 0.0f, 0.0f, -0.53f));
	glTexCoord2d(0, 0); glVertex3f(-0.30f, 0.00f, -0.76f);
	glTexCoord2d(0.5, 1); glVertex3f(0.30f, 0.00f, -0.76f);
	glTexCoord2d(1, 0); glVertex3f(0.0f, 0.0f, -0.53f);

	glColor3ub(155, 128, 255);
	glNormal3dv(FindNormal(0.0f, 0.0f, -0.56f, 0.05f, 0.00f, -0.76f, 0.0f, 0.20f, -0.76f, 1));
	glTexCoord2d(0, 0); glVertex3f(0.0f, 0.0f, -0.56f);
	glTexCoord2d(0.5, 1); glVertex3f(0.05f, 0.00f, -0.76f);
	glTexCoord2d(1, 0); glVertex3f(0.0f, 0.20f, -0.76f);

	glColor3ub(155, 128, 255);
	glNormal3dv(FindNormal(0.0f, 0.0f, -0.56f, -0.05f, 0.00f, -0.76f, 0.0f, 0.20f, -0.76f));
	glTexCoord2d(0, 0); glVertex3f(0.0f, 0.0f, -0.56f);
	glTexCoord2d(0.5, 1); glVertex3f(-0.05f, 0.00f, -0.76f);
	glTexCoord2d(1, 0); glVertex3f(0.0f, 0.20f, -0.76f);

	glColor3ub(155, 128, 255);
	glNormal3dv(FindNormal(0.05f, 0.00f, -0.76f, -0.05f, 0.00f, -0.76f, 0.0f, 0.20f, -0.76f, 1));
	glTexCoord2d(0, 0); glVertex3f(0.05f, 0.00f, -0.76f);
	glTexCoord2d(0.5, 1); glVertex3f(-0.05f, 0.00f, -0.76f);
	glTexCoord2d(1, 0); glVertex3f(0.0f, 0.20f, -0.76f);

	glEnd();
	glPopMatrix();
}

void Bese2(double P1[3], double P2[3], double P3[3], double P4[3], double delta_time)
{
	static double t_max = 0;
	static bool flagReverse = false;

	if (!flagReverse)
	{
		t_max += delta_time / 5; //t_max становится = 1 за 5 секунд
		if (t_max > 1)
		{
			t_max = 1; //после обнуляется
			flagReverse = !flagReverse;
		}
	}
	else
	{
		t_max -= delta_time / 5; //t_max становится = 1 за 5 секунд
		if (t_max < 0)
		{
			t_max = 0; //после обнуляется
			flagReverse = !flagReverse;
		}
	}

	bize(P1, P2, P3, P4);

	Vector3 P_old = bizeWithoutDraw(P1, P2, P3, P4, !flagReverse ? t_max - delta_time : t_max + delta_time);
	Vector3 P = bizeWithoutDraw(P1, P2, P3, P4, t_max);
	Vector3 VecP_P_old = (P - P_old).normolize();

	Vector3 rotateX(VecP_P_old.X(), VecP_P_old.Y(), 0);
	rotateX = rotateX.normolize();

	Vector3 VecPrX = Vector3(1, 0, 0).vectProisvedenie(rotateX);
	double CosX = Vector3(1, 0, 0).ScalarProizv(rotateX);
	double SinAngleZ = VecPrX.Z() / abs(VecPrX.Z());
	double AngleOZ = acos(CosX) * 180 / PI * SinAngleZ;

	double AngleOY = acos(VecP_P_old.Z()) * 180 / PI - 90;

	double A[] = { -0.5,-0.5,-0.5 };
	glPushMatrix();
	glTranslated(P.X(), P.Y(), P.Z());
	glRotated(AngleOZ, 0, 0, 1);
	glRotated(AngleOY, 0, 1, 0);
	//glRotated(CalculateAngleOZ(P_old, P), 0, 0, 1);
	Draw_Cube(A);
	glPopMatrix();

	glColor3d(0, 0, 0);

}


void Ermit(double P1[3], double P2[3], double P3[3], double P4[3])
{
	double vec1[] = { P2[0] - P1[0], P2[1] - P1[1], P2[2] - P1[2] };
	double vec4[] = { P4[0] - P3[0], P4[1] - P3[1], P4[2] - P3[2] };

	glLineWidth(3); //ширина линии
	glBegin(GL_LINE_STRIP);
	for (double t = 0; t <= 1.0001; t += 0.01)
	{
		double P[3];
		P[0] = fE(P1[0], P3[0], vec1[0], vec4[0], t);
		P[1] = fE(P1[1], P3[1], vec1[1], vec4[1], t);
		P[2] = fE(P1[2], P3[2], vec1[2], vec4[2], t);
		glVertex3dv(P); //Рисуем точку P
	}
	glEnd();
	glLineWidth(1); //возвращаем ширину линии = 1

	glBegin(GL_LINES);
	glVertex3dv(P1);
	glVertex3dv(P2);
	glVertex3dv(P4);
	glVertex3dv(P3);
	glEnd();

	glColor3d(1, 0, 1);
	glPointSize(10);
	glBegin(GL_POINTS);
	glVertex3dv(P1);
	glVertex3dv(P3);
	glEnd();
	glColor3d(0, 0, 0);
}

// Тут поверхность 

unsigned long int factorial(int i) 
{
	if (i == 0) return 1;
	else return i * factorial(i - 1);
}

double Bernstein(double u, double n, int index) {
	return (factorial(n) / (factorial(index) * factorial(n - index))) * pow(u, index) * pow(1 - u, n - index);
}

void BezeSurfacePoint(double u, double v, Vector3& vec) {
	Vector3 new_v;
	int n = 3, m = 3;
	for (size_t i = 0; i < points.size(); ++i) {
		for (size_t j = 0; j < points[i].size(); ++j) 
		{
			new_v += points[i][j] * Bernstein(u, n, i) * Bernstein(v, m, j);
		}
	}
	vec = new_v;
}

std::array<Vector3, 4> Triangle(double u, double v, double h) {
	std::array<Vector3, 4> tmp;
	BezeSurfacePoint(u, v, tmp[0]);
	BezeSurfacePoint(u, v + h, tmp[1]);
	BezeSurfacePoint(u + h, v, tmp[2]);
	BezeSurfacePoint(u + h, v + h, tmp[3]);

	return tmp;
}

void Lines_Points() {
	glDisable(GL_LIGHTING);
	glColor3d(0.5, 0.3, 0.7);
	glLineWidth(0.5);
	for (int i = 0; i < points.size(); i++) {
		glBegin(GL_LINE_STRIP);
		for (int j = 0; j < points[i].size(); j++)
			glVertex3dv(points[i][j].toArray());
		glEnd();
	}
	for (int i = 0; i < points.size(); i++) {
		glBegin(GL_LINE_STRIP);
		for (int j = 0; j < points[i].size(); j++)
			glVertex3dv(points[j][i].toArray());
		glEnd();
	}
	glPointSize(10);
	glColor3d(0, 1, 0);
	glBegin(GL_POINTS);
	for (auto& i : points) {
		for (auto& elem : i)
			glVertex3dv(elem.toArray());
	}
	glEnd();
	glEnable(GL_LIGHTING);
}

double MaxCoordX()
{
	double max = points[0][0].X();
	double min = max;
	for (auto i : points) {
		for (auto elem : i)
		{
			if (elem.X() > max)
			{
				max = elem.X();
			}
			if (elem.X() < min)
			{
				min = elem.X();
			}
		}

	}
	return max - min;
}

double MaxCoordY()
{
	double max = points[0][0].Y();
	double min = max;
	for (auto i : points) {
		for (auto elem : i)
		{
			if (elem.Y() > max)
			{
				max = elem.Y();
			}
			if (elem.Y() < min)
			{
				min = elem.Y();
			}
		}
	}
	return max - min;
}

double *TextureCoords(Vector3 vec) 
{
	double coords[3];
	coords[0] = vec.X()/MaxCoordX();
	coords[1] = vec.Y()/MaxCoordY();
	coords[2] = 0;
	return coords;
}

void BezePlane() {
	double h = 0.1;
	for (double u = 0; u < 1.01 - h; u += h) {
		glBegin(GL_TRIANGLES);
		for (double v = 0; v < 1.01 - h; v += h) {
			std::array<Vector3, 4> arr = Triangle(u, v, h);

			glNormal3dv(FindNormal(arr[0].toArray(), arr[1].toArray(), arr[2].toArray(), 1));
			glTexCoord3dv(TextureCoords(arr[0]));
			glVertex3dv(arr[0].toArray());
			glTexCoord3dv(TextureCoords(arr[1]));
			glVertex3dv(arr[1].toArray());
			glTexCoord3dv(TextureCoords(arr[3]));
			glVertex3dv(arr[3].toArray());

			glNormal3dv(FindNormal(arr[0].toArray(), arr[1].toArray(), arr[2].toArray(), 1));
			glTexCoord3dv(TextureCoords(arr[0]));
			glVertex3dv(arr[0].toArray());
			glTexCoord3dv(TextureCoords(arr[2]));
			glVertex3dv(arr[2].toArray());
			glTexCoord3dv(TextureCoords(arr[3]));
			glVertex3dv(arr[3].toArray());
		}
		glEnd();
	}
	Lines_Points();
}

double DeltaTime = 0.01;

void Render(OpenGL *ogl)
{
	glBindTexture(GL_TEXTURE_2D, texId[0]);

	glDisable(GL_TEXTURE_2D);
	glDisable(GL_LIGHTING);

	glEnable(GL_DEPTH_TEST);
	if (textureMode)
		glEnable(GL_TEXTURE_2D);

	if (lightMode)
		glEnable(GL_LIGHTING);

	if (changeTexture)
		glBindTexture(GL_TEXTURE_2D, texId[1]);


	//альфаналожение
	//glEnable(GL_BLEND);
	//glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);


	//настройка материала
	GLfloat amb[] = { 0.3, 0.3, 0.3, 0.3 };
	GLfloat dif[] = { 0.7, 0.7, 0.7, 0.4 };
	GLfloat spec[] = { 0.9, 0.8, 0.4, 0.5 };
	GLfloat sh = 0.1f * 256;


	//фоновая
	glMaterialfv(GL_FRONT, GL_AMBIENT, amb);
	//дифузная
	glMaterialfv(GL_FRONT, GL_DIFFUSE, dif);
	//зеркальная
	glMaterialfv(GL_FRONT, GL_SPECULAR, spec); \
		//размер блика
		glMaterialf(GL_FRONT, GL_SHININESS, sh);

	//чтоб было красиво, без квадратиков (сглаживание освещения)
	glShadeModel(GL_SMOOTH);
	//===================================
	//Прогать тут  

	double P1[] = { 0,0,0 }; //Наши точки
	double P2[] = { -4,6,7 };
	double P3[] = { 5,0,7 };
	double P4[] = { 10,10,0 };

	Bese2(P1, P2, P3, P4, DeltaTime);

	//double A[] = { -0.5,-0.5,-0.5 };
	//Draw_Cube(A);

	BezePlane();

   //Сообщение вверху экрана

	
	glMatrixMode(GL_PROJECTION);	//Делаем активной матрицу проекций. 
	                                //(всек матричные операции, будут ее видоизменять.)
	glPushMatrix();   //сохраняем текущую матрицу проецирования (которая описывает перспективную проекцию) в стек 				    
	glLoadIdentity();	  //Загружаем единичную матрицу
	glOrtho(0, ogl->getWidth(), 0, ogl->getHeight(), 0, 1);	 //врубаем режим ортогональной проекции

	glMatrixMode(GL_MODELVIEW);		//переключаемся на модел-вью матрицу
	glPushMatrix();			  //сохраняем текущую матрицу в стек (положение камеры, фактически)
	glLoadIdentity();		  //сбрасываем ее в дефолт

	glDisable(GL_LIGHTING);



	GuiTextRectangle rec;		   //классик моего авторства для удобной работы с рендером текста.
	rec.setSize(300, 200);
	rec.setPosition(10, ogl->getHeight() - 200 - 10);


	std::stringstream ss;
	ss << "T - вкл/выкл текстур" << std::endl;
	ss << "Q - смена текстур" << std::endl;
	ss << "B - вкл/выкл альфа-наложение" << std::endl;
	ss << "L - вкл/выкл освещение" << std::endl;
	ss << "F - Свет из камеры" << std::endl;
	ss << "G - двигать свет по горизонтали" << std::endl;
	ss << "G+ЛКМ двигать свет по вертекали" << std::endl;
	ss << "CTRL+S сохранить текущее состояние поверхности" << std::endl;
	ss << "CTRL+U загрузить поверхность из файла" << std::endl;
	ss << "Коорд. света: (" << light.pos.X() << ", " << light.pos.Y() << ", " << light.pos.Z() << ")" << std::endl;
	ss << "Коорд. камеры: (" << camera.pos.X() << ", " << camera.pos.Y() << ", " << camera.pos.Z() << ")" << std::endl;
	ss << "Параметры камеры: R="  << camera.camDist << ", fi1=" << camera.fi1 << ", fi2=" << camera.fi2 << std::endl;
	
	rec.setText(ss.str().c_str());
	rec.Draw();

	glMatrixMode(GL_PROJECTION);	  //восстанавливаем матрицы проекции и модел-вью обратьно из стека.
	glPopMatrix();

	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();


}