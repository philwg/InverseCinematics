#define ARMA_USE_LAPACK				// La fonction arma::pinv() s'appuie sur LAPACK ...
#define ARMA_DONT_USE_WRAPPER		// ... mais sans WRAPPER !
#include <armadillo>				// pour les calculs matriciels
#include <GL/glut.h>				// pour l'affichage OpenGL

#define PI 3.14159265359f			// Pour les angles en radians

using namespace std;
using namespace arma;

/** Déclaration des variables ****************************************/

const int N  {10};			// Le nombre de bras du mobile articulé

vec Thetas;					// Les angles associés à chaque bras
vec Longueurs;				// Les longueurs de chaque bras
vec prev_Cible;				// La dernière cible

mat Pivots;					// Les points d'articulation

int width 	{1200};			// La largeur de la fenêtre au départ
int height 	{1200};			// La hauteur de la fenêtre au départ
int xmin 	{0};			// Pour le positionnement du premier Pivot en x ...
int ymin 	{0};			// ... et en y

/** Prototypes des méthodes ******************************************/

void drawObject();									// Dessin du bras articulé
void computePivotCoordinates();						// (Re-)Calcul des coordonnées des Pivots
float norm(vec P);									// Calcul de la norme d'un vecteur
void computeCoordinates(vec new_Cible);				// Calcul des variations et mise à jour
vec get_New_Cible(int x, int y);					// Renvoie une cible valable à partir des entiers
int reverse_Cible_X(vec C);							// Retrouve le x entier à partir des coordonnées
int reverse_Cible_Y(vec C);							// Retrouve le y entier à partir des coordonnées
void myinit();										// Initialisation des données

// Méthodes d'interaction GLUT
void mouse(int button, int state, int x, int y);	// gestion des clics souris
void motion(int x, int y);							// gestion des mouvements de la souris
void parsekey(unsigned char key, int x, int y);		// gestion des touches du clavier
void parsekey_special(int key, int x, int y);		// gestion des touches speciales
void myReshape(int w, int h);						// gestion du redimensionnement
void display();										// gestion de l'affichage

/** Implémentations **************************************************/

// Tracé du bras articulé en empilant les transformations dans OpenGL
void drawObject()
{
	glPointSize(8.0f);
	glLineWidth(5.0f);
	glPushMatrix();
		glTranslatef(0.0f, -0.9f, 0.0f);
		glColor3f(1.0f, 0.0f, 0.0f);
		glBegin(GL_LINES);
			glVertex2f(-1.0f, 0.0f);
			glVertex2f(1.0f, 0.0f);
		glEnd();
		glRotatef(-90.0f, 0.0f, 0.0f, 1.0f);
		glBegin(GL_POINTS);
			glVertex2f(0.0f, 0.0f);
		glEnd();
		for (int i=0; i<N; i++) {
    		glColor3f(0.0f, 1.0f, 0.0f);
    		glRotatef((180.0f/PI)*Thetas(i), 0.0f, 0.0f, 1.0f);
    		glBegin(GL_LINES);
    			glVertex2f(0.0f, 0.0f);
    			glVertex2f(0.0f, Longueurs(i));
    		glEnd();
		    glColor3f(1.0f, 0.0f, 0.0f);
 			glTranslatef(0.0f, Longueurs(i), 0.0f);
    		glBegin(GL_POINTS);
    			glVertex2f(0.0f, 0.0f);
    		glEnd();
  		}
  	glPopMatrix();
}

// Calcul des pivots
void computePivotCoordinates()
{
	float angle;
	angle = 0.0f;
	Pivots(0, 0) = 0.0f;
	Pivots(0, 1) = 0.0f;
	for (int i=1; i<=N; i++) {
		// Attention on cumule les angles pour tenir compte de toute la chaîne
		angle = angle + Thetas(i-1);	
		Pivots(i, 0) = Pivots(i-1, 0) + Longueurs(i-1)*cos(angle);
		Pivots(i, 1) = Pivots(i-1, 1) + Longueurs(i-1)*sin(angle);
	}
}

// Calcul de la norme euclidienne d'un vecteur
float norm(vec V)
{
	return sqrt(V(0)*V(0)+V(1)*V(1));
}

// Méthode de recalcul des paramètres pour atteindre un cible
void computeCoordinates(vec Cible)
{
	mat JacInverse = zeros(2, N);
	mat Jacobien = zeros(N, 2);
	vec lambdas = zeros(N);
	vec E = zeros(2);

	// On évalue la différence de position */
	E(0) = Cible(0) - Pivots(N, 0);
	E(1) = Cible(0) - Pivots(N, 1);
	
	// Tant que la différence est grande, on continue à chercher ...
	while (norm(E) >= 0.01f) {

		// Calcul du Jacobien
		for (int i=0; i<N; i++) {
			Jacobien(i, 0) = -(Pivots(N, 1) - Pivots(i, 1));
			Jacobien(i, 1) =  Pivots(N, 0) - Pivots(i, 0);
		}

		// Calcul de son pseudo Inverse
		JacInverse = pinv(Jacobien);

		// Calcul des variations d'angles en découlant
		float lbd_Max = 0.0f;
		for (int i=0; i<N ; i++) {
			lambdas(i) = JacInverse(0, i)*E(0) + JacInverse(1, i)*E(1);
			// Avec mise à jour de leur maximum
			if (abs(lambdas(i))>lbd_Max) lbd_Max = abs(lambdas(i));
		}

		// si le maximum est supérieur à 2°, on reduit l'ensemble proportionnellement
		if ((lbd_Max*180) >= 2) {
			for (int i=0; i<N; i++) {
				lambdas(i) = (lambdas(i)/lbd_Max)*(2.0f/180.0f);
			}
		}

		// On recalcule les nouveaux angles modifiés et les nouveaux pivots induits
		for (int i=0; i<N; i++) {
			Thetas(i) = Thetas(i) + lambdas(i);
		}
		computePivotCoordinates();

		// On met à jour l'évaluation de l'erreur
		E(0) = Cible(0) - Pivots(N, 0);
		E(1) = Cible(1) - Pivots(N, 1);
	}

	// On conserve la dernière position licite obtenue
	prev_Cible = Cible;
}

/** Méthode de contrôle  : Pour ne pas faire diverger l'algorithme,
/*  il faut que la cible soit à une distance inférieure à la longueur
/*  totale du bras par rapport à son point d'ancrage.
*/
vec get_New_Cible(int x, int y)
{
	vec cible = zeros(2);		
	cible(0) = -1.0f + 2.0f*(float)((x-xmin))/(float)width;
	cible(1) = 1.9f - 2.0f*(float)((y-ymin))/(float)height;
	if (norm(cible)>1.5f) {
		// Si la condition n'est pas remplie, on reste à la dernière position
		cible(0) = prev_Cible(0);
		cible(1) = prev_Cible(1);
	}
	return cible;
}

// On récupère l'abscisse écran entière d'un vecteur
int reverse_Cible_X(vec C)
{
	return round(((C(0)+1.0f)*((float)width)/2.0f)) + xmin;
}

// On récupère l'ordonnée écran entière d'un vecteur
int reverse_Cible_Y(vec C)
{
	return round(((1.9f-C(1))*((float)height)/2.0f)) + ymin;
}

// Gestion de la souris
void mouse(int button, int state, int x, int y)
{
	if (state == GLUT_DOWN) {
		computeCoordinates(get_New_Cible(x, y));
		glutPostRedisplay();
	} 
}

// Gestion des mouvements de la souris
void motion(int x, int y)
{
	computeCoordinates(get_New_Cible(x, y));
	glutPostRedisplay();
}

// Gestion du clavier classique
void parsekey(unsigned char key, int x, int y)
{
	switch (key) {

		case 'q':		//--------------------- pour Quitter
    		exit(0);

    	case 'i':		//--------------------- pour ré-Initialiser
    		myinit();
    		break;

    	default:
    		break;
   	}
	glutPostRedisplay();
}

// Gestion des touches spéciales du clavier
void parsekey_special(int key, int x, int y)
{
	vec delta = zeros(2);
	switch (key) {

  		case GLUT_KEY_UP:		//------------- Déplacement vers le haut
  			delta(1) = 0.01f;
			break;

    	case GLUT_KEY_DOWN:		//------------- Déplacement vers le bas
  			delta(1) = -0.01f;
			break;

		case GLUT_KEY_RIGHT:	//------------- Déplacement à droite
  			delta(0) = 0.01f;
			break;

		case GLUT_KEY_LEFT:	 	//------------- Déplacement à gauche
  			delta(0) = -0.01f;
      		break;

      	default:
      		break;
    }
	computeCoordinates(get_New_Cible(reverse_Cible_X(prev_Cible+delta),
									 reverse_Cible_Y(prev_Cible+delta)));    
	glutPostRedisplay();
}

// Gestion du redimensionnement de la fenêtre
void myReshape(int w, int h)
{
  xmin = 0;
  ymin = 0;
  if (w > h) {
    xmin = (w-h)/2;
    w = h;
  }
  else {
    ymin = (h-w)/2;
    h = w;
  }
  width = w;
  height = h;
  glViewport(xmin, ymin, w, h);
  glutPostRedisplay();
}

// Gestion de l'affichage
void display()
{
  	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
  	glPushMatrix();
		drawObject();
	glPopMatrix();
	glutSwapBuffers(); 
}

// Méthode d'initialisation
void myinit()
{
	prev_Cible = zeros(2);		// Ancienne cible à (0,0)
	Thetas = zeros(N);			// Les angles à 0
	Thetas(0) = PI/2;			// sauf le premier (->alignement vertical)
	glLoadIdentity();
	glClearColor(0.0f, 0.0f, 0.6f, 1.0f);

	// Initialisation des longueurs (ici toutes égales)
	Longueurs = zeros(N);
	for (int i=0; i<N; i++) {
		Longueurs(i) = 1.5f/N;
	}
	
	//Intialisation des Pivots
	Pivots = zeros(N+1, 2);
	computePivotCoordinates();
}


/** Le lanceur de l'application **************************************/

int main(int argc, char *argv[])
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_RGB | GLUT_DOUBLE | GLUT_MULTISAMPLE | GLUT_STENCIL);
	glutInitWindowPosition(200, 0);
	glutInitWindowSize(width, height);
	glutCreateWindow("Inverse Kinematics, N links");

	// Intialisation
	myinit();

	// Enregistrement des fonction de gestion GLUT
	glutDisplayFunc(display);  
	glutMouseFunc(mouse);
	glutMotionFunc(motion);
	glutKeyboardFunc(parsekey);
	glutSpecialFunc(parsekey_special);
	glutReshapeFunc(myReshape);

	// Lancement de la boucle
	glutSwapBuffers();
	glutMainLoop();
	return 0;
}