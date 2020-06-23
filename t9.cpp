//modified by: Jorge Zuniga
//date:
//
//3480 Computer Graphics
//Project
//Author: Gordon Griesel
//Date: 2017

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <X11/Xlib.h>
#include <X11/keysym.h>
#include <X11/Xutil.h>

//Variable types...
typedef double Flt;
typedef Flt Vec[3];
void vecMake(Flt,Flt,Flt,Vec);
#define rnd() ((Flt)rand() / (Flt)RAND_MAX)
#define vecComb(A,a,B,b,c) (c)[0] = (A) * (a)[0] + (B) * (b)[0]; \
                                    (c)[1] = (A) * (a)[1] + (B) * (b)[1]; \
(c)[2] = (A) * (a)[2] + (B) * (b)[2]
//constant definitions...
const int MAXOBJECTS =  300;
const int MAXLIGHTS =  100;
enum {
    TYPE_NONE=0,
    TYPE_DISK,
    TYPE_RING,
    TYPE_SPHERE,
    TYPE_CYLINDER,
    ORTHO,
    PERSPECTIVE,
    SURF_NONE,
    SURF_CHECKER
};

//Ray tracing structures
struct Ray {
    //Ray origin and direction
    Vec o, d;
};

struct Hit {
    //t     = distance to hit point
    //p     = hit point
    //norm  = normal of surface hit
    //color = color of surface hit
    Flt t;
    Vec p, norm, color;
};

//
//---------------------------------------------------
const int MAX_CLIPS = 100; 

struct Clip {
    //structure used for clipping functions of CSG.
    Vec center;
    Vec sCenter;
    Vec normal;
    Flt radius;
    bool inside;
    Clip() {
        radius = 0.0; 
        inside = false;

    }    
    void clear() { 
        radius = 0.0; 
        inside = false;
        vecMake(0.0, 0.0, 0.0, center);
        vecMake(0.0, 0.0, 0.0, sCenter);
        vecMake(0.0, 0.0, 0.0, normal);
        //Clip();
    }    
};
//--------------------------------------------------
//

struct Object {
    int type;
    Vec center;
    Vec norm;
    Flt radius, radius2;
    Vec color;
    int surface;
    bool inside;
    bool specular;
    Vec spec;
    Flt base, apex;
    Clip clip[MAX_CLIPS]; //---------------------
    int nclips; //---------------------
    Object() {
        specular = false;
        inside = false;
        vecMake(1.0, 1.0, 1.0, spec);
        nclips = 0; //--------------------
    }
    void clear_clips() {
        for (int i=0; i<MAX_CLIPS; i++)
            clip[i].clear();
    }
};

struct Light {
    Vec diffuse;
    Vec pos;
    bool shadow;
    Light(){
        shadow = true;
    }
};

class Global {
    public:
        int xres, yres;
        int mode;
        Object object[MAXOBJECTS];
        Light lights[MAXLIGHTS];
        Clip clip[MAX_CLIPS];
        //int nobjects;
        //int nclips;
        int nobjects, nlights, nclips;
        int color;
        int tilt;
        int checker;
        Vec ambient;
        Vec diffuse;
        Vec lightPos;
        Vec from, at, up;
        Flt angle;
        Global() {
            srand((unsigned)time(NULL));
            xres = 640, yres = 480;
            mode = 0;
            nobjects = 0;
            nlights = 0; //--------------------------
            nclips = 0; //--------------------------
            color = 0;
            tilt = 0;
            checker = 0;
            //---------------------------------------------------
            //Multiple light sources creation
            vecMake(.4,.4,.4,ambient);

            float d = 0.9;
            vecMake(d, d, d,lights[nlights].diffuse);
            vecMake(400.0,400.0,500.0,lights[nlights].pos);
            lights[nlights].shadow = false;
            nlights++;

            vecMake(d, d, d,lights[nlights].diffuse);
            vecMake(600.0,1200.0,700.0,lights[nlights].pos);
            lights[nlights].shadow = false;
            nlights++;

            vecMake(d, d, d,lights[nlights].diffuse);
            vecMake(40000.0,200000.0,0.0,lights[nlights].pos); //NEED TO EDIT
            lights[nlights].shadow = false;
            nlights++;

            vecMake(d, d, d,lights[nlights].diffuse);
            vecMake(1000.0,800.0,100.0,lights[nlights].pos);
            lights[nlights].shadow = false;
            nlights++;

            //-------------Light behind the paddle--------------
            vecMake(d, d, d,lights[nlights].diffuse);
            //vecMake(400.0,700.0,-1000.0,lights[nlights].pos);
            vecMake(600.0,700.0,-750.0,lights[nlights].pos);
            lights[nlights].shadow = true;
            nlights++;
            //
            //---------------------------------------------------

            vecMake(4.0, 250.0, 800.0, from);
            vecMake(0.0, 0.0, 0.0, at);
            vecMake(0.0, 1.0, 0.0, up);
            angle = 30.0;
        }
} g;

class X11_wrapper {
    private:
        Display *dpy;
        Window win;
        GC gc;
    public:
        X11_wrapper() {
            //constructor
            if (!(dpy=XOpenDisplay(NULL))) {
                fprintf(stderr, "ERROR: could not open display\n");
                fflush(stderr);
                exit(EXIT_FAILURE);
            }
            int scr = DefaultScreen(dpy);
            win = XCreateSimpleWindow(dpy, RootWindow(dpy, scr), 1, 1,
                    g.xres, g.yres, 0, 0x00ffffff, 0x00000000);
            XStoreName(dpy, win, "3480 ray tracer.  Press M for menu.");
            gc = XCreateGC(dpy, win, 0, NULL);
            XMapWindow(dpy, win);
            XSelectInput(dpy, win, ExposureMask | StructureNotifyMask |
                    PointerMotionMask | ButtonPressMask |
                    ButtonReleaseMask | KeyPressMask);
        }
        ~X11_wrapper() {
            XDestroyWindow(dpy, win);
            XCloseDisplay(dpy);
        }
        void getWindowAttributes(int *width, int *height) {
            XWindowAttributes gwa;
            XGetWindowAttributes(dpy, win, &gwa);
            *width = gwa.width;
            *height = gwa.height;
        }
        XImage *getImage(int width, int height) {
            //XImage *image = XGetImage(dpy, win,
            //            0, 0 , width, height, AllPlanes, ZPixmap);
            return XGetImage(dpy, win, 0, 0, width, height, AllPlanes, ZPixmap);
        }
        void checkResize(XEvent *e) {
            //ConfigureNotify is sent when window size changes.
            if (e->type != ConfigureNotify)
                return;
            XConfigureEvent xce = e->xconfigure;
            g.xres = xce.width;
            g.yres = xce.height;
        }
        void clearScreen() {
            //XClearWindow(dpy, win);
            setColor3i(0,0,0);
            XFillRectangle(dpy, win, gc, 0, 0, g.xres, g.yres);
        }
        bool getXPending() {
            return (XPending(dpy));
        }
        void getXNextEvent(XEvent *e) {
            XNextEvent(dpy, e);
        }
        void setColor3i(int r, int g, int b) {
            unsigned long cref = (r<<16) + (g<<8) + b;
            XSetForeground(dpy, gc, cref);
        }
        unsigned long rgbToLong(Vec rgb) {
            //Convert rgb[3] into an integer
            const float range = 255.999;
            int i;
            unsigned long cref = 0;
            for (i=0; i<3; i++) {
                //Don't let a color get too bright.
                if (rgb[i] > 1.0)
                    rgb[i] = 1.0;
                //Shift left 8 bits
                cref <<= 8;
                //Put next color component in low-order byte
                cref += (int)(rgb[i]*range);
            }
            return cref;
        }
        void drawPixel(int x, int y, Vec rgb) {
            unsigned long cref = rgbToLong(rgb);
            XSetForeground(dpy, gc, cref);
            XDrawPoint(dpy, win, gc, x, y);
        }
        void drawText(int x, int y, const char *text) {
            XDrawString(dpy, win, gc, x, y, text, strlen(text));
        }
} x11;

void init(void);
void checkResize(XEvent *e);
void checkMouse(XEvent *e);
int checkKeys(XEvent *e);
void render(int projection);
//vector functions...
void vecZero(Vec v);
void vecMake(Flt a, Flt b, Flt c, Vec v);
void vecCopy(Vec source, Vec dest);
void vecSub(Vec v0, Vec v1, Vec dest);
void vecNormalize(Vec v);
Flt vecDotProduct(Vec v0, Vec v1);
Flt vecLength(Vec v);


int main(void)
{
    srand((unsigned)time(NULL));
    init();
    x11.clearScreen();
    int done=0;
    while (!done) {
        while (x11.getXPending()) {
            XEvent e;
            x11.getXNextEvent(&e);
            x11.checkResize(&e);
            checkMouse(&e);
            done = checkKeys(&e);
            //render();
        }
    }
    return 0;
}

void takeScreenshot(const char *filename, int reset)
{
    //This function will capture your current X11 window,
    //and save it to a PPM P6 image file.
    //File names are generated sequentially.
    static int picnum = 0;
    int x,y;
    int width, height;
    x11.getWindowAttributes(&width, &height);
    if (reset)
        picnum = 0;
    XImage *image = x11.getImage(width, height);
    //
    //If filename argument is empty, generate a sequential filename...
    char ts[256] = "";
    strcpy(ts, filename);
    if (ts[0] == '\0') {
        sprintf(ts,"./lab5%02i.ppm", picnum);
        picnum++;
    }
    FILE *fpo = fopen(ts, "w");
    if (fpo) {
        fprintf(fpo, "P6\n%i %i\n255\n", width, height);
        for (y=0; y<height; y++) {
            for (x=0; x<width; x++) {
                unsigned long pixel = XGetPixel(image, x, y);
                fputc(((pixel & 0x00ff0000)>>16), fpo);
                fputc(((pixel & 0x0000ff00)>> 8), fpo);
                fputc(((pixel & 0x000000ff)    ), fpo);
            }
        }
        fclose(fpo);
    }
    XFree(image);
}

void init(void)
{
    //Setup some objects
    Object *o;
    g.nobjects = 0;
    Flt r1 = 6.0/255.0;
    Flt g1 = 109.0/255.0;
    Flt b1 = 226.0/255.0;
    //---------------------------------------------------------------
    //floor
    o = &g.object[g.nobjects];
    o->type = TYPE_DISK;
    vecMake(0.0, 0.0, 0.0, o->center);
    vecMake(0.0, 1.0, 0.0, o->norm);
    o->radius = 2000.0;
    o->specular = true;
    vecMake(0.1, 0.1, 0.1, o->spec); // original .2
    vecMake(0.2, 0.2, 0.2, o->spec); // original .6
    vecMake(1.0, 1.0, 1.0, o->color);
    o->surface = SURF_CHECKER;
    vecNormalize(o->norm);
    g.nobjects++;
    //---------------------------------------------------------------
    //Bottom of Paddle
    o = &g.object[g.nobjects];
    o->type = TYPE_DISK;
    vecMake(0.0, 1.0, 0.0, o->center);
    vecMake(0.0, 1.0, 0.0, o->norm);
    o->radius = 100.0;
    o->specular = false;
    vecMake(0.2, 0.2, 0.2, o->spec);
    vecMake(r1, g1, b1, o->color);
    o->surface = SURF_NONE;
    vecNormalize(o->norm);
    g.nobjects++;
    //---------------------------------------------------------------
    //RING
    o = &g.object[g.nobjects];
    o->type = TYPE_RING;
    vecMake(0.0, 40.0, 0.0, o->center);
    vecMake(0.0, 1.0, 0.0, o->norm);
    //o->specular = true;
    vecMake(r1, g1, b1, o->color);
    o->radius = 90;
    o->radius2 = 100;
    o->surface = SURF_NONE;
    g.nobjects++;
    //---------------------------------------------------------------
    //sphere thats going to be clipped for the inside of the paddle
    o = &g.object[g.nobjects];
    o->type = TYPE_SPHERE;
    vecMake(0.0, 40.0, 0.0, o->center);
    o->specular = false;
    vecMake(0.5, 0.5, 0.5, o->spec);
    vecMake(1.0/255.0, 77.0/255.0, 164.0/255.0, o->color);
    o->radius = 90.0;
    o->inside = true;
    o->surface = SURF_NONE;
    g.nobjects++;
    
    //Clip-------------------------------------------
    o->clear_clips();    
    vecMake(0.0, 40.0, 0.0, o->clip[o->nclips].center);
    vecMake(0.0, -1.0, 0.0, o->clip[o->nclips].normal);
    o->nclips++;
    
    //---------------------------------------------------------------
    //sphere sitting on cylinder
    o = &g.object[g.nobjects];
    o->type = TYPE_SPHERE;
    vecMake(0.0, 100.0, 0.0, o->center);
    o->specular = false;
    vecMake(0.5, 0.5, 0.5, o->spec);
    vecMake(r1, g1, b1, o->color);
    o->radius = 50.0;
    o->surface = SURF_NONE;
    g.nobjects++;
    //---------------------------------------------------------------
    //cylinder sitting under the sphere
    o = &g.object[g.nobjects];
    o->type = TYPE_CYLINDER;
    vecMake(0.0, 0.0, 0.0, o->center);
    vecMake(0.0, 0.0, -1.0, o->norm);
    vecMake(r1, g1, b1, o->color);
    o->apex = 100.0;
    o->radius = 50.0;
    o->surface = SURF_NONE;
    g.nobjects++;
    //---------------------------------------------------------------
    //cylinder sitting on floor
    o = &g.object[g.nobjects];
    o->type = TYPE_CYLINDER;
    //o->inside = true;
    vecMake(0.0, 0.0, 0.0, o->center);
    vecMake(0.0, 0.0, -1.0, o->norm);
    vecMake(r1, g1, b1, o->color);
    o->apex = 40.0;
    o->radius = 100.0;
    o->surface = SURF_NONE;
    g.nobjects++;
    //---------------------------------------------------------------
    //cylinder sitting on floor
    o = &g.object[g.nobjects];
    o->type = TYPE_CYLINDER;
    o->inside = true;
    vecMake(0.0, 0.0, 0.0, o->center);
    vecMake(0.0, 0.0, -1.0, o->norm);
    vecMake(r1, g1, b1, o->color);
    o->apex = 40.0;
    o->radius = 90.0;
    o->surface = SURF_NONE;
    g.nobjects++;
    //---------------------------------------------------------------
    //setup light and camera
    //vecMake(600.0, 400.0, 600.0, g.from);
    vecMake(600.0, 325.0, 600.0, g.from);
    vecMake(0.0, 80.0, 0.0, g.at);
    vecMake(0.0, 1.0, 0.0, g.up);
    g.angle = 30.0;
}

Flt vecDotProduct(Vec v0, Vec v1)
{
    return (v0[0]*v1[0] + v0[1]*v1[1] + v0[2]*v1[2]);
}

void vecCrossProduct(Vec v0, Vec v1, Vec dest)
{
    dest[0] = v0[1]*v1[2] - v1[1]*v0[2];
    dest[1] = v0[2]*v1[0] - v1[2]*v0[0];
    dest[2] = v0[0]*v1[1] - v1[0]*v0[1];
}

void vecZero(Vec v)
{
    v[0] = v[1] = v[2] = 0.0;
}

void vecNegate(Vec v)
{
    v[0] = -v[0];
    v[1] = -v[1];
    v[2] = -v[2];
}

void vecMake(Flt a, Flt b, Flt c, Vec v)
{
    v[0] = a;
    v[1] = b;
    v[2] = c;
}

void vecCopy(Vec source, Vec dest)
{
    dest[0] = source[0];
    dest[1] = source[1];
    dest[2] = source[2];
}

Flt vecLength(Vec vec)
{
    return sqrt(vecDotProduct(vec, vec));
}

void vecSub(Vec v0, Vec v1, Vec dest)
{
    dest[0] = v0[0] - v1[0];
    dest[1] = v0[1] - v1[1];
    dest[2] = v0[2] - v1[2];
}

void vecNormalize(Vec vec)
{
    //add code here to normalize a vector.
    Flt len = vecLength(vec);
    if (len == 0.0) {
        vecMake(1,0,0,vec);
        return;
    }
    vec[0] /= len;
    vec[1] /= len;
    vec[2] /= len;
}

void checkMouse(XEvent *e)
{
    if (e->type == ButtonPress) {
        //No mouse in this program.
    }
}

void showMenu()
{
    int y = 20;
    int inc = 16;
    x11.setColor3i(255, 255, 255);
    x11.drawText(10, y, "Menu");
    y += inc;
    x11.setColor3i(255, 255, 0);
    x11.drawText(10, y, "R - Render");
    y += inc;
}

int checkKeys(XEvent *e)
{
    if (e->type == KeyPress) {
        int key = XLookupKeysym(&e->xkey, 0);
        if (key == XK_m) {
            showMenu();
            return 0;
        }
        if (key == XK_a) {
            takeScreenshot("", 0);
            return 0;
        }
        //----------------------------------------
        if (key == XK_r) {
            init();
            render(PERSPECTIVE);
            return 0;
        }
        if (key == XK_p) {
            init();
            render(PERSPECTIVE);
            return 0;
        }
        if (key == XK_Escape) {
            if (g.mode) {
                g.mode = 0;
                x11.clearScreen();
                return 0;
            }
            return 1;
        }
    }
    return 0; 
}

Flt getLength(Vec p1, Vec p2)
{
    Flt d0 = p2[0] - p1[0];
    Flt d1 = p2[1] - p1[1];
    Flt d2 = p2[2] - p1[2];
    Flt len = sqrt(d0*d0 + d1*d1 + d2*d2);
    return len; 
}

Flt getArea(Vec v0, Vec v1, Vec v2)
{
    Flt a = getLength(v1, v0);
    Flt b = getLength(v2, v1);
    Flt c = getLength(v0, v2);
    Flt s = (a+b+c) / 2.0;
    return (sqrt(s * (s-a) * (s-b) * (s-c)));
}

int rayPlaneIntersect(Vec center, Vec normal, Ray *ray, Hit *hit)
{
    //http://mathworld.wolfram.com/Plane.html
    //
    //Where does the ray intersect the plane?
    //
    //plane equation:
    //
    // (p - center) . normal = 0.0
    //
    //where:
    //center = known point on plane
    //normal = normal vector of plane
    //
    //ray equation:
    //
    // O + tD
    //
    //where:
    //   ray origin = O
    //   ray direction = D
    //
    //Substitute ray equation into plane equation,
    //and solve for t
    //
    // (O + t * D - center) . normal = 0.0               (substitution)
    // 
    // (t * D + O - center) . normal = 0.0               (commutative)
    // 
    // t * (D . normal) + (O - center) . normal = 0.0    (distributive)
    // 
    // t * (D . normal) = - (O - center) . normal        (subtraction)
    // 
    // t = - ((O - center) . normal) / (D . normal)      (division)
    // 
    // notes: period is dot product
    //        solve for t
    //        populate hit->t and hit->p[] below...
    //        remember the unary minus sign above
    Vec v0;
    v0[0] = ray->o[0] - center[0];
    v0[1] = ray->o[1] - center[1];
    v0[2] = ray->o[2] - center[2];
    Flt dot1 = vecDotProduct(v0, normal);
    if (dot1 == 0.0)
        return 0;
    Flt dot2 = vecDotProduct(ray->d, normal);
    if (dot2 == 0.0)
        return 0;
    hit->t = -(dot1 / dot2);
    if (hit->t < 0.0)
        return 0.0;
    hit->p[0] = ray->o[0] + hit->t * ray->d[0];
    hit->p[1] = ray->o[1] + hit->t * ray->d[1];
    hit->p[2] = ray->o[2] + hit->t * ray->d[2];
    return 1;
}

int rayDiskIntersect(Object *o, Ray *ray, Hit *hit)
{
    //Does the ray intersect the disk's plane?
    if (rayPlaneIntersect(o->center, o->norm, ray, hit)) {
        //Yes.
        //Check that the hit point is within the disk radius
        //Use the point hit instead of the ray origin
        Flt d0, d1, d2, dist;
        d0 = o->center[0] - hit->p[0];
        d1 = o->center[1] - hit->p[1];
        d2 = o->center[2] - hit->p[2];
        //d2 = 0.0;
        dist = sqrt(d0*d0 + d1*d1 + d2*d2);
        if (dist <= o->radius) {
            for(int i = 0; i < o->nclips; i++) {
                Vec w;
                vecSub(hit->p, o->clip[i].center, w);
                Flt dot; 
                dot = vecDotProduct(o->clip[i].normal, w);
                // printf("dot=%lf\n",dot);
                if(dot>0.0)
                    return 0;
                Vec v;
                vecSub(hit->p, o->clip[i].center, v);
                Flt len = vecLength(v);
                if(len < o->clip[i].radius)
                    return 0;
            }
            //Hit is within radius
            return 1;
        }
    }    
    return 0;
}


void sphereNormal(Vec center,Vec hitPoint, Vec norm)
{
    //Calc normal at hit
    norm[0] = hitPoint[0] - center[0];
    norm[1] = hitPoint[1] - center[1];
    norm[2] = hitPoint[2] - center[2];
    vecNormalize(norm);
}

int raySphereIntersect(Object *o, Ray *ray, Hit *hit)
{
    //printf("raySphereIntersect()...\n");
    //Determine if and where a ray intersects a sphere.
    //
    // sphere equation:
    // (p - c) * (p - c) = r * r
    //
    // where:
    // p = point on sphere surface
    // c = center of sphere
    //
    // ray equation:
    // o + t*d
    //
    // where:
    //   o = ray origin
    //   d = ray direction
    //   t = distance along ray, or scalar
    //
    // substitute ray equation into sphere equation
    //
    // (o + t*d - c) * (o + t*d - c) - r * r = 0
    //
    // we want it in this form:
    // a*t*t + b*t + c = 0
    //
    // (o + d*t - c)
    // (o + d*t - c)
    // -------------
    // o*o + o*d*t - o*c + o*d*t + d*t*d*t - d*t*c - o*c + c*d*t + c*c
    // d*t*d*t + o*o + o*d*t - o*c + o*d*t - d*t*c - o*c + c*d*t + c*c
    // d*t*d*t + 2(o*d*t) - 2(c*d*t) + o*o - o*c - o*c + c*c
    // d*t*d*t + 2(o-c)*d*t + o*o - o*c - o*c + c*c
    // d*t*d*t + 2(o-c)*d*t + (o-c)*(o-c)
    //
    // t*t*d*d + t*2*(o-c)*d + (o-c)*(o-c) - r*r
    //
    // a = dx*dx + dy*dy + dz*dz
    // b = 2(ox-cx)*dx + 2(oy-cy)*dy + 2(oz-cz)*dz
    // c = (ox-cx)*(ox-cx) + (oy-cy)*(oy-cy) + (oz-cz)*(oz-cz) - r*r
    //
    // now put it in quadratic form:
    // t = (-b +/- sqrt(b*b - 4ac)) / 2a
    //
    //
    //1. a, b, and c are given to you just above.
    //2. Create variables named a,b,c, and assign the values you see above.
    //3. Look how a,b,c are used in the quadratic equation.
    //4. Make your code solve for t.
    //5. Remember, a quadratic can have 0, 1, or 2 solutions.
    //
    Flt a = ray->d[0]*ray->d[0] + ray->d[1]*ray->d[1] + ray->d[2]*ray->d[2];
    Flt b = 2.0*(ray->o[0]-o->center[0])*ray->d[0] +
        2.0*(ray->o[1]-o->center[1])*ray->d[1] +
        2.0*(ray->o[2]-o->center[2])*ray->d[2];
    Flt c = (ray->o[0]-o->center[0])*(ray->o[0]-o->center[0]) +
        (ray->o[1]-o->center[1])*(ray->o[1]-o->center[1]) +
        (ray->o[2]-o->center[2])*(ray->o[2]-o->center[2]) -
        o->radius*o->radius;
    Flt t0,t1;
    //discriminant
    Flt disc = b * b - 4.0 * a * c;
    if (disc < 0.0) return 0;
    disc = sqrt(disc);
    t0 = (-b - disc) / (2.0*a);
    t1 = (-b + disc) / (2.0*a);
    if (t0 > 0.0) {
        hit->p[0] = ray->o[0] + ray->d[0] * t0;
        hit->p[1] = ray->o[1] + ray->d[1] * t0;
        hit->p[2] = ray->o[2] + ray->d[2] * t0;
        sphereNormal(hit->p, o->center, hit->norm);
        hit->t = t0;
        //----------------------------------
        if(o->nclips){
            for(int i = 0; i < o->nclips; i++){
                Vec v;
                Vec w;
                vecSub(hit->p, o->clip[i].sCenter, v);
                vecSub(hit->p, o->clip[i].center, w);
                Flt len = vecLength(v);
                Flt dot = vecDotProduct(o->clip[i].normal,w);
                if(len < o->clip[i].radius && !o->clip[i].inside)
                    goto sp_hit2;
                if(len > o->clip[i].radius && o->clip[i].inside)
                    goto sp_hit2;
                if(dot < 0.0)
                    goto sp_hit2;
            }
        }
        return 1;
    }
sp_hit2:
    if (t1 > 0.0) {
        hit->p[0] = ray->o[0] + ray->d[0] * t1;
        hit->p[1] = ray->o[1] + ray->d[1] * t1;
        hit->p[2] = ray->o[2] + ray->d[2] * t1;
        sphereNormal(hit->p, o->center, hit->norm);
        hit->t = t1;
        if(o->nclips){
            for(int i = 0; i < o->nclips; i++){
                Vec v;
                Vec w;
                vecSub(hit->p, o->clip[i].sCenter, v);
                vecSub(hit->p, o->clip[i].center, w);
                Flt len = vecLength(v);
                Flt dot = vecDotProduct(o->clip[i].normal,w);
                if(len < o->clip[i].radius && !o->clip[i].inside)
                    return 0;
                if(len > o->clip[i].radius && o->clip[i].inside)
                    return 0;
                if(dot < 0.0)
                    return 0;
            }
        }
        return 1;
    }
    return 0;
}
void reflect(Vec I, Vec N, Vec R)
{
    //I = incident vector
    //N = the surface normal
    //R = reflected ray
    Flt dot = -vecDotProduct(I, N);
    Flt len = 2.0 * dot;
    R[0] = len * N[0] + I[0];
    R[1] = len * N[1] + I[1];
    R[2] = len * N[2] + I[2];
}

int rayRingIntersect(Object *o, Ray *ray, Hit *hit)
{
    if (rayPlaneIntersect(o->center, o->norm, ray, hit)) {
        Flt d0 = o->center[0] - hit->p[0];
        Flt d1 = o->center[1] - hit->p[1];
        Flt d2 = o->center[2] - hit->p[2];
        Flt dist = sqrt(d0*d0 + d1*d1 + d2*d2);
        if (dist >= o->radius && dist <= o->radius2) {
            vecCopy(o->norm, hit->norm);
            return 1;
        }
    }
    return 0;
}

void cylinderNormal(Vec p, Vec norm)
{
    //Center of cylinder is at the origin
    vecMake(p[0], 0.0, p[2], norm);
    vecNormalize(norm);
}

int rayCylinderIntersect(Object *o, Ray *ray, Hit *hit)
{
    //http://voices.yahoo.com/
    //developing-equation-cone-simplest-case-2522846.html?cat=17
    // x2 + z2 = a*y2
    //----------------------------------
    //for cylinder centered at origin...
    //----------------------------------
    //
    // cylinder equation
    // x2 + z2 = r2
    //
    // where:
    //   x = x component of point on cylinder surface
    //   z = z component of point on cylinder surface
    //   r = radius
    //
    // ray equation:
    // o + t*d
    //
    // where:
    //   o = ray origin
    //   d = ray direction
    //   t = distance along ray, or scalar
    //
    // substitute ray equation into sphere equation
    //
    // (ox+t*dx)2 + (oz+t*dz)2 - r*r = 0
    //
    // where:
    //   ox = x component of ray origin
    //   oz = z component of ray origin
    //   dx = x component of ray direction
    //   dz = z component of ray direction
    //
    // ox + t*dx
    // ox + t*dx
    // --------------------------------
    // ox*ox + 2(ox * t*dx) + t*t*dx*dx
    //
    // add in the z components...
    //
    //ox*ox + 2(ox * t*dx) + t*t*dx*dx + oz*oz + 2(oz * t*dz) + t*t*dz*dz - r*r
    //
    //Goal is to solve for t using the quadratic equation...
    // t = (-b +/- sqrt(b*b - 4ac)) / 2a
    //
    //t*t*dx*dx + t*t*dz*dz + ox*ox + 2(ox * t*dx) + oz*oz + 2(oz * t*dz) - r*r
    //t*t*dx*dx + t*t*dz*dz + 2(t*ox*dx) + 2(t*oz*dz) + ox*ox + oz*oz - r*r
    //t*t*dx*dx + t*t*dz*dz + 2t*ox*dx + 2t*oz*dz + ox*ox + oz*oz - r*r
    //t*t*dx*dx + t*t*dz*dz + t*2*ox*dx + t*2*oz*dz + ox*ox + oz*oz - r*r
    //a = dx*dx + dz*dz
    //b = 2*ox*dx + 2*oz*dz
    //c = ox*ox + oz*oz - r*r
    //
    Ray r;
    vecCopy(ray->o, r.o);
    vecCopy(ray->d, r.d);
    //now put a,b,c into C source code...
    Flt a = r.d[0] * r.d[0] + r.d[2] * r.d[2];
    Flt b = 2.0 * r.o[0] * r.d[0] + 2.0 * r.o[2] * r.d[2];
    Flt c = r.o[0]*r.o[0] + r.o[2]*r.o[2] - o->radius * o->radius;
    //
    Flt t0,t1;
    //disc:  discriminant
    Flt disc = b * b - 4.0 * a * c;
    if (disc < 0.0) return 0;
    disc = sqrt(disc);
    t0 = (-b - disc) / (2.0*a);
    t1 = (-b + disc) / (2.0*a);
    if (t0 > 0.0) {
        //possible hit
        hit->p[0] = r.o[0] + r.d[0] * t0;
        hit->p[1] = r.o[1] + r.d[1] * t0;
        hit->p[2] = r.o[2] + r.d[2] * t0;
        if (hit->p[1] >= 0.0 && hit->p[1] <= o->apex) {
            hit->t = t0;
            cylinderNormal(hit->p, hit->norm);
            return 1;
        }
    }
    if (t1 > 0.0) {
        hit->p[0] = r.o[0] + r.d[0] * t1;
        hit->p[1] = r.o[1] + r.d[1] * t1;
        hit->p[2] = r.o[2] + r.d[2] * t1;
        if (hit->p[1] >= 0.0 && hit->p[1] <= o->apex) {
            hit->t = t1;
            cylinderNormal(hit->p, hit->norm);
            return 1;
        }
    }
    return 0;
}

Flt getShadow(Hit *z, Vec rLight){
    Light *l;
    Flt specHigh = 0.0;
    vecMake(0.0, 0.0, 0.0, rLight);
    bool inShadow;
    Ray ray;
    Flt lightDist;
    for(int h =0; h < g.nlights; h++) {
        l = &g.lights[h];
        Flt dot;
        inShadow = false;
        vecSub(l->pos, z->p, ray.d);
        vecCopy(z->p, ray.o);
        lightDist = vecLength(ray.d);
        vecNormalize(ray.d);
        vecNormalize(z->norm);
        dot = vecDotProduct(z->norm, ray.d);
        //nudge
        ray.o[0] += ray.d[0] + 0.001;   
        ray.o[1] += ray.d[1] + 0.001;   
        ray.o[2] += ray.d[2] + 0.001;   
        Hit hit;
        Object *o;
        for (int i=0; i<g.nobjects; i++) {
            o = &g.object[i];
            switch (o->type) {
                case TYPE_DISK:
                    if (rayDiskIntersect(o, &ray, &hit)){
                        // return (hit.t < lightDist);
                        if(hit.t < lightDist)
                        {
                            if(l->shadow)   
                                inShadow = true;
                        }
                    }
                    break;
                case TYPE_RING:
                    if (rayRingIntersect(o, &ray, &hit)) {
                        // return (hit.t < lightDist);
                        if(hit.t < lightDist)
                        {           

                            if(l->shadow)   
                                inShadow = true;
                        }
                    }
                    break;
                case TYPE_SPHERE:
                    if (raySphereIntersect(o, &ray, &hit)) {
                        // return (hit.t < lightDist);
                        if(hit.t < lightDist)
                        {
                            if(l->shadow)   
                                inShadow = true;
                        }
                    }
                    break;
                case TYPE_CYLINDER:
                    if (rayCylinderIntersect(o, &ray, &hit)) {
                        // return (hit.t < lightDist);
                        if(hit.t < lightDist)
                        {
                            if(l->shadow)   
                                inShadow = true;
                        }
                    }
                    break;
            }
        }
        if (dot < 0) {
            dot =0;
        }
        //specHigh += pow(dot,40)/2;
        specHigh += pow(dot,40)/3.8;
        if (!inShadow ) {
            rLight[0] += (dot * l->diffuse[0])/g.nlights;
            rLight[1] += (dot * l->diffuse[1])/g.nlights;
            rLight[2] += (dot * l->diffuse[2])/g.nlights;
        }

    }
    return specHigh;
}

void trace(Ray *ray, Vec rgb, Flt weight, int level)
{
    if (level > 8) return;
    if (weight < 0.01) return;
    //Trace a ray through the scene.
    int i;
    Hit hit, closehit;
    Object *o;
    int h = -1;
    closehit.t = 9e9;
    for (i=0; i<g.nobjects; i++) {
        o = &g.object[i];
        switch (o->type) {
            case TYPE_DISK:
                if (rayDiskIntersect(o, ray, &hit)) {
                    if (hit.t < closehit.t) {
                        closehit.t = hit.t;
                        vecCopy(hit.p, closehit.p);
                        vecCopy(o->color, closehit.color);
                        vecCopy(o->norm, closehit.norm);
                        h=i;
                    }
                }
                break;
            case TYPE_RING:
                if (rayRingIntersect(o, ray, &hit)) {
                    if (hit.t < closehit.t) {
                        closehit.t = hit.t;
                        vecCopy(hit.p, closehit.p);
                        vecCopy(o->norm, closehit.norm);
                        vecCopy(o->color, closehit.color);
                        h=i;
                    }
                }
                break;
            case TYPE_SPHERE:
                if (raySphereIntersect(o, ray, &hit)) {
                    if (hit.t < closehit.t) {
                        closehit.t = hit.t;
                        vecCopy(hit.p, closehit.p);
                        vecCopy(o->color, closehit.color);
                        sphereNormal(o->center, closehit.p, closehit.norm);
                        if (o->inside)
                            vecNegate(closehit.norm);
                        h=i;
                    }
                }
                break;
            case TYPE_CYLINDER:
                if (rayCylinderIntersect(o, ray, &hit)) {
                    if (hit.t < closehit.t) {
                        closehit.t = hit.t;
                        vecCopy(hit.p, closehit.p);
                        vecCopy(o->color, closehit.color);
                        vecCopy(hit.norm, closehit.norm);
                        if (o->inside)
                            vecNegate(closehit.norm);
                        h=i;
                    }
                }
                break;
        }
    }
    if (h < 0) {
        //ray did not hit an object.
        //if the ray did not hit anything then default to gray instead of black
        rgb[0] = .6;
        rgb[1] = .6;
        rgb[2] = .6;
        return;
    }
    if(closehit.t > 2000) { //2000 dictates the view distance
        rgb[0] = .6; //color of fog
        rgb[1] = .6; //color of fog
        rgb[2] = .6; //color of fog
        return;
    }
    //The ray hit an object.
    //Set the color of the pixel to the color of the object.
    o = &g.object[h];
    if (o->surface == SURF_CHECKER) {
        Flt x = 2000.0 + (closehit.p[0]) / 40.0; //2000 unit shift, the division by 40 determines the size of my squares
        Flt y = 2000.0 + (closehit.p[2]) / 40.0; //2000 unit shift, the division by 40 determines the size of my squares
        Flt cX = x;
        Flt cY = y;
        Flt dX = cX - (int)cX; // gives me decimal values to help determine the white border
        Flt dY = cY - (int)cY; // gives me decimal values to help determine the white border
        if(dX<.05 || dY<.05) {
            closehit.color[0] = 1.0;
            closehit.color[1] = 1.0;
            closehit.color[2] = 1.0;
        } else {
            int ciX = (int) cX;
            int ciY = (int) cY;
            //-------EVEN OR ODD-_------
            x = ciX % 2;
            y = ciY % 2;
            //--------------------------
            if (x==y) {
                closehit.color[0] = 0.62;
                closehit.color[1] = 0.62;
                closehit.color[2] = 0.62;
            } else {
                closehit.color[0] = 0.62;
                closehit.color[1] = 0.62;
                closehit.color[2] = 0.62;
            }
        }
    }
    if (o->specular == true) {
        Vec trgb = {0.0};
        Ray tray;
        vecCopy(closehit.p, tray.o);
        reflect(ray->d, closehit.norm, tray.d); 
        tray.o[0] += tray.d[0] + 0.001; 
        tray.o[1] += tray.d[1] + 0.001; 
        tray.o[2] += tray.d[2] + 0.001; 
        trace(&tray, trgb, weight*0.5, level+1);
        rgb[0] += trgb[0] * o->spec[0] * weight;
        rgb[1] += trgb[1] * o->spec[1] * weight;
        rgb[2] += trgb[2] * o->spec[2] * weight;
    }

    //ambient light 
    rgb[0] += closehit.color[0] *g.ambient[0];
    rgb[1] += closehit.color[1] *g.ambient[1];
    rgb[2] += closehit.color[2] *g.ambient[2];

    Flt specHigh;
    Vec rLight;
    specHigh = getShadow(&closehit, rLight);
    rgb[0] += closehit.color[0] * rLight[0] + specHigh;
    rgb[1] += closehit.color[1] * rLight[1] + specHigh;
    rgb[2] += closehit.color[2] * rLight[2] + specHigh;
    //-----------------
    //
    Flt fog = closehit.t/2000;//this controls opacity, the closer it is to 1, the more it looks like fog
    fog = pow(fog,1.25); //this power function controls the instensity of the fog, meaning how dense/close it is
    rgb[0] = rgb[0]*(1.0-fog)+.6*fog;//link on website helps explain how fog
    rgb[1] = rgb[1]*(1.0-fog)+.6*fog;//interpolation works
    rgb[2] = rgb[2]*(1.0-fog)+.6*fog;
    return;
}

void setupRay(Vec eye, Vec pixel, Ray *ray)
{
    //Make a ray from eye through a screen pixel
    vecCopy(eye, ray->o);
    vecSub(pixel, eye, ray->d);
    vecNormalize(ray->d);
}

/*void render(int projection)
  {
  void castRaysFromCamera();
  castRaysFromCamera();
  return;
  }

*/

void vecAdd(Vec v0, Vec v1, Vec dest) {
    dest[0] = v0[0] + v1[0];
    dest[1] = v0[1] + v1[1];
    dest[2] = v0[2] + v1[2];
}

Flt degreesToRadians(Flt angle) {
    return (angle / 360.0) * (3.14159265 * 2.0);
}

void render(int projection)
{
    printf("castRaysFromCamera()...\n");
    //
    Vec from, at, up;
    Flt angle;
    //
    vecCopy(g.from, from);
    vecCopy(g.at, at);
    vecCopy(g.up, up);
    angle = g.angle;
    //
    Flt fyres = (Flt)g.yres;
    Flt fxres = (Flt)g.xres;
    Flt ty = 1.0 / (fyres - 1.0);
    Flt tx = 1.0 / (fxres - 1.0);
    int px = 1;
    int py = 1;
    vecNormalize(up);
    Flt viewAnglex, aspectRatio;
    Flt frustumheight, frustumwidth;
    Vec rgb, eye, dir, left, out;
    vecSub(at, from, out);
    vecNormalize(out);
    aspectRatio = fxres / fyres;
    //-------------------------------------------------------------------------
    viewAnglex = degreesToRadians(angle * 0.5);
    frustumwidth = tan(viewAnglex);
    frustumheight = frustumwidth / aspectRatio;
    //frustumwidth is actually half the distance across screen
    //compute the left and up vectors...
    vecCrossProduct(out, up, left);
    vecNormalize(left);
    vecCrossProduct(left, out, up);
    //
    Ray ray;
    vecCopy(from, eye);
    //Trace every pixel...
    int istart = 0;
    int iend = g.yres;
    int jstart = 0;
    int jend = g.xres;
    //
    for (int i=istart; i<iend; i++) {
        py = i;
        for (int j=jstart; j<jend; j++) {
            px = j;
            //Start the color at black
            //Start the ray origin at the eye
            vecZero(rgb);
            vecCopy(eye, ray.o);
            vecComb(-frustumheight * (2.0 * (Flt)py*ty - 1.0), up,
                    frustumwidth  * (2.0 * (Flt)px*tx - 1.0), left,
                    dir);
            int sub = 0;
            if(sub == 0) {
                Vec tcol = {0.0};
                for(int k = 0; k < 16; k++){
                    vecComb(-frustumheight * (2.0 *(((Flt)py+rnd()-0.5)*ty) - 1.0), up,
                            frustumwidth  * (2.0 *(((Flt)px+rnd()-0.5)*tx) - 1.0), left,
                            dir);
                    vecAdd(dir,out,ray.d);
                    vecNormalize(ray.d);
                    vecZero(rgb);
                    trace(&ray,rgb,1.0,1);
                    vecAdd(rgb,tcol,tcol);
                }
                rgb[0] = tcol[0] / 16;
                rgb[1] = tcol[1] / 16;
                rgb[2] = tcol[2] / 16;
                x11.drawPixel(j,i,rgb);
            }
            else {
                vecAdd(dir, out, ray.d);
                vecNormalize(ray.d);
                trace(&ray, rgb, 1.0, 1);
                x11.drawPixel(j, i, rgb);   
                /*vecAdd(dir, out, ray.d);
                  vecNormalize(ray.d);
                  trace(&ray, rgb, 1.0, 1);
                  x11.drawPixel(j, i, rgb);*/
            }
        }
    }
}

