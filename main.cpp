
#include <iostream>
#include <gl\glut.h>
#include <math.h>

//----------------------------------------

float z=0;
typedef struct  //warna objek baru
{
    float r;
    float g;
    float b;
} color_t;

typedef struct  //warna objek baru
{
    color_t color;
} gradateColor_t;

typedef struct   //matrix 4 dimensi
{
    float m[4][4];
} matrix3D_t;

typedef struct   //vektor 4 dimensi
{
    float v[4];
} vector3D_t;

typedef struct
{
    float x;
    float y;
    float z;
} point3D_t;

typedef struct
{
    float x;
    float y;
} point2D_t;

typedef struct
{
    float v[3];
} vector2D_t;
//-----------------------------------
color_t bintik = {0.682, 0.678, 0.051};

//--------- matrices and vectors 3D 2nd version----------------//
matrix3D_t createIdentity(void)  //matrix dan vektor 3D ver 2
{
    matrix3D_t u;
    int i,j;
    for (i=0; i<4; i++)
    {
        for(j=0; j<4; j++) u.m[i][j]=0.;
        u.m[i][i]=1.;
    }
    return u;
}

vector3D_t operator * (matrix3D_t a, vector3D_t b)
{
    vector3D_t c;//c=a*b
    for (int i=0; i<4; i++)
    {
        c.v[i]=0;
        for (int j=0; j<4; j++)
            c.v[i]+=a.m[i][j]*b.v[j];
    }
    return c;
}

matrix3D_t operator * (matrix3D_t a, matrix3D_t b)  //matrix PERKALIAN
{
    matrix3D_t c;//c=a*b
    for (int i=0; i<4; i++)
    {
        for (int j=0; j<4; j++)
        {
            c.m[i][j]=0;
            for (int k=0; k<4; k++)
            {
                c.m[i][j]+=a.m[i][k]*b.m[k][j];
            }
        }
    }
    return c;
}

matrix3D_t translationMTX(float dx,float dy,float dz)  //translasi
{
    matrix3D_t trans=createIdentity();
    trans.m[0][3]=dx;
    trans.m[1][3]=dy;
    trans.m[2][3]=dz;
    return trans;
}

matrix3D_t rotationXMTX(float theta)  //rotasi X
{
    matrix3D_t rotate=createIdentity();
    float cs=cos(theta);
    float sn=sin(theta);
    rotate.m[1][1]=cs;
    rotate.m[1][2]=-sn;
    rotate.m[2][1]=sn;
    rotate.m[2][2]=cs;
    return rotate;
}

matrix3D_t rotationYMTX(float theta)  //rotasi Y
{
    matrix3D_t rotate=createIdentity();
    float cs=cos(theta);
    float sn=sin(theta);
    rotate.m[0][0]=cs;
    rotate.m[0][2]=sn;
    rotate.m[2][0]=-sn;
    rotate.m[2][2]=cs;
    return rotate;
}

matrix3D_t rotationZMTX(float theta)  //rotasi Z
{
    matrix3D_t rotate=createIdentity();
    float cs=cos(theta);
    float sn=sin(theta);
    rotate.m[0][0]=cs;
    rotate.m[0][1]=-sn;
    rotate.m[1][0]=sn;
    rotate.m[1][1]=cs;
    return rotate;
}

matrix3D_t scalingMTX(float factorx,float factory,float factorz)  //scaling
{
    matrix3D_t scale=createIdentity();
    scale.m[0][0]=factorx;
    scale.m[1][1]=factory;
    scale.m[2][2]=factorz;
    return scale;
}

matrix3D_t perspectiveMTX(float eyelength)
{
    matrix3D_t perspective=createIdentity();
    perspective.m[3][2]=-1./eyelength;
    return perspective;
}

point3D_t Vector2Point3D(vector3D_t vec)
{
    point3D_t pnt;
    pnt.x=vec.v[0];
    pnt.y=vec.v[1];
    pnt.z=vec.v[2];
    return pnt;
}

point2D_t Vector2Point2D(vector3D_t vec)
{
    point2D_t pnt;
    pnt.x=vec.v[0];
    pnt.y=vec.v[1];
    return pnt;
}

vector3D_t Point2Vector(point3D_t pnt)
{
    vector3D_t vec;
    vec.v[0]=pnt.x;
    vec.v[1]=pnt.y;
    vec.v[2]=pnt.z;
    vec.v[3]=1.;
    return vec;
}

vector3D_t homogenizeVector(vector3D_t vec)
{
    int i;
    for (i=0; i<3; i++)
    {
        vec.v[i]/=vec.v[3];
    }
    vec.v[3]=1.;
    return vec;
}

vector3D_t unitVector(vector3D_t vec)
{
    int i;
    float vec2=0.;
    float vec1,invvec1;
    for (i=0; i<3; i++)
    {
        vec2+=vec.v[i]*vec.v[i];
    }
    vec1=sqrt(vec2);
    if (vec1!=0.)
    {
        invvec1=1./vec1;
        for (i=0; i<3; i++)
        {
            vec.v[i]*=invvec1;
        }
    }
    vec.v[3]=1.;
    return vec;
}
// inner product (dot product) of homogeneous vector
float operator * (vector3D_t a, vector3D_t b)
{
    float c;//c=a*b
    int i;
    c=0;
    for (i=0; i<3; i++)
    {
        c+=a.v[i]*b.v[i];
    }
    return c;
}

// outer product (cross product ) of homogeneous vector
//       i         j         k
//       a0       a1        a2
//       b0       b1        b2
vector3D_t operator ^ (vector3D_t a, vector3D_t b)  //cross product
{
    vector3D_t c; //c=a*b
    c.v[0]=a.v[1]*b.v[2] - a.v[2]*b.v[1];
    c.v[1]=a.v[2]*b.v[0] - a.v[0]*b.v[2];
    c.v[2]=a.v[0]*b.v[1] - a.v[1]*b.v[0];
    c.v[3]=1;
    return c;
}

vector3D_t operator - (vector3D_t v)
{
    vector3D_t c;//c=-v
    c.v[0]=-v.v[0];
    c.v[1]=-v.v[1];
    c.v[2]=-v.v[2];
    c.v[3]=1.;
    return c;
}

vector3D_t operator * (float r, vector3D_t b)
{
    vector3D_t c;//c=r*b
    int i;
    for (i=0; i<3; i++)
    {
        c.v[i]=r*b.v[i];
    }
    c.v[3]=1.;
    return c;
}

vector3D_t operator * (vector3D_t b, float r)
{
    vector3D_t c;//c=r*b
    int i;
    for (i=0; i<3; i++)
    {
        c.v[i]=r*b.v[i];
    }
    c.v[3]=1.;
    return c;
}

float funcPositive(float x)
{
    if (0.<x) return x;
    else return 0.;
}

// x to yth power
float power(float x,float y)
{
    //ln z = y ln x        z = exp (y ln x)
    if (x==0.) return 0;
    return exp(y*log(x));
}

color_t operator + (color_t c1, color_t c2)
{
    color_t col;
    col.r=c1.r+c2.r;
    col.g=c1.g+c2.g;
    col.b=c1.b+c2.b;
    return col;
}

color_t operator * (float r, color_t c)
{
    color_t col;
    col.r=r*c.r;
    col.g=r*c.g;
    col.b=r*c.b;
    return col;
}

color_t operator * (color_t c, float r)
{
    color_t col;
    col.r=r*c.r;
    col.g=r*c.g;
    col.b=r*c.b;
    return col;
}

vector3D_t operator + (vector3D_t a, vector3D_t b)  //matrix PENJUMLAHAN
{
    vector3D_t c;
    for (int i=0; i<4; i++)
    {
        c.v[i]=a.v[i]+b.v[i];
    }
    return c;
}

vector3D_t operator - (vector3D_t a, vector3D_t b)  //matrix PENGURANGAN
{
    vector3D_t c;
    for (int i=0; i<4; i++)
    {
        c.v[i]=a.v[i]-b.v[i];
    }
    return c;
}

//PhongModel color calculation
// LightVector, NormalVector, ViewVector, ColorofObject
color_t PhongModel(vector3D_t Light,vector3D_t Normal,vector3D_t View,color_t col)
{
    float kspe=0.7; // specular reflection coefficient
    float kdif=0.9; // diffuse reflection coefficient
    float kamb=0.0002; // ambient light coefficient
    float tmp,NL,RV;
    color_t ColWhite= {1,1,1};
    vector3D_t ReflectionVector=(2.*(Light*Normal)*Normal)-Light;
    tmp=Normal*Light;
    NL=funcPositive(tmp);
    tmp=ReflectionVector*View;
    RV=funcPositive(tmp);
    return kdif*NL*col+kspe*power(RV,4)*ColWhite+kamb*col;
    //return kdif*NL*col+kamb*col;
}

//----------- End of matrices and vectors 3D 2nd version -----------



//----------- OpenGL drawShape function ver 1 -----------------------
void setColor(float red,float green,float blue)
{
    glColor3f(red, green, blue);
}

void setColor(color_t col)
{
    glColor3f(col.r, col.g, col.b);
}

void drawDot(float x,float y)
{
    glBegin(GL_POINTS);
    glVertex2f(x, y);
    glEnd();
}

void drawLine(float x1, float y1, float x2, float y2)
{
    glLineWidth(4);
    glBegin(GL_LINES);
    glVertex2f(x1, y1);
    glVertex2f(x2, y2);
    glEnd();
}

void drawLine(point2D_t p1,point2D_t p2)
{
    drawLine(p1.x,p1.y,p2.x,p2.y);
}

void drawPolyline(point2D_t pnt[],int n)
{
    int i;
    glBegin(GL_LINE_STRIP);
    for (i=0; i<n; i++)
    {
        glVertex2f(pnt[i].x, pnt[i].y);
    }
    glEnd();
}

void drawPolygon(point2D_t pnt[],int n)
{
    int i;
    glBegin(GL_LINE_LOOP);
    for (i=0; i<n; i++)
    {
        glVertex2f(pnt[i].x, pnt[i].y);
    }
    glEnd();
}

void fillPolygon(point2D_t pnt[],int n,color_t color)
{
    int i;
    setColor(color);
    glBegin(GL_POLYGON);
    for (i=0; i<n; i++)
    {
        glVertex2f(pnt[i].x, pnt[i].y);
    }
    glEnd();
}

void gradatePolygon(point2D_t pnt[],int num,color_t col[])
{
    int i;
    glBegin(GL_POLYGON);
    for (i=0; i<num; i++)
    {
        setColor(col[i]);
        glVertex2f(pnt[i].x, pnt[i].y);
    }
    glEnd();
}
//------------ End of OpenGL drawShape function ver 1 -----------------


//////////////////////////
void drawcharX(float x,float y)
{
    drawLine(x,y,x+10,y+12);
    drawLine(x,y+12,x+10,y);
}

void drawcharY(float x,float y)
{
    drawLine(x+5,y,x+5,y+7);
    drawLine(x,y+12,x+5,y+7);
    drawLine(x+10,y+12,x+5,y+7);
}

void drawcharZ(float x,float y)
{
    drawLine(x,y+12,x+10,y+12);
    drawLine(x+10,y+12,x,y);
    drawLine(x,y,x+10,y);
}
void drawAxes(matrix3D_t view)
{
#define HALFAXIS  220
#define HALFAXIS1 (HALFAXIS-10)
    point3D_t axes[14]=
    {
        {-HALFAXIS,0,0},{HALFAXIS,0,0},{HALFAXIS1,5,0},{HALFAXIS1,0,0},{0,0,0},
        {0,-HALFAXIS,0},{0,HALFAXIS,0},{0,HALFAXIS1,5},{0,HALFAXIS1,0},{0,0,0},
        {0,0,-HALFAXIS},{0,0,HALFAXIS},{5,0,HALFAXIS1},{0,0,HALFAXIS1}
    };
    vector3D_t vec[14];
    point2D_t buff[14];
    int i;
    for (i=0; i<14; i++)
    {
        vec[i]=Point2Vector(axes[i]);
        vec[i]=view*vec[i];
        buff[i]=Vector2Point2D(vec[i]);
    }
    drawPolyline(buff,14);
    drawcharX(buff[1].x,buff[1].y);
    drawcharY(buff[6].x,buff[6].y);
    drawcharZ(buff[11].x-14,buff[11].y);
}



typedef struct
{
    int NumberofVertices; //in the face
    short int pnt[50];
    color_t col;
    int NoObject;
    int z;
} face_t;

typedef struct
{
    int NumberofVertices; //of the object
    point3D_t pnt[1600];
    color_t col[1600];
    int NumberofFaces; //of the object
    face_t fc[1000];
} object3D_t;

void draw3D(object3D_t obyek,matrix3D_t mat)
{
    vector3D_t vec[1600], vecbuff[50];
    vector3D_t vecNormal;
    point2D_t p[50];
    int i,j;
    for(i=0; i<obyek.NumberofVertices; i++)
    {
        vec[i]=Point2Vector(obyek.pnt[i]);
        vec[i]=mat*vec[i];
    }
    //menggambar face invisible
    setColor(1,0,0);
    for(i=0; i<obyek.NumberofFaces; i++)
    {
        for(j=0; j<obyek.fc[i].NumberofVertices; j++)
            vecbuff[j]=vec[obyek.fc[i].pnt[j]];
        vecNormal=(vecbuff[1]-vecbuff[0])^(vecbuff[2]-vecbuff[0]);
        if(vecNormal.v[2]<0)
        {
            for(j=0; j<obyek.fc[i].NumberofVertices; j++)
            {
                p[j]=Vector2Point2D(vecbuff[j]);
            }
            drawPolygon(p,obyek.fc[i].NumberofVertices);
        }
    }
    //menggambar face visible
    setColor(0,1,1);
    for(i=0; i<obyek.NumberofFaces; i++)
    {
        for(j=0; j<obyek.fc[i].NumberofVertices; j++)
            vecbuff[j]=vec[obyek.fc[i].pnt[j]];
        vecNormal=(vecbuff[1]-vecbuff[0])^(vecbuff[2]-vecbuff[0]);
        if(vecNormal.v[2]>=0)
        {
            for(j=0; j<obyek.fc[i].NumberofVertices; j++)
            {
                p[j]=Vector2Point2D(vecbuff[j]);
            }
            drawPolygon(p,obyek.fc[i].NumberofVertices);
        }
    }
}
void draw3Da(object3D_t obyek,matrix3D_t mat)
{
    vector3D_t vec[1600], vecbuff[50];
    //vector3D_t vecNormal;
    point2D_t p[50];
    int i,j;
    for(i=0; i<obyek.NumberofVertices; i++)
    {
        vec[i]=Point2Vector(obyek.pnt[i]);
        vec[i]=mat*vec[i];
    }
    setColor(1,1,1);
    for(i=0; i<obyek.NumberofFaces; i++)
    {
        for(j=0; j<obyek.fc[i].NumberofVertices; j++)
            vecbuff[j]=vec[obyek.fc[i].pnt[j]];
        for(j=0; j<obyek.fc[i].NumberofVertices; j++)
            p[j]=Vector2Point2D(vecbuff[j]);
        drawPolygon(p,obyek.fc[i].NumberofVertices);
    }
}

void draw3Dcg(object3D_t obyek,matrix3D_t mat, color_t col)
{
    vector3D_t vec[1600], vecbuff[50];
    vector3D_t vecNormal;
    vector3D_t lightVector= {0,0,1,1},viewVector= {0,0,1,1};
    color_t colbuff[50];
    point2D_t p[50];
    int i,j;
    for(i=0; i<obyek.NumberofVertices; i++)
    {
        vec[i]=Point2Vector(obyek.pnt[i]);
        vec[i]=mat*vec[i];
    }
    for(i=0; i<obyek.NumberofFaces; i++)
    {
        for(j=0; j<obyek.fc[i].NumberofVertices; j++)
            vecbuff[j]=vec[obyek.fc[i].pnt[j]];
        vecNormal=(vecbuff[1]-vecbuff[0])^(vecbuff[2]-vecbuff[0]);
        if(vecNormal.v[2]<0)
        {
            for(j=0; j<obyek.fc[i].NumberofVertices; j++)
            {
                p[j]=Vector2Point2D(vecbuff[j]);
                vecNormal=unitVector(vecbuff[j]);
                colbuff[j]=PhongModel(lightVector,vecNormal,viewVector,obyek.col[obyek.fc[i].pnt[j]]);
            }
            gradatePolygon(p,obyek.fc[i].NumberofVertices,colbuff);
        }
    }
    for(i=0; i<obyek.NumberofFaces; i++)
    {
        for(j=0; j<obyek.fc[i].NumberofVertices; j++)
            vecbuff[j]=vec[obyek.fc[i].pnt[j]];
        vecNormal=(vecbuff[1]-vecbuff[0])^(vecbuff[2]-vecbuff[0]);
        if(vecNormal.v[2]>=0)
        {
            for(j=0; j<obyek.fc[i].NumberofVertices; j++)
            {
                p[j]=Vector2Point2D(vecbuff[j]);
                vecNormal=unitVector(vecbuff[j]);
                colbuff[j]=PhongModel(lightVector,vecNormal,viewVector,obyek.col[obyek.fc[i].pnt[j]]);
            }
            gradatePolygon(p,obyek.fc[i].NumberofVertices,colbuff);
        }
    }
}

void draw3Dw(object3D_t obyek,matrix3D_t mat,color_t col)
{
    vector3D_t vec[1600], vecbuff[50];
    vector3D_t vecNormal;
    vector3D_t lightVector= {0,0,1,1},viewVector= {0,0,1,1};
    color_t colbuff;
    point2D_t p[50];
    int i,j;
    for(i=0; i<obyek.NumberofVertices; i++)
    {
        vec[i]=Point2Vector(obyek.pnt[i]);
        vec[i]=mat*vec[i];
    }
    for(i=0; i<obyek.NumberofFaces; i++)
    {
        for(j=0; j<obyek.fc[i].NumberofVertices; j++)
            vecbuff[j]=vec[obyek.fc[i].pnt[j]];
        vecNormal=(vecbuff[1]-vecbuff[0])^(vecbuff[2]-vecbuff[0]);
        if(vecNormal.v[2]<0)
        {
            for(j=0; j<obyek.fc[i].NumberofVertices; j++)
            {
                p[j]=Vector2Point2D(vecbuff[j]);
            }
            vecNormal=unitVector(vecNormal);
            colbuff=PhongModel(lightVector,vecNormal,viewVector,col);
            fillPolygon(p,obyek.fc[i].NumberofVertices,colbuff);
        }
    }
    for(i=0; i<obyek.NumberofFaces; i++)
    {
        for(j=0; j<obyek.fc[i].NumberofVertices; j++)
            vecbuff[j]=vec[obyek.fc[i].pnt[j]];
        vecNormal=(vecbuff[1]-vecbuff[0])^(vecbuff[2]-vecbuff[0]);
        if(vecNormal.v[2]>=0)
        {
            for(j=0; j<obyek.fc[i].NumberofVertices; j++)
            {
                p[j]=Vector2Point2D(vecbuff[j]);
            }
            vecNormal=unitVector(vecNormal);
            colbuff=PhongModel(lightVector,vecNormal,viewVector,col);
            fillPolygon(p,obyek.fc[i].NumberofVertices,colbuff);
        }
    }
}

void draw3Dc(object3D_t obyek,matrix3D_t mat)
{
    vector3D_t vec[1600], vecbuff[50];
    vector3D_t vecNormal;
    vector3D_t lightVector= {0,0,1,1},viewVector= {0,0,1,1};
    color_t colbuff;
    point2D_t p[50];
    int i,j;
    for(i=0; i<obyek.NumberofVertices; i++)
    {
        vec[i]=Point2Vector(obyek.pnt[i]);
        vec[i]=mat*vec[i];
    }
    for(i=0; i<obyek.NumberofFaces; i++)
    {
        for(j=0; j<obyek.fc[i].NumberofVertices; j++)
            vecbuff[j]=vec[obyek.fc[i].pnt[j]];
        vecNormal=(vecbuff[1]-vecbuff[0])^(vecbuff[2]-vecbuff[0]);
        if(vecNormal.v[2]<0)
        {
            for(j=0; j<obyek.fc[i].NumberofVertices; j++)
            {
                p[j]=Vector2Point2D(vecbuff[j]);
            }
            vecNormal=unitVector(vecNormal);
            colbuff=PhongModel(lightVector,vecNormal,viewVector,obyek.fc[i].col);
            fillPolygon(p,obyek.fc[i].NumberofVertices,colbuff);
        }
    }
    for(i=0; i<obyek.NumberofFaces; i++)
    {
        for(j=0; j<obyek.fc[i].NumberofVertices; j++)
            vecbuff[j]=vec[obyek.fc[i].pnt[j]];
        vecNormal=(vecbuff[1]-vecbuff[0])^(vecbuff[2]-vecbuff[0]);
        if(vecNormal.v[2]>=0)
        {
            for(j=0; j<obyek.fc[i].NumberofVertices; j++)
            {
                p[j]=Vector2Point2D(vecbuff[j]);
            }
            vecNormal=unitVector(vecNormal);
            colbuff=PhongModel(lightVector,vecNormal,viewVector,obyek.fc[i].col);
            fillPolygon(p,obyek.fc[i].NumberofVertices,colbuff);
        }
    }
}
void createPlane(object3D_t &papan,int m, int n, float dx, float dy)
{
    int i,j,k;
    for(i=0; i<m; i++)
        for(j=0; j<n; j++)
        {
            k=i*n+j;
            papan.pnt[k].x=j*dx;
            papan.pnt[k].y=0;
            papan.pnt[k].z=i*dy;
        }
    papan.NumberofVertices=m*n;
    for(i=0; i<m-1; i++)
        for(j=0; j<n-1; j++)
        {
            k=i*(n-1)+j;
            papan.fc[k].NumberofVertices=4;
            papan.fc[k].pnt[0]=i*n+j;
            papan.fc[k].pnt[1]=(i+1)*n+j;
            papan.fc[k].pnt[2]=(i+1)*n+j+1;
            papan.fc[k].pnt[3]=i*n+j+1;
        }
    papan.NumberofFaces =(m-1)*(n-1);
}
void draw3Dcg(object3D_t obyek,matrix3D_t mat)
{
    vector3D_t vec[1600], vecbuff[50];
    vector3D_t vecNormal;
    vector3D_t lightVector= {0,0,1,1},viewVector= {0,0,1,1};
    color_t colbuff[50];
    point2D_t p[50];
    int i,j;
    for(i=0; i<obyek.NumberofVertices; i++)
    {
        vec[i]=Point2Vector(obyek.pnt[i]);
        vec[i]=mat*vec[i];
    }
    for(i=0; i<obyek.NumberofFaces; i++)
    {
        for(j=0; j<obyek.fc[i].NumberofVertices; j++)
            vecbuff[j]=vec[obyek.fc[i].pnt[j]];
        vecNormal=(vecbuff[1]-vecbuff[0])^(vecbuff[2]-vecbuff[0]);
        if(vecNormal.v[2]<0)
        {
            for(j=0; j<obyek.fc[i].NumberofVertices; j++)
            {
                p[j]=Vector2Point2D(vecbuff[j]);
                vecNormal=unitVector(vecbuff[j]);
                colbuff[j]=PhongModel(lightVector,vecNormal,viewVector,obyek.col[obyek.fc[i].pnt[j]]);
            }
            gradatePolygon(p,obyek.fc[i].NumberofVertices,colbuff);
        }
    }
    for(i=0; i<obyek.NumberofFaces; i++)
    {
        for(j=0; j<obyek.fc[i].NumberofVertices; j++)
            vecbuff[j]=vec[obyek.fc[i].pnt[j]];
        vecNormal=(vecbuff[1]-vecbuff[0])^(vecbuff[2]-vecbuff[0]);
        if(vecNormal.v[2]>=0)
        {
            for(j=0; j<obyek.fc[i].NumberofVertices; j++)
            {
                p[j]=Vector2Point2D(vecbuff[j]);
                vecNormal=unitVector(vecbuff[j]);
                colbuff[j]=PhongModel(lightVector,vecNormal,viewVector,obyek.col[obyek.fc[i].pnt[j]]);
            }
            gradatePolygon(p,obyek.fc[i].NumberofVertices,colbuff);
        }
    }
}

void makeCylinder(object3D_t &silinder, int n, float r, float h, color_t color)
{
    float a=6.28/n;
    int i;
    for(i=0; i<n; i++)
    {
        silinder.pnt[i].x=r*cos(i*a);
        silinder.pnt[i].y=0;
        silinder.pnt[i].z=r*sin(i*a);
        silinder.pnt[n+i].x=r*cos(i*a);
        silinder.pnt[n+i].y=h;
        silinder.pnt[n+i].z=r*sin(i*a);
    }
    silinder.NumberofVertices=2*n;
    for(i=0; i<n; i++)
    {
        silinder.fc[i].NumberofVertices=4;
        silinder.fc[i].pnt[0]=i;
        silinder.fc[i].pnt[1]=n+i;
        silinder.fc[i].pnt[2]=n+i+1;
        silinder.fc[i].pnt[3]=i+1;
        if(i==(n-1))
        {
            silinder.fc[i].pnt[2]=n;
            silinder.fc[i].pnt[3]=0;
        }
    }
    silinder.fc[n].NumberofVertices=n;
    for(i=0; i<n; i++) silinder.fc[n].pnt[i]=i;
    silinder.fc[n+1].NumberofVertices=n;
    for(i=0; i<n; i++) silinder.fc[n+1].pnt[i]=2*n-1-i;
    silinder.NumberofFaces=n+2;
    for(i=0; i<silinder.NumberofFaces; i++) silinder.fc[i].col=color;
    for(i=0; i<silinder.NumberofVertices; i++)
        silinder.col[i]=color;
}


void makeCone(object3D_t &kerucut, int n, float r,float h, color_t color)
{
    float a=6.28/n;
    int i;
    kerucut.pnt[0].x=0;
    kerucut.pnt[0].y=h;
    kerucut.pnt[0].z=0;
    for(i=1; i<=n; i++)
    {
        kerucut.pnt[i].x=r*cos(i*a);
        kerucut.pnt[i].y=0;
        kerucut.pnt[i].z=r*sin(i*a);
    }
    for(i=0; i<n; i++)
    {
        kerucut.fc[i].NumberofVertices=3;
        kerucut.fc[i].pnt[0]=0;
        kerucut.fc[i].pnt[1]=i+2;
        kerucut.fc[i].pnt[2]=i+1;
        if(i==(n-1)) kerucut.fc[i].pnt[1]=1;
    }
    kerucut.fc[n].NumberofVertices=n;
    for(i=0; i<n; i++) kerucut.fc[n].pnt[i]=i+1;
    kerucut.NumberofVertices=n+1;
    kerucut.NumberofFaces=n+1;
    for(i=0; i<kerucut.NumberofFaces; i++) kerucut.fc[i].col=color;
    for(i=0; i<kerucut.NumberofVertices; i++)
        kerucut.col[i]=color;
}

void makeCylinderN(object3D_t &silinder,int m,int n,float r[],float h[],int sw, color_t c)
{
    float a=6.26/n;
    float b=0;
    int i,j;
    silinder.NumberofVertices=(m+1)*n;
    for(i=0; i<=m; i++)
    {
        if(i>0) b=b+h[i-1];
        for(j=0; j<n; j++)
        {
            silinder.pnt[i*n+j].x=r[i]*cos(j*a);
            silinder.pnt[i*n+j].y=b;
            silinder.pnt[i*n+j].z=r[i]*sin(j*a);
        }
    }
    silinder.NumberofFaces=m*n+2;
    for(i=0; i<m; i++)
    {
        for(j=0; j<n; j++)
        {
            silinder.fc[i*n+j].NumberofVertices=4;
            silinder.fc[i*n+j].pnt[0]=i*n+j;
            silinder.fc[i*n+j].pnt[1]=(i+1)*n+j;
            silinder.fc[i*n+j].pnt[2]=(i+1)*n+j+1;
            silinder.fc[i*n+j].pnt[3]=i*n+j+1;
            if(j==(n-1))
            {
                silinder.fc[i*n+j].pnt[2]=i*n+j+1;
                silinder.fc[i*n+j].pnt[3]=(i-1)*n+j+1;
            }
        }
    }
    if(sw==0 || sw==1)
    {
        silinder.fc[m*n].NumberofVertices=n;
        for(i=0; i<n; i++) silinder.fc[m*n].pnt[i]=i;
    }
    if(sw==0 || sw==2)
    {
        silinder.fc[m*n+1].NumberofVertices=n;
        for(i=0; i<n; i++) silinder.fc[m*n+1].pnt[i]=(m+1)*n-1-i;
    }
    for(i=0; i<silinder.NumberofFaces; i++) silinder.fc[i].col=c;
    for(i=0; i<silinder.NumberofVertices; i++)
        silinder.col[i]=c;
}

void makeSphere(object3D_t &sphere,int n,float r, color_t color)
{
    float a=6.28/n;
    float b=6.28/n;
    int i,j;
    sphere.NumberofVertices=(n+1)*n;
    for(i=0; i<=n; i++)
    {
        for(j=0; j<n; j++)
        {
            sphere.pnt[i*n+j].x=r*cos(j*a)*sin(i*b);
            sphere.pnt[i*n+j].y=r*cos(i*b);
            sphere.pnt[i*n+j].z=r*sin(j*a)*sin(i*b);
        }
    }
    sphere.NumberofFaces=n*n+2;
    for(i=0; i<n; i++)
    {
        for(j=0; j<n; j++)
        {
            sphere.fc[i*n+j].NumberofVertices=4;
            sphere.fc[i*n+j].pnt[0]=i*n+j;
            sphere.fc[i*n+j].pnt[1]=(i+1)*n+j;
            sphere.fc[i*n+j].pnt[2]=(i+1)*n+j+1;
            sphere.fc[i*n+j].pnt[3]=i*n+j+1;
            if(j==(n-1))
            {
                sphere.fc[i*n+j].pnt[2]=i*n+j+1;
                sphere.fc[i*n+j].pnt[3]=(i-1)*n+j+1;
            }
        }
    }
    sphere.fc[n*n].NumberofVertices=n;
    for(i=0; i<n; i++) sphere.fc[n*n].pnt[i]=i;
    sphere.fc[n*n+1].NumberofVertices=n;
    for(i=0; i<n; i++) sphere.fc[n*n+1].pnt[i]=(n+1)*n-1-i;
    for(i=0; i<sphere.NumberofFaces; i++) sphere.fc[i].col=color;
    for(i=0; i<sphere.NumberofVertices; i++)
        sphere.col[i]=color;
}

point3D_t interpolate(point3D_t p1,point3D_t p2,float a)
{
    point3D_t p;
    p.x=(1-a)*p1.x+a*p2.x;
    p.y=(1-a)*p1.y+a*p2.y;
    p.z=(1-a)*p1.z+a*p2.z;
    return p;
}

void drawCircle(point2D_t p[], point2D_t p0, float r, int n, color_t color)
{
    //point2D_t p[360];
    float a=6.28/n;
    for(int i=0; i<n; i++)
    {
        p[i].x=p0.x+r*cos(i*a);
        p[i].y=p0.y+r*sin(i*a);
    }
    fillPolygon(p,n,color);

}

void createCircle(float r, int n, float x, float y, color_t color)
{
    point2D_t p[360];
    float a=6.28/n;
    for(int i=0; i<n; i++)
    {
        p[i].x=x+r*(float)cos((float)i*a);
        p[i].y=y+r*(float)sin((float)i*a);
    }
    fillPolygon(p,n,color);
}

void userdraw(void)
{
    //static float tick=0.,tick2=0, stick=0.;
    static float tick=1,stick=0.001,tick1=0.4,tick2=-0.5;

    //,dtick=1,timeTick=0;
    float theta=0.5;
    bool stop=false;
    color_t cokelat = {0.56, 0.24, 0.05};
    color_t hitam = {0,0,0};
    color_t putih = {1,1,1};
    color_t pink = {1, 0.714, 0.757};
    color_t biru = {0,0,1};
    color_t biruLaut = {0.149, 0.725, 0.784};
    color_t kuning = {1, 0.961, 0.424};
    color_t bajuSponge = {0.733, 0.545, 0.122};
    color_t merah= {1,0,0};
    object3D_t objek;
    matrix3D_t tilting=rotationXMTX(0.25)*rotationYMTX(-0.5);

    //lantai
    createPlane(objek,30,30,300,300);
    int k=0;
    for(int x=0; x<20; x++)
        for(int i=0; i<objek.NumberofVertices; i++) objek.col[i]=cokelat;
    matrix3D_t mat=tilting*translationMTX(-200,0,-100);
    draw3Dcg(objek,mat);

    //kota 1 dan 2
    float rbody[4]= {200,200};
    float hbody[3]= {300};
    object3D_t baju;
    object3D_t bajuPutih;
    float rbaju[4]= {200,200};
    float hbaju[3]= {50};
    float rbajuPutih[4]= {200,200};
    float hbajuPutih[3]= {30};

    object3D_t celanaKiri, celanaKanan, kakiKiri, KakiKanan,kakiKiriPutih, KakiKananPutih;
    float rCelana[4]= {30,30,28,28};
    float hCelana[3]= {30,0,20};
    float rKaki[4]= {10,10};
    float hKaki[3]= {50};

    //kaki
    makeCylinderN(kakiKiriPutih,3,30,rKaki,hKaki,0, putih);
    draw3Dc(kakiKiriPutih, tilting*translationMTX(200,-120,500));
    makeCylinderN(KakiKananPutih,3,30,rKaki,hKaki,0, putih);
    draw3Dc(KakiKananPutih, tilting*translationMTX(530,5,750));

    makeCylinderN(kakiKiri,3,30,rKaki,hKaki,0, kuning);
    draw3Dc(kakiKiri, tilting*translationMTX(200,-70,500));
    makeCylinderN(KakiKanan,3,30,rKaki,hKaki,0, kuning);
    draw3Dc(KakiKanan, tilting*translationMTX(530,40,750));

    //sepatu kiri
    int kurangiSepatuy = 300;
    int minX= 90;
    point2D_t sepatu[7] =
    {
        {15-minX,55-kurangiSepatuy}, {15-minX,25-kurangiSepatuy}, {60-minX,25-kurangiSepatuy}, {60-minX,40-kurangiSepatuy}, {60-minX,45-kurangiSepatuy}, {35-minX,45-kurangiSepatuy}, {35-minX,55-kurangiSepatuy}
    };
    fillPolygon(sepatu, 7,hitam);

    //sepatu kanan
    int tambahiSepatuy = 280;
    int minY= -80;
    point2D_t sepatu2[7] =
    {
        {15-minY,55-tambahiSepatuy}, {15-minY,25-tambahiSepatuy}, {60-minY,25-tambahiSepatuy}, {60-minY,40-tambahiSepatuy}, {60-minY,45-tambahiSepatuy}, {35-minY,45-tambahiSepatuy}, {35-minY
        ,55-tambahiSepatuy}
    };
    fillPolygon(sepatu2, 7,hitam);

    //baju
    makeCylinderN(celanaKiri,3,20,rCelana,hCelana,0, cokelat);
    draw3Dcg(celanaKiri, tilting*translationMTX(200,-20,500));
    makeCylinderN(celanaKanan,3,20,rCelana,hCelana,0, cokelat);
    draw3Dcg(celanaKanan, tilting*translationMTX(530,90,750));

    //badan
    makeCylinderN(baju,3,4,rbaju,hbaju,0, cokelat);
    draw3Dcg(baju, tilting*translationMTX(100,-20,180));
    makeCylinderN(bajuPutih,3,4,rbajuPutih,hbajuPutih,0, putih);
    draw3Dcg(bajuPutih, tilting*translationMTX(100,25,180));
    makeCylinderN(objek,3,4,rbody,hbody,0, kuning);
    draw3Dc(objek, tilting);
    //makeCylinderN(bajuPutih,3,4,rbajuPutih,hbajuPutih,0, putih);
    //draw3Dc(bajuPutih, tilting*translationMTX(100,25,180));

    //bintik-bintik
    createCircle(25,30,-60,210,bintik);
    createCircle(12,30,150,210,bintik);
    createCircle(15,30,-70,150,bintik);
    createCircle(20,30,-60,-20,bintik);
    createCircle(15,30,-70,20,bintik);
    createCircle(25,30,150,0,bintik);
    createCircle(15,30,150,50,bintik);

    //mata
    createCircle(50,10,0,150,putih);
    createCircle(50,30,90,158,putih);
    createCircle(20,10,0,150,biru);
    createCircle(20,30,90,158,biru);
    createCircle(10,10,0,150,hitam);
    createCircle(10,30,90,158,hitam);

    drawLine(75,205,68,225);
    drawLine(90,207,90,228);
    drawLine(105,205,113,228);

    drawLine(-18,195,-26,218);
    drawLine(-2,197,-2,222);
    drawLine(15,195,23,218);

    //senyum
    int s=18, g=10;
    point2D_t lidah[4] = {{90,90+g},{90,50},{90+g,50},{60+g,40+g}};
    //point2D_t lidah[4] = {{50,40+g},{50,30},{60+g,30},{60+g,40+g}};
    point2D_t senyum[5] = {{-35,60+s},{30,50+s},{50,55+s},{80,60+s},{120,80+s}};
    drawPolyline(senyum,5);
    fillPolygon(lidah,4,merah);

    /*object3D_t hidung;
    int a = 30;
    float rHidung[10]={30-a+5,40-a,50-a,60-a,70-a,70-a,60-a,50-a};
    float hHidung[9]={10,12,14,25,80,25,14};
    makeCylinderN(hidung,9,20,rHidung,hHidung,0,kuning);
    draw3Dc(hidung, tilting);*/


}

void display(void)
{
    glClear( GL_COLOR_BUFFER_BIT);
    userdraw();

    glutSwapBuffers();
}

int main (int argc, char ** argv)
{
    // insert code here...
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGB);
    glutInitWindowPosition(50,50);
    glutInitWindowSize(800,800);
    glutCreateWindow("Spongebob Squarepants");
    glClearColor(0.494, 0.631, 0.698,0);
    gluOrtho2D(-600,600,-600,600);
    glutIdleFunc(display);
    glutDisplayFunc(display);
    glutMainLoop();
    return 0;
}


