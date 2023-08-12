
/*
 Алгоритм
 Данн прямоугольник с x1 y1 x2 y2 c вырезом x3 y3 x4 y4
 Требуется разбить его на прямоугольные треугольники и затем интерполировать исходную функцию конечными элементами степени три и нарисовать все в гнуплоте
 1.Разбиение прямоугольника с вырезом на маленькие прямоугольники путем деления прямоугольника на 8 маленьких прямоугольников
 2.Разбиваем прямоугольники на треугольники
 3.Высчитываем Pf для каждого треугольника
 4.Создаем мелкую сектку внутри прямоугольника
 5.Считаем интеполиционные значения и реальные значения по сетке
 6.Выводим все два файла
 7.Рисуем
 */

#include <vector>
#include <iostream>
#include <math.h>
#include <cmath>
#include <fstream>

using namespace std;

double fn(double x, double  y){
    return x*x+y*y;
}

class Point{
public:
    double x, y;
    double z;
    Point(){};
    Point(double a, double b){
        x=a; y=b; z=0;
    }
    Point(double a, double b, double c){
        x=a; y=b; z=c;
    };
    void computeZ(){
        z=fn(x, y);
    }
};

double fn(Point A){
    return fn(A.x, A.y);
}

class Triangle{
public:
    Point v1, v2, v3;
    vector<Point> P;
    Triangle(){};
    Triangle(Point a, Point b, Point c){
        v1=a; v2=b; v3=c;
    }
    bool isPointInTriangle();
    
    Point computeFormula(int k, int l, int s, vector<double> x, vector<double> y){
        Point A;
        A.x=((3-s)*x[k]+s*x[l])/3;
        A.y=((3-s)*y[k]+s*y[l])/3;
        return A;
    }
    
    void computeP1_10() {
        Point zero(0, 0, 0);
        P.push_back(zero);
        
        P.push_back(v1); P.push_back(v2); P.push_back(v3);
        
        vector<double> x;
        vector<double> y;
        
        x.push_back(0); x.push_back(v1.x); x.push_back(v2.x); x.push_back(v3.x);
        y.push_back(0); y.push_back(v1.y); y.push_back(v2.y); y.push_back(v3.y);
        
        P.push_back(computeFormula(1, 2, 1, x, y));
        P.push_back(computeFormula(1, 2, 2, x, y));
        P.push_back(computeFormula(2, 3, 1, x, y));
        P.push_back(computeFormula(2, 3, 2, x, y));
        P.push_back(computeFormula(1, 3, 1, x, y));
        P.push_back(computeFormula(1, 3, 2, x, y));
        
        // Вычисление центра тяжести
        double x10 = (v1.x + v2.x + v3.x) / 3.0;
        double y10 = (v1.y + v2.y + v3.y) / 3.0;
        Point center(x10, y10, 0); // Создание вершины центра тяжести
        P.push_back(center); // Добавляем вершину центра тяжести в конец вектора P
    }

    double computeL(Point A, Point B, Point XY){
        double u1=A.x, v1=A.y;
        double u2=B.x, v2=B.y;
        double x=XY.x, y=XY.y;
        return (x-u1)*(v2-v1)-(y-v1)*(u2-u1);
    }
    
    double computePsi(int i, Point XY){
        switch (i) {
            case 1:
                return computeL(P[2], P[3], XY);
                break;
            case 2:
                return computeL(P[1], P[3], XY);
                break;
            case 3:
                return computeL(P[1], P[2], XY);
                break;
            case 4:
                return computeL(P[4], P[8], XY);
                break;
            case 5:
                return computeL(P[5], P[9], XY);
                break;
            case 6:
                return computeL(P[5], P[6], XY);
                break;
            case 7:
                return computeL(P[4], P[7], XY);
                break;
            case 8:
                return computeL(P[7], P[9], XY);
                break;
            case 9:
                return computeL(P[6], P[8], XY);
                break;
            default:
                return 0;
                break;
        }
    }
    
    double multiply3Psi(int a, int b, int c, Point XY){
        return computePsi(a, XY)*computePsi(b, XY)*computePsi(c, XY);
    }
    
    double computeNu(int i, Point XY){
        switch (i) {
            case 1:
                return multiply3Psi(1, 4, 5, XY);
                break;
            case 2:
                return multiply3Psi(2, 6, 7, XY);
                break;
            case 3:
                return multiply3Psi(3, 8, 9, XY);
                break;
            case 4:
                return multiply3Psi(1, 2, 5, XY);
                break;
            case 5:
                return multiply3Psi(1, 2, 7, XY);
                break;
            case 6:
                return multiply3Psi(2, 3, 7, XY);
                break;
            case 7:
                return multiply3Psi(2, 3, 9, XY);
                break;
            case 8:
                return multiply3Psi(1, 3, 5, XY);
                break;
            case 9:
                return multiply3Psi(1, 3, 9, XY);
                break;
            case 10:
                return multiply3Psi(1, 2, 3, XY);
                break;
            default:
                return 0;
                break;
        }
    }
    
    double computeFi(int i, Point XY){
        return computeNu(i, XY)/computeNu(i, P[i]);
    }
    
    double Pf(Point XY){
        double res=0;
        for(int i=1; i<=10; i++){
            res=res+fn(P[i])*computeFi(i, XY);
        }
        return res;
    }
    
};

//Важно что каждый прямоугольник записывается как
//v1      v2
//v4      v3
class Rectangle{
public:
    Point v1, v2, v3, v4;
    Rectangle(){};
    Rectangle(Point a, Point b){
        v1=a; v3=b;
        v2.x=a.x; v2.y=b.y;
        v4.x=b.x; v4.y=a.y;
    }
    Rectangle(double x1, double y1, double x2, double y2){
        Point a(x1, y1), b(x2, y2);
        v1=a; v3=b;
        v2.x=b.x; v2.y=a.y;
        v4.x=a.x; v4.y=b.y;
    }
    Rectangle(Point a, Point b, Point c, Point d){
        v1 = a; v2 = b; v3 = c; v4 = d;
    }
};

vector<Rectangle> splitRectangleToSmallRectagles(Rectangle A, int n);
vector<Triangle> splitRectanglesToTriangles(vector<Rectangle> Rectangle);

vector<Rectangle> splitRectangleToSmallRectagles(Rectangle A, int n){
    vector<Rectangle> Rectagles;
    vector<vector<Point> > Net;
    double hx, hy;
    hx=(A.v2.x-A.v1.x)/n; //шаг по x
    hy=(A.v1.y-A.v4.y)/n; //шаг по y
    for(int i=0; i<=n; i++){
        Net.resize(Net.size()+1);
        for(int j=0; j<=n; j++){
            Point P;
            P.x=A.v4.x+j*hx;
            P.y=A.v4.y+i*hy;
            Net[i].push_back(P);
        }
    }
    for(int i=0; i<=n-1; i++){
        for(int j=0; j<=n-1; j++){
            Point v1, v2, v3, v4;
            v4=Net[i][j];
            v1=Net[i+1][j];
            v2=Net[i+1][j+1];
            v3=Net[i][j+1];
            Rectangle help(v1, v2, v3, v4);
            Rectagles.push_back(help);
        }
    }
    return Rectagles;
}

vector<Triangle> splitRectanglesToTriangles(vector<Rectangle> Rectangle){
    vector<Triangle> Triangles;
    for(int i=0; i<Rectangle.size(); i++){
        Triangle A(Rectangle[i].v1, Rectangle[i].v2, Rectangle[i].v4);
        Triangle B(Rectangle[i].v3, Rectangle[i].v2, Rectangle[i].v4);
        Triangles.push_back(A);
        Triangles.push_back(B);
    }
    return Triangles;
}

void writePointsToFile(vector<Point> Points, const string& filename){
    ofstream file(filename);
    for(int i=0; i<Points.size(); i++){
        file << Points[i].x << " " << Points[i].y << " " << Points[i].z << endl;
    }
}


int main(){
    double x1=0, y1=4, x2=4, y2=0; //исходный прямоугольник
    double x3=1, y3=2, x4=2, y4=1; //Вырезанный прямоугольник
    
    //Создаем вырез
    Rectangle R1(x1, y1, x2, y3);
    Rectangle R2(x1, y3, x3, y4);
    Rectangle R3(x4, y3, x2, y4);
    Rectangle R4(x1, y4, x2, y2);
    
    vector<Rectangle> Rectangles1 = splitRectangleToSmallRectagles(R1, 4);
    vector<Rectangle> Rectangles2 = splitRectangleToSmallRectagles(R2, 4);
    vector<Rectangle> Rectangles3 = splitRectangleToSmallRectagles(R3, 4);
    vector<Rectangle> Rectangles4 = splitRectangleToSmallRectagles(R4, 4);
    
    //Добавляем все 4 прямоугольных разбиения в один список
    vector<Rectangle> Rectangles;
    Rectangles=Rectangles1;
    Rectangles.insert(Rectangles.end(), Rectangles2.begin(), Rectangles2.end());
    Rectangles.insert(Rectangles.end(), Rectangles3.begin(), Rectangles3.end());
    Rectangles.insert(Rectangles.end(), Rectangles4.begin(), Rectangles4.end());
    
    //Треангуляция готова
    vector<Triangle> Triangles = splitRectanglesToTriangles(Rectangles);
    
    //вычисляем z для всей вершин треугольника
    for(int i=0; i<Triangles.size(); i++){
        Triangles[i].v1.computeZ();
        Triangles[i].v2.computeZ();
        Triangles[i].v3.computeZ();
    }
    Triangles[0].computeP1_10();
    writePointsToFile(Triangles[0].P, "test.txt");
    
    return 0;
}
