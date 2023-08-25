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
#include <random>
#include "SolveByGaus.hpp"

using namespace std;

double fn(double x, double  y){
    return sin(x*y);
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

double computeDistance(Point A, Point B){
    return sqrt((A.x-B.x)*(A.x-B.x)+(A.y-B.y)*(A.y-B.y));
}

double fn(Point A){
    return fn(A.x, A.y);
}

class Triangle{
public:
    Point v1, v2, v3;
    
    vector<Point> netPointsInside;
    
    double square;
    
    Point centerOfMass;
    
    double computeSquare(Point A, Point B, Point C){
        double a, b, c, p;
        a=computeDistance(A, B);
        b=computeDistance(C, B);
        c=computeDistance(A, C);
        p=(a+b+c)/2;
        return sqrt(p*(p-a)*(p-b)*(p-c));
    }
    
    void computeSquare(){
        Point A, B, C;
        A=v1; B=v2; C=v3;
        double a, b, c, p;
        a=computeDistance(A, B);
        b=computeDistance(C, B);
        c=computeDistance(A, C);
        p=(a+b+c)/2;
        square = sqrt(p*(p-a)*(p-b)*(p-c));
        
    }
    
    void computeCenterofMass(){
        centerOfMass.x=(v1.x+v2.x+v3.x)/3;
        centerOfMass.y=(v1.y+v2.y+v3.y)/3;
    }
    
    vector<Point> P;
    
    Triangle(){};
    Triangle(Point a, Point b, Point c){
        v1=a; v2=b; v3=c;
        square=computeSquare(a, b, c);
    }
    
    Triangle(double x1, double y1, double x2, double y2, double x3, double y3){
        v1.x = x1; v1.y = y1;
        v2.x = x2; v2.y = y2;
        v3.x = x3; v3.y = y3;
        square=computeSquare(v1, v2, v3);
    }
    
    bool isPointInside(Point A){
        double epsilon=0.000001;
        double a, b, c, S;
        a=computeSquare(A, v1, v2);
        b=computeSquare(A, v1, v3);
        c=computeSquare(A, v3, v2);
        S=a+b+c;
        if(abs(S-square)>epsilon) return false;
        else return true;
    }
    
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
    cout << Points.size() << " Points writen to file" << endl;
}

vector<Point> MakeNet(double x1, double y1, double x2, double y2, int n){
    vector<Point> Net;
    double hx, hy;
    hx=(x2-x1)/n; //шаг по x
    hy=(y1-y2)/n; //шаг по y
    for(int i=0; i<=n; i++){
        for(int j=0; j<=n; j++){
            Point P;
            P.x=x1+j*hx;
            P.y=y2+i*hy;
            Net.push_back(P);
        }
    }
    return Net;
}

vector<Point> generateRandomPointsInsideTriangle(Triangle triangle, int n) {
    vector<Point> randomPoints;
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis(0.0, 1.0);

    for (int i = 0; i < n; ++i) {
        double r1 = dis(gen);
        double r2 = dis(gen);

        if (r1 + r2 <= 1.0) {
            double x = triangle.v1.x + r1 * (triangle.v2.x - triangle.v1.x) +
                       r2 * (triangle.v3.x - triangle.v1.x);
            double y = triangle.v1.y + r1 * (triangle.v2.y - triangle.v1.y) +
                       r2 * (triangle.v3.y - triangle.v1.y);
            randomPoints.push_back({x, y});
        }
        else i--;
    }
    
    return randomPoints;
}

vector<Triangle> Make_Triangulation_Delane(Triangle MasterTriangle, vector<Point> points){
    vector<Triangle> triangles;
    //Найти треугольник которому принадлежит точка
    //Созать 4 новых треугольника по 4м точкам
    //Добавить 4 треугольника в масив
    triangles.push_back(MasterTriangle);
    for(int i=0; i<points.size(); i++){
        for(int j=0; j<triangles.size(); j++){
            if(triangles[j].isPointInside(points[i])){
                Triangle T=triangles[j];
                Triangle A(T.v1, T.v2, points[i]);
                Triangle B(T.v1, points[i], T.v3);
                Triangle C(points[i], T.v2, T.v3);
                triangles.erase(triangles.begin() + j);
                triangles.push_back(A);
                triangles.push_back(B);
                triangles.push_back(C);
                break;
            }
        }
    }
    
    return triangles;
}

double computeIntergralOnTriagle(Triangle T, int accuracy){
    vector<Point> Points = generateRandomPointsInsideTriangle(T, accuracy);
    vector<Triangle> Triangulation = Make_Triangulation_Delane(T, Points);
    for(int i=0; i<Triangulation.size(); i++){
        Triangulation[i].computeCenterofMass();
        Triangulation[i].centerOfMass.computeZ();
        Triangulation[i].computeSquare();
    }
    double Integral=0;
    
    for(int i=0; i<Triangulation.size(); i++){
        Integral += Triangulation[i].square*Triangulation[i].centerOfMass.z;
    }
    return Integral;
}

double multipleTwoFi(int i, int j, Point A, Triangle T){
    
    return T.computeFi(i, A)*T.computeFi(j, A);
}

double multipleFonFi(int i, Point A, Triangle T){
    return fn(A)*T.computeFi(i, A);
}

//Вычисляем левый интегралл в системе, для этого разбиваем треугольник на маленькие и умножаем площадь каждого на центр масс
double computeIntergralOnTriagle(Triangle T, int accuracy, int i, int j){
    vector<Point> Points = generateRandomPointsInsideTriangle(T, accuracy);
    vector<Triangle> Triangulation = Make_Triangulation_Delane(T, Points);
    
    for(int h=0; h<Triangulation.size(); h++){
        Triangulation[h].P=T.P;
    }
    
    for(int t=0; t<Triangulation.size(); t++){
        Triangulation[t].computeCenterofMass();
        Triangulation[t].centerOfMass.z = multipleTwoFi(i, j, Triangulation[t].centerOfMass, Triangulation[t]);
        Triangulation[t].computeSquare();
    }
    
    double Integral=0;
    
    for(int i=0; i<Triangulation.size(); i++){
        Integral += Triangulation[i].square*Triangulation[i].centerOfMass.z;
    }
    
    return Integral;
}

double computeIntergralOnTriagle(Triangle T, int accuracy, int i){
    vector<Point> Points = generateRandomPointsInsideTriangle(T, accuracy);
    vector<Triangle> Triangulation = Make_Triangulation_Delane(T, Points);
    
    for(int h=0; h<Triangulation.size(); h++){
        Triangulation[h].P=T.P;
    }
    
    for(int t=0; t<Triangulation.size(); t++){
        Triangulation[t].computeCenterofMass();
        Triangulation[t].centerOfMass.z = multipleFonFi(i, Triangulation[t].centerOfMass, Triangulation[t]);
        Triangulation[t].computeSquare();
    }
    
    double Integral=0;
    
    for(int i=0; i<Triangulation.size(); i++){
        Integral += Triangulation[i].square*Triangulation[i].centerOfMass.z;
    }
    
    return Integral;
}

vector<vector<double>> computeAlphaMatrix(vector<Triangle> Triangulation, int accuracy){
    vector<vector<double>> alphaMatrix;
    for(int i=1; i<=10; i++){
        vector<double> help;
        for(int j=1; j<=10; j++){
            double a=0;
            for(int l=0; l<Triangulation.size(); l++){
                a+=computeIntergralOnTriagle(Triangulation[l], accuracy, i, j);
            }
            help.push_back(a);
        }
        alphaMatrix.push_back(help);
    }
    return alphaMatrix;
}

vector<double> computeBetaMatrix(vector<Triangle> Triangulation, int accuracy){
    vector<double> betaMatrix;
    for(int i=1; i<=10; i++){
        double a=0;
        for(int l=0; l<Triangulation.size(); l++){
            a+=computeIntergralOnTriagle(Triangulation[l], accuracy, i);
        }
        betaMatrix.push_back(a);
    }
    return betaMatrix;
}


int main(){
    double x1=-10, y1=10, x2=10, y2=-10; //исходный прямоугольник
    double x3=0, y3=5, x4=5, y4=0; //Вырезанный прямоугольник
    int NumPointsInNet=100; //Параметр отвечающий за густоту сетки
    int NumTriangles=5; //Параметр отвечающий за колличестов треугольников
    int IntAccuracy=2000; //точность взятия интеграла
    
    //Проба работы интегралла
    Triangle I1(0, 0, 0, 2, 2, 0);
    Triangle I2(2, 2, 0, 2, 2, 0);
    cout << computeIntergralOnTriagle(I1, 100)+computeIntergralOnTriagle(I2, 100) << endl;
    cout << computeIntergralOnTriagle(I1, 1000)+computeIntergralOnTriagle(I2, 1000) << endl;
    
    //Создаем вырез
    Rectangle R1(x1, y1, x2, y3);
    Rectangle R2(x1, y3, x3, y4);
    Rectangle R3(x4, y3, x2, y4);
    Rectangle R4(x1, y4, x2, y2);
    
    vector<Rectangle> Rectangles1 = splitRectangleToSmallRectagles(R1, NumTriangles);
    vector<Rectangle> Rectangles2 = splitRectangleToSmallRectagles(R2, NumTriangles);
    vector<Rectangle> Rectangles3 = splitRectangleToSmallRectagles(R3, NumTriangles);
    vector<Rectangle> Rectangles4 = splitRectangleToSmallRectagles(R4, NumTriangles);
    
    //Добавляем все 4 прямоугольных разбиения в один список
    vector<Rectangle> Rectangles;
    Rectangles=Rectangles1;
    Rectangles.insert(Rectangles.end(), Rectangles2.begin(), Rectangles2.end());
    Rectangles.insert(Rectangles.end(), Rectangles3.begin(), Rectangles3.end());
    Rectangles.insert(Rectangles.end(), Rectangles4.begin(), Rectangles4.end());
    
    //Треангуляция готова
    vector<Triangle> Triangles = splitRectanglesToTriangles(Rectangles);
    
    Rectangle NewR(x1, y1, x2, y2);
    vector<Rectangle> RectanglesNew = splitRectangleToSmallRectagles(NewR, NumTriangles);
    //vector<Triangle> Triangles = splitRectanglesToTriangles(RectanglesNew);
    cout << Triangles.size() << " Triangles" << endl;
    
    //вычисляем z для всех вершин треугольника
    for(int i=0; i<Triangles.size(); i++){
        Triangles[i].v1.computeZ();
        Triangles[i].v2.computeZ();
        Triangles[i].v3.computeZ();
    }
    
    //Вычисляем для каждого треугольника площадь
    for(int i=0; i<Triangles.size(); i++){
        Triangles[i].square=Triangles[i].computeSquare(Triangles[i].v1, Triangles[i].v2, Triangles[i].v3);
    }
    
    //Создаем сетку для отрисовки
    vector<Point> Net = MakeNet(x1, y1, x2, y2, NumPointsInNet);
    
    //Находим для каждой точки треугольник, которому принадлежит точка
    for(int i=0; i<Net.size(); i++){
        for(int j=0; j<Triangles.size(); j++){
            if(Triangles[j].isPointInside(Net[i])){
                Triangles[j].netPointsInside.push_back(Net[i]);
                break;
            }
        }
    }
    
    //Вычисляем для каждого треугольника P1_10 и площадь
    for(int i=0; i<Triangles.size(); i++){
        Triangles[i].computeP1_10();
    }
    
    //Вычисляем Значения Альфа
    vector<vector<double>> FirstM;
    vector<double> SecondM;
    FirstM = computeAlphaMatrix(Triangles, IntAccuracy);
    SecondM = computeBetaMatrix(Triangles, IntAccuracy);
    
    vector<double> alphas;
    alphas.push_back(0);
    alphas = solveEquations(FirstM, SecondM);
    //Вычисляем интерпалиционные значения
    vector<Point> AproxPoints;
    
    for(int i=0; i<Triangles.size(); i++){
        for(int j=0; j<Triangles[i].netPointsInside.size(); j++){
            Point XY = Triangles[i].netPointsInside[j];
            XY.z=0;
            for(int k=1; k<=10; k++){
                XY.z+=Triangles[i].computeFi(k, XY)*alphas[k];
            }
            XY.z=Triangles[i].Pf(XY);
            AproxPoints.push_back(XY);
        }
    }
    
    
    //Вычисляем реальные значения
    vector<Point> RealPoints;
    
    for(int i=0; i<AproxPoints.size(); i++){
        Point A;
        A=AproxPoints[i];
        A.z=fn(A);
        RealPoints.push_back(A);
    }
    
    
    writePointsToFile(AproxPoints, "aprox_points.txt");
    writePointsToFile(RealPoints, "real_points.txt");
    
    vector<Point> GoraGovnaPoints;
    for(int i=0; i<AproxPoints.size(); i++){
        Point A;
        A.x=AproxPoints[i].x;
        A.y=AproxPoints[i].y;
        A.z=AproxPoints[i].z-RealPoints[i].z;
        GoraGovnaPoints.push_back(A);
    }
    
    writePointsToFile(GoraGovnaPoints, "GoraGovna_points.txt");
    
    
    return 0;
}
