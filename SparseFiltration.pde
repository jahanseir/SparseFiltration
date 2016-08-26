Point[] pointlist;
PFont font;
float radius=50;
String dataFile = "data.txt";
String exportPointsFile = "../exportedPoints.txt";
String filtrationFile = "../filtration.txt";
String imgPath="Images/img-#######.png";
float lastMouseX=0, lastMouseY=0;
float rotateX=0, rotateY=0, rotateZ=0;
float diffX=0, diffY=0;
int inc=10;
float epsilon=0.9;
float scaleFactor=1;
boolean showSimplicies=false, showBalls=false;
boolean showoneSimplicies=false, showtwoSimplicies=false, showthreeSimplicies=false;
boolean sparseRipsMode=false, print=false, saveFrames=false, exportPoints=false;
boolean isMouseReleased=true;
boolean autoRotateX=false, autoRotateY=false, autoIncrease=false, movieMode=false;
SimplicialComplex cechComplex, filteredRipsComplex;
GreedyPermutation gperm;
color regBallsColor=color(128, 255, 0);
color fixBallsColor=color(0, 0, 255);
color smallGrowingBallsColor=color(255, 0, 0);
color triangleColor=color(0, 255, 0);
int secondsNo=10, frameNo=30;
Boolean increasingMovieMode=false;

Ball[] balls;

class Point {
  float x, y, z, dist;
  int index;
  Point(float _x, float _y, float _z) {
    x=_x;
    y=_y;
    z=_z;
    dist=MAX_FLOAT;
  }

  Point(int _index, float _x, float _y, float _z) {
    this(_x, _y, _z);
    index=_index;
  }

  Point(float _x, float _y, float _z, float _dist) {
    this(_x, _y, _z);
    dist=_dist;
  }

  Point(int _index, float _x, float _y, float _z, float _dist) {
    this(_x, _y, _z, _dist);
    index=_index;
  }

  float ComputeEuclideanDist(Point other) {
    return sqrt(pow(x-other.x, 2)+pow(y-other.y, 2)+pow(z-other.z, 2));
  }

  Point Minus(Point second) {
    return new Point(x-second.x, y-second.y, z-second.z);
  }

  Point Plus(Point second) {
    return new Point(x+second.x, y+second.y, z+second.z);
  }

  Point CrossProduct(Point second) {
    return new Point(y*second.z-z*second.y, z*second.x-x*second.z, x*second.y-y*second.x);
  }

  float DotProduct(Point second) {
    return x*second.x+y*second.y+z*second.z;
  }

  Point ScalarProduct(float c) {
    return new Point(c*x, c*y, c*z);
  }

  float Length() {
    return sqrt(x*x+y*y+z*z);
  }
}

class Ball {
  Point center;
  float radius;
  color customColor;
  Ball(Point c, float r) {
    center=c;
    radius=r;
    customColor=regBallsColor;
  }
  Ball(Point c, float r, color _color) {
    this(c, r);
    customColor=_color;
  }
}


class MinEnclosingBall {
  Ball Find(Point[] t, Point[] v) {
    ArrayList<Point> _t=new ArrayList<Point>();
    ArrayList<Point> _v=new ArrayList<Point>();
    for (Point p : t)
      _t.add(p);
    if (v!=null)
      for (Point p : v)
        _v.add(p);
    return FindMinBall(_t, _v);
  }

  Ball FindMinBall(ArrayList<Point> t, ArrayList<Point> v) {
    if (t.size()==0) 
      return FindCircumBall(v);
    int index=(int)random(t.size());
    Point u=t.get(index);
    ArrayList<Point> tempt=new ArrayList<Point>();
    for (int i=0; i<t.size (); i++)
      if (i!=index)
        tempt.add(t.get(i));
    Ball B=FindMinBall(tempt, v);
    if (B==null || B.center.ComputeEuclideanDist(u)>B.radius) {
      tempt=new ArrayList<Point>();
      for (int i=0; i<t.size (); i++)
        if (i!=index)
          tempt.add(t.get(i));
      ArrayList<Point> tempv=new ArrayList<Point>();
      for (int i=0; i<v.size (); i++)
        tempv.add(v.get(i));
      tempv.add(u);
      B=FindMinBall(tempt, tempv);
    }
    return B;
  }

  Ball FindCircumBall(ArrayList<Point> points) {
    switch(points.size()) {
    case 0:
      return null;
    case 1:
      return ComputeOnePointCircle(points.get(0));
    case 2:
      return ComputeTwoPointsCircle(points.get(0), points.get(1));
    case 3:
      return ComputeCircumcircle(points.get(0), points.get(1), points.get(2));
    case 4:
      return ComputeCircumsphere(points.get(0), points.get(1), points.get(2), points.get(3));
    default:
      throw new RuntimeException("An error occured");
    }
  }

  Ball ComputeOnePointCircle(Point p) {
    return new Ball(p, 0);
  }

  Ball ComputeTwoPointsCircle(Point p1, Point p2) {
    return new Ball(new Point((p1.x+p2.x)/2, (p1.y+p2.y)/2, (p1.z+p2.z)/2), p1.ComputeEuclideanDist(p2)/2);
  }

  Ball ComputeCircumcircle(Point p1, Point p2, Point p3) {
    Point p1p2=p1.Minus(p2);
    Point p2p3=p2.Minus(p3);
    Point p1p3=p1.Minus(p3);
    Point p2p1=p2.Minus(p1);
    Point p3p1=p3.Minus(p1);
    Point p3p2=p3.Minus(p2);
    float denom=p1p2.CrossProduct(p2p3).Length();
    if (denom==0) {
      throw new RuntimeException("Degenerate case, 3 colinear points: "+p1.index+","+p2.index+","+p3.index);
    }
    float radius=p1p2.Length()*p2p3.Length()*p3p1.Length()/(2*denom);
    float alpha, beta, gama;
    alpha=pow(p2p3.Length(), 2)*p1p2.DotProduct(p1p3)/(2*pow(denom, 2));
    beta=pow(p1p3.Length(), 2)*p2p1.DotProduct(p2p3)/(2*pow(denom, 2));
    gama=pow(p1p2.Length(), 2)*p3p1.DotProduct(p3p2)/(2*pow(denom, 2));
    return new Ball(p1.ScalarProduct(alpha).Plus(p2.ScalarProduct(beta)).Plus(p3.ScalarProduct(gama)), radius);
  }

  Ball ComputeCircumsphere(Point p1, Point p2, Point p3, Point p4) {
    float[][] a= {
      {
        p1.x, p1.y, p1.z, 1
      }
      , {
        p2.x, p2.y, p2.z, 1
      }
      , {
        p3.x, p3.y, p3.z, 1
      }
      , {
        p4.x, p4.y, p4.z, 1
      }
    };
    float detA=Determinant(a);
    if (detA==0) {
      println("Degenerate case, 4 coplanar points: "+p1.index+","+p2.index+","+p3.index+","+p4.index);
      return null;
    }

    float[][] c= {
      {
        pow(p1.x, 2)+pow(p1.y, 2)+pow(p1.z, 2), p1.x, p1.y, p1.z
      }
      , {
        pow(p2.x, 2)+pow(p2.y, 2)+pow(p2.z, 2), p2.x, p2.y, p2.z
      }
      , {
        pow(p3.x, 2)+pow(p3.y, 2)+pow(p3.z, 2), p3.x, p3.y, p3.z
      }
      , {
        pow(p4.x, 2)+pow(p4.y, 2)+pow(p4.z, 2), p4.x, p4.y, p4.z
      }
    };
    float[][] x= {
      {
        pow(p1.x, 2)+pow(p1.y, 2)+pow(p1.z, 2), p1.y, p1.z, 1
      }
      , {
        pow(p2.x, 2)+pow(p2.y, 2)+pow(p2.z, 2), p2.y, p2.z, 1
      }
      , {
        pow(p3.x, 2)+pow(p3.y, 2)+pow(p3.z, 2), p3.y, p3.z, 1
      }
      , {
        pow(p4.x, 2)+pow(p4.y, 2)+pow(p4.z, 2), p4.y, p4.z, 1
      }
    };
    float[][] y= {
      {
        pow(p1.x, 2)+pow(p1.y, 2)+pow(p1.z, 2), p1.x, p1.z, 1
      }
      , {
        pow(p2.x, 2)+pow(p2.y, 2)+pow(p2.z, 2), p2.x, p2.z, 1
      }
      , {
        pow(p3.x, 2)+pow(p3.y, 2)+pow(p3.z, 2), p3.x, p3.z, 1
      }
      , {
        pow(p4.x, 2)+pow(p4.y, 2)+pow(p4.z, 2), p4.x, p4.z, 1
      }
    };
    float[][] z= {
      {
        pow(p1.x, 2)+pow(p1.y, 2)+pow(p1.z, 2), p1.x, p1.y, 1
      }
      , {
        pow(p2.x, 2)+pow(p2.y, 2)+pow(p2.z, 2), p2.x, p2.y, 1
      }
      , {
        pow(p3.x, 2)+pow(p3.y, 2)+pow(p3.z, 2), p3.x, p3.y, 1
      }
      , {
        pow(p4.x, 2)+pow(p4.y, 2)+pow(p4.z, 2), p4.x, p4.y, 1
      }
    };
    float detX=Determinant(x), detY=-Determinant(y), detZ=Determinant(z), detC=Determinant(c);
    return new Ball(new Point(detX/(2*detA), detY/(2*detA), detZ/(2*detA)), sqrt(pow(detX, 2)+pow(detY, 2)+pow(detZ, 2)-4*detA*detC)/abs(2*detA));
  }
}

float Determinant(float[][] a) {
  int p, h, k, i, j, n=a.length;
  float det=0;
  float[][] temp=new float[n-1][n-1];
  if (n==1) {
    return a[0][0];
  } else if (n==2) {
    det=(a[0][0]*a[1][1]-a[0][1]*a[1][0]);
    return det;
  } else {
    for (p=0; p<n; p++) {
      h = 0;
      k = 0;
      for (i=1; i<n; i++) {
        for (j=0; j<n; j++) {
          if (j==p) {
            continue;
          }
          temp[h][k] = a[i][j];
          k++;
          if (k==n-1) {
            h++;
            k = 0;
          }
        }
      }
      det=det+a[0][p]*pow(-1, p)*Determinant(temp);
    }
    return det;
  }
}

class GreedyPermutation {
  Point[] points;
  GreedyPermutation(Point[] input) {
    points=new Point[input.length];
    for (int i=0; i<input.length; i++)
      points[i]=input[i];
  }
  void Compute() {
    Swap(0, (int)random(points.length));
    for (int i=0; i<points.length-1; i++) {
      int maxIndex = i+1;
      for (int j=i+1; j<points.length; j++) {
        float currentDist=points[j].ComputeEuclideanDist(points[i]);
        points[j].dist=min(currentDist, points[j].dist);
        if (points[maxIndex].dist < points[j].dist)
          maxIndex = j;
      }
      Swap(i+1, maxIndex);
    }
  }

  void Swap(int firstIndex, int secondIndex) {
    Point temp=points[firstIndex];
    points[firstIndex]=points[secondIndex];
    points[secondIndex]=temp;
  }
}

class Simplex {
  Point[] points;
  float birth;
  Simplex(Point[] _points, float _alpha) {
    points=_points;
    birth=_alpha;
  }
}

class SimplicialComplex {
  boolean foundSimplicies=false;
  ArrayList<Simplex> zeroSimplicies;
  ArrayList<Simplex> oneSimplicies;
  ArrayList<Simplex> twoSimplicies;
  ArrayList<Simplex> threeSimplicies;
  SimplicialComplex(Point[] points) {
    zeroSimplicies=new ArrayList<Simplex>();
    for (Point p : points) {
      zeroSimplicies.add(new Simplex(new Point[] {
        p
      }
      , 0));
    }
    oneSimplicies=new ArrayList<Simplex>();
    twoSimplicies=new ArrayList<Simplex>();
    threeSimplicies=new ArrayList<Simplex>();
  }

  void FindCechComplex() {
    MinEnclosingBall eball=new MinEnclosingBall();
    foundSimplicies=true;
    for (int i=0; i<zeroSimplicies.size ()-1; i++)
      for (int j=i+1; j<zeroSimplicies.size (); j++) {
        float dist=zeroSimplicies.get(i).points[0].ComputeEuclideanDist(zeroSimplicies.get(j).points[0]);
        oneSimplicies.add(new Simplex(new Point[] {
          zeroSimplicies.get(i).points[0], zeroSimplicies.get(j).points[0]
        }
        , dist/2));
      }
    for (int i=0; i<pointlist.length-2; i++)
      for (int j=i+1; j<pointlist.length-1; j++)
        for (int k=j+1; k<pointlist.length; k++) {
          Ball b=null;
          float distij=pointlist[i].ComputeEuclideanDist(pointlist[j]);
          float distik=pointlist[i].ComputeEuclideanDist(pointlist[k]);
          float distjk=pointlist[j].ComputeEuclideanDist(pointlist[k]);
          float theta=0;
          float rad;
          if (distij>=distik && distij>=distjk) {
            Point vectIK=pointlist[i].Minus(pointlist[k]);
            Point vectJK=pointlist[j].Minus(pointlist[k]);
            theta=acos(vectIK.DotProduct(vectJK)/(vectIK.Length()*vectJK.Length()));
            rad=distij/2;
          } else if (distik>=distij && distik>=distjk) {
            Point vectIJ=pointlist[i].Minus(pointlist[j]);
            Point vectKJ=pointlist[k].Minus(pointlist[j]);
            theta=acos(vectIJ.DotProduct(vectKJ)/(vectIJ.Length()*vectKJ.Length()));
            rad=distik/2;
          } else {
            Point vectJI=pointlist[j].Minus(pointlist[i]);
            Point vectKI=pointlist[k].Minus(pointlist[i]);
            theta=acos(vectJI.DotProduct(vectKI)/(vectJI.Length()*vectKI.Length()));
            rad=distjk/2;
          }
          if (theta>=PI/2)
            b=new Ball(null, rad);
          else {
            b=eball.ComputeCircumcircle(pointlist[i], pointlist[j], pointlist[k]);
          }
          if (b!=null) {
            twoSimplicies.add(new Simplex(new Point[] {
              pointlist[i], pointlist[j], pointlist[k]
            }
            , b.radius));
          }
        }
  }

  void FindWeightedBalls(float alpha, Ball[] balls) {
    for (Ball b : balls) {
      if (alpha<=b.center.dist*(1+epsilon)/epsilon) {
        b.customColor=regBallsColor;
        b.radius=alpha;
      } else if (alpha<=b.center.dist*pow((1+epsilon), 2)/epsilon) {
        b.customColor=fixBallsColor;
        b.radius=b.center.dist*(1+epsilon)/epsilon;
      } else {
        b.customColor=smallGrowingBallsColor;
        b.radius=0;
      }
    }
  }

  void FindWeightedRips() {
    for (int i=pointlist.length-1; i>=0; i--) {
      ArrayList<Integer> neighbors=new ArrayList<Integer>(); 
      for (int j=i-1; j>=0; j--) {
        float appTime=computeAppearanceTime(pointlist[i], pointlist[j]);
        if (appTime>=0) {
          oneSimplicies.add(new Simplex(new Point[] {
            pointlist[i], pointlist[j]
          }
          , appTime));
          neighbors.add(j);
        }
      }
      for (int k=0; k<neighbors.size (); k++)
        for (int m=k+1; m<neighbors.size (); m++) {
          float appTime=computeAppearanceTime(pointlist[neighbors.get(k)], pointlist[neighbors.get(m)]);
          if (appTime>=0) {
            appTime=max(computeSimplexApprearanceTime(new Point[] {
              pointlist[i], pointlist[neighbors.get(k)], pointlist[neighbors.get(m)]
            }
            ));
            if (appTime>=0) {
              twoSimplicies.add(new Simplex(new Point[] {
                pointlist[i], pointlist[neighbors.get(k)], pointlist[neighbors.get(m)]
              }
              , appTime));
            }
          }
        }
    }
  }

  float[] computeSimplexApprearanceTime(Point[] vertices) {
    float[] appTime=new float[vertices.length*(vertices.length-1)/2];
    int k=0;
    for (int i=0; i<vertices.length-1; i++)
      for (int j=i+1; j<vertices.length; j++) {
        appTime[k]=computeAppearanceTime(vertices[i], vertices[j]);
        k++;
      }
    return appTime;
  }

  float computeAppearanceTime(Point p, Point q)
  {
    float dist=p.ComputeEuclideanDist(q);
    float l1=min(p.dist, q.dist), l2=max(p.dist, q.dist);
    float alpha;

    //r1(a)=a && r2(a)=a
    alpha=dist/2;
    if (alpha<=l1*(1+epsilon)/epsilon) return alpha;

    //r1(a)=l1(1+e)/e && r2(a)=a
    alpha=dist-(l1*(1+epsilon)/epsilon);
    if (alpha<=l2*(1+epsilon)/epsilon) return alpha;

    //r1(a)=l1(1+e)/e && r2(a)=l2(1+e)/e
    if (dist==(l1+l2)*(1+epsilon)/epsilon) return l2*(1+epsilon)/epsilon;

    return -1;
  }

  void ShowOneSimplicies(float alpha, color c) {
    for (Simplex s : oneSimplicies)
      if (s.birth<=alpha) {
        stroke(c);
        strokeWeight(3);
        line(s.points[0].x, s.points[0].y, s.points[0].z, s.points[1].x, s.points[1].y, s.points[1].z);
      }
  }

  void ShowTwoSimplicies(float alpha, color c) {
    for (Simplex s : twoSimplicies)
      if (s.birth<=alpha) {
        noStroke();
        fill(c);
        beginShape();
        vertex(s.points[0].x, s.points[0].y, s.points[0].z);
        vertex(s.points[1].x, s.points[1].y, s.points[1].z);
        vertex(s.points[2].x, s.points[2].y, s.points[2].z);
        endShape();
      }
  }

  void ShowThreeSimplicies(float alpha, color c) {
    for (Simplex s : threeSimplicies)
      if (s.birth<=alpha) {
        noStroke();
        fill(c);
        beginShape();
        vertex(s.points[0].x, s.points[0].y, s.points[0].z);
        vertex(s.points[1].x, s.points[1].y, s.points[1].z);
        vertex(s.points[2].x, s.points[2].y, s.points[2].z);
        endShape();
        beginShape();
        vertex(s.points[2].x, s.points[2].y, s.points[2].z);
        vertex(s.points[1].x, s.points[1].y, s.points[1].z);
        vertex(s.points[3].x, s.points[3].y, s.points[3].z);
        endShape();
        beginShape();
        vertex(s.points[0].x, s.points[0].y, s.points[0].z);
        vertex(s.points[1].x, s.points[1].y, s.points[1].z);
        vertex(s.points[3].x, s.points[3].y, s.points[3].z);
        endShape();
        beginShape();
        vertex(s.points[0].x, s.points[0].y, s.points[0].z);
        vertex(s.points[2].x, s.points[2].y, s.points[2].z);
        vertex(s.points[3].x, s.points[3].y, s.points[3].z);
        endShape();
      }
  }

  void Print(float alpha) {
    PrintWriter output=createWriter(filtrationFile);
    println("One Simplicies:");
    for (Simplex s : zeroSimplicies) {
      output.println(s.birth+"\t"+s.points[0].index);
    }
    for (Simplex s : oneSimplicies) {
      output.println(s.birth+"\t"+s.points[0].index+","+s.points[1].index);
      if (s.birth<=alpha)
        println("edge: "+s.points[0].index+"<->"+s.points[1].index+", birth: "+s.birth);
    }
    println("----------------------");
    println("Two Simplicies:");
    for (Simplex s : twoSimplicies) {
      output.println(s.birth+"\t"+s.points[0].index+","+s.points[1].index+","+s.points[2].index);
      if (s.birth<=alpha)
        println("triangle: "+s.points[0].index+"<->"+s.points[1].index+"<->"+s.points[2].index+", birth: "+s.birth);
    }
    println("----------------------");
    println("Three Simplicies:");
    for (Simplex s : threeSimplicies) {
      output.println(s.birth+"\t"+s.points[0].index+","+s.points[1].index+","+s.points[2].index+","+s.points[3].index);      
      if (s.birth<=alpha)
        println("tetrahedron: "+s.points[0].index+"<->"+s.points[1].index+"<->"+s.points[2].index+"<->"+s.points[3].index+"; birth: "+s.birth);
    }
    output.flush();
    output.close();
  }

  int GetSizeOfOneSimplicies(float alpha) {
    int count=0;
    for (Simplex s : oneSimplicies)
      if (s.birth<=alpha)
        count++;
    return count;
  }

  int GetSizeOfTwoSimplicies(float alpha) {
    int count=0;
    for (Simplex s : twoSimplicies)
      if (s.birth<=alpha)
        count++;
    return count;
  }

  int GetSizeOfThreeSimplicies(float alpha) {
    int count=0;
    for (Simplex s : threeSimplicies)
      if (s.birth<=alpha)
        count++;
    return count;
  }
}

void loadPoints() {
  String[] input = loadStrings(dataFile);
  if (input.length > 0) {
    pointlist=new Point[input.length];
    for (int i = 0; i < input.length; i++) {
      String[] line = splitTokens(input[i], "\t");
      pointlist[i]=new Point((i+1), parseFloat(line[0].trim()), parseFloat(line[1].trim()), parseFloat(line[2].trim()), parseFloat(line[3].trim()));
    }
  }
}

void setup() {
  size(1280, 720, P3D);
  font = createFont("Arial", 14); 
  loadPoints();
  gperm=new GreedyPermutation(pointlist);
  gperm.Compute();
  for (Point p : pointlist)
    println(p.index+"\t"+p.x+"\t"+p.y+"\t"+p.z+"\t"+p.dist);
  balls=new Ball[pointlist.length];
  for (int i=0; i<pointlist.length; i++)
    balls[i]=new Ball(pointlist[i], 0);
  cechComplex=new SimplicialComplex(pointlist);
  cechComplex.FindCechComplex();
  filteredRipsComplex=new SimplicialComplex(pointlist);
  filteredRipsComplex.FindWeightedRips();
}

void draw() {
  boolean stopMovie=false;
  background(255);
  lights();
  textFont(font);       
  fill(0);  
  if (autoIncrease) {
    if (radius==100) {
      radius=0;
    }
    radius++;
  }
  if (movieMode) {
    rotateX=diffX*2*PI/(secondsNo*frameNo);
    diffX++;
    if (diffX>=secondsNo*frameNo)
      stopMovie=true;    
    rotateY=-50;
    if (increasingMovieMode)
      radius+=0.35;
  }
  if (autoRotateX) {
    rotateX=(diffX=diffX+inc)/width*2*PI;
    rotateY=-50;
  }
  if (autoRotateY)
    rotateY=(diffY=diffY+inc)/height*2*PI;
  float alpha=radius;
  SimplicialComplex sc;
  if (sparseRipsMode) {
    sc=filteredRipsComplex;
    sc.FindWeightedBalls(radius, balls);
  } else {
    sc=cechComplex;
    for (Ball b : balls) {
      b.radius=radius;
      b.customColor=regBallsColor;
    }
  }
  if (!movieMode) {
    text("Radius:"+radius, 10, 20, 0);
    text("0-simplicies:"+pointlist.length, 10, 40, 0);
    text("1-simplicies:"+sc.GetSizeOfOneSimplicies(radius), 10, 60, 0);
    text("2-simplicies:"+sc.GetSizeOfTwoSimplicies(radius), 10, 80, 0);
    text("3-simplicies:"+sc.GetSizeOfThreeSimplicies(radius), 10, 100, 0);
  }
  fill(128, 255, 0, 128);
  noStroke();
  translate(width/2, height/2, 0);  
  if (mousePressed) {
    if (isMouseReleased) {
      lastMouseX=mouseX;
      lastMouseY=mouseY;
      isMouseReleased=false;
    }
    diffX+=mouseX-lastMouseX;
    diffY+=mouseY-lastMouseY;
    rotateX=diffX/width*2*PI;
    rotateY=diffY/height*2*PI;
    lastMouseX=mouseX;
    lastMouseY=mouseY;
  }
  rotateX(-rotateY);  
  rotateY(rotateX);  
  rotateZ(rotateZ);
  scale(scaleFactor);

  for (int i = 0; i < pointlist.length; ++i) {
    fill(0, 0, 0);
    pushMatrix();
    translate(pointlist[i].x, pointlist[i].y, pointlist[i].z);
    sphere(2);
    popMatrix();
  }
  if (showBalls)
    for (Ball b : balls) {
      fill(b.customColor);
      pushMatrix();
      translate(b.center.x, b.center.y, b.center.z);
      sphere(b.radius);
      popMatrix();
    }    
  if (showSimplicies) {
    if (showoneSimplicies) 
      sc.ShowOneSimplicies(alpha, color(0, 0, 0));
    if (showtwoSimplicies) 
      sc.ShowTwoSimplicies(alpha, triangleColor);
    if (showthreeSimplicies) 
      sc.ShowThreeSimplicies(alpha, color(256, 0, 0, 128));
  }
  if (print) {
    sc.Print(alpha);
    print=false;
  }
  if (saveFrames)
    saveFrame(imgPath);
  if (stopMovie) {
    saveFrames=movieMode=false;
    diffX=0;
  }
  if (exportPoints) {
    PrintWriter output=createWriter(exportPointsFile);
    for (Point p : pointlist)
      output.print(p.x+","+p.y+","+p.z+";");
    output.flush();
    output.close();
    exportPoints=false;
  }
}

void mouseWheel(MouseEvent e) {
  float temp = scaleFactor-e.getAmount() / 50;
  if (temp<=0.01)
    return;
  scaleFactor=temp;
}

void mouseReleased() {
  isMouseReleased=true;
}

void keyPressed() {  
  if (key=='i') {
    radius++;
  }
  if (key=='d') {
    if (radius>1)
      radius--;
  }
  if (key=='x') {
    autoRotateX=!autoRotateX;
  }
  if (key=='y') {
    autoRotateY=!autoRotateY;
  }
  if (key=='s') {
    if (!showSimplicies)
      showoneSimplicies=showtwoSimplicies=showthreeSimplicies=true;
    showSimplicies=!showSimplicies;
  }
  if (key=='b') {
    showBalls=!showBalls;
  }
  if (key=='1') {
    showoneSimplicies=!showoneSimplicies;
  }
  if (key=='2') {
    showtwoSimplicies=!showtwoSimplicies;
  }
  if (key=='3') {
    showthreeSimplicies=!showthreeSimplicies;
  }
  if (key=='p') {
    print=true;
  }
  if (key=='r') {
    sparseRipsMode=!sparseRipsMode;
  }
  if (key=='f') {
    saveFrames=!saveFrames;
  }
  if (key=='a') {    
    autoIncrease=!autoIncrease;
    if (autoIncrease)
      radius=0;
    else
      radius=50;
  }
  if (key=='m') {
    diffX=0;
    movieMode=!movieMode;
    if (movieMode) {
      if (increasingMovieMode)
        radius=0;
      else
        radius=50;
      saveFrames=true;
    } else {
      radius=50;
      saveFrames=false;
    }
  }
  if (key=='e') {
    exportPoints=true;
  }
}