import processing.core.*; 
import processing.data.*; 
import processing.event.*; 
import processing.opengl.*; 

import java.util.Collections; 
import java.util.Arrays; 
import java.util.Comparator; 
import megamu.mesh.*; 

import java.util.HashMap; 
import java.util.ArrayList; 
import java.io.File; 
import java.io.BufferedReader; 
import java.io.PrintWriter; 
import java.io.InputStream; 
import java.io.OutputStream; 
import java.io.IOException; 

public class watershed extends PApplet {







/* *********************** HIGH RES EXPORT FUNCTIONALITY **************************/

//INCLUDE THESE GLOBAL VARIABLES  
PGraphics render;
PImage img;
String saveFilePath = "../outputs/tattoo_sketch-" + new java.text.SimpleDateFormat("yyyyMMdd-HHmmss").format(new java.util.Date()); //change the XXXX to current project name 
int printWidth = 6;
int printHeight = 6;
int printDpi = 300;
int previewDpi = 72;
boolean renderHighRes = true;
boolean firstFrame = true;
int renderWidth;
int renderHeight;
float scaleFactor = 1;
int seed = 0;

int[] line_palette, background_palette;
ArrayList<PVector> line;
//initialize 
AttractorSystem as;


int INIT_maxChildren = 3 ;

ArrayList<Polygon> polygons;
boolean toggle = true;
ArrayList<Polygon> newPolygons = new ArrayList(); //this is a carrier vessel for updating 
ArrayList<Polygon> deadPolygons = new ArrayList(); //this is used to hold the polygons that die during updating for reconstruction later
int num_beams = 4;
boolean showNormals=false;
int maxFrameLines =2;
int pointCount = 0;
int generation = 0;
int webcast = 4;
float polygonAreaCutoff = 0.001f;
PVector[][] noiseGrid; //set up a noise based flow field for changing tile attributes
Gradient colorGrad;
Polygon poly;


public void setup(){
    
    background(255);
    colorMode(HSB, 360, 100, 100, 100);
    doReset();
}

public ArrayList<PVector> inkscapePathImport(float[][] p, float inputWidth, float inputHeight){
    ArrayList<PVector> out = new ArrayList();
    for(int i=0; i<p.length; i++){
        // println(p[i][0][0]);
        float x = map(p[i][0], 0, inputWidth, 0, renderWidth);
        float y = map(p[i][1], 0, inputHeight, 0 ,renderHeight);
        out.add(new PVector(x, y));
    }
    return out;
}


public void doReset() { //initial setup and used in resetting for high-def export

    int dpi = renderHighRes ? printDpi : previewDpi;
    scaleFactor = dpi / (float)previewDpi;
    renderWidth = printWidth * dpi;
    renderHeight = printHeight * dpi;
    render = createGraphics(renderWidth, renderHeight);
    firstFrame = true;
    noiseSeed(seed);
    randomSeed(seed);
    background_palette = new int[]{color(0xff0f0f0e), color(0xff382a04), color(0xff141524), color(0xff170d1f), color(0xff000000)};
    line_palette = new int[]{color(0xff382a04), color(0xff594a1f), color(0xff073610), color(0xff18361e), color(0xff243618), color(0xff313622), color(0xff473216)};

    float[][][] import_paths = new float[][][]{
        branch5, branch4, branch3, branch2, branch1, arm0, arm1, arm3, arm4, arm5, arm6, arm7, arm8
    };

    ArrayList<ArrayList<PVector>> paths = new ArrayList<ArrayList<PVector>>(); 
    for(float[][] p:import_paths){
        paths.add(inkscapePathImport(p, 3564.00000f, 5014.66650f));
    }
    render.beginDraw();
    // render.background(255);
    // for(int i=0; i<paths.size(); i++){
    //     Ribbon r = new Ribbon(paths.get(i), renderHighRes ? printDpi/previewDpi * 50 : 50, true);
    //     r.vadenWeb(300, 10, new Gradient(line_palette), i<=4 ? true : false );
    // }
    Ribbon r = new Ribbon(paths.get(7), renderHighRes ? printDpi/previewDpi * 50 : 50, true);
    
    Polygon p = new Polygon(paths.get(7), true);
    p.subdivide();
    p.geometricSubdivision.display();

    // r.vadenWeb(300, 10, new Gradient(line_palette), false );

    render.endDraw();

    // as = new AttractorSystem(5);
    // as.addPerlinFlowField(0.005, 4, 0.5, true);
    // as.addPerlinFlowField(0.01, 8, 0.9, false);


}


public void draw(){
    render.beginDraw();
    if(firstFrame){
        firstFrame = false;
        render.colorMode(HSB, 360,100,100,100);
        render.fill(0);
    }

    //ANY LOGIC USED TO DRAW GOES HERE
    // as.calculateAttractorSystem();
    // float[][][] import_paths = new float[][][]{
    //     branch5, branch4, branch3, branch2, branch1, arm0, arm1, arm3, arm4, arm5, arm6, arm7, arm8
    // };

    
    // ArrayList<ArrayList<PVector>> paths = new ArrayList<ArrayList<PVector>>(); 
    // for(float[][] p:import_paths){
    //     paths.add(inkscapePathImport(p, 3564.00000, 5014.66650));
    // }
    // render.beginDraw();
    // render.background(255);
    // for(int i=0; i<paths.size(); i++){
    //     Ribbon r = new Ribbon(paths.get(i), renderHighRes ? printDpi/previewDpi * 50 : 50, true);
    //     r.vadenWeb(300, 10, new Gradient(line_palette), i<=4 ? true : false );
    //     // render.noFill();
    //     // r.display();
    // }
    render.endDraw(); //some settings to display the render object on screen
    int outWidth, outHeight;
    
    float ratio = renderWidth / (float)renderHeight;
    if (ratio > 1) {
        outWidth = width;
        outHeight = (int)(outWidth / ratio);
    } else {
        outHeight = height;
        outWidth = (int)(outHeight * ratio);
    }
    
    background(192);
    image(render, (width-outWidth)/2, (height - outHeight) / 2, outWidth, outHeight);

}


public int lerpColor(int[] arr, float step, int colorMode) {
  int sz = arr.length;
  if (sz == 1 || step <= 0.0f) {
    return arr[0];
  } else if (step >= 1.0f) {
    return arr[sz - 1];
  }
  float scl = step * (sz - 1);
  int i = PApplet.parseInt(scl);
  return lerpColor(arr[i], arr[i + 1], scl - i, colorMode);
}



/* ************************* GENERAL UTILITIES *************************************/

public PVector getTorusPosition (PVector position) {
  PVector pos = position.copy();
  if (pos.x < 0) pos.x = renderWidth + pos.x;
  if (pos.x > renderWidth) pos.x %= renderWidth;
  if (pos.y < 0) pos.y = renderHeight + pos.y;
  if (pos.y > renderHeight) pos.y = pos.y %= renderHeight;
  return pos;
}

/* ************************ INTERSECTION AND COLLISION TESTS ******************************/
public boolean lineLine(float x1, float y1, float x2, float y2, float x3, float y3, float x4, float y4) {

  // calculate the direction of the lines
  float uA = ((x4-x3)*(y1-y3) - (y4-y3)*(x1-x3)) / ((y4-y3)*(x2-x1) - (x4-x3)*(y2-y1));
  float uB = ((x2-x1)*(y1-y3) - (y2-y1)*(x1-x3)) / ((y4-y3)*(x2-x1) - (x4-x3)*(y2-y1));

  // if uA and uB are between 0-1, lines are colliding
  if (uA >= 0 && uA <= 1 && uB >= 0 && uB <= 1) {
    return true;
  }
  return false;
}

public boolean polyPoint(ArrayList<PVector> vertices, float px, float py) {
  boolean collision = false;

  // go through each of the vertices, plus
  // the next vertex in the list
  int next = 0;
  for (int current=0; current<vertices.size(); current++) {

    // get next vertex in list
    // if we've hit the end, wrap around to 0
    next = current+1;
    if (next == vertices.size()) next = 0;

    // get the PVectors at our current position
    // this makes our if statement a little cleaner
    PVector vc = vertices.get(current);    // c for "current"
    PVector vn = vertices.get(next);       // n for "next"

    // compare position, flip 'collision' variable
    // back and forth
    if (((vc.y >= py && vn.y < py) || (vc.y < py && vn.y >= py)) &&
         (px < (vn.x-vc.x)*(py-vc.y) / (vn.y-vc.y)+vc.x)) {
            collision = !collision;
    }
  }
  return collision;
}

public boolean polyLine(ArrayList<PVector> vertices, float x1, float y1, float x2, float y2) {

  // go through each of the vertices, plus the next
  // vertex in the list
  int next = 0;
  for (int current=0; current<vertices.size(); current++) {

    // get next vertex in list
    // if we've hit the end, wrap around to 0
    next = current+1;
    if (next == vertices.size()) next = 0;

    // get the PVectors at our current position
    // extract X/Y coordinates from each
    float x3 = vertices.get(current).x;
    float y3 = vertices.get(current).y;
    float x4 = vertices.get(next).x;
    float y4 = vertices.get(next).y;

    // do a Line/Line comparison
    // if true, return 'true' immediately and
    // stop testing (faster)
    boolean hit = lineLine(x1, y1, x2, y2, x3, y3, x4, y4);
    if (hit) {
      return true;
    }
  }

  // never got a hit
  return false;
}

/************************** EASING FUNCTIONS ********************************/

public float sigmoidEasing(float x){
  return x==0 ? 0 : 1/(1+exp(-x));
}

public float easing(float x){
  return sigmoidEasing(x);
}


/************************* POISSON DISK SAMPLING ****************************/

public boolean isValidPoint(PVector[][] grid, float cellsize,
                     int gwidth, int gheight,
                     PVector p, float radius) {
  /* Make sure the point is on the screen */
  if (p.x < 0 || p.x >= renderWidth || p.y < 0 || p.y >= renderHeight)
    return false;

  /* Check neighboring eight cells */
  int xindex = floor(p.x / cellsize);
  int yindex = floor(p.y / cellsize);
  int i0 = max(xindex - 1, 0);
  int i1 = min(xindex + 1, gwidth - 1);
  int j0 = max(yindex - 1, 0);
  int j1 = min(yindex + 1, gheight - 1);

  for (int i = i0; i <= i1; i++)
    for (int j = j0; j <= j1; j++)
      if (grid[i][j] != null)
        if (dist(grid[i][j].x, grid[i][j].y, p.x, p.y) < radius)
          return false;

  /* If we get here, return true */
  return true;
}

public void insertPoint(PVector[][] grid, float cellsize, PVector point) {
  int xindex = floor(point.x / cellsize);
  int yindex = floor(point.y / cellsize);
  grid[xindex][yindex] = point;
}

public ArrayList<PVector> poissonDiskSampling(float radius, int k) {
  int N = 2;
  /* The final set of points to return */
  ArrayList<PVector> points = new ArrayList<PVector>();
  /* The currently "active" set of points */
  ArrayList<PVector> active = new ArrayList<PVector>();
  /* Initial point p0 */
  PVector p0 = new PVector(random(renderWidth), random(renderHeight));
  PVector[][] grid;
  float cellsize = floor(radius/sqrt(N));

  /* Figure out no. of cells in the grid for our canvas */
  int ncells_width = ceil(renderWidth/cellsize) + 1;
  int ncells_height = ceil(renderHeight/cellsize) + 1;

  /* Allocate the grid an initialize all elements to null */
  grid = new PVector[ncells_width][ncells_height];
  for (int i = 0; i < ncells_width; i++)
    for (int j = 0; j < ncells_height; j++)
      grid[i][j] = null;

  insertPoint(grid, cellsize, p0);
  points.add(p0);
  active.add(p0);

  while (active.size() > 0) {
    int random_index = PApplet.parseInt(random(active.size()));
    PVector p = active.get(random_index);
    
    boolean found = false;
    for (int tries = 0; tries < k; tries++) {
      float theta = random(360);
      float new_radius = random(radius, 2*radius);
      float pnewx = p.x + new_radius * cos(radians(theta));
      float pnewy = p.y + new_radius * sin(radians(theta));
      PVector pnew = new PVector(pnewx, pnewy);
      
      if (!isValidPoint(grid, cellsize,
                        ncells_width, ncells_height,
                        pnew, radius))
        continue;
      
      points.add(pnew);
      insertPoint(grid, cellsize, pnew);
      active.add(pnew);
      found = true;
      break;
    }
    
    /* If no point was found after k tries, remove p */
    if (!found)
      active.remove(random_index);
  }

  return points;
}

/* ************************* BACKGROUND TEXTURES AND OVERLAYS *********************************/
// example color palette setup syntax and usage (parameters in example are safe to use in most cases)
// int[] palette = new int[]{color(#ffe4cd), color(#703642), color(#703642), color(#abcdef), color(#4682b4)};
// watercolorBackgroundTexture(new int[]{color(#ffe4cd), color(#703642), color(#703642), color(#abcdef), color(#4682b4)}, 5000, 5, 50, 0.1, 4); //red rocks scheme

// void watercolorBackgroundTexture(int[] baseColors, int numPoints, int numLayers, float geometryWidth, float colorVar, float jitter){
//   //5000 points and 5 layers each is a good balance (with jitter of 4)

//   Gradient grad = new Gradient(baseColors);

//   // render.blendMode(ADD);
//   render.colorMode(HSB, 360, 100,100,100);
//   render.noStroke();
//   for(int i=0; i<numPoints; i++){ //number of locations we drop geometry at
//       render.pushMatrix();
//       float x = random(-renderWidth/4, renderWidth*5/4);
//       float y = random(-renderHeight/4, renderHeight*5/4);
//       render.translate(x,y);
//       int baseColor = grad.eval(map(x,0,renderWidth,0,1)+randomGaussian()*colorVar, HSB);
//       float w = renderHighRes ? geometryWidth*printDpi/previewDpi : geometryWidth;
//       for(int j=0; j<numLayers; j++){ //number of layers we drop for each of the locations
//         render.rotate(map(j, 0, numLayers, 0, TWO_PI));
//         render.fill(hue(baseColor) + randomGaussian(), saturation(baseColor) + randomGaussian()*8, brightness(baseColor) + randomGaussian()*8, 2+randomGaussian());
//         w = w + randomGaussian()*w/4;
//         // render.rect(0,0,w,w);
//         render.beginShape();
//         for(int k=0; k<8; k++){
//           render.curveVertex(((w/2)+randomGaussian()*w/8)*cos(map(k,0,7,0,TWO_PI)), ((w/2)+randomGaussian()*w/8)*sin(map(k,0,7,0,TWO_PI))); 
//           //the parametric equations used (e.g. (r*cos, r*sin) or (r*tan, r*sin) have a large impact on appearance 
//           //(r*tan, r*tan), (r*sin, r*tan) and visa versa has a cool "woven" appearance that renders well in high resolution 
//           //(r*cos, r*cos) has a more rough, geometric mix with hard lines. this also works well in high res
//           //(r*sin, r*sin) is similar to the above but looks more well mixed in high res than (r*cos, r*cos)
//           //(r*cos, r*tan) and visa versa has a blocky, chaotic pattern. does not render will in high res unless you are looking for some more geometric hard lines 
//         }
//         render.endShape(CLOSE);
//       }
//       render.popMatrix();
//     }

//     render.loadPixels();//then we jitter the resulting pixels to add some grainy vibes 
//     for(int i=0; i<render.pixels.length; i++){
//         color c = render.pixels[i];
//         render.pixels[i] = color(hue(c)+random(-jitter, jitter), saturation(c)+randomGaussian()*jitter, brightness(c) + randomGaussian()*jitter/2);
//     }
//     render.updatePixels();
// }


public void watercolorBackgroundTexture(int[] baseColors, ArrayList<PVector> points, int numLayers, float geometryWidth, float colorVar, float jitter){
  //5000 points and 5 layers each is a good balance (with jitter of 4)

  Gradient grad = new Gradient(baseColors);

  // render.blendMode(ADD);
  render.colorMode(HSB, 360, 100,100,100);
  render.noStroke();
  for(int i=0; i<points.size(); i++){ //number of locations we drop geometry at
      render.pushMatrix();
      float x = points.get(i).x;
      float y = points.get(i).y;
      render.translate(x,y);
      int baseColor = grad.eval(map(x,0,renderWidth,0,1)+randomGaussian()*colorVar, HSB);
      float w = renderHighRes ? geometryWidth*printDpi/previewDpi : geometryWidth;
      for(int j=0; j<numLayers; j++){ //number of layers we drop for each of the locations
        render.rotate(map(j, 0, numLayers, 0, TWO_PI));
        render.fill(hue(baseColor) + randomGaussian(), saturation(baseColor) + randomGaussian()*8, brightness(baseColor) + randomGaussian()*8, 2+randomGaussian());
        w = w + randomGaussian()*w/4;
        // render.rect(0,0,w,w);
        render.beginShape();
        for(int k=0; k<8; k++){
          render.curveVertex(((w/2)+randomGaussian()*w/8)*cos(map(k,0,7,0,TWO_PI)), ((w/2)+randomGaussian()*w/8)*sin(map(k,0,7,0,TWO_PI))); 
          //the parametric equations used (e.g. (r*cos, r*sin) or (r*tan, r*sin) have a large impact on appearance 
          //(r*tan, r*tan), (r*sin, r*tan) and visa versa has a cool "woven" appearance that renders well in high resolution 
          //(r*cos, r*cos) has a more rough, geometric mix with hard lines. this also works well in high res
          //(r*sin, r*sin) is similar to the above but looks more well mixed in high res than (r*cos, r*cos)
          //(r*cos, r*tan) and visa versa has a blocky, chaotic pattern. does not render will in high res unless you are looking for some more geometric hard lines 
        }
        render.endShape(CLOSE);
      }
      render.popMatrix();
    }

    // render.loadPixels();//then we jitter the resulting pixels to add some grainy vibes 
    // for(int i=0; i<render.pixels.length; i++){
    //     color c = render.pixels[i];
    //     render.pixels[i] = color(hue(c)+random(-jitter, jitter), saturation(c)+randomGaussian()*jitter, brightness(c) + randomGaussian()*jitter/2);
    // }
    // render.updatePixels();
}

public void restricted_chaikin_overlay(){
    render.smooth();
    render.noFill();

    float x = random(render.width);
    float y = random(render.height);
    float j = 2.5f;
    PShape s = render.createShape();
    float border = 100;
    s.beginShape();
    for (int i = 0; i < 50000; i++) {
        s.vertex(x, y);
        float qx = random(1) < 0.5f ? -1 : 1;
        float qy = random(1) < 0.5f ? -1 : 1;
        x += qx * j;
        y += qy * j;
        x = constrain(x, -border, render.width + border);
        y = constrain(y, -border, render.height + border);
    }
    s.endShape();
    render.shape(chaikin_open(s, 0.25f, 3), 0, 0);
}

public void gridline(float x1, float y1, float x2, float y2) {
  float tmp;
  /* Swap coordinates if needed so that x1 <= x2 */
  if (x1 > x2) { tmp = x1; x1 = x2; x2 = tmp; tmp = y1; y1 = y2; y2 = tmp; }

  float dx = x2 - x1;
  float dy = y2 - y1;
  float step = 1;

  if (x2 < x1)
    step = -step;

  float sx = x1;
  float sy = y1;
  for (float x = x1+step; x <= x2; x+=step) {
    float y = y1 + step * dy * (x - x1) / dx;
    render.strokeWeight(1 + map(noise(sx, sy), 0, 1, -0.5f, 0.5f));
    render.line(sx, sy, x + map(noise(x, y), 0, 1, -1, 1), y + map(noise(x, y), 0, 1, -1, 1));
    sx = x;
    sy = y;
  }
}

public void canvas_overlay_example1() {
  float spacing = renderHighRes ? 5 * printDpi/previewDpi : 5;
  for (int i = -renderWidth; i < renderHeight + renderWidth; i+=spacing)
    gridline(i, 0, i + renderHeight, renderHeight);
  for (int i = renderHeight + renderWidth; i >= -renderWidth; i-=spacing)
    gridline(i, 0, i - renderHeight, renderHeight);
}
// ************************** Attractor System ****************************************

//global variables needed: 
//float pressure_lag=50; 
//float attractorForce=1000; //1000 is appropriate default for the Lorenz attractor family  

class Attractor{
  PVector pos;
  float force;
  ArrayList<Float> pressure;
  int c;
  float target_number;
  float dt;
  String attribute;
  int alpha;
  boolean disp=true;
  int age=0;
  int pressure_lag = 50;
  int attractorForce = 1000;
  boolean normalize_attractor_force = true; //if normalize_attractor_force == true, kill_Particles should also = true
  float damp = 0.0002f; //applied to particle velocity every frame. set to 1 pill not dampen at all


  
  Attractor(float _x, float _y, float _z){
    pos = new PVector(_x, _y, _z);
    c = color(0);
    float initial_pressure = 5000;
    pressure = new ArrayList<Float>();
    for(int i=0; i<pressure_lag; i++){pressure.add(initial_pressure);}
    force = attractorForce;
    dt = 0.001f;
    }
  

  Attractor(){
      c = color(0,0,100);
  }
  
  public void calc_force_on(Particle p){
    PVector f = PVector.sub(pos, p.pos);
    float d = f.mag();
    p.vel.add( new PVector(0,0));
    
  }
    
  public void update_force(){
    force = (347.3058f + (6147.582f - 347.3058f)/pow(1 + (pressure.get(0)/7004.265f),8.499077f)) + noise(5000,10000)*random(-1,1);
  }

  public void display(){age++;}
  
}

class LorenzAttractor extends Attractor{
  
  LorenzAttractor(float _x, float _y, float _z){
    super(_x, _y, _z);
  }
  
  public void calc_force_on(Particle p){
    PVector f = PVector.sub(pos, p.pos);
    float d = f.mag();
    if(normalize_attractor_force){
      f.normalize();
      float dx = (10*(f.x + f.y))*dt;
      float dy = (-f.x*f.z + 28*f.x-f.y)*dt;
      float dz = (f.x*f.y - 8.0f/3*f.z)*dt;
      // p.vel.add( new PVector(dx, dy, dz).mult(force/(d+0.0001)));
      // float factor = map(d, 0, sqrt(width*height), 0, 1);
      // p.c = color(c);
      // if(random(1)<0.0002){p.DIVISION_CHANCES+=0.1;}
      // if(random(1)<0.00005){p.DEPOSIT_RATE*=max(randomGaussian()+10, 0);}

    }
    else{
      float dx = (10*(f.x + f.y))*dt;
      float dy = (-f.x*f.z + 28*f.x-f.y)*dt;
      float dz = (f.x*f.y - 8.0f/3*f.z)*dt;
      // p.vel.add( new PVector(dx, dy, dz).normalize().mult(force/(d+0.0001)));
    }
  }
}

class PerlinFlowField extends Attractor{
    PVector[][] grid;

    PerlinFlowField(float noiseScale, int noiseOctaves, float noiseGain, boolean scaleNoiseScale){
        this.grid = new PVector[renderWidth][renderHeight];
        noiseDetail(noiseOctaves, noiseGain);
        float xoff = 0;
        for(int x=0; x<renderWidth; x++){
            float yoff = 0;
            for(int y=0; y<renderHeight; y++){
                float n = noise(xoff,yoff);
                //playing around with clmaping the angles produces large effects  
                float angleClamp = map(easing(map(x,0,renderWidth,0,1)),easing(0.f),easing(1.f),1, PI/2); // (PI/10); 
                // float angleClamp = map(x, 0, renderWidth, 1, PI/2);
                // float angle = floor(map(n, 0 , 1, -PI, PI)/angleClamp)*angleClamp; //(randomGaussian()*angleClamp/10+angleClamp);
                float angle = map(n, 0 , 1, -PI, PI); //for a 'normal' flow field
                grid[x][y] = new PVector(angle, n);
                yoff += (renderHighRes && scaleNoiseScale) ? noiseScale/(printDpi/previewDpi) : noiseScale;;
            }
            xoff += (renderHighRes && scaleNoiseScale) ? noiseScale/(printDpi/previewDpi) : noiseScale; //noise works best with step size of 0.005
        } 
    }

    public void calc_force_on(Particle p){
        if(random(1)<0.0005f){p.DEPOSIT_RATE*=max(randomGaussian()+5, 0);}
        p.vel.add( new PVector(cos(grid[constrain(ceil(p.pos.x),0,renderWidth-1)][constrain(ceil(p.pos.y),0,renderHeight-1)].x),sin(grid[constrain(ceil(p.pos.x),0,renderWidth-1)][constrain(ceil(p.pos.y),0,renderHeight-1)].x)).normalize());

    }

    public void display(){
      
    }
}

class AttractorSystem{
  ArrayList<Particle> particles;
  ArrayList<Attractor> attractors;
  int num_attractors = 15;
  int NB_INITIAL_WALKERS = 100;
  boolean scaleNoiseScale = true; //used when scaling pieces with noise into high res
  int attractorForce = 1000; //base force for the attractors
  boolean kill_Particles = true; //if true particles are removed when their age reaches 0
  float damp = 0.0002f; //applied to particle velocity every frame. set to 1 pill not dampen at all
  boolean normalize_attractor_force = true; //if normalize_attractor_force == true, kill_Particles should also = true
  boolean allow_Particle_overlap=false; //will determine if particles die phen they run into each other 
  float STEP_SIZE = 1;
  boolean DISCRETE_DIV_ANGLE = false;
  //noise parameters 
  // int noiseOctaves = 8//defines numbers of octaves. takes int value from 0-8. higher = more high-level structure
  // float noiseGain = 0.5; //defines how much of each octave carries over into the next (0,1). Higher = more small details
  // float noiseScale = 0.005; // 0.005 is 'best' but 0.1 or slightly higher can have more sweeping patterns 
  //initial particle settings
  float initial_TURN_CHANCES = 0.f;
  float initial_TURN_ANGLE = PI/8;
  float initial_DEPOSIT_RATE = 0.01f;
  float  initial_DIVISION_CHANCES = 0.00f;
  float initial_DIVISION_ANGLE = PI / 8;
  float initial_TERMINATION_THRESHOLD = 0.7f;
  float initial_TERMINATION_CHANCES = initial_DIVISION_CHANCES * 0.f;
  boolean initial_DISCRETE_DIV_ANGLE = false;
  float grid[];


  AttractorSystem(){
    particles = new ArrayList<Particle>();
    grid = new float[renderWidth*renderHeight];
    NB_INITIAL_WALKERS *= renderHighRes ? printDpi/previewDpi : 1;
    for (int i = 0; i < NB_INITIAL_WALKERS; i++) {

        float ang = random(0, TWO_PI);
        float x = random(renderWidth);
        float y = random(renderHeight);
        particles.add(
        new Particle(this, 
          new PVector(x, y), 
          ang,
          initial_TURN_CHANCES,
          initial_TURN_ANGLE,
          initial_DEPOSIT_RATE,
          initial_DIVISION_CHANCES,
          initial_DIVISION_ANGLE,
          initial_TERMINATION_THRESHOLD,
          initial_TERMINATION_CHANCES,
          damp
          )
        );
    }
    attractors = new ArrayList<Attractor>();
    // attractors.add(new PerlinFlowField(noiseScale, noiseOctaves, scaleNoiseScale)); 
  }

  AttractorSystem(int n){
    particles = new ArrayList<Particle>();
    grid = new float[renderWidth*renderHeight];
    NB_INITIAL_WALKERS = n;
    NB_INITIAL_WALKERS *= renderHighRes ? printDpi/previewDpi : 1;
    for (int i = 0; i < NB_INITIAL_WALKERS; i++) {

        float ang = random(0, TWO_PI);
        float x = random(renderWidth);
        float y = random(renderHeight);
        particles.add(
        new Particle(this, 
          new PVector(x, y), 
          ang,
          initial_TURN_CHANCES,
          initial_TURN_ANGLE,
          initial_DEPOSIT_RATE,
          initial_DIVISION_CHANCES,
          initial_DIVISION_ANGLE,
          initial_TERMINATION_THRESHOLD,
          initial_TERMINATION_CHANCES,
          damp
          )
        );
    }
    attractors = new ArrayList<Attractor>();
    // attractors.add(new PerlinFlowField(noiseScale, noiseOctaves, scaleNoiseScale)); 
  }

  public void addPerlinFlowField(float _noiseScale, int _noiseOctaves, float _noiseGain, boolean _scaleNoiseScale){
    attractors.add(new PerlinFlowField(_noiseScale, _noiseOctaves, _noiseGain, _scaleNoiseScale));
  }

  public void calculateAttractorSystem(){//logic that updates particles in the system
    ArrayList<Particle> newParticles = new ArrayList<Particle>();

    for (Particle p : particles) {
        // 1. palking step
        //update logic for mapping all attractor forces, here simple addition
        if(attractors.size() > 0){
            for(int j=0; j<attractors.size(); j++){
                attractors.get(j).calc_force_on(p);
            }         
            p.update();

        }
        else{p.update();}
        
        // 2. division step
        float r = random(0, 1);
        if (r < p.DIVISION_CHANCES) {
        float nAngle = p.ang + (DISCRETE_DIV_ANGLE ? round(random(0, 1))*2-1 : random(-1, 1)) * p.DIVISION_ANGLE;
        Particle nParticle = new Particle(new PVector(p.pos.x, p.pos.y), nAngle, p.boxBounds);
        newParticles.add(nParticle);
        }

        p.displayPath(PApplet.parseInt(randomGaussian()*50+200), randomGaussian()*2+8);
    }

    // adds the new particles to the active list
    for (Particle p : newParticles) {
        particles.add(p);
    }

    // checks for dead particles
    if (allow_Particle_overlap==false || attractors.size()==0){
        for (int i = particles.size()-1; i >= 0; i--) {
        Particle p = particles.get(i);
        float r = random(0, 1);
        if (r < particles.get(i).TERMINATION_CHANCES || (particles.get(i).age <= 0 && kill_Particles)) {
            particles.remove(i);
            continue;
        }
        // turn the particle coordinates into an index to sample the environment color
        // to do that compute the "next" particle position
        PVector dir = new PVector(cos(p.ang), sin(p.ang));
        PVector npos = p.pos.copy().add(p.vel.copy().normalize().mult(2*STEP_SIZE));
        npos = getTorusPosition(npos);
        // sample aggregate to check for collision
        int idx = ceil(npos.x) + ceil(npos.y) * renderWidth; 
        int idxLast = PApplet.parseInt(p.pos.x) + PApplet.parseInt(p.pos.y) * renderWidth; //check to ensure particle moved to a new grid location on this step
        if (idx >= renderWidth*renderHeight || idx==idxLast) continue;
        float aggregate = grid[idx];
        // kill the particle if it will run on some aggregate
        if (aggregate >= 1*p.TERMINATION_THRESHOLD) {
            particles.remove(i);
        }
        }
    }
    else{    //just apply termination checks  
        for (int i = particles.size()-1; i >= 0; i--) {
        Particle p = particles.get(i);
        float r = random(0, 1);
        if (r < particles.get(i).TERMINATION_CHANCES || (particles.get(i).age <= 0 && kill_Particles)) {
            particles.remove(i);
            continue;
        }
        }
    }
  }
}

class Particle {
  PVector pos;
  PVector lastPos;
  float ang;
  PVector vel;
  float age;
  float mass;
  int c;
  float TURN_CHANCES;
  float TURN_ANGLE;
  float DEPOSIT_RATE;
  float DIVISION_CHANCES;
  float DIVISION_ANGLE;
  float TERMINATION_THRESHOLD;
  float TERMINATION_CHANCES;
  boolean initial_DISCRETE_DIV_ANGLE = false;
  boolean allowedinSquare = true;
  float[] boxBounds;
  Gradient xGrad =  new Gradient(line_palette);
  Gradient yGrad =  new Gradient(background_palette);
  ArrayList<PVector> path;
  float STEP_SIZE = 1;
  float damp=0.0002f;
  AttractorSystem attractorSystem;

  public Particle (AttractorSystem a, PVector p, float ang, float initial_TURN_CHANCES,float initial_TURN_ANGLE,float initial_DEPOSIT_RATE,float initial_DIVISION_CHANCES,float initial_DIVISION_ANGLE,float initial_TERMINATION_THRESHOLD,float initial_TERMINATION_CHANCES, float damp) 
    {
    attractorSystem = a;
    this.lastPos = p;
    this.pos = p;
    this.ang = ang;
    vel = new PVector(cos(ang),sin(ang),0).normalize();
    mass = 1;
    TURN_CHANCES = initial_TURN_CHANCES;
    TURN_ANGLE = initial_TURN_ANGLE;
    DEPOSIT_RATE = initial_DEPOSIT_RATE;
    DIVISION_CHANCES = initial_DIVISION_CHANCES;
    DIVISION_ANGLE = initial_DIVISION_ANGLE;
    TERMINATION_THRESHOLD = initial_TERMINATION_THRESHOLD;
    TERMINATION_CHANCES = initial_TERMINATION_CHANCES;
    damp = damp;
    age = randomGaussian()*(renderHeight*0.26667f)/2+(renderHeight*0.26667f);
    path = new ArrayList();
  }

  public Particle (PVector p, float ang) {
    this.lastPos = p;
    this.pos = p;
    this.ang = ang;
    vel = new PVector(cos(ang),sin(ang),0).normalize();
    mass = 1;
    // c = pallatte[int(random(pallatte.length))];
    // c = img.get(int(map(pos.x, 0, renderWidth, 0, img.width)), int(map(pos.y, 0, renderHeight, 0, img.height)));
    // c = img.get(int(map(random(renderWidth), 0, renderWidth, 0, img.width)), int(map(random(renderHeight), 0, renderHeight, 0, img.height)));
    if(abs(vel.x) > abs(vel.y)){
      // bright = randomGaussian()*5+45;
      int baseColor = xGrad.eval(map(pos.x,0,renderWidth,0,1)+randomGaussian()*0.1f, HSB);
      c = color(hue(baseColor) + randomGaussian(), saturation(baseColor) + randomGaussian()*8, brightness(baseColor) + randomGaussian()*8);
      }
    else{
      // bright = randomGaussian()*2+20;
      int baseColor = yGrad.eval(map(pos.y,0,renderHeight,0,1)+randomGaussian()*0.1f, HSB);
      c = color(hue(baseColor) + randomGaussian(), saturation(baseColor) + randomGaussian()*8, brightness(baseColor) + randomGaussian()*8);
    }
    age = randomGaussian()*(renderHeight*0.26667f)/2+(renderHeight*0.26667f);
    path = new ArrayList();
  }

  //bounds array is of form [x1, y1, x2, y2]
  public Particle (PVector p, float ang, float[] bounds) {
    this.lastPos = p;
    this.pos = p;
    this.ang = ang;
    vel = new PVector(cos(ang),sin(ang),0).normalize();
    mass = 1;
    // c = color(map(pos.y, 0, renderHeight, 360%240, 360%200),map(pos.x, 0, renderWidth, 70, 100), 100);
    // c = pallatte[int(random(pallatte.length))];

    // c = img.get(int(map(pos.x, 0, renderWidth, 0, img.width)), int(map(pos.y, 0, renderHeight, 0, img.height)));
    if(abs(vel.x) > abs(vel.y)){
      // bright = randomGaussian()*5+45;
      int baseColor = xGrad.eval(map(pos.x,0,renderWidth,0,1)+randomGaussian()*0.1f, HSB);
      c = color(hue(baseColor) + randomGaussian(), saturation(baseColor) + randomGaussian()*8, brightness(baseColor) + randomGaussian()*8);
      }
    else{
      // bright = randomGaussian()*2+20;
      int baseColor = yGrad.eval(map(pos.y,0,renderHeight,0,1)+randomGaussian()*0.1f, HSB);
      c = color(hue(baseColor) + randomGaussian(), saturation(baseColor) + randomGaussian()*8, brightness(baseColor) + randomGaussian()*8);
    }
    // c = color(0, 0, 100);
    age = randomGaussian()*(renderHeight*0.26667f)/2+(renderHeight*0.26667f);
    boxBounds = bounds;
  }
  
  public void update () {
    lastPos = pos.copy();
    // vel = PVector.sub(pos, lastPos);
    // the Particle has a random chances to turn
    if (random(0, 1) < TURN_CHANCES) {
      this.ang+= TURN_ANGLE * (round(random(0, 1)) * 2.f - 1.f);
      vel.add(new PVector(cos(ang), sin(ang)));

    }

    // move along the direction
    pos.add(vel.copy().normalize().mult(STEP_SIZE));
    vel.mult(damp);
    // makes sure that the Particles stays within the window area
    // pos = getTorusPosition(pos);
    path.add(lastPos);
    age--;
  }
  
  public void display () {
    render.colorMode(HSB,360,100,100,100);
    float bright;

    
    // render.stroke(0, int(100. * DEPOSIT_RATE));

    // if(random(1)>map(PVector.sub(new PVector(renderWidth,0),pos).mag(),0.,sqrt(pow(renderWidth*.75,2) + pow(renderHeight*.75,2)), 0,1))
    // {
    // // render.stroke(randomGaussian()*5+45, map(exponentialEasing(map(pos.x,0,renderWidth,0,1)),exponentialEasing(0.),exponentialEasing(1.),0,75)+randomGaussian()*5, map(exponentialEasing(map(pos.y,0,renderHeight,0,1)),exponentialEasing(0.),exponentialEasing(1.),80,100)+randomGaussian()*2, int(100*DEPOSIT_RATE));
    // render.stroke(30,70,bright, int(100*DEPOSIT_RATE)+int(randomGaussian()*2+2));

    // }
    // else if(random(1)<0.005){
    //   // render.stroke(30,bright,40, int(100*DEPOSIT_RATE)+ random(1)<0.0001 ? 10 : 0);
    // }
    // else{
    //   render.stroke(0,0,70, int(100*DEPOSIT_RATE));
    //   }
    render.stroke(0, PApplet.parseInt(100*DEPOSIT_RATE));    
    PVector line = lastPos.copy().sub(pos);
    if (line.mag() < 4*STEP_SIZE) {
      render.line(lastPos.x, lastPos.y, pos.x, pos.y);
      int idx = ceil(pos.x) + ceil(pos.y) * renderWidth;
      if (idx < renderWidth*renderHeight){attractorSystem.grid[idx] += DEPOSIT_RATE;} //update substrate grid when drawing Particle 
    }
  }

  public void displayPath(int path_length, float stroke_weight){
    if(path.size()==path_length){
      render.fill(0,20);
      Ribbon r = new Ribbon(path, stroke_weight, false);
      r.display();
    } 
  }
}


/* *********************** HIGH RES EXPORT FUNCTIONALITY **************************/

//INCLUDE THESE GLOBAL VARIABLES  
// PGraphics render;
// PImage img;
// String saveFilePath = "../outputs/XXXXXXXXX_output-" + new java.text.SimpleDateFormat("yyyyMMdd-HHmmss").format(new java.util.Date()); //change the XXXX to current project name 
// int printWidth = 6;
// int printHeight = 6;
// int printDpi = 300;
// int previewDpi = 72;
// boolean renderHighRes = false;
// boolean firstFrame = true;
// int renderWidth;
// int renderHeight;
// float scaleFactor = 1;
// int seed = 0;

// void setup(){
//     size(750, 750);
//     background(255);
//     colorMode(HSB, 360, 100, 100, 100);
//     doReset();
// }


// void doReset() { //initial setup and used in resetting for high-def export

//     int dpi = renderHighRes ? printDpi : previewDpi;
//     scaleFactor = dpi / (float)previewDpi;
//     renderWidth = printWidth * dpi;
//     renderHeight = printHeight * dpi;
//     render = createGraphics(renderWidth, renderHeight);
//     firstFrame = true;
//     noiseSeed(seed);
//     randomSeed(seed);

//     render.beginDraw();
//     render.colorMode(HSB, 360, 100,100,100);
//     // ANY DRAWING LOGIC TO EXECUTE ON STARTUP/RESET GOES HERE
//     render.endDraw();

// }

// void draw(){
//     render.beginDraw();
//     if(firstFrame){
//         firstFrame = false;
//         render.colorMode(HSB, 360,100,100,100);
//     }
    
//     //ANY LOGIC USED TO DRAW GOES HERE
    
//     render.endDraw(); //some settings to display the render object on screen
//     int outWidth, outHeight;
    
//     float ratio = renderWidth / (float)renderHeight;
//     if (ratio > 1) {
//         outWidth = width;
//         outHeight = (int)(outWidth / ratio);
//     } else {
//         outHeight = height;
//         outWidth = (int)(outHeight * ratio);
//     }
    
//     background(192);
//     image(render, (width-outWidth)/2, (height - outHeight) / 2, outWidth, outHeight);

// }

public void keyPressed() {
  switch (key) {
    case 's':
      render.save(saveFilePath + "-" + "SEED-" + str(seed) + ".png");
      break;
      
    case 'r':
      seed = (int)System.currentTimeMillis();
      renderHighRes = false;
      doReset();
      break;

    case 'R':
      seed = (int)System.currentTimeMillis();
      renderHighRes = true;
      doReset();
      break;
      
    case 'h':
      renderHighRes = true;
      doReset();
      break;

    case 'H':
      renderHighRes = true;
      doReset();
      break;
      
  }
}

/* ********************** CHAIKIN CURVE FUNCTIONS **************************** */

public void chaikin_line(ArrayList<PVector> points, int layers, float layer_var, String mode){
  PShape n = render.createShape();
  n.beginShape();
  for(int i=0; i<points.size(); i++){
    n.vertex(points.get(i).x, points.get(i).y);
  }
  n.endShape();
  
  if(mode == "CLOSE"){
    render.shape(chaikin_close(n, 0.25f, 5), 0, 0);
  
    for(int i=0; i<=layers; i++){
      render.shape(chaikin_close(perturbPShape(n, layer_var, layer_var), .25f, 5),0,0);
    }
  }
  else{
    render.shape(chaikin_open(n, 0.25f, 5), 0, 0);
  
    for(int i=0; i<=layers; i++){
      render.shape(chaikin_open(perturbPShape(n, layer_var, layer_var), .25f, 5),0,0);
    }
  }
}

public PShape perturbPShape(PShape s, float x_std, float y_std){
  PShape n = render.createShape();
  n.beginShape();
  for(int i=0; i<s.getVertexCount(); i++){
    n.vertex(s.getVertex(i).x+randomGaussian()*x_std, s.getVertex(i).y+randomGaussian()*y_std);
  }
  n.endShape();
  return n;
}

public ArrayList<PVector> chaikin_cut(PVector a, PVector b, float ratio) {
  float x, y;
  ArrayList<PVector> n = new ArrayList<PVector>();

  /*
   * If ratio is greater than 0.5 flip it so we avoid cutting across
   * the midpoint of the line.
   */
   if (ratio > 0.5f) ratio = 1 - ratio;

  /* Find point at a given ratio going from A to B */
  x = lerp(a.x, b.x, ratio);
  y = lerp(a.y, b.y, ratio);
  n.add(new PVector(x, y));

  /* Find point at a given ratio going from B to A */
  x = lerp(b.x, a.x, ratio);
  y = lerp(b.y, a.y, ratio);
  n.add(new PVector(x, y));

  return n;
}

public PShape chaikin(PShape shape, float ratio, int iterations, boolean close) {
  // If the number of iterations is zero, return shape as is
  if (iterations == 0)
    return shape;

  PShape next = render.createShape();
  next.beginShape();

  /*
   * Step 1: Figure out how many corners the shape has
   *         depending on whether it's open or closed.
   */
  int num_corners = shape.getVertexCount();
  if (!close)
    num_corners = shape.getVertexCount() - 1;

  /*
   * Step 2: Since we don't have access to edges directly
   *         with a PShape object, do a pairwise iteration
   *         over vertices instead. Same thing.
   */
  for (int i = 0; i < num_corners; i++) {

    // Get the i'th and (i+1)'th vertex to work on that edge.
    PVector a = shape.getVertex(i);
    PVector b = shape.getVertex((i + 1) % shape.getVertexCount());

    // Step 3: Break it using our chaikin_break() function
    ArrayList<PVector> n = chaikin_cut(a, b, ratio);

    /*
     * Now we have to deal with one corner case. In the case
     * of open shapes, the first and last endpoints shouldn't
     * be moved. However, in the case of closed shapes, we
     * cut all edges on both ends.
     */
    if (!close && i == 0) {
      // For the first point of open shapes, ignore vertex A
      next.vertex(a.x, a.y);
      next.vertex(n.get(1).x, n.get(1).y);
    } else if (!close && i == num_corners - 1) {
      // For the last point of open shapes, ignore vertex B
      next.vertex(n.get(0).x, n.get(0).y);
      next.vertex(b.x, b.y);
    } else {
      // For all other cases (i.e. interior edges of open
      // shapes or edges of closed shapes), add both vertices
      // returned by our chaikin_break() method
      next.vertex(n.get(0).x, n.get(0).y);
      next.vertex(n.get(1).x, n.get(1).y);
    }
  }

  if (close)
    next.endShape(CLOSE);
  else
    next.endShape();

  return chaikin(next, ratio, iterations - 1, close);
}

public PShape chaikin_close(PShape original, float ratio, int iterations) {
  return chaikin(original, ratio, iterations, true);
}

public PShape chaikin_open(PShape original, float ratio, int iterations) {
  return chaikin(original, ratio, iterations, false);
}


/*  **************** COLORING ALGORITHMS AND METHODS *****************************/


/*  Example usage of colorstops and gradient class: 
    
    int clrCount = 3;
    Gradient grd;
    ColorStop[] temp = new ColorStop[clrCount];
    float prc;
    float hue = random(1);
    for (int i = 0; i < clrCount; ++i) {
      prc = i == 0 ? 0 : i == clrCount - 1 ? 1 : random(1);
      temp[i] = new ColorStop(HSB, prc, new float[]{hue, random(1), random(1), 1});
    }
    grd = new Gradient(temp);

    for(int i=0; i<renderWidth; i++){
      render.stroke(grd.eval(map(i,0,renderWidth, 0, 1), HSB));
      render.line(i,0,i,renderHeight);
    }
    render.endDraw();
 */



class ColorStop implements Comparable<ColorStop> {
  static final float TOLERANCE = 0.09f;
  float percent; //this is just used for sorting in the gradient class, could be thought of as a heirarchy with arbitrary scale
  int clr;

  ColorStop(int colorMode, float percent, float[] arr) {
    this(colorMode, percent, arr[0], arr[1], arr[2],
      arr.length == 4 ? arr[3] : 1.0f);
  }

  ColorStop(int colorMode, float percent, float x, float y, float z, float w) {
    this(percent, colorMode == HSB ? composeclr(hsbToRgb(x, y, z, w))
      : composeclr(x, y, z, w));
  }

  ColorStop(float percent, int clr) {
    this.percent = constrain(percent, 0.0f, 1.0f);
    this.clr = clr;
  }

  public boolean approxPercent(ColorStop cs, float tolerance) {
    return abs(percent - cs.percent) < tolerance;
  }

  // Mandated by the interface Comparable<ColorStop>.
  // Permits color stops to be sorted by Collections.sort.
  public int compareTo(ColorStop cs) {
    return percent > cs.percent ? 1 : percent < cs.percent ? -1 : 0;
  }
}

class Gradient {
  static final int DEFAULT_COLOR_MODE = RGB;
  ArrayList<ColorStop> colorStops = new ArrayList<ColorStop>();

  Gradient() {
    this(0xff000000, 0xffffffff);
  }

  // Creates equidistant color stops.
  Gradient(int... colors) {
    int sz = colors.length;
    float szf = sz <= 1.0f ? 1.0f : sz - 1.0f;
    for (int i = 0; i < sz; ++i) {
      colorStops.add(new ColorStop(i / szf, colors[i]));
    }
  }

  // Creates equidistant color stops.
  Gradient(int colorMode, float[]... colors) {
    int sz = colors.length;
    float szf = sz <= 1.0f ? 1.0f : sz - 1.0f;
    for (int i = 0; i < sz; ++i) {
      colorStops.add(new ColorStop(colorMode, i / szf, colors[i]));
    }
  }

  Gradient(ColorStop... colorStops) {
    int sz = colorStops.length;
    for (int i = 0; i < sz; ++i) {
      this.colorStops.add(colorStops[i]);
    }
    java.util.Collections.sort(this.colorStops);
    remove();
  }

  Gradient(ArrayList<ColorStop> colorStops) {
    this.colorStops = colorStops;
    java.util.Collections.sort(this.colorStops);
    remove();
  }

  public void add(int colorMode, float percent, float[] arr) {
    add(new ColorStop(colorMode, percent, arr));
  }

  public void add(int colorMode, float percent,
    float x, float y, float z, float w) {
    add(new ColorStop(colorMode, percent, x, y, z, w));
  }

  public void add(final float percent, final int clr) {
    add(new ColorStop(percent, clr));
  }

  public void add(final ColorStop colorStop) {
    for (int sz = colorStops.size(), i = sz - 1; i > 0; --i) {
      ColorStop current = colorStops.get(i);
      if (current.approxPercent(colorStop, ColorStop.TOLERANCE)) {
        println(current, "will be replaced by", colorStop);
        colorStops.remove(current);
      }
    }
    colorStops.add(colorStop);
    java.util.Collections.sort(colorStops);
  }

  public int eval(final float step) {
    return eval(step, DEFAULT_COLOR_MODE);
  }

  public int eval(final float step, final int colorMode) {
    int sz = colorStops.size();

    // Exit from the function early whenever possible.
    if (sz == 0) {
      return 0x00000000;
    } else if (sz == 1 || step < 0.0f) {
      return colorStops.get(0).clr;
    } else if (step >= 1.0f) {
      return colorStops.get(sz - 1).clr;
    }

    ColorStop currStop;
    ColorStop prevStop;
    float currPercent, scaledst;
    for (int i = 0; i < sz; ++i) {
      currStop = colorStops.get(i);
      currPercent = currStop.percent;

      if (step < currPercent) {

        // These can be declared within the for-loop because
        // if step < currPercent, the function will return
        // and no more iterations will be executed.
        float[] originclr = new float[4];
        float[] destclr = new float[4];
        float[] rsltclr = new float[4];

        // If not at the first stop in the gradient (i == 0),
        // then get the previous.
        prevStop = colorStops.get(i - 1 < 0 ? 0 : i - 1);

        scaledst = step - currPercent;
        float denom = prevStop.percent - currPercent;
        if (denom != 0) {
          scaledst /= denom;
        }

        // Assumes that color stops' colors are ints. They could
        // also be float[] arrays, in which case they wouldn't
        // need to be decomposed.
        switch(colorMode) {
        case HSB:
          rgbToHsb(currStop.clr, originclr);
          rgbToHsb(prevStop.clr, destclr);
          smootherStepHsb(originclr, destclr, scaledst, rsltclr);
          return composeclr(hsbToRgb(rsltclr));
        case RGB:
          decomposeclr(currStop.clr, originclr);
          decomposeclr(prevStop.clr, destclr);
          smootherStepRgb(originclr, destclr, scaledst, rsltclr);
          return composeclr(rsltclr);
        }
      }
    }
    return colorStops.get(sz - 1).clr;
  }

  public boolean remove(ColorStop colorStop) {
    return colorStops.remove(colorStop);
  }

  public ColorStop remove(int i) {
    return colorStops.remove(i);
  }

  public int remove() {
    int removed = 0;
    for (int sz = colorStops.size(), i = sz - 1; i > 0; --i) {
      ColorStop current = colorStops.get(i);
      ColorStop prev = colorStops.get(i - 1);
      if (current.approxPercent(prev, ColorStop.TOLERANCE)) {
        println(current, "removed, as it was too close to", prev);
        colorStops.remove(current);
        removed++;
      }
    }
    return removed;
  }
}

//ken perlin's smoother step functions
public float[] smootherStepRgb(float[][] arr, float st, float[] out) {
  int sz = arr.length;
  if (sz == 1 || st < 0) {
    out = java.util.Arrays.copyOf(arr[0], 0);
    return out;
  } else if (st > 1) {
    out = java.util.Arrays.copyOf(arr[sz - 1], 0);
    return out;
  }
  float scl = st * (sz - 1);
  int i = PApplet.parseInt(scl);
  float eval = smootherStep(scl - i);
  out[0] = arr[i][0] + eval * (arr[i + 1][0] - arr[i][0]);
  out[1] = arr[i][1] + eval * (arr[i + 1][1] - arr[i][1]);
  out[2] = arr[i][2] + eval * (arr[i + 1][2] - arr[i][2]);
  out[3] = arr[i][3] + eval * (arr[i + 1][3] - arr[i][3]);
  return out;
}

public float[] smootherStepHsb(float[] a, float[] b, float st, float[] out) {

  // Find difference in hues.
  float huea = a[0];
  float hueb = b[0];
  float delta = hueb - huea;

  // Prefer shortest distance.
  if (delta < -0.5f) {
    hueb += 1.0f;
  } else if (delta > 0.5f) {
    huea += 1.0f;
  }

  float eval = smootherStep(st);

  // The two hues may be outside of 0 .. 1 range,
  // so modulate by 1.
  out[0] = (huea + eval * (hueb - huea)) % 1;
  out[1] = a[1] + eval * (b[1] - a[1]);
  out[2] = a[2] + eval * (b[2] - a[2]);
  out[3] = a[3] + eval * (b[3] - a[3]);
  return out;
}

public float smootherStep(float st) {
  return st * st * st * (st * (st * 6.0f - 15.0f) + 10.0f);
}

public float[] smootherStepRgb(float[] a, float[] b, float st, float[] out) {
  float eval = smootherStep(st);
  out[0] = a[0] + eval * (b[0] - a[0]);
  out[1] = a[1] + eval * (b[1] - a[1]);
  out[2] = a[2] + eval * (b[2] - a[2]);
  out[3] = a[3] + eval * (b[3] - a[3]);
  return out;
}

//general utlities to work with bit representations of color 
public int composeclr(float[] in) {
  return composeclr(in[0], in[1], in[2], in[3]);
}

// Assumes that RGBA are in range 0 .. 1.
public int composeclr(float red, float green, float blue, float alpha) {
  return round(alpha * 255.0f) << 24
    | round(red * 255.0f) << 16
    | round(green * 255.0f) << 8
    | round(blue * 255.0f);
}

public float[] decomposeclr(int clr) {
  return decomposeclr(clr, new float[] { 0.0f, 0.0f, 0.0f, 1.0f });
}

// Assumes that out has 4 elements.
// 1.0 / 255.0 = 0.003921569
public float[] decomposeclr(int clr, float[] out) {
  out[3] = (clr >> 24 & 0xff) * 0.003921569f;
  out[0] = (clr >> 16 & 0xff) * 0.003921569f;
  out[1] = (clr >> 8 & 0xff) * 0.003921569f;
  out[2] = (clr & 0xff) * 0.003921569f;
  return out;
}


//HSB to RBG fxns
public float[] hsbToRgb(float[] in) {
  float[] out = new float[] { 0.0f, 0.0f, 0.0f, 1.0f };
  return hsbToRgb(in[0], in[1], in[2], in[3], out);
}

public float[] hsbToRgb(float[] in, float[] out) {
  if (in.length == 3) {
    return hsbToRgb(in[0], in[1], in[2], 1.0f, out);
  } else if (in.length == 4) {
    return hsbToRgb(in[0], in[1], in[2], in[3], out);
  }
  return out;
}

public float[] hsbToRgb(float hue, float sat, float bri, float alpha) {
  float[] out = new float[] { 0.0f, 0.0f, 0.0f, 1.0f };
  return hsbToRgb(hue, sat, bri, alpha, out);
}

public float[] hsbToRgb(float hue, float sat, float bri, float alpha, float[] out) {
  if (sat == 0.0f) {

    // 0.0 saturation is grayscale, so all values are equal.
    out[0] = out[1] = out[2] = bri;
  } else {

    // Divide color wheel into 6 sectors.
    // Scale up hue to 6, convert to sector index.
    float h = hue * 6.0f;
    int sector = PApplet.parseInt(h);

    // Depending on the sector, three tints will
    // be distributed among R, G, B channels.
    float tint1 = bri * (1.0f - sat);
    float tint2 = bri * (1.0f - sat * (h - sector));
    float tint3 = bri * (1.0f - sat * (1.0f + sector - h));

    switch (sector) {
    case 1:
      out[0] = tint2; out[1] = bri; out[2] = tint1;
      break;
    case 2:
      out[0] = tint1; out[1] = bri; out[2] = tint3;
      break;
    case 3:
      out[0] = tint1; out[1] = tint2; out[2] = bri;
      break;
    case 4:
      out[0] = tint3; out[1] = tint1; out[2] = bri;
      break;
    case 5:
      out[0] = bri; out[1] = tint1; out[2] = tint2;
      break;
    default:
      out[0] = bri; out[1] = tint3; out[2] = tint1;
    }
  }

  out[3] = alpha;
  return out;
}

//RGB to HSB fxns
public float[] rgbToHsb(int clr) {
  return rgbToHsb(clr, new float[] { 0.0f, 0.0f, 0.0f, 1.0f });
}

public float[] rgbToHsb(int clr, float[] out) {
  return rgbToHsb((clr >> 16 & 0xff) * 0.003921569f,
    (clr >> 8 & 0xff) * 0.003921569f,
    (clr & 0xff) * 0.003921569f,
    (clr >> 24 & 0xff) * 0.003921569f, out);
}

public float[] rgbToHsb(float[] in) {
  return rgbToHsb(in, new float[] { 0.0f, 0.0f, 0.0f, 1.0f });
}

public float[] rgbToHsb(float[] in, float[] out) {
  if (in.length == 3) {
    return rgbToHsb(in[0], in[1], in[2], 1.0f, out);
  } else if (in.length == 4) {
    return rgbToHsb(in[0], in[1], in[2], in[3], out);
  }
  return out;
}

public float[] rgbToHsb(float red, float green, float blue, float alpha, float[] out) {

  // Find highest and lowest values.
  float max = max(red, green, blue);
  float min = min(red, green, blue);

  // Find the difference between max and min.
  float delta = max - min;

  // Calculate hue.
  float hue = 0.0f;
  if (delta != 0.0f) {
    if (red == max) {
      hue = (green - blue) / delta;
    } else if (green == max) {
      hue = 2.0f + (blue - red) / delta;
    } else {
      hue = 4.0f + (red - green) / delta;
    }

    hue /= 6.0f;
    if (hue < 0.0f) {
      hue += 1.0f;
    }
  }

  out[0] = hue;
  out[1] = max == 0.0f ? 0.0f : (max - min) / max;
  out[2] = max;
  out[3] = alpha;
  return out;
}

//**********************************AD HOC ******************************************

public ArrayList<PVector> k_nearest_neighbors(PVector p, ArrayList<PVector> points, int k){ //originally used in kenny vaden sketch
    ArrayList<PVector> knn = new ArrayList();
    ArrayList<PVector> tempPoints = new ArrayList(points);

    while(knn.size() < k){
        float minDist = 999999;
        PVector closest = new PVector();
        for(PVector point:tempPoints){
            if (p.equals(point)==false && PVector.dist(p, point) < minDist){
                closest = point; 
                minDist = PVector.dist(p, point);
            }
        }
        knn.add(closest);
        tempPoints.remove(closest);
    }

    return knn;
}

class Circle{ //originally used in kenny vaden sketch
    PVector center;
    float r;
    int c;

    Circle(PVector p, float _r){
        center = p;
        r=_r;
        c = color(0,100,0);
    }

    public boolean isIn(PVector p){
        if(PVector.dist(p, center)<r){
            return true;
        }
        else{return false;}
    }

}

public float compoundTrigFunction(float x, int choice){ //allows compounding of trig functions for interesting lines 
    float coeff1 = 3; //random(2,5);
    float coeff2 = 4; //random(2,5);
    float coeff3 = 5; //random(2,6);
    float coeff4 = 4;// (4,5);
    switch(choice){
      case 0: return cos(x*coeff1+coeff2) - coeff3*sin(randomGaussian()*0.01f + x) + cos(coeff4*x)*pow(sin(pow(x,2)), 2);
      default: return sin(x) + coeff3 * cos(random(2,4)*x);
    }
}


class Ribbon{ //class for drawing a ribbon based on a guide line (as used in flow fields, etc)
    ArrayList<PVector> vertices;

    Ribbon(ArrayList<PVector> guideLine, float stroke_weight, boolean closed){ //should be initialized with ordered set of points 
        // closed bool value set to false if the input line is an open loop (and must be "expanded" into a region)
        //close bool should be set to true if the input list of points forms a closed loop
        vertices = new ArrayList();
        if(closed){
          vertices.addAll(guideLine);
        }
        else{
          ArrayList<PVector> top = new ArrayList();
          ArrayList<PVector> btm = new ArrayList();
          for(int i=1; i<guideLine.size(); i++){
              PVector norm = new PVector(0,0,1).cross(new PVector((guideLine.get(i).x - guideLine.get(i-1).x), (guideLine.get(i).y - guideLine.get(i-1).y))).normalize();
              float theta = norm.heading();
              top.add(new PVector(guideLine.get(i).x + stroke_weight*cos(theta), guideLine.get(i).y + stroke_weight*sin(theta)));
              btm.add(new PVector(guideLine.get(i).x - stroke_weight*cos(theta), guideLine.get(i).y - stroke_weight*sin(theta)));
          }
      
          for(int i=0; i<((top.size() + btm.size())); i++){ // unwrap the top and bottom arrays - first we add all the top points, then start fro the end of hte bottom array to maintain non-self intersection
              vertices.add(
                  i<top.size() ? top.get(i).copy() : btm.get(btm.size()-1-(i-top.size())).copy()
              );
          }
        }
        
    }

    public boolean contains(PVector point){
        return polyPoint(vertices, point.x, point.y);
    }

    public ArrayList<PVector> generatePointsInside(int n){
        ArrayList<PVector> points = new ArrayList();
        int count = 0;
        while(count <= n){
            PVector p = new PVector(random(renderWidth), random(renderHeight));
            if(polyPoint(this.vertices, p.x, p.y)){
                points.add(p);
                count++;
            }
        }
        return points;
    }

    public void display(){
        render.beginShape();
        for(int i=0; i<vertices.size(); i++){
            render.vertex(
                vertices.get(i).x,
                vertices.get(i).y
                );
        }
        render.endShape(CLOSE);
    }

    public void vadenWeb(int n, int _knn, Gradient grad, boolean allow_intersection){
      ArrayList<PVector> points = this.generatePointsInside(n);
        Gradient lineGrad = grad;
        for(PVector p:points){
            ArrayList<PVector> knn = k_nearest_neighbors(p, points, _knn);
            for(PVector k:knn){
              // int baseColor = lineGrad.eval(map(k.x,0,renderWidth,0,1)+randomGaussian()*0.1, HSB);
              // render.stroke(hue(baseColor) + randomGaussian(), saturation(baseColor) + randomGaussian()*8, brightness(baseColor) + randomGaussian()*8);
              if(allow_intersection){
                render.stroke(0);
                render.line(p.x, p.y, k.x, k.y);
              }
              else if(!polyLine(this.vertices, p.x, p.y, k.x, k.y) || polyLine(inkscapePathImport(hole_in_arm_8, 3564.00000f, 5014.66650f), p.x, p.y, k.x, k.y)){
                render.stroke(0);
                render.line(p.x, p.y, k.x, k.y);
              }
            }
      }
    }


}

//**************************** POLYGON AND SUBDIVISION *****************************
public boolean continueSubdivision(float threshold){
  boolean end = false;
  for(Polygon p:polygons){
    if(p.area > (renderWidth *renderHeight)*threshold){
      end = true;
      break;
    }
  }
  return end;
}

class Face{
  PVector p1;
  PVector p2;
  PVector normal;
  PVector center;
  int c = color(0,0,0);
  int nc = color(0,100,100);;
  boolean edge = false;
  int frameLineCount = 0;
  int frameLinesCast = 0;

  Face(float x1, float y1, float x2, float y2, PVector n){
    p1 = new PVector(x1, y1);
    p2 = new PVector(x2, y2);
    center = new PVector((p1.x+p2.x)/2, (p1.y+p2.y)/2);
    normal = n; //used to designate an "outward" direction from the face 
    if((x1==0 && y1==0) || (x1==0 && y1==renderHeight) || (x1==renderWidth && y1==0)) {edge=true;} //allow edges to remain eligible
  }

  Face(float x1, float y1, float x2, float y2){
    p1 = new PVector(x1, y1);
    p2 = new PVector(x2, y2);
    center = new PVector((p1.x+p2.x)/2, (p1.y+p2.y)/2);
    if((x1==0 && y1==0) || (x1==0 && y1==renderHeight) || (x1==renderWidth && y1==0)) {edge=true;}
  }

  public PVector getRandomPoint(float _pct, float std){ //inputs are pct [0,1] and stdev [0,0.5]
    float pct = randomGaussian()*std+_pct; //calc pct based on normal distro
    PVector p = new PVector( //use random pct to pull a point 
        map(pct, 0, 1, p1.x, p2.x),
        map(pct, 0, 1, p1.y, p2.y)
      );
    return p;
  }

  public void display(){
    stroke(c);
    line(p1.x, p1.y, p2.x, p2.y);
    
    if(showNormals){
      stroke(nc);
      int normal_length = 25;
      line(center.x, center.y, center.x+normal_length*normal.x, center.y+normal_length*normal.y); 
      ellipse(center.x+normal_length*normal.x, center.y+normal_length*normal.y, normal_length*0.2f,normal_length*0.2f);
    }
  }

}

class FrameLine extends Face{
  
  FrameLine(float x1, float y1, float x2, float y2){
    super(x1,  y1,  x2,  y2);
    normal = new PVector(0,0,1).cross(new PVector((x2-x1), (y2-y1))).normalize();
    // normal = new PVector(n.x, n.y);
    nc = color(75,100,100);
    c = color(0);
  }
}

class GeometricSubdivision{
  ArrayList<Polygon> polygons;
  ArrayList<Polygon> newPolygons;
  ArrayList<Polygon> deadPolygons;

  GeometricSubdivision(Polygon _p){
    polygons = new ArrayList();
    newPolygons = new ArrayList();
    deadPolygons = new ArrayList();
    this.polygons.add(_p);
  }
  
  public void subdivide(){
    while(continueSubdivision(polygonAreaCutoff)){
        this.newPolygons = new ArrayList();
        for(Polygon p: polygons){
          p.createChildren();
        }
        this.polygons = this.newPolygons;
    }

    for(Polygon p:this.deadPolygons){
      this.polygons.add(p);
    }

    for(Polygon p: this.polygons){
      for(PVector point:p.hull){
        PVector shrinkDirection = PVector.sub(point, p.centroid).normalize();
        float r = 0.9f;
        point.sub(shrinkDirection.mult(PVector.sub(point, p.centroid).mag()*(1-r)));
      }
    }
  }

  public boolean continueSubdivision(float threshold){
    boolean end = false;
    for(Polygon p:this.polygons){
      if(p.area > (renderWidth *renderHeight)*threshold){
        end = true;
        break;
      }
    }
    return end;
  }

  public void display(){
    for(Polygon p: polygons){
        p.display();
    }
  }
}

class Polygon{
  ArrayList<PVector> points;
  ArrayList<Face> edges;
  int c;
  ArrayList<PVector> hull;
  int gen, maxChildren;
  float pct, std, area;
  ArrayList<Polygon> children;
  boolean alive;
  PVector centroid;
  GeometricSubdivision geometricSubdivision;


  Polygon(ArrayList<PVector> _p, boolean ordered){
    if(ordered==false){
      points = _p;
      alive = true;
      hull = new ArrayList();
      edges = new ArrayList();
      children = new ArrayList();
      //TODO: need to solve inheritance mechanism (pass as args, or create and set)
      maxChildren = 5;
      gen=0;
      pct = 0.5f;
      std = 0.1f;
      convexHull();
      area = polygonArea();
      centroid = findCentroid();
      c = color(0,0,100); //colorGrad.eval(map(centroid.x, 0, renderWidth, 0, 1));
      geometricSubdivision = new GeometricSubdivision(this);
    }
    else{
      points = _p;
      alive = true;
      hull = _p;
      edges = new ArrayList();
      for(int i=1; i<hull.size(); i++){ //create FrameLines out of the convex hull
        edges.add(new FrameLine(hull.get(i-1).x, hull.get(i-1).y, hull.get(i).x, hull.get(i).y));
        }
      children = new ArrayList();
      //TODO: need to solve inheritance mechanism (pass as args, or create and set)
      maxChildren = 5;
      gen=0;
      pct = 0.5f;
      std = 0.1f;
      area = polygonArea();
      centroid = findCentroid();
      c = color(0,0,100); //colorGrad.eval(map(centroid.x, 0, renderWidth, 0, 1));
      geometricSubdivision = new GeometricSubdivision(this);
    }
  }

  Polygon(ArrayList<PVector> _p, GeometricSubdivision _g){
      points = _p;
      alive = true;
      hull = new ArrayList();
      edges = new ArrayList();
      children = new ArrayList();
      //TODO: need to solve inheritance mechanism (pass as args, or create and set)
      maxChildren = 5;
      gen=0;
      pct = 0.5f;
      std = 0.1f;
      convexHull();
      area = polygonArea();
      centroid = findCentroid();
      c = color(0,0,100); //colorGrad.eval(map(centroid.x, 0, renderWidth, 0, 1));
      geometricSubdivision = _g;
    }


  public void convexHull(){ //implementing using the Mesh library
    float[][] tempPoints = new float[points.size()][2]; //create 2d array for points to comply with library typing
    for(int i=0; i<points.size(); i++){
      tempPoints[i][0] = points.get(i).x;
      tempPoints[i][1] = points.get(i).y;
    }
    
    Hull tempHull = new Hull(tempPoints); //create hull from points
    int[] temp = tempHull.getExtrema(); // returns hull points as indices of original points array

    for(int i=0; i<temp.length; i++){
      this.hull.add(points.get(temp[i])); //convert indices back to ArrayList of PVectors
    }
    this.hull.add(points.get(temp[0])); //re-add first point for closure

    for(int i=1; i<hull.size(); i++){ //create FrameLines out of the convex hull
      edges.add(new FrameLine(hull.get(i-1).x, hull.get(i-1).y, hull.get(i).x, hull.get(i).y));
    }
  }

  public PVector findCentroid(){
    float x=0;
    float y=0;
    int n = hull.size();
    // Calculate value of shoelace formula
    int j = n - 1;
    for (int i = 0; i < n; i++){
      x += (hull.get(j).x + hull.get(i).x) * (hull.get(j).x * hull.get(i).y - hull.get(i).x*hull.get(j).y);
      y += (hull.get(j).y+hull.get(i).y)*(hull.get(j).x*hull.get(i).y - hull.get(i).x*hull.get(j).y);
      // j is previous vertex to i
      j = i;
    }
    // Return absolute value
    return getTorusPosition(new PVector(Math.abs(x)/(6*this.area),Math.abs(y)/(6*this.area)));
  }

  public void display(){
    render.fill(c,50);
    // render.noStroke();
    render.beginShape();
    for(int i=0; i<hull.size(); i++){
      render.vertex(hull.get(i).x, hull.get(i).y);
    }
    render.endShape(CLOSE);

    // chaikin_line(this.hull, 3, 0.5, "CLOSE");

  }

  public ArrayList<ArrayList<PVector>> split(int numDivisions){
    ArrayList<ArrayList<PVector>> subdivisions = new ArrayList();
    int idx=0; 
    float max_length = 0;
    for(int i=0; i<edges.size(); i++){ //pick the longest edge to start with. this created more well-proportioned subdivisions, less long skinny slices
      if(PVector.sub(edges.get(i).p1, edges.get(i).p2).mag()>max_length){
        max_length = PVector.sub(edges.get(i).p1, edges.get(i).p2).mag();
        idx=i;
      }
    }

    Face f = edges.get(idx); //start with longest face
    PVector np1 = f.getRandomPoint(this.pct, this.std);
    PVector np2 = new PVector();
    boolean goodEdge = false;
    int count=0;
    while(!goodEdge && count < 10000){
      int new_idx = PApplet.parseInt(random(edges.size()-1)); //choose another face randomly (could be improved with some better logic for eligible / preffered faces)
      new_idx = (new_idx>= idx) ? new_idx+1 : new_idx;
      np2 = edges.get(new_idx).getRandomPoint(this.pct, this.std);//select point from new face
      if(!polyLine(this.hull,np1.x, np1.y, np2.x, np2.y)){
        goodEdge = true;
      }
      else{
        // println(count);
      }
      count++;
    }
    // println("esc");
    FrameLine subdivisionEdge = new FrameLine(np1.x, np1.y, np2.x, np2.y);
    edges.add(subdivisionEdge); 
    int firstPointIndex = PApplet.parseInt(random(hull.size())); //first pick random point on the Polygon
    PVector firstPoint = hull.get(firstPointIndex).copy();

    ArrayList<PVector> newPolygon1 = new ArrayList();
    newPolygon1.add(firstPoint); //arbitrarily add this to first new polygon
    ArrayList<PVector> newPolygon2 = new ArrayList();
    for(int i=0; i<hull.size(); i++){
      if(i!=firstPointIndex){ //need to create temp vector from origin of subdivision edge normal vector to each point to test for direction
        if(PVector.dot(PVector.sub(firstPoint,subdivisionEdge.center), subdivisionEdge.normal)*PVector.dot(PVector.sub(hull.get(i), subdivisionEdge.center), subdivisionEdge.normal)>0){//if they are facing the same way as the first point wrt subdivisionEdge plane 
          newPolygon1.add(hull.get(i).copy()); //then we want to add to polygon 1
        }
        else{
          newPolygon2.add(hull.get(i).copy());//otherwise we add to second polygon
        }
      }
    }
    // add both points from the subdividing edge to both new polygons
    newPolygon1.add(subdivisionEdge.p1);
    newPolygon1.add(subdivisionEdge.p2);
    newPolygon2.add(subdivisionEdge.p1);
    newPolygon2.add(subdivisionEdge.p2);
    //now we have to collect all the points 

    subdivisions.add(newPolygon1);
    subdivisions.add(newPolygon2);

    while(subdivisions.size()<numDivisions){
      int childIndex = PApplet.parseInt(random(subdivisions.size()));
      Polygon temp = new Polygon(subdivisions.get(childIndex), false);
      ArrayList<ArrayList<PVector>> newSubdivisionsFromChild = temp.split(2); //call this function to reutrn the subdivisions from the removed set of points
      for(ArrayList<PVector> p:newSubdivisionsFromChild){
        subdivisions.add(p);
      }
      subdivisions.remove(childIndex); //remove the set of poitns that is being subdivided to avoid double counting 
    }
    this.alive = false;
    return subdivisions;

  }

  public void subdivide(){
    geometricSubdivision.subdivide();
  }

  public float polygonArea(){
    area = 0.0f;
    int n = hull.size();
    // Calculate value of shoelace formula
    int j = n - 1;
    for (int i = 0; i < n; i++){
      area += (hull.get(j).x + hull.get(i).x) * (hull.get(j).y - hull.get(i).y);
        
      // j is previous vertex to i
      j = i;
    }

    // Return absolute value
    return Math.abs(area / 2.0f);
    
  }

  public void createChildren(){ 

    if(this.alive && this.area > (renderWidth*renderHeight)*polygonAreaCutoff){
      ArrayList<ArrayList<PVector>> subdivisions = this.split(this.maxChildren);
      int idx = PApplet.parseInt(random(subdivisions.size()));
      for(int i=0; i<subdivisions.size(); i++){
        Polygon p = new Polygon(subdivisions.get(i), this.geometricSubdivision);
        p.pct = this.pct;
        p.std = this.std;
        p.maxChildren = this.maxChildren;
        p.gen = this.gen+1;
        // p.c = (random(1)<map(noiseGrid[constrain(ceil(p.centroid.x),0,renderWidth-1)][constrain(ceil(p.centroid.y),0,renderHeight-1)].y, 0, 1, 0, 0.2))? palette[int(random(palette.length))] : color(hue(c), saturation(c) + 1, brightness(c) - 1); //this could be driven off some noise or other flow field
        p.alive = (idx==i && random(1) < 0.5f && p.gen > 1) ? false : true;
        this.geometricSubdivision.newPolygons.add(p);
        children.add(p);
      }
    }
    else{
      this.geometricSubdivision.deadPolygons.add(this);
    }
  }

  public void displayChildren(){
    for(Polygon c:children){
      if(c.alive){c.display();}
    }
  }
  
}
float[][] branch5 = new float[][]{
{2215.2573f, 3544.0385f},
{2413.5583f, 3607.4816f},
{2609.3978f, 3636.6312000000003f},
{2790.3712f, 3561.5238000000004f},
{2866.9867f, 3576.1536000000006f},
{2993.488f, 3461.2590000000005f},
{2825.3693f, 3613.1379000000006f},
{2685.3289f, 3675.7652000000007f},
{2889.5977f, 3709.670200000001f},
{3097.0305f, 3736.7560000000008f},
{3261.196f, 3719.5314000000008f},
{3400.8477f, 3609.6294000000007f},
{3279.8759f, 3757.683200000001f},
{3233.4437f, 3846.965600000001f},
{3394.9327999999996f, 3956.646200000001f},
{3497.6803999999997f, 4093.0447000000013f},
{3318.9685f, 3977.889300000001f},
{3175.8711f, 3833.736400000001f},
{2996.6650999999997f, 3763.7511000000013f},
{2791.2522999999997f, 3793.6414000000013f},
{2588.5368f, 3742.0634000000014f},
{2363.7653f, 3773.1915000000013f},
{2147.7816000000003f, 3692.9905000000012f},
{2007.0465000000004f, 3539.546000000001f},
{1932.6057000000003f, 3346.5933000000014f},
{1904.8160000000003f, 3106.639500000001f},
{2036.3552000000002f, 2944.601700000001f},
{2109.8113000000003f, 3107.619700000001f},
{2126.7223000000004f, 3311.701700000001f},
{2177.8337f, 3495.397500000001f},
{2215.2573f, 3544.038500000001f},
{2215.2573f, 3544.0385f}
};


float[][] branch4 = new float[][]{
{2135.0079f, 3678.4096f},
{2209.9209f, 3824.1808f},
{2285.0021f, 4012.5571f},
{2159.0264f, 4173.7859f},
{2167.9674f, 4400.9556f},
{2198.4864f, 4434.050200000001f},
{2169.0496999999996f, 4246.3612f},
{2302.9717999999993f, 4100.185600000001f},
{2346.4421999999995f, 4302.913300000001f},
{2317.9576999999995f, 4545.254800000001f},
{2302.8671999999997f, 4773.950000000001f},
{2313.1169999999997f, 4923.3379f},
{2310.0989999999997f, 4751.6601f},
{2361.7443f, 4543.4057f},
{2431.7841f, 4352.596100000001f},
{2506.6393f, 4502.1431f},
{2438.3898f, 4321.7915f},
{2373.6793f, 4121.8543f},
{2408.9476f, 4063.7164f},
{2555.9296f, 4161.3368f},
{2768.6721f, 4237.1829f},
{2856.2007f, 4366.4743f},
{2860.2603f, 4591.5167f},
{2848.6886999999997f, 4789.4578f},
{2964.6766999999995f, 4993.4372f},
{2983.3619999999996f, 5010.304f},
{2901.6200999999996f, 4807.2825f},
{2910.2060999999994f, 4593.0335000000005f},
{3108.3213999999994f, 4598.391100000001f},
{3162.6613999999995f, 4800.761f},
{3257.4327999999996f, 4973.3536f},
{3428.8215999999998f, 5042.9619f},
{3339.2075999999997f, 4971.6224f},
{3190.0606f, 4855.0244f},
{3183.7023999999997f, 4641.7267f},
{3054.1559999999995f, 4504.6259f},
{2934.6866999999993f, 4376.4804f},
{2920.7057999999993f, 4195.3652f},
{2698.400899999999f, 4121.1195f},
{2638.767499999999f, 4056.2342f},
{2823.165499999999f, 4064.0119999999997f},
{3066.199799999999f, 4042.2344f},
{3166.892399999999f, 3984.709f},
{3303.8799999999987f, 4080.0101f},
{3429.4363999999987f, 4253.1208f},
{3448.2928999999986f, 4287.8489f},
{3351.7394999999988f, 4117.5574f},
{3216.966199999999f, 3955.4446999999996f},
{3444.563099999999f, 3852.4867999999997f},
{3545.887999999999f, 3655.0460999999996f},
{3550.097399999999f, 3631.8479999999995f},
{3399.705699999999f, 3835.0750999999996f},
{3177.100399999999f, 3927.0156999999995f},
{2975.441499999999f, 4020.8647999999994f},
{2734.315299999999f, 3984.2053999999994f},
{2536.3115999999986f, 4003.185299999999f},
{2392.492899999999f, 3807.834499999999f},
{2236.947499999999f, 3670.153299999999f},
{2135.0078999999987f, 3678.409599999999f},
{2135.0079f, 3678.4096f}
};

float[][] branch3 = new float[][]{
{2026.7645f, 3650.4156f},
{2058.4782f, 3855.254f},
{1963.5114f, 4014.6502f},
{1830.8584f, 4118.4726f},
{1751.5788f, 4262.731f},
{1667.3677f, 4460.1838f},
{1498.5132f, 4545.986f},
{1354.7655f, 4694.1134999999995f},
{1259.7293f, 4844.825599999999f},
{1275.1397f, 4852.831999999999f},
{1390.4044999999999f, 4707.1757f},
{1418.4454999999998f, 4690.819399999999f},
{1355.1272f, 4863.252799999999f},
{1183.7135999999998f, 4990.4486f},
{1122.5249f, 5141.0756f},
{1263.5455f, 4975.2245f},
{1396.0836f, 4868.7808f},
{1430.5056f, 4695.8819f},
{1556.1943999999999f, 4581.7331f},
{1709.5311f, 4487.2482f},
{1786.2650999999998f, 4314.9829f},
{1907.0652999999998f, 4162.6382f},
{2022.8471999999997f, 4121.754f},
{2017.3375999999996f, 4346.3785f},
{1897.2708999999995f, 4469.912499999999f},
{1750.2583999999995f, 4629.9821999999995f},
{1670.2448999999995f, 4798.404299999999f},
{1622.7664999999995f, 4968.310299999999f},
{1629.3625999999995f, 4992.0181999999995f},
{1693.2162999999996f, 4823.8018999999995f},
{1796.3143999999995f, 4679.223099999999f},
{1900.2487999999996f, 4537.7501999999995f},
{2077.7751f, 4436.188099999999f},
{2061.694f, 4533.757799999999f},
{1931.757f, 4695.217699999999f},
{1846.9538f, 4834.387899999999f},
{1722.5632f, 5009.056999999999f},
{1591.1896000000002f, 5085.6745999999985f},
{1789.0420000000001f, 5000.303799999999f},
{1922.3329f, 4907.609299999999f},
{1959.5646000000002f, 4708.665799999999f},
{2102.6743f, 4526.470999999999f},
{2129.1229000000003f, 4350.985499999999f},
{2110.5723000000003f, 4060.3538999999987f},
{2188.429f, 3874.6422999999986f},
{2136.8051f, 3647.4634999999985f},
{2026.7645f, 3650.4155999999984f},
{2026.7645f, 3650.4156f}
};

float [][] branch2 = new float[][]{
{1994.2528f, 3583.8455f},
{1953.8171f, 3793.1979f},
{1793.4803f, 3973.011f},
{1631.7475f, 4062.6072f},
{1391.1991f, 4085.7639f},
{1270.7037f, 4250.6254f},
{1431.045f, 4112.1166f},
{1640.938f, 4119.1862f},
{1419.9043000000001f, 4210.5455f},
{1357.3653000000002f, 4280.4839f},
{1159.4389f, 4312.5258f},
{935.1975100000001f, 4274.31f},
{748.3725300000001f, 4198.1646f},
{691.4501000000001f, 4207.617f},
{887.5875600000002f, 4314.3871f},
{1083.3904000000002f, 4337.7824f},
{879.8989800000002f, 4396.5615f},
{655.9931800000002f, 4386.6239f},
{422.9270700000002f, 4355.3243999999995f},
{204.19219000000018f, 4357.825199999999f},
{301.2247300000002f, 4345.739299999999f},
{476.6468400000002f, 4393.655899999999f},
{428.6311300000002f, 4563.2532999999985f},
{564.7735900000002f, 4394.293699999998f},
{638.7407800000002f, 4512.042099999999f},
{496.2705200000002f, 4721.381599999999f},
{393.4964300000002f, 4903.624599999999f},
{562.8994800000003f, 4717.1855f},
{665.1466100000002f, 4539.6132f},
{779.6332600000003f, 4438.1852f},
{754.4284200000003f, 4580.2889f},
{749.6279400000003f, 4690.2393999999995f},
{840.6453400000003f, 4490.428099999999f},
{902.7145000000003f, 4479.047099999999f},
{990.9289500000002f, 4553.780899999999f},
{1095.0487000000003f, 4389.601699999999f},
{1308.9034000000004f, 4401.570199999999f},
{1312.7505000000003f, 4510.981199999999f},
{1183.6391000000003f, 4710.172f},
{1174.0237000000004f, 4756.038299999999f},
{1323.7960000000005f, 4547.485199999999f},
{1389.6659000000004f, 4363.3018999999995f},
{1513.0660000000005f, 4193.6977f},
{1727.7316000000005f, 4177.8461f},
{1854.8196000000005f, 4066.2448f},
{2012.1968000000004f, 3931.4455f},
{2075.2385000000004f, 3743.6211999999996f},
{1994.2528000000004f, 3583.8454999999994f},
{1994.2528f, 3583.8455f}
};

float [][] branch1 = new float[][]{
{1795.3476f, 3712.0024f},
{1618.7587f, 3722.9640999999997f},
{1484.5091000000002f, 3681.5955f},
{1331.3434000000002f, 3697.2981999999997f},
{1179.0197000000003f, 3787.7115999999996f},
{1051.5441000000003f, 3861.7744f},
{907.0051100000003f, 3954.6169f},
{834.0606200000003f, 4105.6175f},
{842.2368300000003f, 4131.358200000001f},
{956.4433700000003f, 3974.1943000000006f},
{1090.0128000000004f, 3896.1499000000003f},
{1224.8065000000004f, 3821.8733f},
{1417.7468000000003f, 3756.1372f},
{1563.7613000000003f, 3816.5386000000003f},
{1730.2815000000003f, 3847.4490000000005f},
{1812.9106000000002f, 3831.8100000000004f},
{1972.6428f, 3747.4614000000006f},
{2009.4787000000001f, 3607.4501000000005f},
{1816.3057000000001f, 3681.4654000000005f},
{1795.3476f, 3712.0024000000003f},
{1795.3476f, 3712.0024f}
};

float[][] arm1 = new float [][]{
{382.58446f, 1297.0546f},
{347.47272f, 1470.0165f},
{421.54467999999997f, 1619.6536999999998f},
{338.82122f, 1799.5796999999998f},
{300.57765f, 1967.9969999999998f},
{403.15022f, 2136.9716f},
{586.27181f, 2221.5987999999998f},
{800.3274299999999f, 2196.6695f},
{1012.4735f, 2121.9834f},
{1202.5997f, 2108.4323f},
{1385.8204f, 2099.7744f},
{1521.8300000000002f, 2242.6607f},
{1703.0913f, 2253.5648f},
{1899.4616f, 2244.4955f},
{1959.6401f, 2134.267f},
{1857.0217f, 1978.4272999999998f},
{1678.5159f, 1888.5940999999998f},
{1474.3716000000002f, 1821.0198999999998f},
{1243.1663f, 1813.1894999999997f},
{1043.7624f, 1882.9103999999998f},
{835.2873900000001f, 1934.8496999999998f},
{641.0762800000001f, 1984.8454999999997f},
{497.2113600000001f, 1939.2777999999996f},
{506.83376000000015f, 1764.2727999999997f},
{519.1692900000002f, 1585.5530999999996f},
{369.8731500000001f, 1451.5452999999995f},
{380.82219000000015f, 1301.7785999999996f},
{381.61889000000014f, 1299.3771999999997f},
{382.58446f, 1297.0546f}
};

float[][] arm8 = new float[][]{
{2329.1807f, 2533.7263f},
{2356.5949f, 2516.8855999999996f},
{2413.8574f, 2505.7249999999995f},
{2479.6895f, 2523.6463999999996f},
{2519.0795f, 2529.1107999999995f},
{2566.3952f, 2542.7377999999994f},
{2600.9827999999998f, 2547.865399999999f},
{2633.9656999999997f, 2553.649699999999f},
{2676.2484f, 2558.799599999999f},
{2723.1002f, 2562.086699999999f},
{2724.0238f, 2564.449299999999f},
{2781.9195999999997f, 2578.4037999999987f},
{2837.6639999999998f, 2589.074099999999f},
{2885.9156f, 2595.1544999999987f},
{2986.0458999999996f, 2608.1975999999986f},
{3023.4480999999996f, 2597.6883999999986f},
{3066.8832999999995f, 2579.1362999999988f},
{3097.7619999999997f, 2568.823899999999f},
{3150.854f, 2555.085799999999f},
{3180.3543f, 2531.086999999999f},
{3202.099f, 2513.687099999999f},
{3247.7687f, 2477.1486999999993f},
{3269.7309f, 2450.0313999999994f},
{3293.1113f, 2424.6428999999994f},
{3305.3944f, 2401.9724999999994f},
{3322.5069000000003f, 2368.5056999999993f},
{3331.5368000000003f, 2342.6828999999993f},
{3338.7122000000004f, 2280.3535999999995f},
{3339.9710000000005f, 2220.8060999999993f},
{3337.4506000000006f, 2200.5905999999995f},
{3337.2716000000005f, 2169.9659999999994f},
{3326.8354000000004f, 2132.552899999999f},
{3312.4918000000002f, 2092.646599999999f},
{3320.8022f, 2070.8622999999993f},
{3331.8679f, 2041.2665999999992f},
{3344.4215000000004f, 2015.9423999999992f},
{3358.7273000000005f, 1986.4642999999992f},
{3372.6138000000005f, 1960.238299999999f},
{3378.9760000000006f, 1922.839999999999f},
{3404.4931000000006f, 1862.947199999999f},
{3423.5577000000008f, 1832.187899999999f},
{3465.3985000000007f, 1806.599599999999f},
{3509.9713000000006f, 1790.1172999999992f},
{3540.359200000001f, 1776.1626999999992f},
{3508.544300000001f, 1777.1839999999993f},
{3463.5015000000008f, 1781.5872999999992f},
{3429.093100000001f, 1797.0016999999993f},
{3392.114100000001f, 1834.6280999999994f},
{3344.0616000000014f, 1913.3841999999995f},
{3329.4747000000016f, 1943.8430999999996f},
{3315.5332000000017f, 1976.3640999999996f},
{3277.8908000000015f, 2031.5802999999996f},
{3271.8905000000013f, 2024.3270999999995f},
{3249.2947000000013f, 2005.8575999999996f},
{3218.624200000001f, 1975.5146999999995f},
{3190.131000000001f, 1948.0393999999994f},
{3134.794800000001f, 1925.6582999999994f},
{3111.636200000001f, 1915.2856999999995f},
{3081.4500000000007f, 1910.8996999999995f},
{3042.078200000001f, 1910.2337999999995f},
{3009.8187000000007f, 1919.2255999999995f},
{2982.4608000000007f, 1938.1281999999994f},
{2938.8585000000007f, 1997.9486999999995f},
{2929.338300000001f, 2042.0621999999994f},
{2933.9189000000006f, 2075.2969999999996f},
{2944.1371000000004f, 2103.5206999999996f},
{2960.9492000000005f, 2127.1076999999996f},
{2978.9761000000003f, 2151.1856999999995f},
{3014.3853000000004f, 2177.3028999999997f},
{3042.4235000000003f, 2198.7327999999998f},
{3077.9638000000004f, 2216.4997f},
{3139.5548000000003f, 2217.9094999999998f},
{3171.9846000000002f, 2208.8289999999997f},
{3196.1592f, 2196.513f},
{3218.9091f, 2175.5879999999997f},
{3175.1616f, 2269.8574f},
{3138.2673f, 2309.3399999999997f},
{3068.068f, 2325.8448f},
{3016.9359000000004f, 2333.9323f},
{2955.2285000000006f, 2338.9695f},
{2910.4176000000007f, 2338.6374f},
{2877.9344000000006f, 2333.8642f},
{2805.5949000000005f, 2327.2137f},
{2752.1389000000004f, 2324.1744999999996f},
{2712.4930000000004f, 2311.6886999999997f},
{2655.1718000000005f, 2293.0914999999995f},
{2596.9822000000004f, 2288.7308999999996f},
{2542.1071f, 2279.4975999999997f},
{2483.1736f, 2261.8947999999996f},
{2373.1823f, 2256.8679999999995f},
{2301.0931f, 2248.8322999999996f},
{2185.3902f, 2275.0089999999996f},
{2143.9372999999996f, 2366.0946999999996f},
{2158.4159999999997f, 2478.1858999999995f},
{2209.3068f, 2558.2920999999997f},
{2295.4721f, 2553.6422999999995f},
{2329.1807f, 2533.7262999999994f},
{2329.1807f, 2533.7263f}
};

float[][] arm7 = new float[][]{
{3471.2541f, 936.86536f},
{3370.2335000000003f, 1021.5364f},
{3326.0538f, 1061.5718f},
{3295.106f, 1089.845f},
{3277.5477f, 1166.0358f},
{3290.7588f, 1232.7495000000001f},
{3299.7993f, 1305.1769000000002f},
{3292.9964f, 1352.9651000000001f},
{3272.1526f, 1385.1082000000001f},
{3242.7269f, 1420.1472f},
{3206.2444f, 1446.3559f},
{3163.0803f, 1480.8648f},
{3128.3296f, 1482.409f},
{3096.2326f, 1471.7869f},
{3056.9197999999997f, 1470.0092f},
{2966.6948999999995f, 1524.2285f},
{2904.6805999999997f, 1592.0225f},
{2839.0613999999996f, 1681.3423f},
{2790.6595999999995f, 1737.5978f},
{2751.2330999999995f, 1762.6456f},
{2754.3214999999996f, 1758.5546000000002f},
{2719.1490999999996f, 1778.5512f},
{2638.7389999999996f, 1826.3161f},
{2554.2886999999996f, 1853.8798f},
{2496.8477f, 1869.6968f},
{2416.3469f, 1905.2011f},
{2341.366f, 1936.1592f},
{2341.323f, 1934.441f},
{2301.1835f, 1943.6436f},
{2255.9724f, 1950.4369000000002f},
{2229.1269f, 2056.0936f},
{2302.4985f, 2213.3417f},
{2361.8694f, 2243.5405f},
{2408.2566f, 2194.3824f},
{2425.9922f, 2162.0999f},
{2484.3349000000003f, 2113.27f},
{2528.5317000000005f, 2096.6516f},
{2626.0866000000005f, 2074.7378000000003f},
{2678.5545000000006f, 2057.4709000000003f},
{2715.7458000000006f, 2044.7380000000003f},
{2783.6486000000004f, 2017.8051000000003f},
{2828.0613000000003f, 2002.0788000000002f},
{2877.4749f, 1979.1434000000002f},
{2904.596f, 1955.9471f},
{2942.7061f, 1913.924f},
{2987.4912f, 1839.1761f},
{3029.9151f, 1772.9643999999998f},
{3051.5267000000003f, 1734.5727f},
{3080.6827000000003f, 1706.1752f},
{3125.0251000000003f, 1671.1588f},
{3178.8729000000003f, 1626.1924999999999f},
{3232.9324f, 1595.0180999999998f},
{3285.8227f, 1562.5326999999997f},
{3323.4284000000002f, 1536.1749999999997f},
{3348.4872f, 1503.9372999999998f},
{3376.7576f, 1465.3697999999997f},
{3391.0622f, 1426.7989999999998f},
{3405.9343f, 1392.2341999999996f},
{3414.1978f, 1350.7987999999996f},
{3402.0586f, 1310.6543999999997f},
{3370.0924f, 1276.8765999999996f},
{3332.6469f, 1218.5825999999995f},
{3326.4258f, 1151.4755999999995f},
{3328.5614f, 1104.3558999999996f},
{3382.3714f, 1051.7292999999995f},
{3400.6591f, 1040.7814999999996f},
{3401.7963999999997f, 1037.6385999999995f},
{3461.7153999999996f, 983.3137499999996f},
{3494.1575f, 944.8750499999995f},
{3494.9865f, 852.6067299999995f},
{3483.2122f, 777.9698599999995f},
{3469.1378f, 797.5905999999994f},
{3473.8628f, 928.7132499999994f},
{3471.2536f, 936.8653599999994f},
{3471.2541f, 936.86536f}
};

float[][] arm6 = new float[][]{
{2883.7141f, 691.73371f},
{2906.9888f, 721.11369f},
{2924.3952f, 760.9874f},
{2949.0463999999997f, 814.79033f},
{2974.8275f, 867.81353f},
{2986.4932f, 933.95991f},
{2989.7576f, 972.3244100000001f},
{2987.2038f, 1012.8855000000001f},
{2988.5018999999998f, 1056.2042000000001f},
{2969.2286999999997f, 1097.3642000000002f},
{2942.1265999999996f, 1152.4018f},
{2899.7723999999994f, 1248.2455f},
{2867.2556999999993f, 1287.0214f},
{2852.1956999999993f, 1333.0344f},
{2829.5028999999995f, 1363.2428f},
{2806.6225999999997f, 1396.568f},
{2775.6310999999996f, 1423.3497f},
{2720.5978999999998f, 1481.0578f},
{2654.9239f, 1558.3118f},
{2519.7531999999997f, 1657.2088999999999f},
{2443.6569999999997f, 1713.5561999999998f},
{2384.9902999999995f, 1763.0043999999998f},
{2311.0816999999993f, 1849.6040999999998f},
{2269.830999999999f, 1907.9374999999998f},
{2222.860899999999f, 1906.1887999999997f},
{2184.8060999999993f, 1888.6287999999997f},
{2120.411899999999f, 1749.6402999999998f},
{2142.796699999999f, 1670.417f},
{2210.787299999999f, 1598.4574f},
{2271.794399999999f, 1561.1091f},
{2324.433899999999f, 1536.8698f},
{2360.029699999999f, 1520.5025f},
{2416.9092999999993f, 1486.4764f},
{2453.5047999999992f, 1470.2979f},
{2541.2990999999993f, 1425.1949f},
{2584.4375999999993f, 1404.6308f},
{2648.735399999999f, 1357.783f},
{2684.076799999999f, 1332.4128999999998f},
{2759.559799999999f, 1279.2818999999997f},
{2775.329799999999f, 1247.2946999999997f},
{2875.915399999999f, 1080.5677999999998f},
{2884.732499999999f, 1031.3144999999997f},
{2903.308199999999f, 844.3190199999997f},
{2892.873899999999f, 767.8023499999997f},
{2875.602799999999f, 722.8740799999997f},
{2844.329399999999f, 687.6173199999997f},
{2818.380099999999f, 656.3806499999997f},
{2794.275499999999f, 629.7976099999997f},
{2741.241299999999f, 548.8945699999997f},
{2723.1926999999987f, 514.7782799999997f},
{2699.086399999999f, 477.04475999999966f},
{2698.6200999999987f, 423.64229999999964f},
{2702.5993999999987f, 373.22115999999966f},
{2748.7400999999986f, 295.5305899999997f},
{2771.2698999999984f, 262.8240199999997f},
{2844.7446999999984f, 196.4255999999997f},
{2853.8962999999985f, 207.1587799999997f},
{2846.6416999999983f, 204.0971799999997f},
{2843.5439999999985f, 206.1980099999997f},
{2790.1931999999983f, 270.2395199999997f},
{2731.3516999999983f, 348.7727599999997f},
{2724.9731999999985f, 391.7007699999997f},
{2713.2121999999986f, 454.0327099999997f},
{2752.9157999999984f, 539.3045499999997f},
{2789.1808999999985f, 566.7468999999998f},
{2836.2035999999985f, 613.1134599999998f},
{2863.6360999999984f, 641.9418799999999f},
{2892.7954999999984f, 685.9940399999998f},
{2883.7140999999983f, 691.7337099999999f},
{2883.7141f, 691.73371f}
};

float[][] arm5 = new float[][]{
{2032.3633f, 675.58817f},
{2053.9432f, 642.4028f},
{2065.9561000000003f, 561.4715799999999f},
{2065.7177f, 479.8930499999999f},
{2041.8219000000001f, 444.8747199999999f},
{2010.0268f, 408.13544999999993f},
{1972.6428f, 382.5844599999999f},
{1949.2777f, 343.91306999999995f},
{1944.0349f, 331.90680999999995f},
{1936.7758000000001f, 312.02190999999993f},
{1925.9862f, 300.46876999999995f},
{1928.4057f, 300.06665999999996f},
{1962.4793f, 342.39732f},
{2005.1290999999999f, 387.81624999999997f},
{2042.3273f, 427.60337999999996f},
{2081.1005f, 460.39804999999996f},
{2102.6291f, 500.72347999999994f},
{2130.6318f, 604.89438f},
{2132.0398f, 662.7027499999999f},
{2105.9614f, 733.87727f},
{2062.2236000000003f, 823.02315f},
{2019.0472000000002f, 903.84646f},
{1997.4632000000001f, 949.71767f},
{1983.8404f, 996.58586f},
{1966.4814000000001f, 1049.3052f},
{1953.7535f, 1238.0975f},
{1979.6111f, 1333.8617000000002f},
{2006.0492000000002f, 1403.8808000000001f},
{2023.0611000000001f, 1467.8515000000002f},
{2096.1319000000003f, 1623.4938000000002f},
{2119.2123f, 1700.6901000000003f},
{2135.8788f, 1756.7189000000003f},
{2141.7434f, 1787.8434000000002f},
{2169.0067999999997f, 1892.6749000000002f},
{2188.2545999999998f, 1953.5786000000003f},
{2205.3900999999996f, 2046.5522000000003f},
{2211.9913999999994f, 2095.8163000000004f},
{2228.5972999999994f, 2145.2978000000003f},
{2241.1000999999997f, 2198.4254f},
{2235.9115999999995f, 2250.2745f},
{2216.0418999999993f, 2301.8718f},
{2177.9319999999993f, 2426.1454f},
{2166.9955999999993f, 2460.4915f},
{2123.9393999999993f, 2563.6006f},
{2079.1889999999994f, 2582.2340000000004f},
{2000.4999999999993f, 2548.0825000000004f},
{1970.2407999999994f, 2512.1442000000006f},
{1961.3097999999993f, 2447.2155000000007f},
{1966.7238999999993f, 2408.3684000000007f},
{1957.2701999999992f, 2365.9019000000008f},
{1945.2415999999992f, 2287.8268000000007f},
{1938.299999999999f, 2141.0340000000006f},
{1920.387399999999f, 2086.4850000000006f},
{1904.139699999999f, 2024.9545000000005f},
{1898.6292999999991f, 1983.5123000000006f},
{1898.5196999999991f, 1964.3331000000005f},
{1878.8029999999992f, 1909.0188000000005f},
{1857.5344999999993f, 1859.8776000000005f},
{1839.3870999999992f, 1822.3038000000006f},
{1819.0907999999993f, 1769.6201000000005f},
{1797.5489999999993f, 1694.1050000000005f},
{1776.1063999999992f, 1602.6821000000004f},
{1760.5280999999993f, 1559.5832000000005f},
{1754.7307999999994f, 1451.2250000000006f},
{1755.2887999999994f, 1379.9432000000006f},
{1765.3602999999994f, 1330.4530000000007f},
{1761.1540999999993f, 1303.5692000000006f},
{1788.4230999999993f, 1128.4759000000006f},
{1785.1706999999992f, 1073.4420000000007f},
{1809.0236999999993f, 1032.3423000000007f},
{1829.0050999999992f, 966.0503400000007f},
{1846.4595999999992f, 933.0011700000007f},
{1878.5064999999993f, 882.5668900000007f},
{1907.1536999999994f, 842.5863800000008f},
{1962.4271999999994f, 769.8247100000008f},
{2015.5175999999994f, 695.4856700000008f},
{2032.5112999999994f, 674.0279800000009f},
{2032.3632999999995f, 675.5881700000009f},
{2032.3633f, 675.58817f}
};

float[][] arm4 = new float[][]{
{1564.8638f, 1290.5227f},
{1577.1566f, 1245.9302f},
{1588.1921f, 1186.0118f},
{1593.8055f, 1130.1369f},
{1587.4133f, 1081.8996f},
{1563.6609999999998f, 1047.5819f},
{1539.0647999999999f, 1002.7015999999999f},
{1522.4144999999999f, 959.6469099999998f},
{1507.0095f, 944.3304199999998f},
{1479.8223999999998f, 903.4409799999997f},
{1441.5690999999997f, 870.6199599999998f},
{1439.0539999999996f, 869.0639199999997f},
{1409.1583999999996f, 824.7538499999997f},
{1380.6238999999996f, 788.6791199999997f},
{1356.8293999999996f, 746.5132899999996f},
{1337.1792999999996f, 693.3176399999996f},
{1314.3521999999996f, 648.6458599999996f},
{1311.2200999999995f, 580.8391499999997f},
{1321.6567999999995f, 532.0594199999997f},
{1333.4466999999995f, 489.2009899999997f},
{1352.1094999999996f, 418.9765899999997f},
{1359.7438999999995f, 374.0659399999997f},
{1385.3359999999996f, 330.99039999999974f},
{1412.5954999999994f, 273.11496999999974f},
{1453.2316999999994f, 222.50292999999974f},
{1468.6255999999994f, 168.23304999999974f},
{1477.9175999999993f, 59.00109799999973f},
{1484.5668999999994f, 47.16621399999973f},
{1477.3139999999994f, 192.82657999999975f},
{1444.2460999999994f, 269.5285899999998f},
{1425.3142999999993f, 318.23615999999976f},
{1402.0550999999994f, 367.17048999999975f},
{1374.5818999999995f, 459.95913999999976f},
{1366.5730999999994f, 545.8826599999998f},
{1379.0308999999993f, 616.0163399999998f},
{1397.6768999999993f, 649.5487799999997f},
{1423.5169999999991f, 679.2391299999997f},
{1440.9587999999992f, 695.1671899999997f},
{1486.4703999999992f, 744.8861899999997f},
{1546.5075999999992f, 803.5734599999997f},
{1592.9845999999993f, 836.6716699999997f},
{1655.3898999999992f, 904.8756399999997f},
{1655.040699999999f, 899.3973799999998f},
{1780.4173999999991f, 1130.9569999999999f},
{1782.2430999999992f, 1209.3438999999998f},
{1790.2123999999992f, 1268.0031f},
{1791.614999999999f, 1369.839f},
{1767.9397999999992f, 1409.1065999999998f},
{1702.5731999999991f, 1492.2275f},
{1694.703999999999f, 1549.3385f},
{1664.990299999999f, 1639.5195f},
{1652.095099999999f, 1691.8442f},
{1664.852499999999f, 1750.5612f},
{1682.856299999999f, 1797.1939000000002f},
{1703.174599999999f, 1846.0051000000003f},
{1798.904499999999f, 2010.0464000000004f},
{1863.081399999999f, 2082.9682000000003f},
{1909.320999999999f, 2115.2003000000004f},
{1935.309799999999f, 2160.9631000000004f},
{1992.000299999999f, 2216.5600000000004f},
{2029.130999999999f, 2271.9029000000005f},
{2015.146299999999f, 2301.2842000000005f},
{1991.5153999999989f, 2340.0198000000005f},
{1944.5679999999988f, 2367.3584000000005f},
{1867.0883999999987f, 2369.7194000000004f},
{1656.6154999999987f, 2367.9057000000003f},
{1619.9185999999988f, 2344.0297f},
{1595.841199999999f, 2312.3779f},
{1565.0570999999989f, 2268.9905f},
{1545.5859999999989f, 2246.5933f},
{1501.577499999999f, 2188.7379f},
{1450.707299999999f, 2121.9866f},
{1428.5172999999988f, 2077.5782000000004f},
{1404.810499999999f, 2013.3998000000004f},
{1389.3233999999989f, 1952.0436000000004f},
{1344.2808999999988f, 1853.7582000000004f},
{1330.9649999999988f, 1792.9018000000005f},
{1316.1011999999987f, 1749.1775000000005f},
{1317.3144999999986f, 1685.1769000000004f},
{1337.9129999999986f, 1627.9148000000005f},
{1373.6030999999987f, 1582.7430000000004f},
{1381.9696999999987f, 1566.7300000000005f},
{1377.3040999999987f, 1568.5963000000004f},
{1410.3575999999987f, 1524.9844000000003f},
{1461.1881999999987f, 1426.4284000000002f},
{1502.3438999999987f, 1376.3709000000003f},
{1538.3279999999986f, 1327.4264000000003f},
{1563.0752999999986f, 1291.1793000000002f},
{1564.8637999999987f, 1290.5227000000002f},
{1564.8638f, 1290.5227f}
};

float[][] arm3 = new float[][]{
{2144.1537f, 2089.1787f},
{2114.0557f, 2021.8881f},
{2046.8781999999999f, 1975.8610999999999f},
{1976.3953999999999f, 1941.3164f},
{1895.3535f, 1906.08f},
{1830.8411999999998f, 1854.5301f},
{1772.8944999999999f, 1821.2418f},
{1672.2022f, 1754.7124000000001f},
{1604.6933f, 1704.9867000000002f},
{1548.3616f, 1672.4644f},
{1486.9996999999998f, 1635.2267000000002f},
{1427.7684f, 1590.3055000000002f},
{1380.1716999999999f, 1540.0693f},
{1332.0205999999998f, 1481.498f},
{1255.6014999999998f, 1422.394f},
{1236.5288999999998f, 1339.046f},
{1224.8689999999997f, 1294.1971f},
{1198.2413999999997f, 1206.1593f},
{1238.2604999999996f, 1043.1063f},
{1296.7716999999996f, 939.3373099999999f},
{1333.6792999999996f, 894.5623799999998f},
{1378.4934999999996f, 833.6957999999998f},
{1414.6637999999996f, 770.6750699999999f},
{1454.1738999999995f, 702.8871399999999f},
{1450.5538999999997f, 567.86326f},
{1432.7521999999997f, 514.40519f},
{1404.9476999999997f, 453.19188999999994f},
{1371.6697999999997f, 396.58807999999993f},
{1354.5450999999996f, 324.8727599999999f},
{1332.8455999999996f, 256.0119299999999f},
{1340.0710999999997f, 112.4543599999999f},
{1362.3071999999997f, 91.9626769999999f},
{1363.2939999999996f, 210.0522599999999f},
{1376.7886999999996f, 273.89760999999993f},
{1417.3700999999996f, 345.1889799999999f},
{1462.8065999999997f, 406.9696399999999f},
{1522.8750999999997f, 509.38454999999993f},
{1552.6829999999998f, 609.396f},
{1560.2552999999998f, 702.54805f},
{1542.5469999999998f, 780.01856f},
{1526.7858999999999f, 838.62652f},
{1493.9568f, 914.13715f},
{1447.8102f, 1038.6352f},
{1427.8602999999998f, 1105.8659f},
{1428.6484999999998f, 1181.1731f},
{1482.8379999999997f, 1308.9038f},
{1517.5964999999997f, 1361.8779f},
{1566.5469999999996f, 1411.8881999999999f},
{1612.8160999999996f, 1449.0711f},
{1678.9958999999994f, 1493.7292f},
{1688.0358999999994f, 1502.2943f},
{1736.6581999999994f, 1538.7109f},
{1809.3879999999995f, 1593.3656f},
{1908.2125999999994f, 1654.8400000000001f},
{1979.1966999999993f, 1690.9376000000002f},
{2068.3029999999994f, 1744.9412000000002f},
{2130.1790999999994f, 1766.4037000000003f},
{2182.357499999999f, 1793.5563000000002f},
{2226.550299999999f, 1808.6362000000001f},
{2280.353599999999f, 1847.5087f},
{2331.624699999999f, 1886.3986f},
{2388.593899999999f, 1956.6595f},
{2430.831499999999f, 2022.121f},
{2450.130799999999f, 2095.7699000000002f},
{2435.663499999999f, 2157.0912000000003f},
{2405.0409999999993f, 2226.9551f},
{2380.4411999999993f, 2259.6146000000003f},
{2363.2655999999993f, 2304.4755000000005f},
{2337.266699999999f, 2375.2771000000002f},
{2302.0562999999993f, 2434.5564000000004f},
{2218.222599999999f, 2563.9127000000003f},
{2178.858599999999f, 2637.7700000000004f},
{2162.981099999999f, 2663.1787000000004f},
{2163.673499999999f, 2744.2837000000004f},
{2151.532799999999f, 2821.5692000000004f},
{2177.420999999999f, 2940.1782000000003f},
{2159.365499999999f, 3018.6039f},
{2125.744799999999f, 3090.5215f},
{2093.034199999999f, 3158.0335999999998f},
{2034.8988999999992f, 3182.9936f},
{1922.7770999999993f, 3187.52f},
{1891.1946999999993f, 3140.044f},
{1908.2125999999994f, 2998.2428f},
{1877.1282999999994f, 2921.24f},
{1860.0432999999994f, 2849.3835f},
{1828.4736999999993f, 2721.1673f},
{1834.5588999999993f, 2618.2325f},
{1870.8714999999993f, 2536.7719f},
{1903.6427999999992f, 2485.5849000000003f},
{1933.0935999999992f, 2406.6914f},
{1979.7381999999993f, 2346.4859f},
{2013.5235999999993f, 2293.9896000000003f},
{2066.3618999999994f, 2251.9836000000005f},
{2119.3564999999994f, 2190.6175000000003f},
{2142.8032999999996f, 2132.5034000000005f},
{2144.1536999999994f, 2089.1787000000004f},
{2144.1537f, 2089.1787f}
};

float[][] arm0 = new float[][]{
{85.84822f, 1319.4498f},
{122.94993f, 1374.1015f},
{156.45492f, 1460.8313f},
{196.04480999999998f, 1558.564f},
{216.0084f, 1612.0191f},
{226.3113f, 1679.9102f},
{237.01574f, 1720.6969000000001f},
{261.27719f, 1763.621f},
{287.40491000000003f, 1804.6789f},
{345.26819f, 1844.4376000000002f},
{442.50024f, 1834.5317000000002f},
{495.74199000000004f, 1844.1142000000002f},
{583.1899400000001f, 1840.3095000000003f},
{637.8041800000001f, 1846.0896000000002f},
{643.0702200000001f, 1843.2991000000002f},
{684.1695100000001f, 1861.7696f},
{741.6922500000001f, 1887.7450000000001f},
{787.56411f, 1912.9223000000002f},
{817.4424200000001f, 1928.1191000000001f},
{874.2936000000001f, 1979.1634000000001f},
{915.1853000000001f, 2035.4781f},
{956.52484f, 2084.5659f},
{987.8759100000001f, 2125.3672f},
{1033.9112f, 2179.7983f},
{1081.7274f, 2222.1468f},
{1122.5425f, 2243.3097f},
{1153.6201f, 2259.4013999999997f},
{1232.5324f, 2284.8772999999997f},
{1292.4026000000001f, 2300.8654999999994f},
{1373.5404f, 2309.3412999999996f},
{1422.0944000000002f, 2306.7043999999996f},
{1474.9719000000002f, 2292.4037f},
{1487.5792000000001f, 2293.5681999999997f},
{1534.8739f, 2285.2816999999995f},
{1604.3374000000001f, 2281.6461999999997f},
{1679.6458000000002f, 2271.2528999999995f},
{1769.9632000000001f, 2270.0347999999994f},
{1843.5054000000002f, 2259.9036999999994f},
{1956.4249000000002f, 2287.5120999999995f},
{1986.4758000000002f, 2535.8420999999994f},
{1972.8298000000002f, 2585.3512999999994f},
{1875.1713000000002f, 2640.9569999999994f},
{1826.0583000000001f, 2641.7292999999995f},
{1785.5010000000002f, 2580.9190999999996f},
{1730.6973000000003f, 2546.8269999999998f},
{1675.2345000000003f, 2526.8875f},
{1653.0451000000003f, 2511.1641999999997f},
{1594.6729000000003f, 2492.5114999999996f},
{1485.5037000000002f, 2496.4193999999998f},
{1423.9607000000003f, 2519.4586f},
{1365.2910000000004f, 2529.1299f},
{1250.3980000000004f, 2543.7201f},
{1176.0224000000003f, 2537.0838f},
{1163.0256000000004f, 2535.6589f},
{1102.5358000000003f, 2516.9184999999998f},
{1035.0596000000003f, 2493.7025f},
{974.9400100000003f, 2456.1717f},
{933.7035300000002f, 2415.3635f},
{930.9925800000002f, 2414.6911999999998f},
{865.7710900000002f, 2348.7832f},
{798.4062600000002f, 2257.8329999999996f},
{765.6987200000002f, 2215.7331999999997f},
{743.2693400000002f, 2171.0117999999998f},
{693.3190000000002f, 2091.2771f},
{643.9972500000002f, 2030.3531999999998f},
{603.9973500000002f, 2000.3715999999997f},
{568.7448000000002f, 1980.8084999999996f},
{512.1145900000001f, 1965.6969999999997f},
{432.0317000000001f, 1959.5761999999997f},
{380.4209100000001f, 1953.5447999999997f},
{273.0539900000001f, 1905.9819999999997f},
{232.5430300000001f, 1863.1141999999998f},
{177.0602300000001f, 1780.9826999999998f},
{148.1777400000001f, 1696.5697999999998f},
{140.9649000000001f, 1657.6855999999998f},
{139.87609000000012f, 1581.6256999999998f},
{141.71504000000013f, 1496.4846999999997f},
{127.88078000000013f, 1451.6513999999997f},
{102.47356000000013f, 1383.7150999999997f},
{81.94162400000013f, 1340.1466999999998f},
{52.66278500000013f, 1300.8820999999998f},
{26.05059100000013f, 1251.8832999999997f},
{17.47060900000013f, 1223.4391999999998f},
{85.84822000000014f, 1319.4497999999999f},
{85.84822f, 1319.4498f}
};

float[][] hole_in_arm_8 = new float[][]{
{3057.8763f, 2021.1657f},
{3055.3716999999997f, 2033.4317f},
{3055.3996999999995f, 2045.6158f},
{3057.1749999999993f, 2061.2635f},
{3069.9629999999993f, 2084.0365f},
{3075.5326999999993f, 2090.1854000000003f},
{3078.835899999999f, 2091.8543000000004f},
{3085.8433999999993f, 2097.3260000000005f},
{3085.791399999999f, 2097.7294000000006f},
{3111.9200999999994f, 2113.9742000000006f},
{3115.8522999999996f, 2115.0876000000007f},
{3131.9914999999996f, 2112.5278000000008f},
{3144.8906999999995f, 2111.8692000000005f},
{3169.5273999999995f, 2108.3533000000007f},
{3179.8664999999996f, 2102.027600000001f},
{3186.0398999999998f, 2097.267300000001f},
{3188.9698999999996f, 2092.718700000001f},
{3186.0008999999995f, 2079.9929000000006f},
{3180.0252999999993f, 2065.5576000000005f},
{3166.419799999999f, 2038.0084000000006f},
{3156.819299999999f, 2027.5303000000006f},
{3148.405099999999f, 2021.4347000000005f},
{3129.436399999999f, 2004.7535000000005f},
{3114.1479999999992f, 2001.1575000000005f},
{3087.1801999999993f, 2005.0506000000005f},
{3067.4036999999994f, 2011.4934000000005f},
{3057.8765999999996f, 2021.1657000000005f},
{3057.8763f, 2021.1657f}
};
  public void settings() {  size(750, 750); }
  static public void main(String[] passedArgs) {
    String[] appletArgs = new String[] { "watershed" };
    if (passedArgs != null) {
      PApplet.main(concat(appletArgs, passedArgs));
    } else {
      PApplet.main(appletArgs);
    }
  }
}
