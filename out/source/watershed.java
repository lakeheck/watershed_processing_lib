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
String saveFilePath = "../outputs/watershed-" + new java.text.SimpleDateFormat("yyyyMMdd-HHmmss").format(new java.util.Date()); //change the XXXX to current project name 
int printWidth = 6;
int printHeight = 6;
int printDpi = 300;
int previewDpi = 72;
boolean renderHighRes = false;
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

    ArrayList<PVector> path = inkscapePathImport(p, 3564.00000f, 5014.66650f);
    // printArray(path);
    line = new ArrayList(); //generate a random line 
    int n = 50;
    for(int i=0; i<50; i++){
        line.add(new PVector(map(i, 0, n-1, 0, renderWidth), renderHeight/2 + map(compoundTrigFunction(map(i, 0, n-1, 0, 2*TWO_PI), 0), -3, 4, -50, 50)));
    }

    Ribbon r = new Ribbon(line, renderHighRes ? printDpi/previewDpi * 50 : 50);
    render.beginDraw();

    float[] t = new float[400];
    float scale = 2;
    for(int i=0; i<t.length; i++){
        t[i] = map(i, 0, t.length, -scale*TWO_PI, scale*TWO_PI);
    }
    int sample_size = PApplet.parseInt(t.length*0.5f);
    float rectSize = renderWidth/sample_size; 
    int numRows=20;

    for(int j=0; j < numRows; j++){
        int start = PApplet.parseInt(map(j, 0 , numRows, 0, renderWidth*2/3)); //choose starting point in t for this row
        ArrayList<PVector> tempLine = new ArrayList();

        for(int i=0; i<sample_size*2; i++){
            tempLine.add(new PVector(
                i*rectSize, 
                map(j, 0, numRows, renderHeight*0, renderHeight + randomGaussian()*(renderHighRes ? 10*printDpi/previewDpi : 10)) + (renderHighRes ? 10*printDpi/previewDpi : 10)*compoundTrigFunction(t[PApplet.parseInt((start+i)%t.length)], 0)
            ));
        }

        Ribbon tempRibbon = new Ribbon(tempLine, 20);
        // tempRibbon.vadenWeb(500, 20, new Gradient(line_palette));
        
    }
    render.stroke(0,0,100, 5);
    canvas_overlay_example1();
    render.fill(0);
    render.beginShape();
    for(PVector p: path){
        render.vertex(p.x, p.y); //, 10, 10);
    }
    render.endShape();
    // render.endDraw();

    // Polygon poly = new Polygon(r.vertices, true);
    // poly.subdivide();


    // render.beginDraw();
    // render.background(255);
    // poly.geometricSubdivision.display();


    // r.noFill();
    // r.display();

    
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
      Ribbon r = new Ribbon(path, stroke_weight);
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
    float coeff4 = 4;// random(4,5);
    switch(choice){
      case 0: return cos(x*coeff1+coeff2) - coeff3*sin(randomGaussian()*0.01f + x) + cos(coeff4*x)*pow(sin(pow(x,2)), 2);
      default: return sin(x) + coeff3 * cos(random(2,4)*x);
    }
}


class Ribbon{ //class for drawing a ribbon based on a guide line (as used in flow fields, etc)
    ArrayList<PVector> vertices;

    Ribbon(ArrayList<PVector> guideLine, float stroke_weight){ //should be initialized with ordered set of points 
        vertices = new ArrayList();

        
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

    public void vadenWeb(int n, int _knn, Gradient grad){
      ArrayList<PVector> points = this.generatePointsInside(n);
        Gradient lineGrad = grad;
        for(PVector p:points){
            ArrayList<PVector> knn = k_nearest_neighbors(p, points, _knn);
            for(PVector k:knn){
                    int baseColor = lineGrad.eval(map(k.x,0,renderWidth,0,1)+randomGaussian()*0.1f, HSB);
                    render.stroke(hue(baseColor) + randomGaussian(), saturation(baseColor) + randomGaussian()*8, brightness(baseColor) + randomGaussian()*8);
                    render.line(p.x, p.y, k.x, k.y);
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

    int new_idx = PApplet.parseInt(random(edges.size()-1)); //choose another face randomly (could be improved with some better logic for eligible / preffered faces)
    new_idx = (new_idx>= idx) ? new_idx+1 : new_idx;
    PVector np2 = edges.get(new_idx).getRandomPoint(this.pct, this.std);//select point from new face

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
float[][] p = new float[][]{
{2215.2573f, 3544.0385f},
{2261.5047000000004f, 3548.5978f},
{2290.7342000000003f, 3559.1952f},
{2328.409f, 3573.5438f},
{2382.723f, 3591.8026999999997f},
{2403.0063999999998f, 3601.5912f},
{2436.6746f, 3614.1504999999997f},
{2461.5746999999997f, 3627.9021999999995f},
{2497.0633999999995f, 3644.8167999999996f},
{2526.8517999999995f, 3653.1298999999995f},
{2573.6884999999993f, 3655.3135999999995f},
{2602.6853999999994f, 3639.5784999999996f},
{2638.1402999999996f, 3623.5925999999995f},
{2659.3894999999998f, 3601.4128999999994f},
{2691.3169999999996f, 3578.7973999999995f},
{2733.2663999999995f, 3567.6235999999994f},
{2794.4071999999996f, 3563.5223999999994f},
{2814.7430999999997f, 3540.7060999999994f},
{2837.7985f, 3518.6064999999994f},
{2893.0294f, 3471.456899999999f},
{2924.8514999999998f, 3459.6232999999993f},
{2958.0591999999997f, 3466.7863999999995f},
{2994.7425f, 3469.6443999999997f},
{2998.5045999999998f, 3467.8586999999998f},
{2972.5251999999996f, 3493.3713f},
{2952.4535999999994f, 3522.5946999999996f},
{2931.1429999999996f, 3555.6141999999995f},
{2909.9312999999997f, 3587.5846999999994f},
{2838.3484f, 3608.8631999999993f},
{2778.6634f, 3621.6008999999995f},
{2745.2599f, 3629.9626999999996f},
{2717.7443f, 3641.0411999999997f},
{2688.7536f, 3659.8909999999996f},
{2686.454f, 3675.1917999999996f},
{2717.7249f, 3681.4263999999994f},
{2746.2166f, 3686.7770999999993f},
{2782.4218f, 3695.483199999999f},
{2852.7963f, 3706.1207999999992f},
{2889.2325f, 3710.183499999999f},
{3014.8374f, 3714.985799999999f},
{3051.798f, 3719.343099999999f},
{3082.5386999999996f, 3723.549899999999f},
{3109.7366999999995f, 3741.206799999999f},
{3131.5937999999996f, 3756.792799999999f},
{3162.5689999999995f, 3767.5370999999986f},
{3205.9403999999995f, 3766.8074999999985f},
{3244.4909999999995f, 3729.2212999999983f},
{3282.7612999999997f, 3708.2698999999984f},
{3316.8664999999996f, 3667.4650999999985f},
{3344.6108999999997f, 3634.5044999999986f},
{3389.6215999999995f, 3604.9561999999987f},
{3397.2535999999996f, 3624.344799999999f},
{3367.4373999999993f, 3655.098699999999f},
{3347.6062999999995f, 3678.383999999999f},
{3331.6476999999995f, 3701.211299999999f},
{3310.9209999999994f, 3723.1877999999992f},
{3283.7970999999993f, 3751.023799999999f},
{3247.6023999999993f, 3781.941099999999f},
{3203.5947999999994f, 3811.141899999999f},
{3222.3443999999995f, 3833.079799999999f},
{3250.6821999999993f, 3860.240099999999f},
{3300.419999999999f, 3916.907499999999f},
{3339.869799999999f, 3934.326099999999f},
{3394.932799999999f, 3956.6461999999988f},
{3428.548099999999f, 3979.636399999999f},
{3457.693899999999f, 3996.5220999999988f},
{3484.871699999999f, 4037.4178999999986f},
{3488.147399999999f, 4052.1248999999984f},
{3496.3329999999987f, 4081.0844999999986f},
{3490.577199999999f, 4090.1926999999987f},
{3455.144699999999f, 4044.901199999999f},
{3425.1135999999988f, 4016.880399999999f},
{3359.9689999999987f, 3989.854899999999f},
{3288.4972999999986f, 3963.768199999999f},
{3221.5112999999988f, 3904.639899999999f},
{3195.8158999999987f, 3863.8055999999992f},
{3172.651599999999f, 3829.577099999999f},
{3144.6703999999986f, 3810.019999999999f},
{3105.3445999999985f, 3791.639499999999f},
{3071.4370999999987f, 3773.638899999999f},
{3033.774399999999f, 3763.667099999999f},
{3000.538799999999f, 3762.776099999999f},
{2948.6996999999988f, 3777.321699999999f},
{2919.112299999999f, 3789.475699999999f},
{2875.098399999999f, 3798.1776999999993f},
{2834.864899999999f, 3804.6244999999994f},
{2795.8905999999993f, 3795.6187999999993f},
{2754.140499999999f, 3758.6257999999993f},
{2674.043799999999f, 3715.805799999999f},
{2638.656099999999f, 3719.3226999999993f},
{2587.4529999999986f, 3742.0196999999994f},
{2554.8114999999984f, 3757.7538999999992f},
{2509.6720999999984f, 3773.3116999999993f},
{2428.6887999999985f, 3778.353199999999f},
{2365.9869999999987f, 3771.638699999999f},
{2338.3954999999987f, 3765.0998999999993f},
{2303.0736999999986f, 3741.2954999999993f},
{2248.8500999999987f, 3717.601199999999f},
{2169.3713999999986f, 3696.605499999999f},
{2128.722899999999f, 3681.117799999999f},
{2090.557299999999f, 3620.469699999999f},
{2061.638199999999f, 3578.5747999999994f},
{2032.363299999999f, 3527.2420999999995f},
{2055.140899999999f, 3498.7293999999993f},
{2144.015599999999f, 3507.350999999999f},
{2181.649699999999f, 3526.8323999999993f},
{2212.391199999999f, 3543.963399999999f},
{2213.7662999999993f, 3544.019399999999f},
{2215.2572999999993f, 3544.038399999999f},
{2215.2573f, 3544.0385f}
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
