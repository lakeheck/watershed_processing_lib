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
        branch5, branch4, branch3, branch2, branch1, arm0, arm1_1, arm1_2, arm1_3, arm3, arm4_1, arm4_2, arm4_3, arm4_4, arm5, arm6, arm7, arm8
    };

    float[][][] import_paths_highlight = new float[][][]{
        branch5_highlight, branch4_highlight, branch3_highlight, branch2_highlight, branch1_highlight, arm0_highlight, arm1_1_highlight, arm1_2_highlight, arm1_3_highlight, arm3_highlight, arm4_1_highlight, arm4_2_highlight, arm4_3_highlight, arm4_4_highlight, arm5_highlight, arm6_highlight, arm7_highlight, arm8_highlight
    };
    ArrayList<ArrayList<PVector>> paths = new ArrayList<ArrayList<PVector>>(); 
    ArrayList<ArrayList<PVector>> highlight_paths = new ArrayList<ArrayList<PVector>>(); 

    for(float[][] p:import_paths){
        paths.add(inkscapePathImport(p, 3564.00000f, 5014.66650f));
    }
    for(float[][] p:import_paths_highlight){
        highlight_paths.add(inkscapePathImport(p, 3564.00000f, 5014.66650f));
    }

    render.beginDraw();
    render.background(255);
    for(int i=0; i<paths.size(); i++){
        ArrayList<PVector> tempPoints = new ArrayList();
        Region r = new Region(paths.get(i), renderHighRes ? printDpi/previewDpi * 50 : 50, true);
        if(highlight_paths.get(i).size() > 0){
            Region rh = new Region(highlight_paths.get(i), renderHighRes ? printDpi/previewDpi * (randomGaussian()+7) : randomGaussian()+7, false); //created highlight paths to add
            tempPoints = rh.generatePointsInside(150);
        }
        r.vadenWeb(300, 10, new Gradient(line_palette), i<=4 ? true : false, tempPoints);
    }

    // Region rh = new Region(highlight_paths.get(0), renderHighRes ? printDpi/previewDpi * 10 : 10, false);
    // Region r = new Region(paths.get(7), renderHighRes ? printDpi/previewDpi * 50 : 50, true);
    
    // Polygon p = new Polygon(paths.get(7), true);
    // p.subdivide();
    // p.geometricSubdivision.display();

    // rh.vadenWeb(300, 10, new Gradient(line_palette), false);

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
    //     Region r = new Region(paths.get(i), renderHighRes ? printDpi/previewDpi * 50 : 50, true);
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
      Region r = new Region(path, stroke_weight, false);
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


class Region{ //class for drawing a Region based on a guide line (as used in flow fields, etc)
    ArrayList<PVector> vertices;

    Region(ArrayList<PVector> guideLine, float stroke_weight, boolean closed){ //should be initialized with ordered set of points 
        // closed bool value set to false if the input line is an open loop (and must be "expanded" into a region)
        //close bool should be set to true if the input list of points forms a closed loop
        vertices = new ArrayList();
        if(closed){
          vertices.addAll(guideLine);
        }
        else{
          vertices.addAll(generateBoundsfromFrameLine(guideLine, stroke_weight));
        }
        
    }

    public ArrayList<PVector> generateBoundsfromFrameLine(ArrayList<PVector> guideLine, float stroke_weight){
          ArrayList<PVector> out = new ArrayList();
          ArrayList<PVector> top = new ArrayList();
          ArrayList<PVector> btm = new ArrayList();
          for(int i=1; i<guideLine.size(); i++){
              PVector norm = new PVector(0,0,1).cross(new PVector((guideLine.get(i).x - guideLine.get(i-1).x), (guideLine.get(i).y - guideLine.get(i-1).y))).normalize();
              float theta = norm.heading();
              top.add(new PVector(guideLine.get(i).x + stroke_weight*cos(theta), guideLine.get(i).y + stroke_weight*sin(theta)));
              btm.add(new PVector(guideLine.get(i).x - stroke_weight*cos(theta), guideLine.get(i).y - stroke_weight*sin(theta)));
          }
          for(int i=0; i<((top.size() + btm.size())); i++){ // unwrap the top and bottom arrays - first we add all the top points, then start fro the end of hte bottom array to maintain non-self intersection
              out.add(
                  i<top.size() ? top.get(i).copy() : btm.get(btm.size()-1-(i-top.size())).copy()
              );
          }
          return out;
    }

    public boolean contains(PVector point){
        return polyPoint(vertices, point.x, point.y);
    }

    public ArrayList<PVector> generatePointsInside(int n){
        ArrayList<PVector> points = new ArrayList();
        int count = 0;
        while(count <= n){
            PVector p = new PVector(random(renderWidth), random(renderHeight));
            // println(inkscapeRecoverInitialValue(p, img));
            if(polyPoint(this.vertices, p.x, p.y)){
                points.add(p);
                count++;
            }
        }
        return points;
    }

    public float regionArea(){
        float area = 0.0f;
        int n = vertices.size();
        // Calculate value of shoelace formula
        int j = n - 1;
        for (int i = 0; i < n; i++){
          area += (vertices.get(j).x + vertices.get(i).x) * (vertices.get(j).y - vertices.get(i).y);
            
          // j is previous vertex to i
          j = i;
        }
        return Math.abs(area / 2.0f);

    }

    public ArrayList<PVector> generatePointsInside(ArrayList<PVector> region, int n){
        ArrayList<PVector> points = new ArrayList();
        int count = 0;
        while(count <= n){
            PVector p = new PVector(random(renderWidth), random(renderHeight));
            // println(inkscapeRecoverInitialValue(p, img));
            if(polyPoint(region, p.x, p.y)){
                points.add(p);
                count++;
            }
        }
        return points;
    }

    public ArrayList<PVector> generate_N_points_about_curve(ArrayList<PVector> line, int n, float var){
        ArrayList<PVector> points = new ArrayList();
        int count = 0;
        ArrayList<PVector> region = generateBoundsfromFrameLine(line, 20);
        points.addAll(this.generatePointsInside(region, n));
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

    public void vadenWeb(int n, int _knn, Gradient grad, boolean allow_intersection, ArrayList<PVector> addl_points){
      ArrayList<PVector> points = this.generatePointsInside(n);
      points.addAll(addl_points);
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
    float area = 0.0f;
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
{1871.2624f, 2919.0638f},
{1901.1657f, 3071.1031f},
{1887.3396f, 3154.4330999999997f},
{1918.6283f, 3281.4424f},
{1933.1080000000002f, 3353.6258f},
{2038.5892000000001f, 3555.8673999999996f},
{2170.4837f, 3713.7445f},
{2390.8612000000003f, 3785.5996999999998f},
{2581.1677000000004f, 3729.4687999999996f},
{2634.3890000000006f, 3718.8174999999997f},
{2770.2894000000006f, 3780.0734999999995f},
{2805.3427000000006f, 3821.6109999999994f},
{2796.1611000000007f, 3924.5390999999995f},
{2804.7993000000006f, 3909.6137999999996f},
{2856.3828000000008f, 3813.2387999999996f},
{2923.287500000001f, 3804.3475999999996f},
{2975.398900000001f, 3788.2232999999997f},
{3019.8536000000013f, 3767.5350999999996f},
{3080.0207000000014f, 3789.5039999999995f},
{3235.4610000000016f, 3910.9649999999992f},
{3390.3159000000014f, 4009.4589999999994f},
{3437.3132000000014f, 4024.1525999999994f},
{3491.7357000000015f, 4100.850699999999f},
{3507.4102000000016f, 4097.289499999999f},
{3440.723000000002f, 3999.262899999999f},
{3360.1912000000016f, 3946.8302999999987f},
{3209.3866000000016f, 3826.9823999999985f},
{3161.8537000000015f, 3794.5973999999983f},
{3203.4176000000016f, 3797.4950999999983f},
{3276.4620000000014f, 3761.3660999999984f},
{3353.6770000000015f, 3660.2503999999985f},
{3394.6806000000015f, 3642.1345999999985f},
{3406.9661000000015f, 3596.7507999999984f},
{3366.8151000000016f, 3600.5821999999985f},
{3304.1590000000015f, 3676.5876999999987f},
{3253.0682000000015f, 3745.2598999999987f},
{3235.7886000000017f, 3741.9326999999985f},
{3262.7814000000017f, 3612.9409999999984f},
{3257.322200000002f, 3608.6027999999983f},
{3191.028100000002f, 3742.850499999998f},
{3150.209400000002f, 3754.503799999998f},
{3006.692300000002f, 3722.586299999998f},
{2941.572700000002f, 3722.329499999998f},
{2773.007700000002f, 3692.472099999998f},
{2718.554900000002f, 3674.975899999998f},
{2663.677300000002f, 3674.070799999998f},
{2685.502600000002f, 3645.7223999999983f},
{2796.7363000000023f, 3620.2992999999983f},
{2898.0637000000024f, 3595.540199999998f},
{2937.0070000000023f, 3532.499099999998f},
{2957.7914000000023f, 3489.777799999998f},
{2929.3831000000023f, 3530.218799999998f},
{2853.1106000000023f, 3590.6660999999976f},
{2822.7078000000024f, 3575.2161999999976f},
{2829.5760000000023f, 3533.5237999999977f},
{2892.9605000000024f, 3481.6840999999977f},
{2863.8654000000024f, 3475.8521999999975f},
{2819.8392000000026f, 3514.2489999999975f},
{2788.7717000000025f, 3578.0083999999974f},
{2745.6762000000026f, 3575.7332999999976f},
{2686.1397000000024f, 3596.3020999999976f},
{2645.3470000000025f, 3613.0690999999974f},
{2645.2199000000023f, 3574.5611999999974f},
{2737.025100000002f, 3514.4293999999973f},
{2737.964100000002f, 3499.5036999999975f},
{2640.1180000000018f, 3538.3813999999975f},
{2598.769100000002f, 3602.8196999999973f},
{2520.9894000000018f, 3641.286599999997f},
{2415.2742000000017f, 3628.6087999999972f},
{2236.7567000000017f, 3572.6955999999973f},
{2183.0244000000016f, 3512.2472999999973f},
{2113.485600000002f, 3338.2106999999974f},
{2064.203100000002f, 3204.8667999999975f},
{2041.079900000002f, 3102.9425999999976f},
{1993.875100000002f, 3035.9537999999975f},
{1918.704700000002f, 2941.2591999999977f},
{1876.541000000002f, 2932.2602999999976f}
};


float[][] branch4 = new float[][]{
{2137.8239f, 3753.1378f},
{2285.0020999999997f, 4012.5571f},
{2159.0263999999997f, 4173.7859f},
{2167.9673999999995f, 4400.9556f},
{2198.4863999999993f, 4434.050200000001f},
{2169.049699999999f, 4246.3612f},
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
{3005.5858f, 4989.4783f},
{3024.2711f, 5006.3451f},
{2901.6201f, 4807.282499999999f},
{2910.2061f, 4593.0335f},
{3108.3214f, 4598.3911f},
{3163.9809999999998f, 4817.9164f},
{3257.4327999999996f, 4973.3536f},
{3428.8215999999998f, 5042.9619f},
{3339.2075999999997f, 4971.6224f},
{3204.5766999999996f, 4855.0244f},
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
{2234.3904999999986f, 3742.504499999999f},
{2137.8238999999985f, 3753.137799999999f},
{2137.8239f, 3753.1378f}
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
{1163.9305f, 3774.1964f},
{1134.7179999999998f, 3795.8958f},
{1091.2757f, 3817.3014f},
{1062.3817f, 3820.1209999999996f},
{1027.1248999999998f, 3827.4837999999995f},
{967.4855599999999f, 3834.6287999999995f},
{931.3707599999998f, 3853.1715999999997f},
{914.0526399999998f, 3873.0955999999996f},
{897.3613799999998f, 3900.8826999999997f},
{917.3222199999998f, 3903.16f},
{966.8745599999997f, 3867.6034999999997f},
{1052.2356999999997f, 3861.0721999999996f},
{1097.9479999999996f, 3866.5717999999997f},
{993.0859299999996f, 3884.4172f},
{942.2294499999996f, 3911.4399f},
{925.9819499999996f, 3961.9831999999997f},
{887.9775399999996f, 3972.1679999999997f},
{851.1262499999997f, 3982.2767999999996f},
{844.2602399999997f, 4030.2014999999997f},
{842.2884399999997f, 4066.6388999999995f},
{826.1003299999996f, 4143.6981f},
{832.4873999999996f, 4123.2806f},
{836.1846199999997f, 4119.6531f},
{847.2921599999996f, 4091.1500000000005f},
{880.2619599999996f, 4034.1906000000004f},
{916.2763999999996f, 4015.8682000000003f},
{953.6223299999996f, 3994.4130000000005f},
{998.2613099999995f, 3943.6886000000004f},
{1025.1473999999996f, 3927.6967000000004f},
{1102.7035999999996f, 3913.2831000000006f},
{1144.9274999999996f, 3912.0982000000004f},
{1182.1681999999996f, 3891.6804f},
{1225.1114999999995f, 3850.0343000000003f},
{1245.5813999999996f, 3814.9029f},
{1267.5544999999995f, 3790.9709000000003f},
{1307.5756999999994f, 3781.0787000000005f},
{1379.6054999999994f, 3768.9336000000003f},
{1429.9800999999995f, 3763.3738000000003f},
{1419.7449999999994f, 3790.2403000000004f},
{1373.9168999999995f, 3825.6399000000006f},
{1305.5792999999994f, 3869.8042000000005f},
{1267.0163999999995f, 3910.8515000000007f},
{1252.8934999999994f, 3941.571900000001f},
{1254.8370999999995f, 3950.663500000001f},
{1313.3142999999995f, 3922.644400000001f},
{1361.4319999999996f, 3882.473400000001f},
{1414.9310999999996f, 3842.408200000001f},
{1439.8652999999995f, 3821.7279000000012f},
{1487.4216999999994f, 3787.0789000000013f},
{1554.8782999999994f, 3804.2433000000015f},
{1612.4691999999993f, 3820.3445000000015f},
{1683.0473999999992f, 3833.2166000000016f},
{1732.5591999999992f, 3840.1260000000016f},
{1812.4153999999992f, 3818.287300000002f},
{1872.6200999999992f, 3774.9484000000016f},
{1980.727299999999f, 3712.4639000000016f},
{2043.4780999999991f, 3641.4803000000015f},
{2030.6472999999992f, 3586.2894000000015f},
{2003.8427999999992f, 3507.7986000000014f},
{1977.6482999999992f, 3497.9372000000017f},
{1977.1571999999992f, 3529.657400000002f},
{1975.1271999999992f, 3600.5746000000017f},
{1954.0527999999993f, 3618.0770000000016f},
{1924.6732999999992f, 3623.764600000002f},
{1882.9718999999993f, 3626.664400000002f},
{1838.5028999999993f, 3641.372000000002f},
{1766.9153999999992f, 3641.959200000002f},
{1716.871399999999f, 3629.329500000002f},
{1684.349799999999f, 3628.189700000002f},
{1649.5613999999991f, 3668.624500000002f},
{1673.6452999999992f, 3668.613500000002f},
{1735.9324999999992f, 3669.502200000002f},
{1881.8195999999991f, 3668.624400000002f},
{1846.3410999999992f, 3685.099300000002f},
{1820.4806999999992f, 3697.463500000002f},
{1769.9880999999991f, 3709.857400000002f},
{1745.351599999999f, 3720.710100000002f},
{1626.578999999999f, 3723.621700000002f},
{1476.655499999999f, 3688.314100000002f},
{1445.407899999999f, 3679.630200000002f},
{1356.599299999999f, 3668.624400000002f},
{1280.0905999999989f, 3664.668000000002f},
{1196.921699999999f, 3642.231500000002f},
{1162.6539999999989f, 3654.3161000000023f},
{1157.332199999999f, 3655.4279000000024f},
{1185.850999999999f, 3660.7698000000023f},
{1234.565299999999f, 3670.845300000002f},
{1298.637699999999f, 3698.5668000000023f},
{1346.6038999999992f, 3696.248400000002f},
{1330.791799999999f, 3712.532100000002f},
{1285.504699999999f, 3725.762000000002f},
{1147.475799999999f, 3726.322300000002f},
{882.8595099999991f, 3635.471500000002f},
{833.5937499999991f, 3613.4341000000018f},
{854.7583999999991f, 3632.5258000000017f},
{993.4873899999991f, 3696.5981000000015f},
{1148.868599999999f, 3763.3384000000015f},
{1163.930499999999f, 3774.1963000000014f},
{1163.9305f, 3774.1964f}
};

float[][] arm1_1 = new float [][]{
{1510.6324f, 2241.7276f},
{1724.5534f, 2277.8263f},
{1964.1244f, 2273.5278000000003f},
{1959.6400999999998f, 2134.2670000000003f},
{1857.0216999999998f, 1978.4273000000003f},
{1678.5158999999999f, 1888.5941000000003f},
{1474.3716f, 1821.0199000000002f},
{1243.1662999999999f, 1813.1895000000002f},
{856.9986499999999f, 1958.42f},
{981.0040699999998f, 2103.1317f},
{1355.2021f, 2111.8485f},
{1510.6324f, 2241.7276f},
{1510.6324f, 2241.7276f}
};

float[][] arm1_2 = new float [][]{
{609.17875f, 2007.6571f},
{502.81016000000005f, 1964.4723999999999f},
{308.04271000000006f, 1933.4711f},
{403.15022000000005f, 2136.9716f},
{586.2718100000001f, 2221.5987999999998f},
{752.73766f, 2204.1346f},
{628.07143f, 2003.1061f},
{609.17875f, 2007.6571f}
};

float[][] arm1_3 = new float [][]{
{519.16929f, 1585.5531f},
{394.13460000000003f, 1459.9435f},
{380.82219000000003f, 1301.7786f},
{381.61889f, 1299.3772000000001f},
{382.58446000000004f, 1297.0546000000002f},
{347.47272000000004f, 1470.0165000000002f},
{421.54468f, 1619.6537000000003f},
{321.00487000000004f, 1832.1015000000002f},
{497.48916f, 1851.1015000000002f},
{519.16929f, 1585.5531f},
{519.16929f, 1585.5531f}
};

float[][] arm8 = new float[][]{
{2300.9992f, 2529.6542f},
{2399.5678000000003f, 2504.7842f},
{2559.9898000000003f, 2540.9824f},
{2727.7252000000003f, 2563.9755f},
{2945.101f, 2608.4653f},
{3058.3697f, 2582.8660999999997f},
{3131.3138000000004f, 2562.2897999999996f},
{3275.6252000000004f, 2445.9056999999993f},
{3340.9911f, 2267.2026999999994f},
{3338.509f, 2192.8591999999994f},
{3315.6661f, 2073.257699999999f},
{3338.1353f, 2031.969199999999f},
{3407.4593999999997f, 1853.0375999999992f},
{3492.1854f, 1794.494499999999f},
{3535.2165999999997f, 1773.840599999999f},
{3409.3929f, 1812.6164999999992f},
{3379.1315999999997f, 1855.807999999999f},
{3286.3309999999997f, 2022.0983999999992f},
{3254.5831f, 2011.2451999999992f},
{3152.5946f, 1930.3465999999992f},
{3105.2774999999997f, 1914.1588999999992f},
{2988.0539f, 1933.5157999999992f},
{2940.5702f, 1993.5761999999993f},
{2957.6085000000003f, 2123.845999999999f},
{3074.8776000000003f, 2213.676599999999f},
{3163.3343000000004f, 2210.569999999999f},
{3212.5901000000003f, 2187.1440999999986f},
{3204.4821f, 2215.4619999999986f},
{3120.3299f, 2313.8403999999987f},
{3033.0044000000003f, 2331.669599999999f},
{2820.9338000000002f, 2329.734199999999f},
{2593.0995000000003f, 2287.318999999999f},
{2414.2473000000005f, 2255.1776999999993f},
{2379.6294000000003f, 2257.785799999999f},
{2300.2048000000004f, 2265.686199999999f},
{2236.7491000000005f, 2316.606499999999f},
{2154.2135000000003f, 2466.425599999999f},
{2155.9172000000003f, 2685.3498999999993f},
{2162.5381f, 2718.345099999999f},
{2192.6356f, 2682.500899999999f},
{2293.9333f, 2537.2356999999993f},
{2300.9992f, 2529.6541999999995f},
{2300.9992f, 2529.6542f}
};

float[][] arm7 = new float[][]{
{3471.2541f, 936.86536f},
{3301.1f, 1088.1846f},
{3277.6418f, 1144.1152f},
{3301.9058999999997f, 1313.636f},
{3270.1537999999996f, 1386.5021f},
{3170.1545999999994f, 1476.6858f},
{3093.4939999999992f, 1483.6796f},
{2977.631299999999f, 1526.5184f},
{2897.536999999999f, 1600.0748999999998f},
{2784.7967999999987f, 1740.7341f},
{2607.5853999999986f, 1837.7274f},
{2414.7406999999985f, 1907.6526f},
{2231.9617999999987f, 1987.4942999999998f},
{2213.225799999999f, 2019.0286999999998f},
{2225.718399999999f, 2135.7468f},
{2245.175699999999f, 2266.801f},
{2285.433399999999f, 2258.3276f},
{2372.940299999999f, 2241.5179000000003f},
{2399.599999999999f, 2262.1149000000005f},
{2401.888099999999f, 2223.6247000000003f},
{2484.8046999999992f, 2119.5727f},
{2666.095099999999f, 2059.3085f},
{2872.682899999999f, 1981.0288f},
{2997.290399999999f, 1826.1976f},
{3122.145999999999f, 1670.9501f},
{3290.495399999999f, 1559.0941f},
{3334.8810999999987f, 1521.9686f},
{3410.4398999999985f, 1333.5666999999999f},
{3326.7009999999987f, 1191.2252999999998f},
{3352.578799999999f, 1077.3709f},
{3470.5141999999987f, 975.9119999999999f},
{3499.4017999999987f, 879.5920299999999f},
{3479.346999999999f, 790.9881199999999f},
{3470.230199999999f, 791.7888699999999f},
{3474.192799999999f, 904.9027999999998f},
{3471.2540999999987f, 936.8653599999999f},
{3471.2541f, 936.86536f}
};

float[][] arm6 = new float[][]{
{2977.2074f, 886.1992f},
{2984.2437999999997f, 1059.7493f},
{2956.4273f, 1124.7498f},
{2859.2237f, 1310.6052f},
{2757.4004f, 1440.5815f},
{2696.0602f, 1515.0074f},
{2484.1877f, 1695.8901f},
{2339.942f, 1766.9962f},
{2307.7631f, 1829.7568f},
{2343.3161f, 1906.6868000000002f},
{2323.6685f, 1946.4819000000002f},
{2231.1690000000003f, 1973.2378000000003f},
{2199.1721000000002f, 1951.7591000000004f},
{2164.1784000000002f, 1857.5067000000004f},
{2143.4542f, 1815.6225000000004f},
{2146.7181f, 1666.7449000000004f},
{2264.536f, 1560.3625000000004f},
{2437.147f, 1478.4450000000004f},
{2659.0169f, 1348.4470000000003f},
{2799.9789f, 1203.4160000000004f},
{2826.3695000000002f, 1159.5155000000004f},
{2892.7052000000003f, 996.0605900000004f},
{2897.2533000000003f, 824.5365600000003f},
{2888.5456000000004f, 749.9654100000004f},
{2818.9998000000005f, 644.1005400000004f},
{2753.8021000000003f, 565.1447900000004f},
{2699.5996000000005f, 387.8703800000004f},
{2751.5880000000006f, 290.2984700000004f},
{2834.8032000000007f, 207.1162900000004f},
{2844.122100000001f, 209.1725100000004f},
{2784.939400000001f, 274.5008600000004f},
{2732.2459000000013f, 404.8588800000004f},
{2744.711000000001f, 479.7540700000004f},
{2806.794800000001f, 572.9782900000005f},
{2887.100700000001f, 672.1687400000005f},
{2957.468800000001f, 782.7388200000005f},
{2977.207400000001f, 886.1992000000005f},
{2977.2074f, 886.1992f}
};

float[][] arm5 = new float[][]{
{2042.9607f, 444.15214f},
{1941.7340000000002f, 330.23668999999995f},
{1940.7524f, 314.73488999999995f},
{2095.4726f, 482.01919999999996f},
{2087.5766f, 768.88398f},
{1957.1964999999998f, 1038.7096999999999f},
{1978.1909999999998f, 1367.3640999999998f},
{1996.2301999999997f, 1414.7474999999997f},
{2075.0427999999997f, 1616.3768999999998f},
{2142.9121999999998f, 1817.5188999999998f},
{2162.0694f, 1882.7894f},
{2208.5980999999997f, 2070.8581999999997f},
{2232.3608f, 2151.8079f},
{2218.1593f, 2303.0301999999997f},
{2188.8165f, 2380.2220999999995f},
{2133.4244f, 2544.8163999999997f},
{2048.8088f, 2576.4410999999996f},
{1972.1871999999998f, 2502.9515999999994f},
{1959.9209999999998f, 2375.0422999999996f},
{1930.8312999999998f, 2112.9003999999995f},
{1853.7880999999998f, 1891.4920999999995f},
{1779.0539999999999f, 1611.3322999999996f},
{1728.6715f, 1297.2455999999995f},
{1856.9664f, 986.7398499999995f},
{2001.4018f, 775.1330399999995f},
{2080.6712f, 568.9063199999995f},
{2042.9607000000003f, 444.15213999999946f},
{2042.9607f, 444.15214f}
};

float[][] arm4_1 = new float[][]{
{1380.358f, 1574.4467f},
{1346.4216999999999f, 1609.2798f},
{1322.9515f, 1656.4834f},
{1313.3196999999998f, 1692.2393000000002f},
{1310.0936999999997f, 1726.1014000000002f},
{1310.2934999999998f, 1769.5288000000003f},
{1313.9993999999997f, 1793.3172000000002f},
{1319.3059999999996f, 1804.2033000000001f},
{1331.9742999999996f, 1805.0310000000002f},
{1357.5744999999997f, 1810.1718f},
{1371.1778999999997f, 1812.1826f},
{1396.9805999999996f, 1813.9293f},
{1429.5034999999996f, 1815.1858f},
{1453.6430999999995f, 1821.463f},
{1474.2826999999995f, 1826.4039f},
{1505.8015999999996f, 1838.3628f},
{1544.4182999999996f, 1853.7184000000002f},
{1573.3326999999995f, 1868.9693000000002f},
{1614.5904999999996f, 1874.6530000000002f},
{1642.9478999999997f, 1883.0260000000003f},
{1659.6675999999998f, 1887.6620000000003f},
{1683.7903999999999f, 1901.3780000000002f},
{1689.7205999999999f, 1898.7197f},
{1693.9941999999999f, 1888.7018f},
{1692.2385f, 1872.3705f},
{1684.9878999999999f, 1843.9007f},
{1668.1009f, 1802.1317f},
{1653.8163f, 1768.5091f},
{1638.3099f, 1736.682f},
{1612.3779f, 1723.9335f},
{1590.4824999999998f, 1719.6717f},
{1531.4527999999998f, 1675.2945000000002f},
{1516.7696999999998f, 1665.5986000000003f},
{1480.6348999999998f, 1638.6207000000002f},
{1432.1849999999997f, 1610.4494000000002f},
{1410.4498999999996f, 1600.0563000000002f},
{1401.5132999999996f, 1582.6140000000003f},
{1394.3653999999997f, 1571.0254000000002f},
{1380.3579999999997f, 1574.4467000000002f},
{1380.358f, 1574.4467f}
};

float[][] arm4_2 = new float[][]{
{1532.7873f, 821.50681f},
{1650.2059f, 952.52095f},
{1690.3732f, 1033.5327f},
{1729.3362f, 1114.1217f},
{1752.5846f, 1134.1532f},
{1748.7658999999999f, 1218.3014f},
{1725.3329999999999f, 1328.0924f},
{1706.3082f, 1461.5073f},
{1699.3727f, 1485.8397f},
{1679.1074999999998f, 1511.5289f},
{1645.2104f, 1503.0732f},
{1567.6057999999998f, 1442.1545f},
{1555.7855999999997f, 1433.7097f},
{1505.7780999999998f, 1377.2626f},
{1500.2072999999998f, 1336.721f},
{1525.6625f, 1308.2525f},
{1571.616f, 1112.0259f},
{1567.4095f, 1077.3582000000001f},
{1528.9094f, 1004.3354000000002f},
{1494.8886f, 934.0219200000001f},
{1502.2818f, 883.8232100000001f},
{1532.7873f, 821.5068100000001f},
{1532.7873f, 821.50681f}
};

float[][] arm4_3 = new float[][]{
{1364.2379f, 761.33487f},
{1314.4332000000002f, 656.20919f},
{1307.3264000000001f, 594.82716f},
{1344.1997000000001f, 477.35871000000003f},
{1360.2973000000002f, 453.23979f},
{1372.0977000000003f, 415.22504000000004f},
{1390.0265000000002f, 384.20369000000005f},
{1399.4395000000002f, 379.08022000000005f},
{1404.6342000000002f, 388.68708000000004f},
{1405.5480000000002f, 409.97451f},
{1378.5617000000002f, 498.21079f},
{1389.9547000000002f, 608.40748f},
{1444.6125000000002f, 698.8097399999999f},
{1451.5767f, 749.6107599999999f},
{1439.0126f, 776.1750199999999f},
{1420.9883f, 816.3754399999999f},
{1396.3709f, 822.30401f},
{1378.7024f, 800.0427f},
{1375.1971999999998f, 787.8430099999999f},
{1364.2378999999999f, 761.3348699999999f},
{1364.2379f, 761.33487f}
};

float[][] arm4_4 = new float[][]{
{1436.0914f, 337.32752f},
{1480.1426000000001f, 190.63431f},
{1486.4876000000002f, 126.5111f},
{1483.4293000000002f, 101.30103f},
{1471.4150000000002f, 83.39336899999999f},
{1465.8227000000002f, 82.100111f},
{1466.8606000000002f, 90.595934f},
{1474.1971f, 108.08773f},
{1468.4005000000002f, 186.87892f},
{1423.2029000000002f, 299.46002999999996f},
{1411.3562000000002f, 330.72158999999994f},
{1421.1656000000003f, 351.54121999999995f},
{1431.2596000000003f, 353.3205899999999f},
{1435.9538000000002f, 342.71127999999993f},
{1436.0914000000002f, 337.32751999999994f},
{1436.0914f, 337.32752f}
};

float[][] arm3 = new float[][]{
{1771.5749f, 1797.4881f},
{1576.0361f, 1692.3267f},
{1403.6544000000001f, 1559.1029f},
{1260.2264f, 1427.0197f},
{1201.158f, 1275.1508000000001f},
{1265.5273f, 998.7686500000001f},
{1403.1676f, 822.7786600000001f},
{1461.2669f, 644.3049000000001f},
{1405.4579f, 450.3475900000001f},
{1358.4381f, 73.1661850000001f},
{1378.1859000000002f, 277.88062000000014f},
{1465.0075000000002f, 409.76912000000016f},
{1550.4206000000001f, 587.2010300000002f},
{1541.3472000000002f, 787.9235400000002f},
{1461.6291f, 1003.3881000000002f},
{1429.6233000000002f, 1185.5660000000003f},
{1506.8209000000002f, 1351.4977000000003f},
{1633.4533000000001f, 1466.2577000000003f},
{1742.9677000000001f, 1564.8510000000003f},
{1815.7716f, 1720.1455000000003f},
{1832.3163f, 1827.4903000000004f},
{1771.5749f, 1797.4881000000005f},
{1771.5749f, 1797.4881f}
};

float[][] arm0 = new float[][]{
{1991.6139f, 2520.9184f},
{2060.2612f, 2571.2988f},
{2157.4526f, 2515.3169f},
{2154.4605f, 2597.2720999999997f},
{2162.2392f, 2725.1413999999995f},
{2159.6861f, 2780.9886999999994f},
{2159.3012f, 2835.7123999999994f},
{2181.27f, 2914.5676999999996f},
{2184.3431f, 2979.3994f},
{2163.1951f, 3052.2738f},
{2102.8134f, 3234.3213f},
{2077.3968f, 3181.4596f},
{2025.2337f, 3038.2805000000003f},
{1924.4349f, 2925.6935000000003f},
{1888.6607999999999f, 2930.0371000000005f},
{1869.9981999999998f, 2885.2467000000006f},
{1839.8426999999997f, 2803.1486000000004f},
{1833.3879999999997f, 2677.4714000000004f},
{1769.5340999999996f, 2574.5214000000005f},
{1676.2074999999995f, 2524.1575000000007f},
{1534.8182999999995f, 2490.698900000001f},
{1310.8853999999994f, 2539.067000000001f},
{1121.7695999999994f, 2526.240400000001f},
{964.1817699999995f, 2443.963600000001f},
{812.4504299999994f, 2274.996900000001f},
{672.1194499999995f, 2065.987000000001f},
{484.7815399999995f, 1962.546200000001f},
{275.63038999999947f, 1906.8544000000009f},
{183.03502999999947f, 1790.0585000000008f},
{139.82933999999946f, 1584.5315000000007f},
{126.14346999999947f, 1438.2368000000008f},
{42.76099699999946f, 1285.7744000000007f},
{26.327842999999465f, 1234.0965000000008f},
{66.13150599999946f, 1292.2563000000007f},
{195.56353999999945f, 1563.1277000000007f},
{270.19608999999946f, 1785.8811000000007f},
{350.39938999999947f, 1839.1132000000007f},
{498.2929299999995f, 1843.8705000000007f},
{697.6644399999996f, 1869.0380000000007f},
{870.3013599999996f, 1977.2328000000007f},
{995.9601099999995f, 2133.549200000001f},
{1157.0598999999995f, 2260.404400000001f},
{1377.3040999999996f, 2308.570600000001f},
{1638.3271999999997f, 2273.947100000001f},
{1892.0618999999997f, 2260.656200000001f},
{1948.1276999999998f, 2276.8106000000007f},
{1973.5568999999998f, 2442.0905000000007f},
{1991.6138999999998f, 2520.9184000000005f},
{1991.6139f, 2520.9184f}
};

float[][] hole_in_arm_8 = new float[][]{
{3057.8763f, 2021.1657f},
{3058.0921f, 2062.6998f},
{3062.4386999999997f, 2073.3774f},
{3084.2472f, 2096.2399f},
{3105.6011999999996f, 2110.7893f},
{3126.6453999999994f, 2113.9676999999997f},
{3158.9240999999993f, 2110.9772f},
{3165.0870999999993f, 2109.6274f},
{3182.5798999999993f, 2100.6823999999997f},
{3188.4675999999995f, 2089.0897999999997f},
{3187.2460999999994f, 2083.1368999999995f},
{3168.1805999999992f, 2041.1261999999995f},
{3150.763899999999f, 2022.8579999999995f},
{3127.914399999999f, 2004.6671999999994f},
{3109.252499999999f, 2000.8657999999994f},
{3100.552599999999f, 2002.6722999999993f},
{3072.143699999999f, 2009.1635999999992f},
{3057.876299999999f, 2021.165699999999f},
{3057.8763f, 2021.1657f}
};

//////////////////////////////////////////////////////////////////
float[][] branch5_highlight = new float[][]{

};


float[][] branch4_highlight = new float[][]{

};

float[][] branch3_highlight = new float[][]{

};

float [][] branch2_highlight = new float[][]{

};

float [][] branch1_highlight = new float[][]{

};

float[][] arm1_1_highlight = new float [][]{
{1747.6716f, 2243.1571f},
{1490.7521f, 2099.9373f},
{1297.4938f, 2028.4903f},
{1258.8139999999999f, 2023.3869f},
{1076.6866999999997f, 2028.0908f},
{1027.0048999999997f, 2045.8452f}
};

float[][] arm1_2_highlight = new float [][]{
// {667.00487, 2145.8457},
// {590.63778, 2165.9093},
// {494.80441, 2131.2084999999997},
// {409.1859, 2037.2648999999997},
// {385.81837, 1947.4209999999996}
};

float[][] arm1_3_highlight = new float [][]{

};

float[][] arm8_highlight = new float[][]{
{2199.5628f, 2572.0214f},
{2291.9983f, 2470.0771f},
{2438.1823000000004f, 2449.5326999999997f},
{2634.8963000000003f, 2487.8322f},
{2930.0309f, 2544.384f},
{3125.9327000000003f, 2503.4819f},
{3181.4029000000005f, 2463.8014000000003f},
{3251.3752000000004f, 2323.3581000000004f},
{3269.9459000000006f, 2173.5121000000004f}
};

float[][] arm7_highlight = new float[][]{
{2379.4887f, 2144.3392f},
{2351.6340999999998f, 2181.6618f},
{2332.3828f, 2173.5782f},
{2338.7954f, 2115.8029f},
{2387.7971f, 2086.4909000000002f},
{2416.9233999999997f, 2081.2422f},
{2746.3702f, 1949.5693f},
{2903.8315f, 1807.1542000000002f},
{2927.8599999999997f, 1768.8340000000003f},
{3090.2891999999997f, 1581.3399000000004f},
{3237.5986999999996f, 1515.0200000000004f},
{3271.4947999999995f, 1515.0420000000004f}
};

float[][] arm6_highlight = new float[][]{
{2146.5667f, 1794.3055f},
{2264.1375f, 1635.915f},
{2509.5629999999996f, 1495.5672f},
{2621.5206f, 1433.9978999999998f},
{2916.5333f, 1051.3959f},
{2946.0258f, 955.36616f},
{2957.2205f, 945.0778700000001f}
};

float[][] arm5_highlight = new float[][]{
{2075.2874f, 755.83759f},
{1908.1428f, 1084.3008f},
{1899.7756000000002f, 1370.7526f},
{1932.0897000000002f, 1456.8268f},
{2150.2228000000005f, 1868.6601f},
{2178.5908000000004f, 1960.1401f},
{2181.6645000000003f, 1968.9103f}
};

float[][] arm4_1_highlight = new float[][]{
{1406.0861f, 1598.7549f},
{1370.0986f, 1659.5378999999998f},
{1352.0493000000001f, 1750.5671999999997f},
{1346.7019f, 1798.0218999999997f}
};

float[][] arm4_2_highlight = new float[][]{
{1587.5197f, 1405.7481f},
{1622.3404f, 1238.605f},
{1577.9752f, 1028.164f},
{1511.6716000000001f, 880.14013f}
};

float[][] arm4_3_highlight = new float[][]{
// {1409.0306, 761.43638},
// {1346.1142, 665.2337},
// {1327.6288, 599.29188},
// {1340.2323999999999, 550.04229},
// {1364.9844999999998, 523.51226}
};

float[][] arm4_4_highlight = new float[][]{

};

float[][] arm3_highlight = new float[][]{
// {1715.6716, 1599.8924},
// {1511.5555, 1449.2883},
// {1395.1136999999999, 1186.2008},
// {1412.7830999999999, 960.5289200000001},
// {1495.5701, 661.53683},
// {1494.9434999999999, 620.37293},
// {1359.272, 235.68659000000002},
// {1356.316, 91.67784300000002},
// {1356.5993, 89.73613800000003}
};

float[][] arm0_highlight = new float[][]{
// {225.59297, 1756.7551},
// {316.84755, 1858.7605},
// {518.86972, 1908.6023},
// {649.47843, 1957.8978},
// {797.5637, 2156.6026},
// {832.14535, 2181.6675},
// {1018.1344, 2429.6149},
// {1236.3722, 2503.35},
// {1281.1311, 2503.8901},
// {1521.6865, 2462.8140000000003},
// {1579.6331, 2457.2149000000004},
// {1892.8245, 2644.2736000000004},
// {1936.3987, 2773.2340000000004},
// {1958.6105, 2842.1659000000004},
// {1943.6512, 2853.5120000000006}
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
