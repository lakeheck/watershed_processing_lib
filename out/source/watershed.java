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

    float[][][] import_paths = new float[][][]{
        branch5, branch4, branch3, branch2, branch1, arm0, arm1, arm3, arm4, arm5, arm6, arm7, arm8
    };

    
    ArrayList<ArrayList<PVector>> paths = new ArrayList<ArrayList<PVector>>(); 
    for(float[][] p:import_paths){
        paths.add(inkscapePathImport(p, 3564.00000f, 5014.66650f));
    }
    render.beginDraw();
    for(ArrayList<PVector> path:paths){
        Ribbon r = new Ribbon(path, renderHighRes ? printDpi/previewDpi * 50 : 50, true);
        r.vadenWeb(200, 10, new Gradient(line_palette));
    }
    // printArray(path);
    line = new ArrayList(); //generate a random line 
    int n = 50;
    for(int i=0; i<50; i++){
        line.add(new PVector(map(i, 0, n-1, 0, renderWidth), renderHeight/2 + map(compoundTrigFunction(map(i, 0, n-1, 0, 2*TWO_PI), 0), -3, 4, -50, 50)));
    }

    // Ribbon r = new Ribbon(path, renderHighRes ? printDpi/previewDpi * 50 : 50, true);
    // render.beginDraw();
    // r.vadenWeb(200, 10, new Gradient(line_palette));


    float[] t = new float[400];
    float scale = 2;
    for(int i=0; i<t.length; i++){
        t[i] = map(i, 0, t.length, -scale*TWO_PI, scale*TWO_PI);
    }
    int sample_size = PApplet.parseInt(t.length*0.5f);
    float rectSize = renderWidth/sample_size; 
    int numRows=20;

    // for(int j=0; j < numRows; j++){
    //     int start = int(map(j, 0 , numRows, 0, renderWidth*2/3)); //choose starting point in t for this row
    //     ArrayList<PVector> tempLine = new ArrayList();

    //     for(int i=0; i<sample_size*2; i++){
    //         tempLine.add(new PVector(
    //             i*rectSize, 
    //             map(j, 0, numRows, renderHeight*0, renderHeight + randomGaussian()*(renderHighRes ? 10*printDpi/previewDpi : 10)) + (renderHighRes ? 10*printDpi/previewDpi : 10)*compoundTrigFunction(t[int((start+i)%t.length)], 0)
    //         ));
    //     }

    //     Ribbon tempRibbon = new Ribbon(tempLine, 20, false);
    //     tempRibbon.vadenWeb(500, 20, new Gradient(line_palette));
        
    // }
    render.stroke(0,0,100, 5);
    canvas_overlay_example1();
    // render.fill(0);

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
float[][] branch5 = new float[][]{
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


float[][] branch4 = new float[][]{
{2062.2236f, 3557.1023f},
{2086.4849999999997f, 3598.1602f},
{2132.4361999999996f, 3672.8669999999997f},
{2153.4025999999994f, 3716.7472999999995f},
{2151.4195999999993f, 3712.6903999999995f},
{2172.333199999999f, 3769.8565999999996f},
{2189.419199999999f, 3791.0795f},
{2216.1881999999987f, 3830.5937999999996f},
{2245.9414999999985f, 3872.6539999999995f},
{2281.2048999999984f, 3934.9322999999995f},
{2283.8235999999984f, 4028.8912999999993f},
{2224.872099999998f, 4097.011799999999f},
{2189.471599999998f, 4138.707299999999f},
{2142.193399999998f, 4213.1188999999995f},
{2138.3425999999977f, 4285.6801f},
{2151.4735999999975f, 4331.9459f},
{2168.6189999999974f, 4389.4145f},
{2171.9990999999973f, 4435.3849f},
{2190.995899999997f, 4471.5725f},
{2198.655899999997f, 4429.6213f},
{2185.964099999997f, 4347.0502f},
{2175.143099999997f, 4299.3243999999995f},
{2169.3722999999973f, 4250.3793f},
{2188.633499999997f, 4214.572499999999f},
{2227.1793999999973f, 4168.9884999999995f},
{2267.8216999999972f, 4126.4347f},
{2301.105599999997f, 4100.1856f},
{2299.209899999997f, 4145.8988f},
{2293.016599999997f, 4171.0449f},
{2294.1368999999972f, 4193.0049f},
{2314.169399999997f, 4238.2892999999995f},
{2331.951399999997f, 4275.583f},
{2357.093499999997f, 4331.602599999999f},
{2360.826099999997f, 4367.061599999999f},
{2355.138399999997f, 4425.716799999999f},
{2347.101499999997f, 4479.904999999999f},
{2329.826299999997f, 4529.932299999999f},
{2306.286599999997f, 4572.657799999998f},
{2287.128699999997f, 4635.157499999998f},
{2284.549099999997f, 4730.382399999999f},
{2302.607399999997f, 4771.457799999998f},
{2319.348799999997f, 4817.901099999998f},
{2327.604799999997f, 4856.839999999998f},
{2315.5212999999967f, 4927.834899999998f},
{2321.6344999999965f, 4902.679899999998f},
{2339.7636999999963f, 4864.925699999998f},
{2345.0610999999963f, 4826.341699999998f},
{2330.9657999999963f, 4788.837699999998f},
{2312.780499999996f, 4759.884099999998f},
{2305.520499999996f, 4700.562899999998f},
{2314.522799999996f, 4665.105299999998f},
{2325.366999999996f, 4620.8737999999985f},
{2351.494699999996f, 4563.0194999999985f},
{2378.004999999996f, 4504.087099999999f},
{2395.407099999996f, 4462.410899999999f},
{2417.9121999999957f, 4386.696299999999f},
{2431.7119999999954f, 4351.749799999999f},
{2473.1117999999956f, 4385.8103999999985f},
{2506.3947999999955f, 4471.916999999999f},
{2504.5284999999953f, 4512.630399999998f},
{2521.3248999999955f, 4460.374899999998f},
{2522.2549999999956f, 4412.170499999998f},
{2496.3940999999954f, 4379.441099999998f},
{2436.9300999999955f, 4320.167599999998f},
{2411.2714999999953f, 4283.052199999998f},
{2387.768299999995f, 4232.476799999999f},
{2381.354999999995f, 4180.435099999999f},
{2374.629399999995f, 4142.683299999999f},
{2367.520599999995f, 4100.536499999998f},
{2373.8898999999947f, 4049.7964999999986f},
{2394.192599999995f, 4022.2891999999983f},
{2409.386699999995f, 4066.443099999998f},
{2447.654499999995f, 4109.793299999998f},
{2490.269499999995f, 4132.938199999999f},
{2545.5863999999947f, 4159.906099999998f},
{2588.2275999999947f, 4173.8710999999985f},
{2628.9797999999946f, 4187.177899999999f},
{2726.614099999995f, 4219.626599999999f},
{2768.6720999999948f, 4237.182899999999f},
{2805.049399999995f, 4250.239899999999f},
{2846.1117999999947f, 4268.786099999998f},
{2870.8743999999947f, 4297.932399999999f},
{2859.1189999999947f, 4359.596599999999f},
{2859.1189999999947f, 4396.921899999999f},
{2860.9178999999945f, 4439.603899999999f},
{2857.6028999999944f, 4482.2928999999995f},
{2855.2199999999943f, 4521.750099999999f},
{2860.6034999999943f, 4573.2717999999995f},
{2847.9213999999943f, 4656.332799999999f},
{2842.8377999999943f, 4704.674499999999f},
{2843.6410999999944f, 4780.966199999999f},
{2861.0020999999942f, 4816.912199999999f},
{2888.9791999999943f, 4856.023199999999f},
{2939.3683999999944f, 4891.482299999999f},
{2965.0966999999946f, 4932.671299999999f},
{2965.0336999999945f, 4991.8063999999995f},
{2974.4469999999947f, 5033.935299999999f},
{2978.6274999999946f, 5043.5987f},
{2984.1584999999945f, 4998.7926f},
{2979.7688999999946f, 4922.6255f},
{2946.9552999999946f, 4871.2315f},
{2910.183799999995f, 4827.5138f},
{2886.3846999999946f, 4775.6417f},
{2874.2631999999944f, 4700.4229000000005f},
{2881.5138999999945f, 4641.402700000001f},
{2905.3929999999946f, 4600.405500000001f},
{2927.9521999999947f, 4562.463100000001f},
{2975.044199999995f, 4538.332700000001f},
{3016.2221999999947f, 4556.590400000001f},
{3072.9995999999946f, 4579.750900000001f},
{3112.9307999999946f, 4600.344800000001f},
{3151.768099999995f, 4640.369300000001f},
{3166.2342999999946f, 4742.149400000001f},
{3165.7926999999945f, 4811.403800000001f},
{3184.6860999999944f, 4844.747600000001f},
{3202.7798999999945f, 4878.598900000001f},
{3212.8378999999945f, 4926.716400000001f},
{3239.6326999999947f, 4960.717800000001f},
{3279.2679999999946f, 4978.290500000001f},
{3349.9465999999948f, 4988.528100000001f},
{3383.6222999999945f, 5010.256100000001f},
{3421.8000999999945f, 5036.911900000001f},
{3446.8494999999944f, 5056.683900000001f},
{3425.479599999994f, 5024.011000000001f},
{3397.5061999999944f, 4996.632100000002f},
{3348.0999999999945f, 4976.349100000001f},
{3302.6036999999947f, 4960.319700000001f},
{3269.857599999995f, 4945.696500000001f},
{3215.723399999995f, 4899.827200000001f},
{3181.8206999999948f, 4804.6116f},
{3176.1157999999946f, 4762.6622f},
{3179.4015999999947f, 4671.839f},
{3172.6512999999945f, 4587.281f},
{3142.7599999999948f, 4554.9059f},
{3108.2291999999948f, 4531.342799999999f},
{3071.1153999999947f, 4511.6541f},
{2981.8106999999945f, 4485.0747f},
{2943.2577999999944f, 4455.6926f},
{2931.8420999999944f, 4412.5585f},
{2939.3565999999946f, 4365.9493f},
{2962.6552999999944f, 4314.3777f},
{2972.8772999999946f, 4270.7084f},
{2965.233299999995f, 4227.2222f},
{2920.0890999999947f, 4194.5631f},
{2882.0379999999946f, 4178.064600000001f},
{2834.5727999999945f, 4167.907000000001f},
{2731.0871999999945f, 4131.315900000001f},
{2681.6874999999945f, 4114.313200000001f},
{2634.6071999999945f, 4097.856200000001f},
{2588.3910999999944f, 4067.023300000001f},
{2622.8726999999944f, 4057.0674000000013f},
{2664.8985999999945f, 4051.325000000001f},
{2762.0728999999947f, 4049.7965000000013f},
{2818.4172999999946f, 4064.262700000001f},
{2817.3206999999948f, 4062.1522000000014f},
{2860.9849999999947f, 4077.7905000000014f},
{2931.9076999999947f, 4088.8736000000013f},
{2978.3379999999947f, 4086.403200000001f},
{3068.1571999999946f, 4040.4053000000013f},
{3104.067799999995f, 4010.1099000000013f},
{3154.873999999995f, 3986.5902000000015f},
{3166.8920999999946f, 3984.7090000000017f},
{3169.5102999999945f, 3992.9837000000016f},
{3216.6535999999946f, 4001.9994000000015f},
{3258.4994999999944f, 4031.1338000000014f},
{3303.8796999999945f, 4080.0101000000013f},
{3369.2285999999945f, 4159.512400000001f},
{3396.4999999999945f, 4188.950600000001f},
{3432.5042999999946f, 4263.065300000001f},
{3461.9050999999945f, 4318.409800000001f},
{3454.8571999999945f, 4301.057000000002f},
{3429.3765999999946f, 4247.687700000001f},
{3404.0492999999947f, 4181.376100000001f},
{3382.0658999999946f, 4146.580700000001f},
{3331.7078999999944f, 4097.081900000001f},
{3270.4733999999944f, 4009.6607000000013f},
{3241.7523999999944f, 3980.535900000001f},
{3221.7714999999944f, 3949.728300000001f},
{3262.7073999999943f, 3929.374200000001f},
{3307.0994999999944f, 3903.851300000001f},
{3363.010399999994f, 3874.367500000001f},
{3415.5616999999943f, 3856.9388000000013f},
{3474.9863999999943f, 3840.7747000000013f},
{3503.554199999994f, 3791.8753000000015f},
{3528.234099999994f, 3752.6426000000015f},
{3546.0247999999942f, 3651.7375000000015f},
{3555.235799999994f, 3594.4277000000016f},
{3538.859899999994f, 3659.4680000000017f},
{3517.487699999994f, 3718.9089000000017f},
{3460.617099999994f, 3813.4278000000018f},
{3404.585899999994f, 3834.1362000000017f},
{3330.641799999994f, 3850.474600000002f},
{3225.940599999994f, 3890.0656000000017f},
{3177.100099999994f, 3927.0157000000017f},
{3132.065899999994f, 3944.8200000000015f},
{3074.530299999994f, 3973.6979000000015f},
{3006.752899999994f, 4013.7437000000014f},
{2964.587899999994f, 4021.4973000000014f},
{2908.973399999994f, 4025.4578000000015f},
{2827.0040999999937f, 4020.5072000000014f},
{2781.2566999999935f, 4010.1146000000012f},
{2741.3809999999935f, 3988.5375000000013f},
{2682.7013999999936f, 3985.744700000001f},
{2564.1724999999938f, 4010.7342000000012f},
{2503.1674999999937f, 3978.9258000000013f},
{2417.3554999999938f, 3856.7403000000013f},
{2400.0172999999936f, 3818.3795000000014f},
{2358.4370999999937f, 3776.433600000001f},
{2328.3279999999936f, 3753.146700000001f},
{2285.1426999999935f, 3718.353100000001f},
{2193.0397999999936f, 3614.524400000001f},
{2160.4837999999936f, 3554.146500000001f},
{2112.9869999999937f, 3556.828300000001f},
{2079.0195999999937f, 3564.567400000001f},
{2062.2232999999937f, 3557.102300000001f},
{2062.2236f, 3557.1023f}
};

float[][] branch3 = new float[][]{
{1983.8404f, 3521.6433f},
{2011.4925f, 3579.8854f},
{2025.8249f, 3645.242f},
{2034.6752000000001f, 3685.4999000000003f},
{2051.5463f, 3722.8613000000005f},
{2062.1772f, 3828.6519000000003f},
{2038.4787000000001f, 3900.8179000000005f},
{2024.3529f, 3941.7770000000005f},
{2002.0516f, 3980.3911000000003f},
{1982.2072f, 4001.7771000000002f},
{1927.6826f, 4039.1202000000003f},
{1909.1898f, 4053.5290000000005f},
{1868.0867f, 4094.0512000000003f},
{1812.8929f, 4138.7924f},
{1781.5795f, 4205.0187000000005f},
{1749.7481f, 4267.077200000001f},
{1732.2205000000001f, 4307.479200000001f},
{1715.0981000000002f, 4340.933900000001f},
{1700.4609000000003f, 4387.962200000001f},
{1660.3469000000002f, 4469.987400000001f},
{1620.8523000000002f, 4491.648700000001f},
{1588.6158000000003f, 4499.393300000001f},
{1574.6874000000003f, 4501.854300000001f},
{1533.2164000000002f, 4522.860300000001f},
{1478.3308000000002f, 4551.498400000001f},
{1443.2275000000002f, 4574.8870000000015f},
{1396.7649000000001f, 4629.210300000002f},
{1379.1703000000002f, 4658.199100000002f},
{1361.4869000000003f, 4688.457700000002f},
{1339.1028000000003f, 4718.476600000002f},
{1318.8250000000003f, 4778.906700000001f},
{1291.6480000000004f, 4813.173900000002f},
{1259.7293000000004f, 4844.825600000002f},
{1250.0872000000004f, 4877.209800000002f},
{1280.5087000000003f, 4842.536500000002f},
{1308.5934000000004f, 4816.693800000002f},
{1332.5137000000004f, 4792.570200000002f},
{1360.3250000000005f, 4768.567600000002f},
{1378.6864000000005f, 4726.489600000002f},
{1405.3793000000005f, 4691.419000000002f},
{1418.1959000000006f, 4675.552400000001f},
{1425.8270000000007f, 4665.664100000002f},
{1426.9694000000006f, 4665.464900000002f},
{1413.0872000000006f, 4721.694300000002f},
{1381.0366000000006f, 4788.837700000002f},
{1356.7929000000006f, 4855.672800000002f},
{1337.1246000000006f, 4892.792700000002f},
{1291.1811000000005f, 4933.6181000000015f},
{1256.3454000000004f, 4955.709000000002f},
{1223.2378000000003f, 4964.440200000002f},
{1183.7136000000003f, 4990.4486000000015f},
{1164.7144000000003f, 5031.668500000002f},
{1140.6403000000003f, 5078.457900000001f},
{1119.2652000000003f, 5145.202800000001f},
{1136.5558000000003f, 5128.4980000000005f},
{1161.1834000000003f, 5091.9588f},
{1189.3681000000004f, 5059.345700000001f},
{1214.9389000000003f, 5027.719700000001f},
{1245.4316000000003f, 4989.9457f},
{1286.5167000000004f, 4967.206f},
{1324.3481000000004f, 4948.6033f},
{1373.4581000000003f, 4910.7024f},
{1389.9076000000002f, 4882.9128f},
{1399.6992000000002f, 4844.8256f},
{1403.0461000000003f, 4779.9176f},
{1427.8146000000002f, 4702.277999999999f},
{1443.1965000000002f, 4665.6927f},
{1482.5522000000003f, 4613.6554f},
{1510.6824000000004f, 4589.409299999999f},
{1556.1944000000003f, 4581.7330999999995f},
{1606.5600000000004f, 4570.157399999999f},
{1669.4939000000004f, 4537.989f},
{1695.2361000000003f, 4512.3126999999995f},
{1718.2367000000004f, 4468.8762f},
{1732.8211000000003f, 4418.179f},
{1774.4386000000004f, 4335.6686f},
{1794.3389000000004f, 4301.9352f},
{1818.2796000000005f, 4251.0401f},
{1825.2078000000006f, 4214.0279f},
{1866.0828000000006f, 4176.3939f},
{1918.8989000000006f, 4158.7979000000005f},
{1951.3552000000007f, 4142.273f},
{1980.5413000000005f, 4127.7095f},
{2014.3315000000005f, 4116.1982f},
{2023.0271000000005f, 4124.655299999999f},
{2018.8666000000005f, 4175.6927f},
{2015.4875000000004f, 4214.0192f},
{2014.1845000000003f, 4276.3933f},
{2017.3392000000003f, 4337.2832f},
{2010.1943000000003f, 4371.8029f},
{1982.8041000000003f, 4396.1619f},
{1931.2481000000002f, 4440.6254f},
{1889.9590000000003f, 4475.5768f},
{1813.5389000000002f, 4549.0329f},
{1757.4939000000002f, 4602.4137f},
{1745.0226000000002f, 4682.1246f},
{1728.8350000000003f, 4712.2783f},
{1708.4998000000003f, 4753.629f},
{1681.2337000000002f, 4785.0775f},
{1655.7396000000003f, 4819.5518f},
{1640.1226000000004f, 4865.1209f},
{1626.7533000000003f, 4930.5715f},
{1615.6141000000002f, 4995.0469f},
{1610.4528000000003f, 5019.8928000000005f},
{1641.0179000000003f, 4973.238f},
{1671.7722000000003f, 4917.5315f},
{1690.5795000000003f, 4828.1894f},
{1730.4187000000002f, 4782.5144f},
{1786.0955000000001f, 4717.2051f},
{1798.5558f, 4662.1084f},
{1833.1308000000001f, 4593.0966f},
{1895.4350000000002f, 4543.5389f},
{1927.6356f, 4520.257f},
{1978.2586000000001f, 4509.1061f},
{2023.0577f, 4495.8922f},
{2048.5619f, 4475.581700000001f},
{2079.3740000000003f, 4434.304300000001f},
{2084.7967000000003f, 4429.3932f},
{2097.5098000000003f, 4464.044900000001f},
{2088.3739f, 4505.3955000000005f},
{2031.1887000000002f, 4558.538100000001f},
{1990.3666f, 4598.008400000001f},
{1971.8590000000002f, 4632.013800000001f},
{1953.3532000000002f, 4665.6611f},
{1924.1199000000001f, 4704.8557f},
{1903.6956000000002f, 4738.1021f},
{1884.0694000000003f, 4793.5952f},
{1876.9069000000004f, 4803.271299999999f},
{1851.9023000000004f, 4826.774399999999f},
{1828.9403000000004f, 4869.0871f},
{1794.3037000000004f, 4966.8829f},
{1722.5632000000005f, 5009.057f},
{1671.6380000000006f, 5031.604899999999f},
{1621.9632000000006f, 5058.083299999999f},
{1591.1896000000006f, 5085.674599999999f},
{1641.1797000000006f, 5064.554999999999f},
{1681.4995000000006f, 5043.457399999999f},
{1725.7535000000005f, 5022.124599999999f},
{1798.1214000000004f, 4997.277399999999f},
{1833.7415000000005f, 4989.712199999999f},
{1871.8905000000004f, 4977.332499999999f},
{1904.6005000000005f, 4924.050099999999f},
{1931.5850000000005f, 4885.883499999999f},
{1936.8347000000006f, 4844.408399999999f},
{1942.5441000000005f, 4800.314399999999f},
{1944.5945000000006f, 4760.373799999999f},
{1955.8464000000006f, 4717.919599999999f},
{2000.4201000000005f, 4652.839099999999f},
{2029.7144000000005f, 4615.636099999999f},
{2066.0103000000004f, 4568.983299999999f},
{2103.5879000000004f, 4524.690899999999f},
{2135.0454000000004f, 4470.557399999999f},
{2145.2498000000005f, 4423.388799999999f},
{2138.1539000000007f, 4385.686999999999f},
{2122.1736000000005f, 4332.918599999999f},
{2106.6114000000007f, 4290.815799999999f},
{2095.1665000000007f, 4220.735f},
{2093.9691000000007f, 4125.573899999999f},
{2103.662500000001f, 4076.9091999999996f},
{2118.0356000000006f, 4044.1031f},
{2132.5715000000005f, 4010.3365f},
{2148.0718000000006f, 3971.4132999999997f},
{2166.7344000000007f, 3924.7567f},
{2184.592600000001f, 3884.7246999999998f},
{2203.047600000001f, 3845.6294999999996f},
{2181.9385000000007f, 3730.1178999999997f},
{2163.2840000000006f, 3669.3298999999997f},
{2148.0718000000006f, 3603.7589999999996f},
{2134.2083000000007f, 3554.0279999999993f},
{2119.1447000000007f, 3521.6432999999993f},
{2077.617300000001f, 3513.8195999999994f},
{2015.5669000000007f, 3521.6432999999993f},
{1983.8404000000007f, 3521.6432999999993f},
{1983.8404f, 3521.6433f}
};

float [][] branch2 = new float[][]{
{1994.2528f, 3583.8455f},
{2002.2511f, 3682.2461f},
{1989.3866f, 3711.3772f},
{1982.7363f, 3729.3449f},
{1957.9147f, 3783.4828f},
{1913.1275f, 3854.4865000000004f},
{1879.9450000000002f, 3889.5575000000003f},
{1836.4054f, 3934.088f},
{1747.2903000000001f, 4008.5657f},
{1719.9474f, 4042.4419000000003f},
{1690.0142f, 4071.3785000000003f},
{1652.135f, 4073.5041f},
{1612.9446f, 4047.9558f},
{1565.5198f, 4036.1423f},
{1496.7451f, 4036.7326f},
{1442.6234000000002f, 4062.8603f},
{1366.8088000000002f, 4106.7474999999995f},
{1325.0486000000003f, 4141.2435f},
{1288.2578000000003f, 4172.1485999999995f},
{1274.4391000000003f, 4200.4511999999995f},
{1268.9733000000003f, 4245.104899999999f},
{1282.4813000000004f, 4229.441999999999f},
{1347.4241000000004f, 4168.367399999999f},
{1375.2321000000004f, 4139.907999999999f},
{1417.6800000000003f, 4117.4412999999995f},
{1467.8950000000002f, 4090.3180999999995f},
{1562.8979000000002f, 4069.0411999999997f},
{1610.5873000000001f, 4088.9880999999996f},
{1637.9947000000002f, 4119.7146999999995f},
{1580.7396f, 4120.018999999999f},
{1509.6781f, 4127.676799999999f},
{1463.3838f, 4152.209899999999f},
{1438.7644f, 4186.097699999999f},
{1403.7753f, 4236.588599999999f},
{1382.7172f, 4279.420399999999f},
{1347.4438f, 4288.678499999999f},
{1305.3919f, 4320.384499999999f},
{1255.9968000000001f, 4346.532699999999f},
{1219.0202000000002f, 4342.320499999999f},
{1178.7257000000002f, 4320.921299999999f},
{1115.6041000000002f, 4307.657399999999f},
{1074.4805000000003f, 4291.9944f},
{1021.9179000000004f, 4261.4854f},
{975.0018800000004f, 4256.1491f},
{939.7894000000003f, 4273.054999999999f},
{893.3926600000003f, 4272.6633999999995f},
{847.8579800000003f, 4253.3753f},
{782.3807900000004f, 4218.7955999999995f},
{748.3725300000004f, 4198.164599999999f},
{708.8484400000004f, 4179.951999999999f},
{656.5282400000004f, 4208.524199999999f},
{640.6235800000004f, 4219.883499999999f},
{716.6460100000004f, 4214.027899999999f},
{747.8227300000004f, 4243.109299999999f},
{795.3032100000005f, 4278.933899999999f},
{887.5875600000005f, 4314.387099999999f},
{938.2849700000005f, 4326.076399999999f},
{987.6663400000004f, 4308.282399999999f},
{1027.2764000000004f, 4298.737999999999f},
{1038.6406000000004f, 4307.641799999999f},
{1076.7839000000004f, 4332.592f},
{1080.0658000000003f, 4346.5063f},
{1056.4349000000002f, 4376.3264f},
{1001.9624000000002f, 4411.5531f},
{957.4106500000003f, 4426.883f},
{915.4550000000003f, 4425.9881f},
{866.6801400000003f, 4389.429499999999f},
{820.9361200000003f, 4377.8625999999995f},
{782.2724100000003f, 4390.332399999999f},
{727.3665100000003f, 4393.097099999999f},
{632.1581600000003f, 4383.670499999999f},
{584.4248300000003f, 4371.777899999999f},
{523.1708500000003f, 4353.457699999999f},
{450.3703100000003f, 4361.670799999999f},
{396.8677900000003f, 4341.213699999999f},
{340.3154600000003f, 4319.325599999999f},
{252.2927900000003f, 4331.096f},
{217.2794900000003f, 4349.0921f},
{181.0398800000003f, 4387.4232999999995f},
{216.2287400000003f, 4380.8412f},
{269.0326000000003f, 4362.0788f},
{333.59743000000026f, 4344.1562f},
{378.46992000000023f, 4371.626700000001f},
{450.3224300000002f, 4393.5132f},
{476.6468400000002f, 4393.655900000001f},
{460.96762000000024f, 4432.380900000001f},
{449.7363600000002f, 4492.293300000001f},
{426.4407500000002f, 4545.849500000001f},
{428.6311300000002f, 4563.253300000001f},
{458.6366000000002f, 4528.165600000001f},
{479.44925000000023f, 4476.747300000001f},
{496.42666000000025f, 4426.782100000001f},
{530.0656700000003f, 4396.110300000001f},
{577.4833600000003f, 4396.030300000001f},
{622.2706800000003f, 4425.494100000002f},
{664.2630100000003f, 4451.0556000000015f},
{643.4844800000003f, 4500.973600000001f},
{622.0175400000003f, 4544.247700000001f},
{594.1657300000003f, 4608.942100000001f},
{578.4977600000003f, 4642.494400000001f},
{560.1842800000003f, 4684.700700000001f},
{510.06438000000026f, 4713.644700000002f},
{454.69505000000026f, 4743.924800000002f},
{415.3214200000003f, 4823.657300000002f},
{393.49643000000026f, 4903.624300000001f},
{426.16758000000027f, 4864.269000000001f},
{462.31681000000026f, 4809.577900000001f},
{509.10022000000026f, 4757.031000000001f},
{523.0754600000002f, 4737.855700000001f},
{562.8994800000003f, 4717.185200000001f},
{619.6002000000003f, 4647.001200000001f},
{628.9229700000003f, 4601.987000000001f},
{642.8084500000003f, 4576.184600000001f},
{656.2591000000003f, 4556.143500000001f},
{676.5691500000004f, 4523.508000000001f},
{724.2649900000004f, 4454.418100000001f},
{726.3264200000003f, 4427.884900000001f},
{764.5280700000003f, 4438.486900000001f},
{824.8038500000002f, 4436.283f},
{820.8112700000003f, 4462.1917f},
{798.7930900000002f, 4501.4497f},
{778.1702700000002f, 4533.3777f},
{754.4284200000002f, 4580.2886f},
{738.6634600000002f, 4631.9228f},
{740.7730600000002f, 4691.8936f},
{749.6279400000002f, 4690.239100000001f},
{749.0064700000001f, 4626.515700000001f},
{786.6828900000002f, 4557.1603000000005f},
{819.3058000000002f, 4517.193700000001f},
{849.1508700000002f, 4482.769800000001f},
{880.7094600000001f, 4458.3766000000005f},
{914.2261400000001f, 4480.6171f},
{974.2541400000001f, 4492.2052f},
{993.3839000000002f, 4494.245f},
{991.0229400000002f, 4535.1858f},
{1000.5580000000002f, 4552.0702f},
{1027.1505000000002f, 4508.7147f},
{1052.5373000000002f, 4456.1666000000005f},
{1093.6317000000001f, 4391.322800000001f},
{1096.5843000000002f, 4394.027500000001f},
{1147.9065000000003f, 4407.444000000001f},
{1192.5438000000004f, 4419.316800000001f},
{1243.9186000000004f, 4421.920300000002f},
{1296.8127000000004f, 4408.213600000002f},
{1341.8450000000005f, 4378.258900000002f},
{1346.2184000000004f, 4404.024300000002f},
{1342.4494000000004f, 4460.695500000002f},
{1275.9298000000003f, 4556.436100000002f},
{1244.3385000000003f, 4600.807700000002f},
{1226.3073000000004f, 4654.015200000003f},
{1176.2798000000005f, 4717.142900000003f},
{1162.6835000000005f, 4755.244600000003f},
{1162.4408000000005f, 4792.782100000003f},
{1171.8372000000006f, 4762.584800000003f},
{1207.9040000000007f, 4716.449700000003f},
{1244.2850000000008f, 4677.8749000000025f},
{1277.2336000000007f, 4624.873800000002f},
{1307.5291000000007f, 4584.658000000002f},
{1323.7960000000007f, 4547.484900000002f},
{1356.5627000000006f, 4499.553800000002f},
{1382.1661000000006f, 4460.443100000002f},
{1382.3980000000006f, 4399.764800000003f},
{1392.9530000000007f, 4351.141100000003f},
{1422.1384000000007f, 4299.548500000003f},
{1460.5838000000008f, 4266.665500000003f},
{1494.7011000000007f, 4212.425600000003f},
{1540.6408000000006f, 4190.219900000003f},
{1631.8308000000006f, 4192.809200000002f},
{1718.8307000000007f, 4184.167300000003f},
{1754.6327000000006f, 4159.036200000003f},
{1785.6045000000006f, 4127.888100000003f},
{1825.8707000000006f, 4093.2619000000027f},
{1841.8881000000006f, 4081.7230000000027f},
{1854.8196000000005f, 4066.2445000000025f},
{1893.9240000000004f, 4022.9380000000024f},
{1956.2476000000004f, 3980.6132000000025f},
{2005.3326000000004f, 3935.6562000000026f},
{2026.2678000000003f, 3916.3167000000026f},
{2058.4910000000004f, 3866.9021000000025f},
{2076.2241000000004f, 3822.1772000000024f},
{2087.3476000000005f, 3760.7381000000023f},
{2098.4149000000007f, 3709.8505000000023f},
{2114.4790000000007f, 3646.6828000000023f},
{2110.7465000000007f, 3579.4972000000025f},
{2089.4393000000005f, 3544.1945000000023f},
{2028.4759000000004f, 3560.194200000002f},
{1999.0859000000003f, 3585.2907000000023f},
{1997.4849000000002f, 3584.9216000000024f},
{1994.2528000000002f, 3583.8452000000025f},
{1994.2528f, 3583.8455f}
};

float [][] branch1 = new float[][]{
{1897.9922f, 3040.1468f},
{1897.5349999999999f, 3093.3339f},
{1899.8330999999998f, 3107.3353f},
{1904.2281999999998f, 3167.1629000000003f},
{1914.9341999999997f, 3277.5867000000003f},
{1938.3916999999997f, 3331.3987f},
{1946.7226999999996f, 3362.7791f},
{1953.1697999999997f, 3385.735f},
{1990.4482999999996f, 3480.6939f},
{1996.1086999999995f, 3513.5749f},
{2000.0111999999995f, 3561.3171f},
{1980.5148999999994f, 3606.5394f},
{1952.1138999999994f, 3618.6891f},
{1908.7181999999993f, 3638.5207f},
{1863.3115999999993f, 3664.7333f},
{1833.4197999999992f, 3687.5868f},
{1786.6405999999993f, 3715.5001f},
{1753.0812999999994f, 3722.3768f},
{1710.4117999999994f, 3729.9385f},
{1655.6504999999995f, 3728.8815000000004f},
{1594.0037999999995f, 3716.5519000000004f},
{1528.9085999999995f, 3696.0539000000003f},
{1468.5846999999994f, 3680.9438000000005f},
{1467.1273999999994f, 3682.2476000000006f},
{1403.5176999999994f, 3684.7467000000006f},
{1358.6619999999994f, 3689.5646000000006f},
{1321.8414999999993f, 3697.7901000000006f},
{1274.4785999999992f, 3725.4204000000004f},
{1252.8588999999993f, 3748.3598000000006f},
{1213.0726999999993f, 3773.589200000001f},
{1180.3767999999993f, 3786.399500000001f},
{1132.2396999999994f, 3809.723800000001f},
{1110.1733999999994f, 3822.325600000001f},
{1096.8231999999994f, 3834.432000000001f},
{1060.4977999999994f, 3858.630500000001f},
{1054.2667999999994f, 3860.200000000001f},
{1023.4665999999994f, 3878.1200000000013f},
{987.6143699999994f, 3895.6439000000014f},
{954.0269399999994f, 3918.6792000000014f},
{931.3498499999994f, 3936.6097000000013f},
{904.5533499999993f, 3956.374900000001f},
{860.1402299999993f, 3996.607600000001f},
{850.1472799999993f, 4043.3925000000013f},
{840.6220599999994f, 4085.851300000001f},
{830.5071599999993f, 4112.324200000001f},
{820.3773199999994f, 4150.3676000000005f},
{822.5908199999993f, 4157.419800000001f},
{838.8370099999993f, 4137.955400000001f},
{874.0236699999992f, 4057.341300000001f},
{901.2720999999992f, 4017.611900000001f},
{934.4251199999992f, 3997.318300000001f},
{967.2597999999992f, 3963.6467000000007f},
{1001.4209999999993f, 3932.3358000000007f},
{1028.3123999999993f, 3911.6928000000007f},
{1100.3712999999993f, 3893.6237000000006f},
{1156.1048999999994f, 3878.5633000000007f},
{1177.9035999999994f, 3855.950100000001f},
{1215.0605999999993f, 3827.588000000001f},
{1280.4333999999994f, 3793.091300000001f},
{1366.8092999999994f, 3755.663500000001f},
{1450.1614999999995f, 3753.890800000001f},
{1474.8482999999994f, 3774.360300000001f},
{1531.2327999999993f, 3811.003200000001f},
{1568.0392999999992f, 3817.290500000001f},
{1596.6338999999991f, 3826.266000000001f},
{1632.9823999999992f, 3835.175900000001f},
{1670.208299999999f, 3849.200000000001f},
{1703.900499999999f, 3851.9723000000013f},
{1748.832599999999f, 3845.965300000001f},
{1805.777499999999f, 3837.053800000001f},
{1802.863499999999f, 3834.178800000001f},
{1842.174399999999f, 3824.789100000001f},
{1877.543699999999f, 3805.011000000001f},
{1913.362999999999f, 3781.672400000001f},
{1972.642799999999f, 3747.461400000001f},
{2019.324699999999f, 3690.163100000001f},
{2073.004499999999f, 3634.481400000001f},
{2094.6113999999993f, 3577.915500000001f},
{2086.484999999999f, 3523.509600000001f},
{2069.0819999999994f, 3452.032400000001f},
{2028.4942999999994f, 3389.982800000001f},
{1979.9736999999993f, 3299.061400000001f},
{1971.9495999999992f, 3245.509500000001f},
{1973.9376999999993f, 3256.019300000001f},
{1959.5788999999993f, 3137.1926000000008f},
{1946.2001999999993f, 3080.6192000000005f},
{1910.4776999999992f, 3052.0859000000005f},
{1897.9921999999992f, 3040.1468000000004f},
{1897.9922f, 3040.1468f}
};

float[][] arm1 = new float [][]{
{375.1194f, 1310.1185f},
{350.85794f, 1343.7113f},
{330.82099999999997f, 1377.1331f},
{330.32901999999996f, 1416.4956f},
{353.37647f, 1460.8446999999999f},
{374.6021f, 1491.6475999999998f},
{398.73936000000003f, 1516.8753999999997f},
{430.95067000000006f, 1571.8234999999997f},
{427.88301000000007f, 1620.2502999999997f},
{409.4663300000001f, 1662.9919999999997f},
{392.9759700000001f, 1705.7423999999996f},
{352.9164600000001f, 1792.1244999999997f},
{341.3870100000001f, 1802.7073999999998f},
{336.0500300000001f, 1816.2471999999998f},
{313.0733500000001f, 1855.0607999999997f},
{309.8763500000001f, 1900.8019999999997f},
{308.76485000000014f, 1950.7950999999996f},
{317.30048000000016f, 1999.7080999999996f},
{335.2285200000002f, 2042.9545999999996f},
{360.14957000000015f, 2077.1790999999994f},
{386.17568000000017f, 2116.8431999999993f},
{412.89428000000015f, 2144.3550999999993f},
{456.75238000000013f, 2175.180799999999f},
{514.1273100000001f, 2198.0052999999994f},
{551.01057f, 2212.1961999999994f},
{604.9096400000001f, 2220.8802999999994f},
{674.6359300000001f, 2215.3441999999995f},
{740.9074700000001f, 2204.0596999999993f},
{782.5817800000001f, 2193.9830999999995f},
{869.0699500000001f, 2161.9743999999996f},
{932.8750600000001f, 2144.3282999999997f},
{1001.8627000000001f, 2103.1211999999996f},
{1046.4127f, 2100.3902999999996f},
{1122.6506000000002f, 2098.4547f},
{1173.0142f, 2088.2940999999996f},
{1254.1305f, 2088.3510999999994f},
{1306.5825f, 2094.4618999999993f},
{1383.2992f, 2103.873899999999f},
{1429.1993f, 2125.9593999999993f},
{1456.3545f, 2154.1637999999994f},
{1491.3601999999998f, 2194.4866999999995f},
{1520.7730999999999f, 2226.4671999999996f},
{1550.1777f, 2254.8288999999995f},
{1570.5874f, 2300.2016999999996f},
{1599.2212f, 2350.9796999999994f},
{1643.5182f, 2381.9737999999993f},
{1682.2149f, 2393.9420999999993f},
{1729.8223f, 2402.8631999999993f},
{1774.8186f, 2416.8137999999994f},
{1835.2225f, 2447.8509999999997f},
{1963.3143f, 2462.6269999999995f},
{2090.6757f, 2475.1863999999996f},
{2133.9599999999996f, 2467.0224999999996f},
{2137.2568999999994f, 2400.2101999999995f},
{2155.5367999999994f, 2338.4306999999994f},
{2191.4045999999994f, 2316.4737999999993f},
{2179.2607999999996f, 2279.219899999999f},
{2124.0070999999994f, 2259.993299999999f},
{2039.7152999999994f, 2187.408499999999f},
{1976.3752999999995f, 2166.734199999999f},
{1959.5788999999995f, 2148.071599999999f},
{1905.4571999999996f, 2054.758299999999f},
{1889.6363999999996f, 2022.035499999999f},
{1851.2679999999996f, 1994.191999999999f},
{1813.0279999999996f, 1963.3910999999991f},
{1767.3535999999995f, 1939.049799999999f},
{1728.6664999999994f, 1921.245599999999f},
{1683.2400999999993f, 1899.871499999999f},
{1628.7222999999992f, 1880.1917999999991f},
{1584.5684999999992f, 1860.008999999999f},
{1540.921199999999f, 1847.0478999999991f},
{1481.8148999999992f, 1834.538899999999f},
{1437.368599999999f, 1819.627899999999f},
{1371.7052999999992f, 1814.009999999999f},
{1322.2085999999993f, 1820.9878999999992f},
{1253.5233999999994f, 1815.096899999999f},
{1209.3400999999994f, 1821.4750999999992f},
{1170.0792999999994f, 1838.1007999999993f},
{1130.9569999999994f, 1853.2015999999992f},
{1073.1427999999994f, 1876.644899999999f},
{1035.4608999999994f, 1891.2088999999992f},
{993.8132799999994f, 1907.645799999999f},
{921.2562299999994f, 1942.5142999999991f},
{874.2449499999993f, 1954.9560999999992f},
{825.6072099999993f, 1964.207399999999f},
{742.4124099999993f, 1976.4970999999991f},
{672.5281199999994f, 2000.894099999999f},
{577.3425099999994f, 2005.662399999999f},
{521.0990899999994f, 1981.086699999999f},
{486.90924999999936f, 1937.451599999999f},
{478.96254999999934f, 1885.334299999999f},
{485.67671999999936f, 1804.462199999999f},
{493.9954299999994f, 1769.132299999999f},
{519.1170699999993f, 1642.119999999999f},
{495.83981999999935f, 1586.702399999999f},
{478.80746999999934f, 1543.455799999999f},
{441.97206999999935f, 1510.1565999999991f},
{416.8849699999993f, 1477.555399999999f},
{402.1803599999993f, 1468.320999999999f},
{348.99167999999935f, 1395.966499999999f},
{356.3198599999993f, 1352.903099999999f},
{371.6297999999993f, 1320.5061999999991f},
{375.1193999999993f, 1310.1182999999992f},
{375.1194f, 1310.1185f}
};

float[][] arm8 = new float[][]{
{2153.6673f, 2781.8203f},
{2165.8671f, 2729.8174999999997f},
{2176.5865f, 2695.3262999999997f},
{2190.8509f, 2661.2230999999997f},
{2221.0715f, 2610.5627999999997f},
{2268.4768f, 2569.3567999999996f},
{2323.9208f, 2536.5387999999994f},
{2365.893f, 2512.9904999999994f},
{2459.8259f, 2515.2511999999992f},
{2504.3659f, 2527.0232999999994f},
{2568.0371999999998f, 2542.9637999999995f},
{2610.7295999999997f, 2549.9287999999997f},
{2662.4782999999998f, 2553.7430999999997f},
{2706.6701f, 2559.4493999999995f},
{2726.6699999999996f, 2563.9681999999993f},
{2745.8958f, 2568.8715999999995f},
{2786.8992999999996f, 2579.8418999999994f},
{2840.0774999999994f, 2590.0317999999993f},
{2885.9155999999994f, 2595.154499999999f},
{2988.2427999999995f, 2608.157199999999f},
{3025.4383999999995f, 2596.798399999999f},
{3071.5580999999997f, 2577.066999999999f},
{3146.4903999999997f, 2557.256299999999f},
{3180.3542999999995f, 2531.086999999999f},
{3203.6289999999995f, 2512.6270999999992f},
{3248.9760999999994f, 2475.661699999999f},
{3274.7053999999994f, 2447.5578999999993f},
{3301.9738999999995f, 2407.2507999999993f},
{3318.7997999999993f, 2379.2780999999995f},
{3328.4704999999994f, 2350.4435999999996f},
{3338.7121999999995f, 2280.3535999999995f},
{3339.9818999999993f, 2219.3197999999993f},
{3337.3679999999995f, 2199.6651999999995f},
{3333.7458999999994f, 2151.6316999999995f},
{3319.1417999999994f, 2116.0153999999993f},
{3309.5044999999996f, 2089.0686999999994f},
{3331.8678999999997f, 2041.2665999999992f},
{3344.4215f, 2015.9423999999992f},
{3361.921f, 1980.6652999999992f},
{3370.9417999999996f, 1941.146699999999f},
{3377.1106999999997f, 1895.031799999999f},
{3392.7162f, 1853.552899999999f},
{3437.5054f, 1814.813399999999f},
{3484.2526f, 1806.949499999999f},
{3524.7558f, 1791.918499999999f},
{3506.9269999999997f, 1782.4499999999991f},
{3458.5746999999997f, 1782.2396999999992f},
{3422.5903999999996f, 1802.1535999999992f},
{3364.5952999999995f, 1880.8827999999992f},
{3332.1175999999996f, 1935.611399999999f},
{3316.4820999999997f, 1972.7215999999992f},
{3277.8907999999997f, 2031.5802999999992f},
{3271.8904999999995f, 2024.327099999999f},
{3247.8247999999994f, 2004.527299999999f},
{3217.1349999999993f, 1973.742399999999f},
{3186.1545999999994f, 1944.892099999999f},
{3120.9700999999995f, 1918.769799999999f},
{3092.3507999999997f, 1912.1486999999988f},
{3049.1976999999997f, 1909.221699999999f},
{3004.3738999999996f, 1924.592699999999f},
{2969.6041999999998f, 1952.115299999999f},
{2932.2603f, 2008.505899999999f},
{2932.5525f, 2058.6072999999988f},
{2945.4705999999996f, 2106.811299999999f},
{2970.9601999999995f, 2139.079999999999f},
{2996.8834999999995f, 2165.874899999999f},
{3032.5535999999993f, 2190.617499999999f},
{3079.3998999999994f, 2216.616399999999f},
{3144.1964999999996f, 2217.873499999999f},
{3188.2721999999994f, 2202.4943999999987f},
{3217.5212999999994f, 2174.6251999999986f},
{3206.4551999999994f, 2211.9415999999987f},
{3146.6785999999993f, 2303.725399999999f},
{3105.185299999999f, 2319.422899999999f},
{3037.311399999999f, 2332.112799999999f},
{2959.194299999999f, 2338.518299999999f},
{2906.674499999999f, 2338.344999999999f},
{2859.722499999999f, 2332.129499999999f},
{2773.076399999999f, 2327.936199999999f},
{2720.785299999999f, 2313.605399999999f},
{2652.559799999999f, 2292.718699999999f},
{2590.801699999999f, 2287.938799999999f},
{2539.334799999999f, 2278.394099999999f},
{2480.480499999999f, 2261.293599999999f},
{2317.4395999999992f, 2252.031599999999f},
{2240.1997999999994f, 2238.213099999999f},
{2130.0451999999996f, 2277.921199999999f},
{2049.5745999999995f, 2284.732999999999f},
{1970.5321999999994f, 2334.4831999999988f},
{1900.3953999999994f, 2430.922099999999f},
{1867.2225999999994f, 2558.232399999999f},
{1861.3795999999993f, 2646.046299999999f},
{1877.6501999999994f, 2743.522599999999f},
{1920.0113999999994f, 2789.526899999999f},
{1961.3874999999994f, 2803.5764999999988f},
{2056.6136999999994f, 2816.545999999999f},
{2099.2230999999992f, 2820.816199999999f},
{2143.156899999999f, 2812.710499999999f},
{2153.667299999999f, 2781.820299999999f},
{2153.6673f, 2781.8203f}
};

float[][] arm7 = new float[][]{
{3420.8649f, 974.19067f},
{3365.2918f, 1025.4568f},
{3321.6456f, 1065.3801999999998f},
{3280.1989f, 1138.6715f},
{3290.7588f, 1232.7495f},
{3297.8264f, 1311.7838f},
{3290.6333f, 1364.3835f},
{3247.8828f, 1415.2549999999999f},
{3210.221f, 1443.3428999999999f},
{3179.1709f, 1469.9813f},
{3130.3084f, 1483.1535999999999f},
{3083.5733999999998f, 1470.2914999999998f},
{3023.0544999999997f, 1478.2632999999998f},
{2915.6317999999997f, 1583.4859999999999f},
{2839.0613999999996f, 1681.3422999999998f},
{2790.6595999999995f, 1737.5977999999998f},
{2753.7299999999996f, 1761.1674999999998f},
{2735.2762999999995f, 1770.3209999999997f},
{2681.2041999999997f, 1809.7014999999997f},
{2574.9014999999995f, 1862.9836999999998f},
{2485.7288999999996f, 1891.4622999999997f},
{2414.3756999999996f, 1905.7836999999997f},
{2341.0932f, 1936.7368999999997f},
{2323.4487999999997f, 1938.4220999999995f},
{2253.0354999999995f, 1945.1627999999996f},
{2186.4218999999994f, 1968.0487999999996f},
{2111.7298999999994f, 2014.8264999999997f},
{2052.5945999999994f, 2043.7188999999996f},
{1966.0857999999994f, 2117.7183999999997f},
{1919.5849999999994f, 2199.7313999999997f},
{1903.0166999999994f, 2218.7203999999997f},
{1862.5330999999994f, 2250.7164f},
{1818.1862999999994f, 2278.5247999999997f},
{1771.1272999999994f, 2311.1762f},
{1725.9442999999994f, 2375.1992f},
{1718.8306999999995f, 2426.1454f},
{1746.6830999999995f, 2479.4499f},
{1801.1739999999995f, 2593.8343f},
{1883.1580999999994f, 2676.1696f},
{2015.3864999999994f, 2763.8875000000003f},
{2099.416199999999f, 2725.9408000000003f},
{2117.9087999999992f, 2667.2672000000002f},
{2149.4006999999992f, 2622.7691000000004f},
{2235.7862999999993f, 2534.3888000000006f},
{2297.372999999999f, 2426.145400000001f},
{2319.768199999999f, 2347.7622000000006f},
{2346.862699999999f, 2260.5343000000007f},
{2387.7651999999994f, 2224.5241000000005f},
{2413.5090999999993f, 2180.2457000000004f},
{2456.3051999999993f, 2133.2007000000003f},
{2534.455499999999f, 2091.8176000000003f},
{2612.198699999999f, 2058.1442f},
{2680.013399999999f, 2056.5403f},
{2743.460099999999f, 2031.4394000000002f},
{2807.789299999999f, 2009.0170000000003f},
{2859.039199999999f, 1987.4050000000002f},
{2904.5964999999987f, 1955.9471f},
{2949.0498999999986f, 1904.1369000000002f},
{2991.6316999999985f, 1830.1510000000003f},
{3031.9342999999985f, 1769.1668000000002f},
{3054.1773999999987f, 1730.3344000000002f},
{3090.8378999999986f, 1695.7944000000002f},
{3149.2443999999987f, 1650.4187000000002f},
{3215.1996999999988f, 1607.4033000000002f},
{3285.823199999999f, 1562.5327000000002f},
{3329.3138999999987f, 1530.8836000000001f},
{3359.8392999999987f, 1485.4858000000002f},
{3383.139499999999f, 1442.0169f},
{3407.7007999999987f, 1385.6666f},
{3409.667399999999f, 1332.5137f},
{3381.2513999999987f, 1293.6947f},
{3333.150499999999f, 1220.5377f},
{3325.9127999999987f, 1150.1480000000001f},
{3333.751899999999f, 1097.3371000000002f},
{3382.3718999999987f, 1051.7293000000002f},
{3397.567599999999f, 1040.4342000000001f},
{3467.521599999999f, 990.9870600000002f},
{3460.392099999999f, 959.4735700000001f},
{3420.864899999999f, 974.1906700000001f},
{3420.8649f, 974.19067f}
};

float[][] arm6 = new float[][]{
{2883.7141f, 691.73371f},
{2913.1526000000003f, 731.46233f},
{2947.3865000000005f, 812.23582f},
{2974.8275000000003f, 867.81353f},
{2987.5068f, 936.97436f},
{2997.0602f, 994.3759100000001f},
{2987.7921f, 1058.7665000000002f},
{2967.2879000000003f, 1101.1619000000003f},
{2940.3390000000004f, 1158.0304000000003f},
{2899.7724000000003f, 1248.2455000000004f},
{2864.7178000000004f, 1291.4558000000004f},
{2843.215f, 1344.6129000000003f},
{2812.4623f, 1394.1004000000003f},
{2771.4256f, 1427.0782000000004f},
{2717.5603f, 1486.0411000000004f},
{2647.7372f, 1563.8250000000003f},
{2535.8665f, 1642.0589000000002f},
{2460.8738000000003f, 1701.5490000000002f},
{2456.4506f, 1707.7009000000003f},
{2466.4017000000003f, 1699.1103000000003f},
{2285.7319f, 1884.2755000000002f},
{2256.2234000000003f, 1931.6159000000002f},
{2228.3212000000003f, 1974.5091000000002f},
{2208.9962000000005f, 2083.05f},
{2198.5110000000004f, 2179.4489000000003f},
{2194.7912000000006f, 2272.5697000000005f},
{2195.6268000000005f, 2323.6011000000003f},
{2172.5462000000007f, 2449.7396000000003f},
{2141.8613000000005f, 2485.9483000000005f},
{2073.5176000000006f, 2524.8026000000004f},
{1992.3783000000005f, 2539.8505000000005f},
{1922.8558000000005f, 2493.7081000000003f},
{1846.8855000000005f, 2414.8862000000004f},
{1841.7387000000006f, 2347.5967000000005f},
{1876.4555000000005f, 2287.3650000000007f},
{1910.5539000000006f, 2237.3753000000006f},
{1938.3512000000005f, 2186.7086000000004f},
{1960.2815000000005f, 2132.3083000000006f},
{1984.6330000000005f, 2052.941500000001f},
{2013.2718000000004f, 2002.207500000001f},
{2046.0260000000005f, 1912.5032000000008f},
{2055.7292000000007f, 1834.1105000000007f},
{2071.2246000000005f, 1795.8757000000007f},
{2067.0132000000003f, 1795.6652000000008f},
{2109.1851f, 1732.2583000000009f},
{2137.6266f, 1679.001800000001f},
{2212.9631f, 1597.136700000001f},
{2276.8441f, 1558.3318000000008f},
{2327.8624999999997f, 1534.6415000000009f},
{2387.1829f, 1506.187400000001f},
{2443.7693f, 1476.470200000001f},
{2541.2990999999997f, 1425.1949000000009f},
{2584.4375999999997f, 1404.6308000000008f},
{2648.7353999999996f, 1357.7830000000008f},
{2691.9789999999994f, 1326.7753000000007f},
{2762.6622999999995f, 1276.2999000000007f},
{2808.7297999999996f, 1188.8112000000006f},
{2824.0466999999994f, 1140.4592000000005f},
{2870.6555999999996f, 1087.4633000000006f},
{2883.4390999999996f, 1036.4450000000006f},
{2903.1422f, 842.6152600000006f},
{2891.1800999999996f, 755.7582900000006f},
{2860.9852999999994f, 711.0472200000006f},
{2834.7945999999993f, 673.8434300000006f},
{2798.0464999999995f, 637.5212800000006f},
{2803.6715999999997f, 610.0445500000006f},
{2759.1256f, 570.0927500000006f},
{2719.4822999999997f, 511.0235300000006f},
{2697.4149999999995f, 463.03532000000064f},
{2697.4529999999995f, 395.19069000000064f},
{2739.6782999999996f, 306.0675700000006f},
{2769.5384999999997f, 265.0097200000006f},
{2844.745f, 196.4256000000006f},
{2852.1814999999997f, 201.0011500000006f},
{2846.9678999999996f, 202.7562600000006f},
{2785.0878999999995f, 275.8387500000006f},
{2730.3008999999997f, 352.13885000000056f},
{2718.8201999999997f, 414.14365000000055f},
{2717.3044999999997f, 470.38471000000055f},
{2759.0939f, 545.4360700000005f},
{2827.1382f, 601.6735600000005f},
{2868.8912f, 649.3650900000005f},
{2883.7144f, 691.7337100000004f},
{2883.7141f, 691.73371f}
};

float[][] arm5 = new float[][]{
    {2032.3633f, 675.58817f},
{2058.5995f, 634.70166f},
{2065.9561f, 540.2645299999999f},
{2065.0794f, 476.4379899999999f},
{2031.6503f, 436.0024899999999f},
{1987.5729000000001f, 391.9157899999999f},
{1952.9168000000002f, 347.6455099999999f},
{1942.7825000000003f, 328.4627599999999f},
{1932.3140000000003f, 298.6300899999999f},
{1925.2834000000003f, 299.13799999999986f},
{1960.3672000000001f, 340.77847999999983f},
{2003.1499000000001f, 384.80035999999984f},
{2045.4272f, 429.24109999999985f},
{2089.9903f, 471.3069399999998f},
{2112.6127f, 516.9555899999998f},
{2130.7863f, 625.4320199999997f},
{2125.8817000000004f, 687.7858699999997f},
{2074.1056000000003f, 791.4162999999998f},
{2048.3952000000004f, 847.6025099999997f},
{2009.9084000000005f, 919.8131399999997f},
{1986.1295000000005f, 978.9257099999998f},
{1975.4748000000004f, 1038.4995999999996f},
{1956.2908000000004f, 1117.4435999999996f},
{1954.2830000000004f, 1245.3267999999996f},
{1983.9687000000004f, 1344.0087999999996f},
{2007.2788000000003f, 1412.3190999999997f},
{2028.6860000000004f, 1481.7408999999998f},
{2099.6200000000003f, 1629.6544999999999f},
{2119.5293f, 1704.1763999999998f},
{2137.8764f, 1762.5348f},
{2142.6184000000003f, 1817.9037999999998f},
{2170.7659000000003f, 1896.2217999999998f},
{2196.2243000000003f, 1978.5792f},
{2205.3901000000005f, 2046.5521999999999f},
{2213.3911000000007f, 2101.4150999999997f},
{2232.8784000000005f, 2156.4669f},
{2241.7184000000007f, 2219.9186f},
{2228.3212000000008f, 2276.8441f},
{2185.1345000000006f, 2401.7221999999997f},
{2164.793100000001f, 2448.7688f},
{2146.596400000001f, 2510.0606f},
{2135.741500000001f, 2699.1097999999997f},
{2127.264200000001f, 2784.9503999999997f},
{2136.125300000001f, 2823.9674999999997f},
{2133.141600000001f, 2844.1888999999996f},
{2155.5368000000008f, 2944.9671999999996f},
{2153.9387000000006f, 3011.2524999999996f},
{2120.0778000000005f, 3088.6696999999995f},
{2069.5024000000003f, 3157.1738999999993f},
{1968.9417000000003f, 3211.888999999999f},
{1915.8558000000003f, 3192.0236999999993f},
{1911.0560000000003f, 3124.1286999999993f},
{1912.8699000000004f, 2990.979799999999f},
{1908.6669000000004f, 2928.777699999999f},
{1913.2830000000004f, 2807.5207999999993f},
{1924.6052000000004f, 2707.8199999999993f},
{1953.0788000000005f, 2591.754299999999f},
{1953.1498000000004f, 2604.935199999999f},
{1966.7626000000005f, 2546.780299999999f},
{1968.0316000000005f, 2484.278099999999f},
{1961.9353000000006f, 2410.426899999999f},
{1958.2179000000006f, 2373.463399999999f},
{1945.2421000000006f, 2287.8267999999994f},
{1934.4136000000005f, 2125.6614999999993f},
{1915.8951000000006f, 2069.4677999999994f},
{1900.9845000000007f, 1988.2151999999994f},
{1896.2563000000007f, 1983.3785999999993f},
{1882.1259000000007f, 1917.4022999999993f},
{1858.5320000000006f, 1860.4090999999994f},
{1835.1997000000006f, 1808.4993999999995f},
{1807.0868000000005f, 1750.8305999999995f},
{1788.2328000000005f, 1641.7164999999995f},
{1770.2413000000004f, 1581.6274999999996f},
{1756.3400000000004f, 1482.3194999999996f},
{1756.2990000000004f, 1376.3439999999996f},
{1765.3603000000005f, 1330.4529999999995f},
{1761.1541000000004f, 1303.5691999999995f},
{1788.3680000000004f, 1116.4816999999994f},
{1802.9991000000005f, 1043.8113999999994f},
{1829.0051000000005f, 966.0503399999993f},
{1851.8372000000006f, 918.9391499999994f},
{1888.4408000000005f, 865.9581299999993f},
{1962.4272000000005f, 769.8247099999993f},
{2017.4332000000006f, 692.3845599999993f},
{2032.3633000000007f, 675.5881699999993f},
{2032.3633f, 675.58817f}
};

float[][] arm4 = new float[][]{
{1564.8638f, 1290.5227f},
{1579.3826000000001f, 1220.2122f},
{1590.2010000000002f, 1161.6812f},
{1590.6014000000002f, 1094.6584f},
{1563.6610000000003f, 1047.5819000000001f},
{1532.5840000000003f, 989.4247600000001f},
{1514.6885000000002f, 950.2549500000001f},
{1483.0434000000002f, 908.3759800000001f},
{1438.9037000000003f, 867.1078100000002f},
{1425.4651000000003f, 846.6664100000002f},
{1393.0359000000003f, 802.3039800000001f},
{1364.9452000000003f, 759.7358000000002f},
{1338.3552000000004f, 696.5738100000001f},
{1313.9531000000004f, 641.29861f},
{1313.6839000000004f, 555.3318200000001f},
{1332.5196000000005f, 504.6250200000001f},
{1345.3073000000006f, 435.66353000000015f},
{1365.2514000000006f, 363.4934900000002f},
{1402.9106000000006f, 298.8398700000002f},
{1413.6962000000005f, 221.1524300000002f},
{1436.4692000000005f, 157.54089000000022f},
{1470.9402000000005f, 70.67433300000022f},
{1479.9812000000004f, 99.23408700000022f},
{1463.9717000000005f, 205.3347500000002f},
{1428.6875000000005f, 271.8243200000002f},
{1398.0139000000004f, 332.2867000000002f},
{1370.2206000000003f, 387.5247500000002f},
{1353.8147000000004f, 445.2830900000002f},
{1354.9089000000004f, 559.8796500000002f},
{1388.7226000000003f, 626.5078100000002f},
{1408.8611000000003f, 660.1686900000002f},
{1418.3427000000004f, 676.6703500000002f},
{1463.1288000000004f, 721.8798300000002f},
{1549.0005000000003f, 806.2267100000001f},
{1595.4201000000003f, 838.5981400000002f},
{1658.8349000000003f, 909.8665900000002f},
{1659.4115000000004f, 907.9278400000002f},
{1781.1896000000004f, 1133.4422000000002f},
{1784.1500000000003f, 1254.1305000000002f},
{1795.9361000000004f, 1337.4626000000003f},
{1772.9524000000004f, 1399.6992000000002f},
{1702.5732000000003f, 1492.2275000000002f},
{1694.5692000000004f, 1552.7330000000002f},
{1664.1287000000004f, 1647.9028f},
{1654.1234000000004f, 1705.2871f},
{1667.4693000000004f, 1768.9619f},
{1702.3952000000004f, 1843.4579f},
{1808.9594000000004f, 2022.8192000000001f},
{1873.4806000000003f, 2088.9974f},
{1922.4299000000003f, 2138.5798f},
{1993.3049000000003f, 2218.7518f},
{2045.4925000000003f, 2282.9094f},
{2086.5393000000004f, 2330.9235f},
{2112.9715000000006f, 2371.9308f},
{2131.2754000000004f, 2426.1454f},
{2153.4557000000004f, 2489.0348999999997f},
{2165.6072000000004f, 2559.9761f},
{2173.8955000000005f, 2581.0607f},
{2183.7716000000005f, 2657.7762f},
{2166.9957000000004f, 2718.653f},
{2146.2421000000004f, 2740.7583999999997f},
{2147.2003000000004f, 2937.8259f},
{2119.7802000000006f, 2991.9224f},
{2068.6716000000006f, 3024.8269999999998f},
{1995.2378000000006f, 3031.3477999999996f},
{1916.9472000000005f, 3013.6801999999993f},
{1898.3751000000004f, 2951.7489999999993f},
{1880.0571000000004f, 2918.0877999999993f},
{1857.3345000000004f, 2865.8467999999993f},
{1852.8670000000004f, 2799.9990999999995f},
{1844.1339000000005f, 2728.6983999999993f},
{1835.5842000000005f, 2719.372999999999f},
{1831.5542000000005f, 2665.844799999999f},
{1801.8047000000006f, 2605.927999999999f},
{1766.1129000000005f, 2541.826999999999f},
{1731.8945000000006f, 2485.865899999999f},
{1690.8367000000005f, 2429.877899999999f},
{1644.1217000000006f, 2385.158599999999f},
{1618.3095000000005f, 2341.9779999999987f},
{1593.7909000000004f, 2310.436899999999f},
{1556.7177000000004f, 2258.198199999999f},
{1522.6264000000003f, 2218.273599999999f},
{1452.0262000000002f, 2127.4323999999992f},
{1428.5173000000002f, 2077.578199999999f},
{1403.3087000000003f, 2010.368099999999f},
{1386.8709000000003f, 1942.7688999999991f},
{1343.7904000000003f, 1851.781799999999f},
{1323.3286000000003f, 1783.987499999999f},
{1317.6129000000003f, 1688.0224999999991f},
{1336.2462000000003f, 1631.116199999999f},
{1376.5543000000002f, 1580.551599999999f},
{1376.3709000000003f, 1572.328799999999f},
{1395.8969000000004f, 1549.5947999999992f},
{1436.6236000000004f, 1479.4641999999992f},
{1486.0655000000004f, 1392.6756999999993f},
{1529.1360000000004f, 1340.0249999999994f},
{1564.8638000000005f, 1290.5226999999993f},
{1564.8638f, 1290.5227f}
};

float[][] arm3 = new float[][]{
    {2144.1537f, 2089.1787f},
{2114.0557f, 2021.8881f},
{2046.8781999999999f, 1975.8610999999999f},
{1982.9935999999998f, 1928.1199f},
{1901.9516999999998f, 1892.8835f},
{1837.4393999999998f, 1841.3336f},
{1789.4441999999997f, 1797.3620999999998f},
{1773.3542999999997f, 1800.0735999999997f},
{1759.0338999999997f, 1767.5152999999998f},
{1686.5115999999996f, 1723.4616999999998f},
{1623.1683999999996f, 1675.9543999999999f},
{1568.8637999999996f, 1636.3370999999997f},
{1478.9548999999997f, 1592.5085999999997f},
{1420.0216999999998f, 1556.0073999999997f},
{1340.6240999999998f, 1491.8800999999996f},
{1269.8451999999997f, 1434.8804999999995f},
{1241.1988999999996f, 1369.7726999999995f},
{1234.1548999999995f, 1309.9384999999995f},
{1200.6871999999996f, 1216.1671999999994f},
{1227.8063999999997f, 1062.2736999999995f},
{1296.7716999999998f, 939.3373099999995f},
{1333.6792999999998f, 894.5623799999995f},
{1378.4934999999998f, 833.6957999999995f},
{1414.6637999999998f, 770.6750699999996f},
{1454.1738999999998f, 702.8871399999996f},
{1450.5538999999999f, 567.8632599999996f},
{1432.7522f, 514.4051899999996f},
{1404.9477f, 453.1918899999996f},
{1371.6698f, 396.5880799999996f},
{1354.5450999999998f, 324.8727599999996f},
{1332.8455999999999f, 256.01192999999955f},
{1340.0711f, 112.45435999999955f},
{1362.3072f, 91.96267699999956f},
{1363.2939999999999f, 210.05225999999956f},
{1376.7886999999998f, 273.8976099999996f},
{1417.3700999999999f, 345.18897999999956f},
{1462.8066f, 406.9696399999996f},
{1522.8751f, 509.3845499999996f},
{1552.683f, 609.3959999999996f},
{1560.2553f, 702.5480499999996f},
{1542.547f, 780.0185599999996f},
{1526.7859f, 838.6265199999997f},
{1493.9568000000002f, 914.1371499999997f},
{1447.8102000000001f, 1038.6351999999997f},
{1427.8603f, 1105.8658999999998f},
{1428.6485f, 1181.1730999999997f},
{1482.838f, 1308.9037999999998f},
{1517.5964999999999f, 1361.8778999999997f},
{1566.5469999999998f, 1411.8881999999996f},
{1612.8160999999998f, 1449.0710999999997f},
{1678.9959f, 1493.7291999999998f},
{1688.0358999999999f, 1502.2942999999998f},
{1736.6581999999999f, 1538.7108999999998f},
{1809.388f, 1593.3655999999999f},
{1908.2125999999998f, 1654.84f},
{1979.1966999999997f, 1690.9376f},
{2068.303f, 1744.9412f},
{2130.1791f, 1766.4037f},
{2182.3574999999996f, 1793.5563f},
{2226.5502999999994f, 1808.6362f},
{2280.3535999999995f, 1847.5086999999999f},
{2331.6246999999994f, 1886.3985999999998f},
{2388.5938999999994f, 1956.6594999999998f},
{2430.8314999999993f, 2022.1209999999996f},
{2450.1307999999995f, 2095.7699f},
{2435.6634999999997f, 2157.0912f},
{2405.0409999999997f, 2226.9550999999997f},
{2380.4411999999998f, 2259.6146f},
{2363.2655999999997f, 2304.4755f},
{2337.2666999999997f, 2375.2771f},
{2302.0562999999997f, 2434.5564f},
{2218.2225999999996f, 2563.9127f},
{2178.8585999999996f, 2637.77f},
{2162.9810999999995f, 2663.1787f},
{2163.6734999999994f, 2744.2837f},
{2151.5327999999995f, 2821.5692f},
{2177.4209999999994f, 2940.1782f},
{2159.3654999999994f, 3018.6039f},
{2125.7447999999995f, 3090.5215f},
{2093.0341999999996f, 3158.0335999999998f},
{2034.8988999999997f, 3182.9936f},
{1922.7770999999998f, 3187.52f},
{1891.1946999999998f, 3140.044f},
{1908.2125999999998f, 2998.2428f},
{1877.1282999999999f, 2921.24f},
{1860.0432999999998f, 2849.3835f},
{1828.4736999999998f, 2721.1673f},
{1834.5588999999998f, 2618.2325f},
{1870.8714999999997f, 2536.7719f},
{1903.6427999999996f, 2485.5849000000003f},
{1933.0935999999997f, 2406.6914f},
{1979.7381999999998f, 2346.4859f},
{2013.5235999999998f, 2293.9896000000003f},
{2066.3619f, 2251.9836000000005f},
{2119.3565f, 2190.6175000000003f},
{2142.8033f, 2132.5034000000005f},
{2144.1537f, 2089.1787000000004f},
{2144.1537f, 2089.1787f}
};

float[][] arm0 = new float[][]{
{85.84822f, 1319.4498f},
{122.94993f, 1374.1015f},
{158.02218f, 1463.8568f},
{196.04480999999998f, 1558.564f},
{216.0084f, 1612.0191f},
{226.3113f, 1679.9102f},
{238.072f, 1723.9208f},
{261.27719f, 1763.621f},
{287.40491000000003f, 1804.6789f},
{345.26819f, 1844.4376000000002f},
{444.53371f, 1834.6862f},
{502.05627f, 1843.3482000000001f},
{585.51654f, 1841.2397f},
{643.8616499999999f, 1845.7367000000002f},
{662.9441899999999f, 1853.0447000000001f},
{725.8151799999999f, 1880.6944f},
{768.73304f, 1909.9617f},
{809.99993f, 1924.0449f},
{862.7530499999999f, 1969.9180000000001f},
{915.1852999999999f, 2035.4781f},
{956.5248399999998f, 2084.5659f},
{989.9058499999999f, 2127.2253f},
{1033.9111999999998f, 2179.7983f},
{1081.7273999999998f, 2222.1468f},
{1122.5424999999998f, 2243.3097f},
{1158.0941999999998f, 2261.0478999999996f},
{1232.5323999999998f, 2284.8772999999997f},
{1292.4026f, 2300.8654999999994f},
{1377.3040999999998f, 2308.5705999999996f},
{1436.8771f, 2305.8495999999996f},
{1489.4787f, 2292.1127999999994f},
{1512.7371999999998f, 2287.8395999999993f},
{1595.7366999999997f, 2283.7246999999993f},
{1679.6457999999998f, 2271.2528999999995f},
{1772.0606999999998f, 2270.0033999999996f},
{1851.7425999999998f, 2260.4122999999995f},
{1896.6697f, 2278.1504999999993f},
{1955.8922f, 2302.361499999999f},
{1995.5684f, 2328.622699999999f},
{2049.4953f, 2390.287699999999f},
{2076.0382f, 2422.097099999999f},
{2093.3409f, 2462.754899999999f},
{2110.235f, 2553.526899999999f},
{2111.9897f, 2629.842399999999f},
{2094.9692f, 2727.512799999999f},
{2062.3713f, 2808.911699999999f},
{1958.9693999999997f, 2877.8716999999992f},
{1901.0438999999997f, 2903.9400999999993f},
{1862.7096999999997f, 2902.8928999999994f},
{1855.9623999999997f, 2846.0607999999993f},
{1844.1881999999996f, 2783.9479999999994f},
{1835.7534999999996f, 2721.845599999999f},
{1831.0723999999996f, 2661.4636999999993f},
{1804.7704999999996f, 2599.8662999999992f},
{1744.9583999999995f, 2554.917699999999f},
{1675.2344999999996f, 2526.887499999999f},
{1653.0450999999996f, 2511.164199999999f},
{1594.6728999999996f, 2492.5114999999987f},
{1478.4754999999996f, 2498.0325999999986f},
{1423.9606999999996f, 2519.4585999999986f},
{1365.2909999999997f, 2529.1298999999985f},
{1250.3979999999997f, 2543.7200999999986f},
{1176.0223999999996f, 2537.0837999999985f},
{1163.0255999999997f, 2535.6588999999985f},
{1090.8160999999998f, 2514.9783999999986f},
{1032.6852999999999f, 2491.9769999999985f},
{974.9400099999999f, 2456.1716999999985f},
{933.7035299999999f, 2415.3634999999986f},
{930.9925799999999f, 2414.6911999999984f},
{865.7710899999998f, 2348.7831999999985f},
{796.0562399999999f, 2255.7461999999987f},
{764.0716199999999f, 2214.2908999999986f},
{731.5761399999999f, 2140.6066999999985f},
{691.8814499999999f, 2090.1577999999986f},
{646.9534399999999f, 2032.0217999999986f},
{593.4724799999999f, 1995.0379999999986f},
{558.0134299999999f, 1972.6427999999987f},
{496.5945199999999f, 1962.6141999999986f},
{423.64229999999986f, 1959.5788999999986f},
{367.89364999999987f, 1949.0709999999985f},
{270.8243899999999f, 1903.9522999999986f},
{230.85140999999987f, 1860.6570999999985f},
{177.06022999999988f, 1780.9826999999984f},
{148.1777399999999f, 1696.5697999999984f},
{141.2822599999999f, 1652.6922999999983f},
{139.87999999999988f, 1579.1706999999983f},
{141.71503999999987f, 1496.4846999999984f},
{137.34513999999987f, 1454.4349999999984f},
{104.69472999999988f, 1386.0173999999984f},
{80.18746699999988f, 1338.4297999999983f},
{75.09184899999988f, 1317.2476999999983f},
{85.84821999999988f, 1319.4497999999983f},
{85.84822f, 1319.4498f}
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
