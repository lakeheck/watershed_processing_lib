import processing.core.*; 
import processing.data.*; 
import processing.event.*; 
import processing.opengl.*; 

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

public void setup(){
    
    background(255);
    colorMode(HSB, 360, 100, 100, 100);
    doReset();
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

    
    // line = new ArrayList(); //generate a random line 
    // int n = 50;
    // for(int i=0; i<50; i++){
    //     line.add(new PVector(map(i, 0, n-1, 0, renderWidth), renderHeight/2 + map(compoundTrigFunction(map(i, 0, n-1, 0, 2*TWO_PI)), -3, 4, -50, 50)));
    // }

    // Ribbon r = new Ribbon(line);
    // render.beginDraw();
    // r.display();
    // ArrayList<PVector> points = r.generatePointsInside(500);
    // Gradient lineGrad = new Gradient(line_palette);
    // float colorVar = 0.1;

    // for(PVector p:points){
    //     ArrayList<PVector> knn = k_nearest_neighbors(p, points, 10);
    //     for(PVector k:knn){
    //              int baseColor = lineGrad.eval(map(k.y,0,renderHeight,0,1)+randomGaussian()*colorVar, HSB);
    //              render.stroke(hue(baseColor) + randomGaussian(), saturation(baseColor) + randomGaussian()*8, brightness(baseColor) + randomGaussian()*8);
    //              render.line(p.x, p.y, k.x, k.y);
    //     }
    //     // render.fill(0,100,100);
    //     // render.ellipseMode(CENTER);
    //     // render.ellipse(p.x, p.y, 5, 5); 
    // }
    
    // render.endDraw();

    as = new AttractorSystem();
    as.addPerlinFlowField(0.005f, 4, 0.5f, true);
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
    as.calculateAttractorSystem();
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

public void watercolorBackgroundTexture(int[] baseColors, int numPoints, int numLayers, float geometryWidth, float colorVar, float jitter){
  //5000 points and 5 layers each is a good balance (with jitter of 4)

  Gradient grad = new Gradient(baseColors);

  // render.blendMode(ADD);
  render.colorMode(HSB, 360, 100,100,100);
  render.noStroke();
  for(int i=0; i<numPoints; i++){ //number of locations we drop geometry at
      render.pushMatrix();
      float x = random(-renderWidth/4, renderWidth*5/4);
      float y = random(-renderHeight/4, renderHeight*5/4);
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

    render.loadPixels();//then we jitter the resulting pixels to add some grainy vibes 
    for(int i=0; i<render.pixels.length; i++){
        int c = render.pixels[i];
        render.pixels[i] = color(hue(c)+random(-jitter, jitter), saturation(c)+randomGaussian()*jitter, brightness(c) + randomGaussian()*jitter/2);
    }
    render.updatePixels();
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

        p.displayPath();
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
    pos = getTorusPosition(pos);
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

  public void displayPath(){
    if(path.size()==100){
      render.fill(0,20);
      Ribbon r = new Ribbon(path, 2);
      // render.beginDraw();
      r.display();
      // ArrayList<PVector> points = r.generatePointsInside(500);
      // Gradient lineGrad = new Gradient(line_palette);
      // float colorVar = 0.1;

      // for(PVector p:points){

          // ArrayList<PVector> knn = k_nearest_neighbors(p, points, 10);
          // for(PVector k:knn){
          //         int baseColor = lineGrad.eval(map(k.y,0,renderHeight,0,1)+randomGaussian()*colorVar, HSB);
          //         render.stroke(hue(baseColor) + randomGaussian(), saturation(baseColor) + randomGaussian()*8, brightness(baseColor) + randomGaussian()*8);
          //         render.line(p.x, p.y, k.x, k.y);
          // }
          // render.fill(0,100,100);
          // render.ellipseMode(CENTER);
          // render.ellipse(p.x, p.y, 5, 5); 
        
        // }
    } 
    //   for(int i=1; i<path.size(); i++){
    //     if(PVector.dist(path.get(i), path.get(i-1)) < STEP_SIZE*100){
    //       render.rectMode(CORNERS);
    //       render.fill(c,50);
    //       render.noStroke();
    //       // render.line(path.get(i-1).x, path.get(i-1).y, path.get(i).x, path.get(i).y);
    //       PVector dir = new PVector(0,0,1).cross(new PVector((path.get(i).x-path.get(i-1).x), (path.get(i).y - path.get(i-1).y))).normalize();
    //       PVector start = new PVector(path.get(i-1).x + 10*cos(dir.heading()), path.get(i-1).y + 10*sin(dir.heading())); 
    //       PVector end = new PVector(path.get(i).x - 10*cos(dir.heading()), path.get(i).y - 10*sin(dir.heading())); 
    //       render.rect(start.x, start.y, end.x, end.y);
    //       // render.rect(path.get(i-1).x, path.get(i-1).y + 10*sin(dir.heading()),path.get(i-1).x, path.get(i).y - 10*sin(dir.heading())); 
    //     }
    //   }
    // }
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

public float compoundTrigFunction(float x){ //allows compounding of trig functions for interesting lines 
    return sin(x) + 3*cos(x);
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


}

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
