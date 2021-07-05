


/* ************************* GENERAL UTILITIES *************************************/

PVector getTorusPosition (PVector position) {
  PVector pos = position.copy();
  if (pos.x < 0) pos.x = renderWidth + pos.x;
  if (pos.x > renderWidth) pos.x %= renderWidth;
  if (pos.y < 0) pos.y = renderHeight + pos.y;
  if (pos.y > renderHeight) pos.y = pos.y %= renderHeight;
  return pos;
}

/* ************************ INTERSECTION AND COLLISION TESTS ******************************/
boolean lineLine(float x1, float y1, float x2, float y2, float x3, float y3, float x4, float y4) {

  // calculate the direction of the lines
  float uA = ((x4-x3)*(y1-y3) - (y4-y3)*(x1-x3)) / ((y4-y3)*(x2-x1) - (x4-x3)*(y2-y1));
  float uB = ((x2-x1)*(y1-y3) - (y2-y1)*(x1-x3)) / ((y4-y3)*(x2-x1) - (x4-x3)*(y2-y1));

  // if uA and uB are between 0-1, lines are colliding
  if (uA >= 0 && uA <= 1 && uB >= 0 && uB <= 1) {
    return true;
  }
  return false;
}

boolean polyPoint(ArrayList<PVector> vertices, float px, float py) {
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

boolean polyLine(ArrayList<PVector> vertices, float x1, float y1, float x2, float y2) {

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

float sigmoidEasing(float x){
  return x==0 ? 0 : 1/(1+exp(-x));
}

float easing(float x){
  return sigmoidEasing(x);
}


/************************* POISSON DISK SAMPLING ****************************/

boolean isValidPoint(PVector[][] grid, float cellsize,
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

void insertPoint(PVector[][] grid, float cellsize, PVector point) {
  int xindex = floor(point.x / cellsize);
  int yindex = floor(point.y / cellsize);
  grid[xindex][yindex] = point;
}

ArrayList<PVector> poissonDiskSampling(float radius, int k) {
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
    int random_index = int(random(active.size()));
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


void watercolorBackgroundTexture(int[] baseColors, ArrayList<PVector> points, int numLayers, float geometryWidth, float colorVar, float jitter){
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

void restricted_chaikin_overlay(){
    render.smooth();
    render.noFill();

    float x = random(render.width);
    float y = random(render.height);
    float j = 2.5;
    PShape s = render.createShape();
    float border = 100;
    s.beginShape();
    for (int i = 0; i < 50000; i++) {
        s.vertex(x, y);
        float qx = random(1) < 0.5 ? -1 : 1;
        float qy = random(1) < 0.5 ? -1 : 1;
        x += qx * j;
        y += qy * j;
        x = constrain(x, -border, render.width + border);
        y = constrain(y, -border, render.height + border);
    }
    s.endShape();
    render.shape(chaikin_open(s, 0.25, 3), 0, 0);
}

void gridline(float x1, float y1, float x2, float y2) {
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
    render.strokeWeight(1 + map(noise(sx, sy), 0, 1, -0.5, 0.5));
    render.line(sx, sy, x + map(noise(x, y), 0, 1, -1, 1), y + map(noise(x, y), 0, 1, -1, 1));
    sx = x;
    sy = y;
  }
}

void canvas_overlay_example1() {
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
  color c;
  float target_number;
  float dt;
  String attribute;
  int alpha;
  boolean disp=true;
  int age=0;
  int pressure_lag = 50;
  int attractorForce = 1000;
  boolean normalize_attractor_force = true; //if normalize_attractor_force == true, kill_Particles should also = true
  float damp = 0.0002; //applied to particle velocity every frame. set to 1 pill not dampen at all


  
  Attractor(float _x, float _y, float _z){
    pos = new PVector(_x, _y, _z);
    c = color(0);
    float initial_pressure = 5000;
    pressure = new ArrayList<Float>();
    for(int i=0; i<pressure_lag; i++){pressure.add(initial_pressure);}
    force = attractorForce;
    dt = 0.001;
    }
  

  Attractor(){
      c = color(0,0,100);
  }
  
  void calc_force_on(Particle p){
    PVector f = PVector.sub(pos, p.pos);
    float d = f.mag();
    p.vel.add( new PVector(0,0));
    
  }
    
  void update_force(){
    force = (347.3058 + (6147.582 - 347.3058)/pow(1 + (pressure.get(0)/7004.265),8.499077)) + noise(5000,10000)*random(-1,1);
  }

  void display(){age++;}
  
}

class LorenzAttractor extends Attractor{
  
  LorenzAttractor(float _x, float _y, float _z){
    super(_x, _y, _z);
  }
  
  void calc_force_on(Particle p){
    PVector f = PVector.sub(pos, p.pos);
    float d = f.mag();
    if(normalize_attractor_force){
      f.normalize();
      float dx = (10*(f.x + f.y))*dt;
      float dy = (-f.x*f.z + 28*f.x-f.y)*dt;
      float dz = (f.x*f.y - 8.0/3*f.z)*dt;
      // p.vel.add( new PVector(dx, dy, dz).mult(force/(d+0.0001)));
      // float factor = map(d, 0, sqrt(width*height), 0, 1);
      // p.c = color(c);
      // if(random(1)<0.0002){p.DIVISION_CHANCES+=0.1;}
      // if(random(1)<0.00005){p.DEPOSIT_RATE*=max(randomGaussian()+10, 0);}

    }
    else{
      float dx = (10*(f.x + f.y))*dt;
      float dy = (-f.x*f.z + 28*f.x-f.y)*dt;
      float dz = (f.x*f.y - 8.0/3*f.z)*dt;
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
                float angleClamp = map(easing(map(x,0,renderWidth,0,1)),easing(0.),easing(1.),1, PI/2); // (PI/10); 
                // float angleClamp = map(x, 0, renderWidth, 1, PI/2);
                // float angle = floor(map(n, 0 , 1, -PI, PI)/angleClamp)*angleClamp; //(randomGaussian()*angleClamp/10+angleClamp);
                float angle = map(n, 0 , 1, -PI, PI); //for a 'normal' flow field
                grid[x][y] = new PVector(angle, n);
                yoff += (renderHighRes && scaleNoiseScale) ? noiseScale/(printDpi/previewDpi) : noiseScale;;
            }
            xoff += (renderHighRes && scaleNoiseScale) ? noiseScale/(printDpi/previewDpi) : noiseScale; //noise works best with step size of 0.005
        } 
    }

    void calc_force_on(Particle p){
        if(random(1)<0.0005){p.DEPOSIT_RATE*=max(randomGaussian()+5, 0);}
        p.vel.add( new PVector(cos(grid[constrain(ceil(p.pos.x),0,renderWidth-1)][constrain(ceil(p.pos.y),0,renderHeight-1)].x),sin(grid[constrain(ceil(p.pos.x),0,renderWidth-1)][constrain(ceil(p.pos.y),0,renderHeight-1)].x)).normalize());

    }

    void display(){
      
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
  float damp = 0.0002; //applied to particle velocity every frame. set to 1 pill not dampen at all
  boolean normalize_attractor_force = true; //if normalize_attractor_force == true, kill_Particles should also = true
  boolean allow_Particle_overlap=false; //will determine if particles die phen they run into each other 
  float STEP_SIZE = 1;
  boolean DISCRETE_DIV_ANGLE = false;
  //noise parameters 
  // int noiseOctaves = 8//defines numbers of octaves. takes int value from 0-8. higher = more high-level structure
  // float noiseGain = 0.5; //defines how much of each octave carries over into the next (0,1). Higher = more small details
  // float noiseScale = 0.005; // 0.005 is 'best' but 0.1 or slightly higher can have more sweeping patterns 
  //initial particle settings
  float initial_TURN_CHANCES = 0.;
  float initial_TURN_ANGLE = PI/8;
  float initial_DEPOSIT_RATE = 0.01;
  float  initial_DIVISION_CHANCES = 0.00;
  float initial_DIVISION_ANGLE = PI / 8;
  float initial_TERMINATION_THRESHOLD = 0.7;
  float initial_TERMINATION_CHANCES = initial_DIVISION_CHANCES * 0.;
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

  void addPerlinFlowField(float _noiseScale, int _noiseOctaves, float _noiseGain, boolean _scaleNoiseScale){
    attractors.add(new PerlinFlowField(_noiseScale, _noiseOctaves, _noiseGain, _scaleNoiseScale));
  }

  void calculateAttractorSystem(){//logic that updates particles in the system
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

        p.displayPath(int(randomGaussian()*50+200), randomGaussian()*2+8);
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
        int idxLast = int(p.pos.x) + int(p.pos.y) * renderWidth; //check to ensure particle moved to a new grid location on this step
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
  color c;
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
  float damp=0.0002;
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
    age = randomGaussian()*(renderHeight*0.26667)/2+(renderHeight*0.26667);
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
      int baseColor = xGrad.eval(map(pos.x,0,renderWidth,0,1)+randomGaussian()*0.1, HSB);
      c = color(hue(baseColor) + randomGaussian(), saturation(baseColor) + randomGaussian()*8, brightness(baseColor) + randomGaussian()*8);
      }
    else{
      // bright = randomGaussian()*2+20;
      int baseColor = yGrad.eval(map(pos.y,0,renderHeight,0,1)+randomGaussian()*0.1, HSB);
      c = color(hue(baseColor) + randomGaussian(), saturation(baseColor) + randomGaussian()*8, brightness(baseColor) + randomGaussian()*8);
    }
    age = randomGaussian()*(renderHeight*0.26667)/2+(renderHeight*0.26667);
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
      int baseColor = xGrad.eval(map(pos.x,0,renderWidth,0,1)+randomGaussian()*0.1, HSB);
      c = color(hue(baseColor) + randomGaussian(), saturation(baseColor) + randomGaussian()*8, brightness(baseColor) + randomGaussian()*8);
      }
    else{
      // bright = randomGaussian()*2+20;
      int baseColor = yGrad.eval(map(pos.y,0,renderHeight,0,1)+randomGaussian()*0.1, HSB);
      c = color(hue(baseColor) + randomGaussian(), saturation(baseColor) + randomGaussian()*8, brightness(baseColor) + randomGaussian()*8);
    }
    // c = color(0, 0, 100);
    age = randomGaussian()*(renderHeight*0.26667)/2+(renderHeight*0.26667);
    boxBounds = bounds;
  }
  
  void update () {
    lastPos = pos.copy();
    // vel = PVector.sub(pos, lastPos);
    // the Particle has a random chances to turn
    if (random(0, 1) < TURN_CHANCES) {
      this.ang+= TURN_ANGLE * (round(random(0, 1)) * 2. - 1.);
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
  
  void display () {
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
    render.stroke(0, int(100*DEPOSIT_RATE));    
    PVector line = lastPos.copy().sub(pos);
    if (line.mag() < 4*STEP_SIZE) {
      render.line(lastPos.x, lastPos.y, pos.x, pos.y);
      int idx = ceil(pos.x) + ceil(pos.y) * renderWidth;
      if (idx < renderWidth*renderHeight){attractorSystem.grid[idx] += DEPOSIT_RATE;} //update substrate grid when drawing Particle 
    }
  }

  void displayPath(int path_length, float stroke_weight){
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

void keyPressed() {
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

void chaikin_line(ArrayList<PVector> points, int layers, float layer_var, String mode){
  PShape n = render.createShape();
  n.beginShape();
  for(int i=0; i<points.size(); i++){
    n.vertex(points.get(i).x, points.get(i).y);
  }
  n.endShape();
  
  if(mode == "CLOSE"){
    render.shape(chaikin_close(n, 0.25, 5), 0, 0);
  
    for(int i=0; i<=layers; i++){
      render.shape(chaikin_close(perturbPShape(n, layer_var, layer_var), .25, 5),0,0);
    }
  }
  else{
    render.shape(chaikin_open(n, 0.25, 5), 0, 0);
  
    for(int i=0; i<=layers; i++){
      render.shape(chaikin_open(perturbPShape(n, layer_var, layer_var), .25, 5),0,0);
    }
  }
}

PShape perturbPShape(PShape s, float x_std, float y_std){
  PShape n = render.createShape();
  n.beginShape();
  for(int i=0; i<s.getVertexCount(); i++){
    n.vertex(s.getVertex(i).x+randomGaussian()*x_std, s.getVertex(i).y+randomGaussian()*y_std);
  }
  n.endShape();
  return n;
}

ArrayList<PVector> chaikin_cut(PVector a, PVector b, float ratio) {
  float x, y;
  ArrayList<PVector> n = new ArrayList<PVector>();

  /*
   * If ratio is greater than 0.5 flip it so we avoid cutting across
   * the midpoint of the line.
   */
   if (ratio > 0.5) ratio = 1 - ratio;

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

PShape chaikin(PShape shape, float ratio, int iterations, boolean close) {
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

PShape chaikin_close(PShape original, float ratio, int iterations) {
  return chaikin(original, ratio, iterations, true);
}

PShape chaikin_open(PShape original, float ratio, int iterations) {
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
  static final float TOLERANCE = 0.09;
  float percent; //this is just used for sorting in the gradient class, could be thought of as a heirarchy with arbitrary scale
  int clr;

  ColorStop(int colorMode, float percent, float[] arr) {
    this(colorMode, percent, arr[0], arr[1], arr[2],
      arr.length == 4 ? arr[3] : 1.0);
  }

  ColorStop(int colorMode, float percent, float x, float y, float z, float w) {
    this(percent, colorMode == HSB ? composeclr(hsbToRgb(x, y, z, w))
      : composeclr(x, y, z, w));
  }

  ColorStop(float percent, int clr) {
    this.percent = constrain(percent, 0.0, 1.0);
    this.clr = clr;
  }

  boolean approxPercent(ColorStop cs, float tolerance) {
    return abs(percent - cs.percent) < tolerance;
  }

  // Mandated by the interface Comparable<ColorStop>.
  // Permits color stops to be sorted by Collections.sort.
  int compareTo(ColorStop cs) {
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
    float szf = sz <= 1.0 ? 1.0 : sz - 1.0;
    for (int i = 0; i < sz; ++i) {
      colorStops.add(new ColorStop(i / szf, colors[i]));
    }
  }

  // Creates equidistant color stops.
  Gradient(int colorMode, float[]... colors) {
    int sz = colors.length;
    float szf = sz <= 1.0 ? 1.0 : sz - 1.0;
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

  void add(int colorMode, float percent, float[] arr) {
    add(new ColorStop(colorMode, percent, arr));
  }

  void add(int colorMode, float percent,
    float x, float y, float z, float w) {
    add(new ColorStop(colorMode, percent, x, y, z, w));
  }

  void add(final float percent, final int clr) {
    add(new ColorStop(percent, clr));
  }

  void add(final ColorStop colorStop) {
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

  int eval(final float step) {
    return eval(step, DEFAULT_COLOR_MODE);
  }

  int eval(final float step, final int colorMode) {
    int sz = colorStops.size();

    // Exit from the function early whenever possible.
    if (sz == 0) {
      return 0x00000000;
    } else if (sz == 1 || step < 0.0) {
      return colorStops.get(0).clr;
    } else if (step >= 1.0) {
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

  boolean remove(ColorStop colorStop) {
    return colorStops.remove(colorStop);
  }

  ColorStop remove(int i) {
    return colorStops.remove(i);
  }

  int remove() {
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
float[] smootherStepRgb(float[][] arr, float st, float[] out) {
  int sz = arr.length;
  if (sz == 1 || st < 0) {
    out = java.util.Arrays.copyOf(arr[0], 0);
    return out;
  } else if (st > 1) {
    out = java.util.Arrays.copyOf(arr[sz - 1], 0);
    return out;
  }
  float scl = st * (sz - 1);
  int i = int(scl);
  float eval = smootherStep(scl - i);
  out[0] = arr[i][0] + eval * (arr[i + 1][0] - arr[i][0]);
  out[1] = arr[i][1] + eval * (arr[i + 1][1] - arr[i][1]);
  out[2] = arr[i][2] + eval * (arr[i + 1][2] - arr[i][2]);
  out[3] = arr[i][3] + eval * (arr[i + 1][3] - arr[i][3]);
  return out;
}

float[] smootherStepHsb(float[] a, float[] b, float st, float[] out) {

  // Find difference in hues.
  float huea = a[0];
  float hueb = b[0];
  float delta = hueb - huea;

  // Prefer shortest distance.
  if (delta < -0.5) {
    hueb += 1.0;
  } else if (delta > 0.5) {
    huea += 1.0;
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

float smootherStep(float st) {
  return st * st * st * (st * (st * 6.0 - 15.0) + 10.0);
}

float[] smootherStepRgb(float[] a, float[] b, float st, float[] out) {
  float eval = smootherStep(st);
  out[0] = a[0] + eval * (b[0] - a[0]);
  out[1] = a[1] + eval * (b[1] - a[1]);
  out[2] = a[2] + eval * (b[2] - a[2]);
  out[3] = a[3] + eval * (b[3] - a[3]);
  return out;
}

//general utlities to work with bit representations of color 
int composeclr(float[] in) {
  return composeclr(in[0], in[1], in[2], in[3]);
}

// Assumes that RGBA are in range 0 .. 1.
int composeclr(float red, float green, float blue, float alpha) {
  return round(alpha * 255.0) << 24
    | round(red * 255.0) << 16
    | round(green * 255.0) << 8
    | round(blue * 255.0);
}

float[] decomposeclr(int clr) {
  return decomposeclr(clr, new float[] { 0.0, 0.0, 0.0, 1.0 });
}

// Assumes that out has 4 elements.
// 1.0 / 255.0 = 0.003921569
float[] decomposeclr(int clr, float[] out) {
  out[3] = (clr >> 24 & 0xff) * 0.003921569;
  out[0] = (clr >> 16 & 0xff) * 0.003921569;
  out[1] = (clr >> 8 & 0xff) * 0.003921569;
  out[2] = (clr & 0xff) * 0.003921569;
  return out;
}


//HSB to RBG fxns
float[] hsbToRgb(float[] in) {
  float[] out = new float[] { 0.0, 0.0, 0.0, 1.0 };
  return hsbToRgb(in[0], in[1], in[2], in[3], out);
}

float[] hsbToRgb(float[] in, float[] out) {
  if (in.length == 3) {
    return hsbToRgb(in[0], in[1], in[2], 1.0, out);
  } else if (in.length == 4) {
    return hsbToRgb(in[0], in[1], in[2], in[3], out);
  }
  return out;
}

float[] hsbToRgb(float hue, float sat, float bri, float alpha) {
  float[] out = new float[] { 0.0, 0.0, 0.0, 1.0 };
  return hsbToRgb(hue, sat, bri, alpha, out);
}

float[] hsbToRgb(float hue, float sat, float bri, float alpha, float[] out) {
  if (sat == 0.0) {

    // 0.0 saturation is grayscale, so all values are equal.
    out[0] = out[1] = out[2] = bri;
  } else {

    // Divide color wheel into 6 sectors.
    // Scale up hue to 6, convert to sector index.
    float h = hue * 6.0;
    int sector = int(h);

    // Depending on the sector, three tints will
    // be distributed among R, G, B channels.
    float tint1 = bri * (1.0 - sat);
    float tint2 = bri * (1.0 - sat * (h - sector));
    float tint3 = bri * (1.0 - sat * (1.0 + sector - h));

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
float[] rgbToHsb(int clr) {
  return rgbToHsb(clr, new float[] { 0.0, 0.0, 0.0, 1.0 });
}

float[] rgbToHsb(int clr, float[] out) {
  return rgbToHsb((clr >> 16 & 0xff) * 0.003921569,
    (clr >> 8 & 0xff) * 0.003921569,
    (clr & 0xff) * 0.003921569,
    (clr >> 24 & 0xff) * 0.003921569, out);
}

float[] rgbToHsb(float[] in) {
  return rgbToHsb(in, new float[] { 0.0, 0.0, 0.0, 1.0 });
}

float[] rgbToHsb(float[] in, float[] out) {
  if (in.length == 3) {
    return rgbToHsb(in[0], in[1], in[2], 1.0, out);
  } else if (in.length == 4) {
    return rgbToHsb(in[0], in[1], in[2], in[3], out);
  }
  return out;
}

float[] rgbToHsb(float red, float green, float blue, float alpha, float[] out) {

  // Find highest and lowest values.
  float max = max(red, green, blue);
  float min = min(red, green, blue);

  // Find the difference between max and min.
  float delta = max - min;

  // Calculate hue.
  float hue = 0.0;
  if (delta != 0.0) {
    if (red == max) {
      hue = (green - blue) / delta;
    } else if (green == max) {
      hue = 2.0 + (blue - red) / delta;
    } else {
      hue = 4.0 + (red - green) / delta;
    }

    hue /= 6.0;
    if (hue < 0.0) {
      hue += 1.0;
    }
  }

  out[0] = hue;
  out[1] = max == 0.0 ? 0.0 : (max - min) / max;
  out[2] = max;
  out[3] = alpha;
  return out;
}

//**********************************AD HOC ******************************************

ArrayList<PVector> k_nearest_neighbors(PVector p, ArrayList<PVector> points, int k){ //originally used in kenny vaden sketch
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
    color c;

    Circle(PVector p, float _r){
        center = p;
        r=_r;
        c = color(0,100,0);
    }

    boolean isIn(PVector p){
        if(PVector.dist(p, center)<r){
            return true;
        }
        else{return false;}
    }

}

float compoundTrigFunction(float x, int choice){ //allows compounding of trig functions for interesting lines 
    float coeff1 = 3; //random(2,5);
    float coeff2 = 4; //random(2,5);
    float coeff3 = 5; //random(2,6);
    float coeff4 = 4;// (4,5);
    switch(choice){
      case 0: return cos(x*coeff1+coeff2) - coeff3*sin(randomGaussian()*0.01 + x) + cos(coeff4*x)*pow(sin(pow(x,2)), 2);
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

    ArrayList<PVector> generateBoundsfromFrameLine(ArrayList<PVector> guideLine, float stroke_weight){
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

    boolean contains(PVector point){
        return polyPoint(vertices, point.x, point.y);
    }

    ArrayList<PVector> generatePointsInside(int n){
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

    float regionArea(){
        float area = 0.0;
        int n = vertices.size();
        // Calculate value of shoelace formula
        int j = n - 1;
        for (int i = 0; i < n; i++){
          area += (vertices.get(j).x + vertices.get(i).x) * (vertices.get(j).y - vertices.get(i).y);
            
          // j is previous vertex to i
          j = i;
        }
        return Math.abs(area / 2.0);

    }

    ArrayList<PVector> generatePointsInside(ArrayList<PVector> region, int n){
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

    ArrayList<PVector> generate_N_points_about_curve(ArrayList<PVector> line, int n, float var){
        ArrayList<PVector> points = new ArrayList();
        int count = 0;
        ArrayList<PVector> region = generateBoundsfromFrameLine(line, 20);
        points.addAll(this.generatePointsInside(region, n));
        return points;
    }

    void display(){
        render.beginShape();
        for(int i=0; i<vertices.size(); i++){
            render.vertex(
                vertices.get(i).x,
                vertices.get(i).y
                );
        }
        render.endShape(CLOSE);
    }

    void vadenWeb(int n, int _knn, Gradient grad, boolean allow_intersection){
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
              else if(!polyLine(this.vertices, p.x, p.y, k.x, k.y) || polyLine(inkscapePathImport(hole_in_arm_8, 3564.00000, 5014.66650), p.x, p.y, k.x, k.y)){
                render.stroke(0);
                render.line(p.x, p.y, k.x, k.y);
              }
            }
      }
    }


}

//**************************** POLYGON AND SUBDIVISION *****************************
boolean continueSubdivision(float threshold){
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
  color c = color(0,0,0);
  color nc = color(0,100,100);;
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

  PVector getRandomPoint(float _pct, float std){ //inputs are pct [0,1] and stdev [0,0.5]
    float pct = randomGaussian()*std+_pct; //calc pct based on normal distro
    PVector p = new PVector( //use random pct to pull a point 
        map(pct, 0, 1, p1.x, p2.x),
        map(pct, 0, 1, p1.y, p2.y)
      );
    return p;
  }

  void display(){
    stroke(c);
    line(p1.x, p1.y, p2.x, p2.y);
    
    if(showNormals){
      stroke(nc);
      int normal_length = 25;
      line(center.x, center.y, center.x+normal_length*normal.x, center.y+normal_length*normal.y); 
      ellipse(center.x+normal_length*normal.x, center.y+normal_length*normal.y, normal_length*0.2,normal_length*0.2);
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
  
  void subdivide(){
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
        float r = 0.9;
        point.sub(shrinkDirection.mult(PVector.sub(point, p.centroid).mag()*(1-r)));
      }
    }
  }

  boolean continueSubdivision(float threshold){
    boolean end = false;
    for(Polygon p:this.polygons){
      if(p.area > (renderWidth *renderHeight)*threshold){
        end = true;
        break;
      }
    }
    return end;
  }

  void display(){
    for(Polygon p: polygons){
        p.display();
    }
  }
}

class Polygon{
  ArrayList<PVector> points;
  ArrayList<Face> edges;
  color c;
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
      pct = 0.5;
      std = 0.1;
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
      pct = 0.5;
      std = 0.1;
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
      pct = 0.5;
      std = 0.1;
      convexHull();
      area = polygonArea();
      centroid = findCentroid();
      c = color(0,0,100); //colorGrad.eval(map(centroid.x, 0, renderWidth, 0, 1));
      geometricSubdivision = _g;
    }


  void convexHull(){ //implementing using the Mesh library
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

  PVector findCentroid(){
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

  void display(){
    render.fill(c,50);
    // render.noStroke();
    render.beginShape();
    for(int i=0; i<hull.size(); i++){
      render.vertex(hull.get(i).x, hull.get(i).y);
    }
    render.endShape(CLOSE);

    // chaikin_line(this.hull, 3, 0.5, "CLOSE");

  }

  ArrayList<ArrayList<PVector>> split(int numDivisions){
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
      int new_idx = int(random(edges.size()-1)); //choose another face randomly (could be improved with some better logic for eligible / preffered faces)
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
    int firstPointIndex = int(random(hull.size())); //first pick random point on the Polygon
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
      int childIndex = int(random(subdivisions.size()));
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

  void subdivide(){
    geometricSubdivision.subdivide();
  }

  float polygonArea(){
    float area = 0.0;
    int n = hull.size();
    // Calculate value of shoelace formula
    int j = n - 1;
    for (int i = 0; i < n; i++){
      area += (hull.get(j).x + hull.get(i).x) * (hull.get(j).y - hull.get(i).y);
        
      // j is previous vertex to i
      j = i;
    }

    // Return absolute value
    return Math.abs(area / 2.0);
    
  }

  void createChildren(){ 

    if(this.alive && this.area > (renderWidth*renderHeight)*polygonAreaCutoff){
      ArrayList<ArrayList<PVector>> subdivisions = this.split(this.maxChildren);
      int idx = int(random(subdivisions.size()));
      for(int i=0; i<subdivisions.size(); i++){
        Polygon p = new Polygon(subdivisions.get(i), this.geometricSubdivision);
        p.pct = this.pct;
        p.std = this.std;
        p.maxChildren = this.maxChildren;
        p.gen = this.gen+1;
        // p.c = (random(1)<map(noiseGrid[constrain(ceil(p.centroid.x),0,renderWidth-1)][constrain(ceil(p.centroid.y),0,renderHeight-1)].y, 0, 1, 0, 0.2))? palette[int(random(palette.length))] : color(hue(c), saturation(c) + 1, brightness(c) - 1); //this could be driven off some noise or other flow field
        p.alive = (idx==i && random(1) < 0.5 && p.gen > 1) ? false : true;
        this.geometricSubdivision.newPolygons.add(p);
        children.add(p);
      }
    }
    else{
      this.geometricSubdivision.deadPolygons.add(this);
    }
  }

  void displayChildren(){
    for(Polygon c:children){
      if(c.alive){c.display();}
    }
  }
  
}
