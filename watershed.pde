
import java.util.Collections;
import java.util.Arrays;
import java.util.Comparator;
import megamu.mesh.*;
import processing.pdf.*;

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
float polygonAreaCutoff = 0.001;
PVector[][] noiseGrid; //set up a noise based flow field for changing tile attributes
Gradient colorGrad;
Polygon poly;

PShader shader;

void setup(){
    // size(750, 750, P3D);
    background(255);
    colorMode(HSB, 360, 100, 100, 100);
    size(1000, 1000, P3D);
    noStroke();

    shader = loadShader("distance_field.frag");
    // doReset();
    saveHighRes();
}
void savePDF() {
  println("Saving PDF image...");
  beginRecord(PDF,  "-vector.pdf");
  seededRender();
  endRecord();
  println("Finished");
}

void seededRender() {
  randomSeed(seed);
  noiseSeed(seed);
  render();
}

void render() {
  /* Do your drawing in here */
    shader.set("u_resolution", float(width), float(height));
    shader.set("u_mouse", float(mouseX), float(mouseY));
    shader.set("u_time", millis() / 1000.0);
    // shader(shader);
    rect(0,0,width,height);
    //   if(frameCount==1) 
    save(saveFilePath + "-1.png");
//   if(int(millis()/1000.0) ==15) save(saveFilePath + "-15.png");
}


ArrayList<PVector> inkscapePathImport(float[][] p, float inputWidth, float inputHeight){
    ArrayList<PVector> out = new ArrayList();
    for(int i=0; i<p.length; i++){
        // println(p[i][0][0]);
        float x = map(p[i][0], 0, inputWidth, 0, renderWidth);
        float y = map(p[i][1], 0, inputHeight, 0 ,renderHeight);
        out.add(new PVector(x, y));
    }
    return out;
}



void draw(){
}

void keyPressed(){
    switch(key){
        case 's': saveHighRes();
        break;
    }
}


color lerpColor(color[] arr, float step, int colorMode) {
  int sz = arr.length;
  if (sz == 1 || step <= 0.0) {
    return arr[0];
  } else if (step >= 1.0) {
    return arr[sz - 1];
  }
  float scl = step * (sz - 1);
  int i = int(scl);
  return lerpColor(arr[i], arr[i + 1], scl - i, colorMode);
}

void saveHighRes() {
int dpi = printDpi;
scaleFactor = dpi / (float)previewDpi;
  PGraphics hires = createGraphics(
                        int(width * scaleFactor),
                        int(height * scaleFactor),
                        P3D);
  println("Generating high-resolution image...");
  shader.set("u_resolution", float(hires.width), float(hires.height));
    // shader.set("u_mouse", float(mouseX), float(mouseY));
  shader.set("u_time", millis() / 1000.0);
  hires.beginDraw();
  hires.background(255);
  hires.stroke(0);
  hires.line(0,0,588,588);
  hires.noStroke();
//   hires.scale(scaleFactor);
  hires.filter(shader);
  hires.endDraw();

  hires.save(seed + "-highres.png");
  println("Finished");
}

void doReset() { //initial setup and used in resetting for high-def export
    int dpi = renderHighRes ? printDpi : previewDpi;
    scaleFactor = dpi / (float)previewDpi;
    renderWidth = printWidth * dpi;
    renderHeight = printHeight * dpi;
    render = createGraphics(renderWidth, renderHeight);
    shader=loadShader("distance_field.glsl");
    // shader.set("fraction", 1.0);
    firstFrame = true;
    noiseSeed(seed);
    randomSeed(seed);
    background_palette = new int[]{color(#0f0f0e), color(#382a04), color(#141524), color(#170d1f), color(#000000)};
    line_palette = new int[]{color(#382a04), color(#594a1f), color(#073610), color(#18361e), color(#243618), color(#313622), color(#473216)};
    // render.shader(shader);
    // render.endDraw();




}
