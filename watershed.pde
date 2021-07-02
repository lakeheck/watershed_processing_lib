
import java.util.Collections;
import java.util.Arrays;
import java.util.Comparator;
import megamu.mesh.*;

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


void setup(){
    size(750, 750);
    background(255);
    colorMode(HSB, 360, 100, 100, 100);
    doReset();
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


void doReset() { //initial setup and used in resetting for high-def export

    int dpi = renderHighRes ? printDpi : previewDpi;
    scaleFactor = dpi / (float)previewDpi;
    renderWidth = printWidth * dpi;
    renderHeight = printHeight * dpi;
    render = createGraphics(renderWidth, renderHeight);
    firstFrame = true;
    noiseSeed(seed);
    randomSeed(seed);
    background_palette = new int[]{color(#0f0f0e), color(#382a04), color(#141524), color(#170d1f), color(#000000)};
    line_palette = new int[]{color(#382a04), color(#594a1f), color(#073610), color(#18361e), color(#243618), color(#313622), color(#473216)};

    float[][][] import_paths = new float[][][]{
        branch5, branch4, branch3, branch2, branch1, arm0, arm1, arm3, arm4, arm5, arm6, arm7, arm8
    };

    
    ArrayList<ArrayList<PVector>> paths = new ArrayList<ArrayList<PVector>>(); 
    for(float[][] p:import_paths){
        paths.add(inkscapePathImport(p, 3564.00000, 5014.66650));
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
    int sample_size = int(t.length*0.5);
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


void draw(){
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
