
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

void setup(){
    size(750, 750);
    background(255);
    colorMode(HSB, 360, 100, 100, 100);
    doReset();
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
    as.addPerlinFlowField(0.005, 4, 0.5, true);
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
