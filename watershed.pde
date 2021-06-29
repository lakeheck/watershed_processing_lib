
/* *********************** HIGH RES EXPORT FUNCTIONALITY **************************/

//INCLUDE THESE GLOBAL VARIABLES  
PGraphics render;
PImage img;
String saveFilePath = "../outputs/kenny_vaden_output-" + new java.text.SimpleDateFormat("yyyyMMdd-HHmmss").format(new java.util.Date()); //change the XXXX to current project name 
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
int[] background_palette, line_palette ;

//attractor globals

//poisson disk sampling for box spacing 
int poissonDiskRadius;
ArrayList<PVector> plist;
ArrayList<PVector> lorenzPoisson;
int box_alpha=5;
int layers=0;
int layerFrameSpace = 100;

//initialize
float grid[];
ArrayList<Particle> particles;
ArrayList<Attractor> attractors;
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

 
    as = new AttractorSystem();

    grid = new float[int(renderWidth*renderHeight)];


    // render.beginDraw();
    // render.colorMode(HSB, 360, 100,100,100);
    // render.strokeWeight((0.5 * (renderHighRes ? printDpi/previewDpi : 1 )));
    // render.background(0,0,100);
    // // render.blendMode(ADD);
    // render.ellipseMode(CENTER);
    // // ANY DRAWING LOGIC TO EXECUTE ON STARTUP/RESET GOES HERE


    // int numPanels = 4;
    // float panel_spacing = 10;
    // float panelWidth = renderWidth/numPanels;
    // int numPoints = 600;
    // int k_nearest_neighbors = 30;//number of neighbors to check for 
    // int numCircles = 1;
    // float circleMargin = 0.25;
    // ArrayList<Circle> circles = new ArrayList(); //initialize some ellipses to skew coloring 
    // // render.blendMode(DARKEST);
    // // watercolorBackgroundTexture(background_palette, 5000, 5, 50, 0.1, 4);
    // render.blendMode(BLEND);
    // for(int i=0; i<numCircles; i++){
    //     // circles.add(new Circle(new PVector(random(-renderWidth*circleMargin, renderWidth*circleMargin), random(renderHeight*-circleMargin, renderHeight*circleMargin)), randomGaussian()*renderWidth/100 + renderWidth/20));
    //     Circle c = new Circle(new PVector(random(renderWidth*circleMargin, renderWidth*(1-circleMargin)), random(renderHeight*circleMargin, renderHeight*(1-circleMargin))), randomGaussian()*renderWidth/30 + renderWidth/10);
    //     render.fill(0, 0, 100, 2);
    //     // float noise = randomGaussian()*c.r/10;
    //     // float noise = 0;
    //     // render.beginShape();
    //     for(int j=0; j<50; j++){
    //         float theta = map(j, 0, 49, 0, TWO_PI);
    //         float new_r = (c.r + randomGaussian()*c.r/50);
    //         // render.curveVertex(c.center.x + new_r * cos(theta), c.center.y + new_r * sin(theta));
    //         render.ellipse(c.center.x + randomGaussian()*c.r/30, c.center.y+randomGaussian()*c.r/30, new_r*2, new_r*2);
    //     }
    //     // render.endShape(CLOSE);
    //     circles.add(c);
    // }

    // ArrayList<ArrayList<PVector>> panels = new ArrayList(); 

    // for(int i = 0; i<numPanels; i++){
    //     float startX = map(i, 0, numPanels, 0, renderWidth);
    //     float endX = startX + panelWidth;
    //     float startY = panel_spacing; 
    //     float endY = renderHeight - panel_spacing;
        
    //     startX += ((i==0) ? 2*panel_spacing : panel_spacing);
    //     endX -= ((i==numPanels) ? 2*panel_spacing : panel_spacing);
    //     ArrayList<PVector> tempPoints = new ArrayList();
    //     for(int j=0; j<numPoints; j++){
    //         tempPoints.add(new PVector(random(startX, endX), random(startY, endY)));
    //     }

    //     panels.add(tempPoints);
    // }

    // Gradient lineGrad = new Gradient(line_palette);
    // float colorVar = 0.1;
    // render.blendMode(BLEND);
    // for(ArrayList<PVector> p:panels){
    //     for(int i=0; i<p.size(); i++){
    //         ArrayList<PVector> knn = k_nearest_neighbors(p.get(i), p, k_nearest_neighbors);
    //         for(PVector k:knn){
    //             int baseColor = lineGrad.eval(map(k.y,0,renderHeight,0,1)+randomGaussian()*colorVar, HSB);
    //             render.stroke(hue(baseColor) + randomGaussian(), saturation(baseColor) + randomGaussian()*8, brightness(baseColor) + randomGaussian()*8);

    //             // render.stroke(0);
    //             for(Circle c: circles){
    //                 if(c.isIn(k)){
    //                     render.stroke(c.c);
    //                 }
    //             }
    //             render.line(p.get(i).x, p.get(i).y, k.x, k.y);
    //         }
    //     }
    // }

    // render.endDraw();

}
 
ArrayList<PVector> k_nearest_neighbors(PVector p, ArrayList<PVector> points, int k){
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

class Circle{
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


void draw(){
    render.beginDraw();
    if(firstFrame){
        firstFrame = false;
        render.colorMode(HSB, 360,100,100,100);
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
