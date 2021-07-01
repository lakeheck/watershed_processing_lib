
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

float compoundTrigFunction(float x){
    return sin(x) + 3*cos(x);
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

    
    line = new ArrayList();
    int n = 50;
    for(int i=0; i<50; i++){
        line.add(new PVector(map(i, 0, n-1, 0, renderWidth), renderHeight/2 + map(compoundTrigFunction(map(i, 0, n-1, 0, TWO_PI)), -3, 4, -50, 50)));
    }

    Ribbon r = new Ribbon(line);
    render.beginDraw();
    r.display();
    ArrayList<PVector> points = r.generatePointsInside(50);

    for(PVector p:points){
        render.fill(0,100,100);
        render.ellipseMode(CENTER);
        render.ellipse(p.x, p.y, 5, 5); 
    }
    // ArrayList<PVector> top = new ArrayList();
    // ArrayList<PVector> btm = new ArrayList();

    // for(int i=1; i<line.size(); i++){
    //     PVector norm = new PVector(0,0,1).cross(new PVector((line.get(i).x - line.get(i-1).x), (line.get(i).y - line.get(i-1).y))).normalize();
    //     float theta = norm.heading();
    //     top.add(new PVector(line.get(i).x + 50*cos(theta), line.get(i).y + 50*sin(theta)));
    //     btm.add(new PVector(line.get(i).x - 50*cos(theta), line.get(i).y - 50*sin(theta)));
    // }

    // ArrayList<PVector> verts = new ArrayList();


    // render.beginDraw();
    // render.colorMode(HSB, 360, 100, 100, 100);
    // render.beginShape();
    // for(int i=0; i<((top.size() + btm.size())); i++){
    //     render.vertex(
    //         i<top.size() ? top.get(i).x : btm.get(btm.size() -1 - (i-top.size())).x,
    //         i<top.size() ? top.get(i).y : btm.get(btm.size() -1 - (i-top.size())).y
    //         );
    //     verts.add(
    //         i<top.size() ? top.get(i).copy() : btm.get(btm.size()-1-(i-top.size())).copy()
    //     );
    // }
    // render.endShape(CLOSE);

    // PVector[] points = new PVector[50];
    // for(int i=0; i<500; i++){
    //     PVector p= new PVector(random(renderWidth), random(renderHeight));
    //     if(polyPoint(verts, p.x, p.y)){
    //         render.fill(0,100,100);
    //         render.ellipseMode(CENTER);
    //         render.ellipse(p.x, p.y, 5, 5); 
    //     }
    //     else{            
    //         render.fill(40,100,100);
    //         render.ellipseMode(CENTER);
    //         render.ellipse(p.x, p.y, 5, 5); 
    //     }
    // }
    render.endDraw();

    // as = new AttractorSystem();
    // as.addPerlinFlowField(0.005, 4, 0.5, true);
    // as.addPerlinFlowField(0.01, 8, 0.9, false);


}

class Ribbon{ //class for drawing a ribbon based on a guide line (as used in flow fields, etc)
    ArrayList<PVector> vertices;

    Ribbon(ArrayList<PVector> v){ //should be initialized with ordered set of points 
        vertices = new ArrayList();

        
        ArrayList<PVector> top = new ArrayList();
        ArrayList<PVector> btm = new ArrayList();
        for(int i=1; i<line.size(); i++){
            PVector norm = new PVector(0,0,1).cross(new PVector((line.get(i).x - line.get(i-1).x), (line.get(i).y - line.get(i-1).y))).normalize();
            float theta = norm.heading();
            top.add(new PVector(line.get(i).x + 50*cos(theta), line.get(i).y + 50*sin(theta)));
            btm.add(new PVector(line.get(i).x - 50*cos(theta), line.get(i).y - 50*sin(theta)));
        }
    
        for(int i=0; i<((top.size() + btm.size())); i++){ // unwrap the top and bottom arrays - first we add all the top points, then start fro the end of hte bottom array to maintain non-self intersection
            vertices.add(
                i<top.size() ? top.get(i).copy() : btm.get(btm.size()-1-(i-top.size())).copy()
            );
        }
    }

    boolean contains(PVector point){
        return polyPoint(vertices, point.x, point.y);
    }

    ArrayList<PVector> generatePointsInside(int n){
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
