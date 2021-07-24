#ifdef GL_ES
precision mediump float;
#endif

#define TWO_PI 6.28318530718
#define PI 3.141592653589793
#define HALF_PI 1.5707963267948966

uniform vec2 u_resolution;
uniform vec2 u_mouse;
uniform float u_time;

vec3 rgb2hsb( in vec3 c ){
    vec4 K = vec4(0.0, -1.0 / 3.0, 2.0 / 3.0, -1.0);
    vec4 p = mix(vec4(c.bg, K.wz),
                 vec4(c.gb, K.xy),
                 step(c.b, c.g));
    vec4 q = mix(vec4(p.xyw, c.r),
                 vec4(c.r, p.yzx),
                 step(p.x, c.r));
    float d = q.x - min(q.w, q.y);
    float e = 1.0e-10;
    return vec3(abs(q.z + (q.w - q.y) / (6.0 * d + e)),
                d / (q.x + e),
                q.x);
}

//  Function from Iñigo Quiles
//  https://www.shadertoy.com/view/MsS3Wc
vec3 hsb2rgb( in vec3 c ){
    vec3 rgb = clamp(abs(mod(c.x*6.0+vec3(0.0,4.0,2.0),
                             6.0)-3.0)-1.0,
                     0.0,
                     1.0 );
    rgb = rgb*rgb*(3.0-2.0*rgb); //reassembling rbg using smoothstep()
    return c.z * mix(vec3(1.0), rgb, c.y);
}

float rect(in vec2 st, in vec2 corner, float w, float h){//draw rect with top left corner and width/height
    vec2 end = corner + vec2(w, h);
    vec2 b1 = vec2(st.x<corner.x || st.x > end.x ? 0. : 1., st.y<corner.y || st.y > end.y ? 0. : 1.);
    return b1.x*b1.y;
    
}


float rect(in vec2 st, in vec2 top_left, in vec2 bottom_right){//draw rect with two corners 
    vec2 b1 = vec2(st.x<top_left.x || st.x > bottom_right.x ? 0. : 1., st.y<top_left.y || st.y > bottom_right.y ? 0. : 1.);
    return b1.x*b1.y; 
}

float circle(in vec2 st, in vec2 center, float r){
    return 1.0-step(r, distance(st, center));
}

float rand(float x){
    return fract(sin(12.59585855*PI*x)+4102200.398383);
}

float elasticInOut(float t) {
  return t < 0.5
    ? 0.5 * sin(+13.0 * HALF_PI * 2.0 * t) * pow(2.0, 10.0 * (2.0 * t - 1.0))
    : 0.5 * sin(-13.0 * HALF_PI * ((2.0 * t - 1.0) + 1.0)) * pow(2.0, -10.0 * (2.0 * t - 1.0)) + 1.0;
}

float circle_outline(in vec2 st, in vec2 center, float r){
    return smoothstep(r-0.05, r, distance(st, center)) * (1.-smoothstep(r, r+0.05, distance(st, center)));
}

float random(vec2 st) {
    return fract(sin(dot(st.xy,vec2(12.9898,78.233)))*43758.5453123);
}

float plot(vec2 st, float fxn){//intake y value (prob defined in a separate function) and return based on pixels location relative to the function defined by the input y
    return smoothstep(fxn-0.01, fxn, st.y) // isolate values of y below the smoothstep() line 
    - smoothstep(fxn, fxn+0.01, st.y);  // isolate values of y above the smoothstep () line
    //return 1 for the the union if the set of points in which y is less than the upper bound (pct + 0.02) and greater than the upper bound (pct - 0.02))
}

void main(){
    vec2 st = gl_FragCoord.xy/u_resolution.xy;
    vec3 color = vec3(0.6549, 0.051, 0.051);

    vec2 pos = st - vec2(0.5);

    float r = length(pos)*2.0;
    float a = atan(pos.y,pos.x);

    float f = cos(a*3.);
    // f = abs(cos(a*3.));
    // f = abs(cos(a*2.5))*.5+.3;
    f = abs(cos(a*12.)*sin(a*3.))*.8+.1;
    // f = smoothstep(-.5,1., cos(a*10.))*0.2+0.5;

    color = vec3( 1.-smoothstep(f,f+0.02,r) );

    // color = vec3(1.-plot(st, f));

    gl_FragColor = vec4(color, 1.0);
}