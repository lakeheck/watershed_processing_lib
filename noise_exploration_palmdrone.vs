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


float noise(float p){
	float fl = floor(p);
  float fc = fract(p);
	return mix(rand(fl), rand(fl + 1.0), fc);
}

//below fuctions from https://github.com/ashima/webgl-noise via Patricio Gonzalez

vec3 mod289(vec3 x) { return x - floor(x * (1.0 / 289.0)) * 289.0; }
vec2 mod289(vec2 x) { return x - floor(x * (1.0 / 289.0)) * 289.0; }
vec3 permute(vec3 x) { return mod289(((x*34.0)+1.0)*x); }

// vec3 permute(vec3 x) { return mod(((x*34.0)+1.0)*x, 289.0); }

float snoise(vec2 v) {

    // Precompute values for skewed triangular grid
    const vec4 C = vec4(0.211324865405187,
                        // (3.0-sqrt(3.0))/6.0
                        0.366025403784439,
                        // 0.5*(sqrt(3.0)-1.0)
                        -0.577350269189626,
                        // -1.0 + 2.0 * C.x
                        0.024390243902439);
                        // 1.0 / 41.0

    // First corner (x0)
    vec2 i  = floor(v + dot(v, C.yy));
    vec2 x0 = v - i + dot(i, C.xx);

    // Other two corners (x1, x2)
    vec2 i1 = vec2(0.0);
    i1 = (x0.x > x0.y)? vec2(1.0, 0.0):vec2(0.0, 1.0);
    vec2 x1 = x0.xy + C.xx - i1;
    vec2 x2 = x0.xy + C.zz;

    // Do some permutations to avoid
    // truncation effects in permutation
    i = mod289(i);
    vec3 p = permute(
            permute( i.y + vec3(0.0, i1.y, 1.0))
                + i.x + vec3(0.0, i1.x, 1.0 ));

    vec3 m = max(0.5 - vec3(
                        dot(x0,x0),
                        dot(x1,x1),
                        dot(x2,x2)
                        ), 0.0);

    m = m*m ;
    m = m*m ;

    // Gradients:
    //  41 pts uniformly over a line, mapped onto a diamond
    //  The ring size 17*17 = 289 is close to a multiple
    //      of 41 (41*7 = 287)

    vec3 x = 2.0 * fract(p * C.www) - 1.0;
    vec3 h = abs(x) - 0.5;
    vec3 ox = floor(x + 0.5);
    vec3 a0 = x - ox;

    // Normalise gradients implicitly by scaling m
    // Approximation of: m *= inversesqrt(a0*a0 + h*h);
    m *= 1.79284291400159 - 0.85373472095314 * (a0*a0+h*h);

    // Compute final noise value at P
    vec3 g = vec3(0.0);
    g.x  = a0.x  * x0.x  + h.x  * x0.y;
    g.yz = a0.yz * vec2(x1.x,x2.x) + h.yz * vec2(x1.y,x2.y);
    return 130.0 * dot(m, g);
}
	
float noise(vec2 p){
	vec2 ip = floor(p);
	vec2 u = fract(p);
	u = u*u*(3.0-2.*u);
	
	float res = mix(
		mix(random(ip),random(ip+vec2(1.0,0.0)),u.x),
		mix(random(ip+vec2(0.0,1.0)),random(ip+vec2(1.0,1.0)),u.x),u.y);
	return res*res;
}

float fbm(vec2 x)
{    
    float G = 0.5; //exp2(-H);
    float f = 1.0;
    float a = 1.0;
    float t = 0.0;
    for( int i=0; i<8; i++ )
    {
        t += a*snoise(f*x);
        f = pow(2.0, float(i));
        // f += 2.0;
        a *= G;
    }
    return t;
}


// float fbm(vec2 x)
// {    
//     float G = 0.6; //exp2(-H);
//     float f = 1.0;
//     float a = 1.0;
//     float t = 0.0;
//     for( int i=0; i<8; i++ )
//     {
//         t += a*snoise(f*x);
//         f = pow(3.0, float(i));
//         // f += 2.0;
//         a *= G;
//     }
//     return t;
// }

float fbm(float x)
{    
    float G = 0.5; //exp2(-H);
    float f = 1.0;
    float a = 1.0;
    float t = 0.0;
    for( int i=0; i<8; i++ )
    {
        t += a*snoise(vec2(f*x));
        f *= pow(2.,float(i));
        f*= 2.0;
        a *= G;
    }
    return t;
}


float power_noise(vec2 st, float n){
    return pow(snoise(st),n);
}


vec3 random3(vec3 c) {
	float j = 4096.0*sin(dot(c,vec3(17.0, 59.4, 15.0)));
	vec3 r;
	r.z = fract(512.0*j);
	j *= .125;
	r.x = fract(512.0*j);
	j *= .125;
	r.y = fract(512.0*j);
	return r-0.5;
}

const float F3 =  0.3333333;
const float G3 =  0.1666667;
float snoise(vec3 p) {

	vec3 s = floor(p + dot(p, vec3(F3)));
	vec3 x = p - s + dot(s, vec3(G3));
	 
	vec3 e = step(vec3(0.0), x - x.yzx);
	vec3 i1 = e*(1.0 - e.zxy);
	vec3 i2 = 1.0 - e.zxy*(1.0 - e);
	 	
	vec3 x1 = x - i1 + G3;
	vec3 x2 = x - i2 + 2.0*G3;
	vec3 x3 = x - 1.0 + 3.0*G3;
	 
	vec4 w, d;
	 
	w.x = dot(x, x);
	w.y = dot(x1, x1);
	w.z = dot(x2, x2);
	w.w = dot(x3, x3);
	 
	w = max(0.6 - w, 0.0);
	 
	d.x = dot(random3(s), x);
	d.y = dot(random3(s + i1), x1);
	d.z = dot(random3(s + i2), x2);
	d.w = dot(random3(s + 1.0), x3);
	 
	w *= w;
	w *= w;
	d *= w;
	 
	return dot(d, vec4(52.0));
}

float distanceFunction(vec2 p){
    return 1.;
}

float angleFunction(vec2 p){
    return snoise(p);
}

vec2 domain_warp(vec2 p, float max_distance){
    float angle = angleFunction(p)*TWO_PI;
    float dist = distanceFunction(p)*max_distance;
    return p + vec2(cos(angle), sin(angle)) * dist;
}



float snoiseFractal(vec3 m) {
	return   0.5333333* snoise(m)
				+0.2666667* snoise(2.0*m)
				+0.1333333* snoise(4.0*m)
				+0.0666667* snoise(8.0*m);
}




void main(){
    vec2 st = gl_FragCoord.xy/u_resolution.xy;
    st.x *= u_resolution.x/u_resolution.y;

    st *= 10.;

    vec3 color = vec3(0.3412, 0.0392, 0.0);
    color = vec3(snoise(domain_warp(st, 1.5))+0.4); //, snoise(st));
    color += vec3(0.9686, 0.0275, 0.0275)*abs(fbm(st))*0.9;

    gl_FragColor = vec4(color, 1.0);
}