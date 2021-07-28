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

float rand(vec2 n) { 
	return fract(sin(dot(n, vec2(12.9898, 4.1414))) * 43758.5453);
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



float noise(float p){
	float fl = floor(p);
  float fc = fract(p);
	return mix(rand(fl), rand(fl + 1.0), fc);
}
	
float noise(vec2 p){
	vec2 ip = floor(p);
	vec2 u = fract(p);
	u = u*u*(3.0-2.*u);
	
	float res = mix(
		mix(rand(ip),rand(ip+vec2(1.0,0.0)),u.x),
		mix(rand(ip+vec2(0.0,1.0)),rand(ip+vec2(1.0,1.0)),u.x),u.y);
	return 1.3*res*res*res*res;
}

float fbm(vec2 x)
{    
    float G = 0.5; //exp2(-H);
    float f = 1.0;
    float a = 1.0;
    float t = 0.0;
    for( int i=0; i<8; i++ )
    {
        t += a*noise(f*x);
        f *= 4.0;
        a *= G;
    }
    return t;
}

float fbm(float x)
{    
    float G = 0.5; //exp2(-H);
    float f = 1.0;
    float a = 1.0;
    float t = 0.0;
    for( int i=0; i<8; i++ )
    {
        t += a*noise(f*x);
        f *= 2.0+ u_time*1.;
        a *= G;
    }
    return t;
}


float pattern( in vec2 p, out vec2 q, out vec2 r)
{
    q.x = fbm( p + vec2(0.0,0.0) );
    q.y = fbm( p + vec2(9.2,5.3) );

    r.x = fbm( p + 4.0*q + vec2(1.7,9.2));
    r.y = fbm( p + 4.0*q + vec2(8.3,2.8));

    return fbm( p + 4.0*r);
}

const mat2 mtx = mat2( 0.80,  0.60, -0.60,  0.80 );

float fbm4( vec2 p )
{
    float f = 0.0;

    f += 0.5000*(-1.0+2.0*noise( p )); p = mtx*p*2.02;
    f += 0.2500*(-1.0+2.0*noise( p )); p = mtx*p*2.03;
    f += 0.1250*(-1.0+2.0*noise( p )); p = mtx*p*2.01;
    f += 0.0625*(-1.0+2.0*noise( p ));

    return f/0.9375;
}


float fbm6( vec2 p )
{
    float f = 0.0;

    f += 0.500000*noise( p ); p = mtx*p*2.02;
    f += 0.250000*noise( p ); p = mtx*p*2.03;
    f += 0.125000*noise( p ); p = mtx*p*2.01;
    f += 0.062500*noise( p ); p = mtx*p*2.04;
    f += 0.031250*noise( p ); p = mtx*p*2.01;
    f += 0.015625*noise( p );

    return f/0.96875;
}

vec2 fbm4_2( vec2 p )
{
    return vec2( fbm4(p+vec2(16.0)), fbm4(p+vec2(66.2)) );
}

vec2 fbm6_2( vec2 p )
{
    return vec2( fbm6(p+vec2(33.2)), fbm6(p+vec2(56.7)) );
}


float func( vec2 q, out vec2 o, out vec2 n )
{
    q += 0.05*sin(vec2(0.11,0.13)*u_time + length( q )*4.0);
    
    q *= 0.7 + 0.2*cos(0.05*u_time);

    o = 0.5 + 0.5*fbm4_2( q );
    
    o += 0.02*sin(vec2(0.11,06.136)*u_time*length( o ));

    n = fbm6_2( 4.0*o );

    vec2 p = q + 2.0*n + 1.0;

    float f = 0.5 + 0.5*fbm4( 2.0*p );

    f = mix( f, f*f*f*3.5, f*abs(n.x) );

    f *= 1.0-0.5*pow( 0.5+0.5*sin(86.0*p.x)*sin(8.0*p.y), 8.0 );

    return f;
}

void main(){
    vec2 st = gl_FragCoord.xy/u_resolution;
    vec3 color1 = vec3(1.0, 1.0, 1.0);
    vec3 color2 = vec3(0.1529, 0.1922, 0.1176);
    // vec3 color2 = hsb2rgb(vec3(110./360., .49, .16));
    vec3 color3 = vec3(0.2667, 0.2039, 0.0902);
    float pct;
    vec2 point = vec2(0.7, 0.7);

    vec2 q, r;

    float field = smoothstep(0.9, 0., pow(distance(st,point),distance(st,vec2(0.6)))); //,distance(st,vec2(0.6))));
    
    pct = pattern(st, q ,r);
    vec3 col = vec3(0.2,0.1,0.4);
    col = mix( col, vec3(0.3,0.05,0.05), pct );
    col = mix( col, vec3(1.0, 1.0, 1.0), dot(r,r) );
    col = mix( col, vec3(0.3922, 0.1608, 0.1608), 0.5*q.y*q.y - u_time*0.1);
    col = mix( col, vec3(0.0,0.2,0.4), 0.5*smoothstep(1.2,1.3,abs(r.y)+abs(r.x)+ min(u_time*0.1, 1.5)));
    col *= pct*2.0;
    col += field;
    // color2 += mix(color3, mix(color1, color2, vec3(pct)), pct);
    // color = hsb2rgb(vec3(pct, 0.5, 1.)); 


    gl_FragColor = vec4(col,0.9);
}