
// #extension GL_OES_standard_derivatives : enable

#ifdef GL_ES
precision mediump float;
#endif

varying vec4 v_position;
varying vec4 v_normal;
varying vec2 v_texcoord;
varying vec4 v_color;

const int numOctaves=6;


#if defined(VERTEX)

// attribute vec4 a_position; // myfolder/myfile.obj
attribute vec4 a_position;
attribute vec4 a_normal;
attribute vec2 a_texcoord;
attribute vec4 a_color;

void main(void) {
	v_position = a_position;
	v_normal = a_normal;
	v_texcoord = a_texcoord;
	v_color = a_color;
	gl_Position = v_position;
}

#else // fragment shader

// uniform vec2 u_mouse;
// uniform vec2 u_pos;
// uniform sampler2D u_texture; // https://cdn.jsdelivr.net/gh/actarian/plausible-brdf-shader/textures/mars/4096x2048/diffuse.jpg?repeat=true
// uniform vec2 u_textureResolution;

// float checker(vec2 uv, float repeats) {
// 	float cx = floor(repeats * uv.x);
// 	float cy = floor(repeats * uv.y);
// 	float result = mod(cx + cy, 2.0);
// 	return sign(result);
// }


float rand(float n){return fract(sin(n) * 43758.5453123);}
float rand(vec2 n) { 
	return fract(sin(dot(n, vec2(12.9898, 4.1414))) * 43758.5453);
}

float noise(float p){
	float fl = floor(p);
  float fc = fract(p);
	return mix(rand(fl), rand(fl + 1.0), fc);
}
	

float noise(vec2 p){
	vec2 ip = floor(p);
	vec2 u = fract(p);
	u = u*u*(1.0-.0*u);
	
	float res = mix(
		mix(rand(ip),rand(ip+vec2(1.0,0.0)),u.x),
		mix(rand(ip+vec2(0.0,1.0)),rand(ip+vec2(1.0,1.0)),u.x),u.y);
	return res*res;
}

float fbm(vec2 x)
{    
    float G = 0.5; //exp2(-H);
    float f = 1.0;
    float a = 1.0;
    float t = 0.0;
    for( int i=0; i<numOctaves; i++ )
    {
        t += a*noise(f*x);
        f *= 2.0;
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
    for( int i=0; i<numOctaves; i++ )
    {
        t += a*noise(f*x);
        f *= 2.0;
        a *= G;
    }
    return t;
}


float pattern( in vec2 p )
{
    vec2 q = vec2( fbm( p + vec2(0.0,0.0) ),
                   fbm( p + vec2(0,0) ) );

    return fbm( p + 4.0*q );
}
void main() {
	vec2 p = vec2 (v_texcoord[0], v_texcoord[1]);

	// vec3 ambient = vec3(0.4);
	// vec3 direction = vec3(0.0, 1.0, 1.0);
	// vec3 lightColor = vec3(0.9686, 0.8078, 0.8078);
	// float incidence = max(dot(v_normal.xyz, direction), - 1.0);
	// vec3 light = clamp(ambient + lightColor * incidence, 0.0, 1.0);

	vec3 color = (v_normal.rgb+pattern(p));
	gl_FragColor = vec4(color, 1.0);
}

#endif
