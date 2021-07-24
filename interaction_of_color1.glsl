// Author: Patricio Gonzalez Vivo

// Title: Interaction of color - IV
// Chapter: A color has many faces - the relativity of color

// "To begin the study of how color deceives and how to make use of this,
// the first excercise is
// to make one and the same color look different."
// 															Josef Albers


#ifdef GL_ES
precision mediump float;
#endif

uniform vec2 u_resolution;
uniform vec2 u_mouse;
uniform float u_time;

float rect(in vec2 st, in vec2 size){
	size = 0.25-size*0.25;
    size = vec2(0.24375, 0.2275);
    vec2 uv = smoothstep(size,size,st*(1.0-st)); //binary
	return uv.x*uv.y;
}

void main() {
    vec2 st = gl_FragCoord.xy/u_resolution.xy;
	st.y = 1.-st.y;
    
    vec3 influenced_color = vec3(0.745,0.696,0.529);
    
    vec3 influencing_color_A = vec3(0.114,0.599,0.780);
    vec3 influencing_color_a = vec3(0.065,0.066,0.290);
    vec3 influencing_color_b = vec3(0.865,0.842,0.162);
    vec3 influencing_color_B = vec3(0.980,0.603,0.086);
    
    vec3 color =mix(influencing_color_A,
                         influencing_color_a,
                         step(.3,st.y));
    
    color = mix(color,
               influenced_color,
               rect(abs((st-vec2(.0,.5))*vec2(1.,1.75)),vec2(.025,.09)));

    gl_FragColor = vec4(color,1.0);
}