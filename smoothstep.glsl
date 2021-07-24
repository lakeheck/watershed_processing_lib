#ifdef GL_ES
precision mediump float;
#endif

#define PI 3.14159265359

uniform vec2 u_resolution;
uniform vec2 u_mouse;
uniform float u_time;

float plot(vec2 st, float fxn){//intake y value (prob defined in a separate function) and return based on pixels location relative to the function defined by the input y
    return smoothstep(fxn-0.1, fxn, st.y) // isolate values of y below the smoothstep() line 
    - smoothstep(fxn, fxn+0., st.y);  // isolate values of y above the smoothstep () line
    //return 1 for the the union if the set of points in which y is less than the upper bound (pct + 0.02) and greater than the upper bound (pct - 0.02))
}

void main(){
    vec2 st = gl_FragCoord.xy/u_resolution; //normalize to resolution so that x can be properly interpolated

    float y = smoothstep(0., 1., st.x); //interpolate x coordinate smoothly from 

    vec3 color = vec3(y); //equal values for all color channels based on y's value (i.e. i [0,1] so a greyscale value)

    float pct = plot(st, y); //return 0 if y is outside the line defined in plot, 1 if inside 
    
    //update colors: 
    //if the pixel is on the band defined in plot(), add green to its color >> this is accomplished with pct 
    //otherwise, keep the original color
    color = (1.-pct)*color+pct*vec3(0.,1.,0.);

    gl_FragColor = vec4(color,1.0); //update colors


}