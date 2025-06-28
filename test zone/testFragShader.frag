#version 150
// test script that swaps red and green channels

// pass value into shader from main code
//uniform vec2 resolution;
uniform sampler2D p3d_Texture0;

// input from vertex
in vec2 texcoord;

// out to screen
out vec4 p3d_FragColor;

void main() {
	vec4 color = texture(p3d_Texture0, texcoord);
	p3d_FragColor = color.gbra;
	//p3d_FragColor = vec4(0.0,0.0,1.0,1.0);
}