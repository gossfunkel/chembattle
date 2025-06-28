#version 150


uniform mat4 p3d_ModelViewProjectionMatrix;

// vertex inputs
in vec4 p3d_Vertex;
in vec2 p3d_MultiTexCoord0;

// to fragment
out vec2 texcoord;

void main () {
	gl_Position = p3d_ModelViewProjectionMatrix * p3d_Vertex;
	texcoord = p3d_MultiTexCoord0;
}

/* E/dt = -D*V+∇×M       
Q.x += -b.w*b.x-e.z+w.z;
Q.y += -b.w*b.y-n.z+s.z; */
/*    
V/dt = -D/m*E-D/m*V×M
Q.x += Q.w/2000.*(b.x-Q.x*b.z);
Q.y += Q.w/2000.*(b.y-Q.y*b.z);*/