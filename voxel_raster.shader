shader_type spatial;

uniform sampler2D u_data;

varying vec3 v_cam_pos_local;
varying vec3 v_vertex;
varying flat mat3 v_mvp;

void vertex() {
	v_cam_pos_local = ((inverse(WORLD_MATRIX) * CAMERA_MATRIX) * vec4(0,0,0,1)).xyz;
	v_vertex = VERTEX;
	v_mvp = mat3(MODELVIEW_MATRIX);
}

vec4 get_voxel(vec3 pos) {
	vec2 ts = vec2(16, 16);
	vec2 tc = vec2(4, 4);
	vec3 res = vec3(16);
	pos = floor(clamp(pos, vec3(0), res));
	//pos = floor(mod(pos, res)); // Infinite repeat
	float row = floor(pos.y / tc.x);
	vec2 pixel_pos = vec2(
		pos.x + pos.y * ts.x - row * ts.x * tc.x,
		pos.z + row * ts.y
	);
	vec2 texcoord = pixel_pos / vec2(ts * tc);
	return texture(u_data, texcoord);
}

// For debugging and toying
//float noise(vec2 co) {
//	return fract(sin(dot(co.xy, vec2(12.9898,78.233))) * 43758.5453);
//}
//float bnoise(vec2 co) {
//	return floor(0.5 * noise(co) + 0.5 + 0.01);
//}
//vec4 randcol(vec3 pos) {
//	return vec4(bnoise(pos.xy), bnoise(pos.yz), bnoise(pos.xz), 1.0);
//}
//vec4 get_sphere(vec3 pos) {
//	vec3 center = vec3(8);
//	float radius = 5.0;
//	float d = distance(pos, center);
//	if(floor(d - radius) < 0.0)
//		return vec4(1);
//	return vec4(0);
//}

//vec4 raytrace_dumb(vec3 rpos, vec3 rstep) {
//	vec4 p = vec4(0);
//	for (int i = 0; i < 32; ++i) {
//		p = get_voxel(rpos);
//		if (p.r + p.g + p.b != 0.0) {
//			break;
//		}
//		rpos += rstep;
//	}
//	return p;
//}

vec4 raytrace(vec3 p_rpos, vec3 p_rdir, out vec3 p_normal) {
	vec3 hit_pos = floor(p_rpos);
	vec3 hit_prev_pos = hit_pos;
	
	vec3 istep = sign(p_rdir);
	
	vec3 tdelta = 1.0 / (abs(p_rdir) + 0.000001);
	vec3 tcross = tdelta * mix(p_rpos - floor(p_rpos), ceil(p_rpos) - p_rpos, 0.5 * istep + 0.5);
	
	vec4 col = vec4(0);
	// Iterations should be sqrt(3) * res
	int iterations = 28;
	for (int i = 0; i < iterations; ++i) {
		hit_prev_pos = hit_pos;
		
		// Advance depending on which coordinate is closer to next integer
		// Non-if version is harder to understand but faster (by about 20%)
		float cxy = step(tcross.x, tcross.y);
		float cxz = step(tcross.x, tcross.z);
		float cyz = step(tcross.y, tcross.z);
		vec3 m;
		m.x = cxy * cxz;
		m.y = (1.0 - cxy) * cyz;
		m.z = 1.0 - (m.x + m.y);
		hit_pos += istep * m;
		tcross += tdelta * m;
		
		vec4 v = get_voxel(hit_pos);
		if (v.r + v.g + v.b != 0.0) {
			col = v;
			//col *= exp(-5.0*float(i) / float(iterations)); // Exp fog
			break;
		}
	}
	
	p_normal = hit_prev_pos - hit_pos;
	
	return col;
}

void fragment() {
	vec3 cam_dir = normalize(v_vertex - v_cam_pos_local);
	vec3 pos = v_vertex + vec3(8) - 0.1 * cam_dir;
	vec3 n;
	vec4 col = raytrace(pos, cam_dir, n);
	ALBEDO = col.rgb;
	NORMAL = v_mvp * n;
}
