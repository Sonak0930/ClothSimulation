//
//  main.cpp
//  IK
//
//  Created by Hyun Joon Shin on 2021/06/09.
//

#include <iostream>
#include <JGL/JGL_Window.hpp>
#include "AnimView.hpp"
#include <glm/gtx/quaternion.hpp>
#include <glm/gtc/random.hpp>
#include <Eigen/dense>
#include <stdlib.h>

AnimView* animView;


using namespace std;
using namespace glm;

//damping
const float DEF_K_S = 1;
//air drag
const float K_A = 0.01;
const float K_D = 0.1;
const vec3 g = vec3(0, -980, 0);
const float EPS = 1;
struct Particle {
	vec3 x;
	float m;
	vec3 v;
	vec3 f;
	float r=1;
	bool fixed = false;
	Particle(const vec3& p, float mass = 0.001, const vec3& vel = vec3(0))
		: x(p), m(mass), v(vel), f(0) {}

	void clearForce() { f = vec3(0); }
	void addForce(const vec3& force) { f += force; }
	void integrate(float dt)
	{
		if (fixed) v = vec3(0);

		else
		{
			x = x + v * dt;
			v = v + f / m * dt;

		}
	}

	void draw() const
	{
		drawSphere(x, 1, vec4(1, .7, 0, 1));
	}
};
struct Spring
{
	Particle& a;
	Particle& b;
	float k_s;
	float r;
	
	Spring(Particle& one, Particle& two, float k = DEF_K_S)
		: a(one), b(two), k_s(k), r(length(a.x - b.x)) {}

	void applyForce()
	{
		vec3 deltaX = a.x - b.x;
		vec3 direction = normalize(deltaX);
		vec3 deltaV = a.v - b.v;
		vec3 f = -(k_s * (length(deltaX) - r) + K_D * dot(deltaV, direction)) * direction;
		a.addForce(f);
		b.addForce(-f);

		
	}

	void draw()
	{
		drawCylinder(a.x, b.x, 0.4, vec4(1, .2, 0, 1));
	}

	
};

struct Plain
{
	vec3 p;
	vec3 n;

	void draw()
	{
		drawQuad(p,n, vec2(1000, 1000));

	}

	bool colliding(const Particle& particle)
	{
		if (dot(particle.x-p,n) < EPS && dot(n,particle.v) < EPS)
		{
			return true;
		}

		return false;
	}

	void resolveCollision(Particle& particle, float alpha = 0.8)
	{
		vec3 vn = dot(n,particle.v) * n;
		vec3 vt = particle.v - vn;
		particle.v = vt - alpha * vn;
	}
};

struct Sphere
{
	vec3 p= vec3(0,30,0);
	vec3 n;
	float r = 30.0;
	
	vec3 c = vec3(p.x, p.y/2 + r, p.z);
	void draw()
	{
		drawSphere(p, r, vec4(0.6417,0,0,1.0));

	}

	bool colliding(const Particle& particle)
	{
		
		if (length(particle.x-p) <= (particle.r+r))
		{
			//cout << length(particle.x - p) << "\n";
			return true;
		}

		return false;
	}

	void resolveCollision(Particle& particle, float alpha = 0.8)
	{
		n = particle.x - p;
		n = normalize(n);
		vec3 vn = dot(n, particle.v) * n;
		vec3 vt = particle.v - vn;
		particle.v = vt - alpha * vn;
	}
};

vector<Plain> plains;
vector<Particle> particles;
vector<Spring> springs; 
Sphere sphere;

float randf()
{
	return rand() / (float)RAND_MAX;
}
void init()
{
	particles.clear();
	springs.clear();
	plains.clear();
	plains.push_back({ vec3(0),vec3(0,1,0) });
	sphere = Sphere();
	
	
	const int N = 8;
	const int unitLength = 6;

	
	for (int i = 0; i < N; i++) for (int j = 0; j < N; j++)
	{
		//particles.push_back(Particle({ j*5-17.5, i*5+30,0 }));
		float rx = (randf() * 2 - 1) * 0.5f;
		float ry = (randf() * 2 - 1) * 0.5f;
		float rz = (randf() * 2 - 1) * 0.5f;
		particles.push_back(Particle({ (j - (unitLength - 1) / 2.f) *unitLength+rx, i * unitLength+ 80+ry,rz }));
	}

	for (int i = 0; i < N; i++) for (int j = 0; j < N-1; j++)
	{
		
		springs.push_back({ particles[j + i * N], particles[j + 1 + i * N] });
	}

	for (int i = 0; i < N-1; i++) for (int j = 0; j < N; j++)
	{

		springs.push_back({ particles[j + i * N], particles[j + (1 + i) * N] });
	}

	for (int i = 0; i < N - 1; i++) for (int j = 0; j < N - 1; j++)
	{
		springs.push_back({particles[j + i * N], particles[j+ (i+1)*N+1]});
	
	}
	for (int i = 0; i < N-1; i++) for (int j = 1; j < N; j++)
	{
		springs.push_back({ particles[j + i * N], particles[j +(i+1)*N-1] });
	}

	for (int i = 2; i < N - 2; i++) for (int j = 2; j < N-2; j++)
	{
		//upward
		springs.push_back({ particles[j + i * N], particles[j + (i - 2) * N] });
		//downward
		springs.push_back({ particles[j + i * N], particles[j + (i + 2) * N] });
		//right
		springs.push_back({ particles[j + i * N], particles[j + i  * N+2] });
		//left
		springs.push_back({ particles[j + i * N], particles[j + i * N-2] });

	}
	/*
	particles.push_back(Particle({ 0,50,0 }, 0.001, { 50,0,0 }));
	particles.push_back(Particle({ 10,50,0 }));
	
	springs.push_back(Spring(particles[0], particles[1]));
	*/
}

	
	


float t;
void frame(float dt)
{
	int step = 1000;
	for (int i = 0; i < step; i++)
	{
		for (auto& p : particles) p.clearForce();
		for (auto& p : particles) p.addForce(-K_A* p.v);
		for (auto& p : particles) p.addForce(p.m * g);
	
		for (auto& s : springs) s.applyForce();
		for(auto&p : particles)
			for (auto& plain : plains)
			{
				if (plain.colliding(p))
					plain.resolveCollision(p);
			}

		for (auto& p : particles)
		{
			if (sphere.colliding(p))
				sphere.resolveCollision(p);
		}

			

		for (auto& p : particles) p.integrate(dt/step );
	}
}

void render()
{
	for (auto& p : particles) p.draw();
	for (auto& s : springs) s.draw();
	for (auto& s : plains) s.draw();
	sphere.draw();
}

bool key(int k)
{
	if (k == '1') particles[0].fixed = !particles[0].fixed;
	return false;

	if (k == '0') init();
}

bool press( const glm::vec2& pt2, const glm::vec3& pt3, float d ) {

	return false;
}

bool drag( const glm::vec2& pt2, const glm::vec3& pt3, float d ) {
	
	return false;
}

bool release( const glm::vec2& pt2, const glm::vec3& pt3, float d ) {
	
	return true;
}

int main(int argc, const char * argv[]) {

	
	JGL::Window* window = new JGL::Window(800, 600, "simulation");
	window->alignment(JGL::ALIGN_ALL);
	animView = new AnimView(0, 0, 800, 600);
	animView->renderFunction = render;
	animView->frameFunction = frame;
	animView->initFunction = init;
	animView->pressFunc = press;
	animView->dragFunc = drag;
	animView->releaseFunc = release;
	animView->keyFunc = key;
	init();
	window->show();
	JGL::_JGL::run();
	return 0;
}


