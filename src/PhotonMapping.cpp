/*********************************************************************************
Copyright (C) 2014 Adrian Jarabo (ajarabo@unizar.es)
Copyright (C) 2014 Diego Gutierrez (diegog@unizar.es)
All rights reserved.

This is an educational Ray Tracer developed for the course 'Informatica Grafica'
(Computer Graphics) tought at Universidad de Zaragoza (Spain). As such, it does not
intend to be fast or general, but just to provide an educational tool for undergraduate
students.

This software is provided as is, and any express or implied warranties are disclaimed.
In no event shall copyright holders be liable for any damage.
**********************************************************************************/
#include "PhotonMapping.h"
#include "World.h"
#include "Intersection.h"
#include "Ray.h"
#include "BSDF.h"
#include <random>
#include <iostream>

#define M_PI           3.14159265358979323846  /* pi */

//*********************************************************************
// Compute the photons by tracing the Ray 'r' from the light source
// through the scene, and by storing the intersections with matter
// in the lists 'xx_photons', storing the diffuse (global) and caustic
// photons respectively. For efficiency, both are computed at the same
// time, since computing them separately would result into a lost of
// several samples marked as caustic or diffuse.
// Same goes with the boolean 'direct', that specifies if direct 
// photons (from light to surface) are being stored or not. 
// The initial traced photon has energy defined by the tristimulus
// 'p', that accounts for the emitted power of the light source.
// The function will return true when there are more photons (caustic
// or diffuse) to be shot, and false otherwise.
//---------------------------------------------------------------------
bool PhotonMapping::trace_ray(const Ray& r, const Vector3 &p,
	std::list<Photon> &global_photons, std::list<Photon> &caustic_photons, bool direct)
{

	//Check if max number of shots done...
	if (++m_nb_current_shots > m_max_nb_shots)
	{
		return false;
	}

	// Compute irradiance photon's energy
	Vector3 energy(p);

	Ray photon_ray(r);
	photon_ray.shift();

	bool is_caustic_particle = false;

	//Iterate the path
	while (1)
	{
		// Throw ray and update current_it
		Intersection it;
		world->first_intersection(photon_ray, it);

		if (!it.did_hit())
			break;

		//Check if has hit a delta material...
		if (it.intersected()->material()->is_delta())
		{
			// If delta material, then is caustic...
			// Don't store the photon!
			is_caustic_particle = true;
		}
		else if (photon_ray.get_level() > 0 || direct)
		{
			//If non-delta material, store the photon!
			if (is_caustic_particle)
			{
				//If caustic particle, store in caustics
				if (caustic_photons.size() < m_nb_caustic_photons)
					caustic_photons.push_back(Photon(it.get_position(), photon_ray.get_direction(), energy));
			}
			else
			{
				//If non-caustic particle, store in global
				if (global_photons.size() < m_nb_global_photons)
					global_photons.push_back(Photon(it.get_position(), photon_ray.get_direction(), energy));
			}
			is_caustic_particle = false;
		}

		Real pdf;

		Vector3 surf_albedo = it.intersected()->material()->get_albedo(it);
		Real avg_surf_albedo = surf_albedo.avg();

		Real epsilon2 = static_cast<Real>(rand()) / static_cast<Real>(RAND_MAX);
		while (epsilon2 < 0.)
			epsilon2 = static_cast<Real>(rand()) / static_cast<Real>(RAND_MAX);

		if (epsilon2 > avg_surf_albedo || photon_ray.get_level() > 20)
			break;

		// Random walk's next step
		// Get sampled direction plus pdf, and update attenuation
		it.intersected()->material()->get_outgoing_sample_ray(it, photon_ray, pdf);

		// Shade...
		energy = energy*surf_albedo;
		if (!it.intersected()->material()->is_delta())
			energy *= dot_abs(it.get_normal(), photon_ray.get_direction()) / 3.14159;

		energy = energy / (pdf*avg_surf_albedo);
	}

	if (caustic_photons.size() == m_nb_caustic_photons &&
		global_photons.size() == m_nb_global_photons)
	{
		m_max_nb_shots = m_nb_current_shots - 1;
		return false;
	}

	return true;
}

//*********************************************************************
// TODO: Implement the preprocess step of photon mapping,
// where the photons are traced through the scene. To do it,
// you need to follow these steps for each shoot:
//  1 - Sample a world's light source in the scene to create
//		the initial direct photon from the light source.
//	2 - Trace the photon through the scene storing the inter-
//		sections between the photons and matter. You can use
//		the function 'trace_ray' for this purpose.
//	3 - Finally, once all the photons have been shot, you'll
//		need to build the photon maps, that will be used later
//		for rendering. 
//		NOTE: Careful with function
//---------------------------------------------------------------------
void PhotonMapping::preprocess()
{
	std::vector<LightSource*> lights = world->light_source_list;

	Vector3 Lpos = world->light(0).get_position();
	cout << "\n";
	cout << Lpos.getComponent(0);
	cout << Lpos.getComponent(1);
	cout << Lpos.getComponent(2);
	cout << "\n";
	//Vector3 Ldir = world->light(0).get_incoming_direction;
	//Vector3 Lint = world->light(0).get_incoming_light();
	Vector3 Lint = world->light(0).get_intensities();

	bool moreShots = true;

	double_t radius = 1;

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dis(-radius, radius);

	std::list<Photon> global_photons;
	std::list<Photon> caustic_photons;

	int level = 0;

	while (moreShots){
		Real randomX = Lpos.getComponent(0) + dis(gen);
		Real randomY = Lpos.getComponent(1) + dis(gen);
		Real randomZ = Lpos.getComponent(2) + dis(gen);



		Vector3 pos(randomX, randomY, randomZ); //Punto aleatorio dentro del cubo
		Vector3 dir(Lpos.getComponent(0) - randomX, Lpos.getComponent(1) - randomY, Lpos.getComponent(2) - randomZ);

		if (randomX*randomX + randomY*randomY + randomZ*randomZ <= 1){ //Mira si está dentro de la esfera
			Ray r(pos, dir, level);
			moreShots = trace_ray(r, Lint, global_photons, caustic_photons, false);
		}
	}

	//Guardar en KDTREE
	for (std::list<Photon>::iterator it = global_photons.begin(); it != global_photons.end(); ++it){
		Photon photon = *it;
		m_global_map.store(std::vector<Real>(photon.position.data, photon.position.data + 3), photon);
	}
	cout << "\n";
	cout << global_photons.size();
	cout << "\n";

	for (std::list<Photon>::iterator it = caustic_photons.begin(); it != caustic_photons.end(); ++it){
		Photon photon = *it;
		m_caustics_map.store(std::vector<Real>(photon.position.data, photon.position.data + 3), photon);
	}

	//Balanceo de los arboles
	m_global_map.balance();

	if (!m_caustics_map.is_empty())
		m_caustics_map.balance();

	cout << caustic_photons.size();
	cout << "\n";

}

bool PhotonMapping::insideSphere(Vector3 &sphere, Vector3 &point, double_t &r){

	double_t xdiff = pow(point.getComponent(0) - sphere.getComponent(0), 2);
	double_t ydiff = pow(point.getComponent(1) - sphere.getComponent(1), 2);
	double_t zdiff = pow(point.getComponent(2) - sphere.getComponent(2), 2);
	if ((xdiff + ydiff + zdiff) <= pow(r, 2)){
		return true;
	}
	else{
		return false;
	}
}

//*********************************************************************
// TODO: Implement the function that computes the rendering equation 
// using radiance estimation with photon mapping, using the photon
// maps computed as a proprocess. Note that you will need to handle
// both direct and global illumination, together with recursive the 
// recursive evaluation of delta materials. For an optimal implemen-
// tation you should be able to do it iteratively.
// In principle, the class is prepared to perform radiance estimation
// using k-nearest neighbors ('m_nb_photons') to define the bandwidth
// of the kernel.
//---------------------------------------------------------------------
Vector3 PhotonMapping::shade(Intersection &it0)const
{
	Vector3 L(0);
	Intersection it(it0);

	//Codigo ejemplo busqueda mas cercanos

	//m_global_map.balance();

	Vector3 normal = it.get_normal();
	Real ka = 0.2;
	Real kd = 0.6;
	Real ks = 0.4;
	int n = 50;
	Vector3 albedo = it.intersected()->material()->get_albedo(it);

	//Ambient light
	Vector3 ambientL = world->get_ambient();
	L += ka*albedo*ambientL;

	//Direct light
	if (world->light(0).is_visible(it.get_position())){
		Vector3 directInt = world->light(0).get_incoming_light(it.get_position());
		Vector3 dd = world->light(0).get_incoming_direction(it.get_position());
		Vector3 directDir = dd.operator*(-1);
		Vector3 dV = it.get_ray().get_direction();
		Vector3 directV = dV.operator*(-1);
		Real lambert = 0;
		//Diffuse 
		if (kd > 0){
			lambert = normal.dot(directDir);
			if (lambert > 0){
				L += kd*lambert*directInt*albedo;
			}
		}

		//Specular
		if (ks > 0){
			Real twice = 2 * lambert;
			Vector3 aux = normal.operator*(twice);
			Vector3 LR = aux.operator-(directDir);

			Real dot = directV.dot(LR);
			if (dot > 0){
				Real spec = pow(dot, n);
				L += ks * spec*directInt*albedo;
			}
		}
	}

	//Indirect light

	//Coger n fotones más cercanos
	Vector3 p = it.get_position();
	std::vector<const KDTree<Photon, 3>::Node*> photons;
	Real max_distance = 0; //Radio de la circunferencia
	m_global_map.find(std::vector<Real>(p.data, p.data + 3), m_nb_photons, photons, max_distance);
	double A = M_PI*max_distance*max_distance;
	Vector3 irradiance(0, 0, 0);
	for (const KDTree<Photon, 3>::Node* n : photons){
		Photon ph = n->data();
		double dot = normal.dot(ph.direction);
		if (dot < 0.0) continue;
		irradiance += ph.flux * albedo;
	}
	irradiance = irradiance / 256;
	L+= irradiance * (1 / A) ;

	//**********************************************************************
	// The following piece of code is included here for two reasons: first
	// it works as a 'hello world' code to check that everthing compiles 
	// just fine, and second, to illustrate some of the functions that you 
	// will need when doing the work. Goes without saying: remove the 
	// pieces of code that you won't be using.
	//
	unsigned int debug_mode = 0;

	switch (debug_mode)
	{
	case 1:
		// ----------------------------------------------------------------
		// Display Albedo Only
		L = it.intersected()->material()->get_albedo(it);
		break;
	case 2:
		// ----------------------------------------------------------------
		// Display Normal Buffer
		L = it.get_normal();
		break;
	case 3:
		// ----------------------------------------------------------------
		// Display whether the material is specular (or refractive) 
		L = Vector3(it.intersected()->material()->is_delta());

		//Aqui no tiene sentido el fotonMapping
		break;

	case 4:
		// ----------------------------------------------------------------
		// Display incoming illumination from light(0)
		L = world->light(0).get_incoming_light(it.get_position());
		break;

	case 5:
		// ----------------------------------------------------------------
		// Display incoming direction from light(0)
		L = world->light(0).get_incoming_direction(it.get_position());
		break;

	case 6:
		// ----------------------------------------------------------------
		// Check Visibility from light(0)
		if (world->light(0).is_visible(it.get_position()))
			L = Vector3(1.);
		break;
	}
	// End of exampled code
	//**********************************************************************

	return L;
}