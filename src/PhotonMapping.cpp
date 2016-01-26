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
#include <math.h>

#define M_PI           3.14159265358979323846  /* pi */
#define M_MAX_DEPTH	   5
#define M_MIN_PHOTON   8

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
	bool caustics = false;
	bool globals = false;
	for (LightSource* l : world->light_source_list){
		Vector3 Lpos = l->get_position();
		cout << "\n";
		cout << Lpos.getComponent(0);
		cout << Lpos.getComponent(1);
		cout << Lpos.getComponent(2);
		cout << "\n";
		//Vector3 Ldir = world->light(0).get_incoming_direction;
		//Vector3 Lint = world->light(0).get_incoming_light();
		Vector3 Lint = l->get_intensities();

		bool moreShots = true;

		double_t radius = 1;

		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_real_distribution<Real> dis(-1, 1);
		std::uniform_real_distribution<Real> dis2(0, 2 * M_PI);

		std::list<Photon> global_photons;
		std::list<Photon> caustic_photons;

		int level = 3;

		while (moreShots){

			/*Real randomZ = dis(gen);
			Real r = sqrt(1 - (randomZ*randomZ));
			Real angle = dis2(gen);
			Real randomX = r*cos(angle);
			Real randomY = r*sin(angle);*/


			Real randomX = dis(gen);
			Real randomY = dis(gen);
			Real randomZ = dis(gen);

			Vector3 dir(randomX, randomY, randomZ);

			if (randomX*randomX + randomY*randomY + randomZ*randomZ <= 1){
				Ray ry(Lpos, dir, 0);
				moreShots = trace_ray(ry, Lint, global_photons, caustic_photons, false);
			}
		}

		if (global_photons.size() > 0)
			globals = true;
		if (caustic_photons.size() > 0)
			caustics = true;

		//Guardar en KDTREE
		for (std::list<Photon>::iterator it = global_photons.begin(); it != global_photons.end(); ++it){
			Photon photon = *it;
			m_global_map.store(std::vector<Real>(photon.position.data, photon.position.data + 3), photon);
		}

		for (std::list<Photon>::iterator it = caustic_photons.begin(); it != caustic_photons.end(); ++it){
			Photon photon = *it;
			m_caustics_map.store(std::vector<Real>(photon.position.data, photon.position.data + 3), photon);
		}
	}

	//Balanceo de los arboles
	if (globals)
		m_global_map.balance();

	if (caustics)
		m_caustics_map.balance();
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
Vector3 PhotonMapping::irradianceEstimate(KDTree<Photon, 3> tree, Intersection it, float k,
	Vector3 normal, Vector3 p, Vector3 albedo)const{
	float max_distance, A;
	Vector3 irradiance(0, 0, 0);
	std::vector<const KDTree<Photon, 3>::Node*> list;

	tree.find(std::vector<Real>(p.data, p.data + 3), m_nb_photons, list, max_distance);
	A = M_PI*max_distance*max_distance;

	for (const KDTree<Photon, 3>::Node* n : list){
		Photon ph = n->data();
		double dot = normal.dot(ph.direction);
		if (dot >= 0.0) continue;
		Vector3 vd = (ph.position - it.get_position());
		float wpc = (1 - (vd.length() / (k*max_distance)));
		irradiance += ph.flux * albedo * wpc;
	}
	//Cone filtering
	return irradiance / ((1 - (2 / 3 * k))*A);
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

	Vector3 p = it.get_position();
	Vector3 normal = it.get_normal();

	int n = 50;
	Vector3 albedo = it.intersected()->material()->get_albedo(it);

	Vector3 dV = it.get_ray().get_direction();  //Direccion rayo del punto al ojo
	Vector3 directV = dV.operator*(-1);			//Direccion rayo del punto al ojo corregido
	Vector3 ambientL = world->get_ambient();

	for (LightSource* l : world->light_source_list){

		//Ambient light
		L += albedo*ambientL;

		//Vectors V, L
		Vector3 directInt = l->get_incoming_light(it.get_position()); //Intensidad luz
		Vector3 dd = l->get_incoming_direction(it.get_position());	  //Direccion rayo desde luz
		Vector3 directDir = dd.operator*(-1);						  //Direccion rayo desde luz corregido

		//Direct light
		if (l->is_visible(it.get_position())){

			Real lambert = 0;
			//Diffuse 
			lambert = normal.dot(directDir);
			if (lambert > 0){
				L += lambert*directInt*albedo;
			}

			//Specular
			Real twice = 2 * lambert;
			Vector3 aux = normal.operator*(twice);
			Vector3 LR = aux.operator-(directDir);
			Real dot = directV.dot(LR);
			if (dot > 0){
				Real spec = pow(dot, n);
				L += spec*directInt*albedo;
			}
		}
	}
	//TODO: Revisar y terminar esto
	if (it.intersected()->material()->is_delta()){
		Ray ray;
		Real ignore(0);
		Intersection itdelta;
		it.intersected()->material()->get_outgoing_sample_ray(it, ray, ignore);
		ray.shift();
		world->first_intersection(ray, itdelta);
		int i = 0;
		while (itdelta.did_hit() && itdelta.intersected()->material()->is_delta() && i < M_MAX_DEPTH){
			itdelta.intersected()->material()->get_outgoing_sample_ray(itdelta, ray, ignore);
			world->first_intersection(ray, itdelta);
			i++;
		}

		if (itdelta.did_hit() && i < M_MAX_DEPTH){
			Vector3 color = shade(itdelta);
			L += color;
		}
	}
	else{
		//Es opaco --> Estimacion de radiancia
		//Coger n fotones más cercanos
		float k(1.0); //Constante de filtrado
		L += irradianceEstimate(m_caustics_map, it, k, normal, p, albedo);
		L += irradianceEstimate(m_global_map, it, k, normal, p, albedo);
	}
	return L;
}

